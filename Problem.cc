#include "Problem.h"
#include "Subdomain.h"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <ctime>

using std::vector;

Subdomain subdomain;

Problem::Problem()
{ }

Problem::~Problem()
{ }

void Problem::BuildProblem(Problem_Input* input) 
{
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	// Getting spatial Data from the input
	pin_x = input->pin_x;
	pin_y = input->pin_y;
	refinement = input->refinement;
	z_planes = input->z_planes;
	sp_disc = input->sp_disc;
	bcs = input->bcs;

	// Getting Angular data from the input
	num_polar = input->num_polar;
	num_azim = input->num_azim;
	ang_agg_type = input->ang_agg_type;

	// Build the Quadrature
	quad.BuildQuadrature(input);

	// Compute the Energy Group parameters
	num_groupsets = input->num_groupsets;
	group_per_groupset = input->num_groups/num_groupsets;

	// Getting Partition parameters from the input
	num_SML = input->num_SML;
	num_TpSML = input->num_TpSML;
	partition_type = input->partition_type;
	partition_function = input->partition_function;
	overload = input->overload;
	num_cellsets = input->num_cellsets;
	sched_type = input->sched_type;

	// Build the subdomain 
	subdomain.BuildSubdomain(rank, this);

	// Build the task vector
	int num_tasks = 0;
	// The size is # Anglesets * # Groupsets * # Cellsets on this SML
	num_tasks = quad.num_angleset*num_groupsets*subdomain.total_overload;
	All_Tasks.resize(num_tasks);
	//std::cout << "All_Tasks.size() = " << All_Tasks.size() << std::endl;
	int task_num = 0;
	int octant = -1;
	for (int i = 0; i < subdomain.total_overload; i++)
	{
		for (int j = 0; j < quad.num_angleset; j++)
		{
			for (int k = 0; k < num_groupsets; k++)
			{
				All_Tasks[task_num].cellset_id = subdomain.CellSetIDs[i];
				All_Tasks[task_num].cellset_id_loc = i;
				All_Tasks[task_num].angleset_id = j;
				All_Tasks[task_num].groupset_id = k;
				All_Tasks[task_num].omega = quad.GetOmega(j);
				octant = quad.Anglesets[j].octant;
				All_Tasks[task_num].octant = octant;
				All_Tasks[task_num].depth = subdomain.CellSets[i].GetDepth(subdomain.CellSetIDs[i], octant, this);
				task_num += 1;
			}
		}

	}

	// Now we order the All_Tasks vector by depth, then angleset, then groupset
	std::sort(All_Tasks.begin(), All_Tasks.end(), by_depth());

	// Pre-allocate the DFEM Matrices
	A_tilde.resize(4, std::vector<double>(4, 0.0));
	A.resize(4, std::vector<double>(4, 0.0));
	RHS.resize(4);
	bg.resize(4);

	M.resize(4);
	N.resize(6, vector<vector<Direction > >(4, vector< Direction >(4)));
	L.resize(4, vector<Direction >(4));

	// Resize the source. Piece-wise constant in each cell
	source.resize(4);
	source[0] = 1;

}


void Problem::Sweep()
{
	// Timers
	std::clock_t start = std::clock();
	std::clock_t start_task, start_solve, start_send, start_receive;
	std::clock_t start_cell, start_ang, start_group;
	std::clock_t start_incflux, start_A, start_phi;
	long double duration_task, duration_send, duration_receive, duration_solve;
	long double duration_cell, duration_ang, duration_group;
	long double duration_rhs, duration_incflux, duration_A, duration_phi;

	// PreAllocated incoming and outgoing face vectors;
	std::vector<int> incoming(3,0), outgoing(3,0);

	
//	Loop through task list
//	check to see if task is ready -- wait until it is (mpi poll?)
	std::vector<task>::iterator it = All_Tasks.begin();
	std::vector<task>::iterator it_end = All_Tasks.end();
	int task = 0;
	for (; it != it_end; it++, task++)
	{
		int target = 0;
		start_task = std::clock();
		//start_receive = clock();
		//duration_solve = 0;
		//duration_cell = 0;
		//duration_ang = 0;
		//duration_group = 0;
		//duration_rhs = 0;
		//duration_incflux = 0;
		//duration_A = 0;
		//duration_phi = 0;
		// Number of cells in each direction in this cellset
		int cells_x = subdomain.CellSets[(*it).cellset_id_loc].cells_x;
		int cells_y = subdomain.CellSets[(*it).cellset_id_loc].cells_y;
		int cells_z = subdomain.CellSets[(*it).cellset_id_loc].cells_z;

		int angle_per_angleset = quad.Anglesets[(*it).angleset_id].angle_per_angleset;

		//if (rank == 0){ std::cout << "AS_ID, angle_per_angleset: " << (*it).angleset_id << " " << angle_per_angleset << std::endl; }

		int octant = quad.Anglesets[(*it).angleset_id].octant;

		// Loop through all faces and designate them incoming or outgoing:
		int r = 0;
		int s = 0;
		for (int f = 0; f < 6; f++)
		{
			// Get the neighbors for each face
			Neighbor neighbor = subdomain.CellSets[(*it).cellset_id_loc].neighbors[f];
			// Again check for incoming and get the buffer matrices from faces
			if (dot((*it).omega, neighbor.direction) < 0)
			{
				incoming[r] = f;
				r += 1;
			}
			else
			{
				outgoing[s] = f;
				s += 1;
			}
		}

		// check to see if task has all required incident information
		// First we figure out which faces of the cellset are incoming:
		// Loop over cellset faces
		for (int f = 0; f < 3; f++)
		{
			// Get the neighbors for each face
			Neighbor neighbor = subdomain.CellSets[(*it).cellset_id_loc].neighbors[incoming[f]];

			// Get buffers from boundary conditions or neighbor SMLs
			if (neighbor.SML < 0 || neighbor.SML == rank)
			{
				subdomain.Set_buffer_from_bc(incoming[f]);
			}
			// These are blocking recieves. Since we know the task order of the sweep
			// We have the task wait until it has all its incident fluxes before
			// we continue with the sweep
			else
			{
				target = GetTarget((*it).angleset_id, (*it).groupset_id, (*it).cellset_id);
				if (incoming[f] == 0 || incoming[f] == 1)
				{
					int size = cells_y*cells_z*group_per_groupset*angle_per_angleset * 4;
					MPI_Status status;
					// buffer,size of buffer, data type, target, tag (face), comm
					MPI_Recv(&subdomain.X_buffer[0], size, MPI_DOUBLE, neighbor.SML, target, MPI_COMM_WORLD, &status);
				}
				if (incoming[f] == 2 || incoming[f] == 3)
				{
					int size = cells_x*cells_z*group_per_groupset*angle_per_angleset * 4;
					MPI_Status status;
					// buffer,size of buffer, data type, target, tag (face), comm
					MPI_Recv(&subdomain.Y_buffer[0], size, MPI_DOUBLE, neighbor.SML, target, MPI_COMM_WORLD, &status);
				}
				if (incoming[f] == 4 || incoming[f] == 5)
				{
					int size = cells_x*cells_y*group_per_groupset*angle_per_angleset * 4;
					MPI_Status status;
					// buffer,size of buffer, data type, target, tag (face), comm
					MPI_Recv(&subdomain.Z_buffer[0], size, MPI_DOUBLE, neighbor.SML, target, MPI_COMM_WORLD, &status);
				}
			}
		}
		//duration_receive = (std::clock() - start_task) / (double)CLOCKS_PER_SEC;
		// We march through the cells first in x, then y, then z
		// The order (left to right or right to left) depends on what 
		// octant the angleset is in so we call GetCell to figure out where we are
		
		for (int k = 0; k < cells_z; k++)
		{
			for (int j = 0; j < cells_y; j++)
			{
				for (int i = 0; i < cells_x; i++)
				{
					//start_cell = clock();
					int cell_id = GetCell(i, j, k, cells_x, cells_y, cells_z, octant);

					Cell &my_cell = subdomain.CellSets[(*it).cellset_id_loc].Cells[cell_id];
					std::vector<int> cell_ijk(3, 0);
					my_cell.GetCellijk(cell_id, cells_x, cells_y, cells_z, cell_ijk);

					// Get the cell's DFEM Matrices for building the A matrix
					M = subdomain.CellSets[(*it).cellset_id_loc].Cells[cell_id].M;
					N = subdomain.CellSets[(*it).cellset_id_loc].Cells[cell_id].N;
					L = subdomain.CellSets[(*it).cellset_id_loc].Cells[cell_id].L;

					// Since the source is piece-wise constant, we only need to update the average value
					source[0] = my_cell.GetSource();
					// Loop over angles in the angleset
					for (int m = 0; m < angle_per_angleset; m++)
					{
						//start_ang = clock();
						// Get the direction of the angle
						Direction omega = quad.Anglesets[(*it).angleset_id].Omegas[m];

						// Add in the gradient matrix to the A matrix
						for (int a = 0; a < 4; a++)
						{
							for (int b = 0; b < 4; b++)
							{
								A_tilde[a][b] = dot(omega, L[a][b]);
							}
							RHS[a] = M[a] * source[a];
						}
						// Loop over the incoming faces in this cell
						for (int f = 0; f < 3; f++)
						{
							for (int a = 0; a < 4; a++)
							{
								for (int b = 0; b < 4; b++)
								{
									A_tilde[a][b] += dot(-1 * omega, N[incoming[f]][a][b]);
								}
							}						
						} // faces

						// Loop over groups in this groupset
						for (int g = 0; g < group_per_groupset; g++)
						{
							//start_group = clock();
							// Initialize the b vector
							for (int a = 0; a < 4; a++)
								bg[a] = RHS[a];

							// Retrieve this cells sigma tot (this is in the group loop to simulate
							// multi-group cross sections
							double sigma_t = my_cell.GetSigmaTot();
							// Need to get incoming fluxes for the RHS
							//start_incflux = clock();
							for (int f = 0; f < 3; f++)
							{
								std::vector<double> temp(4, 0);
								subdomain.Get_buffer(cell_ijk[0], cell_ijk[1], cell_ijk[2], g, m, incoming[f], temp);
								//if (rank == 1){ std::cout << "line: " << __LINE__ << std::endl; }
								for (int a = 0; a < 4; a++)
								{
									for (int b = 0; b < 4; b++)
									{
										bg[a] += dot(-1 * omega, N[incoming[f]][a][b])*temp[b];
									}
								}
							}
							//duration_incflux = clock() - start_group;
							//start_A = clock();
							// Add the contribution to the A matrix and the RHS vector
							for (int a = 0; a < 4; a++)
							{
								for (int b = 0; b < 4; b++)
								{
									A[a][b] = A_tilde[a][b] + sigma_t * M[a];	
								}
							}
							//duration_A += clock() - start_group;
							// Solve A^-1*RHS and store it in RHS (4 = number of elements)
							//start_solve = clock();
							GE_no_pivoting(A, bg, 4);
							//duration_solve += (std::clock() - start_group);

							//start_phi = clock();
							// Now we accumulate the fluxes into phi;
							my_cell.phi[(*it).groupset_id*group_per_groupset + g] += bg[0] * quad.Anglesets[(*it).angleset_id].Weights[m];
							//std::cout << bg[0] << " " << quad.Anglesets[(*it).angleset_id].Weights[m] << std::endl;

							// Now we need to translate the cell average to the average on each 
							// face before we push to the down stream neighbers
							// This allows for direct data movement (no cell to cell mapping needed)
							// Again we need to loop over faces and check for incoming/outgoing
							double cell_average = RHS[0];
							for (int f = 0; f < 3; f++)
							{
								Direction facecenter = my_cell.facecenters[outgoing[f]];

								bg[0] = cell_average + facecenter.x*bg[1] + facecenter.y*bg[2] + facecenter.z*bg[3];
								// Now store outgoing fluxes in the buffer arrays
								// Outgoing x faces
								subdomain.Set_buffer(cell_ijk[0], cell_ijk[1], cell_ijk[2], g, m, outgoing[f], bg);
							}
							//duration_phi += clock() - start_group;
							//duration_group += std::clock() - start_group;
						} //groups
						//duration_ang += std::clock() - start_ang;
					} //angles
					//duration_cell += std::clock() - start_cell;
				} // cells in x
			} // cells in y
		} // cells in z
		
		//start_send = clock();

		for (int f = 0; f < 3; f++)
		{
			// figure out target face
			/*if (outgoing[f] == 0 || outgoing[f] == 2 || outgoing[f] == 4)
				target = outgoing[f] + 1;
			else
				target = outgoing[f] - 1;*/
			// Get the neighbors for each face
			Neighbor neighbor = subdomain.CellSets[(*it).cellset_id_loc].neighbors[outgoing[f]];

			if (neighbor.id < 0)
			{
				if (outgoing[f] == 0 || outgoing[f] == 1)
				{
					subdomain.CellSets[(*it).cellset_id_loc].SetBoundaryFlux(outgoing[f], (*it).angleset_id, (*it).groupset_id, subdomain.X_buffer);
				}
				if (outgoing[f] == 2 || outgoing[f] == 3)
				{
					subdomain.CellSets[(*it).cellset_id_loc].SetBoundaryFlux(outgoing[f], (*it).angleset_id, (*it).groupset_id, subdomain.Y_buffer);
				}
				if (outgoing[f] == 4 || outgoing[f] == 5)
				{
					subdomain.CellSets[(*it).cellset_id_loc].SetBoundaryFlux(outgoing[f], (*it).angleset_id, (*it).groupset_id, subdomain.Z_buffer);
				}
			}
			// send an mpi message to neighbor.SML
			else
			{
				target = GetTarget((*it).angleset_id, (*it).groupset_id, neighbor.id);
				if (outgoing[f] == 0 || outgoing[f] == 1)
				{
					int size = cells_y*cells_z*group_per_groupset*angle_per_angleset * 4;
					MPI_Request request;
					// buffer,size of buffer, data type, target, tag (face), comm
					MPI_Isend(&subdomain.X_buffer[0], size, MPI_DOUBLE, neighbor.SML, target, MPI_COMM_WORLD, &request);
					//MPI_Wait(&request, &status);
				}
				if (outgoing[f] == 2 || outgoing[f] == 3)
				{
					int size = cells_x*cells_z*group_per_groupset*angle_per_angleset * 4;
					MPI_Request request;
					// buffer,size of buffer, data type, target, tag (face), comm
					MPI_Isend(&subdomain.Y_buffer[0], size, MPI_DOUBLE, neighbor.SML, target, MPI_COMM_WORLD, &request);
					//MPI_Wait(&request, &status);
				}
				if (outgoing[f] == 4 || outgoing[f] == 5)
				{
					int size = cells_x*cells_y*group_per_groupset*angle_per_angleset * 4;
					MPI_Request request;
					// buffer,size of buffer, data type, target, tag (face), comm
					MPI_Isend(&subdomain.Z_buffer[0], size, MPI_DOUBLE, neighbor.SML, target, MPI_COMM_WORLD, &request);
					//MPI_Wait(&request, &status);

				}
			}

		}
		//duration_send = (std::clock() - start_send) / (double)CLOCKS_PER_SEC;
		duration_task = (std::clock() - start_task) / (double)CLOCKS_PER_SEC;
		
		//duration_cell = duration_cell - duration_ang;
		//duration_ang = duration_ang - duration_group;
		//duration_group = duration_group / (double)CLOCKS_PER_SEC;
		////duration_rhs = duration_rhs / (double)CLOCKS_PER_SEC;
		//duration_phi = (duration_phi - duration_solve) / (double)CLOCKS_PER_SEC;
		//duration_solve = (duration_solve - duration_A) / (double)CLOCKS_PER_SEC;
		//duration_A = duration_A / (double)CLOCKS_PER_SEC;
		//duration_incflux = duration_incflux / (double)CLOCKS_PER_SEC;
		//
		//
		//

		//duration_ang = duration_ang / (double)CLOCKS_PER_SEC;
		//duration_cell = duration_cell / (double)CLOCKS_PER_SEC;
		//if (rank == 0){ 
			//std::cout << "Rank = " << rank << " Task: " << task << " duration = " << duration_task << " seconds." << std::endl;
			//std::cout << "    Get info duration   = " << duration_receive << " seconds." << std::endl;
			//std::cout << "    cell loop duration  = " << duration_cell << " seconds." << std::endl;
			//std::cout << "    angle loop duration = " << duration_ang << " seconds." << std::endl;
			//std::cout << "    group loop duration = " << duration_group << " seconds." << std::endl;
			////std::cout << "        Get rhs duration    = " << duration_rhs << " seconds." << std::endl;
			////std::cout << "        Inc flux duration   = " << duration_incflux << " seconds." << std::endl;
			//std::cout << "        Build A duration    = " << duration_A << " seconds." << std::endl;
			//std::cout << "        solve duration      = " << duration_solve << " seconds." << std::endl;
			//std::cout << "        Get phi duration    = " << duration_phi << " seconds." << std::endl;
			//std::cout << "    send duration       = " << duration_send << " seconds." << std::endl;
		//}
	} // tasks
	long double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;

	//if (rank == 0){ 
		std::cout << "  Rank: " << rank << " Sweep duration: " << duration << " seconds." << std::endl;
	//}
	

	//std::cout << "Phi[0] = " << std::endl;
	//for (int c = 0; c < subdomain.total_overload; c++)
	//{
	//	// Number of cells in each direction in this cellset
	//	int cells_x = subdomain.CellSets[c].cells_x;
	//	int cells_y = subdomain.CellSets[c].cells_y;
	//	int cells_z = subdomain.CellSets[c].cells_z;
	//	for (int k = 0; k < cells_z; k++)
	//	{
	//		for (int j = 0; j < cells_y; j++)
	//		{
	//			for (int i = 0; i < cells_x; i++)
	//			{
	//				int cell_id = GetCell(i, j, k, cells_x, cells_y, cells_z, 0);
	//				std::cout << cell_id << " " << subdomain.CellSets[c].Cells[cell_id].phi[0] << std::endl;
	//			}
	//		}
	//	}
	//}
}

int Problem::GetCell(int i, int j, int k, int cells_x, int cells_y, int cells_z, int octant)
{
	int d = -1;

	if (octant == 0)
		d = i + cells_x*j + cells_x*cells_y*k;
	else if (octant == 1)
		d = (cells_x - i -1) + cells_x*j + cells_x*cells_y*k;
	else if (octant == 2)
		d = (cells_x - i -1) + cells_x*(cells_y - j -1) + cells_x*cells_y*k;
	else if (octant == 3)
		d = i + cells_x*(cells_y - j -1) + cells_x*cells_y*k;
	else if (octant == 4)
		d = i + cells_x*j + cells_x*cells_y*(cells_z - k - 1);
	else if (octant == 5)
		d = (cells_x - i -1) + cells_x*j + cells_x*cells_y*(cells_z - k -1);
	else if (octant == 6)
		d = (cells_x - i -1) + cells_x*(cells_y - j -1) + cells_x*cells_y*(cells_z - k -1);
	else if (octant == 7)
		d = i + cells_x*(cells_y - j -1) + cells_x*cells_y*(cells_z - k -1);

	return d;
}

void Problem::GE_no_pivoting(std::vector<std::vector< double > >& A, std::vector<double>& b, int n)
{
	// Forward elimination
	for (int i = 0; i < n - 1; ++i)
	{
		const std::vector<double>& ai = A[i];
		double bi = b[i];
		double factor = 1.0 / A[i][i];
		for (int j = i + 1; j < n; ++j)
		{
			std::vector<double>& aj = A[j];
			double val = aj[i] * factor;
			b[j] -= val * bi;
			for (int k = i + 1; k < n; ++k)
				aj[k] -= val * ai[k];
		}
	}

	// Back substitution
	for (int i = n - 1; i >= 0; --i)
	{
		const std::vector<double>& ai = A[i];
		double bi = b[i];
		for (int j = i + 1; j < n; ++j)
			bi -= ai[j] * b[j];
		b[i] = bi / ai[i];
	}
}

unsigned int Problem::GetTarget(int as, int gs, int cs)
{
	// The target will be an integer with 10 digits
	// The first 4 digits will be the angleset, the second 3 will be
	// the groupset, and the final 3 digits will be the neighbor cellset
	// This assumes the total number of anglesets is less than 10,000, 
	// and the number of groupsets and cellsets are both less than 1,000

	unsigned int target;
	target = cs + 1000 * gs + pow(1000,2) * as + pow(1000,3);
	//target = pow(10000, 3);
	//std::cout << "CS, GS, AS, target: " << cs << " " << gs << " " << as << " " << target << std::endl;
	return target;
}