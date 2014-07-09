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

	//std::cout << "Rank: " << rank << " and num_tasks: " << num_tasks << std::endl;

	// Now we order the All_Tasks vector by depth, then angleset, then groupset
	std::sort(All_Tasks.begin(), All_Tasks.end(), by_depth());

	//if (rank == 1)
	//{
	//	std::cout << "Task 1's task octants (in order)" << std::endl;
	//	for (int i = 0; i < num_tasks; i++)
	//		std::cout << All_Tasks[i].octant << " " << All_Tasks[i].depth << std::endl;

	//	std::cout << " " << std::endl;
	//}

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
	std::clock_t start = std::clock();
	std::clock_t start_task, start_cell, start_send, start_receive;
	std::clock_t start_solve, start_geo, start_trans, start_ang;
	long double duration_task, duration_get, duration_cell, duration_send, duration_receive, duration_solve;
	long double duration_geo, duration_trans, duration_ang;

	
//	Loop through task list
//	check to see if task is ready -- wait until it is (mpi poll?)
	std::vector<task>::iterator it = All_Tasks.begin();
	std::vector<task>::iterator it_end = All_Tasks.end();
	int i = 0;
	for (; it != it_end; it++, i++)
	{
		start_task = std::clock();
		duration_solve = 0;
		duration_geo = 0;
		duration_ang = 0;
		duration_trans = 0;
		// Number of cells in each direction in this cellset
		int cells_x = subdomain.CellSets[(*it).cellset_id_loc].cells_x;
		int cells_y = subdomain.CellSets[(*it).cellset_id_loc].cells_y;
		int cells_z = subdomain.CellSets[(*it).cellset_id_loc].cells_z;

		int angle_per_angleset = quad.Anglesets[(*it).angleset_id].angle_per_angleset;

		int octant = quad.Anglesets[(*it).angleset_id].octant;
		// check to see if task has all required incident information
		// First we figure out which faces of the cellset are incoming:
		// Loop over cellset faces
		for (int f = 0; f < 6; f++)
		{
			// Get the neighbors for each face
			Neighbor neighbor = subdomain.CellSets[(*it).cellset_id_loc].neighbors[f];
			// Again check for incoming and get the buffer matrices from faces

			if (dot((*it).omega, neighbor.direction) < 0)
			{
				if (neighbor.SML < 0 || neighbor.SML == rank)
				{
					subdomain.Set_buffer_from_bc(f);
				}
				// These are blocking recieves. Since we know the task order of the sweep
				// We have the task wait until it has all its incident fluxes before
				// we continue with the sweep
				else
				{
					if (f == 0 || f == 1)
					{
						int size = cells_y*cells_z*group_per_groupset*angle_per_angleset * 4;
						MPI_Status status;
						// buffer,size of buffer, data type, target, tag (face), comm
						MPI_Recv(&subdomain.X_buffer[0], size, MPI_DOUBLE, neighbor.SML, f, MPI_COMM_WORLD, &status);
					}
					if (f == 2 || f == 3)
					{
						int size = cells_x*cells_z*group_per_groupset*angle_per_angleset * 4;
						MPI_Status status;
						// buffer,size of buffer, data type, target, tag (face), comm
						MPI_Recv(&subdomain.Y_buffer[0], size, MPI_DOUBLE, neighbor.SML, f, MPI_COMM_WORLD, &status);
					}
					if (f == 4 || f == 5)
					{
						int size = cells_x*cells_y*group_per_groupset*angle_per_angleset * 4;
						start_receive = clock();
						MPI_Status status;
						// buffer,size of buffer, data type, target, tag (face), comm
						MPI_Recv(&subdomain.Z_buffer[0], size, MPI_DOUBLE, neighbor.SML, f, MPI_COMM_WORLD, &status);
						duration_receive = (std::clock() - start_task) / (double)CLOCKS_PER_SEC;
						//std::cout << "Rank = " << rank << ", face:   " << f << " and octant: " << octant << ". Receive duration = " << duration_receive << " seconds." << std::endl;
					}
				}
			}
		}

		duration_get = (std::clock() - start_task) / (double)CLOCKS_PER_SEC;
		//if (rank == 0){ std::cout << "    Task: " << i << " get face info duration = " << duration_get << " seconds." << std::endl; }

		// We march through the cells first in x, then y, then z
		// The order (left to right or right to left) depends on what 
		// octant the angleset is in so we call GetCell to figure out where we are
		duration_cell = 0;
		for (int k = 0; k < cells_z; k++)
		{
			start_cell = clock();
			for (int j = 0; j < cells_y; j++)
			{
				for (int i = 0; i < cells_x; i++)
				{
					start_geo = clock();
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

					duration_geo += clock() - start_geo;
					// Loop over angles in the angleset
					for (int m = 0; m < angle_per_angleset; m++)
					{
						start_ang = clock();
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
						// Loop over the faces in this cell
						for (int f = 0; f < 6; f++)
						{
							// Check to see if angle is incoming or out going. 
							// If its incoming, add the surface matrix's contribution to the A matrix
							Direction normal = my_cell.normals[f];
							if (dot(omega, normal) > 0)
							{
								for (int a = 0; a < 4; a++)
								{
									for (int b = 0; b < 4; b++)
									{
										A_tilde[a][b] += dot(-1 * omega, N[f][a][b]);
									}
								}
							}							
						} // faces

						duration_ang += clock() - start_ang;

						// Loop over groups in this groupset
						for (int g = 0; g < group_per_groupset; g++)
						{
							// Initialize the b vector
							for (int a = 0; a < 4; a++)
								bg[a] = RHS[a];

							// Retrieve this cells sigma tot (this is in the group loop to simulate
							// multi-group cross sections
							double sigma_t = my_cell.GetSigmaTot();
							// Need to get incoming fluxes for the RHS
							for (int f = 0; f < 6; f++)
							{
								Direction normal = my_cell.normals[f];
								if (dot(omega, normal) > 0)
								{
									std::vector<double> temp(4, 0);
									subdomain.Get_buffer(cell_ijk[0], cell_ijk[1], cell_ijk[2], g, m, f, temp);
									//if (rank == 1){ std::cout << "line: " << __LINE__ << std::endl; }
									for (int a = 0; a < 4; a++)
									{
										for (int b = 0; b < 4; b++)
										{
											bg[a] += dot(-1 * omega, N[f][a][b])*temp[b];
										}
									}
								}
							}
							// Add the contribution to the A matrix and the RHS vector
							for (int a = 0; a < 4; a++)
							{
								for (int b = 0; b < 4; b++)
								{
									A[a][b] = A_tilde[a][b] + sigma_t * M[a];	
								}
							}
							// Solve A^-1*RHS and store it in RHS (4 = number of elements)
							start_solve = clock();
							GE_no_pivoting(A, bg, 4);
							duration_solve += (std::clock() - start_solve);

							// Now we accumulate the fluxes into phi;
							my_cell.phi[(*it).groupset_id*group_per_groupset + g] += bg[0] * quad.Anglesets[(*it).angleset_id].Weights[m];
							//std::cout << bg[0] << " " << quad.Anglesets[(*it).angleset_id].Weights[m] << std::endl;

							// Now we need to translate the cell average to the average on each 
							// face before we push to the down stream neighbers
							// This allows for direct data movement (no cell to cell mapping needed)
							// Again we need to loop over faces and check for incoming/outgoing
							double cell_average = RHS[0];
							start_trans = clock();
							for (int f = 0; f < 6; f++)
							{
								Direction normal = my_cell.normals[f];
								Direction facecenter = my_cell.facecenters[f];
								if (dot(omega, normal) < 0)
								{
									bg[0] = cell_average + facecenter.x*bg[1] + facecenter.y*bg[2] + facecenter.z*bg[3];
									// Now store outgoing fluxes in the buffer arrays
									// Outgoing x faces
									subdomain.Set_buffer(cell_ijk[0], cell_ijk[1],cell_ijk[2], g, m, f, bg);
								}
							}
							duration_trans += clock() - start_trans;
						} //groups
					} //angles
				} // cells in x
			} // cells in y
			duration_cell += (std::clock() - start_cell) / (double)CLOCKS_PER_SEC;
		} // cells in z
		//if (rank == 0){ std::cout << "    Cell duration = " << duration_cell << " seconds." << std::endl; }

		int target = 0;
		
		for (int f = 0; f < 6; f++)
		{
			// figure out target face
			if (f == 0 || f == 2 || f == 4)
				target = f + 1;
			else
				target = f - 1;
			// Get the neighbors for each face
			Neighbor neighbor = subdomain.CellSets[(*it).cellset_id_loc].neighbors[f];
			// Again check for outgoing and send the buffer matrices to outgoing faces
			if (dot((*it).omega, neighbor.direction) > 0)
			{
				if (neighbor.SML < 0)
				{
					// store the buffer into the boundary information
				}
				// send an mpi message to neighbor.SML
				if (f == 0 || f == 1)
				{
					int size = cells_y*cells_z*group_per_groupset*angle_per_angleset * 4;
					MPI_Request request;
					// buffer,size of buffer, data type, target, tag (face), comm
					MPI_Isend(&subdomain.X_buffer[0], size, MPI_DOUBLE, neighbor.SML, target, MPI_COMM_WORLD, &request);
				}
				if (f == 2 || f == 3)
				{
					int size = cells_x*cells_z*group_per_groupset*angle_per_angleset * 4;
					MPI_Request request;
					// buffer,size of buffer, data type, target, tag (face), comm
					MPI_Isend(&subdomain.Y_buffer[0], size, MPI_DOUBLE, neighbor.SML, target, MPI_COMM_WORLD, &request);
				}
				if (f == 4 || f == 5)
				{
					int size = cells_x*cells_y*group_per_groupset*angle_per_angleset * 4;
					MPI_Request request;
					start_send = clock();
					// buffer,size of buffer, data type, target, tag (face), comm
					MPI_Isend(&subdomain.Z_buffer[0], size, MPI_DOUBLE, neighbor.SML, target, MPI_COMM_WORLD, &request);
					duration_send = (std::clock() - start_send) / (double)CLOCKS_PER_SEC;
					//std::cout << "Rank = " << rank << ", target: " << target << " and octant: " << octant << ". Send    duration = " << duration_send << " seconds." << std::endl;

				}
			}
		}
		
		duration_task = (std::clock() - start_task) / (double)CLOCKS_PER_SEC;
		duration_solve = duration_solve / (double)CLOCKS_PER_SEC;
		duration_ang = duration_ang / (double)CLOCKS_PER_SEC;
		duration_geo = duration_geo / (double)CLOCKS_PER_SEC;
		duration_trans = duration_trans / (double)CLOCKS_PER_SEC;
		if (rank == 0){ 
			std::cout << "Task: " << i << " duration = " << duration_task << " seconds." << std::endl;
			std::cout << "    Task: " << i << " receive duration = " << duration_get << " seconds." << std::endl;
			std::cout << "    Task: " << i << " angle duration   = " << duration_ang << " seconds." << std::endl;
			std::cout << "    Task: " << i << " solve duration   = " << duration_solve << " seconds." << std::endl;
			std::cout << "    Task: " << i << " geo duration     = " << duration_geo << " seconds." << std::endl;
			std::cout << "    Task: " << i << " trans duration   = " << duration_trans << " seconds." << std::endl;
			std::cout << "    Task: " << i << " send duration    = " << duration_send << " seconds." << std::endl;
		}
	} // tasks
	long double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;

	if (rank == 0){ std::cout << "Sweep duration: " << duration << " seconds." << std::endl;}
	

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
