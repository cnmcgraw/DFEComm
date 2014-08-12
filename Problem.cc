#include "Problem.h"
#include "Subdomain.h"
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <deque>

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
	// The size is # Anglesets * # Groupsets * # Cellsets on this SML
	num_tasks = quad.num_angleset*num_groupsets*subdomain.total_overload;
	All_Tasks.resize(num_tasks);
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

	// Allocate the Send Buffers
	subdomain.AllocateBuffers(num_tasks);

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

void Problem::Sweep(std::ofstream &output)
{
	// Timers
	double start_task;
	long double duration_task;

	// Necessary temporary data structures for the sweep
	std::vector<std::vector<int> > incoming(3, std::vector<int>(2,0));
	std::vector<int> outgoing(3,0);
	std::vector<int> cell_ijk(3, 0);
	std::vector<double> temp_solve(4, 0);
	Neighbor neighbor;
	int task_it(0), recv(0), target(0);

	// Request for sends and receives
	MPI_Request sendrequest[3];
	for (int i = 0; i < 3; i++)
		sendrequest[i] = MPI_REQUEST_NULL;

	MPI_Request Request[3 * num_tasks];
	for (int i = 0; i < 3 * num_tasks; i++)
		Request[i] = MPI_REQUEST_NULL;

    // Loop through task list
	std::vector<task>::iterator it = All_Tasks.begin();
	std::vector<task>::iterator it_end = All_Tasks.end();
	for (; it != it_end; it++, task_it++)
	{
		start_task = MPI_Wtime();

		// Get cellset and angleset information
		int cells_x = subdomain.CellSets[(*it).cellset_id_loc].cells_x;
		int cells_y = subdomain.CellSets[(*it).cellset_id_loc].cells_y;
		int cells_z = subdomain.CellSets[(*it).cellset_id_loc].cells_z;
		int angle_per_angleset = quad.Anglesets[(*it).angleset_id].angle_per_angleset;
		int octant = quad.Anglesets[(*it).angleset_id].octant;

		// Loop through all faces and designate them incoming or outgoing:
		int r(0), s(0);
		for (int f = 0; f < 6; f++)
		{
			// Get the neighbors for each face
			neighbor = subdomain.CellSets[(*it).cellset_id_loc].neighbors[f];
			if (dot((*it).omega, neighbor.direction) < 0)
			{
				incoming[r][0] = f;
				incoming[r][1] = neighbor.SML;
				r += 1;
			}
			else
			{
				outgoing[s] = f;
				s += 1;
			}
		}

		// Check to see if task has all required incident information
		bool ready = false;
		std::vector<bool> ready_face(3, false);
		int flag, count;
		MPI_Status status;

		// First probe for any messages waiting to be received
		MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
		// While there are messages waiting to be received, post receives and Probe again
		while (flag)
		{
			MPI_Get_count(&status, MPI_DOUBLE, &count);
			recv = GetPlacement();
			MPI_Recv(&subdomain.Received_buffer[(recv)*subdomain.max_size], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			subdomain.Received_info[recv][0] = status.MPI_TAG;
			subdomain.Received_info[recv][1] = status.MPI_SOURCE;
			subdomain.Received_info[recv][2] = count;

			MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
		}
		// Now look for incoming info from boundary buffers
		for (int f = 0; f < 3; f++)
		{
			if (incoming[f][1] < 0 || incoming[f][1] == rank)
			{
				subdomain.Get_buffer_from_bc(incoming[f][0]);
				ready_face[f] = true;
			}
		}
		ready = (ready_face[0] && ready_face[1] && ready_face[2]);

		// Check in Received buffer for info we already received
		// If we find a matching message, we put its location in the 
		// received buffer back in the queue of empty spots
		if (!ready)
		{
			target = GetTarget((*it).angleset_id, (*it).groupset_id, (*it).cellset_id);
			for (int inf = 0; inf < subdomain.Received_info.size(); inf++)
			{
				// Matching tags?
				if (subdomain.Received_info[inf][0] == target)
				{
					for (int f = 0; f < 3; f++)
					{
						// Matching SML?
						if (subdomain.Received_info[inf][1] == incoming[f][1])
						{
							if (incoming[f][0] == 0 || incoming[f][0] == 1)
							{
								std::vector<double>::iterator buf_it = subdomain.Received_buffer.begin() + inf*subdomain.max_size;
								subdomain.X_buffer.assign(buf_it, buf_it + subdomain.Received_info[inf][2]);
								subdomain.Received_open.push(inf);
								ready_face[f] = true;
							}
							if (incoming[f][0] == 2 || incoming[f][0] == 3)
							{
								std::vector<double>::iterator buf_it = subdomain.Received_buffer.begin() + inf*subdomain.max_size;
								subdomain.Y_buffer.assign(buf_it, buf_it + subdomain.Received_info[inf][2]);
								subdomain.Received_open.push(inf);
								ready_face[f] = true;
							}
							if (incoming[f][0] == 4 || incoming[f][0] == 5)
							{
								std::vector<double>::iterator buf_it = subdomain.Received_buffer.begin() + inf*subdomain.max_size;
								subdomain.Z_buffer.assign(buf_it, buf_it + subdomain.Received_info[inf][2]);
								subdomain.Received_open.push(inf);
								ready_face[f] = true;
							}
						}
					}
				}
			}
		}
		ready = (ready_face[0] && ready_face[1] && ready_face[2]);
		// Keep receiving messages until we have all our incoming info
		while (ready == false)
		{
			// Now we probe for messages
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
			if (flag == true)
			{
				// Check if the message is for this task
				if (status.MPI_TAG == target)
				{
					for (int f = 0; f < 3; f++)
					{
						if (status.MPI_SOURCE == incoming[f][1])
						{
							MPI_Get_count(&status, MPI_DOUBLE, &count);
							if (incoming[f][0] == 0 || incoming[f][0] == 1)
							{
								MPI_Recv(&subdomain.X_buffer[0], count, MPI_DOUBLE, incoming[f][1], target, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								ready_face[f] = true;
							}
							if (incoming[f][0] == 2 || incoming[f][0] == 3)
							{
								MPI_Recv(&subdomain.Y_buffer[0], count, MPI_DOUBLE, incoming[f][1], target, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								ready_face[f] = true;
							}
							if (incoming[f][0] == 4 || incoming[f][0] == 5)
							{
								MPI_Recv(&subdomain.Z_buffer[0], count, MPI_DOUBLE, incoming[f][1], target, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								ready_face[f] = true;
							}
						}
					}
				}
				// if we don't need this messsage yet, store it in Received_buffer and store the info about it in Received_info
				else
				{
					MPI_Get_count(&status, MPI_DOUBLE, &count);
					recv = GetPlacement();
					MPI_Recv(&subdomain.Received_buffer[(recv)*subdomain.max_size], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					subdomain.Received_info[recv][0] = status.MPI_TAG;
					subdomain.Received_info[recv][1] = status.MPI_SOURCE;
					subdomain.Received_info[recv][2] = count;
				}
			}
			ready = (ready_face[0] && ready_face[1] && ready_face[2]);
		}

		// Wait until previously posted Isends have completed, 
		// now that all the Receives have been posted
		MPI_Waitall(3, sendrequest, MPI_STATUSES_IGNORE);


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
						// Get the direction of the angle
						omega = quad.Anglesets[(*it).angleset_id].Omegas[m];

						// Add in the gradient matrix to the A matrix
						for (int a = 0; a < 4; a++)
						{
							for (int b = 0; b < 4; b++)
							{
								A_tilde[a][b] = dot(omega, L[a][b]);
							}
							RHS[a] = M[a] * source[a];
						}
						// Loop over the incoming faces in this cell and 
						// add the surface matrix contribution
						for (int f = 0; f < 3; f++)
						{
							for (int a = 0; a < 4; a++)
							{
								for (int b = 0; b < 4; b++)
								{
									A_tilde[a][b] += dot(-1 * omega, N[incoming[f][0]][a][b]);
								}
							}						
						}

						// Loop over groups in this groupset
						for (int g = 0; g < group_per_groupset; g++)
						{
							// Initialize the b vector
							for (int a = 0; a < 4; a++)
								bg[a] = RHS[a];

							// Retrieve this cells sigma tot (this is in the group loop to simulate
							// multi-group cross sections
							sigma_t = my_cell.GetSigmaTot();

							// Need to get incoming fluxes for the RHS
							for (int f = 0; f < 3; f++)
							{
								subdomain.Get_buffer(cell_ijk[0], cell_ijk[1], cell_ijk[2], g, m, incoming[f][0], temp_solve);
								for (int a = 0; a < 4; a++)
								{
									for (int b = 0; b < 4; b++)
									{
										bg[a] += dot(-1 * omega, N[incoming[f][0]][a][b])*temp_solve[b];
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

							// Solve A^-1*RHS wwith Gaussian Elimination and store it in RHS (4 = number of elements)
							GE_no_pivoting(A, bg, 4);

							// Now we accumulate the fluxes into phi;
							for (int p = 0; p < 4; p++)
								my_cell.phi[(*it).groupset_id*group_per_groupset * 4 + 4 * g + p] += bg[p] * quad.Anglesets[(*it).angleset_id].Weights[m];

							// Now we need to translate the cell average to the average on each 
							// face before we push to the down stream neighbers
							// This allows for direct data movement (no cell to cell mapping needed)
							cell_average = bg[0];
							for (int f = 0; f < 3; f++)
							{
								facecenter = my_cell.facecenters[outgoing[f]];

								bg[0] = cell_average + facecenter.x*bg[1] + facecenter.y*bg[2] + facecenter.z*bg[3];
								// Now store outgoing fluxes in the buffer arrays
								subdomain.Set_buffer(cell_ijk[0], cell_ijk[1], cell_ijk[2], g, m, outgoing[f], task_it, bg);
							}
						} //groups
					} //angles
				} // cells in x
			} // cells in y
		} // cells in z

		for (int f = 0; f < 3; f++)
		{
			// Get the neighbors for each face
			neighbor = subdomain.CellSets[(*it).cellset_id_loc].neighbors[outgoing[f]];

			// If we are on a global boundary, store in a boundary buffer
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
			// Otherwise, send an mpi message to neighbor.SML
			else
			{
				target = GetTarget((*it).angleset_id, (*it).groupset_id, neighbor.id);
				if (outgoing[f] == 0 || outgoing[f] == 1)
				{
					int size = cells_y*cells_z*group_per_groupset*angle_per_angleset * 4;
					MPI_Isend(&subdomain.X_Send_buffer[0], size, MPI_DOUBLE, neighbor.SML, target, MPI_COMM_WORLD, &sendrequest[0]);

				}
				if (outgoing[f] == 2 || outgoing[f] == 3)
				{
					int size = cells_x*cells_z*group_per_groupset*angle_per_angleset * 4;
					MPI_Isend(&subdomain.Y_Send_buffer[0], size, MPI_DOUBLE, neighbor.SML, target, MPI_COMM_WORLD, &sendrequest[1]);

				}
				if (outgoing[f] == 4 || outgoing[f] == 5)
				{
					int size = cells_x*cells_y*group_per_groupset*angle_per_angleset * 4;
					MPI_Isend(&subdomain.Z_Send_buffer[0], size, MPI_DOUBLE, neighbor.SML, target, MPI_COMM_WORLD, &sendrequest[2]);
				}
			}

		}
		duration_task = (MPI_Wtime() - start_task); 
		if (rank == 0){ 
			output << "    Task: " << task_it << " duration = " << duration_task << " seconds." << std::endl;
		}
	} // tasks
	////Printing out Phi
	//if (rank == 0)
	//{
	//	std::cout << "Phi[g = 0] = " << std::endl;
	//	for (int c = 0; c < subdomain.total_overload; c++)
	//	{
	//		// Number of cells in each direction in this cellset
	//		int cells_x = subdomain.CellSets[c].cells_x;
	//		int cells_y = subdomain.CellSets[c].cells_y;
	//		int cells_z = subdomain.CellSets[c].cells_z;
	//		for (int k = 0; k < cells_z; k++)
	//		{
	//			for (int j = 0; j < cells_y; j++)
	//			{
	//				for (int i = 0; i < cells_x; i++)
	//				{
	//					int cell_id = GetCell(i, j, k, cells_x, cells_y, cells_z, 0);
	//					std::cout << cell_id << " " << subdomain.CellSets[c].Cells[cell_id].phi[0] << std::endl;
	//				}
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
	return target;
}

int Problem::GetPlacement()
{
	// If the Receive buffer is full, we need to allocate more space and add its 
	// location to the queue
	if (subdomain.Received_open.empty())
	{
		std::cout << "The queue is empty" << std::endl;
		int size = subdomain.Received_buffer.size();
		int size_info = subdomain.Received_info.size();
		subdomain.Received_buffer.resize(size + subdomain.max_size);
		subdomain.Received_info.resize(size_info + 1, std::vector<int>(3, 0));
		subdomain.Received_open.push(size_info + 1);
	}

	// Now we pick out the first element of the queue
	int location = subdomain.Received_open.front();
	subdomain.Received_open.pop();
	return location;
}

void Problem::ZeroPhi()
{
	// Loop through cellsets in this SMLs subdomain
	for (int i = 0; i < subdomain.CellSets.size(); i++)
	{
		// Loop through cells in this cellset
		for (int j = 0; j < subdomain.CellSets[i].Cells.size(); j++)
		{
			// Loop through groups
			for (int k = 0; k < subdomain.CellSets[i].Cells[j].phi.size(); k++)
					subdomain.CellSets[i].Cells[j].phi[k] = 0;
		}

	}

	for (int i = 0; i < subdomain.Received_buffer.size(); i++)
		subdomain.Received_buffer[i] = 0;

	for (int i = 0; i < subdomain.Received_info.size(); i++)
	{
		for (int j = 0; j < subdomain.Received_info[i].size(); j++)
			subdomain.Received_info[i][j] = 0;
	}
}