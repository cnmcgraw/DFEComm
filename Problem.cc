#include "Problem.h"
#include "Subdomain.h"
#include "Task.h"
#include "ParallelComm.h"
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <deque>

using std::vector;

Subdomain subdomain;
ParallelComm comm;
std::vector<Task> All_Tasks;

// This struct compares tasks:
// First by depth (biggest wins),
// then by direction: first x, then y, then z (positive wins)
// and finally by groupset id (smallest wins)
struct by_depth {
  bool operator()(Task const &a, Task const &b) {
    if (a.depth == b.depth){
      if (a.angleset_id == b.angleset_id)
        return a.groupset_id < b.groupset_id;
      else
      {
        if (a.omega.x == b.omega.x){
          if (a.omega.y == b.omega.y){
            return a.omega.z > b.omega.z;
          }
          else
            return a.omega.y > b.omega.y;
        }
        else
          return a.omega.x > b.omega.x;
      }
    }
    else
      return a.depth > b.depth;
  }
};

Problem::Problem()
{ }

Problem::~Problem()
{ }

void Problem::BuildProblem(Problem_Input* input) 
{
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  // Getting spatial Data from the input
  num_pin_x = input->num_pin_x;
  num_pin_y = input->num_pin_y;
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
  partition_type = input->partition_type;
  partition_function = input->partition_function;
  overload = input->overload;
  num_cellsets = input->num_cellsets;
  sched_type = input->sched_type;
  verbose = input->verbose;

  // Build the subdomain 
  subdomain.BuildSubdomain(rank, this);

  // Build the task vector
  // The size is # Anglesets * # Groupsets * # Cellsets on this SML
  num_tasks = quad.num_angleset*num_groupsets*subdomain.total_overload;
 // All_Tasks.resize(num_tasks);
  int task_num = 0;
  int octant = -1;

  comm.setSizes(subdomain.total_overload, quad.num_angleset, num_groupsets);

  Task_IDs.resize(num_tasks);

  for (int i = 0; i < subdomain.total_overload; i++)
  {
    for (int j = 0; j < quad.num_angleset; j++)
    {   
      for (int k = 0; k < num_groupsets; k++)
      {
        Task task;
        task.BuildTask(rank, this, &(subdomain), i,j,k);
        All_Tasks.push_back(task);
        task_num += 1;
      }
    }

  }

  // Give pointers to the tasks to the comm (so we can figure out neighbors)
  for(int i = 0; i < num_tasks; i++){
    All_Tasks[i].SetIndex(i);
    comm.SetTask(&(All_Tasks[i]));
  }

  double start_sort;
  long double duration_sort;
  start_sort = MPI_Wtime();

  // Now we order the All_Tasks vector by depth, then angleset, then groupset
  std::sort(All_Tasks.begin(), All_Tasks.end(), by_depth());

  duration_sort = (MPI_Wtime() - start_sort); 
  if (rank == 0){
 //   std::cout << "Sort took " << duration_sort << " seconds." << std::endl;
  }

  for(int i = 0 ; i < All_Tasks.size(); i++){
    Task_IDs[All_Tasks[i].task_id] = i;
  }

  for(int i = 0; i < num_tasks; i++){
    All_Tasks[i].AllocateBuffers();
  }
 
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
  double start_task, start_cell, start_wait, start_send;
  long double duration_task(0), duration_innertask(0), duration_cell(0);
  long double duration_wait(0), duration_send(0);
  // First off we need to add all the tasks to the comm:
  for(int i = 0; i < Task_IDs.size(); i++){
    int index = Task_IDs[i];
    comm.addTask(All_Tasks[index].task_id, All_Tasks[index]);
  }
  // Necessary temporary data structures for the sweep
  std::vector<int> cell_ijk(3, 0);
  std::vector<double> temp_solve(4, 0);
  Neighbor neighbor;
  int task_it(0), recv(0), target(0);

  double start_loop, start_angle, start_group;
  long double duration_loop;
  long double total_duration_loop;
  long double duration_angle(0), duration_group(0);
  total_duration_loop = 0;

  // Loop through task list
  std::vector<Task>::iterator it = All_Tasks.begin();
  std::vector<Task>::iterator it_end = All_Tasks.end();
  for (; it != it_end; it++, task_it++)
  {
    start_loop = MPI_Wtime();

    // Get cellset and angleset information
    int cells_x = subdomain.CellSets[(*it).cellset_id_loc].cells_x;
    int cells_y = subdomain.CellSets[(*it).cellset_id_loc].cells_y;
    int cells_z = subdomain.CellSets[(*it).cellset_id_loc].cells_z;
    int cell_per_cellset = (*it).cell_per_cellset;
    int angle_per_angleset = (*it).angle_per_angleset;
    
    int octant = quad.Anglesets[(*it).angleset_id].octant;
    vector<vector<int> > incoming = (*it).incoming;
    vector<vector<int> > outgoing = (*it).outgoing;

    start_wait = MPI_Wtime();
    comm.WaitforReady((*it).task_id);
    duration_wait += MPI_Wtime() - start_wait;

    start_task = MPI_Wtime();
    
    // We march through the cells first in x, then y, then z
    // The order (left to right or right to left) depends on what 
    // octant the angleset is in so we call GetCell to figure out where we are
    for (int k = 0; k < cells_z; k++)
    {
      for (int j = 0; j < cells_y; j++)
      {
        for (int i = 0; i < cells_x; i++)
        {
          start_cell = MPI_Wtime();
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
            start_angle = MPI_Wtime();
            
            // Get the direction of the angle
            omega = quad.Anglesets[(*it).angleset_id].Omegas[m];
            
            // Add in the gradient matrix to the A matrix
            for (int a = 0; a < 4; a++)
            {
              for (int b = 0; b < 4; b++)
              {
   //             A_tilde[a][b] = dot(omega, L[a][b]);
                A_tilde[a][b] = omega.x * L[a][b].x + omega.y * L[a][b].y + omega.z * L[a][b].z;
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
   //               A_tilde[a][b] += dot(-1 * omega, N[incoming[f][0]][a][b]);
                    A_tilde[a][b] +=  (-1 * omega.x * N[incoming[f][0]][a][b].x + -1 * omega.y * N[incoming[f][0]][a][b].y  + -1 * omega.z * N[incoming[f][0]][a][b].z);
                }
              }            
            }
            
            // Loop over groups in this groupset
            for (int g = 0; g < group_per_groupset; g++)
            {
              start_group = MPI_Wtime();
              // Initialize the b vector
              for (int a = 0; a < 4; a++)
                bg[a] = RHS[a];

              // Retrieve this cells sigma tot (this is in the group loop to simulate
              // multi-group cross sections
              sigma_t = my_cell.GetSigmaTot();

              // Need to get incoming fluxes for the RHS
              for (int f = 0; f < 3; f++)
              {
                (*it).Get_buffer(cell_ijk[0], cell_ijk[1], cell_ijk[2], g, m, incoming[f][0], temp_solve);
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

              // Solve A^-1*RHS with Gaussian Elimination and store it in RHS (4 = number of elements)
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
                facecenter = my_cell.facecenters[outgoing[f][0]];

                bg[0] = cell_average + facecenter.x*bg[1] + facecenter.y*bg[2] + facecenter.z*bg[3];
                // Now store outgoing fluxes in the buffer arrays
                (*it).Set_buffer(cell_ijk[0], cell_ijk[1], cell_ijk[2], g, m, outgoing[f][0], task_it, bg);
              }
              duration_group += MPI_Wtime() - start_group;
            } //groups
            duration_angle += MPI_Wtime() - start_angle;
          } //angles
          duration_cell += (MPI_Wtime() - start_cell);
        } // cells in x
      } // cells in y
    } // cells in z
    
    duration_task += (MPI_Wtime() - start_task);
    start_send = MPI_Wtime();
    comm.markComplete((*it).task_id);
    duration_send += MPI_Wtime() - start_send;
     
    duration_loop = (MPI_Wtime() - start_loop);
    total_duration_loop += duration_loop;
  } // tasks
  duration_task -= duration_cell;
  duration_cell -= duration_angle;
  duration_angle -= duration_group;

  if(rank == 0 && verbose){
    printf("Total Duration: %Lf\n", total_duration_loop);
    printf("          Task: %Lf\n", duration_task);
    printf("          Cell: %Lf\n", duration_cell);
    printf("         Angle: %Lf\n", duration_angle);
    printf("         Group: %Lf\n", duration_group);
    printf("          Wait: %Lf\n", duration_wait);
    printf("          Send: %LF\n", duration_send);
    printf("             \n");
  }

  comm.CleanUp();
}

int Problem::GetCell(int i, int j, int k, int cells_x, int cells_y, int cells_z, int octant)
{
  int d = -1;

  if (octant == 0){
    d = i + cells_x*j + cells_x*cells_y*k;
    return d;
  }
  if (octant == 1){
    d = (cells_x - i -1) + cells_x*j + cells_x*cells_y*k;
    return d;
  }
  if (octant == 2){
    d = (cells_x - i -1) + cells_x*(cells_y - j -1) + cells_x*cells_y*k;
    return d;
  }
  if (octant == 3){
    d = i + cells_x*(cells_y - j -1) + cells_x*cells_y*k;
    return d;
  }
  if (octant == 4){
    d = i + cells_x*j + cells_x*cells_y*(cells_z - k - 1);
    return d;
  }
  if (octant == 5){
    d = (cells_x - i -1) + cells_x*j + cells_x*cells_y*(cells_z - k -1);
    return d;
  }
  if (octant == 6){
    d = (cells_x - i -1) + cells_x*(cells_y - j -1) + cells_x*cells_y*(cells_z - k -1);
    return d;
  }
  if (octant == 7){
    d = i + cells_x*(cells_y - j -1) + cells_x*cells_y*(cells_z - k -1);
    return d;
  }


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

  for(int t = 0; t < num_tasks; t++){
    for (int i = 0; i < All_Tasks[t].plane_data.size(); i++)
    {
      for (int j = 0; j < All_Tasks[t].plane_data[i].size(); j++)
        All_Tasks[t].plane_data[i][j] = 0;
    }
  }
}
