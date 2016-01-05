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

  int task_num = 0;
  int octant = -1;

  comm.setSizes(subdomain.total_overload, quad.num_angleset, num_groupsets);

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

  // Now we order the All_Tasks vector by depth, then angleset, then groupset
  std::sort(All_Tasks.begin(), All_Tasks.end(), by_depth());

  
    // Give pointers to the tasks to the comm (so we can figure out neighbors)
  for(int i = 0; i < num_tasks; i++){
    comm.SetTask(&(All_Tasks[i]));
  }

  for(int i = 0; i < num_tasks; i++){
    All_Tasks[i].AllocateBuffers();
  }
 
  // Pre-allocate the DFEM Matrices
  A_tilde.resize(4, std::vector<double>(4, 0.0));
  A.resize(4, std::vector<double>(4, 0.0));
  RHS.resize(4);
  bg.resize(4);

  // Resize the source. Piece-wise constant in each cell
  source.resize(4);
  
}

void Problem::Sweep(std::ofstream &output)
{
  // Timers
  double start_task, start_cell, start_wait, start_send;
  double start_loop, start_angle, start_group;
  long double duration_task(0), duration_cell(0);
  long double duration_wait(0), duration_send(0);
  long double total_duration_loop(0), duration_angle(0), duration_group(0);
  
  // First off we need to add all the tasks to the comm:
  for(int i = 0; i < All_Tasks.size(); i++){
    comm.addTask(All_Tasks[i].task_id, All_Tasks[i]);
  }
  // Necessary temporary data structures for the sweep
  std::vector<int> cell_ijk(3, 0);
  std::vector<double> temp_solve(4, 0);
  Neighbor neighbor;
  int task_it(0), recv(0), target(0);
  
  // Cell Loop
  int cells_x, cells_y, cells_z;
  int cell_per_cellset, angle_per_angleset, octant;
  vector<vector<int> > incoming, outgoing;
  int if1, if2, if3, of1, of2, of3;
  
  // Angle Loop
  int cell_id, cijk0, cijk1, cijk2;
  double ox, oy, oz;
  
  // Group Loop
  int index;
  double wt;

  // Loop through task list
  std::vector<Task>::iterator it = All_Tasks.begin();
  std::vector<Task>::iterator it_end = All_Tasks.end();
  for (; it != it_end; it++, task_it++)
  {
    start_loop = MPI_Wtime();

    // Get cellset and angleset information
    cells_x = subdomain.CellSets[(*it).cellset_id_loc].cells_x;
    cells_y = subdomain.CellSets[(*it).cellset_id_loc].cells_y;
    cells_z = subdomain.CellSets[(*it).cellset_id_loc].cells_z;
    cell_per_cellset = (*it).cell_per_cellset;
    angle_per_angleset = (*it).angle_per_angleset;
    
    octant = quad.Anglesets[(*it).angleset_id].octant;
    incoming = (*it).incoming;
    outgoing = (*it).outgoing;
    int if1 = incoming[0][0], if2 = incoming[1][0], if3 = incoming[2][0];
    int of1 = outgoing[0][0], of2 = outgoing[1][0], of3 = outgoing[2][0];

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
          cell_id = GetCell(i, j, k, cells_x, cells_y, cells_z, octant);
          Cell &my_cell = subdomain.CellSets[(*it).cellset_id_loc].Cells[cell_id];
          my_cell.GetCellijk(cell_id, cells_x, cells_y, cells_z, cell_ijk);
          cijk0 = cell_ijk[0]; 
          cijk1 = cell_ijk[1]; 
          cijk2 = cell_ijk[2];

          // Get the cell's DFEM Matrices for building the A matrix
          std::vector< double > &M = subdomain.CellSets[(*it).cellset_id_loc].Cells[cell_id].M;
          std::vector<std::vector<std::vector< Direction > > > &N = subdomain.CellSets[(*it).cellset_id_loc].Cells[cell_id].N;
          std::vector<std::vector< Direction > > &L = subdomain.CellSets[(*it).cellset_id_loc].Cells[cell_id].L;

          // Since the source is piece-wise constant, we only need to update the average value
          source[0] = my_cell.GetSource();
          
          // Loop over angles in the angleset
          for (int m = 0; m < angle_per_angleset; m++)
          {
            start_angle = MPI_Wtime();
            
            // Get the direction of the angle
            omega = quad.Anglesets[(*it).angleset_id].Omegas[m];
            ox = omega.x;
            oy = omega.y;
            oz = omega.z;
            wt = quad.Anglesets[(*it).angleset_id].Weights[m];
            
            // Add in the gradient matrix to the A matrix
            for (int a = 0; a < 4; a++)
            {
              for (int b = 0; b < 4; b++)
              { 
                // Add the surface matrix contribution for incoming cell faces
                A_tilde[a][b] = ( ox * L[a][b].x      +  oy * L[a][b].y       +  oz * L[a][b].z) 
                              + (-ox * N[if1][a][b].x + -oy * N[if1][a][b].y  + -oz * N[if1][a][b].z)
                              + (-ox * N[if2][a][b].x + -oy * N[if2][a][b].y  + -oz * N[if2][a][b].z)
                              + (-ox * N[if3][a][b].x + -oy * N[if3][a][b].y  + -oz * N[if3][a][b].z);
              }
            }

            // Loop over groups in this groupset
            for (int g = 0; g < group_per_groupset; g++)
            {
              start_group = MPI_Wtime();

              // Initialize the b vector
              for (int a = 0; a < 4; a++)
                bg[a] = M[a] * source[a];

              // Need to get incoming fluxes for the RHS
              (*it).Get_buffer(cijk0, cijk1, cijk2, g, m, if1, temp_solve);
              for (int a = 0; a < 4; a++)
                bg[a] += (-ox * N[if1][a][0].x + -oy * N[if1][a][0].y  + -oz * N[if1][a][0].z)*temp_solve[0]
                      +  (-ox * N[if1][a][1].x + -oy * N[if1][a][1].y  + -oz * N[if1][a][1].z)*temp_solve[1]
                      +  (-ox * N[if1][a][2].x + -oy * N[if1][a][2].y  + -oz * N[if1][a][2].z)*temp_solve[2]
                      +  (-ox * N[if1][a][3].x + -oy * N[if1][a][3].y  + -oz * N[if1][a][3].z)*temp_solve[3];

              (*it).Get_buffer(cijk0, cijk1, cijk2, g, m, if2, temp_solve);
              for (int a = 0; a < 4; a++)
                bg[a] += (-ox * N[if2][a][0].x + -oy * N[if2][a][0].y  + -oz * N[if2][a][0].z)*temp_solve[0]
                      +  (-ox * N[if2][a][1].x + -oy * N[if2][a][1].y  + -oz * N[if2][a][1].z)*temp_solve[1]
                      +  (-ox * N[if2][a][2].x + -oy * N[if2][a][2].y  + -oz * N[if2][a][2].z)*temp_solve[2]
                      +  (-ox * N[if2][a][3].x + -oy * N[if2][a][3].y  + -oz * N[if2][a][3].z)*temp_solve[3];

              (*it).Get_buffer(cijk0, cijk1, cijk2, g, m, if3, temp_solve);
              for (int a = 0; a < 4; a++)
                bg[a] += (-ox * N[if3][a][0].x + -oy * N[if3][a][0].y  + -oz * N[if3][a][0].z)*temp_solve[0]
                      +  (-ox * N[if3][a][1].x + -oy * N[if3][a][1].y  + -oz * N[if3][a][1].z)*temp_solve[1]
                      +  (-ox * N[if3][a][2].x + -oy * N[if3][a][2].y  + -oz * N[if3][a][2].z)*temp_solve[2]
                      +  (-ox * N[if3][a][3].x + -oy * N[if3][a][3].y  + -oz * N[if3][a][3].z)*temp_solve[3];
          
              // Add the contribution to the A matrix and the RHS vector
              for (int a = 0; a < 4; a++)
              {
                for (int b = 0; b < 4; b++)
                  A[a][b] = A_tilde[a][b];  
                // Retrieve this cells sigma tot (this is in the group loop to simulate multi-group cross sections)
                A[a][a] += my_cell.GetSigmaTot() * M[a];
              }

              // Solve A^-1*RHS with Gaussian Elimination and store it in RHS (4 = number of elements)
              GE_no_pivoting(A, bg, 4);

              // Now we accumulate the fluxes into phi;
              index = (*it).groupset_id*group_per_groupset*4+4*g;
              for (int p = 0; p < 4; p++)
                my_cell.phi[index + p] += bg[p] * wt;

              // Now we need to translate the cell average to the average on each 
              // face before we push to the down stream neighbers
              // This allows for direct data movement (no cell to cell mapping needed)
              bg[0] += my_cell.facecenters[of1].x*bg[1] + my_cell.facecenters[of1].y*bg[2] + my_cell.facecenters[of1].z*bg[3];
              (*it).Set_buffer(cijk0, cijk1, cijk2, g, m, of1, task_it, bg);

              bg[0] += my_cell.facecenters[of2].x*bg[1] + my_cell.facecenters[of2].y*bg[2] + my_cell.facecenters[of2].z*bg[3];
              (*it).Set_buffer(cijk0, cijk1, cijk2, g, m, of2, task_it, bg);

              bg[0] += my_cell.facecenters[of3].x*bg[1] + my_cell.facecenters[of3].y*bg[2] + my_cell.facecenters[of3].z*bg[3];
              (*it).Set_buffer(cijk0, cijk1, cijk2, g, m, of3, task_it, bg);
              
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
     
    total_duration_loop += (MPI_Wtime() - start_loop);
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
    ComputePhi();
    std::cout << "=============================" << std::endl;
  }

  //comm.CleanUp();
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

void Problem::ComputePhi()
{

  // First compute the average phi across the problem for all groups
  std::vector<double> phi_average(4,0);
  std::vector<double> phi_stddev(4,0);
  int size = 0;
  for (int i = 0; i < subdomain.CellSets.size(); i++)
  {
    // Loop through cells in this cellset
    for (int j = 0; j < subdomain.CellSets[i].Cells.size(); j++)
    { 
      for (int index = 0; index < num_groupsets*group_per_groupset; index++)
      {
        phi_average[0] += subdomain.CellSets[i].Cells[j].phi[4*index];
        phi_average[1] += subdomain.CellSets[i].Cells[j].phi[4*index+1];
        phi_average[2] += subdomain.CellSets[i].Cells[j].phi[4*index+2];
        phi_average[3] += subdomain.CellSets[i].Cells[j].phi[4*index+3];
        size += 1;
      }
    }
  }
  phi_average[0] /= size;
  phi_average[1] /= size;
  phi_average[2] /= size;
  phi_average[3] /= size;
  
  // Now compute the standard deviation of phi across the problem for all groups
  for (int i = 0; i < subdomain.CellSets.size(); i++)
  {
    // Loop through cells in this cellset
    for (int j = 0; j < subdomain.CellSets[i].Cells.size(); j++)
    { 
      for (int index = 0; index < num_groupsets*group_per_groupset; index++)
      {
        phi_stddev[0] += pow(subdomain.CellSets[i].Cells[j].phi[4*index] - phi_average[0],2);
        phi_stddev[1] += pow(subdomain.CellSets[i].Cells[j].phi[4*index+1]- phi_average[1],2);
        phi_stddev[2] += pow(subdomain.CellSets[i].Cells[j].phi[4*index+2]- phi_average[2],2);
        phi_stddev[3] += pow(subdomain.CellSets[i].Cells[j].phi[4*index+3]- phi_average[3],2);        
      }
    }
  }
  
  phi_stddev[0] = sqrt(phi_stddev[0]/size);
  phi_stddev[1] = sqrt(phi_stddev[1]/size);
  phi_stddev[2] = sqrt(phi_stddev[2]/size);
  phi_stddev[3] = sqrt(phi_stddev[3]/size);  

  std::cout << "Average Phi: " << std::endl;
  std::cout << "       Phi : " << phi_average[0] << std::endl;
  std::cout << "   DPhi/Dx : " << phi_average[1] << std::endl;
  std::cout << "   DPhi/Dy : " << phi_average[2] << std::endl;
  std::cout << "   DPhi/Dz : " << phi_average[3] << std::endl; 
  std::cout << "\n" ;    
  std::cout << "Standard Deviation: " << std::endl;
  std::cout << "       Phi : " << phi_stddev[0] << std::endl;
  std::cout << "   DPhi/Dx : " << phi_stddev[1] << std::endl;
  std::cout << "   DPhi/Dy : " << phi_stddev[2] << std::endl;
  std::cout << "   DPhi/Dz : " << phi_stddev[3] << std::endl; 
  std::cout << "\n" ;  

}

