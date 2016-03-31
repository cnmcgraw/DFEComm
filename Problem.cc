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
#include <string.h>

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
        if (a.omega[0] == b.omega[0]){
          if (a.omega[1] == b.omega[1]){
            return a.omega[2] > b.omega[2];
          }
          else
            return a.omega[1] > b.omega[1];
        }
        else
          return a.omega[0] > b.omega[0];
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
  spider = input->spider;
  if(spider){
    sp_azim.resize(refinement);
    sp_azim = input->sp_azim;
  }
  z_planes = input->z_planes;
  sp_disc = input->sp_disc;
  bcs = input->bcs;
  max_faces = 6;


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
  partition_function.resize(3);
  partition_function = input->partition_function;
  overload.resize(3);
  overload = input->overload;
  num_cellsets.resize(3);
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
  
  // This will be the data structure we'll be writing into.
  // We'll pull from plane_data at the beginning of the task and write into it at the end.
  // This structure assumes that all tasks are the same size.
  interior_data.resize(All_Tasks[0].cells_xy*All_Tasks[0].cells_z*max_faces*All_Tasks[0].group_per_groupset*All_Tasks[0].angle_per_angleset * 4);  
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
  
  std::vector<Task>::iterator it = All_Tasks.begin();
  std::vector<Task>::iterator it_end = All_Tasks.end();
  
  // Necessary temporary data structures for the sweep
  std::vector<int> cell_ijk(2, 0), neighbor_ijk(2,0);
  std::vector<double> temp_solve(4, 0);
  Neighbor neighbor;
  int task_it(0), recv(0), target(0);
  
  // Cell Loop
  int cells_x, cells_y, cells_xy, cells_z;
  int cell_x, cell_y;
  int cell_per_cellset, angle_per_angleset, octant;
  vector<vector<int> > tincoming(3,vector<int>(2,0)), toutgoing(3, std::vector<int>(3,0));
  vector<int> incoming(3,0);
  vector<vector<int> > outgoing(3, vector<int>(3,0));
  int num_inc_faces, num_out_faces;
  double cell_average;
  int it_index = 0;
  int if1, if2, if3, of;
  
  // Angle Loop
  int cell_id, cijk0, cijk1, cijk2;
  omega.resize(All_Tasks[0].angle_per_angleset, std::vector<double>(3,0.0));
  double ox, oy, oz;
  double inf;
 // vector<vector<double> > omN(4,vector<double>(4,0));
  
  // Group Loop
  int index;
  std::vector<double> wt((*it).angle_per_angleset,0.0);
  std::vector<std::vector<double>::iterator> bloc(group_per_groupset);
  
  double sig_tot;

  // Loop through task list
  for (; it != it_end; it++, task_it++)
  {
    start_loop = MPI_Wtime();

    // Get cellset and angleset information
    cells_x = subdomain.CellSets[(*it).cellset_id_loc].cells_x;
    cells_y = subdomain.CellSets[(*it).cellset_id_loc].cells_y;
    cells_xy = subdomain.CellSets[(*it).cellset_id_loc].cells_xy;
    cells_z = subdomain.CellSets[(*it).cellset_id_loc].cells_z;
    cell_per_cellset = (*it).cell_per_cellset;
    angle_per_angleset = (*it).angle_per_angleset;
    
    octant = quad.Anglesets[(*it).angleset_id].octant;
    
    for (int m = 0; m < angle_per_angleset; m++){
      omega[m][0] = quad.Anglesets[(*it).angleset_id].Omegas[m][0];
      omega[m][1] = quad.Anglesets[(*it).angleset_id].Omegas[m][1];
      omega[m][2] = quad.Anglesets[(*it).angleset_id].Omegas[m][2];
      wt[m] = quad.Anglesets[(*it).angleset_id].Weights[m];
    }

    start_wait = MPI_Wtime();
    comm.WaitforReady((*it).task_id);
    duration_wait += MPI_Wtime() - start_wait;

    start_task = MPI_Wtime();
    
    // Need to refill interior_data with zeros
    std::fill(interior_data.begin(), interior_data.end(), 0);
    
    // Need to pull in data to our interior data structure:
    (*it).GetInteriorData(interior_data);

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
          Cell &my_cell = *(subdomain.CellSets[(*it).cellset_id_loc].Cells[cell_id]);
          my_cell.GetCellijk(cell_id, cells_x, cells_y, cells_z, cell_ijk);
          cell_y = (int)(cell_ijk[0] / cells_x);
          cell_x = cell_ijk[0] - cells_x * cell_y;
          sig_tot = my_cell.GetSigmaTot();

          // Get the cell's DFEM Matrices for building the A matrix
          std::vector<std::vector< double > > &M = my_cell.M;
          std::vector<std::vector<std::vector< vector<double> > > > &N = my_cell.N;
          std::vector<std::vector< vector<double> > > &L = my_cell.L;
         

          // Since the source is piece-wise constant, we only need to update the average value
          source[0] = my_cell.GetSource();
          
          // We need the incoming and outgoing face for this cell
          my_cell.GetCellInOut(omega[0], incoming, outgoing);
          
          num_inc_faces = incoming.size();
          num_out_faces = outgoing.size();

          // Loop over angles in the angleset
          for (int m = 0; m < angle_per_angleset; m++)
          {
            start_angle = MPI_Wtime();
            
            ox = omega[m][0];
            oy = omega[m][1];
            oz = omega[m][2];
                
            // Add in the gradient matrix to the A matrix
            for (int a = 0; a < 4; a++)
            {
              for (int b = 0; b < 4; b++)
              { 
                A_tilde[a][b] = ox * L[a][b][0] + oy * L[a][b][1] + oz * L[a][b][2];
              }
              RHS[a] = M[a][a] * source[a];
            }
            // Add the surface matrix contribution for incoming cell faces
            for (int f = 0; f < num_inc_faces; f++)
            {
              inf = incoming[f];
              for (int a = 0; a < 4; a++)
              {
                for (int b = 0; b < 4; b++)
                { 
                  A_tilde[a][b] += -ox * N[ inf ][a][b][0] - oy * N[ inf ][a][b][1] - oz * N[ inf ][a][b][2];
                }
              }
            }

            // Loop over groups in this groupset
            for (int g = 0; g < group_per_groupset; g++)
            {
              start_group = MPI_Wtime();
              
              for(int a = 0; a < 4; a++)
                bg[a] = RHS[a];

              // Initialize b and add incoming fluxes for the RHS
              for (int f = 0; f < num_inc_faces; f++)
              {
                inf = incoming[f];
                
                it_index = (cell_ijk[1]*cells_xy*group_per_groupset*angle_per_angleset * 4*6 + 
                  cell_y*cells_x*group_per_groupset*angle_per_angleset * 4*6 + cell_x*group_per_groupset*angle_per_angleset * 4*6 + inf*group_per_groupset*angle_per_angleset * 4 + g*angle_per_angleset * 4 + m* 4);
                
               // it_index = (*it).Get_buffer_loc(cell_ijk[0], cell_ijk[1], g, m, incoming[f]);
                for (int a = 0; a < 4; a++)
                {
                  for (int b = 0; b < 4; b++)
                  {
                      bg[a] += (-ox * N[ inf ][a][b][0] - oy * N[ inf][a][b][1] - oz * N[ inf ][a][b][2])*interior_data[it_index + b];
                  }
                } 
              }
             
              // Add the contribution to the A matrix and the RHS vector
              for (int a = 0; a < 4; a++)
              {
                for (int b = 0; b < 4; b++)
                  A[a][b] = A_tilde[a][b] + sig_tot * M[a][b];
                // Retrieve this cells sigma tot (this is in the group loop to simulate multi-group cross sections)
              //  A[a][a] += sig_tot * M[a][a];
              }

              // Solve A^-1*RHS with Gaussian Elimination and store it in RHS (4 = number of elements)
              GE_no_pivoting(A, bg, 4);

              // Now we accumulate the fluxes into phi;
              index = (*it).groupset_id*group_per_groupset*4+4*g;
              for (int p = 0; p < 4; p++)
                my_cell.phi[index + p] += bg[p] * wt[m];
                
              

              // Now we need to translate the cell average to the average on each 
              // face before we push to the down stream neighbers
              // This allows for direct data movement (no cell to cell mapping needed)
              cell_average = bg[0];
              for (int f = 0; f < num_out_faces; f++)
              {
                of = outgoing[f][1];
                bg[0] = cell_average + my_cell.facecenters[ of ][0]*bg[1] + my_cell.facecenters[ of ][1]*bg[2] + my_cell.facecenters[ of ][2]*bg[3];
                my_cell.GetCellijk(outgoing[f][0], cells_x, cells_y, cells_z, neighbor_ijk);
                
                
                cell_y = (int)(neighbor_ijk[0] / cells_x);
                cell_x = neighbor_ijk[0] - cells_x * cell_y;
                it_index = (neighbor_ijk[1]*cells_xy*group_per_groupset*angle_per_angleset * 4*6 + 
                  cell_y*cells_x*group_per_groupset*angle_per_angleset * 4*6 + cell_x*group_per_groupset*angle_per_angleset * 4*6 + outgoing[f][2]*group_per_groupset*angle_per_angleset * 4 + g*angle_per_angleset * 4 + m* 4);
                
                
                //it_index = (*it).Get_buffer_loc(neighbor_ijk[0], neighbor_ijk[1], g, m, outgoing[f][2]);
                
                for(int i = 0;i < 4;++i)
                  interior_data[it_index + i] = bg[i];
              }
              duration_group += MPI_Wtime() - start_group;
            } //groups
            duration_angle += MPI_Wtime() - start_angle;
          } //angles
          duration_cell += (MPI_Wtime() - start_cell);
        } // cells in x
      } // cells in y
    } // cells in z
    // Now we need to set the plane_data structure with our interior data 
    
    (*it).SetPlaneData(interior_data);
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
      for (int k = 0; k < subdomain.CellSets[i].Cells[j]->phi.size(); k++)
          subdomain.CellSets[i].Cells[j]->phi[k] = 0;
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
        phi_average[0] += subdomain.CellSets[i].Cells[j]->phi[4*index];
        phi_average[1] += subdomain.CellSets[i].Cells[j]->phi[4*index+1];
        phi_average[2] += subdomain.CellSets[i].Cells[j]->phi[4*index+2];
        phi_average[3] += subdomain.CellSets[i].Cells[j]->phi[4*index+3];
        size += 1;
      }
    }
  //  std::cout << " " << std::endl;
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
        phi_stddev[0] += pow(subdomain.CellSets[i].Cells[j]->phi[4*index] - phi_average[0],2);
        phi_stddev[1] += pow(subdomain.CellSets[i].Cells[j]->phi[4*index+1]- phi_average[1],2);
        phi_stddev[2] += pow(subdomain.CellSets[i].Cells[j]->phi[4*index+2]- phi_average[2],2);
        phi_stddev[3] += pow(subdomain.CellSets[i].Cells[j]->phi[4*index+3]- phi_average[3],2);        
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

