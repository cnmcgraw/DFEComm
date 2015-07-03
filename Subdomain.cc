#include "Subdomain.h"
#include <stdio.h>
#include <iostream>
#include <vector>


Subdomain::Subdomain()
{ }

Subdomain::~Subdomain()
{ }

void Subdomain::BuildSubdomain(int SML_ID, Problem* problem)
{
  ComputeCellSetID(SML_ID, problem);

  CellSets.resize(total_overload);
  for (int i = 0; i < total_overload; i++)
    CellSets[i].BuildCellSet(CellSetIDs[i], problem);
    
  cells_x = 2 * problem->num_pin_x*problem->refinement / problem->num_cellsets[0];
  cells_y = 2 * problem->num_pin_y*problem->refinement / problem->num_cellsets[1];
  cells_z = problem->z_planes / problem->num_cellsets[2];

  // This needs to be the biggest angleset
  angle_per_angleset = problem->quad.Anglesets[0].angle_per_angleset;
  group_per_groupset = problem->group_per_groupset;

  SetBoundaryConditions(problem);
}
void Subdomain::ComputeCellSetID(int SML_ID, Problem* problem)
{
  // The total overload is the total number of 
  // cellsets each SML owns
  total_overload = 1;
  std::vector<int> P;
  P.resize(3);
  // For 2D P_eff in z = 1 and overload in z = 1
  for(int i=0; i< 3; i++)
  {
    total_overload *= problem->overload[i];
    P[i] = problem->num_cellsets[i]/problem->overload[i];
  }

  CellSetIDs.resize(total_overload);
  
  // Compute the i,j,k of the initial cellset
  int k = (int)(SML_ID/(P[0]*P[1]));
  int j = (int)((SML_ID-k*P[0]*P[1])/P[0]);
  int i = SML_ID - j*P[0] - k*P[0]*P[1];

  // If the mesh is blocked in any dimension,
  // we must multiply by the overload in
  // that dimension
  if(problem->partition_function[2] == 1)
    k = k*problem->overload[2];
  if(problem->partition_function[1] == 1)
    j = j*problem->overload[1]; 
  if(problem->partition_function[0] == 1)
    i = i*problem->overload[0];

  std::vector<int> ID_ijk;
  ID_ijk.resize(3);
  int id = 0;
  for(int z = 0; z<problem->overload[2]; z++)
  {
    // Sequential in blocked
    if(problem->partition_function[2] == 1)
      ID_ijk[2] = k + z;
    // Jump by the overload factor for round robin
    else 
      ID_ijk[2] = k + z*problem->overload[2];

    for(int y = 0; y<problem->overload[1]; y++)
    {
      if(problem->partition_function[1] == 1)
        ID_ijk[1] = j + y;
      else 
        ID_ijk[2] = j + y*problem->overload[1];

      for(int x = 0; x<problem->overload[0]; x++, id++)
      {
        if(problem->partition_function[0] == 1)
          ID_ijk[0] = i + x;
        else 
          ID_ijk[2] = i + x*problem->overload[0];

        CellSetIDs[id] = ID_ijk[0] + problem->num_cellsets[0]*ID_ijk[1]
          + problem->num_cellsets[0]*problem->num_cellsets[1]*ID_ijk[2];
      }
    }
  }


}
void Subdomain::SetBoundaryConditions(Problem* problem)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (problem->bcs == 1)
  {
    if (rank == 0){ std::cout << "Reflecting boundary conditions not implemented yet" << std::endl; }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  else
  {
    bc.resize(3);
    for (int i = 0; i < bc.size(); i++)
    {
      bc[i] = 7. / 3.;
    }
  }

}
double Subdomain::GetBoundaryCondition(int Boundary)
{
  return bc[Boundary];
}
void Subdomain::AllocateBuffers(int num_tasks)
{

  buffer.resize(3);
  buffer[0].resize(cells_y*cells_z*group_per_groupset*angle_per_angleset * 4);
  buffer[1].resize(cells_x*cells_z*group_per_groupset*angle_per_angleset * 4);
  buffer[2].resize(cells_x*cells_y*group_per_groupset*angle_per_angleset * 4);

  Send_buffer.resize(3);
  Send_buffer[0].resize(cells_y*cells_z*group_per_groupset*angle_per_angleset * 4);
  Send_buffer[1].resize(cells_x*cells_z*group_per_groupset*angle_per_angleset * 4);
  Send_buffer[2].resize(cells_x*cells_y*group_per_groupset*angle_per_angleset * 4);
  
  int max_size = 0; 
  for(int i = 0; i < 3; i++)
    if(buffer[i].size() > max_size){max_size = buffer[i].size();}

  // This is a naive way of sizing this buffer. It assumes every task receives 3 messages (one from x, y, and z)
  // This is not always the case (some from bc's) so this should be made more efficient
  Received_buffer.resize(3 * num_tasks);
  for(int i = 0; i < num_tasks; i++)
    for(int j = 0; j < 3; j++)
       Received_buffer[i + j].resize(buffer[j].size());

  // Vector of (tag, source, count)'s
  Received_info.resize(3 * num_tasks, std::vector<int>(3,0));

  // This is a queue of open places int the Recieved buffer. At first it will contain the location
  // for every chunk. As the buffer gets filled, the queue will shrink. When chunks in the Received
  // buffer get used, their location gets put back in the queue for reuse. If there are no open locations
  // left in the queue, the Received buffer will be resized.
  Received_open = std::queue<int>();
  for (int i = 0; i < 3 * num_tasks; i++)
    Received_open.push(i);


}
void Subdomain::Set_buffer(int cell_x, int cell_y, int cell_z, int group, int angle, int face, int task, vector<double>& RHS)
{
  if (face == 0 || face == 1)
  {
    for (int i = 0; i < 4; i++)
    {
      Send_buffer[0][cell_y*cells_z*group_per_groupset*angle_per_angleset * 4 + cell_z*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4 + i] = RHS[i];
    }
  }
  else if (face == 2 || face == 3)
  {
    for (int i = 0; i < 4; i++)
    {
      Send_buffer[1][cell_x*cells_z*group_per_groupset*angle_per_angleset * 4 + cell_z*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4 + i] = RHS[i];
    }
  }
  else if (face == 4 || face == 5)
  {
    for (int i = 0; i < 4; i++)
    {
      Send_buffer[2][cell_x*cells_y*group_per_groupset*angle_per_angleset * 4 + cell_y*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4 + i] = RHS[i];
    }
  }
}
void Subdomain::Get_buffer(int cell_x, int cell_y, int cell_z, int group, int angle, int face, vector<double>& RHS)
{
  if (face == 0 || face == 1)
  {
    for (int i = 0; i < 4; i++)
    {
      RHS[i] = buffer[0][cell_y*cells_z*group_per_groupset*angle_per_angleset * 4 + cell_z*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4 + i];
    }
  }
  else if (face == 2 || face == 3)
  {
    for (int i = 0; i < 4; i++)
    {
      RHS[i] = buffer[1][cell_x*cells_z*group_per_groupset*angle_per_angleset * 4 + cell_z*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4 + i];
    }
  }
  else
  {
    for (int i = 0; i < 4; i++)
    {
      RHS[i] = buffer[2][cell_x*cells_y*group_per_groupset*angle_per_angleset * 4 + cell_y*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4 + i];
    }
  }
}
void Subdomain::Get_buffer_from_bc(int face)
{

  if (face == 0 || face == 1)
  {
    // Need to set the cell averages, all the slopes will be zero
    int size = cells_y*cells_z*group_per_groupset*angle_per_angleset;
    for (int i = 0; i < size; i++)
    {
      int j = 4 * i;
      buffer[0][j] = GetBoundaryCondition(0);
      buffer[0][j + 1] = 0;
      buffer[0][j + 2] = 0;
      buffer[0][j + 3] = 0;
      
    }
  }
  else if (face == 2 || face == 3)
  {
    // Need to set the cell averages, all the slopes will be zero
    int size = cells_x*cells_z*group_per_groupset*angle_per_angleset;
    for (int i = 0; i < size; i++)
    {
      int j = 4 * i;
      buffer[1][j] = GetBoundaryCondition(1);
      buffer[1][j + 1] = 0;
      buffer[1][j + 2] = 0;
      buffer[1][j + 3] = 0;
    }
  }
  else
  {
    int size = cells_x*cells_y*group_per_groupset*angle_per_angleset;
    for (int i = 0; i < size; i++)
    {
      int j = 4 * i;
      buffer[2][j] = GetBoundaryCondition(2);
      buffer[2][j + 1] = 0;
      buffer[2][j + 2] = 0;
      buffer[2][j + 3] = 0;
    }
  }
}
