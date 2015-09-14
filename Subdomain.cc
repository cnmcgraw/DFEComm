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
  num_cellsets.resize(3);
  overload.resize(3);
  partition_function.resize(3);

  for(int i = 0; i < 3; i ++){
    num_cellsets[i] = problem->num_cellsets[i];
    overload[i] = problem->overload[i];
    partition_function[i] = problem->partition_function[i];
  }

  ComputeCellSetID(SML_ID);

  CellSets.resize(total_overload);
  for (int i = 0; i < total_overload; i++)
    CellSets[i].BuildCellSet(CellSetIDs[i], problem);
    
  cells_x = 2 * problem->num_pin_x*problem->refinement / problem->num_cellsets[0];
  cells_y = 2 * problem->num_pin_y*problem->refinement / problem->num_cellsets[1];
  cells_z = problem->z_planes / problem->num_cellsets[2];

  // This needs to be the biggest angleset
  angle_per_angleset = problem->quad.Anglesets[0].angle_per_angleset;
  group_per_groupset = problem->group_per_groupset;

}
void Subdomain::ComputeCellSetID(int SML_ID)
{
  // The total overload is the total number of 
  // cellsets each SML owns
  total_overload = 1;
  std::vector<int> P;
  P.resize(3);
  // For 2D P_eff in z = 1 and overload in z = 1
  for(int i=0; i< 3; i++)
  {
    total_overload *= overload[i];
    P[i] = num_cellsets[i]/overload[i];
  }

  CellSetIDs.resize(total_overload);
  
  // Compute the i,j,k of the initial cellset
  int k = (int)(SML_ID/(P[0]*P[1]));
  int j = (int)((SML_ID-k*P[0]*P[1])/P[0]);
  int i = SML_ID - j*P[0] - k*P[0]*P[1];

  // If the mesh is blocked in any dimension,
  // we must multiply by the overload in
  // that dimension
  if(partition_function[2] == 1)
    k = k*overload[2];
  if(partition_function[1] == 1)
    j = j*overload[1]; 
  if(partition_function[0] == 1)
    i = i*overload[0];

  std::vector<int> ID_ijk;
  ID_ijk.resize(3);
  int id = 0;
  for(int z = 0; z< overload[2]; z++)
  {
    // Sequential in blocked
    if(partition_function[2] == 1)
      ID_ijk[2] = k + z;
    // Jump by the overload factor for round robin
    else 
      ID_ijk[2] = k + z*overload[2];

    for(int y = 0; y< overload[1]; y++)
    {
      if(partition_function[1] == 1)
        ID_ijk[1] = j + y;
      else 
        ID_ijk[2] = j + y*overload[1];

      for(int x = 0; x<overload[0]; x++, id++)
      {
        if(partition_function[0] == 1)
          ID_ijk[0] = i + x;
        else 
          ID_ijk[2] = i + x*overload[0];

        CellSetIDs[id] = ID_ijk[0] + num_cellsets[0]*ID_ijk[1]
          + num_cellsets[0]*num_cellsets[1]*ID_ijk[2];
      }
    }
  }
}
int Subdomain::ComputeCellSetIndex(int CS_ID, int SML_ID)
{
  // The total overload is the total number of 
  // cellsets each SML owns
  std::vector<int> P;
  P.resize(3);
  // For 2D P_eff in z = 1 and overload in z = 1
  for(int i=0; i< 3; i++)
  {
    P[i] = num_cellsets[i]/overload[i];
  }
  
  // Compute the i,j,k of the initial cellset
  int k = (int)(SML_ID/(P[0]*P[1]));
  int j = (int)((SML_ID-k*P[0]*P[1])/P[0]);
  int i = SML_ID - j*P[0] - k*P[0]*P[1];

  // If the mesh is blocked in any dimension,
  // we must multiply by the overload in
  // that dimension
  if(partition_function[2] == 1)
    k = k*overload[2];
  if(partition_function[1] == 1)
    j = j*overload[1]; 
  if(partition_function[0] == 1)
    i = i*overload[0];

  int my_k = CS_ID / (num_cellsets[0]*num_cellsets[1]);
  int ij = CS_ID % (num_cellsets[0]*num_cellsets[1]);
  int my_j = ij / num_cellsets[0];
  int my_i = ij % num_cellsets[0];

  int index = (my_i - i) + overload[0]*(my_j - j) + overload[0]*overload[1]*(my_k - k);

  return index;

}

