#include "CellSet.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>

class Problem;

CellSet::CellSet()
{ }

CellSet::~CellSet()
{ }

void CellSet::BuildCellSet(int CS_ID, Problem* problem)
{
  num_cellsets.resize(3);
  for (int i = 0; i < num_cellsets.size(); i++)
  {
    num_cellsets[i] = problem->num_cellsets[i];
  }
  cells_x = 2*problem->num_pin_x*problem->refinement/num_cellsets[0];
  cells_y = 2*problem->num_pin_y*problem->refinement/num_cellsets[1];
  cells_z = problem->z_planes/num_cellsets[2];
  
  cells_xy = cells_x * cells_y;

  num_groupset = problem->num_groupsets;
  num_angleset = problem->quad.num_angleset;

  cells_per_cellset = cells_x*cells_y*cells_z;
  
  Cells.resize(cells_per_cellset);
  for(int i=0; i<cells_per_cellset; i++){
    Cells[i].BuildCell(i, CS_ID, problem);
  }

  // Set global boundary flags
  globalboundary.resize(6);
  BoundaryFlux.resize(6);
  SetGlobalBoundary(CS_ID, problem);
  

  // Figure out the neighbor cellsets (or the boundary this 
  // cellset is on) and the SML that owns that cellset
  // If on a global boundary, the SML will be the negative of
  // this cellsets SML
  neighbors.resize(6);
  GetNeighbors(CS_ID);
  GetNeighborsSML(CS_ID, problem);
}

void CellSet::GetCellSetijk(int CS_ID, int Dx, int Dy, int Dz, std::vector<int>& ijk)
{
  ijk[2] = (int)(CS_ID/(Dx*Dy));
  ijk[1] = (int)((CS_ID-ijk[2]*Dx*Dy)/Dx);
  ijk[0] = CS_ID-ijk[1]*Dx - ijk[2]*Dx*Dy;
}

int CellSet::GetDepth(int CS_ID, int octant, Problem* problem)
{
  int depth = -1;
  // Get the ijk location of each cellset
  std::vector<int> ijk;
  ijk.resize(3);
  GetCellSetijk(CS_ID, problem->num_cellsets[0], problem->num_cellsets[1],
    problem->num_cellsets[2], ijk);

  // Get the depth to the corner element of each octant in cellset space
  if(octant == 0)
  {
    depth = (problem->num_cellsets[0] - ijk[0] - 1) +
      (problem->num_cellsets[1] - ijk[1] - 1) +
      (problem->num_cellsets[2] - ijk[2] - 1);
  }
  else if (octant == 1)
  {
    depth = ijk[0] + (problem->num_cellsets[1] - ijk[1] - 1) +
      (problem->num_cellsets[2] - ijk[2] - 1);
  }
  else if (octant == 2)
  {
    depth = ijk[0] + ijk[1] +
      (problem->num_cellsets[2] - ijk[2] - 1);
  }
  else if (octant == 3)
  {
    depth = (problem->num_cellsets[0] - ijk[0]) + (ijk[1] - 1) +
      (problem->num_cellsets[2] - ijk[2] - 1);
  }
  else if (octant == 4)
  {
    depth = (problem->num_cellsets[0] - ijk[0] - 1) +
      (problem->num_cellsets[1] - ijk[1] - 1) + ijk[2];
  }
  else if (octant == 5)
  {
    depth = ijk[0] + (problem->num_cellsets[1] - ijk[1] - 1) +
      ijk[2];
  }
  else if (octant == 6)
  {
    depth = ijk[0] + ijk[1] + ijk[2];
  }
  else if (octant == 7)
  {
    depth = (problem->num_cellsets[0] - ijk[0]) + (ijk[1] - 1) +
      ijk[2];
  }
  else
  {
    std::cout << "Invalid Octant" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  return depth;

}

void CellSet::SetGlobalBoundary(int CS_ID, Problem* problem)
{
  // Get the globalijk location of each cellset
  std::vector<int> globalijk;
  globalijk.resize(3);
  GetCellSetijk(CS_ID, problem->num_cellsets[0], problem->num_cellsets[1],
    problem->num_cellsets[2], globalijk);


  // boundary[0] = 1 on the -y boundary
  // boundary[0] = 0 on the interior
  if (globalijk[1] == 0)
  {
    globalboundary[0] = 1;
  }
  else
    globalboundary[0] = 0;

  // boundary[1] = 1 on the +x boundary
  // boundary[1] = 0 on the interior
  if (globalijk[0] == problem->num_cellsets[0] - 1)
  {
    globalboundary[1] = 1;
  }
  else
    globalboundary[1] = 0;

  // boundary[2] = 1 on the y boundary
  // boundary[2] = 0 on the interior
  if (globalijk[1] == problem->num_cellsets[1] - 1)
  {
    globalboundary[2] = 1;
  }
  else
    globalboundary[2] = 0;

  // boundary[3] = 1 on the -x boundary
  // boundary[3] = 0 on the interior
  if (globalijk[0] == 0)
  {
    globalboundary[3] = 1;
  }
  else
    globalboundary[3] = 0;

  // boundary[4] = 1 on the -z boundary
  // boundary[4] = 0 on the interior
  if (globalijk[2] == 0)
  {
    globalboundary[4] = 1;
  }
  else
    globalboundary[4] = 0;

  // boundary[5] = 1 on the z boundary
  // boundary[5] = 0 on the interior
  if (globalijk[2] == problem->num_cellsets[2] - 1)
  {
    globalboundary[5] = 1;
  }
  else
    globalboundary[5] = 0;

}

void CellSet::GetNeighbors(int CS_ID)
{
  // Neighbor Directions
  // For arbitrary grids, these will be the face normals
  // Since we're only allowing block grids, the faces
  // are orthogonal along the xyz axes
  for (int i = 0; i < neighbors.size(); i++)
    neighbors[i].direction.resize(3.0);


  neighbors[0].direction[1] = -1;
  neighbors[1].direction[0] = 1;
  neighbors[2].direction[1] = 1;
  neighbors[3].direction[0] = -1;
  neighbors[4].direction[2] = -1;
  neighbors[5].direction[2] = 1;

  // If ID is <0, the CellSet is on a global boundary
  // id = -1 is on -y boundary, id = -2 is on the +x boundary
  // id = -3 is on +y boundary, id = -4 is on the -x boundary
  // id = -5 is on -z boundary, id = -6 is on the +z boundary

  // -Y neighbor
  if (globalboundary[0] == 1)
    neighbors[0].id = -1;
  // Interior
  else
    neighbors[0].id = CS_ID - num_cellsets[0];
    
  // +X neighbor
  if (globalboundary[1] == 1)
    neighbors[1].id = -2;
  // Interior
  else
    neighbors[1].id = CS_ID + 1;

  // +Y neighbor
  if (globalboundary[2] == 1)
    neighbors[2].id = -3;
  // Interior
  else
    neighbors[2].id = CS_ID + num_cellsets[0];
    
  // -X neighbor
  if (globalboundary[3] == 1)
    neighbors[3].id = -4;
  // Interior
  else
    neighbors[3].id = CS_ID - 1;

  // -Z neighbor
  if (globalboundary[4] == 1)
    neighbors[4].id = -5;
  // Interior
  else
    neighbors[4].id = CS_ID - num_cellsets[0] * num_cellsets[1];

  // +Z neighbor
  if (globalboundary[5] == 1)
    neighbors[5].id = -6;
  // Interior
  else
    neighbors[5].id = CS_ID + num_cellsets[0] * num_cellsets[1];


}

void CellSet::GetNeighborsSML(int CS_ID, Problem* problem)
{
  // ijk location of neighbor cellset
  std::vector<int> neighborijk;
  neighborijk.resize(3);

  // Processor ijk of cellset
  std::vector<int> p;
  p.resize(3);
  
  for(int i = 0; i < neighbors.size(); i++)
  {
    // If we're on a global boundary, set the neighbor SML to negative
    if (neighbors[i].id < 0)
    {
      neighbors[i].SML = -1;
    }
    else
    {
      // Find the ijk location of the neighbor CellSet
      GetCellSetijk(neighbors[i].id, problem->num_cellsets[0], problem->num_cellsets[1],
        problem->num_cellsets[2], neighborijk);

      // Using the partition functions, compute the SML_ID
      // that owns that cellset
      for(int j = 0; j < 3; j++)
      {
        // Round Robin
        if(problem->partition_function[j] == 1)
          p[j] = (int)floor(neighborijk[j]/problem->overload[j]);
        // Blocked
        else
          p[j] = neighborijk[j] % (num_cellsets[j]/problem->overload[j]);
      }

      neighbors[i].SML = p[0] + num_cellsets[0]/problem->overload[0]*p[1] +
        (num_cellsets[0]/problem->overload[0])*(num_cellsets[1]/problem->overload[1])*p[2];
    }
    

  }
}
