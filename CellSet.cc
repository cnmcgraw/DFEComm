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

  int size = problem->num_groupsets*problem->group_per_groupset*(problem->num_polar*problem->num_azim) * 4;

  // boundary[0] = 1 on the x boundary
  // boundary[0] = 0 on the interior
  if (globalijk[0] == problem->num_cellsets[0] - 1)
  {
    globalboundary[0] = 1;
    BoundaryFlux[0].resize(cells_y*cells_z*size);
  }
  else
    globalboundary[0] = 0;

  // boundary[1] = 1 on the -x boundary
  // boundary[1] = 0 on the interior
  if (globalijk[0] == 0)
  {
    globalboundary[1] = 1;
    BoundaryFlux[1].resize(cells_y*cells_z*size);
  }
  else
    globalboundary[1] = 0;

  // boundary[2] = 1 on the y boundary
  // boundary[2] = 0 on the interior
  if (globalijk[1] == problem->num_cellsets[1] - 1)
  {
    globalboundary[2] = 1;
    BoundaryFlux[2].resize(cells_x*cells_z*size);
  }
  else
    globalboundary[2] = 0;

  // boundary[3] = 1 on the -y boundary
  // boundary[3] = 0 on the interior
  if (globalijk[1] == 0)
  {
    globalboundary[3] = 1;
    BoundaryFlux[3].resize(cells_x*cells_z*size);
  }
  else
    globalboundary[3] = 0;

  // boundary[4] = 1 on the z boundary
  // boundary[4] = 0 on the interior
  if (globalijk[2] == problem->num_cellsets[2] - 1)
  {
    globalboundary[4] = 1;
    BoundaryFlux[4].resize(cells_x*cells_y*size);
  }
  else
    globalboundary[4] = 0;

  // boundary[5] = 1 on the -z boundary
  // boundary[5] = 0 on the interior
  if (globalijk[2] == 0)
  {
    globalboundary[5] = 1;
    BoundaryFlux[5].resize(cells_x*cells_y*size);
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
  {
    neighbors[i].direction.x = 0;
    neighbors[i].direction.y = 0;
    neighbors[i].direction.z = 0;
  }


  neighbors[0].direction.x = 1;
  neighbors[1].direction.x = -1;
  neighbors[2].direction.y = 1;
  neighbors[3].direction.y = -1;
  if (neighbors.size() == 6){
    neighbors[4].direction.z = 1;
    neighbors[5].direction.z = -1;
  }

  // If ID is <0, the CellSet is on a global boundary
  // id = -1 is on +x boundary, id = -2 is on the -x boundary
  // id = -3 is on +y boundary, id = -4 is on the -y boundary
  // id = -5 is on +z boundary, id = -6 is on the -z boundary

  // +X neighbor
  if (globalboundary[0] == 1)
    neighbors[0].id = -1;
  // Interior
  else
    neighbors[0].id = CS_ID + 1;

  // -X neighbor
  if (globalboundary[1] == 1)
    neighbors[1].id = -2;
  // Interior
  else
    neighbors[1].id = CS_ID - 1;

  // +Y neighbor
  if (globalboundary[2] == 1)
    neighbors[2].id = -3;
  // Interior
  else
    neighbors[2].id = CS_ID + num_cellsets[0];

  // -Y neighbor
  if (globalboundary[3] == 1)
    neighbors[3].id = -4;
  // Interior
  else
    neighbors[3].id = CS_ID - num_cellsets[0];

  // +Z neighbor
  if (globalboundary[4] == 1)
    neighbors[4].id = -5;
  // Interior
  else
    neighbors[4].id = CS_ID + num_cellsets[0] * num_cellsets[1];


  // -Z neighbor
  if (globalboundary[5] == 1)
    neighbors[5].id = -6;
  // Interior
  else
    neighbors[5].id = CS_ID - num_cellsets[0] * num_cellsets[1];


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
    // If the neighbors cellset id is my cellset id 
    //  then the same sml owns it.
    if (neighbors[i].id < 0)
    {
      // If this is a boundary, set the neighbor SML to 
      // this SML's id * -1
      MPI_Comm_rank(MPI_COMM_WORLD, &neighbors[i].SML);
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

void CellSet::SetBoundaryFlux(int boundary, int angleset, int groupset, std::vector<double>& buffer)
{
  // Ordered by groupset then angleset
  for (int i = 0; i < buffer.size(); i++)
    BoundaryFlux[boundary][groupset*num_angleset*buffer.size() + angleset*buffer.size() + i] = buffer[i];
}

