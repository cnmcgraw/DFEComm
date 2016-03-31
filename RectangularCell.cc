#include "RectangularCell.h"
#include <vector>
#include <iostream>
#include <math.h>

using std::vector;


RectangularCell::RectangularCell()
{ }

RectangularCell::~RectangularCell()
{ }

void RectangularCell::BuildCell(int Local_ID, int CS_ID, Problem* problem)
{
  // Compute the Global and Local Cell ID's
  ComputeCellID(Local_ID,CS_ID, problem);

  // Set delta_x, delta_y, and delta_z from problem
  // pin pitch = 4.0, total height = 4 cm
  delta_x = 4.0/(problem->refinement*2.);
  delta_y = 4.0/(problem->refinement*2.);
  delta_z = 4.0/problem->z_planes;

  // Set cross section and source
  sigma_t = 3.0;
  source = 7.0;

  // Get the centroid and vertices of the cell
  GetVertices(CellID);

  // Get the face normals and centers
  GetFaceNormals();
  GetFaceCenters();

  // Set CellSet Boundary Flags
  localboundary.resize(6,0);
  SetLocalBoundary();

  // Figure out global Cell_ID of 6 neighbors
  num_cellsets.resize(3);
  num_cellsets = problem->num_cellsets;
  num_faces = 6;
  max_faces = problem->max_faces;
  neighbors.resize(num_faces);
  GetNeighbors(CS_ID);

  // Build the DFEM Matrices
  ComputeDFEMMatrices();

  // Resize the phi vector
  phi.resize(problem->num_groupsets*problem->group_per_groupset*4);

}

void RectangularCell::SetLocalBoundary()
{
  // Localboundary is filled with zeros. If we're on a boundary, we replace the zero with a 1
  int x,y,z;
  z = localijk[1];
  y = int(localijk[0]/cells_x);
  x = localijk[0] - cells_x * y;
  
  if(y == 0)
    localboundary[0] = 1;
  if(x == cells_x -1)
    localboundary[1] = 1;
  if(y == cells_y - 1)
    localboundary[2] = 1;
  if(x == 0)
    localboundary[3] = 1;
  if(z == 0)
    localboundary[4] = 1;
  if(z == cells_z - 1)
    localboundary[5] = 1;
}

void RectangularCell::GetCellijk(int Cell_ID, int Dx, int Dy, int Dz, std::vector<int>& ijk)
{
  
 // ijk[2] = (int)(Cell_ID/(Dx*Dy));
 // ijk[1] = (int)((Cell_ID-ijk[2]*Dx*Dy)/Dx);
 // ijk[0] = Cell_ID-ijk[1]*Dx - ijk[2]*Dx*Dy;
  
  // This is how we'd calculate z and the xy components
  ijk[1] = (int)(Cell_ID/(Dx*Dy));
  int y = (int)((Cell_ID-ijk[1]*Dx*Dy)/Dx);
  int x = Cell_ID-y*Dx - ijk[1]*Dx*Dy;
  
  ijk[0] = x + Dx*y;
}

void RectangularCell::ComputeCellID(int Local_ID, int CS_ID, Problem* problem)
{
  // These are the number of cells, in each dimension, for the cellset
  cells_x = 2*problem->num_pin_x*problem->refinement/problem->num_cellsets[0];
  cells_y = 2*problem->num_pin_y*problem->refinement/problem->num_cellsets[1];
  cells_z = problem->z_planes/problem->num_cellsets[2];

  cells_per_cellset = cells_x*cells_y*cells_z;

  CellID = CS_ID*cells_per_cellset + Local_ID;
  LocalCellID = Local_ID;

  localijk.resize(2);
  GetCellijk(LocalCellID, cells_x, cells_y, cells_z, localijk);
}

void RectangularCell::GetNeighbors(int CS_ID)
{
  // Neighbor Directions
  // For arbitrary grids, these will be the face normals
  // Since we're only allowing block grids, the faces
  // are orthogonal along the xyz axes
  for (int i = 0; i < neighbors.size(); i++)
    neighbors[i].direction.resize(3,0);

    
  neighbors[0].direction[1] = -1; //-y
  neighbors[1].direction[0] = 1; //+x
  neighbors[2].direction[1] = 1; //+y
  neighbors[3].direction[0] = -1; //-x
  neighbors[4].direction[2] = -1; //-z
  neighbors[5].direction[2] = 1; //+z
  
  // This is the face id for the neighboring face.
  // We can compute this for spiderweb grids, it's just more complicated
  neighbors[0].face = 2; //-y
  neighbors[1].face = 3; //+x
  neighbors[2].face = 0; //+y
  neighbors[3].face = 1; //-x
  neighbors[4].face = 5; //+z
  neighbors[5].face = 4; //-z


    // -Y neighbor
  if(localboundary[0] == 1)
  {
   // neighbors[0].id = CellID - cells_per_cellset*num_cellsets[0] + (cells_y - localijk[1] - 1)*cells_x;
    neighbors[0].id = LocalCellID;
    neighbors[0].cs = 1;
  }
  // Interior
  else
  {
    neighbors[0].id = LocalCellID - cells_x;
    neighbors[0].cs = 0;
  }
  
  // +X neighbor
  if(localboundary[1] == 1)
  {
   // neighbors[1].id = CellID + cells_per_cellset - (cells_x - 1);
    neighbors[1].id = LocalCellID;
    neighbors[1].cs = 1;
  }
  // Interior
  else
  {
    neighbors[1].id = LocalCellID + 1;
    neighbors[1].cs = 0;
  }

  // +Y neighbor
  if(localboundary[2] == 1)
  {
   // neighbors[2].id = CellID + cells_per_cellset*num_cellsets[0] - localijk[1]*cells_x;
    neighbors[2].id = LocalCellID;
    neighbors[2].cs = 1;
  }
  // Interior
  else
  {
    neighbors[2].id = LocalCellID + cells_x;
    neighbors[2].cs = 0;
  }

  // -X neighbor
  if(localboundary[3] == 1)
  {
   // neighbors[3].id = CellID - cells_per_cellset + (cells_x - 1);
    neighbors[3].id = LocalCellID;
    neighbors[3].cs = 1;
  }
  // Interior
  else
  {
    neighbors[3].id = LocalCellID - 1;
    neighbors[3].cs = 0;
  }

  // -Z neighbor
  if(localboundary[4] == 1)
  {
   // neighbors[4].id = CellID - cells_per_cellset*num_cellsets[0]*num_cellsets[1] + (cells_z - localijk[2])*cells_x*cells_y;
    neighbors[4].id = LocalCellID;
    neighbors[4].cs = 1;
  }
  // Interior
  else
  {
    neighbors[4].id = LocalCellID - cells_x*cells_y;
    neighbors[4].cs = 0;
  }
  
  // +Z neighbor
  if(localboundary[5] == 1)
  {
   // neighbors[5].id = CellID + cells_per_cellset*num_cellsets[0]*num_cellsets[1] - localijk[2]*cells_x*cells_y;
    neighbors[5].id = LocalCellID;
    neighbors[5].cs = 1;
  }
  // Interior
  else
  {
    neighbors[5].id = LocalCellID + cells_x*cells_y;
    neighbors[5].cs = 0;
  }
}

void RectangularCell::GetCellInOut(std::vector<double> omega, std::vector<int>& incoming, std::vector<std::vector<int> >& outgoing)
{

    int r(0), s(0);
    for(int f = 0; f < num_faces; f++)
    {
      if (dot(omega, neighbors[f].direction) < 0)
      {
        incoming[r] = f;
        r += 1;
      }
      else
      {
       // If outgoing neighbor is in this cell-set, we'll compute the neighbors id for the face
       // If it is in another cell-set, we set the outgoing face to this cell's id for that face
        if(neighbors[f].cs == 0)
        {
          outgoing[s][0] = neighbors[f].id;
          outgoing[s][1] = f;
          outgoing[s][2] = neighbors[f].face;
          s += 1;
        }
        else
        {
          outgoing[s][0] = LocalCellID;
          outgoing[s][1] = f;
          outgoing[s][2] = f;
          s += 1; 
        }        
      }
    }

}

void RectangularCell::GetFaceNormals()
{
  normals.resize(6, std::vector<double>(3));

  for (int i = 0; i < normals.size(); i++)
  {
    normals[i][0] = 0;
    normals[i][1] = 0;
    normals[i][2] = 0;
  }
  normals[0][1] = -1.;
  normals[1][0] = 1.;
  normals[2][1] = 1.;
  normals[3][0] = -1.;
  normals[4][2] = -1.;
  normals[5][2] = 1.;
}

void RectangularCell::GetFaceCenters()
{
  facecenters.resize(6, std::vector<double>(3));
  face_ids.resize(6, 0);

  for (int i = 0; i < facecenters.size(); i++)
  {
    facecenters[i][0] = 0;
    facecenters[i][1] = 0;
    facecenters[i][2] = 0;
    face_ids[i] = i;
  }
  facecenters[0][1] = -delta_y/2.;
  facecenters[1][0] = delta_x/2.;
  facecenters[2][1] = delta_y/2.;
  facecenters[3][0] = -delta_x/2.;
  facecenters[4][2] = -delta_z/2.;
  facecenters[5][2] = delta_z/2.;
}

void RectangularCell::GetVertices(int CellID)
{
  // First we get the global ijk of the cell using the global CellID
  std::vector<int> cell_ijk(3,0), tmp_ijk(2,0);
  GetCellijk(CellID, cells_x, cells_y, cells_z, tmp_ijk);
  
  cell_ijk[1] = int(tmp_ijk[0]/cells_x);
  cell_ijk[0] = tmp_ijk[0] - cells_x * cell_ijk[1];
  cell_ijk[2] = tmp_ijk[1];
  
  vertices.resize(8, std::vector<double>(3));
  centroid.resize(3, 0.0);
  
  // Centroid of the cell
  centroid[0] = (cell_ijk[0] + 0.5)*delta_x;
  centroid[1] = (cell_ijk[1] + 0.5)*delta_y;
  centroid[2] = (cell_ijk[2] + 0.5)*delta_z;
  

  vertices[0][0] = centroid[0] - delta_x/2;
  vertices[0][1] = centroid[1] - delta_y/2;
  vertices[0][2] = centroid[2] - delta_z/2;
  
  vertices[1][0] = centroid[0] + delta_x/2;
  vertices[1][1] = centroid[1] - delta_y/2;
  vertices[1][2] = centroid[2] - delta_z/2;
  
  vertices[2][0] = centroid[0] + delta_x/2;
  vertices[2][1] = centroid[1] + delta_y/2;
  vertices[2][2] = centroid[2] - delta_z/2;

  vertices[3][0] = centroid[0] - delta_x/2;
  vertices[3][1] = centroid[1] + delta_y/2;
  vertices[3][2] = centroid[2] - delta_z/2;
  
  vertices[4][0] = centroid[0] - delta_x/2;
  vertices[4][1] = centroid[1] - delta_y/2;
  vertices[4][2] = centroid[2] + delta_z/2;
  
  vertices[5][0] = centroid[0] + delta_x/2;
  vertices[5][1] = centroid[1] - delta_y/2;
  vertices[5][2] = centroid[2] + delta_z/2;
  
  vertices[6][0] = centroid[0] + delta_x/2;
  vertices[6][1] = centroid[1] + delta_y/2;
  vertices[6][2] = centroid[2] + delta_z/2;

  vertices[7][0] = centroid[0] - delta_x/2;
  vertices[7][1] = centroid[1] + delta_y/2;
  vertices[7][2] = centroid[2] + delta_z/2;
}

void RectangularCell::ComputeDFEMMatrices()
{
  // M is the mass matrix
  M.resize(4, vector< double >( 4 ) );
  M[0][0] = delta_x*delta_y*delta_z;
  M[1][1] = (1./12.)*delta_x*delta_y*delta_z;
  M[2][2] = (1./12.)*delta_x*delta_y*delta_z;
  M[3][3] = (1./12.)*delta_x*delta_y*delta_z;

  // N is the surface matrix. The first dimension is the face
  // Each face as a 4x4 matrix of (x,y,z) components
  N.resize(6, vector<vector<vector<double > > >(4, vector< vector<double> >(4, vector<double>(3) ) ) );
  // Face 1: e_n = (1, 0, 0)
  N[1][0][0][0] = delta_y*delta_z;
  N[1][0][1][0] = delta_y*delta_z;
  N[1][1][0][0] = delta_y*delta_z;
  N[1][1][1][0] = delta_y*delta_z;
  N[1][2][2][0] = delta_y*delta_z/12.;
  N[1][3][3][0] = delta_y*delta_z/12.;


  // Face 0: e_n = (0, -1, 0)
  N[0][0][0][1] = -delta_x*delta_z;
  N[0][0][2][1] = delta_x*delta_z;
  N[0][1][1][1] = -delta_x*delta_z/12.;
  N[0][2][0][1] = delta_x*delta_z;
  N[0][2][2][1] = -delta_x*delta_z;
  N[0][3][3][1] = delta_x*delta_z/12.;
  
  // Face 3: e_n = (-1, 0, 0)
  for(int i=0; i< N[3].size(); i++)
  {
    for(int j=0; j<N[3].size();j++)
    {
      N[3][i][j][0] = N[1][i][j][0];
      if(i == j)
        N[3][i][j][0] *= -1.;
    }
  }
  
  // Face 2: e_n = (0, -1, 0)
  for(int i=0; i< N[2].size(); i++)
  {
    for(int j=0; j<N[2].size();j++)
    {
      N[2][i][j][1] = N[0][i][j][1];
      if(i == j)
        N[2][i][j][1] *= -1.;
    }
  }

  // Face 4: e_n = (0, 0, -1)
  N[4][0][0][2] = -delta_x*delta_y;
  N[4][0][3][2] = delta_x*delta_y;
  N[4][1][1][2] = -delta_x*delta_y/12.;
  N[4][2][2][2] = -delta_x*delta_y/12.;
  N[4][3][0][2] = delta_x*delta_y;
  N[4][3][3][2] = -delta_x*delta_y;

  // Face 5: e_n = (0, 0, 1)
  for(int i=0; i< N[5].size(); i++)
  {
    for(int j=0; j<N[5].size();j++)
    {
      N[5][i][j][2] = N[4][i][j][2];
      if(i == j)
        N[5][i][j][2] *= -1.;
    }
  }

  // L is the gradient matrix
  // It's mostly zeros, but it's fully built to
  // simulate non-rectangular geometries
  // It is a 4x4 matrix of (x,y,z) components
  L.resize(4, vector< vector<double> >( 4, vector<double>(3) ) );
  L[0][1][0] = delta_y*delta_z;
  L[0][2][1] = delta_x*delta_z;
  L[0][3][2] = delta_x*delta_y;


}
