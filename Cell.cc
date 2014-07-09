#include "Cell.h"
#include <vector>
#include <iostream>
#include <math.h>

using std::vector;


Cell::Cell()
{ }

Cell::~Cell()
{ }

void Cell::BuildCell(int Local_ID, int CS_ID, Problem* problem)
{
	// Compute the Global and Local Cell ID's
	ComputeCellID(Local_ID,CS_ID, problem);

	// Set delta_x, delta_y, and delta_z from problem
	// pin pitch = 12.0, total height = 350 cm
	delta_x = 4.0/(problem->refinement*2);
	delta_y = 4.0/(problem->refinement*2);
	delta_z = 4.0/problem->z_planes;

	// Set cross section and source
	sigma_t = 2.0;
	source = 1.0;

	// Get the face normals and centers
	GetFaceNormals();
	GetFaceCenters();

	// Set CellSet Boundary Flags
	localboundary.resize(3);
	SetLocalBoundary();

	// Figure out global Cell_ID of 6 neighbors
	num_cellsets.resize(3);
	num_cellsets = problem->num_cellsets;
	neighbors.resize(6);
	GetNeighbors(CS_ID);

	// Build the DFEM Matrices
	ComputeDFEMMatrices();

	// Resize the phi vector
	phi.resize(problem->num_groupsets*problem->group_per_groupset);

}

void Cell::SetLocalBoundary()
{
	// boundary[0] = -1 on the -x boundary
	// boundary[0] =  1 on the +x boundary
	// boundary[0] =  0 on the interior
	if(localijk[0] == 0)
		localboundary[0] = -1;
	else if(localijk[0] == cells_x -1)
		localboundary[0] = 1;
	else
		localboundary[0] = 0;

	// boundary[1] = -1 on the -y boundary
	// boundary[1] =  1 on the +y boundary
	// boundary[1] =  0 on the interior
	if(localijk[1] == 0)
		localboundary[1] = -1;
	else if(localijk[1] == cells_y -1)
		localboundary[1] = 1;
	else
		localboundary[1] = 0;

	// boundary[2] = -1 on the -z boundary
	// boundary[2] =  1 on the +z boundary
	// boundary[2] =  0 on the interior
	if(localijk[2] == 0)
		localboundary[2] = -1;
	else if(localijk[2] == cells_z -1)
		localboundary[2] = 1;
	else
		localboundary[2] = 0;
}

void Cell::GetCellijk(int Cell_ID, int Dx, int Dy, int Dz, std::vector<int>& ijk)
{
	ijk[2] = (int)(Cell_ID/(Dx*Dy));
	ijk[1] = (int)((Cell_ID-ijk[2]*Dx*Dy)/Dx);
	ijk[0] = Cell_ID-ijk[1]*Dx - ijk[2]*Dx*Dy;
}

void Cell::ComputeCellID(int Local_ID, int CS_ID, Problem* problem)
{
	cells_x = 2*problem->pin_x*problem->refinement/problem->num_cellsets[0];
	cells_y = 2*problem->pin_y*problem->refinement/problem->num_cellsets[1];
	cells_z = problem->z_planes/problem->num_cellsets[2];

	//std::cout << "cells(x,y,z)= " << cells_x << " " << cells_y << " " << cells_z << std::endl;

	cells_per_cellset = cells_x*cells_y*cells_z;

	//std::cout << "cells_per_cellset= " << cells_per_cellset << std::endl;

	CellID = CS_ID*cells_per_cellset + Local_ID;
	LocalCellID = Local_ID;

	localijk.resize(3);
	GetCellijk(LocalCellID, cells_x, cells_y, cells_z, localijk);
	//std::cout << "Local_ijk: " << localijk[0] << " " << localijk[1] << " " << localijk[2] << std::endl;

	//globalijk.resize(3);
	//// First recomputed Cellsets ijk and store it in global ijk
	//GetCellijk(CS_ID, problem->num_cellsets[0], problem->cellsets[1], 
	//	problem->num_cellsets[2], globalijk);

	//// Now compute cells Global ijk and overwrite what's stored in there
	//globalijk[0] = globalijk[0]*cells_x + localijk[0];
	//globalijk[1] = globalijk[1]*cells_y + localijk[1];
	//globalijk[2] = globalijk[2]*cells_z + localijk[2];

	//std::cout << "Global_ijk: " << globalijk[0] << " " << globalijk[1] << " " << globalijk[2] << std::endl;
}

void Cell::GetNeighbors(int CS_ID)
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
	neighbors[4].direction.z = 1;
	neighbors[5].direction.z = -1;

	
	// If ID is <0, the cell is on the cellset boundary
	// id = -1 is on +x boundary, id = -2 is on the -x boundary
	// id = -3 is on +y boundary, id = -4 is on the -y boundary
	// id = -5 is on +z boundary, id = -6 is on the -z boundary

	// +X neighbor
	if(localboundary[0] == 1)
	{
		neighbors[0].id = CellID + cells_per_cellset - (cells_x - 1);
	}
	// Interior
	else
	{
		neighbors[0].id = CellID + 1;
	}

	// -X neighbor
	if(localboundary[0] == -1)
	{
		neighbors[1].id = CellID - cells_per_cellset + (cells_x - 1);
	}
	// Interior
	else
	{
		neighbors[1].id = CellID - 1;
	}

	// +Y neighbor
	if(localboundary[1] == 1)
	{
		neighbors[2].id = CellID + cells_per_cellset*num_cellsets[0] - localijk[1]*cells_x;
	}
	// Interior
	else
	{
		neighbors[2].id = CellID + cells_x;
	}

	// -Y neighbor
	if(localboundary[1] == -1)
	{
		neighbors[3].id = CellID - cells_per_cellset*num_cellsets[0] + (cells_y - localijk[1] - 1)*cells_x;
	}
	// Interior
	else
	{
		neighbors[3].id = CellID - cells_x;
	}

	// +Z neighbor
	if(localboundary[2] == 1)
	{
		neighbors[4].id = CellID + cells_per_cellset*num_cellsets[0]*num_cellsets[1] - localijk[2]*cells_x*cells_y;
	}
	// Interior
	else
	{
		neighbors[4].id = CellID + cells_x*cells_y;
	}

	// -Z neighbor
	if(localboundary[2] == -1)
	{
		// This needs to be checked....
		neighbors[5].id = CellID - cells_per_cellset*num_cellsets[0]*num_cellsets[1] + (cells_z - localijk[2])*cells_x*cells_y;
	}
	// Interior
	else
	{
		neighbors[5].id = CellID - cells_x*cells_y;
	}
		
}

void Cell::GetFaceNormals()
{
	normals.resize(6);

	for (int i = 0; i < normals.size(); i++)
	{
		normals[i].x = 0;
		normals[i].y = 0;
		normals[i].z = 0;
	}
	normals[0].x = 1;
	normals[1].x = -1;
	normals[2].y = 1;
	normals[3].y = -1;
	normals[4].z = 1;
	normals[5].z = -1;
}

void Cell::GetFaceCenters()
{
	facecenters.resize(6);

	for (int i = 0; i < facecenters.size(); i++)
	{
		facecenters[i].x = 0;
		facecenters[i].y = 0;
		facecenters[i].z = 0;
	}
	facecenters[0].x = delta_x/2;
	facecenters[1].x = -delta_x/2;
	facecenters[2].y = delta_y/2;
	facecenters[3].y = -delta_y/2;
	facecenters[4].z = delta_z/2;
	facecenters[5].z = -delta_z/2;
}

void Cell::ComputeDFEMMatrices()
{
	// M is the mass matrix
	M.resize(4);
	M[0] = 1;
	M[1] = (1/12)*delta_x*delta_y*delta_z;
	M[2] = (1/12)*delta_x*delta_y*delta_z;
	M[3] = (1/12)*delta_x*delta_y*delta_z;

	// N is the surface matrix. The first dimension is the face
	// Each face as a 4x4 matrix of (x,y,z) components
	N.resize(6, vector<vector<Direction > >(4, vector< Direction >(4) ) );
	// Face 0: e_n = (1, 0, 0)
	N[0][0][0].x = delta_y*delta_z;
	N[0][0][1].x = delta_y*delta_z;
	N[0][1][0].x = delta_y*delta_z;
	N[0][1][1].x = delta_y*delta_z;
	N[0][2][2].x = delta_y*delta_z/12;
	N[0][3][3].x = delta_y*delta_z/12;

	// Face 1: e_n = (-1, 0, 0)
	for(int i=0; i< N[1].size(); i++)
	{
		for(int j=0; j<N[1].size();j++)
		{
			N[1][i][j].x = N[0][i][j].x;
			if(i == j)
				N[1][i][j].x *= -1;
		}
	}

	// Face 2: e_n = (0, 1, 0)
	N[2][0][0].y = delta_x*delta_z;
	N[2][0][1].y = delta_x*delta_z;
	N[2][1][1].y = delta_x*delta_z/12;
	N[2][2][0].y = delta_x*delta_z;
	N[2][2][2].y = delta_x*delta_z;
	N[2][3][3].y = delta_x*delta_z/12;
	
	// Face 3: e_n = (0, -1, 0)
	for(int i=0; i< N[3].size(); i++)
	{
		for(int j=0; j<N[3].size();j++)
		{
			N[3][i][j].y = N[2][i][j].y;
			if(i == j)
				N[3][i][j].y *= -1;
		}
	}

	// Face 4: e_n = (0, 0, 1)
	N[4][0][0].z = delta_x*delta_y;
	N[4][0][1].z = delta_x*delta_y;
	N[4][1][1].z = delta_x*delta_y/12;
	N[4][2][2].z = delta_x*delta_y/12;
	N[4][3][0].z = delta_x*delta_y;
	N[4][3][3].z = delta_x*delta_y;

	// Face 5: e_n = (0, 0, -1)
	for(int i=0; i< N[5].size(); i++)
	{
		for(int j=0; j<N[5].size();j++)
		{
			N[5][i][j].z = N[4][i][j].z;
			if(i == j)
				N[5][i][j].z *= -1;
		}
	}

	// L is the gradient matrix
	// It's mostly zeros, but it's fully built to
	// simulate non-rectangular geometries
	// It is a 4x4 matrix of (x,y,z) components
	L.resize(4, vector< Direction >( 4 ) );
	L[0][1].x = delta_y*delta_z;
	L[0][2].y = delta_x*delta_z;
	L[0][3].z = delta_x*delta_y;


}

