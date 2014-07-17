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
		

	cells_x = 2 * problem->pin_x*problem->refinement / problem->num_cellsets[0];
	cells_y = 2 * problem->pin_y*problem->refinement / problem->num_cellsets[1];
	cells_z = problem->z_planes / problem->num_cellsets[2];

	// This needs to be the biggest angleset
	angle_per_angleset = problem->quad.Anglesets[0].angle_per_angleset;
	group_per_groupset = problem->group_per_groupset;

	SetBoundaryConditions(problem);

	X_buffer.resize(cells_y*cells_z*group_per_groupset*angle_per_angleset*4);
	Y_buffer.resize(cells_x*cells_z*group_per_groupset*angle_per_angleset*4);
	Z_buffer.resize(cells_x*cells_y*group_per_groupset*angle_per_angleset*4);
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
	if (problem->bcs == 1)
	{
		std::cout << "Reflecting boundary conditions not implemented yet" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
	{
		bc.resize(3);
		for (int i = 0; i < bc.size(); i++)
		{
			bc[i] = 7. / 3.;
			//bc[i] = 1.0;
		}
	}

}
double Subdomain::GetBoundaryCondition(int Boundary)
{
	return bc[Boundary];
}
void Subdomain::AllocateSendBuffers(int num_tasks)
{
	X_Send_buffer.resize(cells_y*cells_z*group_per_groupset*angle_per_angleset * 4*num_tasks);
	Y_Send_buffer.resize(cells_x*cells_z*group_per_groupset*angle_per_angleset * 4*num_tasks);
	Z_Send_buffer.resize(cells_x*cells_y*group_per_groupset*angle_per_angleset * 4*num_tasks);
}
void Subdomain::Set_buffer(int cell_x, int cell_y, int cell_z, int group, int angle, int face, int task, vector<double>& RHS)
{
	if (face == 0 || face == 1)
	{
		for (int i = 0; i < 4; i++)
		{
			X_Send_buffer[task*X_buffer.size() + cell_y*cells_z*group_per_groupset*angle_per_angleset * 4 + cell_z*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4 + i] = RHS[i];
		}
	}
	else if (face == 2 || face == 3)
	{
		for (int i = 0; i < 4; i++)
		{
			Y_Send_buffer[task*Y_buffer.size() + cell_x*cells_z*group_per_groupset*angle_per_angleset * 4 + cell_z*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4 + i] = RHS[i];
		}
	}
	else if (face == 4 || face == 5)
	{
		for (int i = 0; i < 4; i++)
		{
			Z_Send_buffer[task*Z_buffer.size() + cell_x*cells_y*group_per_groupset*angle_per_angleset * 4 + cell_y*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4 + i] = RHS[i];
		}
	}
}
void Subdomain::Get_buffer(int cell_x, int cell_y, int cell_z, int group, int angle, int face, vector<double>& RHS)
{
	if (face == 0 || face == 1)
	{
		for (int i = 0; i < 4; i++)
		{
			RHS[i] = X_buffer[cell_y*cells_z*group_per_groupset*angle_per_angleset * 4 + cell_z*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4 + i];
		}
	}
	else if (face == 2 || face == 3)
	{
		for (int i = 0; i < 4; i++)
		{
			RHS[i] = Y_buffer[cell_x*cells_z*group_per_groupset*angle_per_angleset * 4 + cell_z*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4 + i];
		}
	}
	else
	{
		for (int i = 0; i < 4; i++)
		{
			RHS[i] = Z_buffer[cell_x*cells_y*group_per_groupset*angle_per_angleset * 4 + cell_y*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4 + i];
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
			X_buffer[j] = GetBoundaryCondition(0);
			X_buffer[j + 1] = 0;
			X_buffer[j + 2] = 0;
			X_buffer[j + 3] = 0;
			
		}
	}
	else if (face == 2 || face == 3)
	{
		// Need to set the cell averages, all the slopes will be zero
		int size = cells_x*cells_z*group_per_groupset*angle_per_angleset;
		for (int i = 0; i < size; i++)
		{
			int j = 4 * i;
			Y_buffer[j] = GetBoundaryCondition(1);
			Y_buffer[j + 1] = 0;
			Y_buffer[j + 2] = 0;
			Y_buffer[j + 3] = 0;
		}
	}
	else
	{
		int size = cells_x*cells_y*group_per_groupset*angle_per_angleset;
		for (int i = 0; i < size; i++)
		{
			int j = 4 * i;
			Z_buffer[j] = GetBoundaryCondition(0);
			Z_buffer[j + 1] = 0;
			Z_buffer[j + 2] = 0;
			Z_buffer[j + 3] = 0;
		}
	}
}
