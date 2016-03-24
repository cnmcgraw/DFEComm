#include "Task.h"
#include <stdio.h>
#include <iostream>
#include <vector>


Task::Task()
{ }

Task::~Task()
{ }

void Task::BuildTask(int SML_ID, Problem* problem, Subdomain* subdomain, int cs, int as, int gs)
{

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  num_cellset = problem->num_cellsets[0] * problem->num_cellsets[1] * problem->num_cellsets[2];
  num_angleset = problem->quad.num_angleset;
  num_groupset = problem->num_groupsets;

  int num_tasks = num_cellset * num_angleset * num_groupset;

  incoming.resize(3, std::vector<int>(2,0));
  outgoing.resize(3, std::vector<int>(3,0));
   // Loop through all faces and designate them incoming or outgoing:
  int r(0), s(0);
  for (int f = 0; f < 6; f++)
  {
    // Get the neighbors for each face
    Neighbor neighbor = subdomain->CellSets[cs].neighbors[f];
    if (dot(problem->quad.GetOmega(as), neighbor.direction) < 0)
    {
      incoming[r][0] = f;
      incoming[r][1] = neighbor.SML;
      r += 1;
    }
    else
    {
      outgoing[s][0] = f;
      outgoing[s][1] = neighbor.SML;
      
      if(neighbor.SML < 0)
        outgoing[s][2] = ComputeTaskID(cs, as, gs);
      else{
        int neighbor_cs_loc = subdomain->ComputeCellSetIndex(neighbor.id, neighbor.SML);
        outgoing[s][2] = ComputeTaskID(neighbor_cs_loc, as, gs);
      }
      s += 1;
    }
  } 


  cellset_id = subdomain->CellSetIDs[cs];
  cellset_id_loc = cs;
  angleset_id = as;
  groupset_id = gs;
  omega.resize(3);
  omega = problem->quad.GetOmega(as);
  octant = problem->quad.Anglesets[as].octant;
  depth = subdomain->CellSets[cs].GetDepth(subdomain->CellSetIDs[cs], octant, problem);

    
  cells_x = 2 * problem->num_pin_x*problem->refinement / problem->num_cellsets[0];
  cells_y = 2 * problem->num_pin_y*problem->refinement / problem->num_cellsets[1];
  cells_xy = cells_x * cells_y;
  cells_z = problem->z_planes / problem->num_cellsets[2];

  cell_per_cellset = cells_x * cells_y * cells_z;
  angle_per_angleset = problem->quad.Anglesets[0].angle_per_angleset;
  group_per_groupset = problem->group_per_groupset;

  SetBoundaryConditions(problem);


  task_id = ComputeTaskID(cellset_id_loc, angleset_id, groupset_id);
}

void Task::AllocateBuffers(){

  plane_data.resize(3);
  // These will be our surface data structures. We'll pull from interior_data for these
  plane_data[0].resize(cells_y*cells_z*group_per_groupset*angle_per_angleset * 4);
  plane_data[1].resize(cells_x*cells_z*group_per_groupset*angle_per_angleset * 4);
  plane_data[2].resize(cells_xy*group_per_groupset*angle_per_angleset * 4);
  

}

int Task::ComputeTaskID(int cs, int as, int gs){
  
  int task_id = cs * (num_groupset * num_angleset) + gs * num_angleset + as;

  return task_id;
}

void Task::SetBoundaryConditions(Problem* problem)
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


double Task::GetBoundaryCondition(int Boundary)
{
  return bc[Boundary];
}


void Task::GetInteriorData(std::vector<double >& interior_data)
{
  int offset(0), index(0), jump(0);
  
    // If we're moving in the +x direction, need plane_data[0] in -x side of interior_data
    if(omega[0] > 0)
    {
      int i = 0;
      jump = cells_x*group_per_groupset*angle_per_angleset*4*6;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (6*k*cells_xy)*group_per_groupset*angle_per_angleset*4;
        for(int j = 0; j < cells_y; j++)
        {
          index = offset +j*jump + 3*group_per_groupset*angle_per_angleset*4;
          for(int gm = 0; gm < group_per_groupset*angle_per_angleset*4; gm++, i++)
            interior_data[index+gm] = plane_data[0][i];
        }
      } 
    }
    // If we're moving in the -x direction, need plane_data[0] in +x side of interior_data
    else
    {
      int i = 0;
      jump = cells_x*group_per_groupset*angle_per_angleset*4*6;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (1 + 6*((cells_x - 1) + k*cells_xy))*group_per_groupset*angle_per_angleset*4;
        for(int j = 0; j < cells_y; j++)
        {
          index = offset +j*jump;
          for(int gm = 0; gm < group_per_groupset*angle_per_angleset*4; gm++, i++)
            interior_data[index+gm] = plane_data[0][i];
        }
      }      
    }
    // If we're moving in the +y direction, need plane_data[0] in -y side of interior_data
    if(omega[1] > 0)
    {
      int i = 0;
      jump = group_per_groupset*angle_per_angleset*4*6;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (0 + 6*k*cells_xy)*group_per_groupset*angle_per_angleset*4;
        for(int j = 0; j < cells_x; j++)
        {
          index =  offset + j*jump;
          for(int gm = 0; gm < group_per_groupset*angle_per_angleset*4; gm++, i++)
            interior_data[index+gm] = plane_data[1][i];
        }
      }      
    }
    // If we're moving in the -y direction, need plane_data[0] in +y side of interior_data
    else
    {
      int i = 0;
      jump = group_per_groupset*angle_per_angleset*4*6;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (2 + 6*(cells_x*(cells_y - 1) + k*cells_xy))*group_per_groupset*angle_per_angleset*4;
        for(int j = 0; j < cells_x; j++)
        {
          index =  offset + j*jump;
          for(int gm = 0; gm < group_per_groupset*angle_per_angleset*4; gm++, i++)
            interior_data[index+gm] = plane_data[1][i];
        }
      }       
    }
    // If we're moving in the +z direction, need plane_data[0] in -z side of interior_data
    if(omega[2] > 0)
    {
      int i = 0;
      jump = group_per_groupset*angle_per_angleset*4*6;
      offset = 4*group_per_groupset*angle_per_angleset*4;
      for (int j = 0; j < cells_xy; j++)
      {
        index =  offset + j*jump;
        for(int gm = 0; gm < group_per_groupset*angle_per_angleset*4; gm++, i++)
            interior_data[index+gm] = plane_data[2][i];
      }      
    }
    // If we're moving in the -z direction, need plane_data[0] in +z side of interior_data
    else
    {
      int i = 0;
      jump = group_per_groupset*angle_per_angleset*4*6;
      offset = (5+6*(cells_z - 1)*cells_xy)*group_per_groupset*angle_per_angleset*4;
      for (int j = 0; j < cells_xy; j++)
      {
        index =  offset + j*jump;
        for(int gm = 0; gm < group_per_groupset*angle_per_angleset*4; gm++, i++)
            interior_data[index+gm] = plane_data[2][i];
      }       
    }
    
  /*  std::cout << "Plane_data[0] = " << std::endl;
    // cells_y*cells_z*group_per_groupset*angle_per_angleset * 4
    for(int j = cells_y -1; j >= 0; j --)
      std::cout << plane_data[0][j*4] << std::endl;
    std::cout << " " << std::endl;
    std::cout << "Plane_data[1] = " << std::endl;
    for(int i = 0; i < cells_x; i ++)
      std::cout << plane_data[1][i*4] << " " ;
    std::cout << " " << std::endl;
      
    std::cout << "Plane_data[2] = " << std::endl;
    // cells_xy*group_per_groupset*angle_per_angleset * 4);
    for(int j = cells_y -1; j >= 0; j --){
      for(int i = 0; i < cells_x; i ++)
        std::cout << plane_data[2][(j*cells_x + i)*4] << " " ;
      std::cout << " " << std::endl;
    }
      
    std::cout << " " << std::endl;
   */ 

}

void Task::SetPlaneData(std::vector<double>& interior_data)
{
   int offset(0), index(0), jump(0);
  
    // If we're moving in the -x direction, need plane_data[0] in -x side of interior_data
    if(omega[0] < 0)
    {
      int i = 0;
      jump = cells_x*group_per_groupset*angle_per_angleset*4*6;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (3 + 6*k*cells_xy)*group_per_groupset*angle_per_angleset*4;
        for(int j = 0; j < cells_y; j++)
        {
          index = offset +j*jump;
          for(int gm = 0; gm < group_per_groupset*angle_per_angleset*4; gm++, i++)
            plane_data[0][i] = interior_data[index+gm] ;
        }
      } 
    }
    // If we're moving in the +x direction, need plane_data[0] in +x side of interior_data
    else
    {
      int i = 0;
      jump = cells_x*group_per_groupset*angle_per_angleset*4*6;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (1 + 6*((cells_x -1) + k*cells_xy))*group_per_groupset*angle_per_angleset*4;
        for(int j = 0; j < cells_y; j++)
        {
          index = offset +j*jump;
          for(int gm = 0; gm < group_per_groupset*angle_per_angleset*4; gm++, i++)
            plane_data[0][i] = interior_data[index+gm];
        }
      }      
    }
    // If we're moving in the -y direction, need plane_data[0] in -y side of interior_data
    if(omega[1] < 0)
    {
      int i = 0;
      jump = group_per_groupset*angle_per_angleset*4*6;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (0 + 6*k*cells_xy)*group_per_groupset*angle_per_angleset*4;
        for(int j = 0; j < cells_x; j++)
        {
          index =  offset + j*jump;
          for(int gm = 0; gm < group_per_groupset*angle_per_angleset*4; gm++, i++)
            plane_data[1][i] = interior_data[index+gm];
        }
      }      
    }
    // If we're moving in the +y direction, need plane_data[0] in +y side of interior_data
    else
    {
      int i = 0;
      jump = group_per_groupset*angle_per_angleset*4*6;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (2 + 6*(cells_x*(cells_y - 1) + k*cells_xy))*group_per_groupset*angle_per_angleset*4;
        for(int j = 0; j < cells_x; j++)
        {
          index =  offset + j*jump;
          for(int gm = 0; gm < group_per_groupset*angle_per_angleset*4; gm++, i++)
            plane_data[1][i] = interior_data[index+gm];
        }
      }       
    }
    // If we're moving in the -z direction, need plane_data[0] in -z side of interior_data
    if(omega[2] < 0)
    {
      int i = 0;
      jump = group_per_groupset*angle_per_angleset*4*6;
      offset = 4*group_per_groupset*angle_per_angleset*4;
      for (int j = 0; j < cells_xy; j++)
      {
        index =  offset + j*jump;
          for(int gm = 0; gm < group_per_groupset*angle_per_angleset*4; gm++, i++)
            plane_data[2][i] = interior_data[index+gm];
      }      
    }
    // If we're moving in the +z direction, need plane_data[0] in +z side of interior_data
    else
    {
      int i = 0;
      jump = group_per_groupset*angle_per_angleset*4*6;
      offset = (5+6*(cells_z - 1)*cells_xy)*group_per_groupset*angle_per_angleset*4;
      for (int j = 0; j < cells_xy; j++)
      {
        index =  offset + j*jump;
          for(int gm = 0; gm < group_per_groupset*angle_per_angleset*4; gm++, i++)
            plane_data[2][i] = interior_data[index+gm];
      }       
    } 

}


// This will be the function for when we have spider web grids
/*void Task::GetInteriorData(std::vector<double >& interior_data)
/{
  int offset(0), index(0);
  
    // If we're moving in the +x direction, need plane_data[0] in -x side of interior_data
    if(omega[0] > 0)
    {
      int i = 0;
      int size = cells_y*group_per_groupset*angle_per_angleset*4;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (2*(cells_x - 1) + (cells_y - 1) + k*cells_xy)*group_per_groupset*angle_per_angleset*4*6 + 3;

        for(int j = 0; j < size; j++, i++)
        {
          index =  offset + 6*j;
          interior_data[index] = plane_data[0][i];
        }
      } 
    }
    // If we're moving in the -x direction, need plane_data[0] in +x side of interior_data
    else
    {
      int i = 0;
      int size = cells_y*group_per_groupset*angle_per_angleset*4;
      for (int k = 0; k < cells_z; k++)
      {
        offset = ((cells_x - 1) + k*cells_xy)*group_per_groupset*angle_per_angleset*4*6 + 1;
        for(int j = 0; j < size; j++, i++)
        {
          index =  offset + 6*j;
          interior_data[index] = plane_data[0][i];
        }
      }      
    }
    // If we're moving in the +y direction, need plane_data[0] in -y side of interior_data
    if(omega[1] > 0)
    {
      int i = 0;
      int size = cells_x*group_per_groupset*angle_per_angleset*4;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (k*cells_xy)*group_per_groupset*angle_per_angleset*4*6 + 2;
        for(int j = 0; j < size; j++, i++)
        {
          index =  offset + 6*j;
          interior_data[index] = plane_data[1][i];
        }
      }      
    }
    // If we're moving in the -y direction, need plane_data[0] in +y side of interior_data
    else
    {
      int i = 0;
      int size = cells_x*group_per_groupset*angle_per_angleset*4;
      for (int k = 0; k < cells_z; k++)
      {
        offset = ((cells_x -1) + (cells_y - 1) + k*cells_xy)*group_per_groupset*angle_per_angleset*4*6 + 0;
        for(int j = 0; j < size; j++, i++)
        {
          index =  offset + 6*j;
          interior_data[index] = plane_data[1][i];
        }
      }       
    }
    // If we're moving in the +z direction, need plane_data[0] in -z side of interior_data
    if(omega[2] > 0)
    {
      int i = 0;
      int size = cells_xy*group_per_groupset*angle_per_angleset*4;
      for (int k = 0; k < size; k++)
      {
        offset = 0 + 4;
        index =  offset + 6*k;
        interior_data[index] = plane_data[2][k];
      }      
    }
    // If we're moving in the -z direction, need plane_data[0] in +z side of interior_data
    else
    {
      int i = 0;
      int size = cells_xy*group_per_groupset*angle_per_angleset*4;
      for (int k = 0; k < size; k++)
      {
        offset = (cells_z - 1)*cells_xy*group_per_groupset*angle_per_angleset*4*6 + 5;
        index =  offset + 6*k;
        interior_data[index] = plane_data[2][k];
      }      
    }
      
}

void Task::SetPlaneData(std::vector<double>& interior_data)
{
   int offset(0), index(0);
  
    // If we're moving in the -x direction, need plane_data[0] from -x side of interior_data
    if(omega[0] < 0)
    {
      int i = 0;
      int size = cells_y*group_per_groupset*angle_per_angleset*4;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (2*cells_x + cells_y + k*cells_xy)*group_per_groupset*angle_per_angleset*4*6 + 3;
        for(int j = 0; j < size; j++, i++)
        {
          index =  offset + 6*j;
          plane_data[0][i] = interior_data[index];
        }
      } 
    }
    // If we're moving in the +x direction, need plane_data[0] in -x side of interior_data
    else
    {
      int i = 0;
      int size = cells_y*group_per_groupset*angle_per_angleset*4;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (cells_x + k*cells_xy)*group_per_groupset*angle_per_angleset*4*6 + 1;
        for(int j = 0; j < size; j++, i++)
        {
          index =  offset + 6*j;
          plane_data[0][i] = interior_data[index];
        }
      }      
    }
    // If we're moving in the -y direction, need plane_data[0] in -y side of interior_data
    if(omega[1] < 0)
    {
      int i = 0;
      int size = cells_x*group_per_groupset*angle_per_angleset*4;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (k*cells_xy)*group_per_groupset*angle_per_angleset*4*6 + 0;
        for(int j = 0; j < size; j++, i++)
        {
          index =  offset + 6*j;
          plane_data[1][i] = interior_data[index];
        }
      }      
    }
    // If we're moving in the +y direction, need plane_data[0] in +y side of interior_data
    else
    {
      int i = 0;
      int size = cells_x*group_per_groupset*angle_per_angleset*4;
      for (int k = 0; k < cells_z; k++)
      {
        offset = (cells_x + cells_y + k*cells_xy)*group_per_groupset*angle_per_angleset*4*6 + 2;
        for(int j = 0; j < size; j++, i++)
        {
          index =  offset + 6*j;
          plane_data[1][i] = interior_data[index];
        }
      }       
    }
    // If we're moving in the -z direction, need plane_data[0] in -z side of interior_data
    if(omega[2] < 0)
    {
      int i = 0;
      int size = cells_xy*group_per_groupset*angle_per_angleset*4;
      for (int k = 0; k < size; k++)
      {
        offset = 0 + 4;
        index =  offset + 6*k;
        plane_data[2][k] = interior_data[index];
      }      
    }
    // If we're moving in the +z direction, need plane_data[0] in +z side of interior_data
    else
    {
      int i = 0;
      int size = cells_xy*group_per_groupset*angle_per_angleset*4;
      for (int k = 0; k < size; k++)
      {
        offset = (cells_z - 1)*cells_xy*group_per_groupset*angle_per_angleset*4*6 + 5;
        index =  offset + 6*k;
        plane_data[2][k] = interior_data[index];
      }      
    } 
  
  
}
 */

void Task::Get_buffer_from_bc(int dim)
{
  if (dim == 0)
  {
    int size = cells_y*cells_z*group_per_groupset*angle_per_angleset;
    for (int i = 0; i < size; i++)
    {
      int j = 4 * i;
      plane_data[0][j] = GetBoundaryCondition(0);
      plane_data[0][j + 1] = 0;
      plane_data[0][j + 2] = 0;
      plane_data[0][j + 3] = 0;
    }
  }
  else if (dim == 1)
  {
    int size = cells_x*cells_z*group_per_groupset*angle_per_angleset;
    for (int i = 0; i < size; i++)
    {
      int j = 4 * i;
      plane_data[1][j] = GetBoundaryCondition(1);
      plane_data[1][j + 1] = 0;
      plane_data[1][j + 2] = 0;
      plane_data[1][j + 3] = 0;
    }
  }
  else
  {
    int size = cells_xy*group_per_groupset*angle_per_angleset;
    for (int i = 0; i < size; i++)
    {
      int j = 4 * i;
      plane_data[2][j] = GetBoundaryCondition(2);
      plane_data[2][j + 1] = 0;
      plane_data[2][j + 2] = 0;
      plane_data[2][j + 3] = 0;
    }
  }

}
