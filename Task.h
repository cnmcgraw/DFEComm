
#ifndef Task_h
#define Task_h 1

#include <vector>
#include <queue>
#include "Quadrature.h"
#include "Problem.h"
#include "Subdomain.h"

using std::vector;

//class Problem;

class Task
{
public:
	Task();
	~Task();

public:

	int cellset_id;
	int cellset_id_loc;
	int angleset_id;
	int groupset_id;
	int depth;
	int octant;
	vector<double> omega;
  vector<vector<int> > incoming, outgoing;


	// Matrices holding buffer space for incoming information
	// They are contiguous vectors of psi's in order of cells, 
	// then groups, and finally angles
  vector<vector<double> > plane_data, Send_buffer;

	// This hold not-yet-needed messages
	vector<vector<double> > Received_buffer;
	vector<vector<int> > Received_info;
	std::queue<int> Received_open;
	
	// Set up the Task
	void BuildTask(int, Problem*, Subdomain*, int, int, int);
	
	// Set the task_id from cs, as, gs
	int ComputeTaskID(int cs, int as, int g);

	// From the input, set the boundary conditions
	void SetBoundaryConditions(Problem*);

	// Function (given a boundary) returns the boundary condition
	double GetBoundaryCondition(int);

	// Allocate the send buffers
	void AllocateBuffers(void);

	// Returns an iterator to the buffer location
  inline std::vector<double>::iterator Get_buffer_loc(int cell_x, int cell_y, int cell_z, int group, int angle, int face)
  {
    if (face == 0 || face == 1)
    {
      return plane_data[0].begin() + (cell_y*cells_z*group_per_groupset*angle_per_angleset * 4 + 
        cell_z*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4);
    }
    else if (face == 2 || face == 3)
    {
      return plane_data[1].begin() + (cell_x*cells_z*group_per_groupset*angle_per_angleset * 4 +
        cell_z*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4);
    }
    else
    {
      return plane_data[2].begin() + (cell_x*cells_y*group_per_groupset*angle_per_angleset * 4 + 
        cell_y*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4);
    }
  }
        
	void Get_buffer_from_bc(int);

	int cells_x;
	int cells_y;
	int cells_z;

  int cell_per_cellset;
	int angle_per_angleset;
	int group_per_groupset;

  int num_cellset;
  int num_angleset;
  int num_groupset;
	
	int task_id;
  int index;
	
	vector<double> bc;


private:

};

#endif
