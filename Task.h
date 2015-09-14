
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
	Direction omega;
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

  // Sets the index in the AllTasks array
  void SetIndex(int i);

	// Set functions for the X, Y, and Z buffers
	void Set_buffer(int, int, int, int, int, int, int, vector<double>&);
	void Get_buffer(int, int, int, int, int, int, vector<double>&);
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
