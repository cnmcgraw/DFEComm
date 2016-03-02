
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

	// Returns an iterator to the buffer location in interior data
//  inline std::vector<double>::iterator Get_buffer_loc(int cell_xy, int cell_z, int group, int angle, int face, std::vector<double>interior_data)
  inline int Get_buffer_loc(int cell_xy, int cell_z, int group, int angle, int face, std::vector<double>interior_data)
  {
    int cell_y = (int)(cell_xy / cells_x);
    int cell_x = cell_xy - cells_x * cell_y;
    int index = (cell_z*cells_xy*group_per_groupset*angle_per_angleset * 4*6 + 
      cell_y*cells_x*group_per_groupset*angle_per_angleset * 4*6 + cell_x*group_per_groupset*angle_per_angleset * 4*6 + face*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4);

    return index;
    
    
 //   return interior_data.begin() + (cell_z*cells_xy*group_per_groupset*angle_per_angleset * 4*6 + 
 //     cell_xy*group_per_groupset*angle_per_angleset * 4*6 + face*group_per_groupset*angle_per_angleset * 4 + group*angle_per_angleset * 4 + angle * 4);  
  }
  
	// Gets the plane data and puts it in interior data
  void GetInteriorData(std::vector<double>& interior_data);
        
  // Puts the interior data into plane data
  void SetPlaneData(std::vector<double>& interior_data);
        
	void Get_buffer_from_bc(int);

        // cells_x and cells_y will be the number of cells in x and y on the boundary.
        // cells_xy will be the total number of cells in the xy plane.
        // In rectangular grids, cells_xy will be the simple product of cells_x and cells_y.
        // In spider web grids, we'll need to sum up the total number of cells for this.
	int cells_x;
	int cells_y;
	int cells_xy;
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
