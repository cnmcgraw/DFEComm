
#ifndef Problem_Input_h
#define Problem_Input_h 1

#include <vector>
#include <string>
#include <mpi.h>


class Problem_Input
{
public:
	Problem_Input();
	~Problem_Input();

public:
	// Data from Input File
	int problem_size;
	// Spatial Data 
	// Number of pins in x and y
	int num_pin_x, num_pin_y;

	// Number of planes in the z dimension (will be 1 for 2D problem)
	int z_planes;

	// The number of cells per pin cell will be (2*refinement)^2
	int refinement;

	// Grid type can be structured or unstructured
	// sp_disc is the spatial discretization. 
	//It can be linear discontinuous or weighted diamond difference
	int sp_disc;

	// Boundary conditions. It can be reflecting or specified.
	// If specified, all boundaries will have the same boundary condition
	int bcs;

	// Angular Data
	// Number of polar and azimuthal angles. Quadrature is product
	int num_polar, num_azim;
	// Angular Aggregation type. Options are single angle, polar, quadrant/octant
	int ang_agg_type;

	// Energy Data
	// Number of groups and groupsets. The number of groups should be able to be
	// divided into the groupsets evenly
	int num_groups, num_groupsets;

	// Partitioning Data
	// Number of Shared Memory Locations
	int num_SML;

	// Partition type: Options are KBA, Hybrid-KBA, Volumetric, or function
	// If function, more parameters are needed
	int partition_type;

	// Partition function can be blocked (1) or round robin (2) in each dimension
	// Overload = number of cellsets per processor in each dimension
	// Num_cellsets = total number of cellsets in each dimension
	// Num_cellsets/Overload should be the number of SML in each dimension
	std::vector< int > partition_function, overload, num_cellsets;
	// Sched_type: options are dependency depth or push to center
	int sched_type;

	// Number of sweeps to perform (for averaging)
	int num_sweeps;

	// Read the input file and store into data structures
	void ProcessInput(std::ifstream&, std::ofstream&);

	// Check to see if the line is a comment (delimiter = #)
	bool CheckComment(std::vector<std::string>);
	
private:
	std::string name;
	std::vector<std::string> inpL;
	int rank;

	// If the user defines a problem_size, need to specify all the parameters
	void DefineProblem();

	// Check to make sure the problem is well defined
	void CheckProblemInput();
	void GetPartitionParameters();

	std::vector<int> FactorSMLCount();

};

#endif
