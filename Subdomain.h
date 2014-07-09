
#ifndef Subdomain_h
#define Subdomain_h 1

#include <vector>
#include "CellSet.h"

using std::vector;

class Problem;

class Subdomain
{
public:
	Subdomain();
	~Subdomain();

public:
	// The total number of cellsets this SML owns
	int total_overload;
	// Vector of boundary conditions
	vector<double> bc;

	// List of cellset IDs this SML owns
	vector<int> CellSetIDs;
	vector< CellSet > CellSets;

	// Matrices holding buffer space for incoming information
	// They are contiguous vectors of psi's in order of cells, 
	// (in the dimensions not in the name)
	// then groups, and finally angles
	vector<double> X_buffer;
	vector<double> Y_buffer;
	vector<double> Z_buffer;

	// Build all the cellsets this SML owns
	void BuildSubdomain(int, Problem*);

	// Compute all the CellSet IDs this SML owns (given the SML ID)
	void ComputeCellSetID(int, Problem*);

	// From the input, set the boundary conditions
	void SetBoundaryConditions(Problem*);

	// Function (given a boundary) returns the boundary condition
	double GetBoundaryCondition( int);

	// Set functions for the X, Y, and Z buffers
	void Set_buffer(int, int, int, int, int, int, vector<double>&);
	void Get_buffer(int, int, int, int, int, int, vector<double>&);
	void Set_buffer_from_bc(int);


private:

	int cells_x;
	int cells_y;
	int cells_z;

	int angle_per_angleset;
	int group_per_groupset;

};

#endif