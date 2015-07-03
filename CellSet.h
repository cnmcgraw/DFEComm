
#ifndef CellSet_h
#define CellSet_h 1

#include <vector>
#include <iostream>
#include "Cell.h"
#include "Problem.h"

//class Problem;

class CellSet
{
public:
	CellSet();
	~CellSet();

public:
	// Vector of cells this CellSet has
	std::vector< Cell > Cells;

	int cells_per_cellset;
	int cells_x, cells_y, cells_z;
	int num_groupset, num_angleset;

	// Vector of neighbor cellsets
	std::vector<Neighbor> neighbors;

	// Vector indicating which global boundary (if any)
	// the cellset is on.
	// -1 indicates the negative (x, y, or z) boundary
	// +1 indicates the positive (x, y, or z) boundary
	// 0 indicates the cellset is in the interior
	std::vector<int> globalboundary;
	std::vector<std::vector<double> > BoundaryFlux;



	// Given an CellSet ID, build its cells
	void BuildCellSet(int, Problem*);

	// Compute the global cell ideas in the CellSet
	void ComputeCellIDs(int, Problem*);

	// Given the CS_ID, and the number of cellset in x,
	// y, and z, return the global ijk of the cellset
	void GetCellSetijk(int, int, int, int, std::vector<int>&);

	// Given the CS_ID, and the octant of the angleset,
	//return the depth to the furthest corner
	int GetDepth(int, int, Problem*);

	// This function sets the globalboundary vector
	void SetGlobalBoundary(int, Problem*);

	// This function gets the upstream direction and 
	// cellset id of each of this cellset neighbors
	void GetNeighbors(int);

	// This function returns the SML of each neighbor
	void GetNeighborsSML(int, Problem*);

	// This function sets the boundary buffer
	void SetBoundaryFlux(int, int, int, std::vector<double>&);



private:
	std::vector<int> num_cellsets;

};

#endif
