
#ifndef Subdomain_h
#define Subdomain_h 1

#include <vector>
#include <queue>
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
  vector<vector<double> > buffer, Send_buffer;

	// This hold not-yet-needed messages
	vector<vector<double> > Received_buffer;
	vector<vector<int> > Received_info;
	std::queue<int> Received_open;
	int max_size;

	// Build all the cellsets this SML owns
	void BuildSubdomain(int, Problem*);

	// Compute all the CellSet IDs this SML owns (given the SML ID)
	void ComputeCellSetID(int);

  // Computes the local CS_ID given the global ID and the SML_ID
  int ComputeCellSetIndex(int CS_ID, int SML_ID);


private:
	int cells_x;
	int cells_y;
	int cells_z;

	int angle_per_angleset;
	int group_per_groupset;

  vector<int> num_cellsets;
  vector<int> overload;
  vector<int> partition_function;

};

#endif
