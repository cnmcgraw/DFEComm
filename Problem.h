
#ifndef Problem_h
#define Problem_h 1

#include <vector>
#include <string>
#include <fstream>
#include "Problem_Input.h"
#include "Quadrature.h"

class Problem
{
public:
	Problem();
	virtual ~Problem();

public:
	int rank;
	// Data Structures from the input
	// Spatial Data
	int num_pin_x, num_pin_y;
	int z_planes;
	int refinement;
	int sp_disc;
	int bcs;
	//Subdomain* subdomain;
	// Angular Data
	int num_polar, num_azim;
	int ang_agg_type;
	Quadrature quad;

	// Energy Data
	int num_groupsets, group_per_groupset;

	// Partitioning Data
	int num_SML;
	int partition_type;
	std::vector< int > partition_function, overload, num_cellsets;
	int sched_type;

	// Solution Data Structures



public:	
	// This function builds the subdomain, quadrature, and energy grid
	void BuildProblem(Problem_Input*);

	// This function performs the sweep
	void Sweep(std::ofstream&);

	// Zeroes out the phi vector before each sweep.
	void ZeroPhi();

private:
	// Matrices to invert during the sweep
	// Their inversion will be replaced with a wait() function
	std::vector< std::vector< double > > A_tilde, A;
	std::vector< double > RHS;
	std::vector< double > bg;
	std::vector< double > M;
	std::vector<std::vector<std::vector< Direction > > > N;
	std::vector<std::vector< Direction > > L;
	std::vector<double> source;
	Direction omega,facecenter;
	double sigma_t, cell_average;
	int num_tasks;

	int GetCell(int, int, int, int, int, int, int);

	void GE_no_pivoting(std::vector<std::vector< double > >&, std::vector<double>&, int);

	// Returns the Target of the messages being sent or received
	unsigned int GetTarget(int, int, int, int);

	// Gets the next open place is the Receive buffer
	int GetPlacement();


};

#endif
