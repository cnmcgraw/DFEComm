#include "Problem_Input.h"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


Problem_Input::Problem_Input()
{ }

Problem_Input::~Problem_Input()
{ }

bool Problem_Input::CheckComment(std::vector<std::string> inputline) {

	if ( inputline[0] == "#") return true;
	if ( inputline.empty() ) return true;
	return false;
}
void Problem_Input::ProcessInput(std::ifstream& input, std::ofstream& fout)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// Build and initialize mesh parameters
	partition_function.resize(3, 1);
	overload.resize(3, 1);
	num_cellsets.resize(3, 1);

	// Problem size indicates if the user has specified
	// tiny, small, medium, or large problem (1,2,3, or 4)
	// If they don't specify it, the problem parameters
	// need to be user input
	problem_size = 0;

	do {
		std::string line;
		do {
			// Read the data line
			inpL.clear();
			if( !line.empty() ) {
				std::string token;
				std::istringstream iss(line);
				while ( getline(iss,token,' ' ) ) {
					if( !token.empty() ) {
						inpL.push_back(token);
					}
				}
				// Check for comments and the end of the input file
				bool isComment = CheckComment(inpL);
				if(!isComment) {
					// Every line should have an ID in the 0th element and a value in the 2nd
					if(inpL.size() < 3 && inpL.size() > 1) {
						std::cout << "No set value for " << inpL[0] << ". " << std::endl;
						break;
					}
					if (inpL[0] == "problem_size")
						problem_size = atoi(inpL[2].c_str());
					
					if (problem_size == 0)
					{	
						//Spatial Data
						if (inpL[0] == "pin_x")
							pin_x = atoi(inpL[2].c_str());
						if (inpL[0] == "pin_y")
							pin_y = atoi(inpL[2].c_str());
						if (inpL[0] == "z_planes")
							z_planes = atoi(inpL[2].c_str());
						if (inpL[0] == "refinement")
							refinement = atoi(inpL[2].c_str());
						bcs = 2;
						if (inpL[0] == "bc")
							bcs = atoi(inpL[2].c_str());

						sp_disc = 1;
						if (inpL[0] == "spatial_discretization")
							sp_disc = atoi(inpL[2].c_str());

						// Angular Data
						if (inpL[0] == "num_polar")
							num_polar = atoi(inpL[2].c_str());
						if (inpL[0] == "num_azim")
							num_azim = atoi(inpL[2].c_str());
						if (inpL[0] == "angle_aggregation")
							ang_agg_type = atoi(inpL[2].c_str());

						// Energy Data
						if (inpL[0] == "num_groups")
							num_groups = atoi(inpL[2].c_str());
						if (inpL[0] == "num_groupsets")
							num_groupsets = atoi(inpL[2].c_str());

						if (inpL[0] == "partition_type")
							partition_type = atoi(inpL[2].c_str());

						GetPartitionParameters();

						if (inpL[0] == "sched_type")
							sched_type = atoi(inpL[2].c_str());
						if (inpL[0] == "num_sweeps")
							num_sweeps = atoi(inpL[2].c_str());
					}

					//Partitioning Data
					if (inpL[0] == "num_SML")
						num_SML = atoi(inpL[2].c_str());
					if (inpL[0] == "num_TpSML")
						num_TpSML = atoi(inpL[2].c_str());
				}
			}
		}while (getline(input,line,'\n') );
	}while( !input.eof() );

	// If the user has defined a problem size, need to define the problem
	if (problem_size != 0)
		DefineProblem();

	// Check for a valid problem (only 1 processor)
	CheckProblemInput();
	//Output(fout);
}
void Problem_Input::GetPartitionParameters()
{
	// Default is KBA (Hybrid for num_SML>1) with overloading of 1
	if(partition_type == 0){
		if (num_SML == 1)
		{
			num_cellsets[0] = (int)sqrt(num_SML);
			num_cellsets[1] = (int)sqrt(num_SML);
			num_cellsets[2] = z_planes;
			overload[0] = 1;
			overload[1] = 1;
			overload[2] = z_planes;
		}
		else
		{
			num_cellsets[0] = (int)sqrt(num_SML / 2);
			num_cellsets[1] = (int)sqrt(num_SML / 2);
			num_cellsets[2] = z_planes / 2;
			overload[0] = 1;
			overload[1] = 1;
			overload[2] = z_planes;
		}
	}
	// KBA
	if(partition_type == 1){
		if(inpL[0] == "overload")
			overload[0] = atoi( inpL[2].c_str() );

		num_cellsets[0] = (int)sqrt(num_SML)/overload[0];
		num_cellsets[1] = (int)sqrt(num_SML)/overload[0];
		num_cellsets[2] = z_planes;
		overload[1] = overload[0];
		overload[2] = z_planes;

	}
	// Hybrid KBA ( Py = 2 (2D), Pz = 2 (3D) )
	if(partition_type == 2){
		if(inpL[0] == "overload")
			overload[0] = atoi( inpL[2].c_str() );

		num_cellsets[0] = (int)sqrt(num_SML/2)/overload[0];
		num_cellsets[1] = (int)sqrt(num_SML/2)/overload[0];
		num_cellsets[2] = z_planes/2;
		overload[1] = overload[0];
		overload[2] = z_planes;
	}
	// Extended Hybrid
	if(partition_type == 3){
		if (rank == 0){
			std::cout << "Extended Hybrid not implemented yet. Using default"
				<< " partitioning (KBA with overload=1)" << std::endl;
		}
		partition_type = 0;
	}
	// Volumetric
	if(partition_type == 4){
		num_cellsets[0] = (int)pow(num_SML,(1/3));
		num_cellsets[1] = (int)pow(num_SML,(1/3));
		num_cellsets[2] = (int)pow(num_SML,(1/3));
	}
	// Volumetric with Red Robin
	if(partition_type == 5){
		partition_function[0] = 2;
		partition_function[1] = 2;
		partition_function[2] = 2;
		if(inpL[0] == "overload")
			overload[0] = atoi( inpL[2].c_str() );
		overload[1] = overload[0];
		overload[2] = overload[0];			
			
		num_cellsets[0] = (int)pow(num_SML/overload[0],(1/3));
		num_cellsets[1] = (int)pow(num_SML/overload[0],(1/3));
		num_cellsets[2] = (int)pow(num_SML/overload[0],(1/3));
	}
	// Need to add partition parameters for predefined partition types
	if(partition_type == 6)
	{
		if(inpL[0] == "func_x")
			partition_function[0] = atoi( inpL[2].c_str() );
		if(inpL[0] == "func_y")
			partition_function[1] = atoi( inpL[2].c_str() );
		if(inpL[0] == "func_z")
			partition_function[2] = atoi( inpL[2].c_str() );
		if(inpL[0] == "overload_x")
			overload[0] = atoi( inpL[2].c_str() );
		if(inpL[0] == "overload_y")
			overload[1] = atoi( inpL[2].c_str() );
		if(inpL[0] == "overload_z")
			overload[2] = atoi( inpL[2].c_str() );
		if(inpL[0] == "num_cellsets_x")
			num_cellsets[0] = atoi( inpL[2].c_str() );
		if(inpL[0] == "num_cellsets_y")
			num_cellsets[1] = atoi( inpL[2].c_str() );
		if(inpL[0] == "num_cellsets_z")
			num_cellsets[2] = atoi( inpL[2].c_str() );
	}

}
void Problem_Input::DefineProblem()
{
	bcs = 2;
	sp_disc = 1;
	sched_type = 1;
	// KBA if only 1 SML, Hybrid if more than 1 SML
	// Factor SML Count will return the optimal 
	// X Y SML layout (as close to sqrt as possible)
	std::vector<int> SML_Counts(3, 0);
	SML_Counts = FactorSMLCount();

	num_cellsets[0] = SML_Counts[0];
	num_cellsets[1] = SML_Counts[1];
	num_cellsets[2] = SML_Counts[2];
	overload[0] = 1;
	overload[1] = 1;
	overload[2] = 1;

	// Tiny Problem
	if (problem_size == 1)
	{
		// Quadrature Data
		num_polar = 4;
		num_azim = 4;
		ang_agg_type = 3; // Octant

		// Energy Data
		num_groups = 10;
		num_groupsets = 1;

		// Spatial Data
		// A cellset = 1 pincell which is 2x2x10 cells
		refinement = 1;
		pin_x = num_cellsets[0];
		pin_y = num_cellsets[1];
		z_planes = 10*num_cellsets[2];
	}
	// Small Problem
	else if (problem_size == 2)
	{
		// Quadrature Data
		num_polar = 4;
		num_azim = 16;
		ang_agg_type = 3; // Octant

		// Energy Data
		num_groups = 50;
		num_groupsets = 1;

		// Spatial Data
		// A cellset = 1 pincell which is 10x10x100 cells
		refinement = 5;
		pin_x = num_cellsets[0];
		pin_y = num_cellsets[1];
		z_planes = 100 * num_cellsets[2];
	}
	// Medium Problem
	else if (problem_size == 3)
	{
		// Quadrature Data
		num_polar = 16;
		num_azim = 16;
		ang_agg_type = 3; // Octant

		// Energy Data
		num_groups = 100;
		num_groupsets = 1;

		// Spatial Data
		// A cellset = 1 pincell which is 10x10x200 cells
		refinement = 5;
		pin_x = num_cellsets[0];
		pin_y = num_cellsets[1];
		z_planes = 100 * num_cellsets[2];
	}
	// Large Problem - This should test the limits of an SML
	else if (problem_size == 4)
	{
		// Quadrature Data
		num_polar = 16;
		num_azim = 64;
		ang_agg_type = 3; // Octant

		// Energy Data
		num_groups = 500;
		num_groupsets = 1;

		// Spatial Data
		// A cellset = 1 pincell which is 10x10x200 cells
		refinement = 5;
		pin_x = num_cellsets[0];
		pin_y = num_cellsets[1];
		z_planes = 200 * num_cellsets[2];
	}
}
void Problem_Input::CheckProblemInput()
{
	// Check that spatial aggregation is integer multiple of 1/4 pin cells
	if(!(2*pin_x/num_cellsets[0] == (int)2*pin_x/num_cellsets[0])){
		if (rank == 0){ std::cout << "Invalid Partition in x" << std::endl; }
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	if(!(2*pin_y/num_cellsets[1] == (int)2*pin_y/num_cellsets[1])){
		if (rank == 0){ std::cout << "Invalid Partition in y" << std::endl; }
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	if(!(z_planes/num_cellsets[2] == (int)z_planes/num_cellsets[2])){
		if (rank == 0){ std::cout << "Invalid Partition in z" << std::endl; }
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	// Check that partitioning gives the correct number of SML
	int p = 1;
	for(int i = 0; i < 3; i++)
		p = p*num_cellsets[i]/overload[i];
	if(!(p == num_SML)){
		if (rank == 0){
			std::cout << "Invalid Partitioning" << std::endl;
			std::cout << "Specified number of SML = " << num_SML << " and P_eff = " << p << std::endl;
		}
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	// Check that the number of angles per angleset is an integer
	if(ang_agg_type == 2 && !(num_polar/2 == (int)num_polar/2)){
		if (rank == 0){ std::cout << "Invalid Polar Aggregation" << std::endl; }
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	if(!(num_azim/4 == (int)num_azim/4)){
		if (rank == 0){ std::cout << "Invalid Number of Azimuthal Angles" << std::endl; }
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	// Check that the number of groups per groupset is an integer
	if(!(num_groups/num_groupsets == (int)num_groups/num_groupsets)){
		if (rank == 0){ std::cout << "Non-integer number of groups per groupset" << std::endl; }
		MPI_Abort(MPI_COMM_WORLD,1);
	}



}
std::vector<int> Problem_Input::FactorSMLCount()
{
	std::vector<int> SMLCounts(3,0);

	// Check if the number of SMLs are greater than 2
	if (num_SML >= 2)
	{
		if (num_SML & 2 != 0)
		{
			if (rank == 0){ std::cout << "Automatic SML Factoring cannot handle odd SML Counts" << std::endl; }
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		// if num_SML > 2 and even, Pz = 2
		SMLCounts[2] = 2;
	}
	else
	{
		SMLCounts[0] = 1;
		SMLCounts[1] = 1;
		SMLCounts[2] = 1;
		return SMLCounts;
	}

	// If the sqrt of the remaining SML is an integer
	// Px = Py = sqrt(num_SML/2)

	if (fmod(sqrt(num_SML / 2), 1) == 0)
	{
		SMLCounts[0] = sqrt(num_SML / 2);
		SMLCounts[1] = sqrt(num_SML / 2);
		return SMLCounts;
	}
	// else find the biggest factor (closest to the sqrt root)
	// Px will be the larger factor
	else
	{
		int root = (int)sqrt(num_SML / 2);
		for (int i = root; i >= 1; i--)
		{
			// While i divides n, Px = num_SML/2/i, Py = i
			if ((num_SML/2)%i == 0)
			{
				SMLCounts[0] = num_SML / (2 * i);
				SMLCounts[1] = i;
				return SMLCounts;
			}
		}
	}
}