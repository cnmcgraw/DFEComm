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
{   
  // Defaults
  partition_function.resize(3, 1);
  overload.resize(3, 1);
  num_cellsets.resize(3, 1);
  problem_size = 0;
  sp_disc = 1;
  num_sweeps = 5;
  partition_type = 0;
  sched_type = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &num_SML);
  bcs = 2;

  partition_bool = false;
}

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
            if (inpL[0] == "num_pin_x")
              num_pin_x = atoi(inpL[2].c_str());
            if (inpL[0] == "num_pin_y")
              num_pin_y = atoi(inpL[2].c_str());
            if (inpL[0] == "z_planes")
              z_planes = atoi(inpL[2].c_str());
            if (inpL[0] == "refinement")
              refinement = atoi(inpL[2].c_str());
            if (inpL[0] == "sp_azim"){
              // If sp_azim is specified, we're doing spiderweb grids
              spider = true;
              // Check if refinement has been defined yet
              if(refinement == 0){
                if (rank == 0){ std::cout << "Refinement must be defined before sp_azim" << std::endl;}
                MPI_Abort(MPI_COMM_WORLD,1);
              }
              // Check that enough sp_azim have been specified
              if( inpL.size() < (refinement + 2) ){
                if (rank == 0){ std::cout << "Not enough sp_azim specified" << std::endl;}
                MPI_Abort(MPI_COMM_WORLD,1);
              }
              sp_azim.resize(refinement,0);
              for(int i = 0; i < refinement; i++){
                sp_azim[i] = atoi(inpL[2+i].c_str());
              }
            }
               
            if (inpL[0] == "bc")
              bcs = atoi(inpL[2].c_str());

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

            if (inpL[0] == "partition_type" || partition_bool){
              if(!partition_bool)
                partition_type = atoi(inpL[2].c_str());
              partition_bool = true;
              GetPartitionParameters();
            }
            
            if (inpL[0] == "sched_type")
              sched_type = atoi(inpL[2].c_str());
            if (inpL[0] == "num_sweeps")
              num_sweeps = atoi(inpL[2].c_str());
          }

          //Partitioning Data
          if (inpL[0] == "num_SML")
            num_SML = atoi(inpL[2].c_str());
        }
      }
    }while (getline(input,line,'\n') );
  }while( !input.eof() );

  // If the user has defined a problem size, need to define the problem
  if (problem_size != 0)
    DefineProblem();

  // Check for a valid problem 
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
      overload[0] = 1;
      overload[1] = 1;
      if(inpL[0] == "overload")
        overload[2] = z_planes / atoi( inpL[2].c_str() );
      else
        overload[2] = z_planes;
        
      num_cellsets[2] = overload[2];
    }
    else
    {
      num_cellsets[0] = (int)sqrt(num_SML / 2);
      num_cellsets[1] = (int)sqrt(num_SML / 2);
      overload[0] = 1;
      overload[1] = 1;
      if(inpL[0] == "overload")
        overload[2] = z_planes / atoi( inpL[2].c_str() );
      else
        overload[2] = z_planes / 2;
        
      num_cellsets[2] = overload[2] * 2;
    }
  }
  // KBA
  if(partition_type == 1){
    if(inpL[0] == "overload")
      overload[0] = atoi( inpL[2].c_str() );

    num_cellsets[0] = (int)sqrt(num_SML) * overload[0];
    num_cellsets[1] = (int)sqrt(num_SML) * overload[0];
    num_cellsets[2] = z_planes;
    overload[1] = overload[0];
    overload[2] = z_planes;

  }
  // Hybrid KBA ( Pz = 2 )
  if(partition_type == 2){
    if(inpL[0] == "overload")
      overload[0] = atoi( inpL[2].c_str() );

    num_cellsets[0] = (int)sqrt(num_SML/2) * overload[0];
    num_cellsets[1] = (int)sqrt(num_SML/2) * overload[0];
    num_cellsets[2] = z_planes;
    overload[1] = overload[0];
    overload[2] = z_planes /2;
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
      
    num_cellsets[0] = (int)pow(num_SML,(1/3)) * overload[0];
    num_cellsets[1] = (int)pow(num_SML,(1/3)) * overload[1];
    num_cellsets[2] = (int)pow(num_SML,(1/3)) * overload[2];
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
    // A process owns 1 pincell which is 2x2x10 cells
    refinement = 1;
    num_pin_x = num_cellsets[0];
    num_pin_y = num_cellsets[1];
    z_planes = 16*num_cellsets[2];
    overload[2] = 16;
    num_cellsets[2] *= 16;
  }
  // Small Problem
  else if (problem_size == 2)
  {
    // Quadrature Data
    num_polar = 4;
    num_azim = 16;
    ang_agg_type = 3; // Octant

    // Energy Data
    num_groups = 10;
    num_groupsets = 1;

    // Spatial Data
    // A process owns 1 pincell which is 10x10x100 cells
    refinement = 3;
    num_pin_x = num_cellsets[0];
    num_pin_y = num_cellsets[1];
    z_planes = 50 * num_cellsets[2];
    overload[2] = 50;
    num_cellsets[2] *= 50;
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
    // A process owns 1 pincell which is 10x10x200 cells
    refinement = 5;
    num_pin_x = num_cellsets[0];
    num_pin_y = num_cellsets[1];
    z_planes = 200 * num_cellsets[2];
    overload[2] = 200;
    num_cellsets[2] *= 200;
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
    // A process owns  1 pincell which is 10x10x200 cells
    refinement = 5;
    num_pin_x = num_cellsets[0];
    num_pin_y = num_cellsets[1];
    z_planes = 200 * num_cellsets[2];
    overload[2] = 200;
    num_cellsets[2] *= 200;
  }
}
void Problem_Input::CheckProblemInput()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Check that spatial aggregation is integer multiple of 1/4 pin cells
  if(!(2*num_pin_x % num_cellsets[0] == 0)){
    if (rank == 0){ std::cout << "Invalid Partition in x" << std::endl; }
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  if(!(2*num_pin_y % num_cellsets[1] == 0)){
    if (rank == 0){ std::cout << "Invalid Partition in y" << std::endl; }
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  if(!(z_planes % num_cellsets[2] == 0)){
    if (rank == 0){ std::cout << "Invalid Partition in z" << std::endl; }
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  
  // Check to make sure we either stay the same or multiply by a factor of 2 
  // for each ring 
  if(spider){
    if(sp_azim[0] % 4 != 0){
      if (rank == 0){ std::cout << "Inner most azimuthal discretization must be a factor of 4" << std::endl; }
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    for(int i = 0; i < refinement - 1; i++){
      double factor = double(sp_azim[i+1]) / double(sp_azim[i]);
      if(!(factor == 1 || fmod(factor, 2.) == 0)){
        if (rank == 0){ std::cout << "Azimuthal discretization must stay the same or multiply by a factor of 2 as we move outward" << std::endl; }
        MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
  }

  // Check that partitioning gives the correct number of SML
  int p = 1;
  for(int i = 0; i < 3; i++)
    p = p*num_cellsets[i]/overload[i];
  
  if(!(p == num_SML)){
    if (rank == 0){
      std::cout << "Invalid Partitioning" << std::endl;
      std::cout << "Specified number of SML = " << num_SML << " and P_eff(x,y,z,p_tot) = (" << num_cellsets[0]/overload[0] << "," << num_cellsets[1]/overload[1] << "," << num_cellsets[2]/overload[2] << "," << p << ")" <<std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  // Check that the number of angles per angleset is an integer
  if(num_polar%2){
    if (rank == 0){ std::cout << "Invalid Polar Aggregation" << std::endl; }
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  if(num_azim%4){
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
    if (num_SML % 2 != 0)
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
    SMLCounts[0] = (int)sqrt(num_SML / 2);
    SMLCounts[1] = (int)sqrt(num_SML / 2);
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
void Problem_Input::SetVerbose(bool verbose_bool)
{
  verbose = verbose_bool;
}
