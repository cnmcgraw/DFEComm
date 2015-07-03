#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <ctime>
#include "Problem_Input.h"
#include "Problem.h"
//#include <mpi.h>

using std::string;
using std::istream;
using std::ifstream;
using std::ostream;
using std::ofstream;
using std::stringstream;

Problem_Input* input_data = 0;
Problem* problem = 0;

int main(int argc, char **argv)
{

  // Begin MPI
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Parse the command line
  string input_file_name, output_file_name;
  std::ifstream input;
  bool problem_size_bool(false), input_file_bool(false);
  int problem_size;
  for (int i = 1; i < argc; i++)
  {
    string arg = (argv)[i];
  // Check for input file
    if (arg == "-f" && i+1 < argc)
  {
    input_file_bool = true;
      input_file_name  = (argv)[i+1];
    input.open(input_file_name.c_str(), std::ios_base::in);
  }
  // Check if problem size is specified
  else if ((arg == "--problemsize" || arg == "-s") && i+1 < argc)
  {
      problem_size_bool = true;
    problem_size = atoi((argv)[i+1]);
  }
  }

  // If both an input file and a problem_size were specified, we abort
  // This could be changed to favor one or the other
  if(problem_size_bool == true && input_file_bool == true)
  {
     if (rank == 0){
      std::cout << "Both a problem size and input file were specified. Please choose only one." << std::endl;
     }
   MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // If neither an input file or a problem_size were specified, we abort
  if(problem_size_bool == false && input_file_bool == false)
  {
     if (rank == 0){
      std::cout << "Please specify an input file or a problem_size." << std::endl;
     }
   MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // If we've gotten here, the use has either specified an input file, or a problem_size (not both)
  std::ofstream output;
  if(input_file_bool == true)
  {
    // If we specified an input file, we need to make sure a valid input file was opened
    if (input_file_name.empty())
    {
     if (rank == 0){ std::cout << "Need to specify input file (Did you forget '-f'?)" << std::endl; }
     MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Derive the output file name
    output_file_name = input_file_name;
    int n = output_file_name.rfind('.');
    if (n != string::npos)
    output_file_name.erase(n, string::npos);

    output_file_name += ".out";

    // Create input and output streams
    if ((input).fail())
    {
      if (rank == 0){
        std::cout << "I/O error trying to open the file " << input_file_name
          << "." << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    output.open(output_file_name.c_str());
    if ((output).fail())
    {
      if (rank == 0){
        std::cout << "I/O error trying to open the file " << output_file_name
          << "." << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Read Input
    input_data = new Problem_Input();
    input_data->ProcessInput(input, output);
  }
  // Problem Size has been specified
  else
  {
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    // Default output file name
    
    std::stringstream ss;
    ss << "Output_" << mpi_size << ".out";
    output_file_name = ss.str();

    output.open(output_file_name.c_str());
    if ((output).fail())
    {
      if (rank == 0){
        std::cout << "I/O error trying to open the file " << output_file_name
          << "." << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Process the inputs from the Command line
    input_data = new Problem_Input();
    MPI_Comm_size(MPI_COMM_WORLD,&input_data->num_SML);
    input_data->problem_size = problem_size;
    input_data->DefineProblem();
    input_data->CheckProblemInput();
  }

  // Build all the data structures, subdomain, quadrature, energy grid
  problem = new Problem();
  problem->BuildProblem(input_data);
  // Perform the sweep
  double start;
  long double duration, total_duration;
  total_duration = 0;
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i = 0; i < input_data->num_sweeps; i++)
  {   
    problem->ZeroPhi();
    start = MPI_Wtime();
    if (rank == 0){ output << "  SWEEP " << i + 1 << std::endl; }
    problem->Sweep(output);
    MPI_Barrier(MPI_COMM_WORLD);
    duration = (MPI_Wtime() - start); // / (double)CLOCKS_PER_SEC;
    if (rank == 0){
      output << "  Sweep " << i + 1 << " took " << duration << " seconds." << std::endl;
      output << " " << std::endl;
    }
    total_duration += duration;
  }
  total_duration = total_duration / input_data->num_sweeps;
  if (rank == 0){
    output << "Average Sweep Time: " << total_duration << " seconds." << std::endl;
  }

  MPI_Finalize();

  delete input_data;
  delete problem;
  input.close();
  output.close();

  return 0;
}
