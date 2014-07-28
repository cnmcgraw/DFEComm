#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
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
  MPI_Init (&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Parse the command line
  string input_file_name, output_file_name;
  std::ifstream input;
  for (int i = 1; i < argc; i++)
  {
    string arg = (argv)[i];
    if (arg == "-f" && i+1 < argc)
	{
      input_file_name  = (argv)[i+1];
	  input.open(input_file_name.c_str(), std::ios_base::in);
	}
  }

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
	  abort();
    }

  // Maybe only have 1 processor do this
  std::ofstream output;
  output.open(output_file_name.c_str());
    if ((output).fail())
    {
		if (rank == 0){
			std::cout << "I/O error trying to open the file " << output_file_name
				<< "." << std::endl;
		}
	  abort();
    }

  // Read Input
  input_data = new Problem_Input();
  input_data->ProcessInput(input, output);

  // Build all the data structures, subdomain, quadrature, energy grid
  problem = new Problem();
  problem->BuildProblem(input_data);
 
  // Perform the sweep
  std::clock_t start;
  long double duration, total_duration;
  total_duration = 0;
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i = 0; i < input_data->num_sweeps; i++)
  {	  
	  problem->ZeroPhi();
	  start = std::clock();
	  if (rank == 0){ output << "  SWEEP " << i + 1 << std::endl; }
	  problem->Sweep(output);
	  MPI_Barrier(MPI_COMM_WORLD);
	  duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
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