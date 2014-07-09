#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
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
	std::cout << "Need to specify input file (Did you forget '-f'?)" << std::endl;
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
      std::cout << "I/O error trying to open the file " << input_file_name
        << "." << std::endl;
	  abort();
    }

  // Maybe only have 1 processor do this
  std::ofstream output;
  output.open(output_file_name.c_str());
    if ((output).fail())
    {
      std::cout << "I/O error trying to open the file " << output_file_name
       << "." << std::endl;
	  abort();
    }

  // Read Input
  input_data = new Problem_Input();
  input_data->ProcessInput(input, output);

  // Build all the data structures, subdomain, quadrature, energy grid
  problem = new Problem();
  problem->BuildProblem(input_data);

  // Perform the sweep
  problem->Sweep();
  

  MPI_Finalize();
  return 0;
}