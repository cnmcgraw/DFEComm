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

void usage(){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0){
    printf("Usage:  [mpirun ...] ./DFEComm [options ...]\n");
    printf("Where options are:\n");
    printf(" -f [Input file]       Input file to run\n");
    printf(" -o [Output file]      Output file to write to\n");
    printf(" -v                    Option for verbose output\n");
    printf(" --problemsize [S]     Preset problem size\n");
    printf(" -s [S]                Preset problem size\n");
    printf("                       Options are 1,2,3, or 4\n");
    printf("                       One z plane/cellset, octant angle aggregation, single groupset\n");
    printf("                       1: Mesh/core = 2x2x10,    Quad = 4x4,   Groups = 10\n");
    printf("                       2: Mesh/core = 10x10x100, Quad = 4x16,  Groups = 50\n");
    printf("                       3: Mesh/core = 10x10x200, Quad = 16x16, Groups = 100\n");
    printf("                       4: Mesh/core = 10x10x200, Quad = 16x64, Groups = 500\n");  
  }
  MPI_Finalize();
  exit(1);

}
int main(int argc, char **argv)
{

  // Begin MPI
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Parse the command line
  string input_file_name, output_file_name;
  std::ifstream input;
  bool problem_size_bool(false), input_file_bool(false), verbose_bool(false);
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
    // Check if output file name is specified
    else if (arg == "-o" && i+1 < argc)
    {
      output_file_name  = (argv)[i+1];   
    }
    // Check if the user wants verbose output
    else if (arg == "-v")
    {
      verbose_bool = true;  
    }
  }

  // If both an input file and a problem_size were specified, we abort
  // This could be changed to favor one or the other
  if(problem_size_bool == true && input_file_bool == true)
  {
   usage();
  }
  // If neither an input file or a problem_size were specified, we abort
  if(problem_size_bool == false && input_file_bool == false)
  {
   usage();
  }
  // If we've gotten here, the use has either specified an input file, or a problem_size (not both)
  std::ofstream output;
  if(input_file_bool == true)
  {
    // If we specified an input file, we need to make sure a valid input file was opened
    if (input_file_name.empty())
    {
     usage();
    }

    // Derive the output file name
    if(output_file_name.empty())
    {
      output_file_name = input_file_name;
      int n = output_file_name.rfind('.');
      if (n != string::npos)
        output_file_name.erase(n, string::npos);

      output_file_name += ".out";
    }

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
    input_data->SetVerbose(verbose_bool);
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
    input_data->SetVerbose(verbose_bool);
  }

  // Build all the data structures, subdomain, quadrature, energy grid
  problem = new Problem();
  problem->BuildProblem(input_data);
    
  // First print problem parameters to output file
  if(rank == 0)
  {
    time_t currentTime;
    struct tm *localTime;

    time( &currentTime );                   // Get the current time
    localTime = localtime( &currentTime );  // Convert the current time to the local time
    int Day    = localTime->tm_mday;
    int Month  = localTime->tm_mon + 1;
    int Year   = localTime->tm_year + 1900;
    int Hour   = localTime->tm_hour;
    int Min    = localTime->tm_min;
    int Sec    = localTime->tm_sec;
    
    int total_cells = 4*problem->refinement*problem->refinement*problem->z_planes*problem->num_pin_x*problem->num_pin_y;
    int cells_per_SML = total_cells/problem->num_SML;
    int unknowns_per_SML = cells_per_SML*problem->num_polar*problem->num_azim*problem->group_per_groupset*problem->group_per_groupset;
    int ang_task_num = 0;
    if(problem->ang_agg_type == 1)
      ang_task_num = problem->num_polar*problem->num_azim;
    else if(problem->ang_agg_type == 2)
      ang_task_num = problem->num_azim;
    else
      ang_task_num = 8;
    int unknowns_per_task = unknowns_per_SML / (ang_task_num * problem->num_groupsets * problem->overload[0]*problem->overload[1]*problem->overload[2]);

    output << "SimpleLD run at " << Day << "/" << Month << "/" << Year << " " << Hour << ":" << Min << ":" << Sec << std::endl;
    output << "   _____ _                 _      _      _____  " << std::endl;
    output << "  / ____(_)               | |    | |    |  __ \\ " << std::endl;
    output << " | (___  _ _ __ ___  _ __ | | ___| |    | |  | |" << std::endl;
    output << "  \\___ \\| | '_ ` _ \\| '_ \\| |/ _ \\ |    | |  | |" << std::endl;
    output << "  ____) | | | | | | | |_) | |  __/ |____| |__| |" << std::endl;
    output << " |_____/|_|_| |_| |_| .__/|_|\\___|______|_____/ " << std::endl;
    output << "                    | |                         " << std::endl;
    output << "                    |_|                         " << std::endl;

    output << "================================================================" << std::endl;
    output << "  Problem parameters: " <<std::endl;
    output << "  Pin_x * Pin_y                                 : " << problem->num_pin_x << "*" << problem->num_pin_y << std::endl;
    output << "  Number of cells per pincell                   : " << 2*problem->refinement << "*" << 2*problem->refinement << "*" << problem->z_planes << std::endl;
    output << "  Quadrature (per octant)                       : " << problem->num_polar/2 << "*" << problem->num_azim/4 << std::endl;
    output << "  Angle Aggregation (1=single 2=polar, 3=octant): " << problem->ang_agg_type << std::endl;
    output << "  Number of groups per groupset                 : " << problem->group_per_groupset << std::endl;
    output << "  Number of groupsets                           : " << problem->num_groupsets << std::endl;
    output << "  Number of SMLs (Px,Py,Pz)                     : " << problem->num_SML << " (" << problem->num_cellsets[0]/problem->overload[0]<< "," << problem->num_cellsets[1]/problem->overload[1]<< "," << problem->num_cellsets[2]/problem->overload[2]<< ")" <<std::endl;
    output << "  Overload (Ox, Oy, Oz)                         : " << "(" << problem->overload[0]<< "," << problem->overload[1]<< "," << problem->overload[2]<< ")" <<std::endl;
    output << "  Number of cells per SML                       : " << cells_per_SML << std::endl;
    output << "  Number of unknowns per SML                    : " << unknowns_per_SML << std::endl;
    output << "  Number of unknowns per task                   : " << unknowns_per_task << std::endl;
    output << "================================================================" << std::endl;
    output << " " << std::endl;
  
  }
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
