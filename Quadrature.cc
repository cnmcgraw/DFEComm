#include "Quadrature.h"
#include <vector>
#include <math.h>
#include <iostream>

Quadrature::Quadrature()
{ }

Quadrature::~Quadrature()
{ }

void Quadrature::BuildQuadrature(Problem_Input* input)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  num_polar = input->num_polar;
  num_azim = input->num_azim;
  num_angles = num_polar*num_azim;

  // Resize all the vectors
  Omegas.resize(num_angles, std::vector<double>(3));
  Weights.resize(num_angles, 1./double(num_angles));

  // Since the answer doesn't matter, make Mu and Theta equally spaced
  double mu, theta;
  int id = 0;
  for(int i=0; i<num_azim; i++)
  {
    for(int j=0; j<num_polar/2; j++, id++)
    {
      mu = j*PI/num_polar - PI/2. + PI/(2.*num_polar);
      theta = i*2.*PI/num_azim + PI/(2.*num_azim);
      Omegas[id][0] = cos(mu)*cos(theta);
      Omegas[id][1] = cos(mu)*sin(theta);
      Omegas[id][2] = sin(mu);
    }
  }
  for (int i = 0; i<num_azim; i++)
  {
    for (int j = num_polar/2; j<num_polar; j++, id++)
    {
      mu = j*PI / num_polar + PI / (2. * num_polar);
      theta = i * 2. * PI / num_azim + PI / (2. * num_azim);
      Omegas[id][0] = cos(mu)*cos(theta);
      Omegas[id][1] = cos(mu)*sin(theta);
      Omegas[id][2] = sin(mu);
    }
  }

  // Single Angle Aggregation
    if(input->ang_agg_type == 1)
  {
    num_angleset = num_angles;
    angle_per_angleset = 1;
    Anglesets.resize(num_angleset);
    for (int i = 0; i < num_angles; i++)
    {
      Anglesets[i].ID = i;
      Anglesets[i].octant = GetOctant(Omegas[i]);
      Anglesets[i].angle_per_angleset = 1;
      Anglesets[i].AngleIDs.resize(angle_per_angleset, i);
      Anglesets[i].Weights.resize(angle_per_angleset, Weights[i]);
      Anglesets[i].Omegas.resize(angle_per_angleset, std::vector<double>(3));
      for(int j = 0; j < angle_per_angleset; j++)
      {
        Anglesets[i].Omegas[j][0] = Omegas[i*angle_per_angleset+j][0];
        Anglesets[i].Omegas[j][1] = Omegas[i*angle_per_angleset+j][1];
        Anglesets[i].Omegas[j][2] = Omegas[i*angle_per_angleset+j][2];
      }
    }

  }
  // Polar Aggregation
  else if(input->ang_agg_type == 2)
  {
    // We need to separate the +z and -z polar angles
    angle_per_angleset = num_polar/2;
    num_angleset = num_azim*2;
    Anglesets.resize(num_angleset);
    for (int i = 0; i < num_azim; i++)
    {
      Anglesets[2 * i].ID = 2 * i;
      Anglesets[2 * i].octant = GetOctant(Omegas[(2 * i + 1)*angle_per_angleset - 1]);
      Anglesets[2 * i].angle_per_angleset = angle_per_angleset;

      Anglesets[2 * i].AngleIDs.resize(angle_per_angleset);
      Anglesets[2 * i].Omegas.resize(angle_per_angleset, std::vector<double>(3));
      Anglesets[2 * i].Weights.resize(angle_per_angleset);

      for (int j = 0; j < num_polar / 2; j++)
      {
        Anglesets[2 * i].AngleIDs[j] = i*num_polar + j;
        Anglesets[2 * i].Omegas[j][0] = Omegas[i*num_polar + j][0];
        Anglesets[2 * i].Omegas[j][1] = Omegas[i*num_polar + j][1];
        Anglesets[2 * i].Omegas[j][2] = Omegas[i*num_polar + j][2];
        Anglesets[2 * i].Weights[j] = Weights[i*num_polar + j];
      }

      Anglesets[2 * i+1].ID = 2 * i + 1;
      Anglesets[2 * i+1].octant = GetOctant(Omegas[(2 * i + 2)*angle_per_angleset - 1]);
      Anglesets[2 * i + 1].angle_per_angleset = angle_per_angleset;

      Anglesets[2 * i+1].AngleIDs.resize(angle_per_angleset);
      Anglesets[2 * i+1].Omegas.resize(angle_per_angleset, std::vector<double>(3));
      Anglesets[2 * i+1].Weights.resize(angle_per_angleset);

      for (int j = num_polar / 2; j < num_polar; j++)
      {
        Anglesets[2 * i].AngleIDs[j  -num_polar / 2] = i*num_polar + j;
        Anglesets[2 * i].Omegas[j - num_polar / 2][0] = Omegas[i*num_polar + j][0];
        Anglesets[2 * i].Omegas[j - num_polar / 2][1] = Omegas[i*num_polar + j][1];
        Anglesets[2 * i].Omegas[j - num_polar / 2][2] = Omegas[i*num_polar + j][2];
        Anglesets[2 * i].Weights[j - num_polar / 2] = Weights[i*num_polar + j];
      }
    }
  }
  //  Octant Aggregation
  else if(input->ang_agg_type == 3)
  {
    angle_per_angleset = num_angles/8;
    num_angleset = 8;
    Anglesets.resize(num_angleset);
    for(int i = 0; i < 8; i++)
    {
      Anglesets[i].ID = i;
      Anglesets[i].angle_per_angleset = angle_per_angleset;
      Anglesets[i].AngleIDs.resize(angle_per_angleset);
      Anglesets[i].Omegas.resize(angle_per_angleset, std::vector<double>(3));
      Anglesets[i].Weights.resize(angle_per_angleset);
      Anglesets[i].octant = GetOctant(Omegas[i*angle_per_angleset]);
      for(int j = 0; j<angle_per_angleset; j++)
      {
        Anglesets[i].AngleIDs[j] = i*angle_per_angleset + j;
        Anglesets[i].Omegas[j][0] = Omegas[i*angle_per_angleset + j][0];
        Anglesets[i].Omegas[j][1] = Omegas[i*angle_per_angleset + j][1];
        Anglesets[i].Omegas[j][2] = Omegas[i*angle_per_angleset + j][2];
        Anglesets[i].Weights[j] = Weights[i*angle_per_angleset + j];
      }
    }

  }
  else
  {
    if (rank == 0){ std::cout << "Invalid Angle Aggregation Type" << std::endl; }
  }
}


void Quadrature::SizeAngleIDs(std::vector<int>& AngleIDs)
{
  AngleIDs.resize(angle_per_angleset);
}

std::vector<double> Quadrature::GetOmega(int AngleSetID)
{
  int my_octant = Anglesets[AngleSetID].octant;
  std::vector<double> my_Omega(3,0);

  if (my_octant == 0)
  {
    my_Omega[0] = 1;
    my_Omega[1] = 1;
    my_Omega[2] = 1;
  }
  else if (my_octant == 1)
  {
    my_Omega[0] = -1;
    my_Omega[1] = 1;
    my_Omega[2] = 1;
  }
  else if (my_octant == 2)
  {
    my_Omega[0] = -1;
    my_Omega[1] = -1;
    my_Omega[2] = 1;
  }
  else if (my_octant == 3)
  {
    my_Omega[0] = 1;
    my_Omega[1] = -1;
    my_Omega[2] = 1;
  }
  else if (my_octant == 4)
  {
    my_Omega[0] = 1;
    my_Omega[1] = 1;
    my_Omega[2] = -1;
  }
  else if (my_octant == 5)
  {
    my_Omega[0] = -1;
    my_Omega[1] = 1;
    my_Omega[2] = -1;
  }
  else if (my_octant == 6)
  {
    my_Omega[0] = -1;
    my_Omega[1] = -1;
    my_Omega[2] = -1;
  }
  else if (my_octant == 7)
  {
    my_Omega[0] = 1;
    my_Omega[1] = -1;
    my_Omega[2] = -1;
  }
  return my_Omega;
}

int Quadrature::GetOctant(std::vector<double> Omega)
{
  int my_octant = 0;

  if (Omega[0] > 0)
  {
    if (Omega[1] > 0)
    {
      if (Omega[2] > 0)
        my_octant = 0;
      else
        my_octant = 4;
    }
    else
    {
      if (Omega[2] > 0)
        my_octant = 3;
      else
        my_octant = 7;
    }
  }
  else
  {
    if (Omega[1] > 0)
    {
      if (Omega[2] > 0)
        my_octant = 1;
      else
        my_octant = 5;
    }
    else
    {
      if (Omega[2] > 0)
        my_octant = 2;
      else
        my_octant = 6;
    }
  }

  return my_octant;

}

 double dot(std::vector<double> a, std::vector<double> b)
 {
   double mydot = 0;
   for(int i = 0; i < a.size(); i++)
      mydot += a[i]*b[i];
   return mydot;
 }

// Direction operator*(double a, Direction omega)
// {
  // omega.x *= a;
  // omega.y *= a;
  // omega.z *= a;

  // return omega;
// }
