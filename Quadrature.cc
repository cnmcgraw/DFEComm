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
	num_polar = input->num_polar;
	num_azim = input->num_azim;
	num_angles = num_polar*num_azim;

	// Resize all the vectors

	Omegas.resize(num_angles);
	Weights.resize(num_angles, 1./num_angles);

	// Since the answer doesn't matter, make Mu and Theta equally spaced
	double mu, theta;
	int id = 0;
	for(int i=0; i<num_azim; i++)
	{
		for(int j=0; j<num_polar/2; j++, id++)
		{
			mu = j*PI/num_polar - PI/2 + PI/(2*num_polar);
			theta = i*2*PI/num_azim + PI/(2*num_azim);
			Omegas[id].x = cos(mu)*cos(theta);
			Omegas[id].y = cos(mu)*sin(theta);
			Omegas[id].z = sin(mu);
		}
	}
	for (int i = 0; i<num_azim; i++)
	{
		for (int j = num_polar/2; j<num_polar; j++, id++)
		{
			mu = j*PI / num_polar + PI / (2 * num_polar);
			theta = i * 2 * PI / num_azim + PI / (2 * num_azim);
			Omegas[id].x = cos(mu)*cos(theta);
			Omegas[id].y = cos(mu)*sin(theta);
			Omegas[id].z = sin(mu);
		}
	}

	// Single Angle Aggregation
    if(input->ang_agg_type == 1)
	{
		num_angleset = num_angles;
		Anglesets.resize(num_angleset);
		for (int i = 0; i < num_angles; i++)
		{
			Anglesets[i].ID = i;
			Anglesets[i].octant = GetOctant(Omegas[i]);
			Anglesets[i].angle_per_angleset = 1;
			Anglesets[i].AngleIDs.resize(angle_per_angleset, i);
			Anglesets[i].Weights.resize(angle_per_angleset, Weights[i]);
			Anglesets[i].Omegas.resize(angle_per_angleset, Omegas[i]);
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
			Anglesets[2 * i].Omegas.resize(angle_per_angleset);
			Anglesets[2 * i].Weights.resize(angle_per_angleset);

			for (int j = 0; j < num_polar / 2; j++)
			{
				Anglesets[2 * i].AngleIDs[j] = i*num_polar + j;
				Anglesets[2 * i].Omegas[j] = Omegas[i*num_polar + j];
				Anglesets[2 * i].Weights[j] = Weights[i*num_polar + j];
			}

			Anglesets[2 * i+1].ID = 2 * i + 1;
			Anglesets[2 * i+1].octant = GetOctant(Omegas[(2 * i + 2)*angle_per_angleset - 1]);
			Anglesets[2 * i + 1].angle_per_angleset = angle_per_angleset;

			for (int j = num_polar / 2; j < num_polar; j++)
			{
				Anglesets[2 * i].AngleIDs[j] = i*num_polar + j;
				Anglesets[2 * i].Omegas[j] = Omegas[i*num_polar + j];
				Anglesets[2 * i].Weights[j] = Weights[i*num_polar + j];
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
			Anglesets[i].Omegas.resize(angle_per_angleset);
			Anglesets[i].Weights.resize(angle_per_angleset);
			Anglesets[i].octant = GetOctant(Omegas[i*angle_per_angleset]);
			for(int j = 0; j<angle_per_angleset; j++)
			{
				Anglesets[i].AngleIDs[j] = i*angle_per_angleset + j;
				Anglesets[i].Omegas[j] = Omegas[i*angle_per_angleset + j];
				Anglesets[i].Weights[j] = Weights[i*angle_per_angleset + j];
			}
		}
		
	}
	else
		std::cout << "Invalid Angle Aggregation Type" << std::endl;
 
}


void Quadrature::SizeAngleIDs(std::vector<int>& AngleIDs)
{
	AngleIDs.resize(angle_per_angleset);
}

Direction Quadrature::GetOmega(int AngleSetID)
{
	int my_octant = Anglesets[AngleSetID].octant;
	Direction my_Omega;

	if (my_octant == 0)
	{
		my_Omega.x = 1;
		my_Omega.y = 1;
		my_Omega.z = 1;
	}
	else if (my_octant == 1)
	{
		my_Omega.x = -1;
		my_Omega.y = 1;
		my_Omega.z = 1;
	}
	else if (my_octant == 2)
	{
		my_Omega.x = -1;
		my_Omega.y = -1;
		my_Omega.z = 1;
	}
	else if (my_octant == 3)
	{
		my_Omega.x = 1;
		my_Omega.y = -1;
		my_Omega.z = 1;
	}
	else if (my_octant == 4)
	{
		my_Omega.x = 1;
		my_Omega.y = 1;
		my_Omega.z = -1;
	}
	else if (my_octant == 5)
	{
		my_Omega.x = -1;
		my_Omega.y = 1;
		my_Omega.z = -1;
	}
	else if (my_octant == 6)
	{
		my_Omega.x = -1;
		my_Omega.y = -1;
		my_Omega.z = -1;
	}
	else if (my_octant == 7)
	{
		my_Omega.x = 1;
		my_Omega.y = -1;
		my_Omega.z = -1;
	}
	return my_Omega;
}

int Quadrature::GetOctant(Direction Omega)
{
	int my_octant = 0;

	if (Omega.x > 0)
	{
		if (Omega.y > 0)
		{
			if (Omega.z > 0)
				my_octant = 0;
			else
				my_octant = 4;
		}
		else
		{
			if (Omega.z > 0)
				my_octant = 3;
			else
				my_octant = 7;
		}
	}
	else
	{
		if (Omega.y > 0)
		{
			if (Omega.z > 0)
				my_octant = 1;
			else
				my_octant = 5;
		}
		else
		{
			if (Omega.z > 0)
				my_octant = 2;
			else
				my_octant = 6;
		}
	}

	return my_octant;
	
}

double dot(Direction a, Direction b)
{
	double mydot = a.x*b.x + a.y*b.y + a.z*b.z;
	return mydot;
}

Direction operator*(double a, Direction omega)
{
	omega.x *= a;
	omega.y *= a;
	omega.z *= a;

	return omega;
}

//double ComputeAngles()
//{
//  // Compute the total number of directions
//  int numAngles = input.num_polar * input.num_azim;
//  qst.resize(numAngles);
//  for(size_t i=0; i<qst.size(); ++i)
//    qst[i].resize(4);  // 0 = weight, 1-3 = direction cosines
//
//  // Order of the Legendre polynomial
//  int nRoots = num_polar;
//  // Temporary storage for the roots and weights (polar data)
//  vector<double> Lxi(nRoots,0.0);
//  vector<double> Lw(nRoots,0.0);
//
//  // Compute roots and weights (polar data)
//  eval_roots_legendre(nRoots, norm, Lxi, Lw);
//
//  // Temporary data for the azimuthal data
//  double azim_angle, azim_weight;
//  double costheta, sintheta = 0.0;
//  double _nazimuth = (double) nazimuth;
//  // Represents the starting angle (delta_phi = half of the angle between azimuthal angles)
//  double delta_phi = pi()/(4.*_nazimuth);
//  // Constant weight per level for all levels
//  azim_weight = 2.0*delta_phi; // = PI/(2.*_nazimuth) !!
//
//  int iord = 0;
//
//  // Starting azimuthal angle for any polar  level
//  azim_angle = delta_phi;
//  // Angle between two azimuthal angles = 2.*delta_phi
//
//  for (int i=0; i<num_azim; ++i)
//  {
//    for (int j=num_polar/2; j<nRoots; ++j)
//    {
//      costheta = Lxi[j];
//      sintheta = std::sqrt(1.0-std::pow(costheta,2));
//      iord = j + (i-1)*num_polar/2;                     // increment direction index
//      qst[iord][MU]  = cos(azim_angle) * sintheta; // cos(angle_azim)*sin(angle_polar)
//      qst[iord][ETA] = sin(azim_angle) * sintheta; // sin(angle_azim)*sin(angle_polar)
//      qst[iord][XI]  = Lxi[j];                     // cos(angle_polar)
//      qst[iord][WT]  = azim_weight * Lw[j];        // global weight = Chebyshev*Gauss_legendre
//    }
//    azim_angle += 2.*delta_phi;
//  }
//}
//
//void Quadrature::eval_roots_legendre()
//{
//  int nRoots = input.num_polar;
//
//  GLRoots.resize(nRoots);
//  GLWeights.resize(nRoots);
//
//  double x1 = -1.0;
//  double x2 = 1.0;
//
//  int m = (nRoots+1)/2;
//  double xm = 0.5 * (x2+x1);
//  double xl = 0.5 * (x2-x1);
//
//  for(int i=1; i<=m; ++i)
//  {
//    double z1, pp, z = cos(PI * (i-0.25) / (nRoots+0.5));
//    do
//    {
//      double p1 = 1.0, p2 = 0.0;
//      for(int j=1; j<=nRoots; j++)
//      {
//        real8 p3 = p2;
//        p2 = p1;
//        p1 = ((2.0*j-1.0)*z*p2 - (j-1.0)*p3) / j;
//      }
//      pp = nRoots * (z*p1 - p2) / (z*z - 1.0);
//      z1 = z;
//      z = z1 - p1/pp;
//    }
//    while (fabs(z - z1) > TOL);
//
//    GLRoots[i-1] = xm - xl * z;
//    GLRoots[nRoots-i] = xm + xl * z;
//    GLWeights[i-1] = norm * xl / ((1.0 - z*z)*pp*pp);
//    GLWeights[nRoots-i] =  GLWeights[i-1];
//  }
//}