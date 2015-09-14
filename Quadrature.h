
#ifndef Quadrature_h
#define Quadrature_h 1

#include <vector>
#include "Problem_Input.h"
#include <iostream>

const double PI = 3.141592653589793238463;

struct Direction{
	double x;
	double y;
	double z;

	Direction operator =(Direction a){
		x = a.x;
		y = a.y;
		z = a.z;
		return *this;
	}

};

double dot(Direction, Direction);

Direction operator*(double, Direction);

struct Angleset
{
	int ID;
	std::vector<int> AngleIDs;
	std::vector<Direction> Omegas;
	std::vector<double> Weights;
	int octant;
	int angle_per_angleset;


};


class Quadrature
{
public:
	Quadrature();
	~Quadrature();

public:

	// vector of Anglesets
	std::vector< Angleset > Anglesets;

	// Number of anglesets
	int num_angleset;

	void BuildQuadrature(Problem_Input*);

	// Why do I need this....
	void SizeAngleIDs(std::vector<int>&);

	// Get the direction of the angleset given the angleset ID
	Direction GetOmega(int);

	// Octant of the angle given a direction
	int GetOctant(Direction);


private:
	int num_polar, num_azim, num_angles;

	int angle_per_angleset;

	std::vector<Direction> Omegas;
	std::vector<double> Weights;
	

	
	


};

#endif
