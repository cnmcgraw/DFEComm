
#ifndef Cell_h
#define Cell_h 1

#include <vector>
#include "Problem.h"
#include "Quadrature.h"

//class Problem;

// This structure is a cell's neighbor
// It includes the neighbors local cell ID
// As well as the direction for which the neighbor
// is upstream.
struct Neighbor{
	int id;
	Direction direction;
	// This will be used for Cellsets neighbors
	int SML;
  int local_id;

	Neighbor operator =(Neighbor a){
		id = a.id;
		direction = a.direction;
		SML = a.SML;
		return *this;
	}
};

class Cell
{
public:
	Cell();
	~Cell();

public:
	// This is the global cell ID
	int CellID;

	// This is the local (to the cellset) cell ID
	int LocalCellID;

	// This is the vector of neighbor cells
	// If the cell is on a cellset boundary, the ID is 
	// negative, and is the ID of the boundary the face is on
	std::vector< Neighbor > neighbors;

	// vector of face normals and face centers (face ID is the vector element)
	std::vector<Direction> normals;
	std::vector<Direction> facecenters;

	// DFEM Matrices (Mass, Surface, and Gradient)
	std::vector< double > M;
	std::vector<std::vector<std::vector< Direction > > > N;
	std::vector<std::vector< Direction > > L;

	// Solution data structures (they're vectors of groupwise numbers
	std::vector<double> phi, current;

	// boundary the cell is on (0 if cell is interior) for each face
	std::vector<int> localboundary;

	// This function builds the cell
	void BuildCell(int, int, Problem*);

	// This function determines which cellset boundary
	// (if any) the cell is on
	void SetLocalBoundary();

	// This function computes the global cell id
	// given the local cell id and the cellset id
	void ComputeCellID(int, int, Problem*);

	// Given the global cell id and
	// the number of cells in the x, y, and z
	// this function computes the global ijk indices of the cell
	void GetCellijk(int, int, int, int, std::vector<int>&);

	// Compute the face normals and face centers for the cell
	void GetFaceNormals();
	void GetFaceCenters();

	// This function computes the cell's 4 (or 6) neighbor cells
	// and the direction those neighbors are upstream
	void GetNeighbors(int);

	// This function computes the Mass, Surface, and Gradient matrix for the cell
	void ComputeDFEMMatrices();

	// Get functions for required geometry data
	double GetSigmaTot(){ return sigma_t; };
	double GetSource(){ return source; };
	double GetDeltaX(){ return delta_x; };
	double GetDeltaY(){ return delta_y; };
	double GetDeltaZ(){ return delta_z; };

	std::vector<int> localijk;

private:
	int cells_per_cellset;
	
	int cells_x, cells_y, cells_z;
	std::vector<int> num_cellsets;

	// Geometry data
	// For 2D problems, delta_z = 1
	double delta_x, delta_y, delta_z;

	double sigma_t;
	double source;
};

#endif
