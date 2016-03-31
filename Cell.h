
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
	std::vector<double> direction;
  int face;
  int cs;
	// This will be used for Cellsets neighbors
	int SML;
  int local_id;

	Neighbor operator =(Neighbor a){
		id = a.id;
		direction = a.direction;
    face = a.face;
    cs = a.cs;
		SML = a.SML;
    local_id = a.local_id;
		return *this;
	}
};

class Cell
{
public:
	Cell(){ };
	~Cell(){ } ;

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
	std::vector<std::vector<double> > normals;
	std::vector<std::vector<double> > facecenters;
  std::vector<int> face_ids;
  
  // vector of vertex locations and the centroid of the cell;
  std::vector<std::vector<double> > vertices;
  std::vector<double> centroid;
  int num_faces, max_faces;

	// DFEM Matrices (Mass, Surface, and Gradient)
	std::vector<std::vector< double > > M;
	std::vector<std::vector<std::vector< std::vector<double> > > > N;
	std::vector<std::vector< std::vector<double> > > L;

	// Solution data structures (they're vectors of groupwise numbers
	std::vector<double> phi, current;

	// boundary the cell is on (0 if cell is interior) for each face
	std::vector<int> localboundary;
  std::vector<int> localijk;
  
  
  

	// This function builds the cell
	virtual void BuildCell(int, int, Problem*) = 0;

	// This function determines which cellset boundary
	// (if any) the cell is on
	virtual void SetLocalBoundary() = 0;

	// This function computes the global cell id
	// given the local cell id and the cellset id
	virtual void ComputeCellID(int, int, Problem*) = 0;

	// Given the global cell id and
	// the number of cells in the x, y, and z
	// this function computes the global ijk indices of the cell
	virtual void GetCellijk(int, int, int, int, std::vector<int>&) = 0;

	// Compute the face normals and face centers for the cell
	virtual void GetFaceNormals() = 0;
	virtual void GetFaceCenters() = 0;
  virtual void GetVertices(int) = 0;

	// This function computes the cell's 6) neighbor cells
	// and the direction those neighbors are upstream
	virtual void GetNeighbors(int) = 0;
  
  // This function returns in the incoming and outgoing faces, given a direction
  virtual void GetCellInOut(std::vector<double>, std::vector<int>& , std::vector<std::vector<int> >&) = 0;

	// This function computes the Mass, Surface, and Gradient matrix for the cell
	virtual void ComputeDFEMMatrices() = 0;

	// Get functions for required geometry data
	double GetSigmaTot(){ return sigma_t; };
	double GetSource(){ return source; };
	double GetDeltaX(){ return delta_x; };
	double GetDeltaY(){ return delta_y; };
	double GetDeltaZ(){ return delta_z; };

protected:
	int cells_per_cellset;
	
	int cells_x, cells_y, cells_z;
	std::vector<int> num_cellsets;

	// Geometry data
	double delta_x, delta_y, delta_z;

	double sigma_t;
	double source;
};

#endif
