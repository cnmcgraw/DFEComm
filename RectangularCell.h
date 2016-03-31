
#ifndef RectangularCell_h
#define RectangularCell_h 1

#include <vector>
#include "Cell.h"


class RectangularCell: public Cell
{
public:
	RectangularCell();
	~RectangularCell();

public:

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
  void GetVertices(int);

	// This function computes the cell's 6) neighbor cells
	// and the direction those neighbors are upstream
	void GetNeighbors(int);
  
  // This function returns in the incoming and outgoing faces, given a direction
  void GetCellInOut(std::vector<double>, std::vector<int>& , std::vector<std::vector<int> >&);

	// This function computes the Mass, Surface, and Gradient matrix for the cell
	void ComputeDFEMMatrices();

};

#endif
