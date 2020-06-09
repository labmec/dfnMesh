/*! 
 *  @brief     Describes a rectangular plane from four corner points.
 *  @details   
 *  @author    Pedro Lima
 *  @date      2019
 */

#ifndef DFNFracPlane_h
#define DFNFracPlane_h

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include <stdio.h>
#include <tuple>
#include "pzmatrix.h"
#include "pzcompel.h"
#include "pzgeoelbc.h"


typedef TPZFMatrix<REAL> Matrix;

/*! 
 *  @brief     Describes a triangular or quadrilateral plane from it's corner points.
 *  @details Enumeration of corner points should follow standard PZ topology, where 
 *  corner nodes are numbered counter-clockwise (clockwise should work as well) from
 *  zero to N. (This condition will automatically be met for triangles, but not 
 *  always for quadrilaterals)
 *  @author    Pedro Lima
 *  @date      2019
 */
class DFNFracPlane
{
  private:
	/// Contains fracture corner points. Matrix 3xn (n is the number of corners)
	Matrix fCornerPoints;
	
	/// Axis that define fracture orientation (Ax0 from node1 to node0, Ax1 from node1 to node2 and Ax2 the normal vector). Matrix 3x3
	Matrix fAxis ;

	/// Area of plane
	double fArea;

	/// Define a default tolerance
	REAL fTolerance = 1.e-4;

	/// If nodes of this plane have been added to a geometric mesh, this vector holds GeoNodes indices 
	TPZManVector<int64_t, 5> fPointsIndex;

  public:
	/// Empty constructor
	DFNFracPlane(){};

	/// Define plane from 3 to 4 corner points. Matrix should be 3xN (3 coordinates for each of the N corner points)
	DFNFracPlane(const Matrix &CornerPoints);

	/// Copy constructor
	DFNFracPlane(const DFNFracPlane &copy);

	/// Assignment operator
	DFNFracPlane &operator=(const DFNFracPlane &copy);

	/// Define corner coordinates
	void SetPlane(Matrix &CornerPoints);

	/// axis(i, j) returns component i of axis j
	REAL axis(int row, int col){return fAxis(row,col);}

	/// Get matrix with axis 0, 1 and 2 on each column
	Matrix axis() const {return fAxis;}

	/// Return corner coordinates
	Matrix GetCornersX() const;

	/// Return area of plane
	double area() const { return fArea; }

	/// Compute area of plane
	REAL ComputeArea();

	/// Return true if the rib intersects the plane
	bool Check_rib(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2);

	/// Return true if the rib intersects the plane
	bool Check_rib(const TPZGeoEl *rib);

	/// Return true if a point is above the fracture plane
   	bool Check_point_above(const TPZVec<REAL> &point) const;

   	/// Check whether the point coordinates are within the plane
   	bool IsPointInPlane(TPZVec<REAL> &point);
	
	/// Computes the intersection point with the plane
	TPZVec<REAL> CalculateIntersection(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2);
   

	/**
	 * @brief Inserts corner nodes in geometric mesh
	 * @param Pointer to geometric mesh
	 * @return Vector with nodes indices in gmesh
	 */
	TPZManVector<int64_t,4> SetPointsInGeomesh(TPZGeoMesh *gmesh);

	/**
	 * @brief Returns index of GeoNode that was created for corner i using DFNFracPlane::SetPointsInGeomesh(gmesh, matID)
	 * @param Local index of corner
	 */
	int64_t CornerIndex(int i) const{
		return fPointsIndex[i];
	}

//deprecated:
	///@brief Creates a geometric element for this plane in pointed mesh
	// int64_t DFNFracPlane::CreateElement(TPZGeoMesh *gmesh)

	// /// Element index
	// int64_t ElementIndex() const {return fFracIndex;}
  private:
	/// Initializes the datastructure of the object
	bool Check_Data_Consistency(Matrix CornerPoints);

};

#endif /* DFNFracPlane */