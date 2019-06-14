/*! 
 *  @brief     Describes a rectangular plane from four corner points.
 *  @details   
 *  @author    Pedro Lima
 *  @date      2019
 */

#ifndef TRSFracPlane_h
#define TRSFracPlane_h

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
 *  @brief     Describes a rectangular plane from four corner points.
 *  @details   
 *  @author    Pedro Lima
 *  @date      2019
 */
class TRSFracPlane
{

  private:
    /// Contains fracture corner points. Matrix 3xn (n is the number of corners)
    Matrix fCornerPoints;
    
    /// Contains fracture axis Ax0,Ax1 and Ax2. Matrix 3x3
    Matrix fAxis ;
    
    /// The coordinates of the center point (should be 3 coordinates x y z )
    TPZManVector<REAL,3> fCenterCo;
    
    // /// Pointer for a geometric mesh
    // TPZGeoMesh *fGMesh;

    // /// Contains fracture edges midpoints. Matrix 3xn (n is the number of edges)
    // Matrix fMidPoints;

    /// Distance between edge midpoints parallel to Ax0
    double fL0;

    /// Distance between edge midpoints parallel to Ax1
    double fL1;

    /// Define a default tolerance
    REAL fTolerance = 1.e-6;
  public:
    /// Empty constructor
    TRSFracPlane(){};

    /// Define plane from 4 corner points. Matrix should be 3x4 (3 coordinates for each of the 4 corner points)
    TRSFracPlane(const Matrix &CornerPoints);

    /// Copy constructor
    TRSFracPlane(const TRSFracPlane &copy);

    /// Assignment operator
    TRSFracPlane &operator=(const TRSFracPlane &copy);

    // /// Define corner coordinates
    // void SetPlane(Matrix &CornerPoints);

    /// Get axis
    REAL axis(int row, int col){return fAxis(row,col);}

    /// Get matrix with axis 0, 1 and 2 on each column
    Matrix axis() const {return fAxis;}

    /// Get center coordinates
    TPZManVector<REAL, 3> center() const {return fCenterCo;}

    /// Get plane's dimension parallel to Ax0
    double L0() {return fL0;}

    /// Get plane's dimension parallel to Ax1
    double L1(){return fL1;}

    /// Return corner coordinates
    Matrix GetCorners() const;

  private:
    /// Initializes the datastructure of the object
    bool Check_Data_Consistency(Matrix CornerPoints);
};

#endif /* TRSFracPlane */