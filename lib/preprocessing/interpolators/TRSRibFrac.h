//
//  RibFrac.hpp
//  RibFrac
//
//  Created by JORGE ORDOÑEZ on 22/5/18.
//  Copyright © 2018 JORGE ORDOÑEZ. All rights reserved.
//

#ifndef TRSRibFrac_h
#define TRSRibFrac_h

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include <stdio.h>
#include <tuple>
#include "pzmatrix.h"
#include "pzcompel.h"
#include "pzgeoelbc.h"
#include "TRSRibs.h"
#include "TRSFace.h"

//#include "tpanic.h"

typedef TPZFMatrix<REAL> Matrix;

class TRSRibFrac
{
private:
    /// Define a default tolerance
    REAL fTolerance= 1.e-6;
    
    /// Map of relevant intersecting ribs
    std::map<int64_t,TRSRibs> fRibs;
    
    /// Map of relevant intersecting ribs
    std::map<int64_t,TRSFace> fFaces;

    /// Contains fracture corner points. Matrix 3xn (n is the number of corner points)
    Matrix fCornerCoordinates;
    
    /// Contains fracture axis Ax0,Ax1 and Ax2. Matrix 3x3
    Matrix fAxis ;
    
    /// The coordinates of the center point (should be 3 coordinates x y z )
    TPZManVector<REAL,3> fCenterCo;
    
    /// Pointer for the geometric mesh
    TPZGeoMesh *fGMesh;

public:
    
    /// Empty constructor
    TRSRibFrac();
    
    /// Define the fracture plane by 4 points
    /// Points should be colinear and define a square
    /// The matrix should be dimension 3x4, each column defining the coordinates
    /// of a point
    TRSRibFrac(const Matrix &data, TPZGeoMesh *gmesh);
    
    /// Copy constructor
    TRSRibFrac(const TRSRibFrac &copy);
    
    /// Assignment operator
    TRSRibFrac &operator=(const TRSRibFrac &copy);

    /// Define the corner coordinates of the fracture
    void SetPlane(Matrix plane);
    
    /// Associate the geometric mesh
    void SetgeoMesh(TPZGeoMesh *gmesh){
        fGMesh=gmesh;
    }
    
    /// Access the geomesh
    TPZGeoMesh* GetgeoMesh(){
        return fGMesh;
    }
    
    /// Return the corner nodes of the fracture
    Matrix GetPlane() const;
    
    /// Modify the default tolerance
    void SetTolerance(REAL tolerance);
    REAL GetTolerance() const;
    
private:
    /// Initializes the datastructure of the object
    bool Check_ConsistencyData(Matrix CornerCoordinates);
    
public:
    /// Return true if a point is above the fracture plane
    bool Check_point_above(const TPZVec<REAL> &point) const;
    
    /// Return true if the rib intersects the plane
    bool Check_rib(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2) const;
    
    ///Return true if the surface needs to be divided
    bool NeedsSurface_Divide(int64_t suface_index, TPZVec<int64_t> interribs) ;
    
    
private:
    
    /// Return true if the intersection point between ribs is within the plane
    bool RibInPlane(TPZVec<REAL> point);
    
    /// Checks the neighbour dimension and return if it is equal
    bool HasEqualDimensionNeighbour(TPZGeoElSide &gelside);
   
public:
    
    /// Insert element in the geometric mesh of lower dimension
    void CreateSkeletonElements(int dimension, int matid);
    
    /// Compute the intersection point with the plane
    TPZVec<REAL> CalculateIntersection(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2);
    
    /// Divide the one dimensional element by the intersection with the plane
    TRSRibs DivideRib(int element_index);
    
   /// Verify if the one dimensional element intersects the plane
    bool CheckElementIntersection(int64_t elindex);
    
   /// Check whether the point coordinates are within the plane
    bool IsPointInPlane(TPZVec<REAL> &point);
    
   /// Access the ribs data structure
    void AddRib(TRSRibs rib);
    
    /// Set a cut face
    void AddFace(TRSFace face);
    
    /// Mao of ribs
    std::map<int64_t ,TRSRibs> GetRibs();
    
    /// Create the children surfaces (not implemented yet)
    void CreateSurfaces(int matID);
    
    TRSRibs *Rib(int index);
    
};

#endif /* TRSRibFrac_h */
