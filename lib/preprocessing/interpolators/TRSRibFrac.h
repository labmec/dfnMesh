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
//#include "tpanic.h"

typedef TPZFMatrix<double> Matrix;

class TRSRibFrac
{
private:
    double ftolerance=0.0;
    
public:
   mutable Matrix fdata;
   mutable Matrix faxis ;
    TPZGeoMesh *fmesh;
    
    TRSRibFrac();
    TRSRibFrac(Matrix data);
    void SetPlane(Matrix plane);
    void SetgeoMesh(TPZGeoMesh *gmesh){
        fmesh=gmesh;
    }
    TPZGeoMesh* GetgeoMesh(){
        return fmesh;
    }
    Matrix GetPlane() const;
    void SetTolerance(double tolerance);
    double GetTolerance() const;
    void Check_ConsistencyData( Matrix plane) ;
    bool Check_point_above(TPZVec<double> point) const;
    bool Check_point_below(TPZVec<double> point) const;
    bool Check_rib( TPZVec<double> p1, TPZVec<double> p2) const;
    bool HasLowerDimensionNeighbour(TPZGeoElSide &gelside);
    void CreateSkeleton(int dimension, int matid);
};

#endif /* TRSRibFrac_h */
