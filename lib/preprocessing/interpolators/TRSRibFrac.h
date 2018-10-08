//
//  InterpolPhil.hpp
//  Interpolator
//
//  Created by JOSE VILLEGAS on 22/5/18.
//  Copyright © 2018 JOSE VILLEGAS. All rights reserved.
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
//#include "tpanic.h"

typedef TPZFMatrix<double> Matrix;

class TRSRibFrac
{
private:
    Matrix fdata;
    Matrix feixos;
    double ftolerance=0.0;
    
public:
    TRSRibFrac();
    TRSRibFrac(Matrix data);
    void SetPlane(Matrix plane);
    Matrix GetPlane();
    void SetTolerance(double tolerance);
    double GetTolerance();
    void Check_ConsistencyData(Matrix plane);
    bool Check_point_above(TPZVec<double> point);
    bool Check_point_below(TPZVec<double> point);
    bool Check_rib(TPZVec<double> p1, TPZVec<double> p2);
};

#endif /* TRSRibFrac_h */
