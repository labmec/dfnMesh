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
//#include "tpanic.h"

typedef TPZFMatrix<double> Matrix;

class TRSRibFrac
{
private:
    double ftolerance=0.0;
    
public:
    Matrix fdata;
    Matrix feixos ;
    
    TRSRibFrac();
    TRSRibFrac(Matrix data);
    void SetPlane(Matrix plane);
    Matrix GetPlane() const;
    void SetTolerance(double tolerance);
    double GetTolerance() const;
    void Check_ConsistencyData( Matrix plane) const;
    bool Check_point_above(TPZVec<double> point) const;
    bool Check_point_below(TPZVec<double> point);
    bool Check_rib(const TPZVec<double> &p1, TPZVec<double> p2) const;
};

#endif /* TRSRibFrac_h */
