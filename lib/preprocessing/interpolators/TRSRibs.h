//
//  RibFrac.hpp
//  RibFrac
//
//  Created by JORGE ORDOÑEZ on 22/5/18.
//  Copyright © 2018 JORGE ORDOÑEZ. All rights reserved.
//

#ifndef TRSRibs_h
#define TRSRibs_h

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

class TRSRibs
{
private:
    TPZVec<TPZVec<double>> fcoords;
    TPZVec<double> fIntersectionPoint;
    TPZVec<TRSRibs> fchild;
    TPZVec<int> fNode_Index;
    double fLength;
    bool fQIsCut=false;

 
    
public:
    TRSRibs();
    TRSRibs(TPZVec<TPZVec<double>> coords);
    TRSRibs(TPZVec<TPZVec<double>> coords,TPZVec<double> IntersectionPoint);
    //@TODO crear destructor;
    //@TODO crear contructor de copia;
    void SetRibsCoords(TPZVec<TPZVec<double>> coords);
    TPZVec<TPZVec<double>> GetRibsCoords();
    void SetIntersectionPoint( TPZVec<double> ipoint);
    TPZVec<double> GetIntersectionPoint();
    void calc_child(TPZVec<TPZVec<double>> coords, TPZVec<double> ipoint);
    void get_child(TPZVec<TRSRibs> &child);
    void SetNodeIndexes(TPZVec<int> indexes);
    TPZVec<int> GetNodeIndexes();
    void calc_flength();
    
};

#endif /* TRSRibs_h */
