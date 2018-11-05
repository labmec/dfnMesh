//
//  InterpolPhil.hpp
//  Interpolator
//
//  Created by JOSE VILLEGAS on 22/5/18.
//  Copyright © 2018 JOSE VILLEGAS. All rights reserved.
//

#ifndef Interpol3RelPerm_h
#define Interpol3RelPerm_h

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include <stdio.h>
#include <tuple>
#include "pzmatrix.h"
//#include "tpanic.h"

typedef TPZFMatrix<double> Matrix;

class Interpol3RelPerm
{
private:
     Matrix fdata;
    typedef std::function<std::tuple<double, double>(double)> functions;
    functions Krw;
    functions Krow;
    functions Krg;
    functions Krog;
public:
    //Set Type
    enum Kralpha {EKrw, EKrow ,EKrg, EKrog};
    enum ModelInterpol {MStoneI, MStoneII};
    void SetData(std::function<std::tuple<double, double>(double)> function, Kralpha alpa);
    double Val(double Sw, double Sg);
    double ValDeriv(double Sw, double Sg);
    
    std::function<double(double)> GetFunction();
    
   functions GetFunctionDeriv();
    
};

#endif /* Interpol3RelPerm_h */
