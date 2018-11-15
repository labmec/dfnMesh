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
#include "TRSLinearInterpolator.h"
//#include "tpanic.h"

typedef TPZFMatrix<double> Matrix;

class Interpol3RelPerm
{
private:
    TRSLinearInterpolator Krw;
    TRSLinearInterpolator Krow;
    TRSLinearInterpolator Krg;
    TRSLinearInterpolator Krog;
    
public:
    //Set Type
    Interpol3RelPerm(){
        
    }
    enum Kralpha {EKrw, EKrow ,EKrg, EKrog};
    enum ModelInterpol {MStoneI, MStoneII};
    
    void SetData(Matrix data, Kralpha alpa);
    double Val(double Sw, double Sg);
    double ValDeriv(double Sw, double Sg);
    
    std::function<double(double)> GetFunction();
    std::function<double(double)> GetFunctionDeriv();
    void ReadData(std::string data);
    
};

#endif /* Interpol3RelPerm_h */
