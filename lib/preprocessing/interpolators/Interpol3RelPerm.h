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
    REAL fswc = -1.0;
    REAL fsorw = -1.0;
    REAL fsorg = -1.0;
public:
    
    
    enum Kralpha {EKrw, EKrow ,EKrg, EKrog};
    enum ModelInterpol {MStoneI, MStoneII};
    ModelInterpol fKroModel;
    
    //Set Type
    Interpol3RelPerm(){
        fKroModel = ModelInterpol::MStoneI;
    }
  
    void SetData(Matrix data, Kralpha alpa);
    void SetParam(REAL swc, REAL sorw, REAL sorg);
    void SetKroModel(ModelInterpol model);
    double Val(double Sw, double Sg);
    Matrix Deriv(double Sw, double Sg);
    
    std::function<double(double)> GetFunction();
    std::function<double(double)> GetFunctionDeriv();
    void ReadData(std::string data, Kralpha alpa);
    
};

#endif /* Interpol3RelPerm_h */
