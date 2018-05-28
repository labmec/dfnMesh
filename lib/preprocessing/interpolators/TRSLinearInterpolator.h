//
//  InterpolPhil.hpp
//  Interpolator
//
//  Created by JOSE VILLEGAS on 22/5/18.
//  Copyright © 2018 JOSE VILLEGAS. All rights reserved.
//

#ifndef TRSLinearInterpolator_h
#define TRSLinearInterpolator_h

#include <stdio.h>
#include <tuple>
#include "TMatrix.h"
#include "tpanic.h"

typedef TMatrix Matrix;
class TRSLinearInterpolator
{
private:
    Matrix fdata;

public:
    //Extrapolation types
    enum Extension {Enone, Linear ,Constant,Slope};
    Extension fextLeft=Enone;
    double fvalLef=0.0;
    Extension fextRight=Enone;
    double fvalRight=0.0;
    
    //constructors
    TRSLinearInterpolator();
    TRSLinearInterpolator(Matrix data);
    
    //set interpolate data
    void SetData(Matrix data);
    
    //Return the interpolate val
    double Val(double x);
   //Return the interpolate val and derivate
    std::tuple<double, double> ValDeriv(double x);
    //Set the extrapolation: left type
    void SetLeftExtension(Extension left, double val=0.);
    //Set the extrapolation: Right type
    void SetRightExtension(Extension right, double val=0.);
    //print the interpolation function
    void Print (std::ostream &out);

};
#endif /* TRSLinearInterpolator_h */
