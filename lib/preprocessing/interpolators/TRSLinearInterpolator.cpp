//
//  InterpolPhil.cpp
//  Interpolator
//
//  Created by JOSE VILLEGAS on 22/5/18.
//  Copyright © 2018 JOSE VILLEGAS. All rights reserved.
//

#include "TRSLinearInterpolator.h"

TRSLinearInterpolator::TRSLinearInterpolator(){
    
}
TRSLinearInterpolator::TRSLinearInterpolator(Matrix data){
    fdata=data;
   
}

void TRSLinearInterpolator::SetData(Matrix data){
    if(data.Rows()>= 2){
        fdata = data;
       
    }
    else{
        DebugStop();
    }
}

double TRSLinearInterpolator::Val(double x){
    std::tuple<double, double> Intp = ValDeriv(x);

    return std::get<0>(Intp);
    
}

std::tuple<double, double> TRSLinearInterpolator::ValDeriv(double x){
    
    int npoints = fdata.Rows();
    double returned=0;
    double deriv = 0;
    
    //Linear Interpolation
    if (x >= fdata(0,0) && x<= fdata(npoints-1,0) ){
        for(int i=0;  i<npoints-1; i++){
            if(x >= fdata(i,0) && x<= fdata(i+1,0)){
                double x1 = fdata(i,0);
                double x2 = fdata(i+1,0);
                double y1 = fdata(i,1);
                double y2 = fdata(i+1,1);
                double a = (x - x2)/(x1 - x2);
                double b = (x - x1)/(x2 - x1);
                returned = a*y1 + b*y2;
                deriv = (y2-y1)/(x2-x1);
                break;
            }
        }
    }
    //
    //Extrapolation left
    if (x < fdata(0,0)){
        double x1 = fdata(0,0);
        double x2 = fdata(1,0);
        double y1 = fdata(0,1);
        double y2 = fdata(1,1);

        switch (fextLeft)
        {
            case Enone:
                returned = fvalLef;
                deriv = 0.0;
                break;
            case Constant:
                returned = fvalLef;
                deriv = 0.0;
                break;
            case Slope:
                returned = y1 - fvalLef*(x1 - x);
                deriv = fvalLef;
                break;
            case Linear:
                returned = y1 - ((y1 - x)*(y2-y1)/(x2-x1));
                deriv = (y2-y1)/(x2-x1);
                break;
        }
    }
    
    // Extrapolation Right
    if (x > fdata(npoints-1,0)){
        double x1 = fdata(npoints-2,0);
        double x2 = fdata(npoints-1,0);
        double y1 = fdata(npoints-2,1);
        double y2 = fdata(npoints-1,1);
        
        switch (fextRight) {
            case Enone:
                returned = fvalRight;
                deriv = 0.0;
                break;
            case Linear:
                returned = y2 + ((x - x2)*(y2-y1)/(x2-x1));
                deriv = (y2-y1)/(x2-x1);
                break;
            case Constant:
                returned = fvalRight;
                deriv =0.0;
                break;
            case Slope:
                returned = y2 + fvalRight*(x-x2);
                deriv = fvalRight;
            default:
                break;
        }
        
    }

    return {returned,deriv};
}
void TRSLinearInterpolator::SetLeftExtension(Extension left, double val){
    fextLeft = left;
    fvalLef =val;
}

void TRSLinearInterpolator::SetRightExtension(Extension right, double val){
    fextRight = right;
    fvalRight = val;
}
void TRSLinearInterpolator::Print (std::ostream &out){
    out<<"El type";
}
    

