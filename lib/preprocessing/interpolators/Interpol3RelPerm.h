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
   // Linear Interpolators functions for Biphasic oil-water system
    TRSLinearInterpolator Krw;
    TRSLinearInterpolator Krow;
    // Linear Interpolators functions for Biphasic oil-gas system
    TRSLinearInterpolator Krg;
    TRSLinearInterpolator Krog;
    
    //EndPoints saturations (if those points are not set, we will take the first points of the data table)
    REAL fswc = -1.0;
    REAL fsorw = -1.0;
    REAL fsorg = -1.0;
    
public:
    
    //Interpolator Type
    enum Kralpha {EKrw, EKrow ,EKrg, EKrog};
    
    //Interpolation Kro Models
    enum ModelInterpol {MStoneI, MStoneII};
    
    ModelInterpol fKroModel;
    
    
    /** @brief Default constructor */
    Interpol3RelPerm(){
        fKroModel = ModelInterpol::MStoneI;
    }
  
    /** @brief Function that reads the interpolation data from a ".txt" file
     * @param data is the name of the file with the interpolation data
    * @param alpa is the name of the file with the interpolation data
     */
    
       void ReadData(std::string data, Kralpha alpa);
    /** @brief Function that sets the matrix((x,f(x)) for linear interpolation or (x,f(x), f'(x)) for Hermite interpolation) with the data to interpolate
     * @param data is a matrix (nx2 for Linear interpolation or nx3 for Hermite interpolation) with one dimensional data to interpolate
     * @param alpa selects the interpolator in which the data matrix is set
     */
    void SetData(Matrix data, Kralpha alpa);
    
    /** @brief Function that sets the interpolator parameters
     * @param swc is the connate water saturation
     * @param sorw is the residual oil saturation on  water phase
     * @param sorg is the residual oil saturation on  gas phase
     */
    void SetParam(REAL swc, REAL sorw, REAL sorg);
    
    /** @brief Function that sets the interpolator model
     * @param model is the interpolator model, you could use StoneI and StoneII model
     */
    void SetKroModel(ModelInterpol model);
    
    /** @brief Function that calculates the value of the interpolated function corresponding to a "Sw" and "Sg" point
     * @param Sw is the water saturantion
     * @param Sg is the gas saturantion
     * @return  kro value evaluated in Sw, Sg.
     */
    double Val(double Sw, double Sg);
    
    /** @brief Function that calculates the value of the  derivative corresponding to "Sw", "Sg" point
     * @param Sw is the water saturantion
     * @param Sg is the gas saturantion
     * @return A matrix with all derivatives.
     (0,0)=dkrw/dsw
     (0,1)=dkro/dsw
     (1,1)=dkro/dsg
     (1,2)=dkrg/dsg
     */
    Matrix Deriv(double Sw, double Sg);
    
    /** @brief Function that transforms the "Interpol3RelPerm" object as an function that returns the interpolated values
     *  @param Function that returns the interpolated value
     *  @return An interpolation function
     */
    std::function<double(double,double)> GetFunction();
    
    /** @brief Function that transforms the "Interpol3RelPerm" object as an function that returns derivatives
     *  @param Function that returns the interpolated value and its derivatives
     *  @return An interpolation function that returns the interpolated value and its derivatives
     */
    std::function<Matrix(double,double)> GetFunctionDeriv();

};

#endif /* Interpol3RelPerm_h */
