//
//  InterpolPhil.hpp
//  Interpolator
//
//  Created by JOSE VILLEGAS on 22/5/18.
//  Copyright © 2018 JOSE VILLEGAS. All rights reserved.
//

#ifndef TRSLinearInterpolator_h
#define TRSLinearInterpolator_h

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include <stdio.h>
#include <tuple>
#include "pzmatrix.h"
//#include "tpanic.h"

typedef TPZFMatrix<double> Matrix;

class TRSLinearInterpolator
{
private:
     Matrix fdata;
public:
   

    //Extrapolation types
    enum Extension {Enone, ELinear ,EConstant,ESlope};
    //Interpolation methods
    enum InterpType {TLinear, THermite};
    //Print format
    enum MatrixOutputFormat {EFormatted, MathemathicaInput};
    
    //Default left extension
    Extension fextLeft=Enone;
    //Default left extension value to extrapolate
    double fvalLef=0.0;
    
     //Default right extension
    Extension fextRight=Enone;
    
    double fvalRight=0.0;
     //Default right extension value to extrapolate
    InterpType fInterType = TLinear;
    
    /** @brief Default constructor */
    TRSLinearInterpolator();
    
    /** @brief Constructor with data to interpolate
     * @param data is a matrix(nx2 for linear interpolation or nx3 for Hermite interpolation) with one dimensional data to interpolate
     */
    TRSLinearInterpolator(Matrix data);
    
    /** @brief Function that sets the matrix((x,f(x)) for linear interpolation or (x,f(x), f'(x)) for Hermite interpolation) with the data to interpolate
     * @param data is a matrix(nx2 for Linear interpolation or nx3 for Hermite interpolation) with one dimensional data to interpolate
     */
    void SetData(Matrix data);
    
    /** @brief Function that returns the data
     * @return  The matrix(nx2 for Linear interpolation or nx3 for Hermite interpolation) with the one dimensional data to interpolate
     */
    Matrix GetData();
    
    /** @brief Function that read the interpolation data from a ".txt" file
     * @param data is the name of the file with the interpolation data
     */
    void ReadData(std::string name);
    
    /** @brief Function that calculates the value of the interpolated function corresponding to a "x" point
     * @param x is the interpolate point
     * @return  The function evaluated at the interpolation point x.
     */
    double Val(double x);
    
    
    /** @brief Function that calculates the value of the interpolated function and its derivative corresponding to "x" point
     * @param x is the point to interpolate
     * @return  The function and its derivative evaluated at the interpolation point x.
     */
    std::tuple<double, double> ValDeriv(double x);
    
    /** @brief Function that sets the left extension type (extrapolation type)
     *  @param left is the extension type, you can use "ELinear", "EConstant" or "ESlope"
     *  @param val is the value corresponding to the extrapolation type (as constant if you select "EConstant" or the slope if you select "ESlope")
     */
    void SetLeftExtension(Extension left, double val=0.);
    
    /** @brief Function that sets the right extension type (extrapolation type)
     *  @param right is the extension type, you can use "ELinear", "EConstant" or "ESlope".
     *  @param val is the value corresponding to the extrapolation type (as constant if you select "EConstant" or the slope if you select "ESlope")
     */
    void SetRightExtension(Extension right, double val=0.);
    
    /** @brief Function that sets the interpolation type
     *  @param type is the interpolation type, you can use "TLinear" or "THermite"
     */
    void SetInterpolationType(InterpType type);
    
    /** @brief Function that calculates the Lagrange Basis value "k" at x point.
     *  @param mat is the matrix that contains the interval to analyze
     *  @param k is the lagrange polynomial number to calculate Lk
     *  @param x is the interpolation point
     *  @return A vector with the Lagrange Basis at x point and its derivatives
     */
    TPZVec<double> LagrangeB(Matrix mat, int k,double x);
    
    /** @brief Function that calculates the Hermite Basis value "k" at x point
     *  @param mat is the matrix that contains the interval to analyze
     *  @param k is the lagrange polynomial number to calculate Lk
     *  @param x is the interpolation point
     *  @return A vector with the Hermite Basis at x point and its derivatives
     */
    TPZVec<double> HermiteB(Matrix mat, int k,double x);
    
    /** @brief Function that transforms the "TRSInterpolator" object as an function that returns the interpolated values
     *  @param Function that returns the interpolated value
     *  @return An interpolation function
     */
    std::function<double(double)> GetFunction();
    
    /** @brief Function that transforms the "TRSInterpolator" object as an function that returns the interpolated values and its derivatives
     *  @param Function that returns the interpolated value and its derivatives
     *  @return An interpolation function that returns the interpolated value and its derivatives
     */
    std::function<std::tuple<double, double>(double)> GetFunctionDeriv();
    
    
    
};
#endif /* TRSLinearInterpolator_h */
