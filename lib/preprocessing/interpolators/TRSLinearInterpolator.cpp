//
//  TRSLinearInterpolator.cpp
//
//  Created by JOSE VILLEGAS on 22/5/18.
//  Copyright © 2018 JOSE VILLEGAS. All rights reserved.
//

#include "TRSLinearInterpolator.h"
/**
 * @brief Default constructor
 */
TRSLinearInterpolator::TRSLinearInterpolator(){
    fdata.Resize(1, 1);
    fdata(0,0)=0;
}

/**
 * @brief Constructor
 * @param data is a matrix that contains the values (x,f(x)) for linear interpolation or (x, f(x), f'(x)) for Hermite Interpolation.
 */
TRSLinearInterpolator::TRSLinearInterpolator(Matrix data){
    fdata=data;
    
}

/** @brief Function that sets the matrix((x,f(x)) for linear interpolation or (x,f(x), f'(x)) for hermite interpolation) with the data to interpolate
 * @param  data is a matrix (nx2 for Linear interpolation or nx3 for Hermite interpolation) with the dimensional data to interpolate
 */
void TRSLinearInterpolator::SetData(Matrix data){
    if(data.Rows()>= 2){
        fdata = data;
    }
    else{
        DebugStop();
    }
}

/** @brief Function that returns a matrix (nx2 for Linear interpolation or nx3 for Hermite interpolation) with one dimensional data to interpolate
 * @return  The matrix (nx2) with one dimensional data to interpolate
 */
Matrix TRSLinearInterpolator::GetData(){
    return fdata;
}

/** @brief Function that reads the interpolation data from a ".txt" file
 * @param data is the name of the file with the interpolation data
 */

/**
 * @brief Interpolation function
 * @param x is a value to interpolate
 * @return Interpolated value if x0<x<xn or extrapolated value if x<x0 or x>xn
 */
double TRSLinearInterpolator::Val(double x){
    std::tuple<double, double> Intp = ValDeriv(x);
    
    return std::get<0>(Intp);
    
}
/**
 * @brief Interpolation function with derivative
 * @param x value to interpolate
 * @return Interpolated value and its derivative if x0<x<xn or extrapolated value and its derivative if x<x0 or x>xn
 */
std::tuple<double, double> TRSLinearInterpolator::ValDeriv(double x){
    
    int npoints = fdata.Rows();
    double returned=0;
    double deriv = 0;
    
    //Linear Interpolation
    
     switch (fInterType)
    {
            
    case TLinear:
    {
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
                return {returned,deriv};
                
                }
            }
        break;
        }
        break;
    }
            
        case THermite:
        {
            if(fdata.Cols()!=3){
                std::cout<<"Error: The interpolation of hermite needs the values ​​of the derivative of the function"<<"\n";
                DebugStop();
            }
            
            if (x >= fdata(0,0) && x<= fdata(npoints-1,0) ){
                for(int i=0;  i<npoints-1; i++){
                    if(x >= fdata(i,0) && x<= fdata(i+1,0)){
                        Matrix Vals(2,3);
                        Vals(0,0) = fdata(i,0);
                        Vals(1,0) = fdata(i+1,0);
                        Vals(0,1) = fdata(i,1);
                        Vals(1,1) = fdata(i+1,1);
                        Vals(0,2) = fdata(i,2);
                        Vals(1,2) = fdata(i+1,2);
                        for(int j = 0; j<2;j++){
                            TPZVec<double> res = HermiteB(Vals, j, x);
                            returned = returned+(Vals(j,1))*res[0] +(Vals(j,2))*res[1];
                            deriv = deriv + (Vals(j,1))*res[2] +(Vals(j,2))*res[3];
                        }
                     return {returned,deriv};
                    }
                }
                break;
            }
            break;
        }
            break;
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
            case EConstant:
                returned = fvalLef;
                deriv = 0.0;
                break;
            case ESlope:
                returned = y1 - fvalLef*(x1 - x);
                deriv = fvalLef;
                break;
            case ELinear:
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
            case ELinear:
                returned = y2 + ((x - x2)*(y2-y1)/(x2-x1));
                deriv = (y2-y1)/(x2-x1);
                break;
            case EConstant:
                returned = fvalRight;
                deriv =0.0;
                break;
            case ESlope:
                returned = y2 + fvalRight*(x-x2);
                deriv = fvalRight;
            default:
                break;
        }
    }
    
    return {returned,deriv};
}

/**
 * @brief Function to define the extrapolation type from the left
 * @param left, you can select: "constant", "linear" (use the same slope of the first two points) or "Slope" to extrapolate with a slope from the first point.
  @param The val parameter sets the value of the slope (in the extrapolation region x < x0). This value is only required in the case of need a linear extrapolation with a given slope, in the other cases this parameter is automatically calculated
 */

void TRSLinearInterpolator::SetLeftExtension(Extension left, double val){
    fextLeft = left;
    fvalLef =val;
}

/**
 * @brief Function to define the extrapolation type from the right
 * @param right, you can select: "constant", "linear" (use the same slope of the last two points) or "Slope" to extrapolate with a slope from the last point.
 @param The parameter val sets the value of the slope (in the extrapolation region x > xn). This value is only required in the case of need a linear extrapolation with a given slope, in the other cases this parameter is automatically calculated
 */
void TRSLinearInterpolator::SetRightExtension(Extension right, double val){
    fextRight = right;
    fvalRight = val;
}

/**
 * @brief This function defines the method to be used to interpolate
 * @param type, you can select: "TLinear" (use linear functions defined in each pair of points to interpolate.) or "THermite" (use hermite polynomials defined by each pair of points)
 */
void TRSLinearInterpolator::SetInterpolationType(InterpType type){
    fInterType = type;
}

/**
 * @brief This Function calculates the Lagrange basis functions
 * @param mat is a matrix that contains the table of interpolation values (x, f (x))
 * @param k is the point value at which the lagrange function is needed
 * @param x is the point at which the lagrange polynomial associated with point k is evaluated
 * @return A vector with the lagrange polynomial value associated with k point and the derivative valueof the polynomial evaluated at point x
 */
TPZVec<double> TRSLinearInterpolator::LagrangeB(Matrix mat, int k,double x){
    
    int rows = mat.Rows();
  
    TPZVec<double> lagran(2);
    double deriv=1.0;
    double val = 1;
    for(int i=0; i<rows; i++){
        if(i!=k){
            val = val*((x-mat(i,0))/(mat(k,0)-mat(i,0)));
            deriv = deriv*(1/(mat(k,0)-mat(i,0)));
        }
    }
    lagran[0]=val;
    lagran[1]=deriv;
    
    return lagran;
}

/**
 * @brief This Function calculates the Hermite basis functions
 * @param mat is a matrix that contains the interpolation table values (x, f (x), f'(x)).
 * @param k is a point value at which the Hermite function is needed
 * @param x is a point at which the Hermite polynomial associated with k point is evaluated
 * @return A vector with the Hermite polynomial value  associated with k point and the derivative value of the polynomial evaluated at x point
 
 */
TPZVec<double> TRSLinearInterpolator::HermiteB(Matrix mat, int k,double x){
    TPZVec<double> BaseLagrange = LagrangeB(mat, k, x);
    TPZVec<double> BaseHermite(4);
    double hermi = (1 - 2*(x-mat(k,0))*(BaseLagrange[1]))*(BaseLagrange[0]*BaseLagrange[0]);
    double hermi2 = (x - mat(k,0))*(BaseLagrange[0]*BaseLagrange[0]);
    
    double hermiD = 2*BaseLagrange[0]*(((1.0 + 2.0*BaseLagrange[1]*(mat(k,0) - x)) * BaseLagrange[1])- (BaseLagrange[1] * BaseLagrange[0]));
    
    double hermi2D = BaseLagrange[0]*(BaseLagrange[0] + 2*(x-mat(k,0))*BaseLagrange[1]);
     BaseHermite[0]=hermi;
     BaseHermite[1]=hermi2;
     BaseHermite[2]=hermiD;
     BaseHermite[3]=hermi2D;
    return BaseHermite;
}

/**
 * @brief Function that generates an interpolation function
 @return  Interpolation function
 */

std::function<double(double)> TRSLinearInterpolator::GetFunction(){
    return [this](double x){
        return this->Val(x);
    };
}

/**
 * @brief Function that generates an interpolation function (and its derivatives)
 @return  Interpolation function (and its derivatives)
 */
std::function<std::tuple<double, double>(double)> TRSLinearInterpolator::GetFunctionDeriv() {
    return [this](double x){
        return this->ValDeriv(x);
    };
}


void TRSLinearInterpolator::ReadData(std::string name){
   
    std::ifstream file;
    file.open(name);
    
    int i=1;
    Matrix data;
    std::string line;
    int n_cols =0;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::istringstream issText(line);
        char l = line[0];
        int nchar = line.size();
        if(l != '/'){
            
            //calcula el numero de columnas
            if(n_cols==0){
                std::string w;
                for(int j=0; j< nchar; j++){
                    w="";
                    if (issText >> w);
                    if(w==""){
                        break;
                    }
                    n_cols ++;
                }
            }
        
            if (n_cols==2) {
                double a, b;
                data.Resize(i, 2);
                if(iss >> a >> b) ;
                data(i-1,0)=a;
                data(i-1,1)=b;
                i=i+1;
            }
            if (n_cols==3) {
                double a, b, c;
                data.Resize(i, 2);
                if(iss >> a >> b >> c) ;
                data(i-1,0)=a;
                data(i-1,1)=b;
                data(i-1,2)=c;
                i=i+1;
                }
            }
         }
        
    
        if(data.Rows()>0){
            std::cout<<"*************************"<<std::endl;
            std::cout<<"Reading file... ok!"<<std::endl;
            std::cout<<"*************************"<<std::endl;
            SetData(data);
            data.Print(std::cout);
        }
        
        
        
        
        
//        char l = line[0];
//        if (l == '{' or l=='}') {
//            std::cout<<"Error: Incorrect Format"<<"\n";
//            DebugStop();
//        }
//        if(l != '/'){
//            // std::string l =line[0];
//            double a, b, c;
//            int val = line.size();
//            int count = 0;
//            std::string w;
//            for(int i =0; i< val; i++){
//                w="";
//                if (issText >> w >> w >> w);
//                if(iss >> a);
//                double valor = a;
//
//            }
//
//
//            if ((iss >> a >> b >> c)) {
//
//                data.Resize(i, 3);
//                data(i-1,0)=a;
//                data(i-1,1)=b;
//                data(i-1,2)=c;
//                i=i+1;
//
//                //break;
//            } // error
//
//            if(abs(c)<1.1E-10){
//                if (data.Cols()!=3){
//                    if(a && b){
//                        iss >> a >> b ;
//                        data.Resize(i, 2);
//                        data(i-1,0)=a;
//                        data(i-1,1)=b;
//                        i=i+1;
//                    }
//
//                }
//                //                std::cout<<"\n"<<"Archivo con problema de lectura: la linea: "<<i<<" ha sido omitida"<<"\n";
//            }
//
//            if(!a or !b){
//                //                std::cout<<"\n"<<"Archivo con problema de lectura: la linea: "<<i<<" ha sido omitida"<<"\n";
//            }
//
//        }
//
//
//        // process pair (a,b)
//    }
//
//    file.close();
//
//    if(data.Rows()>0){
//        std::cout<<"*************************"<<std::endl;
//        std::cout<<"Reading file... ok!"<<std::endl;
//        std::cout<<"*************************"<<std::endl;
//        SetData(data);
//        data.Print(std::cout);
   
    
    // data.Print(std::cout);
    
    
    
}
