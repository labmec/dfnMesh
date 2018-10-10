//
//  InterpolPhil.cpp
//  Interpolator
//
//  Created by JOSE VILLEGAS on 22/5/18.
//  Copyright © 2018 JOSE VILLEGAS. All rights reserved.
//

#include "TRSLinearInterpolator.h"

TRSLinearInterpolator::TRSLinearInterpolator(){
    fdata.Resize(1, 1);
    fdata(0,0)=0;
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

Matrix TRSLinearInterpolator::GetData(){
    return fdata;
}

void TRSLinearInterpolator::ReadData(std::string name){
    std::ifstream file;
    file.open(name);
   
    int i=1;
    Matrix data;
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
       
        char l = line[0];
        if (l == '{' or l=='}') {
            std::cout<<"Error: Incorrect Format"<<"\n";
            DebugStop();
        }
        
        if(l != '/'){
       // std::string l =line[0];
        double a, b, c;
          
        if ((iss >> a >> b >> c)) {
            
            data.Resize(i, 3);
            data(i-1,0)=a;
            data(i-1,1)=b;
            data(i-1,2)=c;
            i=i+1;
          
            //break;
        } // error
            
        if(!c){
            if (data.Cols()!=3){
                if(a && b){
                    iss >> a >> b ;
                    data.Resize(i, 2);
                    data(i-1,0)=a;
                    data(i-1,1)=b;
                    i=i+1;
                    }
          
                }
                std::cout<<"\n"<<"Archivo con problema de lectura: la linea: "<<i<<" ha sido omitida"<<"\n";
          }
            
            if(!a or !b){
                std::cout<<"\n"<<"Archivo con problema de lectura: la linea: "<<i<<" ha sido omitida"<<"\n";
            }
            
        }
       
        
        // process pair (a,b)
    }
    file.close();
   
    if(data.Rows()>0){
        std::cout<<"Reading file... ok!"<<"\n";
        SetData(data);
    }
    
   // data.Print(std::cout);
}
/**
 * @brief Interpolation function
 * @param x value to interpolate
 * @return interpolated value if x0<x<xn or extrapolated value if (x<x0 or x>xn)
 */
double TRSLinearInterpolator::Val(double x){
    std::tuple<double, double> Intp = ValDeriv(x);
    
    return std::get<0>(Intp);
    
}
/**
 * @brief Interpolation function with derivative
 * @param x value to interpolate
 * @return interpolated value and its derivative if x0<x<xn or extrapolated value and its derivative if (x<x0 or x>xn)
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
                break;
                }
            break;
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
                        break;
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
/**
 * @brief function to define the type of extrapolation from the left.
 * @param left You can select: "constant", "linear" (use the same slope of the first two points), or "Slope" to extrapolate with a slope from the first point.
 
  @param val This parameter sets the value of the slope (in the extrapolation region x < x0). This value is only required in the case of needing a linear extrapolation with a given slope, in the other cases this parameter is calculated automatically.
 */

void TRSLinearInterpolator::SetLeftExtension(Extension left, double val){
    fextLeft = left;
    fvalLef =val;
}


/**
 * @brief Function to define the type of extrapolation from the right.
 * @param right  You can select: "constant", "linear" (use the same slope of the last two points), or "Slope" to extrapolate with a slope from the last point.
 @param val This parameter sets the value of the slope (in the extrapolation region x> xn). This value is only required in the case of needing a linear extrapolation with a given slope, in the other cases this parameter is calculated automatically.
 */
void TRSLinearInterpolator::SetRightExtension(Extension right, double val){
    fextRight = right;
    fvalRight = val;
}

/**
 * @brief This function defines the method to be used to interpolate.
 * @param type You can select: "Linear" (Use linear functions defined in each pair of points to interpolate.) or "Hermite" (Use hermite polynomials defined by each pair of points.).
 */
void TRSLinearInterpolator::SetInterpolationType(InterpType type){
    fInterType = type;
}
/**
 * @brief This Function calculates the basic functions of lagrange
 * @param mat Matrix that contains the table of interpolation values. (x, f (x)).
 @param k Value of the point at which the lagrange function is needed.
 @param x Point at which the lagrange polynomial associated with point k is evaluated.
 @return Returns a vector with the value of the lagrange polynomial associated with point k and the value of the derivative of the polynomial evaluated at point x.
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
 * @brief This Function calculates the basic functions of Hermite
 * @param mat Matrix that contains the table of interpolation values. (x, f (x), f'(x)).
 @param k value of the point at which the Hermite function is needed.
 @param x Point at which the Hermite polynomial associated with point k is evaluated.
  @return Returns a vector with the value of the Hermite polynomial associated with point k and the value of the derivative of the polynomial evaluated at point x.
 
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



//void TRSLinearInterpolator::Print(const char *name, std::ostream& out,const MatrixOutputFormat form) const {
//
//}

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
 * @brief Function that generates an interpolation function (and derivatives)
 @return  Interpolation function (and derivatives)
 */

std::function<std::tuple<double, double>(double)> TRSLinearInterpolator::GetFunctionDeriv() {
    return [this](double x){
        return this->ValDeriv(x);
    };
}
