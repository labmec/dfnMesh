//
// RibFrac.cpp
// RibFrac
//
//  Created by JORGE ORDOÑEZ on 22/5/18.
//  Copyright © 2018 JORGE ORDOÑEZ. All rights reserved.
//

#include "TRSRibFrac.h"

TRSRibFrac::TRSRibFrac(){
     feixos.Resize(3, 3);
     fdata.Resize(4, 3);
}

TRSRibFrac::TRSRibFrac(Matrix data){
     feixos.Resize(3, 3);
     Check_ConsistencyData(data);
}

void TRSRibFrac::SetPlane(Matrix plane){
     feixos.Resize(3, 3);
     Check_ConsistencyData(plane);
}

Matrix TRSRibFrac::GetPlane() const{
    return fdata;
}

void TRSRibFrac::SetTolerance(double tolerance){
    ftolerance = tolerance;
}

double TRSRibFrac::GetTolerance() const{
    return ftolerance;
}

void TRSRibFrac::Check_ConsistencyData(Matrix data) const {
  
// Checking vector consistency
    int cols = data.Cols();
    int rows = data.Rows();
    if(cols != 3){
        std::cout<<"Check the input data";
        DebugStop();
    }
    if(rows < 3 or rows > 4){
        std::cout<<"Check the input data";
        DebugStop();
    }
    
//Ax0 computation
   feixos(0,0)=data(1,0)-data(0,0);
    feixos(0,1)=data(1,1)-data(0,1);
    feixos(0,2)=data(1,2)-data(0,2);
    
//Ax1 without normalization
    feixos(1,0)=data(rows-1,0)-data(0,0);
    feixos(1,1)=data(rows-1,1)-data(0,1);
    feixos(1,2)=data(rows-1,2)-data(0,2);
    
//Ax1 normalization
    double normx = feixos(1,0) - (feixos(0,0))*((feixos(0,0)*feixos(1,0)) + (feixos(0,1)*feixos(1,1)) + (feixos(0,2)*feixos(1,2)));
    
    double normy = feixos(1,1)- feixos(0,1)*(feixos(0,0)*feixos(1,0) + feixos(0,1)*feixos(1,1) + feixos(0,2)*feixos(1,2));
    
    double normz = feixos(1,2)- feixos(0,2)*(feixos(0,0)*feixos(1,0) + feixos(0,1)*feixos(1,1) + feixos(0,2)*feixos(1,2));
    
    feixos(1,0)=normx;
    feixos(1,1)=normy;
    feixos(1,2)=normz;
    
//Ax2 computation
    feixos(2,0)=feixos(0,1)*feixos(1,2) - feixos(0,2)*feixos(1,1);
    feixos(2,1)=feixos(0,2)*feixos(1,0) - feixos(0,0)*feixos(1,2);
    feixos(2,2)=feixos(0,0)*feixos(1,1) - feixos(0,1)*feixos(1,0);

    
    if (rows ==4 && cols ==3){
        
//Coplanar verification
        double ver = feixos(2,0)*(data(2,0)-data(0,0))+feixos(2,1)*(data(2,1)-data(0,1))+feixos(2,2)*(data(2,2)-data(0,2));
        
        if(abs(ver) > ftolerance){
            std::cout<<"The points are not coplanar"<<"\n"<<std::endl;
            DebugStop();
        }
    }
    
//Mid points computation
    data.Resize(rows+rows, 3);
    for(int i=0; i<(rows); i++){
        data(rows+i,0)=(data(i,0)+data(i+1,0))/2;
        data(rows+i,1)=(data(i,1)+data(i+1,1))/2;
        data(rows+i,2)=(data(i,2)+data(i+1,2))/2;
    }
    data(rows+rows-1,0)=(data((rows)-1,0)+data(0,0))/2;
    data(rows+rows-1,1)=(data((rows)-1,1)+data(0,1))/2;
    data(rows+rows-1,2)=(data((rows)-1,2)+data(0,2))/2;
    
    data.Resize(rows+rows+1, 3);
    for(int i=0; i< rows; i++){
    
            data(rows+rows,0) += (1.0/rows)*data(i,0);
            data(rows+rows,1) += (1.0/rows)*data(i,1);
            data(rows+rows,2) += (1.0/rows)*data(i,2);
    }
    
    fdata = data;
}

//Checking if the given point is above the plane
bool TRSRibFrac::Check_point_above(TPZVec<double> point) const{
    int pm = fdata.Rows();
    double point_distance = (point[0] - fdata(pm-1,0))*feixos(2,0) + (point[1] - fdata(pm-1,1))*feixos(2,1)+(point[2] - fdata(pm-1,2))*feixos(2,2);
    if (point_distance>0){
        return true;
    }
    else{
        return false;
    }
}

//Checking if the given point is below the plane
bool TRSRibFrac::Check_point_below(TPZVec<double> point){
    bool var = Check_point_above(point);
    if (var == false){
        return true;
    }
    else{
        return false;
    }
}
bool TRSRibFrac::Check_rib(TPZVec<double> p1, TPZVec<double> p2){
    if(Check_point_below(p1)==true && Check_point_below(p2)==false){
        return true;
    }
    else if (Check_point_below(p1)==false && Check_point_below(p2)==true){
        return true;
    }
    else{
        return false;
    }
}
