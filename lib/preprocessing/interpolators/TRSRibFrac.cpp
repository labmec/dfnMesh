//
// RibFrac.cpp
// RibFrac
//
//  Created by JORGE ORDOÑEZ on 22/5/18.
//  Copyright © 2018 JORGE ORDOÑEZ. All rights reserved.
//

#include "TRSRibFrac.h"

TRSRibFrac::TRSRibFrac(){
     faxis.Resize(3, 3);
     fdata.Resize(4, 3);
}

TRSRibFrac::TRSRibFrac(Matrix data){
     faxis.Resize(3, 3);
     Check_ConsistencyData(data);
}

void TRSRibFrac::SetPlane(Matrix plane){
     faxis.Resize(3, 3);
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

/**
 * @brief Consistency of plane points
 * @param Matrix nx3, n is the number of points with the axis x, y , z
 * @return True if the points are coplanar
 * @return False if the points are not coplanar
 */

void TRSRibFrac::Check_ConsistencyData(Matrix data) {
  
// Checking vector consistency
    int cols = data.Cols();
    int rows = data.Rows();
    if(cols != 3){
        std::cout<<"Check the input data";
        DebugStop();
    }
    if(rows < 3 or rows > 4){
        std::cout<<"Check the input data (number of plane points";
        DebugStop();
    }
    
//Ax0 computation
    faxis(0,0)=data(1,0)-data(0,0);
    faxis(0,1)=data(1,1)-data(0,1);
    faxis(0,2)=data(1,2)-data(0,2);
    
//Ax1 without normalization
    faxis(1,0)=data(rows-1,0)-data(0,0);
    faxis(1,1)=data(rows-1,1)-data(0,1);
    faxis(1,2)=data(rows-1,2)-data(0,2);
    
//Ax1 normalization
    
    for (int i=0; i<3; i++) {
     double norm=faxis(1,i) - (faxis(0,i))*((faxis(0,0)*faxis(1,0)) + (faxis(0,1)*faxis(1,1)) + (faxis(0,2)*faxis(1,2)));
          faxis(1,i)=norm;
    };
    
//Ax2 computation
    faxis(2,0)=faxis(0,1)*faxis(1,2) - faxis(0,2)*faxis(1,1);
    faxis(2,1)=faxis(0,2)*faxis(1,0) - faxis(0,0)*faxis(1,2);
    faxis(2,2)=faxis(0,0)*faxis(1,1) - faxis(0,1)*faxis(1,0);

    if (rows ==4 && cols ==3){
        
//Coplanar verification
        double ver = faxis(2,0)*(data(2,0)-data(0,0))+faxis(2,1)*(data(2,1)-data(0,1))+faxis(2,2)*(data(2,2)-data(0,2));
        
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
    
//Center point computation
    
    data.Resize(rows+rows+1, 3);
    for(int i=0; i< rows; i++){
    
            data(rows+rows,0) += (1.0/rows)*data(i,0);
            data(rows+rows,1) += (1.0/rows)*data(i,1);
            data(rows+rows,2) += (1.0/rows)*data(i,2);
    }
    
    fdata = data;
}

/**
 * @brief Check if a point is above or below a plane
 * @param Point
 * @return True if the point is above the plane
 * @return False if the point is below the plane
 */

//Checking if the given point is above the plane
bool TRSRibFrac::Check_point_above(TPZVec<double> point) const{
    int pm = fdata.Rows();
    double point_distance = (point[0] - fdata(pm-1,0))*faxis(2,0) + (point[1] - fdata(pm-1,1))*faxis(2,1)+(point[2] - fdata(pm-1,2))*faxis(2,2);
    if (point_distance>0){
        return true;
    }
    else{
        return false;
    };
};

/**
 * @brief Check two points if they are on both sides of the plan
 * @param Point 1 and Point 2
 * @return True if both are on two sides of the plane
 * @return False otherwise
 */

bool TRSRibFrac::Check_rib(TPZVec<double> p1, TPZVec<double> p2) const {
    if(Check_point_above(p1)==true && Check_point_above(p2)==false){
        return true;
    }
    else if (Check_point_above(p1)==false && Check_point_above(p2)==true){
        return true;
    }
    else{
        return false;
    }
}

/**
 * @brief Check if the neighbour has a lower dimension
 * @param Geo element side
 * @return True if has a lower dimension
 * @return False if has a upper dimension
 */

bool TRSRibFrac::HasLowerDimensionNeighbour(TPZGeoElSide &gelside){
    int dimension = gelside.Dimension();
   
    if (gelside.Element()->Dimension() == dimension){
        return true;
    }

    TPZGeoElSide neighbour = gelside.Neighbour();
    
    while (neighbour != gelside){
        if (neighbour.Element()->Dimension()==dimension){
            return true;
            neighbour = neighbour.Neighbour();
         }
    return false;
    }
}

/**
 * @brief Creates the skeleton mesh
 * @param Dimension
 * @param Material ID number
 * @return Ribs with a selected material ID
 */

void TRSRibFrac::CreateSkeleton(int dimension, int matid){
    if(fmesh){
        int nribs = fRibVec.size();
        int nel = fmesh->NElements();
        for(int iel=0; iel<nel; iel++){
            TPZGeoEl *gel = fmesh->Element(iel);
            int nsides = gel->NSides();
            for(int iside=0; iside<nsides; iside++){
                
                TPZGeoElSide gelside = gel->Neighbour(iside);
//                std::cout<<gelside.Dimension()<<std::endl;
//                gel->Print();
//                gelside.Print(std::cout);
                if (gelside.Dimension()==dimension){
                    bool haskel = HasLowerDimensionNeighbour(gelside);
                    if(haskel==false){
                        TPZGeoElBC(gelside,matid);
                        //crea los ribs;
                        if(gelside.Dimension()==1){
                            int side = gelside.Side();
                            int64_t p1 =  gel->SideNodeIndex(side, 0);
                            int64_t p2 =  gel->SideNodeIndex(side, 1);
                            TPZVec<REAL> pp1(3);
                            TPZVec<REAL> pp2(3);
                            fmesh->NodeVec()[p1].GetCoordinates(pp1);
                            fmesh->NodeVec()[p2].GetCoordinates(pp2);
                            TPZVec<TPZVec<double>> vecpoints(2);
                            vecpoints[0]=pp1;
                            vecpoints[1]=pp2;
                            fRibVec.resize(nribs+1);
                            TRSRibs rib(vecpoints);
                            nribs=fRibVec.size();
                            fRibVec[nribs-1]=TRSRibs(rib);
                        }
                    }
                }
            }
        }
        
    }
    else{
        std::cout<<"XXXX";
    }
}

/**
 * @brief Calculates the intersection point plane-rib
 * @param Point above the plane (node rib)
 * @param Point below the plane (node rib)
 * @return Intersecting point
 */

TPZVec<double> TRSRibFrac::CalculateIntersection(TPZVec<double> pa, TPZVec<double> pb){

    
    double term1 = (pa[0]-pb[0])*faxis(2,0) + (pa[1]-pb[1])*faxis(2,1) + (pa[2]-pb[2])*faxis(2,2);
    int p8_index = fdata.Rows() - 1;
    double term2 = ((fdata(p8_index,0)-pb[0])*(faxis(2,0))) + ((fdata(p8_index,1)-pb[1])*(faxis(2,1))) + ((fdata(p8_index,2)-pb[2])*(faxis(2,2)));
    double alpha = term2/term1;
    TPZVec<double> Pint(3,0.0);
    for (int p=0; p<3; p++){
    Pint[p] = pb[p] + alpha*(pa[p]-pb[p]);
    return Pint;
    }
}
TPZVec<TRSRibs> TRSRibFrac::GetRibsVec(){
    if(fRibVec.size()>0){
        return fRibVec;
    }
}

/**
 * @brief Divide a rib
 * @param Index of the element to divide (for the moment just ribs)
 * @param Point of intersection
 * @return Intersecting point
 */

//void TRSRibFrac::DivideRib(int element_index, TPZVec<double> intersection){
//    TPZGeoEl *gel = fmesh->Element(element_index);
//    if(gel->Dimension()!=1){
//        std::cout<<"Just ribs now!";
//        DebugStop();
//    }
//    TPZVec<TPZGeoEl *> gelsDivide(2);
//    gel->Divide(gelsDivide);
//
//}
