//
// RibFrac.cpp
// RibFrac
//
//  Created by JORGE ORDOÑEZ on 22/5/18.
//  Copyright © 2018 JORGE ORDOÑEZ. All rights reserved.
//

#include "TRSRibFrac.h"

// Constructor
TRSRibFrac::TRSRibFrac(){
}

TRSRibFrac::TRSRibFrac(const Matrix &data, TPZGeoMesh *gmesh){
    if(!Check_ConsistencyData(data)){
        std::cout<<"The input data is not correct";
        DebugStop();
    }
    fCornerCoordinates=data;
    fGMesh=gmesh;
}

TRSRibFrac::TRSRibFrac(const TRSRibFrac &copy){
    fTolerance = copy.GetTolerance();
    fRibs = copy.fRibs;
    
    fCornerCoordinates=copy.GetPlane();
    fAxis=copy.fAxis;
    fCenterCo = copy.fCenterCo;
    fGMesh=copy.fGMesh;
}

TRSRibFrac &TRSRibFrac::operator=(const TRSRibFrac &copy){
    fTolerance = copy.GetTolerance();
    fRibs = copy.fRibs;
    fCornerCoordinates=copy.GetPlane();
    fAxis=copy.fAxis;
    fCenterCo = copy.fCenterCo;
    fGMesh=copy.fGMesh;
}

/**
 * @brief Set the fracture plane
 * @param Fracture plane coordinates (Matrix 3x4)
 * @return The plane set
 */

void TRSRibFrac::SetPlane(Matrix plane){
    if(!Check_ConsistencyData(plane)){
        std::cout<<"The input data is not correct";
        DebugStop();
    }
    fCornerCoordinates=plane;
}

Matrix TRSRibFrac::GetPlane() const{
    return fCornerCoordinates;
}

void TRSRibFrac::SetTolerance(REAL tolerance){
    fTolerance = tolerance;
}

REAL TRSRibFrac::GetTolerance() const{
    return fTolerance;
}

/**
 * @brief Checks the consistency of the input data
 * @param Fracture plane coordinates (Matrix 3x4)
 * @return True if the four points are coplanar and the data is consistent
 * @return False if the points are not coplanar or the data is not consistent
 */

bool TRSRibFrac::Check_ConsistencyData(Matrix CornerCoordinates) {

// Checking vector consistency
    int cols = CornerCoordinates.Cols();
    int rows = CornerCoordinates.Rows();
    Matrix ax(3,3);
    
    if(rows != 3){                                  //Should be 3 (x y z)
        std::cout<<"Check the input data";
        DebugStop();
    }
    if(cols < 3 or cols > 4){             //To form a plane it is needed at least 3 points
        std::cout<<"Check the input data (number of plane points, four is enough)";
        DebugStop();
    }
   Matrix MidPoint;                       //Mid points computation
    if (!(cols ==4 && rows ==3)){
        std::cout<<"Check the input data (number of plane points, four is enough)";
        DebugStop();
    }
    MidPoint.Resize(3,cols);              //3 rows (x y z), mid points

// Mid points computation for the three first points
    for(int i=0; i<(cols-1); i++){
        MidPoint(0,i)=(CornerCoordinates(0,i)+CornerCoordinates(0,i+1))/2;
        MidPoint(1,i)=(CornerCoordinates(1,i)+CornerCoordinates(1,i+1))/2;
        MidPoint(2,i)=(CornerCoordinates(2,i)+CornerCoordinates(2,i+1))/2;
    }
// Mid points computation for the last point
    MidPoint(0,cols-1)=(CornerCoordinates(0,(cols)-1)+CornerCoordinates(0,0))/2;
    MidPoint(1,cols-1)=(CornerCoordinates(1,(cols)-1)+CornerCoordinates(1,0))/2;
    MidPoint(2,cols-1)=(CornerCoordinates(2,(cols)-1)+CornerCoordinates(2,0))/2;

//Ax0 computation
    ax(0,0)=MidPoint(0,cols-1)-MidPoint(0,1);
    ax(1,0)=MidPoint(1,cols-1)-MidPoint(1,1);
    ax(2,0)=MidPoint(2,cols-1)-MidPoint(2,1);

//Ax1 without normalization
    ax(0,1)=MidPoint(0,cols-2)-MidPoint(0,0);
    ax(1,1)=MidPoint(1,cols-2)-MidPoint(1,0);
    ax(2,1)=MidPoint(2,cols-2)-MidPoint(2,0);

//Ax1 normalization
    for (int i=0; i<3; i++) {                 /// i< axis number (x y z)
        double norm=ax(i,1) - (ax(i,0))*((ax(0,0)*ax(0,1)) + (ax(1,0)*ax(1,1)) + (ax(2,0)*ax(2,1)));
        ax(i,1)=norm;
    }

//Ax2 computation
    ax(0,2)=ax(1,0)*ax(2,1) - ax(2,0)*ax(1,1);
    ax(1,2)=ax(2,0)*ax(0,1) - ax(0,0)*ax(2,1);
    ax(2,2)=ax(0,0)*ax(1,1) - ax(1,0)*ax(0,1);

//Coplanar verification
    double ver = ax(0,2)*(CornerCoordinates(0,2)-CornerCoordinates(0,0))+ax(1,2)*(CornerCoordinates(1,2)-CornerCoordinates(1,0))+ax(2,2)*(CornerCoordinates(2,2)-CornerCoordinates(2,0));

    if(abs(ver) > fTolerance){
        std::cout<<"The input points are not coplanar"<<"\n"<<std::endl;
        DebugStop();
    }
        
// After checking the consistency the date is set
    
        fCornerCoordinates=CornerCoordinates;
//        fCornerCoordinates.Print(std::cout);
        fAxis=ax;
        fCenterCo.Resize(3);
// Center point computation
    
        for(int i=0; i< 4; i++){                 /// i< axis number (x y z)
            fCenterCo[0] += (1.0/cols)*fCornerCoordinates(0,i); ///Center point X axis
            fCenterCo[1] += (1.0/cols)*fCornerCoordinates(1,i); ///Center point Y axis
            fCenterCo[2] += (1.0/cols)*fCornerCoordinates(2,i); ///Center point Z axis
        }
    
    return true;
}

/**
 * @brief Checks if a point is above or below the fracture plane
 * @param Point vector with the euclidean coordinates
 * @return True if the point is above the fracture plane
 * @return False if the point is below the fracture plane
 */

bool TRSRibFrac::Check_point_above(const TPZVec<REAL> &point) const{
    
        Matrix fAxis_temp = fAxis;     ///Creating a copy of fAxis matrix
    /// Point distance to the fracture plane computation
        double point_distance = (point[0] - fCenterCo[0])*(fAxis_temp(0,2)) + (point[1] - fCenterCo[1])*(fAxis_temp(1,2)) + (point[2] - fCenterCo[2])*(fAxis_temp(2,2));
        if (point_distance>0){
            return true;
        }
        else{
            return false;
        }
}

/**
 * @brief Checks if a rib is cut by a fracture plane
 * @param Point vector with the euclidean coordinates
 * @param Point vector with the euclidean coordinates
 * @return True if the rib is cut by the fracture plane
 * @return False if the rib is not cut by the fracture plane
 */

bool TRSRibFrac::Check_rib(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2) const{
    
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
  */

void TRSRibFrac::CreateSkeletonElements(int dimension, int matid){
    
    int nel = fGMesh->NElements();
    for(int iel=0; iel<nel; iel++){
        TPZGeoEl *gel = fGMesh->Element(iel);
        int nsides = gel->NSides();
            for(int iside=0; iside<nsides; iside++){
                TPZGeoElSide gelside = gel->Neighbour(iside);
                    if (gelside.Dimension()==dimension){
                        bool haskel = HasLowerDimensionNeighbour(gelside);
                            if(haskel==false){
                                int nel_mesh = fGMesh->NElements();
                                TPZGeoElBC(gelside,matid);
                           
                                    //Setting Ribs
                            
                                    TRSRibs rib(nel_mesh,false);
                                    SetRib(rib);
                }
            }
        }
    }
}

/**
 * @brief Sets the rib idexes
 * @param Ribs
 */

void TRSRibFrac::SetRib(TRSRibs rib){
    int index= rib.ElementIndex();
    fRibs[index]=rib;
}

/**
 * @brief Gives the ribs map
 * @return A map that contains the ribs information
 */

std::map<int64_t ,TRSRibs> TRSRibFrac::GetRibs(){
    return fRibs;
}








//TRSRibFrac::TRSRibFrac(){
//     fAxis.Resize(3, 3);
//     fCornerCoordinates.Resize(3, 4);
//}
//
//TRSRibFrac::TRSRibFrac(Matrix data){
//     fAxis.Resize(3, 3);
//     Check_ConsistencyData();
//}
//
//void TRSRibFrac::SetPlane(Matrix plane){
//     fAxis.Resize(3, 3);
//     Check_ConsistencyData();
//}
//
//Matrix TRSRibFrac::GetPlane() const{
//    return fCornerCoordinates;
//}
//
//void TRSRibFrac::SetTolerance(double tolerance){
//    fTolerance = tolerance;
//}
//
//double TRSRibFrac::GetTolerance() const{
//    return fTolerance;
//}
//
///**
// * @brief Consistency of plane points
// * @param Matrix nx3, n is the number of points with axis x, y , z
// * @return True if the points are coplanar
// * @return False if the points are not coplanar
// */
//
//void TRSRibFrac::Check_ConsistencyData() {
//
//// Checking vector consistency
//    int cols = fCornerCoordinates.Cols();
//    int rows = fCornerCoordinates.Rows();
//    if(rows != 3){                                  //Should be 3 (x y z)
//        std::cout<<"Check the input data";
//        DebugStop();
//    }
//    if(cols < 3 or cols > 4){             //To form a plane it is needed at least 3 points
//        std::cout<<"Check the input data (number of plane points, four is enough)";
//        DebugStop();
//    }
//   Matrix MidPoint;                       //Mid points computation
//    if (cols ==4 && rows ==3){
//    MidPoint.Resize(3,cols);              //3 rows (x y z), mid points and center point columns
//
//    for(int i=0; i<(cols); i++){
//        MidPoint(0,i)=(fCornerCoordinates(0,i)+fCornerCoordinates(0,i+1))/2;
//        MidPoint(1,i)=(fCornerCoordinates(1,i)+fCornerCoordinates(1,i+1))/2;
//        MidPoint(2,i)=(fCornerCoordinates(2,i)+fCornerCoordinates(2,i+1))/2;
//    }
//    MidPoint(0,cols-1)=(fCornerCoordinates(0,(cols)-1)+fCornerCoordinates(0,0))/2;
//    MidPoint(1,cols-1)=(fCornerCoordinates(1,(cols)-1)+fCornerCoordinates(1,0))/2;
//    MidPoint(2,cols-1)=(fCornerCoordinates(2,(cols)-1)+fCornerCoordinates(2,0))/2;
//
//    //Ax0 computation
//    fAxis(0,0)=MidPoint(0,cols-1)-MidPoint(0,1);
//    fAxis(1,0)=MidPoint(1,cols-1)-MidPoint(1,1);
//    fAxis(2,0)=MidPoint(2,cols-1)-MidPoint(2,1);
//
//    //Ax1 without normalization
//    fAxis(0,1)=MidPoint(0,cols-2)-MidPoint(0,0);
//    fAxis(1,1)=MidPoint(1,cols-2)-MidPoint(1,0);
//    fAxis(2,1)=MidPoint(2,cols-2)-MidPoint(2,0);
//
//    //Ax1 normalization
//    for (int i=0; i<3; i++) {
//        double norm=fAxis(i,1) - (fAxis(i,0))*((fAxis(0,0)*fAxis(0,1)) + (fAxis(1,0)*fAxis(1,1)) + (fAxis(2,0)*fAxis(2,1)));
//        fAxis(i,1)=norm;
//    }
//
//    //Ax2 computation
//    fAxis(0,2)=fAxis(1,0)*fAxis(2,1) - fAxis(2,0)*fAxis(1,1);
//    fAxis(1,2)=fAxis(2,0)*fAxis(0,1) - fAxis(0,0)*fAxis(2,1);
//    fAxis(2,2)=fAxis(0,0)*fAxis(1,1) - fAxis(1,0)*fAxis(0,1);
//
//    //Coplanar verification
//    double ver = fAxis(0,2)*(fCornerCoordinates(0,2)-fCornerCoordinates(0,0))+fAxis(1,2)*(fCornerCoordinates(1,2)-fCornerCoordinates(1,0))+fAxis(2,2)*(fCornerCoordinates(2,2)-fCornerCoordinates(2,0));
//
//    if(abs(ver) > fTolerance){
//        std::cout<<"The input points are not coplanar"<<"\n"<<std::endl;
//        DebugStop();
//        }
//
//    // Center point computation
//        int cols = fCornerCoordinates.Cols();
//        for(int i=0; i< 3; i++){                 /// i< axis number (x y z)
//            fCenterCo[0] += (1.0/cols)*fCornerCoordinates(0,i); ///Center point X axis
//            fCenterCo[1] += (1.0/cols)*fCornerCoordinates(1,i); ///Center point Y axis
//            fCenterCo[2] += (1.0/cols)*fCornerCoordinates(2,i); ///Center point Z axis
//        }
//    }
//}
//
///**
// * @brief Check if a point is above or below a plane
// * @param Point (x y z coordinates)
// * @return True if the point is above the plane
// * @return False if the point is below the plane
// */
//
////Checking if the given point is above the plane
//bool TRSRibFrac::Check_point_above(const TPZVec<REAL> &point) const{
//    Matrix fAxis_temp = fAxis;           ///Creating a copy of fAxis matrix
//
//    double point_distance = (point[0] - fCenterCo[0])*(fAxis_temp(0,2)) + (point[1] - fCenterCo[1])*(fAxis_temp(1,2)) + (point[2] - fCenterCo[2])*(fAxis_temp(2,2));
//    if (point_distance>0){
//        return true;
//    }
//    else{
//        return false;
//    }
//}
//
///**
// * @brief Checks two points if they are on both sides of the plane
// * @param Point 1
// * @param Point 2
// * @return True if both are on two sides of the plane
// * @return False otherwise
// */
//
//bool TRSRibFrac::Check_rib(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2) const {
//    if(Check_point_above(p1)==true && Check_point_above(p2)==false){
//        return true;
//    }
//    else if (Check_point_above(p1)==false && Check_point_above(p2)==true){
//        return true;
//    }
//    else{
//        return false;
//    }
//}
//
///**
// * @brief Check if the neighbour has a lower dimension
// * @param Geo element side
// * @return True if has a lower dimension
// * @return False if has a upper dimension
// */
//
//bool TRSRibFrac::HasLowerDimensionNeighbour(TPZGeoElSide &gelside){
//    int dimension = gelside.Dimension();
//
//    if (gelside.Element()->Dimension() == dimension){
//        return true;
//    }
//
//    TPZGeoElSide neighbour = gelside.Neighbour();
//
//    while (neighbour != gelside){
//        if (neighbour.Element()->Dimension()==dimension){
//            return true;
//            neighbour = neighbour.Neighbour();
//         }
//    return false;
//    }
//
//}
//
///**
// * @brief Creates the skeleton mesh
// * @param Dimension
// * @param Material ID number
// * @return Ribs with a selected material ID
// */
//
// void TRSRibFrac::CreateSkeletonElements(int dimension, int matid){
//    if(fGMesh){
//        int nribs = fRibs.size();           /// no se si es lo correcto
//        int nel = fGMesh->NElements();
//        for(int iel=0; iel<nel; iel++){
//            TPZGeoEl *gel = fGMesh->Element(iel);
//            int nsides = gel->NSides();
//            for(int iside=0; iside<nsides; iside++){
//
//                TPZGeoElSide gelside = gel->Neighbour(iside);
////                std::cout<<gelside.Dimension()<<std::endl;
////                gel->Print();
////                gelside.Print(std::cout);
//                if (gelside.Dimension()==dimension){
//                    bool haskel = HasLowerDimensionNeighbour(gelside);
//                    if(haskel==false){
//                        TPZGeoElBC(gelside,matid);
//
//                        //crea los ribs;
//                        if(gelside.Dimension()==1){
//                            int side = gelside.Side();
//                            int64_t p1 =  gel->SideNodeIndex(side, 0);
//                            int64_t p2 =  gel->SideNodeIndex(side, 1);
//                            TPZVec<REAL> pp1(3);
//                            TPZVec<REAL> pp2(3);
//                            fGMesh->NodeVec()[p1].GetCoordinates(pp1);
//                            fGMesh->NodeVec()[p2].GetCoordinates(pp2);
//                            TPZVec<TPZVec<double>> vecpoints(2);
//                            vecpoints[0]=pp1;
//                            vecpoints[1]=pp2;
///// no se que hacer aqui                       fRibs.resize(nribs+1);
///// no se que hacer aqui                       TRSRibs rib(vecpoints);
//                            nribs=fRibs.size();
//                            fRibs[nribs-1]=TRSRibs(rib);
//                        }
//                    }
//                }
//            }
//        }
//
//    }
//    else{
//        std::cout<<"XXXX";
//    }
//}
//
///**
// * @brief Calculates the intersection point plane-rib
// * @param Point above the plane (node rib 1)
// * @param Point below the plane (node rib 2)
// * @return Intersecting point
// */
//
//TPZVec<double> TRSRibFrac::CalculateIntersection(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2){
//
//    bool check = Check_rib(p1, p2);
//    TPZVec<double> Pint(3,1.0);
//    if(check){
//    double term1 = (p1[0]-p2[0])*fAxis(2,0) + (p1[1]-p2[1])*fAxis(2,1) + (p1[2]-p2[2])*fAxis(2,2);
//    double term2 = ((fCenterCo[0]-p2[0])*(fAxis(2,0))) + ((fCenterCo[1]-p2[1])*(fAxis(2,1))) + ((fCenterCo[2]-p2[2])*(fAxis(2,2)));
//    double alpha = term2/term1;
//    TPZVec<double> Pint(3,0.0);
//    for (int p=0; p<3; p++){
//    Pint[p] = p2[p] + alpha*(p1[p]-p2[p]);
//
//    }
//        std::cout<<"Intersection point: "<<std::endl;
//        std::cout<<Pint[0]<<std::endl;
//        std::cout<<Pint[1]<<std::endl;
//        std::cout<<Pint[2]<<std::endl;
//    }
//    else{
//        Pint[0] = -666;
//        Pint[0] = -666;
//        Pint[0] = -666;
//    }
//    return Pint;
//
//}
//TPZVec<TRSRibs> TRSRibFrac::GetRibsVec(){
//    if(fRibs.size()>0){
//        return fRibvec; /// fRibs es una mapa
//    }
//}
//
//
///**
// * @brief Divide a rib
// * @param Index of the element to divide (for the moment just ribs)
// * @param Point of intersection
// * @return Intersecting point
// */
//
////void TRSRibFrac::DivideRib(int element_index, TPZVec<double> intersection){
////    TPZGeoEl *gel = fmesh->Element(element_index);
////    if(gel->Dimension()!=1){
////        std::cout<<"Just ribs now!";
////        DebugStop();
////    }
////    TPZVec<TPZGeoEl *> gelsDivide(2);
////    gel->Divide(gelsDivide);
////
////}
