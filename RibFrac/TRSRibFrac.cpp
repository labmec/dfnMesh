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

// Constructor by using a data matrix and a geomesh
TRSRibFrac::TRSRibFrac(const Matrix &data, TPZGeoMesh *gmesh){
    if(!Check_ConsistencyData(data)){
        std::cout<<"The input data is not correct";
        DebugStop();
    }
    fCornerCoordinates=data;
    fGMesh=gmesh;
}

// Copy constructor
TRSRibFrac::TRSRibFrac(const TRSRibFrac &copy){
    fTolerance = copy.GetTolerance();
    fRibs = copy.fRibs;
    
    fCornerCoordinates=copy.GetPlane();
    fAxis=copy.fAxis;
    fCenterCo = copy.fCenterCo;
    fGMesh=copy.fGMesh;
}

// Assignment operator
TRSRibFrac &TRSRibFrac::operator=(const TRSRibFrac &copy){
    fTolerance = copy.GetTolerance();
    fRibs = copy.fRibs;
    fCornerCoordinates=copy.GetPlane();
    fAxis=copy.fAxis;
    fCenterCo = copy.fCenterCo;
    fGMesh=copy.fGMesh;
    return *this;
}

/**
 * @brief Set the fracture plane
 * @param Fracture plane coordinates (Matrix 3x4)
 */

void TRSRibFrac::SetPlane(Matrix plane){
    if(!Check_ConsistencyData(plane)){
        std::cout<<"The input data is not correct";
        DebugStop();
    }
    fCornerCoordinates=plane;
}

/**
 * @brief Get the fracture plane
 * @return The plane corner coordinates
 */

Matrix TRSRibFrac::GetPlane() const{
    return fCornerCoordinates;
}

/**
 * @brief Set the tolerance for the distance between a point-plane
 * @param Tolerance
 * @return The tolerance
 */

void TRSRibFrac::SetTolerance(REAL tolerance){
    fTolerance = tolerance;
}

/**
 * @brief Get the tolerance
 * @return The tolerance
 */

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
    
    if(rows != 3){                        //Should be 3 (x y z)
        std::cout<<"Check the input data";
        DebugStop();
    }
    if(cols < 3 or cols > 4){//To form a plane it is needed at least 3 points but no more than 4
        std::cout<<"Check the input data (number of plane points, four is enough)";
        DebugStop();
    }
   Matrix MidPoint;                       //Mid points computation
    if (!(cols ==4 && rows ==3)){
        std::cout<<"Check the input data (number of plane points, four is enough)";
        DebugStop();
    }
    MidPoint.Resize(3,cols);              //3 rows (x y z), mid points

// Mid points computation for the first three points
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
    for (int i=0; i<3; i++) {           // i< axis number (x y z)
        double norm=ax(i,1) - (ax(i,0))*((ax(0,0)*ax(0,1)) + (ax(1,0)*ax(1,1)) + (ax(2,0)*ax(2,1)));
        ax(i,1)=norm;
    }

//Ax2 computation
    ax(0,2)=ax(1,0)*ax(2,1) - ax(2,0)*ax(1,1);
    ax(1,2)=ax(2,0)*ax(0,1) - ax(0,0)*ax(2,1);
    ax(2,2)=ax(0,0)*ax(1,1) - ax(1,0)*ax(0,1);

//Coplanar verification
    double ver = ax(0,2)*(CornerCoordinates(0,2)-CornerCoordinates(0,0))+ax(1,2)*(CornerCoordinates(1,2)-CornerCoordinates(1,0))+ax(2,2)*(CornerCoordinates(2,2)-CornerCoordinates(2,0));

//Cheks if the points are coplanar
    if(abs(ver) > fTolerance){
        std::cout<<"The input points are not coplanar"<<"\n"<<std::endl;
        DebugStop();
    }
        
// After checking the consistency the date is set
        fCornerCoordinates=CornerCoordinates;
        fAxis=ax;
        fCenterCo.Resize(3);
// Center point computation
    
        for(int i=0; i< 4; i++){          /// i< axis number (x y z)
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
    
    //Point distance to the fracture plane computation
        double point_distance = (point[0] - fCenterCo[0])*(fAxis.GetVal(0,2)) + (point[1] - fCenterCo[1])*(fAxis.GetVal(1,2)) + (point[2] - fCenterCo[2])*(fAxis.GetVal(2,2));
        if (point_distance>0){
            return true;    //If the point is above de plane
        }
        else{
            return false;   //If the point is below de plane
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
    
        if(Check_point_above(p1) != Check_point_above(p2)){
            return true;    //Rib cut by the plane
        }
        else
        {
            return false;    //Rib does not cut by the plane
        }
}

/**
 * @brief Checks if a rib is cut by a fracture plane
 * @param Surface index (integer)
 * @param Cut ribs vector
 * @return True if the surface needs to be divided
 * @return False if the surface does not need to be divided
 */

bool TRSRibFrac::NeedsSurface_Divide(int64_t suface_index, TPZVec<int64_t> interribs) {
    bool state= false;   //By definition does not need to be divided
    int nribs = interribs.size();
    int count =0;
    for(int i=0; i< nribs; i++){
        int index_anal = interribs[i];
        TRSRibs rib_an = fRibs[index_anal];
        if(rib_an.CutsPlane()==true){
            count++;
        }
    }
    if(count == 2){     //Checks if a surface has 2 ribs cut
        state = true;   //The surface needs to be divided
    }
    if (!(count == 0 or count ==2)){
        DebugStop();//Otherwise the surface does not need to be divided
    }
    return state;       //Returns true or false
}

/**
 * @brief Checks if a rib is within the fracture plane
 * @param Point vector with the euclidean coordinates
 * @return True if the rib is within fracture plane
 * @return False if the rib is without the fracture plane
 */

bool TRSRibFrac::RibInPlane(TPZVec<REAL> point){
    int res;
        res = ((point[0]-fAxis(0,2))*fCenterCo[0])+((point[1]-fAxis(1,2))*fCenterCo[1])+((point[2]-fAxis(2,2))*fCenterCo[2]);

    if (res==0){
        return true;
    }
    else{
        return false;
    }
}

/**
 * @brief Check if the neighbour has a equal dimension
 * @param Geo element side
 * @return True if has a lower dimension
 * @return False if has a upper dimension
 */

bool TRSRibFrac::HasEqualDimensionNeighbour(TPZGeoElSide &gelside){
    
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
                        bool haskel = HasEqualDimensionNeighbour(gelside);
                            if(haskel==false){
                                int nel_mesh = fGMesh->NElements();
                                
                                    TPZGeoElBC(gelside,matid);

                                    //Setting Ribs
                                    TRSRibs rib(nel_mesh,false);
                                    AddRib(rib);
                }
            }
        }
    }
}

/**
 * @brief Calculates the intersection point plane-rib
 * @param Point above the plane (vector)
 * @param Point below the plane (vector)
 * @return Intersecting point
 */

TPZVec<double> TRSRibFrac::CalculateIntersection(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2){
    
    bool check = Check_rib(p1, p2);     //Checks if the rib is cut
    TPZVec<double> Pint;
    if(check){
        double term1 = (p1[0]-p2[0])*fAxis(0,2) + (p1[1]-p2[1])*fAxis(1,2) + (p1[2]-p2[2])*fAxis(2,2);
        double term2 = ((fCenterCo[0]-p2[0])*(fAxis(0,2))) + ((fCenterCo[1]-p2[1])*(fAxis(1,2))) + ((fCenterCo[2]-p2[2])*(fAxis(2,2)));
        double alpha = term2/term1;
        Pint.Resize(3);
        for (int p=0; p<3; p++){
            Pint[p] = p2[p] + alpha*(p1[p]-p2[p]);
        }
//        std::cout<<"Intersection point: "<<std::endl;
//        std::cout<<Pint[0]<<std::endl;
//        std::cout<<Pint[1]<<std::endl;
//        std::cout<<Pint[2]<<std::endl;
    }
    
    return Pint;
}


/**
 * @brief Sets the rib idexes in the map
 * @param Ribs to be set
 */

void TRSRibFrac::AddRib(TRSRibs rib){
    int index= rib.ElementIndex();
    fRibs[index]=rib;
}

/**
 * @brief Add cut faces using indexes
 * @param Face to be set
 */

void TRSRibFrac::AddFace(TRSFace face){
    int index= face.ElementIndex();
    fFaces[index]=face;
}

/**
 * @brief Gives the ribs map
 * @return A map that contains the ribs information
 */

std::map<int64_t ,TRSRibs> TRSRibFrac::GetRibs(){
    return fRibs;
    
    /// Needs to be add to a method for further coding
    
    
//    int nel=fGMesh->NElements();
//    for(int iel=0; iel<nel; iel++){
//        TPZGeoEl *gel = fGMesh->Element(iel);
//        if(gel->Dimension()!=1){continue;}
//        TPZFMatrix<REAL> cooridnates;
//        gel->NodesCoordinates(cooridnates);
//        TPZVec<REAL>p1;
//        TPZVec<REAL>p2;
//        p1[0]=cooridnates(0,0);
//        p1[1]=cooridnates(1,0);
//        p1[2]=cooridnates(2,0);
//        p2[0]=cooridnates(0,1);
//        p2[1]=cooridnates(1,1);
//        p2[2]=cooridnates(2,1);
//        bool check = Check_rib(p1, p2);
//        if(!check){continue;}
//        TPZVec<REAL> point;
//        point = CalculateIntersection(p1, p2);
//        bool check2 = RibInPlane(point);
//        if(!check2){continue;}
//        TRSRibs rib(iel,check2);
//        TPZVec<TPZGeoEl *> gels;
//
//      //  rib.DivideRib(fGMesh, 100);
//        gel->SetMaterialId(100);
//
//
//    }
    
    
}

/**
 * @brief Create cut surfaces
 * @param Material id
 */

void TRSRibFrac::CreateSurfaces(int matID){
    int nel = fGMesh->NElements();
    
    
    for(int iel = 0; iel<nel; iel++){
        int count =0;
        TPZManVector<int64_t,2> cad;
        cad.Resize(2);
        TPZGeoEl *gel = fGMesh->Element(iel);
        int dim = gel->Dimension();
        if (dim == 1){continue;}
        if(gel->MaterialId()==4){continue;}
        for(int iside=4; iside<8; iside++){
            TPZGeoElSide gelside(gel,iside);
            TPZGeoElSide neig = gelside.Neighbour();
            while(neig.Element()->Dimension()!=1){
                neig=neig.Neighbour();
            }
            int rib_index = neig.Element()->Index();
            TRSRibs ribstatus = fRibs[rib_index];
            if(ribstatus.CutsPlane()==true){
                cad[count]=rib_index;
                count++;
                
                //verificar como hacer para setar despues
                
            }
        }
        
        if(count == 2){
            TRSFace face(iel, true);
            std::cout<<"primer rib: "<<cad[0]<<std::endl;
            std::cout<<"segundo rib: "<<cad[1]<<std::endl;
            face.SetRibsInSurface(cad);
            AddFace(face);
            
            gel->SetMaterialId(10);
        }
        if(!(count==0 or count ==2)){
            DebugStop();
        }
    }
  
}

TRSRibs *TRSRibFrac::Rib(int index){
    
    return &fRibs[index];
}










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
