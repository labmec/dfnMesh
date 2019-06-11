//
// RibFrac.cpp
// RibFrac
//
//  Created by JORGE ORDOÑEZ on 22/5/18.
//  Modified By: Pedro Lima
//

#include "TRSRibFrac.h"
#include <math.h>
#include <cstdio>

// Empty Constructor
TRSRibFrac::TRSRibFrac(){
}

// Constructor with corner points and a geomesh
TRSRibFrac::TRSRibFrac(TRSFracPlane &FracPlane, TPZGeoMesh *gmesh){
    fracplane = FracPlane;
    fGMesh=gmesh;
    //printf("constr fracplane.L0 = %.3f \n", fracplane.L0);
}

// Copy constructor
TRSRibFrac::TRSRibFrac(const TRSRibFrac &copy){
    fGMesh=copy.fGMesh;

    fTolerance = copy.GetTolerance();
    fRibs = copy.fRibs;
    fFaces = copy.fFaces;

    fracplane = copy.fracplane;
}

// Assignment operator
TRSRibFrac &TRSRibFrac::operator=(const TRSRibFrac &copy){
    fGMesh=copy.fGMesh;

    fTolerance = copy.GetTolerance();
    fRibs = copy.fRibs;
    fFaces = copy.fFaces;

    fracplane = copy.fracplane;
    return *this;
}

/**
 * @brief Get the fracture plane
 * @return The plane corner coordinates
 */

TRSFracPlane TRSRibFrac::GetPlane() const{
    return fracplane;
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
 * @brief Checks if a point is above or below the fracture plane
 * @param Point vector with the euclidean coordinates
 * @return True if the point is above the fracture plane
 * @return False if the point is below the fracture plane
 */

bool TRSRibFrac::Check_point_above(const TPZVec<REAL> &point) const{
    
    //Point distance to the fracture plane computation
        double point_distance = (point[0] - fracplane.fCenterCo[0])*((fracplane.fAxis).GetVal(0,2)) 
                                + (point[1] - fracplane.fCenterCo[1])*((fracplane.fAxis).GetVal(1,2)) 
                                + (point[2] - fracplane.fCenterCo[2])*((fracplane.fAxis).GetVal(2,2));
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

bool TRSRibFrac::Check_rib(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2) {
        //check for infinite plane
        if(Check_point_above(p1) != Check_point_above(p2)){
            //Rib cut by infinite plane
            //then calculate intersection point and check if it's within plane boundaries
            TPZVec<REAL> intersection = CalculateIntersection(p1, p2);
            return IsPointInPlane(intersection);
        }
        else
        {
            return false;    //Rib is not cut by plane
        }
}


/**
 * @brief Checks if a point is within fracture plane
 * @param Point vector with the euclidean coordinates
 * @return True if the point is within fracture plane
 * @return False if the point is out of fracture plane
 */

bool TRSRibFrac::IsPointInPlane(TPZVec<REAL> &point) 
{
    double  dist = fabs(
                    (point[0] - fracplane.fCenterCo[0])*(fracplane.fAxis).GetVal(0, 0)
                    +(point[1] - fracplane.fCenterCo[1])*(fracplane.fAxis).GetVal(1,0)
                    +(point[2] - fracplane.fCenterCo[2])*(fracplane.fAxis).GetVal(2,0));
    if (dist > fracplane.L0/2){return false;}
            dist = fabs(
                    (point[0] - fracplane.fCenterCo[0])*(fracplane.fAxis).GetVal(0,1)
                    +(point[1] - fracplane.fCenterCo[1])*(fracplane.fAxis).GetVal(1,1)
                    +(point[2] - fracplane.fCenterCo[2])*(fracplane.fAxis).GetVal(2,1));
    if (dist > fracplane.L1/2){return false;}
    // std::cout<<" ___ ";
    return true;
}

/**
 * @brief Checks if a surface needs to be divided
 * @param Surface index (integer)
 * @param Cut ribs vector
 * @return True if the surface needs to be divided
 * @return False if the surface does not need to be divided
 */

bool TRSRibFrac::NeedsSurface_Divide(int64_t suface_index, TPZVec<int64_t> interribs) {
    bool state= false;   //By definition does not need to be divided
    int nribs = interribs.size();
    int nribscut =0;
    for(int i=0; i< nribs; i++){
        int index_anal = interribs[i];
        TRSRibs rib_an = fRibs[index_anal];
        if(rib_an.CutsPlane()==true){
            nribscut++;
        }
    }
    if(nribscut > 0){     //Checks if a surface has ribs cut
        return true;   //The surface needs to be divided
    }
    return false;
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
    return false;    
}

/**
  * @brief Creates the skeleton mesh
  * @param Dimension
  * @param Material ID number
  */

void TRSRibFrac::CreateSkeletonElements(int dimension, int matid)
{
    int nel = fGMesh->NElements();
    for (int iel = 0; iel < nel; iel++)
    {
        TPZGeoEl *gel = fGMesh->Element(iel);
        int nsides = gel->NSides();
        for (int iside = 0; iside < nsides; iside++)
        {
            TPZGeoElSide gelside = gel->Neighbour(iside);
            if (gelside.Dimension() == dimension)
            {
                bool haskel = HasEqualDimensionNeighbour(gelside);
                if (haskel == false)
                {
                    int nel_mesh = fGMesh->NElements();
                    TPZGeoElBC(gelside, matid);
                    switch (dimension)
                    {
                        case 1:{
                            TRSRibs rib(nel_mesh, false);
                            AddRib(rib);
                            break;}
                        case 2:{
                            TRSFace face(nel_mesh, false);
                            AddFace(face);
                            break;}
                        default: {DebugStop();}
                    }
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

TPZVec<double> TRSRibFrac::CalculateIntersection(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2)
{     
    TPZVec<double> Pint;
    
    double term1 = (p1[0]-p2[0])*fracplane.fAxis(0,2) + (p1[1]-p2[1])*fracplane.fAxis(1,2) + (p1[2]-p2[2])*fracplane.fAxis(2,2);
    double term2 = ((fracplane.fCenterCo[0]-p2[0])*(fracplane.fAxis(0,2))) + ((fracplane.fCenterCo[1]-p2[1])*(fracplane.fAxis(1,2))) + ((fracplane.fCenterCo[2]-p2[2])*(fracplane.fAxis(2,2)));
    double alpha = term2/term1;
    Pint.Resize(3);
    for (int p=0; p<3; p++){
        Pint[p] = p2[p] + alpha*(p1[p]-p2[p]);
    }
    //  std::cout<<"Intersection point: "<<std::endl;
    //  std::cout<<Pint[0]<<std::endl;
    //  std::cout<<Pint[1]<<std::endl;
    //  std::cout<<Pint[2]<<std::endl;
    
    
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
//        TPZFMatrix<REAL> coordinates;
//        gel->NodesCoordinates(coordinates);
//        TPZVec<REAL>point;
//        TPZVec<REAL>p2;
//        point[0]=coordinates(0,0);
//        point[1]=coordinates(1,0);
//        point[2]=coordinates(2,0);
//        p2[0]=coordinates(0,1);
//        p2[1]=coordinates(1,1);
//        p2[2]=coordinates(2,1);
//        bool check = Check_rib(point, p2);
//        if(!check){continue;}
//        TPZVec<REAL> point;
//        point = CalculateIntersection(point, p2);
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
        int nribscut =0;
        TPZManVector<int64_t,2> cad;
        cad.Resize(2);
        TPZGeoEl *gel = fGMesh->Element(iel);
        int dim = gel->Dimension();
        if (dim != 2){continue;}
        //40 is MaterialID for fracture plane
        if(gel->MaterialId()==40){continue;}
        for(int iside=4; iside<8; iside++){
            TPZGeoElSide gelside(gel,iside);
            TPZGeoElSide neig = gelside.Neighbour();
            while(neig.Element()->Dimension()!=1){
                neig=neig.Neighbour();
            }
            int rib_index = neig.Element()->Index();
            TRSRibs ribstatus = fRibs[rib_index];
            if(ribstatus.CutsPlane()==true){
                cad[nribscut]=rib_index;
                nribscut++;
            }
        }
        
        if(nribscut > 0){
            TRSFace face(iel, true);
            AddFace(face);
            
            if(nribscut == 2)
            //mid-fracture element
            {
                std::cout<<"first rib: "<<cad[0]<<std::endl;
                std::cout<<"second rib: "<<cad[1]<<std::endl;
                face.SetRibsInSurface(cad);
                gel->SetMaterialId(matID);
            }
            else if(nribscut == 1)
            //end-fracture element
            {
                std::cout<<"single rib cut: "<<cad[0]<<std::endl;
                gel->SetMaterialId(matID+15);
            }
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
