/*! 
 *  @brief     Compares a geomesh with fracture plane to find intersections.
 *  @details   Intersection search is performed after creation of skeleton
 *  elements with TRSRibFrac::CreateSkeletonElements. Fracture plane should
 *  be a TRSFracPlane.
 *  @authors   Jorge Ordoñez
 *  @authors   Pedro Lima
 *  @date      2018-2019
 */

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
    fGMesh = copy.fGMesh;

    fTolerance = copy.GetTolerance();
    fRibs = copy.fRibs;
    fFaces = copy.fFaces;

    fracplane = copy.fracplane;
}

// Assignment operator
TRSRibFrac &TRSRibFrac::operator=(const TRSRibFrac &copy){
    fGMesh = copy.fGMesh;

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
        double point_distance = (point[0] - fracplane.GetCorners()(0,1))*((fracplane.axis()).GetVal(0,2)) 
                                +(point[1] - fracplane.GetCorners()(1,1))*((fracplane.axis()).GetVal(1,2)) 
                                +(point[2] - fracplane.GetCorners()(2,1))*((fracplane.axis()).GetVal(2,2));
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
        TPZVec<REAL> intersection = CalculateIntersection(fracplane, p1, p2);
        return IsPointInPlane(fracplane, intersection);
    }
    else
    {
        return false;    //Rib is not cut by plane
    }
}


/**
 * @brief Checks if a point is within fracture plane
 * @details Enumeration of corner points should follow standard PZ topology, where 
 * corner nodes are numbered counter-clockwise. (This condition will automatically be 
 * met for triangles, but not always for quadrilaterals)
 * @param Point vector with the euclidean coordinates
 * @return True if the point is within fracture plane
 * @return False if the point is out of fracture plane
 */

bool TRSRibFrac::IsPointInPlane(TRSFracPlane &plane, TPZVec<REAL> &point) 
{
    int ncorners = plane.GetCorners().Cols();
    REAL area = 0;
    for(int i = 0; i<ncorners-1; i++){
        //Define vectors from the point to a each one of a pair of corners
        TPZManVector<REAL, 3> ax1(3);
            ax1[0] = plane.GetCorners()(0,i) - point[0];
            ax1[1] = plane.GetCorners()(1,i) - point[1];
            ax1[2] = plane.GetCorners()(2,i) - point[2];
        TPZManVector<REAL, 3> ax2(3);
            ax2[0] = plane.GetCorners()(0,i+1) - point[0];
            ax2[1] = plane.GetCorners()(1,i+1) - point[1];
            ax2[2] = plane.GetCorners()(2,i+1) - point[2];
        //Compute area of trangle outlined by these vectors
        REAL temp = pow(ax1[1]*ax2[2] - ax1[2]*ax2[1],2);
            temp += pow(ax1[2]*ax2[0] - ax1[0]*ax2[2],2);
            temp += pow(ax1[0]*ax2[1] - ax1[1]*ax2[0],2);
                  
        area += sqrt(temp)/2;
    }
    //Then, once more to get initial and last corners
    TPZManVector<REAL, 3> ax1(3);
        ax1[0] = plane.GetCorners()(0,ncorners-1) - point[0];
        ax1[1] = plane.GetCorners()(1,ncorners-1) - point[1];
        ax1[2] = plane.GetCorners()(2,ncorners-1) - point[2];
    TPZManVector<REAL, 3> ax2(3);
        ax2[0] = plane.GetCorners()(0,0) - point[0];
        ax2[1] = plane.GetCorners()(1,0) - point[1];
        ax2[2] = plane.GetCorners()(2,0) - point[2];
    REAL temp = pow(ax1[1]*ax2[2] - ax1[2]*ax2[1],2);
        temp += pow(ax1[2]*ax2[0] - ax1[0]*ax2[2],2);
        temp += pow(ax1[0]*ax2[1] - ax1[1]*ax2[0],2);
                
    area += sqrt(temp)/2;
    // std::cout<<" ___ ";

    //If total computed area is equal to the plane's area, then
    //point is in plane
    return( fabs(area-plane.area()) < fTolerance );
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
        }
        neighbour = neighbour.Neighbour();
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

        //40 is material id for fracture plane, and it should only have an 1D skeleton.
        if(gel->MaterialId() == 40 && dimension == 2){continue;}

        int nsides = gel->NSides();
        for (int iside = 0; iside < nsides; iside++)
        {
            TPZGeoElSide gelside = gel->Neighbour(iside);

            if (gelside.Dimension() != dimension){continue;}
            bool haskel = HasEqualDimensionNeighbour(gelside);
            if (haskel == false)
            {
                int nel_mesh = fGMesh->NElements();
                TPZGeoElBC(gelside, matid);
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
TPZVec<double> TRSRibFrac::CalculateIntersection(TRSFracPlane &plane, const TPZVec<REAL> &p1, const TPZVec<REAL> &p2)
{     
    TPZVec<double> Pint;

    double term1 = ((plane.GetCorners()(0,1)-p2[0])*(plane.axis(0,2))) 
                   +((plane.GetCorners()(1,1)-p2[1])*(plane.axis(1,2)))
                   +((plane.GetCorners()(2,1)-p2[2])*(plane.axis(2,2)));
    double term2 = (p1[0]-p2[0])*plane.axis(0,2) 
                   +(p1[1]-p2[1])*plane.axis(1,2) 
                   +(p1[2]-p2[2])*plane.axis(2,2);
    double alpha = term1/term2;
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
 * @brief Add faces that are cut at the edges of fracture (using indexes)
 * @param Face to be set
 */

void TRSRibFrac::AddEndFace(TRSFace face){
    int index= face.ElementIndex();
    fEndFaces[index]=face;
}

/**
 * @brief Gives the ribs map
 * @return A map that contains the ribs information
 */

std::map<int64_t ,TRSRibs> TRSRibFrac::GetRibs(){
    return fRibs;    
}

/**
 * @brief Create cut surfaces
 * @param Material id
 */

void TRSRibFrac::CreateSurfaces(int matID){
    //CreateSkeletonElements(fGMesh->Dimension(), 4);
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
        
        switch (nribscut)
        {
            case 0: {break;}
            case 2:{ //mid-fracture element
                TRSFace face(iel, true);
                AddFace(face);
                std::cout<<"first rib: "<<cad[0]<<std::endl;
                std::cout<<"second rib: "<<cad[1]<<std::endl;
                face.SetRibsInSurface(cad);
                gel->SetMaterialId(matID);
                break;
            }
            case 1:{ //end-fracture element
                TRSFace face(iel, true);
                AddEndFace(face);
                std::cout<<"single rib cut: "<<cad[0]<<std::endl;
                gel->SetMaterialId(matID+15);
                //Is the fracture skeleton built? Code will break otherwise.
                //TPZVec<REAL> coords = FindEndFracturePoint(face);
                TPZGeoNode ipoint;
                break;
            }
            default: {DebugStop(); break;}
        }
    }
}

TRSRibs *TRSRibFrac::Rib(int index){
    
    return &fRibs[index];
}

// TPZVec<REAL> TRSRibFrac::FindEndFracturePoint(TRSFace face){

// }

