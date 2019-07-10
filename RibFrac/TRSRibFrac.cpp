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

// /**
//  * @brief Add faces that are cut at the edges of fracture (using indexes)
//  * @param Face to be set
//  */

// void TRSRibFrac::AddEndFace(TRSFace face){
//     int index= face.ElementIndex();
//     fEndFaces[index]=face;
// }

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
    int64_t nel = fGMesh->NElements();
    
    
    for(int iel = 0; iel<nel; iel++){
        int nribscut =0;
        TPZManVector<bool,4> ribstatus(4,false);
        TPZManVector<int64_t,2> cad(2);
        TPZGeoEl *gel = fGMesh->Element(iel);
        int dim = gel->Dimension();
        if (dim != 2){continue;}
        //40 is MaterialID for fracture plane
        if(gel->MaterialId()==40){continue;}
        int nsides = gel->NNodes();
        for(int iside = nsides; iside < 2*nsides; iside++){
            //ribstatus[iside-nsides] = false;
            TPZGeoElSide gelside(gel,iside);
            TPZGeoElSide neig = gelside.Neighbour();
            while(neig.Element()->Dimension()!=1){
                neig=neig.Neighbour();
            }
            int rib_index = neig.Element()->Index();
            TRSRibs ribtest = fRibs[rib_index];
            if(ribtest.CutsPlane()==true){
                ribstatus[iside-nsides] = true;
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
                AddFace(face);
                std::cout<<"single rib cut: "<<cad[0]<<std::endl;
                gel->SetMaterialId(matID+15);
                //Is the fracture skeleton built? Code will break otherwise.
                TPZVec<REAL> coords = FindEndFracturePoint(face);

                // int64_t nnodes = fGMesh->NNodes();
                // fGMesh->NodeVec().Resize(nnodes+1);
                // TPZGeoNode ipoint(nnodes, coords, *fGMesh);
                
                // Create geometric element for intersection node
                TPZVec<int64_t> nodeindex(1,0);
                nodeindex[0] = fGMesh->NodeVec().AllocateNewElement();
                fGMesh->NodeVec()[nodeindex[0]].Initialize(coords, *fGMesh);
                int64_t nels = fGMesh->NElements();
                fGMesh->CreateGeoElement(EPoint, nodeindex, 45, nels);
                break;
            }
            default: {DebugStop(); break;}
        }
    }
}

TRSRibs *TRSRibFrac::Rib(int index){
    return &fRibs[index];
}

TPZVec<REAL> TRSRibFrac::FindEndFracturePoint(TRSFace face){
    // Convert TPZGeoEl into TRSFracPlane
    TPZGeoEl *gelface = fGMesh->Element(face.ElementIndex());
    TPZFMatrix<REAL> corners;
    gelface->NodesCoordinates(corners);
    TRSFracPlane faceplane(corners);

    // Check fracplane's ribs for intersection with faceplane
    int nribs = fracplane.GetCorners().Cols();
    for(int irib = 0; irib < nribs; irib++){
        TPZVec<REAL> p1(3);
        TPZVec<REAL> p2(3);
        for(int i = 0; i<3; i++){
            p1[i] = fracplane.GetCorners()(i, irib);
            p2[i] = fracplane.GetCorners()(i, (irib+1)%nribs);
        }
        if(faceplane.Check_rib(p1, p2)){
            return faceplane.CalculateIntersection(p1,p2);
        }
    }
    std::cout << std::endl << "Failed to find intersection point in end-fracture face index: " << face.ElementIndex() << std::endl;
    DebugStop();
    return -1;
}

