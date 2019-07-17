/*! 
 *  @brief     Compares a geomesh with fracture plane to find intersections.
 *  @details   Intersection search is performed after creation of skeleton
 *  elements with TRSFractureMesh::CreateSkeletonElements. Fracture plane should
 *  be a TRSFracPlane.
 *  @authors   Jorge Ordoñez
 *  @authors   Pedro Lima
 *  @date      2018-2019
 */

#include "TRSFractureMesh.h"
#include <math.h>
#include <cstdio>

// Empty Constructor
TRSFractureMesh::TRSFractureMesh(){
}

// Constructor with corner points and a geomesh
TRSFractureMesh::TRSFractureMesh(TRSFracPlane &FracPlane, TPZGeoMesh *gmesh){
    fracplane = FracPlane;
    fGMesh = gmesh;

    // Implement check for previously created skeleton
    // Create skeleton elements
    CreateSkeletonElements(2, 4);
    CreateSkeletonElements(1, 4);

    // Create FracPlane's geometric element into mesh
    int nnodes =  fGMesh->NNodes();
    int ncorners = fracplane.GetCorners().Cols();
    fGMesh->NodeVec().Resize(nnodes + ncorners);
    TPZVec<TPZGeoNode> corners(ncorners);
    TPZVec<int64_t> CornerIndexes(ncorners);
    for (int i = 0; i < ncorners; i++)
    {
        TPZVec<REAL> nodeX(3, 0);
        nodeX[0] = fracplane.GetCorners()(0,i);
        nodeX[1] = fracplane.GetCorners()(1,i);
        nodeX[2] = fracplane.GetCorners()(2,i);

        corners[i].SetCoord(nodeX);
        fGMesh->NodeVec()[nnodes + i] = corners[i];
        CornerIndexes[i] = nnodes + i;
    }
    MElementType elemtype;
    switch (ncorners)
    {
        case 3: elemtype = ETriangle; break;
        case 4: elemtype = EQuadrilateral; break;
        default: DebugStop();
    }
    int64_t nels = fGMesh->NElements();
    fGMesh->CreateGeoElement(elemtype, CornerIndexes, 40, nels);
    

}

// Copy constructor
TRSFractureMesh::TRSFractureMesh(const TRSFractureMesh &copy){
    fGMesh = copy.fGMesh;

    fTolerance = copy.GetTolerance();
    fRibs = copy.fRibs;
    fFaces = copy.fFaces;

    fracplane = copy.fracplane;
}

// Assignment operator
TRSFractureMesh &TRSFractureMesh::operator=(const TRSFractureMesh &copy){
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

TRSFracPlane TRSFractureMesh::GetPlane() const{
    return fracplane;
}

/**
 * @brief Set the tolerance for the distance between a point-plane
 * @param Tolerance
 * @return The tolerance
 */

void TRSFractureMesh::SetTolerance(REAL tolerance){
    fTolerance = tolerance;
}

/**
 * @brief Get the tolerance
 * @return The tolerance
 */

REAL TRSFractureMesh::GetTolerance() const{
    return fTolerance;
}



// /**
//  * @brief Checks if a surface needs to be divided
//  * @param Surface index (integer)
//  * @param Cut ribs vector
//  * @return True if the surface needs to be divided
//  * @return False if the surface does not need to be divided
//  */

// bool TRSFractureMesh::NeedsSurface_Divide(int64_t suface_index, TPZVec<int64_t> interribs) {
//     bool state= false;   //By definition does not need to be divided
//     int nribs = interribs.size();
//     int nribscut =0;
//     for(int i=0; i< nribs; i++){
//         int index_anal = interribs[i];
//         TRSRibs rib_an = fRibs[index_anal];
//         if(rib_an.CutsPlane()==true){
//             nribscut++;
//         }
//     }
//     if(nribscut > 0){     //Checks if a surface has ribs cut
//         return true;   //The surface needs to be divided
//     }
//     return false;
// }

/**
 * @brief Check if the neighbour has a equal dimension
 * @param Geo element side
 * @return True if has a lower dimension
 * @return False if has a upper dimension
 */

bool TRSFractureMesh::HasEqualDimensionNeighbour(TPZGeoElSide &gelside){
    
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

void TRSFractureMesh::CreateSkeletonElements(int dimension, int matid)
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

void TRSFractureMesh::AddRib(TRSRibs rib){
    int index= rib.ElementIndex();
    fRibs[index]=rib;
}

/**
 * @brief Add cut faces using indexes
 * @param Face to be set
 */

void TRSFractureMesh::AddFace(TRSFace face){
    int index= face.ElementIndex();
    fFaces[index]=face;
}

/**
 * @brief Add faces that are cut at the edges of fracture (using indexes)
 * @param Face to be set
 */

void TRSFractureMesh::AddEndFace(TRSFace face){
    int index= face.ElementIndex();
    fEndFaces[index]=face;
}

/**
 * @brief Gives the ribs map
 * @return A map that contains the ribs information
 */

std::map<int64_t ,TRSRibs> TRSFractureMesh::GetRibs(){
    return fRibs;    
}

/**
 * @brief Create cut surfaces
 * @param Material id
 */

void TRSFractureMesh::SplitFaces(int matID){
    //CreateSkeletonElements(1, 4);
    //CreateSkeletonElements(2, 4);
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
            TPZGeoElSide gelside(gel,iside);
            TPZGeoElSide neig = gelside.Neighbour();
            while(neig.Element()->Dimension()!=1){
                neig=neig.Neighbour();
            }
            int rib_index = neig.Element()->Index();
            TRSRibs ribtest = fRibs[rib_index];
            if(ribtest.IsCut()==true){
                ribstatus[iside-nsides] = true;
                // store the rest of the ribs
                cad[nribscut]=rib_index;
                nribscut++;
            }
        }
        
        switch (nribscut)
        {
            case 0: {break;}
            case 2:{ //mid-fracture element
            // Create TRSFace
                TRSFace face(iel, true);
                AddFace(face);
                std::cout<<"first rib: "<<cad[0]<<std::endl;
                std::cout<<"second rib: "<<cad[1]<<std::endl;
                face.SetRibsCut(cad);
                gel->SetMaterialId(matID);
                
            // Connect intersection points
                TPZVec<int64_t> ipoints(2);
                ipoints[0] = Rib(cad[0])->IntersectionIndex();
                ipoints[1] = Rib(cad[1])->IntersectionIndex();
                int64_t nels = fGMesh->NElements();
                fGMesh->CreateGeoElement(EOned, ipoints, 46, nels);

                break;
            }
            case 1:{ //end-fracture element
                TRSFace face(iel, true);
                AddFace(face);
                std::cout<<"single rib cut: "<<cad[0]<<std::endl;
                gel->SetMaterialId(matID+15);
                //Is the fracture skeleton built? Code won't work otherwise.
                TPZVec<REAL> coords = FindEndFracturePoint(face);
                
                // Create geometric element for intersection node
                TPZVec<int64_t> nodeindex(1,0);
                nodeindex[0] = fGMesh->NodeVec().AllocateNewElement();
                fGMesh->NodeVec()[nodeindex[0]].Initialize(coords, *fGMesh);
                int64_t nels = fGMesh->NElements();
                fGMesh->CreateGeoElement(EPoint, nodeindex, 45, nels);

                // Connect intersection points
                nels++;
                TPZVec<int64_t> ipoints(2);
                ipoints[0] = Rib(cad[0])->IntersectionIndex();
                ipoints[1] = nodeindex[0];
                fGMesh->CreateGeoElement(EOned, ipoints, 46, nels);
                break;
            }
            default: {std::cout<<"\nNo more than 2 ribs should've been cut\n";DebugStop();}
        }
    }
}

TRSRibs *TRSFractureMesh::Rib(int index){
    return &fRibs[index];
}

TPZVec<REAL> TRSFractureMesh::FindEndFracturePoint(TRSFace face){
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
    std::cout<<"\n TRSFractureMesh::FindEndFracturePoint\n";
    std::cout << "\nFailed to find intersection point in end-fracture face index: " << face.ElementIndex() << std::endl;
    DebugStop();
    return -1;
}

void TRSFractureMesh::SplitRibs(int matID){
    //search gmesh for cut ribs
    int64_t Nels = fGMesh->NElements();
    for (int iel = 0; iel < Nels; iel++){
        TPZGeoEl *gel = fGMesh->Element(iel);
        int nSides = gel->NSides();
        //skip all elements that aren't ribs
        if (gel->Dimension() != 1){continue;}
        for (int side = 0; side < nSides; side++){
            if (gel->NSideNodes(side) == 2){

                // Get rib's vertexes
                int64_t p1 = gel->SideNodeIndex(side, 0);
                int64_t p2 = gel->SideNodeIndex(side, 1);
                TPZVec<REAL> pp1(3);
                TPZVec<REAL> pp2(3);
                fGMesh->NodeVec()[p1].GetCoordinates(pp1);
                fGMesh->NodeVec()[p2].GetCoordinates(pp2);

                // Check rib
                bool resul = fracplane.Check_rib(pp1, pp2);

                // Split rib
                if (resul == true){
                    TRSRibs rib(iel, true);
                    AddRib(rib);
                    TPZVec<REAL> ipoint = fracplane.CalculateIntersection(pp1, pp2);
                    //5O is the material of children ribs
                    Rib(iel)->DivideRib(fGMesh, ipoint, matID);
                    gel->SetMaterialId(12);

                    // std::cout<<"Element: "<<iel<<" Side: "<<side<<" Rib status: "<<resul<<" Fracture : 1"<<"\n";
                }
            }
        }
    }
}

