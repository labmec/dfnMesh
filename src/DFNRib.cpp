/*! 
 *	DFNRib.cpp
 *  @authors   Pedro Lima
 *  @authors   Jorge Ordoñez
 *  @date      2018-2019
 */

#include "DFNRib.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzbiharmonic.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"
#include "TPZRefPatternDataBase.h"
#include "DFNFracture.h"
#include <math.h>

/// Default constructor takes a pointer to geometric element
DFNRib::DFNRib(TPZGeoEl *gel, DFNFracture *Fracture) :
    fGeoEl(gel),
    fFracture(Fracture),
    fStatus({0,0,0})
{};

// Copy constructor
DFNRib::DFNRib(const DFNRib &copy){
    this->operator=(copy);
}

// Assignment operator
DFNRib &DFNRib::operator=(const DFNRib &copy){
    fGeoEl = copy.fGeoEl;
    fStatus = copy.fStatus;
    fCoord = copy.fCoord;
    return *this;
}

inline REAL VectorNorm(TPZManVector<REAL,3> &vector);
inline REAL Distance(TPZManVector<REAL,3> &vector1, TPZManVector<REAL,3> &vector2);


/// Get real intersection coordinates
TPZManVector<REAL, 3> DFNRib::RealCoord(){
    TPZManVector<REAL, 3> coord(3, 0);
    TPZGeoMesh *gmesh = fGeoEl->Mesh();

    if(fStatus[2] && fIntersectionIndex >= 0){
        gmesh->NodeVec()[fIntersectionIndex].GetCoordinates(coord);
    }
    else if(fStatus[2]){coord = fCoord;}
    else if(fStatus[0]){fGeoEl->NodePtr(0)->GetCoordinates(coord);}
    else if(fStatus[1]){fGeoEl->NodePtr(1)->GetCoordinates(coord);}
    
    return coord;
}


void DFNRib::Refine(){
    if(!this->NeedsRefinement()) return;
    if(!fGeoEl){
        std::cout<<"No gel associated to the Rib\n";
        DebugStop();
    }

    // set new node in gmesh if it hasn't been set yet
    TPZGeoMesh *gmesh = fGeoEl->Mesh();
    if(fIntersectionIndex < 0){
        fIntersectionIndex = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[fIntersectionIndex].Initialize(fCoord,*gmesh);
    }
    
    // Set refinement pattern
    this->CreateRefPattern();

    // set children
        TPZManVector<int64_t,2> cornerindices(2,0);
        TPZGeoEl *child;
    // first child
        cornerindices[0] = fGeoEl->NodeIndex(0);
        cornerindices[1] = fIntersectionIndex;
        int64_t elindex = gmesh->NElements();
        child = gmesh->CreateGeoElement(EOned, cornerindices, fGeoEl->MaterialId(), elindex);
        fGeoEl->SetSubElement(0,child);
        child->SetFatherIndex(fGeoEl->Index());
    // second child
        cornerindices[0] = fIntersectionIndex;
        cornerindices[1] = fGeoEl->NodeIndex(1);
        elindex++;
        child = gmesh->CreateGeoElement(EOned, cornerindices, fGeoEl->MaterialId(), elindex);
        fGeoEl->SetSubElement(1,child);
        child->SetFatherIndex(fGeoEl->Index());
}



// ab




void DFNRib::CreateRefPattern(){
        // refinement mesh
        TPZGeoMesh refPatternMesh;
        int refnnodes = 2+fStatus[2];
        TPZGeoMesh *gmesh = fGeoEl->Mesh();

        // set nodes
        refPatternMesh.NodeVec().Resize(refnnodes);
        TPZManVector<REAL,3> coord(3);
        for(int i = 0; i<3; i++){
            if(i < 2) {gmesh->NodeVec()[fGeoEl->NodeIndex(i)].GetCoordinates(coord);}
            else      {gmesh->NodeVec()[fIntersectionIndex].GetCoordinates(coord);}

            refPatternMesh.NodeVec()[i].Initialize(coord,refPatternMesh);
        }
        
        // insert father
        TPZManVector<int64_t,2> cornerindices({0,1});
        int64_t elindex = 0;
        refPatternMesh.CreateGeoElement(EOned, cornerindices, DFNMaterial::Erefined, elindex);
        
        // insert children
        cornerindices = {0,2};
        elindex = 1;
        refPatternMesh.CreateGeoElement(EOned, cornerindices, DFNMaterial::Erefined, elindex);
        cornerindices = {2,1};
        elindex = 2;
        refPatternMesh.CreateGeoElement(EOned, cornerindices, DFNMaterial::Erefined, elindex);
                    
        // define refPattern
        refPatternMesh.BuildConnectivity(); 
        TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(refPatternMesh);
        fGeoEl->SetRefPattern(refpat);
    }

/**
 * @brief Check geometry of intersection against a tolerance, snaps intersection 
 * to closest side(s) if necessary and modifies affected neighbours.
 * @return True if any optimization has been made.
 * @param tolDist: Minimum acceptable length
*/
bool DFNRib::Optimize(REAL tolDist){
    if(!this->NeedsRefinement()) return false;
    
    TPZManVector<REAL,3> coord(3,0);
    fGeoEl->NodePtr(0)->GetCoordinates(coord);
    REAL dist0 = Distance(fCoord,coord);

    fGeoEl->NodePtr(1)->GetCoordinates(coord);
    REAL dist1 = Distance(fCoord,coord);

    int64_t closestnode = (dist0 < dist1 ? 0 : 1);
    REAL dist = MIN(dist0,dist1);

    // if(dist<tolDist){
    if(fStatus[0] || fStatus[1] || dist<tolDist){ // this condition seems to be more robust
        fIntersectionIndex = fGeoEl->NodeIndex(closestnode);
        fStatus[closestnode] = 1;
        fStatus[2] = 0;
        UpdateNeighbours(2);            // clear origin
        UpdateNeighbours(closestnode);  // flag destination
        UpdateMaterial(); //@todo I might find this to be unnecessary
        return true;
    }
    return false;
}





void DFNRib::UpdateNeighbours(int iside){
    // If neighbour doesn't have a DFNElement associated to it, create
    // Go through every neighbour and match the iside entry on their status vector
    // And run optimization

    TPZGeoElSide gelside(fGeoEl,iside);
    TPZGeoElSide neighbour = gelside.Neighbour();
    TPZGeoEl *gel;
    for(/*void*/; neighbour != gelside; neighbour = neighbour.Neighbour()){
        gel = neighbour.Element();
        if(gel->HasSubElement()) continue;
        int neig_side = neighbour.Side();

        switch(gel->Dimension()){
            case 0: {continue;}
            case 1:{
                // check if DFNRib exists
                DFNRib *rib_ptr = fFracture->Rib(gel->Index());
                if(!rib_ptr){
                    DFNRib neig_rib(gel,fFracture);
                    fFracture->AddRib(neig_rib);
                    rib_ptr = fFracture->Rib(gel->Index());
                }
                // skip if StatusVec entry already match
                if(fStatus[iside] == rib_ptr->StatusVec()[neig_side]){continue;}
                // else, match StatusVec entry for iside and optimize
                rib_ptr->StatusVec()[neig_side] = fStatus[iside];
                rib_ptr->Optimize();
                break;
            }
            case 2:{
                // check if DFNFace exists
                DFNFace *face_ptr = fFracture->Face(gel->Index());
                if(!face_ptr){
                    DFNFace neig_face(gel,fFracture);
                    fFracture->AddFace(neig_face);
                    face_ptr = fFracture->Face(gel->Index());
                }
                // skip if StatusVec entry already match
                if(fStatus[iside] == face_ptr->StatusVec()[neig_side]){continue;}
                // else, match StatusVec entry for iside and optimize
                face_ptr->StatusVec()[neig_side] = fStatus[iside];
                face_ptr->Optimize();
                break;
            }
            case 3: {continue;}
        }
    }
    return;
}






bool DFNRib::UpdateMaterial(){
    if(fGeoEl->MaterialId() == DFNMaterial::Efracture) {return false;}
    if(fStatus[0] && fStatus[1]){
        fGeoEl->SetMaterialId(DFNMaterial::Efracture);
        return true;
    }
    return false;
}




inline REAL VectorNorm(TPZManVector<REAL,3> &vector){
    int n = vector.NElements();
    REAL temp = 0.0;
    for (int i = 0; i < n; i++){
        temp += vector[i]*vector[i];
    }
    return sqrt(temp);
}

inline REAL Distance(TPZManVector<REAL,3> &vector1, TPZManVector<REAL,3> &vector2){
    int n = vector1.NElements();
    TPZManVector<REAL,3> difference(n,0);
    for(int i = 0; i<n;i++){
        difference[i] = vector1[i] - vector2[i];
    }
    return VectorNorm(difference);
}