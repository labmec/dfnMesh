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
DFNRib::DFNRib(TPZGeoEl *gel, DFNFracture *Fracture, int status) :
    fGeoEl(gel),
    fStatus(status),
    fFracture(Fracture),
    fIntersectionIndex(-1)
{fCoord.resize(0);};

// Copy constructor
DFNRib::DFNRib(const DFNRib &copy){
    this->operator=(copy);
}

// Assignment operator
DFNRib &DFNRib::operator=(const DFNRib &copy){
    fGeoEl = copy.fGeoEl;
    fStatus = copy.fStatus;
    fCoord = copy.fCoord;
    fIntersectionIndex = copy.fIntersectionIndex;
    fFracture = copy.fFracture;
    return *this;
}

inline REAL VectorNorm(TPZManVector<REAL,3> &vector);
inline REAL Distance(TPZManVector<REAL,3> &vector1, TPZManVector<REAL,3> &vector2);


/// Get real intersection coordinates
TPZManVector<REAL, 3> DFNRib::RealCoord(){
    TPZManVector<REAL, 3> coord(3, 0);
    TPZGeoMesh *gmesh = fGeoEl->Mesh();

    if(fStatus == 2 && fIntersectionIndex < 0)
        {coord = this->fCoord;}
    else
        {gmesh->NodeVec()[fIntersectionIndex].GetCoordinates(coord);}
    
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
    // @TODO discuss with your advisor if this is really necessary
    // where is the refinement pattern used?
    this->CreateRefPattern();

    // set children
        TPZManVector<int64_t,2> cornerindices(2,0);
        TPZGeoEl *child;
        int64_t elindex = gmesh->NElements();
        for(int i=0; i<2; i++){
            cornerindices[i%2] = fGeoEl->NodeIndex(i%2);
            cornerindices[(i+1)%2] = fIntersectionIndex;
            child = gmesh->CreateGeoElement(EOned, cornerindices, fGeoEl->MaterialId(), elindex);
            fGeoEl->SetSubElement(i,child);
            child->SetFatherIndex(fGeoEl->Index());
            elindex++;
        }
}







void DFNRib::CreateRefPattern(){
        // refinement mesh
        TPZGeoMesh refPatternMesh;
        int refnnodes = fGeoEl->NCornerNodes() + this->NeedsRefinement();
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





void DFNRib::SnapIntersection_force(){
    REAL BIG_NUMBER = 1.e12;
    this->SnapIntersection_try(BIG_NUMBER);
}



bool DFNRib::SnapIntersection_try(REAL tolDist){
    // If there's an intersection at the 1D side of this DFNRib, check if 
    // that intersection violates the tolerable distance and, if necessary, snap it down to a lower-dimensional side
    int64_t closestnode = -1;
    if(this->NeedsSnap(closestnode,tolDist)){
        fIntersectionIndex = fGeoEl->NodeIndex(closestnode);
        fStatus = closestnode;
        UpdateNeighbours(closestnode);
        return true;
    }
    return false;
}

bool DFNRib::NeedsSnap(int64_t& closestnode, REAL tolDist){
    if(!this->NeedsRefinement()) return false;
    
    TPZManVector<REAL,3> coord(3,0);
    fGeoEl->NodePtr(0)->GetCoordinates(coord);
    REAL dist0 = Distance(fCoord,coord);

    fGeoEl->NodePtr(1)->GetCoordinates(coord);
    REAL dist1 = Distance(fCoord,coord);

    closestnode = (dist0 < dist1 ? 0 : 1);
    REAL dist = MIN(dist0,dist1);

    if(dist<tolDist)
        return true;
    else 
        return false;
}



void DFNRib::UpdateNeighbours(int iside){

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
                    break;
                    // DFNRib neig_rib(gel,fFracture);
                    // fFracture->AddRib(neig_rib);
                    // rib_ptr = fFracture->Rib(gel->Index());
                }
                // // skip if StatusVec entry already match
                // if(fStatus[iside] == rib_ptr->StatusVec()[neig_side]){continue;}
                // // else, match StatusVec entry for iside and optimize
                // rib_ptr->StatusVec()[neig_side] = fStatus[iside];
                // rib_ptr->SnapIntersection_try();
                // if(rib_ptr->NeedsRefinement()) {rib_ptr->SnapIntersection_try();}
                if(rib_ptr->NeedsRefinement()) {rib_ptr->SnapIntersection_force();}
                break;
            }
            case 2:{
                // check if DFNFace exists
                DFNFace *neig_face = fFracture->Face(gel->Index());
                if(!neig_face) break;
                if(!neig_face->UpdateStatusVec()) break;
		        neig_face->UpdateRefMesh();
                break;
            }
            case 3: {continue;}
        }
    }
    return;
}





// // Might be useful at some point... but not yet
// bool DFNRib::UpdateMaterial(){
//     DebugStop();
//     if(fGeoEl->MaterialId() == DFNMaterial::Efracture) {return false;}
//     if(fStatus[0] && fStatus[1]){
//         fGeoEl->SetMaterialId(DFNMaterial::Efracture);
//         return true;
//     }
//     return false;
// }




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



void DFNRib::Print(std::ostream &out) const
{
    out<<"\nGeometric Element # "<<fGeoEl->Index();
	out<<"\nSide intersected:"<<fStatus;
	out<<"\nIntersection Coord : "<< fCoord;
	out<<"\nIntersection Node index: " << fIntersectionIndex << "\n";
}