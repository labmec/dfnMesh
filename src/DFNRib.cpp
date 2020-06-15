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
    fIntersectionIndex = copy.fIntersectionIndex;
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

                fIntersectionIndex = 2;
                if(qsi[0] < -0.95 || qsi[0] > 0.95){
                    refnnodes = 2;
                    if(qsi[0] < -0.95) fIntersectionIndex = 0;
                    else fIntersectionIndex = 1; 
                }
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
                    
            refPatternMesh.BuildConnectivity(); 
        // define refPattern
        refPatternMesh.BuildConnectivity(); 
        TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(refPatternMesh);
        fGeoEl->SetRefPattern(refpat);
    }


    // set new node in gmesh 
    int64_t nnodes = gmesh->NNodes();
    MElementType etype;
    switch(fIntersectionIndex){
        case 0:{
            fIntersectionIndex = gel->NodeIndex(0); 
            etype = EPoint;
            break;}
        case 1:{
            fIntersectionIndex = gel->NodeIndex(1);
            etype = EPoint;
            break;}
        case 2:{
            gmesh->NodeVec().Resize(nnodes+1);
            gmesh->NodeVec()[nnodes].Initialize(fCoord,*gmesh);
            fIntersectionIndex = nnodes;
            etype = EOned;
            break;
        }
        default: DebugStop();
    }

    // set children
        TPZManVector<int64_t,2> cornerindices(2);
        TPZGeoEl *child;
    // first child
        cornerindices[0] = gel->NodeIndex(0);
        cornerindices[1] = (etype == EPoint ? gel->NodeIndex(1) : fIntersectionIndex);
        int64_t elindex = gmesh->NElements();
        child = gmesh->CreateGeoElement(EOned, cornerindices, matID, elindex);
        gel->SetSubElement(0,child);
        child->SetFatherIndex(gel->Index());
    // second child
        if(etype == EOned){
            cornerindices[1] = gel->NodeIndex(1);
        }else{
            cornerindices.Resize(1);
        }
        cornerindices[0] = fIntersectionIndex;
        elindex++;
        child = gmesh->CreateGeoElement(etype, cornerindices, matID, elindex);
        gel->SetSubElement(1,child);
        child->SetFatherIndex(gel->Index());
    
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