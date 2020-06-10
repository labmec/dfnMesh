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


/// Default constructor takes a pointer to geometric element
DFNRib::DFNRib(TPZGeoEl *gel) :
    fGeoEl(gel),
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

/// Get real intersection coordinates
TPZManVector<REAL, 3> DFNRib::RealCoord(){
    TPZManVector<REAL, 3> coord(3, 0);
    TPZGeoMesh *gmesh = fGeoEl->Mesh();

    if(this->IsCut()){
        if(fIntersectionIndex >= 0){
            gmesh->NodeVec()[fIntersectionIndex].GetCoordinates(coord);
        }else{
            coord = fCoord;
        }
    }else{
        if(fStatus[0] == 1){
            fGeoEl->NodePtr(0)->GetCoordinates(coord);
        }else{
            fGeoEl->NodePtr(1)->GetCoordinates(coord);
        }
    }
    return coord;
}

/// Divide the given rib and generate the subelements of material id matID
void DFNRib::RefineRib(int matID){
    int iel_index = this->Index();          //Element index
    TPZGeoMesh *gmesh = fGeoEl->Mesh();
    if(!gmesh->Element(iel_index)){     // If the element does not exist the code is going to break
        std::cout<<"No gel associated to the Rib\n";
        DebugStop();
    }
    TPZGeoEl *gel = gmesh->Element(iel_index);
    if(matID == -1) matID = gel->MaterialId();
    // Set refinement pattern
    {
        // dummy mesh
            TPZGeoMesh refPatternMesh;
            int refnnodes = 3;
            // check if point should be snapped to vertex
                TPZManVector<REAL,3> qsi(1,2);
                gel->ComputeXInverse(fCoord,qsi,1E-5);

                fIntersectionIndex = 2;
                if(qsi[0] < -0.95 || qsi[0] > 0.95){
                    refnnodes = 2;
                    if(qsi[0] < -0.95) fIntersectionIndex = 0;
                    else fIntersectionIndex = 1; 
                }
        // set nodes
            refPatternMesh.NodeVec().Resize(refnnodes);
            TPZManVector<REAL,3> coord(3);
            for(int i = 0; i<2; i++){
                gmesh->NodeVec()[gel->NodeIndex(i)].GetCoordinates(coord);
                refPatternMesh.NodeVec()[i].Initialize(coord,refPatternMesh);
            }
            if(fIntersectionIndex == 2) refPatternMesh.NodeVec()[2].Initialize(fCoord,refPatternMesh);
        // insert father
            TPZManVector<int64_t,2> cornerindices(2);
            cornerindices[0] = 0;
            cornerindices[1] = 1;
            int64_t elindex = 0;
            refPatternMesh.CreateGeoElement(EOned, cornerindices, matID, elindex);
        // insert children
            cornerindices[0] = 0;
            cornerindices[1] = (fIntersectionIndex == 0 ? 1 : fIntersectionIndex);
            elindex = 1;
            refPatternMesh.CreateGeoElement(EOned, cornerindices, matID, elindex);
            elindex = 2;
            if(fIntersectionIndex == 2){
                cornerindices[0] = 2;
                cornerindices[1] = 1;
                refPatternMesh.CreateGeoElement(EOned, cornerindices, matID, elindex);
            }else{
                TPZManVector<int64_t,1> pointindex(1,fIntersectionIndex);
                refPatternMesh.CreateGeoElement(EPoint, pointindex, matID, elindex);
            }
                        
            refPatternMesh.BuildConnectivity(); 
            // define refPattern
            TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(refPatternMesh);
            gel->SetRefPattern(refpat);
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


/**
 * @brief Check geometry of intersection against a tolerance, snaps intersection 
 * to closest side(s) if necessary and modifies affected neighbours.
 * @return True if any optimization has been made.
 * @param fracture: A pointer to a DFNFracture object
 * @param tol: Minimum acceptable length
*/
bool DFNRib::Optimize(DFNFracture *fracture, REAL tol){
    //@todo
    return false;
}
