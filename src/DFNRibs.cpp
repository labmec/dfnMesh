/*! 
 *	DFNRibs.cpp
 *  @authors   Pedro Lima
 *  @authors   Jorge Ordoñez
 *  @date      2018-2019
 */

#include "DFNRibs.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzbiharmonic.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"
#include "TPZRefPatternDataBase.h"

// Constructor
DFNRibs::DFNRibs(){
    fIsCut=false;
}

// Constructor using an index and a asking if the plane is cut or not
DFNRibs::DFNRibs(int64_t index, bool IsCut){
    fRibIndex=index;
    fIsCut=IsCut;
}

// Copy constructor
DFNRibs::DFNRibs(const DFNRibs &copy){
    fRibIndex = copy.fRibIndex;
    fIsCut = copy.fIsCut;
    fSubElements = copy.fSubElements;
    fIntersectionIndex = copy.fIntersectionIndex;
    fFather = copy.fFather;
}

// Assignment operator
DFNRibs &DFNRibs::operator=(const DFNRibs &copy){
    fRibIndex = copy.fRibIndex;
    fIsCut = copy.fIsCut;
    fSubElements = copy.fSubElements;
    fIntersectionIndex = copy.fIntersectionIndex;
    fFather = copy.fFather;
    return *this;
}

/**
 * @brief Set the index for an element if the fracture plane is cut by this element
 * @param elindex Element index
 * @param IsCut: flag if the element is cut by the fracture plane
 */

void DFNRibs::SetElementIndex(int64_t elindex, bool IsCut){
    fRibIndex=elindex;
    fIsCut=IsCut;
}

/**
* @return An index (geomesh associated) elements vector
*/

int64_t DFNRibs::ElementIndex() const{
    return fRibIndex;
    
}

/**
 * @brief Set if the cut is within the plane
 * @param True or False if the plane is cut
 */

void DFNRibs::SetIsCut(bool IsCut) {
    fIsCut = IsCut;
}
/**
 * @return If the intersection point is within the fracture plane
 */

bool DFNRibs::IsCut() const{
    return fIsCut;
}

/**
 * @return An index (geomesh associated) subelements vector
 */

TPZVec<int64_t> DFNRibs::SubElements() const{
    return fSubElements;
}

/**
 * @brief Define the divided rib
 * @param Subelements vector
 */

void DFNRibs::SetChildren(const TPZVec<int64_t> &subels){
    fSubElements = subels;
}

/**
 * @brief Divide a 1D element into 2 1D elements (for now)
 * @param Mesh
 * @param Intersection point between the fracture plane and the rib
 * @param Material id to assign for the new ribs
 */

void DFNRibs::DivideRib(TPZGeoMesh *gmesh,TPZVec<REAL> IntersectionCoord, int matID){
    int iel_index = fRibIndex;          //Element index
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
                gel->ComputeXInverse(IntersectionCoord,qsi,1E-5);

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
            if(fIntersectionIndex == 2) refPatternMesh.NodeVec()[2].Initialize(IntersectionCoord,refPatternMesh);
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
            gmesh->NodeVec()[nnodes].Initialize(IntersectionCoord,*gmesh);
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
    
}



