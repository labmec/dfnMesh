/*! 
 *	TRSRibs.cpp
 *  @authors   Jorge Ordoñez
 *  @authors   Pedro Lima
 *  @date      2018-2019
 */

#include "TRSRibs.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzbiharmonic.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"
#include "TPZRefPatternDataBase.h"

// Constructor
TRSRibs::TRSRibs(){
    fIsCut=false;
}

// Constructor using an index and a asking if the plane is cut or not
TRSRibs::TRSRibs(int64_t index, bool IsCut){
    fRibIndex=index;
    fIsCut=IsCut;
}

// Copy constructor
TRSRibs::TRSRibs(const TRSRibs &copy){
    fRibIndex = copy.fRibIndex;
    fIsCut = copy.fIsCut;
    fSubElements = copy.fSubElements;
    fIntersectionIndex = copy.fIntersectionIndex;
    fFather = copy.fFather;
}

// Assignment operator
TRSRibs &TRSRibs::operator=(const TRSRibs &copy){
    fRibIndex = copy.fRibIndex;
    fIsCut = copy.fIsCut;
    fSubElements = copy.fSubElements;
    fIntersectionIndex = copy.fIntersectionIndex;
    fFather = copy.fFather;
    return *this;
}

/**
 * @brief Set the index for an element if the fracture plane is cut by this element
 * @param Element index
 * @param Yes or no if the element is cutting the fracture plane
 */

void TRSRibs::SetElementIndex(int64_t elindex, bool IsCut){
    fRibIndex=elindex;
    fIsCut=IsCut;
}

/**
* @return An index (geomesh associated) elements vector
*/

int64_t TRSRibs::ElementIndex() const{
    return fRibIndex;
    
}

/**
 * @brief Set if the cut is within the plane
 * @param True or False if the plane is cut
 */

void TRSRibs::SetIsCut(bool IsCut) {
    fIsCut = IsCut;
}
/**
 * @return If the intersection point is within the fracture plane
 */

bool TRSRibs::IsCut() const{
    return fIsCut;
}

/**
 * @return An index (geomesh associated) subelements vector
 */

TPZVec<int64_t> TRSRibs::SubElements() const{
    return fSubElements;
}

/**
 * @brief Define the divided rib
 * @param Subelements vector
 */

void TRSRibs::SetChildren(const TPZVec<int64_t> &subels){
    fSubElements = subels;
}

/**
 * @brief Divide a 1D element into 2 1D elements (for now)
 * @param Mesh
 * @param Intersection point between the fracture plane and the rib
 * @param Material id to assign for the new ribs
 */

void TRSRibs::DivideRib(TPZGeoMesh *gmesh,TPZVec<REAL> IntersectionCoord, int matID){
    int iel_index = fRibIndex;          //Element index
    if(!gmesh->Element(iel_index)){     // If the element does not exist the code is going to break
        std::cout<<"No gel associated to the Rib\n";
        DebugStop();
    }
    TPZGeoEl *gel = gmesh->Element(iel_index);
    int64_t nelements = gmesh->NElements();
    int64_t nnodes = gmesh->NNodes();
    // gmesh->NodeVec().Resize(nnodes+1);       //Adding an extra node
    // gmesh->NodeVec()[nnodes].Initialize(IntersectionCoord, *gmesh);
    
    // // create GeoEl for intersection point
    // TPZVec<int64_t> vnnodes(1,nnodes);
    // gmesh->CreateGeoElement(EPoint,vnnodes,45,nelements);     // can't get 1D neighbours from nodes... so had to create GeoEl EPoint for all intersection points and track GeoEl index instead of nodeindex
    // fIntersectionIndex = nelements;

    // // Create children ribs
    // // TPZManVector<TPZGeoEl *, 2> subelements(2);
	// fSubElements.Resize(2);
    // TPZVec<int64_t> cornerindexes(2);
    // TPZManVector<TPZGeoEl *, 2> children(2); 
    // //child 1
    //     cornerindexes[0] = gel->NodeIndex(0);    //Setting the new rib node as node 0
    //     cornerindexes[1] = nnodes;   //Setting the new rib node as node 1
    //     nelements++;
    //     children[0] = gmesh->CreateGeoElement(EOned, cornerindexes, matID, nelements);
    // // child 2
    //     cornerindexes[0] = gel->NodeIndex(1);    //Setting the new rib node as node 0
    //     cornerindexes[1] = nnodes;   //Setting the new rib node as node 1
    //     nelements++;
    //     children[1] = gmesh->CreateGeoElement(EOned, cornerindexes, matID, nelements);
	// // SetChildren(SubElements);
    //     nnodes = gmesh->NNodes();
    //     gel->Divide(children);
    //     nnodes = gmesh->NNodes();
    // fSubElements[0] = nelements-1;              //Setting a new subelement
    // fSubElements[1] = nelements;              //Setting a second new subelement

    // TPZAutoPointer<TPZRefPattern> ptr = gel->GetRefPattern();
    // gel->SetRefPattern(ptr);

    TPZManVector<TPZGeoEl *,2> children(2);
    // if(!gRefDBase.GetUniformRefPattern(EOned))
    // {
    //     gRefDBase.InitializeUniformRefPattern(EOned);
    // }
    gel->Divide(children);
    // gmesh->NodeVec()[nnodes].Initialize(IntersectionCoord, *gmesh);
    gmesh->NodeVec()[nnodes].SetCoord(IntersectionCoord);


    // create GeoEl for intersection point
    TPZVec<int64_t> vnnodes(1,nnodes);
    gmesh->CreateGeoElement(EPoint,vnnodes,45,nelements);     
    fIntersectionIndex = nelements;

    fSubElements.resize(2);
    fSubElements[0] = children[0]->Index();
    fSubElements[1] = children[1]->Index();

}



