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
    gmesh->NodeVec().Resize(nnodes+1);       //Adding an extra node
    gmesh->NodeVec()[nnodes].Initialize(IntersectionCoord, *gmesh);
    
    // create GeoEl for intersection point
    TPZVec<int64_t> vnnodes(1,nnodes);
    gmesh->CreateGeoElement(EPoint,vnnodes,45,nelements);     // can't get 1D neighbours from nodes... so had to create GeoEl EPoint for all intersection points and track GeoEl index instead of nodeindex
    fIntersectionIndex = nelements;

    // Create children ribs
	fSubElements.Resize(2);
    TPZVec<int64_t> cornerindexes(2);
    //child 1
    cornerindexes[0] = gel->NodeIndex(0);    //Setting the new rib node as node 0
    cornerindexes[1] = nnodes;   //Setting the new rib node as node 1
    nelements++;
    gmesh->CreateGeoElement(EOned, cornerindexes, matID, nelements);
    // child 2
    cornerindexes[0] = gel->NodeIndex(1);    //Setting the new rib node as node 0
    cornerindexes[1] = nnodes;   //Setting the new rib node as node 1
    nelements++;
    gmesh->CreateGeoElement(EOned, cornerindexes, matID, nelements);
    
	// TPZVec<int64_t> SubElements = {nelements-1,nelements};
	// SetChildren(SubElements);

    fSubElements[0] = nelements-1;              //Setting a new subelement
    fSubElements[1] = nelements;              //Setting a second new subelement
	
	// Create a TRSRib objects for children?
}









/// Information of methods used before for further coding



//
//TRSRibs::TRSRibs(TPZVec<TPZVec<double>> fathercoords,TPZVec<double> IntersectionPoint){
//    if(fathercoords.size()!=2){
//        std::cout<<"Se necesitan dos puntos para definir un rib";
//        DebugStop();
//    }
//   //Creates Father
//    SetRibsCoords(fathercoords);
//    calc_flength();
////
////    std::cout<<fathercoords[0][0]<<std::endl;
////    std::cout<<fathercoords[0][1]<<std::endl;
////    std::cout<<fathercoords[0][2]<<std::endl;
//    calc_child(fathercoords, IntersectionPoint);

//
///**
// * @brief
// * @param
// * @return
// * @return
// */
//
//void TRSRibs::SetRibsCoords(TPZVec<TPZVec<double>> coords){
//    if(coords.size()!=2){
//        std::cout<<"Se necesitan dos puntos para definir un rib";
//        DebugStop();
//    }
//    fcoords=coords;
//}
//TPZVec<TPZVec<double>> TRSRibs::GetRibsCoords(){
//
//    return fcoords;
//}
//
///**
// * @brief
// * @param
// * @return
// * @return
// */
//
//void TRSRibs::SetIntersectionPoint( TPZVec<double> ipoint){
//    if (ipoint.size()!=3){
//        std::cout<<"Se necesitan 3 coordenadas para definir un punto";
//        DebugStop();
//    }
//    fIntersectionPoint = ipoint;
//}
//
///**
// * @brief
// * @param
// * @return
// * @return
// */
//
//TPZVec<double> TRSRibs::GetIntersectionPoint(){
//
//    return fIntersectionPoint;
//}
//
///**
// * @brief
// * @param
// * @return
// * @return
// */
//
//void TRSRibs::calc_child(TPZVec<TPZVec<double>> coords, TPZVec<double> ipoint){
//    //@TODO VERIFICAR CONSISTENCIA DE LA INTERSECCION
//    if(fQIsCut == true){
//        DebugStop();
//    }
//    this->fQIsCut = true;
//    fchild.resize(2);
//    TPZVec<TPZVec<double>> neo_child;
//    neo_child.resize(2);
//    neo_child[0] =coords[0];
//    neo_child[1] = ipoint;
//    TRSRibs child1(neo_child);
//    child1.calc_flength();
//    fchild[0] = child1;
//    neo_child[0] = ipoint;
//    neo_child[1] = coords[1];
//    TRSRibs child2(neo_child);
//    fchild[1] = child2;
//}
//
///**
// * @brief
// * @param
// * @return
// * @return
// */
//
//void TRSRibs::SetNodeIndexes(TPZVec<int> indexes){
//    if (indexes.size()!=3){
//        std::cout<<"Se necesitan 3 coordenadas para definir un punto";
//        DebugStop();
//    }
//    fNode_Index = indexes;
//}
//
///**
// * @brief
// * @param
// * @return
// * @return
// */
//
//TPZVec<int> TRSRibs::GetNodeIndexes(){
//
//   return  fNode_Index;
//}
//
///**
// * @brief
// * @param
// * @return
// * @return
// */
//
//void TRSRibs::calc_flength(){
//    double norm=0.0;
//
//    norm += (fcoords[0][0]-fcoords[1][0])*(fcoords[0][0]-fcoords[1][0]);
//    norm += (fcoords[0][1]-fcoords[1][1])*(fcoords[0][1]-fcoords[1][1]);
//    norm += (fcoords[0][2]-fcoords[1][2])*(fcoords[0][2]-fcoords[1][2]);
//    norm = sqrt(norm);
//    fLength=norm;
//}
//
///**
// * @brief
// * @param
// * @return
// * @return
// */
//
//void TRSRibs::get_child(TPZVec<TRSRibs> &child){
//    if (fchild.size() !=2){
//        std::cout<<"No child"<<std::endl;
//    }
//    else{
//    child = fchild;
//    }
//}
