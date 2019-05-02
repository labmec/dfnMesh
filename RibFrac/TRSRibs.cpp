//
// RibFrac.cpp
// RibFrac
//
//  Created by JORGE ORDOÑEZ on 22/5/18.
//  Copyright © 2018 JORGE ORDOÑEZ. All rights reserved.
//

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
    fCutisInplane=false;
}

// Constructor using an index and a asking if the plane is cut or not
TRSRibs::TRSRibs(int64_t index, bool inplane){
    fRibIndex=index;
    fCutisInplane=inplane;
}

// Copy constructor
TRSRibs::TRSRibs(const TRSRibs &copy){
    fRibIndex=copy.fRibIndex;
    fCutisInplane=copy.fCutisInplane;
    fSubElements = copy.fSubElements;
}

// Assignment operator
TRSRibs &TRSRibs::operator=(const TRSRibs &copy){
    fRibIndex=copy.fRibIndex;
    fCutisInplane=copy.fCutisInplane;
    if(copy.fSubElements.size()==2){
    fSubElements = copy.fSubElements;
    }
    return *this;
}

/**
 * @brief Set the index for an element if the fracture plane is cut by this element
 * @param Element index
 * @param Yes or no if the element is cutting the fracture plane
 */

void TRSRibs::SetElementIndex(int64_t elindex, bool cutsplane){
    fRibIndex=elindex;
    fCutisInplane=cutsplane;
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

void TRSRibs::SetCutsPlane(bool is) {
    fCutisInplane=is;
}
/**
 * @return If the intersection point is within the fracture plane
 */

bool TRSRibs::CutsPlane() const{
    return fCutisInplane;
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

void TRSRibs::DefineRibDivide(const TPZVec<int64_t> &subels){
    fSubElements = subels;
}

/**
 * @brief Divide a 1D element into 2 1D elements (for now)
 * @param Mesh
 * @param Intersection point between the fracture plane and the rib
 * @param Material id to assign for the new ribs
 */

void TRSRibs::DivideRib(TPZGeoMesh *gmesh,TPZVec<REAL> interpoint, int matid){
    int iel_index = fRibIndex;          //Element index
    if(!gmesh->Element(iel_index)){     // If the element does not                          exist the code is going to break
        std::cout<<"No gel asociado a Rib"<<std::endl;
        DebugStop();
    }
    TPZGeoEl *gel = gmesh->Element(iel_index);
    int64_t nelements = gmesh->NElements();
    
    int nNods= gmesh->NNodes();
    
    TPZVec<TPZGeoNode> Node(1);
    Node[0].SetCoord(interpoint);
    gmesh->NodeVec().Resize(nNods+1);       //Adding an extra node
    gmesh->NodeVec()[nNods]=Node[0];
   
    TPZVec<int64_t> cornerindexes(2);
    cornerindexes[0] = gel->NodeIndex(0);   //Setting the new rib node as node 0
    cornerindexes[1] = nNods;               //Setting the new rib node as node 1
    gmesh->CreateGeoElement(EOned, cornerindexes, matid, nelements);
    fSubElements.Resize(2);
    fSubElements[0]=nelements;              //Setting a new subelement
    
    cornerindexes[0] = nNods;              //Setting a second new rib node as node 0
    cornerindexes[1] = gel->NodeIndex(1);  //Setting a second new rib node as node 1
    nelements++;
     gmesh->CreateGeoElement(EOned, cornerindexes, matid, nelements);
    fSubElements[1]=nelements;              //Setting a second new subelement
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
