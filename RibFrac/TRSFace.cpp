//
//  FaceFrac.cpp
//  reservoirlib
//
//  Created by Jorge Paúl Ordóñez Andrade on 29/10/2018.
//

#include "TRSFace.h"

// Empty constructor
TRSFace::TRSFace(){
    fFaceIndex = -1;
    fCutisInplane = false;
}

//Constructor
TRSFace::TRSFace(int64_t index, bool inplane){
    fFaceIndex = index;
    fCutisInplane = inplane;
}

/// Copy constructor
TRSFace::TRSFace(const TRSFace &copy){
    fFaceIndex = copy.fFaceIndex;
    fCutisInplane = copy.fCutisInplane;
    fFaceIndex = copy.fFaceIndex;
    fSubFaces = copy.fSubFaces;
}

/// Assignment operator
TRSFace &TRSFace::operator=(const TRSFace &copy){
    fFaceIndex = copy.fFaceIndex;
    fCutisInplane = copy.fCutisInplane;
    fFaceIndex = copy.fFaceIndex;
    fSubFaces = copy.fSubFaces;
    return *this;
}

/// Define the element index and whether it cuts the plane
void TRSFace::SetElementIndex(int64_t elindex, bool cutsplane){
    fFaceIndex = elindex;
    fCutisInplane = cutsplane;
}

/// Element index
int64_t TRSFace::ElementIndex() const{
    return fFaceIndex;
}

/// Intersects the plane or not
bool TRSFace::CutsPlane() const{
    return fCutisInplane;
}

/// Return the subelement indices
TPZVec<int64_t> TRSFace::SubElements() const{
    return fSubFaces;
}

/// Set the subelement indices
void TRSFace::DefineRibDivide(const TPZVec<int64_t> &subels){
    fSubFaces = subels;
}

///
///Divide the given rib and generate the subelements
void TRSFace::DivideSurface(TPZGeoMesh *gmesh, int matid){
    //@TODO Check consistence
//    if(!(gmesh->Element(fFaceIndex))){DebugStop();}
//    TPZGeoEl *gel = gmesh->Element(fFaceIndex);
//    int no_intr_P1 = fRibInter[0].second;
//    int no_intr_P2 = fRibInter[1].second;
//    int P1 = fRibInter[0].first;
//    int P2 = fRibInter[1].first;
//    TPZGeoEl *gel1 = gmesh->Element(P1);
//    TPZGeoEl *gel2 = gmesh->Element(P2);
//    TPZVec<int> status(4,0);
//    for(int iside=4; iside<8; iside++){
//        TPZGeoElSide gelside(gel,iside);
//        TPZGeoElSide neig = gelside.Neighbour();
//        int index_neig = neig.Element()->Index();
//        if(index_neig==P1){
//            status[iside]=1;
//        }
//        if(index_neig==P2){
//            status[iside]=2;
//        }
//    }
//    int cases=-1;
//    if(status[0]!=0 && status[2]!=0){
//        cases=0;
//    }
//    if(status[1]!=0 && status[3]!=0){
//         cases=1;
//    }
//    if (cases < 0){
//        
//    }
    
  
    
    
}

void TRSFace::SetRibsInSurface(TPZManVector<int64_t,2> ribsinsurface){
    fRibs=ribsinsurface;
}
TPZManVector<int64_t,2> TRSFace::RibsInSurface(){
    return fRibs;
}





//bool TRSFace::DivideSurface(TPZGeoMesh *gmesh){
//   // si una cara tiene dos ribs cortados en la misma cara entonces debe ser dividido
//
//}
