//
// RibFrac.cpp
// RibFrac
//
//  Created by JORGE ORDOÑEZ on 22/5/18.
//  Copyright © 2018 JORGE ORDOÑEZ. All rights reserved.
//

#include "TRSRibs.h"
TRSRibs::TRSRibs(){
    fLength = -1;
}

/**
 * @brief
 * @param
 * @return
 * @return
 */

TRSRibs::TRSRibs(TPZVec<TPZVec<double>> coords){
    if(coords.size()!=2){
        std::cout<<"Se necesitan dos puntos para definir un rib";
        DebugStop();
    }
    SetRibsCoords(coords);
    calc_flength();
}

/**
 * @brief
 * @param
 * @return
 * @return
 */

TRSRibs::TRSRibs(TPZVec<TPZVec<double>> fathercoords,TPZVec<double> IntersectionPoint){
    if(fathercoords.size()!=2){
        std::cout<<"Se necesitan dos puntos para definir un rib";
        DebugStop();
    }
   //Creates Father
    SetRibsCoords(fathercoords);
    calc_flength();
//
//    std::cout<<fathercoords[0][0]<<std::endl;
//    std::cout<<fathercoords[0][1]<<std::endl;
//    std::cout<<fathercoords[0][2]<<std::endl;
    calc_child(fathercoords, IntersectionPoint);
    
    
    
    //here
    
}

/**
 * @brief
 * @param
 * @return
 * @return
 */

void TRSRibs::SetRibsCoords(TPZVec<TPZVec<double>> coords){
    if(coords.size()!=2){
        std::cout<<"Se necesitan dos puntos para definir un rib";
        DebugStop();
    }
    fcoords=coords;
}
TPZVec<TPZVec<double>> TRSRibs::GetRibsCoords(){
   
    return fcoords;
}

/**
 * @brief
 * @param
 * @return
 * @return
 */

void TRSRibs::SetIntersectionPoint( TPZVec<double> ipoint){
    if (ipoint.size()!=3){
        std::cout<<"Se necesitan 3 coordenadas para definir un punto";
        DebugStop();
    }
    fIntersectionPoint = ipoint;
}

/**
 * @brief
 * @param
 * @return
 * @return
 */

TPZVec<double> TRSRibs::GetIntersectionPoint(){
   
    return fIntersectionPoint;
}

/**
 * @brief
 * @param
 * @return
 * @return
 */

void TRSRibs::calc_child(TPZVec<TPZVec<double>> coords, TPZVec<double> ipoint){
    //@TODO VERIFICAR CONSISTENCIA DE LA INTERSECCION
    if(fQIsCut == true){
        DebugStop();
    }
    this->fQIsCut = true;
    fchild.resize(2);
    TPZVec<TPZVec<double>> neo_child;
    neo_child.resize(2);
    neo_child[0] =coords[0];
    neo_child[1] = ipoint;
    TRSRibs child1(neo_child);
    child1.calc_flength();
    fchild[0] = child1;
    neo_child[0] = ipoint;
    neo_child[1] = coords[1];
    TRSRibs child2(neo_child);
    fchild[1] = child2;
}

/**
 * @brief
 * @param
 * @return
 * @return
 */

void TRSRibs::SetNodeIndexes(TPZVec<int> indexes){
    if (indexes.size()!=3){
        std::cout<<"Se necesitan 3 coordenadas para definir un punto";
        DebugStop();
    }
    fNode_Index = indexes;
}

/**
 * @brief
 * @param
 * @return
 * @return
 */

TPZVec<int> TRSRibs::GetNodeIndexes(){
   
   return  fNode_Index;
}

/**
 * @brief
 * @param
 * @return
 * @return
 */

void TRSRibs::calc_flength(){
    double norm=0.0;

    norm += (fcoords[0][0]-fcoords[1][0])*(fcoords[0][0]-fcoords[1][0]);
    norm += (fcoords[0][1]-fcoords[1][1])*(fcoords[0][1]-fcoords[1][1]);
    norm += (fcoords[0][2]-fcoords[1][2])*(fcoords[0][2]-fcoords[1][2]);
    norm = sqrt(norm);
    fLength=norm;
}

/**
 * @brief
 * @param
 * @return
 * @return
 */

void TRSRibs::get_child(TPZVec<TRSRibs> &child){
    if (fchild.size() !=2){
        std::cout<<"No child"<<std::endl;
    }
    else{
    child = fchild;
    }
}
