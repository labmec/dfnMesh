//
//  RibFrac.hpp
//  RibFrac
//
//  Created by JORGE ORDOÑEZ on 22/5/18.
//  Copyright © 2018 JORGE ORDOÑEZ. All rights reserved.
//

#ifndef TRSRibs_h
#define TRSRibs_h

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include <stdio.h>
#include <tuple>
#include "pzmatrix.h"
#include "pzcompel.h"
#include "pzgeoelbc.h"
//#include "tpanic.h"

typedef TPZFMatrix<REAL> Matrix;

/// TRSRibs class describes the intersection of ribs with a plane
class TRSRibs
{
private:
    
    /// The index of the one dimensional element
    int64_t fRibIndex;

    /// Indicates whether the intersection point is in the plane
    bool fCutisInplane;
    
    /// Indices of the subelements
    TPZManVector<int64_t,2> fSubElements;
    
public:
    /// Empty constructor
    TRSRibs();
    
    /// Rib defined by the element index, indicating whether the
    /// intersection point is within the fracture surface
    TRSRibs(int64_t index, bool inplane);
    
    /// Copy constructor
    TRSRibs(const TRSRibs &copy);
    
    /// Assignment operator
    TRSRibs &operator=(const TRSRibs &copy);
    
    /// Define the element index and whether it cuts the plane
    void SetElementIndex(int64_t elindex, bool cutsplane);
    
    /// Element index
    int64_t ElementIndex() const;
    
    /// Intersects the plane or not
    void SetCutsPlane(bool is) ;
    
    /// Intersects the plane or not
    bool CutsPlane() const;
    
    /// Return the subelement indices
    TPZVec<int64_t> SubElements() const;
    
    /// Set the subelement indices
    void DefineRibDivide(const TPZVec<int64_t> &subels);
    
    ///Divide the given rib and generate the subelements
    void DivideRib(TPZGeoMesh *gmesh, TPZVec<REAL> interpoint,int matid);
};

#endif /* TRSRibs_h */
