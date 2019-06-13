//
//  FaceFrac.hpp
//  reservoirlib
//
//  Created by Jorge Paúl Ordóñez Andrade on 29/10/2018.
//

#ifndef TRSFace_hpp
#define TRSFace_hpp

#include <stdio.h>
#include "pzmanvector.h"
#include "pzmatrix.h"
#include "pzcompel.h"
#include "pzgeoelbc.h"

/// TRSFace class describes a plane and whether it's cut by a plane. It also carries a method to split itself into smaller sub faces.
class TRSFace
{
    
private:
    /// Face element index whithin it's gmesh file
    int64_t fFaceIndex;

    /// Indicates whether it intersects a plane
    bool fCutisInplane;

//  TPZManVector<int64_t,2> fRibIndexes;
    
    /// Vector of it's ribs
    TPZManVector<int64_t,2> fRibs;

    /// Vector of sub faces
    TPZManVector<int64_t,4> fSubFaces;
   
    
public:
    
    /// Empty constructor
    TRSFace();
    
    /// Constructor
    TRSFace(int64_t index, bool inplane);
    
    /// Copy constructor
    TRSFace(const TRSFace &copy);
    
    /// Assignment operator
    TRSFace &operator=(const TRSFace &copy);
    
    /// Define the element index and whether it cuts the plane
    void SetElementIndex(int64_t elindex, bool cutsplane);
    
    /// Element index
    int64_t ElementIndex() const;
    
    /// Intersects the plane or not
    bool CutsPlane() const;
    
    /// Return the subelement indices
    TPZVec<int64_t> SubElements() const;
    
    /// Set the subelement indices
    void DefineRibDivide(const TPZVec<int64_t> &subels);
    
    /// Divide the given surface and generate subelements
    void DivideSurface(TPZGeoMesh *gmesh, int matid);
    
    /// Set ribs in the surface
    void SetRibsInSurface(TPZManVector<int64_t,2> ribsinsurface);
    TPZManVector<int64_t,2> RibsInSurface();
    
};
#endif /* TRSFace_h */
