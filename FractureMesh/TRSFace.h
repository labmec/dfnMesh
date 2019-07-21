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
    bool fIsCut;

    /// Indicates which ribs are cut
    TPZManVector<bool,4> fRibStatus;
    
    /// Vector of its ribs
    TPZManVector<int64_t,4> fRibs;

    /// Vector of sub faces
    TPZManVector<int64_t,5> fSubFaces;
    
    /// Index of intersection point (only for FractureMesh::EndFaces)
    int64_t fIntersection = -1;
    
public:
    
    /// Empty constructor
    TRSFace();
    
    /// Constructor
    TRSFace(int64_t index, bool iscut);
    
    /// Copy constructor
    TRSFace(const TRSFace &copy);
    
    /// Assignment operator
    TRSFace &operator=(const TRSFace &copy);
    
    /// Define the element index and whether it cuts the plane
    void SetElementIndex(int64_t elindex, bool iscut);
    
    /// Element index
    int64_t ElementIndex() const;
    
    /// Intersects the plane or not
    bool IsCut() const;
    
    /// Return the subelement indices
    TPZVec<int64_t> SubElements() const;
    
    /// Set the subelement indices
    void DefineRibDivide(const TPZVec<int64_t> &subels);
    
    /// Divide the given surface and generate subelements
    void DivideSurface(TPZGeoMesh *gmesh, int matid);
    
    /// Set ribs cut
    void SetRibsCut(TPZManVector<int64_t,2> ribsinsurface);
    TPZManVector<int64_t,2> RibsInSurface();
    
    /// Access intersection index data (only for FractureMesh::EndFaces)
    void SetIntersectionIndex(int64_t index) {fIntersection = index;}

    /// Returns index for intersection point in face (only for FractureMesh::EndFaces)
    int64_t IntersectionIndex() {return fIntersection;}


};
#endif /* TRSFace_h */
