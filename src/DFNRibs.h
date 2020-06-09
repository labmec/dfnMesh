/*! 
 *	DFNRibs.hpp
 *  @authors   Pedro Lima
 *  @authors   Jorge Ordoñez
 *  @date      2018-2019
 */

#ifndef DFNRibs_h
#define DFNRibs_h

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

/// DFNRibs class describes ribs and whether it's cut by a plane. It also carries a method to split itself into two ribs at a given point.
class DFNRibs
{
private:
    
    /// The index of the one dimensional element
    int64_t fRibIndex;

    /// Indicates whether it intersects the plane
    bool fIsCut;
    
    /// Indices of the subelements
    // TPZManVector<int64_t,2> fSubElements;
    
    // /// Intersection point index
    // int64_t fIntersectionIndex = -1;
    /// Index of geometric element EPoint at intersection
    int64_t fIntersectionIndex = -1;

    /// Father Element
    int64_t fFather = -1;

public:
    /// Empty constructor
    DFNRibs();
    
    /// Rib defined by the element index, indicating whether the
    /// intersection point is within the fracture surface
    DFNRibs(int64_t index, bool inplane);
    
    /// Copy constructor
    DFNRibs(const DFNRibs &copy);
    
    /// Assignment operator
    DFNRibs &operator=(const DFNRibs &copy);
    
    /**
     * @brief Set the index for an element if the fracture plane is cut by this element
     * @param Element index
     * @param Yes or no if the element is cutting the fracture plane
     */
    void SetElementIndex(int64_t elindex, bool IsCut);
    
    /// Element index
    int64_t ElementIndex() const;
    
    /// Intersects the plane or not
    void SetIsCut(bool is) ;
    
    /// Intersects the plane or not
    bool IsCut() const;
    
    /// Return the subelement indices
    TPZVec<int64_t> SubElements() const;
    
    /// Returns intersection EPoint element index
    int64_t IntersectionIndex() const {return fIntersectionIndex;}

    /// Set the subelement indices
    void SetChildren(const TPZVec<int64_t> &subels);
    
    ///Divide the given rib and generate the subelements
    void DivideRib(TPZGeoMesh *gmesh, TPZVec<REAL> interpoint,int matID);

};

#endif /* DFNRibs_h */
