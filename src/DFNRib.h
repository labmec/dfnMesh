/*! 
 *	DFNRib.hpp
 *  @authors   Pedro Lima
 *  @authors   Jorge Ordoñez
 *  @date      2018-2019
 */

#ifndef DFNRib_h
#define DFNRib_h

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
class DFNFracture;

/// DFNRib class describes ribs and whether it's cut by a plane. It also carries a method to split itself into two ribs at a given point.
class DFNRib
{
private:
    
    /// pointer to its original geometric element
	TPZGeoEl *fGeoEl = nullptr;

	/// Contains the side of this element that has been intersected
	int fStatus = -1;

    /// Anticipated coordinates of an intersection point initialy found for this DFNRib and a DFNPolygon. 
    /// The point may be moved afterwards during DFNRib::SnapIntersection()
	TPZManVector<REAL, 3> fCoord;
	
	/// Index of intersection node if any correspond to it in the geometric mesh
	int64_t fIntersectionIndex = -1;

    /// Pointer to the fracture
    DFNFracture *fFracture = nullptr;

    /// Offbound ribs are edges of an intersected polyhedral volume, and intersected by the plane that contains a DFNPolygon
    /// however, their intersection node is not within DFNPolygon's bounds
    bool fOffbound = false;
    // bool fAdjacent = false;
    // bool fAttached = false;
    // bool fVirtual = false;
    // bool fExternal = false;
    // bool fOutter = false;
    // bool fProwling = false;

public:
    /// Empty constructor
    DFNRib(): fGeoEl(nullptr), fStatus(-1){fCoord.resize(0);};
    
    /// Default constructor takes a pointer to geometric element
    DFNRib(TPZGeoEl *gel, DFNFracture *Fracture, int status = -1);
    
    /// Copy constructor
    DFNRib(const DFNRib &copy);
    
    /// Assignment operator
    DFNRib &operator=(const DFNRib &copy);
    
    /// Destructor
    ~DFNRib(){};

    /// Index of the associated geometric element
    int64_t Index() const {return fGeoEl->Index();}

    /// Element pointer
    TPZGeoEl* GeoEl() const {return fGeoEl;}
    
    /// Set index of intersection node
    void SetIntersectionIndex(int64_t index){fIntersectionIndex = index;}
    
    /// Get index of intersection node
    int64_t IntersectionIndex() const {return fIntersectionIndex;}
    
    /// Set anticipated intersection coordinates (real coordinates are defined after DFNRib::SnapIntersection_try())
    void SetIntersectionCoord(TPZManVector<REAL, 3> coord){
        fCoord = coord;
    }

    /// Get anticipated intersection coordinates
    TPZManVector<REAL, 3> AntCoord(){return fCoord;}

    /// Get real intersection coordinates (after SnapIntersection)
    TPZManVector<REAL, 3> RealCoord();

    /// Give face a pointer to which fracture is cutting it
    void SetFracture(DFNFracture *Fracture){fFracture = Fracture;}

    /**
     * @brief Flag if rib was found to be intersected at any side
    */
    inline bool IsIntersected() const{
        return fStatus != -1;
    }

    /// Check if this rib should be refined
    inline bool NeedsRefinement() const {return (fStatus == 2);}

    /// Set Status (the side where intersection happens)
    void SetStatus(const int status){
        if(status < 0 || status > 2) DebugStop();
        fStatus = status;
    }

    /// Get side of this rib that has been intersected
    const int Status() const {return fStatus;}

    /// Check if should be refined and generate the subelements
    void Refine();

    /// Sets material id
    void SetMaterialId(int matID){
        fGeoEl->SetMaterialId(matID);
    };

    /// After intersection snap, recursively snap intersections at neighbours through side iside
    void UpdateNeighbours(int iside);

    /// Check if rib will be incorporated onto fracture and correct its material id as such.
    /// Returns true if changes have been made
    bool UpdateMaterial();
    
    /**
     * @brief Check geometry of intersection against a tolerance, snaps intersection 
     * to closest side(s) if necessary and recursevily modifies affected neighbours.
     * @return True if any snap happened.
     * @param tolDist: Minimum acceptable length
    */
    bool SnapIntersection_try(REAL tolDist = 1e-4);

    /// Forces projection (snap) of the intersection node to the closest lower dimensional side (nothing happens if insersection is already at zero dim side)
    void SnapIntersection_force();

    /**
     * @brief Check geometry of intersection against tolerance, to test if intersection should be snapped
     * @return False if no tolerance is violated or if intersection has already been snapped.
     * @param tolDist: Minimum acceptable distance
    */
    bool NeedsSnap(int64_t& closestnode, REAL tolDist = 1e-4);

    /** @brief Print method for logging */
    void Print(std::ostream& out = std::cout) const;

    private:
        /// Creates refinement pattern based on status vector and intersection node index
        void CreateRefPattern();

};

#endif /* DFNRib_h */
