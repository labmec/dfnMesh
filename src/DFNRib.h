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

	/// Contains the side of this element that has been intersected. By default, a fracture intersects a rib at the 1D side, but the intersection may be snap to one of the nodes.
	int fIntersectionSide = 2;

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
    DFNRib(): fGeoEl(nullptr), fIntersectionSide(-1){fCoord.resize(0);};
    
    /// Default constructor takes a pointer to geometric element
    DFNRib(TPZGeoEl *gel, DFNFracture *Fracture);
    
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
    
    /// Sets the GeoEl
    void SetGeoEl(TPZGeoEl* geoEl) {fGeoEl = geoEl;}
    
    /// Set index of intersection node
    void SetIntersectionIndex(int64_t index){fIntersectionIndex = index;}
    
    /// Get index of intersection node
    int64_t IntersectionIndex() const {return fIntersectionIndex;}
    
    /// Set anticipated intersection coordinates (real coordinates are defined after DFNRib::SnapIntersection_try())
    void SetIntersectionCoord(TPZManVector<REAL, 3>& coord){fCoord = coord;}

    /// Set flag for offbound rib
    void FlagOffbound(bool flag){fOffbound = flag;}

    /// Get anticipated intersection coordinates
    const TPZManVector<REAL, 3>& AntCoord(){return fCoord;}

    /// Get real intersection coordinates (after SnapIntersection)
    TPZManVector<REAL, 3> RealCoord() const;

    /// Give face a pointer to which fracture is cutting it
    void SetFracture(DFNFracture *Fracture){fFracture = Fracture;}

    /**
     * @brief Flag if rib was found to be intersected at any side
    */
    inline bool IsIntersected() const{
        return fIntersectionSide != -1;
    }

    /// Check if intersection node could be snapped (to actually check if snap will happen given a tolerance, use DFNRib::NeedsSnap)
    inline bool CanBeSnapped() const {return (fIntersectionSide == 2);}
    /// Check if this rib should be refined
    inline bool NeedsRefinement() const {return (fIntersectionSide == 2);}

    /// Set IntersectionSide (the side where intersection happens)
    void SetIntersectionSide(const int sideindex){
        if(sideindex < 0 || sideindex > 2) DebugStop();
        fIntersectionSide = sideindex;
    }

    /// Get side of this rib that has been intersected
    int IntersectionSide() const {return fIntersectionSide;}

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

    /// Forces projection (snap) of the intersection node to the closest lower dimensional side (nothing happens if insersection is already at zero dim side).
    /// @param closestnode If filled, snapped will be forced to node of local index == closestnode. If left -1, code will compute closest node.
    void SnapIntersection_force(int closestnode = -1);

    /**
     * @brief Check geometry of intersection against tolerance, to test if intersection should be snapped
     * @return False if no tolerance is violated or if intersection has already been snapped.
     * @param tolDist: Minimum acceptable distance
    */
    bool NeedsSnap(int& closestnode, REAL tolDist = 1e-4);

    /** @brief Adds a pointer of this rib into the corresponding position of its neighbour faces ribvectors*/
    void AppendToNeighbourFaces();

    bool IsOffbound() const{return fOffbound;}

    /** @brief Print method for logging */
    void Print(std::ostream& out = std::cout) const;

    private:
        /// Creates refinement pattern based on how a fracture intersects this rib
        void CreateRefPattern();

};

inline std::ostream& operator<<(std::ostream &out, const DFNRib& rib){
    rib.Print(out);
    return out;
}

#endif /* DFNRib_h */
