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
	TPZGeoEl *fGeoEl;

	/// status vector describes how the sides of this element have been intersected
    // @TODO for a mere human reader, the word "describes" does not explain
    // this data structure (I suggest to just delete it... ;-)
	TPZManVector<int, 3> fStatus;

	/// Anticipated coordinates of intersection (may not be real coordinates)
    // @TODO for a mere human reader (not superman) the statement "may not be real
    // coordinates" is a worrysome statement
    // @TODO change the name to fIntersectionPoint (for instance)
	TPZManVector<REAL, 3> fCoord;
	
	/// Index of intersection node
    // @TODO does that mean the index of a node in the geometric mesh?
	int64_t fIntersectionIndex = -1;

    /// Pointer to the fracture
    DFNFracture *fFracture = nullptr;

public:
    /// Empty constructor
    DFNRib(): fGeoEl(nullptr), fStatus({0,0,0}){};
    
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
    
    /// Set index of intersection node
    void SetIntersectionIndex(int64_t index){fIntersectionIndex = index;}
    
    /// Get index of intersection node
    int64_t IntersectionIndex() const {return fIntersectionIndex;}
    
    /// Set anticipated intersection coordinates (real coordinates are defined by DFNRib::Optimize())
    void SetIntersectionCoord(TPZManVector<REAL, 3> coord){
        fCoord = coord;
    }

    /// Get anticipated intersection coordinates
    TPZManVector<REAL, 3> AntCoord(){return fCoord;}

    /// Get real intersection coordinates (after optimization)
    TPZManVector<REAL, 3> RealCoord();

    /// Give face a pointer to which fracture is cutting it
    void SetFracture(DFNFracture *Fracture){fFracture = Fracture;}

    /**
     * @brief Flag if rib was found to be intersected at any side
     * @return True if any element of status vector is true
     * @TODO GET RID OF FSTATUS!! Maybe status can be substituted by an enum?
     *  I dont understand how a rib can have no intersection? An object rib is created to manage
     *  the intersection!
    */
    inline bool IsIntersected() const{
        return (fStatus[0] || fStatus[1] || fStatus[2]);
    }

    /// Check if this rib should be refined
    // @TODO BEST DOCUMENTATION EVER!!!
    inline bool NeedsRefinement() const {return fStatus[2];}

    /// Set Status Vector (topology of intersection)
    void SetStatusVec(const TPZVec<int> &status){
        if(status.size() != 3) DebugStop();
        fStatus = status;
    }
    void SetStatusVec(const std::initializer_list<int>& status){
        if(status.size() != 3) DebugStop();
        fStatus = status;
    }

    /// Reference to status vector
    const TPZVec<int> &StatusVec() const {return fStatus;}
    
    /**
     * @brief Check geometry of intersection against a tolerance, snaps intersection 
     * to closest side(s) if necessary and modifies affected neighbours.
     * @return True if any optimization has been made.
     * @param tol: Minimum acceptable length
    */
    // @TODO Suggestion : change the name to SnapIfNecessary or
    // something else. It is not an optimization
    // Maybe split in two methods : NeedsSnap and projectnode
    // then ForceProjection does not need to use "gambiarra" with
    // bignumber
    bool Optimize(REAL tolDist = 1e-4);

    /// Check if should be refined and generate the subelements
    void Refine();

    /// Sets material id
    void SetMaterialId(int matID){
        fGeoEl->SetMaterialId(matID);
    };

    /// After optimization, update neighbours through side iside
    // @TODO BEST COMMENT EVER = Just kidding, the comment doesnt
    // explain anything...
    void UpdateNeighbours(int iside);

    /// Check if rib will be incorporated onto fracture and correct its material id as such.
    /// Returns true if changes have been made
    bool UpdateMaterial();

    /// Forces projection (optimization) of the intersection to the closest lower dimensional side
    // @TODO does it call force projection recursively?
    void ForceProjection();

    private:
        /// Creates refinement pattern based on status vector and intersection node index
        void CreateRefPattern();

};

#endif /* DFNRib_h */
