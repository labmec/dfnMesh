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
	TPZManVector<int, 3> fStatus;

	/// Anticipated coordinates of intersection (may not be real coordinates)
	TPZManVector<REAL, 3> fCoord;
	
	/// Index of intersection node
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

    /// Element index
    int64_t Index() const {return fGeoEl->Index();}
    
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
    */
    inline bool IsIntersected() const{
        return (fStatus[0] || fStatus[1] || fStatus[2]);
    }

    /// Check if this rib should be refined
    inline bool NeedsRefinement() {return fStatus[2];}

    /// Set Status Vector (topology of intersection)
    void SetStatusVec(TPZManVector<int, 3> status){
        if(status.size() != 3) DebugStop();
        fStatus = status;
    }
    void SetStatusVec(const std::initializer_list<int>& status){
        if(status.size() != 3) DebugStop();
        fStatus = status;
    }

    /// Reference to status vector
    TPZManVector<int, 3> &StatusVec(){return fStatus;}
    
    /**
     * @brief Check geometry of intersection against a tolerance, snaps intersection 
     * to closest side(s) if necessary and modifies affected neighbours.
     * @return True if any optimization has been made.
     * @param tol: Minimum acceptable length
    */
    bool Optimize(REAL tolDist = 1e-4);

    /// Check if should be refined and generate the subelements
    void Refine();

    /// Sets material id
    void SetMaterialId(int matID){
        fGeoEl->SetMaterialId(matID);
    };

    /// After optimization, update neighbours through side iside
    void UpdateNeighbours(int iside);

    private:
        /// Creates refinement pattern based on status vector and intersection node index
        void CreateRefPattern();

};

#endif /* DFNRib_h */
