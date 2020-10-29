/*! 
 *	DFNFace.hpp
 *  @authors   Pedro Lima
 *  @authors   Philippe Devloo
 *  @date      2018-2020
 */

#ifndef DFNFace_hpp
#define DFNFace_hpp

#include <stdio.h>
#include "pzmanvector.h"
#include "pzmatrix.h"
#include "pzcompel.h"
#include "pzgeoelbc.h"
#include "DFNRib.h"

class DFNFracture;

/// DFNFace class describes a plane and how it's cut by a fracture. It also carries a method to split itself into smaller sub faces.
// @TODO I dont understand the meaning of Face based on this description>
// It represents the intersection of the fracture polygon with a given volume?
class DFNFace
{
    
private:
	/// pointer to its original geometric element
	TPZGeoEl *fGeoEl;

	/** A status vector describes the topology of the intersection
     *  each vector entry corresponds to a side, which may contain 
     *  an intersection node with the fracture
     * @TODO phil : what is the meaning if the status vector is entirely ZERO
     */
	TPZManVector<int> fStatus;

	/// Anticipated coordinates of an in-plane intersection node (if any has beeen found)
    /// This vector is only filled if this face intersects the limits of the fracture polygon
	TPZManVector<REAL, 3> fCoord;
	
	/// Index of in-plane intersection node
	int64_t fIntersectionIndex = -1;
	
	/// Pointer to a fracture mesh
	DFNFracture *fFracture = nullptr;

	/// Vector with pointers to respective DFNRibs (if there's not a rib at that edge, it'll contain a nullptr)
	TPZManVector<DFNRib*, 4> fRibs;

    /// Refinement mesh is used to create a refinement pattern and decide how to optimize it
    TPZGeoMesh fRefMesh;

public:
    /// Empty constructor
    DFNFace() : fGeoEl(nullptr), fStatus(9,0), fRibs(0), fCoord(0){};
    
    /// Default constructor takes a pointer to geometric element and a fracture
    DFNFace(TPZGeoEl *gel, DFNFracture *fracture);
    
    /// Copy constructor
    DFNFace(const DFNFace &copy);
    
    /// Assignment operator
    DFNFace &operator=(const DFNFace &copy);

    /// Destructor
    ~DFNFace(){};
    
    /** @brief Print method for logging */
    void Print(std::ostream& out = std::cout, bool print_refmesh = false) const;
    
    /// Pointer to geometric element
    TPZGeoEl* GeoEl() const {return fGeoEl;}
    
    /// Element index
    int64_t Index() const {return fGeoEl->Index();}

    /**
     * @brief Return pointer to i-th DFNrib.
     * @returns nullptr if rib isn't in Fracture's rib map
    */
    DFNRib *Rib(int i);
    
    /// Set index of in-plane intersection node
    void SetIntersectionIndex(int64_t index){fIntersectionIndex = index;}

    /// Get index of in-plane intersection node
    int64_t IntersectionIndex() const {return fIntersectionIndex;}
    
    /// Set anticipated in-plane intersection coordinates (real coordinates are defined by DFNRib::SnapIntersection_try())
    void SetIntersectionCoord(TPZManVector<REAL, 3> coord){
        fCoord = coord;
    }

    /// Get in-plane intersection coordinates
    TPZManVector<REAL, 3> IntersectionCoord() const {return fCoord;}
    
    /// Set ribs of face
    void SetRibs(TPZVec<DFNRib*> &RibVec) {fRibs = RibVec;}

    /// Get ribs of face
    TPZVec<DFNRib*> GetRibs() const {return fRibs;}

    /// Give face a pointer to which fracture is cutting it
    void SetFracture(DFNFracture *Fracture){fFracture = Fracture;}

    /**
     * @brief Check if element should be refined
     * @return False if only one node or only two consecutive nodes have been intersected
    */
    bool NeedsRefinement() const{
        int nsides = fGeoEl->NSides();
        if(fStatus[nsides-1]) return true;
        int nnodes = fGeoEl->NCornerNodes();
        int nedges = nnodes;
        int cutedges = 0;
        int cutnodes = 0;
        bool consecutive_nodes = false;
        for(int i=0; i<nedges; i++){
            cutnodes += fStatus[i];
            cutedges += fStatus[i+nnodes];
            if(fStatus[i]&&fStatus[(i+1)%nnodes]){
                // check anterior and posterior edges for intersection
                consecutive_nodes = (fStatus[i-1+nnodes]==0)&&(fStatus[(i+1)%nnodes+nnodes]==0);
                //                             anterior edge                   posterior edge
            }
        }
        if(cutnodes==3) consecutive_nodes = false;

        if(consecutive_nodes && cutedges<=1) return false;
        if(cutnodes == 1 && cutedges == 0) return false;
        return true;
    }

    /// Return true if this face is in contact at all with the fracture
    bool IsIntersected() const{
        int nsides = fGeoEl->NSides();
        for(int i=0; i<nsides; i++){
            if(fStatus[i]) return true;
        }
        return false;
    }

    /**
     * @brief Checks if face is intersected by one of the fracture edges
     * @note A face on the fracture boundary will only have one node/edge intersected
     * @TODO Change the name? I associate boundary with boundary of the domain IntersectsFractureBoundary?
    */
    bool IsOnBoundary();

    /**
     * @brief Searches for coordinates of a possible in-plane intersection point.
     * @note 1: Assumes user has run IsOnBoundary() and it returned true;
     * @note 2: Won't create a node in mesh. This must be done later.
     * @TODO this could be done later
    */
    bool FindInPlanePoint();
    
    /// Reference to status vector
    const TPZVec<int> &StatusVec(){return fStatus;}
    
    /**
     * @brief Create/update face status vector from ribs' status
     * @returns True if any changes have been made
     */
    bool UpdateStatusVec();

    /**
     * @brief Uses status vector, in-plane point coordinates and geometric element topology
     * to create/update the refinement mesh that describes how this face will be refined
     * @returns True if any changes have been made
     * @TODO as there will be no inplane point, this can be considerably simplified
    */
    void UpdateRefMesh();

    void UpdateMaterial();

    /**
     * @brief Returns the split pattern that should be used to split this face
     * @param Status vector that indicates which sides are cut
     * @return Integer that indicates which split pattern to use. (check documentation)
     * @TODO I believe this can be simplified if only the ribs are divided!! Please update the code
     */
    int GetSplitPattern(const TPZManVector<int> &status) const;

    /**
     * @brief Check geometry of intersection against a tolerance, snaps intersection 
     * to closest side(s) if necessary and modifies affected neighbours.
     * @return True if any optimization has been made.
     * @param tolDist: Minimum acceptable distance
     * @param tolAngle_cos: Cosine of minimum acceptable angle
     * @todo
    */
    bool SnapIntersection_try(REAL tolDist = 1e-4, REAL tolAngle_cos = 0.9);

    /**
     * @brief Force snap intersection nodes down to closest nodes
    */
    bool SnapIntersection_force();

    /**
     * @brief Check geometry of intersection against tolerances, to test if intersections should be snapped
     * @return False if no tolerance is violated or if intersections have already been snapped.
     * @param tolDist: Minimum acceptable distance
     * @param tolAngle_cos: Cosine of minimum acceptable angle
    */
    bool NeedsSnap(REAL tolDist = 1e-4, REAL tolAngle_cos = 0.9);

    /// Return false if all intersections have already been snapped
    bool CanBeSnapped();

    /// Check if should be refined and generate the subelements of material id matID
    void Refine();

    /// After optimization, update neighbours through side iside
    void UpdateNeighbours(int iside);

    /**
     * @brief If intersection of fracture with this face results in a 1D element, get its index.
     * @return Index of 1D element, -1 if non existant or if line has not been created yet
     */ 
    int64_t LineInFace();

    
    /// @brief Returns the first 1D side that contains a DFNRib 
    int FirstRibSide(){
        int nedges = fGeoEl->NSides(1);
        for(int iedge=0; iedge<nedges; iedge++){
            if(fRibs[iedge]) return iedge + fGeoEl->NSides(0);
        }
        DebugStop();
        return -1;
    }

    /// @brief Returns the other 1D side with a DFNRib
    int OtherRibSide(int inletside);

    /// @brief Make children of this element inherit its polyhedral indices
    /// @note: Does nothing if mesh is not 3D
    void InheritPolyhedra();

private: 

    /**
     * @brief Determines split case and fill nodes of children and indices of the intersection points
    */
    void FillChildrenAndNewNodes(TPZManVector<TPZManVector<int64_t,4>,6> &child, TPZManVector<TPZManVector<REAL,3>> &newnode);

};
#endif /* DFNFace_h */
