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
#include "DFNNamespace.h"
#include "DFNRib.h"

class DFNFracture;

/// DFNFace class describes a mesh 2D element (face) and how it's cut by a fracture. It also carries a method to split itself into smaller sub faces.
class DFNFace
{
    
private:
	/// pointer to its original geometric element
	TPZGeoEl *fGeoEl;

	/** A status vector describes the topology of the intersection
     *  each vector entry corresponds to a side, which may contain 
     *  an intersection node with the fracture.
     * @details The vector usually starts with 1 entries on its 1D sides (4~7 for quads and 3~5 for triangles), then snapping algorithms may move them to the nodes. If a face starts with nodes set to 1, then that means a tolerance was verified at the rib level before the face was created.
     * For example: 
     * {1,0,1,0,0,0,0,0,0} is a quadrilateral face intersected by a fracture, whose intersection nodes have been snapped to nodes 0 and 2;
     * {0,0,0,0,1,0,1,0,0} is a quadrilateral face intersected by a fracture passing through sides 4 and 6 (the 1st and 3rd edges), where no snap was performed;
     * {1,0,0,0,1,0,0} is a triangular face intersected by a fracture passing through sides 0 and 2. Meaning one intersection was snapped to node 0 and the other wasn't snapped;
     * {0,0,0,0,0,0,0} is an uninitialized Status vector of a triangle;
     * @attention As it currently stands, our methodology doesn't have a Status Vector with more than two entries == 1; If you find this ever to be violated, it's a bug;
     * @note @phil if you are seeing an entirely zero StatusVector, you probably printed it before calling DFNFace::UpdateStatusVector, where faces build their status vector by analyzing their ribs.
     */
	TPZManVector<int,9> fStatus;

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
    DFNFace(TPZGeoEl *gel, DFNFracture *fracture, TPZVec<DFNRib *> &ribvec);
    
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
    DFNRib *Rib(int i) const;
    
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

    void AddRib(DFNRib* ribptr, int side);

    /// Get ribs of face
    TPZVec<DFNRib*> GetRibs() const {return fRibs;}

    /// Give face a pointer to which fracture is cutting it
    void SetFracture(DFNFracture *Fracture){fFracture = Fracture;}

    /**
     * @brief Check if element should be refined
     * @return False if only one node or only two consecutive nodes have been intersected
    */
    bool NeedsRefinement() const;

    /// Return true if this face is in contact at all with the fracture
    bool IsIntersected() const{
        int nsides = fGeoEl->NSides();
        for(int i=0; i<nsides; i++){
            if(fStatus[i]) return true;
        }
        return false;
    }

    /// Check if intersections were snapped to the same node
    bool AllSnapsWentToSameNode() const;

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
    bool SnapIntersection_try(REAL tolDist = DFN::default_tolDist, REAL tolAngle_cos = DFN::default_tolCos);

    /**
     * @brief Force snap intersection nodes down to closest nodes
    */
    bool SnapIntersection_force(const int edge_index);

    /**
     * @brief Check geometry of intersection against tolerances, to test if intersections should be snapped
     * @return False if no tolerance is violated or if intersections have already been snapped.
     * @param edge_index: (out) Fill with index of the rib that needs snap. If doesn't need snap, will get attributed -1
     * @param tolDist: Minimum acceptable distance
     * @param tolAngle_cos: Cosine of minimum acceptable angle
    */
    bool NeedsSnap(int& edge_index, REAL tolDist = DFN::default_tolDist, REAL tolAngle_cos = DFN::default_tolCos);

    /// Return false if all intersections have already been snapped
    bool CanBeSnapped() const;

    /// Check if should be refined and generate the subelements of material id matID
    void Refine();

    /// After optimization, update neighbours through side iside
    void UpdateNeighbours(int iside);

    /**
     * @brief If intersection of fracture with this face results in a 1D element, get its index.
     * @return Index of 1D element, -1 if non existant or if line has not been created yet
     */ 
    int64_t LineInFace() const;
    
    /**
     * @brief If the intersections of this face were snapped to 2 consecutive nodes, get the local side index for the 1D side that connects those nodes.
     * @return Side Index, if there exists an assimilated side, get -1 if non existent;
     */ 
    int CheckAssimilatedSide() const;

    
    /// @brief Returns the first 1D side that contains a DFNRib 
    int FirstRibSide() const;
    /// @brief Returns the second 1D side that contains a DFNRib 
    int SecondRibSide() const;

    /// @brief Returns the other 1D side with a DFNRib
    int OtherRibSide(int inletside) const;

    /// @brief Make children of this element inherit its polyhedral indices
    /// @note: Does nothing if mesh is not 3D
    // void InheritPolyhedra();

    int NIntersectedRibs() const;

    /** @return Number of ribs whose intersection point is within the fracture limits*/
    int NInboundRibs() const;

    /// Return number of intersected edges whose intersection weren't yet snapped
    int NSnappableRibs(int& first_rib_localindex) const;

    void SketchStatusVec(std::ostream& out = std::cout) const;
private: 

    /**
     * @brief Determines split case and fill nodes of children and indices of the intersection points
     * @details 1. Children's normal vector is guaranteed to match their father's;
    */
    void FillChildrenAndNewNodes(TPZManVector<TPZManVector<int64_t,4>,6> &child, TPZManVector<TPZManVector<REAL,3>> &newnode);

    /// Given an index of a snappable node in the RefMesh, get the local rib index that contains it as intersection node
    /// This perfoms floating point operations: at most 2 vector<double, 3> subtractions and the norm of the resulting vector
    int ComputeRibContainingNode(const int64_t snapNode) const;

};
};
#endif /* DFNFace_h */
