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

// Forward declarations
class DFNFracture;
enum ESplitPattern : int;

/// DFNFace class describes a mesh 2D element (face) and how it's cut by a fracture. It also carries a method to split itself into smaller sub faces.
class DFNFace
{
    
private:
	/// pointer to its original geometric element
	TPZGeoEl *fGeoEl;

	/** A status vector describes the topology of the intersection. 
     *  Each vector entry corresponds to a side, which may contain 
     *  an intersection node with the fracture.
     * @details The vector usually starts with zero everywhere and a pair of ones on two of its 1D sides (4~7 for quads and 3~5 for triangles), 
     * then snapping algorithms may move them to the nodes.\n 
     * If a face starts with nodes set to 1, then that means a tolerance was verified at the rib level before the face was created.
     * @example For example: 
     * {1,0,1,0,0,0,0,0,0} is a quadrilateral face, intersected by a fracture, whose intersection nodes have been snapped to nodes 0 and 2;
     * {0,0,0,0,1,0,1,0,0} is a quadrilateral face intersected by a fracture passing through sides 4 and 6 (the 1st and 3rd edges), where no snap was performed;
     * {1,0,0,0,1,0,0} is a triangular face intersected by a fracture passing through sides 0 and 2. Meaning one intersection was snapped to node 0 and the other wasn't snapped;
     * {0,0,0,0,0,0,0} is an uninitialized Status vector of a triangle;
     * @attention As it currently stands, our methodology doesn't have a Status Vector with more than two entries == 1; If you find this ever to be violated, it's a bug;
     * @note @phil if you are seeing an entirely zero StatusVector, you probably printed it before calling DFNFace::UpdateStatusVector, where faces build their status vector by analyzing their ribs.
     */
	TPZManVector<int,9> fStatus;

	/// Anticipated coordinates of an in-plane intersection node (if any has beeen found)
    /// This vector is only filled if this face intersects the limits of the fracture polygon
	// TPZManVector<REAL, 3> fCoord;
	
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
    DFNFace() : fGeoEl(nullptr), fStatus(9,0), fRibs(0){};
    
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
    
    /// Sets the geoel
    void SetGeoEl(TPZGeoEl* geoEl) {fGeoEl = geoEl;}
    
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
    /// @deprecated (But could be brought back to support 2D DFNs)
    // void SetIntersectionCoord(TPZManVector<REAL, 3> coord){
    //     fCoord = coord;
    // }

    /// Get in-plane intersection coordinates
    /// @deprecated (But could be brought back to support 2D DFNs)
    // TPZManVector<REAL, 3> IntersectionCoord() const {return fCoord;}
    
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
    /**
     * @brief (STATIC) Check if element should be refined
     * @return False if only one node or only two consecutive nodes have been intersected
    */
    static bool NeedsRefinement(const TPZVec<int>& statusvec);

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
    // bool FindInPlanePoint();
    
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
     * @brief Returns the split pattern that should be used to split this face.
     * @details This function is the actual definition of the SplitPatterns, you can see
     *        the comments in the enum definition for some examples
     * @param Status vector that indicates which sides are cut
     */
    static ESplitPattern GetSplitPattern(const TPZManVector<int> &status);

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

    /// Given an index of a snappable node in the RefMesh, get the local rib index that contains it as intersection node. 
    /// This perfoms floating point operations: at most 2 vector<double, 3> subtractions and the norm of the resulting vector
    /// @deprecated I found a way to do it with no flop
    int ComputeRibContainingNode(const int64_t snapNode) const;

};

inline std::ostream& operator<<(std::ostream &out, const DFNFace& face){
    face.Print(out);
    return out;
}

/** @brief Enumeration of possible SplitPatterns to be used when refining faces to 
 * conform to their intersection with the fracture.
 * 
 * @details Split patterns are statically defined, since they describe a small set of 
 * possible topologies for the refinement patterns that will be (dinamically) created.
 * Split patterns are TOPOLOGY, Refinement Patterns are GEOMETRY.
 * 
 * @example 1: For example: ESplitPattern::Quad_2_OppositeEdges is the topological 
 * construction of a fracture crossing through a 2D skeleton mesh quadrilateral element 
 * and intersecting 2 of its edges, without triggering any snap of nodes (no tolerance 
 * is violated). There's no imposition of what permutation of this construction we're 
 * handling, or on the coordinates of the intersection nodes, other than that they 
 * exist in the space of 2 opposing edges and be not snapped to a closest node.
 * 
 * @example 2: Another example: ESplitPattern::Triang_1_Edge_1_Node is the topological 
 * construction of a fracture crossing through a 2D skeleton mesh triangle element and 
 * intersecting an edge without triggering a snap, and another edge whose intersection 
 * was snapped to its closest node. 
*/
enum ESplitPattern : int{
    Uninitialized           = 0,
    None                    = 0,
    
    // Quadrilaterals ==============================================================
    Quad_2_OppositeEdges   =  1, // Quadrilateral with 2 opposite refined edges
    Quad_2_AdjacentEdges   =  2, // Quadrilateral with 2 adjacent refined edges
    Quad_2_OppositeNodes   =  3, // Quadrilateral with 2 opposite intersected corners
    Quad_1_Edge            =  4, // Quadrilateral with a single refined edge
    // Quad_1_Edge_1_MidFace  =  5, // (Deprecated) Quadrilateral with a refined edge and a mid-face intersection node 
    // Quad_1_Node_1_MidFace  =  6, // (Deprecated) Quadrilateral with an intersected corner and a mid-face intersection node
    Quad_1_Edge_1_Node     =  7, // Quadrilateral with an intersected corner and a refined edge 
    // Quad_1_MidFace         =  8, // (Deprecated) Quadrilateral with a single mid-face intersection node 
    // Quad_2_MidFace         =  9, // (Deprecated) Quadrilateral with 2 mid-face intersection nodes

    // Triangles ===================================================================
    Triang_2_AdjacentEdges  = 10, // Triangle with 2 adjacent refined edges
    Triang_1_Edge_1_Node    = 11, // Triangle with an intersected corner and a refined edge
    // Triang_1_Edge_1_MidFace= 12, // (Deprecated) Triangle with a refined edge and a mid-face intersection node
    // Triang_1_Node_1_MidFace= 13, // (Deprecated) Triangle with an intersected corner and a mid-face intersection node
    Triang_1_Edge           = 14  // Triangle with a single refined edge
    // Triang_1_MidFace       = 15, // (Deprecated) Triangle with a single mid-face intersection node
    // Trang_2_MidFace        = 16, // (Deprecated) Triangle with 2 mid-face intersection nodes
};



#endif /* DFNFace_h */
