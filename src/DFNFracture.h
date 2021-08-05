/*! 
 *  @authors   Pedro Lima
 *  @date      2018-2020
 */

#ifndef DFNFracture_h
#define DFNFracture_h

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include <stdio.h>
#include <tuple>
#include "pzmatrix.h"
#include "pzcompel.h"
#include "pzgeoelbc.h"

#include "DFNRib.h"
#include "DFNFace.h"
#include "DFNVolume.h"
#include "DFNPolygon.h"
#include "DFNMesh.h"
#include "DFNNamespace.h"
#include <gmsh.h>

class DFNMesh;
class DFNPolyhedron;
typedef TPZFMatrix<REAL> Matrix;




/** 
 *  @brief     Describes a surface mesh for a fracture and all ribs & faces that are intersected by it.
 *  @details   Intersection search is performed after creation of skeleton
 *  elements with DFNMesh::CreateSkeletonElements. Fracture plane should
 *  be a DFNPolygon.
 */
class DFNFracture
{
private:
	
	/// Pointer for the complete DFN mesh
	DFNMesh *fdfnMesh = nullptr;

	/// Map of ribs affected by this fracture {gel_index, DFNRib}
	std::map<int64_t, DFNRib> fRibs;

	/// Map of faces affected by this fracture {gel_index, DFNFace}
	std::map<int64_t, DFNFace> fFaces;

	/// A planar convex polygon that defines the outline of the fracture
	DFNPolygon fPolygon;
	
    /// A material id for this fracture elements. You want to keep it different than
    /// all other mesh matids, otherwise it'll cost you to find the fracture boundary 
    /// conditions. Also, if you want to merge fractures, just match their matids.
    int fmatid = DFNMaterial::Efracture;

    /// Index of this fracture at fdfnMesh fracture vector
    int fIndex = -1;

    /// Directive for determining if this fracture limits are truncated, extended or recovered
    FracLimit fLimit = Eextended;

	/// Set (of indices) of 2D geo elements on fracture surface
	std::set<int64_t> fSurfaceFaces;
	
	/// Set (of indices) of 1D geo elements on fracture surface. 
    /// @comment This structure is used (and built) twice. First at the recovery of fracture limits, then later at the search for fracture-fracture instersections.
	std::set<int64_t> fSurfaceEdges;

    TPZStack<int> f_frac_frac_intersections;
    
    TPZVec<std::pair<TPZGeoNode*,REAL> > fDelimitingNodes;

    int fmatid_BC = fmatid-1;

public:

#if PZ_LOG
    /// During debug, I experimented with a log file for each fracture. It's not as 
    /// useful as one would think, but I left it here if you want to use it
    // log4cxx::LoggerPtr fLogger = nullptr;

    /** @brief Creates a logger for this object and fills pointer to this->fLogger
      * @note A method to create a separate logger + appender for each DFNFracture. (as opposed to the main logger which logs everything from the DFNMesh)*/
    // log4cxx::LoggerPtr CreateLogger(std::string filename="default", std::string layout_convpattern = "default");
#endif // PZ_LOG


    /// Empty constructor
    DFNFracture();

    ///Destructor
    ~DFNFracture(){};
    
    /**
     * @brief Constructor from a DFNPolygon
     */
    DFNFracture(DFNPolygon &Polygon, DFNMesh *dfnMesh, FracLimit limithandling = Eextended, int matid = DFNMaterial::Efracture);
    
    /// Copy constructor
    DFNFracture(const DFNFracture &copy);
    
    /// Assignment operator
    DFNFracture &operator=(const DFNFracture &copy);
    
    
    /// Return the corner nodes of the fracture
    DFNPolygon &Polygon();

    /// Get pointer to the DFNMesh
    DFNMesh* dfnMesh() const{return fdfnMesh;}

    /// Get reference to the set of 2D elements at the surface of this fracture
	std::set<int64_t>& Surface(){return fSurfaceFaces;}

    int MaterialId() const{return fmatid;}
    void SetMaterialId(int matid){fmatid = matid;}
    int Index() const{return fIndex;}

    /// Number of 2D elements at the surface of this fracture
    int NSurfElements() const{return fSurfaceFaces.size();}
    
    /// return the indices of all polyhedra intersected by the fracture
    std::set<int> IdentifyIntersectedPolyhedra();

    /** @brief Given an intersected polyhedron, check if it could lead to undesirable features due to the fracture intersecting it. Use only the SnapRibs.
     * this method should return false when all the snap ribs loop around a TPZGeoEl
     * @return False if volume has N_SnapRibs == 0
     * @return False if there are no pairs of neighbour SnapRibs
     * @return False if all pairs of neighbour SnapRibs are co-linear (same EldestAncestor)
     * @return False if VolumeSnapRibs perfectly cover (no more, no less) their eldest ancestors AND the set of rib elders form a closed loop of 3 or 4 edges (which means we've snapped onto a mesh face which will latter be incorporated during DFNFracture::MeshFractureSurface)
     * @return True(?) if there are 2 non-colinear neighbour SnapRibs (?). Not really, take the following construction: Triangulate a sphere, cut through it with a plane perfectly splitting it in 2 halfs. Let all intersections be snapped. That is not a problem volume, but this condition would return true.
     */
    bool IsProblemVolume(const std::set<int64_t>& AllSnapRibs, const DFNPolyhedron& IntersectedVolume) const;
    
    // return the indices of the geometric elements that belong to the discretized fracture
    // through snapping AND were in the original mesh
    /** @brief Fill a set with all SnapRibs due to this Fracture
     *  @pre Should be called after DFNFracture::SnapIntersections_faces. The set will often be empty if you fail to do so. It will ALWAYS be empty if called before DFNFracture::SnapIntersections_ribs
     *  @details SnapRibs are the 1D elements that previously existed in the mesh and this DFNFracture will try to incorporate to its surface. Not to be confused with DFNRibs (which are the 1D elements that were intersected).
     *  @note This new 'SnapRibs' set is a subset of DFNFracture::fSurfaceEdges which is temporarily assembled during the limit recovery of this fracture. So there's some work being done twice, which we could optimize later.
    */
    std::set<int64_t> IdentifySnapRibs();
    
    /** facet set of 2d elements with the same normal direction
     * divide facets of polyhedra that have non-colinear snap-ribs
     * loop over polyhedra
     * identify the snap-ribs the belong to the polyhedra
     * verify if the snap-ribs require meshing of a facet (hard)
     * divide the facet (hard part)
     * @deprecated
    */
    void MeshSnapPlanes();

    /**
     * @brief Checks if this fracture is overlapping (partially, or completely) planar subsets
     * of the shell of every DFNPolyhedron that it intersects. 
     * @details If any trivial case is found, code moves on and the incorporation of any 
     * overlapped TPZGeoEl is done at DFNFracture::FindPolygon() (as it used to). If any case is
     * found that FindPolygon can't handle, then the DFNPolyhedron is split into tetrahedra and 
     * pyramids by Gmsh, and we recursively update the set of faces and ribs intersected until 
     * no non-trivial overlap is left.
     * @note This function is somewhat costly, if you know it won't be necessary, you can just
     * not call it
    */
    void CheckSnapInducedOverlap();
    
    // /** @briefUpdates DFNFracture's DFNFace list, swapping elements of a given polyhedron for their children
    //  * @param polyh : A polyhedron from where to get refined DFNFaces
    //  * @ATTENTION this method will fail if you call it after DFNPolyhedron::SwapForChildren. SwapForChildren is called during DFNMesh::UpdatePolyhedra. A bug from this could be very hard to catch.
    // */
    // void SwapRefinedDFNFaces(DFNPolyhedron& polyh);

    /** @brief Removes refined DFNFaces from the fFaces datastructure
    */
    void RemoveRefinedDFNFaces(const int vol_index);

private:

    
    void Initialize(DFNPolygon &Polygon, DFNMesh *dfnMesh, FracLimit limithandling = Eextended, int matid = DFNMaterial::Efracture);
        
    /// Checks neighbour's dimension and returns true if it is equal
    bool HasEqualDimensionNeighbour(TPZGeoElSide &gelside);
    
    /// Finds intersection point of fracture boundaries and geometric mesh faces
    bool FindEndFracturePoint(DFNFace &face, TPZManVector<REAL,3> &ipoint);
    

    // void ImportElementsFromGMSH(TPZGeoMesh * gmesh, int dimension, std::set<int64_t> &oldnodes, TPZStack<int64_t> &newelements);

    /**
     * @name Surface Meshing auxiliar methods
     * @{
    */
    /** @brief Mesh a convex polygon from a list of sequentialy connected edges. If not simple, calls on Gmsh
     * @param polygon a loop of edges that don't necessarily occupy the same plane
    */
    void MeshPolygon(TPZStack<int64_t>& polygon);
    /** @brief Projects a non-planar polygon onto its best fitting plane and uses Gmsh to mesh it
     * @param orientedpolygon an oriented loop of edges that don't necessarily occupy the same plane
    */
    void MeshPolygon_GMSH(TPZStack<int64_t>& orientedpolygon,std::set<int64_t>& nodes, TPZStack<int64_t>& newelements, bool isplane=false);
    /** @brief Recursively finds the next face and gets its in-plane edge to build a SubPolygon (subset of the Fracture DFNPolygon contained in a polyhedron)
     *  @param Polygon_per_face a structure to store the two subpolygons per face
    */
    void BuildSubPolygon(TPZVec<std::array<int, 2>>& Polygon_per_face,
                        std::pair<int64_t,int> currentface_orient,
                        int inlet_side,
                        TPZStack<int64_t>& subpolygon);

    /// @brief Given a side and an oriented face, get a neighbour that shares the same polyhedron
    std::pair<int64_t,int> PolyhNeighbour(std::pair<int64_t,int>& currentface_orient, int currentside, int& neigside);

    /** @brief Setup the edges that form the subpolygon as an oriented loop in Gmsh fashion*/
    void SetLoopOrientation(TPZStack<int64_t>& edgelist);
    
    /// set a subpolygon index for a face in the structure Polygon_per_face
    void SetPolygonIndex(std::pair<int64_t,int> face_orient, int polyg_index,TPZVec<std::array<int, 2>>& Polygon_per_face);
    /// get the subpolygon index for a face from the structure Polygon_per_face
    int GetPolygonIndex(std::pair<int64_t,int> face_orient,const TPZVec<std::array<int, 2>>& Polygon_per_face);
    /// Builds and fills a list with this fracture outer loop of edges
    void GetOuterLoop(std::vector<int> &outerLoop);

    /** @brief Insert 2D elements in the set of 2D elements in surface */
    void InsertFaceInSurface(int64_t elindex);
    void InsertFaceInSurface(const TPZVec<int64_t>& el_indices){
        for(int64_t index : el_indices){InsertFaceInSurface(index);}
    }
    /** @brief Insert 1D elements in the set of 1D elements in surface */
    void InsertEdgeInSurface(int64_t elindex){fSurfaceEdges.insert(elindex);}
    void InsertEdgeInSurface(const TPZVec<int64_t>& el_indices){
        for(int64_t index : el_indices){InsertEdgeInSurface(index);}
    }
    /// Removes negative integers from a stack
    void ClearNegativeEntries(TPZStack<int64_t>& subpolygon);

    /** @} */

    // Special setup of fracture mat id when working in 2D
    void SetFracMaterial_2D();

    /** @brief For every child of every face in fFaces, check if it's above or below the fracture plane. It changes the material id of the elements and corrects the real fracture datastructures accordingly
     * @attention This method was written to be called by a virtual fracture. One that is auxiliar orthogonal to a real fracture.
    */
    void SortFacesAboveBelow(int id_above, int id_below, DFNFracture& realfracture);
    void RemoveFromSurface(TPZGeoEl* gel);
    void AddToSurface(TPZGeoEl* gel);

    /** @brief Check if face is above or below fracture surface
     * @param use_face_centroid: true -> check using the centroid of the face; false -> check using centroid of the edges of the face.
    */
    bool CheckFaceAbove(TPZGeoEl* face, bool use_face_centroid);


    void ResetSurfaceMaterial(const int matid);

public:

    /// @brief Check if 2D element on surface has at most 1 other neighbour (through edges) on this fracture surface
    bool CheckIsLegalSurfaceElement(const int64_t elindex) const;

    void PlotVTK(const std::string exportname, bool putGraphicalElements = true);
    void SetupGraphicsFractureIntersections(TPZStack<int>& fracfrac_int);
    void SetupGraphicsFractureBC();

    /// @brief from a set of 1D elements find if they form a lineloop of an existing 2D element in the mesh
    TPZGeoEl* FindPolygon(TPZStack<int64_t>& polygon);

    /// Triangulates fracture surface from outline
    void MeshFractureSurface();

    /** @brief Following directive DFNFracture::fLimit, modify the extended 
     * limits of this fracture to better represent the limits (boundary edges) 
     * of the polygon that defined this fracture*/
    void RecoverFractureLimits();

    /// Access the ribs data structure
    DFNRib* AddRib(DFNRib &rib);
    
    DFNFace* AddFace(DFNFace &face);

    /// Pointer to rib corresponding to geometric element with index
    // return NULL if the geometric element is not intersected
    DFNRib *Rib(int64_t index);

    /// Pointer to face of index 'index'
    DFNFace *Face(int64_t index);
    
    /** @brief Search for intersected 1D elements and create DFNRib objects for them;
     *  @details Loop over all 1D elements in the mesh; use DFNPolygon::fNodesAbove to check if an edge has nodes on opposite sides of the plane that contains the DFNPolygon; and check if the intersection is within the DFNPolygon bounds.
     **/
    void CreateRibs();
    /// Find and intersect ribs within the element indices in a set.
    /// @attention Unlike DFNFracture::CreateRibs(), this method won't check if intersection point is within fPolygon limits. Since ribs are already limited to a set, intersection search is done with the unbounded plane that contains fPolygon
    void FindRibs(const std::set<int64_t>& ribset);
    
    /// verify proximity of rib intersection node
    /// Coalesce intersected ribs
    void SnapIntersections_ribs(REAL tolDist = -1.0);
    
    /// Set Refinement Patterns and create sub elements
    void RefineRibs();

    /// Find intersected faces
    void CreateFaces();
    /// Coalesce intersected faces
    void SnapIntersections_faces(REAL tolDist = -1.0, REAL tolAngle = -1.0);
    /// Set Refinement Patterns and create sub elements
    void RefineFaces();

    /// Remove refined faces from the surface gelindex set and replace them by their youngest children
    void UpdateFractureSurface();


    /**
     * @name Handling frature limits
     * @{
    */
    /** @brief Identify Ribs, Faces and Polyhedra that are affected by the limits of the fracture
     * @note Should be called after DFNFracture::CreateFaces()
     * @details This method is a preparation for the actual recovery of fracture limits, which happens during DFNFracture::RecoverFractureLimits
    */
    void IsolateFractureLimits();
    /** @brief Find Ribs that lie out of the bounds of this fracture but are still influenced by its limits */
    void FindOffboundRibs();
    /** @brief Find Faces that lie out of the bounds of this fracture but are still influenced by its limits */
    void FindOffboundFaces();
    /** @brief Fill a set with every 1D element that is contained in the surface of this fracture*/
    void GetEdgesInSurface(std::set<int64_t>& edges);
    /** @brief Creates a fracture (orthogonal to this) that contains the specified edge and whose normal vector is oriented towards the interior of this
     * @param orthfracture [out] The reference to the created fracture
     * @param edgeindex Index of the edge of the original fracture that is calling this method
     * @details its polygon plane is a triangle with the nodes of the specified edge and a node above the real fracture
    */
    void CreateOrthogonalFracture(DFNFracture& orthfracture, const int edgeindex);

    /** @brief Clears material ids that might have been changed by other fractures and updates fracture surface list of elements
     * @note Because limit recovery is done after all fractures have been inserted, some cleaning is necessary*/
    void CleanUp();

    /** @} */

    /** @brief Print method for logging */
    void Print(std::ostream& out = std::cout) const;
    
    /**
     * Modifies the matids of the TPZGeoElBCs on the boundary. This lets one easily identify what elements
     * are in the boundary for future assesment of boundary conditions
     */
    void ModifyBoundaryMatID(const int matidboundary);
    
    /**
     * For a given TPZGeoEl, checks if any of its nodes are one of the delimiting nodes of the user given DFNPolygon
     * TODO: This method uses geometrical distances for now and can be optimized. Also, there is hardcoded tolerance that may not be robust
     */
    void CheckIfNodesAreDelimitingNodes(TPZGeoEl* gel);
                
        
    /** @brief This is a draft to find the boundary conditions for fractures. I had the idea and quickly wrote down... but have not tested yet.
     * Supposedly, the way it's written, it should allow a user to essentially merge 2 or more fractures into one.
     * Meaning, they would have a continuous boundary 
    */
    void ExportFractureBC(int matid, std::ofstream& out);
    
    /**
     * Function to try to quickly find a GeoEl with the required matid that is the boundary
     */
    TPZGeoEl* findBCGeoElWithMatID(const int matid) const;

    /**
     * Checks is a certain node is a delimiting node from DFNPolygon
     */
    const bool IsDelimitingNode(TPZGeoNode *node) const;

    /**
     * Loops over the GeoEls closed list to find one that start on a delimiting node.
     */
    TPZGeoEl* FindGeoStartingOnDelimitingNode(TPZGeoEl* gel, const int matid) const;
    
    /**
     * Function to export boundary conditions that are different in each side of the fracture. May not work with some cases that have fracture intesections
     */
    void ExportFractureBCDifferentTags(int matid, std::ofstream& out);

    /** @brief Find edges that define the intersection between this and another DFNFracture, by solving a shortest
     * path in graph problem to best approximate a segment.
     * @param OtherFrac The other fracture the code should search for intersection
     * @param Segment 2 coordinates representing an initial and a final point that graph should approximate
     * @param EdgeList [output] A stack to fill with the result of the search
     * @details This code performs 2 geometrical searches to get, within the set of nodes of the 
     * CommonFaces, the mesh nodes that are closest to the coordinates informed in the Segment
    */
    bool FindFractureIntersection_NonTrivial(const DFNFracture& OtherFrac, 
                                            // const std::set<int64_t>& CommonFaces, 
                                            const Segment& Segment,
                                            TPZStack<int64_t>& EdgeList);

    /** @brief Find edges that define the intersection between this and another DFNFracture, by searching for
     * 1D sides of the faces in this fracture's surface through where a neighbour can be found with the matid
     * of the OtherFrac.
     * @param OtherFrac The other fracture the code should search for intersection
     * @param EdgeList [output] A stack to fill with the result of the search
    */
    void FindFractureIntersection_Trivial(const DFNFracture& OtherFrac, TPZStack<int64_t>& EdgeList);
};

inline std::ostream& operator<<(std::ostream &out, const DFNFracture& fracture){
    fracture.Print(out);
    return out;
}

namespace DFN{
    static FracLimit StringToFracLimit(const std::string& name){
        if(name[0] != 'E'){DebugStop();}
        if(name == "Etruncated"){	return FracLimit::Etruncated;}
        if(name == "Eextended"){	return FracLimit::Eextended;}
        if(name == "Erecovered"){	return FracLimit::Erecovered;}
        PZError << "\nUnrecognized FracLimit directive type: \""<<name<<"\"\n";
        DebugStop();
        return FracLimit::Eextended;
    }
}

#endif /* DFNFracture_h */

