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
    const DFNPolygon& Polygon() const {return fPolygon;}

    /// Get pointer to the DFNMesh
    DFNMesh* dfnMesh() const{return fdfnMesh;}

    /// Get reference to the set of 2D elements at the surface of this fracture
	std::set<int64_t>& Surface(){return fSurfaceFaces;}
    const std::set<int64_t>& Surface() const {return fSurfaceFaces;}

    /// Get reference to the set of 1D elements at the surface of this fracture
    std::set<int64_t>& SurfaceEdges(){return fSurfaceEdges;}
    const std::set<int64_t>& SurfaceEdges() const {return fSurfaceEdges;}
    
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

    /// @brief Debugging purpose function for plotting intersection between fractures
    void PlotVTK_SharedSurface(const std::string& filepath, DFNFracture& otherfrac, const Segment& seg);
    
    void Initialize(DFNPolygon &Polygon, DFNMesh *dfnMesh, FracLimit limithandling = Eextended, int matid = DFNMaterial::Efracture);
        
    /// Checks neighbour's dimension and returns true if it is equal
    bool HasEqualDimensionNeighbour(TPZGeoElSide &gelside);
    
    /// Plot VTK for a SubPolygon
    /// @param relativefolderpath Relative path for output folder from directory ./LOG
    void PlotVTK_SubPolygon(const TPZVec<int64_t>& subpolygon, const int volumeindex, const std::string relativefolderpath) const;
    
    /**
     * @name Surface Meshing auxiliar methods
     * @{
    */
    /** @brief Mesh a convex polygon from a list of sequentialy connected edges. If not simple, calls on Gmsh
     * @param polygon a loop of edges that don't necessarily occupy the same plane
     * @param polyhindex passed only for debugging purposes in case one needs to plot the polyhedron
     * @param newelements newly created elements based on the intersection
    */
    void MeshPolygon(TPZStack<int64_t>& polygon, const int polyhindex, TPZStack<int64_t>& newelements);
    
    /**
     * @brief Check if there are any duplicated edge elements in the subpolygon. If so, adds the local indices of
     * the repeated entries to the return set
     * @param subpolygon a loop of edges elements that don't necessarily occupy the same plane
     * @param locDuplicateIndices local indices of duplicate edges in subpolygon
     * @return true if there are duplicate edges
     */
    const bool CheckForDuplicateEdges(const TPZStack<int64_t>& subpolygon,
                                      std::set<int>& locDuplicateIndices) const;
    
    /**
     * @brief Clears duplicate edges in subpolygon based on the local indices in locDuplicateIndices
     * @param locDuplicateIndices local indices of duplicate edges in subpolygon
     * @param subpolygon a loop of edges elements that don't necessarily occupy the same plane
     */
    void ClearDuplicateEdges(const std::set<int>& locDuplicateIndices, TPZStack<int64_t>& subpolygon) const;

    /**
     * @brief Checks if subpolygon forms a closed loop
     * @param subpolygon a loop of edges elements that don't necessarily occupy the same plane
     * @param badVolumes list of volume polyhedras that have small angles and need to be refined into simplexes during rollback
     */
    const bool CheckIfPolygonIsClosedLoop(const TPZStack<int64_t>& subpolygon) const;

    
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

    /**
     * @brief Checks if a subpolygon can be considered planar "enough". Tolerance when this method was created was 45 degrees. 
     */
    const bool CheckSubPolygonPlanarity(TPZStack<int64_t>& subpolygon, const int polyhindex) const;
    
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

    /** @brief Sort 2D elements on the surface of realfracture as being either above or below an orthogonal plane. Removes from surface those that are found to be 'below'
     * @attention This method was written to be called by a virtual fracture. One that is auxiliar and orthogonal to a real fracture.
    */
    void SortFacesAboveBelow(int id_above, int id_below, DFNFracture& realfracture);
    void RemoveFromSurface(TPZGeoEl* gel);
    void RemoveFromSurface(const TPZVec<int64_t>& indices);
    void AddToSurface(TPZGeoEl* gel);
    void AddToSurface(const std::set<int64_t>& indices);
    /// @brief Adds elements to the fracture surface by index. If an element was already on surface, removes it instead.
    void AddOrRemoveFromSurface(const std::set<int64_t>& indices);

    /** @brief Check if face is above or below fracture surface
     * @param use_face_centroid: true -> check using the centroid of the face; false -> check using centroid of the edges of the face.
    */
    bool CheckFaceAbove(TPZGeoEl* face, bool use_face_centroid);


    void ResetSurfaceMaterial(const int matid);
    
    void CreateSetsOfContiguousEdges(std::set<int64_t> common_edges,
                                     TPZManVector<std::set<int64_t>,2> &contiguousEdges) const;
    
    void ComputeStartAndEndOfSetOfEdges(const std::set<int64_t>& common_edges,
                                        const Segment& Segment,
                                        int64_t& start, int64_t& end) const;

public:

    /// @brief Check if 2D element on surface has at most 1 other neighbour (through edges) on this fracture surface
    bool CheckIsLegalSurfaceElement(const int64_t elindex) const;

    void PlotVTK(const int surface_matid, const std::string exportname,
                 bool putGraphicalElements = true, bool plotIntersections = true);
    void SetupGraphicsFractureIntersections(TPZStack<int>& fracfrac_int);
    void SetupGraphicsFractureBC();

    /// @brief from a set of 1D elements find if they form a lineloop of an existing 2D element in the mesh
    TPZGeoEl* FindPolygon(const TPZStack<int64_t>& polygon) const;
    
    /**
     * @brief Triangulates fracture surface from outline
     * @param badVolumes list of volume polyhedras that have small angles and need to be refined into simplexes during rollback
     */
    void MeshFractureSurface(TPZStack<int> &badVolumes);

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
    void GetEdgesInSurface(){GetEdgesInSurface(fSurfaceEdges);}
    /** @brief Creates a fracture (orthogonal to this) that contains the specified edge and whose normal vector is oriented towards the interior of this
     * @param orthfracture [out] The reference to the created fracture
     * @param edgeindex Index of the edge of the original fracture that is calling this method
     * @details its polygon plane is a triangle with the nodes of the specified edge and a node above the real fracture
    */
    void CreateOrthogonalFracture(DFNFracture& orthfracture, const int edgeindex);

    /** @brief Clears material ids that might have been changed by other fractures and updates fracture surface list of elements
     * @note Because limit recovery is done after all fractures have been inserted, some cleaning is necessary*/
    void CleanUp(int surface_matid);

    /** @} */

    /** @brief Print method for logging */
    void Print(std::ostream& out = std::cout) const;

    /** @brief This is a draft to find the boundary conditions for fractures. I had the idea and quickly wrote down... but have not tested yet.
     * Supposedly, the way it's written, it should allow a user to essentially merge 2 or more fractures into one.
     * Meaning, they would have a continuous boundary 
    */
    void ExportFractureBC(int matid, std::ofstream& out);

    /** @brief Find edges that define the intersection between this and another DFNFracture, by solving a shortest
     * path in graph problem to best approximate a segment.
     * @param OtherFrac The other fracture the code should search for intersection
     * @param EdgeList [output] A stack to fill with the result of the search
     * @details This code performs 2 geometrical searches to get, within the set of nodes of the 
     * CommonFaces, the mesh nodes that are closest to the coordinates informed in the Segment
    */
    bool FindFractureIntersection(DFNFracture& OtherFrac, 
                                  TPZStack<int64_t>& EdgeList);

    /** @brief Find edges that define the intersection between this and another DFNFracture, by searching for
     * 1D sides of the faces in this fracture's surface through where a neighbour can be found with the matid
     * of the OtherFrac.
     * @param OtherFrac The other fracture the code should search for intersection
     * @param EdgeList [output] A stack to fill with the result of the search
    */
    void FindFractureIntersection_Trivial(const DFNFracture& OtherFrac, TPZStack<int64_t>& EdgeList);
    
    /**
     * @brief returns data structure to state before starting the cut by current fracture
     * @param gmeshbackup geometric mesh before performing operations regarding the cut of the fracture
     */
    void RollBack(TPZGeoMesh *gmeshBackup);
    
    /**
     * @brief Search for a 2D element which is looped by this subpolygon based only on element connectivity (no flop).
     * @details If any is found, its youngest children are incorporated to the fracture surface.
     * @param subpolygon [in] vector with oriented closed loop of edges
     * @param polyhindex [in] index of volume that contains the subpolygon
     * @param subpolygMesh [out] elements that cover subpolygon area
     * @returns true: if a subset of the volume shell was incorporated to the fracture surface.
     */
    bool TryFaceIncorporate_Topology(const TPZStack<int64_t>& subpolygon,
                                     const int polyhindex,
                                           TPZStack<int64_t>& subpolygMesh);
    /**
     * @brief Compute internal dihedral angles between a volume shell and the surface mesh of the subpolygon created inside it. And compare it to a threshold.
     * @details If all angles violate the threshold, deletes subpolygMesh and incorporate a subset of the volume shell to the fracture surface
     * @details If not all, but at least one, angles violates the threshold, tags this volume as a badVolume
     * @details If no angle violates the threshold, nothing happens
     * @param subpolygon vector with oriented closed loop of edges
     * @param polyhindex index of volume that contains the subpolygon
     * @param subpolygMesh elements created to cover subpolygon area (will get deleted in the case of an incorporation)
     * @param badVolumes list (to be filled) with polyhedra that need to be refined into simplexes
     * @returns true: if a subset of the volume shell was incorporated to the fracture surface.
     */
    bool TryFaceIncorporate_Geometry(const TPZStack<int64_t>& subpolygon,
                                    const int polyhindex,
                                    const TPZStack<int64_t>& subpolygMesh,
                                        TPZStack<int>& badVolumes);
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

