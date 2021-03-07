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
	
    /// A material id for this fracture elements
    int fmatid = DFNMaterial::Efracture;

    /// Index of this fracture at fdfnMesh fracture vector
    int fIndex = -1;

    /// Directive for determining if this fracture limits are truncated, extended or recovered
    FracLimit fLimit = Eextended;

	/// Set (of indices) of 2D geo elements on fracture surface
	std::set<int64_t> fSurfaceFaces;
	
	/// Set (of indices) of 1D geo elements on fracture surface. 
    /// @comment Right now, I'm thinking of filling this datastructure only temporarily, to use in DFNFracture::RecoverFractureLimits recovery. This may change later...
	std::set<int64_t> fSurfaceEdges;

public:

#ifdef LOG4CXX
    log4cxx::LoggerPtr fLogger = nullptr;

    /** @brief Creates a logger for this object and fills pointer to this->fLogger
      * @note A method to create a separate logger + appender for each DFNFracture. (as opposed to the main logger which logs everything from the DFNMesh)*/
    log4cxx::LoggerPtr CreateLogger(std::string filename="default", std::string layout_convpattern = "default");
#endif // LOG4CXX


    /// Empty constructor
    DFNFracture();

    ///Destructor
    ~DFNFracture(){};
    
    /**
     * @brief Constructor from a DFNPolygon
     */
    DFNFracture(DFNPolygon &Polygon, DFNMesh *dfnMesh, FracLimit limithandling = Eextended);
    
    /// Copy constructor
    DFNFracture(const DFNFracture &copy);
    
    /// Assignment operator
    DFNFracture &operator=(const DFNFracture &copy);
    
    
    /// Return the corner nodes of the fracture
    DFNPolygon &Polygon();
    
    DFNMesh* dfnMesh() const{return fdfnMesh;}


	std::set<int64_t>& Surface(){return fSurfaceFaces;}

    int MaterialId() const{return fmatid;}
    int Index() const{return fIndex;}

    int NSurfElements() const{return fSurfaceFaces.size();}
    
private:

    void Initialize(DFNPolygon &Polygon, DFNMesh *dfnMesh, FracLimit limithandling = Eextended);
        
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

    /** @brief Search for quadrilaterals that violate our general criterion for refinement*/
    void SearchForSpecialQuadrilaterals();
    /** @brief A special quadrilateral is one that should be refined because of a special exception.
     * @details It has 3 or more neighbour ribs (some of them NOT occupying its 1D sides) whose intersection nodes were snapped toward this quadrilateral.
     * If this quadrilateral is not refined, it'll (sometimes) induce a non-convex polyhedron with overlapped faces that GMSH can't mesh.
    */
    bool IsSpecialQuadrilateral(TPZGeoEl* gel);

public:

    /// @brief Check if there is a common neighbour to 3 geoelsides of dimension dim
    /// @param dim: Filter by dimension. Set -1 to skip filter
    TPZGeoEl* FindCommonNeighbour(TPZGeoElSide& gelside1, TPZGeoElSide& gelside2, TPZGeoElSide& gelside3, int dim = -1);
    /// @brief from a set of 1D elements find if they form a lineloop of an existing 2D element in the mesh
    TPZGeoEl* FindPolygon(TPZStack<int64_t>& polygon);

    /// Triangulates fracture surface from outline
    void MeshFractureSurface();

    /** @brief Following directive DFNFracture::fLimit, modify the extended limits of this fracture to better represent the limits (boundary edges) of the polygon that defined this fracture
     * @todo After debugging, this method should be private and called in the end of DFNFracture::MeshFractureSurface
    */
    void RecoverFractureLimits();

    /// Access the ribs data structure
    DFNRib* AddRib(DFNRib &rib);
    
    DFNFace* AddFace(DFNFace &face);

    /// Insert new volume in data structure
    void AddVolume(DFNVolume volume);

    /// Pointer to rib corresponding to geometric element with index
    // return NULL if the geometric element is not intersected
    DFNRib *Rib(int64_t index);

    /// Pointer to face of index 'index'
    DFNFace *Face(int64_t index);
    
    /** @brief Search for intersected 1D elements and create DFNRib objects for them;
     *  @details Loop over all 1D elements in the mesh; use DFNPolygon::fNodesAbove to check if an edge has nodes on opposite sides of the plane that contains the DFNPolygon; and check if the intersection is within the DFNPolygon bounds.
     **/
    void FindRibs();
    /// Find and intersect ribs within the element indices in a set.
    /// @attention Unlike DFNFracture::FindRibs(), this method won't check if intersection point is within fPolygon limits. Since ribs are already limited to a set, intersection search is done with the unbounded plane that contains fPolygon
    void FindRibs(const std::set<int64_t>& ribset);
    
    /// verify proximity of rib intersection node
    /// Coalesce intersected ribs
    // @TODO change the name of this method
    void SnapIntersections_ribs(REAL tolDist);
    
    /// Set Refinement Patterns and create sub elements
    void RefineRibs();

    /// Find intersected faces
    void FindFaces();
    /// Coalesce intersected faces
    void SnapIntersections_faces(REAL tolDist = 1e-4, REAL tolAngle = 0.1);
    /// Set Refinement Patterns and create sub elements
    void RefineFaces();

    /// Remove refined faces from the surface gelindex set and replace them by their youngest children
    void UpdateFractureSurface();


    /**
     * @name Handling frature limits
     * @{
    */
    /** @brief Identify Ribs, Faces and Polyhedra that are affected by the limits of the fracture
     * @note Should be called after DFNFracture::FindFaces()
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

    /** @brief This is a draft to find the boundary conditions for fractures. I had the idea and quickly wrote down... but have not tested yet.
     * Supposedly, the way it's written, it should allow a user to essentially merge 2 or more fractures into one.
     * Meaning, they would have a continuous boundary 
    */
    void ExportFractureBC(int matid, std::ofstream& out);
};


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

