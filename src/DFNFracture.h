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
	
	/// Map of elements on fracture surface {El_index,El_pointer}
	std::map<int64_t, TPZGeoEl *> fSurface;

public:
    
    /// Empty constructor
    DFNFracture();

    ///Destructor
    ~DFNFracture(){};
    
    /**
     * @brief Constructor from a DFNPolygon
     */
    DFNFracture(DFNPolygon &Polygon, DFNMesh *dfnMesh);
    
    /// Copy constructor
    DFNFracture(const DFNFracture &copy);
    
    /// Assignment operator
    DFNFracture &operator=(const DFNFracture &copy);
    
    
    /// Return the corner nodes of the fracture
    DFNPolygon &Polygon();
    
    DFNMesh* dfnMesh() const{return fdfnMesh;}
    
private:
        
    /// Checks neighbour's dimension and returns true if it is equal
    bool HasEqualDimensionNeighbour(TPZGeoElSide &gelside);
    
    /// Finds intersection point of fracture boundaries and geometric mesh faces
    bool FindEndFracturePoint(DFNFace &face, TPZManVector<REAL,3> &ipoint);
    

    /**
     * @brief Read dim-dimensional geometric elements from a gmsh::model into a TPZGeoMesh, 
     * and imported elements are pushed to the back of TPZGeoMesh::ElementVector 
     * (Must be called between the pair gmsh::initialize and gmsh::finalize of 
     * the model from which elements should be read).
     * @param gmsh: Pointer to geometric mesh where elements should be inserted.
     * @param dimension of elements to be imported
     * @param oldnodes: a set of old nodes that don't require importing
     * @param newelements: a vector with the indices of imported elements
     * @note If GMsh has created any new nodes, those will be inserted into TPZGeoMesh aswell
    */
    void ImportElementsFromGMSH(TPZGeoMesh * gmesh, int dimension, std::set<int64_t> &oldnodes, TPZVec<int64_t> &newelements);

    /**
     * @name Surface Meshing auxiliar methods
     * @{
    */
    /** @brief Projects a non-planar polygon onto its best fitting plane and uses Gmsh to mesh it
     * @param polygon an oriented loop of edges that don't necessarily occupy the same plane
    */
    void MeshPolygon(TPZStack<int64_t>& polygon);
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

    /** @brief Insert elements in the map of elements in surface */
    void InsertElementsInSurface(TPZVec<int64_t> &idvec){
        for(auto id : idvec)
            {InsertElementsInSurface(fRibs.begin()->second.GeoEl()->Mesh()->Element(id));}
    }
    void InsertElementsInSurface(TPZVec<TPZGeoEl*> &elvec){
        for(auto el : elvec)
            {InsertElementsInSurface(el);}
    }
    void InsertElementsInSurface(TPZGeoEl* el){fSurface.insert({el->Index(),el});}
    /// Removes negative integers from a stack
    void ClearNegativeEntries(TPZStack<int64_t>& subpolygon);

    /** @} */

    // Special setup of fracture mat id when working in 2D
    void SetFracMaterial_2D();
public:

    /// @brief Check if there is a common neighbour to 3 geoelsides of dimension dim
    /// @param dim: Filter by dimension. Set -1 to skip filter
    TPZGeoEl* FindCommonNeighbour(TPZGeoElSide& gelside1, TPZGeoElSide& gelside2, TPZGeoElSide& gelside3, int dim = -1);
    /// @brief from a set of 1D elements find if they form a lineloop of an existing 2D element in the mesh
    TPZGeoEl* FindPolygon(TPZStack<int64_t>& polygon);
    /// Triangulates fracture surface from outline
    void MeshFractureSurface();
   /// Access the ribs data structure
    void AddRib(DFNRib &rib);
    
    void AddFace(DFNFace &face);

    /// Insert new volume in data structure
    void AddVolume(DFNVolume volume);

    /// Pointer to rib corresponding to geometric element with index
    // return NULL if the geometric element is not intersected
    DFNRib *Rib(int64_t index);

    /// Pointer to face of index 'index'
    DFNFace *Face(int64_t index);
    
    /// Find intersected ribs, create DFNRib objects
    void FindRibs();
    
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


    /** @brief Print method for logging */
    void Print(std::ostream& out = std::cout) const;

};

#endif /* DFNFracture_h */

