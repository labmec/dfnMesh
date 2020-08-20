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
	DFNMesh *fdfnMesh;
	
	/// Map of ribs affected by this fracture
    // @TODO what is the value of the first key?
	std::map<int64_t, DFNRib> fRibs;

	/// Map of faces affected by this fracture
    // @TODO what is the value of the first key?
	std::map<int64_t, DFNFace> fFaces;

	/// A planar convex polygon that defines the outline of the fracture
	DFNPolygon fPolygon;
	
	/// Map of elements on fracture surface
    // @TODO what is the value of the first key?
	std::map<int64_t, TPZGeoEl *> fSurface;

	/// A set of constraints to the surface mesh of the fracture
    // @TODO I have no idea what constraint means? I see a map of geometric elements
    // @TODO what is the value of the first key?
	std::map<int64_t, TPZGeoEl *> fOutline;

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

public:
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

    /// Triangulates fracture surface from outline
    void MeshFractureSurface();
    /// Assemble the set of constraints that outlines the fracture surface
    void AssembleOutline();
    /// Assemble subpolygons (as ordered line loops of 1D members of the Outline)
    void GetSubPolygons();
    /// Builds and fills a list with this fracture outer loop of edges
    void GetOuterLoop(std::vector<int> &outerLoop);
    /// Find faces that should be incorporated to fracture surface
    void GetFacesInSurface(std::vector<TPZGeoEl*> &faces);

    /**
     * @brief Insert elements in the map of elements in surface
    */
    void InsertElementsInSurface(TPZVec<int64_t> &idvec){
        for(auto id : idvec){
            InsertElementsInSurface(fRibs.begin()->second.GeoEl()->Mesh()->Element(id));
        }
    }
    void InsertElementsInSurface(TPZVec<TPZGeoEl*> &elvec){
        for(auto el : elvec){
            InsertElementsInSurface(el);
        }
    }
    void InsertElementsInSurface(TPZGeoEl* el){
        fSurface.insert({el->Index(),el});
    }
};

#endif /* DFNFracture_h */

