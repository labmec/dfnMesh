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
#include "DFNFracPlane.h"
#include "DFNMesh.h"
#include <gmsh.h>

class DFNMesh;
typedef TPZFMatrix<REAL> Matrix;

/** 
 *  @brief     Describes a surface mesh for a fracture and all ribs & faces that are intersected by it.
 *  @details   Intersection search is performed after creation of skeleton
 *  elements with DFNFracture::CreateSkeletonElements. Fracture plane should
 *  be a DFNFracPlane.
 */
class DFNFracture
{
private:
	
	/// Pointer for the complete DFN mesh
	DFNMesh *fdfnMesh;
	
	/// Map of ribs affected by this fracture
	std::map<int64_t, DFNRib> fRibs;

	/// Map of faces affected by this fracture
	std::map<int64_t, DFNFace> fFaces;

	/// A planar convex polygon that indicates an insertion region for a fracture
	DFNFracPlane *fFracplane;
	
	/// Map of elements on fracture surface
	std::map<int64_t, TPZGeoEl *> fSurface;

	/// A set of constraints to the surface mesh of the fracture
	std::map<int64_t, TPZGeoEl *> fOutline;

public:
    
    /// Empty constructor
    DFNFracture();
    
    /**
     * Define the fracture plane from 3 to 4 points
     * Points should be coplanar
     * The matrix should be dimension 3xN, each column defining the coordinates
     * of a point
     *  
     */
    DFNFracture(DFNFracPlane *FracPlane, DFNMesh *dfnMesh);
    
    /// Copy constructor
    DFNFracture(const DFNFracture &copy);
    
    /// Assignment operator
    DFNFracture &operator=(const DFNFracture &copy);
    
    
    /// Return the corner nodes of the fracture
    DFNFracPlane *FracPlane();
    
    
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
     * @note If GMsh has created any new nodes, those will be inserted into TPZGeoMesh aswell
    */
    void ImportElementsFromGMSH(TPZGeoMesh * gmesh, int dimension);

public:
    /// Access the ribs data structure
    void AddRib(DFNRib &rib);
    
    void AddFace(DFNFace &face);

    /// Insert new volume in data structure
    void AddVolume(DFNVolume volume);

    /// Pointer to rib of index 'index'
    DFNRib *Rib(int64_t index);

    /// Pointer to face of index 'index'
    DFNFace *Face(int64_t index);
    
    /// Find intersected ribs
    void FindRibs();
    /// Coalesce intersected ribs
    void OptimizeRibs(REAL tolDist);
    /// Set Refinement Patterns and create sub elements
    void RefineRibs();

    /// Find intersected faces
    void FindFaces();
    /// Coalesce intersected faces
    void OptimizeFaces(REAL tolDist = 1e-4, REAL tolAngle = 0.1);
    /// Set Refinement Patterns and create sub elements
    void RefineFaces();

    /// Triangulates fracture surface from outline
    void MeshFractureSurface();
    /// Assemble the set of constraints that outlines the fracture surface
    void AssembleOutline();
    /// Builds and fills a list with this fracture outer loop of edges
    void GetOuterLoop(std::vector<int> &outerLoop);
    /// Find faces that should be incorporated to fracture surface
    void GetFacesInSurface(std::vector<TPZGeoEl*> &faces);
};

#endif /* DFNFracture_h */

