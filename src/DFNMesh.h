/*! 
 *  @authors   Pedro Lima
 *  @date      2020.06
 */

#ifndef DFNMesh_h
#define DFNMesh_h

#include "pzgmesh.h"
#include "DFNNamespace.h"
#include "DFNFracture.h"
#include "DFNVolume.h"
#include "DFNRolodex.h"
#include <gmsh.h>
#include <array>
/**
 * @brief Describes a fractured MHM geomesh
 */
class DFNMesh{
private:
    // A set of fractures
    std::list<DFNFracture *> fFractures;
    
    // A set of volumes
    std::map<int64_t, DFNVolume> fVolumes;
    
    // Minimum acceptable distance/length for point insertion
    REAL fTolDist = 1e-5;
    
    // Minimum acceptable angle
    REAL fTolAngle = 0.2;
    REAL fTolAngle_cos = 0.98006657784; // == cos(0.2)

    // Pointer to geometric mesh
    TPZGeoMesh *fGMesh;


    
    
    /// For each edge, a vector of sorted faces around that edge
    TPZVec<TRolodex> fSortedFaces;


    /** 
     * For each 2D skeleton element the left and right polyhedral index
     * - fPolyh_per_2D_el[i][0] is the index of the polyhedron on the positive side of element of index i
     * - fPolyh_per_2D_el[i][1] is the index of the polyhedron on the negative side of element of index i
    */
	TPZVec<std::array<int,2>> fPolyh_per_2D_el;

public:
    // // Local enumeration of material indices
    // enum class Ematerial {intact = 1, fracture = 2, surface = 2, refined = 3, transition = 3};
    
    /// Material tags
    // const int MatIntact = 1;
    // const int MatFracture = 2;
    // const int MatRefined = 3;
    
    /// Constructor
    DFNMesh(TPZGeoMesh *gmesh, REAL tolerableLength = 1e-5, REAL tolerableAngle = 0.2);

    /// Empty constructor
    DFNMesh(): fGMesh(nullptr){fTolAngle_cos = std::cos(fTolAngle);}
    /// Destructor
    ~DFNMesh(){};
    /// Copy constructor
    DFNMesh(const DFNMesh &copy);
    /// Assignment operator
    DFNMesh &operator=(const DFNMesh &copy);

    /// Add new fracture
    void AddFracture(DFNFracture *fracture){fFractures.push_back(fracture);}
    
    /// Pointer to volume of index 'index'
    DFNVolume *Volume(int64_t index){return &fVolumes[index];}
    
    /// Pointer to geometric mesh
    TPZGeoMesh *Mesh(){return fGMesh;}

    /// Set tolerances
    void SetTolerances(REAL tolerableLength, REAL tolerableAngle){
        fTolDist = tolerableLength;
        fTolAngle = tolerableAngle;
        fTolAngle_cos = std::cos(tolerableAngle);
    }

    /// Get minimal tolerable length
    const REAL TolDist(){return fTolDist;}
    /// Get minimal tolerable angle (rad)
    const REAL TolAngle(){return fTolAngle;}
    /// Get cosine of minimal tolerable angle
    const REAL TolAngle_cos(){return fTolAngle_cos;}

    /// Return reference to list of fractures
    std::list<DFNFracture *>& FractureList(){return fFractures;}
    
    /** 
     * @brief Insert intersection elements of lower dimension in the geometric mesh.
     * @note matid = -1 will use the material id of first that is found. Set it if you want to force
     */
    void CreateSkeletonElements(int dimension, int matid = -1);
    
    /// Setup datastructure for fractured volumes (including finding fracture elements enclosed by them)
    void CreateVolumes();
    
    /// Exports a .geo file for this mesh
    void ExportGMshCAD(std::string filename);
    
    /// Uses gmsh API to tetrahedralize a DFNVolume
    void Tetrahedralize(DFNVolume *volume);
    
    /// Find the volumetrical element that encloses a 2D element
    bool FindEnclosingVolume(TPZGeoEl *ifracface);
    /**
     * @brief Tetrahedralize fractured volumes and refine intact volumes
     * @todo refine intact volumes and maybe we could pass a target measure for subelements size as parameter
     */
    void GenerateSubMesh();

    /** @brief Print method for logging */
    void Print(std::ostream& out = std::cout) const;
    void PrintRolodexes(std::ostream& out = std::cout) const;
    void PrintPolyhedra(std::ostream& out = std::cout) const;

    /**
     * @brief Prints DFN Geometric Mesh. 
     * @param pzmesh : File name for geomesh txt. Feed "skip" to skip
     * @param vtkmesh : File name for geomesh vtk. Feed "skip" to skip
     * @param MaterialIDs...
     * @todo this method is unfinished
     */
    void PrintVTK(std::string pzmesh = "pzmesh.txt"
               ,std::string vtkmesh = "vtkmesh.vtk"
               ,int fracture = 2
               ,int transition = 3
               ,int intact = 1);
    /**
     * @brief Prints DFN Geometric Mesh and material ids are renumbered for VTK colorful print of refinement of 2D elements :) 
     * @param pzmesh : File name for geomesh txt. Feed "skip" to skip
     * @param vtkmesh : File name for geomesh vtk. Feed "skip" to skip
     */
    void PrintVTKColorful(std::string pzmesh = "pzmesh.txt"
                       ,std::string vtkmesh = "vtkmesh.vtk");
    
    /**
     * @brief Imports d-dimensional elements from a GMsh::model to a TPZGeoMesh. Imported 
     * elements are pushed to the back of TPZGeoMesh::ElementVector.
     * @note (1) Must be called between the pair gmsh::initialize and gmsh::finalize of the
     * model from which new elements should be read.
     * @note (2) If GMsh has created any new nodes, those will be inserted into TPZGeoMesh aswell
     * @param gmesh: Pointer to TPZGeoMesh
     * @param dimension: Dimension of elements that should be imported
     * @param oldnodes: A set of indices of the old nodes that were used to define the geometry in 
     * GMsh (so that new nodes may be identified)
     */
    void ImportElementsFromGMSH(TPZGeoMesh * gmesh, int dimension, std::set<int64_t> &oldnodes);
    
    int Dimension(){return fGMesh->Dimension();}
    
    
    
    
    
    
    
    /**
     * @brief Sort faces around each 1D element of the mesh
     * @note Fills SortedFaces data
    */
    void SortFacesAroundEdges();
    /**
     * @brief Assembles mesh's polyhedral volumes as lists of the faces that outline them, and fills DFNMesh::fPolyh_by_2D_el
     */
    void BuildPolyhedra();
    /**
     * @brief Identify and match polyh-index of neighbour faces that enclose the same polyhedron
     * @param initial_face_orient Index of geoEl of current face paired with orientation {index,orientation}
     * @param IsConvex Method will flag this as false if a dihedral angle is found to be bigger than 180 degrees (non-convex polyhedron)
     * @param polyhedron: A stack to list every face (and corresponding orientation) that share this polyhedron 
     * @note This method is called by DFNMesh::BuildPolyhedra to assemble all polyhedral volumes in the mesh
    */
    template<int Talloc>
    void BuildVolume(std::pair<int64_t,int> initial_face_orient, bool& IsConvex, TPZStack<std::pair<int64_t,int>,Talloc>& polyhedron);
    /// return a vector of indices for edges occupying 1D sides of a face
    TPZManVector<int64_t,4> GetEdgeIndices(int64_t face_index);
    /// set a polyhedral index for a face in the structure fPolyh_per_2D_el
    void SetPolyhedralIndex(std::pair<int64_t,int> face_orient, int polyh_index);
    /// get the polyhedral index for a face from the structure fPolyh_per_2D_el
    int GetPolyhedralIndex(std::pair<int64_t,int> face_orient);
    /**
     * @brief Given a set of faces that enclose a volume, call on Gmsh to generate 3D mesh
    */
    template<int Talloc>
    void MeshPolyhedron(TPZStack<std::pair<int64_t,int>,Talloc>& polyhedron);

private:
    
    /**
     * @brief Navigate through neighbours of first level, than second level, and so on, until an element of a specific material id is found
     * @returns Index of eldest ancestor of such element
     */
    int64_t SearchIndirectNeighbours(TPZGeoEl* gel);
    /**
     * @brief Goes through all neighbours of gel and identifies if any of them has material id different of that from surface elements
     * @returns Index of eldest ancestor of first found neighbour that matches the material id criteria (macro element)
     * @note Later I might modify this to fill a vector/list with all neighbours that match the material id criteria
     */
    int64_t FindAdjacentMacroEl(TPZGeoEl* gel);
    /**
     * @brief Pushes all neighbours of a geometric element onto the back of a list
     * @note Not all neighbours are pushed to the list, but rather some criteria are specified so that it gets only those that are candidates of higher interest. Currently these would be 2D elements that are neighbours through the edges of gel.
     * @param gel: Pointer to geometric element
     * @param candidate_queue: Reference to current list of candidates 
     */
    void QueueNeighbours(TPZGeoEl* gel,   std::list<int64_t> &candidate_queue);
    
    /**
     * @brief Deletes face + ribs (if left isolated) + nodes (if left isolated)
     * @param face: 2D element to be deleted
     */
    void DeleteElementAndRibs(TPZGeoEl *face);
    /**
     * @brief Deletes gel + children + isolated ribs + unused nodes
     * @note It will assume element has been found not to belong to the domain of interest and will not verify
     * @param gel: pointer to the geometric element
     */
    void CropExternalElement(TPZGeoEl *gel);
    /**
     * @brief Deletes gel and all elements that share the same eldest ancestor, then deletes the ancestor
     * @param gel: Any member of the family
     */
    void DeleteFamily(TPZGeoEl *gel);
    /**
     * @brief Check if the neighbour has equal dimension
     * @param geliside GeoElement side
     */
    
    bool HasEqualDimensionNeighbour(TPZGeoElSide &gelside);
    
    /// Set all material ids to 1
    void ClearMaterials();
    
    // @todo
    void RestoreMaterials(TPZGeoMesh *originalmesh);
};













#endif /* DFNMesh_h */
