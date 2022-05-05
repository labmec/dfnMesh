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
#include "DFNPolyhedron.h"
#include "DFNRolodex.h"
#include <gmsh.h>
#include <array>
/**
 * @class DFNMesh
 * @brief Describes a fractured MHM geomesh
 */
class DFNMesh{
private:
    // A set of fractures
    std::vector<DFNFracture *> fFractures;
    
    // Minimum acceptable distance/length for point insertion
    REAL fTolDist = DFN::default_tolDist;
    
    // Minimum acceptable angle
    REAL fTolAngle = DFN::default_tolAngle;
    REAL fTolAngle_cos = DFN::default_tolCos; // == cos(0.2)

    // Pointer to geometric mesh
    TPZGeoMesh *fGMesh = nullptr;

    /// For each edge, a vector of sorted faces around that edge. For efficiency, this structure should live through the life of every DFNMesh object
    TPZVec<TRolodex> fSortedFaces;

    /** 
     * For each 2D skeleton element the left and right polyhedral index
     * - fPolyh_per_face[i][0] is the index of the polyhedron on the positive side of element of index i
     * - fPolyh_per_face[i][1] is the index of the polyhedron on the negative side of element of index i
    */
	TPZVec<std::array<int,2>> fPolyh_per_face;
    
    // A set of polyhedral volumes
    TPZStack<DFNPolyhedron*,20> fPolyhedra;
	
	/// For each dimension (0,1,2,3), it stores the physicalTag (Material ID in PZ) and Name
	TPZManVector<std::map<int,std::string>,4> m_dim_physical_tag_and_name;

public:
    
    /// Constructor
    DFNMesh(TPZGeoMesh *gmesh, TPZManVector<std::map<int,std::string>,4>& dim_physical_tag_and_name, REAL tolerableLength = DFN::default_tolDist, REAL tolerableAngle = DFN::default_tolAngle, int prerefine = 0);

    /// Empty constructor
    DFNMesh(): fGMesh(nullptr){}
    /// Destructor
    ~DFNMesh();
    /// Copy constructor
    DFNMesh(const DFNMesh &copy);
    /// Assignment operator
    DFNMesh &operator=(const DFNMesh &copy);

    /// Pointer to geometric mesh
    TPZGeoMesh *Mesh(){return fGMesh;}

    /// Get a reference to fPolyh_per_face
    /// [for each 2D element an index for each of the polyhedral volume that it encloses]
    TPZVec<std::array<int,2>>& GetPolyhedraPerFace(){return fPolyh_per_face;}

    /// Set tolerances
    void SetTolerances(REAL tolerableLength, REAL tolerableAngle){
        fTolDist = std::max(tolerableLength, gDFN_SmallNumber);
        fTolAngle = std::max(std::min(M_PI,tolerableAngle), gDFN_SmallNumber);
        fTolAngle_cos = std::cos(fTolAngle);
        std::cout<<"\nTolerable distance:  "<<fTolDist;
        std::cout<<"\nTolerable angle:     "<<fTolAngle;
        std::cout<<"\n      cos(angle):    "<<fTolAngle_cos;
        std::cout<<std::endl;
    }

    /// Get minimal tolerable length
    const REAL TolDist(){return fTolDist;}
    /// Get minimal tolerable angle (rad)
    const REAL TolAngle(){return fTolAngle;}
    /// Get cosine of minimal tolerable angle
    const REAL TolAngle_cos(){return fTolAngle_cos;}

    /// Return reference to list of fractures
    std::vector<DFNFracture *>& FractureList(){return fFractures;}
    
    /** 
     * @brief Insert intersection elements of lower dimension in the geometric mesh.
     * @note matid = -1 will use the material id of first neighbour that is found. Set it if you want to force any matID
     */
    void CreateSkeletonElements(int dimension, int matid = -1);

    /** @brief Standard command to create a DFNFracture */
    DFNFracture* CreateFracture(DFNPolygon &Polygon, FracLimit limithandling = Eextended, int materialid = DFNMaterial::Efracture);
    
    
    /// Exports a .geo file for this mesh
    void ExportGMshCAD(std::string filename);
    void ExportGMshCAD_nodes(std::ofstream& out);
    void ExportGMshCAD_edges(std::ofstream& out);
    void ExportGMshCAD_faces(std::ofstream& out);
    void ExportGMshCAD_volumes(std::ofstream& out);
    void ExportGMshCAD_fractures(std::ofstream& out);
    void ExportGMshCAD_boundaryconditions(std::ofstream& out);
    void ExportGMshCAD_fractureIntersections(std::ofstream& out);

    void ExportDetailedGraphics(const std::string ColorPreset = "Rainbow Uniform");
    
    // /// Uses gmsh API to tetrahedralize a DFNVolume
    // void Tetrahedralize(DFNVolume *volume);
    
    // /// Find the volumetrical element that encloses a 2D element
    // bool FindEnclosingVolume(TPZGeoEl *ifracface);
    /**
    //  * @brief Tetrahedralize fractured volumes and refine intact volumes
    //  * @todo refine intact volumes and maybe we could pass a target measure for subelements size as parameter
    //  */
    // void GenerateSubMesh();

    /** @brief Print method for logging */
    void Print(std::ostream& out = std::cout, char* name = nullptr) const;
    void PrintRolodexes(std::ostream& out = std::cout) const;
    void PrintPolyhedra(std::ostream& out = std::cout) const;

    /**
     * @brief Prints Geometric Mesh to a VTK file, a PZGeoMesh.txt or both
     * @param pzmesh : File name for geomesh txt. Feed "skip" to skip
     * @param vtkmesh : File name for geomesh vtk. Feed "skip" to skip
     */
    void PrintVTK(std::string vtkmesh = "LOG/vtkmesh.vtk"
                ,std::string pzmesh = "LOG/pzmesh.txt");
    /**
     * @brief Prints DFN Geometric Mesh and material ids are renumbered for VTK colorful print of refinement of 2D elements :) 
     * @param pzmesh : File name for geomesh txt. Feed "skip" to skip
     * @param vtkmesh : File name for geomesh vtk. Feed "skip" to skip
     */
    void PrintVTKColorful(std::string vtkmesh = "LOG/vtkmesh.vtk"
                        ,std::string pzmesh = "skip");
    
    
    inline int Dimension() const{return fGMesh->Dimension();}
    
    
    
    
    
    
    
    /**
     * @brief Assembles mesh's polyhedral volumes as lists of the faces that outline them, and fills DFNMesh::fPolyh_by_2D_el
     * @attention Assumes boundary polyhedron is already in DFNMesh::fPolyh_by_2D_el
     */
    void UpdatePolyhedra();
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
    void BuildVolume(std::pair<int64_t,int> initial_face_orient, bool& IsConvex, TPZStack<std::pair<int64_t,int>,Talloc>& polyhedron, int& coarseindex);
    /// set a polyhedral index for a face in the structure fPolyh_per_face
    void SetPolyhedralIndex(std::pair<int64_t,int> face_orient, int polyh_index);
    /// get the polyhedral index for a face from the structure fPolyh_per_face
    int GetPolyhedralIndex(const std::pair<int64_t,int>& face_orient) const;
    /// Get the polyhedral index for a face from the structure fPolyh_per_face (alias I usually use with gdb)
    int GetPolyhedralIndex(int64_t faceindex, int orientation) const;
    /**
     * @brief Given a set of faces that enclose a volume, call on Gmsh to generate 3D mesh
    */
    void MeshPolyhedron(const TPZVec<std::pair<int64_t,int>>& polyhedron, int coarseindex, TPZStack<int64_t>& newgels);

    /** @brief Break quadrilaterals down to 2 triangles each in a stack of oriented faces*/
    // template<int Talloc>
    void RefineQuads(const TPZVec<std::pair<int64_t,int>>& polyhedron);

    /**
     * @brief Reference to polyhedra stack
    */
    const TPZStack<DFNPolyhedron*,20>& Polyhedra() const{return fPolyhedra;}
    // template<int Talloc>

    /** @brief adds geoels for graphics that illustrate the tolerance*/
    void PlotTolerance(TPZManVector<int64_t>&);

    // /// Runs DFNMesh::BuildPolyhedra without convexity verification to isolate the boundary
    // void BuildPolyhedra_firstrun();

    /** @brief Initializes the polyhedral volumes datastructure. 
     * @attention It isolates the boundary as the polyhedron of index 0.*/
    void InitializePolyhedra();

    void ExpandPolyhPerFace()
        {fPolyh_per_face.Resize(fGMesh->NElements(),{-1,-1});}

    /// Create a new polyhedron at the end of the polyhedra vector of this mesh
    DFNPolyhedron* CreatePolyhedron(const TPZVec<std::pair<int64_t,int>>& shell,int64_t coarseindex = -1, bool isConvex = true);

    /// return a vector of indices for edges occupying 1D sides of a face
    TPZVec<int64_t> GetEdgeIndices(int64_t face_index){return DFN::GetEdgeIndices(fGMesh->Element(face_index));}

    DFNPolyhedron& Polyhedron(int i){return *fPolyhedra[i];}

    /// Get number of fractures in the fracture vector of this DFN
    int NFractures() const{return fFractures.size();}

    /// Get number of fractures that have at least 1 surface element
    int RealNFractures() const;

    /// Get number of polyhedra for this dfnMesh. Should always be at least 2, unless they aren't initialized
    int NPolyhedra() const {return fPolyhedra.size();}

    /** @brief For a specific face, pass its polyhedral index to its children*/
    void InheritPolyhedra(TPZGeoEl* father);
    
    /// Set all material ids to matid
    void ClearMaterials(int matid = DFNMaterial::Eintact);
    void ClearMaterials(const int matid, TPZVec<int>& backup);
    // Restores material ids from a backup vector
    void RestoreMaterials(TPZVec<int>& backup);

    /// @brief At any point in the code, dump a colored VTK file for graphical debugging 
    /// @warning This method can break the code downstream. It's meant for graphical debugging only.
    /// @note Nothing should break if you leave polygon_representation == false, though;
    /// @param polygon_representation Set true to insert elements that represent the original definition for each fracture. (this will likely break the code downstream and should only be used if you plan to stop the execution soon. It is only meant for graphical debugging);
    /// @param clearmaterials Will set every element's material id in the mesh to 1, before changing again fractures to DFNMaterial::Efracture and also coloring refined faces
    /// @param filename A string with the relative path to the file where to dump the vtk graphics.
    void DumpVTK(bool polygon_representation = false, bool clearmaterials = true, std::string filename = "LOG/vtkmesh.vtk");

    void PrintSummary(){
        std::ofstream out("LOG/dfn.summary.log");
        this->Print(out);
    }

    // DebugStop but dump some more information
    void DFN_DebugStop();

    /// PreRefines the coarse mesh n-times using uniform refpatterns. Can be used to improve mesh quality
    void PreRefine(int n = 1);

private:
    /// Gather oriented faces around a 3D element to define a shell, then create a new polyhedron in the polyhedra vec DFNMesh::fPolyhedra
    int CreateGelPolyhedron(TPZGeoEl* vol, int coarseindex);

    /**
     * @brief Navigate through neighbours of first level, than second level, and so on, until an element of a specific material id is found
     * @returns Index of eldest ancestor of such element
     */
    int64_t SearchIndirectNeighbours(TPZGeoEl* gel);
    /**
     * @brief Pushes all neighbours of a geometric element onto the back of a list
     * @note Not all neighbours are pushed to the list, but rather some criteria are specified so that it gets only those that are candidates of higher interest. Currently these would be 2D elements that are neighbours through the edges of gel.
     * @param gel: Pointer to geometric element
     * @param candidate_queue: Reference to current list of candidates 
     */
    void QueueNeighbours(TPZGeoEl* gel,   std::list<int64_t> &candidate_queue);
    
    /**
     * @brief Check if the neighbour has equal dimension
     * @param geliside GeoElement side
     */
    
    bool HasEqualDimensionNeighbour(TPZGeoElSide &gelside);
    


    /** @brief set polyh index -1 for every face in a stack*/
    void ClearPolyhIndex(TPZVec<std::pair<int64_t,int>>& facestack);
    
    /**
     * @brief prints the rolodex rol and the element with index indexNotFoundCard
     */
    void PrintProblematicRolodex(const int &indexNotFoundCard, TRolodex &rol);
    
    
public:
    /** @brief For every face without a polyh index inherit their father's*/
    void InheritPolyhedra();

    /**
     * @brief Prints the polyhedron associated with element with index coarseindex
     */
    void PlotVolumesByCoarseIndex(const int64_t coarseindex, const std::string dirpath) const;

    
    /**
     * @brief Rollback to the state before the last fracture was introduced
     * @param gmeshBackup geometric mesh before actual fracture cutting process
     * @param badVolumes list with polyhedra volumes that need to be refined into simplexes
     */
    void RollBackLastFracture(TPZGeoMesh *gmeshBackup, TPZStack<int>& badVolumes);

    /** @brief Print and plot detailed info about a Rolodex
     * @param AxleIndex : Index of 1D element which is the axle of the Rolodex
     */
    void PrintRolodexBugReport(const int64_t AxleIndex);

    /// Plot all DFNPolygons to the same file (to plot to different files use DFNPolygon::PlotVTK())
    void PlotAllPolygons(const std::string filepath) const;
};













#endif /* DFNMesh_h */
