/*! 
 *  @authors   Pedro Lima
 *  @date      2020.06
 */

#ifndef DFNMesh_h
#define DFNMesh_h

#include "pzgmesh.h"
#include "DFNFracture.h"
#include "DFNVolume.h"
#include <gmsh.h>

    // A 2D element sorted around an edge (a rolodex card)
    struct TRolodexCard{
        int64_t fgelindex;
        REAL fangle_to_reference;
        // int forientation; //1 or -1

        /// empty constructor and destructor
        TRolodexCard(int64_t id=-1):fgelindex(-1),fangle_to_reference(0.0){};
        ~TRolodexCard(){};
        /// copy constructor and attribution operator
        TRolodexCard& operator=(const TRolodexCard& copy){
            fangle_to_reference = copy.fangle_to_reference;
            fgelindex = copy.fgelindex;
        }
        TRolodexCard(const TRolodexCard& copy){
            this->operator=(copy);
        }
        /// less than operator
        bool operator<(const TRolodexCard& other){
            return fgelindex < other.fgelindex;
        }
        bool operator==(const TRolodexCard& other){
            return fgelindex == other.fgelindex;
        }
    };
    /// A set of 2D elements around an edge, sorted by an angle to a reference
    // The reference is the first card = fcards.begin()
    struct TRolodex{
        std::vector<TRolodexCard> fcards;
        int64_t fedgeindex;

        /// default constructor and destructor
        TRolodex(): fedgeindex(-1){fcards.resize(0);};
        ~TRolodex(){};
        /// copy constructor and attribution operator
        TRolodex& operator=(const TRolodex& copy){
            fcards = copy.fcards;
            fedgeindex = copy.fedgeindex;
        }
        TRolodex(const TRolodex& copy){
            this->operator=(copy);
        }
        /**
         * @brief Get the next face forward or backwards according to direction and fill the angle between faces
         * @param current_index = index of current element
         * @param angle to be filled with angle between current and next face
         * @param direction 1 or -1
        */
        // TRolodexCard& NextFace(int64_t current_index, REAL& angle, int direction=1){
        //     return fcards[0];
        // }
        TRolodexCard& NextFace(int64_t current_index, REAL& angle, int direction=1){
            TRolodexCard current_card;
            current_card.fgelindex = current_index;
            int ncards = fcards.size();
            auto it = std::find(fcards.begin(),
                                fcards.end(),
                                current_card);
            if(it == fcards.end()) DebugStop();
            int jcard = it - fcards.begin();
            
            if(direction == 1){
                return fcards[(jcard+1)%ncards];
            }
            if(direction == -1){
                return fcards[(jcard-1+ncards)%ncards];
            }
            DebugStop();
        }

        /**
         * @brief Get a reference to a card using an index
         * @param index of the geometric element of that card
         * @param position to fill with the position where this card appears in the card vector
         * @note This methods involve a search (std::find) through a vector
        */
        TRolodexCard& Card(int64_t index, int & position){
            TRolodexCard dummycard;
            dummycard
            auto it = std::find(fcards.begin(),
                                fcards.end(),
                                index);
            if(it == fcards.end()) DebugStop();
            position = it - fcards.begin();
            TRolodexCard& ref = *it;
            return ref;
        }
    };






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
    /**
     * @brief Prints DFN Geometric Mesh. 
     * @param pzmesh : File name for geomesh txt. Feed "skip" to skip
     * @param vtkmesh : File name for geomesh vtk. Feed "skip" to skip
     * @param MaterialIDs...
     * @todo this method is unfinished
     */
    void Print(std::string pzmesh = "pzmesh.txt"
               ,std::string vtkmesh = "vtkmesh.vtk"
               ,int fracture = 2
               ,int transition = 3
               ,int intact = 1);
    /**
     * @brief Prints DFN Geometric Mesh and material ids are renumbered for VTK colorful print of refinement of 2D elements :) 
     * @param pzmesh : File name for geomesh txt. Feed "skip" to skip
     * @param vtkmesh : File name for geomesh vtk. Feed "skip" to skip
     */
    void PrintColorful(std::string pzmesh = "pzmesh.txt"
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
    
    
    
    
    
    
    //@todo: Talk to Phil about this FaceTracker usage... I'm not sure if this is what he had in mind
    // @phil If you dont know what its for, then it is useless
    
    //Maps how many polyhedra this face is part of (at most 2)
    // @TODO Since when a class variable does not begin with "f"?
    std::map<int64_t,int> FaceTracker;
    // @TODO why are these method inline?
    // @TODO documentation of this method??
    void InitializeFaceTracker(){
        if(this->Dimension()<3) return;
        if(this->FaceTracker.size() > 0){std::cout<<"\n\n"<<__PRETTY_FUNCTION__<<"\nTried to initialize an already initialized structure\n\n"; return;}
        TPZGeoEl* gel = nullptr;
        int64_t nels = this->Mesh()->NElements();
        for(int64_t iel=0; iel<nels; iel++){
            gel = fGMesh->Element(iel);
            if(!gel) continue;
            if(gel->Dimension() != 2) continue;
            int n3Dneighbours = 0;
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            TPZGeoElSide neig = gelside.Neighbour();
            while(neig!=gelside){
                if(neig.Element()->Dimension()==3)   n3Dneighbours++;
                neig = neig.Neighbour();
            }
            if(n3Dneighbours > 2) DebugStop(); // A face cannot be shared by more than 2 volumes in 3D
            FaceTracker[iel] = n3Dneighbours;
        }
    }
    // @TODO No idea what this means
    void UpdateFaceTracker(){
        for(TPZGeoEl* el : fGMesh->ElementVec()){
            if(!el) continue;
            if(el->Dimension() != 2) continue;
            TPZGeoEl* elder = nullptr;
            if(el->Father()){
                elder = el->EldestAncestor();
            }else{
                elder = el;
            }
            auto end = FaceTracker.end();
            auto eldertrack = FaceTracker.find(elder->Index());
            if( eldertrack == end){
                FaceTracker[el->Index()] = 2; //fracture surface element
            }else{
                FaceTracker[el->Index()] = eldertrack->second; //child inherits father npolyh
            }
        }
    }
    
    
    
    
    
    /**
     * @brief Assembles mesh's polyhedral volumes as lists of the faces that outline them
     */
    void GetPolyhedra();
    void GetPolyhedra2();

    /**
     * @brief Sort faces around each 1D element of the mesh
     * @note Fills SortedFaces data
    */
    void SortFacesAroundEdges();

private:
    // /**
    //  *  @brief Navigate children tree to access most extreme branches
    //  *  @param gel: Pointer to geometric element of ancestor
    //  *  @param outfile: ofstream in which to write accessed data
    //  */
    // void PrintYoungestChildren(TPZGeoEl *gel, std::ofstream &outfile);
    
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


namespace DFN{
/**
 * @brief Tests if a 2D element is an interface or boundary for 3D coarse elements in the context of DFN meshing
 */
bool IsInterface(TPZGeoEl* gel);
}

// @TODO there should be lengthy explanation here documenting the why and how
// of these material ids
enum DFNMaterial{
    Eintact = 1, 
    Efracture = 2, 
    // Esurface = 2, 
    Erefined = 3, 
    // Etransition = 3
};


// Set Material ID for element and its children
static void SetMaterialIDChildren(int id, TPZGeoEl* gel){
    gel->SetMaterialId(id);
    if(gel->HasSubElement()){
        int nchildren = gel->NSubElements();
        for(int i=0; i<nchildren; i++){
            SetMaterialIDChildren(id,gel->SubElement(i));
        }
    }
}








namespace DFN{
// 2*3.1415...
const float _2PI = 6.2831853071795865;
// A small number for geometric tolerances
static const double gSmallNumber = 1.e-3;



template<typename Ttype>
float DotProduct_f(TPZManVector<Ttype,3> &vec1, TPZManVector<Ttype,3> &vec2){
    int size1 = vec1.size();
    int size2 = vec2.size();
    if(size1 != size2){throw std::bad_exception();}
    float dot = 0.;
    for(int j=0; j<size1; j++){
        dot += vec1[j]*vec2[j];
    }
    return dot;
}

template<typename Ttype>
float Norm_f(TPZManVector<Ttype, 3> &vec){
    float norm = 0.;
    for(int j=0, size=vec.size(); j<size; j++){
        norm += vec[j]*vec[j];
    }
    return std::sqrt(norm);
}

template<typename Ttype>
TPZManVector<Ttype,3> CrossProduct_f(TPZManVector<Ttype,3> &vec1, TPZManVector<Ttype,3> &vec2){
    if(vec1.size() != 3){throw std::bad_exception();}
    if(vec2.size() != 3){throw std::bad_exception();}
    
    TPZManVector<REAL,3> result(3,0.);
    result[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    result[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    result[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
    return result;
}
template <class T, int NumExtAlloc1, int NumExtAlloc2>
TPZManVector<T,3> operator-(TPZManVector<T,NumExtAlloc1>& v1,TPZManVector<T,NumExtAlloc2>& v2){
    int64_t size1 = v1.size();
    int64_t size2 = v2.size();
    if(size1 != size2) throw std::bad_exception();
    TPZManVector<T,3> result(size1);
    for(int64_t i = 0; i<size1; i++){
        result[i] = v1[i] - v2[i];
    }
    return result;
}


/**
 * @brief Returns the oriented dihedral angle between gel and neighbour
 * @note:1 Make sure neighbour is an actual neighbour, otherwise this method will spit nonsense
 * @note:2 Returned angle is in interval [0, 2pi)
 * @param sideorientation decides if method returns angle or 2*pi-angle, as a right-handed axis 
 * orientation, with the right thumb place over the shared 1D side, and considering the first element node distribution.
 * If thumb orientation matches the orientation of gelside, use 1, else, use -1.
 */
static float DihedralAngle(TPZGeoElSide &gelside, TPZGeoElSide &neighbour, int sideorientation = 1){

    // Consistency checks
    if(gelside.Element()->Dimension() != 2)     DebugStop();
    if(gelside.Dimension() != 1)                DebugStop();
    if(neighbour.Element()->Dimension() !=2)    DebugStop();
    if(neighbour.Dimension() != 1)              DebugStop();

    if(gelside == neighbour) return 0.;
    // {
    //     switch(sideorientation){
    //         case  1: return 0.;
    //         case -1: return DFN::_2PI;
    //     }
    // }

    if(!gelside.NeighbourExists(neighbour))     DebugStop();
    
    TPZGeoEl* gel = gelside.Element();
    TPZGeoMesh* gmesh = gel->Mesh();
    const int side = gelside.Side();
    TPZManVector<double,3> sharednode0(3,0);
    TPZManVector<double,3> sharednode1(3,0);
    gmesh->NodeVec()[gelside.SideNodeIndex(0)].GetCoordinates(sharednode0);
    gmesh->NodeVec()[gelside.SideNodeIndex(1)].GetCoordinates(sharednode1);
    
    TPZManVector<REAL,3> oppositenode_gel(3,0);
    TPZManVector<REAL,3> oppositenode_neig(3,0);
    gel->Node((gelside.Side()+2)%gel->NNodes()).GetCoordinates(oppositenode_gel);
    neighbour.Element()->Node((neighbour.Side()+2)%neighbour.Element()->NNodes()).GetCoordinates(oppositenode_neig);
    TPZManVector<REAL,3> tangentvec_gel = oppositenode_gel - sharednode0;
    TPZManVector<REAL,3> tangentvec_neig = oppositenode_neig - sharednode0;
    TPZManVector<REAL,3> tangentvec_edge(3);
    switch(sideorientation){
        case -1:{tangentvec_edge = sharednode1 - sharednode0; break;}
        case  1:{tangentvec_edge = sharednode0 - sharednode1; break;}
        default: DebugStop();
    }
    
    TPZManVector<REAL,3> normalvec_gel = CrossProduct_f(tangentvec_gel,tangentvec_edge);
    TPZManVector<REAL,3> normalvec_neig = CrossProduct_f(tangentvec_neig,tangentvec_edge);;
    TPZManVector<REAL,3> aux = CrossProduct_f(normalvec_neig,normalvec_gel);
    float x = Norm_f(tangentvec_edge)*DotProduct_f(normalvec_neig,normalvec_gel);
    float y = DotProduct_f(tangentvec_edge,aux);
    float angle = atan2(y,x);
    
    return (angle >= 0? angle : angle + _2PI);
}

/**
 * @brief Get a vector from node 0 to node 1 of a 1D side
 */
static void GetSideVector(TPZGeoElSide &gelside, TPZManVector<REAL,3>& vector){
    if(gelside.Dimension() != 1) DebugStop();
    int node0 = gelside.SideNodeLocIndex(0);
    int node1 = gelside.SideNodeLocIndex(1);
    
    TPZManVector<REAL,3> coord0(3,0);
    TPZManVector<REAL,3> coord1(3,0);
    gelside.Element()->Node(node0).GetCoordinates(coord0);
    gelside.Element()->Node(node1).GetCoordinates(coord1);
    
    vector = coord1 - coord0;
}

/**
 * @brief Check if the side that connects 2 neighbours has the same orientation in each element
 * @note currently exclusive to 1D sides
 */
static bool OrientationMatch(TPZGeoElSide &neig1, TPZGeoElSide &neig2){
    if(neig1.Dimension() != 1) DebugStop();
    if(!neig1.NeighbourExists(neig2)) DebugStop();
    return (neig1.SideNodeIndex(0) == neig2.SideNodeIndex(0));
}

/**
 * @brief Computes the cossine of the angle at a corner of a 2D element
*/
static REAL CornerAngle_cos(TPZGeoEl *gel, int corner){
	int ncorners = gel->NCornerNodes();
	if(corner >= ncorners) DebugStop();

	int nsides = gel->NSides();

	TPZManVector<REAL,3> point_corner(3);
	TPZManVector<REAL,3> point_anterior(3);
	TPZManVector<REAL,3> point_posterior(3);

	gel->Node(corner).GetCoordinates(point_corner);
	gel->Node((corner+1)%ncorners).GetCoordinates(point_posterior);
	gel->Node((corner-1+ncorners)%ncorners).GetCoordinates(point_anterior);

	TPZManVector<REAL,3> vec1 = point_posterior - point_corner;
	TPZManVector<REAL,3> vec2 = point_anterior  - point_corner;

	REAL cosine = DotProduct_f(vec1,vec2)/(Norm_f(vec1)*Norm_f(vec2));
	return cosine;
}

    /**
     * @brief Get a pointer to an element that is superposed in a lower dimensional side of a geometric element
     * @return nullptr if there is no element in that side
    */
    static TPZGeoEl* GetSkeletonNeighbour(TPZGeoEl* gel, int side){
        if(gel->SideDimension(side) == gel->Dimension()) return nullptr;
        TPZGeoElSide gelside(gel,side);
        TPZGeoElSide neig;
        for(neig = gelside.Neighbour(); neig!=gelside; neig=neig.Neighbour()){
            if(neig.Element()->Dimension() == gel->SideDimension(side)){
                return neig.Element();
            }
        }
        return nullptr;
    }
}

#endif /* DFNMesh_h */
