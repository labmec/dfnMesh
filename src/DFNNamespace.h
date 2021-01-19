/*! 
 *  @authors   Pedro Lima
 *  @date      2020.09
 *  @brief     Namespace DFN declarations
 */

#ifndef DFNNamespace_h
#define DFNNamespace_h


#include "pzgmesh.h"


/// 2*3.1415...
#define gDFN_2PI 6.2831853071795865
#define gDFN_SmallNumber 1.e-4
#define gNoMaterial GMESHNOMATERIAL
#define gDFN_NoIndex -999999
/// gmsh doesn't like zero indexed entities
#define gmshshift int(1)

/**
 * @brief Contains methods and variables of the namespace DFN
*/
namespace DFN{

    /// 2*3.1415...
    static const float _2PI = 6.2831853071795865;
    // A small number for geometric tolerances
    // static const double gSmallNumber = 1.e-4;

    /**
     * @brief Tests if a 2D element is an interface or boundary for 3D coarse elements in the context of DFN meshing
     */
    bool IsInterface(TPZGeoEl* gel);

    template<typename TReturnType,typename Ttype1,typename Ttype2>
    TReturnType DotProduct(TPZManVector<Ttype1,3> &vec1, TPZManVector<Ttype2,3> &vec2);

    template<typename Ttype1, typename Ttype2>
    float DotProduct_f(TPZManVector<Ttype1,3> &vec1, TPZManVector<Ttype2,3> &vec2);

    /** @brief Returns the norm of a vector with template precision*/
    template<typename Ttype>
    Ttype Norm(TPZManVector<Ttype, 3> &vec);

    /** @brief Returns the norm of a vector with float precision*/
    template<typename Ttype>
    float Norm_f(TPZManVector<Ttype, 3> &vec);

    /** 
     * @brief Vector cross product with template return type
     * @param ReturnType CrossProduct<ReturnType>(vec1,vec2)
    */
    template<typename TReturnType, typename T2>
    TPZVec<TReturnType> CrossProduct(TPZManVector<T2,3> &vec1, TPZManVector<T2,3> &vec2);

    template<typename Ttype>
    TPZVec<float> CrossProduct_f(TPZManVector<Ttype,3> &vec1, TPZManVector<Ttype,3> &vec2);

    template <typename T>
    TPZVec<T> operator-(TPZVec<T>& v1,TPZVec<T>& v2);

    template <typename T>
    TPZVec<T> operator+(TPZVec<T>& v1,TPZVec<T>& v2);

    /// @brief Get an orthogonal projection of x onto an arbitrary plane represented by a normal vector and a reference point in the plane
    /// @warning It assumes the given normal vector has norm == 1 and won't verify (to keep it efficient)
    template<typename Ttype>
    TPZManVector<Ttype, 3> GetProjectedX(TPZManVector<Ttype, 3> &x,TPZManVector<Ttype,3>& inplane_point,TPZManVector<Ttype,3>& normal);

    /**
     * @brief Returns the oriented dihedral angle between gel and neighbour
     * @note:1 Make sure neighbour is an actual neighbour, otherwise this method will spit nonsense
     * @note:2 Returned angle is in interval [0, 2pi)
     * @param sideorientation decides if method returns angle or 2*pi-angle, as a right-handed axis 
     * orientation, with the right thumb place over the shared 1D side, and considering the first element node distribution.
     * If thumb orientation matches the orientation of gelside, use 1, else, use -1.
     */
    float DihedralAngle(TPZGeoElSide &gelside, TPZGeoElSide &neighbour, int sideorientation = 1);

    /**
     * @brief Get a vector from node 0 to node 1 of a 1D side
     */
    void GetSideVector(TPZGeoElSide &gelside, TPZManVector<REAL,3>& vector);

    /**
     * @brief Check if the side that connects 2 neighbours has the same orientation in each element
     * @note currently exclusive to 1D sides
     */
    bool OrientationMatch(TPZGeoElSide &neig1, TPZGeoElSide &neig2);

    /**
     * @brief Computes the cossine of the angle at a corner of a 2D element
    */
    REAL CornerAngle_cos(TPZGeoEl *gel, int corner);

    /**
     * @brief Get a pointer to an element that is superposed in a lower dimensional side of a geometric element
     * @return nullptr if there is no element in that side
    */
    TPZGeoEl* GetSkeletonNeighbour(TPZGeoEl* gel, int side);


    /// builds a loop of oriented 1D elements occupying the 1D sides of a 2D el
    /// @param shift: indices will get shifted by a constant 
    template<class Tcontainer>
    void GetLineLoop(TPZGeoEl* face_el, Tcontainer& lineloop, const int shift=0);

    /**
     * @brief Generates the best fitting plane to approximate a point cloud in R3 using least squares
     * @param pointcloud A 3xN matrix of coordinates for N points
     * @param centroid: Reference to a vector to fill with centroid coordinates 
     * @param normal: Reference to a vector to fill with normal vector 
    */
    void BestFitPlane(TPZFMatrix<REAL>& pointcloud, TPZManVector<REAL,3>& centroid, TPZManVector<REAL,3>& normal);

    /** @brief Takes a simple oriented lineloop with 3 or 4 edges and create a geometric element
     * @param lineloop an oriented loop of 3 or 4 edges
    */  
    template<class Tcontainer>
    static TPZGeoEl* MeshSimplePolygon(TPZGeoMesh* gmesh, Tcontainer lineloop, int matid);

    /** @brief Check if a set of mesh nodes are coplanar
     * @warning While using with DFNFracture subpolygons, this method is only robust for sets
     * of 4 points, since it fails in the case of colinearity of nodes 0, 1 and 2. I could make it better
     * by picking the next node-triple whenever this colinearity is found, but there's no need for it
     * right now. I'll leave it as a maybe-to-do...
    */
    template<class TContainer>
    bool AreCoPlanar(TPZGeoMesh* gmesh, TContainer nodeindices, REAL tolerance = gDFN_SmallNumber);

    /** @brief Set material ids for a set of element indices and skip negative entries*/
    template<class TContainer>
    void BulkSetMaterialId(TPZGeoMesh* gmesh, TContainer elementindices, int matid);

    /** @brief Return true if polygon has at least 3 valid edges
     * @note Assumes any negative index represents an edge that was colapsed down to zero dimension
    */
    bool IsValidPolygon(TPZStack<int64_t> polygon);

    /** @brief returns the intersection of 2 sets*/
    template<typename Ttype>
    std::set<Ttype> set_intersection(std::set<Ttype>& set1, std::set<Ttype>& set2);

    /** Outward/inward (true/false) orientation of faces in NeoPZ's tetrahedron */
    static const bool TPZTetrahedron_faceorient[4] = {0,1,1,0};
    /** Outward/inward (true/false) orientation of faces in NeoPZ's hexahedron */
    static const bool TPZCube_faceorient[6] =        {0,1,1,0,0,1};
    /** Outward/inward (true/false) orientation of faces in NeoPZ's prism */
    static const bool TPZPrism_faceorient[5] =       {0,1,1,0,1};
    /** Outward/inward (true/false) orientation of faces in NeoPZ's pyramid */
    static const bool TPZPyramid_faceorient[5] =     {0,1,1,0,0};
    /** @brief Check if a 2D side of a 3D element is oriented outward according to NeoPZ topolgy */
    bool PZOutwardPointingFace(TPZGeoElSide faceside);

    /// return a vector of indices for edges occupying 1D sides of a face
    TPZManVector<int64_t,4> GetEdgeIndices(TPZGeoEl* face);

    /// Creates a refinement pattern using a father element and a list of children, than associates that refpattern to father element
    void CreateRefPattern(TPZGeoEl* father, TPZVec<TPZGeoEl*> children);

    /// Return 1 if subel's normal vector matches its father's, and -1 otherwise
    int SubElOrientation(TPZGeoEl* father, int ichild);

    // Get a vector normal to a 2D element that matches its orientation
    void ElementOrientation(TPZGeoEl* gel, TPZManVector<REAL,3>& orientvec);

    /// Given a volume and a face occupying a side of a volume, check if the face is oriented inwards or outwards the volume
    int SkeletonOrientation(TPZGeoElSide volside, TPZGeoEl* face);

    void SetEdgesMaterialId(TPZGeoEl* gel, int matid);

    // Get the sign of a number
    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }
} /*namespace DFN*/



/// The robustness of these material ids is dependant upon the method DFNMesh::ClearMaterials()
/// Leaving Erefined as the biggest material id is important for better graphics, but isn't mandatory for the code to work. 
/// The choice of ids is, otherwise, arbitrary.
enum DFNMaterial{
    Eintact = 1, 
    Efracture = 2, 
    // Esurface = 2, 
    Erefined = 3, //@todo default used to be 3, I've changed it for a test 
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

/** @brief a better overload of operator<< for std::pair*/
template <typename T1, typename T2>
std::ostream& 
operator<<(      std::ostream&      out, 
           const std::pair<T1,T2>&  pair )
{
    out << pair.first << "|";
    if(std::is_convertible<T2,int>::value){
        out << (pair.second < 0 ? "":" ");
    }
    out << pair.second;
    // out << std::showpos << pair.second;
    // out << std::noshowpos;
	return out;
}

template<typename T>
void Print(const std::set<T>& s, std::ostream& out = std::cout, std::string name = "no-name"){
    out << "\nstd::set : "<<name<<"\n{ ";
    for(auto& el : s){
        out << el << ", ";
    }
    out << "\b\b }\n"<<std::endl;
}

// template <class T>
// int printf(TPZVec<T> vec){
//     for(auto& el : vec){
//         printf(el);
//     }
//     return 1;
// }

// template <typename T1, typename T2>
// std::string to_string(std::pair<T1,T2> pair){
//     std::string buf = "";
//     buf += std::to_string(pair.first);
//     buf += '|';
//     buf += std::to_string(pair.second);
//     return buf;
// }
// template <typename T1, typename T2>
// int printf(std::pair<T1,T2> pair){
//     return std::printf(to_string(pair));
// }

// // TPZGeoElSideIndex &operator= (const TPZGeoElSideIndex &A );
// TPZGeoElSide operator++(TPZGeoElSide& gelside,int){
//     TPZGeoElSide pre = gelside;
//     gelside = gelside.Neighbour();
//     return pre;
// }


// static const std::string RefQuad_2triang_1 = 
// "4 3\n"
// "\n"
// "161601 noname\n"
// "\n"
// "-1 -1  0\n"
// " 1 -1  0\n"
// " 1  1  0\n"
// "-1  1  0\n"
// "\n"
// "3 1 0 1 2 3\n"
// "2 1 0 1 2\n"
// "2 1 0 2 3";
// static const std::string RefQuad_2triang_2 = 
// "4 3\n"
// "\n"
// "161602 noname\n"
// "\n"
// "-1 -1  0\n"
// " 1 -1  0\n"
// " 1  1  0\n"
// "-1  1  0\n"
// "\n"
// "3 1 0 1 2 3\n"
// "2 1 1 2 3\n"
// "2 1 1 3 0";

#include "DFNNamespaceTemp.cpp"
#endif /* DFNNamespace_h */
