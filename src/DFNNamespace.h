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
#define gDFN_SmallNumber DFN::gSmallNumber
/// gmsh doesn't like zero indexed entities
#define gmshshift int(1)

/**
 * @brief Contains methods and variables of the namespace DFN
*/
namespace DFN{

    /// 2*3.1415...
    static const float _2PI = 6.2831853071795865;
    // A small number for geometric tolerances
    static const double gSmallNumber = 1.e-3;

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
    TPZGeoEl* GetSkeletonNeighbour(TPZGeoEl* gel, int side);


    /// builds a loop of oriented 1D elements occupying the 1D sides of a 2D el
    /// @param shift: indices will get shifted by a constant 
    static void GetLineLoop(TPZGeoEl* face_el, std::vector<int>& lineloop, const int shift=0){
        if(face_el->Dimension() != 2) DebugStop();
        int nsides = face_el->NSides();
        TPZManVector<int,4> lineloop_debug(4,-999999999);
        lineloop.resize(face_el->NSides(1));
        for(int iside = face_el->NSides(0); iside<nsides-1; iside++){
            TPZGeoElSide gelside(face_el,iside);
            TPZGeoElSide neig;
            for(neig = gelside.Neighbour(); neig != gelside; neig = neig.Neighbour()){
                if(neig.Element()->Dimension()!=1) continue;
                int orientation = OrientationMatch(gelside,neig)?1:-1;
                lineloop[iside-face_el->NSides(1)] = orientation*(neig.Element()->Index()+shift);
                lineloop_debug[iside-face_el->NSides(1)] = orientation*neig.Element()->Index();
                break;
            }
        }
        return;
    }

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

} /*namespace DFN*/



/// The robustness of these material ids is dependant upon the method DFNMesh::ClearMaterials()
/// Leaving Erefined as the biggest material id is important for better graphics, but isn't mandatory for the code to work. 
/// The choice of ids is, otherwise, arbitrary.
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

template <typename T1, typename T2>
std::ostream& 
operator<<(std::ostream& out, 
           const std::pair<T1,T2>& pair )
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

#include "DFNNamespaceTemp.cpp"
#endif /* DFNNamespace_h */
