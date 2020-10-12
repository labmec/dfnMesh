/*! 
 *  @authors   Pedro Lima
 *  @date      2020.09
 */

#ifndef DFNNamespace_h
#define DFNNamespace_h


#include "pzgmesh.h"



/**
 * @brief Contains methods and variables of the namespace DFN
*/
namespace DFN{
// 2*3.1415...
const float _2PI = 6.2831853071795865;
// A small number for geometric tolerances
static const double gSmallNumber = 1.e-3;

/**
 * @brief Tests if a 2D element is an interface or boundary for 3D coarse elements in the context of DFN meshing
 */
bool IsInterface(TPZGeoEl* gel);







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

/** @brief Returns the norm of a vector with float precision*/
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
    static void BestFitPlane(TPZFMatrix<REAL>& pointcloud, TPZManVector<REAL,3>& centroid, TPZManVector<REAL,3>& normal){
        // Coherence checks
        int npoints = pointcloud.Cols();
        if(npoints < 3) DebugStop();
        if(pointcloud.Rows() != 3) DebugStop();
        
        // compute centroid
        centroid.Resize(3,0.);
        centroid.Fill(0.);
        for(int i=0; i<3; i++){
            for(int j=0; j<npoints; j++){
                centroid[i] += pointcloud(i,j);
            }
            centroid[i] /= npoints;
        }

        // A = pointcloud - centroid
        TPZFMatrix<REAL> A(pointcloud);
        for(int i=0; i<3; i++){
            for(int j=0; j<npoints; j++){
                A(i,j) -= centroid[i];
            }
        }

        // Singular Value Decomposition (SVD): A = U*SIGMA*VT
        #ifdef USING_MKL
        {   // using MKL dgvesd()
            // U (3x3)
            TPZFMatrix<REAL> U(3,3,0.);
            // S (3x1) = Diagonal of Sigma with non zero elements
            TPZFMatrix<REAL> S(3,1,0.);
            // VT (npointsxnpoints)
            TPZFMatrix<REAL> VT(npoints,npoints,0.);
            A.SingularValueDecomposition(U,S,VT,'S','N');
            // Normal vector is the column of U corresponding to the least eigenvalue of A
            normal[0] = U(0,2);
            normal[1] = U(1,2);
            normal[2] = U(2,2);
        }
        #else
            PZError << "\nThis method needs NeoPZ compiled with MKL\n"<<__PRETTY_FUNCTION__;
            DebugStop();
        #endif // USING_MKL
        
        // float norm = Norm_f(normal);
        // if(std::fabs(1.-norm) > DFN::gSmallNumber) DebugStop();
    }
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



#endif /* DFNNamespace_h */
