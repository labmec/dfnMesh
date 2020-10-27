/*! 
 *  @authors   Pedro Lima
 *  @date      2020.10
 */



#include "DFNNamespace.h"



namespace DFN{
    

    /**
     * @brief Get a pointer to an element that is superposed in a lower dimensional side of a geometric element
     * @return nullptr if there is no element in that side
    */
    TPZGeoEl* GetSkeletonNeighbour(TPZGeoEl* gel, int side){
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


    /**
     * @brief Generates the best fitting plane to approximate a point cloud in R3 using least squares
     * @param pointcloud A 3xN matrix of coordinates for N points
     * @param centroid: Reference to a vector to fill with centroid coordinates 
     * @param normal: Reference to a vector to fill with normal vector 
    */
    void BestFitPlane(TPZFMatrix<REAL>& pointcloud, TPZManVector<REAL,3>& centroid, TPZManVector<REAL,3>& normal){
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
        {
            // using MKL dgvesd()
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

        #ifdef PZDEBUG
        {
            float norm = Norm_f(normal);
            if(std::fabs(1.-norm) > gDFN_SmallNumber) DebugStop();
        }
        #endif //PZDEBUG
    }


    /**
     * @brief Returns the oriented dihedral angle between gel and neighbour
     * @note:1 Make sure neighbour is an actual neighbour, otherwise this method will spit nonsense
     * @note:2 Returned angle is in interval [0, 2pi)
     * @param sideorientation decides if method returns angle or 2*pi-angle, as a right-handed axis 
     * orientation, with the right thumb place over the shared 1D side, and considering the first element node distribution.
     * If thumb orientation matches the orientation of gelside, use 1, else, use -1.
     */
    float DihedralAngle(TPZGeoElSide &gelside, TPZGeoElSide &neighbour, int sideorientation){

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
        
        TPZManVector<float,3> normalvec_gel = CrossProduct<float>(tangentvec_gel,tangentvec_edge);
        TPZManVector<float,3> normalvec_neig = CrossProduct<float>(tangentvec_neig,tangentvec_edge);;
        TPZManVector<float,3> aux = CrossProduct<float>(normalvec_neig,normalvec_gel);
        float x = Norm_f(tangentvec_edge)*DotProduct_f(normalvec_neig,normalvec_gel);
        float y = DotProduct_f(tangentvec_edge,aux);
        float angle = atan2(y,x);
        
        return (angle >= 0? angle : angle + _2PI);
    }

    /** @brief Return true if polygon has at least 3 valid edges
     * @note Assumes any negative index represents an edge that was colapsed down to zero dimension
    */
    bool IsValidPolygon(TPZStack<int64_t> polygon){
        int nedges = 0;
        for(auto& edge : polygon){
            nedges += int(edge>-1);
        }
        return nedges > 2;
    }

    /** @brief Alternate name of dfn::IsValidPolygon for readability
     * @returns False
    */
    bool IsDegeneratePolygon(TPZStack<int64_t> polygon){
        return !IsValidPolygon(polygon);
    }

} /* namespace DFN*/