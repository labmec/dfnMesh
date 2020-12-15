/*! 
 *  @authors   Pedro Lima
 *  @date      2020.10
 */



#include "DFNNamespace.h"
#include "TPZRefPattern.h"



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

    /** @brief Check if a 2D side of a 3D element is oriented outward according to NeoPZ topolgy */
    bool PZOutwardPointingFace(TPZGeoElSide faceside){
        // consistency check
        TPZGeoEl* gel = faceside.Element();
        if(faceside.Dimension() != 2 || gel->Dimension() < 3) 
            {PZError<<__PRETTY_FUNCTION__<<" Invalid_Argument\n"; DebugStop();}
        
        int nlowdimsides = gel->NSides(0) + gel->NSides(1);
        int iface = faceside.Side() - nlowdimsides;

        switch(gel->Type()){
            case ETetraedro: return TPZTetrahedron_faceorient[iface];
            case ECube:      return TPZCube_faceorient[iface];
            case EPrisma:    return TPZPrism_faceorient[iface];
            case EPiramide:  return TPZPyramid_faceorient[iface];
            default: DebugStop();
        }

        DebugStop();
        return false;
    }

    TPZManVector<int64_t,4> GetEdgeIndices(TPZGeoEl* face){
        // consistency checks
        if(!face) 							DebugStop();
        if(face->Dimension() != 2) 			DebugStop();
        // if(face->HasSubElement()) 			DebugStop();

        int nedges = face->NSides(1);
        TPZManVector<int64_t,4> output(nedges,-1);
        for(int iside = nedges; iside<face->NSides()-1; iside++){
            output[iside-nedges] = DFN::GetSkeletonNeighbour(face,iside)->Index();
        }

        return output;
    }

    /**
     * @brief Computes the cossine of the angle at a corner of a 2D element
    */
    REAL CornerAngle_cos(TPZGeoEl *gel, int corner){
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
     * @brief Check if the side that connects 2 neighbours has the same orientation in each element
     * @note currently exclusive to 1D sides
     */
    bool OrientationMatch(TPZGeoElSide &neig1, TPZGeoElSide &neig2){
        if(neig1.Dimension() != 1) DebugStop();
        if(!neig1.NeighbourExists(neig2)) DebugStop();
        return (neig1.SideNodeIndex(0) == neig2.SideNodeIndex(0));
    }



    /**
     * @brief Get a vector from node 0 to node 1 of a 1D side
     */
    void GetSideVector(TPZGeoElSide &gelside, TPZManVector<REAL,3>& vector){
        if(gelside.Dimension() != 1) DebugStop();
        int node0 = gelside.SideNodeLocIndex(0);
        int node1 = gelside.SideNodeLocIndex(1);
        
        TPZManVector<REAL,3> coord0(3,0);
        TPZManVector<REAL,3> coord1(3,0);
        gelside.Element()->Node(node0).GetCoordinates(coord0);
        gelside.Element()->Node(node1).GetCoordinates(coord1);
        
        vector = coord1 - coord0;
    }






    void CreateRefPattern(TPZGeoEl* father, TPZVec<TPZGeoEl*> children){
        TPZGeoMesh refmesh;
        TPZGeoMesh* gmesh;
        int nfathernodes = father->NCornerNodes();
        TPZManVector<int64_t,4> nodeindices(nfathernodes,-1);
        TPZManVector<REAL,3> coord(3,0.);
        std::map<int64_t,int64_t> orig_to_copy;
        // std::map<int64_t,int64_t> copy_to_orig;
        // Copy nodes
		for(int i=0; i<nfathernodes; i++){
            nodeindices[i] = refmesh.NodeVec().AllocateNewElement();
			father->Node(i).GetCoordinates(coord);
			refmesh.NodeVec()[i].Initialize(coord,refmesh);

            // map
            orig_to_copy[father->NodeIndex(i)] = i;
            // copy_to_orig[i] = father->NodeIndex(i);
		}
		// Copy father
		MElementType etype = father->Type();
		int64_t index = -1;
		refmesh.CreateGeoElement(etype,nodeindices,1,index);
        
        // Copy children
        int nchildren = children.size();
        for(int ichild=0; ichild < nchildren; ichild++){
            TPZGeoEl* child = children[ichild];
            int nchildnodes = child->NCornerNodes();
            nodeindices.resize(nchildnodes);
            // gather children nodes that may be other than those of their father
            for(int inode=0; inode < nchildnodes; inode++){
                auto itr = orig_to_copy.find(child->NodeIndex(inode));
                int64_t copyindex;
                if(itr == orig_to_copy.end()){
                    copyindex = orig_to_copy.insert({child->NodeIndex(inode),refmesh.NodeVec().AllocateNewElement()}).second;
                    child->Node(copyindex).GetCoordinates(coord);
			        refmesh.NodeVec()[copyindex].Initialize(coord,refmesh);
                }else{
                    copyindex = itr->second;
                }
                
                nodeindices[inode] = copyindex;
            }
            etype = child->Type();
            index = -1;
            refmesh.CreateGeoElement(etype,nodeindices,1,index);
        }

        // create refpattern from refmesh
        TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(refmesh);
	    father->SetRefPattern(refpat);
    }

    void ElementOrientation(TPZGeoEl* gel, TPZManVector<REAL,3>& orientvec){
        // @todo This could be done using gram-schimidt orthogonalization
        // gel->Jacobian()
        switch(gel->Dimension()){
            case  1:{
                TPZGeoElSide gelside(gel,2);
                DFN::GetSideVector(gelside,orientvec); break;
            }
            case  2:{
                TPZGeoElSide edgeside1(gel,gel->FirstSide(1));
                TPZGeoElSide edgeside2(gel,gel->FirstSide(1)+1);
                TPZManVector<REAL,3> edgevec1(3,0.);
                TPZManVector<REAL,3> edgevec2(3,0.);
                DFN::GetSideVector(edgeside1,edgevec1);
                DFN::GetSideVector(edgeside2,edgevec2);
                orientvec = CrossProduct<REAL>(edgevec1,edgevec2);
                REAL norm = Norm<REAL>(orientvec);
                orientvec[0] /= norm;
                orientvec[1] /= norm;
                orientvec[2] /= norm;
                break;
            }
            default: PZError << "\nTried to get orientation vector of element of dimension ="<<gel->Dimension()<<"\n";DebugStop();
        }
    }

    int SubElOrientation(TPZGeoEl* father, int ichild){
        if(!father->HasSubElement()) DebugStop();

        TPZGeoEl* child = father->SubElement(ichild);
        if(father->Dimension() != child->Dimension()){PZError<<"\n I don't see how this could be considered consistent\n";DebugStop();}

        TPZGeoElSide fatherside(father,father->NSides()-1);
        TPZGeoElSide childside(child,child->NSides()-1);
        TPZManVector<REAL,3> qsi(fatherside.Dimension(),0.);
        int orientation=0;
        switch(father->Dimension()){
            case  2:{
                // TPZManVector<REAL,3> father_normal(3,0.);
                // fatherside.Normal(qsi,father_normal);
                // TPZManVector<REAL,3> child_normal(3,0.);
                // childside.Normal(qsi,child_normal);
                TPZManVector<REAL,3> father_normal(3,0.);
                TPZManVector<REAL,3> child_normal(3,0.);
                ElementOrientation(father,father_normal);
                ElementOrientation(child,child_normal);
                orientation = int(DFN::DotProduct<float>(father_normal,child_normal));
                break;
            }
            default:{PZError<<"\nCurrently implemented for 2D els only\n";DebugStop();}
        }

        if(orientation != 1 && orientation != -1) DebugStop();
        return orientation;
    }



    int SkeletonOrientation(TPZGeoElSide volside, TPZGeoEl* face){
        TPZGeoElSide faceside(face,face->NSides()-1);
        // Get transformation from volume to skeleton face
        TPZTransform<REAL> transf = volside.NeighbourSideTransform(faceside);
        // Compute matrix determinant
        TPZFMatrix<REAL> inverse;
        REAL det = 0.;
        transf.Mult().DeterminantInverse(det,inverse);
        // positive determinant means skeleton's orientation matches volume's faceside orientation
        int neig_side_orientation = sgn(det);
        // get 2D side orientation
        int sideorientation = (PZOutwardPointingFace(volside)?1:-1);
        // return their product
        return sideorientation * neig_side_orientation;
    }

} /* namespace DFN*/