/*! 
 *  @authors   Pedro Lima
 *  @date      2020.10
 *  @brief     Namespace DFN definitions
 */



#include "DFNNamespace.h"
#include "TPZRefPattern.h"

#if PZ_LOG
static TPZLogger logger("dfn.mesh");
#endif


namespace DFN{
    

    void PrintGeoEl(TPZGeoEl* gel) {

        if (!gel) {
            DebugStop();
        }
        const int index = gel->Index();
        std::string filename = "gel_index_" + std::to_string(index) + ".vtk";
        std::cout << "\n =====> Printing Element to " << filename << std::endl;
        TPZGeoMesh *gmesh = gel->Mesh();
        
        TPZGeoMesh bogusMesh;
        bogusMesh.ElementVec().Resize(gmesh->NElements());
        for (int i = 0; i < bogusMesh.NElements(); i++) {
            bogusMesh.ElementVec()[i] = nullptr;
        }
        bogusMesh.NodeVec() = gmesh->NodeVec();
        
        // Adding el to print to bogus mesh
        TPZGeoEl* copiedEl = gel->Clone(bogusMesh);
        copiedEl->SetMaterialId(2);
        
        std::ofstream out(filename);
        TPZVTKGeoMesh::PrintGMeshVTK(&bogusMesh, out);
    }

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
#ifdef PZ_USING_LAPACK
        {
            // using LAPACK dgvesd()
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
#else // PZ_USING_LAPACK
            PZError << "\nThis method needs NeoPZ compiled with LAPACK\n"<<__PRETTY_FUNCTION__;
            DebugStop();
#endif // PZ_USING_LAPACK

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
    REAL DihedralAngle(TPZGeoElSide &gelside, TPZGeoElSide &neighbour, int sideorientation){

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
     * @note 1: Assumes any negative index represents an edge that was colapsed down to zero dimension
     * @note 2: Removes all negative entries from the polygon
    */
    bool IsValidPolygon(TPZStack<int64_t>& polygon,TPZGeoMesh* gmesh){
        // clear collapsed edges from polygon lineloop
        DFN::ClearNegativeEntries(polygon);
        int nedges = polygon.size();
        if(nedges <= 1) return false;
        // check for co-linearity of edges through shared eldest ancestor
        std::set<int64_t> ancestors;
        for(const auto index : polygon){
            TPZGeoEl* elder = gmesh->Element(index)->EldestAncestor();
            ancestors.insert(elder->Index());
        }
        return ancestors.size() > 1;
    }

    /** @brief Alternate name of dfn::IsValidPolygon for readability
     * @returns False
    */
    bool IsDegeneratePolygon(TPZStack<int64_t>& polygon, TPZGeoMesh* gmesh){
        return !IsValidPolygon(polygon,gmesh);
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






    void CreateRefPattern(TPZGeoEl* father, TPZVec<TPZGeoEl*>& children){
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
			refmesh.NodeVec()[nodeindices[i]].Initialize(coord,refmesh);

            // map
            orig_to_copy[father->NodeIndex(i)] = nodeindices[i];
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
                    copyindex = refmesh.NodeVec().AllocateNewElement();
                    orig_to_copy.insert({child->NodeIndex(inode),copyindex});
                    child->Node(inode).GetCoordinates(coord);
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
        // std::stringstream stream;
        // refmesh.Print(stream);
        // LOG4CXX_DEBUG(logger,stream.str());

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

    /// @todo Rewrite this using the determinant of the matrix of the transformation between parametric space of father-child, so it won't be limited to 2D elements
    int SubElOrientation(TPZGeoEl* father, int ichild){
        if(!father->HasSubElement()) DebugStop();

        TPZGeoEl* child = father->SubElement(ichild);
        int dim = father->Dimension();
        if(dim != child->Dimension()){PZError<<"\n I don't see how this could be considered consistent\n";DebugStop();}
        
        

        TPZGeoElSide fatherside(father,father->NSides()-1);
        TPZGeoElSide childside(child,child->NSides()-1);
        TPZManVector<REAL,3> qsi(fatherside.Dimension(),0.);
        int orientation=0;
        switch(dim){
            case  2:{
                // TPZManVector<REAL,3> father_normal(3,0.);
                // TPZManVector<REAL,3> child_normal(3,0.);
                // ElementOrientation(father,father_normal);
                // ElementOrientation(child,child_normal);
                // orientation = int(DFN::DotProduct<float>(father_normal,child_normal));
                TPZTransform<REAL> transf = father->GetTransform(child->NSides()-1,ichild);
                REAL det = 1.;
                // TPZFMatrix<REAL> inverse(dim,dim);
                // transf.Mult().DeterminantInverse(det,inverse);
                det = transf.Mult()(0,0)*transf.Mult()(1,1) - transf.Mult()(0,1)*transf.Mult()(1,0);
                orientation = DFN::sgn(det);
                break;
            }
            default:{PZError<<"\nCurrently implemented for 2D els only\n";DebugStop();}
        }

        if(orientation != 1 && orientation != -1) DebugStop();
        return orientation;
    }



    int SkeletonOrientation(TPZGeoElSide& volside, TPZGeoEl* face){
        TPZGeoElSide faceside(face,face->NSides()-1);
        // Get transformation from volume to skeleton face
        TPZTransform<REAL> transf = volside.NeighbourSideTransform(faceside);
        if(transf.Mult().Rows() != 2 || transf.Mult().Cols() != 2) {PZError<<"\n[Error] Shouldn't an affine transformation between 2D parametric spaces be 2x2?\n";DebugStop();}
        // Compute matrix determinant
        REAL det = 0.;
        // TPZFMatrix<REAL> test = {{0,1},{-1,-1}};
        // TPZFMatrix<REAL> inverse;
        // test.InitializePivot();
        // test.DeterminantInverse(det,inverse);
        // transf.Mult().DeterminantInverse(det,inverse); //(@note Can't use this because LUDecomposition sometimes changes the sign of the determinant)
        det = transf.Mult()(0,0)*transf.Mult()(1,1) - transf.Mult()(0,1)*transf.Mult()(1,0);
        // positive determinant means skeleton's orientation matches volume's faceside orientation
        int neig_side_orientation = sgn(det);
        // get 2D side orientation
        int sideorientation = (PZOutwardPointingFace(volside)?1:-1);
        // return their product
        return sideorientation * neig_side_orientation;
    }

    void SetEdgesMaterialId(TPZGeoEl* gel, int matid){
        TPZManVector<int64_t,4> edgeindices = DFN::GetEdgeIndices(gel);
        TPZGeoMesh* gmesh = gel->Mesh();
        for(int64_t index : edgeindices){
            gmesh->Element(index)->SetMaterialId(matid);
        }
    }




    void ImportElementsFromGMSH(TPZGeoMesh * gmesh, int dimension, std::set<int64_t> &oldnodes, TPZStack<int64_t> &newelements){
        // GMsh does not accept zero index entities
        constexpr int shift = 1;

        // First, if GMsh has created new nodes, these should be inserted in PZGeoMesh
        // create a map <node,point>
        std::map<int,int> mapGMshToPZ;

        for(int64_t pznode : oldnodes){
            std::vector<size_t> node_identifiers;
            std::vector<double> coord;
            std::vector<double> parametricCoord;
            gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord,0,pznode+shift,true);
            int gmshnode = (int) node_identifiers[0];
            // insert with hint (since oldnodes is an already sorted set, these nodes will all go in the end)
            mapGMshToPZ.insert(mapGMshToPZ.end(),{gmshnode,pznode+shift});
        }

        // add new nodes into PZGeoMesh
        {
            // get all nodes from GMsh
                std::vector<size_t> node_identifiers;
                std::vector<double> coord;
                std::vector<double> parametricCoord;
                gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord);
            // iterate over node_identifiers
            int nnodes = node_identifiers.size();
            for(int i = 0; i < nnodes; i++){
                int gmshnode = node_identifiers[i];
                // check if it is contained in the map
                if(mapGMshToPZ.find(gmshnode) == mapGMshToPZ.end()){
                    // New node -> add to PZGeoMesh
                    int pznode = (int) gmesh->NodeVec().AllocateNewElement();
                    TPZManVector<REAL,3> newnodeX(3);
                    newnodeX[0] = coord[3*i];
                    newnodeX[1] = coord[3*i+1];
                    newnodeX[2] = coord[3*i+2];
                    gmesh->NodeVec()[pznode].Initialize(newnodeX,*gmesh);
                    // int pznode = (int) gmesh->NNodes();
                    // gmesh->NodeVec().resize(pznode+1);
                    // insert it in map
                    mapGMshToPZ.insert({gmshnode,pznode+shift});
                }

            }
        }
        

        
        int64_t nels = gmesh->NElements();
        std::vector<std::pair<int, int> > dim_to_physical_groups;
        gmsh::model::getPhysicalGroups(dim_to_physical_groups,dimension);
        
        /// inserting the elements
        for (auto group: dim_to_physical_groups) {
        
            int dim = group.first;
            // only want elements of a given dimension
            if(dim != dimension) continue;
            int physical_identifier = group.second;
        
            std::vector< int > entities;
            gmsh::model::getEntitiesForPhysicalGroup(dim, physical_identifier, entities);

            for (auto tag: entities) {
            // std::cout<<"______________________test - tag = "<<tag;
            
                std::vector<int> group_element_types;
                std::vector<std::vector<std::size_t> > group_element_identifiers;
                std::vector<std::vector<std::size_t> > group_node_identifiers;
                gmsh::model::mesh::getElements(group_element_types,group_element_identifiers,group_node_identifiers, dim, tag);
                int n_types = group_element_types.size();
                for (int itype = 0; itype < n_types; itype++){
                    int el_type = group_element_types[itype];
                    int n_nodes = TPZGeoMeshBuilder::GetNumberofNodes(el_type);
                    std::vector<int> node_identifiers(n_nodes);
                    int n_elements = group_element_identifiers[itype].size();
                    for (int iel = 0; iel < n_elements; iel++) {
                        // int el_identifier = group_element_identifiers[itype][iel]+nels;
                        int el_identifier = gmesh->CreateUniqueElementId()+gmshshift;
                        // std::cout<<"\n"<<el_identifier<<"\n";

                        for (int inode = 0; inode < n_nodes; inode++) {
                            // Use mapGMshToPZ to translate from GMsh node index to PZ nodeindex
                            node_identifiers[inode] = mapGMshToPZ[group_node_identifiers[itype][iel*n_nodes+inode]];
                        }
                        TPZGeoEl* newel = TPZGeoMeshBuilder::InsertElement(gmesh, physical_identifier, el_type, el_identifier, node_identifiers);
                        newelements.push_back(newel->Index());
                        #ifdef PZDEBUG
                            if(newel->Dimension() != dimension){ 
                                PZError << "\nGmsh tried to group a " << newel->Dimension() << "D element as " << dimension << "D.\n";
                                gmesh->Element(el_identifier)->Print(PZError);
                                DebugStop();
                            }
                        #endif // PZDEBUG
                        // int64_t ntest = gmesh->NElements();
                        // std::cout<<"nelements = "<<ntest<<"\n";
                    }
                }
            }
        }
        gmesh->BuildConnectivity();

#if PZDEBUG
//         for(int64_t index : newelements){
//             TPZGeoEl* gel = gmesh->Element(index);
//             TPZVec<REAL> qsi(gel->Dimension(),0.); TPZFNMatrix<9,REAL> jac; TPZFNMatrix<9,REAL> axes; REAL detjac; TPZFNMatrix<9,REAL> jacinv;
//             gel->Jacobian(qsi,jac,axes,detjac,jacinv);
//             if(detjac < gDFN_SmallNumber){
//                 DFN::PrintGeoEl(gel);
//                 LOGPZ_FATAL(logger,"Gmsh created element with zero or negative det jacobian."
//                                     << "\n\tdetjac = " << detjac
//                                     << "\n\tElement = \n" << gel
//                 );
//                 DebugStop();
//             }
//         }
#endif
    }

    TPZGeoEl* FindCommonNeighbour(TPZGeoElSide& gelside1, TPZGeoElSide& gelside2, TPZGeoElSide& gelside3, int dim){
        TPZGeoMesh* gmesh = gelside1.Element()->Mesh();
        std::set<int64_t> neighbours1;
        std::set<int64_t> neighbours2;
        std::set<int64_t> neighbours3;

        TPZGeoElSide neig;
        for(neig = gelside1.Neighbour(); neig != gelside1; neig = neig.Neighbour()){
            TPZGeoEl* neig_el = neig.Element();
            if(dim > -1 && neig_el->Dimension() != dim) continue;
            neighbours1.insert(neig_el->Index());
        }
        if(neighbours1.size() < 1) return nullptr;
        for(neig = gelside2.Neighbour(); neig != gelside2; neig = neig.Neighbour()){
            TPZGeoEl* neig_el = neig.Element();
            if(dim > -1 && neig_el->Dimension() != dim) continue;
            neighbours2.insert(neig_el->Index());
        }
        if(neighbours2.size() < 1) return nullptr;
        for(neig = gelside3.Neighbour(); neig != gelside3; neig = neig.Neighbour()){
            TPZGeoEl* neig_el = neig.Element();
            if(dim > -1 && neig_el->Dimension() != dim) continue;
            neighbours3.insert(neig_el->Index());
        }
        if(neighbours3.size() < 1) return nullptr;

        std::set<int64_t> common = DFN::set_intersection(neighbours1,neighbours3);
        /** @warning: you may feel tempted to use:
         *  if(common.size() == 1) return gmesh->Element(*(common.begin()));
         *  but a common neighbour of 2 faces is not a condition for an existing face. It has to be neighbour of 3.
        */
        if(common.size() < 1) return nullptr;
        common = DFN::set_intersection(common,neighbours2);

        if(common.size() > 1) DebugStop(); // I don't think this could possibly happen, but if it ever does, I've left a weaker imposition rather than DebugStop() commented below
        // {
            // // in this case, what you probably want is 
            // for(auto& iel : common){
            //     if(gmesh->Element(*itr)->HasSubElement()) continue;
            //     return gmesh->Element(*itr);
            // }
            // DebugStop();
            // // or maybe just bet on the highest index candidate
            // return gmesh->Element(*(common.rbegin()));
        // }
        if(common.size() < 1) return nullptr;

        return gmesh->Element(*(common.begin()));

    }



    std::set<int64_t> YoungestChildren(const std::set<int64_t>& parents, TPZGeoMesh* gmesh){
        TPZStack<TPZGeoEl*> childrenstack;
        for(const auto index : parents){
            TPZGeoEl* parentgel = gmesh->Element(index);
            if(parentgel->HasSubElement()){
                parentgel->YoungestChildren(childrenstack);
            }else{
                childrenstack.push_back(parentgel);
            }
        }
        std::set<int64_t> childrenset;
        for(const TPZGeoEl* child : childrenstack){
            childrenset.insert(child->Index());
        }
        return childrenset; // move semantics makes this O(1)
    }

    template<typename Ttype>
    bool Is2PIorZero(Ttype angle, Ttype tolerance){
        PZError << "\nUninplemented type\n";
        DebugStop();
        return false;
    }

    // float
    template<>
    bool Is2PIorZero<float>(float angle, float tolerance){
        return abs(float(DFN::_2PI) - angle) < tolerance || abs(angle) < tolerance;
    }
    // double
    template<>
    bool Is2PIorZero<double>(double angle, double tolerance){
        return abs(DFN::_2PI - angle) < tolerance || abs(angle) < tolerance;
    }

    void SketchStatusVec(const TPZManVector<int,9>& statusvec, std::ostream& out){
		std::stringstream sout;
		sout << "\nStatusVec Sketch:";
		switch(statusvec.size()){
			case pztopology::TPZTriangle::NSides:{
				sout << "\n  " << (statusvec[2]?"x":" ");
				sout << "\n | \\";
				sout << "\n" << (statusvec[5]?"x":" ") <<"|  \\" << (statusvec[4]?"x":" ");
				sout << "\n |___\\";
				sout << "\n" << (statusvec[0]?"x":" ") << "  " << (statusvec[3]?"x":" ") << "  " << (statusvec[1]?"x":" ");
				break;
			}
			case pztopology::TPZQuadrilateral::NSides:{
				sout << "\n" << (statusvec[3]?"x":" ") << " __" << (statusvec[6]?"x":"_") << "__ " << (statusvec[2]?"x":" ");
				sout << "\n |     | ";
				sout << "\n" << (statusvec[7]?"x":" ") << "|     |" << (statusvec[5]?"x":" ");
				sout << "\n |_____| ";
				sout << "\n" << (statusvec[0]?"x":" ") << "   " << (statusvec[4]?"x":" ") << "   " << (statusvec[1]?"x":" ");
				break;
			}
			default: DebugStop();
		}

		out << sout.str();
    }

    bool AreCoPlanar(TPZGeoMesh* gmesh,const std::set<int64_t>& nodeindices, REAL tolerance){
        int nnodes = nodeindices.size();
        if(nnodes < 4) return true; // this return may introduce a bug for sets of 3 nodes that are colinear, I'll fix it if the code ever breaks because of it.
        TPZManVector<TPZManVector<REAL,3>,3> p(3,{0.,0.,0.});

        auto it = nodeindices.begin();
        // In case we stumble on a triple of colinear points, we can keep on trying until a non-colinear triple is found or all possible triples were tried
        for(int attempt=1; attempt <= nnodes; attempt++){
            if(attempt > 1) it++;

            // This is a little hard to read, but it's done this way to work with a template container. It would be easier to read if I'd limited the function to vectors, but I need it to work primarily with std::sets
            auto auxit = it;
            // p[0]
            gmesh->NodeVec()[ *it  ].GetCoordinates(p[0]);
            // p[1]
            if(++auxit == nodeindices.end()) auxit = nodeindices.begin();
            gmesh->NodeVec()[*auxit].GetCoordinates(p[1]);
            // p[2]
            if(++auxit == nodeindices.end()) auxit = nodeindices.begin();
            gmesh->NodeVec()[*auxit].GetCoordinates(p[2]);

            
            TPZManVector<REAL,3> vec0 = p[0] - p[1]; // < Vector connecting p[1] to p[0]
            TPZManVector<REAL,3> vec2 = p[2] - p[1]; // < Vector connecting p[1] to p[2]
            TPZManVector<REAL,3> vecj {0., 0., 0.};  // < Vector connecting p[1] to p[j]

            // Get a normal vector to the p0,p1,p2 plane
            TPZManVector<REAL,3> normal = CrossProduct<REAL>(vec0,vec2);
            REAL norm = DFN::Norm<REAL>(normal);
            // Test if p0, p1 and p2 are colinear by check if vec0 and vec2 are linearly independent
            if(fabs(norm) < tolerance){
                continue; // get the next triple of points and try again
            }
            normal[0] /= norm;
            normal[1] /= norm;
            normal[2] /= norm;

            // Test every other j-th node for coplanarity
            TPZManVector<REAL,3> pj {0., 0., 0.};
            if(++auxit == nodeindices.end()) auxit = nodeindices.begin();
            for(/*void*/; auxit != it; (++auxit == nodeindices.end() ? auxit = nodeindices.begin() : auxit)){
                int64_t jnode = *auxit;
                gmesh->NodeVec()[jnode].GetCoordinates(pj);
                vecj = pj - p[1];
                REAL orth_dist = DotProduct<REAL>(normal,vecj);
                if(fabs(orth_dist) > tolerance) return false; 
            }
            return true;
        }
        /* Surely, colinear points are indeed co-planar, but that's not what this function was intended for
        if you want to adapt it to return true for co-linear points, maybe a bool argument to give the
        choice to the user would be the better way go. It should of course default to not supporting it.
        */
        PZError << "\n" << __PRETTY_FUNCTION__;
        PZError << "\n\t Failed due to co-linear points";
        DebugStop();
        return false;
    }

    void PlotJacobian(TPZGeoMesh* gmesh, std::string filename /*="LOG/jacplot.vtk"*/){

        int64_t nels = gmesh->NElements();
        const REAL notcomputed = -999.0;
        TPZVec<REAL> elData(nels,notcomputed);

        for(auto gel : gmesh->ElementVec()){
            if(!gel) continue;
            if(gel->HasSubElement()) continue;

            TPZVec<REAL> qsi(gel->Dimension(),0.); TPZFMatrix<double> jac; TPZFMatrix<double> axes; REAL detjac = notcomputed; TPZFMatrix<double> jacinv;
            gel->Jacobian(qsi,jac,axes,detjac,jacinv);
            elData[gel->Index()] = detjac;
        }

        std::ofstream file(filename);
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,file, elData);
    }
} /* namespace DFN*/
