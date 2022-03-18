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
    // Return angle in radians from a cosine and a sine within the range [0,2pi)
    REAL Arc(const REAL cos, const REAL sin){
#ifdef PZDEBUG
        constexpr REAL fuzzyOne = 1.+gDFN_SmallNumber;
        if(cos < -fuzzyOne || cos > fuzzyOne || sin < -fuzzyOne || sin > fuzzyOne)
            DebugStop();
#endif // PZDEBUG
        const REAL arccos = std::acos(cos);
        constexpr REAL _2pi = 2.*M_PI;
        return sin > 0. ? arccos : _2pi - arccos;
    }

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
    TPZGeoEl* GetSkeletonNeighbour(const TPZGeoEl* gel, const int side){
        if(gel->SideDimension(side) == gel->Dimension()) return nullptr;
        const TPZGeoElSide gelside(const_cast<TPZGeoEl*>(gel),side);
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
    void BestFitPlane(const TPZFMatrix<REAL>& pointcloud, TPZManVector<REAL,3>& centroid, TPZManVector<REAL,3>& normal){
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
            float norm = Norm<REAL>(normal);
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
    REAL DihedralAngle(const TPZGeoElSide &gelside, const TPZGeoElSide &neighbour, const int sideorientation){

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
        TPZManVector<REAL,3> sharednode0(3,0);
        TPZManVector<REAL,3> sharednode1(3,0);
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
        
        TPZManVector<REAL,3> normalvec_gel = CrossProduct<REAL>(tangentvec_gel,tangentvec_edge);
        TPZManVector<REAL,3> normalvec_neig = CrossProduct<REAL>(tangentvec_neig,tangentvec_edge);;
        TPZManVector<REAL,3> aux = CrossProduct<REAL>(normalvec_neig,normalvec_gel);
        REAL x = Norm<REAL>(tangentvec_edge)*DotProduct<REAL>(normalvec_neig,normalvec_gel);
        REAL y = DotProduct<REAL>(tangentvec_edge,aux);
        REAL angle = atan2(y,x);
        
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
        return ancestors.size() > 2;
    }

    /** @brief Alternate name of dfn::IsValidPolygon for readability
     * @returns False
    */
    bool IsDegeneratePolygon(TPZStack<int64_t>& polygon, TPZGeoMesh* gmesh){
        return !IsValidPolygon(polygon,gmesh);
    }

    /** @brief Check if a 2D side of a 3D element is oriented outward according to NeoPZ topolgy */
    bool PZOutwardPointingFace(const TPZGeoElSide faceside){
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
    REAL CornerAngle_cos(TPZGeoEl *gel, const int corner){
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

        REAL cosine = DotProduct<REAL>(vec1,vec2)/(Norm<REAL>(vec1)*Norm<REAL>(vec2));
        return cosine;
    }

    /**
     * @brief Check if the side that connects 2 neighbours has the same orientation in each element
     * @note currently exclusive to 1D sides
     */
    bool OrientationMatch(const TPZGeoElSide &neig1, const TPZGeoElSide &neig2){
        if(neig1.Dimension() != 1) DebugStop();
        if(!neig1.NeighbourExists(neig2)) DebugStop();
        return (neig1.SideNodeIndex(0) == neig2.SideNodeIndex(0));
    }



    /**
     * @brief Get a vector from node 0 to node 1 of a 1D side
     */
    void GetSideVector(const TPZGeoElSide &gelside, TPZManVector<REAL,3>& vector){
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


    int SubElOrientation(TPZGeoEl* father, const int ichild){
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

    void SetEdgesMaterialId(const TPZGeoEl* gel, const int matid){
        for(int iside=gel->FirstSide(1); iside<gel->NSides(); iside++){
            TPZGeoEl* edge = DFN::GetSkeletonNeighbour(gel,iside);
            if(!edge) continue;
            edge->SetMaterialId(matid);
        }
    }




    void ImportElementsFromGMSH(TPZGeoMesh * gmesh, const int dimension, const std::set<int64_t> &oldnodes, TPZStack<int64_t> &newelements){
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
        // DFN::TestPositiveJacobian(gmesh,newelements);
        // DFN::TestInternalAngles(gmesh,newelements,2.*M_PI/180.);
#endif // PZDEBUG
    }

    TPZGeoEl* FindCommonNeighbour(const TPZGeoElSide gelside1, const TPZGeoElSide gelside2, const TPZGeoElSide gelside3, const int dim){
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
    bool Is2PIorZero(const Ttype angle, const Ttype tolerance){
        PZError << "\nUninplemented type\n";
        DebugStop();
        return false;
    }

    // float
    template<>
    bool Is2PIorZero<float>(const float angle, const float tolerance){
        return abs(float(DFN::_2PI) - angle) < tolerance || abs(angle) < tolerance;
    }
    // double
    template<>
    bool Is2PIorZero<double>(const double angle, const double tolerance){
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

    bool AreCoPlanar(TPZGeoMesh* gmesh,const std::set<int64_t>& nodeindices, const REAL tolerance){
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



    void PlotNeighbours(const std::string filepath,
                        const TPZGeoElSide geoelside, 
                        const int filterDimension, 
                        const bool UnrefinedOnly, 
                        const bool orientationMatch){
        
        // Consistency
        if(geoelside.Element() == nullptr) DebugStop();

        TPZGeoMesh* gmesh = geoelside.Element()->Mesh();
        TPZGeoMesh auxMesh;
        auxMesh.ElementVec().Resize(gmesh->NElements());
        for (int i = 0; i < auxMesh.NElements(); i++) {
            auxMesh.ElementVec()[i] = nullptr;
        }
        auxMesh.NodeVec() = gmesh->NodeVec();

        TPZGeoElSide neig = geoelside.Neighbour();
        for(/*void*/; neig != geoelside; neig++){
            if(filterDimension > 0 && neig.Element()->Dimension() != filterDimension) continue;
            if(UnrefinedOnly && neig.Element()->HasSubElement()) continue;
            
            TPZGeoEl* newel = neig.Element()->Clone(auxMesh);
            if(orientationMatch){
                int orientation = DFN::OrientationMatch(geoelside,neig);
                newel->SetMaterialId(orientation);
            }
        }
        std::ofstream out(filepath);
        TPZVTKGeoMesh::PrintGMeshVTK(&auxMesh, out, true, true);
    }

    void PlotVTK_elementList(const std::string filepath, const TPZVec<int64_t>& newgels, TPZGeoMesh* gmesh){
        if(newgels.size() == 0) return;

        TPZGeoMesh graphicMesh;
        graphicMesh.ElementVec().Resize(gmesh->NElements());
#ifdef MACOSX
        for (int i = 0; i < graphicMesh.NElements(); i++) {
            graphicMesh.ElementVec()[i] = nullptr;
        }
#else
        std::fill(graphicMesh.ElementVec().begin(),graphicMesh.ElementVec().end(),nullptr);
#endif
        graphicMesh.NodeVec() = gmesh->NodeVec();
        for(int64_t index : newgels){
            TPZGeoEl* newel = gmesh->Element(index)->Clone(graphicMesh);
            if(newel->HasSubElement()) newel->SetSubElement(0,nullptr);
        }

        std::ofstream out(filepath);
        TPZVTKGeoMesh::PrintGMeshVTK(&graphicMesh, out, true, true);
    }
    void PlotVTK_SideList(const std::string filepath, const TPZVec<TPZGeoElSide>& sidelist){
        if(sidelist.size() == 0) return;
        TPZGeoMesh* gmesh = sidelist[0].Element()->Mesh();
        TPZGeoMesh graphicMesh;
        graphicMesh.ElementVec().Resize(gmesh->NElements());
#ifdef MACOSX
        for (int i = 0; i < graphicMesh.NElements(); i++) {
            graphicMesh.ElementVec()[i] = nullptr;
        }
#else
        std::fill(graphicMesh.ElementVec().begin(),graphicMesh.ElementVec().end(),nullptr);
#endif
        graphicMesh.NodeVec() = gmesh->NodeVec();
        
        const int64_t nels = gmesh->NElements()+sidelist.size();
        TPZVec<REAL> elData(nels,-1);

        for(TPZGeoElSide geoside : sidelist){
            if(!geoside.Element()){
                int discard = graphicMesh.ElementVec().AllocateNewElement(); // To make sure elData has the same size as graphicMesh.ElementVec
                graphicMesh.ElementVec()[discard] = nullptr;
                continue;
            }
            TPZGeoEl* newel = geoside.Element()->Clone(graphicMesh);
            newel->SetMaterialId(0);
            newel->ResetSubElements();
            
            const int nsides = newel->NSides();
            for(int i=0;i<nsides;i++) newel->SetSideDefined(i);

            TPZGeoEl* newbc = newel->CreateBCGeoEl(geoside.Side(),1);
            elData[newbc->Index()] = geoside.Side();
        }

        std::ofstream out(filepath);
        TPZVTKGeoMesh::PrintGMeshVTK(&graphicMesh, out, elData);
    }

    /// @brief Check if internal angles of a planar quadrilateral SubPolygon violate a threshold.
    /// @details We're trying to avoid a quadrilateral with consecutive parallel edges, so this functions warns the code of an angle that gets too close to 180deg
    /// @param tol_angle_cos cos(t),  t being the maximum tolerable angle (convexity of polyhedral volumes guarantees angles between 0 and 180deg)
    /// @param splitThisAngle [output] local index of the first angle that violates the threshhold
    bool CheckSubPolygonInternalAngles(TPZGeoMesh* gmesh,const TPZVec<int64_t>& subpolygon, const REAL tol_angle_cos, int& splitThisAngle){
        const int nelements = subpolygon.size();
        bool paralleledgesQ = false;
        for (int64_t in = 0 ; in < nelements ; in++) {
            TPZManVector<REAL,3> c0(3,-1.),c1(3,-1.),vec0(3,-1.),vec1(3,-1.);
            {
                int orientation = sgn(subpolygon[in]);
                const int64_t index = std::abs(subpolygon[in]);
                TPZGeoEl* gel = gmesh->Element(index);
                if(gel->Dimension() != 1) DebugStop();
                gel->NodePtr(0)->GetCoordinates(c0);
                gel->NodePtr(1)->GetCoordinates(c1);
                vec0 = c1 - c0;
                vec0[0] *= -orientation;
                vec0[1] *= -orientation;
                vec0[2] *= -orientation;
            }

            {
                int orientation = sgn(subpolygon[(in+1)%nelements]);
                const int64_t index = std::abs(subpolygon[(in+1)%nelements]);
                TPZGeoEl* gel = gmesh->Element(index);
                if(gel->Dimension() != 1) DebugStop();
                gel->NodePtr(0)->GetCoordinates(c0);
                gel->NodePtr(1)->GetCoordinates(c1);
                vec1 = c1 - c0;
                vec0[0] *= orientation;
                vec0[1] *= orientation;
                vec0[2] *= orientation;
            }
            
            const REAL norm0 = DFN::Norm<REAL>(vec0);
            const REAL norm1 = DFN::Norm<REAL>(vec1);
            REAL cosangle = DFN::DotProduct<REAL>(vec0,vec1) / norm0 / norm1;
            

            if(cosangle < tol_angle_cos){
                if (paralleledgesQ){
                    DebugStop(); // 3 parallel edges!
                }
                splitThisAngle = (in+1)%nelements;
                paralleledgesQ = true;
                return paralleledgesQ;
            }
        }
        return paralleledgesQ;
    }

    // Fill a vector with the global node indices of a subpolygon in the order they appear locally
    void GetSubPolygonGlobalNodeIndices(TPZGeoMesh* gmesh,const TPZVec<int64_t>& subpolygon, TPZVec<int64_t>& globNodes){
        const int nelements = subpolygon.size();
        globNodes.Resize(nelements,-1);
        for (int64_t in = 0 ; in < nelements ; in++) {
            const int64_t index = std::abs(subpolygon[in]);
            TPZGeoEl* edge = gmesh->Element(index);
            if(subpolygon[in] < 0){
                globNodes[in] = edge->NodeIndex(1);
            }else{
                globNodes[in] = edge->NodeIndex(0);
            }
        }        
    }


    /** @brief Takes a simple oriented lineloop with 3 or 4 edges and create a geometric element
     * @param lineloop an oriented loop of 3 or 4 edges
    */  
    // template<class Tcontainer>
    void MeshSimplePolygon(TPZGeoMesh* gmesh,const TPZVec<int64_t>& lineloop, int matid, TPZStack<int64_t>& newelements){
        int nelements = lineloop.size();
        if(nelements < 3 || nelements > 4) DebugStop();
        newelements.clear();
        
        // associate an elementside with each node in lineloop
        std::map<int64_t,TPZGeoElSide> nodeconnectivity;
        for(auto elindex : lineloop)
        {
            TPZGeoEl *gel = gmesh->Element(std::abs(elindex));
            int ncorner = gel->NCornerNodes();
            for (int ic=0; ic<ncorner; ic++) {
                int64_t nodeindex = gel->NodeIndex(ic);
                TPZGeoElSide gelside(gel,ic);
                nodeconnectivity[nodeindex] = gelside;
            }
        }
        // parallelEdgesQ will be true if the polygon has 4 elements and one of the internal angles
        // is larger than 170 degrees
        bool parallelEdgesQ = false;
        int midnode = -1;
        
        TPZManVector<int64_t,4> globNodes(nelements,-1);
        
        if(nelements == 4){
            GetSubPolygonGlobalNodeIndices(gmesh,lineloop,globNodes);
            constexpr REAL cos170 = -0.984807753;
            // midnode is the node that is aligned
            parallelEdgesQ = CheckSubPolygonInternalAngles(gmesh,lineloop,cos170,midnode);
            if (parallelEdgesQ && midnode == -1) DebugStop();
        }

        // we are generating 2 triangles from 4 sides
        if (parallelEdgesQ){
            TPZManVector<int64_t,3> cornerindices0(3,-1),cornerindices1(3,-1);
            cornerindices0[0] = globNodes[midnode];
            cornerindices0[1] = globNodes[(midnode+1)%nelements];
            cornerindices0[2] = globNodes[(midnode+2)%nelements];

            cornerindices1[0] = globNodes[(midnode+2)%nelements];
            cornerindices1[1] = globNodes[(midnode+3)%nelements];
            cornerindices1[2] = globNodes[midnode];
            
            MElementType eltype = MElementType::ETriangle;
            
            int64_t elindex = -1;
            TPZGeoEl* new_el0 = gmesh->CreateGeoElement(eltype,cornerindices0,matid,elindex);
            newelements.push_back(new_el0->Index());
            TPZGeoEl* new_el1 = gmesh->CreateGeoElement(eltype,cornerindices1,matid,elindex);
            newelements.push_back(new_el1->Index());
            // gmesh->BuildConnectivity();

            constexpr int desloc = 2;
            const int nedges = nelements;
            for(int iedge=0; iedge<2; iedge++){
                {
                    int edge_index_in_subpolygon = (midnode+iedge)%nedges;
                    TPZGeoEl* edge = gmesh->Element(std::abs(lineloop[edge_index_in_subpolygon]));
                    TPZGeoElSide edgeside(edge,2);
                    TPZGeoElSide faceside(new_el0,iedge+3);
                    faceside.SetConnectivity(edgeside);
                    // @maybeToDo set Node connectivity
                }
                {
                    int edge_index_in_subpolygon = (midnode+iedge+desloc)%nedges;
                    TPZGeoEl* edge = gmesh->Element(std::abs(lineloop[edge_index_in_subpolygon]));
                    TPZGeoElSide edgeside(edge,2);
                    TPZGeoElSide faceside(new_el1,iedge+3);
                    faceside.SetConnectivity(edgeside);
                    // @maybeToDo set Node connectivity
                }
            }
            // set the connectivity along the line that connects both triangles
            TPZGeoElSide faceside0(new_el0,5);
            TPZGeoElSide faceside1(new_el1,5);
            faceside0.SetConnectivity(faceside1);
        }
        else{ // create a triangle or quadrilateral
            TPZManVector<int64_t,4> cornerindices(nelements,-1);
            int i=0;
            for(int64_t edge : lineloop){
                int64_t index = std::abs(edge);
                if(edge < 0){
                    cornerindices[i] = gmesh->Element(index)->NodeIndex(1);
                }else{
                    cornerindices[i] = gmesh->Element(index)->NodeIndex(0);
                }
                i++;
            }
            int64_t index = -1;
            MElementType eltype;
            switch (nelements){
    //            case  2: eltype = MElementType::EInterfaceLinear; break;
                case  3: eltype = MElementType::ETriangle; break;
                case  4: eltype = MElementType::EQuadrilateral; break;
                default: DebugStop();
            }
            TPZGeoEl* new_el = gmesh->CreateGeoElement(eltype,cornerindices,matid,index);
            newelements.push_back(new_el->Index());

            const int nedges = nelements;
            for(int iedge=0; iedge<nelements; iedge++){
                TPZGeoEl* edge = gmesh->Element(std::abs(lineloop[iedge]));
                TPZGeoElSide edgeside(edge,2);
                TPZGeoElSide faceside(new_el,iedge+nedges);
                faceside.SetConnectivity(edgeside);
                // @maybeToDo set Node connectivity
            }
        }
        for(auto el : newelements)
        {
            TPZGeoEl *gel = gmesh->Element(el);
            int nc = gel->NCornerNodes();
            for (int ic = 0; ic<nc; ic++) {
                TPZGeoElSide gelside(gel,ic);
                int64_t nodeindex = gel->NodeIndex(ic);
#ifdef PZDEBUG
                if(gelside.Neighbour()) DebugStop();
                if(nodeconnectivity.find(nodeindex) == nodeconnectivity.end()) DebugStop();
#endif
                gelside.SetConnectivity(nodeconnectivity[nodeindex]);
            }
            TPZGeoElSide gelside(gel);
            gelside.SetConnectivity(gelside);
        }
    }


    TPZGeoElSide TestInternalAngles(TPZGeoMesh* gmesh, const TPZVec<int64_t>& el_indices, const REAL tol){
        for(const int64_t index : el_indices){
            TPZGeoEl* gel = gmesh->Element(index);
            const int d = gel->Dimension() - 2;
            const int nsides = gel->NSides();
            for(int iside=0; iside<nsides; iside++){
                if(gel->SideDimension(iside) != d) continue;
                const TPZGeoElSide subfacet{gel,iside};
                REAL angle = DFN::InternalAngle(subfacet);
                if(angle < tol)
                    return subfacet;
            }
        }
        return {nullptr,-1};
    }

    bool TestPositiveJacobian(TPZGeoMesh* gmesh, const TPZVec<int64_t>& el_indices){
        for(const int64_t index : el_indices){
            TPZGeoEl* gel = gmesh->Element(index);
            TPZVec<REAL> qsi(gel->Dimension(),0.); TPZFNMatrix<9,REAL> jac; TPZFNMatrix<9,REAL> axes; REAL detjac; TPZFNMatrix<9,REAL> jacinv;
            gel->Jacobian(qsi,jac,axes,detjac,jacinv);
            if(detjac < 1e-4){
                // DFN::PrintGeoEl(gel);
                LOGPZ_FATAL(logger,"Gmsh created element with zero or negative det jacobian."
                                    << "\n\tdetjac = " << detjac
                                    << "\n\tElement = \n" << gel
                );
                DebugStop();
                return false;
            }
        }
        return true;
    }


    REAL InternalAngle(const TPZGeoElSide subfacet){
        // Consistency
        bool consistentInput = true;
        consistentInput = consistentInput && subfacet.Element();
        consistentInput = consistentInput && (subfacet.Dimension() == subfacet.Element()->Dimension()-2);
        consistentInput = consistentInput && (subfacet.Element()->Dimension() > 1);
        if(!consistentInput){
            throw std::invalid_argument("I don't know how to compute internal angle around side "
                                        +std::to_string(subfacet.Side())
                                        +" of element "
                                        +std::to_string(subfacet.Element()->Index()));
            REAL nonsense = -10.*M_PI;
            return nonsense;
        }
        const int d = subfacet.Element()->Dimension();
        // @definition: facet = (d-1)-dimensional side of a d-dimensional element
        // @definition: subfacet = (d-2)-dimensional side of a d-dimensional element
        // A subfacet side is shared by 2 facets
        TPZStack<TPZGeoElSide> facets;
        subfacet.Element()->AllHigherDimensionSides(subfacet.Side(),2,facets);
        if(facets.size() != 2) DebugStop();

        // Compute normals of facets
        TPZManVector<REAL,3> midpointQsi(d-1,0.);
        facets[0].CenterPoint(midpointQsi);
        TPZManVector<REAL,3> normal0(3,0.);
        facets[0].Normal(midpointQsi,normal0);
        facets[1].CenterPoint(midpointQsi);
        TPZManVector<REAL,3> normal1(3,0.);
        facets[1].Normal(midpointQsi,normal1);
        
        // @note: TPZGeoElSide::Normal() always builds an Outward pointing unitary vector

        // Compute cosine of supplementary angle
        const REAL cos_sup = DFN::DotProduct<REAL>(normal0,normal1);
        // Compute sine of supplementary angle
        TPZManVector<REAL,3> cross = DFN::CrossProduct<REAL>(normal0,normal1);
        const REAL sin_sup = DFN::Norm<REAL>(cross);
        // Compute supplementary angle
        // const REAL sup_angle = DFN::Arc(cos_sup,sin_sup);
        
        // Compute internal angle
        const REAL cos = -cos_sup;
        const REAL sin =  sin_sup;
        REAL internal_angle = DFN::Arc(cos,sin);
        // REAL internal_angle = M_PI - sup_angle;

#ifdef PZDEBUG
//         // If you'd like to account for non-convex elements, it goes something like this
//         TPZManVector<double,3> qsi(d,0.); TPZFMatrix<double> jac(d,d); TPZFMatrix<double> axes(3,d); REAL detjac = -999.0; TPZFMatrix<double> jacinv(d,d);
//         subfacet.Element()->Jacobian(qsi,jac,axes,detjac,jacinv);
//         if(detjac < 0.)
//             {internal_angle = DFN::_2PI - internal_angle;}
#endif // PZDEBUG

        return internal_angle;
    }
} /* namespace DFN*/


void DFN::GmshConfig(){
    // gmsh::logger::start();
    // std::vector<std::string> log;
    // gmsh::logger::get(log);
    // // gmsh::logger::stop();
    // gmsh::option::
    gmsh::option::setNumber("General.Verbosity",1);
}
