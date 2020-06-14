/*! 
 *  @brief     Contains implementation of methods for class DFNMesh.
 *  @authors   Pedro Lima
 *  @date      2020.06
 */

#include "DFNMesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZGeoMeshBuilder.h"



DFNMesh::DFNMesh(TPZGeoMesh *gmesh){
	fGMesh = gmesh;
	for (int i = 1; i < fGMesh->Dimension(); i++){
		CreateSkeletonElements(i,-1);
	}
	
}


void DFNMesh::Print(std::string pzmesh
                    ,std::string vtkmesh
                    ,int fracture
                    ,int transition
                    ,int intact)
{
	TPZGeoMesh *gmesh = this->fGMesh;
	if(pzmesh != "skip"){
		std::ofstream meshprint(pzmesh);
		gmesh->Print(meshprint);
	}
	if(vtkmesh != "skip"){
		std::ofstream out1(vtkmesh);
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out1, true, true);
	}
}








void PrintColorful(std::string pzmesh = "pzmesh.txt",std::string vtkmesh = "vtkmesh.vtk"){
	
	// for()
	// TPZVTKGeoMesh::PrintGMeshVTK()
}











void DFNMesh::GenerateSubMesh(){
    TPZGeoMesh *gmesh = this->Mesh();
	int dim = gmesh->Dimension();
	if(dim != 3) return; //@todo temporary while DFNVolume isn't generalized to 2D meshes
	
	gmsh::initialize();
	//Loop over list of fractured volumes
	for (auto itr = fVolumes.begin(); itr != fVolumes.end(); itr++){
    	DFNVolume *ivolume = &itr->second;
		// if(ivolume->ElementIndex() != 1) continue;
		// Use GMsh to tetrahedralize volumes
    	Tetrahedralize(ivolume);
	}
	gmsh::finalize();	
}



/**
 * @brief Deletes face + ribs (if left isolated) + nodes (if left isolated)
 * @param face: 2D element to be deleted
*/
void DFNMesh::DeleteElementAndRibs(TPZGeoEl *face){
	TPZGeoMesh *gmesh = face->Mesh();
	// queue ribs
	int nribs = face->NNodes();
	// int nribs = (face->Type() == MElementType::EQuadrilateral ? 4 : 3);
	TPZManVector<int64_t, 4> ribs(nribs,-1);
	for(int irib = 0; irib < nribs; irib++){
		TPZGeoElSide faceside(face,irib+nribs);
		for(TPZGeoElSide neighbour = faceside.Neighbour(); neighbour != faceside; neighbour = neighbour.Neighbour()){
			if(neighbour.Element()->Dimension() == 1){
				ribs[irib] = neighbour.Element()->Index(); 
				break;
			}
		}
	}
	// delete face
	gmesh->DeleteElement(face);
	// then, delete face's ribs and nodes if necessary
	for(int irib = 0; irib < nribs; irib++){
		if(ribs[irib] == -1) continue;
		TPZGeoEl *ribgel = gmesh->Element(ribs[irib]);
		TPZManVector<int64_t,2> nodes(2,-1);
		for(int iside = 0; iside < 3 ; iside++){
			TPZGeoElSide ribgelside(ribgel,iside);
			TPZGeoElSide neighbour = ribgelside.Neighbour();
			bool delete_Q = true;
			// @todo maybe add exception for 0D neighbour
			if(iside < 2){
				if(neighbour == ribgelside){
					nodes[iside] = ribgel->NodeIndex(iside);
				}
			}else{
				while(neighbour != ribgelside){
					if(neighbour.Element()->Dimension() == 2){
						delete_Q = false; 
						break;
					}
					neighbour = neighbour.Neighbour();
				}
				if(delete_Q){gmesh->DeleteElement(ribgel,ribs[irib]);}
			}
		}
		// delete node
		for(int inode : nodes){
			if(inode < 0) continue;
			gmesh->NodeVec().SetFree(inode);
			gmesh->NodeVec()[inode].SetNodeId(-1);
			// delete &gmesh->NodeVec()[inode];
		}
	}
}









/**
 * @brief Deletes gel + children + isolated ribs + unused nodes
 * @note It will assume element has been found not to belong to the domain of interest and will not verify
 * @param gel: pointer to the geometric element
*/
void DFNMesh::CropExternalElement(TPZGeoEl *gel){
	TPZGeoMesh *gmesh = gel->Mesh();
	// If an element is external to the domain, then its eldest ancestor and all the refinement tree are also external
	TPZGeoEl *elder = gel;
	if(gel->Father()){elder = gel->EldestAncestor();}
	// Start from youngest children and go up the tree
	while(elder->HasSubElement()){
		TPZStack<TPZGeoEl*> youngestChildren;
		elder->YoungestChildren(youngestChildren);
		for(auto child : youngestChildren){
			DeleteElementAndRibs(child);
		}
	}
	DeleteElementAndRibs(elder);
}









/**
 * @brief Deletes gel and all elements that share the same eldest ancestor, then deletes the ancestor
 * @param gel: Any member of the family
*/
void DFNMesh::DeleteFamily(TPZGeoEl *gel){
	TPZGeoMesh *gmesh = gel->Mesh();
	if(!gel->Father()){
		gmesh->DeleteElement(gel);
		return;
	}
	TPZGeoEl *elder = gel->EldestAncestor();
	while(elder->HasSubElement()){ //this looks redundant, but I'll need it to delete ribs and nodes eventually
		TPZStack<TPZGeoEl*> youngestChildren;
		elder->YoungestChildren(youngestChildren);
		for(auto child : youngestChildren){
			gmesh->DeleteElement(child);
		}	
	}
	gmesh->DeleteElement(elder);
}






void DFNMesh::ImportElementsFromGMSH(TPZGeoMesh * gmesh, int dimension, std::set<int64_t> &oldnodes){
    // GMsh does not accept zero index entities
    const int shift = 1;

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

		#ifdef PZDEBUG
			if(dimension == 3){
				physical_identifier++; // to differ cut volumes from cut faces
				//std::cout<<"\n@comment: For better graphics, volumetrical transition elements have material id shifted\n";
			}
		#endif
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
                    int el_identifier = group_element_identifiers[itype][iel]-1+nels;
					// std::cout<<"\n"<<el_identifier<<"\n";

                    for (int inode = 0; inode < n_nodes; inode++) {
                        // node_identifiers[inode] = group_node_identifiers[itype][iel*n_nodes+inode];
                        // Use mapGMshToPZ to translate from GMsh node index to PZ nodeindex
                        node_identifiers[inode] = mapGMshToPZ[group_node_identifiers[itype][iel*n_nodes+inode]];
                    }
                    TPZGeoMeshBuilder::InsertElement(gmesh, physical_identifier, el_type, el_identifier, node_identifiers);
					int64_t ntest = gmesh->NElements();
					// std::cout<<"nelements = "<<ntest<<"\n";
                }
            }
        }
    }
    gmesh->BuildConnectivity();
}





























int64_t DFNMesh::SearchIndirectNeighbours(TPZGeoEl* gel){
    
    std::list<int64_t> candidate_queue;
	std::set<int64_t> verified;
    candidate_queue.push_back(gel->Index());
    for (auto index : candidate_queue) {
		if(verified.find(index) != verified.end()) continue;
        TPZGeoEl *currentgel = gel->Mesh()->Element(index);
		int64_t macroElindex = FindAdjacentMacroEl(currentgel);
        if (macroElindex == -1) {
            QueueNeighbours(currentgel, candidate_queue);
        }
        if (macroElindex >= 0) {
            return macroElindex;;
        }   
		verified.insert(index);
    }
	// If this point is reached, current element is surrounded by elements that have been deleted and, therefore, is not in the domain (and should also be deleted).
	return -1;
}







int64_t DFNMesh::FindAdjacentMacroEl(TPZGeoEl* gel){
	int surfaceMaterial = (*fFractures.begin())->GetSurfaceMaterial();
    int nsides = gel->NSides();
        
	for(int iside = nsides-2; iside >= 0; iside--){
		TPZGeoElSide gelside(gel, iside);
		if(gelside.Dimension() < 1) break;
		TPZGeoElSide neig = gelside.Neighbour();
		for( ; neig != gelside; neig = neig.Neighbour()){
			if(neig.Element()->Dimension() != 2) continue;
			int mat = neig.Element()->MaterialId();
			// if (mat != surfaceMaterial){
			if (mat < surfaceMaterial){
				return neig.Element()->EldestAncestor()->Index();
			}
		}
	}
    
    return -1;
}




void DFNMesh::QueueNeighbours(TPZGeoEl* gel, std::list<int64_t> &candidate_queue){
    
    int nsides = gel->NSides();
	int ncorners = gel->NCornerNodes();
    for(int iside = nsides-2; iside > ncorners-1; iside--){
        TPZGeoElSide gelside(gel, iside);
        TPZGeoElSide neig = gelside.Neighbour();
        for( ; neig != gelside; neig = neig.Neighbour()){
			// only need 2D neighbours
			if(neig.Element()->Dimension() != 2) continue;
            candidate_queue.push_back(neig.Element()->Index());
        }
    }
}

















bool DFNMesh::FindEnclosingVolume(TPZGeoEl *ifracface){
	bool gBigPlanes_Q = true; //placeholder
	if(ifracface->Dimension()!=2) DebugStop();
	int64_t ifracfaceindex = ifracface->Index();
	int surfaceMaterial = (*fFractures.begin())->GetSurfaceMaterial();
	TPZGeoMesh *gmesh = ifracface->Mesh();
    // get coordinates of geometric center of face
    TPZVec<REAL> faceCenter(3);
    {
        TPZGeoElSide geliside(ifracface, ifracface->NSides()-1);
        geliside.CenterX(faceCenter);
    }

    // map of indices for volumes that could contain the face
    std::map<REAL, int64_t> candidates;
    {
		int64_t macroElindex = SearchIndirectNeighbours(ifracface);
		if(macroElindex == -1){
			if(!gBigPlanes_Q){
				std::cout<<"\n "<<__PRETTY_FUNCTION__<<" found no enclosing volume for element #"<<ifracface->Index()<<"\n";
				DebugStop();
			}
			CropExternalElement(ifracface);
			// gmesh->DeleteElement(ifracface,ifracface->Index());
			return false;
		}
		TPZGeoEl *macroEl = ifracface->Mesh()->Element(macroElindex);
		// get macroEl's center coordinates
		TPZManVector<REAL,3> macroElCenter(3);
		TPZGeoElSide macroElfaceside(macroEl, macroEl->NSides()-1);
		macroElfaceside.CenterX(macroElCenter);
		// construct vector from center of macroEl to center of ifracface
		TPZManVector<REAL,3> v1(3,0);
			v1[0] = faceCenter[0] - macroElCenter[0];
			v1[1] = faceCenter[1] - macroElCenter[1];
			v1[2] = faceCenter[2] - macroElCenter[2];
		// Normalize v1
		REAL norm = sqrtl(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
			v1[0] = v1[0]/norm;
			v1[1] = v1[1]/norm;
			v1[2] = v1[2]/norm;
		// iterate over volumetric neighbours through macroEl's face
		TPZGeoElSide ivolume = macroElfaceside.Neighbour();
		for( ; ivolume != macroElfaceside; ivolume = ivolume.Neighbour()){
			if(ivolume.Element()->Dimension() != 3){continue;}
			// get coordinates for center of volume
			TPZManVector<REAL,3> volumeCenter(3);
			{
				TPZGeoElSide gelsidevolume (ivolume.Element(),ivolume.Element()->NSides()-1);
				gelsidevolume.CenterX(volumeCenter);
			}
			// construct vector from center of ifracface to center of volume
			TPZManVector<REAL,3> v2(3,0);
				v2[0] = volumeCenter[0] - macroElCenter[0];
				v2[1] = volumeCenter[1] - macroElCenter[1];
				v2[2] = volumeCenter[2] - macroElCenter[2];
			// Normalize v2
			norm = sqrtl(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
				v2[0] = v2[0]/norm;
				v2[1] = v2[1]/norm;
				v2[2] = v2[2]/norm;
			// if dot product between the vectors constructed for centers is
			// positive, that volume is a candidate
			REAL dot = 0;
			for(int ico = 0; ico < 3; ico++){dot += v1[ico]*v2[ico];}
			if(dot>0){
				candidates[dot] = ivolume.Element()->Index();
			}
		}
	}
	// return best candidate 
    if(candidates.size() > 0){
        // reverse iterator (rbegin) gives biggest key in map
        int64_t volumeindex = candidates.rbegin()->second;
		// For planes that cut the boundary of the domain (gBigPlanes_Q == true), one more verification is required
		bool candidate_is_encloser = true;
		if(gBigPlanes_Q){
			// verify if centroid of ifracface is inside best candidate parametric domain
			TPZGeoEl *candidate_gel = gmesh->Element(volumeindex);
			TPZManVector<REAL,3> qsi(3,2.);
			candidate_is_encloser = candidate_gel->ComputeXInverse(faceCenter,qsi,1e-3);
		}
		if(candidate_is_encloser){
			fVolumes[volumeindex].SetFaceInVolume(ifracface->Index());
			// std::cout<<"Face #"<<ifracface->Index()<<" \t in volume #"<<volumeindex<<"\n";
			return true;
		}
    }
	if(!gBigPlanes_Q){
    	std::cout<<"\n "<<__PRETTY_FUNCTION__<<" found no enclosing volume for element #"<<ifracface->Index()<<"\n";
		DebugStop();
	}
	// DeleteFamily(ifracface);
	CropExternalElement(ifracface);
	// gmesh->DeleteElement(ifracface,ifracface->Index());
    return false;
}











/**
 * 	@brief Creates a .geo for the mesh
 */ 
void DFNMesh::ExportGMshCAD(std::string filename){
	// gmsh doesn't accept index zero elements
	const int shift = 1;
	// materials
    int mtransition = 18;
    int msurface = 40;
    int mintact = 1;
	std::ofstream outfile(filename);

    TPZGeoMesh *pzgmesh = this->Mesh();
	// Number of elements in the mesh
    int64_t nels = pzgmesh->NElements();
	// Dimension of the mesh
	int mesh_dim = pzgmesh->Dimension();

    // Giving fGMesh another name for readability's sake
    pzgmesh->BuildConnectivity();
    CreateSkeletonElements(1,DFNMaterial::Erefined);
    // Title
    outfile<<"//  Geo file generated by Discrete Fracture Network methods \n"
            <<"// Fracture #1 \n\n";
    
    // write nodes
    outfile<< "// POINTS DEFINITION \n\n";
    int64_t nnodes = pzgmesh->NNodes();
    // @ToDo Do we need physical groups for points too?
    for (int64_t inode = 0; inode < nnodes; inode++){
		if(pzgmesh->NodeVec()[inode].Id() < 0) continue;
        TPZManVector<REAL, 3> co(3,0.);
        pzgmesh->NodeVec()[inode].GetCoordinates(co);
        outfile << "Point(" << inode+shift << ") = {" << co[0] << ',' << co[1] << ',' << co[2] << "};\n";
    }
    
    // write edges
    outfile << "\n\n// LINES DEFINITION \n\n";
    {
        // declare lists to define physical groups
        std::list<int64_t> groupSurface;
        std::list<int64_t> groupTransition;
        std::list<int64_t> groupIntact;
        // iterate over all 1D elements
        for (int64_t iel = 0; iel < nels; iel++){
            TPZGeoEl *gel = pzgmesh->Element(iel);
			if(!gel) continue;
            if(gel->Dimension() != 1) continue;
            if(gel->HasSubElement()) continue;
            // Only 2 types of lines are needed:
			// fracture lines
			// i.e. lines with material id of fracture
			// @todo if(gel->MaterialId() == msurface){
			bool entered_Q = false;
			if(gel->MaterialId() >= msurface){
				outfile << "Line(" << iel+shift << ") = {" << gel->NodeIndex(0)+shift << ',' << gel->NodeIndex(1)+shift << "};\n";
				entered_Q = true;
			}
			// edges of macro elements
			// i.e. lines whose elder is an edge of an element with no father
			if(entered_Q == false){
				TPZGeoEl *elder = gel;
				if(gel->Father()) {elder = gel->EldestAncestor();}
				TPZGeoElSide elderside(elder,2);
				TPZGeoElSide neighbour = elderside.Neighbour();
				while(neighbour != elderside){
					// if(neighbour.Element()->Dimension() == 3){ // @todo[see next line] this condition is great but won't work if we're in 2D
					if(!neighbour.Element()->Father()){ // @todo this should fix it
            			outfile << "Line(" << iel+shift << ") = {" << gel->NodeIndex(0)+shift << ',' << gel->NodeIndex(1)+shift << "};\n";
						entered_Q = true;
						break;
					}
					neighbour = neighbour.Neighbour();
				}
			// if(entered_Q) continue;
			}
			
            // list it according to material
			if(!entered_Q) continue;
			if(mesh_dim == 2){
				if(gel->MaterialId() >= mtransition && gel->MaterialId() < msurface){groupTransition.push_back(iel+shift);}
				else if(gel->MaterialId() >= msurface){groupSurface.push_back(iel+shift);}
				else groupIntact.push_back(iel+shift);
				// @todo
				// if(gel->MaterialId() == mtransition) groupTransition.push_back(iel+shift);
				// else if (gel->MaterialId() == mtransition) groupTransition.push_back(iel+shift);
				// else groupIntact.push_back(iel+shift);
			}
        }
        // write physical groups
		if(mesh_dim == 2){
			outfile<<"\nPhysical Curve("<<msurface<<") = {";
			for(auto itr = groupSurface.begin(); itr != groupSurface.end();/*Increment in next line*/){
				outfile<<*itr<<(++itr!=groupSurface.end()? "," : "};\n");
			}
		}
    }
    // write faces
    outfile << "\n\n// FACES DEFINITION \n\n";
    {
        // declare lists to define physical groups
        std::list<int64_t> groupSurface;
        std::list<int64_t> groupTransition;
        std::list<int64_t> groupIntact;
        // iterate over all 2D elements
        for (int64_t iel = 0; iel < nels; iel++){
            TPZGeoEl *gel = pzgmesh->Element(iel);
			if(!gel) continue;
            if(gel->Dimension() != 2) continue;

			// 2 types of faces should enter
			int type = 0;

			// faces of the fracture surface
			do{//pseudo do-while just so I can use break to skip code
				if(mesh_dim==2) break;
				if(gel->HasSubElement()) break;
				// @todo if(gel->MaterialId() != msurface) break;
				if(gel->MaterialId() < msurface) break;
				TPZGeoEl *elder = gel;
				if(gel->Father()) elder = gel->EldestAncestor();
				int nsides = elder->NSides();
				TPZGeoElSide elderside(elder,nsides-1);
				TPZGeoElSide neighbour = elderside.Neighbour();
				if(elderside == neighbour) type = 1;
			} while (false); //pseudo do-while just so I can use break to skip code

			// faces of coarse elements
			do{//pseudo do-while just so I can use break to skip code
				if(type == 1) break;
				if(mesh_dim==2 && !gel->Father()){type = 2; break;}
				if(gel->Father()) break;
				// if(gel->MaterialId() == mintact) {type = 2; break;}
				int nsides = gel->NSides();
				TPZGeoElSide gelside(gel,nsides-1);
				TPZGeoElSide neighbour = gelside.Neighbour();
				while(neighbour != gelside){
					if(neighbour.Element()->Dimension() == 3){
						type = 2;
						break;
					}
					neighbour = neighbour.Neighbour();
				} 
			} while (false);//pseudo do-while just so I can use break to skip code

            if(!type) continue;
			bool planesurface_Q = false;
            // curve loop _______________________________________
            outfile << "Curve Loop(" << iel+shift << ") = {";
            int nnodes = gel->NCornerNodes();
            int nedges = nnodes; //for readability 
            TPZManVector<int64_t,4> facenodevec(nnodes);
            gel->GetNodeIndices(facenodevec);
			int64_t reference_node = facenodevec[0];
            // curve loops require a proper orientation of lines
            for(int iside = nedges; iside < 2*nedges; iside++){
                TPZGeoElSide gelside(gel,iside);
                TPZGeoElSide side = gelside.Neighbour();
                // find line element
                while(side.Element()->Dimension() != 1) {side = side.Neighbour();}
				TPZGeoEl *irib = side.Element();
				TPZStack<TPZGeoEl*> children;
				if(irib->HasSubElement()){
					irib->YoungestChildren(children);
				}else{
					children.Push(irib);
				}
				int nchildren = children.NElements();
				if(nchildren > 1) planesurface_Q = 1;
				int children_added = 0;
                int64_t index = 0;
				//@todo for an unfortunate quick improvisation I had ribs refined into a rib and a point when fractures pass through the node, but this will be fixed soon
				{// @todo temporary exception
					TPZStack<TPZGeoEl*> aux_children;
					for(auto child : children){
						if(child->Dimension() != 0){
							aux_children.Push(child);
						}else{nchildren--;}
					}
					children = aux_children;
				}
				TPZManVector<bool> added_list(nchildren,false);
				while(children_added < nchildren){
					// find next child and check orientation
					int orientation = 1;
					for(int ichild = 0; ichild<nchildren; ichild++){
						if(added_list[ichild]) continue;
						TPZGeoEl *child = children[ichild];
						if(child->NodeIndex(0) == reference_node){
							orientation = 1;
							reference_node = child->NodeIndex(1); // next node becomes reference
							index = child->Index();
							added_list[ichild] = true;
							break;
						}
						else if(child->NodeIndex(1) == reference_node){
							orientation = -1;
							reference_node = child->NodeIndex(0); // next node becomes reference
							index = child->Index();
							added_list[ichild] = true;
							break;
						}
					}
					// add child
					children_added++;
					index = orientation*(index+shift);
                	outfile << index <<(iside + children_added < 2*nedges-1 + nchildren ? "," : "};\n");
				}
            } // close curve loop _______________________________________
            // surface
			if(planesurface_Q) outfile << "Plane ";
            outfile << "Surface("<<iel+shift<<") = {"<<iel+shift<<"};\n";
            if(gel->MaterialId() >= mtransition && gel->MaterialId() < msurface){groupTransition.push_back(iel+shift);}
            else if(gel->MaterialId() >= msurface){groupSurface.push_back(iel+shift);}
            else groupIntact.push_back(iel+shift);
        }
        // write physical groups
        outfile<<"\nPhysical Surface("<<mtransition<<") = {";
        for(auto itr = groupTransition.begin(); itr != groupTransition.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupTransition.end()? "," : "};\n");
        }
        if(mesh_dim == 3) outfile<<"\nPhysical Surface("<<msurface<<") = {";
        for(auto itr = groupSurface.begin(); itr != groupSurface.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupSurface.end()? "," : "};\n");
        }
        outfile<<"\nPhysical Surface("<<mintact<<") = {";
        for(auto itr = groupIntact.begin(); itr != groupIntact.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupIntact.end()? "," : "};\n");
        }
    }

	// Line {x} In Surface {y}
	{
		outfile << "\n\n// EMBEDDED CURVES \n\n";
		for(int iel = 0; iel < nels; iel++){
			TPZGeoEl *gel = pzgmesh->Element(iel);
			if(!gel) continue;
			if(gel->Dimension() != 2) continue;
			if(gel->Father()) continue;
			if(!gel->HasSubElement()) continue;
			if(mesh_dim == 3){
				int nsides = gel->NSides();
				TPZGeoElSide gelside(gel,nsides-1);
				TPZGeoElSide neighbour = gelside.Neighbour();
				if(neighbour == gelside) continue;
			}

			TPZStack<TPZGeoEl*> children;
			gel->YoungestChildren(children);
			// queue all possible ribs using youngest children
			std::set<TPZGeoEl*> candidate_ribs;
			for(TPZGeoEl* child : children){
				int nribs = child->NNodes();
				for(int cside = nribs; cside < 2*nribs; cside++){
					TPZGeoElSide childside(child,cside);
					TPZGeoElSide neig = childside.Neighbour();
					for(/*void*/; neig != childside; neig = neig.Neighbour()){
						if(neig.Element()->Dimension() != 1) continue;
						// @todo if(neig.Element()->MaterialId() != msurface) continue;
						if(neig.Element()->MaterialId() < msurface) continue;
						candidate_ribs.insert(neig.Element());
					}
				}
			}
			outfile << "Curve{";
			bool should_enter = true;
			auto end = candidate_ribs.end();
			for(auto it = candidate_ribs.begin(); it != end;/*void*/){
				TPZGeoEl *rib = *it;
				TPZGeoEl *elder = (rib->Father()? rib->EldestAncestor() : rib);
				TPZGeoElSide elderside(elder,2);
				TPZGeoElSide neig = elderside.Neighbour();
				while(neig != elderside){
					if(neig.Element()->Dimension() == mesh_dim && !neig.Element()->Father()){
						should_enter = false;
						break;
					}
					neig = neig.Neighbour();
				}
				if(!should_enter) {++it;continue;}
				outfile<<rib->Index()+shift<<(++it != end?",":"}");
			}
			outfile << " In Surface{"<<iel+shift<<"};\n";
		}

	}
	

    // write volumes
    if(mesh_dim == 3){
    	outfile << "\n\n// VOLUMES DEFINITION \n\n";
        // declare lists to define physical groups
        std::list<int64_t> groupSurface;
        std::list<int64_t> groupTransition;
        std::list<int64_t> groupIntact;
        // iterate over all 3D elements
        for (int64_t iel = 0; iel < nels; iel++){
            TPZGeoEl *gel = pzgmesh->Element(iel);
			if(!gel) continue;
            if(gel->Dimension() != 3) continue;
            // if(gel->HasSubElement()) continue; // redundant for now, but leaving it here just in case

            // Surface loop
            // gmsh doesn't accept zero index elements
            outfile << "Surface Loop(" << iel+shift << ") = {";

            // iterate over 2D sides to look for faces that close the surface loop
            int nnodes = gel->NCornerNodes();
            int nsides = gel->NSides();
            bool volumeIsCut = false;
            for(int iside = nnodes; iside < nsides-1; iside++){
                if(gel->SideDimension(iside) != 2) continue;
                // find face element
                TPZGeoElSide gelside(gel,iside);
                TPZGeoElSide side = gelside.Neighbour();
                while(side.Element()->Dimension() != 2) {side = side.Neighbour();}
				outfile << side.Element()->Index()+shift << (iside < nsides-2? "," : "};\n");
				if(side.Element()->HasSubElement()) volumeIsCut = true;
            }

            // volume
            outfile << "Volume("<< iel+shift << ") = {"<< iel+shift <<"};\n"; /* gmsh doesn't accept zero index elements */
            if(volumeIsCut){groupTransition.push_back(iel);}
            else groupIntact.push_back(iel);

            if(volumeIsCut){
                int nsurfaces = Volume(iel)->GetFacesInVolume().size();
                TPZManVector<int64_t,6> enclosedSurfaces(nsurfaces);
                enclosedSurfaces = Volume(iel)->GetFacesInVolume();

                // @ToDo may need PrintYoungestChildren here, depending on which elements FindEnclosingVolume is called on
                // -------------------------------------------------------------------------------------
                outfile << "Surface{";
                for(int i = 0; i<nsurfaces; i++){
                    // TPZGeoEl *surface = fGMesh->Element(enclosedSurfaces[i]);
                    // if(surface->HasSubElement()){
                    //     PrintYoungestChildren(surface,outfile);
                    // }
                    // else{
                    //     outfile << enclosedSurfaces[i] << (i<nsurfaces-1?",":"} ");
                    // }
                    outfile << enclosedSurfaces[i]+shift << (i<nsurfaces-1?",":"} ");
                }
                // -------------------------------------------------------------------------------------
                outfile << "In Volume{"<< iel+shift <<"};\n"; /* gmsh doesn't accept zero index elements */
            }
        }
        // write physical groups
        outfile<<"\nPhysical Volume("<<mtransition<<") = {";
        for(auto itr = groupTransition.begin(); itr != groupTransition.end();/*Increment in next line*/){
            outfile<<*itr+shift<<(++itr!=groupTransition.end()? "," : "};\n"); /* gmsh doesn't accept zero index elements */
        }

		if(groupIntact.size() != 0){
			outfile<<"\nPhysical Volume("<<mintact<<") = {";
			for(auto itr = groupIntact.begin(); itr != groupIntact.end();/*Increment in next line*/){
				outfile<<*itr+shift<<(++itr!=groupIntact.end()? "," : "};\n"); /* gmsh doesn't accept zero index elements */
			}
    		outfile<<"\nTransfinite Volume {Physical Volume("<<mintact<<")};\n";
		}
    }
    
	// outfile<<"\nTransfinite Curve {:} = 2;\n";
    outfile<<"Transfinite Surface {Physical Surface("<<mintact<<")};\n";
    outfile<<"Recombine Surface {Physical Surface("<<mintact<<")};\n";
    // outfile<<"Recombine Surface {Physical Surface("<<mtransition<<")};\n";

}



















// /**
//  *  @brief Navigate children tree to access most extreme branches
//  *  @param gel pointer to geometric element of eldest ancestor
//  *  @param outfile ofstream in which to write accessed data
//  */
// void DFNMesh::PrintYoungestChildren(TPZGeoEl *gel, std::ofstream &outfile){
    
//     int nchildren = gel->NSubElements();
//     for(int i = 0; i<nchildren; i++){
//         TPZGeoEl *ichild = gel->SubElement(i);
//         if(ichild->HasSubElement()){
//             PrintYoungestChildren(ichild,outfile);
//         }else{
//             outfile << ichild->Index() +1;//+shift
//         }
// 		outfile << (i < nchildren-1? "," : "");
//     }
// }




























void DFNMesh::CreateVolumes(){
    TPZGeoMesh *gmesh = this->Mesh();
	
	int dim = gmesh->Dimension();
	if(dim != 3) return; //@todo temporary while DFNVolume isn't generalized to 2D meshes
    // map all volumes that are cut
    int64_t nels = gmesh->NElements();
	for (int64_t iel = 0; iel < nels; iel++){
        TPZGeoEl *gel = gmesh->Element(iel);
        if(gel->Dimension() != dim){continue;}
        int nsides = gel->NSides();
        // int ncorners = gel->NCornerNodes();
		// int nfaces = (int) (nsides+1-2*ncorners)/2 //from Euler's characteristic
        for (int iside = nsides-2; iside > 0; iside--){
            TPZGeoElSide gelside(gel,iside);
            if (gelside.Dimension() != 2){break;}
            TPZGeoElSide neighbour = gelside.Neighbour();
			while(neighbour.Element()->Dimension() != 2 && neighbour != gelside){
				neighbour = neighbour.Neighbour();
			}
			if(neighbour == gelside) continue;
			TPZGeoEl *sideface = neighbour.Element();
			if(sideface->HasSubElement()){
                DFNVolume volume(iel,true);
                fVolumes.insert({iel,volume});
                break;
            }            
        }
    }
    
    // gmesh->BuildConnectivity(); //@todo remove this after test
	// search through each 2D element of the triangulated fractures surfaces to find their enclosing volume
	int surfaceMaterial = (*fFractures.begin())->GetSurfaceMaterial();
	for(int64_t iel = 0; iel < nels; iel++){
		TPZGeoEl *gel = gmesh->Element(iel);
		if(!gel) continue;
		// if(gel->MaterialID() != surfaceMaterial){continue;}
		// During development, elements at fracture surface have material id over 40
		if(gel->MaterialId() <= surfaceMaterial) continue;
		if(gel->Dimension() != 2) continue;
		if(gel->HasSubElement()) continue;
		// Find volume that encloses that element
		FindEnclosingVolume(gel);
	}

}































/**
 * 	@brief Uses GMsh API to tetrahedralize a DFNVolume
 */ 
void DFNMesh::Tetrahedralize(DFNVolume *volume){
	// GMsh doesn't like zero index entities
    const int shift = 1;
	int surfaceMaterial = 40;
	int transitionMaterial = 18;

	TPZGeoMesh *gmesh = fGMesh;
	int mesh_dim = gmesh->Dimension();
	int64_t ivol = volume->ElementIndex();
	TPZGeoEl *volGel = gmesh->Element(ivol);
	// std::cout<<"\n\n _________   volume # "<<ivol<<"\n";
	// List faces that form the volume shell
	std::vector<int> surfaceloop;
	for(int nsides = volGel->NSides(),
					 iside = nsides-2; iside >= 0; iside--){
		TPZGeoElSide gelside(volGel,iside);
		if(gelside.Dimension() != 2) break;
		TPZGeoElSide neig = gelside.Neighbour();
		while(neig.Element()->Dimension() != 2){ 
			neig = neig.Neighbour();
			if(neig == gelside){
				PZError << "\n\n Error at "<<__PRETTY_FUNCTION__<<"\n"<<"There shouldn't be a volume without skeleton\n\n";
				DebugStop();
			}
		}
		surfaceloop.push_back(neig.Element()->Index());
	}
	// List faces that are enclosed in the volume
	TPZManVector<int64_t> enclosedFaces = volume->GetFacesInVolume();

	// List all lines
	std::set<int64_t> lines;
	// List edges from volume shell
	for(int iface = 0,
			nfaces = surfaceloop.size(); iface < nfaces; iface++){
		TPZGeoEl *gel = gmesh->Element(surfaceloop[iface]);
		for(int nsides = gel->NSides(),
				nnodes = gel->NCornerNodes(),
				iside = nsides-2; iside >= nnodes; iside--){
			TPZGeoElSide gelside(gel,iside);
			if(gelside.Dimension() != 1) break;
			TPZGeoElSide neig = gelside.Neighbour();
			while(neig.Element()->Dimension() != 1) neig = neig.Neighbour();
			TPZStack<TPZGeoEl*> unrefined_lines;
			if(neig.Element()->HasSubElement()){
				neig.Element()->YoungestChildren(unrefined_lines);
			}else{
				unrefined_lines.Push(neig.Element());
			}
			for(auto line : unrefined_lines){
				lines.insert(line->Index());
			}
		}
	}
	// List edges from faces in volume
	for(int iface = 0,
			nfaces = enclosedFaces.size(); iface < nfaces; iface++){
		TPZGeoEl *gel = gmesh->Element(enclosedFaces[iface]);
		for(int nsides = gel->NSides(),
				nnodes = gel->NCornerNodes(),
				iside = nsides-2; iside >= nnodes; iside--){
			TPZGeoElSide gelside(gel,iside);
			if(gelside.Dimension() != 1) break;
			TPZGeoElSide neig = gelside.Neighbour();
			while(neig.Element()->Dimension() != 1) neig = neig.Neighbour();
			TPZStack<TPZGeoEl*> unrefined_lines;
			if(neig.Element()->HasSubElement()){
				neig.Element()->YoungestChildren(unrefined_lines);
			}else{
				unrefined_lines.Push(neig.Element());
			}
			for(auto line : unrefined_lines){
				lines.insert(line->Index());
			}
		}
	}

	// List nodes
	std::set<int64_t> nodes;
	for(int64_t line : lines){
		TPZGeoEl *gel = gmesh->Element(line);
		nodes.insert(gel->NodeIndex(0));
		nodes.insert(gel->NodeIndex(1));
	}
	
	// gmsh::initialize();
	std::string modelname = "model"+std::to_string(ivol);
	gmsh::model::add(modelname);
	std::string mshfilename("LOG/testAPI_volume");
	mshfilename += std::to_string(ivol);
	// gmsh::model::add(mshfilename);
	mshfilename += ".msh";
	gmsh::option::setNumber("Mesh.Algorithm3D",1);  // (1: Delaunay, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT) Default value: 1
	// Insert nodes ____________________________________
	for(int64_t inode : nodes){
		TPZManVector<REAL,3> coord(3);
		gmesh->NodeVec()[inode].GetCoordinates(coord);
		gmsh::model::geo::addPoint(coord[0],coord[1],coord[2],0.,inode+shift);
	}
	// std::cout<<"\n\n\n";
	// Insert lines ____________________________________
	for(int64_t iline : lines){
		TPZGeoEl *gel = gmesh->Element(iline);
		int64_t node0 = gel->NodeIndex(0)+shift;
		int64_t node1 = gel->NodeIndex(1)+shift;
		gmsh::model::geo::addLine(node0,node1,iline+shift);
		gmsh::model::geo::mesh::setTransfiniteCurve(iline+shift,2);
	}
	// Insert faces ____________________________________
	{
	std::vector<int> wiretag(1);
	// Faces in volume shell
	for(int64_t faceindex : surfaceloop){
		TPZGeoEl *face = gmesh->Element(faceindex);
		int nnodes = face->NCornerNodes();
		int nedges = nnodes;
		TPZManVector<int64_t,4> facenodevec(nnodes);
		face->GetNodeIndices(facenodevec);
		// line loop ________________________________________________________________________
		std::vector<int> lineloop;
		int64_t reference_node = face->NodeIndex(0);
		bool planesurface_Q = false;
		for(int iside = nnodes; iside < nnodes+nedges; iside++){
			TPZGeoElSide gelside(face,iside);
			TPZGeoElSide side = gelside.Neighbour();
			// find line element
			while(side.Element()->Dimension() != 1) {side = side.Neighbour();}
			TPZGeoEl *irib = side.Element();
			TPZStack<TPZGeoEl*> children;
			if(irib->HasSubElement()){
				irib->YoungestChildren(children);
			}else{
				children.Push(irib);
			}
			int nchildren = children.NElements();
			if(nchildren > 1) planesurface_Q = 1;
			int children_added = 0;
			int64_t index = 0;
			TPZManVector<bool> added_list(nchildren,false);
			while(children_added < nchildren){
				// find next child and check orientation
				int orientation = 1;
				for(int ichild = 0; ichild<nchildren; ichild++){
					if(added_list[ichild]) continue;
					TPZGeoEl *child = children[ichild];
					if(child->NodeIndex(0) == reference_node){
						orientation = 1;
						reference_node = child->NodeIndex(1); // next node becomes reference
						index = child->Index();
						added_list[ichild] = true;
						break;
					}
					else if(child->NodeIndex(1) == reference_node){
						orientation = -1;
						reference_node = child->NodeIndex(0); // next node becomes reference
						index = child->Index();
						added_list[ichild] = true;
						break;
					}
				}
				// add child
				children_added++;
				index = orientation*(index+shift);
				lineloop.push_back(index);
			}
		}
		// insert curve loop ________________________________________________________________
		wiretag[0] = gmsh::model::geo::addCurveLoop(lineloop,face->Index()+shift);
		// insert surface
		if(planesurface_Q){gmsh::model::geo::addPlaneSurface(wiretag,wiretag[0]);
		}else{gmsh::model::geo::addSurfaceFilling(wiretag,wiretag[0]);}
		
		// @todo gmsh::model::geo::mesh::setTransfiniteSurface(wiretag[0]);
		// if(face->Type() == EQuadrilateral) gmsh::model::geo::mesh::setRecombine(2,wiretag[0]);
	}
	// Enclosed faces
	for(int64_t faceindex : enclosedFaces){
		TPZGeoEl *face = gmesh->Element(faceindex);
		int nnodes = face->NCornerNodes();
		int nedges = nnodes;
		TPZManVector<int64_t,4> facenodevec(nnodes);
		face->GetNodeIndices(facenodevec);
		// line loop
		std::vector<int> lineloop(nedges);
		for(int iside = nnodes; iside < nnodes+nedges; iside++){
			TPZGeoElSide gelside(face,iside);
			TPZGeoElSide neig = gelside.Neighbour();
			// find line element
			while(neig.Element()->Dimension()!=1){neig = neig.Neighbour();}
			// find first node of line at the face
			int inode = 0;
			while(facenodevec[inode] != neig.SideNodeIndex(0)) inode++;
			// check orientation by comparing second node of line with next node of face
			if(neig.SideNodeIndex(1) == facenodevec[(inode+1)%nnodes]){
				lineloop[iside-nnodes] = neig.Element()->Index() + shift;
			}else{
				lineloop[iside-nnodes] = - (neig.Element()->Index() + shift);
			}
		}
		// insert curve loop
		wiretag[0] = gmsh::model::geo::addCurveLoop(lineloop,face->Index()+shift);
		// insert surface
		gmsh::model::geo::addSurfaceFilling(wiretag,wiretag[0]);	
		gmsh::model::geo::mesh::setTransfiniteSurface(wiretag[0]);
		// if(face->Type() == EQuadrilateral) gmsh::model::geo::mesh::setRecombine(2,wiretag[0]);
	}
	}
	
	{// begin Lines-in-Surface ________________________________________________
	// "lines of fracture surface whose elder does not have a mesh_dim neighbour without father"
	for(int iface : surfaceloop){
		TPZGeoEl *gel = gmesh->Element(iface);
		if(!gel->HasSubElement()) continue;
		TPZStack<TPZGeoEl *> children;
		gel->YoungestChildren(children);
		// queue all possible lines by checking 1D neighbours of children
		std::set<TPZGeoEl *> candidate_ribs;
		for(auto child : children){
			int nribs = child->NNodes();
			for(int cside = nribs; cside < 2*nribs; cside++){
				TPZGeoElSide childside(child,cside);
				TPZGeoElSide neig = childside.Neighbour();
				for(/*void*/; neig != childside; neig = neig.Neighbour()){
					if(neig.Element()->Dimension() != 1) continue;
					// @todo if(neig.Element()->MaterialId() != surfaceMaterial) continue;
					if(neig.Element()->MaterialId() < surfaceMaterial) continue;
					candidate_ribs.insert(neig.Element());
				}
			}
		}
		bool should_enter = true;
		std::vector<int> lines_in_surface;
		auto end = candidate_ribs.end();
		for(auto it = candidate_ribs.begin(); it != end;it++){
			TPZGeoEl *rib = *it;
			TPZGeoEl *elder = (rib->Father()? rib->EldestAncestor() : rib);
			TPZGeoElSide elderside(elder,2);
			TPZGeoElSide neig = elderside.Neighbour();
			while(neig != elderside){
				if(neig.Element()->Dimension() == mesh_dim && !neig.Element()->Father()){
					should_enter = false;
					break;
				}
				neig = neig.Neighbour();
			}
			if(!should_enter){continue;}
			lines_in_surface.push_back(rib->Index()+shift);
		}
		gmsh::model::geo::synchronize(); // synchronize is required before embedding
		gmsh::model::mesh::embed(1,lines_in_surface,2,iface+shift);
	}
	}// end Lines-in-Surface __________________________________________________


	// Insert volumes ____________________________________
	std::vector<int> shelltag(1);
	// Shift surfaceloop indices
	for(int nsurfaces = surfaceloop.size(),
			i = 0; i < nsurfaces; i++){
		surfaceloop[i] += shift;
	}
	shelltag[0] = gmsh::model::geo::addSurfaceLoop(surfaceloop,ivol+shift);
	gmsh::model::geo::addVolume(shelltag,shelltag[0]);
	
	// Surfaces in Volume ________________________________
		gmsh::model::geo::synchronize(); // synchronize is required before embedding
		int nfacesenclosed = enclosedFaces.size();
		std::vector<int> facesinvolume(nfacesenclosed);
		{
			for(int i = 0; i<nfacesenclosed; i++){
				facesinvolume[i] = enclosedFaces[i] + shift;
			}
			gmsh::model::mesh::embed(2, facesinvolume, 3, ivol+shift);
		}
	// Physical groups ____________________________
		gmsh::model::addPhysicalGroup(2,facesinvolume,surfaceMaterial);
		gmsh::model::addPhysicalGroup(2,surfaceloop,transitionMaterial);
		gmsh::model::addPhysicalGroup(3,shelltag,transitionMaterial);
	
	// synchronize before meshing
		gmsh::model::geo::synchronize();
	// mesh
		gmsh::model::mesh::generate(3);
	gmsh::write(mshfilename);
	// gmsh::write("LOG/testAPI_volume.msh");
	// import meshed volume back into PZ geoMesh
		ImportElementsFromGMSH(gmesh,3,nodes);
	gmsh::model::remove();
	// gmsh::clear();
	// gmsh::finalize();
}















/**
 * @brief Check if the neighbour has equal dimension
 * @param geliside GeoElement side
 */

bool DFNMesh::HasEqualDimensionNeighbour(TPZGeoElSide &gelside){
    
    int dimension = gelside.Dimension();

    if (gelside.Element()->Dimension() == dimension){
        return true;
    }

    TPZGeoElSide neighbour = gelside.Neighbour();

    while (neighbour != gelside){
        if (neighbour.Element()->Dimension()==dimension){
            return true;
        }
        neighbour = neighbour.Neighbour();
    }
    return false;    
}










 /**
  * @brief Creates the skeleton mesh
  * @param dimension
  */

void DFNMesh::CreateSkeletonElements(int dimension, int matid)
{
    int nel = fGMesh->NElements();
    for (int iel = 0; iel < nel; iel++)
    {
        TPZGeoEl *gel = fGMesh->Element(iel);
        if(!gel) continue;
        // Elements can't have a skeleton of higher dimension than itself
        if(gel->Dimension() <= dimension) continue;
        
        int nsides = gel->NSides();
        int ncorners = gel->NCornerNodes();
        // iterating from higher-dimensional sides to lower-dimensional should narrow the search
        for (int iside = nsides-2; iside >= ncorners; iside--)
        {
            TPZGeoElSide gelside = gel->Neighbour(iside);

            if (gelside.Dimension() != dimension){continue;}
            bool haskel = HasEqualDimensionNeighbour(gelside);
            if (haskel == false)
            {
                if(matid == -1) matid = gel->MaterialId();
                // TPZGeoElBC(gelside, matid);
                if(gel->MaterialId() == DFNMaterial::Efracture) TPZGeoElBC(gelside, DFNMaterial::Efracture);
                else TPZGeoElBC(gelside, matid);
            }
        }
    }
}








