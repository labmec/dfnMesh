/*! 
 *  @brief     Contains implementation of methods for class DFNMesh.
 *  @authors   Pedro Lima
 *  @date      2020.06
 */

#include "DFNMesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZGeoMeshBuilder.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("dfn.mesh"));
#endif

DFNMesh::DFNMesh(TPZGeoMesh *gmesh, REAL tolDist, REAL tolAngle){
	fGMesh = gmesh;
    fTolDist = tolDist;
    fTolAngle = tolAngle;
    fTolAngle_cos = std::cos(tolAngle);

	// create a copy of the mesh because materials will be reset afterwards @todo
	this->ClearMaterials();
	for (int i = 1; i < fGMesh->Dimension(); i++){
		CreateSkeletonElements(i,1);
	}
	
	fSortedFaces.clear();
	fPolyh_per_face.clear();
	fPolyhedra.clear();
	fFractures.clear();
	fVolumes.clear();

	if(fGMesh->Dimension() == 3){
		InitializePolyhedra();
		// BuildPolyhedra_firstrun();
	}
}


/// Copy constructor
DFNMesh::DFNMesh(const DFNMesh &copy){
    this->operator=(copy);
}

/// Assignment operator
DFNMesh &DFNMesh::operator=(const DFNMesh &copy){
    fFractures = 		copy.fFractures;
    fTolDist = 			copy.fTolDist;
    fTolAngle = 		copy.fTolAngle;
    fTolAngle_cos = 	copy.fTolAngle_cos;
    fGMesh = 			copy.fGMesh;
	fSortedFaces = 		copy.fSortedFaces;
    fPolyh_per_face = 	copy.fPolyh_per_face;
	fPolyhedra = 		copy.fPolyhedra;
	fVolumes = 			copy.fVolumes;
    return *this;
}

void DFNMesh::PrintVTK(std::string pzmesh
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








void DFNMesh::PrintVTKColorful(std::string pzmesh,std::string vtkmesh){
	// mesh.txt doesn't gain much... so print it normal first
	this->PrintVTK(pzmesh,"skip");
	TPZGeoMesh *gmesh = this->fGMesh;
	int64_t nels = gmesh->NElements();
	int64_t iel;
	TPZGeoEl *gel;
	// shift material ids
	for(iel = 0; iel < nels; iel++){
		gel = gmesh->Element(iel);
		if(!gel) continue;
		if(gel->Dimension() != 2) continue;
		if(gel->HasSubElement()) continue;
		if(!gel->Father()) continue;
		if(gel->MaterialId() == DFNMaterial::Efracture) continue;
		if(gel->MaterialId() == DFNMaterial::Eintact) continue;
		int subindex = gel->WhichSubel();
		int matid = gel->MaterialId();
		gel->SetMaterialId(DFNMaterial::Erefined+subindex);
	}
	// print vtk only, since txt has already been print
	this->PrintVTK("skip",vtkmesh);

	// then, restore original mat ids
	for(iel = 0; iel < nels; iel++){
		gel = gmesh->Element(iel);
		if(!gel) continue;
		if(gel->Dimension() != 2) continue;
		if(gel->HasSubElement()) continue;
		if(!gel->Father()) continue;
		int subindex = gel->WhichSubel();
		int matid = gel->MaterialId();
		gel->SetMaterialId(matid-subindex);
	}
}


void InwardNormal(const TPZGeoElSide& gelside, TPZManVector<REAL,3>& qsi, TPZManVector<REAL,3>& normal){
	gelside.Normal(qsi,normal);
	for(int dim=normal.size(), i=0; i<dim; i++){
		normal[i] *= -1;
	}
}




/** @brief adds geoels for graphics that illustrate the tolerance
 * @warning only implemented for 2D right now
*/
void DFNMesh::PlotTolerance(TPZManVector<int64_t>& indices){
	using namespace DFN;
	TPZManVector<int64_t,8> nodeindices(8,-1);
	TPZGeoEl* gel = nullptr;
	TPZManVector<REAL,3> orig_coord(3,0.);
	TPZManVector<REAL,3> clone_coord(3,0.);

	TPZStack<TPZManVector<REAL,3>,4> Side_inNormals(4,{0.,0.,0.});

	for(int64_t iel : indices){
		gel = fGMesh->Element(iel);
		if(gel->Dimension() != 2) DebugStop();
		int sidedim = gel->Dimension()-1;

		int nnodes = gel->NCornerNodes();
		int nedges = gel->NSides(1);
		nodeindices.Fill(-1);
		nodeindices.Resize(nnodes);
		
		// Get inward normal vector for every edge and scale them by the tolerable distance
		Side_inNormals.resize(nedges);
		for(int iside=nnodes, iedge=0; iside<gel->NSides()-1; iside++, iedge++){
			TPZGeoElSide gelside(gel,iside);
			TPZManVector<REAL,3> qsi(gelside.Dimension(),0.);
			InwardNormal(gelside,qsi,Side_inNormals[iedge]);
			for(int idim=0; idim<3; idim++)
				{Side_inNormals[iedge][idim] *= fTolDist;}
			Side_inNormals[iedge][2] += gDFN_SmallNumber;
		}

		// Clone nodes and move inward by the tolerable distance 
		for(int i=0; i<nnodes; i++){
			nodeindices[i] = fGMesh->NodeVec().AllocateNewElement();
			gel->Node(i).GetCoordinates(orig_coord);
			// desloc from posterior and anterior edges
			for(int idim=0; idim<3; idim++){
				clone_coord[idim] = orig_coord[idim] 
									+Side_inNormals[i][idim]
									+Side_inNormals[(i+nedges-1)%nedges][idim];
			}
			fGMesh->NodeVec()[nodeindices[i]].Initialize(clone_coord,*fGMesh);
		}

		// Create clone geoel
		MElementType etype = gel->Type();
		int64_t index = -1;
		fGMesh->CreateGeoElement(etype,nodeindices,-1,index);
	}
}





// void DFNMesh::GenerateSubMesh(){
//     TPZGeoMesh *gmesh = this->Mesh();
// 	int dim = gmesh->Dimension();
// 	if(dim != 3) return; //@todo temporary while DFNVolume isn't generalized to 2D meshes
	
// 	//Loop over list of fractured volumes
// 	for (auto itr = fVolumes.begin(); itr != fVolumes.end(); itr++){
//     	DFNVolume *ivolume = &itr->second;
// 		// if(ivolume->ElementIndex() == 0) continue;
// 		// Use GMsh to tetrahedralize volumes
// 		gmsh::logger::start();

//     	Tetrahedralize(ivolume);
		
// 		std::vector<std::string> logvec;
// 		gmsh::logger::get(logvec);
// 		for(std::string str : logvec){
// 			std::cout << str << std::endl;
// 		}
// 		gmsh::logger::stop();
// 	}
// }



// /**
//  * @brief Deletes face + ribs (if left isolated) + nodes (if left isolated)
//  * @param face: 2D element to be deleted
// */
// void DFNMesh::DeleteElementAndRibs(TPZGeoEl *face){
// 	TPZGeoMesh *gmesh = face->Mesh();
// 	// queue ribs
// 	int nribs = face->NCornerNodes();
// 	// int nribs = (face->Type() == MElementType::EQuadrilateral ? 4 : 3);
// 	TPZManVector<int64_t, 4> ribs(nribs,-1);
// 	for(int irib = 0; irib < nribs; irib++){
// 		TPZGeoElSide faceside(face,irib+nribs);
// 		for(TPZGeoElSide neighbour = faceside.Neighbour(); neighbour != faceside; neighbour = neighbour.Neighbour()){
// 			if(neighbour.Element()->Dimension() == 1){
// 				ribs[irib] = neighbour.Element()->Index(); 
// 				break;
// 			}
// 		}
// 	}
// 	// delete face
// 	gmesh->DeleteElement(face);
// 	// then, delete face's ribs and nodes if necessary
// 	for(int irib = 0; irib < nribs; irib++){
// 		if(ribs[irib] == -1) continue;
// 		TPZGeoEl *ribgel = gmesh->Element(ribs[irib]);
// 		TPZManVector<int64_t,2> nodes(2,-1);
// 		for(int iside = 0; iside < 3 ; iside++){
// 			TPZGeoElSide ribgelside(ribgel,iside);
// 			TPZGeoElSide neighbour = ribgelside.Neighbour();
// 			bool delete_Q = true;
// 			// @todo maybe add exception for 0D neighbour
// 			if(iside < 2){
// 				if(neighbour == ribgelside){
// 					nodes[iside] = ribgel->NodeIndex(iside);
// 				}
// 			}else{
// 				while(neighbour != ribgelside){
// 					if(neighbour.Element()->Dimension() == 2){
// 						delete_Q = false; 
// 						break;
// 					}
// 					neighbour = neighbour.Neighbour();
// 				}
// 				if(delete_Q){gmesh->DeleteElement(ribgel,ribs[irib]);}
// 			}
// 		}
// 		// delete node
// 		for(int inode : nodes){
// 			if(inode < 0) continue;
// 			gmesh->NodeVec().SetFree(inode);
// 			gmesh->NodeVec()[inode].SetNodeId(-1);
// 			// delete &gmesh->NodeVec()[inode];
// 		}
// 	}
// }









// /**
//  * @brief Deletes gel + children + isolated ribs + unused nodes
//  * @note It will assume element has been found not to belong to the domain of interest and will not verify
//  * @param gel: pointer to the geometric element
// */
// void DFNMesh::CropExternalElement(TPZGeoEl *gel){
// 	// return; //debugging
// 	TPZGeoMesh *gmesh = gel->Mesh();
// 	// If an element is external to the domain, then its eldest ancestor and all the refinement tree are also external
// 	TPZGeoEl *elder = gel;
// 	if(gel->Father()){elder = gel->EldestAncestor();}
// 	// Start from youngest children and go up the tree
// 	while(elder->HasSubElement()){
// 		TPZStack<TPZGeoEl*> youngestChildren;
// 		elder->YoungestChildren(youngestChildren);
// 		for(auto child : youngestChildren){
// 			DeleteElementAndRibs(child);
// 		}
// 	}
// 	DeleteElementAndRibs(elder);
// }









// /**
//  * @brief Deletes gel and all elements that share the same eldest ancestor, then deletes the ancestor
//  * @param gel: Any member of the family
// */
// void DFNMesh::DeleteFamily(TPZGeoEl *gel){
// 	TPZGeoMesh *gmesh = gel->Mesh();
// 	if(!gel->Father()){
// 		gmesh->DeleteElement(gel);
// 		return;
// 	}
// 	TPZGeoEl *elder = gel->EldestAncestor();
// 	while(elder->HasSubElement()){ //this looks redundant, but I'll need it to delete ribs and nodes eventually
// 		TPZStack<TPZGeoEl*> youngestChildren;
// 		elder->YoungestChildren(youngestChildren);
// 		for(auto child : youngestChildren){
// 			gmesh->DeleteElement(child);
// 		}	
// 	}
// 	gmesh->DeleteElement(elder);
// }






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
                    // int matid = DFNMaterial::Erefined;
					// TPZGeoMeshBuilder::InsertElement(gmesh, matid, el_type, el_identifier, node_identifiers);
					int64_t ntest = gmesh->NElements();
					// std::cout<<"nelements = "<<ntest<<"\n";
                }
            }
        }
    }
    gmesh->BuildConnectivity();
}





























// int64_t DFNMesh::SearchIndirectNeighbours(TPZGeoEl* gel){
    
//     std::list<int64_t> candidate_queue;
// 	std::set<int64_t> verified;
//     candidate_queue.push_back(gel->Index());
//     for (auto index : candidate_queue) {
// 		if(verified.find(index) != verified.end()) continue;
//         TPZGeoEl *currentgel = gel->Mesh()->Element(index);
// 		int64_t macroElindex = FindAdjacentMacroEl(currentgel);
//         if (macroElindex == -1) {
//             QueueNeighbours(currentgel, candidate_queue);
//         }
//         if (macroElindex >= 0) {
//             return macroElindex;;
//         }   
// 		verified.insert(index);
//     }
// 	// If this point is reached, current element is surrounded by elements that have been deleted and, therefore, is not in the domain (and should also be deleted).
// 	return -1;
// }







// int64_t DFNMesh::FindAdjacentMacroEl(TPZGeoEl* gel){
//     int nsides = gel->NSides();
        
// 	for(int iside = nsides-2; iside >= 0; iside--){
// 		TPZGeoElSide gelside(gel, iside);
// 		if(gelside.Dimension() < 1) break;
// 		TPZGeoElSide neig = gelside.Neighbour();
// 		for(/*void*/; neig != gelside; neig = neig.Neighbour()){
// 			if(neig.Element()->Dimension() != 2) continue;
// 			if (DFN::IsInterface(neig.Element())){
// 				if(neig.Element()->Father())
// 					{return neig.Element()->EldestAncestor()->Index();}
// 				else
// 					{return neig.Element()->Index();}
// 			}
// 		}
// 	}
    
//     return -1;
// }




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










bool DFN::IsInterface(TPZGeoEl* gel){
	// gel can only be an interface if it is (d-1)-dimensional
	if(gel->Dimension() != gel->Mesh()->Dimension()-1) {return false;}
 
	TPZGeoEl* elder = nullptr;
	if(gel->Father()){elder = gel->EldestAncestor();}
	else{elder = gel;}
 
	// Return true if eldest ancestor has an orphan neighbour through side nsides-1
	TPZGeoElSide gelside(elder,elder->NSides()-1);
	TPZGeoElSide neig = gelside.Neighbour();
	for(/*void*/; neig != gelside; neig = neig.Neighbour()){
		// if(neig.Element()->Dimension() == 3) {return true;}
		if(!neig.Element()->Father()) {return true;}
	}
 
	return false;
}






// bool DFNMesh::FindEnclosingVolume(TPZGeoEl *ifracface){
// 	bool gBigPlanes_Q = true; //placeholder
// 	if(ifracface->Dimension()!=2) DebugStop();
// 	if(DFN::IsInterface(ifracface)) return false;
// 	int64_t ifracfaceindex = ifracface->Index();
// 	TPZGeoMesh *gmesh = ifracface->Mesh();
//     // get coordinates of geometric center of face
//     TPZVec<REAL> faceCenter(3);
//     {
//         TPZGeoElSide geliside(ifracface, ifracface->NSides()-1);
//         geliside.CenterX(faceCenter);
//     }

//     // map of indices for volumes that could contain the face
//     std::map<REAL, int64_t> candidates;
//     {
// 		int64_t macroElindex = SearchIndirectNeighbours(ifracface);
// 		if(macroElindex == -1){
// 			if(!gBigPlanes_Q){
// 				std::cout<<"\n "<<__PRETTY_FUNCTION__<<" found no enclosing volume for element #"<<ifracface->Index()<<"\n";
// 				DebugStop();
// 			}
// 			CropExternalElement(ifracface);
// 			// gmesh->DeleteElement(ifracface,ifracface->Index());
// 			return false;
// 		}
// 		TPZGeoEl *macroEl = ifracface->Mesh()->Element(macroElindex);
// 		// get macroEl's center coordinates
// 		TPZManVector<REAL,3> macroElCenter(3);
// 		TPZGeoElSide macroElfaceside(macroEl, macroEl->NSides()-1);
// 		macroElfaceside.CenterX(macroElCenter);
// 		// construct vector from center of macroEl to center of ifracface
// 		TPZManVector<REAL,3> v1(3,0);
// 			v1[0] = faceCenter[0] - macroElCenter[0];
// 			v1[1] = faceCenter[1] - macroElCenter[1];
// 			v1[2] = faceCenter[2] - macroElCenter[2];
// 		// Normalize v1
// 		REAL norm = sqrtl(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
// 			v1[0] = v1[0]/norm;
// 			v1[1] = v1[1]/norm;
// 			v1[2] = v1[2]/norm;
// 		// iterate over volumetric neighbours through macroEl's face
// 		TPZGeoElSide ivolume = macroElfaceside.Neighbour();
// 		for( ; ivolume != macroElfaceside; ivolume = ivolume.Neighbour()){
// 			if(ivolume.Element()->Dimension() != 3){continue;}
// 			// get coordinates for center of volume
// 			TPZManVector<REAL,3> volumeCenter(3);
// 			{
// 				TPZGeoElSide gelsidevolume (ivolume.Element(),ivolume.Element()->NSides()-1);
// 				gelsidevolume.CenterX(volumeCenter);
// 			}
// 			// construct vector from center of ifracface to center of volume
// 			TPZManVector<REAL,3> v2(3,0);
// 				v2[0] = volumeCenter[0] - macroElCenter[0];
// 				v2[1] = volumeCenter[1] - macroElCenter[1];
// 				v2[2] = volumeCenter[2] - macroElCenter[2];
// 			// Normalize v2
// 			norm = sqrtl(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
// 				v2[0] = v2[0]/norm;
// 				v2[1] = v2[1]/norm;
// 				v2[2] = v2[2]/norm;
// 			// if dot product between the vectors constructed for centers is
// 			// positive, that volume is a candidate
// 			REAL dot = 0;
// 			for(int ico = 0; ico < 3; ico++){dot += v1[ico]*v2[ico];}
// 			if(dot>0){
// 				candidates[dot] = ivolume.Element()->Index();
// 			}
// 		}
// 	}
// 	// return best candidate 
//     if(candidates.size() > 0){
//         // reverse iterator (rbegin) gives biggest key in map
//         int64_t volumeindex = candidates.rbegin()->second;
// 		// For planes that cut the boundary of the domain (gBigPlanes_Q == true), one more verification is required
// 		bool candidate_is_encloser = true;
// 		if(gBigPlanes_Q){
// 			// verify if centroid of ifracface is inside best candidate parametric domain
// 			TPZGeoEl *candidate_gel = gmesh->Element(volumeindex);
// 			TPZManVector<REAL,3> qsi(3,2.);
// 			candidate_is_encloser = candidate_gel->ComputeXInverse(faceCenter,qsi,1e-3);
// 		}
// 		if(candidate_is_encloser){
// 			fVolumes[volumeindex].SetFaceInVolume(ifracface->Index());
// 			// std::cout<<"Face #"<<ifracface->Index()<<" \t in volume #"<<volumeindex<<"\n";
// 			return true;
// 		}
//     }
// 	if(!gBigPlanes_Q){
//     	std::cout<<"\n "<<__PRETTY_FUNCTION__<<" found no enclosing volume for element #"<<ifracface->Index()<<"\n";
// 		DebugStop();
// 	}
// 	// DeleteFamily(ifracface);
// 	CropExternalElement(ifracface);
// 	// gmesh->DeleteElement(ifracface,ifracface->Index());
//     return false;
// }











/**
 * 	@brief Creates a .geo for the mesh
 */ 
void DFNMesh::ExportGMshCAD(std::string filename){
	DebugStop(); // @todo New Surface meshing has left this method faulty for 3D meshes. I'll have to refactor.
	// gmsh doesn't accept index zero elements
	const int shift = 1;
	// materials
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
			bool entered_Q = false;
			if(gel->MaterialId() == DFNMaterial::Efracture){
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
				if(gel->MaterialId() == DFNMaterial::Erefined){groupTransition.push_back(iel+shift);}
				else if(gel->MaterialId() == DFNMaterial::Efracture){groupSurface.push_back(iel+shift);}
				else groupIntact.push_back(iel+shift);
			}
        }
        // write physical groups
		if(mesh_dim == 2){
			outfile<<"\nPhysical Curve(\"Fractures\","<<DFNMaterial::Efracture<<") = {";
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
				if(gel->MaterialId() != DFNMaterial::Efracture) break;
				// if(gel->MaterialId() < DFNMaterial::Efracture) break;
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
				// if(gel->MaterialId() == DFNMaterial::Eintact) {type = 2; break;}
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
				// //@todo for an unfortunate quick improvisation I had ribs refined into a rib and a point when fractures pass through the node, but this will be fixed soon
				// {// @todo temporary exception
				// 	TPZStack<TPZGeoEl*> aux_children;
				// 	for(auto child : children){
				// 		if(child->Dimension() != 0){
				// 			aux_children.Push(child);
				// 		}else{nchildren--;}
				// 	}
				// 	children = aux_children;
				// }
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
            if(gel->MaterialId() == DFNMaterial::Erefined){groupTransition.push_back(iel+shift);}
            else if(gel->MaterialId() == DFNMaterial::Efracture){groupSurface.push_back(iel+shift);}
            else groupIntact.push_back(iel+shift);
        }
        // write physical groups
        outfile<<"\nPhysical Surface("<<DFNMaterial::Erefined<<") = {";
        for(auto itr = groupTransition.begin(); itr != groupTransition.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupTransition.end()? "," : "};\n");
        }
        if(mesh_dim == 3) outfile<<"\nPhysical Surface(\"Fractures\","<<DFNMaterial::Efracture<<") = {";
        for(auto itr = groupSurface.begin(); itr != groupSurface.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupSurface.end()? "," : "};\n");
        }
        outfile<<"\nPhysical Surface("<<DFNMaterial::Eintact<<") = {";
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
				int nribs = child->NCornerNodes();
				for(int cside = nribs; cside < 2*nribs; cside++){
					TPZGeoElSide childside(child,cside);
					TPZGeoElSide neig = childside.Neighbour();
					for(/*void*/; neig != childside; neig = neig.Neighbour()){
						if(neig.Element()->Dimension() != 1) continue;
						if(neig.Element()->MaterialId() != DFNMaterial::Efracture) continue;
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
                while(side.Element()->Dimension() != 2) {
					side = side.Neighbour();
					if(side == gelside){
						PZError << "3D element without 2D skeleton:\n\t index = " << gel->Index() << " side = "<<iside;
						DebugStop();
					}
				}
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

                outfile << "Surface{";
                for(int i = 0; i<nsurfaces; i++){
                    outfile << enclosedSurfaces[i]+shift << (i<nsurfaces-1?",":"} ");
                }
                outfile << "In Volume{"<< iel+shift <<"};\n"; /* gmsh doesn't accept zero index elements */
            }
        }
        // write physical groups
		if(groupTransition.size() != 0){
			outfile<<"\nrefinedvol[] = {";
			for(auto itr = groupTransition.begin(); itr != groupTransition.end();/*Increment in next line*/){
				outfile<<*itr+shift<<(++itr!=groupTransition.end()? "," : "};\n"); /* gmsh doesn't accept zero index elements */
			}
		}
		if(groupIntact.size() != 0){
			outfile<<"\nintactvol[] = {";
			for(auto itr = groupIntact.begin(); itr != groupIntact.end();/*Increment in next line*/){
				outfile<<*itr+shift<<(++itr!=groupIntact.end()? "," : "};\n"); /* gmsh doesn't accept zero index elements */
			}
    		outfile<<"\nTransfinite Volume {intactvol[]};\n";
		}
		for(int64_t iel=0;iel<fGMesh->NElements(); iel++){
			if(!fGMesh->Element(iel)) continue;
			if(fGMesh->Element(iel)->Dimension() == 3){
				outfile<<"\nPhysical Volume("<<iel+shift<<") = {"<<iel+shift<<"};";
			}
		}
    }
    
	// outfile<<"\nTransfinite Curve {:} = 2;";
    outfile<<"\nTransfinite Surface {Physical Surface("<<DFNMaterial::Eintact<<")};";
    outfile<<"\nRecombine Surface {Physical Surface("<<DFNMaterial::Eintact<<")};";
    // outfile<<"Recombine Surface {Physical Surface("<<DFNMaterial::Erefined<<")};\n";

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




























// void DFNMesh::CreateVolumes(){
//     TPZGeoMesh *gmesh = this->Mesh();
	
// 	int dim = gmesh->Dimension();
// 	if(dim != 3) return; //@todo temporary while DFNVolume isn't generalized to 2D meshes
//     // map all volumes that are cut
//     int64_t nels = gmesh->NElements();
// 	for (int64_t iel = 0; iel < nels; iel++){
//         TPZGeoEl *gel = gmesh->Element(iel);
//         if(gel->Dimension() != dim){continue;}
//         int nsides = gel->NSides();
//         // int ncorners = gel->NCornerNodes();
// 		// int nfaces = (int) (nsides+1-2*ncorners)/2 //from Euler's characteristic
//         for (int iside = nsides-2; iside > 0; iside--){
//             TPZGeoElSide gelside(gel,iside);
//             if (gelside.Dimension() != 2){break;}
//             TPZGeoElSide neighbour = gelside.Neighbour();
// 			while(neighbour.Element()->Dimension() != 2 && neighbour != gelside){
// 				neighbour = neighbour.Neighbour();
// 			}
// 			if(neighbour == gelside) continue;
// 			TPZGeoEl *sideface = neighbour.Element();
//             // @phil I dont understand the logic of this. Please insert documentation
// 			if(sideface->HasSubElement()){
//                 DFNVolume volume(iel,true);
//                 fVolumes.insert({iel,volume});
//                 break;
//             }            
//         }
//     }
    
//     // gmesh->BuildConnectivity(); //@todo remove this after test
// 	// search through each 2D element of the triangulated fractures surfaces to find their enclosing volume
// 	for(int64_t iel = 0; iel < nels; iel++){
// 		TPZGeoEl *gel = gmesh->Element(iel);
// 		if(!gel) continue;
// 		if(gel->MaterialId() != DFNMaterial::Efracture){continue;}
// 		if(gel->Dimension() != 2) continue;
// 		if(gel->HasSubElement()) continue;
// 		// Find volume that encloses that element
// 		FindEnclosingVolume(gel);
// 	}

// }































// /**
//  * 	@brief Uses GMsh API to tetrahedralize a DFNVolume
//  */ 
// void DFNMesh::Tetrahedralize(DFNVolume *volume){
// 	// GMsh doesn't like zero index entities
//     const int shift = 1;

// 	TPZGeoMesh *gmesh = fGMesh;
// 	int mesh_dim = gmesh->Dimension();
// 	int64_t ivol = volume->ElementIndex();
// 	TPZGeoEl *volGel = gmesh->Element(ivol);
// 	// std::cout<<"\n\n _________   volume # "<<ivol<<"\n";
// 	// List faces that form the volume shell
// 	std::vector<int> surfaceloop;
// 	for(int nsides = volGel->NSides(),
// 					 iside = nsides-2; iside >= 0; iside--){
// 		TPZGeoElSide gelside(volGel,iside);
// 		if(gelside.Dimension() != 2) break;
// 		TPZGeoElSide neig = gelside.Neighbour();
// 		while(neig.Element()->Dimension() != 2){ 
// 			neig = neig.Neighbour();
// 			if(neig == gelside){
// 				PZError << "\n\n Error at "<<__PRETTY_FUNCTION__<<"\n"<<"There shouldn't be a volume without skeleton\n\n";
// 				DebugStop();
// 			}
// 		}
// 		surfaceloop.push_back(neig.Element()->Index());
// 	}
// 	// List faces that are enclosed in the volume
// 	TPZManVector<int64_t> enclosedFaces = volume->GetFacesInVolume();

// 	// List all lines
// 	std::set<int64_t> lines;
// 	// List edges from volume shell
// 	for(int iface = 0,
// 			nfaces = surfaceloop.size(); iface < nfaces; iface++){
// 		TPZGeoEl *gel = gmesh->Element(surfaceloop[iface]);
// 		for(int nsides = gel->NSides(),
// 				nnodes = gel->NCornerNodes(),
// 				iside = nsides-2; iside >= nnodes; iside--){
// 			TPZGeoElSide gelside(gel,iside);
// 			if(gelside.Dimension() != 1) break;
// 			TPZGeoElSide neig = gelside.Neighbour();
// 			while(neig.Element()->Dimension() != 1) neig = neig.Neighbour();
// 			TPZStack<TPZGeoEl*> unrefined_lines;
// 			if(neig.Element()->HasSubElement()){
// 				neig.Element()->YoungestChildren(unrefined_lines);
// 			}else{
// 				unrefined_lines.Push(neig.Element());
// 			}
// 			for(auto line : unrefined_lines){
// 				lines.insert(line->Index());
// 			}
// 		}
// 	}
// 	// List edges from faces in volume
// 	for(int iface = 0,
// 			nfaces = enclosedFaces.size(); iface < nfaces; iface++){
// 		TPZGeoEl *gel = gmesh->Element(enclosedFaces[iface]);
// 		for(int nsides = gel->NSides(),
// 				nnodes = gel->NCornerNodes(),
// 				iside = nsides-2; iside >= nnodes; iside--){
// 			TPZGeoElSide gelside(gel,iside);
// 			if(gelside.Dimension() != 1) break;
// 			TPZGeoElSide neig = gelside.Neighbour();
// 			while(neig.Element()->Dimension() != 1) neig = neig.Neighbour();
// 			TPZStack<TPZGeoEl*> unrefined_lines;
// 			if(neig.Element()->HasSubElement()){
// 				neig.Element()->YoungestChildren(unrefined_lines);
// 			}else{
// 				unrefined_lines.Push(neig.Element());
// 			}
// 			for(auto line : unrefined_lines){
// 				lines.insert(line->Index());
// 			}
// 		}
// 	}

// 	// List nodes
// 	std::set<int64_t> nodes;
// 	for(int64_t line : lines){
// 		TPZGeoEl *gel = gmesh->Element(line);
// 		nodes.insert(gel->NodeIndex(0));
// 		nodes.insert(gel->NodeIndex(1));
// 	}
	
// 	// gmsh::initialize();
// 	std::string modelname = "model"+std::to_string(ivol);
// 	gmsh::model::add(modelname);
// 	gmsh::model::setCurrent(modelname);
// 	std::string mshfilename("LOG/testAPI_volume");
// 	mshfilename += std::to_string(ivol);
// 	// gmsh::model::add(mshfilename);
// 	mshfilename += ".msh";
// 	gmsh::option::setNumber("Mesh.Algorithm3D",1);  // (1: Delaunay, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT) Default value: 1
// 	// Insert nodes ____________________________________
// 	for(int64_t inode : nodes){
// 		TPZManVector<REAL,3> coord(3);
// 		gmesh->NodeVec()[inode].GetCoordinates(coord);
// 		gmsh::model::geo::addPoint(coord[0],coord[1],coord[2],0.,inode+shift);
// 	}
// 	// std::cout<<"\n\n\n";
// 	// Insert lines ____________________________________
// 	for(int64_t iline : lines){
// 		TPZGeoEl *gel = gmesh->Element(iline);
// 		int64_t node0 = gel->NodeIndex(0)+shift;
// 		int64_t node1 = gel->NodeIndex(1)+shift;
// 		gmsh::model::geo::addLine(node0,node1,iline+shift);
// 		gmsh::model::geo::mesh::setTransfiniteCurve(iline+shift,2);
// 	}
// 	// Insert faces ____________________________________
// 	{
// 	std::vector<int> wiretag(1);
// 	// Faces in volume shell
// 	for(int64_t faceindex : surfaceloop){
// 		TPZGeoEl *face = gmesh->Element(faceindex);
// 		int nnodes = face->NCornerNodes();
// 		int nedges = nnodes;
// 		TPZManVector<int64_t,4> facenodevec(nnodes);
// 		face->GetNodeIndices(facenodevec);
// 		// line loop ________________________________________________________________________
// 		std::vector<int> lineloop;
// 		int64_t reference_node = face->NodeIndex(0);
// 		bool planesurface_Q = false;
// 		for(int iside = nnodes; iside < nnodes+nedges; iside++){
// 			TPZGeoElSide gelside(face,iside);
// 			TPZGeoElSide side = gelside.Neighbour();
// 			// find line element
// 			while(side.Element()->Dimension() != 1) {side = side.Neighbour();}
// 			TPZGeoEl *irib = side.Element();
// 			TPZStack<TPZGeoEl*> children;
// 			if(irib->HasSubElement()){
// 				irib->YoungestChildren(children);
// 			}else{
// 				children.Push(irib);
// 			}
// 			int nchildren = children.NElements();
// 			if(nchildren > 1) planesurface_Q = 1;
// 			int children_added = 0;
// 			int64_t index = 0;
// 			TPZManVector<bool> added_list(nchildren,false);
// 			while(children_added < nchildren){
// 				// find next child and check orientation
// 				int orientation = 1;
// 				for(int ichild = 0; ichild<nchildren; ichild++){
// 					if(added_list[ichild]) continue;
// 					TPZGeoEl *child = children[ichild];
// 					if(child->NodeIndex(0) == reference_node){
// 						orientation = 1;
// 						reference_node = child->NodeIndex(1); // next node becomes reference
// 						index = child->Index();
// 						added_list[ichild] = true;
// 						break;
// 					}
// 					else if(child->NodeIndex(1) == reference_node){
// 						orientation = -1;
// 						reference_node = child->NodeIndex(0); // next node becomes reference
// 						index = child->Index();
// 						added_list[ichild] = true;
// 						break;
// 					}
// 				}
// 				// add child
// 				children_added++;
// 				index = orientation*(index+shift);
// 				lineloop.push_back(index);
// 			}
// 		}
// 		// insert curve loop ________________________________________________________________
// 		wiretag[0] = gmsh::model::geo::addCurveLoop(lineloop,face->Index()+shift);
// 		// insert surface
// 		if(planesurface_Q){gmsh::model::geo::addPlaneSurface(wiretag,wiretag[0]);
// 		}else{gmsh::model::geo::addSurfaceFilling(wiretag,wiretag[0]);}
		
// 		// @todo gmsh::model::geo::mesh::setTransfiniteSurface(wiretag[0]);
// 		// if(face->Type() == EQuadrilateral) gmsh::model::geo::mesh::setRecombine(2,wiretag[0]);
// 	}
// 	// Enclosed faces
// 	for(int64_t faceindex : enclosedFaces){
// 		TPZGeoEl *face = gmesh->Element(faceindex);
// 		int nnodes = face->NCornerNodes();
// 		int nedges = nnodes;
// 		TPZManVector<int64_t,4> facenodevec(nnodes);
// 		face->GetNodeIndices(facenodevec);
// 		// line loop
// 		std::vector<int> lineloop(nedges);
// 		for(int iside = nnodes; iside < nnodes+nedges; iside++){
// 			TPZGeoElSide gelside(face,iside);
// 			TPZGeoElSide neig = gelside.Neighbour();
// 			// find line element
// 			while(neig.Element()->Dimension()!=1){neig = neig.Neighbour();}
// 			// find first node of line at the face
// 			int inode = 0;
// 			while(facenodevec[inode] != neig.SideNodeIndex(0)) inode++;
// 			// check orientation by comparing second node of line with next node of face
// 			if(neig.SideNodeIndex(1) == facenodevec[(inode+1)%nnodes]){
// 				lineloop[iside-nnodes] = neig.Element()->Index() + shift;
// 			}else{
// 				lineloop[iside-nnodes] = - (neig.Element()->Index() + shift);
// 			}
// 		}
// 		// insert curve loop
// 		wiretag[0] = gmsh::model::geo::addCurveLoop(lineloop,face->Index()+shift);
// 		// insert surface
// 		gmsh::model::geo::addSurfaceFilling(wiretag,wiretag[0]);	
// 		gmsh::model::geo::mesh::setTransfiniteSurface(wiretag[0]);
// 		// if(face->Type() == EQuadrilateral) gmsh::model::geo::mesh::setRecombine(2,wiretag[0]);
// 	}
// 	}
	
// 	{// begin Lines-in-Surface ________________________________________________
// 	// "lines of fracture surface whose elder does not have a mesh_dim neighbour without father"
// 	for(int iface : surfaceloop){
// 		TPZGeoEl *gel = gmesh->Element(iface);
// 		if(!gel->HasSubElement()) continue;
// 		TPZStack<TPZGeoEl *> children;
// 		gel->YoungestChildren(children);
// 		// queue all possible lines by checking 1D neighbours of children
// 		std::set<TPZGeoEl *> candidate_ribs;
// 		for(auto child : children){
// 			int nribs = child->NCornerNodes();
// 			for(int cside = nribs; cside < 2*nribs; cside++){
// 				TPZGeoElSide childside(child,cside);
// 				TPZGeoElSide neig = childside.Neighbour();
// 				for(/*void*/; neig != childside; neig = neig.Neighbour()){
// 					if(neig.Element()->Dimension() != 1) continue;
// 					if(neig.Element()->MaterialId() != DFNMaterial::Efracture) continue;
// 					candidate_ribs.insert(neig.Element());
// 				}
// 			}
// 		}
// 		bool should_enter = true;
// 		std::vector<int> lines_in_surface;
// 		auto end = candidate_ribs.end();
// 		for(auto it = candidate_ribs.begin(); it != end;it++){
// 			TPZGeoEl *rib = *it;
// 			TPZGeoEl *elder = (rib->Father()? rib->EldestAncestor() : rib);
// 			TPZGeoElSide elderside(elder,2);
// 			TPZGeoElSide neig = elderside.Neighbour();
// 			while(neig != elderside){
// 				if(neig.Element()->Dimension() == mesh_dim && !neig.Element()->Father()){
// 					should_enter = false;
// 					break;
// 				}
// 				neig = neig.Neighbour();
// 			}
// 			if(!should_enter){continue;}
// 			lines_in_surface.push_back(rib->Index()+shift);
// 		}
// 		gmsh::model::geo::synchronize(); // synchronize is required before embedding
// 		gmsh::model::mesh::embed(1,lines_in_surface,2,iface+shift);
// 	}
// 	}// end Lines-in-Surface __________________________________________________


// 	// Insert volumes ____________________________________
// 	std::vector<int> shelltag(1);
// 	// Shift surfaceloop indices
// 	for(int nsurfaces = surfaceloop.size(),
// 			i = 0; i < nsurfaces; i++){
// 		surfaceloop[i] += shift;
// 	}
// 	shelltag[0] = gmsh::model::geo::addSurfaceLoop(surfaceloop,ivol+shift);
// 	gmsh::model::geo::addVolume(shelltag,shelltag[0]);
	
// 	// Surfaces in Volume ________________________________
// 		gmsh::model::geo::synchronize(); // synchronize is required before embedding
// 		int nfacesenclosed = enclosedFaces.size();
// 		std::vector<int> facesinvolume(nfacesenclosed);
// 		{
// 			for(int i = 0; i<nfacesenclosed; i++){
// 				facesinvolume[i] = enclosedFaces[i] + shift;
// 			}
// 			gmsh::model::mesh::embed(2, facesinvolume, 3, ivol+shift);
// 		}
// 	// Physical groups ____________________________
// 		gmsh::model::addPhysicalGroup(2,facesinvolume,DFNMaterial::Efracture);
// 		gmsh::model::addPhysicalGroup(2,surfaceloop,DFNMaterial::Erefined);
// 		gmsh::model::addPhysicalGroup(3,shelltag,DFNMaterial::Erefined);
	
// 	// synchronize before meshing
// 		gmsh::model::geo::synchronize();
// 	// mesh
// 		gmsh::model::mesh::generate(3);
// 	gmsh::write(mshfilename);
// 	// gmsh::write("LOG/testAPI_volume.msh");
// 	// import meshed volume back into PZ geoMesh
// 		ImportElementsFromGMSH(gmesh,3,nodes);
// 	gmsh::model::remove();
// 	gmsh::clear();
// 	// gmsh::finalize();
// }















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








void DFNMesh::ClearMaterials(){
	for(auto el : fGMesh->ElementVec()){
		if(!el) continue;
		el->SetMaterialId(DFNMaterial::Eintact);
	}
}

void DFNMesh::RestoreMaterials(TPZGeoMesh *originalmesh){
	// @todo
	DebugStop();
}

// template<int Talloc>
// bool DFNMesh::IsConvexPolyhedron(TPZStack<std::pair<int64_t,int>,Talloc>& polyhedron){
// 	for(auto iface_orient : polyhedron){
// 		if(iface_orient.first < 0) DebugStop();
// 		if(iface_orient.second == 0) DebugStop();
// 		// Edges occupying the one dimensional sides
// 		TPZManVector<int64_t,4> edges = GetEdgeIndices(iface_orient.first);
// 		// Neighbour cards through the edges
// 		TPZManVector<TRolodexCard,4> cards(edges.size());
// 		// List the neighbour cards
// 		for(int i=0; i<edges.size(); i++){
// 			int64_t iedge = edges[i];
// 			TRolodex& rolodex = fSortedFaces[iedge];
// 			cards[i] = rolodex.Card(iface_orient.first);
// 		}

// 		// Determine orientation of neighbour cards
// 		TPZManVector<std::pair<TRolodexCard, int>,4> facingcards(edges.size());
// 		for(int i=0; i<edges.size(); i++){
// 			REAL angle_debug = 0.0;
// 			int64_t iedge = edges[i];
// 			TRolodex& rolodex = fSortedFaces[iedge];
// 			std::pair<TRolodexCard, int> current_card = {cards[i],iface_orient.second};
// 			facingcards[i] = rolodex.FacingCard(current_card,angle_debug);


// 			REAL angle = iface_orient.second*current_card.first.forientation*(facingcards[i].first.fangle_to_reference - current_card.first.fangle_to_reference);
// 			if(angle < 0.) {angle = angle +DFN::_2PI;}
// 			REAL test = std::abs(angle - angle_debug);
// 			if( test > 1e-5) DebugStop();
// 			if(angle > M_PI+gDFN_SmallNumber){ return false; }
			
// 		}
// 	}
// 	return true;
// }

void DFNMesh::InitializePolyhedra(){
	// Clear data structure
	fPolyhedra.clear();
	fPolyh_per_face.clear();
	fPolyh_per_face.Resize(fGMesh->NElements(),{-1,-1});
	int ipolyh=0; 											///< index of polyhedron
	TPZStack<std::pair<int64_t,int>,20> shell(20,{-1,0});	///< oriented faces that enclose the polyhedron


	// Gather mesh boundary first
	shell.clear();
	DFNPolyhedron polyhedron;
	TPZGeoElSide gelside;
	TPZGeoElSide neig;
	for(TPZGeoEl* gel : fGMesh->ElementVec()){
		if(gel->Dimension() != 2) continue;
		// if(gel->MaterialId() != DFNMaterial::Eskeleton) continue;

		gelside = {gel,gel->NSides()-1};
		int nneighbours = 0;
		TPZGeoElSide volneig;
		for(neig=gelside.Neighbour(); neig != gelside; neig = neig.Neighbour()){
			if(neig.Element()->Dimension() < 3) continue;
			nneighbours += neig.Element()->Dimension() > 2;
			volneig = neig; //< catch a volume in case there's a boundary condition 
		}
		if(nneighbours > 1) continue;
		// Faces don't have guaranteed positive orientation, since they were created using CreateBCGeoEl(). @see  DFN::CreateSkeletonElements()
		int orient = DFN::PZOutwardPointingFace(volneig) ? 1 : -1; 
		shell.push_back({gel->Index(),orient});
		SetPolyhedralIndex({gel->Index(),orient},0);
	}
	polyhedron.Initialize(this,ipolyh,shell);
	// polyhedron.Print();
	fPolyhedra.push_back(polyhedron);
	shell.Fill({-1,0});
	shell.clear();

	// Then initialize the rest of the polyhedra
	// BuildPolyhedra();
	UpdatePolyhedra();
	// @maybeTODO - Looping over 3D elements and matching a polyhedral index to their shell would be more efficient. But I have to pick my battles
}

void DFNMesh::UpdatePolyhedra(){
	// Start by sorting faces around edges and filling the this->fSortedFaces datastructure
	this->SortFacesAroundEdges();
	std::cout<<"Updating polyhedral volumes\r"<<std::flush;
	fPolyh_per_face.Resize(fGMesh->NElements(),{-1,-1});
	// Refined faces pass down their polyh index to their subelements
	InheritPolyhedra();
	int polyh_index = fPolyhedra.size();

	
	TPZStack<std::pair<int64_t,int>,20> polyhedron(20,{-1,0});
	// loop over 2D skeleton elements
	for(TPZGeoEl* initial_face : fGMesh->ElementVec()){
		if(!initial_face) continue;
		if(initial_face->Dimension() != 2) continue;
		if(initial_face->HasSubElement()) continue;
		// if(initial_face->MaterialId() != DFNMaterial::Eskeleton && 
		// 	initial_face->MaterialId() != DFNMaterial::Efracture) continue;
		int64_t initial_id = initial_face->Index();
		
		// look for polyhedron on each orientation of the initial_face
		for(int ipolyh_local=0; ipolyh_local<2; ipolyh_local++){
			// if the first has been found, continue to the next
			if(fPolyh_per_face[initial_id][ipolyh_local] != -1) continue;
			polyhedron.Fill({-1,0});
			polyhedron.clear();
			int orientation = ipolyh_local?-1:1;
			std::pair<int64_t,int> initial_face_orient = {initial_id,orientation};
			polyhedron.push_back({initial_id,orientation});
			polyh_index = fPolyhedra.size();
			SetPolyhedralIndex(initial_face_orient,polyh_index);
			bool IsConvex = true;
			BuildVolume(initial_face_orient,IsConvex,polyhedron);

			#ifdef PZDEBUG
				// std::cout << "Polyh#"<< std::setw(4) << fPolyhedra.size();;
				// std::cout << ":  " << polyhedron << std::endl;
			#endif //PZDEBUG

			if(!IsConvex) {
				MeshPolyhedron(polyhedron);
				ClearPolyhIndex(polyhedron);
				this->SortFacesAroundEdges();
				--ipolyh_local;
			}else{
				DFNPolyhedron new_polyhedron(this,polyh_index,polyhedron);
				fPolyhedra.push_back(new_polyhedron);
			}
		}
	}
	std::cout<<"                             \r"<<std::flush;
}

void DFNMesh::ClearPolyhIndex(TPZVec<std::pair<int64_t,int>>& facestack){
	for(auto& faceorient : facestack){
		SetPolyhedralIndex(faceorient, -1);
	}
}

// void DFNMesh::BuildPolyhedra(){
// 	// Start by sorting faces around edges and filling the this->fSortedFaces datastructure
// 	this->SortFacesAroundEdges();
// 	fPolyh_per_face.Resize(fGMesh->NElements(),{-1,-1});
// 	int polyh_counter = fPolyhedra.size();
	
// 	TPZStack<std::pair<int64_t,int>,20> polyhedron(20,{-1,0});
// 	// loop over 2D skeleton elements
// 	for(TPZGeoEl* initial_face : fGMesh->ElementVec()){
// 		if(!initial_face) continue;
// 		if(initial_face->Dimension() != 2) continue;
// 		if(initial_face->HasSubElement()) continue;
// 		int64_t initial_id = initial_face->Index();

// 		// look for polyhedron on each orientation of the initial_face
// 		for(int ipolyh_local=0; ipolyh_local<2; ipolyh_local++){
// 			// if the first has been found, continue to the next
// 			if(fPolyh_per_face[initial_id][ipolyh_local] != -1) continue;
// 			polyhedron.Fill({-1,0});
// 			polyhedron.clear();
// 			int orientation = ipolyh_local?-1:1;
// 			std::pair<int64_t,int> initial_face_orientation = {initial_id,orientation};
// 			polyhedron.push_back({initial_id,orientation});
// 			SetPolyhedralIndex(initial_face_orientation,polyh_counter);
// 			bool IsConvex = true;
// 			BuildVolume(initial_face_orientation,IsConvex,polyhedron);

// 			#ifdef PZDEBUG
// 				// std::cout << "Polyh#"<< std::setw(4) << polyh_counter;
// 				// std::cout << ":  " << polyhedron << std::endl;
// 			#endif //PZDEBUG

// 			if(!IsConvex) MeshPolyhedron(polyhedron);
// 			else{
// 				DFNPolyhedron new_polyhedron(this,polyh_counter,polyhedron);
// 				fPolyhedra.push_back(new_polyhedron);
// 			}
// 			polyh_counter++;
// 		}
// 	}
// }





void DFNMesh::SetPolyhedralIndex(std::pair<int64_t,int> face_orient, int polyh_index){
	switch(face_orient.second){
		case  1: fPolyh_per_face[face_orient.first][0] = polyh_index; break;
		case -1: fPolyh_per_face[face_orient.first][1] = polyh_index; break;
		default: DebugStop();
	}
}
int DFNMesh::GetPolyhedralIndex(std::pair<int64_t,int> face_orient){
	int polyh_index = -1;
	switch(face_orient.second){
		case  1: polyh_index = fPolyh_per_face[face_orient.first][0]; break;
		case -1: polyh_index = fPolyh_per_face[face_orient.first][1]; break;
		default: DebugStop();
	}
	return polyh_index;
}

template<int Talloc>
void DFNMesh::BuildVolume(std::pair<int64_t,int> initial_face_orient, bool& IsConvex, TPZStack<std::pair<int64_t,int>, Talloc>& polyhedron){
	int polyh_index = GetPolyhedralIndex(initial_face_orient);
	if(polyh_index == -1) DebugStop();
	// Edges occupying the one dimensional sides
	TPZManVector<int64_t,4> edges = GetEdgeIndices(initial_face_orient.first);
	// Neighbour cards through the edges
	TPZManVector<TRolodexCard,4> cards(edges.size());
	// List the neighbour cards
	for(int i=0; i<edges.size(); i++){
		int64_t iedge = edges[i];
		TRolodex& rolodex = fSortedFaces[iedge];
		cards[i] = rolodex.Card(initial_face_orient.first);
	}

	// Determine orientation of neighbour cards
	TPZManVector<std::pair<TRolodexCard, int>,4> facingcards(edges.size());
	for(int i=0; i<edges.size(); i++){
		int64_t iedge = edges[i];
		TRolodex& rolodex = fSortedFaces[iedge];
		std::pair<TRolodexCard, int> current_card = {cards[i],initial_face_orient.second};
		REAL angle = 0.0;
		facingcards[i] = rolodex.FacingCard(current_card,angle);
		if(angle > M_PI+gDFN_SmallNumber){ IsConvex = false; }
	}
	// queue neighbour cards to verify
	std::list<std::pair<int64_t, int>> to_verify;
	for(int i=0; i<edges.size(); i++){
		int64_t iedge = edges[i];
		std::pair<int64_t,int> faceorient = {facingcards[i].first.fgelindex,facingcards[i].second};
		int nextface_polyindex = GetPolyhedralIndex(faceorient);
		if(nextface_polyindex != polyh_index){
			to_verify.push_back(faceorient);
			SetPolyhedralIndex(faceorient,polyh_index);
			polyhedron.push_back(faceorient);
		}
	}
	// recursively try to add neighbours of cards that have been queued 
	for(auto orientedface : to_verify){
		BuildVolume(orientedface,IsConvex,polyhedron);
	}

}

// template<int Talloc>
// void DFNMesh::BuildVolume(std::pair<int64_t,int> initial_face_orient, bool& IsConvex, TPZStack<std::pair<int64_t,int>, Talloc>& polyhedron){
// 	int polyh_index = GetPolyhedralIndex(initial_face_orient);
// 	if(polyh_index == -1) DebugStop();
// 	// Edges occupying the one dimensional sides
// 	TPZManVector<int64_t,4> edges = GetEdgeIndices(initial_face_orient.first);
// 	// Neighbour cards through the edges
// 	TPZManVector<TRolodexCard,4> cards(edges.size());
// 	// List the neighbour cards
// 	for(int i=0; i<edges.size(); i++){
// 		int64_t iedge = edges[i];
// 		TRolodex& rolodex = fSortedFaces[iedge];
// 		cards[i] = rolodex.Card(initial_face_orient.first);
// 	}

// 	// Determine orientation of neighbour cards
// 	TPZManVector<std::pair<TRolodexCard, int>,4> facingcards(edges.size());
// 	for(int i=0; i<edges.size(); i++){
// 		int64_t iedge = edges[i];
// 		TRolodex& rolodex = fSortedFaces[iedge];
// 		std::pair<TRolodexCard, int> current_card = {cards[i],initial_face_orient.second};
// 		REAL angle = 0.0;
// 		facingcards[i] = rolodex.FacingCard(current_card,angle);
// 		if(angle > M_PI+gDFN_SmallNumber){ IsConvex = false; }
// 	}
// 	// queue neighbour cards to verify
// 	std::list<std::pair<int64_t, int>> to_verify;
// 	for(int i=0; i<edges.size(); i++){
// 		int64_t iedge = edges[i];
// 		std::pair<int64_t,int> faceorient = {facingcards[i].first.fgelindex,facingcards[i].second};
// 		int nextface_polyindex = GetPolyhedralIndex(faceorient);
// 		if(nextface_polyindex == -1) {
// 			to_verify.push_back(faceorient);
// 			SetPolyhedralIndex(faceorient,polyh_index);
// 			polyhedron.push_back(faceorient);
// 		}else if(nextface_polyindex != polyh_index){
// 			DebugStop();
// 		}
// 	}
// 	// recursively verify cards that have been queued 
// 	for(auto orientedface : to_verify){
// 		BuildVolume(orientedface,IsConvex,polyhedron);
// 	}
	
// }




template<int Talloc>
void DFNMesh::MeshPolyhedron(TPZStack<std::pair<int64_t,int>,Talloc>& polyhedron){
	// GMsh doesn't like zero index entities
    const int shift = 1;

	// copy faces of polyhedron into an std::vector for gmsh
	std::vector<int> surfaceloop;
	int nfaces = polyhedron.size();
	surfaceloop.resize(nfaces);
	for(int i=0; i<nfaces; i++){
		surfaceloop[i] = polyhedron[i].first;
	}

	// List all lines
	std::set<int64_t> lines;
	// List edges from volume shell
	for(int iface = 0; iface < nfaces; iface++){
		TPZManVector<int64_t,4> face_lines = GetEdgeIndices(surfaceloop[iface]);
		for(auto iline : face_lines){
			lines.insert(iline);
		}
	}

	// List nodes
	std::set<int64_t> nodes;
	for(int64_t line : lines){
		TPZGeoEl *gel = fGMesh->Element(line);
		nodes.insert(gel->NodeIndex(0));
		nodes.insert(gel->NodeIndex(1));
	}

	// gmsh::initialize();
	std::string modelname = "model_polyh";
	gmsh::model::add(modelname);
	gmsh::model::setCurrent(modelname);
	std::string mshfilename = "LOG/gmshAPI_polyh.msh";
	gmsh::option::setNumber("Mesh.Algorithm3D",1);  // (1: Delaunay, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT) Default value: 1
	// Insert nodes ____________________________________
	{TPZManVector<REAL,3> coord(3);
	for(int64_t inode : nodes){
		fGMesh->NodeVec()[inode].GetCoordinates(coord);
		REAL meshsize = 0.;
		gmsh::model::geo::addPoint(coord[0],coord[1],coord[2],meshsize,inode+shift);
	}}
	// Insert lines ____________________________________
	for(int64_t iline : lines){
		TPZGeoEl *gel = fGMesh->Element(iline);
		int64_t node0 = gel->NodeIndex(0)+shift;
		int64_t node1 = gel->NodeIndex(1)+shift;
		gmsh::model::geo::addLine(node0,node1,iline+shift);
		gmsh::model::geo::mesh::setTransfiniteCurve(iline+shift,2);
	}
	// Insert faces ____________________________________
	{
		// wiretag is a dummy vector with the shifted index of the face/curve-loop
		std::vector<int> wiretag(1,-1);
		std::vector<int> lineloop;
		lineloop.reserve(4);
		for(int64_t iface : surfaceloop){
			TPZGeoEl* face_el = fGMesh->Element(iface);
			DFN::GetLineLoop(face_el,lineloop,shift);
			wiretag[0] = iface+shift;
			gmsh::model::geo::addCurveLoop(lineloop,wiretag[0]);
			gmsh::model::geo::addSurfaceFilling(wiretag,wiretag[0]);
			gmsh::model::geo::mesh::setTransfiniteSurface(wiretag[0]);
			// gmsh::model::geo::mesh::setRecombine(2,wiretag[0]);
		}
	}

	// Set volume _______________________________________
	{
		std::vector<int> shelltag(1,1);
		// shift indices in surface loop
		for(int i=0; i < nfaces; i++){
			surfaceloop[i] += shift;
		}
		gmsh::model::geo::addSurfaceLoop(surfaceloop,shelltag[0]);
		gmsh::model::geo::addVolume(shelltag,shelltag[0]);
		gmsh::model::addPhysicalGroup(3,shelltag,DFNMaterial::Eintact);
	}

	// synchronize before meshing
	gmsh::model::geo::synchronize();
	// mesh
	gmsh::model::mesh::generate(3);
	#ifdef PZDEBUG
		gmsh::write(mshfilename);
	#endif //PZDEBUG
	// import meshed volume back into PZ geoMesh
	std::set<int64_t>& old_nodes = nodes;
	ImportElementsFromGMSH(fGMesh,3,old_nodes);
	gmsh::model::remove();
	gmsh::clear();
	// gmsh::finalize();
	CreateSkeletonElements(1);
	CreateSkeletonElements(2);
}






// /// Appends all rolodex-neighbours of a face to a polyhedron (if they've not been appended already)
// /// If any dihedral angle is found higher than pi radians, change flag to non-convex polyhedron
// void DFNMesh::AppendNeighboursToPolyhedron(TPZGeoEl* current_face, std::vector<int>& polyhedron, const int polyh_index, bool& convexPolyh){
// 	// consistency checks
// 	if(!current_face) 							DebugStop();
// 	if(current_face->Dimension() != 2) 			DebugStop();
// 	if(current_face->HasSubElement()) 			DebugStop();
	
// 	int64_t current_index = current_face->Index();
// 	int currentface_orientation = (fPolyh_per_face[current_index][0]==polyh_index ? 1:-1);
// 	int nedges = current_face->NSides(1);
// 	bool newface=false;
// 	// Loop over each edge of the current 2D_el to append neighbours to the polyhedron. 
// 	for(int iside=nedges; iside<current_face->NSides()-1; iside++){
// 		TPZGeoEl* edge = DFN::GetSkeletonNeighbour(current_face,iside);
// 		TPZGeoElSide edgeside(edge,2);
// 		TPZGeoElSide gelside(current_face,iside);
// 		int direction = currentface_orientation * (DFN::OrientationMatch(edgeside,gelside)? 1:-1);
// 		REAL angle = 0.0;
// 		TRolodex& rolodex = fSortedFaces[edge->Index()];
// 		TRolodexCard& nextcard = rolodex.NextCard(current_index,angle,direction);
// 		int nextcard_orientation = -direction*nextcard.forientation;
		
// 		//@todo Phil did it better
// 		newface = std::find(polyhedron.begin(), polyhedron.end(), nextcard.fgelindex) == polyhedron.end();
// 		if(newface) polyhedron.push_back(nextcard.fgelindex);
// 		else continue;
		
// 		if(nextcard_orientation > 0)
// 			{fPolyh_per_face[nextcard.fgelindex][0] = polyh_index;}
// 		else
// 			{fPolyh_per_face[nextcard.fgelindex][1] = polyh_index;}
// 		// check if polyh is not convex
// 		if(angle > M_PI+gDFN_SmallNumber){convexPolyh = false;}
// 	}//for each edge

// }


// void DFNMesh::GetPolyhedra2(){
// 	// Start by sorting faces around edges and filling the this->fSortedFaces datastructure
// 	this->SortFacesAroundEdges();
// 	fPolyh_per_face.Resize(fGMesh->NElements(),{-1,-1});
// 	TPZGeoEl* initial_face = nullptr;
//     // For each polyhedra, a vector of face indexes 
// 	// @todo this will probably only be used for debug. We can decide after this method is finished
// 	TPZStack<TPZAutoPointer<std::vector<int>>> polyhedra_vec;
// 	int polyh_counter = 0;

// 	std::cout<<"\n\n Print polyhedra:\n\n";
// 	// loop over 2D skeleton elements
// 	for(TPZGeoEl* el : fGMesh->ElementVec()){
// 		if(!el) continue;
// 		if(el->Dimension() != 2) continue;
// 		if(el->HasSubElement()) continue;
// 		//@todo filter skeleton + fracture material IDs

// 		initial_face = el;
// 		// number of polyhedra that have to be found for this face
// 		const int npolyh_total = 2;

// 		// flank defines if we're looking for a polyhedron on the positive or the negative side of the face
// 		int initial_face_flank = (fPolyh_per_face[el->Index()][0] >= 0 ? -1 : 1);
// 		// number of polyhedra found on either flank of this face
// 		int npolyh_found = int(fPolyh_per_face[el->Index()][0] >= 0) + int(fPolyh_per_face[el->Index()][1] >= 0);

// 		// for each 2D skeleton we get 2 polyhedra
// 		while(npolyh_found < npolyh_total){
// 			// a vector of faces that enclose a polyhedron
// 			TPZAutoPointer<std::vector<int>> polyhedron = new std::vector<int>;
// 			// a set of faces whose neighbours have already been queued and added to the polyhedron (to avoid double-checking)
// 			std::set<int64_t> verified;
// 			// match face and polyhedron in datastructure
// 			polyhedron->push_back(initial_face->Index());
// 			if(initial_face_flank > 0)
// 				{fPolyh_per_face[initial_face->Index()][0] = polyh_counter;}
// 			else
// 				{fPolyh_per_face[initial_face->Index()][1] = polyh_counter;}
			
// 			// flags
// 			bool convexPolyh = true;
// 			// std::cout<<"\n\nPolyhedron #"<<polyh_counter;
// 			// for every face found to belong to the polyhedron, try to append its neighbours
// 			for(int i=0; i < polyhedron->size(); i++){
// 				int64_t current_index = polyhedron->operator[](i);
// 				// try to insert in verified, if it's not inserted than neighbours of this element have already been queued and we can move down the queue 
// 				bool newface = verified.insert(current_index).second;
// 				if(!newface) continue;
// 				TPZGeoEl* current_face = fGMesh->Element(current_index);
// 				AppendNeighboursToPolyhedron(current_face,*polyhedron,polyh_counter,convexPolyh);
// 			}//for each face in the polyhedron
// 			polyhedra_vec.push_back(polyhedron);
// 			polyh_counter++;
// 			npolyh_found++;
// 			// revert initial orientation to look for polyhedron on the other side of initial_face
// 			initial_face_flank *= -1;
// 			// if(!convexPolyh) MeshPolyhedron(polyhedron); //@todo
// 		}//while(npolyh_found < npolyh_total)
// 	}//for(TPZGeoEl* el : fGMesh->ElementVec())


// 	// print result for debug________________________________________
// 	int counter = 0;
// 	for(auto polyh : polyhedra_vec){
// 		std::cout<<"\n\nPolyhedron #"<<counter;
// 		int nfaces = (*polyh).size();
// 		for(int i=0; i<nfaces; i++){
// 			std::cout<<std::endl<<(*polyh)[i];
// 		}
// 		counter++;
// 	}
// 	// print result for debug________________________________________
// }




void DFNMesh::SortFacesAroundEdges(){
	// std::cout<<"\n\n Print rolodexes:\n\n";
	// fSortedFaces.clear();
	fSortedFaces.Resize(fGMesh->ElementVec().NElements());
	std::cout<<"Sorting Rolodexes\r"<<std::flush;

	// Loop over 1D elements
	for(TPZGeoEl* gel : fGMesh->ElementVec()){
		if(!gel) continue;
		if(gel->Dimension() != 1) continue;
		if(gel->HasSubElement()) continue;


		std::map<REAL,TPZGeoElSide> facemap;
		TRolodex& rolodex = fSortedFaces[gel->Index()];

		TPZGeoElSide edgeside(gel,2);
		// @todo std::set<int> mats = {DFNMaterial::Efracture, DFNMaterial::Eskeleton}; int nfaces = edgeside.NNeighbours(2,mats);
		int nfaces = edgeside.NNeighbours(2);
		int ncards = rolodex.NCards();
		// skip if there aren't any new faces on the rolodex
		if(nfaces <= ncards) continue;

		TPZGeoElSide neig = edgeside.Neighbour();
		TPZGeoElSide reference_el;
		int reference_orientation;
		for(/*void*/; neig != edgeside; neig = neig.Neighbour()){
			if(neig.Element()->Dimension() != 2) continue;
			if(neig.Element()->HasSubElement()) continue;
			if(facemap.size() == 0) {
				reference_el = neig;
				reference_orientation = (DFN::OrientationMatch(reference_el,edgeside)?1:-1);
			}
			REAL angle = DFN::DihedralAngle(reference_el,neig,reference_orientation);
			facemap.insert({angle,neig});
		}
		int j = 0;
		// std::cout<<"\n Edge # "<<gel->Index();
		rolodex.fedgeindex = gel->Index();
		rolodex.fcards.clear();
		rolodex.fcards.resize(facemap.size());
		for(auto iterator : facemap){
			TPZGeoElSide& faceside = iterator.second;
			TRolodexCard& facecard = rolodex.fcards[j];
			facecard.fgelindex = faceside.Element()->Index();
			facecard.fSide = faceside.Side();
			facecard.fangle_to_reference = iterator.first;
			facecard.forientation = (DFN::OrientationMatch(edgeside,faceside)?1:-1);
			facecard.fposition = j;
			j++;
			// std::cout<<"\n"<<faceside.Element()->Index()<<"\t"<< iterator.first;
		}
		// std::cout<<"\n";
	}
	std::cout<<"                 \r"<<std::flush;
}

// void DFNMesh::IsolateBoundaryPolyh2(){
// 	if(fPolyhedra.size() == 0) DebugStop();
// 	int aux_id = -1;
// 	DFNPolyhedron aux_polyh;
// 	for(auto& polyh : fPolyhedra){
// 		if(polyh.NElements() > pztopology::TPZCube::NumSides(2)){
// 			if(polyh.Index() == 0) return;
// 			aux_id = polyh.Index();

// 			aux_polyh = *fPolyhedra.begin();
// 			*fPolyhedra.begin() = polyh;
// 			polyh = aux_polyh;
			
// 			fPolyhedra.begin()->SwapIndex(0);
// 			polyh.SwapIndex(aux_id);
// 			return;
// 		}
// 	}
// 	PZError << "\nCouldn't find the boundary of the domain?\n";
// 	DebugStop();
// }

void DFNMesh::InheritPolyhedra(){
	fPolyh_per_face.Resize(fGMesh->NElements(),{-1,-1});

	for(TPZGeoEl* father : fGMesh->ElementVec()){
		if(!father) continue;
		if(!father->HasSubElement()) continue;

		TPZGeoEl* child = nullptr;
		int nchildren = father->NSubElements();
		for(int i=0; i<2; i++){
			int orient = i?-1:1;
			int father_polyhindex = GetPolyhedralIndex({father->Index(),orient});
			if(father_polyhindex<0) continue;
			for(int jchild=0; jchild<nchildren; jchild++){
				child = father->SubElement(jchild);
				if(GetPolyhedralIndex({child->Index(),orient}) < 0){
					SetPolyhedralIndex({child->Index(),orient},father_polyhindex);
				}
			}
			SetPolyhedralIndex({father->Index(),orient},-1); //< refined face doesn't need a polyhedral index
			DFNPolyhedron& polyhedron = fPolyhedra[father_polyhindex];
			polyhedron.SwapForChildren(father);
		}
	}
}

void DFNMesh::Print(std::ostream & out, char* name) const
{
	out << "\n\t DFN MESH INFORMATION:\n";
	if(name) out <<"\"" << name <<"\"\n";
	out << "\n\n\t\t GEOMETRIC TPZGeoMesh INFORMATIONS:\n\n";
	out << "TITLE-> " << fGMesh->Name() << "\n\n";
	out << "Number of nodes               = " << fGMesh->NodeVec().NElements() << "\n";
	out << "Number of elements            = " << fGMesh->ElementVec().NElements()-fGMesh->ElementVec().NFreeElements() << "\n";
	out << "Number of free elements       = " << fGMesh->ElementVec().NFreeElements() << "\n\n";

	out << "Number of fractures           = " << fFractures.size() << "\n";
	out << "Number of DFNVolumes          = " << fVolumes.size() << "\n";
	out << "Tolerable distance            = " << fTolDist << "\n";
	out << "Tolerable angle            	  = " << fTolAngle << "\n";
	out << "Tolerable cos(angle)          = " << fTolAngle_cos << "\n";
	
	out << "\n\tFRACTURES INFO:\n";

	int ifrac = -1;
	for(auto frac : fFractures){
		out << "\n Fracture #" << ++ifrac << "\n";
		frac->Print(out);
	}

	PrintRolodexes(out);
	PrintPolyhedra(out);

}

void DFNMesh::PrintRolodexes(std::ostream & out) const{
	if(fSortedFaces.size() == 0) {
		out << "\n- No rolodexes in this mesh\n"; 
		return;
	}
	out << "\n\nROLODEXES:_____________\n";
	for(TPZGeoEl* edge : fGMesh->ElementVec()){
		if(!edge) continue;
		if(edge->Dimension() != 1) continue;
		if(edge->HasSubElement()) continue;
		if(edge->Index() >= fSortedFaces.size()) break;
		fSortedFaces[edge->Index()].Print(out);
	}
}
void DFNMesh::PrintPolyhedra(std::ostream & out) const{
	if(fPolyhedra.size() == 0) {
		out << "\n- No polyhedra in this mesh\n"; 
		return;
	}
	out <<"\n\nPOLYHEDRA BY FACE:_________________\n";
	{
		out<<"                 (+)  |  (-)";
		// int npolyhedra=0;
		for(TPZGeoEl* gel : fGMesh->ElementVec()){
			if(!gel) continue;
			if(gel->Dimension() != 2) continue;
			if(gel->HasSubElement()) continue;
			if(gel->Index() >= fPolyh_per_face.size()) break;
			
			out << "\nFace #" << std::setw(5)<<std::right << gel->Index() << ":    "; 
			out << std::setw(3)<<std::right << fPolyh_per_face[gel->Index()][0] << "   | ";
			out << std::setw(3)<<std::right << fPolyh_per_face[gel->Index()][1];
			// npolyhedra = std::max({npolyhedra,fPolyh_per_face[gel->Index()][0],fPolyh_per_face[gel->Index()][1]});
		}
		// out<<"\n\nNumber of polyhedra found: "<<npolyhedra+1;
		out<<"\n\nNumber of polyhedra found: "<<fPolyhedra.size();
		for(auto& polyh : fPolyhedra){
			polyh.Print(out);
		}
	}
	out.flush();
}