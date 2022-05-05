/*! 
 *  @brief     Contains implementation of methods for class DFNMesh.
 *  @authors   Pedro Lima
 *  @date      2020.06
 */

#include "DFNMesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZGeoMeshBuilder.h"
#include <filesystem>
#include "pzlog.h"

#if PZ_LOG
static TPZLogger logger("dfn.mesh");
#endif

DFNMesh::DFNMesh(TPZGeoMesh *gmesh, TPZManVector<std::map<int,std::string>,4>& dim_physical_tag_and_name, REAL tolDist, REAL tolAngle, int prerefine){
	fGMesh = gmesh;
	m_dim_physical_tag_and_name = dim_physical_tag_and_name;
    SetTolerances(tolDist,tolAngle);

	// Option to apply a uniform refinement to the coarse mesh before inserting fractures
	PreRefine(prerefine);

	// create a copy of the mesh because materials will be reset afterwards @todo
	// this->ClearMaterials();
	for (int i = 1; i < fGMesh->Dimension(); i++){
		CreateSkeletonElements(i,1);
	}
	
	fSortedFaces.clear();
	fPolyh_per_face.clear();
	fPolyhedra.clear();
	fFractures.clear();
	// fVolumes.clear();
	if(fGMesh->Dimension() == 3){
		InitializePolyhedra();
		// try{InitializePolyhedra();}
		// catch(...){
		// 	DumpVTK();
		// 	PZError << "\nFailed to initialize polyhedra. VTK representation written to vtkmesh.vtk\n";
		// 	throw std::bad_exception();
		// }
		
#if PZ_LOG
		if(logger.isDebugEnabled())
		{
			std::stringstream sout;
			sout << "[Polyhedra initialization]\n";
			this->PrintPolyhedra(sout);
			sout << "\n\n";
			LOGPZ_DEBUG(logger,sout.str());
		}
#endif // PZ_LOG
	}
	else{PZError << "\n2D meshes are currently unsupported\n"; DebugStop();}
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
	// fVolumes = 			copy.fVolumes;
    return *this;
}


DFNMesh::~DFNMesh(){
	for(auto fracture : fFractures){
		delete fracture;
	}
	for(auto volume : fPolyhedra) delete volume;
}


void DFNMesh::PrintVTK(std::string vtkmesh
                    ,std::string pzmesh)
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








void DFNMesh::PrintVTKColorful(std::string vtkmesh,std::string pzmesh){
	// mesh.txt doesn't gain much from coloring... so print it normal first
	this->PrintVTK("skip",pzmesh);
	TPZGeoMesh *gmesh = this->fGMesh;
	int64_t nels = gmesh->NElements();
	int64_t iel;
	TPZGeoEl *gel;

	// Gather set of material ids for DFNFractures
	std::set<int> frac_material;
	for(auto frac : fFractures){
		frac_material.insert(frac->MaterialId());
	}

	// shift material ids
	for(iel = 0; iel < nels; iel++){
		gel = gmesh->Element(iel);
		if(!gel) continue;
		if(gel->Dimension() != 2) continue;
		if(gel->HasSubElement()) continue;
		if(!gel->Father()) continue;
		// if there will be a unique material id for each fracture, we might want to create an std::set<int> with all the material ids of fractures and do the next if as a binary search
		if(frac_material.find(gel->MaterialId()) != frac_material.end()) continue;
		// if(gel->MaterialId() == DFNMaterial::Efracture) continue;
		// if(gel->MaterialId() == DFNMaterial::Eintact) continue;
		int subindex = gel->WhichSubel();
		int matid = gel->MaterialId();
		gel->SetMaterialId(DFNMaterial::Erefined+subindex);
	}
	// print vtk only, since txt has already been print
	this->PrintVTK(vtkmesh,"skip");

	// then, restore original mat ids
	// for(iel = 0; iel < nels; iel++){
	// 	gel = gmesh->Element(iel);
	// 	if(!gel) continue;
	// 	if(gel->Dimension() != 2) continue;
	// 	if(gel->HasSubElement()) continue;
	// 	if(!gel->Father()) continue;
	// 	int subindex = gel->WhichSubel();
	// 	int matid = gel->MaterialId();
	// 	gel->SetMaterialId(matid-subindex);
	// }
}


void InwardNormal(const TPZGeoElSide& gelside, TPZManVector<REAL,3>& qsi, TPZManVector<REAL,3>& normal){
	gelside.Normal(qsi,normal);
	int dim=normal.size();
	for(int i=0; i<dim; i++){
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


DFNFracture* DFNMesh::CreateFracture(DFNPolygon &Polygon, FracLimit limithandling, int materialid){
	// For the export graphics code, we're starting with this limit of max_intersections and reserving material ids within every 1000 for frac_frac_intersections. It's just for graphics
	// int frac_tag = this->NFractures() + 1;
	// materialid = frac_tag*max_intersections;
	DFNFracture* fracture = new DFNFracture(Polygon,this,limithandling, materialid);

#if PZ_LOG
    LOGPZ_INFO(logger, "[Start][Fracture " << fracture->Index() << "]");
#endif // PZ_LOG

	std::cout<<"\nFracture #"<<fracture->Index();
	fFractures.push_back(fracture);
	return fracture;
}



// void DFNMesh::ImportElementsFromGMSH(TPZGeoMesh * gmesh, int dimension, std::set<int64_t> &oldnodes, TPZVec<TPZGeoEl*>& newgels){
//     // GMsh does not accept zero index entities
//     const int shift = 1;

//     // First, if GMsh has created new nodes, these should be inserted in PZGeoMesh
//     // create a map <node,point>
//     std::map<int,int> mapGMshToPZ;

//     for(int64_t pznode : oldnodes){
// 		std::vector<size_t> node_identifiers;
//         std::vector<double> coord;
//         std::vector<double> parametricCoord;
//         gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord,0,pznode+shift,true);
//         int gmshnode = (int) node_identifiers[0];
// 		// insert with hint (since oldnodes is an already sorted set, these nodes will all go in the end)
//         mapGMshToPZ.insert(mapGMshToPZ.end(),{gmshnode,pznode+shift});
// 	}

//     // add new nodes into PZGeoMesh
//     {
//         // get all nodes from GMsh
//             std::vector<size_t> node_identifiers;
//             std::vector<double> coord;
//             std::vector<double> parametricCoord;
//             gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord);
//         // iterate over node_identifiers
//         int nnodes = node_identifiers.size();
//         for(int i = 0; i < nnodes; i++){
//             int gmshnode = node_identifiers[i];
//             // check if it is contained in the map
//             if(mapGMshToPZ.find(gmshnode) == mapGMshToPZ.end()){
//                 // New node -> add to PZGeoMesh
//                 int pznode = (int) gmesh->NodeVec().AllocateNewElement();
//                 TPZManVector<REAL,3> newnodeX(3);
//                 newnodeX[0] = coord[3*i];
//                 newnodeX[1] = coord[3*i+1];
//                 newnodeX[2] = coord[3*i+2];
//                 gmesh->NodeVec()[pznode].Initialize(newnodeX,*gmesh);
//                 // int pznode = (int) gmesh->NNodes();
//                 // gmesh->NodeVec().resize(pznode+1);
//                 // insert it in map
//                 mapGMshToPZ.insert({gmshnode,pznode+shift});
//             }

//         }
//     }
    

    
//     int64_t nels = gmesh->NElements();
//     std::vector<std::pair<int, int> > dim_to_physical_groups;
//     gmsh::model::getPhysicalGroups(dim_to_physical_groups,dimension);
// 	int nnew = 0;
//     /// inserting the elements
//     for (auto group: dim_to_physical_groups) {
       
//         int dim = group.first;
//         // only want elements of a given dimension
//         if(dim != dimension) continue;
//         int physical_identifier = group.second; 
       
//         std::vector< int > entities;
//         gmsh::model::getEntitiesForPhysicalGroup(dim, physical_identifier, entities);

// 		for (auto tag: entities) {
// 		// std::cout<<"______________________test - tag = "<<tag;
           
//             std::vector<int> group_element_types;
//             std::vector<std::vector<std::size_t> > group_element_identifiers;
//             std::vector<std::vector<std::size_t> > group_node_identifiers;
//             gmsh::model::mesh::getElements(group_element_types,group_element_identifiers,group_node_identifiers, dim, tag);
//             int n_types = group_element_types.size();
//             for (int itype = 0; itype < n_types; itype++){
//                 int el_type = group_element_types[itype];
//                 int n_nodes = TPZGeoMeshBuilder::GetNumberofNodes(el_type);
//                 std::vector<int> node_identifiers(n_nodes);
//                 int n_elements = group_element_identifiers[itype].size();
//                 for (int iel = 0; iel < n_elements; iel++) {
//                     int el_identifier = group_element_identifiers[itype][iel]+nels-gmshshift;
// 					// std::cout<<"\n"<<el_identifier<<"\n";

//                     for (int inode = 0; inode < n_nodes; inode++) {
//                         // node_identifiers[inode] = group_node_identifiers[itype][iel*n_nodes+inode];
//                         // Use mapGMshToPZ to translate from GMsh node index to PZ nodeindex
//                         node_identifiers[inode] = mapGMshToPZ[group_node_identifiers[itype][iel*n_nodes+inode]];
//                     }
//                     TPZGeoMeshBuilder::InsertElement(gmesh, physical_identifier, el_type, el_identifier, node_identifiers);
// 					nnew++;
// 					newgels.resize(nnew);
// 					newgels[nnew-1] = gmesh->Element(el_identifier);
// 				}
//             }
//         }
//     }
//     gmesh->BuildConnectivity();
// }





























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









int DFNMesh::RealNFractures() const{
	int counter=0;
	for(auto frac : fFractures){
		counter += (frac->NSurfElements() > 0);
	}
	return counter;
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








void DFNMesh::ExportGMshCAD_nodes(std::ofstream& out){
	const float meshsize = 0.0;
    // write nodes
    out<< "// POINTS DEFINITION \n\n";
    out<< "h = " << meshsize << ";\n\n";
    int64_t nnodes = fGMesh->NNodes();
    for (int64_t inode = 0; inode < nnodes; inode++){
		if(fGMesh->NodeVec()[inode].Id() < 0) continue;
        TPZManVector<REAL, 3> co(3,0.);
        fGMesh->NodeVec()[inode].GetCoordinates(co);
        out << "Point(" << inode+gmshshift << ") = {" 
		    << co[0] << ',' 
			<< co[1] << ',' 
			<< co[2] 
			<< ", h"
			<<"};\n";
    }
}





void DFNMesh::ExportGMshCAD_edges(std::ofstream& out){
	out << "\n\n// LINES DEFINITION \n\n";
	// declare lists to define physical groups
	// std::list<int64_t> groupFracture;

	int64_t nels = fGMesh->NElements();

	// iterate over all 1D elements
	for (int64_t iel = 0; iel < nels; iel++){
		TPZGeoEl *gel = fGMesh->Element(iel);
		if(!gel) continue;
		if(gel->Dimension() != 1) continue;
		if(gel->HasSubElement()) continue;
		
		out << "Line(" << iel+gmshshift << ") = {" << gel->NodeIndex(0)+gmshshift << ',' << gel->NodeIndex(1)+gmshshift << "};\n";
		// list it according to material
		// if(gel->MaterialId() == DFNMaterial::Efracture){groupFracture.push_back(iel+gmshshift);}
	}
	// write physical groups
	// if(fGMesh->Dimension() == 2){
	// 	out<<"\nPhysical Curve(\"Fractures\","<<DFNMaterial::Efracture<<") = {";
	// 	for(auto itr = groupFracture.begin(); itr != groupFracture.end();/*Increment in next line*/){
	// 		out<<*itr<<(++itr!=groupFracture.end()? "," : "};\n");
	// 	}
	// }
}










void DFNMesh::ExportGMshCAD_faces(std::ofstream& out){
	std::stringstream stream;
    stream << "\n\n// FACES DEFINITION \n\n";
	int64_t nels = fGMesh->NElements();
	int mesh_dim = fGMesh->Dimension();
	// // declare lists to define physical groups
	// std::list<int64_t> groupFracture;
	// std::list<int64_t> groupTransition;
	// std::list<int64_t> groupIntact;

	TPZManVector<int64_t,4> lineloop;
	// iterate over all 2D elements
	for (int64_t iel = 0; iel < nels; iel++){
		TPZGeoEl *gel = fGMesh->Element(iel);
		if(!gel) continue;
		if(gel->Dimension() != 2) continue;
		if(gel->HasSubElement()) continue;

		// curve loop _______________________________________
		stream << "Curve Loop(" << iel+gmshshift << ") = {";
		DFN::GetLineLoop(gel,lineloop,gmshshift);
		int nedges = gel->NSides(1);
		for(int i=0; i<nedges; i++){
			stream << lineloop[i] << ",";
		}
		stream.seekp(stream.str().length()-1);
		stream << "};\n";
		// surface
		stream << "Surface("<<iel+gmshshift<<") = {"<<iel+gmshshift<<"};\n";

		// if(gel->MaterialId() == DFNMaterial::Efracture){groupFracture.push_back(iel+gmshshift);}
		// if(mesh_dim == 2){
		// 	if(gel->MaterialId() == DFNMaterial::Erefined){groupTransition.push_back(iel+gmshshift);}
		// 	else groupIntact.push_back(iel+gmshshift);
		// }
	}



	
	// if(groupTransition.size() > 0){
	// 	stream<<"\nPhysical Surface("<<DFNMaterial::Erefined<<") = {";
	// 	for(auto index : groupTransition){
	// 		stream<<index<<",";
	// 	}
	// 	stream.seekp(stream.str().length()-1);
	// 	stream << "};\n";
	// }
	// if(groupFracture.size() > 0){
	// 	if(mesh_dim == 3) stream<<"\nPhysical Surface(\"Fractures\","<<DFNMaterial::Efracture<<") = {";
	// 	for(auto index : groupFracture){
	// 		stream<<index<<",";
	// 	}
	// 	stream.seekp(stream.str().length()-1);
	// 	stream << "};\n";
	// }
	// if(groupIntact.size() > 0){
	// 	stream<<"\nPhysical Surface("<<DFNMaterial::Eintact<<") = {";
	// 	for(auto index : groupIntact){
	// 		stream<<index<<",";
	// 	}
	// 	stream.seekp(stream.str().length()-1);
	// 	stream << "};\n";
	// }
	out << stream.str() << std::endl;
}












void DFNMesh::ExportGMshCAD_volumes(std::ofstream& out){
	/** @note: There's a seemingly strange behaviour which should be remembered in the future.
	 * If element of index zero is a coarse volume, the entry for its coarse grouping in the exported .geo is likely going to be
	 * Physical Volume("c1",1) = {2};
	 * This is the desired behaviour. This is the third time I've had to re-explain it to myself so I might aswell write it down.
	 * c1 = coarse element of index zero + gmshshift
	 * 1  = coarse element of index zero + gmshshift
	 * 2  = the polyhedral index + gmshshift
	 * polyhedral index is not zero because zero is reserved for the 'boundary polyhedron'
	 * Keeping gmshshift for every index, regardless of its necessity or lack of it, makes the code less error-prone since we know that everything needs to gmshshift back later
	 **/


	// We'll need a multimap to define physical tags for coarse elements
	// maps <coarse index, polyh index>
	std::multimap<int64_t,int> coarse_groups;

	std::stringstream stream;
	stream << "\n\n// VOLUMES DEFINITION \n\n";
	// declare lists to define physical groups
	std::list<int64_t> groupTransition;
	std::list<int64_t> groupIntact;
	int nels = fPolyhedra.size();
	for(int ipolyh=1; ipolyh < nels; ipolyh++){
		DFNPolyhedron& polyh = Polyhedron(ipolyh);
		if(polyh.IsRefined()) continue;

		stream << "Surface Loop(" << polyh.Index()+gmshshift << ") = {";
		for(auto& face_orient : polyh.Shell()){
			stream << face_orient.first+gmshshift << ",";
		}
		stream.seekp(stream.str().length()-1);
		stream << "};\n";
		stream << "Volume(" << polyh.Index()+gmshshift << ") = {" << polyh.Index()+gmshshift << "};\n";

		coarse_groups.insert({polyh.CoarseIndex(),ipolyh});
	}

	// Give physical tags to differenciate coarse element groups
	stream << "\n\n// COARSE ELEMENTS GROUPING\n";
	int countEls = 0;
	if(coarse_groups.size()){
		int current_coarse = coarse_groups.begin()->first;
		std::stringstream groupname;
		groupname << "\"c" << countEls++ << "\",";
		stream 	<< "\nPhysical Volume("
				<< groupname.str()
				<< current_coarse+gmshshift
				<<") = {";
		for(auto& itr : coarse_groups){
			if(itr.first != current_coarse){
				current_coarse = itr.first;
				groupname.str(std::string());
				groupname << "\"c" << countEls++ << "\",";
				stream.seekp(stream.str().length()-1);
				stream << "};";
				stream 	<< "\nPhysical Volume("
						<< groupname.str()
						<< current_coarse+gmshshift
						<<") = {";
			}
			stream << itr.second +gmshshift << ",";
		}
		stream.seekp(stream.str().length()-1);
		stream << "};\n";
	}
	out << stream.str() << std::endl;
}











/**
 * 	@brief Creates a .geo for the mesh
 */ 
void DFNMesh::ExportGMshCAD(std::string filename){
	std::cout<<" -Exporting Gmsh CAD file\r"<<std::flush;
	// fGMesh->BuildConnectivity();
	// Surface cleanup (necessary if limit recovery is done after all fractures have been inserted to the mesh)
	for(auto frac : FractureList()){
		frac->CleanUp(frac->MaterialId());
	}
	// gmsh doesn't accept index zero elements
	// const int shift = 1;
	std::ofstream out(filename);

    // Title
    out<<"//  Geo file generated by DFNMesh project\n//\thttps://github.com/labmec/dfnMesh";

	// write nodes
    ExportGMshCAD_nodes(out);
    // write edges
    ExportGMshCAD_edges(out);
    // write faces
    ExportGMshCAD_faces(out);
    // write volumes
    ExportGMshCAD_volumes(out);
	// write fractures
    ExportGMshCAD_fractures(out);
	// boundary conditions
	ExportGMshCAD_boundaryconditions(out);
	ExportGMshCAD_fractureIntersections(out);
    
    
	out << "\n// OPTIONS\n";
    // out<<"\nTransfinite Surface {Physical Surface("<<DFNMaterial::Eintact<<")};";
    // out<<"\nRecombine Surface {Physical Surface("<<DFNMaterial::Eintact<<")};";
    // out<<"Recombine Surface {Physical Surface("<<DFNMaterial::Erefined<<")};\n";
	// out << "\nCoherence;";
	out << "\nCoherence Mesh;";
	out << "\nTransfinite Curve {:} = 2;";
	out << "\nTransfinite Surface{:};";
	out << "\nTransfinite Volume{:};";
	out << "\nRecombine Surface{:};";
	out << "\nRecombine Volume{:};";
	out.flush();
	std::cout<<"                         \r"<<std::flush;
}

void DFNMesh::ExportGMshCAD_fractures(std::ofstream& out){
	if(this->RealNFractures() == 0) return;
	std::stringstream stream;

	stream << "\n\n // FRACTURES\n";
	// write physical groups for fractures
	int ifrac = 0;
	std::multimap<int,int> fracturegroups;
	// loop over fractures and list the elements on their surface
	for(DFNFracture* frac : this->fFractures){
		ifrac++;
		if(frac->Surface().size() == 0) continue;
		stream << "\nfrac" << frac->Index() << "[] = {";
		for(int64_t gelindex : frac->Surface()){
			stream << gelindex+gmshshift << ",";
		}
		stream.seekp(stream.str().length()-1);
		stream << "};";
		fracturegroups.insert({frac->Index(),frac->MaterialId()});
		// stream<<"Physical Surface(\"Fracture "<< ifrac << "\","<<frac->MaterialId()<<") = frac" << ifrac <<"[];\n";
	}
	// If fractures share a material id, we can (must?) merge their listings into the same physical group
	int igroup = fracturegroups.begin()->first;
    int index = fracturegroups.begin()->second;
	stream 	<< "\n\nPhysical Surface(\"Fracture"
            << igroup
			<< "\","
			<< index
			<<") = {";
	for(auto& itr : fracturegroups){
		if(itr.first != igroup){
			igroup = itr.first;
            index = itr.second;
			stream.seekp(stream.str().length()-2);
			stream << "};";
			stream 	<< "\nPhysical Surface(\"Fracture"
					<< igroup
					<< "\","
					<< index
					<<") = {";
		}
		stream << "frac" << igroup  << "[], ";
	}
	stream.seekp(stream.str().length()-2);
	stream << "};\n";
	
	
	out << stream.str() << std::endl;
}
void DFNMesh::ExportGMshCAD_boundaryconditions(std::ofstream& out){
	std::multimap<int, int64_t> physicalgroups;
    if(fPolyhedra.size() == 0)
    {
        std::cout << __PRETTY_FUNCTION__ << "does not contain a boundary polyhedron\n";
        return;
    }
	// Use the boundary polyhedron to get boundary conditions
	DFNPolyhedron& boundary = Polyhedron(0);
	for(auto orientedface : boundary.Shell()){
		int64_t index = orientedface.first;
		int matid = fGMesh->Element(index)->MaterialId();
		physicalgroups.insert({matid,index});
	}
	/** @note I was initially tempted to have this. It is a non-robust temporary patch, though. 
	  * It can fail when a volume is broken to tetrahedra for some reason.
	  * (see DFNFracture::CheckSnapInducedOverlap and DFNMesh::UpdatePolyhedra for examples)
	  * The actual robust solution is the one above, using the boundary polyhedron.
	  * for(TPZGeoEl* gel : Mesh()->ElementVec()){
	  * 	if(!gel) continue;
	  * 	if(gel->Dimension() != 2) continue;
	  * 	TPZGeoEl* ancestor = (gel->Father() ? gel->EldestAncestor() : gel);
	  * 	TPZGeoElSide gelside(ancestor, ancestor->NSides()-1);
	  * 	if(gelside.NNeighbours(3) != 1) continue;
	  * 	physicalgroups.insert({gel->MaterialId(),gel->Index()});
	  * }
	  */

	// Give physical tags to differenciate coarse element groups
	std::stringstream stream;
	stream << "\n\n// BOUNDARY CONDITIONS\n";
	if(physicalgroups.size()){
		int current_bc = physicalgroups.begin()->first;
		std::stringstream bcname;
		std::map<int,std::string>::const_iterator it = m_dim_physical_tag_and_name[2].find(current_bc);
		if(it != m_dim_physical_tag_and_name[2].end())
			bcname << "\"" <<  it->second << "\",";
		else
			bcname << "\"bc" << current_bc << "\",";
		
			
		stream 	<< "\nPhysical Surface("
				<< bcname.str()
				<< current_bc
				<<") = {";
		for(auto& itr : physicalgroups){
			if(itr.first != current_bc){
				current_bc = itr.first;
				stream.seekp(stream.str().length()-1);
				it = m_dim_physical_tag_and_name[2].find(current_bc);
				bcname.str(std::string());
				if(it != m_dim_physical_tag_and_name[2].end())
					bcname << "\"" <<  it->second << "\",";
				else
					bcname << "\"bc" << current_bc << "\",";
				
				stream << "};";
				stream 	<< "\nPhysical Surface("
						<< bcname.str()
						<< current_bc
						<<") = {";
			}
			stream << itr.second +gmshshift << ",";
		}
		stream.seekp(stream.str().length()-1);
		stream << "};\n";
	}
	out << stream.str() << std::endl;
	
	for(DFNFracture* fracture : fFractures){
		int bcindex = 10;
		const int bcFracMatId = fracture->MaterialId() + 1;
		fracture->ExportFractureBC(bcindex+fracture->Index(),out);
		out << "Physical Curve(\"BCfrac" 
			<< fracture->Index() 
			<< "\", " << bcFracMatId
			<< ") = {BCfrac"
			<< fracture->Index()
			<< "[]};\n";
	}
}





void DFNMesh::ExportGMshCAD_fractureIntersections(std::ofstream& out){
	/// A geometrical intersection between 2 bounded planes in R^3 is a line segment, 
	/// so we represent it by the coordinates of its nodes
	Segment int_segment;
	out << "\n// INTER-FRACTURE INTERSECTIONS\n";
	const int nfrac = this->NFractures();

	// Build the set with every 1D element at the surface of each fracture
    for(auto fracture : fFractures){
		fracture->GetEdgesInSurface();
	}

	// Test every pair of fractures for intersection
	for(int jfrac = 0; jfrac<nfrac; jfrac++){
		for(int kfrac = jfrac+1; kfrac<nfrac; kfrac++){
			// const DFNPolygon& jpolygon = fFractures[jfrac]->Polygon();
			// const DFNPolygon& kpolygon = fFractures[kfrac]->Polygon();
			// bool geom_intersection_Q = jpolygon.ComputePolygonIntersection(kpolygon,int_segment);
			// if(!geom_intersection_Q) continue;
			
			TPZStack<int64_t> intersection_edges;
			fFractures[jfrac]->FindFractureIntersection(*fFractures[kfrac],intersection_edges);
			if(intersection_edges.size() == 0) continue;
			std::stringstream stream;
			// Setup Physical Groups for each intersection
			stream << "\nfracIntersection_" << jfrac << '_' << kfrac << "[] = { ";
			for(int64_t iedge : intersection_edges){
				stream << iedge+gmshshift << ',';
			}
			stream.seekp(stream.str().length()-1);
			stream << "};";
			stream << "\nfracIntersection_" << kfrac << '_' << jfrac << "[] = fracIntersection_" << jfrac << '_' << kfrac << ';';
			stream <<"\nPhysical Curve(\"fracIntersection_" << jfrac << '_' << kfrac<<"\") = "
								  << "fracIntersection_" << jfrac << '_' << kfrac<<"[];";
			out << stream.str() << std::endl;
		}
	}
	out << std::endl;
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
	if(dimension == 0) return;
	std::cout << " -Creating skeleton elements\r" << std::flush;
    int nel = fGMesh->NElements();
    for (int iel = 0; iel < nel; iel++)
    {
        TPZGeoEl *gel = fGMesh->Element(iel);
        if(!gel) continue;
        // Elements can't have a skeleton of higher dimension than itself
        if(gel->Dimension() <= dimension) continue;
		if(gel->HasSubElement()) continue;
        
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
                // if(gel->MaterialId() == DFNMaterial::Efracture) TPZGeoElBC(gelside, DFNMaterial::Efracture);
                TPZGeoElBC(gelside, matid);
            }
        }
    }
	std::cout << "                            \r" << std::flush;
}








void DFNMesh::ClearMaterials(int matid){
	std::cout << "\n\n[WARNING] Boundary conditions were erased to print VTK graphics.\n\n" << std::flush;
	for(auto el : fGMesh->ElementVec()){
		if(!el) continue;
		el->SetMaterialId(matid);
	}
}
void DFNMesh::ClearMaterials(const int matid, TPZVec<int>& backup){
	backup.Resize(fGMesh->NElements(),-1);
	for(auto el : fGMesh->ElementVec()){
		if(!el) continue;
		backup[el->Index()] = el->MaterialId();
		el->SetMaterialId(matid);
	}
}

void DFNMesh::RestoreMaterials(TPZVec<int>& backup){
	backup.Resize(fGMesh->NElements(),-1);
	for(auto el : fGMesh->ElementVec()){
		if(!el) continue;
		el->SetMaterialId(backup[el->Index()]);
	}
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










DFNPolyhedron* DFNMesh::CreatePolyhedron(const TPZVec<std::pair<int64_t,int>>& shell,int64_t coarseindex, bool isConvex){
	int ipolyh = fPolyhedra.size();
	fPolyhedra.resize(ipolyh+1);
	DFNPolyhedron* newpolyhedron = new DFNPolyhedron;
	fPolyhedra[ipolyh] = newpolyhedron;
	newpolyhedron->Initialize(this,ipolyh,shell,coarseindex,isConvex);
#if PZ_LOG
	if(logger.isDebugEnabled()){
		std::stringstream stream;
		stream << "[Adding Polyhedron]";
		newpolyhedron->Print(stream);
		LOGPZ_DEBUG(logger,stream.str());
	}
#endif // PZ_LOG
	return newpolyhedron;
}










void DFNMesh::InitializePolyhedra(){
	// Clear data structure
	fPolyhedra.clear();
	fPolyh_per_face.clear();
	fPolyh_per_face.Resize(fGMesh->NElements(),{-1,-1});
	int ipolyh=0; 											///< index of polyhedron
	TPZStack<std::pair<int64_t,int>,20> shell(20,{-1,0});	///< oriented faces that enclose the polyhedron

	std::cout<<" -Initializing polyhedral volumes\r"<<std::flush;

	// Gather boundary faces first
	shell.clear();
	TPZGeoElSide gelside;
	TPZGeoElSide neig;
	for(TPZGeoEl* gel : fGMesh->ElementVec()){
		if(gel->Dimension() != 2) continue;
		if(gel->HasSubElement()) continue;
		// if(gel->Type() != EQuadrilateral) continue; //testing
		// if(gel->MaterialId() != DFNMaterial::Eskeleton) continue;

		gelside = {gel,gel->NSides()-1};
		int nneighbours = 0;
		TPZGeoElSide volneig;
		for(neig=gelside.Neighbour(); neig != gelside; neig = neig.Neighbour()){
			if(neig.Element()->Dimension() < 3) continue;
			nneighbours += neig.Element()->Dimension() > 2;
			volneig = neig; //< catch a volume in case there's a boundary condition 
		}
		if(nneighbours != 1) continue;
		// Faces don't have guaranteed positive orientation, since they were created using CreateBCGeoEl(). @see  DFN::CreateSkeletonElements()
		int orient = DFN::SkeletonOrientation(volneig,gel);
		shell.push_back({gel->Index(),orient});
		SetPolyhedralIndex({gel->Index(),orient},0);
	}
	CreatePolyhedron(shell,-1,false);
#ifdef PZDEBUG
	shell.Fill({-1,0});
#endif //PZDEBUG
	shell.clear();

	// Then initialize the rest of the polyhedra
	// Looping over 3D elements and matching a polyhedral index to their shell
	for(TPZGeoEl* vol : fGMesh->ElementVec()){
		if(!vol) continue;
		if(vol->Dimension() != 3) continue;
		if(vol->HasSubElement()) continue;
		TPZGeoEl* ancestor = vol->EldestAncestor();
		CreateGelPolyhedron(vol,ancestor->Index());
	}

	std::cout<<"                                 \r"<<std::flush;
}

int DFNMesh::CreateGelPolyhedron(TPZGeoEl* vol, int coarseindex){
	if(vol->Dimension() != 3) return -1;
	TPZStack<std::pair<int64_t,int>,6> shell(6,{-1,0});	///< oriented faces that enclose the polyhedron
	shell.clear();
	int ipolyh = fPolyhedra.size();
	
	#if PZDEBUG
		// TPZManVector<REAL,3> qsi(3,0.33333333); Matrix jac; Matrix axes; REAL detjac; Matrix jacinv;
		// vol->Jacobian(qsi,jac, axes, detjac, jacinv);
		// LOGPZ_DEBUG(logger, "detjac vol " << vol->Index() << " = " << detjac)
		// if(detjac < 0.) DebugStop();
	#endif //PZDEBUG
	// Loop over 2D sides
	for(int iside=vol->FirstSide(2); iside < vol->NSides()-1; iside++){
		// Get oriented face in that side
		TPZGeoElSide volside(vol,iside);
		TPZGeoEl* face = DFN::GetSkeletonNeighbour(vol,iside);
		int orient = -DFN::SkeletonOrientation(volside,face);
		std::pair<int64_t,int> face_orient = {face->Index(),orient};
		shell.push_back(face_orient);
		// Set polyhedral index to that oriented face
		int currentface_polyh_index = GetPolyhedralIndex(face_orient);
		// if(currentface_polyh_index > -1 && currentface_polyh_index != ipolyh) DFN_DebugStop();
		SetPolyhedralIndex(face_orient,ipolyh);
	}
	// DFNPolyhedron polyhedron(this,ipolyh,shell,coarseindex);
	// polyhedron.Print();
	// fPolyhedra.push_back(polyhedron);
	CreatePolyhedron(shell,coarseindex);
	shell.Fill({-1,0});
	return ipolyh;
}

void DFNMesh::UpdatePolyhedra(){
#if PZ_LOG
	if(logger.isInfoEnabled()) LOGPZ_INFO(logger,"[Start][Update Polyhedra]");
#endif // PZ_LOG
	// Start by sorting faces around edges and filling the this->fSortedFaces datastructure
	this->SortFacesAroundEdges();
	std::cout<<"-Updating polyhedral volumes\r"<<std::flush;
	fPolyh_per_face.Resize(fGMesh->NElements(),{-1,-1});
	// Refined faces pass down their polyh index to their subelements
	InheritPolyhedra();
	int polyh_index = fPolyhedra.size();

	
	TPZStack<std::pair<int64_t,int>,20> polyhedron(20,{-1,0});
	// loop over 2D skeleton elements
	for(int64_t iel=0; iel < fGMesh->NElements(); iel++){
		TPZGeoEl* initial_face = fGMesh->Element(iel);
		if(!initial_face) continue;
		if(initial_face->Dimension() != 2) continue;
		int64_t initial_id = initial_face->Index();
		
		// look for polyhedron on each orientation of the initial_face
		for(int ipolyh_local=0; ipolyh_local<2; ipolyh_local++){
			if(initial_face->HasSubElement()) break;
			int orientation = ipolyh_local?-1:1;
			// if the first has been found, continue to the next
			if(GetPolyhedralIndex({initial_id,orientation}) != -1) continue;
			polyhedron.Fill({-1,0});
			polyhedron.clear();
			std::pair<int64_t,int> initial_face_orient = {initial_id,orientation};
			polyhedron.push_back({initial_id,orientation});
			polyh_index = fPolyhedra.size();
			SetPolyhedralIndex(initial_face_orient,polyh_index);
			bool IsConvex = true;
			int coarseindex = -1;

			BuildVolume(initial_face_orient,IsConvex,polyhedron,coarseindex);
			DFNPolyhedron* newvolume = CreatePolyhedron(polyhedron,coarseindex,IsConvex);

			#ifdef PZDEBUG
				if(coarseindex < 0){
					this->DumpVTK(true);
					PZError << "\nFailed to attribute a coarse element index to: ";
					newvolume->Print(PZError,true);
					newvolume->PlotVTK("./LOG/OrphanVolume.vtk");
					DebugStop();
				}
			#endif //PZDEBUG

			if(!IsConvex) {
				ClearPolyhIndex(polyhedron);
				newvolume->Refine();
				this->SortFacesAroundEdges();
				--ipolyh_local;
			}
		}
	}
	std::cout<<"                               \r"<<std::flush;
#if PZ_LOG
	if(logger.isInfoEnabled()) LOGPZ_INFO(logger,"[End][Update Polyhedra]");
	if(logger.isDebugEnabled()){
		std::stringstream stream;
		stream << "[Result][Update Polyhedra]\n";
		PrintPolyhedra(stream);
		LOGPZ_DEBUG(logger,stream.str());
	}
#endif // PZ_LOG
}

void DFNMesh::ClearPolyhIndex(TPZVec<std::pair<int64_t,int>>& facestack){
	for(auto& faceorient : facestack){
		SetPolyhedralIndex(faceorient, -1);
	}
}





void DFNMesh::SetPolyhedralIndex(std::pair<int64_t,int> face_orient, int polyh_index){
	switch(face_orient.second){
		case  1: fPolyh_per_face[face_orient.first][0] = polyh_index; break;
		case -1: fPolyh_per_face[face_orient.first][1] = polyh_index; break;
		default: DebugStop();
	}
}
int DFNMesh::GetPolyhedralIndex(int64_t faceindex, int orientation) const{
	return GetPolyhedralIndex(std::make_pair(faceindex,orientation));
}
int DFNMesh::GetPolyhedralIndex(const std::pair<int64_t,int>& face_orient) const{
	int polyh_index = -1;
	if(face_orient.first >= fPolyh_per_face.size()) return -1;
	switch(face_orient.second){
		case  1: polyh_index = fPolyh_per_face[face_orient.first][0]; break;
		case -1: polyh_index = fPolyh_per_face[face_orient.first][1]; break;
		default: DebugStop();
	}
	return polyh_index;
}

template<int Talloc>
void DFNMesh::BuildVolume(std::pair<int64_t,int> initial_face_orient, bool& IsConvex, TPZStack<std::pair<int64_t,int>, Talloc>& polyhedron, int& coarseindex){
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
		try{
            cards[i] = rolodex.Card(initial_face_orient.first);
        }
		catch(...){
			PrintRolodexBugReport(iedge);
            PrintProblematicRolodex(initial_face_orient.first,rolodex);
            DFN_DebugStop();
        }
	}

	// Determine orientation of neighbour cards
	TPZManVector<std::pair<TRolodexCard, int>,4> facingcards(edges.size());
	for(int i=0; i<edges.size(); i++){
		int64_t iedge = edges[i];
		TRolodex& rolodex = fSortedFaces[iedge];
		std::pair<TRolodexCard, int> current_card = {cards[i],initial_face_orient.second};
		REAL angle = 0.0;
		facingcards[i] = rolodex.FacingCard(current_card,angle);
		if(angle > M_PI+gDFN_SmallNumber)
			{IsConvex = false; }
		if(DFN::Is2PIorZero(angle)){
			PZError << "\nCollapsed polyhedron:\n"
					<< polyhedron
					<< "\nCollapse happens between faces:\n"
					<< initial_face_orient << "\n"
					<< std::make_pair(facingcards[i].first.fgelindex , facingcards[i].second)
					<< "\nAround rolodex:\n"
					<< rolodex;
			PrintRolodexBugReport(iedge);
			DFN_DebugStop();
		}
	}
	// queue neighbour cards to verify
	TPZStack<std::pair<int64_t, int>> to_verify;
	for(int i=0; i<edges.size(); i++){
		int64_t iedge = edges[i];
		std::pair<int64_t,int> faceorient = {facingcards[i].first.fgelindex,facingcards[i].second};
		int nextface_polyindex = GetPolyhedralIndex(faceorient);
		if(nextface_polyindex != polyh_index){
			if(nextface_polyindex > 0) {coarseindex = Polyhedron(nextface_polyindex).CoarseIndex();}
			to_verify.push_back(faceorient);
			SetPolyhedralIndex(faceorient,polyh_index);
			polyhedron.push_back(faceorient);
		}
	}
	// recursively try to add neighbours of cards that have been queued 
	for(auto orientedface : to_verify){
		BuildVolume(orientedface,IsConvex,polyhedron,coarseindex);
	}

}

/** Check neighbours to see if there are overlapping elements in the plane of this element
 * We're looking for a pair of triangles overlapping a quadrilateral
*/
TPZStack<TPZGeoEl*,2> GetOverlappedElements(TPZGeoEl* gel){
	// DebugStop(); // Under construction
	if(gel->Dimension() != 2){PZError << "\nCurrently limited to 2D elements\n"; DebugStop();}
	TPZStack<TPZGeoEl*,2> overlapped;
	
	TPZGeoMesh* gmesh = gel->Mesh();
	int nsides = gel->NSides();
	TPZGeoElSide neig(nullptr,0);
	int dim = gel->Dimension();

	// For every 2 consecutive sides, see if there's a common neighbour between them
	for(int iside = gel->NSides(0); iside < nsides-1; iside++){
		TPZGeoElSide gelside1(gel,iside);
		TPZGeoElSide gelside2(gel,(iside+1)%gel->NSides(1)+gel->NSides(1));
		std::set<int64_t> neighbours1;
		std::set<int64_t> neighbours2;
		// gather neighbours for gelside1
		for(neig = gelside1.Neighbour(); neig != gelside1; ++neig){
			TPZGeoEl* neig_el = neig.Element();
			if(neig_el->Dimension() != dim) continue;
			// neighbours1.insert(neig_el);
			neighbours1.insert(neig_el->Index());
		}
		if(neighbours1.size() < 1) continue;
		// gather neighbours for gelside2
		for(neig = gelside2.Neighbour(); neig != gelside2; ++neig){
			TPZGeoEl* neig_el = neig.Element();
			if(neig_el->Dimension() != dim) continue;
			// neighbours2.insert(neig_el);
			neighbours2.insert(neig_el->Index());
		}
		if(neighbours2.size() < 1) continue;

		// get intersection
		std::set<int64_t> common = DFN::set_intersection(neighbours1,neighbours2);
		// a non-empty intersection means an overlapped element
		if(common.size() == 0) continue;
		if(common.size() >  1) {Print(common) ; DebugStop();}
		
		// overlapped.push_back(*(common.begin()));
		overlapped.push_back(gmesh->Element(*(common.begin())));
		// iside++; // assuming the overlap is 2 triangles, we can skip a side if an overlap was found
	}

	return overlapped;
}



// template<int Talloc>
// void DFNMesh::RefineQuads(TPZStack<std::pair<int64_t,int>,Talloc>& polyhedron){
void DFNMesh::RefineQuads(const TPZVec<std::pair<int64_t,int>>& polyhedron){
	// DebugStop(); // Under construction
	// TPZStack<std::pair<int64_t,int>,Talloc> aux;
	for(auto& orient_face : polyhedron){
		TPZGeoEl* gel = fGMesh->Element(orient_face.first);
		if(gel->HasSubElement()) DFN_DebugStop();
		if(gel->Type() == MElementType::ETriangle) {continue;}
		
		// Get index of the polyhedron on the other side of this gel to pass it to the children
		int permut_orient = -orient_face.second;
		int otherpolyh = GetPolyhedralIndex({gel->Index(),permut_orient});

		// Get the elements overlapped by gmsh
		TPZStack<TPZGeoEl*,2> children = GetOverlappedElements(gel);
		// Create a refpattern
		DFN::CreateRefPattern(gel,children);
		// Setup refinement tree
		for(int i=0; i<children.size(); i++){
			gel->SetSubElement(i,children[i]);
			children[i]->SetMaterialId(gel->MaterialId());
			children[i]->SetFather(gel);
			// It's important that children inherit the polyhedron on the other side
			int child_orient = permut_orient*DFN::SubElOrientation(gel,i);
			SetPolyhedralIndex({children[i]->Index(),child_orient},otherpolyh);
		}
		InheritPolyhedra(gel);
	}
	
}


// template<int Talloc>
// void DFNMesh::MeshPolyhedron(TPZStack<std::pair<int64_t,int>,Talloc>& polyhedron, int coarseindex){
void DFNMesh::MeshPolyhedron(const TPZVec<std::pair<int64_t,int>>& polyhedron, int coarseindex, TPZStack<int64_t>& newgels){
	// GMsh doesn't like zero index entities
    constexpr int shift = 1;

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
	constexpr char mshfilename[] = "LOG/gmshAPI_LastVolumeMeshed.msh";
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
		// gmsh::model::addPhysicalGroup(2,surfaceloop,DFNMaterial::Erefined);
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
	DFN::ImportElementsFromGMSH(fGMesh,3,old_nodes,newgels);
	gmsh::model::remove();
	gmsh::clear();
	// gmsh::finalize();

	// New volumetrical elements should have skeleton elements
	CreateSkeletonElements(1,DFNMaterial::Eintact);
	CreateSkeletonElements(2,DFNMaterial::Eintact);
	fPolyh_per_face.Resize(fGMesh->NElements(),{-1,-1});

	// Refine quadrilaterals down to 2 triangles so not to depend on pyramids to mesh the polyhedron
	// I've written a code that does this at the end in order to give gmsh the freedom to optimize the volumetrical mesh as it sees fit
	RefineQuads(polyhedron);

	// Make sure these new 3D elements have corresponding polyhedra, and coarse index is inherited
	for(int64_t index : newgels){
		TPZGeoEl* vol = fGMesh->Element(index);
		int newid = CreateGelPolyhedron(vol,coarseindex);
	}
}










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
        const int dimOfNeigh = 2;
		int nfaces = edgeside.NNeighbours(dimOfNeigh);
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
            std::map<REAL,TPZGeoElSide>::iterator it = facemap.find(angle);
             if (it != facemap.end()) {
                 rolodex.Initialize(gel->Index(),facemap);
                 PrintProblematicRolodex(neig.Element()->Index(), rolodex);
                 DebugStop(); // two faces with same angle in rolodex!
             }
			facemap.insert({angle,neig});
		}
		rolodex.Initialize(gel->Index(),facemap);
	}
	std::cout<<"                 \r"<<std::flush;
}


void DFNMesh::InheritPolyhedra(){
	fPolyh_per_face.Resize(fGMesh->NElements(),{-1,-1});

	for(TPZGeoEl* father : fGMesh->ElementVec()){
		if(!father) continue;
		if(father->Dimension() != 2) continue;
		if(!father->HasSubElement()) continue;
		InheritPolyhedra(father);
	}
}

void DFNMesh::InheritPolyhedra(TPZGeoEl* father){
	if(!father->HasSubElement()) return;
	TPZGeoEl* child = nullptr;
	int nchildren = father->NSubElements();

	// Gather children orientation
	TPZManVector<int,4> childorient(nchildren,0);
	for(int jchild=0; jchild<nchildren; jchild++){
		childorient[jchild] = DFN::SubElOrientation(father,jchild);
	}

	for(int i=0; i<2; i++){
		int fatherorient = i?-1:1;
		int father_polyhindex = GetPolyhedralIndex({father->Index(),fatherorient});
		if(father_polyhindex<0) continue;
		for(int jchild=0; jchild<nchildren; jchild++){
			child = father->SubElement(jchild);
			int orient = childorient[jchild]*fatherorient;
			int current_child_polyh = GetPolyhedralIndex({child->Index(),orient});
			if(current_child_polyh < 0){
				SetPolyhedralIndex({child->Index(),orient},father_polyhindex);
			}else if(father_polyhindex != current_child_polyh){
				PZError << "\nFather element is trying to overwrite its child's polyh index\n"
						<<	"Father index " << father->Index()
						<<	"Father polyh " << father_polyhindex
						<<	"Child index " << child->Index()
						<<	"Child polyh " << current_child_polyh;
				DFN_DebugStop();
			}
		}
		SetPolyhedralIndex({father->Index(),fatherorient},-1); //< refined face doesn't need a polyhedral index
		DFNPolyhedron& polyhedron = Polyhedron(father_polyhindex);
		polyhedron.SwapForChildren(father);
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
	// out << "Number of DFNVolumes          = " << fVolumes.size() << "\n";
	out << "Tolerable distance            = " << fTolDist << "\n";
	out << "Tolerable angle            	  = " << fTolAngle << "\n";
	out << "Tolerable cos(angle)          = " << fTolAngle_cos << "\n";
	
	out << "\n\tFRACTURES INFO:\n";

	for(auto frac : fFractures){
		frac->Print(out);
	}

	PrintRolodexes(out);
	PrintPolyhedra(out);

}

void DFNMesh::PrintRolodexes(std::ostream & out) const{
	out << "\n\nROLODEXES:\n";
	if(fSortedFaces.size() == 0) {
		out << "\n\"No rolodexes initialized in this mesh\"\n"; 
		return;
	}
	for(TPZGeoEl* edge : fGMesh->ElementVec()){
		if(!edge) continue;
		if(edge->Dimension() != 1) continue;
		if(edge->HasSubElement()) continue;
		if(edge->Index() >= fSortedFaces.size()) break;
		fSortedFaces[edge->Index()].Print(out);
	}
}
void DFNMesh::PrintPolyhedra(std::ostream & out) const{
	out <<"\n\nPOLYHEDRA BY FACE:\n";
	if(fPolyhedra.size() == 0) {
		out << "\n\"No polyhedra initialized in this DFNMesh\"\n"; 
		return;
	}
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
			polyh->Print(out,false);
		}
	}
	out.flush();
}





void DFNMesh::DumpVTK(bool polygon_representation, bool clearmaterials, std::string filename){
	if(clearmaterials) ClearMaterials();
	std::cout<<" -Dumping VTK graphics\r"<<std::flush;
	for(auto frac : FractureList()){
		if(clearmaterials) frac->CleanUp(frac->MaterialId());
		if(polygon_representation) frac->Polygon().InsertGeomRepresentation(fGMesh,100+frac->Index(),1);
	}
	if(clearmaterials) 
		PrintVTKColorful(filename,"skip");
	else
		PrintVTK(filename,"skip");
	std::cout<<"                      \r"<<std::flush;
}




void DFNMesh::DFN_DebugStop(){
	fGMesh->BuildConnectivity();
	PlotAllPolygons("./LOG/AllDFNPolygons.vtk");
	PrintSummary();
	std::ofstream pzmesh("LOG/pzmesh.txt");
	fGMesh->Print(pzmesh);
	DumpVTK(true,true);
	DebugStop();
}



void DFNMesh::PreRefine(int n){
	if(n == 0) return;
	if(fFractures.size() > 0){
		PZError << "\n__PRETTY_FUNCTION__\nWas called on a mesh that might have already been refined for a fracture.\n";
		DebugStop();
	}

	gRefDBase.InitializeUniformRefPattern(MElementType::EOned);
	gRefDBase.InitializeUniformRefPattern(MElementType::EQuadrilateral);
	gRefDBase.InitializeUniformRefPattern(MElementType::ECube);

	gRefDBase.InitializeUniformRefPattern(MElementType::ETriangle);
	gRefDBase.InitializeUniformRefPattern(MElementType::ETetraedro);
	gRefDBase.InitializeUniformRefPattern(MElementType::EPrisma);
	gRefDBase.InitializeUniformRefPattern(MElementType::EPiramide);

	int maxdim = fGMesh->Dimension();
	TPZManVector<TPZGeoEl*,8> children;
	for(int i=0; i<n; i++){
		std::cout << "Pre-refining mesh x" << i+1 << std::endl;
		int64_t nels = fGMesh->NElements();
		for(int dim=1; dim<=maxdim; dim++){
			for(int64_t iel =0; iel < nels ; iel++){
				TPZGeoEl* gel = fGMesh->Element(iel);
				if(!gel) continue;
				if(gel->Dimension() != dim) continue;
				if(gel->HasSubElement()) continue;
				
				gel->Divide(children);
			}
		}
		// CreateSkeletonElements(maxdim-1);
		// CreateSkeletonElements(maxdim-2);
		// UpdatePolyhedra();
		// DFN::PlotJacobian(fGMesh,"LOG/jacplot." + std::to_string(i+1) + ".vtk");
	}

}

void CreateFilterScript(DFNMesh& dfn, std::ofstream& filter,std::string filename, const std::string ColorPreset = "Rainbow Uniform");
// void PlotAllPolygons(DFNMesh& dfn, std::string filename);


void DFNMesh::ExportDetailedGraphics(const std::string ColorPreset){
	TPZGeoMesh* gmesh = this->Mesh();
#ifdef PZDEBUG
	std::cout << "\n[WARNING] " << __PRETTY_FUNCTION__ << " inserts graphical elements that may leave the TPZGeoMesh inconsistent. It was meant to be called at the end of your script.";
#endif // PZDEBUG
	TPZVec<int> matid_backup;
	ClearMaterials(GMESHNOMATERIAL,matid_backup);
	
	// Standard control
	const std::string filename = "Fracture";
	const std::string dirname = "./graphics";
	std::filesystem::create_directory(dirname);
	

	// Plot each of the fractures
	PlotAllPolygons(dirname + "/allPolygons.vtk");
	for(auto frac : fFractures){
		std::string exportname = dirname + '/' + filename + "." + std::to_string(frac->Index()) + ".vtk";
		frac->PlotVTK(frac->Index(),exportname);
	}

	// Plot complete graphics
	// PrintVTK(dirname+'/'+"Complete"+'.'+std::to_string(NFractures())+".vtk","skip");

	// Create a filter script
	std::ofstream filter("graphics/GraphicsScript.py");
	CreateFilterScript(*this,filter,filename,ColorPreset);


	RestoreMaterials(matid_backup);
}

void CreateFilterScript(DFNMesh& dfn, std::ofstream& filter,std::string filename,const std::string ColorPreset){
	const std::string cwd = std::filesystem::current_path();

	// TUTORIAL
	filter  << "########################################################################\n"
			<< "#  _   _                 _                                              \n"
			<< "# | | | | _____      __ | |_ ___    _   _ ___  ___   _ __ ___   ___   _ \n"
			<< "# | |_| |/ _ \\ \\ /\\ / / | __/ _ \\  | | | / __|/ _ \\ | '_ ` _ \\ / _ \\ (_)\n"
			<< "# |  _  | (_) \\ V  V /  | || (_) | | |_| \\__ \\  __/ | | | | | |  __/  _ \n"
			<< "# |_| |_|\\___/ \\_/\\_/    \\__\\___/   \\__,_|___/\\___| |_| |_| |_|\\___| (_)\n"
			<< "#                                                                       \n"
			<< "# 1. Clear paraview (File > Disconnect > Yes);\n"
			<< "# 2. Open python shell (View > Python Shell);\n"
			<< "# 3. Select 'Run Script' in the lower right corner of the Python Shell;\n"
			<< "# 4. Search and execute this python script;\n"
			<< "########################################################################\n\n\n";

	std::printf("\n-Creating python filter script\n");
	filter <<   "#### import the simple module from the paraview\n"
				"from paraview.simple import *\n"
				"#### disable automatic camera reset on 'Show'\n"
				"paraview.simple._DisableFirstRenderCameraReset()\n";
	
	filter << "\n# create a new \'Legacy VTK Reader\' for CoarseMesh file\n"
			<< "coarseMeshvtk = LegacyVTKReader(FileNames=[\'" << cwd << "/graphics/CoarseMesh.vtk'])\n"
			<< "RenameSource(\'CoarseMesh\', coarseMeshvtk)\n";
	// Setup layout and view to create filters
	filter << "\n\n";
	filter << 	"# get active view\n"
				"renderView1 = GetActiveViewOrCreate('RenderView')\n"
				"# get layout\n"
				"layout1 = GetLayout()\n";

	
	filter <<	"\n# find source of Coarse Mesh\n"
				"CoarseMesh = FindSource(\'"<< "CoarseMesh.vtk" << "\')\n";
	
	filter  << "CoarseDisplay = Show("<<"CoarseMesh"<<", renderView1,'UnstructuredGridRepresentation')\n"
			<< "CoarseDisplay.SetRepresentationType('Wireframe')\n"
			<< "CoarseDisplay.AmbientColor = [0.443137255, 0.450980392, 0.6]\n"
			<< "CoarseDisplay.DiffuseColor = [0.443137255, 0.450980392, 0.6]\n";
	

	filter << "\n# create a new \'Legacy VTK Reader\' for DFNPolygons file\n"
			<< "DFNPolygonsvtk = LegacyVTKReader(FileNames=[\'" << cwd << "/graphics/allPolygons.vtk'])\n"
			<< "RenameSource(\'AllDFNPolygons\', DFNPolygonsvtk)\n";
	// Setup layout and view to create filters
	filter << "\n\n";
	filter << 	"# get active view\n"
				"renderView1 = GetActiveViewOrCreate('RenderView')\n"
				"# get layout\n"
				"layout1 = GetLayout()\n";

	
	filter <<	"\n# find source of allPolygons\n"
				"DFNPolygons = FindSource(\'"<< "AllDFNPolygons" << "\')\n";
	
	filter  << "allPolygonsDisplay = Show("<<"DFNPolygons"<<", renderView1,'UnstructuredGridRepresentation')\n"
			<< "allPolygonsDisplay.SetRepresentationType('Wireframe')\n"
			<< "allPolygonsDisplay.AmbientColor = [0.0, 0.0, 0.0]\n"
			<< "allPolygonsDisplay.DiffuseColor = [0.0, 0.0, 0.0]\n";
	






	filter << "# create a new \'Legacy VTK Reader\' for Fracture list\n"
			<< "fracturevtk = LegacyVTKReader(";
	
	std::stringstream fracturelist;
	// std::stringstream timesteps;
	fracturelist << "FileNames=[";
	// timesteps << "TimestepValues=[";
	for(auto frac : dfn.FractureList()){
		fracturelist << "\n\t\'" << cwd << "/graphics/Fracture." << frac->Index() << ".vtk\',";
		// timesteps << ' ' << frac->Index() << ',';
	}
	fracturelist.seekp(fracturelist.str().length()-1); ///< removes spare comma
	fracturelist << ']';
	// timesteps.seekp(timesteps.str().length()-1); ///< removes spare comma
	// timesteps << ']';

	filter << fracturelist.str();
	// filter << ",\n" << timesteps.str();
	filter << "\n)\n";
	filter << "RenameSource(\'Fracture\', fracturevtk)\n";

	filter <<	"\n# find source of fractures\n"
				"Fracture = FindSource(\'"<< filename << "\')\n";
	// Update current layout and view to create filters
	filter << "\n\n";
	filter << 	"# get active view\n"
				"renderView1 = GetActiveViewOrCreate('RenderView')\n"
				"# get layout\n"
				"layout1 = GetLayout()\n";

	// FILTERS
	// @{

	int nfrac = dfn.NFractures();

	filter << "\n\n";
	std::string rootfilter = "Fracture";
	std::string frac = "Fracture";
	std::string Polygon = "Polygon";
	std::string Surface = "Surface";
	std::string Boundary = "Boundary";
	std::string Intersection = "Intersections";
	std::string Shrink = "Shrink";
	std::string two_dim = "two_dim";
	std::string one_dim = "one_dim";
	filter  << "\n# (2D filter)" << "\n"
			<< two_dim <<" = Threshold(Input=" <<rootfilter<<")\n"
			<< two_dim <<".Scalars = ['CELLS', 'Dimension']\n"
			<< two_dim <<".ThresholdRange = [2, 2]\n"
			<< "RenameSource('2D', "<<two_dim<<")\n";
	
	filter  << "\n# Shrink (" <<two_dim<< ")\n"
			<< Shrink <<" = Shrink(Input=" <<two_dim<<")\n"
			<< Shrink <<".ShrinkFactor = 0.9\n"
			<< "RenameSource('Shrink', "<<Shrink<<")\n";


	filter  << "\n\n# (" << Surface << ")\n"
			<< Surface << " = Threshold(Input="<<Shrink<<")\n"
			<< Surface << ".Scalars = ['CELLS', 'material']\n"
			<< Surface << ".ThresholdRange = [0, "<<nfrac-1<<"]\n"
			<< "RenameSource('" << Surface << "', " << Surface << ")\n"
			<< "# show data in view\n"
			<< Surface << "Display = Show(" << Surface << ", renderView1,'UnstructuredGridRepresentation')\n"
			<< "# set representation type\n"
			<< Surface << "Display.SetRepresentationType('Surface')\n";
	filter  << "# set scalar coloring using an separate color/opacity maps\n"
			<< "ColorBy(" << Surface << "Display, ('CELLS', 'material'), separate=True)\n"
			<< "# Hide colormap bar\n"
			<< Surface << "Display.SetScalarBarVisibility(renderView1, False)\n"
			<< "# rescale color and/or opacity maps used to include current data range\n"
			<< Surface << "Display.RescaleTransferFunctionToDataRange(True, False)\n"
			<< "# get separate color transfer function/color map for 'material'\n"
			<< Surface << "Display_materialLUT = GetColorTransferFunction('material', " << Surface << "Display, separate=True)\n"
			<< "# get separate opacity transfer function/opacity map for 'material'\n"
			<< Surface << "Display_materialPWF = GetOpacityTransferFunction('material', " << Surface << "Display, separate=True)\n"
			<< "# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.\n"
			<< Surface << "Display_materialLUT.ApplyPreset('"<<ColorPreset<<"', True)\n"
			<< "# Rescale transfer function\n"
			<< Surface << "Display_materialLUT.RescaleTransferFunction(0.0, " << nfrac-1 << ".0)\n"
			<< "# Rescale transfer function\n"
			<< Surface << "Display_materialPWF.RescaleTransferFunction(0.0, " << nfrac-1 << ".0)\n\n";
	


	filter  << "\n# " << two_dim <<" (Polygon)" << "\n"
			<< Polygon <<" = Threshold(Input=" <<two_dim<<")\n"
			<< Polygon <<".Scalars = ['CELLS', 'material']\n"
			<< Polygon <<".ThresholdRange = [" << -(nfrac+1) <<','<< -1 <<"]\n"
			<< "RenameSource('Polygon', "<<Polygon<<")\n";
	filter  << Polygon <<"Display = Show("<<Polygon<<", renderView1,'UnstructuredGridRepresentation')\n"
			<< Polygon<<"Display.SetRepresentationType('Wireframe')\n"
			<< Polygon<<"Display.AmbientColor = [0.0, 0.0, 0.0]\n"
			<< Polygon<<"Display.DiffuseColor = [0.0, 0.0, 0.0]\n";
			// << "ColorBy("<<Polygon<<"Display, ('CELLS', 'material'))\n"
			// << Polygon<<"Display.RescaleTransferFunctionToDataRange(True, False)\n"
			// << Polygon<<"Display.SetScalarBarVisibility(renderView1, False)\n";
	
	filter  << "\n# (1D filter)" << "\n"
			<< one_dim <<" = Threshold(Input=" <<rootfilter<<")\n"
			<< one_dim <<".Scalars = ['CELLS', 'Dimension']\n"
			<< one_dim <<".ThresholdRange = [1, 1]\n"
			<< "RenameSource('1D', "<<one_dim<<")\n";
	
	


	filter  << "\n\n# (" << Boundary << ")\n"
			<< Boundary << " = Threshold(Input="<<one_dim<<")\n"
			<< Boundary << ".Scalars = ['CELLS', 'material']\n"
			<< Boundary << ".ThresholdRange = ["<<nfrac<<", "<<2*nfrac-1<<"]\n"
			<< "RenameSource('" << Boundary << "', " << Boundary << ")\n"
			<< "# show data in view\n"
			<< Boundary << "Display = Show(" << Boundary << ", renderView1,'UnstructuredGridRepresentation')\n"
			<< "# set representation type\n"
			<< Boundary << "Display.SetRepresentationType('Surface With Edges')\n";
	filter  << "# set scalar coloring using an separate color/opacity maps\n"
			<< "ColorBy(" << Boundary << "Display, ('CELLS', 'material'), separate=True)\n"
			<< "# Hide colormap bar\n"
			<< Boundary << "Display.SetScalarBarVisibility(renderView1, False)\n"
			// << Boundary<<"Display.AmbientColor = [1.0, 0.0, 0.0]\n"
			// << Boundary<<"Display.DiffuseColor = [1.0, 0.0, 0.0]\n"
			<< "# Set Line Width\n"
			<< Boundary<<"Display.LineWidth = 4.0\n"
			<< "# rescale color and/or opacity maps used to include current data range\n"
			<< Boundary << "Display.RescaleTransferFunctionToDataRange(True, False)\n"
			<< "# get separate color transfer function/color map for 'material'\n"
			<< Boundary << "Display_materialLUT = GetColorTransferFunction('material', " << Boundary << "Display, separate=True)\n"
			<< "# get separate opacity transfer function/opacity map for 'material'\n"
			<< Boundary << "Display_materialPWF = GetOpacityTransferFunction('material', " << Boundary << "Display, separate=True)\n"
			<< "# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.\n"
			<< Boundary << "Display_materialLUT.ApplyPreset('"<<ColorPreset<<"', True)\n"
			<< "# Rescale transfer function\n"
			<< Boundary << "Display_materialLUT.RescaleTransferFunction("<<nfrac<<".0, "<<2*nfrac-1<<".0)\n"
			<< "# Rescale transfer function\n"
			<< Boundary << "Display_materialPWF.RescaleTransferFunction("<<nfrac<<".0, "<<2*nfrac-1<<".0)\n\n";


	int TriangleNumber = nfrac*(nfrac+1)/2; ///< Triangle Number bounds the range for frac-to-frac intersections (wikipedia explains it better than I, but you can see this as the answer for the question "how many unique entries in a square symmetric matrix")
	int max_intersections = TriangleNumber;

	filter  << "\n\n# (" << Intersection << ")\n"
			<< Intersection << " = Threshold(Input="<<one_dim<<")\n"
			<< Intersection << ".Scalars = ['CELLS', 'material']\n"
			<< Intersection << ".ThresholdRange = [" << 2*nfrac<<','<<2*nfrac + max_intersections << "]\n"
			<< "RenameSource('" << Intersection << "', " << Intersection << ")\n"
			<< "# show data in view\n"
			<< Intersection << "Display = Show(" << Intersection << ", renderView1,'UnstructuredGridRepresentation')\n"
			<< "# set representation type\n"
			<< Intersection << "Display.SetRepresentationType('Surface With Edges')\n";
	filter  << "# set scalar coloring using an separate color/opacity maps\n"
			<< "ColorBy(" << Intersection << "Display, ('CELLS', 'material'), separate=True)\n"
			<< "# Hide colormap bar\n"
			<< Intersection << "Display.SetScalarBarVisibility(renderView1, False)\n"
			// << Intersection <<"Display.AmbientColor = [0.0, 1.0, 0.0]\n"
			// << Intersection <<"Display.DiffuseColor = [0.0, 1.0, 0.0]\n"
			<< "# Set Line Width\n"
			<< Intersection<<"Display.LineWidth = 4.0\n"
			<< "# rescale color and/or opacity maps used to include current data range\n"
			<< Intersection << "Display.RescaleTransferFunctionToDataRange(True, False)\n"
			<< "# get separate color transfer function/color map for 'material'\n"
			<< Intersection << "Display_materialLUT = GetColorTransferFunction('material', " << Intersection << "Display, separate=True)\n"
			<< "# get separate opacity transfer function/opacity map for 'material'\n"
			<< Intersection << "Display_materialPWF = GetOpacityTransferFunction('material', " << Intersection << "Display, separate=True)\n"
			<< "# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.\n"
			<< Intersection << "Display_materialLUT.ApplyPreset('"<<ColorPreset<<"', True)\n"
			<< "# Rescale transfer function\n"
			<< Intersection << "Display_materialLUT.RescaleTransferFunction(" << 2*nfrac<<','<<2*nfrac + max_intersections << ".0)\n"
			<< "# Rescale transfer function\n"
			<< Intersection << "Display_materialPWF.RescaleTransferFunction(" << 2*nfrac<<','<<2*nfrac + max_intersections << ".0)\n\n";
	// }@ FILTERS


	// TIMESTEP CONFIGURATION
	filter  << "\n# get animation scene\n"
			<< "animationScene1 = GetAnimationScene()\n"
			<< "# get the time-keeper\n"
			<< "timeKeeper1 = GetTimeKeeper()\n"
			<< "# Properties modified on animationScene1\n"
			<< "animationScene1.PlayMode = 'Snap To TimeSteps'\n";


	// {@ FINAL SETUP
	filter  << '\n' << '\n';
	filter  << "Hide(" << Polygon << ", renderView1)\n";

	filter  << "\n\nrenderView1.Update()\n"
			<< "# current camera placement for renderView1\n"
			<< "renderView1.InteractionMode = '3D'\n"
			<< "renderView1.CameraPosition = [3.47, 4.0, 10000.0]\n"
			<< "renderView1.CameraFocalPoint = [3.47, 4.0, 0.0]\n"
			<< "renderView1.CameraParallelScale = 7.1307240904601645\n"
			<< "\n# reset view to fit data\n"
			<< "renderView1.ResetCamera()\n";
	// }@ FINAL SETUP

	filter.flush();
}

void DFNMesh::PrintProblematicRolodex(const int &indexNotFoundCard, TRolodex &rol) {
    
    std::cout << "\n =====> Printing problematic rolodex to ProblemRolodex.vtk" << std::endl;
    
    if (indexNotFoundCard < 0) {
        DebugStop();
    }
    
    TPZGeoMesh bogusMesh;
    bogusMesh.ElementVec().Resize(Mesh()->NElements());
    for (int i = 0; i < bogusMesh.NElements(); i++) {
        bogusMesh.ElementVec()[i] = nullptr;
    }
    bogusMesh.NodeVec() = Mesh()->NodeVec();
    
    // Not found card
    TPZGeoEl* notfoundcard = Mesh()->Element(indexNotFoundCard);
    if (!notfoundcard) {
        DebugStop();
    }
    TPZGeoEl* copiedNotFoundCard = notfoundcard->Clone(bogusMesh);
    copiedNotFoundCard->SetMaterialId(3);

    // Axle
    TPZGeoEl* axle = Mesh()->Element(rol.fedgeindex);
    if (!axle) {
        DebugStop();
    }
    TPZGeoEl* newaxle = axle->Clone(bogusMesh);
    newaxle->SetMaterialId(2);
    
    // All the preexisting cards in rolodex
    for (auto& card : rol.fcards) {
        const int index = card.fgelindex;
        if (index < 0) {
            DebugStop();
        }
        if (index == indexNotFoundCard) {
            std::cout << "\n====> Card is in rolodex! Everything should be fine!" << std::endl;
        }
        TPZGeoEl *gelInRolodex = Mesh()->Element(index);
        if (!gelInRolodex) {
            DebugStop();
        }
        TPZGeoEl* newel = gelInRolodex->Clone(bogusMesh);
    }
    
    std::ofstream out("ProblemRolodex.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&bogusMesh, out, true, true);
}

void DFNMesh::RollBackLastFracture(TPZGeoMesh *gmeshBackup, TPZStack<int>& badVolumes) {
    LOGPZ_DEBUG(logger,"\nBad volumes found during cut of fracture"
                << NFractures() << "\nThe volume indexes are = "
                << badVolumes);
    for (int i = 0; i < fFractures.size()-1; i++) {
        fFractures[i]->RollBack(gmeshBackup);
    }
    
    delete fGMesh;
    
    DFNFracture *frac = fFractures[fFractures.size()-1];
    delete frac;
    fFractures[fFractures.size()-1] = nullptr;
    fFractures.resize(fFractures.size()-1);
    
    fGMesh = gmeshBackup;
    
    // refine polyhedra
    for (auto& polindex : badVolumes){
        DFNPolyhedron& volume = Polyhedron(polindex);
        volume.Refine();
    }
}



void DFNMesh::PrintRolodexBugReport(const int64_t AxleIndex){
	// Consistency
	if(AxleIndex < 0 || AxleIndex > fGMesh->NElements()) throw std::invalid_argument("AxleIndex does not point to an edge element in the mesh.\n");
	TPZGeoEl* edge = fGMesh->Element(AxleIndex);
	if(edge->Dimension() != 1) throw std::invalid_argument("AxleIndex does not point to an edge element in the mesh.\n");

	const std::string dirname = "./LOG/RolodexBugReport";
	std::filesystem::create_directories(dirname);

	// Polygons are usually of interest when debugging
	PlotAllPolygons(dirname+"/AllDFNPolygons.vtk");
	// You could alternativelly print them as timesteps to ParaView
	// for(auto frac : fFractures) frac->Polygon().PlotVTK(dirname+"DFNPolygon."+std::to_string(frac->Index())+".vtk");

	// Plot current Rolodex (if existent)
	if(AxleIndex < fSortedFaces.size() && fSortedFaces[AxleIndex].NCards() != 0){
		TRolodex& rolodex = fSortedFaces[AxleIndex];
		rolodex.PlotVTK(dirname+"/CurrentRolodex.vtk",fGMesh);
	}

	// Plot Neighbours of edge to see if they match the Rolodex
	TPZGeoElSide edgeside = {edge,2};
	constexpr int twoD = 2, UnrefinedOnly = true, orientationMatch = true;
	DFN::PlotNeighbours(dirname+"/twoDNeighbours.vtk",edgeside,twoD,UnrefinedOnly,orientationMatch);
	DFN::PlotNeighbours(dirname+"/threeDNeighbours.vtk",edgeside,3,UnrefinedOnly,orientationMatch);

	// Plot polyhedral volumes
	std::set<int> toPlot;
	TPZGeoElSide neig = edgeside.Neighbour();
	for(/*void*/; neig!= edgeside; ++neig){
		if(neig.Element()->Dimension() != 2) continue;
		if(neig.Element()->HasSubElement()) continue;

		int polyhindex = GetPolyhedralIndex(neig.Element()->Index(),+1);
		if(!(polyhindex < 0)) toPlot.insert(polyhindex);
		polyhindex = GetPolyhedralIndex(neig.Element()->Index(),-1);
		if(!(polyhindex < 0)) toPlot.insert(polyhindex);
	}
	for(int polyhindex : toPlot){
		if(polyhindex >= NPolyhedra()) continue;
		DFNPolyhedron& polyh = Polyhedron(polyhindex);
		polyh.PlotVTK(dirname+"/Poly_"+std::to_string(polyhindex)+".vtk");
	}
}

void DFNMesh::PlotVolumesByCoarseIndex(const int64_t coarseindex, const std::string dirpath) const {
    
    for (auto& pol : fPolyhedra) {
        if(pol->CoarseIndex() != coarseindex) continue;
        std::string filename = dirpath + "poly_" + std::to_string(pol->Index()) + "." + std::to_string(NFractures()) + ".vtk";
        pol->PlotVTK(filename);
        
    }
}


void DFNMesh::PlotAllPolygons(const std::string filepath) const{
    constexpr int type = 7; // Arbitrary polygon
    constexpr int eldimension = 2;

	std::ofstream file(filepath);
    std::stringstream node, connectivity, materialstream, indexstream, dimensionstream, typestream;

    //Header
    file << "# vtk DataFile Version 3.0\n";
    file << "DFNPolygon VTK Visualization\n";
    file << "ASCII\n\n";

    file << "DATASET UNSTRUCTURED_GRID\n";

	int totalNNodes = 0;

	
	for(DFNFracture* fracture : fFractures){
		const DFNPolygon& polygon = fracture->Polygon();

		const int elNnodes = polygon.NCornerNodes();
		connectivity << elNnodes;
		TPZManVector<REAL,3> cornerX(3,0.0);
		for(int64_t inode = 0; inode < elNnodes; inode++) {
			polygon.iCornerX(inode,cornerX);
			for (int c = 0; c < 3; c++) {
				REAL coord = cornerX[c];
				node << coord << " ";
			}
			node << '\n';
			connectivity << ' ' << totalNNodes + inode;
		}
    	connectivity << '\n';

		materialstream << fracture->MaterialId() << '\n';
		indexstream << fracture->Index() << '\n';
		// dimensionstream << eldimension << '\n';
		typestream << type << '\n';
		totalNNodes += elNnodes;
	}
    node << '\n';
	connectivity << '\n';

    file << "POINTS " << totalNNodes << " float\n";
    file << node.str();
    int64_t size = totalNNodes+NFractures();
    file << "CELLS " << NFractures() << ' ' << size << '\n';
    file << connectivity.str() << '\n';


    file << "CELL_TYPES " << NFractures() << '\n' << typestream.str() << '\n';
    file << "CELL_DATA "<< NFractures() <<'\n';
    file << "FIELD FieldData 2\n";
    file << "material 1 " << NFractures() << " int\n" << materialstream.str();
    file << "elIndex 1 " << NFractures() << " int\n" << indexstream.str();
    // file << "Dimension 1 " << NFractures() << " int\n" << dimensionstream.str();
    file.close();
}
