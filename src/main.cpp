
//includes
	#ifdef HAVE_CONFIG_H
	#include <pz_config.h>
	#endif
	#include "pzgmesh.h"
	#include "pzgengrid.h"
	#include "TPZExtendGridDimension.h"
	#include "TPZVTKGeoMesh.h"

	// #include "pzlog.h"

	#include <stdio.h>
	#include <math.h>
	#include <iostream>
	#include <fstream>
	#include <string>
	#include <sstream>
	#include <cstdio>

	#include <set>
	#include <unordered_set>
	#include <queue>
	#include <map>
	#include <vector>

	#include "DFNFractureMesh.h"
	#include "DFNRibs.h"
	#include "DFNFace.h"
	#include "DFNVolume.h"

	#include "TPZRefPatternDataBase.h"

	#include <gmsh.h>

	#define fTolerance 10e-5
//includes


void ReadFractureFromFile(std::string filename, TPZFMatrix<REAL> &plane);
TPZGeoMesh* ReadExampleFromFile(std::string filename, TPZManVector<TPZFMatrix<REAL> > &planevector);



struct DFNMesh{
	// private:
		std::list<DFNFractureMesh *> fFractures;
		std::map<int64_t, DFNVolume> fVolumes;
	// public:
		/// Pointer to volume of index 'index'
		DFNVolume *Volume(int64_t index){return &fVolumes[index];}
		/// Uses GMsh to mesh volumes cut by fracture plane
		void CreateVolumes();
		/// Uses gmsh API to tetrahedralize volumes
    	void Tetrahedralize(DFNVolume *volume);
		/// Uses gmsh API to tetrahedralize a DFNVolume
		void Tetrahedralize2(DFNVolume *volume);
    	/// Find the volumetrical element that encloses a 2D element
    	bool FindEnclosingVolume(TPZGeoEl *ifracface);
	// private:
		/**
    	 *  @brief Navigate children tree to access most extreme branches
    	 *  @param gel: Pointer to geometric element of eldest ancestor
    	 *  @param outfile: ofstream in which to write accessed data
    	 */
    	void PrintYoungestChildren(TPZGeoEl *gel, std::ofstream &outfile);
		
		/**
		 * @brief Navigate through neighbours of first level, than second level, and so on, until an element of a specific material id is found
		 * @returns Index of eldest ancestor of such element
		*/
		int64_t SearchIndirectNeighbours(TPZGeoEl* gel);
		/**
		 * @brief Goes through all neighbours of gel and identifies if any of them has material id different of that from surface elements
		 * @returns Index of eldest ancestor of first found neighbour that matches the material id criteria (macro element)
		 * @note Later I might modify this to fill a vector/list with all neighbours that match the material id criteria
		*/
		int64_t FindAdjacentMacroEl(TPZGeoEl* gel);
		/**
		 * @brief Pushes all neighbours of a geometric element onto the back of a list
		 * @note Not all neighbours are pushed to the list, but rather some criteria is specified so that it gets only those that are candidates of higher interest. Currently these would be 2D elements that are neighbours through the edges of gel.
		 * @param gel: Pointer to geometric element
		 * @param candidate_queue: Reference to current list of candidates 
		*/
		void QueueNeighbours(TPZGeoEl* gel,   std::list<int64_t> &candidate_queue);
};















//MATERIAL ID MAP
// 1 gmesh (default)
// 4 skeleton
// 12 ribs that will be divided
// 18 children ribs
// 20 mid-fracture cut faces
// 35 end-fracture cut faces
// 40 Fracture plane
// 45 Intersection points in end-faces

using namespace std;

int main(){
	gRefDBase.InitializeUniformRefPattern(EOned);
	TPZManVector< TPZFMatrix<REAL>> planevector;
	TPZGeoMesh *gmesh;
	gmesh = ReadExampleFromFile("examples/exampleOctagon.txt",planevector);

	int surfaceMaterial = 40;
	int transitionMaterial = 20;
	DFNMesh dfn;
	for(int iplane = 0,
		    nfractures = planevector.size(); iplane < nfractures; iplane++){
		DFNFracPlane *fracplane = new DFNFracPlane(planevector[iplane]);
		DFNFractureMesh *fracmesh = new DFNFractureMesh(*fracplane,gmesh,surfaceMaterial);
		// Find and split intersected ribs
		fracmesh->SplitRibs(transitionMaterial);
		// Find and split intersected faces
		fracmesh->SplitFaces(transitionMaterial);
		std::cout<<"\n";
		// Mesh fracture surface
		fracmesh->SplitFracturePlane();
		dfn.fFractures.push_back(fracmesh);
		//Print result
		std::ofstream meshprint("meshprint.txt");
		std::ofstream out1("./TestSurfaces.vtk");
		gmesh->Print(meshprint);
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out1, true);
	}
	// Mesh transition volumes
	dfn.CreateVolumes();


	return 0;
}














void DFNMesh::CreateVolumes(){
    TPZGeoMesh *gmesh = (*fFractures.begin())->GetGeoMesh();
	

    // map all volumes that are cut
    int64_t nels = gmesh->NElements();
	for (int64_t iel = 0; iel < nels; iel++){
        TPZGeoEl *gel = gmesh->Element(iel);
        if(gel->Dimension() != 3){continue;}
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
    
    
	// search through each 2D element of the triangulated fractures surfaces to find their enclosing volume
	int surfaceMaterial = (*fFractures.begin())->GetSurfaceMaterial();
	for(int64_t iel = 0; iel < nels; iel++){
		TPZGeoEl *gel = gmesh->Element(iel);
		// if(gel->MaterialID() != surfaceMaterial){continue;}
		// During development, elements at fracture surface have material id over 40
		if(gel->MaterialId() <= surfaceMaterial) continue;
		if(gel->Dimension() != 2) continue;
		if(gel->HasSubElement()) continue;
		// Find volume that encloses that element
		FindEnclosingVolume(gel);
	}
	
	// debug & verification__________________________________________________________________________
	// for(auto ivolume = fVolumes.begin(); ivolume != fVolumes.end(); ivolume++){
	// 	TPZManVector<int64_t> faces = ivolume->second.GetFacesInVolume();
	// 	for(int i = 0; i<faces.size(); i++){
	// 		std::cout<<"vol # "<<ivolume->first<<" \t face # "<<faces[i]<<"\n";
	// 	}
	// }
	// debug & verification__________________________________________________________________________
	
	
	Tetrahedralize(&fVolumes.begin()->second);
	
	// //Loop over list of volumes cut
	// for (auto itr = fVolumes.begin(); itr != fVolumes.end(); itr++){
    // 	DFNVolume *ivolume = &itr->second;
	// 	// Use GMsh to tetrahedralize volumes
    // 	Tetrahedralize(ivolume);
	// }
	

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
}







int64_t DFNMesh::FindAdjacentMacroEl(TPZGeoEl* gel){
	int surfaceMaterial = (*fFractures.begin())->GetSurfaceMaterial();
    int nsides = gel->NSides();
        
	for(int iside = nsides-2; iside >= 0; iside--){
		TPZGeoElSide gelside(gel, iside);
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
	if(ifracface->Dimension()!=2) DebugStop();
	int surfaceMaterial = (*fFractures.begin())->GetSurfaceMaterial();
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
        fVolumes[volumeindex].SetFaceInVolume(ifracface->Index());
		// std::cout<<"Face #"<<ifracface->Index()<<" \tin volume #"<<volumeindex<<"\n";
        return true;
    }
    std::cout<<"\n DFNFractureMesh::FindEnclosingVolume found no enclosing volume for element #"<<ifracface->Index()<<"\n";
    DebugStop();
    return false;
}










/**
 * 	@brief Uses GMsh API to tetrahedralize a DFNVolume
 */ 
void DFNMesh::Tetrahedralize2(DFNVolume *volume){
	// GMsh does not accept zero index entities
    const int shift = 1;

	TPZGeoMesh *gmesh = (*fFractures.begin())->GetGeoMesh();
	int64_t ivol = volume->ElementIndex();
	TPZGeoEl *volGel = gmesh->Element(ivol);
	TPZManVector<int64_t> enclosedFaces = volume->GetFacesInVolume();

	// Initialize GMsh
	gmsh::initialize();
		gmsh::option::setNumber("Mesh.Algorithm3D",1); // (1: Delaunay, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT) Default value: 1
		// Insert nodes ______________________________
		{
			std::set<int64_t> nodes;
			// corners of volume
			for(int icorner = 0,
					ncorners = volGel->NCornerNodes(); icorner < ncorners; icorner++){
				int64_t inode = volGel->NodeIndex(icorner);
				nodes.insert(inode);
			}
			// nodes of skeleton children faces
			for(int nsides = volGel->NSides(),
					iside = nsides-2; iside >= 0; iside--){
				TPZGeoElSide gelside(volGel,iside);
				if(gelside.Dimension() != 2) break;
				TPZStack<TPZGeoEl *> unrefinedSons;
				gelside.Element()->GetAllSiblings(unrefinedSons);
				for(TPZGeoEl *child : unrefinedSons){
					for(int icorner = 0,
					ncorners = child->NCornerNodes(); icorner < ncorners; icorner++){
						int64_t inode = child->NodeIndex(icorner);
						nodes.insert(inode);
					}
				}
			}
			// nodes of enclosed faces
			for(int iface = 0,
					nfaces = enclosedFaces.size(); iface < nfaces; iface++){
				TPZGeoEl *gel = gmesh->Element(enclosedFaces[iface]);
				for(int icorner = 0,
				ncorners = gel->NCornerNodes(); icorner < ncorners; icorner++){
					int64_t inode = gel->NodeIndex(icorner);
					nodes.insert(inode);
				}
			}
			// Pass nodes to GMsh
			for(int64_t inode : nodes){
				TPZManVector<REAL,3> coord(3);
				gmesh->NodeVec()[inode].GetCoordinates(coord);
				gmsh::model::geo::addPoint(coord[0],coord[1],coord[2],0.,inode+shift);
			}
		}


		// Insert Lines ______________________________
		{
			std::set<int64_t> lines;			
			// edges of skeleton children faces
			for(int nsides = volGel->NSides(),
					iside = nsides-2; iside >= 0; iside--){
				TPZStack<TPZGeoEl *> unrefinedSons;{
					TPZGeoElSide vol_gelside(volGel,iside);
					if(vol_gelside.Dimension() != 2) break;
					vol_gelside.Element()->GetAllSiblings(unrefinedSons);
				}
				for(TPZGeoEl *child : unrefinedSons){
					for(int iside_child = child->NCornerNodes(),
							nsides_child = child->NSides(); iside_child < nsides_child-1; iside_child++){
						TPZGeoElSide gelside(child,iside_child);
						TPZGeoElSide neig = gelside.Neighbour();
						while(neig.Dimension() != 1){
							neig = neig.Neighbour();
						}
						int64_t iline = neig.Element()->Index();
						lines.insert(iline);
					}
				}
			}
			// edges of enclosed faces
			for(int iface = 0,
					nfaces = enclosedFaces.size(); iface < nfaces; iface++){
				TPZGeoEl *gel = gmesh->Element(enclosedFaces[iface]);
				for(int iside = gel->NCornerNodes(),
						nsides = gel->NSides(); iside < nsides-1; iside++){
					TPZGeoElSide gelside(gel,iside);
					TPZGeoElSide neig = gelside.Neighbour();
					while(neig.Dimension() != 1){
						neig = neig.Neighbour();
					}
					int64_t iline = neig.Element()->Index();
					lines.insert(iline);
				}
			}
			// Pass lines to GMsh
			for(int64_t iline : lines){
				TPZGeoEl *gel = gmesh->Element(iline);
				int64_t node0 = gel->NodeIndex(0)+shift;
				int64_t node1 = gel->NodeIndex(1)+shift;
				gmsh::model::geo::addLine(node0,node1,iline+shift);
				gmsh::model::geo::mesh::setTransfiniteCurve(iline+shift,2);
			}
		}
		// Insert Faces ______________________________
		{
			std::set<int64_t> faces;
			// skeleton children faces
			// enclosed faces
		}
		// Insert Volume _____________________________
		{
			// Surface loop of skeleton youngest children
		}
		// Surfaces in Volume ________________________
		{
			// Construct an std::vector
			gmsh::model::geo::synchronize();
			// gmsh::model::mesh::embed(int dim, const std::vector<int> &tags, int inDim, int inTag);
		}
		// Physical groups ____________________________
		{

		}
		// synchronize before meshing
			gmsh::model::geo::synchronize();
		// mesh
			gmsh::model::mesh::generate(3);
			gmsh::model::mesh::optimize("Netgen");
		// write (for testing)
			// gmsh::write("testAPI.msh");
		// import meshed plane back into PZ geoMesh
			// ImportElementsFromGMSH(gmesh,3);
	// close GMsh
	gmsh::finalize();
}











/**
 * 	@brief Uses GMsh to tetrahedralize a DFNVolume
 */ 
void DFNMesh::Tetrahedralize(DFNVolume *volume){
    int mtransition = 20;
    int msurface = 40;
    int mintact = 1;
	std::ofstream outfile("fracture1.geo");

    TPZGeoMesh *pzgmesh = (*fFractures.begin())->GetGeoMesh();
	int64_t ivol = volume->ElementIndex();

    // Giving fGMesh another name for readability's sake
    pzgmesh->BuildConnectivity();
    (*fFractures.begin())->CreateSkeletonElements(1,mtransition);
    // Title
    outfile<<"//  Geo file generated by Discrete Fracture Network methods \n"
            <<"// Fracture #1 \n\n";
    
    // write nodes
    outfile<< "// POINTS DEFINITION \n\n";
    int64_t nnodes = pzgmesh->NNodes();
    // @ToDo Do we need physical groups for points too?
    for (int64_t inode = 0; inode < nnodes; inode++){
        TPZManVector<REAL, 3> co(3,0.);
        pzgmesh->NodeVec()[inode].GetCoordinates(co);
        outfile << "Point(" << inode << ") = {" << co[0] << ',' << co[1] << ',' << co[2] << "};\n";
    }
    
    // write edges
    int64_t nels = pzgmesh->NElements();
    outfile << "\n\n// LINES DEFINITION \n\n";
    {
        // declare lists to define physical groups
        std::list<int64_t> groupSurface;
        std::list<int64_t> groupTransition;
        std::list<int64_t> groupIntact;
        // iterate over all 1D elements
        for (int64_t iel = 0; iel < nels; iel++){
            TPZGeoEl *gel = pzgmesh->Element(iel);
            if(gel->Dimension() != 1) continue;
            if(gel->HasSubElement()) continue;
            // if(gel->MaterialId() == fTransitionMaterial) continue;
            outfile << "Line(" << iel << ") = {" << gel->NodeIndex(0) << ',' << gel->NodeIndex(1) << "};\n";
    // @ToDo this is kind of a mess, but only for debugging
            // list it according to material
            if(gel->MaterialId() >= mtransition && gel->MaterialId() < msurface){groupTransition.push_back(iel);}
            else if(gel->MaterialId() >= msurface){groupSurface.push_back(iel);}
            else groupIntact.push_back(iel);
        }
        // write physical groups
        outfile<<"\nPhysical Curve("<<mtransition<<") = {";
        for(auto itr = groupTransition.begin(); itr != groupTransition.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupTransition.end()? "," : "};\n");
        }
        outfile<<"\nPhysical Curve("<<msurface<<") = {";
        for(auto itr = groupSurface.begin(); itr != groupSurface.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupSurface.end()? "," : "};\n");
        }
        outfile<<"\nPhysical Curve("<<mintact<<") = {";
        for(auto itr = groupIntact.begin(); itr != groupIntact.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupIntact.end()? "," : "};\n");
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
            if(gel->Dimension() != 2) continue;
            if(gel->HasSubElement()) continue;
            // if(gel->MaterialId() == fTransitionMaterial) continue;
            
            int nnodes = gel->NCornerNodes();
            int nedges = nnodes; //for readability 
            TPZManVector<int64_t,4> facenodevec(nnodes);
            gel->GetNodeIndices(facenodevec);
            // line loop
            outfile << "Line Loop(" << iel << ") = {";
            // line loops require a proper orientation of lines
            for(int iside = nedges; iside<2*nedges; iside++){
                TPZGeoElSide gelside(gel,iside);
                TPZGeoElSide side = gelside.Neighbour();
                // find line element
                while(side.Element()->Dimension() != 1) {side = side.Neighbour();}
                // find first node of line
                int inode = 0;
                while(facenodevec[inode] != side.SideNodeIndex(0)) ++inode;
                // check orientation by comparing second node of line with next node of face
                int64_t index=0;
                    if(side.SideNodeIndex(1)==facenodevec[(inode+1)%nedges]){
                        index = side.Element()->Index();
                    }
                    else{
                        index = -side.Element()->Index();
                    }
                outfile << index <<(iside < 2*nedges-1? "," : "};\n");
            }
            // surface
            outfile << "Surface("<<iel<<") = {"<<iel<<"};\n";
            // @ToDo this is kind of a mess, but only for debugging
            if(gel->MaterialId() >= mtransition && gel->MaterialId() < msurface){groupTransition.push_back(iel);}
            else if(gel->MaterialId() >= msurface){groupSurface.push_back(iel);}
            else groupIntact.push_back(iel);
        }
        // write physical groups
        outfile<<"\nPhysical Surface("<<mtransition<<") = {";
        for(auto itr = groupTransition.begin(); itr != groupTransition.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupTransition.end()? "," : "};\n");
        }
        outfile<<"\nPhysical Surface("<<msurface<<") = {";
        for(auto itr = groupSurface.begin(); itr != groupSurface.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupSurface.end()? "," : "};\n");
        }
        outfile<<"\nPhysical Surface("<<mintact<<") = {";
        for(auto itr = groupIntact.begin(); itr != groupIntact.end();/*Increment in next line*/){
            outfile<<*itr<<(++itr!=groupIntact.end()? "," : "};\n");
        }
    }

    // write volumes
    outfile << "\n\n// VOLUMES DEFINITION \n\n";
    {
        // declare lists to define physical groups
        std::list<int64_t> groupSurface;
        std::list<int64_t> groupTransition;
        std::list<int64_t> groupIntact;
        // iterate over all 3D elements
        for (int64_t iel = 0; iel < nels; iel++){
            TPZGeoEl *gel = pzgmesh->Element(iel);
            if(gel->Dimension() != 3) continue;
            if(gel->HasSubElement()) continue;

            // Surface loop
            // gmsh doesn't accept zero index elements
            outfile << "Surface Loop(" << iel+1 << ") = {";

            // iterate over 2D sides to look for faces that close the surface loop
            int nnodes = gel->NCornerNodes();
            int nsides = gel->NSides();
            bool volumeIsCut = false;
            for(int iside = nnodes; iside < nsides-1; iside++){
                if(gel->SideDimension(iside) != 2) continue;
                TPZGeoElSide gelside(gel,iside);
                // if(gelside.Dimension() < 2) DebugStop();
                // find face element
                TPZGeoElSide side = gelside.Neighbour();
                while(side.Element()->Dimension() != 2) {side = side.Neighbour();}
                // DFNFace *iface = Face(side.Element()->Index());
                // if face is not cut, add it to the loop, else, add its children
                if(side.Element()->HasSubElement() == false){
                    outfile << side.Element()->Index() << (iside < nsides-2? "," : "};\n");
                }
                else{
                    volumeIsCut = true;
                    TPZGeoEl *sidegel = side.Element();
                    PrintYoungestChildren(sidegel,outfile);
                    outfile << (iside < nsides-2? "," : "};\n");
                }
            }

            // volume
            outfile << "Volume("<< iel+1 << ") = {"<< iel+1 <<"};\n"; /* gmsh doesn't accept zero index elements */
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
                    outfile << enclosedSurfaces[i] << (i<nsurfaces-1?",":"} ");
                }
                // -------------------------------------------------------------------------------------
                outfile << "In Volume{"<< iel+1 <<"};\n"; /* gmsh doesn't accept zero index elements */
            }
        }
        // write physical groups
        outfile<<"\nPhysical Volume("<<mtransition<<") = {";
        for(auto itr = groupTransition.begin(); itr != groupTransition.end();/*Increment in next line*/){
            outfile<<*itr+1<<(++itr!=groupTransition.end()? "," : "};\n"); /* gmsh doesn't accept zero index elements */
        }
        outfile<<"\nPhysical Volume("<<mintact<<") = {";
        for(auto itr = groupIntact.begin(); itr != groupIntact.end();/*Increment in next line*/){
            outfile<<*itr+1<<(++itr!=groupIntact.end()? "," : "};\n"); /* gmsh doesn't accept zero index elements */
        }
    }
    
    outfile<<"\nTransfinite Surface {Physical Surface("<<mintact<<")};\n";
    outfile<<"Recombine Surface {Physical Surface("<<mintact<<")};\n";
    outfile<<"\nTransfinite Volume {Physical Volume("<<mintact<<")};\n";

}



















/**
 *  @brief Navigate children tree to access most extreme branches
 *  @param gel pointer to geometric element of eldest ancestor
 *  @param outfile ofstream in which to write accessed data
 */
void DFNMesh::PrintYoungestChildren(TPZGeoEl *gel, std::ofstream &outfile){
    
    int nchildren = gel->NSubElements();
    for(int i = 0; i<nchildren; i++){
        TPZGeoEl *ichild = gel->SubElement(i);
        if(ichild->HasSubElement()){
            PrintYoungestChildren(ichild,outfile);
        }else{
            outfile << ichild->Index();
        }
		outfile << (i < nchildren-1? "," : "");
    }
}























void ReadFractureFromFile(std::string filename, TPZFMatrix<REAL> &plane){
	//    Reading coordinates of a plane from txt file
	string value;
	int npoints = 0;
	string line;
	// count number of corners
	while (npoints == 0){
		ifstream plane_file(filename);
		if (!plane_file){
			std::cout << "Error reading file" << std::endl;
			DebugStop();
		}
		getline(plane_file, line);
		std::stringstream ss(line);
		while (getline(ss, value, ' ')){
			while (value.length() == 0){
				getline(ss, value, ' ');
			}
			npoints++;
		}
	}
	// then read points
	plane.Resize(3,npoints);
	{ //just a scope
		int i = 0;
		int j = 0;
		std::cout << "Fracture plane defined as: \n";
		ifstream plane_file(filename);
		while (getline(plane_file, line)){
			std::stringstream ss(line);
			while (getline(ss, value, ' ')){
				while (value.length() == 0){
					getline(ss, value, ' ');
				}
				plane(i, j) = std::stod(value);
				std::cout << plane(i, j) << (j<npoints-1?" , ":"\n");
				j++;
			}
			j = 0;
			i++;
		}
		std::cout << std::endl;
	}
}


TPZGeoMesh* ReadExampleFromFile(std::string filename, TPZManVector<TPZFMatrix<REAL> > &planevector){
	/*_______________________________________________________________
						FILE FORMAT 

		Domain 								//dimensions of the domain
		Lx Ly Lz (double)

		Mesh			
		2D ELTYPE (string)
		Nx Ny Nz (int)						// number of divisions

		NumberOfFractures [N] (int)

		Fracture 0	[ncorners]				// coordinates for j corner nodes
		[x0] [x1] ... [xj]
		[y0] [y1] ... [yj]
		[z0] [z1] ... [zj]

		Fracture 1	[ncorners]
		...
		Fracture N	[ncorners]
		_______________________________________________________________
						EXAMPLE
		Domain
		2.0 2.0 2.0

		Mesh
		EQuadrilateral
		2 2 2

		NumberOfFractures 2

		Fracture 0 4
		0.3 1.3 1.3 0.3
		0.3 0.3 1.4 1.4
		0.7 0.7 0.5 0.5
		
		Fracture 1 4
		0.6 1.6 1.6 0.6
		0.55 0.7 0.55 0.40
		1.3 1.3 0.25 0.25
	 */

	string line, word;
	// const string Domain = "Domain";
	// const string Mesh("Mesh"), Fractures("NumberOfFractures");
	int i, j, nfractures;
	REAL Lx, Ly, Lz;
	MElementType eltype;
	int nx, ny, nz;

	// Read file
	ifstream plane_file(filename);
	if (!plane_file){
		std::cout << "Error reading file" << std::endl;
		DebugStop();
	}
	// Go through it line by line
	while (getline(plane_file, line)){
		{
			std::stringstream ss(line);
			getline(ss, word, ' ');
		}
		if(word == "Domain"){
			getline(plane_file, line);
			std::stringstream ss(line);
			getline(ss, word, ' ');
			while (word.length() == 0){getline(ss, word, ' ');}
			Lx = std::stod(word);
			getline(ss, word, ' ');
			while (word.length() == 0){getline(ss, word, ' ');}
			Ly = std::stod(word);
			getline(ss, word, ' ');
			while (word.length() == 0){getline(ss, word, ' ');}
			Lz = std::stod(word);
		}

		{
			getline(plane_file, line);
			while(line.length() == 0){getline(plane_file, line);}
			std::stringstream ss(line);
			getline(ss, word, ' ');
		}
		if(word == "Mesh"){
			getline(plane_file, line);
			{
				std::stringstream ss(line);
				getline(ss, word, ' ');
				while (word.length() == 0){getline(ss, word, ' ');}
				if(word == "EQuadrilateral"){
					eltype = EQuadrilateral;
				}else if(word == "ETriangle"){
					eltype = ETriangle;
				}else{ std::cout<<"\nError reading file\n"; DebugStop();}
			}
			getline(plane_file, line);
			{
				std::stringstream ss(line);
				getline(ss, word, ' ');
				while (word.length() == 0){getline(ss, word, ' ');}
				nx = std::stoi(word);
				getline(ss, word, ' ');
				while (word.length() == 0){getline(ss, word, ' ');}
				ny = std::stoi(word);
				getline(ss, word, ' ');
				while (word.length() == 0){getline(ss, word, ' ');}
				nz = std::stoi(word);
			}
		}

		{
			getline(plane_file, line);
			while(line.length() == 0){getline(plane_file, line);}
			std::stringstream ss(line);
			getline(ss, word, ' ');
			if(word == "NumberOfFractures") {
				getline(ss, word, ' ');
				while (word.length() == 0){getline(ss, word, ' ');}
				nfractures = std::stoi(word);
			}
		}
		planevector.resize(nfractures);
		for(int ifrac = 0; ifrac < nfractures; ifrac++){
			string aux;
			while(aux != "Fracture"){
				getline(plane_file, line);
				std::stringstream ss(line);
				getline(ss, aux, ' ');
			}
			int fracid;
			int ncorners;
			{
				std::stringstream ss(line);
				getline(ss, word, ' ');
				while (word.length() == 0 || word == "Fracture"){getline(ss, word, ' ');}
				fracid = std::stoi(word);
				getline(ss, word, ' ');
				while (word.length() == 0 || word == "Fracture"){getline(ss, word, ' ');}
				ncorners = std::stoi(word);
			}
			planevector[ifrac].Resize(3,ncorners);
			std::cout<<"\nCorners of fracture #"<<fracid<<":\n";
			for(int i=0; i<3; i++){
				getline(plane_file, line);
				std::stringstream ss(line);
				int j = 0;
				while (getline(ss, word, ' ')){
					while (word.length() == 0){getline(ss, word, ' ');}
					planevector[ifrac](i, j) = std::stod(word);
					std::cout << planevector[ifrac](i, j) << (j<ncorners-1?", \t":"\n");
					j++;
				}
			}
			std::cout<<"\n";
		}
	}




	// Creating the Geo mesh
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	TPZManVector<REAL, 3> x0(3, 0.), x1(3, 0.);
	x1[0] = Lx;
	x1[1] = Ly;
	x1[2] = 0.;
	TPZManVector<int, 2> ndiv(2);
	ndiv[0] = nx;
	ndiv[1] = ny;
	TPZGenGrid gengrid(ndiv, x0, x1);
	gengrid.SetElementType(eltype);
	gmesh->SetDimension(2);
	gengrid.Read(gmesh);

	// Mesh 3D
	Lz = Lz/nz;
	TPZExtendGridDimension extend(gmesh,Lz);
	extend.SetElType(1);
	TPZGeoMesh *gmesh3d = extend.ExtendedMesh(nz);
	gmesh = gmesh3d;
	return gmesh;
}
