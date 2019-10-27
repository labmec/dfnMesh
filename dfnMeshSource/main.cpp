
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
	#include <map>
	#include <vector>

	#include "DFNFractureMesh.h"
	#include "DFNRibs.h"
	#include "DFNFace.h"

	#include "TPZRefPatternDataBase.h"

	#include <gmsh.h>
//includes


void ReadFractureFromFile(std::string filename, TPZFMatrix<REAL> &plane);
TPZGeoMesh* ReadExampleFromFile(std::string filename, TPZManVector<TPZFMatrix<REAL> > &planevector);



struct DFNMesh{
	// private:
		std::list<DFNFractureMesh> fFractures;
		std::map<int64_t, DFNVolume> fVolumes;
	// public:
		/// Pointer to volume of index 'index'
		DFNVolume *Volume(int64_t index){return &fVolumes[index];}
		/// Uses GMsh to mesh volumes cut by fracture plane
		void CreateVolumes();
		/// Uses gmsh API to tetrahedralize volumes
    	void Tetrahedralize(DFNVolume *volume);
    	/// Find the volumetrical element that encloses a 2D element
    	bool FindEnclosingVolume(TPZGeoEl *ifracface);
	// private:
		/**
    	 *  @brief Navigate children tree to access most extreme branches
    	 *  @param gel: Pointer to geometric element of eldest ancestor
    	 *  @param outfile: ofstream in which to write accessed data
    	 */
    	void PrintYoungestChildren(TPZGeoEl *gel, std::ofstream &outfile);
};

//MATERIAL ID MAP
// 1 gmesh (default)
// 4 skeleton
// 12 ribs that will be divided
// 19 children ribs
// 20 mid-fracture cut faces
// 35 end-fracture cut faces
// 40 Fracture plane
// 45 Intersection points in end-faces

using namespace std;

int main(){
	gRefDBase.InitializeUniformRefPattern(EOned);
	TPZManVector< TPZFMatrix<REAL>> planevector(2);
	TPZGeoMesh *gmesh;
	gmesh = ReadExampleFromFile("example.txt",planevector);



// FIRST PLANE
	Matrix plane(planevector[0]);
	//  Construction of fracplane and FractureMesh
	DFNFracPlane fracplane(plane);
	DFNFractureMesh fracmesh(fracplane, gmesh, 40);
	// Find and split intersected ribs
	fracmesh.SplitRibs(19);
	// Find and split intersected faces
	fracmesh.SplitFaces(18);
	// triangulation of fracture plane
	fracmesh.SplitFracturePlane();
	//Print result
	std::ofstream meshprint1("meshprint1.txt");
	std::ofstream out1("./TestSurfaces.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out1, true);
	gmesh->Print(meshprint1);



// SECOND PLANE
	plane = planevector[1];
	// Construction of fracplane and FractureMesh
	DFNFracPlane fracplane2(plane);
	DFNFractureMesh fracmesh2(fracplane2, gmesh, 40);
	// Find and split intersected ribs
	fracmesh2.SplitRibs(19);
	// gmesh->BuildConnectivity();
	// Find and split intersected faces
	fracmesh2.SplitFaces(18);
	// triangulation of fracture plane
	fracmesh2.SplitFracturePlane();
	

	//Print result
	std::ofstream meshprint2("meshprint2.txt");
	std::ofstream out2("./TestSurfaces2.vtk");
	gmesh->Print(meshprint2);
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2, true);

	// Mesh transition volumes
	DFNMesh dfn;
	dfn.fFractures.push_back(fracmesh);
	dfn.fFractures.push_back(fracmesh2);
	// dfn.CreateVolumes();
	return 0;
}














// void DFNMesh::CreateVolumes(){
// 	DFNFractureMesh *fracmesh;
//     TPZAutoPointer<TPZGeoMesh> gmesh = fracmesh->GetGeoMesh();
	

//     // map all volumes that are cut
//     int64_t nels = gmesh->NElements();
// 	for (int64_t iel = 0; iel < nels; iel++){
//         TPZGeoEl *gel = gmesh->Element(iel);
//         if(gel->Dimension() != 3){continue;}
//         int nsides = gel->NSides();
//         // int ncorners = gel->NCornerNodes();
// 		// int nfaces = (int) (nsides+1-2*ncorners)/2 //from Euler's characteristic
// 		// 2D sides wont start before index 9 for any 3D element
//         for (int iside = nsides-2; iside > 0; iside--){
//             TPZGeoElSide gelside(gel,iside);
//             if (gelside.Dimension() != 2){break;}
//             TPZGeoElSide neighbour = gelside.Neighbour();
// 			while(neighbour.Element()->Dimension() != 2 && neighbour != gelside){
// 				neighbour = neighbour.Neighbour();
// 			}
// 			if(neighbour == gelside) continue;
// 			TPZGeoEl *sideface = neighbour.Element();
// 			if(sideface->HasSubElement()){
//                 DFNVolume volume(iel,true);
//                 fVolumes.insert({iel,volume});
//                 break;
//             }            
//         }
//     }
    
//     //Loop over list of fracmeshes
// 	for(auto fracmesh = this->fFractures.begin(); fracmesh != fFractures.end(); fracmesh++){
// 		// search through each 2D element of the triangulated fracture surface to find their enclosing volume
// 		// iterate over fracplane's elements created at SplitFracturePlane
// 		for(int64_t iel = 0; iel < nels; iel++){
// 			TPZGeoEl *gel = gmesh->Element(iel);
// 			// if(gel->MaterialID() != fSurfaceMaterial){continue;}
// 			// During development, elements at fracture surface have material id over 40
// 			if(gel->MaterialId() <= fracmesh->GetSurfaceMaterial()) continue;
// 			if(gel->Dimension() != 2) continue;
// 			if(gel->HasSubElement()) continue;
// 			// Find volume that encloses that element
// 			FindEnclosingVolume(gel);
// 		}
// 	}

// 	//Loop over list of volumes cut
// 	for (auto itr = fVolumes.begin(); itr != fVolumes.end(); itr++){
//     	DFNVolume *ivolume = &itr->second;
// 		// Use GMsh to tetrahedralize volumes
//     	Tetrahedralize(ivolume);
// 	}
	

// }














// // @ToDo: Rewrite this using GMsh API
// void DFNMesh::Tetrahedralize(DFNVolume *volume){
//     int mtransition = 19;
//     int msurface = 40;
//     int mintact = 1;

//     // Giving fGMesh another name for readability's sake
//     TPZAutoPointer<TPZGeoMesh> pzgmesh = fFractures.begin()->GetGeoMesh();
// 	// @ToDo: check if building connectivity is necessary at this point
//     pzgmesh->BuildConnectivity();
//     CreateSkeletonElements(1,mtransition);
//     // Title
//     outfile<<"//  Geo file generated by Discrete Fracture Network methods \n"
//             <<"// Fracture #1 \n\n";
    
//     // write nodes
//     outfile<< "// POINTS DEFINITION \n\n";
//     int64_t nnodes = pzgmesh->NNodes();
//     // @ToDo Do we need physical groups for points too?
//     for (int64_t inode = 0; inode < nnodes; inode++){
//         TPZManVector<REAL, 3> co(3,0.);
//         pzgmesh->NodeVec()[inode].GetCoordinates(co);
//         outfile << "Point(" << inode << ") = {" << co[0] << ',' << co[1] << ',' << co[2] << "};\n";
//     }
    
//     // write edges
//     int64_t nels = pzgmesh->NElements();
//     outfile << "\n\n// LINES DEFINITION \n\n";
//     {
//         // declare lists to define physical groups
//         std::list<int64_t> groupSurface;
//         std::list<int64_t> groupTransition;
//         std::list<int64_t> groupIntact;
//         // iterate over all 1D elements
//         for (int64_t iel = 0; iel < nels; iel++){
//             TPZGeoEl *gel = pzgmesh->Element(iel);
//             if(gel->Dimension() != 1) continue;
//             if(gel->HasSubElement()) continue;
//             // if(gel->MaterialId() == fTransitionMaterial) continue;
//             outfile << "Line(" << iel << ") = {" << gel->NodeIndex(0) << ',' << gel->NodeIndex(1) << "};\n";
//     // @ToDo this is kind of a mess, but only for debugging
//             // list it according to material
//             if(gel->MaterialId() >= mtransition && gel->MaterialId() < msurface){groupTransition.push_back(iel);}
//             else if(gel->MaterialId() >= msurface){groupSurface.push_back(iel);}
//             else groupIntact.push_back(iel);
//         }
//         // write physical groups
//         outfile<<"\nPhysical Curve("<<mtransition<<") = {";
//         for(auto itr = groupTransition.begin(); itr != groupTransition.end();/*Increment in next line*/){
//             outfile<<*itr<<(++itr!=groupTransition.end()? "," : "};\n");
//         }
//         outfile<<"\nPhysical Curve("<<msurface<<") = {";
//         for(auto itr = groupSurface.begin(); itr != groupSurface.end();/*Increment in next line*/){
//             outfile<<*itr<<(++itr!=groupSurface.end()? "," : "};\n");
//         }
//         outfile<<"\nPhysical Curve("<<mintact<<") = {";
//         for(auto itr = groupIntact.begin(); itr != groupIntact.end();/*Increment in next line*/){
//             outfile<<*itr<<(++itr!=groupIntact.end()? "," : "};\n");
//         }
//     }
//     // write faces
//     outfile << "\n\n// FACES DEFINITION \n\n";
//     {
//         // declare lists to define physical groups
//         std::list<int64_t> groupSurface;
//         std::list<int64_t> groupTransition;
//         std::list<int64_t> groupIntact;
//         // iterate over all 2D elements
//         for (int64_t iel = 0; iel < nels; iel++){
//             TPZGeoEl *gel = pzgmesh->Element(iel);
//             if(gel->Dimension() != 2) continue;
//             if(gel->HasSubElement()) continue;
//             // if(gel->MaterialId() == fTransitionMaterial) continue;
            
//             int nnodes = gel->NCornerNodes();
//             int nedges = nnodes; //for readability 
//             TPZManVector<int64_t,4> facenodevec(nnodes);
//             gel->GetNodeIndices(facenodevec);
//             // line loop
//             outfile << "Line Loop(" << iel << ") = {";
//             // line loops require a proper orientation of lines
//             for(int iside = nedges; iside<2*nedges; iside++){
//                 TPZGeoElSide gelside(gel,iside);
//                 TPZGeoElSide side = gelside.Neighbour();
//                 // find line element
//                 while(side.Element()->Dimension() != 1) {side = side.Neighbour();}
//                 // find first node of line
//                 int inode = 0;
//                 while(facenodevec[inode] != side.SideNodeIndex(0)) ++inode;
//                 // check orientation by comparing second node of line with next node of face
//                 int64_t index=0;
//                     if(side.SideNodeIndex(1)==facenodevec[(inode+1)%nedges]){
//                         index = side.Element()->Index();
//                     }
//                     else{
//                         index = -side.Element()->Index();
//                     }
//                 outfile << index <<(iside < 2*nedges-1? "," : "};\n");
//             }
//             // surface
//             outfile << "Surface("<<iel<<") = {"<<iel<<"};\n";
//             // @ToDo this is kind of a mess, but only for debugging
//             if(gel->MaterialId() >= mtransition && gel->MaterialId() < msurface){groupTransition.push_back(iel);}
//             else if(gel->MaterialId() >= msurface){groupSurface.push_back(iel);}
//             else groupIntact.push_back(iel);
//         }
//         // write physical groups
//         outfile<<"\nPhysical Surface("<<mtransition<<") = {";
//         for(auto itr = groupTransition.begin(); itr != groupTransition.end();/*Increment in next line*/){
//             outfile<<*itr<<(++itr!=groupTransition.end()? "," : "};\n");
//         }
//         outfile<<"\nPhysical Surface("<<msurface<<") = {";
//         for(auto itr = groupSurface.begin(); itr != groupSurface.end();/*Increment in next line*/){
//             outfile<<*itr<<(++itr!=groupSurface.end()? "," : "};\n");
//         }
//         outfile<<"\nPhysical Surface("<<mintact<<") = {";
//         for(auto itr = groupIntact.begin(); itr != groupIntact.end();/*Increment in next line*/){
//             outfile<<*itr<<(++itr!=groupIntact.end()? "," : "};\n");
//         }
//     }

//     // write volumes
//     outfile << "\n\n// VOLUMES DEFINITION \n\n";
//     {
//         // declare lists to define physical groups
//         std::list<int64_t> groupSurface;
//         std::list<int64_t> groupTransition;
//         std::list<int64_t> groupIntact;
//         // iterate over all 3D elements
//         for (int64_t iel = 0; iel < nels; iel++){
//             TPZGeoEl *gel = pzgmesh->Element(iel);
//             if(gel->Dimension() != 3) continue;
//             if(gel->HasSubElement()) continue;

//             // Surface loop
//             // gmsh doesn't accept zero index elements
//             outfile << "Surface Loop(" << iel+1 << ") = {";

//             // iterate over 2D sides to look for faces that close the surface loop
//             int nnodes = gel->NCornerNodes();
//             int nsides = gel->NSides();
//             bool volumeIsCut = false;
//             for(int iside = nnodes; iside < nsides-1; iside++){
//                 if(gel->SideDimension(iside) != 2) continue;
//                 TPZGeoElSide gelside(gel,iside);
//                 // if(gelside.Dimension() < 2) DebugStop();
//                 // find face element
//                 TPZGeoElSide side = gelside.Neighbour();
//                 while(side.Element()->Dimension() != 2) {side = side.Neighbour();}
//                 DFNFace *iface = Face(side.Element()->Index());
//                 // if face is not cut, add it to the loop, else, add its children
//                 if(side.Element()->HasSubElement() == false){
//                     outfile << side.Element()->Index() << (iside < nsides-2? "," : "};\n");
//                 }
//                 else{
//                     volumeIsCut = true;
//                     TPZGeoEl *sidegel = side.Element();
//                     PrintYoungestChildren(sidegel,outfile);
//                     outfile << (iside < nsides-2? "," : "};\n");
//                 }
//             }

//             // volume
//             outfile << "Volume("<< iel+1 << ") = {"<< iel+1 <<"};\n"; /* gmsh doesn't accept zero index elements */
//             if(volumeIsCut){groupTransition.push_back(iel);}
//             else groupIntact.push_back(iel);

//             if(volumeIsCut){
//                 int nsurfaces = Volume(iel)->GetFacesInVolume().size();
//                 TPZManVector<int64_t,6> enclosedSurfaces(nsurfaces);
//                 enclosedSurfaces = Volume(iel)->GetFacesInVolume();

//                 outfile << "Surface{";
//                 for(int i = 0; i<nsurfaces; i++){
//                     // TPZGeoEl *surface = fGMesh->Element(enclosedSurfaces[i]);
//                     // if(surface->HasSubElement()){
//                     //     PrintYoungestChildren(surface,outfile);
//                     // }
//                     // else{
//                     //     outfile << enclosedSurfaces[i] << (i<nsurfaces-1?",":"} ");
//                     // }
//                     outfile << enclosedSurfaces[i] << (i<nsurfaces-1?",":"} ");
//                 }
//                 outfile << "In Volume{"<< iel+1 <<"};\n"; /* gmsh doesn't accept zero index elements */
//             }
//         }
//         // write physical groups
//         outfile<<"\nPhysical Volume("<<mtransition<<") = {";
//         for(auto itr = groupTransition.begin(); itr != groupTransition.end();/*Increment in next line*/){
//             outfile<<*itr+1<<(++itr!=groupTransition.end()? "," : "};\n"); /* gmsh doesn't accept zero index elements */
//         }
//         outfile<<"\nPhysical Volume("<<mintact<<") = {";
//         for(auto itr = groupIntact.begin(); itr != groupIntact.end();/*Increment in next line*/){
//             outfile<<*itr+1<<(++itr!=groupIntact.end()? "," : "};\n"); /* gmsh doesn't accept zero index elements */
//         }
//     }
    
//     outfile<<"\nTransfinite Surface {Physical Surface("<<mintact<<")};\n";
//     outfile<<"Recombine Surface {Physical Surface("<<mintact<<")};\n";
//     outfile<<"\nTransfinite Volume {Physical Volume("<<mintact<<")};\n";


// }


































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
