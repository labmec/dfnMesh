
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


void ReadFractureFromFile(std::string filename, TPZFMatrix<REAL> &plane);
void ReadExampleFromFile(std::string filename, TPZManVector<TPZFMatrix<REAL> *> &planevector, TPZGeoMesh *gmesh);
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

int main()
{
	TPZGeoMesh *gmesh = new TPZGeoMesh;

	TPZManVector< TPZFMatrix<REAL> *> planevector(2);
	planevector[0] = new Matrix;
	planevector[1] = new Matrix;

	ReadExampleFromFile("example.txt",planevector,gmesh);
	


	// example fractures
	// 2.9084405 1.6236484 0.091559516 1.3763516
	// 2.5516489 2.2694336 2.4483511 2.7305664
	// 1.3832619 2.8898022 1.6167381 0.11019779

	// 7x7x7 (0.5)
	// 2.8 1.2 1.2 2.8
	// 2.8 2.8 1.2 1.2
	// 1.3 1.3 1.3 1.3


	Matrix plane = *planevector[0];
	// ReadFractureFromFile("fracture.txt", plane);
	//  Construction of fracplane and FractureMesh
	DFNFracPlane fracplane(* planevector[0]);
	DFNFractureMesh fracmesh(fracplane, gmesh, 40);

	// Find and split intersected ribs
	gRefDBase.InitializeUniformRefPattern(EOned);
	fracmesh.SplitRibs(19);

	// Find and split intersected faces
	fracmesh.SplitFaces(18);

	// triangulation of fracture plane
	fracmesh.SplitFracturePlane();
	

	// // Mesh transition volumes
	// fracmesh.CreateVolumes();


	//Print result
	std::ofstream meshprint1("meshprint1.txt");
	std::ofstream out1("./TestSurfaces.vtk");
	gmesh->Print(meshprint1);
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out1, true);






// SECOND PLANE
	plane = *planevector[1];
	// ReadFractureFromFile("fracture2.txt", plane);
	// Construction of fracplane and FractureMesh
	DFNFracPlane fracplane2(* planevector[1]);
	DFNFractureMesh fracmesh2(fracplane2, gmesh, 40);

	// Find and split intersected ribs
	fracmesh2.SplitRibs(19);
	// gmesh->BuildConnectivity();
	// Find and split intersected faces
	fracmesh2.SplitFaces(18);
	// triangulation of fracture plane
	fracmesh2.SplitFracturePlane();

	gmesh->BuildConnectivity();
	

	//Print result
	std::ofstream meshprint2("meshprint2.txt");
	std::ofstream out2("./TestSurfaces2.vtk");
	gmesh->Print(meshprint2);
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2, true);

	// Mesh transition volumes
	fracmesh2.CreateVolumes();
	return 0;
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


void ReadExampleFromFile(std::string filename, TPZManVector<TPZFMatrix<REAL> *> &planevector, TPZGeoMesh *gmesh){
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
	const string Domain = "Domain";
	const string Mesh("Mesh"), Fractures("NumberOfFractures");
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
			planevector[ifrac]->Resize(3,ncorners);
			std::cout<<"\nCorners of fracture #"<<fracid<<":\n";
			for(int i=0; i<3; i++){
				getline(plane_file, line);
				std::stringstream ss(line);
				int j = 0;
				while (getline(ss, word, ' ')){
					while (word.length() == 0){getline(ss, word, ' ');}
					planevector[ifrac]->operator()(i, j) = std::stod(word);
					std::cout << planevector[ifrac]->operator()(i, j) << (j<ncorners-1?", \t":"\n");
					j++;
				}
			}
			std::cout<<"\n";
		}
	}




	// Creating the Geo mesh

	TPZManVector<REAL, 3> x0(3, 0.), x1(3, 0.);
	x1[0] = Lx;
	x1[1] = Ly;
	x1[2] = 0.;
	TPZManVector<int, 2> ndiv(2);
	ndiv[0] = nx;
	ndiv[1] = ny;
	TPZGenGrid gengrid(ndiv, x0, x1);
	gengrid.SetElementType(eltype);
	// TPZGeoMesh *gmesh = new TPZGeoMesh;
	gmesh->SetDimension(2);
	gengrid.Read(gmesh);

	// Mesh 3D

	TPZExtendGridDimension extend(gmesh,Lz);
	extend.SetElType(1);
	TPZGeoMesh *gmesh3d = extend.ExtendedMesh(nz);
	gmesh = gmesh3d;

}
