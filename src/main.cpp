
//includes
	#ifdef HAVE_CONFIG_H
		#include <pz_config.h>
	#endif
	
	#include "TPZGenGrid2D.h"
	#include "TPZExtendGridDimension.h"
	#include "TPZRefPatternDataBase.h"
	#include "TPZGmshReader.h"
	#include "pzlog.h"
	#include "TPZTimer.h"
	
	#include <stdio.h>
	#include <math.h>
	#include <iostream>
	#include <fstream>
	#include <string>
	#include <sstream>
	#include <cstdio>
	#include <set>
	#include <unordered_set>
	#include <map>
	#include <vector>

	#include "DFNFracture.h"
	#include "DFNMesh.h"

	#include <gmsh.h>
//includes

/**
 * @brief Calls ReadExampleFromFile but for a single fracture
 * @param filename: path to the file that defines the fracture
 * @param planevector: Matrix to fill with corners of the fracture
*/
void ReadFracture(std::string filename, TPZFMatrix<REAL> &plane);

/**
 * @brief Define which example to run. See example file sintax in function definition
 * @param filename: path to the file that defines the example
 * @param planevector: vector to fill with corners of the fractures
 * @param mshfile: [optional] path to .msh file (if applicable)
 * @returns pointer to geometric mesh created/read
*/
TPZGeoMesh* ReadExampleFromFile(std::string filename, TPZManVector<TPZFMatrix<REAL> > &planevector, std::string mshfile = "no-msh-file");






/**
 * @brief information and assumptions
*/
void PrintPreamble(){
	std::string neopzversion = "/commit/85f5651"; // https://github.com/labmec/neopz/commit/...
	std::string gmshversion = "4.5.6";
	std::cout<<"\n";
	std::cout<<"\nNeoPZ assumed version: " << neopzversion;
	std::cout<<"\nGMsh assumed version: " << gmshversion << "\n\n";
	std::cout<<"Runing...\n\n";

}

TPZGeoMesh* ReadInput(int argc, char* argv[], TPZManVector< TPZFMatrix<REAL>> &planevector, REAL &toldist, REAL &tolangle);









//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

int main(int argc, char* argv[]){
	#ifdef LOG4CXX
    	InitializePZLOG();
	#endif
	TPZTimer time("DFNMesh");
	time.start();
	PrintPreamble();
    /// this data structure defines the fractures which will cut the mesh
    // each matrix is dimension (3xn) where n is the number of vertices
	TPZManVector< TPZFMatrix<REAL>> planevector;
	TPZGeoMesh *gmesh = nullptr;
	REAL tol_dist = 1.e-4;
	REAL tol_angle = 1.e-3; 
	gmesh = ReadInput(argc,argv,planevector,tol_dist,tol_angle);
	gmsh::initialize();
	
    /// Constructor of DFNMesh initializes the skeleton mesh
	DFNMesh dfn(gmesh);
    /// this will initialize a strange data structure allowing for each dim-1 element to know how many volumes are
    // connected
	dfn.InitializeFaceTracker();
	// Loop over fractures and refine mesh around them
	DFNFracture *fracture = nullptr;
	for(int iplane = 0, nfractures = planevector.size(); iplane < nfractures; iplane++){
        // a polygon represents a set of points in a plane
        // @TODO : wouldnt it be wiser to pass the geometric mesh so that the complete
        // data structure can be initialized?
		DFNPolygon polygon(planevector[iplane]);
        
        // Initialize the basic data of fracture
		fracture = new DFNFracture(polygon,&dfn);
        
        // @TODO shouldn't these methods be called in the constructor??
		// @reply they're here temporarily for debug fase. We'll move them wherever's more intuitive in the end.

	// Find and split intersected ribs
		fracture->FindRibs();
		fracture->OptimizeRibs(tol_dist);
		fracture->RefineRibs();
	// Find and split intersected faces
		fracture->FindFaces();
		fracture->RefineFaces();
		dfn.GetPolyhedra();
	// Mesh fracture surface
		fracture->AssembleOutline();
		// fracture->GetSubPolygons();
		if(gmesh->Dimension() == 3){
			fracture->MeshFractureSurface();
		}
	//insert fracture
		dfn.AddFracture(fracture);
	}
	// Mesh transition volumes
		dfn.CreateVolumes();
		// dfn.ExportGMshCAD("dfnExport.geo"); // this is optional, I've been mostly using it for graphical debugging purposes
		// dfn.GenerateSubMesh();
	time.stop();
	std::cout<<"\n\n"<<time;
	//Print graphics
		for(auto frac : dfn.FractureList()){
			frac->Polygon().InsertGeoEl(gmesh);
		}
		dfn.PrintColorful();

	std::cout<<"\n ...the end.\n\n";

	gmsh::finalize();
	return 0;
}




TPZGeoMesh* ReadInput(int argc, char* argv[], TPZManVector< TPZFMatrix<REAL>> &planevector, REAL &toldist, REAL &tolangle){
	TPZGeoMesh* gmesh = nullptr;
	switch(argc){
		case 0:
		case 1: gmesh = ReadExampleFromFile("examples/2D-mult-fracture.txt",planevector); 
				break;
		case 2: gmesh = ReadExampleFromFile(argv[1],planevector);
				break;
		case 5: tolangle = std::stod(argv[4]);
		case 4: toldist = std::stod(argv[3]);
		case 3: gmesh = ReadExampleFromFile(argv[1],planevector,argv[2]);
				break;
		default: PZError << "\n\n Invalid parameters \n\n"; DebugStop();
	}
	return gmesh;
}




void ReadFracture(std::string filename, TPZFMatrix<REAL> &plane){
	TPZManVector<TPZFMatrix<REAL>> planevector;
	ReadExampleFromFile(filename, planevector);
	plane = planevector[0];
}


TPZGeoMesh* ReadExampleFromFile(std::string filename, TPZManVector<TPZFMatrix<REAL> > &planevector, std::string mshfile){
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
	bool create_mesh_Q = mshfile == "no-msh-file";
	// const string Domain = "Domain";
	// const string Mesh("Mesh"), Fractures("NumberOfFractures");
	int i, j, nfractures;
	REAL Lx, Ly, Lz;
	MMeshType eltype;
	int nx, ny, nz;

	// Read file
	ifstream plane_file(filename);
	if (!plane_file){
		std::cout << "\nCouldn't find file " << filename << std::endl;
		DebugStop();
	}
	// Go through it line by line
	while (getline(plane_file, line)){
		if(create_mesh_Q){
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
						eltype = MMeshType::EQuadrilateral;
					}else if(word == "ETriangle" || word == "ETriangular"){
						eltype = MMeshType::ETriangular;
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
			getline(plane_file, line);
		}

		{
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
					std::cout << std::setw(14) << std::setprecision(6) << std::right << planevector[ifrac](i, j) << (j<ncorners-1?",":"\n");
					j++;
				}
			}
			std::cout<<"\n";
		}
	}




	// Creating the Geo mesh
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	if(create_mesh_Q){
		TPZManVector<REAL, 3> x0(3, 0.), x1(3, 0.);
		x1[0] = Lx;
		x1[1] = Ly;
		x1[2] = 0.;
		TPZManVector<int, 2> ndiv(2);
		ndiv[0] = nx;
		ndiv[1] = ny;
		TPZGenGrid2D gengrid(ndiv, x0, x1);
		gengrid.SetElementType(eltype);
		gengrid.SetRefpatternElements(true);
		gengrid.Read(gmesh);
		gmesh->SetDimension(2);

		// Mesh 3D
		if(nz != 0){
			Lz = Lz/nz;
			TPZExtendGridDimension extend(gmesh,Lz);
			extend.SetElType(1);
			TPZGeoMesh *gmesh3d = extend.ExtendedMesh(nz);
			gmesh = gmesh3d;
		}
	}else{
		TPZGmshReader reader;
		gmesh = reader.GeometricGmshMesh4(mshfile, gmesh);
	}
	return gmesh;
}
