
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
 * @brief Calls SetupExampleFromFile but for a single fracture
 * @param filename: path to the file that defines the fracture
 * @param polyg_stack: Matrix to fill with corners of the fracture
*/
void ReadFracture(std::string filename, TPZFMatrix<REAL> &plane);

/**
 * @brief Define which example to run. See example file sintax in function definition
 * @param filename: path to the file that defines the example
 * @param polyg_stack: vector to fill with corners of the fractures
 * @param mshfile: [optional] path to .msh file (if applicable)
 * @returns pointer to geometric mesh created/read
*/
TPZGeoMesh* SetupExampleFromFile(std::string filename, TPZStack<TPZFMatrix<REAL> > &polyg_stack, std::string mshfile, REAL& toldist, REAL& tolangle);


void ReadFile(	const std::string&				filename, 
					TPZStack<TPZFMatrix<REAL>>& polygonmatrices, 
					std::string& 				mshfile,
					TPZManVector<REAL,3>& 		x0,
					TPZManVector<REAL,3>& 		xf,
					TPZManVector<int,3>& 		nels,
					MMeshType&					eltype,
					TPZManVector<REAL,2>&		tol
					);



/**
 * @brief information and assumptions
*/
void PrintPreamble(){
	std::string neopzversion = "/commit/29373e1"; // https://github.com/labmec/neopz/commit/...
	std::string gmshversion = "4.5.6";
	std::cout<<"\n";
	std::cout<<"\nNeoPZ assumed version: " << neopzversion;
	std::cout<<"\nGMsh assumed version: " << gmshversion << "\n\n";
	std::cout<<"Runing...\n\n";

}

TPZGeoMesh* ReadInput(int argc, char* argv[], TPZStack<TPZFMatrix<REAL>> &polyg_stack, REAL &toldist, REAL &tolangle);


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("dfn.fracture"));
#endif







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
	PrintPreamble();
    /// this data structure defines the fractures which will cut the mesh
    // each matrix is dimension (3xn) where n is the number of vertices
	TPZStack<TPZFMatrix<REAL>> polyg_stack;
	TPZGeoMesh *gmesh = nullptr;
	REAL tol_dist = 1.e-4;
	REAL tol_angle = 1.e-3; 
	gmesh = ReadInput(argc,argv,polyg_stack,tol_dist,tol_angle);
	gmsh::initialize();
	
    /// Constructor of DFNMesh initializes the skeleton mesh
	time.start();
	DFNMesh dfn(gmesh);
	// dfn.PrintPolyhedra();
	dfn.SetTolerances(tol_dist,tol_angle);
	// Loop over fractures and refine mesh around them
	DFNFracture *fracture = nullptr;
	for(int iplane = 0, nfractures = polyg_stack.size(); iplane < nfractures; iplane++){
		gmesh->BuildConnectivity();//@todo There is a proper buildconnectivity missing... this is a temporary patching until I find out where it's actually supposed to be
        // a polygon represents a set of points in a plane
		DFNPolygon polygon(polyg_stack[iplane], gmesh);
        // Initialize the basic data of fracture
		fracture = new DFNFracture(polygon,&dfn);
		dfn.AddFracture(fracture);
        
	// Find and split intersected ribs
		fracture->FindRibs();
		fracture->SnapIntersections_ribs(tol_dist);
	// Build the DFNFace objects and split intersected faces if necessary
		fracture->FindFaces();
		fracture->IsolateFractureLimits();
		fracture->SnapIntersections_faces(tol_dist,tol_angle);
		// fracture->Polygon().PlotNodesAbove_n_Below(gmesh);

#ifdef LOG4CXX
        if(logger->isDebugEnabled()){
            std::stringstream sout;
            fracture->Print(sout);
            LOGPZ_DEBUG(logger,sout.str());
        }
#endif
		fracture->RefineRibs();
		fracture->RefineFaces();
	// Mesh fracture surface
		if(gmesh->Dimension() == 3){
			fracture->MeshFractureSurface();
			dfn.UpdatePolyhedra();
		}

		std::ofstream logtest("LOG/dfnprint.log");
		dfn.Print(logtest,argv[1]);
	}
	// Mesh intersected volumes
    // dfn.ExportGMshCAD("dfnExport.geo"); // this is optional, I've been mostly using it for graphical debugging purposes
		// dfn.GenerateSubMesh();
	time.stop();
	std::cout<<"\n\n"<<time<<" ms"<<std::endl;
	//Print graphics
		for(auto frac : dfn.FractureList()){
			frac->Polygon().InsertGeoEl(gmesh);
		}
		dfn.PrintVTK("pzmesh.txt","skip");
		dfn.PrintVTKColorful();
	std::cout<<"\n ...the end.\n\n";

	gmsh::finalize();
	return 0;
}



// Takes program input and creates a mesh, matrices with the point coordinates, and writes tolerances
TPZGeoMesh* ReadInput(int argc, char* argv[], TPZStack<TPZFMatrix<REAL>> &polyg_stack, REAL &toldist, REAL &tolangle){
	TPZGeoMesh* gmesh = nullptr;
	std::string default_example("examples/two-hex-and-a-frac.txt");
	std::string example = default_example;
	std::string mshfile = "no-msh-file";
	for(int iarg=1; iarg < argc; ++iarg){
		std::string aux = argv[iarg];
		if(argv[iarg][0] != '-'){example = argv[iarg];}
		else if(aux == "-m"){mshfile = argv[++iarg];}
		else if(aux == "-f"){example = argv[++iarg];}
		else if(aux == "-td"){toldist = std::stod(argv[++iarg]);}
		else if(aux == "-ta"){tolangle = std::stod(argv[++iarg]);}
		else if(aux == "-tc"){tolangle = std::acos(std::stod(argv[++iarg]));}
		else{
			PZError << "\nUnrecognized arguments passed\n\t\""<<argv[iarg]<<"\" \""<<argv[++iarg]<<"\"\n\n"; 
			DebugStop();
		}
	}
	gmesh = SetupExampleFromFile(example,polyg_stack,mshfile,toldist,tolangle);
	return gmesh;
}






TPZGeoMesh* SetupExampleFromFile(std::string filename, TPZStack<TPZFMatrix<REAL> > &polyg_stack, std::string mshfile, REAL& toldist, REAL& tolangle){


	MMeshType eltype;
	TPZStack<Matrix> planestack;
	TPZManVector<REAL, 3> x0(3, 0.), x1(3, 0.);
	TPZManVector<int,3> nels(3,0);
	TPZManVector<REAL,2> tol = {toldist,tolangle};
	
	ReadFile(filename,polyg_stack,mshfile,x0,x1,nels,eltype,tol);
	toldist = tol[0];
	tolangle = tol[1];


	// Creating the Geo mesh
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	if(mshfile == "no-msh-file"){
		TPZManVector<int, 2> ndiv(2);
		ndiv[0] = nels[0];
		ndiv[1] = nels[1];
		TPZGenGrid2D gengrid(ndiv, x0, x1);
		gengrid.SetElementType(eltype);
		gengrid.SetRefpatternElements(true);
		gengrid.Read(gmesh);
		gmesh->SetDimension(2);

		// Mesh 3D
		if(nels[2] != 0){
			REAL Lz = (x1[2]-x0[2])/nels[2];
			TPZExtendGridDimension extend(gmesh,Lz);
			extend.SetElType(1);
			TPZGeoMesh *gmesh3d = extend.ExtendedMesh(nels[2]);
			gmesh = gmesh3d;
		}
	}else{
		TPZGmshReader reader;
		gmesh = reader.GeometricGmshMesh4(mshfile, gmesh);
	}
	return gmesh;
}





void ReadFile(	const std::string&			filename, 
				TPZStack<TPZFMatrix<REAL>>& polygonmatrices, 
				std::string& 				mshfile,
				TPZManVector<REAL,3>& 		x0,
				TPZManVector<REAL,3>& 		xf,
				TPZManVector<int,3>& 		nels,
				MMeshType&					eltype,
				TPZManVector<REAL,2>&		tol
				)
{
		/*_______________________________________________________________
						FILE FORMAT 

		Domain 								//dimensions of the domain
		Lx Ly Lz (double)

		Mesh			
		2D ELTYPE (string)
		Nx Ny Nz (int)						// number of divisions

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
	int i, j, nfractures, ifrac=0;
	REAL Lx, Ly, Lz;
	TPZManVector<REAL,3> L(3,-1.);
	x0.Fill(0.);
	// Read file
	ifstream file(filename);
	if (!file){
		std::cout << "\nCouldn't find file " << filename << std::endl;
		DebugStop();
	}
	// Go through it line by line
	while (getline(file, line)){
		while(line.length() == 0){getline(file, line);}
		std::stringstream ss(line);
		getline(ss, word, ' ');
		while (word.length() == 0){getline(ss, word, ' ');}
		if(word == "Domain"){
			getline(file,line);
			ss.clear();
			ss.str(line);

			for(int i=0; i<3; i++){
				getline(ss,word, ' ');
				while (word.length() == 0){getline(ss, word, ' ');}
				L[i] = std::stod(word);
			}
		}
		else if(word == "Origin"){
			getline(file,line);
			ss.clear();
			ss.str(line);
			for(int i=0; i<3; i++){
				getline(ss,word, ' ');
				while (word.length() == 0){getline(ss, word, ' ');}
				x0[i] = std::stod(word);
			}
		}
		else if(word == "Mesh"){
			getline(file,line);
			ss.clear();
			ss.str(line);
			getline(ss,word, ' ');
			if(word[0] == 'E'){
				if(word == "EQuadrilateral") eltype = MMeshType::EQuadrilateral;
				else if(word == "ETriangle" || word == "ETriangular") eltype = MMeshType::ETriangular;
				else {PZError << "\nUnrecognized mesh type\n"<<word<<std::endl; DebugStop();}

				getline(file,line);
				ss.clear();
				ss.str(line);
				for(int i=0; i<3; i++){
					getline(ss,word, ' ');
					while (word.length() == 0){getline(ss, word, ' ');}
					nels[i] = std::stoi(word);
				}
			}else if(word[0] == '\"' || word[0] == '\''){
				mshfile = word.substr(1, word.size()-2);
			}else{PZError<<"\nUnrecognized mesh info syntax\n"; DebugStop();}
		}
		else if(word == "Fracture" || word == "fracture"){
			getline(ss,word, ' ');
			ifrac = std::stoi(word);
			getline(ss,word, ' ');
			int ncorners = std::stoi(word);
			int size = MAX(ifrac+1,polygonmatrices.size());
			polygonmatrices.resize(size);
			polygonmatrices[ifrac].Resize(3,ncorners);
			for(int i=0; i<3; i++){
				getline(file, line);
				while(line.length() == 0){getline(file, line);}
				ss.clear();
				ss.str(line);
				int j = 0;
				while (getline(ss, word, ' ')){
					while (word.length() == 0){getline(ss, word, ' ');}
					polygonmatrices[ifrac](i, j) = std::stod(word);
					// std::cout << std::setw(14) << std::setprecision(6) << std::right << polygonmatrices[ifrac](i, j) << (j<ncorners-1?",":"\n");
					j++;
				}
			}
		}
		else if(word == "toldist" || word == "tolDist" || word == "TolDist" || word == "TolerableDistance"){
			getline(ss, word, ' ');
			tol[0] = std::stod(word);
		}
		else if(word == "tolangle" || word == "tolAngle" || word == "TolAngle" || word == "TolerableAngle"){
			getline(ss, word, ' ');
			tol[1] = std::stod(word);
		}
		else if(word == "tolcos" || word == "tolCos" || word == "TolCos" || word == "TolerableCosine"){
			getline(ss, word, ' ');
			tol[1] = std::acos(std::stod(word));
		}

	}
	for(int i=0; i<3; i++) {xf[i] = x0[i] + L[i];}
	std::cout<<std::endl;

	
}