
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

MMeshType StringToMeshtype(const std::string& name){
	if(name[0] != 'E'){}
    if(name == "EQuadrilateral"){	return MMeshType::EQuadrilateral;}
    if(name == "ETriangular"){		return MMeshType::ETriangular;}
    if(name == "EHexahedral"){		return MMeshType::EHexahedral;}
    if(name == "ETetrahedral"){		return MMeshType::ETetrahedral;}
    if(name == "EPyramidal"){		return MMeshType::EPyramidal;}
    if(name == "EPrismatic"){		return MMeshType::EPrismatic;}
    if(name == "EHexaPyrMixed"){	return MMeshType::EHexaPyrMixed;}
    if(name == "ENoType"){			return MMeshType::ENoType;}
	DebugStop();
	return MMeshType::ENoType;
}

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
TPZGeoMesh* SetupExampleFromFile(std::string filename, TPZStack<TPZFMatrix<REAL> > &polyg_stack, std::string mshfile, REAL& toldist, REAL& tolangle,TPZManVector<int>& matid,TPZManVector<FracLimit>& limit_directives);


void ReadFile(	const std::string			& filename, 
				TPZStack<TPZFMatrix<REAL>>	& polygonmatrices, 
				std::string 				& mshfile,
				TPZManVector<REAL,3> 		& x0,
				TPZManVector<REAL,3> 		& xf,
				TPZManVector<int,3> 		& nels,
				MMeshType					& eltype,
				TPZManVector<REAL,2>		& tol,
				TPZManVector<int>			& matid,
				TPZManVector<FracLimit>		& limit_directives
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
	std::cout<<"Running...\n";

}

TPZGeoMesh* ReadInput(int argc, char* argv[], TPZStack<TPZFMatrix<REAL>> &polyg_stack, REAL &toldist, REAL &tolangle,TPZManVector<int>& matid,TPZManVector<FracLimit>& limit_directives);


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
	TPZManVector<int> matid;
	TPZManVector<FracLimit> limit_directives;
	gmesh = ReadInput(argc,argv,polyg_stack,tol_dist,tol_angle,matid,limit_directives);
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
		// fracture = new DFNFracture(polygon,&dfn,FracLimit::Etruncated);
		fracture = new DFNFracture(polygon,&dfn,limit_directives[iplane]);
		dfn.AddFracture(fracture);
        
		// Find intersected ribs and impose tolerance
		fracture->FindRibs();
		fracture->SnapIntersections_ribs(tol_dist);
		// Find intersected faces
		fracture->FindFaces();
		// 
		// fracture->IsolateFractureLimits();
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
#ifdef PZDEBUG
		std::ofstream logtest("LOG/dfnprint.log");
		dfn.Print(logtest,argv[1]);
#endif //PZDEBUG
	}
	// Recover Limits
	for(auto frac : dfn.FractureList()){
		frac->RecoverFractureLimits();
	}

	// Mesh intersected volumes
    // dfn.ExportGMshCAD("dfnExport.geo"); // this is optional, I've been mostly using it for graphical debugging purposes
		// dfn.GenerateSubMesh();
	
	if(polyg_stack.size() == 0){std::cout<<"\nNo fractures were recognized.\n";}
	time.stop();
	std::cout<<"\nTotal running time:\n"<<time<<" ms"<<std::endl;
	//Print graphics
	for(auto frac : dfn.FractureList()){
		frac->Polygon().InsertGeomRepresentation(gmesh);
	}
	// dfn.PrintVTK("pzmesh.txt","skip");
	// dfn.PrintVTK();
	dfn.PrintVTKColorful();
	std::cout<<"\n ...the end.\n\n\a";

	gmsh::finalize();
	return 0;
}



// Takes program input and creates a mesh, matrices with the point coordinates, and writes tolerances
TPZGeoMesh* ReadInput(int argc, char* argv[], TPZStack<TPZFMatrix<REAL>> &polyg_stack, REAL &toldist, REAL &tolangle,TPZManVector<int>& matid,TPZManVector<FracLimit>& limit_directives){
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
			PZError << "\nUnrecognized arguments passed:\n\t\""<<argv[iarg]<<"\" \""<<argv[++iarg]<<"\"\n\n"; 
			DebugStop();
		}
	}
	std::cout<<"input file: "<<example<<"\n";
	gmesh = SetupExampleFromFile(example,polyg_stack,mshfile,toldist,tolangle,matid,limit_directives);
	return gmesh;
}


#include "TPZGenGrid3D.h"



TPZGeoMesh* SetupExampleFromFile(std::string filename, TPZStack<TPZFMatrix<REAL> > &polyg_stack, std::string mshfile, REAL& toldist, REAL& tolangle,TPZManVector<int>& matid,TPZManVector<FracLimit>& limit_directives){


	MMeshType eltype;
	TPZStack<Matrix> planestack;
	TPZManVector<REAL, 3> x0(3, 0.), x1(3, 0.);
	TPZManVector<int,3> nels(3,0);
	TPZManVector<REAL,2> tol = {toldist,tolangle};
	
	ReadFile(filename,polyg_stack,mshfile,x0,x1,nels,eltype,tol,matid,limit_directives);
	toldist = tol[0];
	tolangle = tol[1];


	// Creating the Geo mesh
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	if(eltype != MMeshType::ENoType){
		if(!nels[2]){ // 2D mesh
			TPZManVector<int, 2> ndiv(2);
			ndiv[0] = nels[0];
			ndiv[1] = nels[1];
			TPZGenGrid2D gengrid(ndiv, x0, x1);
			gengrid.SetElementType(eltype);
			gengrid.SetRefpatternElements(true);
			gengrid.Read(gmesh);
			gmesh->SetDimension(2);

		}else{ // 3D mesh
			TPZGenGrid3D gengrid(x0,x1,nels,eltype);
			gmesh = gengrid.BuildVolumetricElements(1);
		}
	}else{ // mesh file
		TPZGmshReader reader;
		gmesh = reader.GeometricGmshMesh4(mshfile, gmesh);
	}
	return gmesh;
}





void ReadFile(	const std::string			& filename, 
				TPZStack<TPZFMatrix<REAL>>	& polygonmatrices, 
				std::string 				& mshfile,
				TPZManVector<REAL,3> 		& x0,
				TPZManVector<REAL,3> 		& xf,
				TPZManVector<int,3> 		& nels,
				MMeshType					& eltype,
				TPZManVector<REAL,2>		& tol,
				TPZManVector<int>			& matid,
				TPZManVector<FracLimit>		& limit_directives
				)
{
		/*_______________________________________________________________
						FILE FORMAT 

		Domain 								// dimensions of the domain
		Lx Ly Lz (double)

		Origin								// position of x0 to generate TPZGenGrid
		x0 y0 z0 (double)

		Mesh								
		2D ELTYPE (string)					// to create a TPZGenGrid
		Nx Ny Nz (int)						// number of divisions

		Mesh
		"/path/to/msh/file.msh"				// to read a .msh file (either quotes or double quotes will work and are required)

		TolDist d (double)					// (optional) Tolerable distance

		TolAngle a (double)					// (optional) Tolerable angle (use TolCos to prescribe cosine)

		Fracture 0							// Fracture keyword and index
		Limit directive						// (optional) Directive for limit handling (options are limited to enums available in FracLimit)
		[x0] [x1] ... [xj]					// Corner coordinates matrix
		[y0] [y1] ... [yj]
		[z0] [z1] ... [zj]

		Fracture 1
		...
		Fracture N
		_______________________________________________________________
						EXAMPLE
		Origin
		-1.0 -1.0 -1.0

		Domain
		2.0 2.0 2.0

		Mesh
		EHexahedral
		1 1 1

		tolDist 0.0001

		Fracture 0
		Limit Eextended
		2.1213 -1.2247 -2.1213  1.2247
		2.1213  1.2247 -2.1213 -1.2247
		0.000  -2.4495  0.000   2.4495
	 */
	string line, word;
	eltype = MMeshType::ENoType;
	// const string Domain = "Domain";
	// const string Mesh("Mesh"), Fractures("NumberOfFractures");
	int i, j, nfractures, ifrac=0;
	REAL Lx, Ly, Lz;
	TPZManVector<REAL,3> L(3,-1.);
	x0.Fill(0.);
	// Read file
	ifstream file(filename);
	if (!file){
		std::cout << "\nCouldn't find file \"" << filename << "\""<< std::endl;
		DebugStop();
	}
	// Go through it line by line
	while (getline(file, line)){
		if(line.length() == 0) continue;
		std::stringstream ss(line);
		getline(ss, word, ' ');
		while (word.length() == 0){getline(ss, word, ' ');}
		if(word[0]=='#') continue; // comment
		// Domain dimensions
		if(word == "Domain"){
			getline(file,line);
			ss.clear();
			ss.str(line);

			for(int i=0; i<3; i++){
				getline(ss,word, ' ');
				while (word.length() == 0){getline(ss, word, ' ');}
				L[i] = std::stod(word);
			}
			continue;
		}
		// Grid origin position
		else if(word == "Origin"){
			getline(file,line);
			ss.clear();
			ss.str(line);
			for(int i=0; i<3; i++){
				getline(ss,word, ' ');
				while (word.length() == 0){getline(ss, word, ' ');}
				x0[i] = std::stod(word);
			}
			continue;
		}
		// Mesh/grid definition
		else if(word == "Mesh"){
			getline(file,line);
			ss.clear();
			ss.str(line);
			getline(ss,word, ' ');
			if(word[0] == 'E'){
				// if(word == "EQuadrilateral") eltype = MMeshType::EQuadrilateral;
				// else if(word == "ETriangle" || word == "ETriangular") eltype = MMeshType::ETriangular;
				// else {PZError << "\nUnrecognized mesh type\n"<<word<<std::endl; DebugStop();}
				eltype = StringToMeshtype(word);

				getline(file,line);
				ss.clear();
				ss.str(line);
				for(int i=0; i<3; i++){
					getline(ss,word, ' ');
					while (word.length() == 0){getline(ss, word, ' ');}
					nels[i] = std::stoi(word);
				}
			}else if(word[0] == '\"' || word[0] == '\''){
				// mshfile = word.substr(1, word.size()-2);
				mshfile = line.substr(1, line.size()-2);
			}else{
				PZError<<"\nUnrecognized mesh info syntax\n\""<<line<<"\""; 
				PZError<<"\nMeshtypes should start with an uppercase E, and .msh file path should be written between quotes or double-quotes.\n";
				DebugStop();
			}
			continue;
		}
		// Fracture reading
		else if(word == "Fracture" || word == "fracture"){
			int ncorners = 1;
			while (getline(ss, word, ' ') && word.length() == 0){/*void*/};
			ifrac = std::stoi(word);
			int npolygons = MAX(ifrac+1,polygonmatrices.size());
			matid.Resize(npolygons,DFNMaterial::Efracture);
			limit_directives.Resize(npolygons,FracLimit::Eextended);
			while (getline(file, line)){
			// while (polygonmatrices[ifrac].Rows() !=3 && getline(file, line)){
				ss.clear();
				ss.str(line);
				while (getline(ss, word, ' ') && word.length() == 0){/*void*/};
				if(word.length() == 0) continue; // to skip empty line with spaces
				// Material id
				if(word[0] == 'm' || word[0] == 'M'){
					while (getline(ss, word, ' ') && word.length() == 0){/*void*/};
					matid[ifrac] = std::stoi(word);
				}
				// Limit handling
				else if(word == "limit" || word == "Limit"){
					while (getline(ss, word, ' ') && word.length() == 0){/*void*/};
					limit_directives[ifrac] = DFN::StringToFracLimit(word);
				}
				// Corner coordinates
				else{
					// Count number of corners from the number of columns in the first line of the matrix in the stringstream
					while (getline(ss, word, ' ')){ncorners += word.length() != 0;}
					ss.clear();
					ss.seekg(0);
					if(ncorners < 3){
						PZError<<"\nInvalid input file.\t Invalid number of corners.\n"; 
						PZError<<"\tFile:"<<filename<<"\n\tFracture index: "<<ifrac<<"\n\tNumber of corners: "<<ncorners<<"\n\n"; 
						DebugStop();
					}
					polygonmatrices.resize(npolygons);
					if(polygonmatrices[ifrac].Cols()>2){
						PZError<<"\nInvalid input file.\t Do you have fractures with repeated indices?\n"; 
						PZError<<"\tFile:"<<filename<<"\n\tFracture index:"<<ifrac<<"\n\n"; 
						DebugStop();
					}
					polygonmatrices[ifrac].Resize(3,ncorners);
					for(int i=0; i<3; i++){
						ss.clear();
						ss.str(line);
						int j = 0;
						while (getline(ss, word, ' ')){
							while (word.length() == 0){getline(ss, word, ' ');}
							polygonmatrices[ifrac](i, j) = std::stod(word);
							// std::cout << std::setw(14) << std::setprecision(6) << std::right << polygonmatrices[ifrac](i, j) << (j<ncorners-1?",":"\n");
							j++;
						}
						if(j!=ncorners){
							PZError<<"\nInvalid input file.\n\t Maybe coordinate matrix has inconsistant number of components?\n";
							PZError<<"\tFile:"<<filename<<"\n\tFracture index:"<<ifrac<<"\n\n"; 
							DebugStop();
						}
						if(i==2) break;
						getline(file, line);
						while(line.length() == 0){getline(file, line);}
					}
					break;
				}
			}
			continue;
		}
		else if(word == "toldist" || word == "tolDist" || word == "TolDist" || word == "TolerableDistance"){
			while (getline(ss, word, ' ') && word.length() == 0){/*void*/};
			tol[0] = std::stod(word);
			continue;
		}
		else if(word == "tolangle" || word == "tolAngle" || word == "TolAngle" || word == "TolerableAngle"){
			while (getline(ss, word, ' ') && word.length() == 0){/*void*/};
			getline(ss, word, ' ');
			tol[1] = std::stod(word);
			continue;
		}
		else if(word == "tolcos" || word == "tolCos" || word == "TolCos" || word == "TolerableCosine"){
			while (getline(ss, word, ' ') && word.length() == 0){/*void*/};
			getline(ss, word, ' ');
			tol[1] = std::acos(std::stod(word));
			continue;
		}

	}
	for(int i=0; i<3; i++) {xf[i] = x0[i] + L[i];}

	// Test if informed element type matches grid definition
	if(eltype != MMeshType::ENoType){
		int dim = 3;
		if(nels[2] == 0 || L[2] == 0.0 ||
		nels[1] == 0 || L[1] == 0.0 ||
		nels[0] == 0 || L[0] == 0.0)
		{
			dim = 2;
		}
		if(dim != MMeshType_Dimension(eltype)) {
			PZError << "\n\nInvalid mesh type ("<<eltype<<") for a "<<dim<<"D domain.\n\n";
			DebugStop();
		}
	}
	
}