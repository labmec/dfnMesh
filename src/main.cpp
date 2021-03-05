
//includes
	#ifdef HAVE_CONFIG_H
		#include <pz_config.h>
	#endif
	
	#include "TPZGenGrid2D.h"
	#include "TPZGenGrid3D.h"
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
	#include "util/json.hpp"
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
void ReadFileJSON(const std::string			& filename, 
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
// TPZGeoMesh* ReadInputJSON(int argc, char* argv[], TPZStack<TPZFMatrix<REAL>> &polyg_stack, REAL &toldist, REAL &tolangle,TPZManVector<int>& matid,TPZManVector<FracLimit>& limit_directives);


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("dfn.mesh"));
#endif




// void ScriptForBug2(TPZGeoMesh* gmesh){
// 	std::vector<int> ref0 = {20,0,1,2};
// 	std::vector<int> ref1 = {21,3,4,5};
// 	std::vector<int> ref2 = {22,6,7,8};
// 	std::vector<int> ref3 = {23,9,10,11};
// 	std::vector<int> ref4 = {19,12,13,14,15,16,17};

// 	std::vector<std::vector<int>> refs = {ref0, ref1, ref2, ref3, ref4};

// 	for(auto& ref : refs){
// 		TPZGeoEl* father = gmesh->Element(ref[0]);
// 		int nchildren = ref.size()-1;
// 		TPZManVector<TPZGeoEl*,6> children(nchildren,nullptr);
// 		for(int i=1; i<=nchildren; i++){
// 			children[i-1] = gmesh->Element(ref[i]);
// 		}
// 		DFN::CreateRefPattern(father,children);
// 		for(int i=0; i<nchildren;i++) father->SetSubElement(i,children[i]);
// 		children.Fill(nullptr);
// 		children.clear();
// 	}
// }


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
	REAL tol_angle = 1.e-6; 
	TPZManVector<int> matid;
	TPZManVector<FracLimit> limit_directives;
	gmesh = ReadInput(argc,argv,polyg_stack,tol_dist,tol_angle,matid,limit_directives);
	gmsh::initialize();
	
    /// Constructor of DFNMesh initializes the skeleton mesh
	// ScriptForBug2(gmesh);
	time.start();
	DFNMesh dfn(gmesh);
	dfn.SetTolerances(tol_dist,tol_angle);

	// dfn.PrintPolyhedra();


    // Loop over fractures and refine mesh around them
	for(int iplane = 0, nfractures = polyg_stack.size(); iplane < nfractures; iplane++){
        // a polygon represents a set of points in a plane
        // poly_stack[iplane] is a matrix 3xn where n is the number of points 
		DFNPolygon polygon(polyg_stack[iplane], gmesh);
        // Initialize the basic data of fracture
		// fracture = new DFNFracture(polygon,&dfn,FracLimit::Etruncated);
        // initialize an empty DFNFracture object
		DFNFracture *fracture = dfn.CreateFracture(polygon,limit_directives[iplane]);
        
		// Find intersected ribs and create a corresponding DFNRib object (administered by DFNFracture)
        // @predro document what this method does
		fracture->FindRibs();
        // For the rib object whose intersection point is closer than toldist to an endnode,
        // snap the intersection point to the nearest node
		fracture->SnapIntersections_ribs(tol_dist);
		// create DFNFace objects - these are faces that have intersected ribs (always 2)
		fracture->FindFaces();
		fracture->SnapIntersections_faces(tol_dist,tol_angle);

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
			// dfn.InheritPolyhedra();
			fracture->MeshFractureSurface();
			// dfn.DumpVTK();
			dfn.UpdatePolyhedra();
		}
#ifdef PZDEBUG
		std::ofstream logtest("LOG/dfn.summary.txt");
		dfn.Print(logtest,argv[1]);
#endif //PZDEBUG
	}
	// Recover Limits
	for(auto frac : dfn.FractureList()){
		frac->RecoverFractureLimits();
	}
	


	// Generate submesh
    dfn.ExportGMshCAD("dfnExport.geo");
	
	if(polyg_stack.size() == 0){std::cout<<"\nNo fractures were recognized.\n";}
	time.stop();
	std::cout<<"\nTotal running time:\n"<<time<<" ms"<<std::endl;
	//Print graphics
	dfn.DumpVTK();
	std::cout<<"\n ...the end.\n\n";

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





TPZGeoMesh* SetupExampleFromFile(std::string filename, TPZStack<TPZFMatrix<REAL> > &polyg_stack, std::string mshfile, REAL& toldist, REAL& tolangle,TPZManVector<int>& matid,TPZManVector<FracLimit>& limit_directives){


	MMeshType eltype = MMeshType::ENoType;
	TPZStack<Matrix> planestack;
	TPZManVector<REAL, 3> x0(3, 0.), x1(3, 0.);
	TPZManVector<int,3> nels(3,0);
	TPZManVector<REAL,2> tol = {toldist,tolangle};

	// Choose extension and read file
	std::string extension;
	try{extension = filename.substr(filename.find_last_of('.'));}
	catch(std::out_of_range){extension = '?';}
	if(extension == ".txt")	
		{ReadFile(filename,polyg_stack,mshfile,x0,x1,nels,eltype,tol,matid,limit_directives);}
	else if(extension == ".json" || extension == ".jsonc")
		{ReadFileJSON(filename,polyg_stack,mshfile,x0,x1,nels,eltype,tol,matid,limit_directives);}
	else{
		PZError << "\nUnrecognized file extension:"
				<< "\nFile = " << filename
				<< "\nExtension = " << extension << std::endl;
		DebugStop();
	}
	
	// Get tolerances
	toldist = tol[0];
	tolangle = tol[1];


	// Creating the Geo mesh
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	if(eltype != MMeshType::ENoType){
		if(!nels[2]){ // 2D mesh
            if(nels[0]==0 || nels[1] == 0) DebugStop();
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
		MatID m (int)						// (optional) Material id of this fracture. If not provided, defaults to DFNMaterial::Efracture
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







void ReadFileJSON(const std::string			& filename, 
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
						JSON OPTIONS 
	{												// <Don't forget opening brace>
		"PZGenGrid": {								// To create a PZGenGrid
			"Origin": [x, y, z], 		<double>	// (optional) for the origin point of the grid. Assumes [0,0,0] if not provided
			// "x0": 								// (alias for "Origin")
			// "minX": 								// (alias for "Origin")
			"Endpoint": [x, y, z], 		<double>	// (optional) for the end point of the grid. Uses "Dimensions" + "Origin" if not provided
			// "xf": 								// (alias for "Endpoint")
			// "maxX": 								// (alias for "Endpoint")
			"Dimensions": [Lx, Ly, Lz],	<double>	// (optional) for the dimensions of the grid. Uses "Endpoint" - "Origin" if not provided
			"MMeshType": "type",		<string>	// Mesh type (EHexahedral, ETetrahedral, ...)
			"Nels": [Nx, Ny, Nz],		<int>		// Number of subdivisions of each direction of the grid
			// "nelDiv": 							// (alias for "Nels")
		},
		"Mesh": "path",					<string>	// Path for .msh file
		"TolDist": value,				<double>	// Tolerable distance (optional). Assumes 1e-4 if not provided
		"TolAngle": value,				<double>	// Tolerable angle (optional). Assumes 1e-4 if not provided
		"Fractures":[
			{
				"Index": i,				<int>		// Fracture index
				"Limit": "directive",	<string>	// Limit handling directive for this fracture {Etruncated, Eextended, Erecovered}
				"MaterialID": id, 		<int>		// Material id. Assumes DFNMaterial::Efracture if not provided
				"Nodes":[				<double>	// Node coordinates for this fracture
					[x0, y0, z0],
					[    ...   ],
					[xN, yN, zN]
				]
			}
			,{
				<Next fracture>
			}
		]
	}												// <Don't forget closing brace>
	______________________________________________________________
						EXAMPLE
	{
		"PZGenGrid":{
			"minX": [0.0, 0.0, 0.0],
			"Dimensions": [2.0, 2.0, 2.0],
			"MMeshType": "EHexahedral",
			"Nels": [1,1,1]
		},
		// "Mesh" : "examples/cube.msh",
		"TolDist": 1e-4,
		"TolAngle": 1e-2,
		"Fractures":[
			{ //Fracture 0
				"Index": 0,
				"Limit": "Etruncated",
				"MaterialID": 5,
				"Nodes":[
					[0.95, 2.50,-0.50],
					[0.95,-0.50,-0.50],
					[0.95,-0.50, 2.50],
					[0.95, 2.50, 2.50]
				]
			}
			,{//Fracture 1
				"Index": 1,
				"Limit": "Etruncated",
				"MaterialID": 5,
				"Nodes":[
					[0.45, 2.50,-0.50],
					[0.45,-0.50,-0.50],
					[0.45,-0.50, 2.50],
					[0.45, 2.50, 2.50]
				]
			}
		]
	}
	 */
	using json = nlohmann::json;


	// Read file
	std::ifstream file(filename);
	if (!file){
		std::cout << "\nCouldn't find file \"" << filename << "\""<< std::endl;
		DebugStop();
	}

	// Parse json
	json input;
	// file >> input;
	input = json::parse(file,nullptr,true,true); // to ignore comments in json file


	// Coarse Mesh
	if(input.find("PZGenGrid") != input.end()){
		if(input.find("Mesh") != input.end()){
			PZError << "\nInput file has a PZGenGrid and a .msh file. Choose one.\n"; DebugStop();
		}
		auto& gengrid = input["PZGenGrid"];

		// minX
		x0.Resize(3,0.0);
		json minX;
		if(gengrid.find("Origin")!= gengrid.end())
			{minX = gengrid["Origin"];}
		else if(gengrid.find("x0") != gengrid.end())
			{minX = gengrid["x0"];}
		else if(gengrid.find("minX") != gengrid.end())
			{minX = gengrid["minX"];}
		else
			{minX = {0.0,0.0,0.0};}
		for(int i=0; i<3; i++){x0[i] = minX[i];}
		
		// maxX
		xf.Resize(3,1.0);
		json maxX;
		if(gengrid.find("Endpoint")!= gengrid.end())
			{maxX = gengrid["Endpoint"];}
		else if(gengrid.find("xf") != gengrid.end())
			{maxX = gengrid["xf"];}
		else if(gengrid.find("maxX") != gengrid.end())
			{maxX = gengrid["maxX"];}
		else
			{maxX = {1.0,1.0,1.0};}
		for(int i=0; i<3; i++){xf[i] = maxX[i];}

		// Dimensions
		if(gengrid.find("Dimensions") != gengrid.end()){
			for(int i=0; i<3; i++){xf[i] = x0[i] + (double)gengrid["Dimensions"][i];}
		}

		// Number of grid divisions
		if(gengrid.find("Nels")!= gengrid.end())
			for(int i=0; i<3; i++){nels[i] = gengrid["Nels"][i];}
		else if(gengrid.find("nelDiv")!= gengrid.end())
			for(int i=0; i<3; i++){nels[i] = gengrid["nelDiv"][i];}

		// MMeshType
		if(gengrid.find("MMeshType") != gengrid.end())
			{eltype = StringToMeshtype((std::string)gengrid["MMeshType"]);}
		else
			{eltype = MMeshType::ENoType;}
	}
	else if(input.find("Mesh") != input.end()){
		mshfile = (std::string)input["Mesh"];
	}
	else{
		PZError << "\nMissing coarse mesh definition.\n"; DebugStop();
	}

	// Tolerances
	if(input.find("TolDist") != input.end())
		{tol[0] = input["TolDist"];}
	else if(input.find("toldist") != input.end())
		{tol[0] = input["toldist"];}
	else if(input.find("Tolerable Distance") != input.end())
		{tol[0] = input["Tolerable Distance"];}
	
	if(input.find("TolAngle") != input.end())
		{tol[1] = input["TolAngle"];}
	else if(input.find("tolangle") != input.end())
		{tol[1] = input["tolangle"];}
	else if(input.find("Tolerable Angle") != input.end())
		{tol[1] = input["Tolerable Angle"];}
	
	// Fractures
	int nfractures = input["Fractures"].size();
	matid.Resize(nfractures,DFNMaterial::Efracture);
	limit_directives.Resize(nfractures,FracLimit::Etruncated);
	polygonmatrices.Resize(nfractures);
	for(auto& fracture : input["Fractures"]){
		int i = fracture["Index"];
		if(polygonmatrices[i].Cols() != 0){
			PZError << "\nInput file has fractures with repeated indices. Index = " << i << "\n"; DebugStop();
		}
		if(fracture.find("MaterialID") != fracture.end()){
			matid[i] = (int)fracture["MaterialID"];
		}
		if(fracture.find("Limit") != fracture.end()){
			limit_directives[i] = DFN::StringToFracLimit((std::string)fracture["Limit"]);
		}
		int npoints = fracture["Nodes"].size();
		polygonmatrices[i].Resize(3,npoints);
		for(int j=0; j<npoints; j++){
			for(int k=0; k<3; k++){
				polygonmatrices[i](k,j) = (REAL)fracture["Nodes"][j][k];
			}
		}
		// polygonmatrices[i].Print(std::cout);
		// std::cout.flush();
	}

}
