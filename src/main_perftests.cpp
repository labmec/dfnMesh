
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
	#include "json.hpp"
	#include "filereader.h"
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
	// Loop over fractures and refine mesh around them
	DFNFracture *fracture = nullptr;

	// dfn.PrintPolyhedra();


	for(int iplane = 0, nfractures = polyg_stack.size(); iplane < nfractures; iplane++){
        // a polygon represents a set of points in a plane
		DFNPolygon polygon(polyg_stack[iplane], gmesh);
        // Initialize the basic data of fracture
		// fracture = new DFNFracture(polygon,&dfn,FracLimit::Etruncated);
		fracture = dfn.CreateFracture(polygon,limit_directives[iplane]);
        
		// Find intersected ribs and impose tolerance
		fracture->CreateRibs();
		fracture->SnapIntersections_ribs(tol_dist);
		// Find intersected faces
		fracture->CreateFaces();
		fracture->SnapIntersections_faces(tol_dist,tol_angle);

#ifdef LOG4CXX
        // if(logger->isDebugEnabled()){
        //     std::stringstream sout;
        //     fracture->Print(sout);
        //     LOGPZ_DEBUG(logger,sout.str());
        // }
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
// 		std::ofstream logtest("LOG/dfn.summary.log");
// 		dfn.Print(logtest,argv[1]);
#endif //PZDEBUG
	}
	// Recover Limits
	for(auto frac : dfn.FractureList()){
		frac->RecoverFractureLimits();
	}
	


	// Generate submesh
    // dfn.ExportGMshCAD("dfnExport.geo");
	
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
	std::string default_example("/home/projetos/dfnMesh/examples/time_consuming_test.json");
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


