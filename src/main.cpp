
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
	#include <filesystem>
//includes





/**
 * @brief information and assumptions
*/
void PrintPreamble(){
	std::string neopzversion = "/commit/bb94d514"; // https://github.com/labmec/neopz/commit/...
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
	std::string configpath = PROJECT_ROOT "/src/util/DFNlog4cxx.cfg";
	log4cxx::PropertyConfigurator::configure(configpath);
#endif // LOG4CXX
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
		fracture->CreateRibs();
        // For the rib object whose intersection point is closer than toldist to an endnode,
        // snap the intersection point to the nearest node
		fracture->SnapIntersections_ribs(tol_dist);
		// create DFNFace objects - these are faces that have intersected ribs (always 2)
		fracture->CreateFaces();
        // this method will snap the ribs with small angles to coincide with
		fracture->SnapIntersections_faces(tol_dist,tol_angle);

		fracture->Handle_SnapInducedOverlap();
		// std::set<int64_t> SnapRibs;
		// fracture->IdentifySnapRibs(SnapRibs);
		// for(int64_t index : SnapRibs){
		// 	gmesh->Element(index)->SetMaterialId(20);
		// }

#ifdef LOG4CXX
        if(logger->isDebugEnabled()){
            std::stringstream sout;
            fracture->Print(sout);
            LOGPZ_DEBUG(logger,sout.str());
        }
#endif
        // we decided that the ribs can be cut. Apply the refinement to the geometric elements
		fracture->RefineRibs();
        // apply the refinement to the faces
        // @pedro - shouldnt we verify if the snap created a fracture line across a face
        // maybe we need the notion of snap-ribs - ribs that belong to the fracture but
        // are the result of a snap operation
        // set of snap-ribs is suggested data structure of DFNFracture
        // problems occurr when snap-ribs force the division of a face that was not intersected
        // by the fracture plane "in the first place"
        // In order to "know" if a fracture line will cross a face, we need to know the
        // for each set of coplanar faces belonging to a same polyhedra a set of neighbouring snap-ribs
        
        // For each fracture, we need to identify the set of polyhedra - 
        // For each intersected polyhedra, we need to verify if the snap-ribs will induce a refinement
        // of one of its facets
        // this method must be implemented in DFNFracture as it affects objects out of the scope
        // of DFNFace and DFNRib
        
        
		fracture->RefineFaces();
		// dfn.PrintVTKColorful();
	// Mesh fracture surface
		if(gmesh->Dimension() == 3){
			// dfn.InheritPolyhedra();
            // divide the fracture in simple geometries using the mesh created in RefineFaces
			fracture->MeshFractureSurface();
			// dfn.DumpVTK();
            // this is where the code can crash : a face that should be internal can be aligned
            // with an existing face
			dfn.UpdatePolyhedra();
		}
#ifdef PZDEBUG
        {
            std::ofstream logtest("LOG/dfn.summary.txt");
            dfn.Print(logtest,argv[1]);
			dfn.DumpVTK(false,true,"LOG/gmesh.vtk");
        }
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
	dfn.DumpVTK(true);
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






