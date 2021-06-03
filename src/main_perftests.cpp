
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
 * @brief information and assumptions
*/
void PrintPreamble(){
	std::string neopzversion = "/commit/eb88bc8"; // https://github.com/labmec/neopz/commit/...
	std::string neopzbranch = "MHM_HybridMixed"; // https://github.com/labmec/neopz/commit/...
	std::string gmshversion = "4.7.0";
	std::cout<<"\nDFNMesh PERFORMANCE TESTS\n";
	std::cout<<"\nNeoPZ assumed version: " << neopzversion << "\tBranch: " << neopzbranch;
	std::cout<<"\nGMsh assumed version: " << gmshversion << "\n\n";
	std::cout<<"Running...\n";

}



#if PZ_LOG
static TPZLogger logger("dfn.mesh");
#endif




void Test1();
void Test2();


//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

int main(int argc, char* argv[]){
#if PZ_LOG
	std::string configpath = PROJECT_ROOT "/src/util/DFNlog4cxx.cfg";
	log4cxx::PropertyConfigurator::configure(configpath);
#endif // PZ_LOG
	PrintPreamble();
	gmsh::initialize();
	Test1();

	std::cout<<"\n ...the end.\n\n";
	gmsh::finalize();
	return 0;
}







namespace test1{

	TPZGeoMesh* CreateMesh(int reflevel = 1){
		MMeshType eltype = MMeshType::ENoType;
		TPZManVector<int,3> nels(3,reflevel);
		TPZManVector<REAL,3> x0(3,-1.), x1(3,1.);


		enum EConfig {Ecube = 1, Etetrahedron = 2};

		EConfig config = Ecube;

		switch(config){
			case Ecube:{
				eltype = MMeshType::EHexahedral;
				x0 = {-1,-1,-1};
				x1 = { 1, 1, 1};
				break;
			}
			case Etetrahedron: DebugStop(); break;
			default : DebugStop();
		}

		std::cout << " -Creating Mesh\r" << std::flush;
		// Creating the Geo mesh
		TPZGeoMesh *gmesh = new TPZGeoMesh;
		if(!nels[2]){ // 2D mesh
			DebugStop();
		}else{ // 3D mesh
			TPZGenGrid3D gengrid(x0,x1,nels,eltype);
			gmesh = gengrid.BuildVolumetricElements(1);
		}
		std::cout << "               \r" << std::flush;
		return gmesh;
	}

}
/// @brief Constant fracture with a gradually refining coarse mesh
void Test1(){
	TPZFMatrix<REAL> frac_corners = {{ 1.05, 1.15,-1.20},
									 {-1.10, 1.05, 1.05},
									 { 1.05,-1.20, 1.15}};
	std::ofstream timelog("LOG/TEST1.txt");
	TPZTimer time("test1_timer");
	time.start();
	for(int reflevel=64; reflevel <= 128; reflevel *= 2){
		std::cout << "\nreflvl = " << reflevel << std::endl;
		std::unique_ptr<TPZGeoMesh> gmesh(test1::CreateMesh(reflevel));
		time.reset();
		DFNMesh dfn(gmesh.get(),2.,3.);

		DFNPolygon polygon(frac_corners,gmesh.get());
		std::unique_ptr<DFNFracture> fracture(dfn.CreateFracture(polygon,FracLimit::Etruncated));

		fracture->CreateRibs();
		fracture->SnapIntersections_ribs(dfn.TolDist());
		fracture->CreateFaces();
		fracture->SnapIntersections_faces(dfn.TolDist(),dfn.TolAngle());
		fracture->RefineRibs();
		fracture->RefineFaces();
		fracture->MeshFractureSurface();
		dfn.UpdatePolyhedra();
		fracture->RecoverFractureLimits();
		
		time.stop();
		timelog << "reflvl = " << std::setw(4) << std::right << reflevel;
		timelog << " [time] " << time << std::endl;
	}

}