
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
	std::string neopzversion = "/commit/eb88bc8"; // https://github.com/labmec/neopz/commit/...
	std::string gmshversion = "4.5.6";
	std::cout<<"\n";
	std::cout<<"\nNeoPZ assumed version: " << neopzversion;
	std::cout<<"\nGMsh assumed version: " << gmshversion << "\n\n";
	std::cout<<"Running...\n";

}

TPZGeoMesh* ReadInput(int argc, char* argv[], TPZStack<TPZFMatrix<REAL>> &polyg_stack, REAL &toldist, REAL &tolangle,TPZManVector<int>& matid,TPZManVector<FracLimit>& limit_directives, int& prerefine);


#if PZ_LOG
	static TPZLogger logger("dfn.mesh");
#endif // PZ_LOG





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
	TPZLogger::InitializePZLOG(configpath);
#endif // PZ_LOG
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
	int prerefine = 0;
	gmesh = ReadInput(argc,argv,polyg_stack,tol_dist,tol_angle,matid,limit_directives,prerefine);
	gmsh::initialize();
	DFN::GmshConfig();
	
	// ScriptForBug2(gmesh);
	time.start();
    /// Constructor of DFNMesh initializes the skeleton mesh
	DFNMesh dfn(gmesh,tol_dist,tol_angle,prerefine);
    dfn.ExportGMshCAD("dfnExportCoarse.geo");

	// std::ofstream out1("graphics/CoarseMesh.vtk");
	// TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out1, true, true);
	// dfn.InheritPolyhedra();
	// dfn.PrintPolyhedra();

    

    // Loop over fractures and refine mesh around them
    TPZGeoMesh *gmeshbackup = nullptr;

    // Uncomment to print polygons
//    std::set<int64_t> polygels;
//    std::ofstream polygfile("allPolygons.vtk");
//    for(int iplane = 0, nfractures = polyg_stack.size(); iplane < nfractures; iplane++){
////        if(iplane != 11) continue;
//        DFNPolygon polygon(polyg_stack[iplane], dfn.Mesh());
//        DFNFracture *fracture = dfn.CreateFracture(polygon,limit_directives[iplane],matid[iplane]);
//        TPZVec<TPZGeoEl*> graphical_elements = polygon.InsertGeomRepresentation(dfn.Mesh(), -(fracture->Index()+1), 1);
//        for(TPZGeoEl* gel : graphical_elements){polygels.insert(gel->Index());}
//        // PlotAllPolygons(dirname + "allPolygons.vtk");
//    }
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh,polygels,polygfile);
        
    int fracrollingback = -1; // used for debugging purposes
    int count = 0;
    
	for(int iplane = 0, nfractures = polyg_stack.size(); iplane < nfractures; iplane++){
        
        std::cout << "\n\n\t\t-----------------  Beginning fracture " << iplane << "  -----------------\n" << std::endl;
        
        // At the beggining of each step create a backup copy of
        // -> DFNMesh with its respective fGMesh and other attributes
        // -> All the DFNFractures so far with its respective faces and ribs
        // These have to be hard copies? Means, when there is a pointer, we should do a new
        gmeshbackup = new TPZGeoMesh(*dfn.Mesh());
        
        // a polygon represents a set of points in a plane
        // poly_stack[iplane] is a matrix 3xn where n is the number of points 
		DFNPolygon polygon(polyg_stack[iplane], dfn.Mesh());
        // Initialize the basic data of fracture
        // initialize an empty DFNFracture object
		DFNFracture *fracture = dfn.CreateFracture(polygon,limit_directives[iplane],matid[iplane]);
        
		// Find intersected ribs and create a corresponding DFNRib object (administered by DFNFracture)
		fracture->CreateRibs();
		// create DFNFace objects - these are faces that have intersected ribs
		fracture->CreateFaces();
        // For the rib object whose intersection point is closer than toldist to an endnode,
        // snap the intersection point to the nearest node
		fracture->SnapIntersections_ribs(tol_dist);
        // this method will snap the ribs with small angles to coincide with
		fracture->SnapIntersections_faces(tol_dist,tol_angle);

		// Search and fix possibly problematic overlaps of fracture surface and existing mesh elements
//		fracture->CheckSnapInducedOverlap();

		LOGPZ_DEBUG(logger, *fracture);
        // we decided that the ribs can be cut. Apply the refinement to the geometric elements
		fracture->RefineRibs();
        // apply the refinement to the faces
		fracture->RefineFaces();
		// dfn.PrintVTKColorful();
	// Mesh fracture surface
        
        //dfn.Polyhedra()[133].PrintVTK();
		if(dfn.Mesh()->Dimension() == 3){
			// divide the fracture in simple geometries using the mesh created in RefineFaces
            TPZStack<int> badVolumes;
			fracture->MeshFractureSurface(badVolumes);
            // now verify angles and create possible list of polyhedra that needs remeshing
            if (badVolumes.size() > 0){
                iplane--;                
                dfn.RollBackLastFracture(gmeshbackup,badVolumes);
                if (iplane != fracrollingback) {
                    fracrollingback = iplane;
                    count = 0;
                }
                else{
                    count++;
                }
                
                std::cout << "\n====> Found volumes with bad angles. Rolled back fracture and refined polyhedra into simplexes..." << std::endl;
                std::cout << badVolumes.size() << " bad volume" << (badVolumes.size()>1?'s':' ') << std::endl;
                std::cout << "This is the " << count << " time rollback is applied to fracture number " << iplane+1 << std::endl;
                continue;
            }
            else{
                std::cout << "\n*** All polyhedra are consistent ***" << std::endl;
            }
            
			dfn.UpdatePolyhedra();
		}
        
        
#ifdef PZDEBUG
        // {
        //     std::ofstream logtest("LOG/dfn.summary.log");
        //     dfn.Print(logtest,argv[1]);
		// 	// dfn.DumpVTK(false,false,"LOG/vtkmesh.vtk");
        // }
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
	// dfn.FractureList()[0]->PlotVTK("LOG/vtkmesh.0.vtk");
	dfn.ExportDetailedGraphics();
	dfn.DumpVTK(true,true);
	dfn.PrintSummary();
	dfn.PrintVTK("skip","LOG/pzmesh.txt");
	std::cout<<"\n ...the end.\n\n";

	gmsh::finalize();
	return 0;
}



// Takes program input and creates a mesh, matrices with the point coordinates, and writes tolerances
TPZGeoMesh* ReadInput(int argc, char* argv[], TPZStack<TPZFMatrix<REAL>> &polyg_stack, REAL &toldist, REAL &tolangle,TPZManVector<int>& matid,TPZManVector<FracLimit>& limit_directives, int& prerefine){
	TPZGeoMesh* gmesh = nullptr;
	std::string default_example("examples/two-hex-and-a-frac.json");
	std::string example = default_example;
	std::string mshfile = "no-msh-file";
	for(int iarg=1; iarg < argc; ++iarg){
		std::string aux = argv[iarg];
		try{
			if(argv[iarg][0] != '-'){example = argv[iarg];}
			else if(aux == "-m"){mshfile = argv[++iarg];}
			else if(aux == "-f"){example = argv[++iarg];}
			else if(aux == "-td"){toldist = std::stod(argv[++iarg]);}
			else if(aux == "-ta"){tolangle = std::stod(argv[++iarg]);}
			else if(aux == "-tc"){tolangle = std::acos(std::stod(argv[++iarg]));}
			else if(aux == "-r"){prerefine = std::stoi(argv[++iarg]);}
			else{
				throw std::bad_exception();
			}
		}catch(...){
			PZError << "\nUnrecognized argument passed:\n\t\""<< argv[iarg] << "\"\n" << std::endl; 
			DebugStop();
		}
	}
	std::cout<<"input file: "<<example<<"\n";
	gmesh = SetupExampleFromFile(example,polyg_stack,mshfile,toldist,tolangle,matid,limit_directives,prerefine);
	return gmesh;
}






