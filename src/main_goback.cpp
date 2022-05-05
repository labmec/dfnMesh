
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

    #include "dfnrawdata.h"
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

using namespace std;

TPZGeoMesh* ReadInput(int argc, char* argv[], map<int, DFNRawData>& dfnrawdata, REAL &toldist, REAL &tolangle, int& prerefine, TPZManVector<std::map<int,std::string>,4>& dim_physical_tag_and_name);
void CheckForOverlappingFractures(DFNMesh& dfn);

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
int main(int argc, char* argv[]){
#if PZ_LOG
	std::string configpath = PROJECT_ROOT "/src/util/DFNlog4cxx.cfg";
	TPZLogger::InitializePZLOG(configpath);
#endif // PZ_LOG
	TPZTimer time("DFNMesh");
	PrintPreamble();
	TPZGeoMesh *gmesh = nullptr;
	REAL tol_dist = 1.e-4;
	REAL tol_angle = 1.e-6; 
    map<int, DFNRawData> map_dfnrawdata;
	TPZManVector<std::map<int,std::string>,4> dim_physical_tag_and_name;
	int prerefine = 0;
	gmesh = ReadInput(argc,argv,map_dfnrawdata,tol_dist,tol_angle,prerefine,dim_physical_tag_and_name);
	gmsh::initialize();
	DFN::GmshConfig();
	
	
	// ScriptForBug2(gmesh);
	time.start();
    /// Constructor of DFNMesh initializes the skeleton mesh
	DFNMesh dfn(gmesh,dim_physical_tag_and_name,tol_dist,tol_angle,prerefine);
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
    
//	for(int iplane = 0, nfractures = nfracfromfile; iplane < nfractures; iplane++){
    int iplane = 0;
    std::cout << "\n\n==========>  Total number of fractures " << map_dfnrawdata.size() << std::endl;
    map<int, DFNRawData>::iterator it_dfnrawdata = map_dfnrawdata.begin();
    for ( ; it_dfnrawdata != map_dfnrawdata.end() ; it_dfnrawdata++) {
                        
        std::cout << "\n\n\t\t-----------------  Beginning fracture " << iplane << "  -----------------\n" << std::endl;
        if (iplane == 0) it_dfnrawdata = map_dfnrawdata.begin(); // in case it rolls back the first fracture
        DFNRawData &dfnrawdata = it_dfnrawdata->second;
        
        // At the beggining of each step create a backup copy of
        // -> DFNMesh with its respective fGMesh and other attributes
        // -> All the DFNFractures so far with its respective faces and ribs
        // These have to be hard copies? Means, when there is a pointer, we should do a new
        gmeshbackup = new TPZGeoMesh(*dfn.Mesh());
        
        // a polygon represents a set of points in a plane
        // poly_stack[iplane] is a matrix 3xn where n is the number of points 
		DFNPolygon polygon(dfnrawdata.fpolygonmatrices, dfn.Mesh());
        // Initialize the basic data of fracture
        // initialize an empty DFNFracture object
		DFNFracture *fracture = dfn.CreateFracture(polygon,dfnrawdata.flimit_directives,dfnrawdata.fmatid);
        
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
//                iplane--;
                if (iplane > 0) {
                    it_dfnrawdata--;
                }
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
                std::cout << "This is the " << count << " time rollback is applied to fracture index " << it_dfnrawdata->first << std::endl;
                continue;
            }
            else{
                std::cout << "\n*** All polyhedra are consistent ***" << std::endl;
            }
            
			dfn.UpdatePolyhedra();
		}
        
        iplane++;
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
        std::cout << "\n\n\t\t-----------------  Recovering limits fracture " << frac->Index() << "  -----------------\n" << std::endl;
		frac->RecoverFractureLimits();
	}
    
    // Check for overlapping fracures. This can be used for debugging or for simply taking away the fractures
    // that overlap in case we don't have a way to treat this particular case in iMRS yet.
    CheckForOverlappingFractures(dfn);

	// Generate submesh
    dfn.ExportGMshCAD("dfnExport.geo");

	if(map_dfnrawdata.size() == 0){std::cout<<"\nNo fractures were recognized.\n";}
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
TPZGeoMesh* ReadInput(int argc, char* argv[], map<int, DFNRawData>& dfnrawdata, REAL &toldist, REAL &tolangle, int& prerefine, TPZManVector<std::map<int,std::string>,4>& dim_physical_tag_and_name) {
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
	gmesh = SetupExampleFromFile(example,dfnrawdata,mshfile,toldist,tolangle,prerefine,dim_physical_tag_and_name);
	return gmesh;
}

void CheckForOverlappingFractures(DFNMesh& dfn) {
    cout << "\n======================= Searching for overlapping fractures =======================" << endl;
    const int nfrac = dfn.NFractures();
    int noverlapfracs = 0;
    for (int i = 0; i < nfrac; ++i) {
        DFNFracture* fracI = dfn.FractureList()[i];
        std::set<int64_t>& fracIset = fracI->Surface();
        for (int j = i+1; j < nfrac; ++j) {
            DFNFracture* fracJ = dfn.FractureList()[j];
            std::set<int64_t>& fracJset = fracJ->Surface();
            std::vector<int64_t> intersection(fracIset.size()+fracJset.size());
            std::vector<int64_t>::iterator it;
            it = std::set_intersection(fracIset.begin(),
                                       fracIset.end(),
                                       fracJset.begin(),
                                       fracJset.end(),
                                       intersection.begin());
            intersection.resize(it-intersection.begin());
            if (intersection.size() > 0) {
                noverlapfracs++;
                cout << "--------> Found overlapping fractures" << endl;
                cout << "Fracture index " <<  fracI->Index() << " and id " << fracI->MaterialId()
                    << " overlaps with fracture index " << fracJ->Index() << " and id " << fracJ->MaterialId() << endl;
                cout << "Intersection set has size = " << intersection.size() << " and contains:" << endl;
                for (auto& inter : intersection) {
                    cout << inter << " ";
                }
                cout << endl << endl;
            }
        }
    }
    if (noverlapfracs == 0) {
        cout << "----> There are no overlapping fractures!" << endl;
    }
    else {
        cout << "\n----> Found " << noverlapfracs << " overlapping fracture(s)\n" << endl;
    }
}




