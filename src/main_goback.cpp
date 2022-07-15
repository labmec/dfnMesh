
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

using namespace std;
namespace fs = std::filesystem;

void ParseInputArguments(int argc, char* argv[], std::string &inputfiledir, std::string &meshfile, REAL &toldist, REAL &tolangle, int& prerefine);

TPZGeoMesh* ReadInput(std::string &pathToJSon, std::string &MeshFileName, map<int, DFNRawData>& dfnrawdata, REAL &toldist, REAL &tolangle, int& prerefine, TPZManVector<std::map<int,std::string>,4>& dim_physical_tag_and_name);
void CheckForOverlappingFractures(DFNMesh& dfn);
void CreateOutputFolders(std::string& outputFolder);
void CopyJsonAndMeshToOutputFolder(std::string& pathToJson,std::string& outputFolder,std::string& meshFile);

#if PZ_LOG
	static TPZLogger logger("dfn.mesh");
#endif // PZ_LOG


void getOutputFileNames(std::string inputfiledir , std::string& outputFolder, std::string& coarseOutputName,
                        std::string& fineOutputName, std::string& pathToJson, std::string& meshFile);
bool fileExists(const fs::path& p, fs::file_status s = fs::file_status{});


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
    std::string inputfiledir;

    std::string coarseOutputName,fineOutputName,outputFolder,pathToJson,meshFile;
    ParseInputArguments(argc,  argv, inputfiledir, meshFile, tol_dist, tol_angle, prerefine);
    // if the only argument is the json filename or directory
    // method that parses the input argument and identifies the json file
    // outputFolder :
    getOutputFileNames(inputfiledir,outputFolder,coarseOutputName,fineOutputName,pathToJson,meshFile);
    
    gmesh = ReadInput(pathToJson, meshFile ,map_dfnrawdata,tol_dist,tol_angle,prerefine,dim_physical_tag_and_name);
	gmsh::initialize();
	DFN::GmshConfig();
    
    // Creating output folder
    CreateOutputFolders(outputFolder);
    CopyJsonAndMeshToOutputFolder(pathToJson,outputFolder,meshFile);
	
//    dfn.ExportGMshCAD(coarseOutputName);
	// ScriptForBug2(gmesh);
	time.start();
    /// Constructor of DFNMesh initializes the skeleton mesh
	DFNMesh dfn(gmesh,dim_physical_tag_and_name,tol_dist,tol_angle,prerefine);
    

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
    dfn.ExportGMshCAD(fineOutputName);

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

void ParseInputArguments(int argc, char* argv[], std::string &inputfiledir, std::string &meshfile, REAL &toldist, REAL &tolangle, int& prerefine)
{
    std::string default_example("examples/two-hex-and-a-frac.json");
    inputfiledir = default_example;
    meshfile = "no-msh-file";
    for(int iarg=1; iarg < argc; ++iarg){
        std::string aux = argv[iarg];
        try{
            if(argv[iarg][0] != '-'){inputfiledir = argv[iarg];}
            else if(aux == "-m"){meshfile = argv[++iarg];}
            else if(aux == "-f"){inputfiledir = argv[++iarg];}
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
}

// Takes program input and creates a mesh, matrices with the point coordinates, and writes tolerances
TPZGeoMesh* ReadInput(std::string &pathToJSon, std::string &mshfile , map<int, DFNRawData>& dfnrawdata, REAL &toldist, REAL &tolangle, int& prerefine, TPZManVector<std::map<int,std::string>,4>& dim_physical_tag_and_name) {
	TPZGeoMesh* gmesh = nullptr;
	std::cout<<"input file: "<<pathToJSon<<"\n";
	gmesh = SetupExampleFromFile(pathToJSon,dfnrawdata,mshfile,toldist,tolangle,prerefine,dim_physical_tag_and_name);
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

void getOutputFileNames(std::string ArgName, std::string& outputFolder, std::string& coarseOutputName,
                        std::string& fineOutputName, std::string& JsonFilename, std::string& meshFile){
    
    // Checking if a folder or a json was provided. If folder, append the json file in the folder to the string ArgName
    int njson = 0;
    std::string basemeshpath(INPUTMESHES);
    std::string outFileName;
    std::string dirname;
    if(ArgName.find(".json") == std::string::npos){
        dirname = basemeshpath + "/" + ArgName;
        if(!fileExists(dirname)){
            cout << "\n\n=====> ERROR! Folder " << dirname << " does not exist" << endl;
            DebugStop();
        }
        for (auto const& dir_entry : std::filesystem::directory_iterator{dirname}){
            std::string filename = dir_entry.path().string();
//            std::cout << "dir entry:" << dir_entry.path().string() << endl;
            auto pos = filename.find(".json");
            if(pos == std::string::npos) continue;
//            std::cout << "position : " << pos << std::endl;
//            std::cout << "substr :" << filename.substr(pos) << std::endl;
            if(filename.substr(pos) == ".json"){
                if(njson == 0){
                    outFileName = filename.substr(filename.find_last_of("/") + 1);
                    auto example_pos = dirname.rfind("/examples/");
                    auto dirname_noexample = dirname.substr(example_pos+10,std::string::npos);
                    JsonFilename = "../examples/" + dirname_noexample + "/" + outFileName;
                    outFileName = outFileName.substr(0,outFileName.find("."));
                    njson++;
                }
                else {
                    cout << "\n\n=====> ERROR! There are two json files in the provided input folder" << endl;
                    DebugStop();
                }
            }
        }
    }else{
        JsonFilename = ArgName;
//        ArgName = ArgName.substr(0, ArgName.find_last_of("."));
        auto lastslashpos = ArgName.find_last_of("/");
        dirname = ArgName.substr(0,lastslashpos);
        // everything after the last "/"
        outFileName = ArgName.substr(lastslashpos+1,ArgName.length()-lastslashpos);
        outFileName = outFileName.substr(0, outFileName.find_last_of("."));
        dirname = dirname + "/" + outFileName;
    }
    auto example_pos = dirname.rfind("/examples/");
//    std::cout << "dirname size " << dirname.size() << std::endl;
    auto dirname_noexample = dirname.substr(example_pos+10,std::string::npos);
    outputFolder = "Outputs/" + dirname_noexample + "/";
    fineOutputName = outputFolder + outFileName + "_fine.geo";
    coarseOutputName = outputFolder + outFileName + "_coarse.geo";
    
    // Getting mesh if available
    nlohmann::json input;
    auto filename = basemeshpath+"/"+JsonFilename;
    std::cout << "opening file :" << filename << std::endl;
    std::ifstream file(basemeshpath+"/"+JsonFilename);
    if(!file) DebugStop();
    input = nlohmann::json::parse(file,nullptr,true,true); // to ignore comments in json file
    meshFile.clear();
    if(input.find("Mesh") != input.end()){
        meshFile = (std::string)input["Mesh"];
        auto lastslash = JsonFilename.find_last_of("/");
        auto meshdirname = JsonFilename.substr(0,lastslash+1);
//        meshFile = meshFile.substr(meshFile.find("examples/") + 9,meshFile.length());
        meshFile = basemeshpath + "/" +meshdirname + meshFile;
    }
    std::cout << "outputFolder : " << outputFolder << std::endl;
    std::cout << "coarseOutputName : " << coarseOutputName << std::endl;
    std::cout << "fineOutputName : " << fineOutputName << std::endl;
    std::cout << "JsonFilename : " << JsonFilename << std::endl;
    std::cout << "meshFile : " << meshFile << std::endl;
}

void CreateOutputFolders(std::string& outputFolder) {

    std::string folders = outputFolder;
    char c = '/';
    std::string folderToCreate = "";
    int nfolders = 0;
    while (folders.find("/") != std::string::npos) {
        if(nfolders == 0) folderToCreate = folders.substr(0,folders.find("/"));
        else folderToCreate = folderToCreate + "/" + folders.substr(0,folders.find("/"));
        folders = folders.substr(folders.find("/")+1);
        if(!fileExists(folderToCreate)){
            if (!fs::create_directory(folderToCreate))
                DebugStop();
            else
                cout << "Directory created with name " << folderToCreate << endl;
        }
        nfolders++;
    }
}

bool fileExists(const fs::path& p, fs::file_status s) {
    if(fs::status_known(s) ? fs::exists(s) : fs::exists(p))
        return true;
    else
        return false;
}

void CopyJsonAndMeshToOutputFolder(std::string& pathToJson,std::string& outputFolder,std::string& meshFile){
    std::string basemeshpath(INPUTMESHES);
    std::string outFileName = outputFolder.substr(outputFolder.find_last_of("/",outputFolder.length()-2)+1);
    auto slashpos = meshFile.find_last_of("/");
    auto dirname = meshFile.substr(0,slashpos+1);
    auto extensionpos = meshFile.rfind(".msh");
    auto rootname = meshFile.substr(slashpos+1,extensionpos-slashpos-1);
    
    auto JsonSlash = pathToJson.find_last_of("/");
    auto JsonExt = pathToJson.rfind(".json");
    auto JsonRootname = pathToJson.substr(JsonSlash+1,JsonExt-JsonSlash-1);
    outFileName = meshFile.substr(slashpos+1,std::string::npos);
    auto geoinputname = dirname + rootname + ".geo";
    auto geooutputname = outputFolder + rootname + ".geo";
    auto geooutputnameJson = outputFolder + JsonRootname + "_coarse.geo";
    auto meshoutputnameJson = outputFolder + JsonRootname + "_coarse.msh";
    
    fs::copy(basemeshpath + "/" + pathToJson, outputFolder, fs::copy_options::update_existing);
    if(meshFile.size() && meshFile != "none"){
        fs::copy(meshFile, outputFolder + outFileName, fs::copy_options::update_existing);
    }
    fs::copy(geoinputname,geooutputname,fs::copy_options::update_existing);
    fs::copy(geoinputname,geooutputnameJson,fs::copy_options::update_existing);
    fs::copy(meshFile,meshoutputnameJson,fs::copy_options::update_existing);
}
