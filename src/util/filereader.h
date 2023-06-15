#include "TPZGmshReader.h"
#include "pzlog.h"

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

#include "json.hpp"

#include "dfnrawdata.h"
#include "dfn_config.h"
#include <iomanip> // std::setprecision.
#include <libInterpolate/Interpolate.hpp>

/**
 * @brief Define which example to run. See example file sintax in function definition
 * @param filename: path to the json file that defines the example
 * @param dfnrawdata : data containing fracture geometries - output read from json file
 * @param mshfile: [optional] path to .msh file (if applicable) ??? read from json file
 * @param toldist: distance tolerance - read from json file
 * @param tolangle : angle tolerance - read from jscon file
 * @param dim_physical_tag_and_name : name of the physical tags read from gmsh - read from gmsh file
 * @returns pointer to geometric mesh created/read
*/
TPZGeoMesh* SetupExampleFromFile(std::string filename, std::map<int, DFNRawData>& dfnrawdata, std::string mshfile, REAL& toldist, REAL& tolangle, TPZManVector<std::map<int,std::string>,4>& dim_physical_tag_and_name, REAL& meshEdgesBaseSize);
void ModifyTopeAndBase2(TPZGeoMesh * gmesh ,int nlayers);
void ReadData(std::string name, bool print_table_Q, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);
TPZGeoMesh* GenerateUnisimMesh(int nlayers);
void AllForceCoplanarity(TPZGeoMesh *gmesh);
bool ForceCoplanarity(TPZVec<REAL> &c1, TPZVec<REAL> &c2, TPZVec<REAL> &c3, TPZVec<REAL> &c4);
double determinante(TPZVec<TPZVec<REAL>> matriz);
bool IsCoplanar(TPZVec<REAL> &c1, TPZVec<REAL> &c2, TPZVec<REAL> &c3, TPZVec<REAL> &c4);
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
                std::map<int, DFNRawData>& dfnrawdata,
				std::string 				& mshfile,
				TPZManVector<REAL,3> 		& x0,
				TPZManVector<REAL,3> 		& xf,
				TPZManVector<int,3> 		& nels,
				MMeshType					& eltype,
				TPZManVector<REAL,2>		& tol,
                REAL                        & meshEdgesBaseSize
				);


// filename : name of file that can be opened
// map_dfnrawdata : polygon data and properties of each fracture
// meshfile : name of gmsh file that will be used for intersection - as read in the json file
// x0 : no idea - coordinates of what?
// xf : no idea - coordinates of what?
// nels : number of elements ??
// tol : distance and angle tolerance
// prefine : number of refinements??
void ReadFile(	const std::string			& filename,
				TPZStack<TPZFMatrix<REAL>>	& polygonmatrices, 
				std::string 				& mshfile,
				TPZManVector<REAL,3> 		& x0,
				TPZManVector<REAL,3> 		& xf,
				TPZManVector<int,3> 		& nels,
				MMeshType					& eltype,
				TPZManVector<REAL,2>		& tol,
				TPZManVector<int>			& matid,
				TPZManVector<FracLimit>		& limit_directives,
				int							& prerefine
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
	using namespace std;

	std::string line, word;
	eltype = MMeshType::ENoType;
	// const string Domain = "Domain";
	// const string Mesh("Mesh"), Fractures("NumberOfFractures");
	int i, j, nfractures, ifrac=0;
	REAL Lx, Ly, Lz;
	TPZManVector<REAL,3> L(3,-1.);
	x0.Fill(0.);
	// Read file
	std::ifstream file(filename);
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
		else if(word == "PreRefine" || word == "prerefine"){
			while (getline(ss, word, ' ') && word.length() == 0){/*void*/};
			getline(ss, word, ' ');
			prerefine = std::stoi(word);
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

// filename : name of file that can be opened
// map_dfnrawdata : polygon data and properties of each fracture
// meshfile : name of gmsh file that will be used for intersection - as read in the json file
// x0 : no idea - coordinates of what?
// xf : no idea - coordinates of what?
// nels : number of elements ??
// tol : distance and angle tolerance
// prefine : number of refinements??
void ReadFileJSON(const std::string			& filename,
                std::map<int, DFNRawData>&  map_dfnrawdata,
				std::string 				& mshfile,
				TPZManVector<REAL,3> 		& x0,
				TPZManVector<REAL,3> 		& xf,
				TPZManVector<int,3> 		& nels,
				MMeshType					& eltype,
				TPZManVector<REAL,2>		& tol,
				int                         & prerefine,
                REAL                        & meshEdgesBaseSize
				)
{
		/*_______________________________________________________________
						JSON OPTIONS 
	{												// <Don't forget opening brace>
		"$schema": "path-schema.json"	<string>	// I have set up a json schema for the user to get autocomplete. It's in ${workspaceFolder}/examples/dfn_schema.json. Just add the $schema line and be happy. Make sure you've got your relative path right.
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
		"PreRefine": value,				<int>		// Number of uniform refinements to apply to the coarse mesh before introducing fractures. Defaults to zero.
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
        auto lastslash = filename.find_last_of("/");
        auto dirname = filename.substr(0,lastslash+1);
        mshfile = dirname+mshfile;
	}
	else{
		PZError << "\nMissing coarse mesh definition.\n"; DebugStop();
	}

	// Tolerances
    tol.resize(2);
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
	
	// An option to apply uniform refinement in the mesh before fractures
	if(input.find("PreRefine") != input.end())
		{prerefine = input["PreRefine"];}
	else if(input.find("prerefine") != input.end())
		{prerefine = input["prerefine"];}
	
	//Init matid for fractures
	const int fracInitMatId = input["FractureInitMatId"];
	int actualMatId = fracInitMatId;
    
    // Mesh size
    if(input.find("MeshEdgesBaseSize") != input.end())
        {meshEdgesBaseSize = (REAL)input["MeshEdgesBaseSize"];}
	
	// Fractures
	int nfractures = input["Fractures"].size();
	for(auto& fracture : input["Fractures"]){
		int i = fracture["Index"];
        
        if(map_dfnrawdata.find(i) != map_dfnrawdata.end()){ // repeated indexes!
            DebugStop();
        }
        else{
            DFNRawData mydata;
            if(fracture.find("MaterialID") != fracture.end()){
				DebugStop(); // Not using this structure anymore
                mydata.fmatid = (int)fracture["MaterialID"];
            }else if(fracture.find("MatID") != fracture.end()){
				DebugStop(); // Not using this structure anymore
                mydata.fmatid = (int)fracture["MatID"];
            }else{
				mydata.fmatid = actualMatId;
				actualMatId += 2;
			}
            if(fracture.find("Limit") != fracture.end()){
                mydata.flimit_directives = DFN::StringToFracLimit((std::string)fracture["Limit"]);
            }
            if(fracture.find("NRefFracBorder") != fracture.end()){
                mydata.fnrefborder = (int)fracture["NRefFracBorder"];
            }
            if(fracture.find("SizeEdgesTouchFracBorder") != fracture.end()){
                if(meshEdgesBaseSize < 0) DebugStop();
                mydata.fSizeOfEdgesTouchFracBorder = (REAL)fracture["SizeEdgesTouchFracBorder"];
            }
            int npoints = fracture["Nodes"].size();
            mydata.fpolygonmatrices.Resize(3,npoints);
            for(int j=0; j<npoints; j++){
                for(int k=0; k<3; k++){
                    mydata.fpolygonmatrices(k,j) = (REAL)fracture["Nodes"][j][k];
                }
            }
            map_dfnrawdata[i] = mydata; // this makes a copy.
        }
        
	}

}

// filename : a name that with INPUTMESHES directory leads to a json file
// dnfrawdata : no idea - output of this function
// mshfile : I dont know. It is not used - overwritten by the .json file
// toldist : distance tolerance input/output should be in the json file
TPZGeoMesh* SetupExampleFromFile(std::string filename, std::map<int, DFNRawData>& dfnrawdata, std::string mshfile, REAL& toldist, REAL& tolangle, int& prerefine, TPZManVector<std::map<int,std::string>,4>& dim_physical_tag_and_name, REAL& meshEdgesBaseSize){

	std::string basemeshpath(INPUTMESHES);
	filename = basemeshpath + "/" + filename;

	MMeshType eltype = MMeshType::ENoType;
	TPZStack<Matrix> planestack;
	TPZManVector<REAL, 3> x0(3, 0.), x1(3, 0.);
	TPZManVector<int,3> nels(3,0);
	TPZManVector<REAL,2> tol = {toldist,tolangle};

	// Choose extension and read file
	std::string extension;
	try{extension = filename.substr(filename.find_last_of('.'));}
	catch(std::out_of_range){extension = '?';}
    if(extension == ".txt"){
        DebugStop(); // DEPRECATED!
//        {ReadFile(filename,polyg_stack,mshfile,x0,x1,nels,eltype,tol,matid,limit_directives,prerefine);}
    }
    else if(extension == ".json" || extension == ".jsonc") {
		ReadFileJSON(filename,dfnrawdata,mshfile,x0,x1,nels,eltype,tol,prerefine,meshEdgesBaseSize);
    }
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
        std::string fullfilename = mshfile;
//        mshfile    std::string    "/Users/jose/Documents/GitHub/dfnMesh/examples/../examples/ResultsJose/UNISIM_Test/fl_case1.msh"
//#ifdef MACOSX
//        fullfilename = "../" + mshfile;
//#endif
//        fullfilename = basemeshpath + "/" + fullfilename;
//		gmesh = reader.GeometricGmshMesh(fullfilename, gmesh);
        
        //
        
        gmesh=GenerateUnisimMesh(2);
//        TPZPersistenceManager::OpenRead("/Users/jose/Documents/GitHub/dfnMesh/examples/ResultsJose/UNISIM_Test/test_coarse.txt");
//        TPZSavable *restore = TPZPersistenceManager::ReadFromFile();
//        gmesh = dynamic_cast<TPZGeoMesh *>(restore);
        
        //
        
		dim_physical_tag_and_name = reader.GetDimPhysicalTagName();
#ifdef PZDEBUG
		std::cout << "------------------- Materials from supplied coarse mesh ---------------" << std::endl;
		for(int i = 0 ; i < reader.GetDimPhysicalTagName().size() ; i++){
			std::cout << "=========> dimension " << i << std::endl;
			std::map<int,std::string>& stringPhysTag = reader.GetDimPhysicalTagName()[i];
			for(auto el : stringPhysTag){
				std::cout << "phystag: " << el.first;
				std::cout << "\t\tname: " << el.second << std::endl;
			}
		}
#endif
	}
	
	

	std::ofstream out1("graphics/CoarseMesh.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out1, true, true);
	return gmesh;
}
TPZGeoMesh * GenerateUnisimMesh(int nlayers){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    /*
     2 4 "inlet"
     2 5 "outlet"
     2 6 "noflux"
     3 3 "k33"
     3 10 "k31"
     */
    dim_name_and_physical_tagCoarse[2]["k33"] = 1;
    dim_name_and_physical_tagCoarse[2]["k31"] = 2;
    dim_name_and_physical_tagCoarse[1]["inlet"] = 3;
    dim_name_and_physical_tagCoarse[1]["outlet"] = 4;
    dim_name_and_physical_tagCoarse[1]["noflux"] = 5;
    
    // Creating gmsh reader
    TPZGmshReader  GeometryFine;
    TPZGeoMesh *gmesh2D;
    REAL l = 1.0;
    GeometryFine.SetCharacteristiclength(l);
   // std::string filename("/Users/jose/Documents/GitHub/dfnMesh/examples/ResultsJose/UNISIM_Test/unisim_2D.msh");
    
    std::string filename("/home/jose/GitHub/dfnMesh/examples/ResultsJose/UNISIM_Test/unisim_2D.msh");
    // Reading mesh
    GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagCoarse);
    gmesh2D = GeometryFine.GeometricGmshMesh(filename,nullptr,false);
    gmesh2D->BuildConnectivity();
    std::cout<<"Dim: "<<gmesh2D->Dimension()<<std::endl;
    std::cout<<"Nels: "<<gmesh2D->NElements()<<std::endl;
    
    std::ofstream file2("COARSEfromDFN2D.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh2D, file2);
    
    
    
   
    int topID= 5;
    int baseID = 5;
    int w=200;
    
    TPZGeoMesh * returnedMesh = nullptr;
    
    
    TPZExtendGridDimension extend(gmesh2D, w);
    extend.SetElType(1);
    returnedMesh = extend.ExtendedMesh( nlayers,topID,baseID);
    
    ModifyTopeAndBase2(returnedMesh ,nlayers);
    std::ofstream file("COARSEfromDFNaNTES.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(returnedMesh, file);
    
//    AllForceCoplanarity(returnedMesh);
//    AllForceCoplanarity(returnedMesh);
//    std::ofstream file20("COARSEfromDFN.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(returnedMesh, file20);
//    coarse3D=returnedMesh;
    
//    std::ofstream file("COARSEfromDFN.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(returnedMesh, file);
    return returnedMesh;
    
}
void ModifyTopeAndBase2(TPZGeoMesh * gmesh ,int nlayers){
//    std::string filename2 = "Reservoir/base_unisimMOD.txt";
    
//     std::string filename1 = "topeMOD.txt";
//    std::string filename2 = "baseMOD.txt";
//
//    std::string filename1 = "baseMOD.txt";
//   std::string filename2 = "topeMOD.txt";
//
//    std::vector<double> x, y, z, x1,y1,z1;
//    ReadData(filename1, true, x, y, z);
//    ReadData(filename2, true, x1, y1, z1);
//
//    //
//    double sum=0.0;
//    for (auto val:z) {
//        sum += val;
//    }
//    double val_base= sum / z1.size();
////    double val_tope= sum / z.size();
//    sum=0.0;
//    for (auto val:z1) {
//        sum += val;
//    }
////    double val_base= sum / z1.size();
//    double val_tope= sum / z.size();
//
//    _2D::ThinPlateSplineInterpolator <double> interpTope;
//    _2D::ThinPlateSplineInterpolator <double> interpBase;
//    int nCoordinates = gmesh->NodeVec().NElements();
//
//    interpTope.setData(x,y,z);
//    interpBase.setData(x1,y1,z1);
    
//    interpBase.setData(x,y,z);
//    interpTope.setData(x1,y1,z1);
    
//    double delt = 200/nlayers;
//    double maxval =200*nlayers;
    
//    for(int inode=0; inode<nCoordinates; inode++){
//        TPZGeoNode node = gmesh->NodeVec()[inode];
//        TPZVec<REAL> co(3);
//        node.GetCoordinates(co);
//        double topeinterpol_val =interpTope(co[0],co[1]);
//        double baseinterpol_val = interpBase(co[0],co[1]);
//        double dist = std::abs(topeinterpol_val - baseinterpol_val);
////        topeinterpol_val = baseinterpol_val - dist;
////        baseinterpol_val = topeinterpol_val + dist;
////        double topeinterpol_val =interpBase(co[0],co[1]);
////        double baseinterpol_val = interpTope(co[0],co[1]);
//        double zcord = co[2];
//        if (topeinterpol_val==0) {
//            topeinterpol_val = val_tope;
//            if (co[0]>1000.00) {
//                topeinterpol_val -= 120;
//            }
//        }
//        if (baseinterpol_val==0) {
//            baseinterpol_val = val_base;
//            if (co[0]>1000.00) {
//                baseinterpol_val = val_base-80;
//            }
//
//        }
//        if(zcord==0.0){
////            co[2] = baseinterpol_val;
//            co[2] = topeinterpol_val;
//            gmesh->NodeVec()[inode].SetCoord(co);
//        }
//        else if(zcord==maxval){
//            co[2] = baseinterpol_val;
////            co[2] = topeinterpol_val;
//            gmesh->NodeVec()[inode].SetCoord(co);
//        }
//        else{
//            double nivel = zcord/200;
//            std::cout<<" funciona"<<std::endl;
//            co[2] =  topeinterpol_val + (nivel)*(baseinterpol_val - topeinterpol_val)/(nlayers);
////            co[2] =  baseinterpol_val + (nivel)*(topeinterpol_val - baseinterpol_val)/(nlayers);
//            gmesh->NodeVec()[inode].SetCoord(co);
//        }
//
//    }
        std::string filename1 = "topeMOD.txt";
        std::string filename2 = "baseMOD.txt";
        std::vector<double> x, y, z, x1,y1,z1;
        ReadData(filename1, true, x, y, z);
        ReadData(filename2, true, x1, y1, z1);
    
        _2D::ThinPlateSplineInterpolator <double> interpTope;
        _2D::ThinPlateSplineInterpolator <double> interpBase;
    
        interpTope.setData(x,y,z);
        interpBase.setData(x1,y1,z1);
    
        int nCoordinates = gmesh->NodeVec().NElements();
        double sum=0.0;
        for (auto val:z) {
            sum += val;
        }
        double val_tope= sum / z.size();
        sum=0.0;
        for (auto val:z1) {
            sum += val;
        }
        double val_base= sum / z1.size();
    //    val_base = 1000;
    //    val_tope = 5000;
    //
    //    val_tope = 3000;
    //    val_base = 3000;
        int npointsPerLayer = nCoordinates/(nlayers+1);
        double valinter=0.0;
        for (int ilay = 1; ilay <= nlayers+1; ilay++) {
            for (int ipoint = (ilay-1)*npointsPerLayer; ipoint<(ilay)*npointsPerLayer; ipoint++) {
                TPZGeoNode node = gmesh->NodeVec()[ipoint];
                
                TPZVec<REAL> co(3);
                node.GetCoordinates(co);
                double topeinterpol =interpTope(co[0],co[1]);
                double baseinterpol = interpBase(co[0],co[1]);
                if (topeinterpol==0) {
                    topeinterpol = val_tope;
                    if (co[0]>1000.00) {
                        topeinterpol -= 120;
                    }
                }
                if (baseinterpol==0) {
    
                    baseinterpol = val_base;
                    if (co[0]>1000.00) {
                       baseinterpol = val_base-80;
                    }
    
                }
    
                if (ilay==1) {
                    valinter=topeinterpol;
    //                valinter = 3500;
                    co[2]=valinter;
                    gmesh->NodeVec()[ipoint].SetCoord(co);
                }
                if (ilay==nlayers+1) {
                    valinter = baseinterpol;
    //                valinter = 2850;
                    co[2]=valinter;
                    gmesh->NodeVec()[ipoint].SetCoord(co);
                }
                if (ilay>1   && ilay < nlayers+1) {
                    valinter = topeinterpol + (ilay-1)*(baseinterpol - topeinterpol)/(nlayers);
                    co[2]=valinter;
                    gmesh->NodeVec()[ipoint].SetCoord(co);
                }
            }
        }
}
void ReadData(std::string name, bool print_table_Q, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z){
    
    bool modpoints = true;
    std::ifstream file;
    std::string basemeshpath("/home/jose/GitHub/iMRS/iMRS/FracMeshes/dfnimrs/unisim_meshes/Reservoir_props/");
    basemeshpath = basemeshpath  + name;
    file.open(basemeshpath);
    int i=1;
    
    
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::istringstream issText(line);
        char l = line[0];
        if(l != '/'){
            i=i+1;
            int val = i%15;
            if(val ==0){
                double a, b, c;
                if(iss >> a >> b >> c) ;
                if (modpoints) {
                    x.push_back(a - 350808.47);
                    y.push_back(b - 7.51376238e6);
                    z.push_back(c);
                }
                else{
                x.push_back(a);
                y.push_back(b);
                z.push_back(c);
                }
            };
        };
    };
    
    if(x.size() == 0){
        std::cout<<"No data read."<<std::endl;
        
        DebugStop();
    }
    if(print_table_Q){
        std::cout<<"*************************"<<std::endl;
        std::cout<<"Reading file... ok!"<<std::endl;
        std::cout<<"*************************"<<std::endl;
        std::cout<<x.size()<<std::endl;
        std::cout<<y.size()<<std::endl;
        std::cout<<z.size()<<std::endl;
    }
    
}
void AllForceCoplanarity(TPZGeoMesh *gmesh){
    std::setprecision(15);
    int nCoordinates = gmesh->NodeVec().NElements();
    std::vector<int> verify(nCoordinates,-1);
//    TPZGeoNode node = gmesh->NodeVec()[ipoint];
//    TPZVec<REAL> co(3);
//    node.GetCoordinates(co);
    int nmods = 0;
    int nels =gmesh->NElements();
    std::set<std::pair<TPZGeoElSide, std::set<int>>> infNonCoplanar;
    int coplanaarfaces =0;
    int Noncoplanaarfaces =0;
    for(int iel = 0; iel< nels; iel++){
        TPZGeoEl *gel = gmesh->Element(iel);
        int nnodesEl = gel->NNodes();
        if(nnodesEl!=8){
            continue;
        }
        
        TPZFMatrix<REAL> cooridnates;
        int init=gel->FirstSide(2)-1;
        int end = gel->FirstSide(3);
        int nside = gel->NSides();
        if(gel->Index()==4676){
            int ok=0;
        }
        for(int iside= init; iside < end; iside++){
            int nNodesSide = gel->NSideNodes(iside);
            if(nNodesSide<3){
                continue;
            }
            TPZGeoNode *node1 =  gel->SideNodePtr(iside, 0);
            TPZGeoNode *node2 =  gel->SideNodePtr(iside, 1);
            TPZGeoNode *node3 =  gel->SideNodePtr(iside, 2);
            TPZGeoNode *node4 =  gel->SideNodePtr(iside, 3);

            TPZVec<REAL> co1(3),co2(3),co3(3),co4(3);
            node1->GetCoordinates(co1);
            node2->GetCoordinates(co2);
            node3->GetCoordinates(co3);
            node4->GetCoordinates(co4);

            int id1 = node1->Id();
            int id2 = node2->Id();
            int id3 = node3->Id();
            int id4 = node4->Id();

            bool is_coplanar = IsCoplanar(co1, co2, co3, co4);
            if(!is_coplanar){
                std::set<int> indexes;
                indexes.insert(id1);
                indexes.insert(id2);
                indexes.insert(id3);
                indexes.insert(id4);
                TPZGeoElSide gelside(gel,iside);
                std::pair<TPZGeoElSide, std::set<int>> pair = std::make_pair(gelside, indexes);
                infNonCoplanar.insert(pair);
                Noncoplanaarfaces++;
            }
            coplanaarfaces++;

        }
    }
    int ok=0;
    std::map<int, std::set<TPZGeoElSide>> idCommonGeoEls;
//    for(auto it = infNonCoplanar.begin(); it!=infNonCoplanar.end(); it++){
//        std::pair<TPZGeoElSide, std::set<int>> puntos_partida = *it;
//        for(auto it2 = it; it2!=infNonCoplanar.end(); it2++){
//            std::pair<TPZGeoElSide, std::set<int>> puntos_llegada = *it2;
//            std::set<int> intersection;
//            std::set_intersection(puntos_partida.second.begin(), puntos_partida.second.end(),
//                                  puntos_llegada.second.begin(), puntos_llegada.second.end(),
//                                      std::inserter(intersection, intersection.begin()));
//            int ncommon =intersection.size();
//            if(ncommon==4 || ncommon ==0){
//                continue;
//            }
//            if(ncommon>2){
//                DebugStop();
//            }
//            std::cout<<"ncommon: "<<ncommon<< std::endl;
//            for(auto itFin = intersection.begin(); itFin!=intersection.end(); itFin++){
//                idCommonGeoEls[*itFin].insert(puntos_partida.first);
//                idCommonGeoEls[*itFin].insert(puntos_llegada.first);
//            }
//
//        }
    
    //TestRapido
    
    for(auto it = infNonCoplanar.begin(); it!=infNonCoplanar.end(); it++){
        std::pair<TPZGeoElSide, std::set<int>> puntos_partida = *it;
        for(auto it2 = it; it2!=infNonCoplanar.end(); it2++){
            std::pair<TPZGeoElSide, std::set<int>> puntos_llegada = *it2;
            std::set<int> intersection;
            std::set_intersection(puntos_partida.second.begin(), puntos_partida.second.end(),
                                  puntos_llegada.second.begin(), puntos_llegada.second.end(),
                                      std::inserter(intersection, intersection.begin()));
            int ncommon =intersection.size();
            if(ncommon==4 || ncommon ==0){
                continue;
            }
            if(ncommon>2){
                DebugStop();
            }
            std::cout<<"ncommon: "<<ncommon<< std::endl;
            for(auto itFin = intersection.begin(); itFin!=intersection.end(); itFin++){
                idCommonGeoEls[*itFin].insert(puntos_partida.first);
                idCommonGeoEls[*itFin].insert(puntos_llegada.first);
//                verify[*itFin]=1;
            }
            
        }
    }
        //testRaapido
    int contador =0;
    int nose =0;
    int sise =0;
        for(auto it = infNonCoplanar.begin(); it!=infNonCoplanar.end(); it++){
            std::pair<TPZGeoElSide, std::set<int>> puntos_partida = *it;
            TPZGeoElSide gelsideTest = puntos_partida.first;
            TPZGeoEl *gel = gelsideTest.Element();
            int iside = gelsideTest.Side();
            TPZGeoNode *node1 =  gel->SideNodePtr(iside, 0);
            TPZGeoNode *node2 =  gel->SideNodePtr(iside, 1);
            TPZGeoNode *node3 =  gel->SideNodePtr(iside, 2);
            TPZGeoNode *node4 =  gel->SideNodePtr(iside, 3);
            
            TPZVec<REAL> co1(3),co2(3),co3(3),co4(3);
            node1->GetCoordinates(co1);
            node2->GetCoordinates(co2);
            node3->GetCoordinates(co3);
            node4->GetCoordinates(co4);
            
            int id1 = node1->Id();
            int id2 = node2->Id();
            int id3 = node3->Id();
            int id4 = node4->Id();
            
            if(verify[id1]==-1){
                verify[id1]=1;
                verify[id2]=1;
                verify[id3]=1;
                verify[id4]=1;
                ForceCoplanarity(co2,co3,co4,co1);
                node1->SetCoord(co1);
                sise++;
                continue;
            }
            if(verify[id2]==-1){
                verify[id1]=1;
                verify[id2]=1;
                verify[id3]=1;
                verify[id4]=1;
                ForceCoplanarity(co1,co3,co4,co2);
                node2->SetCoord(co2);
                sise++;
                continue;
            }
            if(verify[id3]==-1){
                verify[id1]=1;
                verify[id2]=1;
                verify[id3]=1;
                verify[id4]=1;
                ForceCoplanarity(co1,co2,co4,co3);
                node3->SetCoord(co3);
                sise++;
                continue;
            }
            if(verify[id4]==-1){
                verify[id1]=1;
                verify[id2]=1;
                verify[id3]=1;
                verify[id4]=1;
                ForceCoplanarity(co1,co2,co3,co4);
                node4->SetCoord(co4);
                sise++;
                continue;
            }
            std::cout<<"no se que hacer =("<<std::endl;
            nose++;
          
      }
    ok=0;
    
}
bool ForceCoplanarity(TPZVec<REAL> &c1, TPZVec<REAL> &c2, TPZVec<REAL> &c3, TPZVec<REAL> &c4){
    TPZVec<TPZVec<REAL>> v = {{c2[0] - c1[0], c3[0] - c1[0], c4[0] - c1[0]},
                                              {c2[1] - c1[1], c3[1] - c1[1], c4[1] - c1[1]},
                                              {c2[2] - c1[2], c3[2] - c1[2], c4[2] - c1[2]}};

    REAL det = determinante(v);
    REAL cordant = c1[2];
        // Si los puntos ya son coplanares
        if (std::abs(det) < 1e-12) {
            std::cout << "Los puntos ya son coplanares" << std::endl;
            
            return false;
        }

        // Si no son coplanares, calculamos la nueva componente z de c1
    REAL z = c1[2] - det / ((c2[0] - c1[0]) * (c3[1] - c1[1]) - (c2[1] - c1[1]) * (c3[0] - c1[0]));

        c1[2] = z;
    
    std::cout<<" cord ant: "<<cordant <<" cord new: "<< z <<std::endl;
    
    return true;
   
}
// Cálculo del determinante de una matriz 3x3
double determinante(TPZVec<TPZVec<REAL>> matriz) {
    double sum = matriz[0][0] * ((matriz[1][1] * matriz[2][2]) - (matriz[2][1] * matriz[1][2])) -
                 matriz[0][1] * ((matriz[1][0] * matriz[2][2]) - (matriz[1][2] * matriz[2][0])) +
                 matriz[0][2] * ((matriz[1][0] * matriz[2][1]) - (matriz[1][1] * matriz[2][0]));

    return sum;
}
bool IsCoplanar(TPZVec<REAL> &c1, TPZVec<REAL> &c2, TPZVec<REAL> &c3, TPZVec<REAL> &c4){
    TPZVec<TPZVec<REAL>> v = {{c2[0] - c1[0], c3[0] - c1[0], c4[0] - c1[0]},
                                              {c2[1] - c1[1], c3[1] - c1[1], c4[1] - c1[1]},
                                              {c2[2] - c1[2], c3[2] - c1[2], c4[2] - c1[2]}};

    REAL det = determinante(v);
    REAL cordant = c1[2];
        // Si los puntos ya son coplanares
        if (std::abs(det) < 1e-12) {
            std::cout << "Los puntos ya son coplanares" << std::endl;
            
            return true;
        }
    
    return false;
    
}
