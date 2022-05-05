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


/**
 * @brief Define which example to run. See example file sintax in function definition
 * @param filename: path to the file that defines the example
 * @param polyg_stack: vector to fill with corners of the fractures
 * @param mshfile: [optional] path to .msh file (if applicable)
 * @returns pointer to geometric mesh created/read
*/
TPZGeoMesh* SetupExampleFromFile(std::string filename, std::map<int, DFNRawData>& dfnrawdata, std::string mshfile, REAL& toldist, REAL& tolangle, TPZManVector<std::map<int,std::string>,4>& dim_physical_tag_and_name);


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
				TPZManVector<REAL,2>		& tol
				);


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


void ReadFileJSON(const std::string			& filename,
                std::map<int, DFNRawData>&  map_dfnrawdata,
				std::string 				& mshfile,
				TPZManVector<REAL,3> 		& x0,
				TPZManVector<REAL,3> 		& xf,
				TPZManVector<int,3> 		& nels,
				MMeshType					& eltype,
				TPZManVector<REAL,2>		& tol,
				int                         & prerefine
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
	
	// An option to apply uniform refinement in the mesh before fractures
	if(input.find("PreRefine") != input.end())
		{prerefine = input["PreRefine"];}
	else if(input.find("prerefine") != input.end())
		{prerefine = input["prerefine"];}
	
	//Init matid for fractures
	const int fracInitMatId = input["FractureInitMatId"];
	int actualMatId = fracInitMatId;
	
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


TPZGeoMesh* SetupExampleFromFile(std::string filename, std::map<int, DFNRawData>& dfnrawdata, std::string mshfile, REAL& toldist, REAL& tolangle, int& prerefine, TPZManVector<std::map<int,std::string>,4>& dim_physical_tag_and_name){

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
	else if(extension == ".json" || extension == ".jsonc")
		{ReadFileJSON(filename,dfnrawdata,mshfile,x0,x1,nels,eltype,tol,prerefine);}
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
		gmesh = reader.GeometricGmshMesh(mshfile, gmesh);
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

