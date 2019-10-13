
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgmesh.h"
#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZVTKGeoMesh.h"

// #include "pzlog.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>

#include <set>
#include <map>
#include <vector>

#include "DFNFractureMesh.h"
#include "DFNRibs.h"
#include "DFNFace.h"



#include "TPZRefPatternDataBase.h"

#include <gmsh.h>

//MATERIAL ID MAP
// 1 gmesh (default)
// 4 skeleton
// 12 ribs that will be divided
// 19 children ribs
// 20 mid-fracture cut faces
// 35 end-fracture cut faces
// 40 Fracture plane
// 45 Intersection points in end-faces

using namespace std;

int main()
{

	// gmsh::initialize();
    // gmsh::option::setNumber("Mesh.Algorithm", 1); // (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
    
    //     // points
	// 	gmsh::model::geo::addPoint(-1.,-1.,0.,0.,1);
	// 	gmsh::model::geo::addPoint(1.,-1.,0.,0.,2);
	// 	gmsh::model::geo::addPoint(1.,1.,0.,0.,3);
	// 	gmsh::model::geo::addPoint(-1.,1.,0.,0.,4);
	// 	gmsh::model::geo::addPoint(-0.5,-1.,0.,0.,5);
	// 	gmsh::model::geo::addPoint(0.5,1.,0.,0.,6);
    //     // lines
	// 	gmsh::model::geo::addLine(1, 5, 1);
  	// 	gmsh::model::geo::addLine(5, 2, 2);
  	// 	gmsh::model::geo::addLine(2, 3, 3);
  	// 	gmsh::model::geo::addLine(3, 6, 4);
  	// 	gmsh::model::geo::addLine(6, 4, 5);
  	// 	gmsh::model::geo::addLine(4, 1, 6);
  	// 	gmsh::model::geo::addLine(5, 6, 7);
	// 	// surfaces
	// 	gmsh::model::geo::addCurveLoop({1, 2, 3, 4, 5, 6}, 1);
  	// 	gmsh::model::geo::addPlaneSurface({1}, 1);
		    
    // 	// line in surface
    //     gmsh::model::geo::synchronize();
    //     gmsh::model::mesh::embed(1,{7},2,1);
    // // // PHYSICAL GROUPS
    // //     // physical curve
    // //     std::vector<int> * allLines;
    // //     std::vector<int> * auxVector;
    // //     if(curvesInSurface.size() > edgeloopvector.size()){
    // //         auxVector = &edgeloopvector;
    // //     }else{
    // //         auxVector = &curvesInSurface;
    // //     }
    // //     allLines->insert(allLines->end(), auxVector->begin(), auxVector->end() );
    // //     gmsh::model::addPhysicalGroup(1,*allLines,fSurfaceMaterial);
    // //     // physical surface
    // //     gmsh::model::addPhysicalGroup(2,{surfaceIndex},fSurfaceMaterial);

    // // // synchronize
    // //     gmsh::model::geo::synchronize();
    // // mesh
    //     gmsh::model::mesh::generate(2);
    // // write (for testing)
    //     gmsh::write("testAPI.msh");
	// gmsh::finalize();

	// Creating the Geo mesh

	int dimel = 2;
	TPZManVector<REAL, 3> x0(3, 0.), x1(3, 2);
	x1[2] = 0.;
	TPZManVector<int, 2> nelx(2, dimel);
	TPZGenGrid gengrid(nelx, x0, x1);
	gengrid.SetElementType(EQuadrilateral);
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	gmesh->SetDimension(2);
	gengrid.Read(gmesh);

	// Mesh 3D

	TPZExtendGridDimension extend(gmesh,1.);
	extend.SetElType(1);
	TPZGeoMesh *gmesh3d = extend.ExtendedMesh(dimel);
	gmesh = gmesh3d;
	//std::ofstream out3("3DMESH.vtk");
	//TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out3, true);

	//    Reading coordinates of a plane from txt file
	string value;
	int npoints = 0;
	string line;
	// count number of corners
	while (npoints == 0){
		ifstream plane_file("fracture.txt");
		if (!plane_file){
			std::cout << "Error reading file" << std::endl;
			DebugStop();
		}
		getline(plane_file, line);
		std::stringstream ss(line);
		while (getline(ss, value, ' ')){
			while (value.length() == 0){
				getline(ss, value, ' ');
			}
			npoints++;
		}
	}
	// then read points
	Matrix plane(3, npoints);
{ //just a scope
	int i = 0;
	int j = 0;
	std::cout << "Fracture plane defined as: \n";
	ifstream plane_file("fracture.txt");
	while (getline(plane_file, line)){
		std::stringstream ss(line);
		while (getline(ss, value, ' ')){
			while (value.length() == 0){
				getline(ss, value, ' ');
			}
			plane(i, j) = std::stod(value);
			std::cout << plane(i, j) << (j<npoints-1?" , ":"\n");
			j++;
		}
		j = 0;
		i++;
	}
	std::cout << std::endl;
}







	// example fractures
	// 2.9084405 1.6236484 0.091559516 1.3763516
	// 2.5516489 2.2694336 2.4483511 2.7305664
	// 1.3832619 2.8898022 1.6167381 0.11019779

	// 7x7x7 (0.5)
	// 2.8 1.2 1.2 2.8
	// 2.8 2.8 1.2 1.2
	// 1.3 1.3 1.3 1.3





	//  Construction of fracplane and FractureMesh
	DFNFracPlane fracplane(plane);
	DFNFractureMesh fracmesh(fracplane, gmesh, 40);

	// Find and split intersected ribs
	gRefDBase.InitializeUniformRefPattern(EOned);
	fracmesh.SplitRibs(19);

	// Find and split intersected faces
	fracmesh.SplitFaces(18);
	// Split edge of fracture
	// fracmesh.SplitFractureEdge();

	// triangulation of fracture plane
	fracmesh.SplitFracturePlane();
	gmesh->BuildConnectivity();
	fracmesh.CreateSkeletonElements(1,19);

	// // Mesh transition volumes
	// fracmesh.CreateVolumes();


	//Print result
	std::ofstream meshprint1("meshprint1.txt");
	std::ofstream out1("./TestSurfaces.vtk");
	gmesh->Print(meshprint1);
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out1, true);























/*
// SECOND PLANE


	//    Reading coordinates of a plane from txt file
	// count number of points
	npoints = 0;
	// count number of corners
	while (npoints == 0){
		ifstream plane_file("fracture2.txt");
		if (!plane_file){
			std::cout << "Error reading file" << std::endl;
			DebugStop();
		}
		getline(plane_file, line);
		std::stringstream ss(line);
		while (getline(ss, value, ' ')){
			while (value.length() == 0){
				getline(ss, value, ' ');
			}
			npoints++;
		}
	}
	// then read points
	plane.Resize(3, npoints);
{ //just a scope
	int i = 0;
	int j = 0;
	std::cout << "\n Second Fracture plane defined as: \n";
	ifstream plane_file("fracture2.txt");
	while (getline(plane_file, line)){
		std::stringstream ss(line);
		while (getline(ss, value, ' ')){
			while (value.length() == 0){
				getline(ss, value, ' ');
			}
			plane(i, j) = std::stod(value);
			std::cout << plane(i, j) << (j<npoints-1 ?" , ":"\n");
			j++;
		}
		j = 0;
		i++;
	}
	std::cout << std::endl;
}




	//  Construction of fracplane and FractureMesh
	DFNFracPlane fracplane2(plane);
	DFNFractureMesh fracmesh2(fracplane2, gmesh, 40);

	// Find and split intersected ribs
	fracmesh2.SplitRibs(19);
	gmesh->BuildConnectivity();
	// Find and split intersected faces
	fracmesh2.SplitFaces(18);
	// Split edge of fracture
	fracmesh2.SplitFractureEdge();

	// triangulation of fracture plane
	// fracmesh2.SplitFracturePlane();

	// Mesh transition volumes
	// fracmesh2.CreateVolumes();
	gmesh->BuildConnectivity();
	fracmesh2.CreateSkeletonElements(1,19);

	//Print result
	std::ofstream meshprint2("meshprint2.txt");
	std::ofstream out2("./TestSurfaces2.vtk");
	gmesh->Print(meshprint2);
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2, true);
*/
	return 0;
}

