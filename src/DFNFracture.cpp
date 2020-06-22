/*! 
 *  @brief     Contains implementation of class DFNFracture
 *  @authors   Pedro Lima
 *  @date      2019
 */

#include "DFNFracture.h"
#include "DFNMesh.h"
#include <math.h>
#include <cstdio>
// #include <unordered_set>
#include <algorithm>
#include "TPZRefPatternDataBase.h"
#include "TPZGeoMeshBuilder.h"

// Empty Constructor
DFNFracture::DFNFracture(){
}

// Constructor with corner points, a geomesh and material ID
DFNFracture::DFNFracture(DFNFracPlane &FracPlane, DFNMesh *dfnMesh){
    fFracplane = FracPlane;
    fdfnMesh = dfnMesh;

    // Set corner nodes of fracture into mesh
    if(fdfnMesh->Dimension() == 3){
        TPZManVector<int64_t,4> nodeindices = fFracplane.SetPointsInGeomesh(fdfnMesh->Mesh());
    }
}

// Copy constructor
DFNFracture::DFNFracture(const DFNFracture &copy){
    this->operator=(copy);
}

// Assignment operator
DFNFracture &DFNFracture::operator=(const DFNFracture &copy){
    fdfnMesh = copy.fdfnMesh;
    fTolerance = copy.GetTolerance();
    fRibs = copy.fRibs;
	fFaces = copy.fFaces;
    fFracplane = copy.fFracplane;
    return *this;
}

/**
 * @brief Get the fracture plane
 * @return The plane corner coordinates
 */

DFNFracPlane &DFNFracture::FracPlane() {
    return fFracplane;
}

/**
 * @brief Set the tolerance for the distance between a point-plane
 * @param Tolerance
 */

void DFNFracture::SetTolerance(REAL tolerance){
    fTolerance = tolerance;
}

/**
 * @brief Get the tolerance
 * @return The tolerance
 */

REAL DFNFracture::GetTolerance() const{
    return fTolerance;
}









void DFNFracture::AddFace(DFNFace &face){
    int index= face.Index();
    fFaces.emplace(index,face);
}
void DFNFracture::AddRib(DFNRib &rib){
    int index= rib.Index();
    fRibs.emplace(index,rib);
}

// /**
//  * @brief Add cut faces using indexes
//  * @param Face to be set
//  */

// bool DFNFracture::AddMidFace(DFNFace &face){
//     // iterate over ribs to find intersected ones
//     int nribscut = 0;
//     TPZManVector<int64_t,2> CutRibsIndex(2);
//     TPZVec<int64_t> rib_index = face.GetRibs();
//     int nribs = rib_index.size();
//     for (int irib = 0; irib < nribs; irib++)
//     {
//         DFNRib *ribtest = Rib(rib_index[irib]);
//         if(ribtest->IsIntersected() == true){
//             CutRibsIndex[nribscut]=rib_index[irib];
//             nribscut++;
//         }
//     }
//     // std::cout<<"first rib: "<<CutRibsIndex[0]<<std::endl;
//     // std::cout<<"second rib: "<<CutRibsIndex[1]<<std::endl;
//     TPZGeoMesh *gmesh = fdfnMesh->Mesh();
//     // Connect intersection points
//     TPZVec<int64_t> ipoints(2);
//     ipoints[0] = Rib(CutRibsIndex[0])->IntersectionIndex();
//     ipoints[1] = Rib(CutRibsIndex[1])->IntersectionIndex();
//     int64_t nels = gmesh->NElements();
//     int matid = DFNMaterial::Efracture;
//     this->fSurface[nels] = gmesh->CreateGeoElement(EOned, ipoints, matid, nels);
    
//     int index= face.Index();
//     fFaces[index]=face;
//     return true;
// }

// /**
//  * @brief Add faces that are cut at the edges of fracture (using indexes)
//  * @param Face to be set
//  */

// bool DFNFracture::AddEndFace(DFNFace &face){
//     // Create geometric element for intersection node beetween EndFace and fracture edge
//         TPZManVector<REAL,3> coords(3);
//         bool ipoint_exists = FindEndFracturePoint(face,coords);
//         if(ipoint_exists == false) return false;
//         TPZManVector<int64_t> ipoints(2);
//         ipoints[0] = fdfnMesh->Mesh()->NodeVec().AllocateNewElement();
//         fdfnMesh->Mesh()->NodeVec()[ipoints[0]].Initialize(coords, *fdfnMesh->Mesh());
//         face.SetIntersectionIndex(ipoints[0]);
//     //iterate over ribs to connect intersection points
//     TPZVec<int64_t> rib_index = face.GetRibs();
//     int nribs = rib_index.size();
//     for (int irib = 0; irib < nribs; irib++)
//     {
//         DFNRib *ribtest = Rib(rib_index[irib]);
//         if(ribtest->IsIntersected() == true){
//             // std::cout<<"single rib cut: "<<rib_index[irib]<<std::endl;
//             // Connect intersection points
//                 int64_t nels = fdfnMesh->Mesh()->NElements();
//                 ipoints[1] = Rib(rib_index[irib])->IntersectionIndex();
//                 this->fSurface[nels] = fdfnMesh->Mesh()->CreateGeoElement(EOned, ipoints, DFNMaterial::Efracture, nels);
//                 break;
//         }
//     }
//     int index= face.Index();
//     fFaces[index]=face;
//     return true;
// }









DFNRib * DFNFracture::Rib(int64_t index){
    auto candidate = fRibs.find(index);
    if(candidate != fRibs.end()){
        return &candidate->second;
    }
    return nullptr;
}
DFNFace * DFNFracture::Face(int64_t index){
    auto candidate = fFaces.find(index);
    if(candidate != fFaces.end()){
        return &candidate->second;
    }
    return nullptr;
}


void DFNFracture::FindFaces(){
    TPZGeoMesh *gmesh = fdfnMesh->Mesh();
    TPZGeoEl *gel;
 
    // iterate over 2D elements and check their 1D neighbours for intersections
    int64_t nel = gmesh->NElements();
    for(int iel=0; iel<nel; iel++){
        gel = gmesh->Element(iel);
        if(!gel)continue;
        if(gel->Dimension() != 2) continue;
        if(gel->HasSubElement() ) continue;
 
        int nnodes = gel->NCornerNodes();
        int nedges = nnodes;
 
        // gather ribs
        TPZManVector<int64_t,4>  rib_index(nedges,-1);
        for(int iedge = 0; iedge < nedges; iedge++){
            TPZGeoElSide gelside(gel,iedge+nnodes);
            TPZGeoElSide neig = gelside.Neighbour();
            for(/*void*/;neig != gelside; neig = neig.Neighbour()){
                if(neig.Element()->Dimension() != 1) continue;
                rib_index[iedge] = neig.Element()->Index();
            }
        }

        // build face
        DFNFace candidate(gel,this);
        candidate.SetRibs(rib_index);
        candidate.UpdateStatusVec();
        if(!candidate.IsIntersected()) continue;
        if(candidate.IsOnBoundary()){
            candidate.FindInPlanePoint();
        }
        gel->SetMaterialId(DFNMaterial::Erefined);
        candidate.UpdateRefMesh();
        AddFace(candidate);
    }
}

// /**
//  * @brief Create cut surfaces
//  * @param Material id
//  * @note This method was getting too long. Moved some of it into AddEndFace and AddMidFace.
//  */

// void FindFacesOld(){
//     TPZGeoMesh *gmesh = fdfnMesh->Mesh();
//     TPZGeoEl *gel;

//     // iterate over all 2D elements and check their 1D neighbours for intersections
//     int64_t nel = gmesh->NElements();
//     for(int iel = 0; iel<nel; iel++){
//         gel = gmesh->Element(iel);
//         if(!gel) continue;
//         if(gel->Dimension() != 2) continue;
//         if(gel->HasSubElement()) continue;
        
//         int nribscut = 0;
//         int nedges = gel->NCornerNodes();

//         // vector with status for each node and rib of face
//         TPZManVector<int> sidestatus(nedges*2+1,0);
//         // TPZManVector<int64_t,2> CutRibsIndex(2);
//         // vector with indices of father ribs that outline the face
//         TPZManVector<int64_t,4> rib_index(nedges,-1);

//         // iterate over ribs to check for intersection
//         for(int iside = 0; iside < nedges; iside++){
//             TPZGeoElSide gelside(gel,iside+nedges);
//             TPZGeoElSide neig;
//             // TPZGeoElSide neig = gelside.Neighbour();
//             // while(neig.Element()->Dimension()!=1 || !neig.Element()->HasSubElement()){ 
//             //     neig=neig.Neighbour();
//             // }
//             // rib_index[iside] = neig.Element()->Index();
//             for(neig = gelside.Neighbour(); neig!= gelside; neig = neig.Neighbour()){
//                 if(neig.Element()->Dimension() == 1 && neig.Element()->HasSubElement()) {
//                     rib_index[iside] = neig.Element()->Index();
//                     break;
//                 }
//                 if(neig.Element()->Dimension() == 1) rib_index[iside] = neig.Element()->Index();
//             }
//             DFNRib *ribtest = Rib(rib_index[iside]);
//             if(ribtest->IsIntersected()==true){
//             //check if ribtest was divided into two ribs, or a rib and a point
//                 // get node where ribtest is cut
//                 int64_t cutnode = ribtest->IntersectionIndex();
//                 TPZVec<int64_t> ribtestNodes(2);
//                 // gmesh->Element(ribtest->Index())->GetNodeIndices(ribtestNodes);
//                 ribtestNodes[0] = gel->SideNodeIndex(iside+nedges,0);
//                 ribtestNodes[1] = gel->SideNodeIndex(iside+nedges,1);
//                 // check if intersection is a corner node and assign status accordingly
//                 if(cutnode == ribtestNodes[0]){
//                     sidestatus[iside] = true;
//                 }
//                 else if(cutnode == ribtestNodes[1]){
//                     sidestatus[(iside+1)%nedges] = true;
//                 }
//                 else{
//                     sidestatus[iside+nedges] = true;
//                 }
                
//                 nribscut++;
//             }
//         }

//         // if there are ribs cut, create a face object
//         if(nribscut == 0){continue;}

//         int edgescut = 0, cornerscut = 0;
//         for(int i = 0; i<nedges; i++){
//             if(sidestatus[i]) cornerscut++;
//             if(sidestatus[i+nedges]) edgescut++;
//         }
//         if(edgescut == 0 && cornerscut == 1) nribscut = 1;
//         // Create DFNFace
//         DFNFace face(gel, this);
//         face.SetRibs(rib_index); 
//         // During development, elements at fracture surface have material id over DFNMaterial::Efracture
//         // if(gel->MaterialId() != DFNMaterial::Efracture) {gel->SetMaterialId(matID);}
//         face.StatusVec() = sidestatus;
//         face.SetFracture(this);
//         // if(nribscut == 1) {gel->SetMaterialId(matID+17);} // this is here for graphical debugging only... comment it on release

//         // Add face to map
//         bool face_ok = false;
//         switch (nribscut){
//             case  2: face_ok = AddMidFace(face);break;
//             case  1: face_ok = AddEndFace(face);break;
//             default: std::cout<<"\n\n"<<__PRETTY_FUNCTION__<<"\n Face # "<<iel<<"\nNo more than 2 ribs should've been cut\n\n";DebugStop();
//         }
//         if(face_ok){ 
//             if(gel->MaterialId() != DFNMaterial::Efracture) {gel->SetMaterialId(DFNMaterial::Erefined);}
//             Face(iel)->Refine(gel->MaterialId());
//         }

//         // //Print result
// 		// 	std::ofstream out1("./TestSurfaces.vtk");
// 		// 	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out1, true);
//     }
//     gmesh->BuildConnectivity();
//     fdfnMesh->CreateSkeletonElements(1, DFNMaterial::Erefined);
// }





















void DFNFracture::FindRibs(){
    //search gmesh for intersected ribs
    int64_t Nels = fdfnMesh->Mesh()->NElements();
    TPZManVector<int64_t, 2> inode(2,0);
    for (int iel = 0; iel < Nels; iel++){
        TPZGeoEl *gel = fdfnMesh->Mesh()->Element(iel);
        //skip all elements that aren't ribs
        if (gel->Dimension() != 1){continue;}
        // skip all elements that have been cut by a previous fracture
        if(gel->HasSubElement()){continue;}

        // Check rib
        TPZManVector<REAL,3> intpoint(3,0);
        bool resul = fFracplane.Check_rib(gel, &intpoint);

        // Add rib
        if (resul == true){
            DFNRib rib(gel, this);
            rib.StatusVec()[2] = 1; //StatusVec = {0,0,1}
            rib.SetIntersectionCoord(intpoint);
            if(gel->MaterialId() != DFNMaterial::Efracture) {gel->SetMaterialId(DFNMaterial::Erefined);}
            AddRib(rib);
        }
    }
}


void DFNFracture::RefineRibs(){
    for(auto itr = fRibs.begin(); itr!=fRibs.end(); itr++){
        DFNRib* rib = &itr->second;
        rib->Refine();
    }
}
void DFNFracture::RefineFaces(){
    for(auto itr = fFaces.begin(); itr!=fFaces.end(); itr++){
        DFNFace *face = &itr->second;
        face->Refine();
    }
}





void DFNFracture::OptimizeRibs(REAL tolDist){
    for(auto itr = fRibs.begin(); itr!=fRibs.end(); itr++){
        DFNRib* rib = &itr->second;
        rib->Optimize(tolDist);
    }
}
void DFNFracture::OptimizeFaces(REAL tolDist, REAL tolAngle){
    for(auto itr = fFaces.begin(); itr!=fFaces.end(); itr++){
        DFNFace* face = &itr->second;
        face->Optimize(tolDist, tolAngle);
    }
}










void DFNFracture::SplitFractureEdge(std::list<int> &fracEdgeLoop){ 
    
    // Compute fracture edges' lenghts
    int nedges = fFracplane.GetCornersX().Cols();
    TPZVec<REAL> edgelength(nedges,0);
    Matrix fraccorners(3,nedges);
    fraccorners = fFracplane.GetCornersX();
    for(int i = 0; i < nedges; i++){
        edgelength[i] = sqrtl(pow(fraccorners(0,i)-fraccorners(0,(i+1)%nedges),2)
                              +pow(fraccorners(1,i)-fraccorners(1,(i+1)%nedges),2)
                              +pow(fraccorners(2,i)-fraccorners(2,(i+1)%nedges),2));
    }
    
	// vector of pointers to maps
    // @todo refactor this to tpzautopointer<map> to prevent memory leak
	TPZManVector<std::map<REAL, int64_t>* > edgemap(nedges);
    for(int i = 0; i < nedges; i++){
        edgemap[i] = new std::map<REAL, int64_t>;
    }
    
    // iterate over all endfaces and map it to the fracture-edge that cuts it
    for (auto it = fFaces.begin(); it != fFaces.end(); it++){
        // get intersection node index and coordinates
        DFNFace *iface = &it->second;
        int64_t ipointindex = iface->IntersectionIndex();
        TPZVec<REAL> ipointcoord(3);
        fdfnMesh->Mesh()->NodeVec()[ipointindex].GetCoordinates(ipointcoord);

        // iterate over edges to check if ipoint belongs to it
        for(int iedge = 0; iedge < nedges; iedge++){
			//vectors from ipoint to iedge's nodes
			TPZManVector<REAL, 3> v1(3);
				v1[0] = fraccorners(0,iedge) - ipointcoord[0];
				v1[1] = fraccorners(1,iedge) - ipointcoord[1];
				v1[2] = fraccorners(2,iedge) - ipointcoord[2];
			TPZManVector<REAL, 3> v2(3);
				v2[0] = fraccorners(0,(iedge+1)%nedges) - ipointcoord[0];
				v2[1] = fraccorners(1,(iedge+1)%nedges) - ipointcoord[1];
				v2[2] = fraccorners(2,(iedge+1)%nedges) - ipointcoord[2];
			// square of cross product
			REAL temp = pow(v1[1]*v2[2] - v1[2]*v2[1],2);
				temp += pow(v1[2]*v2[0] - v1[0]*v2[2],2);
				temp += pow(v1[0]*v2[1] - v1[1]*v2[0],2);
				temp = sqrtl(temp);
            // check if point is in edge by calculating if it's normal distance to the edge is zero
            REAL dist = temp/edgelength[iedge];
			if(dist<fTolerance){
				// compute local 1D coordinate (alpha)
				REAL norm = sqrtl(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
				REAL alpha = norm/edgelength[iedge];
				if(alpha > 1+fTolerance){std::cout<<"\nProblem with alpha\n";DebugStop();}
				// map intersection indexes from smallest alpha to biggest
                // map <alpha, index>
				edgemap[iedge]->insert({alpha, ipointindex});
			}
        }
    }
    
    // GMsh does not accept zero index entities
    const int shift = 1;

    bool warning_message_Q = true;
    //Once intersections on fracture-edges have been properly ordered and mapped by edge
	//iterate over edges to split them
	for (int iedge = 0; iedge < nedges; iedge++)
	{
		int64_t nels = fdfnMesh->Mesh()->NElements();
		TPZManVector<int64_t,2> inodes(2);     //index of nodes to be connected
        int icorner = iedge; //for readability
        
        if(edgemap[iedge]->size() == 0){
            #ifdef PZDEBUG
                if(warning_message_Q){
                    std::cout<<"\n Warning: Is there an edge of a fracture that doesn't cut any element? \n\n";
                    warning_message_Q = false;
                }
            #endif //PZDEBUG
            // DebugStop();
            inodes[0] = fFracplane.CornerIndex(icorner);
            inodes[1] = fFracplane.CornerIndex((icorner+1)%nedges);
            this->fSurface[nels] = fdfnMesh->Mesh()->CreateGeoElement(EOned, inodes, DFNMaterial::Efracture+6, nels);
            fracEdgeLoop.push_back((int) nels + shift);
            continue;
        }


		// connect first end-face intersection to iedge's first node
		auto it = edgemap[iedge]->begin();
        inodes[0] = fFracplane.CornerIndex(icorner);
		inodes[1] = it->second;
		this->fSurface[nels] = fdfnMesh->Mesh()->CreateGeoElement(EOned, inodes, DFNMaterial::Efracture+6, nels);
        fracEdgeLoop.push_back((int) nels + shift);
        
        // iterate over iedge's map
        while(++it != edgemap[iedge]->end()){
            nels++;
            inodes[0] = inodes[1];
            inodes[1] = it->second;
            this->fSurface[nels] = fdfnMesh->Mesh()->CreateGeoElement(EOned, inodes, DFNMaterial::Efracture+6, nels);
            fracEdgeLoop.push_back((int) nels + shift);
        }

		// connect last end-intersection to edge last node
		nels++;
		inodes[0] = inodes[1];
        inodes[1] = fFracplane.CornerIndex((icorner+1)%nedges);
		this->fSurface[nels] = fdfnMesh->Mesh()->CreateGeoElement(EOned, inodes, DFNMaterial::Efracture+6, nels);
        fracEdgeLoop.push_back((int) nels + shift);
    }
	
fdfnMesh->Mesh()->BuildConnectivity();  
}












void DFNFracture::MeshFractureSurface(){
    // GMsh does not accept zero index entities
    const int shift = 1;
    // First construct the edges of the fracture surface
    std::list<int> fracEdgeLoop;
    SplitFractureEdge(fracEdgeLoop);

    // initialize GMsh
    gmsh::initialize();
    gmsh::option::setNumber("Mesh.Algorithm", 5); // (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
    // INSERT POINTS
        // iterate over ribs and get their ipoints
        std::set<int64_t> pointset;
        
        for(auto itr = fRibs.begin(); itr != fRibs.end(); itr++){
            DFNRib *irib = &itr->second;
            if (irib->IsIntersected() == false){continue;}
            int64_t inode = irib->IntersectionIndex();
            pointset.insert(inode);
        }
        for(auto point : pointset){
            TPZManVector<REAL, 3> coord(3);
            fdfnMesh->Mesh()->NodeVec()[point].GetCoordinates(coord);

            gmsh::model::geo::addPoint(coord[0],coord[1],coord[2],0.,point+shift);

        }
        // iterate over endFaces and get their ipoints
        for(auto itr = fFaces.begin(); itr != fFaces.end(); itr++){
            DFNFace *iface = &itr->second;
            if (iface->IsIntersected() == false){continue;}
            int64_t inode = iface->IntersectionIndex();
            TPZManVector<REAL, 3> coord(3);
            fdfnMesh->Mesh()->NodeVec()[inode].GetCoordinates(coord);

            gmsh::model::geo::addPoint(coord[0],coord[1],coord[2],0.,inode+shift);
        }
        // Corners of fracture plane
        {
            int ncorners = fFracplane.GetCornersX().Cols();
            for(int i = 0; i<ncorners; i++){
                int64_t inode = fFracplane.CornerIndex(i);
                TPZManVector<REAL, 3> coord(3);
                fdfnMesh->Mesh()->NodeVec()[inode].GetCoordinates(coord);

                gmsh::model::geo::addPoint(coord[0],coord[1],coord[2],0.,inode+shift);
            }
        }
    // INSERT LINES
        // move loop list into a vector
        std::vector<int> edgeloopvector{std::make_move_iterator(std::begin(fracEdgeLoop)),
                                        std::make_move_iterator(std::end(fracEdgeLoop))};
        std::vector<int> curvesInSurface;
        for(auto iter = fSurface.begin(); iter != fSurface.end(); iter++){
            int64_t iel = iter->first+shift;
            TPZGeoEl *gel = iter->second;
            if(gel->Dimension() != 1) continue;
            int64_t node0 = gel->NodeIndex(0)+shift;
            int64_t node1 = gel->NodeIndex(1)+shift;

            gmsh::model::geo::addLine(node0,node1,iel);
            gmsh::model::geo::mesh::setTransfiniteCurve(iel,2); // to reduce number of nodes created by GMsh
            bool lineIsInEdge = (std::find(edgeloopvector.begin(),edgeloopvector.end(),iel) != edgeloopvector.end());
            if(lineIsInEdge == false){
                curvesInSurface.push_back(iel);
            }
        }
    
    // CURVE LOOP
        int surfaceIndex = 1;
        gmsh::model::geo::addCurveLoop(edgeloopvector,surfaceIndex);

    // SURFACE
        gmsh::model::geo::addPlaneSurface({surfaceIndex},surfaceIndex);
    
    // lines in surface
        // @comment GMsh requires synchronize before embedding geometric entities
        gmsh::model::geo::synchronize();
        gmsh::model::mesh::embed(1,curvesInSurface,2,surfaceIndex);
    // PHYSICAL GROUPS
        // physical curve
        int nlines = curvesInSurface.size() + edgeloopvector.size();
        if(curvesInSurface.size() > edgeloopvector.size()){
            curvesInSurface.reserve(nlines);
            curvesInSurface.insert(curvesInSurface.end(), edgeloopvector.begin(), edgeloopvector.end() );

            gmsh::model::addPhysicalGroup(1,curvesInSurface,DFNMaterial::Efracture+6);
        }else{
            edgeloopvector.reserve(nlines);
            edgeloopvector.insert(edgeloopvector.end(), curvesInSurface.begin(), curvesInSurface.end() );

            gmsh::model::addPhysicalGroup(1,edgeloopvector,DFNMaterial::Efracture+6);
        }
        // physical surface
        gmsh::model::addPhysicalGroup(2,{surfaceIndex},DFNMaterial::Efracture+7);

    // synchronize before meshing
        gmsh::model::geo::synchronize();
    // mesh
        gmsh::model::mesh::generate(2);
        // gmsh::model::mesh::optimize("Netgen");
    // write (for testing)
        // gmsh::write("testAPI.msh");
    // import meshed plane back into PZ geoMesh
        ImportElementsFromGMSH(fdfnMesh->Mesh(),2);
    // close GMsh
    gmsh::finalize();
    
    fdfnMesh->CreateSkeletonElements(1, DFNMaterial::Efracture);
}







































void DFNFracture::ImportElementsFromGMSH(TPZGeoMesh * gmesh, int dimension){
    // GMsh does not accept zero index entities
    const int shift = 1;

    // First, if GMsh has created new nodes, these should be inserted in PZGeoMesh
    // create a map <node,point>
    std::map<int,int> mapGMshToPZ;

    // iterate over ribs and get intersection nodes
    for(auto itr = fRibs.begin(); itr != fRibs.end(); itr++){
        DFNRib *irib = &itr->second;
        if (irib->IsIntersected() == false) continue;
        int pznode = (int) irib->IntersectionIndex() +shift;

        std::vector<size_t> node_identifiers;
        std::vector<double> coord;
        std::vector<double> parametricCoord;
        gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord,0,pznode,true);
        int gmshnode = (int) node_identifiers[0];
        mapGMshToPZ.insert({gmshnode,pznode});
    }
    // iterate over endFaces and get their ipoints
    for(auto itr = fFaces.begin(); itr != fFaces.end(); itr++){
        DFNFace *iface = &itr->second;
        if (iface->IsIntersected() == false) continue;
        int pznode = (int) iface->IntersectionIndex() +shift;
        
        std::vector<size_t> node_identifiers;
        std::vector<double> coord;
        std::vector<double> parametricCoord;
        gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord,0,pznode,true);
        int gmshnode = (int) node_identifiers[0];
        mapGMshToPZ.insert({gmshnode,pznode});
    }
    // Corners of fracture plane
    {
        int ncorners = fFracplane.GetCornersX().Cols();
        for(int i = 0; i<ncorners; i++){
            int pznode = (int) fFracplane.CornerIndex(i) +shift;
            
            std::vector<size_t> node_identifiers;
            std::vector<double> coord;
            std::vector<double> parametricCoord;
            gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord,0,pznode,true);
            int gmshnode = (int) node_identifiers[0];
            mapGMshToPZ.insert({gmshnode,pznode});
        }
    }

    // add new nodes into PZGeoMesh
    {
        // get all nodes from GMsh
            std::vector<size_t> node_identifiers;
            std::vector<double> coord;
            std::vector<double> parametricCoord;
            gmsh::model::mesh::getNodes(node_identifiers, coord, parametricCoord);
        // iterate over node_identifiers
        int nnodes = node_identifiers.size();
        for(int i = 0; i < nnodes; i++){
            int gmshnode = node_identifiers[i];
            // check if it is contained in the map
            if(mapGMshToPZ.find(gmshnode) == mapGMshToPZ.end()){
                // New node -> add to PZGeoMesh
                int pznode = (int) gmesh->NodeVec().AllocateNewElement();
                TPZManVector<REAL,3> newnodeX(3);
                newnodeX[0] = coord[3*i];
                newnodeX[1] = coord[3*i+1];
                newnodeX[2] = coord[3*i+2];
                gmesh->NodeVec()[pznode].Initialize(newnodeX,*gmesh);
                // int pznode = (int) gmesh->NNodes();
                // gmesh->NodeVec().resize(pznode+1);
                // insert it in map
                mapGMshToPZ.insert({gmshnode,pznode+shift});
            }

        }
    }
    
    // remember to use mapGMshToPZ to translate from GMsh node index to PZ nodeindex




    
    int64_t nels = gmesh->NElements();
    std::vector<std::pair<int, int> > dim_to_physical_groups;
    gmsh::model::getPhysicalGroups(dim_to_physical_groups);
   
    std::vector<std::pair<int, int> > entities_0d;
    std::vector<std::pair<int, int> > entities_1d;
    std::vector<std::pair<int, int> > entities_2d;
    std::vector<std::pair<int, int> > entities_3d;
    gmsh::model::getEntities(entities_0d,0);
    gmsh::model::getEntities(entities_1d,1);
    gmsh::model::getEntities(entities_2d,2);
    gmsh::model::getEntities(entities_3d,3);
   
    /// inserting the elements
    for (auto group: dim_to_physical_groups) {
       
        int dim = group.first;
        // only want elements of a given dimension
        if(dim != dimension) continue;
        int physical_identifier = group.second;
       
        std::vector< int > entities;
        gmsh::model::getEntitiesForPhysicalGroup(dim, physical_identifier, entities);
       
        for (auto tag: entities) {
           
            std::vector<int> group_element_types;
            std::vector<std::vector<std::size_t> > group_element_identifiers;
            std::vector<std::vector<std::size_t> > group_node_identifiers;
            gmsh::model::mesh::getElements(group_element_types,group_element_identifiers,group_node_identifiers, dim, tag);
            int n_types = group_element_types.size();
            for (int itype = 0; itype < n_types; itype++){
                int el_type = group_element_types[itype];
                int n_nodes = TPZGeoMeshBuilder::GetNumberofNodes(el_type);
                std::vector<int> node_identifiers(n_nodes);
                int n_elements = group_element_identifiers[itype].size();
                for (int iel = 0; iel < n_elements; iel++) {
                    int el_identifier = group_element_identifiers[itype][iel]-1+nels;
                    for (int inode = 0; inode < n_nodes; inode++) {
                        // node_identifiers[inode] = group_node_identifiers[itype][iel*n_nodes+inode];
                        // Use mapGMshToPZ to translate from GMsh node index to PZ nodeindex
                        node_identifiers[inode] = mapGMshToPZ[group_node_identifiers[itype][iel*n_nodes+inode]];
                    }
                    TPZGeoMeshBuilder::InsertElement(gmesh, physical_identifier, el_type, el_identifier, node_identifiers);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
}


















void DFNFracture::AssembleOutline(){
    // Build skeleton of 1D elements between new subelements
    fdfnMesh->CreateSkeletonElements(1);
    DFNFace *face = nullptr;
    TPZGeoEl *face_gel = nullptr;

    for(auto itr = fFaces.begin(); itr!=fFaces.end(); itr++){
        face = &itr->second;
        face_gel = face->GeoEl();
        int n_intersection_points = 0;
        for(int istate : face->StatusVec()){
            n_intersection_points += istate;
        }
        if(n_intersection_points < 2) continue;

        int nsides = face_gel->NSides();
        int nnodes = face_gel->NCornerNodes();
        TPZStack<int64_t> framenodes;
        for(int i=0; i<nsides; i++){
            if(face->StatusVec()[i]){
                if(i<nnodes){
                    framenodes.push_back(face_gel->NodeIndex(i));
                }else if(i<nsides-1){
                    framenodes.push_back(face->Rib(i-nnodes)->IntersectionIndex());
                }else if(i == nsides-1){
                    framenodes.push_back(face->IntersectionIndex());
                }
            }
        }

        // search for side of child that contains foutline element (frame_el)
        TPZGeoEl *frame_el = nullptr;
        int64_t frame_el_index = -1;
        int nchildren = face_gel->NSubElements();
        if(nchildren == 0){nchildren++;}
        for(int ichild=0; ichild<nchildren; ichild++){
            TPZGeoEl* child;
            if(face_gel->HasSubElement()){
                child = face_gel->SubElement(ichild);
            }else{
                child = face_gel;
            }
            int child_nnodes = child->NCornerNodes();
            int ichildside = -1;
            for(int inode=0; inode<child_nnodes; inode++){
                if(framenodes[0]==child->NodeIndex(inode)){
                    if(framenodes[1]==child->NodeIndex((inode+1)%child_nnodes)){
                        ichildside = inode+child_nnodes;
                        break;
                    }else if(framenodes[1]==child->NodeIndex((inode+child_nnodes-1)%child_nnodes)){
                        ichildside = (inode+child_nnodes-1)%child_nnodes + child_nnodes;
                        break;
                    }
                }
            }
            if(ichildside < 0){continue;}
            TPZGeoElSide gelside(child,ichildside);
            TPZGeoElSide neig = gelside.Neighbour();
            for(/*void*/; neig!=gelside; neig = neig.Neighbour()){
                if(neig.Element()->Dimension() != 1) continue;
                frame_el = neig.Element();
                frame_el_index = frame_el->Index();
                break;
            }
            if(!frame_el) DebugStop();
            break;
        }
        fOutline.insert({frame_el_index,frame_el});
        frame_el->SetMaterialId(DFNMaterial::Efracture);
    }
}