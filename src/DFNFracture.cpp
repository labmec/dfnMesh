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
DFNFracture::DFNFracture(DFNFracPlane *FracPlane, DFNMesh *dfnMesh){
    fFracplane = FracPlane;
    fdfnMesh = dfnMesh;

    // Set corner nodes of fracture into mesh
    if(fdfnMesh->Dimension() == 3){
        TPZManVector<int64_t,4> nodeindices = fFracplane->SetPointsInGeomesh(fdfnMesh->Mesh());
    }
}

// Copy constructor
DFNFracture::DFNFracture(const DFNFracture &copy){
    this->operator=(copy);
}

// Assignment operator
DFNFracture &DFNFracture::operator=(const DFNFracture &copy){
    fdfnMesh = copy.fdfnMesh;
    fRibs = copy.fRibs;
	fFaces = copy.fFaces;
    fFracplane = copy.fFracplane;
    return *this;
}

/**
 * @brief Get the fracture plane
 * @return The plane corner coordinates
 */

DFNFracPlane *DFNFracture::FracPlane() {
    return fFracplane;
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
        // gel->SetMaterialId(DFNMaterial::Erefined);
        candidate.UpdateRefMesh();
        AddFace(candidate);
    }
}























void DFNFracture::FindRibs(){
    //search gmesh for intersected ribs
    int64_t Nels = fdfnMesh->Mesh()->NElements();
    TPZManVector<int64_t, 2> inode(2,0);
    for (int iel = 0; iel < Nels; iel++){
        TPZGeoEl *gel = fdfnMesh->Mesh()->Element(iel);
        if(!gel) continue;
        //skip all elements that aren't ribs
        if (gel->Dimension() != 1){continue;}
        // skip all elements that have been cut by a previous fracture
        if(gel->HasSubElement()){continue;}

        // Check rib
        TPZManVector<REAL,3> intpoint(3,0);
        bool resul = fFracplane->Check_rib(gel, &intpoint);

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













void DFNFracture::GetOuterLoop(std::vector<int> &loop){
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    // Compute fracture edges' lenghts
    int nedges = fFracplane->GetCornersX().Cols();
    TPZVec<REAL> edgelength(nedges,0);
    Matrix fraccorners(3,nedges);
    fraccorners = fFracplane->GetCornersX();
    for(int i = 0; i < nedges; i++){
        edgelength[i] = sqrtl(pow(fraccorners(0,i)-fraccorners(0,(i+1)%nedges),2)
                              +pow(fraccorners(1,i)-fraccorners(1,(i+1)%nedges),2)
                              +pow(fraccorners(2,i)-fraccorners(2,(i+1)%nedges),2));
    }
    
	// vector of pointers to maps
    // @todo refactor this to tpzautopointer<map> to prevent memory leak
	// TPZManVector<std::map<REAL, int64_t>* > edgemap(nedges);
	TPZManVector<TPZAutoPointer<std::map<REAL, int64_t>> > edgemap(nedges);
    for(int i = 0; i < nedges; i++){
        edgemap[i] = new std::map<REAL, int64_t>;
    }
    // The set of points that have already been checked
    std::set<int64_t> checked;

    // iterate over all endfaces and map it to the fracture-edge that cuts it
    for (auto it = fFaces.begin(); it != fFaces.end(); it++){
        // get intersection node index and coordinates
        DFNFace *iface = &it->second;
        // @todo iface->IsOnBoundary2
        TPZManVector<REAL,3> ipointcoord = iface->IntersectionCoord();
        if(ipointcoord.size()<2) continue;
        int64_t ipointindex = iface->IntersectionIndex();

        // check if point was already checked
        auto aux = checked.insert(ipointindex);
        bool already_checked = !aux.second;
        if(already_checked) continue;

        // @todo if(ipointindex == any of the fracplane.CornerIndex) continue;
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
			if(dist<gDFN_SmallNumber){
				// compute local 1D coordinate (alpha)
				REAL norm = sqrtl(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
				REAL alpha = norm/edgelength[iedge];
				if(alpha > 1+gDFN_SmallNumber){std::cout<<"\nMisattribution of point to edge\n";DebugStop();}
				// map intersection indexes from smallest alpha to biggest
                // map <alpha, index>
				edgemap[iedge]->insert({alpha, ipointindex});
			}
        }
    }
    bool warning_message_Q = true;
    //Once intersections on fracture-edges have been properly ordered and mapped by edge
	//iterate over edges to split them
	for (int iedge = 0; iedge < nedges; iedge++)
	{
		int64_t nels = gmesh->NElements();
		TPZManVector<int64_t,2> inodes(2);     //index of nodes to be connected
        int icorner = iedge; //for readability
        int forshow = DFNMaterial::Efracture;

        if(edgemap[iedge]->size() == 0){
            #ifdef PZDEBUG
                if(warning_message_Q){
                    PZError<<"\n Warning: Is there an edge of a fracture that doesn't cut any element? \n\n";
                    warning_message_Q = false;
                }
            #endif //PZDEBUG
            // DebugStop();
            inodes[0] = fFracplane->CornerIndex(icorner);
            inodes[1] = fFracplane->CornerIndex((icorner+1)%nedges);
            gmesh->CreateGeoElement(EOned, inodes, forshow, nels);
            loop.push_back((int) nels);
            continue;
        }


		// connect first end-face intersection to iedge's first node
		auto it = edgemap[iedge]->begin();
        inodes[0] = fFracplane->CornerIndex(icorner);
		inodes[1] = it->second;
		gmesh->CreateGeoElement(EOned, inodes, forshow, nels);
        loop.push_back((int) nels);
        
        // iterate over iedge's map
        while(++it != edgemap[iedge]->end()){
            nels++;
            inodes[0] = inodes[1];
            inodes[1] = it->second;
            gmesh->CreateGeoElement(EOned, inodes, forshow, nels);
            loop.push_back((int) nels);
        }

		// connect last end-intersection to edge last node
		nels++;
		inodes[0] = inodes[1];
        inodes[1] = fFracplane->CornerIndex((icorner+1)%nedges);
		gmesh->CreateGeoElement(EOned, inodes, forshow, nels);
        loop.push_back((int) nels);
    }
	
    gmesh->BuildConnectivity();
    // correct duplicates
    int nlines = loop.size();
    TPZGeoEl* gel;
    for(int iline=0; iline<nlines; iline++){
        int iedge = loop[iline];
        gel = gmesh->Element(iedge);
        if(!gel || gel->Dimension() != 1) DebugStop();
        TPZGeoElSide gelside(gel,2);
        TPZGeoElSide neig = gelside.Neighbour();
        bool newedge_Q = true;
        for(/*void*/; neig!=gelside; neig = neig.Neighbour()){
            if(neig.Element()->Dimension() != 1) continue;
            // if duplicate, replace in loop and delete duplicate
            newedge_Q = false;
            loop[iline] = neig.Element()->Index(); 
            gmesh->DeleteElement(gel,iedge);
        }
        // if it's new, add it to fOutline
        if(newedge_Q){fOutline.insert({iedge,gel});}
    }

    // fix orientation
    for(int iline=1; iline<nlines; iline++){
        int iline_index = loop[iline];
        gel = gmesh->Element(iline_index);
        int anterior_index = loop[(iline+nlines-1)%nlines];
        TPZGeoEl *anterior = gmesh->Element(anterior_index);
        if(gel->NodeIndex(0) != anterior->NodeIndex(1)){
            loop[iline] *= -1;
        }
    }
}













void DFNFracture::GetFacesInSurface(std::vector<TPZGeoEl*> &faces){
    faces.reserve(20);
    std::map<int64_t, int> candidates;
    TPZGeoEl *cand_face;
    int nedges;
    TPZGeoEl *edge;
    TPZGeoMesh* gmesh = fdfnMesh->Mesh();
    // loop through edges in outline and check their 2D neighbours
    for(auto itr = fOutline.begin(); itr != fOutline.end(); itr++){
        edge = itr->second;
        TPZGeoElSide gelside(edge,2);
        TPZGeoElSide neig = gelside.Neighbour();
        for(/*void*/;neig!=gelside; neig = neig.Neighbour()){
            cand_face = neig.Element();
            if(cand_face->Dimension() != 2) continue;
            nedges = cand_face->NCornerNodes();
            int64_t index = cand_face->Index();
            if(++candidates[index] < nedges) continue;
            faces.push_back(gmesh->Element(index));
            SetMaterialIDChildren(DFNMaterial::Efracture,cand_face);
        }
    }
}














void DFNFracture::MeshFractureSurface(){
    // GMsh does not accept zero index entities
    const int shift = 1;
    // First construct the edges of the fracture surface
    std::vector<int> outerLoop;
    GetOuterLoop(outerLoop);
    std::vector<TPZGeoEl*> facesInSurface;
    GetFacesInSurface(facesInSurface);
    return;
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
            int ncorners = fFracplane->GetCornersX().Cols();
            for(int i = 0; i<ncorners; i++){
                int64_t inode = fFracplane->CornerIndex(i);
                TPZManVector<REAL, 3> coord(3);
                fdfnMesh->Mesh()->NodeVec()[inode].GetCoordinates(coord);

                gmsh::model::geo::addPoint(coord[0],coord[1],coord[2],0.,inode+shift);
            }
        }
    // INSERT LINES
        std::vector<int> curvesInSurface;
        for(auto iter = fSurface.begin(); iter != fSurface.end(); iter++){
            int64_t iel = iter->first+shift;
            TPZGeoEl *gel = iter->second;
            if(gel->Dimension() != 1) continue;
            int64_t node0 = gel->NodeIndex(0)+shift;
            int64_t node1 = gel->NodeIndex(1)+shift;

            gmsh::model::geo::addLine(node0,node1,iel);
            gmsh::model::geo::mesh::setTransfiniteCurve(iel,2); // to reduce number of nodes created by GMsh
            bool lineIsInEdge = (std::find(outerLoop.begin(),outerLoop.end(),iel) != outerLoop.end());
            if(lineIsInEdge == false){
                curvesInSurface.push_back(iel);
            }
        }
    
    // CURVE LOOP
        int surfaceIndex = 1;
        gmsh::model::geo::addCurveLoop(outerLoop,surfaceIndex);

    // SURFACE
        gmsh::model::geo::addPlaneSurface({surfaceIndex},surfaceIndex);
    
    // lines in surface
        // @comment GMsh requires synchronize before embedding geometric entities
        gmsh::model::geo::synchronize();
        gmsh::model::mesh::embed(1,curvesInSurface,2,surfaceIndex);
    // PHYSICAL GROUPS
        // physical curve
        int nlines = curvesInSurface.size() + outerLoop.size();
        if(curvesInSurface.size() > outerLoop.size()){
            curvesInSurface.reserve(nlines);
            curvesInSurface.insert(curvesInSurface.end(), outerLoop.begin(), outerLoop.end() );

            gmsh::model::addPhysicalGroup(1,curvesInSurface,DFNMaterial::Efracture);
        }else{
            outerLoop.reserve(nlines);
            outerLoop.insert(outerLoop.end(), curvesInSurface.begin(), curvesInSurface.end() );

            gmsh::model::addPhysicalGroup(1,outerLoop,DFNMaterial::Efracture);
        }
        // physical surface
        gmsh::model::addPhysicalGroup(2,{surfaceIndex},DFNMaterial::Efracture);

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
        int ncorners = fFracplane->GetCornersX().Cols();
        for(int i = 0; i<ncorners; i++){
            int pznode = (int) fFracplane->CornerIndex(i) +shift;
            
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
        TPZGeoEl *frame_el = nullptr;
        int64_t frame_el_index = -1;
        int nchildren = face_gel->NSubElements();
        if(nchildren == 0){nchildren++;}
        // queue all possible lines by checking 1D neighbours of children
		std::set<TPZGeoEl *> candidate_ribs;
		for(int ichild=0; ichild<nchildren; ichild++){
            TPZGeoEl* child;
            if(face_gel->HasSubElement()){
                child = face_gel->SubElement(ichild);
            }else{
                child = face_gel;
            }
			int nribs = child->NCornerNodes();
			for(int cside = nribs; cside < 2*nribs; cside++){
				TPZGeoElSide childside(child,cside);
				TPZGeoElSide neig = childside.Neighbour();
				for(/*void*/; neig != childside; neig = neig.Neighbour()){
					if(neig.Element()->Dimension() != 1) continue;
					candidate_ribs.insert(neig.Element());
                    break;
				}
			}
		}
        bool enters_Q = false;
        int nframenodes = framenodes.size();
        int counter = 0;
        for(auto candidate : candidate_ribs){
            for(int inode=0; inode<nframenodes; inode++){
                if(framenodes[inode]==candidate->NodeIndex(0)){
                    for(int nextnode=0; nextnode<nframenodes; nextnode++){
                        if(framenodes[nextnode]==candidate->NodeIndex(1)){
                            enters_Q = true;
                            break;
                        }
                    }
                }
                if(enters_Q) break;
            }
            if(!enters_Q) continue;
            frame_el = candidate;
            frame_el_index = candidate->Index();
            fOutline.insert({frame_el_index,frame_el});
            frame_el->SetMaterialId(DFNMaterial::Efracture);
            enters_Q = false;
            counter++;
            if(counter==2) break; //2 ribs is already an exception... 3 should never happen
            // if(counter==framenodes.size()-1) break;
        }
    }
}