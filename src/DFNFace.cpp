/*! 
 *	DFNFace.cpp
 *  @authors   Pedro Lima
 *  @authors   Philippe Devloo
 *  @date      2018-2020
 */

#include "DFNFace.h"
#include "DFNFracture.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"


//Constructor
DFNFace::DFNFace(TPZGeoEl *gel, DFNFracture *fracture) :
	fGeoEl(gel),
	fFracture(fracture)
{
	int nsides = gel->NSides();
	fStatus.Resize(nsides,0);
}

/// Copy constructor
DFNFace::DFNFace(const DFNFace &copy){
    this->operator=(copy);
}

/// Assignment operator
DFNFace &DFNFace::operator=(const DFNFace &copy){
    fGeoEl = copy.fGeoEl;
	fStatus = copy.fStatus;
	fCoord = copy.fCoord;
	fIntersectionIndex = copy.fIntersectionIndex;
	fFracture = copy.fFracture;
	fRibs = copy.fRibs;
    fRefMesh = copy.fRefMesh;
    return *this;
}













/// Check if should be refined and generate the subelements of material id matID
void DFNFace::Refine(){
	TPZGeoEl *face = fGeoEl;
	int matid = fGeoEl->MaterialId();
	TPZGeoMesh *gmesh = face->Mesh();
	// check if face exists in mesh
   	if(!face){DebugStop();}
	
	// Vector of face's node indices
	TPZManVector<int64_t,4> node;
	face->GetNodeIndices(node);
	
	// Determine pattern of refinement
	int splitcase = GetSplitPattern(fStatus);
	
	// Vector of children elements
	TPZManVector<TPZManVector<int64_t,4>,6> child;
	// Maps geometric mesh index to refinement mesh index
	std::map<int64_t, int64_t> gmesh_to_rmesh;
	for(int i = 0; i<node.size(); i++) gmesh_to_rmesh.insert({node[i], i});

	// Divide surface according to split pattern (these algorithms make no sense without documentation)
	int64_t nodeA, nodeB;
	switch(splitcase){
		case 1:{
			child.resize(2);
			int i = 0;
			while(fStatus[i+4]==false){i++;}
			// Intersection node at rib i
			DFNRib *ribA = fFracture->Rib(fRibs[i]);
			nodeA = ribA->IntersectionIndex();
			// Intersection node at rib opposite to rib i
			DFNRib *ribB = fFracture->Rib(fRibs[i+2]);
			nodeB = ribB->IntersectionIndex();
			
			child[0].Resize(4);
				child[0][0] = nodeB;
				child[0][1] = nodeA;
				child[0][2] = node[i+1];
				child[0][3] = node[i+2];
			child[1].resize(4);
				child[1][0] = nodeA;
				child[1][1] = nodeB;
				child[1][2] = node[(i+3)%4];
				child[1][3] = node[i];
			break;}
		case 2:{
			child.resize(4);
			int i = 0;
			while(fStatus[i+4]==false || fStatus[(i+3)%4+4]==false){i++;}
			// Intersection node at rib clockwise adjacent to rib i
			DFNRib *ribA = fFracture->Rib(fRibs[(i+3)%4]);
			nodeA = ribA->IntersectionIndex();
			// Intersection node at rib i
			DFNRib *ribB = fFracture->Rib(fRibs[i]);
			nodeB = ribB->IntersectionIndex();

			child[0].Resize(3);
				child[0][0] = nodeA;
				child[0][1] = nodeB;
				child[0][2] = node[i];
			child[1].resize(3);
				child[1][0] = nodeB;
				child[1][1] = node[(i+1)%4];
				child[1][2] = node[(i+2)%4];
			child[2].resize(3);
				child[2][0] = nodeA;
				child[2][1] = nodeB;
				child[2][2] = node[(i+2)%4];
			child[3].resize(3);
				child[3][0] = nodeA;
				child[3][1] = node[(i+2)%4];
				child[3][2] = node[(i+3)%4];
			break;}
		case 3:{
			child.resize(2);
			int i = fStatus[1];
			child[0].Resize(3);
				child[0][0] = node[i];
				child[0][1] = node[i+1];
				child[0][2] = node[i+2];
			child[1].resize(3);
				child[1][0] = node[i];
				child[1][1] = node[i+2];
				child[1][2] = node[(i+3)%4];
			break;}
		case 4:
		case 7:{
			child.resize(3);
			int i = 0;
			while(fStatus[i+4]==false) i++;
			// Intersection node at rib i
			DFNRib *ribA = fFracture->Rib(fRibs[i]);
			nodeA = ribA->IntersectionIndex();

			child[0].resize(3);
				child[0][0] = nodeA;
				child[0][1] = node[(i+1)%4];
				child[0][2] = node[(i+2)%4];
			child[1].resize(3);
				child[1][0] = nodeA;
				child[1][1] = node[(i+2)%4];
				child[1][2] = node[(i+3)%4];
			child[2].resize(3);
				child[2][0] = nodeA;
				child[2][1] = node[(i+3)%4];
				child[2][2] = node[i];
			break;}
		case 5:{
			child.resize(5);
			int i = 0;
			while(fStatus[i+4]==false){i++;}
			// Intersection node at rib i
			DFNRib *ribA = fFracture->Rib(fRibs[i]);
			nodeA = ribA->IntersectionIndex();
			// In-plane itersection node
			nodeB = fIntersectionIndex;

			child[0].Resize(3);
				child[0][0] = nodeB;
				child[0][1] = nodeA;
				child[0][2] = node[(i+1)%4];
			child[1].resize(3);
				child[1][0] = nodeB;
				child[1][1] = node[(i+1)%4];
				child[1][2] = node[(i+2)%4];
			child[2].resize(3);
				child[2][0] = nodeB;
				child[2][1] = node[(i+2)%4];
				child[2][2] = node[(i+3)%4];
			child[3].resize(3);
				child[3][0] = nodeB;
				child[3][1] = node[(i+3)%4];
				child[3][2] = node[i];
			child[4].resize(3);
				child[4][0] = nodeB;
				child[4][1] = node[i];
				child[4][2] = nodeA;
			break;}
		case 6:
		// case 7 == case 4
		case 8:{ // case 8 == case 6
			// In-plane itersection node
			nodeA = fIntersectionIndex;
			child.resize(4);
			for(int i=0; i<4; i++) {
				child[i].resize(3);
				child[i][0] = nodeA;
				child[i][1] = node[i];
				child[i][2] = node[(i+1)%4];
			}
			break;}
		case 9:{
			// @ToDo
			break;}
		case 10:{
			child.resize(2);
			int i = 0;
			while(fStatus[i+3]==false || fStatus[(i+2)%3+3]==false){i++;}
			// Intersection node at rib right before (counter-clockwise) to rib i
			DFNRib *ribA = fFracture->Rib(fRibs[(i+2)%3]);
			nodeA = ribA->IntersectionIndex();
			// Intersection node at rib opposite to rib i
			DFNRib *ribB = fFracture->Rib(fRibs[i]);
			nodeB = ribB->IntersectionIndex();
			
			child[0].Resize(3);
				child[0][0] = nodeB;
				child[0][1] = nodeA;
				child[0][2] = node[i];
			child[1].resize(4);
				child[1][0] = nodeA;
				child[1][1] = nodeB;
				child[1][2] = node[(i+1)%3];
				child[1][3] = node[(i+2)%3];
			break;}
		case 11:
		case 14:{ //case 11 == case 14
			child.resize(2);
			int i=0;
			while(fStatus[i+3] == false) i++;
			// Intersection node at rib i
			DFNRib *ribA = fFracture->Rib(fRibs[i]);
			nodeA = ribA->IntersectionIndex(); 

			child[0].resize(3);
				child[0][0] = nodeA;
				child[0][1] = node[(i+1)%3];
				child[0][2] = node[(i+2)%3];
			child[1].resize(3);
				child[1][0] = nodeA;
				child[1][1] = node[(i+2)%3];
				child[1][2] = node[i];
			break;}
		case 12:{
			child.resize(4);
			int i = 0;
			while(fStatus[i+3]==false){i++;}
			// Intersection node at rib i
			DFNRib *ribA = fFracture->Rib(fRibs[i]);
			nodeA = ribA->IntersectionIndex();
			// In-plane itersection node
			nodeB = fIntersectionIndex;

			child[0].Resize(3);
				child[0][0] = nodeB;
				child[0][1] = nodeA;
				child[0][2] = node[(i+1)%3];
			child[1].resize(3);
				child[1][0] = nodeB;
				child[1][1] = node[(i+1)%3];
				child[1][2] = node[(i+2)%3];
			child[2].resize(3);
				child[2][0] = nodeB;
				child[2][1] = node[(i+2)%3];
				child[2][2] = node[i];
			child[3].resize(3);
				child[3][0] = nodeB;
				child[3][1] = node[i];
				child[3][2] = nodeA;
			break;}
		case 13: // case 13 == case 15
		// case 14 == case 11
		case 15:{
			// In-plane itersection node
			nodeA = fIntersectionIndex;
			child.resize(3);
			for(int i=0; i<3; i++) {
				child[i].resize(3);
				child[i][0] = nodeA;
				child[i][1] = node[i];
				child[i][2] = node[(i+1)%3];
			}
			break;}
		case 16:{
			// @ToDo
			break;}
		default: DebugStop();
	}

	int nchildren = child.size();
	TPZManVector<int64_t,6> childrenIndices(nchildren,0);


// defining Refinement Pattern
//----------------------------------------------------------------------------
	{
		gmesh_to_rmesh.insert({nodeA,node.size()});
		gmesh_to_rmesh.insert({nodeB,node.size()+1});
		// set mesh to define refinement pattern
		TPZGeoMesh refPatternMesh;
		// count number of nodes for refinement pattern
		int refNNodes = face->NCornerNodes();
		if(this->fIntersectionIndex >= 0) refNNodes++;
		int ncornersfather = face->NCornerNodes();
		for(int irib=ncornersfather;irib<fStatus.size();irib++){
			if(this->fStatus[irib]) refNNodes++;
		}
	
		// insert nodes
		refPatternMesh.NodeVec().Resize(refNNodes);
		for(auto itr = gmesh_to_rmesh.begin(); itr != gmesh_to_rmesh.end(); itr++){
			int64_t refnode = itr->second;
			if(refnode == refNNodes) continue;
			TPZManVector<REAL,3> coord(3);
			int64_t meshnode = itr->first;
			gmesh->NodeVec()[meshnode].GetCoordinates(coord);
			// refPatternMesh.NodeVec().AllocateNewElement();
			refPatternMesh.NodeVec()[refnode].Initialize(coord,refPatternMesh);
		}
	
		
		// insert father
		{
			MElementType elemtype = face->Type();
			TPZManVector<int64_t,4> cornerindices(ncornersfather);
			for(int i = 0; i<ncornersfather; i++) cornerindices[i] = i;
			int64_t index = 0;
			refPatternMesh.CreateGeoElement(elemtype, cornerindices, matid, index);
		}
		// insert children
		for (int i = 0; i < nchildren; i++)
		{
			int ncorners = child[i].size();
			MElementType elemtype;
			switch (ncorners){
				case 3: elemtype = ETriangle; break;
				case 4: elemtype = EQuadrilateral; break;
			}
			int64_t index = i+1;
			// int64_t index = i;
			TPZManVector<int64_t,4> refchild(ncorners);
			for(int k = 0; k<ncorners; k++) refchild[k] = gmesh_to_rmesh[child[i][k]];
			// for colorful printing use matid+3*i or something like that
			refPatternMesh.CreateGeoElement(elemtype, refchild, matid, index);
		}
		refPatternMesh.BuildConnectivity(); 
    	    //Print result
			// if(true);
				// std::ofstream out2("./TestRefMesh.vtk");
				// TPZVTKGeoMesh::PrintGMeshVTK(&refPatternMesh, out2, true);
		// define refPattern
		TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(refPatternMesh);
		face->SetRefPattern(refpat);
	}
//----------------------------------------------------------------------------

	for (int i = 0; i < nchildren; i++)
	{
		int ncorners = child[i].size();
		MElementType elemtype;
		switch (ncorners){
			case 3: elemtype = ETriangle; break;
			case 4: elemtype = EQuadrilateral; break;
		}
		int64_t index = gmesh->NElements();
		// for colorful printing use matid+3*i or something like that
		TPZGeoEl *newface = gmesh->CreateGeoElement(elemtype, child[i], matid, index);
		childrenIndices[i] = index;

		
		face->SetSubElement(i,newface);
		// Tell the child who its father is
		newface->SetFatherIndex(this->Index());
	}
	// create skeleton?
}
















/**
 * @brief Returns the split pattern that should be used to split this face
 * @param Status vector (boolean) that indicates which ribs and/or nodes are cut
 * @return Integer that indicates which split pattern to use. (check documentation)
 */
int DFNFace::GetSplitPattern(TPZManVector<int> &status){
    // Count number of ribs and nodes cut
    int ribscut = 0;
	int nodescut = 0;
	int n = (int) status.size()/2;
	for (int i = 0; i < n; i++){
		ribscut += status[i+n];
		nodescut += status[i];
	}

	// Get split case
    // Check documentation for consult on what each splitcase means
	int splitcase;
	if (n == 4){ //quadrilateral
		switch(ribscut){
			case 0:
				switch(nodescut){
					case 0: splitcase = 8;break; //or 9
					case 1: splitcase = 6;break;
					case 2: splitcase = 3;break; 
				}
				break;
			case 1:
				switch(nodescut){
					case 0: if(fIntersectionIndex == -1){splitcase = 4;break;}
							else{splitcase = 5;break;}
					case 1: splitcase = 7;break;
				}
				break;
			case 2:{
				int i=0;
				while(status[i+4]==false){i++;}
				if(status[i+5]==true||status[(i+3)%4+4]==true){splitcase = 2;}
				else splitcase = 1;
				break;}
			default: DebugStop();break;
		}
	}
	else{ //triangle
		switch(ribscut){
			case 0:
				switch(nodescut){
					case 0: splitcase = 15;break; //or 16
					case 1: splitcase = 13;break;
				}
				break;
			case 1:
				switch(nodescut){
					case 0: if(fIntersectionIndex == -1){splitcase = 14;break;}
							else{splitcase = 12;break;}
					case 1: splitcase = 11;break;
				}
				break;
			case 2:
				splitcase = 10;
				break;
			default: DebugStop();break;
		}
	}
	return splitcase;
}

// Naive vector comparison
template <class T>
bool operator!=(const class TPZManVector<T>& v1,const class TPZManVector<T>& v2){
	int64_t size1 = v1.size();
	int64_t size2 = v2.size();
	if(size1 != size2) return true;
	for(int64_t i = 0; i<size1; i++){
		if(v1[i]!=v2[i]) return true;
	}
	return false;
}

bool DFNFace::UpdateStatusVec(){
	if(fRibs.size()<1 || fRibs[0]<0) DebugStop();
	
	TPZGeoMesh *gmesh = fGeoEl->Mesh();
	int nnodes = fGeoEl->NCornerNodes();
	int nribs = nnodes;
	DFNRib *rib;
	TPZGeoEl *rib_gel;
	
	TPZManVector<int> backup = fStatus;
	
	int orientation = 0;
	for(int i=0; i<nribs; i++){
		rib = fFracture->Rib(fRibs[i]);
		if(rib){
			rib_gel = gmesh->Element(fRibs[i]);
			orientation = (rib_gel->NodeIndex(0) == fGeoEl->NodeIndex(i)) ? 0 : 1;
			fStatus[i] = rib->StatusVec()[orientation];
			// fStatus[(i+1)%nnodes] = rib->StatusVec()[(orientation+1)%2];
			fStatus[i+nnodes] = rib->StatusVec()[2];
		}
		else{
			fStatus[i] = 0;
			// fStatus[(i+1)%nnodes] = 0;
			fStatus[i+nnodes] = 0;
		}
	}

	return backup != fStatus;
}