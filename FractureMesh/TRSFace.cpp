/*! 
 *	TRSFace.cpp
 *  @authors   Pedro Lima
 *  @authors   Jorge Ordoñez
 *  @date      2018-2019
 */

#include "TRSFace.h"
#include "TRSFractureMesh.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"

// Empty constructor
TRSFace::TRSFace(){
    fFaceIndex = -1;
    fIsCut = false;
}

//Constructor
TRSFace::TRSFace(int64_t index, bool iscut){
    fFaceIndex = index;
    fIsCut = iscut;
}

/// Copy constructor
TRSFace::TRSFace(const TRSFace &copy){
    fFaceIndex = copy.fFaceIndex;
    fIsCut = copy.fIsCut;
    fStatus = copy.fStatus;
    fRibs = copy.fRibs;
    fSubFaces = copy.fSubFaces;
    fIntersection = copy.fIntersection;
	fFracMesh = copy.fFracMesh;
	fFather = copy.fFather;
}

/// Assignment operator
TRSFace &TRSFace::operator=(const TRSFace &copy){
    fFaceIndex = copy.fFaceIndex;
    fIsCut = copy.fIsCut;
    fStatus = copy.fStatus;
    fRibs = copy.fRibs;
    fSubFaces = copy.fSubFaces;
    fIntersection = copy.fIntersection;
	fFracMesh = copy.fFracMesh;
	fFather = copy.fFather;
    return *this;
}

/// Define the element index and whether it cuts the plane
void TRSFace::SetElementIndex(int64_t elindex, bool iscut){
    fFaceIndex = elindex;
    fIsCut = iscut;
}

/// Element index
int64_t TRSFace::ElementIndex() const{
    return fFaceIndex;
}

/// Intersects the plane or not
bool TRSFace::IsCut() const{
    return fIsCut;
}

/// Return the subelement indices
TPZVec<int64_t> TRSFace::SubElements() const{
    return fSubFaces;
}

/// Set the subelement indices
void TRSFace::SetChildren(TPZVec<int64_t> subels){
    fSubFaces = subels;
}











/// Divide the given this surface, generate subelements and return vector with indices
void TRSFace::DivideSurface(int matid){
	TPZGeoMesh *gmesh = fFracMesh->GetGeoMesh();
	TPZGeoEl *face = gmesh->Element(fFaceIndex);
	// check if face exists in mesh
   	if(!face){DebugStop();}
	
	// Vector of face's node indices
	TPZVec<int64_t> node;
	face->GetNodeIndices(node);
	
	// Determine pattern of refinement
	int splitcase = GetSplitPattern(fStatus);
	
	// Vector of children elements
	TPZVec<TPZVec<int64_t>> child;
	// Vector of refinement elements (nodes numbered as in master element)
	std::map<int64_t, int64_t> refnode;
	for(int i = 0; i<node.size(); i++) refnode.insert({node[i], i});

	// Divide surface according to split pattern (these algorithms make no sense without documentation)
	int64_t nodeA, nodeB;
	switch(splitcase){
		case 1:{
			child.resize(2);
			int i = 0;
			while(fStatus[i+4]==false){i++;}
			// Intersection node at rib i
			TRSRibs *ribA = fFracMesh->Rib(fRibs[i]);
			nodeA = gmesh->Element(ribA->IntersectionIndex())->NodeIndex(0);
			// Intersection node at rib opposite to rib i
			TRSRibs *ribB = fFracMesh->Rib(fRibs[i+2]);
			nodeB = gmesh->Element(ribB->IntersectionIndex())->NodeIndex(0);
			
			TPZVec<int64_t> child1(4);
				child1[0] = nodeB;
				child1[1] = nodeA;
				child1[2] = node[i+1];
				child1[3] = node[i+2];
			TPZVec<int64_t> child2(4);
				child2[0] = nodeA;
				child2[1] = nodeB;
				child2[2] = node[(i+3)%4];
				child2[3] = node[i];
			
			child[0] = child1;
			child[1] = child2;
			break;}
		case 2:{
			child.resize(4);
			int i = 0;
			while(fStatus[i+4]==false || fStatus[(i+3)%4+4]==false){i++;}
			// Intersection node at rib clockwise adjacent to rib i
			TRSRibs *ribA = fFracMesh->Rib(fRibs[(i+3)%4]);
			nodeA = gmesh->Element(ribA->IntersectionIndex())->NodeIndex(0);
			// Intersection node at rib i
			TRSRibs *ribB = fFracMesh->Rib(fRibs[i]);
			nodeB = gmesh->Element(ribB->IntersectionIndex())->NodeIndex(0);

			TPZVec<int64_t> child1(3);
				child1[0] = nodeA;
				child1[1] = nodeB;
				child1[2] = node[i];
			TPZVec<int64_t> child2(3);
				child2[0] = nodeB;
				child2[1] = node[(i+1)%4];
				child2[2] = node[(i+2)%4];
			TPZVec<int64_t> child3(3);
				child3[0] = nodeA;
				child3[1] = nodeB;
				child3[2] = node[(i+2)%4];
			TPZVec<int64_t> child4(3);
				child4[0] = nodeA;
				child4[1] = node[(i+2)%4];
				child4[2] = node[(i+3)%4];

			child[0] = child1;
			child[1] = child2;
			child[2] = child3;
			child[3] = child4;
			break;}
		case 3:{
			break;}
		case 4:{
			break;}
		case 5:{
			child.resize(5);
			int i = 0;
			while(fStatus[i+4]==false){i++;}
			// Intersection node at rib i
			TRSRibs *ribA = fFracMesh->Rib(fRibs[i]);
			nodeA = gmesh->Element(ribA->IntersectionIndex())->NodeIndex(0);
			// In-plane itersection node
			nodeB = gmesh->Element(fIntersection)->NodeIndex(0);

			TPZVec<int64_t> child1(3);
				child1[0] = nodeB;
				child1[1] = nodeA;
				child1[2] = node[(i+1)%4];
			TPZVec<int64_t> child2(3);
				child2[0] = nodeB;
				child2[1] = node[(i+1)%4];
				child2[2] = node[(i+2)%4];
			TPZVec<int64_t> child3(3);
				child3[0] = nodeB;
				child3[1] = node[(i+2)%4];
				child3[2] = node[(i+3)%4];
			TPZVec<int64_t> child4(3);
				child4[0] = nodeB;
				child4[1] = node[(i+3)%4];
				child4[2] = node[i];
			TPZVec<int64_t> child5(3);
				child5[0] = nodeB;
				child5[1] = node[i];
				child5[2] = nodeA;

			child[0] = child1;
			child[1] = child2;
			child[2] = child3;
			child[3] = child4;
			child[4] = child5;
			break;}
		case 6:{
			break;}
		case 7:{
			break;}
		case 8:{
			break;}
		case 9:{
			break;}
		case 10:{
			child.resize(2);
			int i = 0;
			while(fStatus[i+3]==false || fStatus[(i+2)%3+3]==false){i++;}
			// Intersection node at rib right before (counter-clockwise) to rib i
			TRSRibs *ribA = fFracMesh->Rib(fRibs[(i+2)%3]);
			nodeA = gmesh->Element(ribA->IntersectionIndex())->NodeIndex(0);
			// Intersection node at rib opposite to rib i
			TRSRibs *ribB = fFracMesh->Rib(fRibs[i]);
			nodeB = gmesh->Element(ribB->IntersectionIndex())->NodeIndex(0);
			
			TPZVec<int64_t> child1(3);
				child1[0] = nodeB;
				child1[1] = nodeA;
				child1[2] = node[i];
			TPZVec<int64_t> child2(4);
				child2[0] = nodeA;
				child2[1] = nodeB;
				child2[2] = node[(i+1)%3];
				child2[3] = node[(i+2)%3];

			child[0] = child1;
			child[1] = child2;
			break;}
		case 11:{
			break;}
		case 12:{
			child.resize(4);
			int i = 0;
			while(fStatus[i+3]==false){i++;}
			// Intersection node at rib i
			TRSRibs *ribA = fFracMesh->Rib(fRibs[i]);
			nodeA = gmesh->Element(ribA->IntersectionIndex())->NodeIndex(0);
			// In-plane itersection node
			nodeB = gmesh->Element(fIntersection)->NodeIndex(0);

			TPZVec<int64_t> child1(3);
				child1[0] = nodeB;
				child1[1] = nodeA;
				child1[2] = node[(i+1)%3];
			TPZVec<int64_t> child2(3);
				child2[0] = nodeB;
				child2[1] = node[(i+1)%3];
				child2[2] = node[(i+2)%3];
			TPZVec<int64_t> child3(3);
				child3[0] = nodeB;
				child3[1] = node[(i+2)%3];
				child3[2] = node[i];
			TPZVec<int64_t> child4(3);
				child4[0] = nodeB;
				child4[1] = node[i];
				child4[2] = nodeA;

			child[0] = child1;
			child[1] = child2;
			child[2] = child3;
			child[3] = child4;
			// DebugStop();
			break;}
		case 13:{
			break;}
		case 14:{
			break;}
		case 15:{
			break;}
		case 16:{
			break;}
		default: DebugStop();
	}

	int nchildren = child.size();
	TPZManVector<int64_t,6> childrenIndices(nchildren,0);
	childrenIndices.Shrink();

// defining Refinement Pattern
//----------------------------------------------------------------------------
	
	refnode.insert({nodeA,node.size()});
	refnode.insert({nodeB,node.size()+1});
	// set mesh to define refinement pattern
	TPZGeoMesh refPatternMesh;
	// count number of nodes for refinement pattern
	int refNNodes = face->NCornerNodes();
	if(this->fIntersection >= 0) refNNodes++;
	int ncornersfather = face->NCornerNodes();
	for(int irib=ncornersfather;irib<fStatus.size();irib++){
		if(this->fStatus[irib]) refNNodes++;
	}

	// insert nodes
	refPatternMesh.NodeVec().Resize(refNNodes);
	for(auto itr = refnode.begin(); itr != refnode.end(); itr++){
		if(itr->second == refNNodes) break;
		TPZManVector<REAL,3> coord(3);
		int64_t meshnode = itr->first;
		int64_t refnode = itr->second;
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
		for(int k = 0; k<ncorners; k++) refchild[k] = refnode[child[i][k]];
		// for colorful printing use matid+3*i or something like that
		refPatternMesh.CreateGeoElement(elemtype, refchild, matid, index);
		childrenIndices[i] = index;
	}
	refPatternMesh.BuildConnectivity();
	// int problematicface = 117;
	// if(fFaceIndex == problematicface){
	// 	std::ofstream meshprint3("meshprint3.txt");
	// 	refPatternMesh.Print(meshprint3);
	// } 
	// define refPattern
	TPZRefPattern *refpat = new TPZRefPattern(refPatternMesh);
	TPZAutoPointer<TPZRefPattern> patternPointer(refpat);
	face->SetRefPattern(patternPointer);
	
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
		newface->SetFather(fFaceIndex);
	}
	// SetChildren(childrenIndices);
	// create skeleton?
}
















/**
 * @brief Returns the split pattern that should be used to split this face
 * @param Status vector (boolean) that indicates which ribs and/or nodes are cut
 * @return Integer that indicates which split pattern to use. (check documentation)
 */
int TRSFace::GetSplitPattern(TPZVec<bool> &status){
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
					case 0: if(fIntersection == -1){splitcase = 4;break;}
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
					case 0: if(fIntersection == -1){splitcase = 14;break;}
							else{splitcase = 12;break;}
					case 1: splitcase = 11;break;
				}
				break;
			case 2:
				splitcase = 10;
				break;
		}
	}
	return splitcase;
}
  