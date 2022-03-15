/*! 
 *  @brief     Describes a rectangular polygon from four corner points.
 *  @details   
 *  @author    Pedro Lima
 *  @date      2019
 */


#include "DFNPolygon.h"
#include <math.h>
#include "DFNMesh.h"
#include "pzvec_extras.h"



#if PZ_LOG
    static TPZLogger logger("dfn.mesh");
#endif


//Constructor
DFNPolygon::DFNPolygon(const Matrix &CornerPoints, const TPZGeoMesh* gmesh) 
:   fPointsIndex(CornerPoints.Cols(),-1), 
    fArea(-1.),
    fAxis(3,3,0.)
{
  //If data is consistent, fAxis was computed during consistency check
  fCornerPoints = CornerPoints;
  ComputeAxis();
  Check_Data_Consistency();
   ComputeArea();
  // initialize the fNodesAbove data structure
  SortNodes(gmesh);
}

DFNPolygon::DFNPolygon(const Matrix &CornerPoints) {
    fCornerPoints = CornerPoints;
}

// Copy constructor
DFNPolygon::DFNPolygon(const DFNPolygon &copy){
    this->operator=(copy);
}



 DFNPolygon &DFNPolygon::operator=(const DFNPolygon &copy)
 {
	fCornerPoints = copy.fCornerPoints;
	fAxis = copy.fAxis;
	fArea = copy.area();
	fPointsIndex = copy.fPointsIndex;
    fNodesAbove = copy.fNodesAbove;
    fMax_component = copy.fMax_component;
	return *this;
 }

void DFNPolygon::ComputeAxis()
{
    int cols = fCornerPoints.Cols();
    int rows = fCornerPoints.Rows();
    
    if(rows != 3){    //Should be 3 (x y z)
        std::cout<<"Check the input data. Poorly defined fracture\n";
        DebugStop();
    }
    if(cols < 3){
        std::cout<<"Check the input data (number of corner points for a fracture seems to be less than 3)\n";
        DebugStop();
    }

    //Ax0
    fAxis(0,0) = fCornerPoints(0,0) - fCornerPoints(0,1);
    fAxis(1,0) = fCornerPoints(1,0) - fCornerPoints(1,1);
    fAxis(2,0) = fCornerPoints(2,0) - fCornerPoints(2,1);

    //Ax1
    fAxis(0,1) = fCornerPoints(0,2) - fCornerPoints(0,1);
    fAxis(1,1) = fCornerPoints(1,2) - fCornerPoints(1,1);
    fAxis(2,1) = fCornerPoints(2,2) - fCornerPoints(2,1);

    //Ax2
    fAxis(0,2) = fAxis(1,0)*fAxis(2,1) - fAxis(2,0)*fAxis(1,1);
    fAxis(1,2) = fAxis(2,0)*fAxis(0,1) - fAxis(0,0)*fAxis(2,1);
    fAxis(2,2) = fAxis(0,0)*fAxis(1,1) - fAxis(1,0)*fAxis(0,1);
    
    //Ax2 norm
    REAL norm = sqrt(fAxis(0,2)*fAxis(0,2)+fAxis(1,2)*fAxis(1,2)+fAxis(2,2)*fAxis(2,2));
    
    if(norm < 1.e-8){
        PZError << "\n[FATAL] User defined polygon seems to be inconsistent."
                << " Your first 3 corners might be co-linear nodes.\n"
                << "Corner nodes:\n"
                << fCornerPoints << '\n';
        DebugStop();
    }
    
    for (int i = 0; i < 3; i++) // i< axis number (x y z)
    {
        fAxis(i, 2) = fAxis(i, 2)/norm;
    }
    Matrix orth(3,3,0.);
    Matrix transf(3,3,0.);
    fAxis.GramSchmidt(orth,transf);
    fAxis = orth;

    // Get max absolute component of normal vector
    {
        REAL abs_x = fabs(fAxis(0,2));
        REAL abs_y = fabs(fAxis(1,2));
        REAL abs_z = fabs(fAxis(2,2));

        fMax_component = 2;
        if (abs_x > abs_y) {
            if (abs_x > abs_z) fMax_component = 0;
        }
        else if (abs_y > abs_z) fMax_component = 1;
    }
}



bool DFNPolygon::Check_Data_Consistency() const
{
  
  int cols = fCornerPoints.Cols();
  //Coplanarity verification for polygon
  if(cols > 3){
    for(int ic = 3; ic<cols; ic++) {
      //scalar product between Ax2 and <P(ic)-P1> should be zero
      
      // 1) Compute the normalized vector <P(ic)-P1>
      TPZManVector<REAL,3> vecToTest(3,0.);
      REAL normvec = 0.;
      for (int idim = 0; idim < 3; idim++) {
        vecToTest[idim] = fCornerPoints.g(idim,ic) - fCornerPoints.g(idim,1);
        normvec += vecToTest[idim]*vecToTest[idim];
      }
      normvec = sqrt(normvec);
      for (int idim = 0; idim < 3; idim++)
        vecToTest[idim] /= normvec;
      
#ifdef PZDEBUG
      const REAL newnorm = sqrt(vecToTest[0]*vecToTest[0] + vecToTest[1]*vecToTest[1] + vecToTest[2]*vecToTest[2]);
      if((std::abs(newnorm) - 1.) > ZeroTolerance()){
        std::cout<<"Fracture corner points are not coplanar"<<"\n"<<std::endl;
        DebugStop();
      }
#endif

      // 2) Perform dot product
      REAL ver = fAxis.GetVal(0,2)*vecToTest[0]
      +fAxis.GetVal(1,2)*vecToTest[1]
      +fAxis.g(2,2)*vecToTest[2];
      
      // 3) Checks if points are coplanar
      if(std::abs(ver) > gDFN_SmallNumber){
        std::cout<<"Fracture corner points are not coplanar"<<"\n"<<std::endl;
        DebugStop();
        return false;
      }
    }
  }
  
  return true;
}

/**
 * @brief Get polygon's corner points
 * @return polygon corner coordinates
 */
const Matrix& DFNPolygon::GetCornersX() const{
    return fCornerPoints;
}

void DFNPolygon::SetCornersX(const Matrix &CornerPoints){
    if(fCornerPoints.Cols() != 0) LOGPZ_WARN(logger,"Corners of a DFNPolygon were reset. You need to re-sort nodes above and below.\n");
    fCornerPoints.Resize(CornerPoints.Rows(),CornerPoints.Cols());
    fCornerPoints = CornerPoints;
    ComputeAxis();
    Check_Data_Consistency();
    ComputeArea();
    fNodesAbove.clear();
}


/**
 * @brief Checks if a point is above or below the fracture polygon
 * @param Point vector with the euclidean coordinates
 * @return True if the point is above the fracture polygon
 * @return False if the point is below the fracture polygon
 */

bool DFNPolygon::Compute_PointAbove(const TPZVec<REAL> &point) const{
    
    //Point distance to the fracture polygon computation
        double point_distance = (point[0] - GetCornersX().g(0,1))*(axis().GetVal(0,2)) 
                                +(point[1] - GetCornersX().g(1,1))*(axis().GetVal(1,2)) 
                                +(point[2] - GetCornersX().g(2,1))*(axis().GetVal(2,2));
        if (point_distance>0){
            return true;    //If the point is above the polygon
        }
        else{
            return false;   //If the point is below the polygon
        }
}
bool DFNPolygon::IsCutByPlane(TPZGeoEl *gel, TPZManVector<REAL,3> &intersection) const{
    if(gel->Dimension() != 1) {PZError<<"\n\n This ain't no rib \n\n"; DebugStop();}
    // Get rib's vertices
    TPZManVector<int64_t,2> inode(2,0);
    gel->GetNodeIndices(inode);
    TPZManVector<REAL,3> node0(3,0);
    TPZManVector<REAL,3> node1(3,0);
    gel->NodePtr(0)->GetCoordinates(node0);
    gel->NodePtr(1)->GetCoordinates(node1);
    
    //check for infinite plane
    if(IsPointAbove(inode[0]) != IsPointAbove(inode[1])){
        //Rib cut by infinite plane
        //then calculate intersection point and check if it's within polygon boundaries
        intersection = CalculateIntersection(node0, node1);
        return true;
    }else{
        return false;    //Rib is not cut by polygon
    }
}

bool DFNPolygon::IsCutByPolygon(TPZGeoEl *gel, TPZManVector<REAL,3> &intersection) const{
    return IsCutByPlane(gel,intersection) && IsPointInPolygon(intersection);
}


bool DFNPolygon::Check_pair(const TPZVec<REAL>& p1, const TPZVec<REAL>& p2, TPZManVector<REAL,3> &intersection) const{
    //check for infinite plane
    if(Compute_PointAbove(p1) != Compute_PointAbove(p2)){
        //then calculate intersection point and check if it's within polygon boundaries
        intersection = CalculateIntersection(p1,p2);
        return IsPointInPolygon(intersection);
    }else{
        return false;    //Rib is not cut by polygon
    }
}

/**
 * @brief Calculates the intersection point polygon-rib
 * @param Point above the polygon (vector)
 * @param Point below the polygon (vector)
 * @return Intersecting point
 */
TPZManVector<double, 3> DFNPolygon::CalculateIntersection(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2) const
{     
    TPZManVector<double,3> Pint(3,0.);

    double orthdist_p2 = ((fCornerPoints.g(0,1)-p2[0])*(fAxis.g(0,2))) 
                        +((fCornerPoints.g(1,1)-p2[1])*(fAxis.g(1,2)))
                        +((fCornerPoints.g(2,1)-p2[2])*(fAxis.g(2,2)));
    double rib_normal_component = (p1[0]-p2[0])*fAxis.g(0,2) 
                                 +(p1[1]-p2[1])*fAxis.g(1,2) 
                                 +(p1[2]-p2[2])*fAxis.g(2,2);
    double alpha = orthdist_p2/rib_normal_component;
    for (int p=0; p<3; p++){
        Pint[p] = p2[p] + alpha*(p1[p]-p2[p]);
    }
    //  std::cout<<"Intersection point: "<<std::endl;
    //  std::cout<<Pint[0]<<std::endl;
    //  std::cout<<Pint[1]<<std::endl;
    //  std::cout<<Pint[2]<<std::endl;
        
    return Pint;
}

/**
 * @brief Checks if a coplanar point is within fracture polygon
 * @details Enumeration of corner points should follow standard PZ topology, where 
 * corner nodes are numbered counter-clockwise. (This condition will automatically be 
 * met for triangles, but not always for quadrilaterals)
 * @param Point vector with the euclidean coordinates
 * @return True if the point is within fracture polygon
 * @return False if the point is out of fracture polygon
 */


bool DFNPolygon::IsPointInPolygon(const TPZVec<REAL> &point) const
{
    const bool isUseAreaMethod = true;
    if (isUseAreaMethod) {
        int ncorners = fCornerPoints.Cols();
        REAL area = 0;
        for(int i = 0; i<ncorners; i++){
            //Define vectors from the point to a each one of a pair of corners
            TPZManVector<REAL, 3> ax1(3);
                ax1[0] = fCornerPoints.g(0,i) - point[0];
                ax1[1] = fCornerPoints.g(1,i) - point[1];
                ax1[2] = fCornerPoints.g(2,i) - point[2];
            TPZManVector<REAL, 3> ax2(3);
                ax2[0] = fCornerPoints.g(0,(i+1)%ncorners) - point[0];
                ax2[1] = fCornerPoints.g(1,(i+1)%ncorners) - point[1];
                ax2[2] = fCornerPoints.g(2,(i+1)%ncorners) - point[2];
            //Compute area of trangle outlined by these vectors
            REAL temp = pow(ax1[1]*ax2[2] - ax1[2]*ax2[1],2);
                temp += pow(ax1[2]*ax2[0] - ax1[0]*ax2[2],2);
                temp += pow(ax1[0]*ax2[1] - ax1[1]*ax2[0],2);
                      
            area += sqrtl(temp)/2;
        }
        
        // std::cout<<" ___ ";

        //If total computed area is equal to the polygon's area, then
        //point is in polygon
        return( fabs(area-fArea) < gDFN_SmallNumber );
    }
    else {
        // Projecting all points to 2D by ignoring the max absolute component of the normal vector of this polygon
        int x = (fMax_component+1)%3;   // Projected x
        int y = (x+1)%3;                // Projected y
        
        int N_intersections = 0;
        int ncorners = fCornerPoints.Cols();
        
        constexpr REAL ignore = -999.;
        TPZManVector<REAL, 3> pA(3,ignore); // First node of the i-th edge
        TPZManVector<REAL, 3> pB(3,ignore); // Second node of the i-th edge
                                            // Set a ray leaving the point and going to infinity along x axis
                                            // Count how many edges of the projected polygon are intersected by the ray
        for(int i = 0; i<ncorners; i++){
            pA[x] = fCornerPoints.g(x,i);
            pA[y] = fCornerPoints.g(y,i);
            // pA[fMax_component] = ignore;
            pB[x] = fCornerPoints.g(x,(i+1)%ncorners);
            pB[y] = fCornerPoints.g(y,(i+1)%ncorners);
            // pB[fMax_component] = ignore;
            
            // Consecutive vertices on opposite sides of the point
            if((point[y] > pB[y]) != (point[y] > pA[y])){ // Implies pA[y]-pB[y] != 0
                                                          // Parametrize intersection
                REAL alpha = (point[y] - pB[y])/(pA[y]-pB[y]);
                // If intersection point is to the right, it's hit by the ray
                N_intersections += (pB[x] + alpha*(pA[x]-pB[x]) > point[x]);
            }
        }
        // An odd number of intersections, means point is inside polygon
        return (N_intersections%2);
    }
}


const TPZVec<int64_t> &DFNPolygon::SetPointsInGeomesh(TPZGeoMesh *gmesh){
	int64_t nels = gmesh->NElements();
	int ncorners = fCornerPoints.Cols();
	fPointsIndex.resize(ncorners);
	int64_t nnodes = gmesh->NNodes();
	gmesh->NodeVec().Resize(nnodes + ncorners);
	for(int ipoint = 0; ipoint < ncorners; ipoint++){
		TPZManVector<REAL,3> coords(3);
		for(int ico = 0; ico<3; ico++){
			coords[ico] = fCornerPoints(ico,ipoint);
		}
		gmesh->NodeVec()[nnodes+ipoint].Initialize(coords,*gmesh);
		fPointsIndex[ipoint] = nnodes+ipoint;
	}
	return fPointsIndex;
}




/**
 * @brief Computes area of polygon
 * @details Enumeration of corner points should follow standard PZ topology, where 
 * corner nodes are numbered counter-clockwise (clockwise should work as well) from
 * zero to N. (This condition will automatically be met for triangles, but not 
 * always for quadrilaterals)
 */
REAL DFNPolygon::ComputeArea(){
	int npoints = fCornerPoints.Cols();
    if(npoints<3) return 0;
	
    TPZManVector<REAL,3> normal(3);
	for(int i = 0; i<3; i++){
		normal[i] = fAxis(i,2);
	}

    // select largest abs coordinate to ignore for projection
    REAL abs_x = fabs(normal[0]);
    REAL abs_y = fabs(normal[1]);
    REAL abs_z = fabs(normal[2]);

    // int  ignore;           // coord to ignore: 0=x, 1=y, 2=z
    // ignore = 2;                         // ignore z-coord
    // if (abs_x > abs_y) {
    //     if (abs_x > abs_z) ignore = 0;  // ignore x-coord
    // }
    // else if (abs_y > abs_z) ignore = 1; // ignore y-coord
    int ignore = fMax_component;

    // compute area of the 2D projection
    REAL area = 0;
    for(int i=1; i<=npoints; i++){
        area += (   fCornerPoints((ignore+1)%3,i%npoints)
                    *(  fCornerPoints((ignore+2)%3,(i+1)%npoints)
                       -fCornerPoints((ignore+2)%3,(i-1)%npoints) )      );
    }

    // scale to get area before projection
    REAL norm = sqrt( abs_x*abs_x + abs_y*abs_y + abs_z*abs_z); // norm of normal vector
    area *= (norm / (2 * normal[ignore]));
	fArea = fabs(area);
    return fArea;
}



// /**
//  * @brief Creates a geometric element for this polygon in pointed mesh
//  * @param Pointer to geometric mesh
//  * @return Index for newly created element in gmesh
//  */
// int64_t DFNPolygon::CreateElement(TPZGeoMesh *gmesh){
// 	// number of nodes for gmesh
// 	int nnodes =  gmesh->NNodes();
// 	// nomber of corners for polygon
// 	int ncorners = fCornerPoints.Cols();
// 	// add corner as nodes in gmesh
// 	gmesh->NodeVec().Resize(nnodes + ncorners);
// 	TPZVec<int64_t> CornerIndexes(ncorners);
// 	for (int i = 0; i < ncorners; i++)
// 	{
// 		TPZVec<REAL> nodeX(3, 0);
// 		nodeX[0] = fCornerPoints(0,i);
// 		nodeX[1] = fCornerPoints(1,i);
// 		nodeX[2] = fCornerPoints(2,i);

// 		gmesh->NodeVec()[nnodes + i].Initialize(nodeX, *gmesh);
// 		CornerIndexes[i] = nnodes + i;
// 	}
// 	// set element type
// 	MElementType elemtype;
// 	switch (ncorners){
// 		case 3: elemtype = ETriangle; break;
// 		case 4: elemtype = EQuadrilateral; break;
// 		default: DebugStop();
// 	}
// 	// create geometric element
// 	int64_t polygonindex = gmesh->NElements();
// 	gmesh->CreateGeoElement(elemtype, CornerIndexes, 40, polygonindex);

// 	return polygonindex;
// }



TPZManVector<REAL, 3> DFNPolygon::GetProjectedX(const TPZManVector<REAL, 3> &x)const{
    TPZManVector<REAL, 3> projection(3,0);
    REAL alpha = fAxis.g(0,2)*(x[0]-fCornerPoints.g(0,1))
                +fAxis.g(1,2)*(x[1]-fCornerPoints.g(1,1))
                +fAxis.g(2,2)*(x[2]-fCornerPoints.g(2,1));
    projection[0] = x[0] - alpha*fAxis.g(0,2);
    projection[1] = x[1] - alpha*fAxis.g(1,2);
    projection[2] = x[2] - alpha*fAxis.g(2,2);
    return projection;
}



TPZGeoEl* DFNPolygon::InsertGeoEl(TPZGeoMesh* gmesh, int matid, TPZVec<int64_t>* nodes){
    int nnodes = this->NCornerNodes();
    if(nnodes > 4) return nullptr; //pz only supports triangles and quadrilaterals as 2D topologies

    MElementType eltype = (nnodes == 4? EQuadrilateral : ETriangle);
    int64_t elindex = -1;
    TPZManVector<int64_t> newnodes(nnodes,-1);
    if(nodes){
        if(nodes->size() != nnodes) DebugStop();
    }else{
        for(int inode=0; inode<nnodes; inode++){
            int64_t nodeindex = gmesh->NodeVec().AllocateNewElement();
            newnodes[inode] = nodeindex;
            TPZManVector<REAL, 3> coord(3,0);
            coord[0] = fCornerPoints(0,inode);
            coord[1] = fCornerPoints(1,inode);
            coord[2] = fCornerPoints(2,inode);
            gmesh->NodeVec()[nodeindex].Initialize(coord,*gmesh);
        }
        nodes = &newnodes;
    }
    return gmesh->CreateGeoElement(eltype,*nodes,matid,elindex);
}

void DFNPolygon::ComputeCentroid(TPZVec<REAL>& centroid){
    centroid.resize(3);
    centroid.Fill(0.);
    int nnodes = this->NCornerNodes();
    for(int inode=0; inode<nnodes; ++inode){
        centroid[0] += fCornerPoints(0,inode);
        centroid[1] += fCornerPoints(1,inode);
        centroid[2] += fCornerPoints(2,inode);
    }
    centroid[0] /= nnodes;
    centroid[1] /= nnodes;
    centroid[2] /= nnodes;
}



TPZVec<TPZGeoEl*> DFNPolygon::InsertGeomRepresentation(TPZGeoMesh* gmesh, int matid, int exclusive_dim){
    int nnodes = this->NCornerNodes();
    
    // Create nodes
    int needscentroid = (nnodes > 4 && exclusive_dim!=1);
    TPZManVector<int64_t> polyg_nodes(nnodes+needscentroid,0);
    TPZManVector<REAL, 3> coord(3,0);
    int64_t nodeindex=-1;
    for(int inode=0; inode < nnodes; ++inode){
        nodeindex = gmesh->NodeVec().AllocateNewElement();
        polyg_nodes[inode] = nodeindex;
        coord[0] = fCornerPoints(0,inode);
        coord[1] = fCornerPoints(1,inode);
        coord[2] = fCornerPoints(2,inode);
        gmesh->NodeVec()[nodeindex].Initialize(coord,*gmesh);
    }

    int nels2D = nnodes*(exclusive_dim!=1 && nnodes>4)  // 2D elements for non-quad and non-triang
                +(nnodes <= 4);                     // 2D element if quad or triangle
    int nels1D = nnodes*(exclusive_dim!=2);             // Boundary elements
    int nels = nels1D + nels2D;
    TPZStack<TPZGeoEl*> els;
    
    
    TPZManVector<int64_t,3> el_nodes(3,-1);
    if(exclusive_dim != 1){ // 2D elements
        if(needscentroid){
            // Create centroid
            nodeindex = gmesh->NodeVec().AllocateNewElement();
            polyg_nodes[nnodes] = nodeindex;
            this->ComputeCentroid(coord);
            gmesh->NodeVec()[nodeindex].Initialize(coord,*gmesh);
            el_nodes[0] = polyg_nodes[nnodes]; ///< centroid
            // Create 2D elements
            for(int iel=0; iel<nels2D; ++iel){
                el_nodes[1] = polyg_nodes[iel];
                el_nodes[2] = polyg_nodes[(iel+1)%nnodes];
                int64_t dummyindex=-1;
                els.push_back(gmesh->CreateGeoElement(ETriangle,el_nodes,matid,dummyindex,0));
            }
        }else{
            els.push_back(InsertGeoEl(gmesh,matid));
        }
    }

    if(exclusive_dim != 2){// 1D elements at boundary
        el_nodes.resize(2);
        for(int inode=0; inode<nnodes; inode++){
            el_nodes[0] = polyg_nodes[inode];
            el_nodes[1] = polyg_nodes[(inode+1)%nnodes];
            int64_t dummyindex = -1;
            els.push_back(gmesh->CreateGeoElement(MElementType::EOned,el_nodes,matid,dummyindex,0));
        }
    }
    gmesh->BuildConnectivity();
    return els;
}


void DFNPolygon::Print(std::ostream & out) const
{
	out << "\nDFNPolygon:\n";
	out << "Number of corners             = " << NCornerNodes() << "\n";
	out << "Area                          = " << fArea << "\n";
	out << "Corner coordinates:\n";
    fCornerPoints.Print("CornersX",out,EFixedColumn);
	out << "Axes:\n";
    fAxis.Print("Axes",out,EFixedColumn);
	out << "GeoNode indices:\n\t\t";
    out << fPointsIndex;
    out << std::endl;
}



void DFNPolygon::SortNodes(const TPZGeoMesh* gmesh){
    LOGPZ_INFO(logger,"[Start][Sorting nodes to either side of the plane]");
    LOGPZ_DEBUG(logger, "Nodes above before = \n\t{" << fNodesAbove << '}')


    std::cout << " -Sorting nodes to either side of the plane\r" << std::flush;
    int64_t i = fNodesAbove.size();
    const int64_t nnodes = gmesh->NNodes();
    fNodesAbove.Resize(nnodes,false);
    TPZManVector<REAL,3> coord(3,0.0);
    for(/*void*/; i<nnodes; i++){
        TPZGeoNode& node = gmesh->NodeVec()[i];
        if(node.Id() < 0) continue;
        node.GetCoordinates(coord);
        fNodesAbove[i] = Compute_PointAbove(coord);
    }
    std::cout << "                                           \r" << std::flush;

    LOGPZ_DEBUG(logger, "Nodes above after = \n\t{" << fNodesAbove << '}');
    LOGPZ_INFO(logger,"[End][Sorting nodes to either side of the plane]");

}

void DFNPolygon::PlotNodesAbove_n_Below(TPZGeoMesh* gmesh){
    SortNodes(gmesh);

    TPZManVector<int64_t,1> nodeindices(1,-1);
    const int64_t nnodes = gmesh->NNodes();
    for(int64_t i=0; i<nnodes; i++){
        if(gmesh->NodeVec()[i].Id() < 0) continue;
        nodeindices[0] = i;
        int64_t index;
        gmesh->CreateGeoElement(EPoint,nodeindices,fNodesAbove[i],index);
    }
}





REAL DFNPolygon::EdgeLength(int edgeindex) const{
    // consistency checks
    // if(edgeindex < 0) DebugStop();
    // if(edgeindex >= NEdges()) DebugStop();
    TPZManVector<REAL,3> dif(3,0.);
    GetEdgeVector(edgeindex,dif);
    return DFN::Norm<REAL>(dif);
}



void DFNPolygon::GetEdgeVector(int edgeindex, TPZVec<REAL>& edgevector) const{
    // consistency checks
    if(edgeindex < 0) DebugStop();
    if(edgeindex >= NEdges()) DebugStop();

    int node0 = edgeindex;
    int node1 = (edgeindex+1)%NEdges();

    edgevector[0] = fCornerPoints.g(0,node0) - fCornerPoints.g(0,node1);
    edgevector[1] = fCornerPoints.g(1,node0) - fCornerPoints.g(1,node1);
    edgevector[2] = fCornerPoints.g(2,node0) - fCornerPoints.g(2,node1);

}


bool DFNPolygon::ComputePolygonIntersection(const DFNPolygon& otherpolyg, Segment& segment) const{

    int n_int = 0; //< Number of intersection points found
    segment.clear();
    
    const bool isUsePlaneCrossMethod = false;
    if (isUsePlaneCrossMethod) {
        const DFNPolygon& polyg1 = *this;
        const DFNPolygon& polyg2 = otherpolyg;
        
        TPZManVector<REAL,3> p1n,p2n,p3n,p3np2n,p1np3n,pinter;
        polyg1.GetNormal(p1n);
        polyg2.GetNormal(p2n);
        Cross(p1n, p2n, p3n);
        REAL squaredLength = 0.;
        for (int i = 0; i < 3; ++i) squaredLength += p3n[i]*p3n[i];
        Cross(p3n, p2n, p3np2n);
        Cross(p1n, p3n, p1np3n);
        REAL p1ap1n = 0., p2ap2n = 0.;
        for (int i = 0; i < 3; ++i) {
            p1ap1n -= polyg1.PointAndCoor(0,i)*p1n[i];
            p2ap2n -= polyg2.PointAndCoor(0,i)*p2n[i];
        }
        
        for (int i = 0; i < 3; ++i) {
            pinter[i] = p3np2n[i]*p1ap1n + p1np3n[i]*p2ap2n;
        }
        
        pinter /= squaredLength;
        
        // Now we need to find start and end based on intersection of the
        // line defined by the point (pinter) and the normal (p3n)
                
    }
    else{
        TPZManVector<REAL,3> p1(3,0.), p2(3,0.), intpoint(3,0.);
        TPZManVector<const DFNPolygon*,2> polyg {this, &otherpolyg};
        // polyg[0] = *this;
        // polyg[1] = otherpolyg;

        for(int A=0; A<2; A++){
            int B = !A;
            // Test edges of DFNPolygonA being cut by DFNPolygonB, then
            // test edges of DFNPolygonB being cut by DFNPolygonA
            int nnodesA = polyg[A]->NCornerNodes();
            for(int inode=0; inode<nnodesA; inode++){
                // Get 2 consecutive corners
                // @suggestion: This code is affected by machine precision. A further extension to this implementation would be to displace p1 and p2 from polyg[A] centroid by some geometrical tolerance (using a displacement vector from the centroid to each node). Fractures that perfectly end on another may not have their intersection properly computed, and this small displacement would fix that.
                polyg[A]->iCornerX(inode,p1);
                polyg[A]->iCornerX((inode+1)%nnodesA,p2);
                
                // Check if corners lie in opposite sides, and if the
                // intersection point is within bounds of the polygon
                bool intersects_Q = polyg[B]->Check_pair(p1,p2,intpoint);
                if(!intersects_Q) continue;

                segment.push_back(intpoint);
                n_int++;
                if(n_int == 2) return true;
            }
        }

        return false;
    }
}


void DFNPolygon::PlotVTK(const std::string filepath, const int materialID, const int64_t index) const{
    std::ofstream file(filepath);

    std::stringstream node, connectivity;

    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "DFNPolygon VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    const int elNnodes = this->NCornerNodes();
    file << "POINTS " << elNnodes << " float\n";

    connectivity << elNnodes;
    TPZManVector<REAL,3> cornerX(3,0.0);
    for(int64_t inode = 0; inode < elNnodes; inode++) {
        iCornerX(inode,cornerX);
        for (int c = 0; c < 3; c++) {
            REAL coord = cornerX[c];
            node << coord << " ";
        }
        node << '\n';
        connectivity << " " << inode;
    }
    connectivity << '\n';
    node << '\n';
    file << node.str();
    int64_t size = elNnodes+1;
    file << "CELLS 1 " << size << '\n';
    file << connectivity.str() << std::endl;

    constexpr int type = 7; // Arbitrary polygon
    constexpr int eldimension = 2; 

    file << "CELL_TYPES 1\n7\n";
    file << "CELL_DATA 1\n";
    file << "FIELD FieldData 3\n";
    file << "material 1 1 int\n" << materialID << '\n';
    file << "elIndex 1 1 int\n" << index << '\n';
    file << "Dimension 1 1 int\n2";
    file.close();
}

using namespace DFN;

const REAL DFNPolygon::GetWorstAngleCos() {
    REAL area = 0.;
    TPZManVector<REAL,3> centroid(3,0.);
    ComputeCentroid(centroid);
    const int npts = fCornerPoints.Cols();
    for (int i = 0; i < npts; i++) {
        area += SubTriangleArea(centroid, i);
    }
    
    TPZManVector<REAL,3> crossprod;
    TPZStack<TPZManVector<REAL,3>> allnormals;
    const REAL eps = 1.e-2;
    for (int i = 0; i < npts; i++) {
        for (int j = i+1; j < npts-2; j++) {
            for (int k = j+1; k < npts; k++) {
                crossprod = GetCrossProduct(i,j,k);
                const REAL norm = DFN::Norm<REAL>(crossprod);
                if (norm/2. < eps*area)
                    continue;
                                
                crossprod /= norm;
                allnormals.Push(crossprod);
            }
        }
    }
    
    REAL worstAngleCos = 1.;
    for (int i = 0; i < allnormals.size(); i++) {
        for (int j = i+1; j < allnormals.size(); j++) {
            const REAL cosangle = DFN::DotProduct<REAL>(allnormals[i], allnormals[j]);
            if (cosangle < worstAngleCos){
                worstAngleCos = cosangle;
            }
        }
    }
    return worstAngleCos;
}


const REAL DFNPolygon::SubTriangleArea(const TPZVec<REAL>& centroid, const int i) const {
    const int npts = fCornerPoints.Cols();
    const int inext = (i+1)%npts;
    TPZManVector<REAL,3> vec0 = TPZVec<REAL>({fCornerPoints(0,i),fCornerPoints(1,i),fCornerPoints(2,i)});
    vec0 -= centroid;
    TPZManVector<REAL,3> vec1 = {fCornerPoints(0,inext),fCornerPoints(1,inext),fCornerPoints(2,inext)};
    vec1 -= centroid;
    
    const REAL area = DFN::Norm<REAL>(DFN::CrossProduct<REAL>(vec0, vec1))/2.;
    
    return area;
}

TPZManVector<REAL,3> DFNPolygon::GetCrossProduct(const int i, const int j, const int k) {
    TPZManVector<REAL,3> vec0 = TPZVec<REAL>({fCornerPoints(0,i),fCornerPoints(1,i),fCornerPoints(2,i)});
    for (int idim = 0; idim < 3; idim++) vec0[idim] = vec0[idim] - fCornerPoints(idim,j);
    TPZManVector<REAL,3> vec1 = {fCornerPoints(0,k),fCornerPoints(1,k),fCornerPoints(2,k)};
    for (int idim = 0; idim < 3; idim++) vec1[idim] = vec1[idim] - fCornerPoints(idim,j);

    TPZManVector<REAL,3> crossprod = DFN::CrossProduct<REAL>(vec0, vec1);
    
    return crossprod;
}
