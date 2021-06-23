/*! 
 *  @brief     Describes a rectangular polygon from four corner points.
 *  @details   
 *  @author    Pedro Lima
 *  @date      2019
 */


#include "DFNPolygon.h"
#include <math.h>
#include "DFNMesh.h"


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
    if(fCornerPoints.Cols() != 0) LOG_DFN_WARN("Corners of a DFNPolygon were reset. You need to re-sort nodes above and below.\n");
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

    int  ignore;           // coord to ignore: 0=x, 1=y, 2=z
    ignore = 2;                         // ignore z-coord
    if (abs_x > abs_y) {
        if (abs_x > abs_z) ignore = 0;  // ignore x-coord
    }
    else if (abs_y > abs_z) ignore = 1; // ignore y-coord

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
/**
 * @brief Inserts a few geometrical elements to graphically represent this polygon. If NCorners <= 4 it'll be only one element.
*/
TPZVec<TPZGeoEl*> DFNPolygon::InsertGeomRepresentation(TPZGeoMesh* gmesh, int matid){
    int nnodes = this->NCornerNodes();
    if(nnodes <= 4) return {InsertGeoEl(gmesh,matid)};
    
    // Create nodes
    TPZManVector<int64_t> polyg_nodes(nnodes+1,0);
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
    // Create centroid
    nodeindex = gmesh->NodeVec().AllocateNewElement();
    polyg_nodes[nnodes] = nodeindex;
    this->ComputeCentroid(coord);
    gmesh->NodeVec()[nodeindex].Initialize(coord,*gmesh);
    // Create elements
    int nels = nnodes;
    TPZVec<TPZGeoEl*> els(nels,nullptr);
    TPZManVector<int64_t> el_nodes(3,-1);
    el_nodes[0] = polyg_nodes[nnodes]; ///< centroid
    for(int iel=0; iel<nels; ++iel){
        el_nodes[1] = polyg_nodes[iel];
        el_nodes[2] = polyg_nodes[(iel+1)%nnodes];
        int64_t dummyindex=-1;
        els[iel] = gmesh->CreateGeoElement(ETriangle,el_nodes,matid,dummyindex,0);
    }
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
    const DFNPolygon& polygA = *this;
    const DFNPolygon& polygB = otherpolyg;

    int n_int = 0; //< Number of intersection points found
    segment.clear();

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
