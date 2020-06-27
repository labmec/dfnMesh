/*! 
 *  @brief     Describes a rectangular plane from four corner points.
 *  @details   
 *  @author    Pedro Lima
 *  @date      2019
 */


#include "DFNFracPlane.h"
#include <math.h>

//Constructor
DFNFracPlane::DFNFracPlane(const Matrix &CornerPoints)
{
    if( !Check_Data_Consistency(CornerPoints) )	{
		std::cout<<"Error at DFNFracPlane: Bad input data \n";
		DebugStop();
	}
	//If data is consistent, fAxis was computed during consistency check
	fCornerPoints = CornerPoints;
	ComputeArea();
}

// Copy constructor
DFNFracPlane::DFNFracPlane(const DFNFracPlane &copy){
    this->operator=(copy);
}



 DFNFracPlane &DFNFracPlane::operator=(const DFNFracPlane &copy)
 {
	fCornerPoints = copy.GetCornersX();
	fAxis = copy.fAxis;
	fArea = copy.area();
	fTolerance = copy.fTolerance;
	fPointsIndex = copy.fPointsIndex;
	return *this;
 }

/**
 * @brief Checks the consistency of the input data
 * @param Fracture plane coordinates (Matrix 3x4)
 * @return True if the four points are coplanar and the data is consistent
 * @return False if the points are not coplanar or the data is not consistent
 */
bool DFNFracPlane::Check_Data_Consistency(Matrix CornerPoints)
{
	// Checking vector consistency
	int cols = CornerPoints.Cols();
	int rows = CornerPoints.Rows();
	Matrix ax(3,3);
	
	if(rows != 3){    //Should be 3 (x y z)
		std::cout<<"Check the input data";
		DebugStop();
	}
	if(cols < 3 || cols > 10){//To form a plane it is needed at least 3 points but no more than 4
		std::cout<<"Check the input data (number of corner points, four is enough)";
		DebugStop();
	}

	//Ax0 without normalization
	ax(0,0) = CornerPoints(0,0) - CornerPoints(0,1);
	ax(1,0) = CornerPoints(1,0) - CornerPoints(1,1);
	ax(2,0) = CornerPoints(2,0) - CornerPoints(2,1);

	//Ax1 without normalization
	ax(0,1) = CornerPoints(0,2) - CornerPoints(0,1);
	ax(1,1) = CornerPoints(1,2) - CornerPoints(1,1);
	ax(2,1) = CornerPoints(2,2) - CornerPoints(2,1);

	//Ax2 computation
	ax(0,2) = ax(1,0)*ax(2,1) - ax(2,0)*ax(1,1);
	ax(1,2) = ax(2,0)*ax(0,1) - ax(0,0)*ax(2,1);
	ax(2,2) = ax(0,0)*ax(1,1) - ax(1,0)*ax(0,1);
	
	//Ax2 normalization
	REAL norm = sqrt(ax(0,2)*ax(0,2)+ax(1,2)*ax(1,2)+ax(2,2)*ax(2,2));
	for (int i = 0; i < 3; i++) // i< axis number (x y z)
	{ 
		ax(i, 2) = ax(i, 2)/norm;
	}

	//Coplanarity verification for quadrilateral plane
	if(cols == 4){
		//scalar product between Ax2 and <P3-P1> should be zero
		REAL ver = ax(0,2)*(CornerPoints(0,3)-CornerPoints(0,1))
				  +ax(1,2)*(CornerPoints(1,3)-CornerPoints(1,1))
				  +ax(2,2)*(CornerPoints(2,3)-CornerPoints(2,1));
		//Checks if points are coplanar
		if(std::abs(ver) > fTolerance){
			std::cout<<"Fracture corner points are not coplanar"<<"\n"<<std::endl;
			// DebugStop();
		}
	}

	// After checking the consistency the axis can be set
	fAxis = ax; // Set those in here to avoid re-computation in the constructor

	return true;
}

/**
 * @brief Get plane's corner points
 * @return Plane corner coordinates
 */
Matrix DFNFracPlane::GetCornersX() const{
    return fCornerPoints;
}

// /**
//  * @brief Computes area of plane
//  * @details Enumeration of corner points should follow standard PZ topology, where 
//  * corner nodes are numbered counter-clockwise (clockwise should work as well) from
//  * zero to N. (This condition will automatically be met for triangles, but not 
//  * always for quadrilaterals)
//  */
// void DFNFracPlane::ComputeArea(){
// 	int npoints = fCornerPoints.Cols();
// 	switch(npoints){
// 		case 3:{ //triangle
// 			//Area equals half the norm of the cross product between two edges
// 			REAL temp = pow(fAxis(1,0)*fAxis(2,1) - fAxis(2,0)*fAxis(1,1),2);
// 				temp += pow(fAxis(2,0)*fAxis(0,1) - fAxis(0,0)*fAxis(2,1),2);
// 				temp += pow(fAxis(0,0)*fAxis(1,1) - fAxis(1,0)*fAxis(0,1),2);
						
// 			fArea = sqrtl(temp)/2;
// 			break;
// 		}
// 		case 4:{ //quadrilateral
// 			//Area equals the sum of two triangles that fill the quadrilateral
// 			//Compute area of first triangle
// 			REAL temp1 = pow(fAxis(1,0)*fAxis(2,1) - fAxis(2,0)*fAxis(1,1),2);
// 				temp1 += pow(fAxis(2,0)*fAxis(0,1) - fAxis(0,0)*fAxis(2,1),2);
// 				temp1 += pow(fAxis(0,0)*fAxis(1,1) - fAxis(1,0)*fAxis(0,1),2);
				
// 				temp1 = sqrtl(temp1)/2;
			
// 			//Define vectors that correspond to the other two edges (as opposed to Ax0 and Ax1)
// 			TPZManVector<REAL, 3> ax3(3);
// 					ax3[0] = fCornerPoints(0,1) - fCornerPoints(0,2);
// 					ax3[1] = fCornerPoints(1,1) - fCornerPoints(1,2);
// 					ax3[2] = fCornerPoints(2,1) - fCornerPoints(2,2);
// 			TPZManVector<REAL, 3> ax4(3);
// 					ax4[0] = fCornerPoints(0,3) - fCornerPoints(0,2);
// 					ax4[1] = fCornerPoints(1,3) - fCornerPoints(1,2);
// 					ax4[2] = fCornerPoints(2,3) - fCornerPoints(2,2);
// 			//Compute area of second triangle
// 			REAL temp2 = pow(ax3[1]*ax4[2] - ax3[2]*ax4[1],2);
// 				temp2 += pow(ax3[2]*ax4[0] - ax3[0]*ax4[2],2);
// 				temp2 += pow(ax3[0]*ax4[1] - ax3[1]*ax4[0],2);
				
// 				temp2 = sqrtl(temp2)/2;

// 			fArea = temp1 + temp2;
// 			break;
// 		}
// 		default:{std::cout<<"Plane is not triangle nor quadrilateral? \n";

// 			DebugStop();
// 		}
// 	}
// }

/**
 * @brief Checks if a point is above or below the fracture plane
 * @param Point vector with the euclidean coordinates
 * @return True if the point is above the fracture plane
 * @return False if the point is below the fracture plane
 */

bool DFNFracPlane::Check_point_above(const TPZVec<REAL> &point) const{
    
    //Point distance to the fracture plane computation
        double point_distance = (point[0] - GetCornersX()(0,1))*(axis().GetVal(0,2)) 
                                +(point[1] - GetCornersX()(1,1))*(axis().GetVal(1,2)) 
                                +(point[2] - GetCornersX()(2,1))*(axis().GetVal(2,2));
        if (point_distance>0){
            return true;    //If the point is above de plane
        }
        else{
            return false;   //If the point is below de plane
        }
}

bool DFNFracPlane::Check_rib(TPZGeoEl *gel, TPZManVector<REAL,3> *intersection){
    if(gel->Dimension() != 1) {PZError<<"\n\n This ain't no rib \n\n"; DebugStop();}
    // Get rib's vertices
    TPZManVector<int64_t,2> inode(2,0);
    gel->GetNodeIndices(inode);
    TPZManVector<REAL,3> node0(3,0);
    TPZManVector<REAL,3> node1(3,0);
    gel->NodePtr(0)->GetCoordinates(node0);
    gel->NodePtr(1)->GetCoordinates(node1);
    return Check_rib(node0,node1,intersection);
}

bool DFNFracPlane::Check_rib(const TPZManVector<REAL,3> &p1, const TPZManVector<REAL,3> &p2, TPZManVector<REAL,3> *intersection) {
    //check for infinite plane
    if(Check_point_above(p1) != Check_point_above(p2)){
        //Rib cut by infinite plane
        //then calculate intersection point and check if it's within plane boundaries
        *intersection = CalculateIntersection(p1, p2);
        return IsPointInPlane(*intersection);
    }
    else
    {
        return false;    //Rib is not cut by plane
    }
}

/**
 * @brief Calculates the intersection point plane-rib
 * @param Point above the plane (vector)
 * @param Point below the plane (vector)
 * @return Intersecting point
 */
TPZManVector<double, 3> DFNFracPlane::CalculateIntersection(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2)
{     
    TPZVec<double> Pint;

    double term1 = ((fCornerPoints(0,1)-p2[0])*(fAxis(0,2))) 
                   +((fCornerPoints(1,1)-p2[1])*(fAxis(1,2)))
                   +((fCornerPoints(2,1)-p2[2])*(fAxis(2,2)));
    double term2 = (p1[0]-p2[0])*fAxis(0,2) 
                   +(p1[1]-p2[1])*fAxis(1,2) 
                   +(p1[2]-p2[2])*fAxis(2,2);
    double alpha = term1/term2;
    Pint.Resize(3);
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
 * @brief Checks if a coplanar point is within fracture plane
 * @details Enumeration of corner points should follow standard PZ topology, where 
 * corner nodes are numbered counter-clockwise. (This condition will automatically be 
 * met for triangles, but not always for quadrilaterals)
 * @param Point vector with the euclidean coordinates
 * @return True if the point is within fracture plane
 * @return False if the point is out of fracture plane
 */

bool DFNFracPlane::IsPointInPlane(TPZVec<REAL> &point) 
{
    int ncorners = fCornerPoints.Cols();
    REAL area = 0;
    for(int i = 0; i<ncorners; i++){
        //Define vectors from the point to a each one of a pair of corners
        TPZManVector<REAL, 3> ax1(3);
            ax1[0] = fCornerPoints(0,i) - point[0];
            ax1[1] = fCornerPoints(1,i) - point[1];
            ax1[2] = fCornerPoints(2,i) - point[2];
        TPZManVector<REAL, 3> ax2(3);
            ax2[0] = fCornerPoints(0,(i+1)%ncorners) - point[0];
            ax2[1] = fCornerPoints(1,(i+1)%ncorners) - point[1];
            ax2[2] = fCornerPoints(2,(i+1)%ncorners) - point[2];
        //Compute area of trangle outlined by these vectors
        REAL temp = pow(ax1[1]*ax2[2] - ax1[2]*ax2[1],2);
            temp += pow(ax1[2]*ax2[0] - ax1[0]*ax2[2],2);
            temp += pow(ax1[0]*ax2[1] - ax1[1]*ax2[0],2);
                  
        area += sqrtl(temp)/2;
    }
	
    // std::cout<<" ___ ";

    //If total computed area is equal to the plane's area, then
    //point is in plane
    return( fabs(area-fArea) < fTolerance );
}


TPZManVector<int64_t,4> DFNFracPlane::SetPointsInGeomesh(TPZGeoMesh *gmesh){
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
 * @brief Computes area of plane
 * @details Enumeration of corner points should follow standard PZ topology, where 
 * corner nodes are numbered counter-clockwise (clockwise should work as well) from
 * zero to N. (This condition will automatically be met for triangles, but not 
 * always for quadrilaterals)
 */
REAL DFNFracPlane::ComputeArea(){
	int n = fCornerPoints.Cols();
	TPZManVector<REAL,3> normal(3);
	for(int i = 0; i<3; i++){
		normal[i] = fAxis(i,2);
	}

    REAL area = 0;
    REAL an, ax, ay, az; // abs value of normal and its coords
    int  coord;           // coord to ignore: 1=x, 2=y, 3=z
    int  i, j, k;         // loop indices

    if (n < 3) return 0;  // a degenerate polygon

    // select largest abs coordinate to ignore for projection
    ax = fabs(normal[0]);    		// abs x-coord
    ay = fabs(normal[1]);    		// abs y-coord
    az = fabs(normal[2]);    		// abs z-coord

    coord = 3;                    // ignore z-coord
    if (ax > ay) {
        if (ax > az) coord = 1;   // ignore x-coord
    }
    else if (ay > az) coord = 2;  // ignore y-coord

    // compute area of the 2D projection
    switch (coord) {
      case 1:
        for (i=1, j=2, k=0; i<n; i++, j++, k++){
			j = j%n;
            area += (fCornerPoints(1,i) * (fCornerPoints(2,j) - fCornerPoints(2,k)));
		}
        break;
      case 2:
        for (i=1, j=2, k=0; i<n; i++, j++, k++){
			j = j%n;
            area += (fCornerPoints(2,i) * (fCornerPoints(0,j) - fCornerPoints(0,k)));
		}
		break;
      case 3:
        for (i=1, j=2, k=0; i<n; i++, j++, k++){
			j = j%n;
            area += (fCornerPoints(0,i) 
					*(fCornerPoints(1,j) 
					- fCornerPoints(1,k)));
		}
        break;
    }
    switch (coord) {    // wrap-around term
      case 1:
		area += (fCornerPoints(1,0) * (fCornerPoints(2,1) - fCornerPoints(2,n-1)));
        break;
      case 2:
		area += (fCornerPoints(2,0) * (fCornerPoints(0,1) - fCornerPoints(0,n-1)));
        break;
      case 3:
		area += (fCornerPoints(0,0) * (fCornerPoints(1,1) - fCornerPoints(1,n-1)));
        break;
    }

    // scale to get area before projection
    an = sqrt( ax*ax + ay*ay + az*az); // length of normal vector
    switch (coord) {
      case 1:
        area *= (an / (2 * normal[0]));
        break;
      case 2:
        area *= (an / (2 * normal[1]));
        break;
      case 3:
        area *= (an / (2 * normal[2]));
    }
	fArea = fabs(area);
    return fArea;
}



// /**
//  * @brief Creates a geometric element for this plane in pointed mesh
//  * @param Pointer to geometric mesh
//  * @return Index for newly created element in gmesh
//  */
// int64_t DFNFracPlane::CreateElement(TPZGeoMesh *gmesh){
// 	// number of nodes for gmesh
// 	int nnodes =  gmesh->NNodes();
// 	// nomber of corners for fracplane
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
// 	int64_t fracplaneindex = gmesh->NElements();
// 	gmesh->CreateGeoElement(elemtype, CornerIndexes, 40, fracplaneindex);

// 	return fracplaneindex;
// }



TPZManVector<REAL, 3> DFNFracPlane::GetProjectedX(TPZManVector<REAL, 3> &x){
    TPZManVector<REAL, 3> projection(3,0);
    REAL alpha = fAxis(0,2)*(x[0]-fCornerPoints(0,1))
                +fAxis(1,2)*(x[1]-fCornerPoints(1,1))
                +fAxis(2,2)*(x[2]-fCornerPoints(2,1));
    projection[0] = x[0] - alpha*fAxis(0,2);
    projection[1] = x[1] - alpha*fAxis(1,2);
    projection[2] = x[2] - alpha*fAxis(2,2);
    return projection;
}