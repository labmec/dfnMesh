/*! 
 *  @brief     Describes a rectangular plane from four corner points.
 *  @details   
 *  @author    Pedro Lima
 *  @date      2019
 */


#include "TRSFracPlane.h"
#include <math.h>

//Constructor
TRSFracPlane::TRSFracPlane(const Matrix &CornerPoints)
{
    if( !Check_Data_Consistency(CornerPoints) )	{
		std::cout<<"Error at TRSFracPlane: Bad input data \n";
		DebugStop();
	}
	//If data is consistent, fAxis was computed during consistency check
	fCornerPoints = CornerPoints;
    // Center point computation 
	fCenterCo.Resize(3);
	for(int i=0; i< 4; i++){          /// i< axis number (x y z)
		fCenterCo[0] += (1./4)*fCornerPoints(0,i); ///Center point X axis
		fCenterCo[1] += (1./4)*fCornerPoints(1,i); ///Center point Y axis
		fCenterCo[2] += (1./4)*fCornerPoints(2,i); ///Center point Z axis
	}

    //L0 and L1 computation
	fL0 = fabs(
				(fCornerPoints(0,0)-fCornerPoints(0,1))*fAxis(0,0)
				+(fCornerPoints(1,0)-fCornerPoints(1,1))*fAxis(1,0)
				+(fCornerPoints(2,0)-fCornerPoints(2,1))*fAxis(2,0)
			);
	fL1 = fabs(
				(fCornerPoints(0,0)-fCornerPoints(0,3))*fAxis(0,1)
				+(fCornerPoints(1,0)-fCornerPoints(1,3))*fAxis(1,1)
				+(fCornerPoints(2,0)-fCornerPoints(2,3))*fAxis(2,1)
			);
}

// Copy constructor
TRSFracPlane::TRSFracPlane(const TRSFracPlane &copy){
	fCornerPoints = copy.GetCorners();
	fAxis = copy.fAxis;
	fCenterCo = copy.fCenterCo;
	fL0 = copy.fL0;
	fL1 = copy.fL1;
}



 TRSFracPlane &TRSFracPlane::operator=(const TRSFracPlane &copy)
 {
	fCornerPoints = copy.GetCorners();
	fAxis = copy.fAxis;
	fCenterCo = copy.fCenterCo;
	fL0 = copy.fL0;
	fL1 = copy.fL1;
	return *this;
 }

/**
 * @brief Checks the consistency of the input data
 * @param Fracture plane coordinates (Matrix 3x4)
 * @return True if the four points are coplanar and the data is consistent
 * @return False if the points are not coplanar or the data is not consistent
 */
bool TRSFracPlane::Check_Data_Consistency(Matrix CornerPoints)
{
	// Checking vector consistency
	int cols = CornerPoints.Cols();
	int rows = CornerPoints.Rows();
	Matrix ax(3,3);
	
	if(rows != 3){    //Should be 3 (x y z)
		std::cout<<"Check the input data";
		DebugStop();
	}
	if(cols < 3 or cols > 4){//To form a plane it is needed at least 3 points but no more than 4
		std::cout<<"Check the input data (number of corner points, four is enough)";
		DebugStop();
	}
		if (!(cols ==4 && rows ==3)){
			std::cout<<"Check the input data (number of plane points, four is enough)";
			DebugStop();
		}

    // Mid points computation for the first three points
		Matrix MidPoints;
		MidPoints.Resize(3,cols);              //3 rows (x y z), mid points
		for(int i=0; i<(cols-1); i++){
			MidPoints(0,i)=(CornerPoints(0,i)+CornerPoints(0,i+1))/2;
			MidPoints(1,i)=(CornerPoints(1,i)+CornerPoints(1,i+1))/2;
			MidPoints(2,i)=(CornerPoints(2,i)+CornerPoints(2,i+1))/2;
		}
    // Mid points computation for the last point
		MidPoints(0,cols-1)=(CornerPoints(0,(cols)-1)+CornerPoints(0,0))/2;
		MidPoints(1,cols-1)=(CornerPoints(1,(cols)-1)+CornerPoints(1,0))/2;
		MidPoints(2,cols-1)=(CornerPoints(2,(cols)-1)+CornerPoints(2,0))/2;

    //Ax0 without normalization
		ax(0,0)=MidPoints(0,cols-1)-MidPoints(0,1);
		ax(1,0)=MidPoints(1,cols-1)-MidPoints(1,1);
		ax(2,0)=MidPoints(2,cols-1)-MidPoints(2,1);

    //Ax1 without normalization
		ax(0,1)=MidPoints(0,cols-2)-MidPoints(0,0);
		ax(1,1)=MidPoints(1,cols-2)-MidPoints(1,0);
		ax(2,1)=MidPoints(2,cols-2)-MidPoints(2,0);

    //Ax0 and Ax1 normalization
		double norm = sqrt(ax(0,0)*ax(0,0)+ax(1,0)*ax(1,0)+ax(2,0)*ax(2,0));
		for (int i = 0; i < 3; i++) // i< axis number (x y z)
		{ 
			ax(i, 0) = ax(i, 0)/norm;
		}
		norm = sqrt(ax(0,1)*ax(0,1)+ax(1,1)*ax(1,1)+ax(2,1)*ax(2,1));
		for (int i = 0; i < 3; i++) // i< axis number (x y z)
		{ 
			ax(i, 1) = ax(i, 1)/norm;
		}

    //Ax2 computation
        ax(0,2)=ax(1,0)*ax(2,1) - ax(2,0)*ax(1,1);
        ax(1,2)=ax(2,0)*ax(0,1) - ax(0,0)*ax(2,1);
        ax(2,2)=ax(0,0)*ax(1,1) - ax(1,0)*ax(0,1);

    //Coplanar verification
		//scalar product between Ax2 and <P3-P1> should be zero
		double ver = ax(0,2)*(CornerPoints(0,3)-CornerPoints(0,1))
					+ax(1,2)*(CornerPoints(1,3)-CornerPoints(1,1))
					+ax(2,2)*(CornerPoints(2,3)-CornerPoints(2,1));
		//Checks if the points are coplanar
		if(abs(ver) > fTolerance){
			std::cout<<"The input points are not coplanar"<<"\n"<<std::endl;
			DebugStop();
		}
		
    // After checking the consistency the axis can be set
        fAxis = ax; // Set those in here to avoid re-computation in the constructor
	
	// fMidPoints.Resize(3,cols);
	// fMidPoints = MidPoints;
	return true;
}

/**
 * @brief Get plane's corner points
 * @return Plane corner coordinates
 */
Matrix TRSFracPlane::GetCorners() const{
    return fCornerPoints;
}

