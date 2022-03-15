/*! 
 *  @brief     Describes a planar polygon from it's corner points.
 *  @details   
 *  @author    Pedro Lima
 *  @date      2019.05
 */

#ifndef DFNPolygon_h
#define DFNPolygon_h

#include "pzfmatrix.h"
#include "pzstack.h"


typedef TPZFMatrix<REAL> Matrix;
/// A geometrical intersection between 2 bounded planes in R^3 is a line segment, 
/// so we can represent it by the coordinates of its nodes
typedef TPZStack<TPZManVector<REAL,3>,2> Segment;

/*! 
 *  @brief     Describes a planar convex polygon from it's corner points. (not to be confused with DFNPolyhedron)
 *  @details Enumeration of corner points should follow standard PZ topology, where 
 *  corner nodes are numbered counter-clockwise (clockwise should work as well) from
 *  zero to N. (This condition will automatically be met for triangles, but not 
 *  always for quadrilaterals and higher order polygons)
 */
class DFNPolygon
{
  private:
	/// Contains polygon corner coordinates. Matrix 3xn (n is the number of corners)
	Matrix fCornerPoints;
	
	/// Axis that define polygon orientation (Ax0 from node1 to node0, Ax1 from node1 to node2 and Ax2 the normal vector). Matrix 3x3
	Matrix fAxis ;

	/// Area of polygon
	double fArea = -1.0;

	/// If nodes of this polygon have been added to a geometric mesh, this vector holds GeoNodes indices
	TPZManVector<int64_t> fPointsIndex;

    /// Tracks which nodes are above this Polygon - for each node of the geometric mesh, indicates whether the node is above the plane
    TPZVec<bool> fNodesAbove;
	
	/// Index of biggest absolute component of normal vector (for costless projection)
	int fMax_component = -1;
	
  public:
	/// Empty constructor
	DFNPolygon(){};
	/// Empty destructor
	~DFNPolygon(){};

	/// Default constructor. Matrix should be 3xN (3 coordinates for each of the N corner points)
    // the first 3 points cannot be colinear
	DFNPolygon(const Matrix &CornerPoints, const TPZGeoMesh* gmesh);
    
    DFNPolygon(const Matrix &CornerPoints);

	/// Copy constructor
	DFNPolygon(const DFNPolygon &copy);

	/// Assignment operator
	DFNPolygon &operator=(const DFNPolygon &copy);
    
	/// Return number of corners of polygon
	int NCornerNodes() const{return fCornerPoints.Cols();}
	int NEdges() const{return fCornerPoints.Cols();}
    
    /// compute the direction of the axes based on the first three nodes
    void ComputeAxis();

	/// Define corner coordinates
	void SetCornersX(const Matrix &CornerPoints);

	/// axis(i, j) returns component i of axis j
	REAL axis(int row, int col){return fAxis(row,col);}

	/// Get matrix with axis 0, 1 and 2 on each column
	Matrix axis() const {return fAxis;}

	/// Return corner coordinates
	const Matrix& GetCornersX() const;

	/// Fill corner coordinates
	void iCornerX(int icorner, TPZManVector<REAL,3>& coord) const{
		coord[0] = fCornerPoints.g(0,icorner);
		coord[1] = fCornerPoints.g(1,icorner);
		coord[2] = fCornerPoints.g(2,icorner);
	}
    
    const REAL PointAndCoor(int icorner, int coor) const{
        return fCornerPoints.g(coor,icorner);
    }

	/// Return area of polygon
	double area() const { return fArea; }

	/// Compute area of polygon
	REAL ComputeArea();
 
    /// Sort nodes as above or below this polygon
    void SortNodes(const TPZGeoMesh* gmesh);

	/// Same as DFNPolygon::SortNodes but only for possible new nodes in gmesh. Starts from the index in the gmesh::NodeVector equal to the current size of this polygon fNodesAbove vector.
	/// If you call this method with an uninitialized DFNPolygon::fNodesAbove, it'll just call on DFPolygon::SortNodes
	// void SortNodes_update(const TPZGeoMesh* gmesh);

	/// Check if the segment that conects 2 coordinates has intersection with this polygon
	bool Check_pair(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2, TPZManVector<REAL,3> &intersection) const;

	/** @brief Return true if the rib intersects the polygon by also checking if intersection point is within bounds of the polygon
	 * @attention If you want to skip check on polygon limits, use DFNPolygon::IsCutByPlane()
	 * @param intersection [out] If intersected, fills vector with intersection coordinates
	*/
	bool IsCutByPolygon(TPZGeoEl *rib, TPZManVector<REAL,3> &intersection) const;
	/** @brief Return true if the rib intersects the plane that contains the polygon
	 * @attention This doesn't check if intersection is in the bounds of the polygon. If you want to check on polygon limits, use DFNPolygon::IsCutByPolygon()
	 * @param intersection [out] If intersected, fills vector with intersection coordinates
	*/
	bool IsCutByPlane(TPZGeoEl *rib, TPZManVector<REAL,3> &intersection) const;

	/// Return true if a point is above the polygon plane
   	bool Compute_PointAbove(const TPZVec<REAL> &point) const;

   	/// Check whether the point coordinates are within the polygon region
   	bool IsPointInPolygon(const TPZVec<REAL> &point) const;
	
	/// Computes the intersection point with the polygon
	TPZManVector<REAL,3> CalculateIntersection(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2) const;
   

	/**
	 * @brief Inserts corner nodes in geometric mesh
	 * @param Pointer to geometric mesh
	 * @return Vector with nodes indices in gmesh
	 */
	const TPZVec<int64_t> &SetPointsInGeomesh(TPZGeoMesh *gmesh);

	/**
	 * @brief Returns index of GeoNode that was created for corner i using DFNPolygon::SetPointsInGeomesh(gmesh, matID)
	 * @param Local index of corner
	 */
	int64_t CornerIndex(int i) const{
		DebugStop(); // Deprecated
		return fPointsIndex[i];
	}

	/**
	 * @brief Takes a set of coordinates in 3D and returns its projection onto the polygon
	*/
	TPZManVector<REAL, 3> GetProjectedX(const TPZManVector<REAL, 3> &point) const;

	/// Fill a 3D vector with the components of the normal direction of this polygon
	void GetNormal(TPZManVector<REAL,3>& normal_vec) const{
		normal_vec.resize(3);
		normal_vec[0] = fAxis(0,2);
		normal_vec[1] = fAxis(1,2);
		normal_vec[2] = fAxis(2,2);
	}

	/**
	 * @brief Inserts a geometric element + nodes for this polygon.
	 * @note 1. It only works for triangles and quadrilaterals
	 * @note 2. I'm leaving this here for graphical debugging. Not for the actual flow of the code
	 * @param nodes: Optional pointer to vector of existing nodes
	*/
	TPZGeoEl* InsertGeoEl(TPZGeoMesh* gmesh, int matid = 100, TPZVec<int64_t>* nodes = nullptr);

	/**
	 * @brief Inserts a few geometric elements to graphically represent this polygon. 
	 * @param gmesh - Pointer to geometric mesh where to insert elements; 
	 * @param matid - Material id to use for graphical element(s); 
	 * @param exclusive_dim - Exclusive dimension of graphical elements: 
	 * If given 2, inserts 2D elements to cover surface of polygon; 
	 * If given 1, inserts only 1D elements for the polygon boundary; 
	 * If given 0, do both 1D and 2D.
	*/
	TPZVec<TPZGeoEl*> InsertGeomRepresentation(TPZGeoMesh* gmesh, int matid = 100, int exclusive_dim = 0);

	void ComputeCentroid(TPZVec<REAL>& centroid);

	/// Plot polygon graphic to VTK file. You can set a material id and an element index to it and they'll be exported to VTK file
	void PlotVTK(const std::string filepath, const int materialID = 0, const int64_t index = 0) const;

	/** @brief Print method for logging */
    void Print(std::ostream& out = std::cout) const;

	void PlotNodesAbove_n_Below(TPZGeoMesh* gmesh);

	/** @brief Check if node is above plane by checking fNodesAbove */
	bool IsPointAbove(int64_t index) const{return fNodesAbove[index];}

	/// Compute the length of an edge
	REAL EdgeLength(int edgeindex) const;
	/// Compute a vector from the first to last nodes of an edge
	void GetEdgeVector(int edgeindex, TPZVec<REAL>& edgevector) const;

	/// Compute an intersection segment between 2 DFNPolygons (if there is any)
	bool ComputePolygonIntersection(const DFNPolygon& otherpolyg, Segment& intersection_segment) const;

    /// For a subpolygon, returns the cosine of the "worst angle". The bigger is the angles, the less planar is the subpolygon
    const REAL GetWorstAngleCos();
    
  private:
	/// Checks consistency and initializes the datastructure of the object
	bool Check_Data_Consistency() const;
    
    /// Calculates the area of a subtriangle given by a centroid and two consecutives indexes (i and i+1)
    const REAL SubTriangleArea(const TPZVec<REAL>& centroid, const int i) const;
    
    /// Returns the cross products of of the nodes in index i, j, and k
    TPZManVector<REAL,3> GetCrossProduct(const int i, const int j, const int k);

};

inline std::ostream& operator<<(std::ostream &out, const DFNPolygon& polygon){
    polygon.Print(out);
    return out;
}

#endif /* DFNPolygon */
