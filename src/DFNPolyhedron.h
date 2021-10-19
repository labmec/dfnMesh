/*! 
 *  @brief     Describes a convex polyhedral volume as a set of oriented faces that enclose the shell around it.
 *  @details   
 *  @authors   Pedro Lima, Philippe Devloo
 *  @date      2020
 */

#ifndef DFNPolyhedron_h
#define DFNPolyhedron_h

#include "pzgmesh.h"
// #include "DFNMesh.h"
// #include "DFNFracture.h"
#include "DFNFace.h"

// class DFNFace;
class DFNMesh;
// class DFNFracture;

/*! 
 *  @brief Describes a convex polyhedral volume as a set of oriented faces that enclose the shell around it. (not to be confused with DFNPolygon)
 */
// @pedro - what is the difference/relation between a DFNPolyhedron and a DFNVolume?
// looked at DFNVolume that is deprecated...
class DFNPolyhedron
{
    private:
		// Pointer to the dfn this polyhedron belongs
		DFNMesh* fDFN = nullptr;

		// A set of oriented faces that form the polyhedron. {index, orientation}
		std::map<int64_t, int> fShell;

		// Index given to this polyhedron by DFNMesh
		int fIndex = -1;

		// Index of the coarse element where this polyhedron is contained
		int fCoarseIndex = -1;

		/** A flag for convexity. Polyhedra can only be non-convex if they have been refined into convex subsets.
		 * @note We weren't going to need this, but I noticed creating instances of this class for non-convex is the more elegant way to keep a continuous inheritance of fCoarseIndex as these polyhedra get refined. That's the only reason I'm keeping it.
		 */
		int fIsConvex = true;
		
    public:
		/// Empty constructor
		DFNPolyhedron(){};
		/// Empty destructor
		~DFNPolyhedron(){};

		// Initialize a polyhedron's datastructure
		void Initialize(DFNMesh* dfn, int id, const TPZVec<std::pair<int64_t,int>>& oriented_shell, const int coarseindex, const bool isConvex = true)
		{
			fDFN = dfn;
			fIndex = id;
			fShell.clear();
			for(auto& orientedface : oriented_shell){
				fShell.insert(orientedface);
			}
			fCoarseIndex = coarseindex;
			fIsConvex = isConvex;
		}
		
		/// Constructor from a container
		DFNPolyhedron(
					DFNMesh* dfn, 
					int id, 
					const TPZVec<std::pair<int64_t,int>>& oriented_shell, 
					const int coarseindex,
					const bool isConvex = true)
		{
			this->Initialize(dfn, id, oriented_shell, coarseindex, isConvex);
		}
		
		/// Copy constructor
		DFNPolyhedron(const DFNPolyhedron &copy);
		
		/// Assignment operator
		DFNPolyhedron &operator=(const DFNPolyhedron &copy);

		/** @brief Print method for logging */
		void Print(std::ostream& out = std::cout, bool detailed = true) const;

		// Index of polyhedron
		int Index() const{return fIndex;}

		// Index of coarse element to which this polyhedron is a subset
		int CoarseIndex() const{return fCoarseIndex;}

		// Get number of faces around this polyhedron
		int NElements() const{return fShell.size();}

		// Check flag of this polyhedron for convexity
		bool IsConvex() const{return fIsConvex;}

		// Set a new index for this polyhedron
		void SwapIndex(const int newid);

		// Reference to set of oriented faces that form this polyhedron
		std::map<int64_t, int>& Shell(){return fShell;}

		// Lists DFNFaces that are part of this polyhedron
		void ListDFNFaces(DFNFracture* fracture, TPZStack<DFNFace*> facelist);

		/** @brief Remove faces from this polyhedron*/
		void RemoveFaces(const TPZVec<std::pair<int64_t,int>>& facestack);

		/** @brief Checks if this polyhedron was split into smaller polyhedra*/
		bool IsRefined()const;

		/** @brief Remove father from shell and add its subelements */
		void SwapForChildren(TPZGeoEl* father);

		/** @brief Get a set with the indices of the 1D skeleton elements around the faces of this polyhedron */
		std::set<int64_t> GetEdges() const;
		/** @brief Get the indices of edges around this volume that are elements of a set. This returns the intersection between DFNPolyhedron::GetEdges and the argument 'SuperSet' */
		std::set<int64_t> GetEdges_InSet(const std::set<int64_t>& SuperSet) const;

		/** @brief Checks if this polyhedron was intersected by a fracture
		 * @todo(maybe) - This method could be made more restrictive if we wished to do so. I'll leave it returning true if any face of fShell is intersected, but we could potentially want it to return true only if the intersection would actually lead to the polyhedron being split. This, of course, would cost some extra computations, and I don't really see a need for this more restrictive version yet.
		*/
		bool IsIntersected(DFNFracture& fracture)const;

		/** @brief Returns true if any of the faces in this polyhedron's shell contains only one Inbound rib*/
		bool IntersectsFracLimit(DFNFracture& fracture)const;

		/** @brief Call on gmsh to refine itself 
		 * @details Any quadrilateral in the shell will get refined to at least 2 triangles
		*/
		void Refine();
    
        /**
         * Prints the shells of the polyhedron to vtk
         */
        void PlotVTK(const std::string filepath = "") const;

		/// Plots this volume and its neighbour volumes
		void PlotVTK_NeighbourVolumes(const std::string filepath = "") const;

		/// Checks if this Polyhedron is a tetrahedron
		bool IsTetrahedron() const;

		/** @brief Check if the 2D element of a given index is on the shell of this volume. Checks both possible orientations.
		 * @param faceindex Index of 2D geoelement
		*/
		bool IsBoundedBy(const int64_t faceindex) const;
		/** @brief Check if the oriented 2D element of a given index is on the shell of this volume
		 * @param faceindex Index of 2D geoelement
		 * @param orientation Element orientation (positive means element normal vector points toward the interior of the volume)
		*/
		bool IsBoundedBy(const int64_t faceindex, const int orientation) const;
		/** @brief Check if the oriented 2D element of a given index is on the shell of this volume
		 * @details faceindex: Index of 2D geoelement
		 * @details orientation: Element orientation (positive means element normal vector points toward the interior of the volume)
		 * @param orientedFace: Paired {elIndex, orientation}.
		*/
		bool IsBoundedBy(const std::pair<int64_t, int> orientedFace) const;
		/** @brief Same as DFNPolyhedron::IsBoundedBy(), but uses the father of the element whose index is input.
		 * @param childindex Index of 2D geoelement whose father is to be tested
		*/
		bool IsBoundedByFather(const int64_t childindex) const;

		/** @brief Given a closed loop of edges on this volume's shell, get a set of faces which are exclusively to a side of this loop
		 * @param delimiter a set of edge-sides of faces {2D_gel | 1D_side}. 1D sides should form a closed loop to partition the shell, and the 2D els should all lie on the same subset of the shell you want me to assemble.
		 * @warning I will not verify if your delimiter is coherent. If it's not, the returned set will be the complete shell instead of a subset
		*/
		std::set<int64_t> GetShellSubset(const std::set<TPZGeoElSide>& delimiter) const;

		/** @brief Computes total surface area of the elements in the shell of this volume*/
		REAL ComputeArea() const;

#ifdef PZDEBUG
		void CoherentRefinementTest(const TPZVec<int64_t>& newgels);
#endif // PZDEBUG
};

inline std::ostream& operator<<(std::ostream &out, const DFNPolyhedron& polyhedron){
    polyhedron.Print(out);
    return out;
}

#endif /* DFNPolyhedron_h */

