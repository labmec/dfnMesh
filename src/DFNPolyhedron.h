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

class DFNMesh;
// class DFNFracture;

/*! 
 *  @brief      Describes a convex polyhedral volume as a set of oriented faces that enclose the shell around it.
 */
class DFNPolyhedron
{
    private:
		// A set of oriented faces that form the polyhedron. {index, orientation}
		std::map<int64_t, int> fShell;

		// Index given to this polyhedron by DFNMesh
		int fIndex = -1;

		// Pointer to the dfn this polyhedron belongs
		DFNMesh* fDFN = nullptr;

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
		void Initialize(DFNMesh* dfn, int id, const TPZVec<std::pair<int64_t,int>>& oriented_shell, const int coarseindex)
		{
			fDFN = dfn;
			fIndex = id;
			fShell.clear();
			for(auto& orientedface : oriented_shell){
				fShell.insert(orientedface);
			}
			fCoarseIndex = coarseindex;
		}
		
		/// Constructor from a container
		DFNPolyhedron(
					DFNMesh* dfn, 
					int id, 
					const TPZVec<std::pair<int64_t,int>>& oriented_shell, 
					const int coarseindex)
		{
			this->Initialize(dfn,id,oriented_shell,coarseindex);
		}
		
		/// Copy constructor
		DFNPolyhedron(const DFNPolyhedron &copy);
		
		/// Assignment operator
		DFNPolyhedron &operator=(const DFNPolyhedron &copy);

		/** @brief Print method for logging */
		void Print(std::ostream& out = std::cout, bool detailed = true);

		// Index of polyhedron
		int Index(){return fIndex;}

		// Index of coarse element to which this polyhedron is a subset
		int CoarseIndex(){return fCoarseIndex;}

		// Get number of faces around this polyhedron
		int NElements(){return fShell.size();}

		// Set a new index for this polyhedron
		void SwapIndex(const int newid);

		// Reference to set of oriented faces that form this polyhedron
		std::map<int64_t, int>& Shell(){return fShell;}

		// Lists DFNFaces that are part of this polyhedron
		void ListDFNFaces(DFNFracture* fracture, TPZStack<DFNFace*> facelist);

		/** @brief Remove faces from this polyhedron*/
		void RemoveFaces(const TPZVec<std::pair<int64_t,int>>& facestack);

		/** @brief Checks if this polyhedron was split in 2 smaller polyhedra*/
		bool IsRefined();

		/** @brief Remove father from shell and add its subelements */
		void SwapForChildren(TPZGeoEl* father);

		/** @brief Fills a vector with the indices of the 1D skeleton elements around the faces of this polyhedron */
		void GetEdges(TPZVec<TPZGeoEl*>& edgelist);

		/** @brief Checks if this polyhedron was intersected by a fracture
		 * @todo(maybe) - This method could be made more restrictive if we wished to do so. I'll leave it returning true if any face of fShell is intersected, but we could potentially want it to return true only if the intersection would actually lead to the polyhedron being split. This, of course, would cost some extra computations, and I don't really see a need for this more restrictive version yet.
		*/
		bool IsIntersected(DFNFracture& fracture);

		/** @brief Returns true if any of the faces in this polyhedron's shell contains only one Inbound rib*/
		bool IntersectsFracLimit(DFNFracture& fracture);
};


#endif /* DFNPolyhedron_h */

