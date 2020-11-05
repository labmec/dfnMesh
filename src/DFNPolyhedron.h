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
    public:
		/// Empty constructor
		DFNPolyhedron(){};
		/// Empty destructor
		~DFNPolyhedron(){};

		// Initialize a polyhedron's datastructure
		template<int t_NAlloc>
		void Initialize(DFNMesh* dfn, int id, const TPZStack<std::pair<int64_t,int>,t_NAlloc>& oriented_shell)
		{
			fDFN = dfn;
			fIndex = id;
			fShell.clear();
			for(auto& orientedface : oriented_shell){
				fShell.insert(orientedface);
			}
		}
		
		/// Constructor from a container
		template<int t_NAlloc>
		DFNPolyhedron(
					DFNMesh* dfn, 
					int id, 
					const TPZStack<std::pair<int64_t,int>,t_NAlloc>& oriented_shell)
		{
			this->Initialize(dfn,id,oriented_shell);
		}
		
		/// Copy constructor
		DFNPolyhedron(const DFNPolyhedron &copy);
		
		/// Assignment operator
		DFNPolyhedron &operator=(const DFNPolyhedron &copy);

		/** @brief Print method for logging */
		void Print(std::ostream& out = std::cout) const;

		// Index of polyhedron
		int Index(){return fIndex;}

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
};


#endif /* DFNPolyhedron_h */

