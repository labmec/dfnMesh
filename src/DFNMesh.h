/*! 
 *  @authors   Pedro Lima
 *  @date      2020.06
 */

#ifndef DFNMesh_h
#define DFNMesh_h

#include "pzgmesh.h"
#include "DFNFracture.h"
#include "DFNVolume.h"
#include <gmsh.h>


/**
 * @brief Describes a fractured MHM geomesh
*/
class DFNMesh{
	private:
        // A set of fractures
		std::list<DFNFracture *> fFractures;
		
        // A set of volumes
        std::map<int64_t, DFNVolume> fVolumes;
		
        // Minimum acceptable distance/length for point insertion
        REAL fToleranceLength = 1e-5;

        // Minimum acceptable angle/aspect-ratio
		REAL fToleranceRatio = 0.2;

        // Pointer to geometric mesh
		TPZGeoMesh *fGMesh;

	public:
        // // Local enumeration of material indices
        // enum class Ematerial {intact = 1, fracture = 2, surface = 2, refined = 3, transition = 3};

        /// Material tags
        // const int MatIntact = 1;
        // const int MatFracture = 2;
        // const int MatRefined = 3;

		/// Constructor
		DFNMesh(TPZGeoMesh *gmesh);

        /// Add new fracture
        void AddFracture(DFNFracture *fracture){fFractures.push_back(fracture);}

		/// Pointer to volume of index 'index'
		DFNVolume *Volume(int64_t index){return &fVolumes[index];}

        /// Pointer to geometric mesh
        TPZGeoMesh *Mesh(){return fGMesh;}

		/** 
         * @brief Insert intersection elements of lower dimension in the geometric mesh.
         * @note matid = -1 will use the material id of first that is found. Set it if you want to force
         */
    	void CreateSkeletonElements(int dimension, int matid = -1);
		
		/// Setup datastructure for fractured volumes (including finding fracture elements enclosed by them)
		void CreateVolumes();
		
		/// Exports a .geo file for this mesh
    	void ExportGMshCAD(std::string filename);
		
		/// Uses gmsh API to tetrahedralize a DFNVolume
		void Tetrahedralize(DFNVolume *volume);
    	
		/// Find the volumetrical element that encloses a 2D element
    	bool FindEnclosingVolume(TPZGeoEl *ifracface);
		/**
		 * @brief Tetrahedralize fractured volumes and refine intact volumes
		 * @todo refine intact volumes and maybe we could pass a target measure for subelements size as parameter
		*/
		void GenerateSubMesh();
        /**
         * @brief Prints DFN Geometric Mesh. 
         * @param pzmesh : File name for geomesh txt. Feed "skip" to skip
         * @param vtkmesh : File name for geomesh vtk. Feed "skip" to skip
         * @param MaterialIDs...
         * @todo this method is unfinished
        */
        void Print(std::string pzmesh = "pzmesh.txt"
                  ,std::string vtkmesh = "vtkmesh.vtk"
                  ,int fracture = 3
                  ,int transition = 2
                  ,int intact = 1);
        /**
         * @brief Prints DFN Geometric Mesh and material ids are renumbered for VTK colorful print of refinement of 2D elements :) 
         * @param pzmesh : File name for geomesh txt. Feed "skip" to skip
         * @param vtkmesh : File name for geomesh vtk. Feed "skip" to skip
        */
        void PrintColorful(std::string pzmesh = "pzmesh.txt"
							,std::string vtkmesh = "vtkmesh.vtk");
		
        /**
         * @brief Imports d-dimensional elements from a GMsh::model to a TPZGeoMesh. Imported 
         * elements are pushed to the back of TPZGeoMesh::ElementVector.
         * @note (1) Must be called between the pair gmsh::initialize and gmsh::finalize of the
         * model from which new elements should be read.
         * @note (2) If GMsh has created any new nodes, those will be inserted into TPZGeoMesh aswell
         * @param gmesh: Pointer to TPZGeoMesh
         * @param dimension: Dimension of elements that should be imported
         * @param oldnodes: A set of indices of the old nodes that were used to define the geometry in 
         * GMsh (so that new nodes may be identified)
         */
        void ImportElementsFromGMSH(TPZGeoMesh * gmesh, int dimension, std::set<int64_t> &oldnodes);

        int Dimension(){return fGMesh->Dimension();}

	private:
		// /**
    	//  *  @brief Navigate children tree to access most extreme branches
    	//  *  @param gel: Pointer to geometric element of ancestor
    	//  *  @param outfile: ofstream in which to write accessed data
    	//  */
    	// void PrintYoungestChildren(TPZGeoEl *gel, std::ofstream &outfile);
		
		/**
		 * @brief Navigate through neighbours of first level, than second level, and so on, until an element of a specific material id is found
		 * @returns Index of eldest ancestor of such element
		*/
		int64_t SearchIndirectNeighbours(TPZGeoEl* gel);
		/**
		 * @brief Goes through all neighbours of gel and identifies if any of them has material id different of that from surface elements
		 * @returns Index of eldest ancestor of first found neighbour that matches the material id criteria (macro element)
		 * @note Later I might modify this to fill a vector/list with all neighbours that match the material id criteria
		*/
		int64_t FindAdjacentMacroEl(TPZGeoEl* gel);
		/**
		 * @brief Pushes all neighbours of a geometric element onto the back of a list
		 * @note Not all neighbours are pushed to the list, but rather some criteria are specified so that it gets only those that are candidates of higher interest. Currently these would be 2D elements that are neighbours through the edges of gel.
		 * @param gel: Pointer to geometric element
		 * @param candidate_queue: Reference to current list of candidates 
		*/
		void QueueNeighbours(TPZGeoEl* gel,   std::list<int64_t> &candidate_queue);
        
        /**
         * @brief Deletes face + ribs (if left isolated) + nodes (if left isolated)
         * @param face: 2D element to be deleted
        */
        void DeleteElementAndRibs(TPZGeoEl *face);
        /**
         * @brief Deletes gel + children + isolated ribs + unused nodes
         * @note It will assume element has been found not to belong to the domain of interest and will not verify
         * @param gel: pointer to the geometric element
        */
        void CropExternalElement(TPZGeoEl *gel);
        /**
         * @brief Deletes gel and all elements that share the same eldest ancestor, then deletes the ancestor
         * @param gel: Any member of the family
        */
        void DeleteFamily(TPZGeoEl *gel);
        /**
         * @brief Check if the neighbour has equal dimension
         * @param geliside GeoElement side
         */

        bool HasEqualDimensionNeighbour(TPZGeoElSide &gelside);

        /// Set all material ids to 1
        void ClearMaterials();

        // @todo
        void RestoreMaterials(TPZGeoMesh *originalmesh);
};



enum DFNMaterial{
    Eintact = 1, 
    Efracture = 2, 
    // Esurface = 2, 
    Erefined = 3, 
    // Etransition = 3
};

// A small number for geometric tolerances
static const double gDFN_SmallNumber = 1.e-3;

#endif /* DFNMesh_h */