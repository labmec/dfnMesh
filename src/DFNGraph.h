#include "DFNNamespace.h"
#include "micropather.h"

using namespace micropather;

class DFNGraph : public micropather::Graph
{
    private:
        
        // MPVector<int64_t> fpath;
        
        MicroPather* fpather;

        // Stores the distance between nodes
        TPZFMatrix<float> fdist;

        // Least cost
        float fleast = 1e-12;

        /// Mesh to graph node index map
        std::map<int64_t,int> fMeshToGraphNode;
        /// Graph to mesh node index map
        std::vector<int64_t> fGraphToMeshNode;

        // A matrix to cache the connectivities by edge for efficiency
        TPZFMatrix<int64_t> fedges;

        /// Pointer to the DFN
        DFNMesh* fDFN = nullptr;

    public:
        /// Default constructor
        DFNGraph(DFNMesh* dfnmesh, const std::set<int64_t>& edges, int nnodes = 0){
            fpather = new MicroPather(this,20);
            fdist.Resize(nnodes,nnodes);
            fedges.Resize(nnodes,nnodes);
            fDFN = dfnmesh;
            
            ComputeCostMatrix(dfnmesh->Mesh(), edges);
            CacheEdgeConnectivity(dfnmesh->Mesh(), edges);
            
            fleast = ComputeMinimumLength();
        }
        /// Empty constructor
        DFNGraph()
           :fpather(0),
            fdist(0.,1,1)
        {}

        /// Destructor
        ~DFNGraph(){
            delete fpather;
        }

        virtual float LeastCostEstimate( int nodeStart, int nodeEnd ){
            return (float)fleast;
            // return 0.6666667f;

            // using namespace DFN;
            // int graphstart = *(int*)nodeStart;
            // int64_t meshstart = GraphToMeshNode(graphstart);
            // int graphend = *(int*)nodeEnd;
            // int64_t meshend = GraphToMeshNode(graphend);
            // TPZManVector<REAL,3> coordS(3,0.);
            // TPZManVector<REAL,3> coordE(3,0.);

            // TPZGeoMesh* fGMesh;
            // fGMesh->NodeVec()[meshstart].GetCoordinates(coordS);
            // fGMesh->NodeVec()[meshend].GetCoordinates(coordE);

            // TPZManVector<REAL,3> dif = coordS - coordE;
            // return Norm(dif);
        }


        int NNodes() const{return fdist.Rows();}

        /// Compute minimal distance in the weighted graph by taking the minimal non-zero entry in fdist
        float ComputeMinimumLength() const{
            int nnodes = NNodes();
            const float BIG_NUMBER = __FLT_MAX__;
            float minimal = BIG_NUMBER;
            for(int i=0; i<nnodes; i++){
                for(int j=i+1; j<nnodes; j++){
                    if(IsZero(fdist.g(i,j))) continue;
                    minimal = std::min(fdist.g(i,j),minimal);
                }
            }
            return minimal;
        }

        /// Convert mesh node index to graph node index. If not indexed yet, creates one
        int MeshToGraphNode(const int64_t index){
            int nnodes = fMeshToGraphNode.size();
            // insert only inserts if key is not found. Having the function like this, we save 2 binary searches per edge
            auto test = fMeshToGraphNode.insert(std::make_pair(index,nnodes));
            return test.first->second;
        }
        /// Convert graph node index to mesh node index. Simple vector element access
        int GraphToMeshNode(const int index) const{
            return fGraphToMeshNode[index];
        }

        /// @brief Solve shortest path between 2 nodes
        const bool ComputeShortestPath(int64_t start, int64_t end, TPZStack<int64_t>& Path){
            float totalCost = 0.;
            // solution gets the nodes of the shortest path
            micropather::MPVector<int> solution;

            int graphstart = MeshToGraphNode(start);
            int graphend = MeshToGraphNode(end);

            // Test gets the success result as an enum
            int test = fpather->Solve(graphstart,graphend,&solution,&totalCost);
            if(test == MicroPather::NO_SOLUTION){
                PZError << "\nCouldn't solve for shortest path while analysing graph for inter-fracture intersection.\n"
                        << "\nstart_localindex = " << graphstart
                        << "\nend_localindex   = " << graphend
                        << "\nstart_globalindex = " << start
                        << "\nend_globalindex = " << end;
                PZError << "\nConnectivity Matrix for the graph:\n";
                fedges.Print("GraphConnectivity", PZError, MatrixOutputFormat::EFormatted);
                return false;
//                fDFN->DFN_DebugStop();
            }
            
            // if(test == MicroPather::SOLVED)
            ConvertNodePathToEdgePath(solution,Path);
            return true;
        }

    private:
        void ComputeCostMatrix(TPZGeoMesh* gmesh, const std::set<int64_t>& edges){
            if(fdist.Rows() == 0){
                int estimative = edges.size()+2; // Number of nodes can never be bigger then twice the number of edges
                fdist.Resize(estimative,estimative);
            }
            fdist.Zero();
            for(auto index : edges){
                TPZGeoEl* edge = gmesh->Element(index);
                if(edge->Dimension() != 1) DebugStop();
                int i = MeshToGraphNode(edge->NodeIndex(0));
                int j = MeshToGraphNode(edge->NodeIndex(1));
                fdist(i,j) = edge->SideArea(2);
                fdist(j,i) = fdist(i,j); // Symmetry
            }
            int real_nnodes = fMeshToGraphNode.size();
            fdist.Resize(real_nnodes,real_nnodes);
            #if PZ_LOG
                // fdist.Print("cost",std::cout,EFixedColumn);
            #endif // PZ_LOG
        }
        /** 
         * Return the exact cost from the given node to all its neighboring nodes. This
         * may be called multiple times, or cached by the solver. It *must* return the same
         * exact values for every call to MicroPather::Solve(). It should generally be a simple,
         * fast function with no callbacks into the pather.
        */	
        virtual void AdjacentCost(int node, MP_VECTOR< micropather::StateCost > *adjacent ){
            int nnodes = this->NNodes();
            int inode = node;
            for(int jneig=0; jneig<nnodes; jneig++){
                if(IsZero(fdist(inode,jneig))) continue;
                micropather::StateCost travelcost = {jneig,fdist(inode,jneig)};
                adjacent->push_back(travelcost);
            }
        }
        /**
            This function is only used in DEBUG mode - normally you print out some concise info (like "(1,2)") 
            without an ending newline.
        */
        virtual void PrintStateInfo( int node ){
            DebugStop();
        }

        /** @brief Store edge connectivity but using graph node indices
         * this is only used in the end to convert a path in local nodes to edge indices
        */
        void CacheEdgeConnectivity(TPZGeoMesh* gmesh, const std::set<int64_t>& edges){
            if(fedges.Rows() == 0){
                int nnodes = NNodes();
                fedges.Resize(nnodes,nnodes);
            }
            fedges.Zero();
            for(int64_t index : edges){
                TPZGeoEl* edge = gmesh->Element(index);
                int i = MeshToGraphNode(edge->NodeIndex(0));
                int j = MeshToGraphNode(edge->NodeIndex(1));
                fedges(i,j) = index;
                fedges(j,i) = index; // Symmetry
            }
        }

        /** @brief converts a path given in graph nodes to a path listed as edges*/
        void ConvertNodePathToEdgePath(const micropather::MPVector<int>& NodePath, TPZStack<int64_t>& EdgePath) const{
            const int nnodes = NodePath.size();
            if(nnodes < 2) DebugStop();

            for(int inode=0; inode<nnodes-1; inode++){
                int i = NodePath[inode];
                int j = NodePath[inode+1];
                EdgePath.push_back(fedges.g(i,j));
            }
        }
};
