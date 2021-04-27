#include "DFNNamespace.h"
#include "micropather.h"

using namespace micropather;

class DFNGraph : public micropather::Graph
{
    private:
        
        MPVector<int64_t> fpath;
        
        MicroPather* fpather;

        // Stores the distance between nodes
        TPZFMatrix<float> fdist;

        // Least cost
        float fleast = 1e-12;

        /// Mesh to graph node index map
        std::map<int64_t,int> fMeshToGraphNode;
        /// Graph to mesh node index map
        std::vector<int64_t> fGraphToMeshNode;

    public:
        /// Default constructor
        DFNGraph(TPZGeoMesh* gmesh, const std::set<int64_t>& edges){
            fpather = new MicroPather(this,20);
            ComputeCostMatrix(gmesh, edges);
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

        virtual float LeastCostEstimate( void* nodeStart, void* nodeEnd ){
            return 0.6666667f;
            // return (float)fleast;

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
            const float BIG_NUMBER = FLT_MAX;
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
        void ComputeShortestPath(int64_t start, int64_t end, TPZManVector<int64_t>& Path){
            float totalCost = 0.;
            micropather::MPVector<void*> solution;

            int graphstart = MeshToGraphNode(start);
            int graphend = MeshToGraphNode(end);

            int test = fpather->Solve(&graphstart,&graphend,&solution,&totalCost);
            int nedges = solution.size();
            Path.Resize(nedges,0);
            for(int i = 0; i<nedges; i++)
                {Path[i] = int64_t(solution[i]);}
        }

    private:
        void ComputeCostMatrix(TPZGeoMesh* gmesh, const std::set<int64_t>& edges){
            int estimative = edges.size()+2; // Number of nodes can never be bigger then twice the number of edges
            fdist.Resize(estimative,estimative);
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
            fdist.Print("cost1",std::cout,EFixedColumn);
            fdist.Resize(real_nnodes,real_nnodes);
            fdist.Print("cost2",std::cout,EFixedColumn);
        }
        /** 
         * Return the exact cost from the given node to all its neighboring nodes. This
         * may be called multiple times, or cached by the solver. It *must* return the same
         * exact values for every call to MicroPather::Solve(). It should generally be a simple,
         * fast function with no callbacks into the pather.
        */	
        virtual void AdjacentCost(void* node, MP_VECTOR< micropather::StateCost > *adjacent ){
            int nnodes = this->NNodes();
            int inode = *(int*)node;
            for(int jneig=0; jneig<nnodes; jneig++){
                if(IsZero(fdist(inode,jneig))) continue;
                micropather::StateCost travelcost = {&jneig,fdist(inode,jneig)};
                adjacent->push_back(travelcost);
            }
        }
        /**
            This function is only used in DEBUG mode - it dumps output to stdout. Since void* 
            aren't really human readable, normally you print out some concise info (like "(1,2)") 
            without an ending newline.
        */
        virtual void PrintStateInfo( void* node ){
            int x = 0;   
            int y = 1;   
            printf( "(%d,%d)", x, y );
        }
};