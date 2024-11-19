#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "cnagraph.h"
#include "cnatree.h"


namespace py = pybind11;

typedef std::tuple<int, int, int, int> genotype;
typedef std::map<int, genotype> snvGenotype;
typedef std::map<int, snvGenotype> genotypes;


//wrapper to convert CnaTree:CnaEdgeSetSet to python friendly format (python doesn't seem to like sets)
const std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> genCNATrees(const IntPairSet& L,
                                                    int root_x, int root_y){

    CnaTree::CnaEdgeSetSet scriptS = CnaGraph::getCnaTrees(L, root_x, root_y);
    std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> scriptT;

    for (auto tree : scriptS) {
        std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> tr;

        for (auto edge : tree) {
            // Extract parent and child values
            int parent_x = edge.first._x;
            int parent_y = edge.first._y;
            int child_x = edge.second._x;
            int child_y = edge.second._y;

     
            // // Create pairs of integers for parent and child
            std::pair<int, int> par(parent_x, parent_y);
            std::pair<int, int> chi(child_x, child_y);

            // Create a pair of pairs for the edge
            std::pair<std::pair<int, int>, std::pair<int, int>> edge_pair(par, chi);

            // Insert the pair of pairs into the vector
            tr.push_back(edge_pair);
        }
  
        // Insert the vector of pairs of pairs into the outer vector
        scriptT.push_back(tr);
    }


    return scriptT;

 }
//  

 std::vector<std::vector<std::pair<std::tuple<int, int,int, int>, std::tuple<int, int, int, int>>>> genGenoTrees(std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> cnatree){
        
        CnaTree::CnaEdgeSet tree;
//        std::cout << "Number of edges: " << cnatree.size() << std::endl;

        for(auto edge: cnatree){
            std::pair<int, int> par = edge.first; 
            std::pair<int, int> child = edge.second; 
//            std::cout << par.first << "," << par.second << "->" << child.first << "," << child.second <<std::endl;
            CnaTree::CnaGenotype pGeno(par.first, par.second);
            CnaTree::CnaGenotype cGeno(child.first, child.second);
            CnaTree:: CnaEdge ed = std::make_pair(pGeno, cGeno);
            tree.insert(ed);
        }

        CnaTree cntree(tree);


        

//	std::cout << "About to enumerate..." << std::endl;

        GenotypeTree::GenotypeEdgeSetSet result;
        cntree.enumerateGenotypeTrees(result);
 


         std::vector<std::vector<std::pair<std::tuple<int, int,int, int>, std::tuple<int, int, int, int>>>> scriptT;

        for(auto gtree: result){
            std::vector<std::pair<std::tuple<int, int,int, int>, std::tuple<int, int, int, int>>> geno_tree;
            for(auto e: gtree){
                GenotypeTree::Genotype parent= e.first;
                GenotypeTree::Genotype child= e.second;
                geno_tree.push_back(std::make_pair( 
                    std::make_tuple(parent._x, parent._y, parent._xbar, parent._ybar),
                    std::make_tuple(child._x, child._y, child._xbar, child._ybar) 
                    ));


            }
            scriptT.push_back(geno_tree);

        }
   
//	std::cout << "Done" << std::endl;
        

        return scriptT;


 }





PYBIND11_MODULE(clonelib, m) {

    m.def("get_cna_trees", &genCNATrees, "Get CNA trees as a Python list of lists of tuples of int pair tuples (tree edges)");
    m.def("get_genotype_trees", &genGenoTrees, "Get genotype trees as a Python list of lists of tuples of int 4 tuples (tree edges)");

}
