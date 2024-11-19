/*
 * genotypegraph.h
 *
 *  Created on: 11-jan-2016
 *      Author: M. El-Kebir
 */

#ifndef GENOTYPEGRAPH_H
#define GENOTYPEGRAPH_H

#include <lemon/adaptors.h>
#include <vector>
#include "utils.h"
#include "genotypetree.h"

/// This class models a directed graph containing all genotype trees
class GenotypeGraph
{
public:
  /// Constructor
  ///
  /// @param max_x Maximum number of maternal copies
  /// @param max_y Maximum number of paternal copies
  /// @param root_x Root number of maternal copies
  /// @param root_y Root number of paternal copies
  GenotypeGraph(int max_x, int max_y, int root_x, int root_y);
  
  /// Mutation type
  typedef enum { MUTATION, AMPLIFICATION, DELETION } EdgeType;
 
  /// Edge type map
  typedef Digraph::ArcMap<EdgeType> EdgeTypeArcMap;
  
  /// Write DOT visualization of genotype graph
  ///
  /// @param out Output stream
  void writeDOT(std::ostream& out) const;
  
  /// Enumerate genotype trees and return the number of enumerated state trees
  ///
  /// @param L Copy number states that must be present
  int enumerate(const IntPairSet& L);
  
private:
  typedef std::vector<NodeMatrix> Node3Matrix;
  typedef std::vector<Node3Matrix> Node4Matrix;
  
  typedef std::list<Arc> ArcList;
  typedef ArcList::const_iterator ArcListIt;
  typedef ArcList::iterator ArcListNonConstIt;
  typedef ArcList::const_reverse_iterator ArcListRevIt;
  
  /// Initialize state graph
  void init();
  
  /// Initialize enumeration by adding root node to T and setting the frontier F
  ///
  /// @param subG Subgraph of state graph
  /// @param T Partial state tree
  /// @param F Frontier
  void init(SubDigraph& subG,
            SubDigraph& T,
            ArcList& F);
  
  /// Report enumerated state tree
  ///
  /// @param L
  /// @param LL
  /// @param T Enumerated state tree
  /// @param V
  bool finalize(const IntPairSet& L,
                const BoolMatrix& LL,
                const SubDigraph& T,
                const IntMatrix& V,
                const Arc& mutationEdge);
  
  /// Write partial state tree in DOT format
  ///
  /// @param T Partial state tree
  /// @param out Output stream
  void writeDOT(const SubDigraph& T,
                std::ostream& out) const;
  
  /// Write partial state tree and frontier in DOT format
  ///
  /// @param T Partial state tree
  /// @param F Frontier
  /// @param out Output stream
  void writeDOT(const SubDigraph& T,
                const ArcList& F,
                std::ostream& out) const;
  
  /// Add edge to G
  ///
  /// @param v Source node
  /// @param type Edge type (mutation, amplification, deletion)
  /// @param x Number of maternal copies of target node
  /// @param y Number of paternal copies of target node
  /// @param xbar Number of mutated maternal copies of target node
  /// @param ybar Number of mutated paternal copies of target node
  void addEdge(Node v, EdgeType type, int x, int y, int xbar, int ybar)
  {
    Node w = _toNode[x][y][xbar][ybar];
    Arc vw = _G.addArc(v, w);
    _type[vw] = type;
  }
  
  /// Return whether proposed partial genotype tree is valid
  ///
  /// @param L
  /// @param LL
  /// @param T Proposed partial state tree
  /// @param V
  /// @param verticesCopyStates Copy number states of vertices
  /// @param leavesCopyStates Copy number states of leaves
  bool isValid(const IntPairSet& L,
               const BoolMatrix& LL,
               const SubDigraph& T,
               const IntMatrix& V,
               const IntPairSet& verticesCopyStates,
               const IntPairSet& leavesCopyStates) const;
  
  /// Return whether partial genotype tree is valid
  ///
  /// @param T Partial genotype tree
  bool isValid(const SubDigraph& T) const;
  
  /// Return whether partial genotype tree remains valid upon edge of specified edge
  ///
  /// @param T Partial genotype tree
  /// @param a Arc
  bool isValid(const SubDigraph& T,
               Arc a);
  
  /// Grow partial genotype tree
  ///
  /// @param L Set of required copy number states
  /// @param LL Set of required copy number states in matrix form
  /// @param G Graph
  /// @param T Current partial state tree
  /// @param F Frontier of edges to expand into
  /// @param V Nodes of G in matrix form
  /// @param mutationEdge Mutation edge
  /// @param verticesCopyStates Copy number states of vertices
  /// @param leavesCopyStates Copy number states of leaves
  void grow(const IntPairSet& L,
            const BoolMatrix& LL,
            SubDigraph& G,
            SubDigraph& T,
            ArcList& F,
            IntMatrix& V,
            Arc& mutationEdge,
            IntPairSet& verticesCopyStates,
            IntPairSet& leavesCopyStates);

private:
  typedef Digraph::NodeMap<GenotypeTree::Genotype> GenotypeNodeMap;
  
private:
  /// Maximum maternal copy number
  const int _maxCopyNumberX;
  /// Maximum paternal copy number
  const int _maxCopyNumberY;
  /// Root maternal copy number
  const int _rootX;
  /// Root paternal copy number
  const int _rootY;
  /// Genotype graph
  Digraph _G;
  /// Maternal copy number node map
  IntNodeMap _x;
  /// Paternal copy number node map
  IntNodeMap _y;
  /// Mutated maternal copy number node map
  IntNodeMap _xbar;
  /// Mutated paternal copy number node map
  IntNodeMap _ybar;
  /// Mutation type edge map
  EdgeTypeArcMap _type;
  /// Genotype to node map
  Node4Matrix _toNode;
  /// Root node
  Node _root;
  /// Enumerated genotype trees
  GenotypeTree::GenotypeEdgeSetSet _result;
  
public:
  static const GenotypeTree::GenotypeEdgeSetSet& getStateTrees(const IntPairSet& L,
                                                               int root_x, int root_y);
  
  static void writeGenotypeTrees(std::ostream& out);
  
  static void readGenotypeTrees(std::istream& in);
  
  static void partition(const int root_x,
                        const int root_y,
                        const GenotypeTree::GenotypeEdgeSet& S,
                        const GenotypeTree::Genotype& vertex,
                        bool& preMutationFlag,
                        GenotypeTree::Genotype& mutation,
                        GenotypeTree::GenotypeSet& preMutation,
                        GenotypeTree::GenotypeSet& postMutation);
  
  static void partition(const int root_x,
                        const int root_y,
                        const GenotypeTree::GenotypeEdgeSet& S,
                        GenotypeTree::GenotypeSet& preMutation,
                        GenotypeTree::GenotypeSet& postMutation,
                        GenotypeTree::Genotype& mutation);

  static void print(const GenotypeTree::GenotypeEdgeSet& S,
                    std::ostream& out)
  {
    for (GenotypeTree::GenotypeEdgeSetIt it = S.begin(); it != S.end(); ++it)
    {
      out << "( (" << it->first._x << "," << it->first._y << "," << it->first._xbar << "," << it->first._ybar << ") , ("
          << it->second._x << "," << it->second._y << "," << it->second._xbar << "," << it->second._ybar << ") )  ";
    }
  }
  
  static void print(const GenotypeTree::GenotypeEdgeSetSet& setS,
                    std::ostream& out)
  {
    for (GenotypeTree::GenotypeEdgeSetSetIt it1 = setS.begin(); it1 != setS.end(); ++it1)
    {
      const GenotypeTree::GenotypeEdgeSet& S = *it1;
      print(S, std::cout);
    }
  }
  
public:
  typedef std::map<IntPairSet, GenotypeTree::GenotypeEdgeSetSet> Dictionary;
  
public:
  static Dictionary _dict;
};
  
std::ostream& operator<<(std::ostream& out,
                         const GenotypeGraph::Dictionary& dict);

std::istream& operator>>(std::istream& in,
                         GenotypeGraph::Dictionary& dict);

#endif // GENOTYPEGRAPH_H
