//
// Created by Mohammed El-Kebir on 8/13/23.
//

#ifndef CNAGRAPH_H
#define CNAGRAPH_H

#include <lemon/adaptors.h>
#include <vector>
#include "utils.h"
#include "cnatree.h"

class CnaGraph
{
public:
  /// Constructor
  ///
  /// @param max_x Maximum number of maternal copies
  /// @param max_y Maximum number of paternal copies
  /// @param root_x Root number of maternal copies
  /// @param root_y Root number of paternal copies
  CnaGraph(int max_x, int max_y, int root_x, int root_y);

  /// Mutation type
  typedef enum { AMPLIFICATION, DELETION } EdgeType;

  /// Edge type map
  typedef Digraph::ArcMap<EdgeType> EdgeTypeArcMap;

  /// Write DOT visualization of CNA graph
  ///
  /// @param out Output stream
  void writeDOT(std::ostream& out) const;

  /// Enumerate CNA trees and return the number of enumerated trees
  ///
  /// @param L Copy number states that must be present
  int enumerate(const IntPairSet& L);

private:
  typedef std::list<Arc> ArcList;
  typedef ArcList::const_iterator ArcListIt;
  typedef ArcList::iterator ArcListNonConstIt;
  typedef ArcList::const_reverse_iterator ArcListRevIt;

  /// Initialize CNA graph
  void init();

  /// Initialize enumeration by adding root node to T and setting the frontier F
  ///
  /// @param subG Subgraph of CNA graph
  /// @param T Partial CNA tree
  /// @param F Frontier
  void init(SubDigraph& subG,
            SubDigraph& T,
            ArcList& F);

  /// Report enumerated CNA tree
  ///
  /// @param L Set of required copy number states
  /// @param LL Set of required copy number states in matrix form
  /// @param T Enumerated CNA tree
  /// @param V Nodes of G in matrix form
  void finalize(const IntPairSet& L,
                const BoolMatrix& LL,
                const SubDigraph& T,
                const IntMatrix& V);

  /// Write partial CNA tree in DOT format
  ///
  /// @param T Partial CNA tree
  /// @param out Output stream
  void writeDOT(const SubDigraph& T,
                std::ostream& out) const;

  /// Write partial CNA tree and frontier in DOT format
  ///
  /// @param T Partial CNA tree
  /// @param F Frontier
  /// @param out Output stream
  void writeDOT(const SubDigraph& T,
                const ArcList& F,
                std::ostream& out) const;

  /// Add edge to G
  ///
  /// @param v Source node
  /// @param type Edge type (amplification, deletion)
  /// @param x Number of maternal copies of target node
  /// @param y Number of paternal copies of target node
  void addEdge(Node v, EdgeType type, int x, int y)
  {
    Node w = _toNode[x][y];
    Arc vw = _G.addArc(v, w);
    _type[vw] = type;
  }

  /// Return whether proposed partial CNA tree is valid
  ///
  /// @param L Set of required copy number states
  /// @param LL Set of required copy number states in matrix form
  /// @param T Proposed partial CNA tree
  /// @param V Nodes of G in matrix form
  /// @param verticesCopyStates Copy number states of vertices
  /// @param leavesCopyStates Copy number states of leaves
  bool isValid(const IntPairSet& L,
               const BoolMatrix& LL,
               const SubDigraph& T,
               const IntMatrix& V,
               const IntPairSet& verticesCopyStates,
               const IntPairSet& leavesCopyStates) const;

  /// Return whether partial CNA tree is valid
  ///
  /// @param T Partial CNA tree
  bool isValid(const SubDigraph& T) const;

  /// Return whether partial CNA tree remains valid upon addition of specified edge
  ///
  /// @param T Partial CNA tree
  /// @param a Arc
  bool isValid(const SubDigraph& T,
               Arc a);

  /// Grow partial CNA tree
  ///
  /// @param L Set of required copy number states
  /// @param LL Set of required copy number states in matrix form
  /// @param G Graph
  /// @param T Current partial state tree
  /// @param F Frontier of edges to expand into
  /// @param V Nodes of G in matrix form
  /// @param verticesCopyStates Copy number states of vertices
  /// @param leavesCopyStates Copy number states of leaves
  void grow(const IntPairSet& L,
            const BoolMatrix& LL,
            SubDigraph& G,
            SubDigraph& T,
            ArcList& F,
            IntMatrix& V,
            IntPairSet& verticesCopyStates,
            IntPairSet& leavesCopyStates);

private:
  typedef Digraph::NodeMap<CnaTree::CnaGenotype> CnaPairNodeMap;

private:
  /// Maximum maternal copy number
  const int _maxCopyNumberX;
  /// Maximum paternal copy number
  const int _maxCopyNumberY;
  /// Root maternal copy number
  const int _rootX;
  /// Root paternal copy number
  const int _rootY;
  /// CNA graph
  Digraph _G;
  /// Maternal copy number node map
  IntNodeMap _x;
  /// Paternal copy number node map
  IntNodeMap _y;
  /// Mutation type edge map
  EdgeTypeArcMap _type;
  /// Copy number to node map
  NodeMatrix _toNode;
  /// Root node
  Node _root;
  /// Enumerated CNA trees
  CnaTree::CnaEdgeSetSet _result;

public:
  /// @param L Set of required copy number states
  /// @param root_x Root number of maternal copies
  /// @param root_y Root number of paternal copies
  static const CnaTree::CnaEdgeSetSet& getCnaTrees(const IntPairSet& L,
                                                   int root_x,
                                                   int root_y);

  static void writeCnaTrees(std::ostream& out);

  static void readCnaTrees(std::istream& in);

  static void print(const CnaTree::CnaEdgeSet& S,
                    std::ostream& out)
  {
    for (CnaTree::CnaEdgeSetIt it = S.begin(); it != S.end(); ++it)
    {
      out << "( (" << it->first._x << "," << it->first._y << ") , ("
          << it->second._x << "," << it->second._y << ") )  ";
    }
  }

  static void print(const CnaTree::CnaEdgeSetSet& setS,
                    std::ostream& out)
  {
    for (CnaTree::CnaEdgeSetSetIt it1 = setS.begin(); it1 != setS.end(); ++it1)
    {
      const CnaTree::CnaEdgeSet& S = *it1;
      print(S, std::cout);
    }
  }

public:
  typedef std::map<IntPairSet, CnaTree::CnaEdgeSetSet> Dictionary;
  static Dictionary _dict;

  static CnaTree sampleCnaTree();

private:
  static std::vector<IntPairSet> _keys;
};

#endif //CNAGRAPH_H
