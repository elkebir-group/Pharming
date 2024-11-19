//
// Created by Mohammed El-Kebir on 8/14/23.
//

#ifndef BASETREE_H
#define BASETREE_H

#include "utils.h"

/// This class models a base tree
class BaseTree
{
public:
  /// Default constructor
  BaseTree();

  /// Constructor
  ///
  /// @param k Number of states
  BaseTree(int k);

    // Delete the copy constructor and copy assignment operator
//    BaseTree(const BaseTree&) = delete;
//    BaseTree& operator=(const BaseTree&) = delete;
//
//    // Provide move constructor and move assignment operator
//    BaseTree(BaseTree&&) noexcept = default;
//    BaseTree& operator=(BaseTree&&) noexcept = default;
  /// Copy constructor
  ///
  /// @param other Other state tree
//  BaseTree(const BaseTree& other);

  /// Assignment operator
  ///
  /// @param other Other state tree
//  BaseTree& operator=(const BaseTree& other);

  /// Constructor
  ///
  /// @param pi Vector of parental states
  BaseTree(const IntVector& pi);

  /// Return number of distinct states
  int k() const
  {
    return _k;
  }

  /// Return parent state; -1 is returned if i is the root state,
  /// -2 is returned if state is absent from the state tree
  ///
  /// @param i State
  int parent(int i) const
  {
    assert(0 <= i && i < _k);
    Node v_i = _stateToNode[i];
    Arc a = InArcIt(_T, v_i);
    if (a == lemon::INVALID)
    {
      if (v_i == _root)
      {
        return -1;
      }
      else
      {
        return -2;
      }
    }
    else
    {
      return _nodeToState[_T.source(a)];
    }
  }

  /// Return number of vertices
  int numVertices() const
  {
    return lemon::countNodes(_T);
  }

  /// Return whether specified state is present
  ///
  /// @param i State
  bool isPresent(int i) const
  {
    assert(0 <= i && i < _k);
    return parent(i) != -2;
  }

  /// Return whether specified states are in a parent-child relationship
  ///
  /// @param i Parent state
  /// @param j Child state
  bool isParent(int i, int j) const
  {
    assert(0 <= i && i < _k);
    assert(1 <= j && j < _k);

    Node v_i = _stateToNode[i];
    Node v_j = _stateToNode[j];

    Arc a = InArcIt(_T, v_j);
    if (a == lemon::INVALID)
    {
      return false;
    }
    else
    {
      return v_i == _T.source(a);
    }
  }

  /// Return whether specified states are siblings
  ///
  /// @param i State
  /// @param j State
  bool areSiblings(int i, int j) const
  {
    assert(1 <= i && i < _k);
    assert(1 <= j && j < _k);

    Node v_i = _stateToNode[i];
    Node v_j = _stateToNode[j];

    Arc a1 = InArcIt(_T, v_i);
    if (a1 == lemon::INVALID)
      return false;

    Node v_pi_i = _T.source(a1);

    Arc a2 = InArcIt(_T, v_j);
    if (a2 == lemon::INVALID)
      return false;
    Node v_pi_j = _T.source(a2);

    return v_pi_i == v_pi_j;
  }

  /// Return whether specified states are in an ancestor-descendant relationship
  ///
  /// @param i Ancestor state
  /// @param j Descendant state
  bool isAncestor(int i, int j) const
  {
    assert(0 <= i && i < _k);
    assert(0 <= j && j < _k);

    Node v_i = _stateToNode[i];

    Node v = _stateToNode[j];
    while (v != v_i && v != _root)
    {
      v = _T.source(InArcIt(_T, v));
    }

    return v == v_i;
  }

  /// Return whether specified states are in a descendant-ancestor relationship
  ///
  /// @param i Descendant state
  /// @param j Ancestor state
  bool isDescendant(int i, int j) const
  {
    return isAncestor(j, i);
  }

  /// Return whether specified states occur in distinct branches
  ///
  /// @param i State
  /// @param j State
  bool isIncomparable(int i, int j) const
  {
    return !isAncestor(i, j) && !isAncestor(j, i);
  }

  /// Return set of descendant states of specified state
  ///
  /// @param i State
  const IntSet& D(int i) const
  {
    assert(0 <= i && i < _k);
    return _D[_stateToNode[i]];
  }

  /// Return label of specified state
  ///
  /// @param i State
  const std::string& label(int i) const
  {
    assert(0 <= i && i < _k);
    return _label[_stateToNode[i]];
  }

  /// Set label of specified state
  ///
  /// @param i State
  /// @param l Label
  void setLabel(int i, const std::string& l)
  {
    _label[_stateToNode[i]] = l;
  }

  /// Return tree
  const Digraph& T() const
  {
    return _T;
  }

  /// Return the state of the specified node
  ///
  /// @param v_i Node
  int state(Node v_i) const
  {
    assert(v_i != lemon::INVALID);
    return _nodeToState[v_i];
  }

  /// Return the node of the specified state
  ///
  /// @param i State
  Node node(int i) const
  {
    assert(0 <= i && i < _k);
    return _stateToNode[i];
  }

  /// Return root state
  int rootState() const
  {
    return _nodeToState[_root];
  }

  /// Write edge list
  ///
  /// @param out Output stream
  void writeEdgeList(std::ostream& out) const;

  /// Write state tree in DOT format
  ///
  /// @param out Output stream
  virtual void writeDOT(std::ostream& out) const;

  /// Return the number of samples
  int getNrSamples() const
  {
    assert(_root != lemon::INVALID);
    return _s[_root].size();
  }

  /// Get mixture proportions
  ///
  /// @param i State
  const DoubleVector& getMixtureProportion(int i) const
  {
    assert(0 <= i && i < _k);
    Node v_i = _stateToNode[i];

    return _s[v_i];
  }

  /// Get mixture proportion
  ///
  /// @param i State
  /// @param p Sample
  double getMixtureProportion(int i, int p) const
  {
    Node v_i = node(i);
    assert(0 <= p && p < _s[v_i].size());
    return _s[v_i][p];
  }

  /// Set mixture proportion
  ///
  /// @param i State
  /// @param s Proportion
  void setMixtureProportion(int i, double s)
  {
    assert(0 <= i && i < _k);
    Node v_i = _stateToNode[i];

    _s[v_i].push_back(s);
  }

  /// Increment mixture proportion
  ///
  /// @param p Sample
  /// @param i State
  /// @param s Proportion delta
  void incrementMixtureProportion(int p, int i, double s)
  {
    assert(0 <= i && i < _k);
    Node v_i = _stateToNode[i];

    _s[v_i][p] += s;
  }

  /// Reset mixture proportions to 0
  ///
  /// @param nrSamples Number of samples
  void resetMixtureProportions(int nrSamples)
  {
    for (int j = 0; j < _k; ++j)
    {
      Node v_j = _stateToNode[j];
      if (v_j != lemon::INVALID)
      {
        _s[v_j] = DoubleVector(nrSamples, 0);
      }
    }
  }

  /// Clear mixture proportions
  void clearMixtureProportions()
  {
    for (NodeIt v_i(_T); v_i != lemon::INVALID; ++v_i)
    {
      _s[v_i].clear();
    }
  }

  /// TODO: Sample mixture proportions
  ///
  /// @param m Number of samples
  void sampleMixtureProportions(int m);

protected:
  /// Number of states
  int _k;
  /// Genotype tree
  Digraph _T;
  /// Root node
  Node _root;
  /// Node labeling
  StringNodeMap _label;
  /// State to node mapping
  NodeVector _stateToNode;
  /// Node to state mapping
  IntNodeMap _nodeToState;
  /// Descendant state set node map
  IntSetNodeMap _D;
  /// Mixture proportions;
  DoubleVectorNodeMap _s;

  /// Initialize descendant sets of subtree rooted at the specified node
  ///
  /// @param v_i node
  void initD(Node v_i);

  /// Initialize state tree using specified parental state vector
  ///
  /// @param pi Vector of parental states
  void init(const IntVector& pi);
};

#endif //BASETREE_H
