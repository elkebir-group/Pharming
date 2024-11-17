//
// Created by Mohammed El-Kebir on 8/14/23.
//

#include "basetree.h"
#include <lemon/connectivity.h>
#include <iomanip>

BaseTree::BaseTree()
   : _k(0)
   , _T()
   , _root(lemon::INVALID)
   , _label(_T)
   , _stateToNode(_k, lemon::INVALID)
   , _nodeToState(_T)
   , _D(_T)
   , _s(_T)
{
}

BaseTree::BaseTree(int k)
   : _k(k)
   , _T()
   , _root(lemon::INVALID)
   , _label(_T)
   , _stateToNode(_k, lemon::INVALID)
   , _nodeToState(_T)
   , _D(_T)
   , _s(_T)
{
  IntVector pi(_k, -1);
  for (int i = 1; i < _k; ++i)
  {
    pi[i] = i - 1;
  }

  init(pi);
}

BaseTree::BaseTree(const IntVector& pi)
   : _k(pi.size())
   , _T()
   , _root(lemon::INVALID)
   , _label(_T, "")
   , _stateToNode(_k, lemon::INVALID)
   , _nodeToState(_T)
   , _D(_T)
   , _s(_T)
{
  init(pi);
}

//BaseTree::BaseTree(const BaseTree& other)
//   : _k(other._k)
//   , _T()
//   , _root(lemon::INVALID)
//   , _label(_T)
//   , _stateToNode(_k, lemon::INVALID)
//   , _nodeToState(_T)
//   , _D(_T)
//   , _s(_T)
//{
//  lemon::digraphCopy(other._T, _T)
//     .node(other._root, _root)
//     .nodeMap(other._nodeToState, _nodeToState)
//     .nodeMap(other._D, _D)
//     .nodeMap(other._label, _label)
//     .run();
//
//  for (NodeIt v_i(_T); v_i != lemon::INVALID; ++v_i)
//  {
//    _stateToNode[_nodeToState[v_i]] = v_i;
//  }
//}

//BaseTree& BaseTree::operator=(const BaseTree& other)
//{
//  if (this != &other)
//  {
//    _k = other._k;
//    _stateToNode = NodeVector(_k, lemon::INVALID);
//
//    lemon::digraphCopy(other._T, _T)
//       .node(other._root, _root)
//       .nodeMap(other._nodeToState, _nodeToState)
//       .nodeMap(other._label, _label)
//       .nodeMap(other._D, _D)
//       .nodeMap(other._s, _s)
//       .run();
//
//    for (NodeIt v_i(_T); v_i != lemon::INVALID; ++v_i)
//    {
//      _stateToNode[_nodeToState[v_i]] = v_i;
//    }
//  }
//  return *this;
//}

void BaseTree::writeEdgeList(std::ostream& out) const
{
  bool first = true;
  for (ArcIt a_ij(_T); a_ij != lemon::INVALID; ++a_ij)
  {
    if (first)
    {
      first = false;
    }
    else
    {
      out << " ; ";
    }

    Node v_i = _T.source(a_ij);
    Node v_j = _T.target(a_ij);

    out << _label[v_i] << " -> " << _label[v_j];
  }
}

void BaseTree::writeDOT(std::ostream& out) const
{
  out << "digraph T {" << std::endl;

  for (NodeIt v_i(_T); v_i != lemon::INVALID; ++v_i)
  {
    int i = _nodeToState[v_i];
    if (!isPresent(i)) continue;

    out << "\t" << i << " [label=\"" << _label[v_i];
    for (double s : _s[v_i])
    {
      out << std::setprecision(2) << "\\n" << s;
    }
    out << "\"]" << std::endl;
  }

  for (ArcIt a_ij(_T); a_ij != lemon::INVALID; ++a_ij)
  {
    int i = _nodeToState[_T.source(a_ij)];
    int j = _nodeToState[_T.target(a_ij)];

    out << "\t" << i << " -> " << j << std::endl;
  }

  out << "}" << std::endl;
}


void BaseTree::init(const IntVector& pi)
{
  // state 0 should be the root
  assert(pi[0] == -1);
  char buf[1024];

  _root = _T.addNode();
  _stateToNode[0] = _root;
  _nodeToState[_root] = 0;
  snprintf(buf, 1024, "%d", 0);
  _label[_root] = buf;

  // init nodes of _T
  for (int i = 1; i < _k; ++i)
  {
    Node v_i = _T.addNode();
    _stateToNode[i] = v_i;
    _nodeToState[v_i] = i;
    snprintf(buf, 1024, "%d", i);
    _label[v_i] = buf;
  }

  // init edges of _T
  for (int i = 1; i < _k; ++i)
  {
    int pi_i = pi[i];

    //assert(0 <= pi_i < _k);
    if (0 <= pi_i && pi_i < _k)
    {
      _T.addArc(_stateToNode[pi_i], _stateToNode[i]);
    }
  }

  initD(_root);

  assert(lemon::dag(_T));
}

void BaseTree::initD(Node v_i)
{
  IntSet& D_i = _D[v_i];
  D_i.clear();

  D_i.insert(_nodeToState[v_i]);

  for (OutArcIt a(_T, v_i); a != lemon::INVALID; ++a)
  {
    Node v_j = _T.target(a);

    initD(v_j);
    D_i.insert(_D[v_j].begin(), _D[v_j].end());
  }
}