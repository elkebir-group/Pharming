//
// Created by Mohammed El-Kebir on 8/14/23.
//

#include <lemon/connectivity.h>
#include "cnatree.h"
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/uniform_int.hpp>

CnaTree::CnaTree()
    : BaseTree(), _cnState(_T)
{
}

CnaTree::CnaTree(const CnaTree::CnaEdgeSet& edges)
    : BaseTree(), _cnState(_T)
{
  std::map<CnaTree::CnaGenotype, int> nodes;

  _k = edges.size() + 1;

  // 1. identify vertices
  for (const CnaTree::CnaEdge& edge: edges)
  {
    int state = nodes.size() + 1;
    nodes[edge.second] = state;
  }

  // 2. identify root
  CnaTree::CnaGenotype root(-1, -1);
  for (const CnaTree::CnaEdge& edge: edges)
  {
    if (nodes.count(edge.first) == 0)
    {
      root = edge.first;
      break;
    }
  }
  nodes[root] = 0;
  assert(root._x != -1 && root._y != -1);

  _stateToNode = NodeVector(nodes.size(), lemon::INVALID);

  // 3. add nodes
  char buf[1024];
  for (const auto& node: nodes)
  {
    assert(_stateToNode.size() > node.second);
    Node v = _T.addNode();
    _nodeToState[v] = node.second;
    _stateToNode[node.second] = v;
    _cnState[v] = node.first;
    snprintf(buf, 1024, "(%d,%d)", _cnState[v]._x, _cnState[v]._y);
    _label[v] = buf;
  }

  _root = _stateToNode[0];

  // 4. add edges
  for (const CnaTree::CnaEdge& edge: edges)
  {
    _T.addArc(_stateToNode[nodes[edge.first]], _stateToNode[nodes[edge.second]]);
  }

  initD(_root);
  assert(lemon::dag(_T));
}

CnaTree::CnaTree(int k)
    : BaseTree(k), _cnState(_T)
{
}

CnaTree::CnaTree(const IntVector& pi)
    : BaseTree(pi), _cnState(_T)
{
}

CnaTree::CnaTree(const CnaTree& other)
   : BaseTree()
   , _cnState(_T)
{
  _k = other._k;
  _stateToNode = NodeVector(_k, lemon::INVALID);
  lemon::digraphCopy(other._T, _T)
     .node(other._root, _root)
     .nodeMap(other._nodeToState, _nodeToState)
     .nodeMap(other._D, _D)
     .nodeMap(other._label, _label)
     .nodeMap(other._cnState, _cnState)
     .run();

  for (NodeIt v_i(_T); v_i != lemon::INVALID; ++v_i)
  {
    _stateToNode[_nodeToState[v_i]] = v_i;
  }
}

bool operator<(const CnaTree::CnaGenotype& lhs, const CnaTree::CnaGenotype& rhs)
{
  return lhs._x < rhs._x || (lhs._x == rhs._x && lhs._y < rhs._y);
}

CnaTree& CnaTree::operator=(const CnaTree& other)
{
  if (this != &other)
  {
    _k = other._k;
    _stateToNode = NodeVector(_k, lemon::INVALID);

    lemon::digraphCopy(other._T, _T)
        .node(other._root, _root)
        .nodeMap(other._nodeToState, _nodeToState)
        .nodeMap(other._label, _label)
        .nodeMap(other._D, _D)
        .nodeMap(other._s, _s)
        .nodeMap(other._cnState, _cnState)
        .run();

    for (NodeIt v_i(_T); v_i != lemon::INVALID; ++v_i)
    {
      _stateToNode[_nodeToState[v_i]] = v_i;
    }
  }
  return *this;
}

std::ostream& operator<<(std::ostream& out, const CnaTree& T)
{
  out << T.k() << " # nodes" << std::endl;
  // output vertices
  for (int i = 0; i < T.k(); ++i)
  {
    Node v_i = T.node(i);
    out << i << " " << T._label[v_i] << " "
        << T._cnState[v_i]._x
        << " " << T._cnState[v_i]._y;
    for (double s: T._s[v_i])
    {
      out << " " << s;
    }
    out << std::endl;
  }

  for (ArcIt a(T._T); a != lemon::INVALID; ++a)
  {
    out << T.state(T._T.source(a)) << " " << T.state(T._T.target(a)) << std::endl;
  }

  return out;
}

std::istream& operator>>(std::istream& in, CnaTree& T)
{
  std::string line;
  getline(in, line);
  std::stringstream ss(line);

  T._T.clear();


  int k = -1;
  ss >> k;
  if (k < 0)
  {
    throw std::runtime_error("Expected positive number of states");
  }

//  IntVector states(k);
  StringVector labels(k);
  IntVector x(k);
  IntVector y(k);
  DoubleMatrix prop(k);

  for (int i = 0; i < k; ++i)
  {
    getline(in, line);
    StringVector s;
    boost::split(s, line, boost::is_any_of(" \t"));

    if (s.size() < 4)
    {
      throw std::runtime_error("Expected at least four space-separated tokens");
    }

    int state = boost::lexical_cast<int>(s[0]);

    labels[state] = s[1];
    x[state] = boost::lexical_cast<int>(s[2]);
    y[state] = boost::lexical_cast<int>(s[3]);

    for (int idx = 4; idx < s.size(); ++idx)
    {
      prop[state].push_back(boost::lexical_cast<double>(s[idx]));
    }
  }

  // parse edge list
  IntVector pi(k, -1);
  for (int i = 0; i < k - 1; ++i)
  {
    getline(in, line);
    StringVector s;
    boost::split(s, line, boost::is_any_of(" \t"));

    if (s.size() != 2)
    {
      throw std::runtime_error("Expected two integer tokens");
    }

    int source = boost::lexical_cast<int>(s[0]);
    int target = boost::lexical_cast<int>(s[1]);

    if (!(0 <= source && source < k))
    {
      throw std::runtime_error("Invalid source node");
    }

    if (!(0 <= target && target < k))
    {
      throw std::runtime_error("Invalid target node");
    }

    pi[target] = source;
  }

  T = CnaTree(pi);

  for (int i = 0; i < T.k(); ++i)
  {
    Node v_i = T.node(i);
    T._cnState[v_i]._x = x[i];
    T._cnState[v_i]._y = y[i];
    T._label[v_i] = labels[i];
    T._s[v_i] = prop[i];
  }

  return in;
}

std::ostream& operator<<(std::ostream& out, const CnaTreeVector& vecT)
{
  out << vecT.size() << " #CNA trees" << std::endl;
  for (const CnaTree& T: vecT)
  {
    out << T;
  }

  return out;
}

std::istream& operator>>(std::istream& in, CnaTreeVector& vecT)
{
  g_lineNumber = 0;

  std::string line;
  getline(in, line);
  std::stringstream ss(line);

  int nrStateTrees = -1;
  ss >> nrStateTrees;

  if (nrStateTrees < 1)
  {
    throw std::runtime_error(getLineNumber() + "Error: incorrect number of state trees");
  }

  for (int i = 0; i < nrStateTrees; ++i)
  {
    CnaTree T;
    in >> T;
    vecT.push_back(T);
  }

  return in;
}

void CnaTree::enumerateGenotypeTrees(GenotypeTree::GenotypeEdgeSetSet& result) const
{
  for (int i = 0; i < _k; ++i)
  {
    Node mutationNode = node(i);
    GenotypeTree::Genotype mutationNodeGenotype(_cnState[mutationNode]._x, _cnState[mutationNode]._y, 0, 0);

    NodeVector mutationNodeChildren;
    for (OutArcIt a(_T, mutationNode); a != lemon::INVALID; ++a)
    {
      mutationNodeChildren.push_back(_T.target(a));
    }

    GenotypeTree::GenotypeEdgeSet fixedTree;
    fixPreMutation(_root, mutationNode, fixedTree);

    if (mutationNodeGenotype._x > 0)
    {
      GenotypeTree::GenotypeEdgeSet tmp = fixedTree;
      splitEnumerate(mutationNode, mutationNodeGenotype, mutationNodeChildren, mutationNode, true, true, tmp, result);
    }
    if (mutationNodeGenotype._y > 0)
    {
      GenotypeTree::GenotypeEdgeSet tmp = fixedTree;
      splitEnumerate(mutationNode, mutationNodeGenotype, mutationNodeChildren, mutationNode, false, true, tmp, result);
    }
  }
}

bool CnaTree::next(BoolVector& boolVector)
{
  // rightmost 0
  size_t i = 0;
  for (; i < boolVector.size(); ++i)
  {
    if (!boolVector[i]) break;
  }

  if (i == boolVector.size())
  {
    return false;
  }

  boolVector[i] = true;

  for (size_t j = 0; j < i; ++j)
  {
    boolVector[j] = false;
  }

  return true;
}

void CnaTree::fixPreMutation(Node u, Node mutationNode, GenotypeTree::GenotypeEdgeSet& tree) const
{
  if (u != mutationNode)
  {
    GenotypeTree::Genotype genotype_u(_cnState[u]._x, _cnState[u]._y, 0, 0);
    for (OutArcIt a(_T, u); a != lemon::INVALID; ++a)
    {
      Node v = _T.target(a);
      GenotypeTree::Genotype genotype_v(_cnState[v]._x, _cnState[v]._y, 0, 0);
      tree.insert({genotype_u, genotype_v});
      fixPreMutation(v, mutationNode, tree);
    }
  }
}

void CnaTree::splitEnumerate(Node u, GenotypeTree::Genotype genotype_u,
                             NodeVector children,
                             const Node mutationNode, const bool mut_x,
                             bool preMutation,
                             GenotypeTree::GenotypeEdgeSet& currentTree,
                             GenotypeTree::GenotypeEdgeSetSet& result) const
{
  const auto& cnState_u = _cnState[u];
  assert(genotype_u._x == cnState_u._x);
  assert(genotype_u._y == cnState_u._y);
  assert(0 <= genotype_u._xbar && genotype_u._xbar <= genotype_u._x);
  assert(0 <= genotype_u._ybar && genotype_u._ybar <= genotype_u._y);
  assert(genotype_u._xbar == 0 || genotype_u._ybar == 0);

  if (currentTree.size() == _k)
  {
    assert(result.count(currentTree) == 0);
    result.insert(currentTree);
    currentTree.clear();
    return;
  }

  if (u == mutationNode && preMutation)
  {
    GenotypeTree::GenotypeEdgeSet cpyCurrentTree = currentTree;
    BoolVector splitVector(children.size(), false);
    do
    {
      GenotypeTree::Genotype genotype_v = genotype_u;
      if (cnState_u._x > 0 && mut_x)
      {
        genotype_v._xbar = 1;
      }
      else
      {
        genotype_v._ybar = 1;
      }
      currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
      NodeVector childrenMut;
      NodeVector childrenNonMut;
      for (size_t idx = 0; idx < children.size(); ++idx)
      {
        if (splitVector[idx])
        {
          childrenMut.push_back(children[idx]);
        }
        else
        {
          childrenNonMut.push_back(children[idx]);
        }
      }
      splitEnumerate(u, genotype_u, childrenNonMut, mutationNode, mut_x, false, currentTree, result);
      splitEnumerate(u, genotype_v, childrenMut, mutationNode, mut_x, false, currentTree, result);
      currentTree = cpyCurrentTree;
    } while (next(splitVector));

//    // split
//    if (cnState_u._x > 0 && mut_x)
//    {
//      GenotypeTree::Genotype genotype_v = genotype_u;
//      genotype_v._xbar = 1;
//      currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
//      GenotypeTree::GenotypeEdgeSet cpyCurrentTree = currentTree;
//      splitEnumerate(u, genotype_v, mutationNode, mut_x, false, currentTree, result);
//      if (OutArcIt(_T, u) != lemon::INVALID)
//      {
//        currentTree = cpyCurrentTree;
//        splitEnumerate(u, genotype_u, mutationNode, mut_x, false, currentTree, result);
//      }
//
//      // TODO: partition children of u into two sets, those that will be children of the mutated node and those that will not
//    }
//    if (cnState_u._y > 0 && !mut_x)
//    {
//      GenotypeTree::Genotype genotype_v = genotype_u;
//      genotype_v._ybar = 1;
//      currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
//      GenotypeTree::GenotypeEdgeSet cpyCurrentTree = currentTree;
//      splitEnumerate(u, genotype_v, mutationNode, mut_x, false, currentTree, result);
//      if (OutArcIt(_T, u) != lemon::INVALID)
//      {
//        currentTree = cpyCurrentTree;
//        splitEnumerate(u, genotype_u, mutationNode, mut_x, false, currentTree, result);
//      }
//    }
  }
  else
  {
    for (Node v : children)
    {
      NodeVector children_v;
      for (OutArcIt a(_T, v); a != lemon::INVALID; ++a)
      {
        children_v.push_back(_T.target(a));
      }

      const auto& cnState_v = _cnState[v];
      if (genotype_u._xbar == 0 && genotype_u._ybar == 0)
      {
        GenotypeTree::Genotype genotype_v;
        genotype_v._x = cnState_v._x;
        genotype_v._y = cnState_v._y;
        genotype_v._xbar = genotype_v._ybar = 0;

        currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
        splitEnumerate(v, genotype_v, children_v, mutationNode, mut_x, preMutation, currentTree, result);
      }
      else if (genotype_u._xbar > 0)
      {
        assert(genotype_u._ybar == 0);
        GenotypeTree::Genotype genotype_v;
        genotype_v._x = cnState_v._x;
        genotype_v._y = cnState_v._y;
        genotype_v._ybar = 0;

        if (genotype_u._x == cnState_v._x)
        {
          genotype_v._xbar = genotype_u._xbar;
          currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
          splitEnumerate(v, genotype_v, children_v, mutationNode, mut_x, preMutation, currentTree, result);
        }
        else if (genotype_u._x < cnState_v._x)
        {
          // amplification
          if (genotype_u._xbar == genotype_u._x)
          {
            // all amplified copies must be mutated
            genotype_v._xbar = genotype_v._x;
            currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
            splitEnumerate(v, genotype_v, children_v, mutationNode, mut_x, preMutation, currentTree, result);
          }
          else
          {
            int delta = genotype_v._x - genotype_u._x;
            GenotypeTree::GenotypeEdgeSet cpyCurrentTree = currentTree;
            for (int l = 0; l <= delta; ++l)
            {
              genotype_v._xbar = genotype_u._xbar + l;
              currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
              splitEnumerate(v, genotype_v, children_v, mutationNode, mut_x, preMutation, currentTree, result);
              if (l < delta)
              {
                currentTree = cpyCurrentTree;
              }
            }
          }
        }
        else
        {
          // deletion
          if (genotype_u._xbar == genotype_u._x)
          {
            // all deleted copies must be mutated
            genotype_v._xbar = genotype_v._x;
            currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
            splitEnumerate(v, genotype_v, children_v, mutationNode, mut_x, preMutation, currentTree, result);
          }
          else
          {
            int genotype_u_non_mut_x = genotype_u._x - genotype_u._xbar;

            GenotypeTree::GenotypeEdgeSet cpyCurrentTree = currentTree;
            for (int xbar_v = 0; xbar_v <= genotype_v._x; ++xbar_v)
            {
              int xxbar_v = genotype_v._x - xbar_v;
              if (xxbar_v > genotype_u_non_mut_x) continue;
              if (xbar_v > genotype_u._xbar) continue;

              genotype_v._xbar = xbar_v;
              currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
              splitEnumerate(v, genotype_v, children_v, mutationNode, mut_x, preMutation, currentTree, result);
              if (xbar_v < genotype_v._x)
              {
                currentTree = cpyCurrentTree;
              }
            }
          }
        }
      }
      else
      {
        assert(genotype_u._ybar > 0);
        assert(genotype_u._xbar == 0);

        GenotypeTree::Genotype genotype_v;
        genotype_v._x = cnState_v._x;
        genotype_v._y = cnState_v._y;
        genotype_v._xbar = 0;

        if (genotype_u._y == cnState_v._y)
        {
          genotype_v._ybar = genotype_u._ybar;
          currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
          splitEnumerate(v, genotype_v, children_v, mutationNode, mut_x, preMutation, currentTree, result);
        }
        else if (genotype_u._y < cnState_v._y)
        {
          // amplification
          if (genotype_u._ybar == genotype_u._y)
          {
            // all amplified copies must be mutated
            genotype_v._ybar = genotype_v._y;
            currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
            splitEnumerate(v, genotype_v, children_v, mutationNode, mut_x, preMutation, currentTree, result);
          }
          else
          {
            GenotypeTree::GenotypeEdgeSet cpyCurrentTree = currentTree;

            int delta = genotype_v._y - genotype_u._y;
            for (int l = 0; l <= delta; ++l)
            {
              genotype_v._ybar = genotype_u._ybar + l;
              currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
              splitEnumerate(v, genotype_v, children_v, mutationNode, mut_x, preMutation, currentTree, result);
              if (l < delta)
              {
                currentTree = cpyCurrentTree;
              }
            }
          }
        }
        else
        {
          // deletion
          if (genotype_u._ybar == genotype_u._y)
          {
            // all deleted copies must be mutated
            genotype_v._ybar = genotype_v._y;
            currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
            splitEnumerate(v, genotype_v, children_v, mutationNode, mut_x, preMutation, currentTree, result);
          }
          else
          {
            GenotypeTree::GenotypeEdgeSet cpyCurrentTree = currentTree;
            int genotype_u_non_mut_y = genotype_u._y - genotype_u._ybar;

            for (int ybar_v = 0; ybar_v <= genotype_v._y; ++ybar_v)
            {
              int yybar_v = genotype_v._y - ybar_v;
              if (yybar_v > genotype_u_non_mut_y) continue;
              if (ybar_v > genotype_u._ybar) continue;

              genotype_v._ybar = ybar_v;
              currentTree.insert(GenotypeTree::GenotypeEdge(genotype_u, genotype_v));
              splitEnumerate(v, genotype_v, children_v, mutationNode, mut_x, preMutation, currentTree, result);
              if (ybar_v < genotype_v._y)
              {
                currentTree = cpyCurrentTree;
              }
            }
          }
        }
      }
    }
  }
}