//
// Created by Mohammed El-Kebir on 8/13/23.
//

#include <lemon/bfs.h>
#include "cnagraph.h"
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/uniform_int.hpp>

CnaGraph::Dictionary CnaGraph::_dict = CnaGraph::Dictionary();
std::vector<IntPairSet> CnaGraph::_keys;

const CnaTree::CnaEdgeSetSet& CnaGraph::getCnaTrees(const IntPairSet& L,
                                                    int root_x, int root_y)
{
  int max_x = 0;
  int max_y = 0;
  for (const IntPair& xy : L)
  {
    max_x = std::max(max_x, xy.first);
    max_y = std::max(max_y, xy.second);
  }

  CnaGraph G(max_x, max_y, root_x, root_y);

  if (_dict.find(L) == _dict.end())
  {
    G.enumerate(L);
    _dict[L] = G._result;
    for (size_t i = 0; i < G._result.size(); ++i)
    {
      _keys.push_back(L);
    }
  }

  return _dict[L];
}

CnaGraph::CnaGraph(int max_x, int max_y, int root_x, int root_y)
  : _maxCopyNumberX(std::max(max_x, 1))
  , _maxCopyNumberY(std::max(max_y, 1))
  , _rootX(root_x)
  , _rootY(root_y)
  , _G()
  , _x(_G, 0)
  , _y(_G, 0)
  , _type(_G)
  , _toNode(_maxCopyNumberX + 1,
            NodeVector(_maxCopyNumberY + 1, lemon::INVALID))
  , _root(lemon::INVALID)
  , _result()
{
  init();
//  writeDOT(std::cout);
}

void CnaGraph::init()
{
  assert(0 <= _rootX && _rootX <= _maxCopyNumberX);
  assert(0 <= _rootY && _rootY <= _maxCopyNumberY);

  // vertices
  for (int x = 0; x <= _maxCopyNumberX; ++x)
  {
    for (int y = 0; y <= _maxCopyNumberY; ++y)
    {
      Node v_xy = _G.addNode();
      _toNode[x][y] = v_xy;
      _x[v_xy] = x;
      _y[v_xy] = y;
    }
  }

  _root = _toNode[_rootX][_rootY];

  // edges
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    const int x = _x[v];
    const int y = _y[v];

    assert(0 <= x && x <= _maxCopyNumberX);
    assert(0 <= y && y <= _maxCopyNumberY);

    // amplification edges
    if (0 < x && x < _maxCopyNumberX)
    {
      addEdge(v, AMPLIFICATION, x+1, y);
    }
    if (0 < y && y < _maxCopyNumberY)
    {
      addEdge(v, AMPLIFICATION, x, y+1);
    }

    // deletion edges
    if (x > 0)
    {
      addEdge(v, DELETION, x-1, y);
    }
    if (y > 0)
    {
      // delete nonmutated y copy
      addEdge(v, DELETION, x, y-1);
    }
  }
}

void CnaGraph::writeDOT(std::ostream& out) const
{
  IntNodeMap level(_G, 0);
  lemon::bfs(_G).distMap(level).run(_root);

  out << "digraph G {" << std::endl;

  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    out << "\t" << _G.id(v) << " [label=\"("
        << _x[v] << "," << _y[v] << ")\\n" << level[v] << "\"]" << std::endl;
  }

  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node v = _G.source(a);
    Node w = _G.target(a);

    out << "\t" << _G.id(v) << " -> " << _G.id(w) << " [color=";
    switch (_type[a])
    {
      case AMPLIFICATION:
        out << "green";
        break;
      case DELETION:
        out << "red";
        break;
    }
    out << "]" << std::endl;
  }
  out << "}" << std::endl;
}

void CnaGraph::init(SubDigraph& subG,
                    SubDigraph& T,
                    ArcList& F)
{
  T.enable(_root);
  F.clear();

  for (OutArcIt a(_G, _root); a != lemon::INVALID; ++a)
  {
    F.push_back(a);
  }
}

int CnaGraph::enumerate(const IntPairSet& L)
{
  BoolNodeMap filterNodesT(_G, false);
  BoolArcMap filterArcsT(_G, false);
  SubDigraph T(_G, filterNodesT, filterArcsT);

  BoolNodeMap filterNodesG(_G, true);
  BoolArcMap filterArcsG(_G, true);
  SubDigraph subG(_G, filterNodesG, filterArcsG);

  ArcList F;
  init(subG, T, F);

  IntMatrix V(_maxCopyNumberX + 1,
              IntVector(_maxCopyNumberY + 1, 0));
  V[_rootX][_rootY] = 1;

  BoolMatrix LL(_maxCopyNumberX + 1,
                BoolVector(_maxCopyNumberY + 1, false));
  for (const IntPair& xy : L)
  {
    LL[xy.first][xy.second] = true;
  }

  _result.clear();

  IntPairSet leavesCopyStates;
  leavesCopyStates.insert(IntPair(_rootX, _rootY));
  IntPairSet verticesCopyStates;
  verticesCopyStates.insert(IntPair(_rootX, _rootY));

  grow(L, LL, subG, T, F, V, verticesCopyStates, leavesCopyStates);

  return _result.size();
}

bool CnaGraph::isValid(const IntPairSet& L,
                       const BoolMatrix& LL,
                       const SubDigraph& T,
                       const IntMatrix& V,
                       const IntPairSet& verticesCopyStates,
                       const IntPairSet& leavesCopyStates) const
{
  // check #1: all copy states are in T
  IntPairSet diff;

  std::set_difference(L.begin(), L.end(),
                      verticesCopyStates.begin(), verticesCopyStates.end(),
                      std::inserter(diff, diff.begin()));

  if (!diff.empty())
  {
    return false;
  }

  // check #2: all leaves are in L
  std::set_difference(leavesCopyStates.begin(), leavesCopyStates.end(),
                      L.begin(), L.end(),
                      std::inserter(diff, diff.begin()));
  if (!diff.empty())
  {
    return false;
  }

  return true;
}

bool CnaGraph::isValid(const SubDigraph& T) const
{
  if (!T.status(_root))
    return false;

  // inf sites on copy states
  // this boils down to checking that |V_(x,y)| = 1
  IntMatrix V(_maxCopyNumberX + 1, IntVector(_maxCopyNumberY + 1, 0));
  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    int x = _x[v];
    int y = _y[v];

    if (V[x][y] == 0)
    {
      V[x][y] = 1;
    }
    else
    {
      return false;
    }
  }

  return true;
}

bool CnaGraph::isValid(const SubDigraph& T, Arc a)
{
  Node w = T.target(a);

  assert(T.status(T.source(a)));
  assert(!T.status(a));
  assert(!T.status(w));

  T.enable(a);
  T.enable(w);

  bool res = isValid(T);

  T.disable(a);
  T.disable(w);

  return res;
}

void CnaGraph::writeDOT(const SubDigraph& T, std::ostream& out) const
{
  out << "digraph S {" << std::endl;

  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    out << "\t" << T.id(v) << " [label=\"("
        << _x[v] << "," << _y[v] << ")\"]" << std::endl;
  }

  for (SubArcIt a(T); a != lemon::INVALID; ++a)
  {
    out << "\t" << T.id(T.source(a)) << " -> " << T.id(T.target(a)) << std::endl;
  }

  out << "}" << std::endl;
}

void CnaGraph::writeDOT(const SubDigraph& T, const ArcList& F, std::ostream& out) const
{
  out << "digraph S {" << std::endl;

  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    out << "\t" << T.id(v) << " [label=\"("
        << _x[v] << "," << _y[v] << ")\"]" << std::endl;
  }

  for (ArcListIt it = F.begin(); it != F.end(); ++it)
  {
    Arc st = *it;
    Node t = _G.target(st);

    out << "\t" << T.id(t) << " [label=\"("
        << _x[t] << "," << _y[t] << ")\",color=red]" << std::endl;
  }

  for (SubArcIt a(T); a != lemon::INVALID; ++a)
  {
    out << "\t" << T.id(T.source(a)) << " -> " << T.id(T.target(a)) << std::endl;
  }

  for (ArcListIt it = F.begin(); it != F.end(); ++it)
  {
    Arc st = *it;
    out << "\t" << _G.id(_G.source(st)) << " -> " << _G.id(_G.target(st)) << " [color=red]" << std::endl;
  }

  out << "}" << std::endl;
}

void CnaGraph::finalize(const IntPairSet& L,
                        const BoolMatrix& LL,
                        const SubDigraph& T,
                        const IntMatrix& V)
{
//  static int idx = 0;
//  ++idx;

  // convert to pairs
  CnaPairNodeMap toPair(_G, CnaTree::CnaGenotype(0, 0));

  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    int x_v = _x[v];
    int y_v = _y[v];

    toPair[v]._x = x_v;
    toPair[v]._y = y_v;
  }

  CnaTree::CnaEdgeSet S;

  // short-circuit inner nodes not in LL
  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    int x_v = _x[v];
    int y_v = _y[v];

    if (LL[x_v][y_v])
    {
      // find parent
      SubInArcIt a(T, v);
      while (a != lemon::INVALID)
      {
        Node u = T.source(a);

        int x_u = _x[u];
        int y_u = _y[u];

        if (LL[x_u][y_u] || u == _root)
        {
          S.insert(std::make_pair(toPair[u], toPair[v]));
          break;
        }
        else
        {
          a = SubInArcIt(T, u);
        }
      }
    }
  }

  _result.insert(S);
}

void CnaGraph::grow(const IntPairSet& L,
                    const BoolMatrix& LL,
                    SubDigraph& G,
                    SubDigraph& T,
                    ArcList& F,
                    IntMatrix& V,
                    IntPairSet& verticesCopyStates,
                    IntPairSet& leavesCopyStates)
{
  if (isValid(L, LL, T, V, verticesCopyStates, leavesCopyStates))
  {
    // report
    finalize(L, LL, T, V);
//    writeDOT(T, std::cout);
  }
  else
  {
    ArcList FF;
    while (!F.empty())
    {
      assert(!F.empty());

      // pick first edge (u,v) from frontier and add to T
      Arc uv = F.back();
      F.pop_back();

      Node u = G.source(uv);
      Node v = G.target(uv);

      // Get target (x_v, y_v)
      int x_v = _x[v];
      int y_v = _y[v];

      assert(V[x_v][y_v] == 0);

      // update V
      ++V[x_v][y_v];
      {
        IntPair xy_u(_x[u], _y[u]);
        leavesCopyStates.erase(xy_u);
      }

      leavesCopyStates.insert(IntPair(x_v, y_v));
      verticesCopyStates.insert(IntPair(x_v, y_v));

      assert(T.status(u));
      assert(!T.status(v));
      assert(!T.status(uv));

      // add uv to T
      T.enable(v);
      T.enable(uv);

      ArcList newF = F;

      // remove arcs from F
      for (ArcListNonConstIt it = newF.begin(); it != newF.end();)
      {
        Arc st = *it;
        Node s = T.source(st);
        Node t = T.target(st);

        int x_t = _x[t];
        int y_t = _y[t];

        //
        // condition 1: s is a parent of v
        // condition 2: t is not a child of v and has the same copy states as v
        if (v == t || (v != s && x_v == x_t && y_v == y_t))
        {
          assert(T.status(s));
          it = newF.erase(it);
        }
        else
        {
          ++it;
        }
      }

      // push each arc vw where w not in V(T) onto F
      for (SubOutArcIt vw(G, v); vw != lemon::INVALID; ++vw)
      {
        Node w = G.target(vw);
        if (T.status(w))
          continue;

        int x_w = _x[w];
        int y_w = _y[w];

        if (V[x_w][y_w] == 0)
        {
          newF.push_back(vw);
        }
      }

//      static int idx = 0;
//      std::cout << ++idx << std::endl;
//      writeDOT(T, newF, std::cout);

      grow(L, LL, G, T, newF, V, verticesCopyStates, leavesCopyStates);

      G.disable(uv);

      T.disable(uv);
      T.disable(v);
      --V[x_v][y_v];
      if (V[x_v][y_v] == 0)
      {
        leavesCopyStates.erase(IntPair(x_v, y_v));
        verticesCopyStates.erase(IntPair(x_v, y_v));
      }

      FF.push_back(uv);
    }

    for (ArcListRevIt it = FF.rbegin(); it != FF.rend(); ++it)
    {
      Arc a = *it;
      assert(!G.status(a));

      F.push_back(*it);
      G.enable(a);
    }
  }
}

std::istream& operator>>(std::istream& in, CnaGraph::Dictionary& dict)
{
  g_lineNumber = 0;
  std::string line;
  getline(in, line);
  std::stringstream ss(line);

  CnaGraph::Dictionary newDict;

  int nrCopyNumberStates = -1;
  ss >> nrCopyNumberStates;
  if (nrCopyNumberStates < 0)
  {
    throw std::runtime_error(getLineNumber() + "Error: invalid number of copy number states");
  }

  for (int i = 0; i < nrCopyNumberStates; ++i)
  {
    IntPairSet L;
    StringVector s;

    getline(in, line);
    boost::split(s, line, boost::is_any_of(" "));

    for (const std::string& xy_str : s)
    {
      int x = -1;
      int y = -1;
      if (sscanf(xy_str.c_str(), "%d,%d", &x, &y) != 2 || x < 0 || y < 0)
      {
        throw std::runtime_error(getLineNumber() + "Error: invalid copy number '" + xy_str + "'");
      }
      L.insert(IntPair(x, y));
    }

    if (newDict.count(L) > 0)
    {
      throw std::runtime_error(getLineNumber() + "Error: copy number states '" + line + "' already present");
    }

    int nrStateTrees = -1;
    getline(in, line);
    ss.str(line);
    ss >> nrStateTrees;
    if (nrStateTrees < 0)
    {
      throw std::runtime_error(getLineNumber() + "Error: invalid number of CNA trees");
    }

    for (int j = 0; j < nrStateTrees; ++j)
    {
      getline(in, line);
      boost::split(s, line, boost::is_any_of(" "));

      CnaTree::CnaStateVector xyVector;
      for (const std::string& xy_str : s)
      {
        int x = -1;
        int y = -1;
        if (sscanf(xy_str.c_str(), "%d,%d", &x, &y) != 2 || x < 0 || y < 0)
        {
          throw std::runtime_error(getLineNumber() + "Error: invalid copy number '" + xy_str + "'");
        }
        xyVector.push_back(CnaTree::CnaGenotype(x, y));
      }

      if (xyVector.size() % 2 != 0)
      {
        throw std::runtime_error(getLineNumber() + "Error: odd number of xy triples '" + line + "'");
      }

      CnaTree::CnaEdgeSet cnaTree;
      for (int jj = 0; jj < xyVector.size() / 2; ++jj)
      {
        cnaTree.insert(CnaTree::CnaEdge(xyVector[2 * jj], xyVector[2 * jj + 1]));
      }

      if (newDict[L].count(cnaTree) > 0)
      {
        throw std::runtime_error(getLineNumber() + "Error: duplicate CNA tree '" + line + "'");
      }
      newDict[L].insert(cnaTree);
    }
  }

  dict = newDict;

  return in;
}

std::ostream& operator<<(std::ostream& out, const CnaGraph::Dictionary& dict)
{
  out << dict.size() << " #copynumber_states" << std::endl;
  for (const auto& kv : dict)
  {
    bool first = true;
    for (const IntPair& xy : kv.first)
    {
      if (first)
        first = false;
      else
        out << " ";
      out << xy.first << "," << xy.second;
    }
    out << std::endl;

    out << kv.second.size() << " #CNA trees" << std::endl;
    for (const CnaTree::CnaEdgeSet& S : kv.second)
    {
      first = true;
      for (const CnaTree::CnaEdge& edge : S)
      {
        if (first)
          first = false;
        else
          out << " ";

        out << edge.first._x << "," << edge.first._y << " ";
        out << edge.second._x << "," << edge.second._y;
      }
      out << std::endl;
    }
  }

  return out;
}

void CnaGraph::writeCnaTrees(std::ostream& out)
{
  out << _dict;
}

void CnaGraph::readCnaTrees(std::istream& in)
{
  Dictionary newDict;
  in >> newDict;

  _dict.insert(newDict.begin(), newDict.end());
  _keys.clear();
  for (const auto kv : _dict)
  {
    for (size_t i = 0; i < kv.second.size(); ++i)
    {
      _keys.push_back(kv.first);
    }
  }
}

CnaTree CnaGraph::sampleCnaTree()
{
  boost::random::uniform_int_distribution<> intDist(0, _keys.size() - 1);
  int sample = intDist(g_rng);
  const auto& L = _keys[sample];
  assert(_dict.count(L) > 0);

  boost::random::uniform_int_distribution<> intDist2(0, _dict[L].size() - 1);

  int sample2 = intDist2(g_rng);
  for (const auto& cnaTree : _dict[L])
  {
    if (sample2-- == 0)
    {
      return CnaTree(cnaTree);
    }
  }
  assert(false);
  return CnaTree();
}