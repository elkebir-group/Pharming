/*
 * genotypegraph.cpp
 *
 *  Created on: 11-jan-2016
 *      Author: M. El-Kebir
 */

#include <lemon/bfs.h>
#include "genotypegraph.h"

GenotypeGraph::Dictionary GenotypeGraph::_dict = GenotypeGraph::Dictionary();

const GenotypeTree::GenotypeEdgeSetSet& GenotypeGraph::getStateTrees(const IntPairSet& L,
                                                                     int root_x, int root_y)
{
  int max_x = 0;
  int max_y = 0;
  for (const IntPair& xy : L)
  {
    max_x = std::max(max_x, xy.first);
    max_y = std::max(max_y, xy.second);
  }
  
  GenotypeGraph G(max_x, max_y, root_x, root_y);
//  G.writeDOT(std::cout);
  
  if (_dict.find(L) == _dict.end())
  {
    G.enumerate(L);
    _dict[L] = G._result;
  }
  
//  std::cout << _dict << std::endl;
  
  return _dict[L];
}
  
GenotypeGraph::GenotypeGraph(int max_x, int max_y,
                             int root_x, int root_y)
  : _maxCopyNumberX(std::max(max_x, root_x))
  , _maxCopyNumberY(std::max(max_y, root_y))
  , _rootX(root_x)
  , _rootY(root_y)
  , _G()
  , _x(_G, 0)
  , _y(_G, 0)
  , _xbar(_G, 0)
  , _ybar(_G, 0)
  , _type(_G)
  , _toNode(_maxCopyNumberX + 1,
            Node3Matrix(_maxCopyNumberY + 1,
                        NodeMatrix(_maxCopyNumberX + 1,
                                   NodeVector(_maxCopyNumberY + 1, lemon::INVALID))))
  , _root(lemon::INVALID)
  , _result()
{
  init();
//  writeDOT(std::cout);
}
  
void GenotypeGraph::init()
{
  // vertices
  for (int x = 0; x <= _maxCopyNumberX; ++x)
  {
    for (int y = 0; y <= _maxCopyNumberY; ++y)
    {
      for (int xbar = 0; xbar <= x; ++xbar)
      {
        Node v_xy_xbar0 = _G.addNode();
        _toNode[x][y][xbar][0] = v_xy_xbar0;
        _x[v_xy_xbar0] = x;
        _y[v_xy_xbar0] = y;
        _xbar[v_xy_xbar0] = xbar;
        _ybar[v_xy_xbar0] = 0;
      }
  
      for (int ybar = 1; ybar <= y; ++ybar)
      {
        Node v_xy_0ybar = _G.addNode();
        _toNode[x][y][0][ybar] = v_xy_0ybar;
        _x[v_xy_0ybar] = x;
        _y[v_xy_0ybar] = y;
        _xbar[v_xy_0ybar] = 0;
        _ybar[v_xy_0ybar] = ybar;
      }
    }
  }
  
  _root = _toNode[_rootX][_rootY][0][0];
  
  // edges
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    const int x = _x[v];
    const int y = _y[v];
    const int xbar = _xbar[v];
    const int ybar = _ybar[v];
    
    assert(0 <= x && x <= _maxCopyNumberX);
    assert(0 <= y && y <= _maxCopyNumberY);
    assert(0 <= xbar && xbar <= x);
    assert(0 <= ybar && ybar <= y);
  
    // mutation edges
    if (xbar == 0 && ybar == 0 && x > 0)
    {
      // insert mutation in x copy
      addEdge(v, MUTATION, x, y, xbar + 1, ybar);
    }
    if (xbar == 0 && ybar == 0 && y > 0)
    {
      // insert mutation in y copy
      addEdge(v, MUTATION, x, y, xbar, ybar + 1);
    }
    // amplification edges
    if (0 < x && x < _maxCopyNumberX && xbar < x)
    {
      // amplify nonmutated x copy
      addEdge(v, AMPLIFICATION, x+1, y, xbar, ybar);
    }
    if (0 < x && x < _maxCopyNumberX && xbar > 0)
    {
      // amplify mutated x copy
      addEdge(v, AMPLIFICATION, x+1, y, xbar+1, ybar);
    }
    if (0 < y && y < _maxCopyNumberY && ybar < y)
    {
      // amplify nonmutated y copy
      addEdge(v, AMPLIFICATION, x, y+1, xbar, ybar);
    }
    if (0 < y && y < _maxCopyNumberY && ybar > 0)
    {
      // amplify mutated y copy
      addEdge(v, AMPLIFICATION, x, y+1, xbar, ybar+1);
    }
    // deletion edges
    if (x > xbar)
    {
      // delete nonmutated x copy
      addEdge(v, DELETION, x-1, y, xbar, ybar);
    }
    if (xbar > 0)
    {
      // delete mutated x copy
      addEdge(v, DELETION, x-1, y, xbar-1, ybar);
    }
    if (y > ybar)
    {
      // delete nonmutated y copy
      addEdge(v, DELETION, x, y-1, xbar, ybar);
    }
    if (ybar > 0)
    {
      // delete mutated y copy
      addEdge(v, DELETION, x, y-1, xbar, ybar-1);
    }
  }
}
  
void GenotypeGraph::writeDOT(std::ostream& out) const
{
  IntNodeMap level(_G, 0);
  lemon::bfs(_G).distMap(level).run(_root);
  
  out << "digraph G {" << std::endl;
  
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    out << "\t" << _G.id(v) << " [label=\"("
        << _x[v] << "," << _y[v] << ","
        << _xbar[v] << "," << _ybar[v] << ")\\n" << level[v] << "\"]" << std::endl;
  }
  
  for (ArcIt a(_G); a != lemon::INVALID; ++a)
  {
    Node v = _G.source(a);
    Node w = _G.target(a);
    
    out << "\t" << _G.id(v) << " -> " << _G.id(w) << " [color=";
    switch (_type[a])
    {
      case MUTATION:
        out << "black";
        break;
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
  
void GenotypeGraph::init(SubDigraph& subG,
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
  
int GenotypeGraph::enumerate(const IntPairSet& L)
{
  BoolNodeMap filterNodesT(_G, false);
  BoolArcMap filterArcsT(_G, false);
  SubDigraph T(_G, filterNodesT, filterArcsT);
  
  BoolNodeMap filterNodesG(_G, true);
  BoolArcMap filterArcsG(_G, true);
  SubDigraph subG(_G, filterNodesG, filterArcsG);
  
  ArcList F;
  init(subG, T, F);
  
  Arc a = lemon::INVALID;
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

  grow(L, LL, subG, T, F, V, a, verticesCopyStates, leavesCopyStates);
  
  return _result.size();
}
  
bool GenotypeGraph::isValid(const IntPairSet& L,
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
  
  if (diff.size() > 0)
  {
    return false;
  }
  
  // check #2: all leaves are in L
  std::set_difference(leavesCopyStates.begin(), leavesCopyStates.end(),
                      L.begin(), L.end(),
                      std::inserter(diff, diff.begin()));
  if (diff.size() > 0)
  {
    return false;
  }
    
  return true;
}
  
bool GenotypeGraph::isValid(const SubDigraph& T) const
{
  if (!T.status(_root))
    return false;
  
  // at most one mutation edge
  int mutCount = 0;
  Arc mutationEdge = lemon::INVALID;
  for (SubArcIt a_cidj(T); a_cidj != lemon::INVALID; ++a_cidj)
  {
    if (_type[a_cidj] == MUTATION)
    {
      ++mutCount;
      mutationEdge = a_cidj;
    }
  }
  
  if (mutCount > 1)
    return false;
  
  int mut_x = _x[T.source(mutationEdge)];
  int mut_y = _y[T.source(mutationEdge)];
  
  // inf sites on copy states
  // this boils down to checking that |V_(x,y)| = 1
  IntMatrix V(_maxCopyNumberX + 1, IntVector(_maxCopyNumberY + 1, 0));
  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    int x = _x[v];
    int y = _y[v];
    
    if (x == mut_x && y == mut_y)
      continue;
    
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
  
bool GenotypeGraph::isValid(const SubDigraph& T,
                            Arc a)
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
  
void GenotypeGraph::writeDOT(const SubDigraph& T, std::ostream& out) const
{
  out << "digraph S {" << std::endl;

  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    out << "\t" << T.id(v) << " [label=\"("
        << _x[v] << "," << _y[v] << ","
        << _xbar[v] << "," << _ybar[v] << ")\"]" << std::endl;
  }
  
  for (SubArcIt a(T); a != lemon::INVALID; ++a)
  {
    out << "\t" << T.id(T.source(a)) << " -> " << T.id(T.target(a)) << std::endl;
  }
  
  out << "}" << std::endl;
}
  
void GenotypeGraph::writeDOT(const SubDigraph& T, const ArcList& F, std::ostream& out) const
{
  out << "digraph S {" << std::endl;
  
  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    out << "\t" << T.id(v) << " [label=\"("
        << _x[v] << "," << _y[v] << ","
        << _xbar[v] << "," << _ybar[v] << ")\"]" << std::endl;
  }
  
  for (ArcListIt it = F.begin(); it != F.end(); ++it)
  {
    Arc st = *it;
    Node t = _G.target(st);
    
    out << "\t" << T.id(t) << " [label=\"("
        << _x[t] << "," << _y[t] << ","
        << _xbar[t] << "," << _ybar[t] << ")\",color=red]" << std::endl;
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
  
bool GenotypeGraph::finalize(const IntPairSet& L,
                             const BoolMatrix& LL,
                             const SubDigraph& T,
                             const IntMatrix& V,
                             const Arc& mutationEdge)
{
  // convert to triples
  GenotypeNodeMap toGenotype(_G);

  for (SubNodeIt v(T); v != lemon::INVALID; ++v)
  {
    int x_v = _x[v];
    int y_v = _y[v];
    int xbar_v = _xbar[v];
    int ybar_v = _ybar[v];

    toGenotype[v]._x = x_v;
    toGenotype[v]._y = y_v;
    toGenotype[v]._xbar = xbar_v;
    toGenotype[v]._ybar = ybar_v;
  }

  // mutation edge must be in LL (mutation is introduced in an observed copy state)
  if (mutationEdge == lemon::INVALID ||
      !LL[_x[T.source(mutationEdge)]][_y[T.source(mutationEdge)]])
  {
    return false;
  }
  
  GenotypeTree::GenotypeEdgeSet S;
  
  // shortcircuit inner nodes not in LL
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
        
        if (LL[x_u][y_u] || (u == _root))
        {
          S.insert(std::make_pair(toGenotype[u], toGenotype[v]));
          break;
        }
        else
        {
          a = SubInArcIt(T, u);
        }
      }
    }
  }
  
  bool ok = false;
  for (const auto st : S)
  {
    ok |= (st.second._xbar > 0 || st.second._ybar > 0);
  }
  
  if (ok)
  {
    _result.insert(S);
    return false;
  }
  
  return false;
}
  
void GenotypeGraph::grow(const IntPairSet& L,
                         const BoolMatrix& LL,
                         SubDigraph& G,
                         SubDigraph& T,
                         ArcList& F,
                         IntMatrix& V,
                         Arc& mutationEdge,
                         IntPairSet& verticesCopyStates,
                         IntPairSet& leavesCopyStates)
{
  if (isValid(L, LL, T, V, verticesCopyStates, leavesCopyStates) && mutationEdge != lemon::INVALID)
  {
    // report
    finalize(L, LL, T, V, mutationEdge);
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
      
      // mutation edge => V[x_v][y_v] == 0
      assert(mutationEdge == lemon::INVALID || V[x_v][y_v] == 0);
      assert(V[x_v][y_v] <= 1);

      // update V
      ++V[x_v][y_v];
      if (_type[uv] == MUTATION)
      {
        mutationEdge = uv;
      }
      else
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
        // condition 3: remove other mutation edges if picked edge (u,v) was a mutation edge
        if (v == t || (v != s && x_v == x_t && y_v == y_t) || (mutationEdge == uv && V[x_t][y_t] == 1))
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
        
        if ((mutationEdge == lemon::INVALID && x_v == x_w && y_v == y_w)
            || V[x_w][y_w] == 0)
        {
          newF.push_back(vw);
        }
      }
      
//      static int idx = 0;
//      std::cout << ++idx << std::endl;
//      writeDOT(T, newF, std::cout);
      
      grow(L, LL, G, T, newF, V, mutationEdge, verticesCopyStates, leavesCopyStates);
      
      G.disable(uv);
      if (uv == mutationEdge)
      {
        mutationEdge = lemon::INVALID;
      }
      
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

void GenotypeGraph::partition(const int root_x,
                              const int root_y,
                              const GenotypeTree::GenotypeEdgeSet& S,
                              const GenotypeTree::Genotype& vertex,
                              bool& preMutationFlag,
                              GenotypeTree::Genotype& mutation,
                              GenotypeTree::GenotypeSet& preMutation,
                              GenotypeTree::GenotypeSet& postMutation)
{
  if (preMutationFlag)
  {
    preMutation.insert(vertex);
  }
  else
  {
    postMutation.insert(vertex);
  }
  
  for (const GenotypeTree::GenotypeEdge& st : S)
  {
    if (st.first == vertex)
    {
      if (preMutationFlag && (st.second._xbar == 1 || st.second._ybar == 1))
      {
        bool newPreMutationFlag = false;
        partition(root_x, root_y, S, st.second, newPreMutationFlag, mutation, preMutation, postMutation);
        mutation = st.second;
      }
      else
      {
        partition(root_x, root_y, S, st.second, preMutationFlag, mutation, preMutation, postMutation);
      }
    }
  }
}

void GenotypeGraph::partition(const int root_x,
                              const int root_y,
                              const GenotypeTree::GenotypeEdgeSet& S,
                              GenotypeTree::GenotypeSet& preMutation,
                              GenotypeTree::GenotypeSet& postMutation,
                              GenotypeTree::Genotype& mutation)
{
  GenotypeTree::Genotype root(root_x, root_y, 0, 0);
  bool preMutationFlag = true;
  mutation = GenotypeTree::Genotype();
  partition(root_x, root_y, S, root, preMutationFlag, mutation, preMutation, postMutation);
}

std::istream& operator>>(std::istream& in,
                         GenotypeGraph::Dictionary& dict)
{
  g_lineNumber = 0;
  std::string line;
  getline(in, line);
  std::stringstream ss(line);
  
  GenotypeGraph::Dictionary newDict;

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
      throw std::runtime_error(getLineNumber() + "Error: invalid number of genotype trees");
    }
    
    for (int j = 0; j < nrStateTrees; ++j)
    {
      getline(in, line);
      boost::split(s, line, boost::is_any_of(" "));
      
      GenotypeTree::GenotypeVector genotypeVector;
      for (const std::string& genotype_str : s)
      {
        int x = -1;
        int y = -1;
        int xbar = -1;
        int ybar = -1;
        if (sscanf(genotype_str.c_str(), "%d,%d,%d,%d", &x, &y, &xbar, &ybar) != 4 || x < 0 || y < 0 || xbar < 0 || ybar < 0)
        {
          throw std::runtime_error(getLineNumber() + "Error: invalid copy number '" + genotype_str + "'");
        }
        genotypeVector.push_back(GenotypeTree::Genotype(x, y, xbar, ybar));
      }
      
      if (genotypeVector.size() % 2 != 0)
      {
        throw std::runtime_error(getLineNumber() + "Error: odd number of genotypes '" + line + "'");
      }
      
      GenotypeTree::GenotypeEdgeSet genotypeTree;
      for (int jj = 0; jj < genotypeVector.size() / 2; ++jj)
      {
        genotypeTree.insert(GenotypeTree::GenotypeEdge(genotypeVector[2 * jj], genotypeVector[2 * jj + 1]));
      }
      
      if (newDict[L].count(genotypeTree) > 0)
      {
        throw std::runtime_error(getLineNumber() + "Error: duplicate state tree '" + line + "'");
      }
      newDict[L].insert(genotypeTree);
    }
  }
  
  dict = newDict;
  
  return in;
}

std::ostream& operator<<(std::ostream& out,
                         const GenotypeGraph::Dictionary& dict)
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
    
    out << kv.second.size() << " #state_trees" << std::endl;
    for (const GenotypeTree::GenotypeEdgeSet& S : kv.second)
    {
      first = true;
      for (const GenotypeTree::GenotypeEdge& edge : S)
      {
        if (first)
          first = false;
        else
          out << " ";
        
        out << edge.first._x << "," << edge.first._y << "," << edge.first._xbar << "," << edge.first._ybar << " ";
        out << edge.second._x << "," << edge.second._y << "," << edge.second._xbar << "," << edge.second._ybar;
      }
      out << std::endl;
    }
  }
  
  return out;
}

void GenotypeGraph::writeGenotypeTrees(std::ostream& out)
{
  out << _dict;
}

void GenotypeGraph::readGenotypeTrees(std::istream& in)
{
  Dictionary newDict;
  in >> newDict;
  
  _dict.insert(newDict.begin(), newDict.end());
}
