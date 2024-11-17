/*
 *  genotypetree.cpp
 *
 *   Created on: 28-sep-2015
 *       Author: M. El-Kebir
 */

#include "genotypetree.h"

#include <lemon/connectivity.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <iomanip>
#include <boost/random/gamma_distribution.hpp>

GenotypeTree::GenotypeTree(const GenotypeTree::GenotypeEdgeSet& edges)
    : BaseTree(), _cnState(_T)
{
  std::map<GenotypeTree::Genotype, int> nodes;
  _k = edges.size() + 1;

  // 1. identify vertices
  for (const GenotypeTree::GenotypeEdge & edge: edges)
  {
    nodes[edge.second] = nodes.size() + 1;
  }

  // 2. identify root
  GenotypeTree::Genotype root(-1, -1, -1, -1);
  for (const GenotypeTree::GenotypeEdge& edge: edges)
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
    Node v = _T.addNode();
    _nodeToState[v] = node.second;
    _stateToNode[node.second] = v;
    _cnState[v] = node.first;
    snprintf(buf, 1024, "(%d,%d,%d,%d)", _cnState[v]._x, _cnState[v]._y, _cnState[v]._xbar, _cnState[v]._ybar);
    _label[v] = buf;
  }

  _root = _stateToNode[0];

  // 4. add edges
  for (const GenotypeTree::GenotypeEdge& edge: edges)
  {
    _T.addArc(_stateToNode[nodes[edge.first]], _stateToNode[nodes[edge.second]]);
  }

  initD(_root);
  assert(lemon::dag(_T));
}


GenotypeTree::GenotypeTree()
  : BaseTree()
  , _cnState(_T)
{
}

GenotypeTree::GenotypeTree(int k)
  : BaseTree(k)
  , _cnState(_T)
{
  IntVector pi(_k, -1);
  for (int i = 1; i < _k; ++i)
  {
    pi[i] = i - 1;
  }
  
  init(pi);
}

//void GenotypeTree::writeSummary(const ReadMatrix& R,
//                                const GenotypeTreeVector& sol,
//                                std::ostream& out)
//{
//  const int n = R.getNrCharacters();
//  const int m = R.getNrSamples();
//
//  out << "SNV_index\t"
//      << "SNV_label\t"
//      << "sample_index\t"
//      << "sample_label\t"
//      << "ref\t"
//      << "alt\t"
//      << "VAF\t"
//      << "sanger_n_mut\t"
//      << "sanger_CCF\t"
//      << "decifer_CCF\t"
//      << "decifer_DCF\t"
//      << "decifer_cluster_CCF\t"
//      << "decifer_cluster_DCF\t"
////      << "decifer_cluster_index\t"
////      << "decifer_state_tree_index\t"
//      << "nr_nonzero_z"
//      << std::endl;
//
//  for (int i = 0; i < n; ++i)
//  {
//    for (int p = 0; p < m; ++p)
//    {
//      int ref = R.getRef(p, i);
//      int var = R.getVar(p, i);
//      double vaf = R.getVAF(p, i);
//      double purity = 1;
//      auto sanger = ReadMatrix::computeSangerCCF(ref, var,
//                                                                purity, R.getCopyNumberStates(p, i));
//      int sanger_n_mut = sanger.first;
//      double sangerCCF = sanger.second;
//
//      IntSet zSet;
//      IntPairSet mutTumorSet;
//      IntPairSet nonMutTumorSet;
//      const GenotypeTree& T_i = sol[i];
//      for (int j = 0; j != T_i.numVertices(); ++j)
//      {
//        if (T_i.isPresent(j))
//        {
//          int x_j = T_i.x(j);
//          int y_j = T_i.y(j);
//          int z_j = std::max(T_i.xbar(j), T_i.ybar(j));
//          if (z_j > 0)
//          {
//            zSet.insert(z_j);
//          }
//          if (x_j != 1 || y_j != 1)
//          {
//            if (z_j > 0)
//            {
//              mutTumorSet.insert(IntPair(x_j, y_j));
//            }
//            else
//            {
//              nonMutTumorSet.insert(IntPair(x_j, y_j));
//            }
//          }
//        }
//      }
//
//      double deciferCCF = T_i.maxLikelihoodCCF(p, vaf);
//      double deciferDCF = T_i.maxLikelihoodDCF(p, vaf);
//      double deciferClusterDCF = T_i.dcf(p);
//      double deciferClusterCCF = T_i.cf(p);
//
//      out << i << "\t" << R.indexToCharacter(i) << "\t"
//          << p << "\t" << R.indexToSample(p) << "\t"
//          << ref << "\t" << var << "\t"
//          << R.getVAF(p, i) << "\t"
//          << sanger_n_mut << "\t" << sangerCCF << "\t"
//          << deciferCCF << "\t" << deciferDCF << "\t"
//          << deciferClusterCCF << "\t" << deciferClusterDCF << "\t"
//          << zSet.size()
//          << std::endl;
//    }
//  }
//}

Node GenotypeTree::computeMaxLikelihoodStateProportions(DoubleNodeMap& s,
                                                        int p,
                                                        double vaf) const
{
  // 1. identify mutation vertex
  int star = -1;
  int x_star = -1;
  int y_star = -1;
  double mu_star = 0;
  Node v_star = lemon::INVALID;
  
  for (NodeIt v_j(_T); v_j != lemon::INVALID; ++v_j)
  {
    const int j = _nodeToState[v_j];
    if (isPresent(j))
    {
      int i = parent(j);
      if (i == -1 || !isPresent(i)) continue;
      if ((xbar(i) == 0 && xbar(j) == 1) || (ybar(i) == 0 && ybar(j) == 1))
      {
        star = j;
        x_star = x(i);
        y_star = y(i);
        v_star = v_j;
        mu_star = _s[_stateToNode[i]][p] + _s[v_j][p];
        break;
      }
    }
  }
  
  assert(star != -1);
  
  double C = mu_star * (x_star + y_star);
  double D = 0;
  for (NodeIt v_j(_T); v_j != lemon::INVALID; ++v_j)
  {
    const int j = _nodeToState[v_j];
    if (isPresent(j) && (x(j) != x_star || y(j) != y_star))
    {
      C += _s[v_j][p] * (x(j) + y(j));
      D += _s[v_j][p] * (xbar(j) + ybar(j));
    }
  }
  
  for (NodeIt v_j(_T); v_j != lemon::INVALID; ++v_j)
  {
    const int j = _nodeToState[v_j];
    if (x_star == x(j) && y_star == y(j))
    {
      if (j == star)
      {
        s[v_j] = std::min(mu_star, vaf * C - D);
      }
      else
      {
        s[v_j] = std::max(0., mu_star - vaf * C + D);
      }
    }
    else
    {
      s[v_j] = _s[v_j][p];
    }
  }
  
  return v_star;
}

double GenotypeTree::maxLikelihoodCCF(int p,
                                      double vaf) const
{
  DoubleNodeMap s(_T, 0.);
  computeMaxLikelihoodStateProportions(s, p, vaf);
  
  double ccf = 0;
  for (NodeIt v_j(_T); v_j != lemon::INVALID; ++v_j)
  {
    const int j = _nodeToState[v_j];
    if (isPresent(j) && (xbar(j) > 0 || ybar(j) > 0))
    {
      ccf += s[v_j];
    }
  }
  
  return ccf;
}

double GenotypeTree::maxLikelihoodDCF(int p,
                                      double vaf) const
{
  DoubleNodeMap s(_T, 0.);
  Node v_star = computeMaxLikelihoodStateProportions(s, p, vaf);
  int star = _nodeToState[v_star];
  
  double dcf = 0;
  for (NodeIt v_j(_T); v_j != lemon::INVALID; ++v_j)
  {
    const int j = _nodeToState[v_j];
//    if (isPresent(j))
//      std::cout << "s_{" << x(j) << "," << y(j) << "," << z(j) << "} = " << s[v_j] << " " << _cnState[v_j]._s[p] << std::endl;
    if (isPresent(j) && isDescendant(j, star))
    {
      dcf += s[v_j];
    }
  }
  
  return dcf;
}
  
GenotypeTree::GenotypeTree(const IntVector& pi)
  : BaseTree(pi)
  , _cnState(_T)
{
  init(pi);
}
  
GenotypeTree::GenotypeTree(const GenotypeTree& other)
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
  
GenotypeTree& GenotypeTree::operator=(const GenotypeTree& other)
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

void GenotypeTree::writeDOT(const int m,
                            std::ostream& out) const
{
  out << "digraph S {" << std::endl;
  out << "\tlabel=\"VAF =";
  for (int p = 0; p < m; ++p)
  {
    out << std::setprecision(2) << " " << vaf(p);
  }
  out << "\\nCF =";
  for (int p = 0; p < m; ++p)
  {
    out << std::setprecision(2) << " " << cf(p);
  }
  out << "\\nDCF =";
  for (int p = 0; p < m; ++p)
  {
    out << std::setprecision(2) << " " << dcf(p);
  }
  out << "\"" << std::endl;
  
  for (NodeIt v_i(_T); v_i != lemon::INVALID; ++v_i)
  {
    int i = _nodeToState[v_i];
    if (!isPresent(i)) continue;

    const auto& cnState_v_i = _cnState[v_i];
    
    out << "\t" << i << " [label=\"("
        << cnState_v_i._x << ","
        << cnState_v_i._y << ","
        << cnState_v_i._xbar << ","
        << cnState_v_i._ybar << ")";
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

void GenotypeTree::sampleMixtureProportions(double f, double purity)
{
  const int mutState = mutationState();
  const IntSet& D_mutState = D(mutState);
  boost::random::gamma_distribution<> gamma_dist(1, 1);
  
  NodeSet postMutNodes;
  NodeSet preMutNodes;
  double sum_pre = 0;
  double sum_post = 0;
  for (NodeIt v_i(_T); v_i != lemon::INVALID; ++v_i)
  {
    const int i = _nodeToState[v_i];
    if (!isPresent(i)) continue;

    _s[v_i].push_back(gamma_dist(g_rng));
    if (D_mutState.count(i) == 1) {
      postMutNodes.insert(v_i);
      sum_post += _s[v_i].back();
    }
    else
    {
      preMutNodes.insert(v_i);
      sum_pre += _s[v_i].back();
    }
  }
  
  // zero out nodes with small weights
  for (Node v_i : postMutNodes)
  {
    if ((_s[v_i].back() / sum_post) < 0.05)
    {
      sum_post -= _s[v_i].back();
      _s[v_i].back() = 0;
    }
  }
  for (Node v_i : preMutNodes)
  {
    if ((_s[v_i].back() / sum_pre) < 0.05)
    {
      sum_pre -= _s[v_i].back();
      _s[v_i].back() = 0;
    }
  }
  
  for (Node v_i : postMutNodes)
  {
    _s[v_i].back() /= sum_post;
    _s[v_i].back() *= f;
    if (!g_tol.nonZero(_s[v_i].back()))
    {
      _s[v_i].back() = 0;
    }
  }
  
  for (Node v_i : preMutNodes)
  {
    _s[v_i].back() /= sum_pre;
    _s[v_i].back() *= (1 - f - (1 - purity));
    if (!g_tol.nonZero(_s[v_i].back()))
    {
      _s[v_i].back() = 0;
    }
  }
  
  _s[_root].back() += 1 - purity;
  if (!g_tol.nonZero(_s[_root].back()))
  {
    _s[_root].back() = 0;
  }
}

void GenotypeTree::writeClusterDOT(double gamma,
                                   int t,
                                   std::ostream& out) const
{
  const int m = getNrSamples();
  
  out << "\tsubgraph cluster_" << t << " {" << std::endl;
  out << "\t\tlabel=\"gamma = " << gamma;
  
  out << "\\nVAF = ";
  for (int p = 0; p < m; ++p)
  {
    out << std::setprecision(2) << " " << vaf(p);
  }
  out << "\\nCF =";
  for (int p = 0; p < m; ++p)
  {
    out << std::setprecision(2) << " " << cf(p);
  }
  out << "\\nDCF =";
  for (int p = 0; p < m; ++p)
  {
    out << std::setprecision(2) << " " << dcf(p);
  }
  out << "\"" << std::endl;
  
  for (NodeIt v_i(_T); v_i != lemon::INVALID; ++v_i)
  {
    int i = _nodeToState[v_i];
    if (!isPresent(i)) continue;

    const auto& cnState_v_i = _cnState[v_i];
    
    out << "\t\t" << i + t * 100 << " [label=\"("
        << cnState_v_i._x << ","
        << cnState_v_i._y << ","
        << cnState_v_i._xbar << ","
        << cnState_v_i._ybar << ")";
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
    
    out << "\t\t" << i + t * 100 << " -> " << j + t * 100 << std::endl;
  }
  
  out << "\t}" << std::endl;
}

std::ostream& operator<<(std::ostream& out, const GenotypeTree& T)
{
  // output pi
  out << -1;
  for (int i = 1; i < T.k(); ++i)
  {
    out << " " << T.parent(i);
  }
  out << std::endl;
  
  // output vertices
  for (int i = 0; i < T.k(); ++i)
  {
    Node v_i = T.node(i);
    out << T._cnState[v_i]._x
        << " " << T._cnState[v_i]._y
        << " " << T._cnState[v_i]._xbar
        << " " << T._cnState[v_i]._ybar;
    for (double s : T._s[v_i])
    {
      out << " " << s;
    }
    out << std::endl;
  }
  
  // output labels
  out << T.label(0);
  for (int i = 1; i < T.k(); ++i)
  {
    out << " " << T.label(i);
  }
  out << std::endl;

  return out;
}
  
std::istream& operator>>(std::istream& in, GenotypeTree& T)
{
  std::string line;
  getline(in, line);
  std::stringstream ss(line);
  
  StringVector s;
  boost::split(s, line, boost::is_any_of(" \t"));
  
  IntVector pi(s.size(), -1);
  for (int i = 0; i < s.size(); ++i)
  {
    ss >> pi[i];
  }

  T = GenotypeTree(pi);
  
  for (int i = 0; i < T.k(); ++i)
  {
    getline(in, line);
    ss.clear();
    ss.str(line);
    
    Node v_i = T.node(i);
    GenotypeTree::Genotype& cnState_i = T._cnState[v_i];
    
    boost::split(s, line, boost::is_any_of(" "));
    
    cnState_i._x = boost::lexical_cast<int>(s[0]);
    cnState_i._y = boost::lexical_cast<int>(s[1]);
    cnState_i._xbar = boost::lexical_cast<int>(s[2]);
    cnState_i._ybar = boost::lexical_cast<int>(s[3]);

    for (int j = 4; j < s.size(); ++j)
    {
      T._s[v_i].push_back(boost::lexical_cast<double>(s[j]));
    }
  }
  
  getline(in, line);
  ss.clear();
  ss.str(line);
  for (int i = 0; i < T.k(); ++i)
  {
    ss >> T._label[T._stateToNode[i]];
  }
  
  return in;
}

std::ostream& operator<<(std::ostream& out, const GenotypeTreeVector& vecT)
{
  out << vecT.size() << " #genotype trees" << std::endl;
  for (const GenotypeTree& T : vecT)
  {
    out << T;
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, GenotypeTreeVector& vecT)
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
    GenotypeTree T;
    in >> T;
    vecT.push_back(T);
  }
  
  return in;
}

bool operator<(const GenotypeTree::Genotype& lhs, const GenotypeTree::Genotype& rhs)
{
  return lhs._x < rhs._x
    || (lhs._x == rhs._x && lhs._y < rhs._y)
    || (lhs._x == rhs._x && lhs._y == rhs._y && lhs._xbar < rhs._xbar)
    || (lhs._x == rhs._x && lhs._y == rhs._y && lhs._xbar == rhs._xbar && lhs._ybar < rhs._ybar);
}
