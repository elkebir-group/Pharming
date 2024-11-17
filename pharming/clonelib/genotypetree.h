/*
 *  genotypetree.h
 *
 *   Created on: 28-sep-2015
 *       Author: M. El-Kebir
 */

#ifndef GENOTYPETREE_H
#define GENOTYPETREE_H

#include "basetree.h"
#include "utils.h"
#include "readmatrix.h"

// Forward declaration
class GenotypeTree;

/// Vector of state trees
typedef std::vector<GenotypeTree> GenotypeTreeVector;

/// This class models a genotype tree
class GenotypeTree : public BaseTree
{
public:
  /// This struct models an SNV-specific copy number state
  struct Genotype
  {
    /// Default constructor
    Genotype()
      : _x(-1)
      , _y(-1)
      , _xbar(-1)
      , _ybar(-1)
    {
    }

    /// Constructor
    Genotype(int x, int y, int xbar, int ybar)
      : _x(x)
      , _y(y)
      , _xbar(xbar)
      , _ybar(ybar)
    {
    }

    bool operator==(const Genotype& other) const
    {
      return _x == other._x && _y == other._y && _xbar == other._xbar && _ybar == other._ybar;
    }

    bool operator!=(const Genotype& other) const
    {
      return !(_x == other._x && _y == other._y && _xbar == other._xbar && _ybar == other._ybar);
    }
    double vaf() const
    {
      int totalCopies = _x + _y;
      if (totalCopies == 0) {
          return 0.0; // Handle division by zero case
      }
      return (_xbar + _ybar) / (totalCopies * 1.0);
    }

      double mutCopies() const
      {
          return (_xbar + _ybar);
      }
    /// Number of maternal copies
    int _x;
    /// Number of paternal copies
    int _y;
    /// Number of mutated maternal copies
    int _xbar;
    /// Number of mutated paternal copies
    int _ybar;
  };

public:
  typedef std::set<Genotype> GenotypeSet;
  typedef GenotypeSet::const_iterator GenotypeSetIt;
  typedef std::vector<Genotype> GenotypeVector;

  typedef std::pair<Genotype, Genotype> GenotypeEdge;
  typedef std::set<GenotypeEdge> GenotypeEdgeSet;
  typedef GenotypeEdgeSet::const_iterator GenotypeEdgeSetIt;

  typedef std::set<GenotypeEdgeSet> GenotypeEdgeSetSet;
  typedef GenotypeEdgeSetSet::const_iterator  GenotypeEdgeSetSetIt;
  typedef GenotypeEdgeSetSet::iterator GenotypeEdgeSetSetNonConstIt;

public:
  /// Default constructor
  GenotypeTree();
  
  /// Constructor
  ///
  /// @param k Number of states
  GenotypeTree(int k);
  
  /// Copy constructor
  ///
  /// @param other Other genotype tree
  GenotypeTree(const GenotypeTree& other);
  
  /// Assignment operator
  ///
  /// @param other Other genotype tree
  GenotypeTree& operator=(const GenotypeTree& other);
  
  /// Constructor
  ///
  /// @param pi Vector of parental states
  GenotypeTree(const IntVector& pi);

  /// Constructor
  ///
  /// \param edges Edge set
  GenotypeTree(const GenotypeEdgeSet& edges);
  
  /// Write summary
  ///
  /// @param R Read matrix
  /// @param sol State tree vector
  /// @param out Output stream
//  static void writeSummary(const ReadMatrix& R,
//                           const GenotypeTreeVector& sol,
//                           std::ostream& out);

  /// Reset mixture proportions using specified dcf
  ///
  /// \param p Sample
  /// \param dcf DCF value
  void resetMixtureProportions(int p, double dcf)
  {
    // Identify the mutation node
    const Node v_i = mutationNode(_root);
    const int i = _nodeToState[v_i];

    // Identify the parent of the mutation node
    const Node v_pi_i = _stateToNode[parent(i)];

    // Compute mixture proportion of CN state of mutation node
    const double marginal_sum = _s[v_i][p] + _s[v_pi_i][p];
    
    double sum = 0;
    for (int j = 0; j < _k; ++j)
    {
      if (i == j) continue;
      Node v_j = _stateToNode[j];
      if (v_j == lemon::INVALID)
        continue;
      
      if (isAncestor(i, j))
      {
        sum += _s[v_j][p];
      }
    }
    
    _s[v_i][p] = dcf - sum;
    _s[v_pi_i][p] = marginal_sum - _s[v_i][p];
  }
  
  /// Set copy number state
  ///
  /// @param i State
  /// @param x Number of maternal copies
  /// @param y Number of paternal copies
  /// @param xbar Number of mutated maternal copies
  /// @param ybar Number of mutated paternal copies
  void setCnState(int i, int x, int y, int xbar, int ybar)
  {
    assert(0 <= i && i < _k);
    Node v_i = _stateToNode[i];
    _cnState[v_i]._x = x;
    _cnState[v_i]._y = y;
    _cnState[v_i]._xbar = xbar;
    _cnState[v_i]._ybar = ybar;
  }
  
  /// Get copy number state
  ///
  /// @param i State
  /// @param x Number of maternal copies
  /// @param y Number of paternal copies
  /// @param xbar Number of mutated maternal copies
  /// @param ybar Number of mutated paternal copies
  void getCnState(int i, int& x, int& y, int& xbar, int& ybar) const
  {
    assert(0 <= i && i < _k);
    Node v_i = _stateToNode[i];
    x = _cnState[v_i]._x;
    y = _cnState[v_i]._y;
    xbar = _cnState[v_i]._xbar;
    ybar = _cnState[v_i]._ybar;
  }
  
  void writeClusterDOT(double gamma,
                       int t,
                       std::ostream& out) const;
  
  /// Write state tree in DOT format
  ///
  /// @param m Number of samples
  /// @param out Output stream
  void writeDOT(const int m,
                std::ostream& out) const;
  
  /// Return number of maternal copies of specified state
  ///
  /// @param i State
  int x(int i) const
  {
    Node v_i = node(i);
    return _cnState[v_i]._x;
  }
  
  /// Return number of paternal copies of specified state
  ///
  /// @param i State
  int y(int i) const
  {
    Node v_i = node(i);
    return _cnState[v_i]._y;
  }

  /// Return number of mutated copies of specified state
  ///
  /// @param i State
  int z(int i) const
  {
    Node v_i = node(i);
    return std::max(_cnState[v_i]._xbar, _cnState[v_i]._ybar);
  }

  /// Return number of mutated maternal copies of specified state
  ///
  /// @param i State
  int xbar(int i) const
  {
    Node v_i = node(i);
    return _cnState[v_i]._xbar;
  }

  /// Return number of mutated paternal copies of specified state
  ///
  /// @param i State
  int ybar(int i) const
  {
    Node v_i = node(i);
    return _cnState[v_i]._ybar;
  }
  
  /// Return the unique state whose incoming edge is a mutation edge
  int mutationState() const
  {
    Node mut_node = mutationNode(_root);
    if (mut_node == lemon::INVALID)
    {
      return -1;
    }
    else
    {
      return _nodeToState[mut_node];
    }
  }
  
  /// Return variant allele frequency of the specified sample
  ///
  /// @param p Sample
  double vaf(int p) const
  {
    double numerator = 0;
    double denominator = 0;
    for (NodeIt v_i(_T); v_i != lemon::INVALID; ++v_i)
    {
      const int i = _nodeToState[v_i];
      const auto& cnState_v_i = _cnState[v_i];

      if (isPresent(i))
      {
        numerator += (cnState_v_i._xbar + cnState_v_i._ybar) * _s[v_i][p];
        denominator += (cnState_v_i._x + cnState_v_i._y) * _s[v_i][p];
      }
    }
    return numerator / denominator;
  }
  
  /// Return descendant cell fraction of the specified sample
  ///
  /// @param p Sample
  double dcf(int p) const
  {
    double res = 0;
    
    int i = mutationState();
    assert(i != -1);
    
    for (int j : D(i))
    {
      Node v_j = _stateToNode[j];
      res += _s[v_j][p];
    }
    
    return res;
  }
  
  /// Return cell fraction of the specified sample
  ///
  /// @param p Sample
  double cf(int p) const
  {
    double res = 0;
    
    for (NodeIt v_i(_T); v_i != lemon::INVALID; ++v_i)
    {
      const int i = _nodeToState[v_i];
      const auto& cnState_v_i = _cnState[v_i];
      if (isPresent(i) && (cnState_v_i._xbar > 0 || cnState_v_i._ybar > 0))
      {
        res += _s[v_i][p];
      }
    }

    return res;
  }
  
  /// Return maximum likelihood CCF given vaf
  ///
  /// @param p Sample
  /// @param vaf
  double maxLikelihoodCCF(int p,
                          double vaf) const;
  
  /// Return maximum likelihood DCF given vaf
  ///
  /// @param p Sample
  /// @param vaf
  double maxLikelihoodDCF(int p,
                          double vaf) const;
  
  /// Generate mixture proportions for specified sample and dcf
  ///
  /// @param f DCF
  /// @param purity Purity
  void sampleMixtureProportions(double f, double purity);
  
  /// Check whether mixture proportions are ok
  ///
  /// @param m Number of samples
  bool okMixtureProportions(int m) const
  {
    for (NodeIt v_i(_T); v_i != lemon::INVALID; ++v_i)
    {
      const int i = _nodeToState[v_i];
      const Genotype& cnState_i = _cnState[v_i];
      if (isPresent(i))
      {
        if (_s[v_i].size() != m)
        {
          return false;
        }
        
        bool ok = false;
        for (int p = 0; p < m; ++p)
        {
          ok |= g_tol.nonZero(_s[v_i][p]);
        }
        if (!ok)
        {
          return false;
        }
      }
    }
    
    return true;
  }

private:
  typedef Digraph::NodeMap<Genotype> CopyNumberStateNodeMap;
  
  Node computeMaxLikelihoodStateProportions(DoubleNodeMap& s,
                                            int p,
                                            double vaf) const;
  
private:
  /// Copy number state node map
  CopyNumberStateNodeMap _cnState;

  /// Identify mutation node recursively from v_i
  ///
  /// @param v_i Node
  Node mutationNode(Node v_i) const
  {
    if (_cnState[v_i]._xbar > 0 || _cnState[v_i]._ybar > 0)
    {
      return v_i;
    }
    else
    {
      for (OutArcIt a_ij(_T, v_i); a_ij != lemon::INVALID; ++a_ij)
      {
        Node v_j = _T.target(a_ij);
        Node mut_node = mutationNode(v_j);
        if (mut_node != lemon::INVALID)
        {
          return mut_node;
        }
      }
      
      return lemon::INVALID;
    }
  }
  
  friend std::ostream& operator<<(std::ostream& out, const GenotypeTree& S);
  friend std::istream& operator>>(std::istream& in, GenotypeTree& S);
  friend std::ostream& operator<<(std::ostream& out, const GenotypeTreeVector& vecS);
  friend std::istream& operator>>(std::istream& in, GenotypeTreeVector& vecS);
};

/// Output state tree
///
/// @param out Output stream
/// @param T Genotype tree
std::ostream& operator<<(std::ostream& out, const GenotypeTree& T);

/// Input state tree
///
/// @param in Input stream
/// @param T Genotype tree
std::istream& operator>>(std::istream& in, GenotypeTree& T);

/// Output state tree vector
///
/// @param out Output stream
/// @param vecS Genotype tree vector
std::ostream& operator<<(std::ostream& out, const GenotypeTreeVector& vecS);

/// Input state tree vector
///
/// @param in Input stream
/// @param vecS Genotype tree vector
std::istream& operator>>(std::istream& in, GenotypeTreeVector& vecS);

bool operator<(const GenotypeTree::Genotype& lhs, const GenotypeTree::Genotype& rhs);

#endif // GENOTYPETREE_H
