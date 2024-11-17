//
// Created by Mohammed El-Kebir on 8/14/23.
//

#ifndef CNATREE_H
#define CNATREE_H

#include "utils.h"
#include "basetree.h"
#include "genotypetree.h"

// Forward declaration
class CnaTree;

/// Vector of state trees
typedef std::vector<CnaTree> CnaTreeVector;

/// This class models a CNA tree
class CnaTree : public BaseTree {
public:
    /// This struct models a copy number state
    struct CnaGenotype {
        /// Default constructor
        CnaGenotype()
                : _x(-1), _y(-1) {
        }

        CnaGenotype(int x, int y)
                : _x(x), _y(y) {
        }

        bool operator==(const CnaGenotype &other) const {
            return _x == other._x && _y == other._y;
        }

        bool operator!=(const CnaGenotype &other) const {
            return !(_x == other._x && _y == other._y);
        }


        /// Number of maternal copies
        int _x;
        /// Number of paternal copies
        int _y;
    };

    typedef std::set<CnaGenotype> CnaStateSet;
    typedef CnaStateSet::const_iterator CnaStateSetIt;
    typedef std::vector<CnaGenotype> CnaStateVector;

    typedef std::pair<CnaGenotype, CnaGenotype> CnaEdge;
    typedef std::set<CnaEdge> CnaEdgeSet;
    typedef CnaEdgeSet::const_iterator CnaEdgeSetIt;

    typedef std::set<CnaEdgeSet> CnaEdgeSetSet;
    typedef CnaEdgeSetSet::const_iterator CnaEdgeSetSetIt;
    typedef CnaEdgeSetSet::iterator CnaEdgeSetSetNonConstIt;

public:
    /// Default constructor
    CnaTree();

    /// Constructor
    ///
    /// \param edges Edge set
    CnaTree(const CnaEdgeSet &edges);

    /// Constructor
    ///
    /// @param k Number of states
    CnaTree(int k);

    /// Copy constructor
    ///
    /// @param other Other CNA tree
    CnaTree(const CnaTree &other);

    /// Assignment operator
    ///
    /// @param other Other CNA tree
    CnaTree &operator=(const CnaTree &other);

    /// Constructor
    ///
    /// @param pi Vector of parental states
    CnaTree(const IntVector &pi);

    /// Set copy number state
    ///
    /// @param i State
    /// @param x Number of maternal copies
    /// @param y Number of paternal copies
    void setCnState(int i, int x, int y) {
        assert(0 <= i && i < _k);
        Node v_i = _stateToNode[i];
        _cnState[v_i]._x = x;
        _cnState[v_i]._y = y;
    }

    /// Get copy number state
    ///
    /// @param i State
    /// @param x Number of maternal copies
    /// @param y Number of paternal copies
    void getCnState(int i, int &x, int &y) const {
        assert(0 <= i && i < _k);
        Node v_i = _stateToNode[i];
        x = _cnState[v_i]._x;
        y = _cnState[v_i]._y;
    }

    /// Return number of maternal copies of specified state
    ///
    /// @param i State
    int x(int i) const {
        Node v_i = node(i);
        return _cnState[v_i]._x;
    }

    /// Return number of paternal copies of specified state
    ///
    /// @param i State
    int y(int i) const {
        Node v_i = node(i);
        return _cnState[v_i]._y;
    }

    /// Returns whether CNA tree is truncal
    bool truncal() const {
        return _k > 1 && lemon::countOutArcs(_T, _root) == 1;
    }

    bool hasLoss() const {
        return hasLoss_helper(_root);
    }

    bool hasLoss_helper(Node source) const {
        for (OutArcIt a(_T, source); a != lemon::INVALID; ++a) {
            Node target = _T.target(a);
            int x_source = _cnState[source]._x;
            int y_source = _cnState[source]._y;
            int x_target = _cnState[target]._x;
            int y_target = _cnState[target]._y;
            if (x_target < x_source | y_target < y_source) {
                return true;
            } else {
                bool childLoss = hasLoss_helper(target);
                if (childLoss) {
                    return true;
                }
            }
        }
        return false; //if we have arrived this far with no loss
    }

    void enumerateGenotypeTrees(GenotypeTree::GenotypeEdgeSetSet &result) const;

    void splitEnumerate(Node u, GenotypeTree::Genotype genotype_u,
                        NodeVector children,
                        const Node mutationNode, const bool mut_x,
                        bool preMutation,
                        GenotypeTree::GenotypeEdgeSet &currentTree,
                        GenotypeTree::GenotypeEdgeSetSet &result) const;

private:
    typedef Digraph::NodeMap <CnaGenotype> CopyNumberStateNodeMap;

    void fixPreMutation(Node u, Node mutationNode, GenotypeTree::GenotypeEdgeSet &tree) const;

    static bool next(BoolVector &boolVector);

    /// Copy number state node map
    CopyNumberStateNodeMap _cnState;

    friend std::ostream &operator<<(std::ostream &out, const CnaTree &S);

    friend std::istream &operator>>(std::istream &in, CnaTree &S);

    friend std::ostream &operator<<(std::ostream &out, const CnaTreeVector &vecS);

    friend std::istream &operator>>(std::istream &in, CnaTreeVector &vecS);
};

bool operator<(const CnaTree::CnaGenotype &lhs, const CnaTree::CnaGenotype &rhs);

#endif //CNATREE_H
