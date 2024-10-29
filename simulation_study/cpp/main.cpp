#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <cassert>
#include <iterator>
#include <map>
//#include <unordered_set>

/*
 * Output:
 * - cell clustering TP, FP, TN, FN
 * - cell ancestral TP, FP, TN, FN
 * - cell incomparable TP, FP, TN, FN
 * - genotype_hamming
 * - n_inf_cell_clusters
 * - n_inf_mut_clusters
 */

//// Define a custom hash function for pairs of uint64_t integers
//struct PairHash {
//    size_t operator()(const std::pair<uint64_t, uint64_t>& pair) const {
//        // Combine the hash values of the two elements
//        return std::hash<uint64_t>()(pair.first) ^ (std::hash<uint64_t>()(pair.second) << 1);
//    }
//};
//
//// Define a custom equality predicate for pairs of uint64_t integers
//struct PairEqual {
//    bool operator()(const std::pair<uint64_t, uint64_t>& left, const std::pair<uint64_t, uint64_t>& right) const {
//        return left.first == right.first && left.second == right.second;
//    }
//};

// Define ClusteringPairs as an unordered_set of pairs of uint64_t integers
//using ClusteringPairs = std::unordered_set<std::pair<uint64_t, uint64_t>, PairHash, PairEqual>;
typedef std::set <uint64_t> IntSet;
typedef std::vector <IntSet> IntSetVector;

typedef std::vector <uint64_t> IntVector;
typedef std::unordered_map <std::string, IntVector> Clustering;
typedef std::vector <std::pair <uint64_t, Clustering>> ClusteringVector;

typedef std::pair <uint64_t, uint64_t> IntPair;
typedef std::set <IntPair> ClusteringPairs;
typedef std::vector <ClusteringPairs> ClusteringPairsVector;

typedef std::vector <std::string> StringVector;
typedef std::vector<std::set<uint64_t>> Tree;


struct geno_presence{
    uint64_t x;
    uint64_t y;
    int pres;

  


    bool operator!=(const geno_presence& other) const {
        return x != other.x || y != other.y || pres != other.pres;
    }

    bool operator==(const geno_presence& other) const {
        return x == other.x && y == other.y && pres == other.pres;
    }

    bool operator<(const geno_presence& other) const {

        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
      
        return pres < other.pres;
    }

    bool operator>(const geno_presence& other) const {

        if (x != other.x) return x > other.x;
        if (y != other.y) return y > other.y;
        return pres > other.pres;
    }


};

struct geno{
    uint64_t x;
    uint64_t y;
    uint64_t xbar;
    uint64_t ybar;

    int mut_copies() const {
        return std::max(xbar, ybar);
    }

     int presence() const {
        return xbar > 0 || ybar > 0;
    }


    bool operator!=(const geno& other) const {
        return x != other.x || y != other.y || xbar != other.xbar || ybar != other.ybar;
    }

    bool operator==(const geno& other) const {
        return x == other.x && y == other.y && mut_copies() ==  other.mut_copies();
    }

    bool operator<(const geno& other) const {

        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        if (xbar != other.xbar) return xbar < other.xbar;
        return ybar < other.ybar;
    }

    bool operator>(const geno& other) const {

        if (x != other.x) return x > other.x;
        if (y != other.y) return y > other.y;
        if (xbar != other.xbar) return xbar > other.xbar;
        return ybar > other.ybar;
    }


};

struct CNA{
    uint64_t x;
    uint64_t y;

    bool operator==(const CNA& other) const {
        return x == other.x && y == other.y;
    }

    bool operator!=(const CNA& other) const {
        return x != other.x || y != other.y;
    }


    bool operator<(const CNA& other) const {

        if (x != other.x) return x < other.x;
        return y < other.y;

    }

};

typedef std::map<int, std::vector<geno>> Genotypes;
typedef std::map<int, std::map<int,CNA>> CNAGeno;
typedef std::set<std::pair<geno, geno>> SNVTree;
typedef std::set<std::pair<geno_presence, geno_presence>> SNVTreePresence;
typedef std::set<std::pair<CNA,CNA>> CNATree;
//StringVector split(const std::string& input, const std::string& regex) {
//    // passing -1 as the submatch index parameter performs splitting
//    std::regex re(regex);
//    std::sregex_token_iterator
//            first{input.begin(), input.end(), re, -1},
//            last;
//    return {first, last};
//}

StringVector split (const std::string &s, char delim) {
    StringVector result;
    std::stringstream ss (s);
    std::string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

ClusteringPairs getIncomparableNodePairs(const Tree& T, const ClusteringPairs& ancPairs, const IntSet nodes) {
    // IntSet nodes;
    // for (uint64_t source = 0; source < T.size(); ++source) {
    //     nodes.insert(source);
    //     for (uint64_t target : T[source]) {
    //         nodes.insert(target);
    //     }
    // }

    ClusteringPairs incPairs;
    for (auto it1 = nodes.begin(); it1 != nodes.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != nodes.end(); ++it2) {
            if (ancPairs.count({*it1, *it2}) == 0 && ancPairs.count({*it2, *it1}) == 0) {
                incPairs.insert({*it1, *it2});
            }
        }
    }

    return incPairs;
}

ClusteringPairs getAncestralNodePairs(const Tree& T) {
    ClusteringPairs pairs;
    for (uint64_t source = 0; source < T.size(); ++source) {
        for (auto target : T[source]) {
            pairs.insert({source, target});
        }
    }
    // std::cout << pairs.size() << std::endl;
    bool change = true;
    while (change) {
        change = false;
        for (auto it1 = pairs.begin(); !change && it1 != pairs.end(); ++it1) {
            for (auto it2 = pairs.begin(); !change && it2 != pairs.end(); ++it2) {
                if (it1->second == it2->first && pairs.count({it1->first, it2->second}) == 0) {
                    change = true;
                    pairs.insert({it1->first, it2->second});
                    break;
                }
            }
        }
    }

    return pairs;
}

IntSetVector cellClusterGenotypes(const Tree& T) {
    IntSet nodes;
    for (uint64_t source = 0; source < T.size(); ++source) {
        nodes.insert(source);
        for (uint64_t target : T[source]) {
            nodes.insert(target);
        }
    }

    IntSetVector genotypes(nodes.size());

    auto ancPairs = getAncestralNodePairs(T);
    for (uint64_t target = 0; target < nodes.size(); ++target)
    {
        genotypes[target].insert(target);
        for (const auto& ancPair : ancPairs) {
            if (ancPair.second == target) {
                genotypes[target].insert(ancPair.first);
            }
        }
    }

    return genotypes;
}

std::pair<IntSet, Tree> parseTree(std::istream &in) {
    Tree adjList;

    std::string line;
    IntSet nodes;
    getline(in, line);
    std::stringstream ss(line);

    uint64_t nrEdges = 0;
    ss >> nrEdges;

    for (uint64_t idx = 0; idx < nrEdges; ++idx) {
        getline(in, line);

        StringVector s = split(line, ' ');

        assert(s.size() == 2);

        char *end;
        uint64_t source = std::strtoull(s[0].c_str(), &end, 10);
        uint64_t target = std::strtoull(s[1].c_str(), &end, 10);
        nodes.insert(source);
        nodes.insert(target);

        uint64_t n = adjList.size();
    
        if (source >= n) {
            for (uint64_t i = 0; i <= source - n; ++i) {
                adjList.push_back(std::set<uint64_t>());
            }
        }

        n = adjList.size();
        if (target >= n) {
            for (uint64_t i = 0; i <= target - n; ++i) {
                adjList.push_back(std::set<uint64_t>());
            }
        }

        adjList[source].insert(target);
       
    }

    return {nodes, adjList};
}

std::map<uint64_t,uint64_t> reverseMap(const IntSetVector clustMap){
    std::map<uint64_t,uint64_t> pairs;
    for(int q=0; q < clustMap.size(); q++){
        for(auto i: clustMap[q]){
            pairs.insert({i,q});
        }
    }

    return pairs;
}

std::map<uint64_t, IntPair> makeCellMap(const IntSetVector gt,const IntSetVector inf ){
    std::map<uint64_t, IntPair> cellMap;
    std::map<uint64_t, uint64_t>  gtPairs = reverseMap(gt);
    std::map<uint64_t, uint64_t>  infPairs = reverseMap(inf);
    for (const auto& pair : gtPairs) {

        cellMap.insert({pair.first, {pair.second,  infPairs[pair.first] }});

    }


    return cellMap;

}
void computeGenotypeSimilarity(const std::pair<Genotypes, CNAGeno > gt_genos, const std::pair<Genotypes, CNAGeno> inf_genos,
                               const std::map<uint64_t, IntPair> cellMap,
                               const std::set<uint64_t> snvs,
                               uint64_t m,
                               uint64_t &pres,  uint64_t &allspec, uint64_t &sim,uint64_t &cna, uint64_t &correct_sim, uint64_t &correct_presence )
{
    // Initialize accumulators to zero

    Genotypes gt = gt_genos.first;
    Genotypes inf =inf_genos.first;

    CNAGeno cngt = gt_genos.second;
    CNAGeno cninf = inf_genos.second;


    for(size_t i = 0; i < cellMap.size(); ++i){
        IntPair nodes = cellMap.at(i);

        for(auto j: snvs){
            if(snvs.count(j) > 0){
                geno g1 = gt.at(nodes.first)[j];
                geno g2 = inf.at(nodes.second)[j];

                // Compare mutation copies
                pres += (g1.mut_copies() > 0) == (g2.mut_copies() > 0);

                // Compare xbar and ybar values
                allspec += (g1.xbar == g2.xbar) && (g1.ybar == g2.ybar);

                // Compare mutation copies
                sim += (g1.mut_copies() == g2.mut_copies());

                correct_presence += (g1.x == g2.x) && (g1.y == g2.y) && (g1.presence() == g2.presence());
                correct_sim += (g1.x == g2.x) && (g1.y == g2.y) && (g1.mut_copies() == g2.mut_copies());
            }
        }
        for (const auto& pair : cngt.at(nodes.first)) {
            int seg = pair.first;
            CNA c1 = pair.second;
            if(cninf.at(nodes.second).count(seg) > 0){
                CNA c2 = cninf.at(nodes.second).at(seg);
                cna += (c1.x == c2.x) && (c1.y == c2.y);
            }

        }
    }
}

std::pair<Genotypes , CNAGeno> parseGenotypes(std::istream &in, const IntSet nodes, const int m) {
    Genotypes genotypes;
    CNAGeno cnas;


    std::string line;

    for(auto n: nodes){
        genotypes[n] = std::vector<geno>(m);
    }

    // skip the first line, header
    getline(in, line);

    while (in.good()) {
        getline(in, line);
        if (line.empty()) break;

        StringVector s = split(line, ',');

        geno g;
        CNA cn;
        char *end;
        uint64_t n = std::strtoull(s[0].c_str(), &end, 10);
        uint64_t snv = std::strtoull(s[1].c_str(), &end, 10);
        g.x = std::strtoull(s[2].c_str(), &end, 10);
        g.y = std::strtoull(s[3].c_str(), &end, 10);
        g.xbar = std::strtoull(s[4].c_str(), &end, 10);
        g.ybar = std::strtoull(s[5].c_str(), &end, 10);
        int seg = std::strtoull(s[6].c_str(), &end, 10);
        cn.x  = g.x;
        cn.y = g.y;


        genotypes[n][snv]= g;

        if(cnas.count(n) ==1){
            if(cnas[n].count(seg) < 1){
                cnas[n].insert({seg, cn});
            }
        }else{
            cnas[n].insert( {seg, cn});
        }

    }

    return {genotypes, cnas};
}

std::pair<uint64_t, IntSetVector> parseClustering(std::istream &in, int elementIdx, int clusterIdx) {
    uint64_t count = 0;
    IntSetVector res;
    std::string line;

    // skip the first line, header
    getline(in, line);

    while (in.good()) {
        getline(in, line);
        if (line.empty()) break;

        StringVector s = split(line, ',');

        assert(s.size() > elementIdx && s.size() > clusterIdx);

        char *end;
        uint64_t cluster = std::strtoull(s[clusterIdx].c_str(), &end, 10);
//        int indx = s[elementIdx].find('_');
        uint64_t element = std::strtoull(s[elementIdx].c_str() + s[elementIdx].find('_') + 1, &end, 10);

        uint64_t n = res.size();
        if (cluster >= n) {
            for (uint64_t i = 0; i <= cluster - n; ++i) {
                res.push_back(IntSet());
            }
        }

        res[cluster].insert(element);
    
        ++count;
    }

    return {count, res};
}

void performComparison(const uint64_t n,
                       const ClusteringPairs &GT, const ClusteringPairs &inferred,
                       uint64_t &TP, uint64_t &TN, uint64_t &FN, uint64_t &FP) {
    std::set <std::pair<uint64_t, uint64_t>> S;


    std::set_difference(GT.begin(), GT.end(),
                        inferred.begin(), inferred.end(),
                        std::inserter(S, S.begin()));
    FN = S.size();

    S.clear();

    std::set_difference(inferred.begin(), inferred.end(),
                        GT.begin(), GT.end(),
                        std::inserter(S, S.begin()));
    FP = S.size();

    S.clear();

    std::set_intersection(GT.begin(), GT.end(),
                          inferred.begin(), inferred.end(),
                          std::inserter(S, S.begin()));
    TP = S.size();

    TN = n - FN - FP - TP;
}

ClusteringPairs makePairs(const IntSetVector &clustering) {
    ClusteringPairs pairs;

    for (const auto &set: clustering) {
        for (auto it1 = set.begin(); it1 != set.end(); ++it1) {
            for (auto it2 = std::next(it1); it2 != set.end(); ++it2) {
                if(*it1 < *it2){
                    pairs.insert(std::make_pair(*it1, *it2));
                }else{
                    pairs.insert(std::make_pair(*it2, *it1));
                }

            }
        }
    }

    return pairs;
}

ClusteringPairs expandClustering(const ClusteringPairs& unexpPairs, const IntSetVector& clustering, const bool ancestral) {
   
    ClusteringPairs res;
    // std::cout << clustering.size() << std::endl;
    for (const auto& unexpPair : unexpPairs) {
        uint64_t sourceIdx = unexpPair.first;
        uint64_t targetIdx = unexpPair.second;
//         std::cout <<"Source size:" << clustering[sourceIdx].size() << std::endl;
         if(clustering.size() <= sourceIdx || clustering.size() <= targetIdx ) {
             continue;
         }
 
        if(!clustering[sourceIdx].empty() && !clustering[targetIdx].empty()){
                for (uint64_t item1 : clustering[sourceIdx]) {
                        for (uint64_t item2 : clustering[targetIdx]) {
//
                                if(ancestral or item1 < item2){
                                    res.insert({item1, item2});
                                }else{
                                    res.insert({item2, item1});
                                }

                    }
                 
                }
    
        }
    }
    // std::cout << res.size() << std::endl;
    return res;
       

}

uint64_t getNrNonEmptyClusters(const IntSetVector& clustering) {
    uint64_t count = 0;
    for (const auto& part : clustering) {
        if (part.size() > 0) ++count;
    }
    return count;
}

//double computeGenotypeSimilarity(uint64_t nrCells,
//                                 uint64_t nrMutations,
//                                 const IntSetVector& gtGenotypes,
//                                 const IntSetVector& infGenotypes,
//                                 const IntSetVector& gtCellClustering,
//                                 const IntSetVector& infCellClustering,
//                                 const IntSetVector& gtMutClustering,
//                                 const IntSetVector& infMutClustering) {
//    uint64_t count = 0;
//    IntVector cellToClusterGT(nrCells, 0);
//    for (uint64_t clusterIdx = 0; clusterIdx < gtCellClustering.size(); ++clusterIdx) {
//        for (uint64_t cell : gtCellClustering[clusterIdx]) {
//            cellToClusterGT[cell] = clusterIdx;
//        }
//    }
//
//    IntVector cellToClusterInf(nrCells , 0);
//    for (uint64_t clusterIdx = 0; clusterIdx < infCellClustering.size(); ++clusterIdx) {
//        for (uint64_t cell : infCellClustering[clusterIdx]) {
//            cellToClusterInf[cell] = clusterIdx;
//        }
//    }
//
//    for (uint64_t cell = 0; cell < nrCells; ++cell) {
//        IntSet cellGenotypeGT;
//        uint64_t clusterIdxGT = cellToClusterGT[cell];
//        for (uint64_t mutClusterIdx : gtGenotypes[clusterIdxGT]) {
//            for (uint64_t mut : gtMutClustering[mutClusterIdx]) {
//                cellGenotypeGT.insert(mut);
//            }
//        }
//
//        IntSet cellGenotypeInf;
//        uint64_t clusterIdxInf = cellToClusterInf[cell];
//        for (uint64_t mutClusterIdx : infGenotypes[clusterIdxInf]) {
//            for (uint64_t mut : infMutClustering[mutClusterIdx]) {
//                cellGenotypeInf.insert(mut);
//            }
//        }
//
//        IntSet S;
//        std::set_intersection(cellGenotypeGT.begin(), cellGenotypeGT.end(),
//                              cellGenotypeInf.begin(), cellGenotypeInf.end(),
//                              std::inserter(S, S.begin()));
//
//        count += S.size();
//    }
//
//    return count / (double) (nrCells * nrMutations);
//}
SNVTree make_snv_tree(  ClusteringPairs edges,    const Genotypes G, int j){
    SNVTree ts;


    for(auto e: edges){
        geno g1 = G.at(e.first)[j];

        geno g2 =G.at(e.second)[j];
            if(g1 != g2){
                ts.insert({g1, g2});
            }

        }

    return ts;
}


SNVTreePresence make_snv_tree_presence(  ClusteringPairs edges,    const Genotypes G, int j){
    SNVTreePresence ts;


    for(auto e: edges){
        geno g1 = G.at(e.first)[j];
        geno_presence g1p;
        g1p.pres = g1.presence();
        g1p.x = g1.x;
        g1p.y = g1.y;
       

        geno g2 =G.at(e.second)[j];
        geno_presence g2p;
        g2p.pres = g2.presence();
        g2p.x = g2.x;
        g2p.y = g2.y;
        

            if(g1p != g2p){
                ts.insert({g1p, g2p});
            }

        }

    return ts;
}
ClusteringPairs get_edges(Tree T){
    ClusteringPairs edges;
    for (uint64_t source = 0; source < T.size(); ++source) {
        if(T[source].size() > 0){
            for (auto target : T[source]) {
                edges.insert({source, target});
            }
        }

    }

    return edges;
}


CNATree make_cna_tree(const ClusteringPairs edges,  const CNAGeno G, int ell){
    CNATree S;

    for(auto e: edges){
        CNA g1 = G.at(e.first).at(ell);
        if (G.at(e.second).count(ell) > 0) {

            CNA g2 = G.at(e.second).at(ell);
            if (g1 != g2) {
                S.insert({g1, g2});
            }
        }


    }

    return S;
}

std::pair<double, double> genotype_tree_accuracy(Tree GT, Tree INF, Genotypes gtGenos, Genotypes infGenos, std::set<uint64_t> snvs, uint64_t m){
    int correct = 0;
    int correct_tree_presence = 0;
    ClusteringPairs gtEdges = get_edges(GT);
    ClusteringPairs infEdges = get_edges(INF);
    for(auto j: snvs){
        if(snvs.count(j) > 0){
            SNVTree gtsnv = make_snv_tree(gtEdges, gtGenos, j);
            SNVTree infsnv =make_snv_tree(infEdges, infGenos, j);
            SNVTreePresence gtsnvp = make_snv_tree_presence(gtEdges, gtGenos, j);
            SNVTreePresence infsnvp =make_snv_tree_presence(infEdges, infGenos, j);
            correct += (gtsnv ==infsnv);
            correct_tree_presence +=  (gtsnvp == infsnvp);
        }


    }
    double perc_snv_tree = (double) correct / snvs.size();
    double perc_snv_tree_presence = (double)  correct_tree_presence / snvs.size();
    return {perc_snv_tree, perc_snv_tree_presence};
}

uint64_t cna_tree_accuracy(Tree GT, Tree INF, CNAGeno gtGenos, CNAGeno infGenos, std::set<int> segments){
    int correct = 0;
    ClusteringPairs gtEdges = get_edges(GT);
    ClusteringPairs infEdges = get_edges(INF);
    for(int ell: segments){
        CNATree gtS = make_cna_tree(gtEdges, gtGenos, ell);
        CNATree infS =make_cna_tree(infEdges, infGenos, ell);
        correct += (gtS == infS);

    }
    return correct;
}

int main(int argc, char **argv) {
    if (argc != 9) {
        std::cerr << "Usage: " << argv[0] << " <GT_TREE> <GT_CELL_CLUSTERING> <GT_MUT_CLUSTERING> <GT_GENOTYPES> <INF_TREE> <INF_CELL_CLUSTERING> <INF_MUT_CLUSTERING> <INF_GENOTYPES>" << std::endl;
        std::cout << "anc_cell_TP,anc_cell_TN,anc_cell_FN,anc_cell_FP,anc_cell_recall,anc_cell_precision"
                  << ",inc_cell_TP,inc_cell_TN,inc_cell_FN,inc_cell_FP,inc_cell_recall,inc_cell_precision"
                  << ",cell_TP,cell_TN,cell_FN,cell_FP,cell_recall,cell_precision,cell_RI,cell_ARI"
                  << ",cell_GT_non_empty" << ",cell_inf_non_empty"
                  << ",anc_mut_TP,anc_mut_TN,anc_mut_FN,anc_mut_FP,anc_mut_recall,anc_mut_precision"
                  << ",inc_mut_TP,inc_mut_TN,inc_mut_FN,inc_mut_FP,inc_mut_recall,inc_mut_precision"
                  << ",mut_TP,mut_TN,mut_FN,mut_FP,mut_recall,mut_precision,mut_RI,mut_ARI"
                  << ",mut_GT_non_empty" << ",mut_inf_non_empty"
                  << ",genotype_pres_similarity" << ",genotype_allele_spec_similarity"<< ",genotype_similarity" <<",cna_genotype_similarity"
                  << ",snv_tree_accuracy" << ",cna_tree_accuracy" << ",inf_segments" << ",inf_mut"
                  << ",full_genotype_presence_similarity" << ",full_genotype_similarity" << ",snv_tree_presence_accuracy"
                  << std::endl;
        return 1;
    }

        // << "," <<  (double) correct_presence / double(n * snvs.size())
        //         << "," <<  (double) correct_sim / double(n * snvs.size())
        //         << "," <<  snv_tree_acc.second 
    std::ifstream inGtTree(argv[1]);
    auto gtTree = parseTree(inGtTree);
    inGtTree.close();

    auto ancPairsGtTree = getAncestralNodePairs(gtTree.second);
  
    auto incPairsGtTree = getIncomparableNodePairs(gtTree.second, ancPairsGtTree, gtTree.first);

//    auto genotypesGtTree = cellClusterGenotypes(gtTree.second);
    
   
    
    std::ifstream inInfTree(argv[5]);
    auto infTree = parseTree(inInfTree);


    auto ancPairsInfTree = getAncestralNodePairs(infTree.second);

    auto incPairsInfTree = getIncomparableNodePairs(infTree.second, ancPairsInfTree, infTree.first);

 
//    auto genotypesInfTree = cellClusterGenotypes(infTree.second);



    // SNV genotype similarity
    // precision / recall for tree metrics
    // number of nonempty cell clusters
    // number of nonempty mutation clusters
    // ancestral pair TP/TN/FN/FP
    // incomparable pair TP/TN/FN/FP



    // cell clustering
    std::ifstream inGtCell(argv[2]);
    auto gt_cell_clustering = parseClustering(inGtCell, 0, 1);
    inGtCell.close();

    std::ifstream inInfCell(argv[6]);
    auto inf_cell_clustering = parseClustering(inInfCell, 0, 1);
    inInfCell.close();



    if (gt_cell_clustering.first != inf_cell_clustering.first) {
        std::cerr << "Error: different number of elements in cell clustering." << std::endl;
    }

    {

        auto anc_pairs_gt = expandClustering(ancPairsGtTree, gt_cell_clustering.second, true);
   
        auto anc_pairs_inf = expandClustering(ancPairsInfTree, inf_cell_clustering.second, true);
    
        uint64_t n = gt_cell_clustering.first;
        uint64_t anc_cell_TP = 0, anc_cell_TN = 0, anc_cell_FN = 0, anc_cell_FP = 0;
        performComparison(n * (n - 1), anc_pairs_gt, anc_pairs_inf,
                          anc_cell_TP, anc_cell_TN, anc_cell_FN, anc_cell_FP);

        std::cout << anc_cell_TP << "," << anc_cell_TN << "," << anc_cell_FN << "," << anc_cell_FP
                  << "," << ((anc_cell_TP + anc_cell_FN) == 0 ? 1. : (double) anc_cell_TP / (double) (anc_cell_TP + anc_cell_FN))
                  << "," << ((anc_cell_TP + anc_cell_FP) == 0 ? 1. : (double) anc_cell_TP / (double) (anc_cell_TP + anc_cell_FP)) << std::flush;
 
    }


    {
//
        auto inc_pairs_gt = expandClustering(incPairsGtTree, gt_cell_clustering.second, false);
        // for(auto pr: inc_pairs_gt){
        //     std::cout << pr.first << ":" << pr.second << std::endl;
        // }

        auto inc_pairs_inf = expandClustering(incPairsInfTree, inf_cell_clustering.second, false);

        uint64_t n = gt_cell_clustering.first;

        uint64_t inc_cell_TP = 0, inc_cell_TN = 0, inc_cell_FN = 0, inc_cell_FP = 0;
        performComparison(n * (n - 1) / 2, inc_pairs_gt, inc_pairs_inf,
                          inc_cell_TP, inc_cell_TN, inc_cell_FN, inc_cell_FP);

        std::cout << "," << inc_cell_TP << "," << inc_cell_TN << "," << inc_cell_FN << "," << inc_cell_FP
                  << "," << ((inc_cell_TP + inc_cell_FN) == 0 ? 1. : (double) inc_cell_TP / (double) (inc_cell_TP + inc_cell_FN))
                  << "," << ((inc_cell_TP + inc_cell_FP) == 0 ? 1. : (double) inc_cell_TP / (double) (inc_cell_TP + inc_cell_FP)) << std::flush;
    }

    {
        

        auto pairs_gt_cell_clustering = makePairs(gt_cell_clustering.second);
        auto pairs_inf_cell_clustering = makePairs(inf_cell_clustering.second);

        uint64_t cell_TP = 0, cell_TN = 0, cell_FN = 0, cell_FP = 0;
        uint64_t n = gt_cell_clustering.first;
        performComparison(n * (n - 1) / 2, pairs_gt_cell_clustering, pairs_inf_cell_clustering,
                          cell_TP, cell_TN, cell_FN, cell_FP);

        std::cout << "," << cell_TP << "," << cell_TN << "," << cell_FN << "," << cell_FP
                  << "," << ((cell_TP + cell_FN) == 0 ? 1. : (double) cell_TP / (double) (cell_TP + cell_FN))
                  << "," << ((cell_TP + cell_FP) == 0 ? 1. : (double) cell_TP / (double) (cell_TP + cell_FP))
                  << "," << (double) (cell_TP + cell_TN) / (double) (cell_TP + cell_FN + cell_FP + cell_TN)
                  << "," << 2.0 * (cell_TP * cell_TN - cell_FN * cell_FP) /
                            (double) ((cell_TP + cell_FN) * (cell_FN + cell_TN) +
                                      (cell_TP + cell_FP) * (cell_FP + cell_TN))
                  << "," << getNrNonEmptyClusters(gt_cell_clustering.second)
                  << "," << getNrNonEmptyClusters(inf_cell_clustering.second)
                  << std::flush;
    }



    // mutation clustering
    std::ifstream inGtMut(argv[3]);
    auto gt_mut_clustering = parseClustering(inGtMut, 0, 1);
    inGtMut.close();

    std::ifstream inInfMut(argv[7]);
    auto inf_mut_clustering = parseClustering(inInfMut, 0, 1);

    inInfMut.close();
    std::vector<std::set<uint64_t>> psi_inverse = inf_mut_clustering.second;
    std::set <uint64_t> snvs;

    for(auto clust: psi_inverse){
        for(auto j: clust){
            snvs.insert(j);
        }
    }

    if (gt_mut_clustering.first != inf_mut_clustering.first) {
        std::cerr << "Error: different number of elements in mutation clustering." << std::endl;
    }

    {
        auto anc_pairs_gt = expandClustering(ancPairsGtTree, gt_mut_clustering.second, true);
        auto anc_pairs_inf = expandClustering(ancPairsInfTree, inf_mut_clustering.second, true);

        uint64_t n = gt_mut_clustering.first;
        uint64_t anc_mut_TP = 0, anc_mut_TN = 0, anc_mut_FN = 0, anc_mut_FP = 0;
        performComparison(n * (n - 1), anc_pairs_gt, anc_pairs_inf,
                          anc_mut_TP, anc_mut_TN, anc_mut_FN, anc_mut_FP);

        std::cout << "," << anc_mut_TP << "," << anc_mut_TN << "," << anc_mut_FN << "," << anc_mut_FP
                  << "," << ((anc_mut_TP + anc_mut_FN) == 0 ? 1. : (double) anc_mut_TP / (double) (anc_mut_TP + anc_mut_FN))
                  << "," << ((anc_mut_TP + anc_mut_FP) == 0 ? 1. : (double) anc_mut_TP / (double) (anc_mut_TP + anc_mut_FP)) << std::flush;
    }

    {
        auto inc_pairs_gt = expandClustering(incPairsGtTree, gt_mut_clustering.second, false);
        auto inc_pairs_inf = expandClustering(incPairsInfTree, inf_mut_clustering.second, false);

        uint64_t n = gt_mut_clustering.first;
        uint64_t inc_mut_TP = 0, inc_mut_TN = 0, inc_mut_FN = 0, inc_mut_FP = 0;
        performComparison(n * (n - 1) / 2, inc_pairs_gt, inc_pairs_inf,
                          inc_mut_TP, inc_mut_TN, inc_mut_FN, inc_mut_FP);

        std::cout << "," << inc_mut_TP << "," << inc_mut_TN << "," << inc_mut_FN << "," << inc_mut_FP
                  << "," << ((inc_mut_TP + inc_mut_FN) == 0 ? 1. : (double) inc_mut_TP / (double) (inc_mut_TP + inc_mut_FN))
                  << "," << ((inc_mut_TP + inc_mut_FP) == 0 ? 1. : (double) inc_mut_TP / (double) (inc_mut_TP + inc_mut_FP)) << std::flush;
    }

    {
        auto pairs_gt_mut_clustering = makePairs(gt_mut_clustering.second);
        auto pairs_inf_mut_clustering = makePairs(inf_mut_clustering.second);

        uint64_t n = gt_mut_clustering.first;
        uint64_t mut_TP = 0, mut_TN = 0, mut_FN = 0, mut_FP = 0;
        performComparison(n * (n - 1) / 2, pairs_gt_mut_clustering, pairs_inf_mut_clustering,
                          mut_TP, mut_TN, mut_FN, mut_FP);

        std::cout << "," << mut_TP << "," << mut_TN << "," << mut_FN << "," << mut_FP
                  << "," << (double) mut_TP / (double) (mut_TP + mut_FN)
                  << "," << (double) mut_TP / (double) (mut_TP + mut_FP)
                  << "," << (double) (mut_TP + mut_TN) / (double) (mut_TP + mut_FN + mut_FP + mut_TN)
                  << "," << 2.0 * (mut_TP * mut_TN - mut_FN * mut_FP) /
                            (double) ((mut_TP + mut_FN) * (mut_FN + mut_TN) +
                                      (mut_TP + mut_FP) * (mut_FP + mut_TN))
                  << "," << getNrNonEmptyClusters(gt_mut_clustering.second)
                  << "," << getNrNonEmptyClusters(inf_mut_clustering.second)
                  << std::flush;
    }

  
    uint64_t pres = 0;
    uint64_t allspec = 0;
    uint64_t sim =0;
    uint64_t cna =0;
    uint64_t correct_presence =0;
    uint64_t correct_sim =0;
    // uint64_t m = 5000;  // Directly initialize m with 5000
    uint64_t m = *snvs.rbegin() + 1;
    // uint64_t m =  gt_mut_clustering.first;


    std::ifstream gtGenos(argv[4]);

    auto gtGenotypes = parseGenotypes(gtGenos, gtTree.first, m);
    gtGenos.close();
    

    CNAGeno cnGT = gtGenotypes.second;
    std::set<int> segments;
    for(const auto pair: cnGT){
        for(const auto pr: pair.second){
            segments.insert(pr.first);
        }
        break;
    }


    std::map<uint64_t, IntPair> cellMap = makeCellMap(gt_cell_clustering.second, inf_cell_clustering.second);
    uint64_t n = cellMap.size();

    std::ifstream infGenos(argv[8]);
    auto infGenotypes = parseGenotypes(infGenos, infTree.first, m);
    infGenos.close();

    CNAGeno cnINF = infGenotypes.second;
    std::set<int> infsegments;
    for(const auto pair: cnINF){
        for(const auto pr: pair.second){
            infsegments.insert(pr.first);
        }
        break;
    }

    std::pair<double, double>  snv_tree_acc;
    snv_tree_acc =genotype_tree_accuracy(gtTree.second, infTree.second, gtGenotypes.first, infGenotypes.first, snvs, m);


    computeGenotypeSimilarity(gtGenotypes, infGenotypes, cellMap, snvs, m, pres, allspec, sim, cna, correct_sim, correct_presence);
    uint64_t cna_tree_correct = cna_tree_accuracy(gtTree.second, infTree.second, gtGenotypes.second, infGenotypes.second, infsegments);

    std::cout << "," << (double) pres / double(n * snvs.size())
              <<"," << (double) allspec / double(n * snvs.size())
              <<"," << (double) sim / double(n * snvs.size())
              <<"," << (double) cna / double(n * segments.size())
              << "," << snv_tree_acc.first 
              << "," << (double) cna_tree_correct / double(segments.size())
              << "," <<  infsegments.size() << "," << snvs.size()
               << "," <<  (double) correct_presence / double(n * snvs.size())
                << "," <<  (double) correct_sim / double(n * snvs.size())
                << "," <<  snv_tree_acc.second 
              << std::flush;
//    std::cout << "," << computeGenotypeSimilarity(gt_cell_clustering.first, gt_mut_clustering.first,
//                                                  cellClusterGenotypes(gtTree.second), cellClusterGenotypes(infTree.second),
//                                                  gt_cell_clustering.second, inf_cell_clustering.second,
//                                                  gt_mut_clustering.second, inf_mut_clustering.second)
//              << std::endl;

    return 0;
}
