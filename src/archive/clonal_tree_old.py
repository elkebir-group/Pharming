
import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter
from itertools import product, chain, combinations
from scipy.stats import binom
from scipy.special import logsumexp
from itertools import chain
import pickle 
import pygraphviz as pgv
# import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from sklearn.metrics.cluster import adjusted_rand_score
from genotype import genotype


SET3_HEX = [
    '#8dd3c7',
    '#ffffb3',
    '#bebada',
    '#fb8072',
    '#80b1d3',
    '#fdb462',
    '#b3de69',
    '#fccde5',
    '#d9d9d9',
    '#bc80bd',
    '#ccebc5',
    '#ffed6f'
]


class ClonalTree:
    """
    A class to model a clonal tree with associated SNV/CNA genotypes and cell clustering

    ...

    Attributes
    ----------
    key : int
        the unique key for the clonal tree
    
    tree : networkx DiGraph
        the clonal tree graph
    
    cell_mapping : dict
        a dictionary with node ids as keys and np.arrays with the indices of cells attached to that node
    
    mut_mapping : dict
        a dictionary with node ids as keys and np.arrays with the indices of gained SNVs assigned to the 
        incoming edge of that node
    
    
    mut_loss_mapping : dict 
        a dictionary with node ids as keys and np.arrays with the indices of lost SNVs assigned to the 
        incoming edge of that node
    
    
    snv_genotypes: dict
        a dictionary with node ids as keys and an np.array of length M of the number of mutated copies for each SNV
    
    loglikelihood : float
        the loglikelihood of the clonal tree for a given set of data


    
   

    Methods
    -------

    get_tip_cells(t)
        returns the set of cells in node t

    get_tip_muts(t)
        returns the set of SNVs on the incoming edge of node t

    set_key(key)
        updates the key identifier for the clonal tree

    has_loss()
        returns a boolean if any SNVs are lost

    get_seed_by_node(node, lamb, tau, anc_muts)
        returns a Seed for the given node provided it meets criteria specified by lamb, tau

    inorder_traversal(root)
        performs an inorder traversal of the tree for relabeling purpuses

    relabel()
        relabels the nodes consecutively from an inorder traversal and updates the mapping dictionaries

    merge(tree, l)
        merges the clonal tree with a given tree by replacing the the leaf node l of the clonal tree
        with the root node of the given tree

    snv_genotypes()
        returns a dictionary with nodes as keys containing the snv genotype vectors

    presence_by_node(m, cells, node)
        returns an expanded snv_genotype matrix (cells x m) for a given set of cells for a specific node

    event_by_node(m, cells, node, bin_mapping)
        returns an expanded cna_genotype matrix (cells x m) for a given set of cells for a specific node

    cna_genotype( n,  nbins)
        returns a vector of length nbins containing the CNA genotype vector for node n

    save(path)
        a helper method to save a clonal tree to disk via pickle

    rdr_likelihood_by_node(n, cnn_hmm, cells=None)
        computes the likelihood of the bin count data for a specified node
    
    tree_png(fname, chrom_map_file: str=None)
        saves the png drawing of the clonal tree to disk
    
    tree_dot(fname, chrom_map_file: str=None)
        saves the dot string for the clonal tree to disk
    
    get_leaves()
        returns the node ids of a all leaf nodes
    
    compute_likelihood(data)
        computes the likelihood of the given data for the current clonal tree, genotypes and cell clustering
    
    compute_likelihood_by_node(n, data)
        computes the likelihood of the given data for a specified node  
        given the current clonal tree, genotypes and cell clustering
    
    get_loglikelihood()
        accesses the current values of the loglikelihood, variannt read count loglikelihood and binned 
        read count loglikelihood
    
    generate_results(cell_lookup, mut_lookup)
        converts all the internal mapping dictionaries to output dataframes 
    
    compute_variant_likelihood_by_node_without_events(node, like0, like1, bin_mapping=None)
        computes the loglikelihood of the variant read count data for the specified node of the clonal tree 
        without making use of the CNA genotypes 
    
    compute_variant_likelihood_by_node_with_events(node, like0, like1, bin_mapping=None)
        computes the loglikelihood of the variant read count data for the specified node of the clonal tree 
        making use of the CNA genotypes 

    find_root()
        returns the node id of the root of the clonal tree
    
    get_all_cells()
        returns an np.array of all cell indices currently assigned in the clonal tree
    
    get_all_muts()
        returns an np.array of all SNV indices currently assigned in the clonal tree


    """

    def __init__(self, key, tree, mut_mapping, mut_loss_mapping=None,cell_mapping=None, mutated_copies=None ):
        self.tree: nx.DiGraph = tree
        
        self.root= self.find_root()
        if cell_mapping is None:
            self.cell_mapping = {}
        else:    
            self.cell_mapping = cell_mapping

        self.mut_mapping = mut_mapping

        self.set_genotypes()

        self.mutated_copies = mutated_copies

        if mut_loss_mapping is None:
            self.mut_loss_mapping = {}
        else:
            self.mut_loss_mapping = mut_loss_mapping


        self.key = key
        self.loglikelihood = None

        # self.node_attrs = {k for node_data in self.tree.nodes.data() for k in node_data[1].keys()}
        # self.edge_attrs =  {k for edge_data in self.tree.edges.data() for k in edge_data[1].keys()}
    
    def set_genotypes(self):
        # muts = self.get_all_muts()
        genotypes = {v: {self.key: {} } for v in self.clones()}
        
        for v, muts in self.mut_mapping.items():
  
            x,y = self.tree.nodes[v]["genotype"]
            desc_nodes = self.get_descendants(v)
            not_pres_nodes = set(self.clones()) -{v} - desc_nodes
            # ancestors_nodes = self.get_ancestors(v)
            
            
       
            for m in muts:
        
                genotypes[v][self.key][m] = genotype(x,y,1,0)
                for a in not_pres_nodes:
                    x_1, y_1  = self.tree.nodes[a]["genotype"]
                    genotypes[a][self.key][m] = genotype(x_1, y_1,0,0)
                for d in desc_nodes:
                    x_1, y_1  = self.tree.nodes[d]["genotype"]
                    genotypes[d][self.key][m] = genotype(x_1, y_1,1,0)
        self.genotypes = genotypes  


    def get_ancestors(self, v):

        def find_ancestors(tree, node, ancestors=None):
            if ancestors is None:
                ancestors = set()

            predecessors = list(tree.predecessors(node))

            for parent in predecessors:
                ancestors.add(parent)
                find_ancestors(tree, parent, ancestors)

            return ancestors
        return find_ancestors(self.tree, v)
    
    def get_descendants(self, v):
        def find_descendants(tree, node, descendants=None):
            if descendants is None:
                descendants = set()

            successors = list(tree.successors(node))

            for child in successors:
                descendants.add(child)
                find_descendants(tree, child, descendants)

            return descendants
        return find_descendants(self.tree,v)
    
    
    def clones(self):
        return list(self.tree.nodes)
    
    def edges(self):
        return list(self.tree.edges)
    
    def get_clone_cells(self, n):
        if n in self.cell_mapping:
            return np.concatenate([self.cell_mapping[t][k] for k in self.cell_mapping[t]]).tolist()
        else:
            return []

    def parent(self, v):
        if v ==self.root:
            return None 
        else:
            return list(self.tree.predecessors(v))[0]
    
    def children(self, v):
        return list(self.tree.neighbors[v])
    
    def get_clone_snvs(self,n):
        clone = nx.shortest_path(self.tree, source=self.root, target=n)
        snvs = []
        for c in clone:
            snvs += self.mut_mapping[c]
        return snvs
        
    def set_key(self, key):
        self.key = key

    def has_loss(self):
        return len(self.mut_loss_mapping) > 0
    
    def set_cell_mapping(self, cell_mapping):
        self.cell_mapping  = cell_mapping

    def __str__(self) -> str:
        all_cells = self.get_all_cells()
        all_snvs  = self.get_all_muts()
        mystr =f"\nClonal Tree {self.key}"
        mystr += "\n--------------------------------"

        mystr +=f"\n{len(list(self.tree.nodes))} Nodes and {len(list(self.tree.edges))} Edges"
        mystr += f"\n{len(all_snvs)} SNVs"
        mystr+= f"\n{len(all_cells)} cells\n"


        node_order = list(nx.dfs_preorder_nodes(self.tree, source=self.root))
        mystr += "\nParent->Child\tEdge Type"
        mystr += "\n--------------------------------"
        for u in node_order:
        
            for v in self.tree.neighbors(u):
                etype= ""
                # if "event" in self.tree.edge_attrs:
                #     etype = self.tree[u][v]["event"]
                # else:
                #     etype = "NA"
           
                mystr+= f"\nNode {u}->Node {v}\t{etype}"
               
    
        mystr += "\n\nNode\tSNVs\tCells\tGenotype"
        mystr += "\n--------------------------------"
        for n in  node_order:
            genotype = self.tree.nodes[n]['genotype']
            if n in self.cell_mapping:
                ncells = len(self.cell_mapping[n])
            else:
                ncells = 0
            if n in self.mut_mapping:
                nmuts = len(self.mut_mapping[n])
            else:
                nmuts = 0
            mystr += f"\n{n}\t{nmuts}\t{ncells}\t{genotype}"
        mystr += "\n"
        mystr += f"\nLog-likelihood: {self.loglikelihood}\n"
        return mystr
  




 
   

    def snv_genotypes(self, m=None):
        if m is None:
            m = len(self.get_all_muts())
        
        y_dict = {}
        for node in self.mut_mapping:
            y = np.zeros(m, dtype=int)
            ancestral_present_muts = self.get_ancestral_muts(node).astype(int)
            present_muts = np.concatenate(
                [ancestral_present_muts, self.mut_mapping[node]])
            if node in self.mut_loss_mapping:
                present_muts = np.setdiff1d(
                    present_muts, self.mut_loss_mapping[node])
            y[present_muts.astype(int)] = 1
            y_dict[node] = y

        return y_dict
    
    def cell_cluster_genotypes(self, n=None):
        if n is None:
            n = len(self.get_all_cells())
        c_dict = {}
        for node in self.cell_mapping:
            c = np.zeros(n, dtype=int)
            present_clade = list(nx.dfs_preorder_nodes(self.tree, node))
            present_cells = np.concatenate([self.cell_mapping[i][0] for i in present_clade])
   
            c[present_cells] = 1
            c_dict[node] = c

        return c_dict

    def presence_by_node(self, m, cells, node):

        presence = np.zeros(shape=(len(cells), m), dtype=int)
        if len(cells) == 0:
            return presence

        ancestral_present_muts = self.get_ancestral_muts(node).astype(int)
        present_muts = np.concatenate(
            [ancestral_present_muts, self.mut_mapping[node]])
        if node in self.mut_loss_mapping:
            present_muts = np.setdiff1d(
                present_muts, self.mut_loss_mapping[node])

        presence[:, present_muts] = 1

        return presence  # n x m binary matrix cell y_ij=1 if mutation j  is harbored in cell i





    def draw(self, fname, cmap='Set3'):

        mut_count = {n : len(self.mut_mapping[n]) for n in self.mut_mapping}
        cell_count = {n : len(self.cell_mapping[n]) for n in self.cell_mapping}
        labels = {}
        # color_values, colormap = self.create_color_map(cmap)
        for n in self.tree:
                labels[n] = str(n)
                if n in self.cell_mapping:
                    if cell_count[n] > 0:
                        labels[n] += "\nCells:" + str(cell_count[n])
                # SNV
                if n in self.mut_mapping:
                    if mut_count[n] > 0:
                        labels[n] += "\nSNVs:" + str(mut_count[n])
        like_label = f"Segment {self.key}\n"
        tree = pgv.AGraph(strict=False, directed=False)
        tree.node_attr['style']='filled'
        if self.loglikelihood is not None:
            total_like = np.round(self.loglikelihood)
            like_label += f"Log Likelihood: {total_like}"
            tree.graph_attr["label"] = like_label
 
        # colormap = cm.get_cmap(cmap)
        for n in self.tree:

            tree.add_node(n, label=labels[n])
      
            node_attr = tree.get_node(n)
            try:
                x,y = self.tree.nodes[n]["genotype"]
                color_value = x+y
            except:
                color_value = None
        
            if color_value is not None:
               
                
                # color = colormap(color_value)
                # hex_color = mcolors.rgb2hex(color)
                node_attr.attr['fillcolor'] =SET3_HEX[color_value]
                # node_attr['fillcolor'] = hex_color
    
        tree.add_edges_from(list(self.tree.edges))
        tree.layout("dot")
        tree.draw(fname)
  


    def save_text(self, path):
        
        
        leafs = [n for n in self.tree.nodes if len(list(self.tree.successors(n))) ==0]
      
                    
        with open(path, "w+") as file:
            file.write(f"{len(list(self.tree.edges))} #edges\n")
            for u,v in list(self.tree.edges):
                file.write(f"{u} {v}\n")
            file.write(f"{len(leafs)} #leaves\n")
            for l in leafs:
                file.write(f"{l}\n")
            
     
            

    # def tree_png(self, fname, chrom_map_file: str = None):
    #     #self.node_events = self.relabel()
    #     self.node_events = None
    #     bin2chrom = None
    #     if chrom_map_file is not None:
    #         bin2chrom = pd.read_csv(chrom_map_file, names=[
    #                                 "chrom", "arm", "start", "end"])
    #     dt = DrawClonalTree(self, bin2chrom)
    #     dt.savePNG(fname)

    # def tree_dot(self, fname, chrom_map_file: str = None):
    #     bin2chrom = None
    #     if chrom_map_file is not None:
    #         bin2chrom = pd.read_csv(chrom_map_file, names=[
    #                                 "chrom", "arm", "start", "end"])
    #     dt = DrawClonalTree(self, bin2chrom)
    #     dt.saveDOT(fname)

    def get_node_muts(self, t):
        if t in self.mut_mapping:
            return self.mut_mapping[t]
        return []
    
    
    def get_node_cells(self, t):
        if t in self.cell_mapping:
            return self.cell_mapping[t]
        return []

    def get_leaves(self):
        leaves = [l for l in list(self.tree.nodes())
                  if self.tree.out_degree(l) == 0]
        return leaves

    def get_all_cells(self):
        all_cells = []
        for n in self.cell_mapping:
            all_cells += self.cell_mapping[n]
        return all_cells

    def get_cell_clusters(self):
        n_cell = len(self.get_all_cells())
        clusters = np.zeros(n_cell, dtype=int)
        for cluster, cells in self.cell_mapping.items():
            if len(cells) > 0:
                clusters[cells] = cluster
        return clusters

    def get_all_muts(self):
        muts = list(chain.from_iterable([mut for n,mut in self.mut_mapping.items()]))
        muts.sort()
        return muts

    def get_mut_clusters(self, n_mut=None):
        if n_mut is None:
            muts = self.get_all_muts()
            clusters = np.zeros(len(muts), dtype=int)
            mut_to_index= {m: i for i,m in enumerate(muts)}
        else:
            clusters  = np.full(shape=n_mut, fill_value=-1)
        for cluster, muts in self.mut_mapping.items():
            mut_indices = [mut_to_index[m] for m in muts]
            if len(muts) > 0:
                clusters[mut_indices] = cluster
        return clusters

    @staticmethod
    def recall(gt_pairs, inf_pairs) -> float:
        if sum(gt_pairs.values()) == 0:
            return 1
        return  sum((gt_pairs & inf_pairs).values())\
                           / sum(gt_pairs.values())



    @staticmethod
    def ari(v1, v2) -> float:
        return adjusted_rand_score(v1, v2)
    
    def compute_cell_ari(self, obj) -> float:
        gt_cell = self.get_cell_clusters()
        pred_cell = obj.get_cell_clusters()
        return self.ari(gt_cell, pred_cell)

    
    def compute_mut_ari(self,obj) -> float:
         gt_mut = self.get_mut_clusters()
         pred_mut = obj.get_mut_clusters()
         return self.ari(gt_mut, pred_mut)





    def get_cell_cluster(self, c):
        for node in self.cell_mapping:
            if c in self.cell_mapping[node][0]:
                return node

    def get_snv_cluster(self, s):
        for node in self.mut_mapping:
            if s in self.mut_mapping[node]:
                return node  

    @staticmethod
    def is_incomparable(graph: nx.DiGraph, u, v) -> bool:
        for path in nx.all_simple_paths(graph, source=0, target=v):
            if u in path:
                return False
        for path in nx.all_simple_paths(graph, source=0, target=u):
            if v in path:
                return False
        return True

    
    def ancestor_pair_recall(self, obj, include_loss=False):
        if type(obj) != ClonalTree:
            raise TypeError("Comparable must be a ClonalTree") 
        else:
           
            ancestral = self.get_ancestor_pairs(include_loss)
            
            return self.recall(ancestral, obj.get_ancestor_pairs(include_loss))
    
    def ancestor_cell_pair_recall(self, obj):
        if type(obj) != ClonalTree:
            raise TypeError("Comparable must be a ClonalTree") 
        else:
           
            ancestral = self.get_cell_ancestor_pairs()
            return self.recall(ancestral, obj.get_cell_ancestor_pairs())

    def incomp_pair_recall(self, obj, include_loss=False):
        if type(obj) != ClonalTree:
            raise TypeError("Comparable must be a ClonalTree") 
        else:
           
            ancestral = self.get_incomparable_pairs(include_loss)
            
            return self.recall(ancestral, obj.get_incomparable_pairs(include_loss))
    
    
    def clustered_pair_recall(self, obj, feature="cell") -> float:
        if feature == "cell":
            mapping= self.cell_mapping 
            obj_mapping = obj.cell_mapping
        else:
            mapping = self.mut_mapping
            obj_mapping = obj.mut_mapping
        
        clustered = self.get_cluster_pairs(mapping)
        recall = self.recall(clustered, obj.get_cluster_pairs(obj_mapping))
        return recall 
    
    def incomp_cell_pair_recall(self, obj):
        if type(obj) != ClonalTree:
            raise TypeError("Comparable must be a ClonalTree") 
        else:
           
            ancestral = self.get_cell_incomparable_pairs()
            return self.recall(ancestral, obj.get_cell_incomparable_pairs())

    def get_clone_proportions(self):
        clone_prop = {}
        ncells = sum([len(val) for key,  val in self.cell_mapping.items()])
        for n in self.tree:
            if n in self.cell_mapping:
                clone_prop[n] = len(self.cell_mapping[n])/ncells
            else:
                clone_prop[n] = 0
        return pd.Series(clone_prop).sort_index()
        
    def get_incomparable_pairs(self, include_loss: bool=True) -> Counter:
        mut_mapping = self.update_mapping(self.tree,self.mut_mapping)
        pairs = Counter()
        if include_loss:
            raise NotImplementedError("not impletmented")
            # for u, v in combinations(self.tree.nodes, 2):
            #     if include_loss and u in self.mut_loss_mapping:
            #         u_loss = self.mut_loss_mapping[u]
            #     else:
            #         u_loss = []
            #     if include_loss and v in self.mut_loss_mapping:
            #         v_loss = self.mut_loss_mapping[v]
            #     else:
            #         v_loss = []

            #     if self.is_incomparable(self.tree, u, v):
            #         for mut1 in chain(
            #             product(self.mut_mapping[u], [1]),
            #             product(u_loss, [0])
            #         ):
            #             for mut2 in chain(
            #                 product(self.mut_mapping[v], [1]),
            #                 product(v_loss, [0])
            #             ):
            #                 if mut1 < mut2:
            #                     pairs[(mut1, mut2)] += 1
            #                 else:
            #                     pairs[(mut2, mut1)] += 1
        else:
            for u, v in combinations(self.tree.nodes, 2):
                if self.is_incomparable(self.tree, u, v):
                    pairs.update((min(a, b), max(a, b)) for a, b in product(mut_mapping[u], mut_mapping[v]))
        return pairs

    @staticmethod
    def update_mapping(tree, mapping):
        full_mapping = mapping.copy()
        for n in tree:
            if n not in full_mapping:
                full_mapping[n] = []
        return full_mapping

    def score_cells(self, obj) -> dict:
       
        if type(obj) != ClonalTree:
            raise TypeError("Comparable must be a ClonalTree") 
        else:
             
             scores = {
                 'feature' : 'cell',
                 'ari' : self.compute_cell_ari(obj),
                 'anc_pair_recall' : self.ancestor_cell_pair_recall(obj),
                 'incomp_pair_recall': self.incomp_cell_pair_recall(obj),
                 'clustered_pair_recall' : self.clustered_pair_recall(obj, feature="cell"),
                 'gt_nodes' : len(self.cell_mapping),
                 'n_assigned': len(self.get_all_cells()),
                 'inf_nodes' : len(obj.cell_mapping)
             }
             return scores 
    
    def score_snvs(self, obj) -> dict:
        if type(obj) != ClonalTree:
            raise TypeError("Comparable must be a ClonalTree")
        else:
             
             scores = {
                 'feature' : 'SNV',
                 'ari' : self.compute_mut_ari(obj),
                 'anc_pair_recall' : self.ancestor_pair_recall(obj),
                 'incomp_pair_recall': self.incomp_pair_recall(obj),
                 'clustered_pair_recall' : self.clustered_pair_recall(obj, feature="snv"),
                 'gt_nodes' : len(self.mut_mapping),
                'n_assigned': len(self.get_all_muts()),
                'inf_nodes' : len(obj.mut_mapping)
             }
             return scores 
    
    def score(self, obj):
        if type(obj) != ClonalTree:
            raise TypeError("Comparable must be a ClonalTree")
        else:
            cell_scores = self.score_cells(obj)
            snv_scores = self.score_snvs(obj)   

            return pd.concat([pd.DataFrame(cell_scores, index=[self.key]), pd.DataFrame(snv_scores,index=[self.key])])     
           

    def get_cluster_pairs(self, mapping) -> Counter:
        mapping = self.update_mapping(self.tree, mapping)
        pairs = Counter()
        for node in self.tree.nodes:
                pairs.update(combinations(sorted(mapping[node]), 2))
        return pairs

        # if include_loss:
        #     raise NotImplementedError("not impletmented")
        #     for node in self.tree.nodes:
        #         if include_loss and node in self.mut_loss_mapping:
        #             node_loss = self.mut_loss_mapping[node]
        #         else:
        #             node_loss = []
        #         for mut1, mut2 in combinations(
        #             chain(
        #                 product(self.mut_mapping[node], [1]),
        #                 product(node_loss, [0])
        #             ),
        #             2
        #         ):
        #             if mut1 < mut2:
        #                 pairs[(mut1, mut2)] += 1
        #             else:
        #                 pairs[(mut2, mut1)] += 1
        # else:
    def get_cell_ancestor_pairs(self) -> Counter:
        mapping = self.update_mapping(self.tree,self.cell_mapping)
        pairs = Counter()
        for node in self.tree.nodes:
            for children in nx.dfs_successors(self.tree, source=node).values():
                for child in children:
                    pairs.update(product(mapping[node], mapping[child]))
        return pairs

    def get_cell_cluster_pairs(self) -> Counter:
        pairs = Counter()
        for node in self.tree.nodes:
            pairs.update(combinations(sorted(self.cell_mapping[node][0]), 2))
        return pairs

    def get_cell_incomparable_pairs(self) -> Counter:
        cell_mapping= self.update_mapping(self.tree, self.cell_mapping)
        pairs = Counter()
        for u, v in combinations(self.tree.nodes, 2):
            if self.is_incomparable(self.tree, u, v):
                pairs.update((min(a, b), max(a, b)) for a, b in product(cell_mapping[u], cell_mapping[v]))
           
        return pairs

    def get_ancestor_pairs(self, include_loss: bool=True) -> Counter:
        mut_mapping = self.update_mapping(self.tree, self.mut_mapping)
        pairs = Counter()
        if include_loss:
            raise NotImplementedError("not impletmented")
            # for node in self.tree.nodes:
            #     for children in nx.dfs_successors(self.tree, source=node).values():
            #         for child in children:
            #             if include_loss and node in self.mut_loss_mapping:
            #                 node_loss = self.mut_loss_mapping[node]
            #             else:
            #                 node_loss = []

            #             for mut1 in chain(
            #                 product(self.mut_mapping[node], [1]),
            #                 product(node_loss, [0])
            #             ):
            #                 if include_loss and child in self.mut_loss_mapping:
            #                     child_loss = self.mut_loss_mapping[child]
            #                 else:
            #                     child_loss = []
            #                 for mut2 in chain(
            #                     product(self.mut_mapping[child], [1]),
            #                     product(child_loss, [0])
            #                 ):
            #                     pairs[(mut1, mut2)] += 1
        else:
            for node in self.tree.nodes:
                for children in nx.dfs_successors(self.tree, source=node).values():
                    for child in children:
                        pairs.update(product(mut_mapping[node], mut_mapping[child]))
        return pairs

        # return pairs
    @staticmethod
    def likelihood_function( a,d,y,c, alpha):
        
        # if d ==0:
        #     return 1e-10
        
        vaf = (1/c)*y
        # elif y ==0:

        #     val =  binom.pmf(a,d,alpha)
        #     return val
          
        # else:
            # vaf = np.arange(1, c)/c
            # vaf  = 1/c
            # vaf = 0.5
        adjusted_vaf =  vaf*(1- alpha) + (1-vaf)*(alpha/3)
        val = binom.logpmf(a,d,adjusted_vaf)

        # val[np.isnan(val)] = np.log(1e-10)
            # return binom.pmf(a,d,adjusted_vaf)
        return val.sum(axis=1)
    
    def compute_map_likelihood(self, data, seg, alpha):
        # print(self.tree_to_string(self.tree))
        clone_order = list(nx.dfs_preorder_nodes(self.tree, source=self.root))
        # clone_order = [c for c in clone_order if c != self.root]
        cells_by_cn = data.cells_by_cn(seg)
        total_cn_states = {}
        like_list =[]
        seg_snvs = data.seg_to_snvs[seg]
        m = len(seg_snvs)
        cell_mapping = {}
        cell_mapping[self.root] = []
        for n in clone_order:
            x,y = self.tree.nodes[n]["genotype"]
            total_cn = x+y 
            total_cn_states[n]= total_cn
            cna_geno = np.full(m, total_cn, dtype=int).reshape(1,-1)
            clone = nx.shortest_path(self.tree, source=self.root, target=n)
            # if self.tree[clone[-2]][clone[-1]]["event"]== "mutation":
            snvs = []
            for c in clone:
                if c in self.mut_mapping:
                    snvs += self.mut_mapping[c]
                
            y_vec = np.zeros(shape=data.M, dtype=int)
           
            y_vec[snvs] = 1
            y_vec= y_vec[seg_snvs]
            y_vec = y_vec.reshape(1,-1)

            # cell_by_snv_like = np.zeros((data.var.shape[0],self.m))
            # for i in range(self.data.var.shape[0]):
            #     for j in range(self.m):
            #         a = self.data.var[i,j]
            #         d = self.data.total[i,j]
            #         c = cna_geno[:,j][0]
            #         alpha= self.alpha 
            #         y = y_vec[:,j][0]
                    
            #         out = likelihood_function(a,d,y,c,alpha)
            #         if np.log(out) > 0:
            #             print(f"a:{a} d:{d} y:{y} c:{c}: prob: {out} logprob:{np.log(out)}")
            #         cell_by_snv_like[i,j] =out
            # print(cell_by_snv_like)
            cell_like = self.likelihood_function(data.var[:, seg_snvs], data.total[:,seg_snvs], y_vec, cna_geno, alpha)
         
            # assert(np.array_equal(cell_by_snv_like, cell_by_snv_like2))
            # print(cell_by_snv_like2)
            # cell_by_snv_like = np.log(cell_by_snv_like)
            # print(cell_by_snv_like)
            # print(f"node: {n} :{cell_by_snv_like.sum(axis=1)}")
            like_list.append(cell_like)
     

        #rows are nodes and columns are cells
        cell_likes = np.vstack(like_list)

        likelihood =0
        cell_assign = {}
 
        #now only consider clone assignments for cells with the correct CN
        for s in cells_by_cn:
            # cn = self.total_cn_by_sample[s]
            clone_filter = [total_cn_states[n]==s for n in clone_order]

            clones = [c for i,c in enumerate(clone_order) if clone_filter[i]]
            for c in clones:
                cell_mapping[c] = []
            cells = cells_by_cn[s]
            sample_likes = cell_likes[clone_filter][:, cells]

            like = np.max(sample_likes, axis=0).sum()
            likelihood += like
            map_clone = np.argmax(sample_likes, axis=0)
            for c, m in zip(cells, map_clone):
                cell_assign[c]= clones[m]
   
                cell_mapping[clones[m]].append(c)
          
            

      
      
        # self.print_verb(f"Log Likelihood: {self.likelihood}")
        # if self.verbose:
        #     for n in self.tree:

        #         print(f"Node{n}: {len(cell_mapping[n])} cells assigned")

        # cell_mapping = {n: np.array(cell_mapping[n], dtype=int) for n in self.tree}
        # self.print_verb(self)
        self.cell_mapping = cell_mapping
        return  likelihood
    
    def compute_marginal_likelihood(self, data,seg, alpha):
            seg_snvs = data.seg_to_snvs[seg]
            all_nodes = []
            total_cn = {}
            prob = {}
            for n in self.tree:
                x,y = self.tree.nodes[n]["genotype"]
                total_cn[n] = x+y 
                if x+y in prob:
                    prob[x+y] += 1
                else:
                    prob[x+y] =1
            prob = {key: 1/val for key, val in prob.items()}
            

            for n in self.tree:
         
                vaf = np.zeros(data.M)
        
                cn = total_cn[n]
                snvs = self.get_ancestral_muts(n)
                # for s in snvs:
                vaf[snvs]= [self.mutated_copies[s] for s in snvs]
                vaf= vaf/cn
                # vaf[snvs] = self.mutated_copies[snvs]/cn

                adj_vaf =vaf*(1- alpha) + (1-vaf)*(alpha/3)
                adj_vaf = adj_vaf[seg_snvs].reshape(1,-1)
          
                #should have dim = cells x snvs
                part1 = binom.pmf(data.var[:, seg_snvs], data.total[:,seg_snvs], adj_vaf)
                part2 =(data.copy_numbers[:, seg] == cn)*prob[cn]
                part2 = part2.reshape(-1,1)
                node_like= np.log(part1*part2)
                node_like = node_like.sum(axis=1)

            
                all_nodes.append(node_like)
         
            loglikelihood = logsumexp(all_nodes, axis=0).sum()
            
            return loglikelihood
    
    def compute_likelihood(self, data, seg, alpha=0.001, attachment="marginal", use_existing=False):
        if attachment =="marginal":
            self.loglikelihood = self.compute_marginal_likelihood(data,seg, alpha)
        elif attachment =="map" and not use_existing:
            self.loglikelihood = self.compute_map_likelihood(data, seg, alpha)
        else:
            if len(self.cell_mapping) > 0:
                raise NotImplementedError("not impletmented")
                self.loglikelihood = sum([self.compute_likelihood_by_node(n) for n in self.cell_mapping])
            else:
                self.loglikelihood = np.NINF
        return self.loglikelihood

        # seg_snvs = data.seg_to_snvs[seg]
      
        # # self.loglikelihood_dict = {"total": 0, "variant": 0, "bin": 0}
        # self.node_likelihood = {}
        # if attachment =="marginal":
            # all_nodes = []
            # total_cn = {}
            # prob = {}
            # for n in self.tree:
            #     x,y = self.tree.nodes[n]["genotype"]
            #     total_cn[n] = x+y 
            #     if x+y in prob:
            #         prob[x+y] += 1
            #     else:
            #         prob[x+y] =1
            # prob = {key: 1/val for key, val in prob.items()}
            

            # for n in self.tree:
         
            #     vaf = np.zeros(data.M)
        
            #     cn = total_cn[n]
            #     snvs = self.get_ancestral_muts(n)
            #     # for s in snvs:
            #     vaf[snvs]= [self.mutated_copies[s] for s in snvs]
            #     vaf= vaf/cn
            #     # vaf[snvs] = self.mutated_copies[snvs]/cn

            #     adj_vaf =vaf*(1- alpha) + (1-vaf)*(alpha/3)
            #     adj_vaf = adj_vaf[seg_snvs].reshape(1,-1)
          
            #     #should have dim = cells x snvs
            #     part1 = binom.pmf(data.var[:, seg_snvs], data.total[:,seg_snvs], adj_vaf)
            #     part2 =(data.copy_numbers[:, seg] == cn)*prob[cn]
            #     part2 = part2.reshape(-1,1)
            #     node_like= np.log(part1*part2)
            #     node_like = node_like.sum(axis=1)

            
            #     all_nodes.append(node_like)
         
            # self.loglikelihood = logsumexp(all_nodes, axis=0).sum()
            
            # return self.loglikelihood

         


                


                # for key in self.loglikelihood_dict:
                #  self.loglikelihood_dict[key] += self.node_likelihood[n][key]

    
        # elif attachment == 'map' and use_existing:
        #     for n in self.cell_mapping:
        #         if len(self.get_tip_cells(n)) > 0:
        #             self.node_likelihood[n] = self.compute_likelihood_by_node(n, data)
        # else:
        #     #MAP assign cells to clone and compute likelihood 
        #     return None
     



        # self.loglikelihood = self.loglikelihood_dict['total']
        # self.variant_likelihood = self.loglikelihood_dict['variant']
        # self.bin_count_likelihood = self.loglikelihood_dict['bin']
        # n= len(self.get_all_cells())
        # m= len(self.get_all_muts())
        # self.norm_loglikelihood = self.loglikelihood/(n*m)

        # return self.loglikelihood

    def compute_likelihood_by_node(self, node, data):

        like0 = data.like0
        if self.cna_genotype_mode:
            like1 = data.like1_dict 
            
        else:
            like1 = data.like1_marg
      
       
        node_like_dict = {}
        node_likelihood = self.compute_variant_likelihood(
            node, like0, like1, None)
        node_like_dict["variant"] = node_likelihood
        node_like_dict["total"] = node_likelihood

        if self.cna_genotype_mode:
            cnn_hmm = data.cna_hmm
            
            bin_node_likelihood, _ = self.rdr_likelihood_by_node(node, cnn_hmm)
            node_like_dict["bin"] = bin_node_likelihood
            node_like_dict['total'] += bin_node_likelihood
        else:
            node_like_dict["bin"] =0
            if self.use_rd:
                rd_node_likelihood, _ = self.read_depth_likelihood_by_node(node, data.read_depth)
                node_like_dict["bin"] = rd_node_likelihood
            node_like_dict['total'] += node_like_dict["bin"]

        return node_like_dict

    def get_loglikelihood(self):
        return self.loglikelihood


    def get_muts(self,node):
        if node not in self.mut_mapping:
            return []
        else:
            return self.mut_mapping[node]
    

    def get_cells(self,node):
        if node not in self.cell_mapping:
            return []
        else:
            return self.cell_mapping[node]
    
    def get_desc_cells(self, node):
        '''
        node: int node in the tree

        returns a list of cell indices that are either assigned to that node or are descendents of given node
        '''

        cells = []
        path = nx.dfs_preorder_nodes(self.tree, node)
        for p in path:
            cells += self.get_cells(p)
        return cells 

    # def relabel(self, label_map):
    #     self.tree = nx.relabel(label_map)
    #     for old, new in label_map.items():
    #         if new in self.ce
    #             raise ValueError("proposed new node label already exists")
    #         cells = self.cell_mapping[old]


    def edges(self):
        return list(self.tree.edges)
    # def compute_variant_likelihood_by_node_without_events(self, node, like0, like1, bin_mapping=None):

    #     m = like0.shape[1]
    #     cells = self.get_tip_cells(node)
    #     like0 = like0[cells, :]
    #     like1 = like1[cells, :]

    #     y = self.presence_by_node(m, cells, node)     

    #     loglikelihood = np.multiply(
    #         (1-y), like0).sum() + np.multiply(y, like1).sum()

    #     return loglikelihood




    def find_root(self):
        for n in self.tree:
            if self.tree.in_degree(n) == 0:
                return n

    def get_ancestral_muts(self, node):
       
        path = list(nx.shortest_path(self.tree, self.root, node))

     
        present_muts =list(chain.from_iterable([self.mut_mapping[p] for p in path if p in self.mut_mapping]))
       
        lost_muts = list(chain.from_iterable([self.mut_loss_mapping[p]
                        for p in path if p in self.mut_loss_mapping]))
        present_muts = list(set(present_muts) - set(lost_muts))
        # if len(lost_muts) > 0:
        #     lost_muts = list(chain.from_iterable((lost_muts)))
   
        return present_muts
    @staticmethod
    def mapping_to_dataframe(mapping, id_name):
  
        pred_list = []
        if len(mapping)==0:
            return pd.DataFrame(columns=[id_name, "cluster"])
        for k in mapping:
                temp = pd.DataFrame(mapping[k], columns=[id_name])
                temp['cluster'] =k
                pred_list.append(temp)

        pred = pd.concat(pred_list)
        pred= pred.sort_values(by=[id_name])
        

        return pred
    
    def generate_mut_dataframe( self,lookup):
    #convert mapping to series in order of mutations
        pred_df= self.mapping_to_dataframe(self.mut_mapping, "mutation_id")

        pred_df["mutation"] = lookup[pred_df['mutation_id']].values

        pred_df = pred_df.drop(['mutation_id'], axis=1)
        return pred_df
    
    def generate_cell_dataframe( self,lookup):
    #convert mapping to series in order of cells
        pred_df= self.mapping_to_dataframe(self.cell_mapping, "cell_id")

        pred_df["cell"] = lookup[pred_df['cell_id']].values

        pred_df = pred_df.drop(['cell_id'], axis=1)
        return pred_df 
    

    def generate_results(self, cell_lookup, mut_lookup):
        pcell = self.generate_cell_dataframe(cell_lookup)
        pmut = self.generate_mut_dataframe(mut_lookup)
        # ploss = generate_mut_dataframe(self.mut_loss_mapping, mut_lookup)


        return pcell, pmut
    
    def save(self, path):
        with open(path, "wb") as file:
              pickle.dump(self, file)

    

    def save_text(self, path):
        
        
        leafs = [n for n in self.tree.nodes if len(list(self.tree.successors(n))) ==0]
      
                    
        with open(path, "w+") as file:
            file.write(f"{len(list(self.tree.edges))} #edges\n")
            for u,v in list(self.tree.edges):
                file.write(f"{u} {v}\n")
            file.write(f"{len(leafs)} #leaves\n")
            for l in leafs:
                file.write(f"{l}\n")
            
def load(fname):
    with open(fname, "rb") as file:
        ct = pickle.load(file)
    return ct 
    # def save_results(self, cell_lookup, mut_lookup, pcell_fname, pmut_fname, ploss_fname, pevents_fname):
    #     pcell, pmut, ploss, pevents = self.generate_results(
    #         cell_lookup, mut_lookup)
    #     pcell.to_csv(pcell_fname, index=False)
    #     pmut.to_csv(pmut_fname, index=False)
    #     ploss.to_csv(ploss_fname, index=False)
    #     pevents.to_csv(pevents_fname)
# class SegmentTree(ClonalTree):
#     def __init__(self, tree, cell_mapping, mut_mapping, cna_genotypes, mut_loss_mapping=None, key=0):
#         super().__init__(tree, {}, {}, {}, {}, key, type)
    
    


# @dataclass
# class Clone:

#     cells: np.array
#     muts: np.array
#     id: int= None
#     cna_genotype: list=None


     
#     def __str__(self):

#         outstring = f"Cells: {len(self.cells)} Muts: {len(self.muts)}" # Ancestral Muts: {len(self.ancestral_muts)} "
#         return outstring

#     def __eq__(self, object):

#         # ancestral_muts_same = np.array_equal(
#         #     np.sort(self.ancestral_muts), np.sort(object.ancestral_muts))

#         if type(object) is type(self):
#             return np.array_equal(self.cells, object.cells) \
#                 and np.array_equal(self.muts, object.muts) #\
#         #         and ancestral_muts_same
#         else:
#             return False

#     def set_id(self, id):
#         self.id = id
    
#     def get_id(self):
#         return self.id


    # def strip(self, var):
    #     var_counts_by_snv= var[np.ix_(self.cells, self.muts)].sum(axis=0)
    #     bad_snvs = self.muts[var_counts_by_snv==0]
    #     self.muts = np.setdiff1d(self.muts, bad_snvs)
        
    #     var_counts_by_cells = var[np.ix_(self.cells,self.muts)].sum(axis=1)
    #     bad_cells = self.cells[var_counts_by_cells ==0]
    #     self.cells = np.setdiff1d(self.cells, bad_cells)

    # def count_obs(self,total):
    #     nobs =np.count_nonzero(total[np.ix_(self.cells,self.muts)])
    #     return nobs
    

