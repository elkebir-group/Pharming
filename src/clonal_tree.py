
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

        self.mutated_copies = mutated_copies

        if mut_loss_mapping is None:
            self.mut_loss_mapping = {}
        else:
            self.mut_loss_mapping = mut_loss_mapping


        self.key = key
        self.loglikelihood = None

        # self.node_attrs = {k for node_data in self.tree.nodes.data() for k in node_data[1].keys()}
        # self.edge_attrs =  {k for edge_data in self.tree.edges.data() for k in edge_data[1].keys()}
    def get_clone_cells(self, n):
        if n in self.cell_mapping:
            return np.concatenate([self.cell_mapping[t][k] for k in self.cell_mapping[t]]).tolist()
        else:
            return []

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
                # if "event" in self.edge_attrs:
                etype = self.tree[u][v]["event"]
           
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

    def create_color_map(self, cmap):
        genos = []
        for n in self.tree:
            geno = self.tree.nodes[n]["genotype"]
            genos.append(sum(geno))
     
        color_values = set(genos)
        colormap = cm.get_cmap(cmap, len(color_values))
        return list(color_values), colormap

# Define the colormap

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
 
        colormap = cm.get_cmap(cmap)
        for n in self.tree:

            tree.add_node(n, label=labels[n])
      
            node_attr = tree.get_node(n)
            try:
                x,y = self.tree.nodes[n]["genotype"]
                color_value = x+y
            except:
                color_value = None
        
            if color_value is not None:
               
                
                color = colormap(color_value)
                hex_color = mcolors.rgb2hex(color)
                node_attr.attr['fillcolor'] =hex_color
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
            if len(cells[0]) > 0:
                clusters[cells[0]] = cluster
        return clusters

    def get_all_muts(self):
        muts = np.concatenate([self.mut_mapping[i] for i in self.mut_mapping])
        muts = np.sort(muts)
        return muts

    def get_mut_clusters(self, n_mut=None):
        if n_mut is None:
            n_mut = len(self.get_all_muts())
            clusters = np.zeros(n_mut, dtype=int)
        else:
            clusters  = np.full(shape=n_mut, fill_value=-1)
        for cluster, muts in self.mut_mapping.items():
            if len(muts) > 0:
                clusters[muts] = cluster
        return clusters










    def get_cell_cluster(self, c):
        for node in self.cell_mapping:
            if c in self.cell_mapping[node][0]:
                return node

    def get_snv_cluster(self, s):
        for node in self.mut_mapping:
            if s in self.mut_mapping[node]:
                return node  



    def get_cell_ancestor_pairs(self) -> Counter:
        pairs = Counter()
        for node in self.tree.nodes:
            for children in nx.dfs_successors(self.tree, source=node).values():
                for child in children:
                    pairs.update(product(self.cell_mapping[node][0], self.cell_mapping[child][0]))
        return pairs

    def get_cell_cluster_pairs(self) -> Counter:
        pairs = Counter()
        for node in self.tree.nodes:
            pairs.update(combinations(sorted(self.cell_mapping[node][0]), 2))
        return pairs

    def get_cell_incomparable_pairs(self) -> Counter:
        pairs = Counter()
        for u, v in combinations(self.tree.nodes, 2):
            if self.is_incomparable(self.tree, u, v):
                pairs.update((min(a, b), max(a, b)) for a, b in product(self.cell_mapping[u][0], self.cell_mapping[v][0]))
                # for cell1 in self.cell_mapping[u][0]:
                #     for cell2 in self.cell_mapping[v][0]:
                #         if cell1 < cell2:
                #             pairs[(cell1, cell2)] += 1
                #         else:
                #             pairs[(cell2, cell1)] += 1
        return pairs

    def get_ancestor_pairs(self, include_loss: bool=True) -> Counter:
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
                        pairs.update(product(self.mut_mapping[node], self.mut_mapping[child]))
        return pairs

        return pairs

    def compute_likelihood(self, data, seg, alpha=0.001, attachment="marginal", use_existing=False):

        seg_snvs = data.seg_to_snvs[seg]
      
        # self.loglikelihood_dict = {"total": 0, "variant": 0, "bin": 0}
        self.node_likelihood = {}
        if attachment =="marginal":
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
         
            self.loglikelihood = logsumexp(all_nodes, axis=0).sum()
            
            return self.loglikelihood

         


                


                # for key in self.loglikelihood_dict:
                #  self.loglikelihood_dict[key] += self.node_likelihood[n][key]

    
        elif attachment == 'map' and use_existing:
            for n in self.cell_mapping:
                if len(self.get_tip_cells(n)) > 0:
                    self.node_likelihood[n] = self.compute_likelihood_by_node(n, data)
        else:
            #MAP assign cells to clone and compute likelihood 
            return None
     



        self.loglikelihood = self.loglikelihood_dict['total']
        self.variant_likelihood = self.loglikelihood_dict['variant']
        self.bin_count_likelihood = self.loglikelihood_dict['bin']
        n= len(self.get_all_cells())
        m= len(self.get_all_muts())
        self.norm_loglikelihood = self.loglikelihood/(n*m)

        return self.loglikelihood

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
    

