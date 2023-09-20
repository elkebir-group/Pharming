
import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter, defaultdict 
from itertools import product, chain, combinations
import pickle 
import pygraphviz as pgv
# import matplotlib.pyplot as pl
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from sklearn.metrics.cluster import adjusted_rand_score



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
    '#ffed6f']


            
def load(fname):
    with open(fname, "rb") as file:
        ct = pickle.load(file)
    return ct 

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
    
    genotypes: dict of dict of dict of genotype datastructs
        genotypes[v][s][m]
    
    cell_mapping : dict
        a dictionary with node ids as keys and np.arrays with the indices of cells attached to that node
    
    mut_mapping : dict
        a dictionary with node ids as keys and np.arrays with the indices of gained SNVs on the incoming edge
    
    
    mut_loss_mapping : dict 
        a dictionary with node ids as keys and np.arrays with the indices of lost SNVs assigned to the 
        incoming edge of that node

    score : float
    
    




    Methods
    -------


    set_key(key)
        updates the key identifier for the clonal tree

    has_loss()
        returns a boolean if any SNVs are lost


    save(path)
        a helper method to save a clonal tree to disk via pickle

    
    get_leaves()
        returns a list with the node ids of all leaf nodes
    

    
    generate_results(cell_lookup, mut_lookup)
        converts all the internal mapping dictionaries to output dataframes 
    


    find_root()
        returns the node id of the root of the clonal tree
    
    get_all_cells()
        returns an np.array of all cell indices currently assigned in the clonal tree
    
    get_all_muts()
        returns an np.array of all SNV indices currently assigned in the clonal tree



    Genotypes[v][s][m] = genotype
    """

    def __init__(self, tree, genotypes: dict, cell_mapping:dict=None, key=0, cost=np.Inf ):
        self.tree: nx.DiGraph = tree
        self.genotypes = genotypes 
        
        self.root= self.find_root()
        

        if cell_mapping is None:
            self.cell_mapping = {}
            
        else:    
            self.cell_mapping = cell_mapping

     
        self.mut2seg = {}
        for v in self.genotypes:
            for s in self.genotypes[v]:
                for m in self.genotypes[v][s]:
                    self.mut2seg[m]=s
        # self.mut2seg = {m: s  for v in self.genotypes for s, muts in self.genotypes[v].items() for m in muts}
        
        self.seg2muts = {}
        for m,s in self.mut2seg.items():
            if s in self.seg2muts:
                self.seg2muts[s].append(m)
            else:
                self.seg2muts[s] = [m]


        self.mut_mapping, self.mut_loss_mapping = self.get_mut_mapping()
        self.psi = self.get_psi()
        self.phi = self.get_phi()
        self.key = key
        self.cost = cost

    
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
            # genotype = self.tree.nodes[n]['genotype']
            genotype = ""
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
        mystr += f"\nJ(T,phi): {self.cost}\n"
        return mystr
    
    def find_root(self):
        for n in self.tree:
            if self.tree.in_degree(n) == 0:
                return n
    def get_cost(self):
        return self.cost 

    def set_cost(self, cost):
        self.cost  = cost 

    def edges(self):
        return list(self.tree.edges)
    
    def parent(self, v):
        if v ==self.root:
            return None 
        else:
            return list(self.tree.predecessors(v))[0]
    
    def children(self, v):
        return list(self.tree.neighbors[v])
    
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
    
    def get_phi(self):
         self.phi = {i : k  for k, cells in self.cell_mapping.items() for i in cells}
         return self.phi
    
    def get_psi(self):
         self.psi = {m : k  for k, snvs in self.mut_mapping.items()  for m in snvs}
         return self.psi

    def phi_to_cell_mapping(self, phi):
        self.phi = phi 
        self.cell_mapping = {v: [] for v in self.tree}
        for i, v in self.phi.items():
            self.cell_mapping[v].append(i)
      
    
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
    
    def get_cell_assigment(self, i):
        for node in self.cell_mapping:
            if i in self.cell_mapping[node]:
                return node

    def get_snv_assignment(self, m):
        for node in self.mut_mapping:
            if m in self.mut_mapping[node]:
                return node  

    def get_clone_snvs(self,n):
        clone = nx.shortest_path(self.tree, source=self.root, target=n)
        snvs = []
        for c in clone:
            snvs += self.mut_mapping[c]
        return snvs
        
    def set_key(self, key):
        self.key = key
    
    def get_key(self):
        return self.key 

    def has_loss(self):
        return len(self.mut_loss_mapping) > 0
    
    def set_cell_mapping(self, cell_mapping):
        self.cell_mapping  = cell_mapping
        # self.psi = self.get_psi()
        self.phi = self.get_phi()
        # print("Warning, cost is not automatically updated, ensure cost is manually set.")


    def prune_tree(self, node_to_remove):
        parent = next(self.tree.predecessors(node_to_remove))

        # Get the children of the node to be removed
        children = list(self.tree.successors(node_to_remove))

        # Remove the node to be removed from the graph
        self.tree.remove_node(node_to_remove)

        # Reattach the children to the parent
        for child in children:
            self.tree.add_edge(parent, child)

        if node_to_remove in self.cell_mapping:
            del self.cell_mapping[node_to_remove]
        
        if node_to_remove in self.genotypes:
            del self.genotypes[node_to_remove]

    def get_latent_vafs(self, v, s=None):
        vafs = {}
        if s is None:
            segs = list(self.seg2muts.keys())
        else:
            segs  = [s]
        for s in segs:
            for m in self.seg2muts[s]:
                vafs[m] =self.genotypes[v][s][m].vaf
                
        
        return vafs 

    def get_mut_mapping(self):
        gained= []
        lost = []
        muts = []
        mut_mapping = {v: [] for v in self.tree}
        mut_loss_mapping = {v: [] for v in self.tree}
        for v in self.tree:
            for s in self.seg2muts:
                for m, geno in self.genotypes[v][s].items():
                    if m not in muts:
                        muts.append(m)
                    # if m in [523, 792, 451,831]:
                    #     print(m)
                    if geno.z > 0 and m not in gained:
                        mut_mapping[v].append(m)
                        gained.append(m)
                    if m in gained and geno.z ==0 and m not in lost:
                        mut_loss_mapping[v].append(m)
                        lost.append(m)
        missing= set(muts)- (set(gained).union(set(lost)))
        # print(missing)
        for m in missing:
            mut_mapping[self.root].append(m)
        return mut_mapping, mut_loss_mapping 
   
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



    def get_ancestral_muts(self, node):
       
        path = list(nx.shortest_path(self.tree, self.root, node))
        present_muts =list(chain.from_iterable([self.mut_mapping[p] for p in path if p in self.mut_mapping]))
        lost_muts = list(chain.from_iterable([self.mut_loss_mapping[p]
                        for p in path if p in self.mut_loss_mapping]))
        present_muts = list(set(present_muts) - set(lost_muts))

        return present_muts
    

    def node_snv_cost(self, v, cells, data):

        snvs = self.get_all_muts()
        latent_vafs = self.get_latent_vafs(v)
        lvafs = np.array(list(latent_vafs.values())).reshape(1,-1)
        snvs = list(latent_vafs.keys())
        #dims = cells by snvs
        obs_vafs = data.obs_vafs(cells, snvs)

        cell_scores = np.nansum(np.abs(obs_vafs - lvafs), axis=1)


        return cell_scores
    
    def pooled_node_snv_cost(self, v, cells, data):

 
        latent_vafs = self.get_latent_vafs(v)
        lvafs = np.array(list(latent_vafs.values())).reshape(1,-1)
        snvs = list(latent_vafs.keys())
        #dims = cells by snvs
        obs_vafs = data.compute_vafs(cells, snvs)

        cell_scores = np.nansum(np.abs(obs_vafs - lvafs))


        return cell_scores

    
    def assign_cells(self, data):

        nodes = np.array(self.tree)
        cell_scores = np.vstack([self.node_snv_cost(v, data.cells, data) for v in nodes])
        assignments = np.argmin(cell_scores, axis=0)
        
        self.phi = {i: v for i,v in zip(data.cells, nodes[assignments])}
        self.phi_to_cell_mapping(self.phi)
        # self.cell_mapping = defaultdict(list)
        # for i, v in self.phi.items():
        #     self.cell_mapping[v].append(i)
        # self.cell_mapping = dict(self.cell_mapping)



    def compute_costs(self, data, lamb=0):
        if len(self.cell_mapping) ==0:
            self.assign_cells(data)

        self.node_cost = {}
        self.cost = 0
        for v in self.tree:
                cells = self.cell_mapping[v]
                if len(cells) ==0:
                    continue
                cell_scores = self.node_snv_cost(v, cells, data)
                self.node_cost[v] = cell_scores.sum()
        
        self.cost = sum([score for _, score in self.node_cost.items()])

        return self.cost 
    
    def compute_pooled_costs(self, data, lamb=0):
        if len(self.cell_mapping) ==0:
            self.assign_cells(data)
        self.node_cost = {}
        self.cost = 0
        for v in self.tree:
                cells = self.cell_mapping[v]
                if len(cells) ==0:
                    continue
                cell_scores = self.pooled_node_snv_cost(v, cells, data)
                self.node_cost[v] = cell_scores.sum()
        
        self.cost = sum([score for v, score in self.node_cost.items()])

        return self.cost 




   #-------------------------- Save Methods ---------------------------------------#
    def draw(self, fname, mapping=None, cmap='Set3'):

        mut_count = {n : len(self.mut_mapping[n]) for n in self.mut_mapping}
        cell_count = {n : len(self.cell_mapping[n]) for n in self.cell_mapping}
        labels = {}
        # color_values, colormap = self.create_color_map(cmap)
        for n in self.tree:
                if mapping is not None:
                    labels[n] = str(mapping[n])
                else:
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
        if self.cost is not None:
            total_like = np.round(self.cost)
            tree.graph_attr["label"] = f"Objective: {total_like}"
 
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
        '''
        writes the tree to text file
        '''
        
        leafs = [n for n in self.tree.nodes if len(list(self.tree.successors(n))) ==0]
      
                    
        with open(path, "w+") as file:
            file.write(f"{len(list(self.tree.edges))} #edges\n")
            for u,v in list(self.tree.edges):
                file.write(f"{u} {v}\n")
            file.write(f"{len(leafs)} #leaves\n")
            for l in leafs:
                file.write(f"{l}\n")

    def save(self, path):
        with open(path, "wb") as file:
              pickle.dump(self, file)
    
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

    #--------------------- End Save Functions ------------------------------------#
  

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

    def get_cell_assignments(self):
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

    def get_mut_assignments(self, n_mut=None):
        '''

        '''
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
    #----------------------------- Performance Metric Functions ------------------------------------#
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






    

    

    

    
