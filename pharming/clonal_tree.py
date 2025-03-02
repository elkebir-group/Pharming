
import numpy as np
import pandas as pd
import networkx as nx
from itertools import  chain
import pickle 
import pygraphviz as pgv
from .cell_mapping import CellAssign
from scipy.stats import binom
from collections import defaultdict
from scipy.special import logsumexp


EPSILON = -1e5

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
    genotypes: dict (node ids as keys) of dict (snv index as keys) of tuples (x,y,w,z)
        genotypes[v][m]
    cell_mapping : dict of lists
        a dictionary with node ids as keys and np.arrays with lists of cells attached to that node
    mut_mapping : dict of lists
        a dictionary with node ids as keys and lists with the indices of gained SNVs on the incoming edge
    mut_loss_mapping : dict of lists 
        a dictionary with node ids as keys and lists with the indices of lost SNVs assigned to the 
        incoming edge of that node
    cost : float current cost of cell assignment
    root : int id of root node 
    

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



    """

    def __init__(self, tree, genotypes: dict, seg_to_muts:dict=None,  rho:dict=None, k=None):
        self.tree: nx.DiGraph = tree
        self.root= self.find_root()
        # self.preorder_nodes = list(nx.dfs_preorder_nodes(self.tree, self.root))
        self.genotypes = genotypes 
        if seg_to_muts is None:
            self.seg_to_muts = {}
        else:
            self.seg_to_muts = seg_to_muts
        
        self.mut_to_seg = {}
        for seg, muts in self.seg_to_muts.items():
            for m in muts:
                self.mut_to_seg[m] = seg
        
      


        self.mut_mapping, self.mut_loss_mapping = None, None

        self.psi = {}

  
        self.cost = np.inf
        self.snv_cost = np.inf
        self.cna_cost = np.inf

        if rho is not None:
            for ell in rho:
                for q in rho[ell]:
                    rho[ell][q] = sorted(rho[ell][q], key=lambda t: self.snv_tree_has_loss(t))
                
        self.rho = rho
        if k is None and self.rho is not None:
            for ell in self.rho:
                self.k = len(self.rho[ell])
                break 
        else:
            self.k = None 
    

    
    def __str__(self) -> str:
        # all_cells = self.get_all_cells()
        all_snvs  = self.get_all_muts()
        mystr =f"\nClonal Tree T"
        mystr += "\n--------------------------------"

        mystr +=f"\n{len(list(self.tree.nodes))} Nodes and {len(list(self.tree.edges))} Edges"
        mystr += f"\n{len(all_snvs)} SNVs"
        return mystr
        # mystr+= f"\n{len(all_cells)} cells\n"


        # node_order = list(nx.dfs_preorder_nodes(self.tree, source=self.root))
        # mystr += "\nParent->Child\tEdge Type"
        # mystr += "\n--------------------------------"
        # for u in node_order:
        
        #     for v in self.tree.neighbors(u):
        #         etype= ""
        #         # if "event" in self.tree.edge_attrs:
        #         #     etype = self.tree[u][v]["event"]
        #         # else:
        #         #     etype = "NA"
           
        #         mystr+= f"\nNode {u}->Node {v}\t{etype}"
               
        # mystr += "\n\nNode\tSNVs\tGenotype"
        # mystr += "\n--------------------------------"
        # for n in  node_order:
        #     # genotype = self.tree.nodes[n]['genotype']
        #     genotype = ""
        #     # if n in self.cell_mapping:
        #     #     ncells = len(self.cell_mapping[n])
        #     # else:
        #     #     ncells = 0
        #     if n in self.mut_mapping:
        #         nmuts = len(self.mut_mapping[n])
        #     else:
        #         nmuts = 0
        #     mystr += f"\n{n}\t{nmuts}\t{genotype}"
        # mystr += "\n"
        # mystr += f"\nJ(T,phi): {self.cost}\n"
        # return mystr
    
    def check_genotypes(self):
        return all(len(self.genotypes[v]) >0 for v in self.tree)
    
    def preorder(self, source=None):

        # if source is None:
        #     return self.preorder_nodes

        if source is None:
            return list(nx.dfs_preorder_nodes(self.tree, source=self.root))
   
    
        return list(nx.dfs_preorder_nodes(self.tree, source=source))

    def get_tree(self):
        return self.tree.copy()

    def get_genotypes(self):
        return self.genotypes

    def get_seg_to_muts(self):
        return self.seg_to_muts.copy()
    
    def get_segments(self):
        return set(self.seg_to_muts.keys())
    
    def compute_dcfs(self, ca):
   
        self.update_mappings()

        dcfs = {}
        counts = ca.get_cell_count()
        ncells = ca.n
        for n in self.tree:
            
            if len(self.mut_mapping[n]) > 0:
          
                dcfs[n] = sum([counts[u] for u in self.preorder(source=n)])/ncells 
        return dcfs 
    


    def posterior_dcf(self, j:int, dcf:float,a: int,d:int ,cn_prop:dict ):

        # if j == 290:
        #     print(f"{j}: {dcf}")
        if d==0:
            posterior = EPSILON
        else:
            parent, u, p_geno, u_geno = self.get_split_nodes(j)
            m_star = u_geno[2] + u_geno[3]
            #compute fractional copy number as weighted sum 
            F = sum([ (cn[0]+ cn[1])*cn_prop[cn] for cn in cn_prop])
            desc_geno = self.get_desc_genotypes(u, j)
            # cn_prop = {x+y: 1.0 if x+y==cn else 0.0 for x,y,m in self.gamma}
            #v = self.dcf(dcf, cn)
            #copy code so that posterior 
            v = (dcf*m_star)/F 
            for geno in desc_geno:
                z = geno[2] + geno[3]
                cn_state = (geno[0], geno[1])

                v += (1/F)* (z-m_star)*cn_prop[cn_state]
                # + (1/F)*sum([(geno.w-m_star)*cn_prop[cn_state] for x,y,m in self.desc_genotypes])
    
            # posterior = max(beta.logpdf(v, a+1, d-a + 1),EPSILON)
            posterior = max(binom.logpmf(a, d,v), EPSILON)
            # if posterior > 0:
            #     print(f"{posterior}: {dcf} : {v} : {a/d}")
        return posterior



    def relabel(self, mapping=None):
        """
        dict mapping: dictionary containing old labels as keys, new labels as values
                        may be partial

        Updates the nodes labels. If mapping is None, the nodes are labeled sequential
        via a preorder traversal
    
        """
        if mapping is None:
            mapping = {n: i for i,n in enumerate(self.preorder())}
        
        self.tree = nx.relabel_nodes(self.tree, mapping)
        genotypes = {}
        #TODO: check if new mapping will result in a collison 
        for u in mapping:
 
            genotypes[mapping[u]] = self.genotypes[u]
        
        self.genotypes = genotypes
        self.update_mappings()
        
        
    def has_loss(self):
        self.update_mappings()
        return any(len(val) for node, val in self.mut_loss_mapping.items())
 
    def mutation_cluster_tree(self):
        if self.mut_mapping is None:
            self.update_mappings()

        cluster_nodes = [n for n in self.tree if len(self.mut_mapping[n]) > 0]
        # for key, val in self.rho.items():
        #     cluster_nodes = list(val.keys())
        #     break


        def find_closest_ancestor(u, nodes_list):
        # Initialize variables to keep track of the closest ancestor and its distance
            closest_ancestor = None
            min_distance =  np.inf  # Initialize to infinity
            
            # Calculate the shortest path from u to each node in the list
            for node in nodes_list:
                # Calculate the shortest path length from u to the current node
                if u !=node and nx.has_path(self.tree, source=node, target=u ):
                    distance = nx.shortest_path_length(self.tree, source=node, target=u)
                
                    # Update the closest ancestor and its distance if the current node is closer
                    if distance < min_distance:
                        closest_ancestor = node
                        min_distance = distance
            return closest_ancestor
    

        T_m = nx.DiGraph()
        for n in cluster_nodes:
            T_m.add_node(n)
            parent = find_closest_ancestor(n, cluster_nodes)
            if parent is not None:
                T_m.add_edge(parent, n)
        return T_m
      


    
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
    
    # def get_cna_genos(self):
    #     self.cna_genotypes = {}
        

    #     for seg, snvs in self.seg_to_muts.items():
    #         m = snvs[0]
    #         self.cna_genotypes[seg] = {}
    #         for v in self.tree:
    #             g = self.genotypes[v][m]
    #             self.cna_genotypes[seg][v] = (g[0],g[1])
         
    #     return self.cna_genotypes
    
    def get_cna_genos(self):
        self.cna_genotypes = {}

        # Pre-calculate the mutation for each segment
        # segment_mutation = {seg: snvs[0] for seg, snvs in self.seg_to_muts.items() if len(snvs) >0}

        for seg,snvs in self.seg_to_muts.items():
            if len(snvs) > 0:
                m =snvs[0]
                self.cna_genotypes[seg] = {}
                for v in self.tree:
                    g = self.genotypes[v][m]
                    self.cna_genotypes[seg][v] = (g[0], g[1])

        return self.cna_genotypes

    def get_cna_tree(self, seg):
        self.get_cna_genos()
        S = nx.DiGraph()

        for u,v in self.tree.edges:
            u_cna = self.cna_genotypes[seg][u] 
            v_cna = self.cna_genotypes[seg][v]
            if u_cna != v_cna:
                S.add_edge(u_cna, v_cna)
        if len(S) ==0:
            S.add_node(self.cna_genotypes[seg][self.root])
        return S
    
    def get_non_infinite_alleles_cna_tree(self, seg):
        S = nx.DiGraph()
        state_dict ={}
        for u,v in self.tree.edges:
            u_cna = self.cna_genotypes[seg][u] 
            v_cna = self.cna_genotypes[seg][v]
            if u_cna != v_cna:
                S.add_edge(u, v)
                if u not in state_dict:
                    state_dict[u] = u_cna
                if v not in state_dict:
                    state_dict[v] = v_cna
        return S, state_dict
                

            
            

                    

    
    def parent(self, v):
        if v ==self.root:
            return None 
        else:
            return list(self.tree.predecessors(v))[0]
    
    def children(self, v):
        return list(self.tree.neighbors(v))
    
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
    
    def get_desc_genotypes(self, v, j):
        
        #TODO: fix convert to a set so that it works on clonal trees and not just genotype trees
        return [self.genotypes[u][j] for u in nx.dfs_preorder_nodes(self.tree, v) if u != v]
    
    @staticmethod
    def cna_eq(g1, g2):
        return (g1[0]== g2[0] and g1[1]==g2[1])

    @staticmethod
    def mut_copies(g):
        return g[2] + g[3]
    
    @staticmethod
    def vaf(g):
        return (g[2] + g[3]) / (g[0] + g[1])

    def get_split_nodes(self, j):
        
        for u in self.preorder():
            if u != self.root:
                parent= self.parent(u)
                p_geno =self.genotypes[parent][j]
                u_geno = self.genotypes[u][j]
                u_z = self.mut_copies(u_geno)
                if self.cna_eq(p_geno, u_geno) and u_z > 0:
                    return parent, u, p_geno, u_geno

    # def trim(self, cellAssign):
    #     '''
    #     remove nodes that only have cell assignments and pushees the cell assignments up the tree
    #     '''
        
    #     if self.mut_mapping is None or self.mut_loss_mapping is None:
    #         self.update_mappings()
    #     cell_mapping = cellAssign.cell_mapping
    #     to_del = []
    #     for k in cell_mapping:
    #         if len(self.mut_mapping[k])==0 and len(self.mut_loss_mapping[k])==0:
    #             if len(cell_mapping[k]) > 0 and self.is_leaf(k):
    #                 path =nx.shortest_path(self.tree, source=self.root, target=k)
    #                 path = path[::-1]
    #                 for p in path:
    #                     if len(self.mut_mapping[p]) >0 or len(self.mut_loss_mapping[p]) > 0:
    #                         cell_mapping[p] += cell_mapping[k] 
    #                         to_del.append(k)
    #                         self.tree.remove_node(k)
    #                         break
    #     for k in  to_del:
    #         del cell_mapping[k]
    #     cellAssign.update_phi()
    #     cellAssign.update_clones(self.clones())
    
    def loss_recall(self, snvs: list) -> float:
        """Calculate the SNV loss recall for a given set of ground truth lost SNVs.
        
        Parameters
        ----------
        snvs : list
            A list of SNVs that are the ground truth lost SNVs
        
        Returns
        -------
        float
            The recall of the lost SNVs
        """
        self.update_mappings()
        TP = 0
        FN = 0
        lost_snvs = self.get_lost_snvs()
        for j in snvs:
            if j in lost_snvs:
                TP += 1
            else:
                FN += 1
        
        if TP + FN == 0:
            return 0.0  # Handle division by zero if no positives are found
        
        return TP / (TP + FN)


    def get_lost_snvs(self):
        self.update_mappings()
        return list(chain.from_iterable(self.mut_loss_mapping.values()))
    
    def loss_precision(self, snvs: list) -> float:
        """Calculate the SNV loss precision for a given set of ground truth lost SNVs.
        
        Parameters
        ----------
        snvs : list
            A list of SNVs that are the ground truth lost SNVs
        
        Returns
        -------
        float
            The precision of the lost SNVs
        
        """
        self.update_mappings()
        TP = 0
        FP = 0
        lost_snvs = self.get_lost_snvs()
        for j in lost_snvs:
            if j in snvs:
                TP += 1
            else:
                FP += 1  # This counts False Positives
        
        if TP + FP == 0:
            return 0.0  # Handle division by zero if no positives are found
        
        return TP / (TP + FP)


    def fp_leaves(self, snvs:list) -> dict:
        self.update_mappings()
        fps = []
        for u, lost_snvs in self.mut_loss_mapping.items():
            is_leaf = self.is_leaf(u)
            for j in lost_snvs:
                fps.append([u, j, is_leaf, j in snvs])
    
        return pd.DataFrame(fps, columns=["node", "snv", "is_leaf", "TP"])


    # def toJson(self, fname, segment_csv=None, snv_csv=None):
    #     json = convertToJson(self, segment_csv, snv_csv)
    #     with open(fname, "w") as f:
    #         f.write(str(json))
    
    def filter_snvs(self,  snvs):
      

        snvs = set(snvs)
        # snvs = set(self.seg_to_muts[seg])
        keys = set(self.genotypes[self.root].keys())
        to_del  = keys - snvs 


        for v in self.genotypes:
            # self.mut_mapping[v] = [j for j in self.mut_mapping[v] if j in snvs]

            for s in to_del:
                del self.genotypes[v][s]
        
        for s in to_del:
            seg = self.mut_to_seg[s]
            self.seg_to_muts[seg].remove(s)
            del self.mut_to_seg[s]
        
    def filter_segments(self,segments ):

        snvs = [j for ell in segments for j in self.seg_to_muts[ell] if ell in self.seg_to_muts]
        self.filter_snvs(snvs)
        to_remove = set(self.seg_to_muts.keys()) - set(segments)
        for ell in to_remove:
            del self.seg_to_muts[ell]
            
           
            

    
    def add_snv(self,j, seg, node, snv_tree):
        '''
        update the genotypes, mut_mapping and psi to include snv j
        '''
        if node not in self.tree:
            raise KeyError("Node id is not in the clonal tree.")

        if seg not in self.seg_to_muts:
            raise KeyError("Segment is not in the clonal tree.")

        self.mut_to_seg[j] = seg 
        self.seg_to_muts[seg].append(j)
        # self.mut_mapping[node].append(j)
        self.psi[j] = node 


        def match_snv_genotype(tree, cn_state):
            for u in tree:
                # geno = genotype(*u)
                if self.mut_copies(u) > 0 and u[0] == cn_state[0] and u[1]==cn_state[1]:
                    return u
        
        def update_desc_genotype(u, snv_tree):
            descendants.append(u)
       
            cn_state = self.cna_genotypes[seg][u]
            self.genotypes[u][j] =match_snv_genotype(snv_tree, cn_state)

            for v in self.tree.successors(u):
                update_desc_genotype(v, snv_tree)
        
        descendants = []
        self.get_cna_genos()
        update_desc_genotype(node, snv_tree)

        not_present = set(self.clones())  - set(descendants)
        for u in not_present:
            cn_tup =self.cna_genotypes[seg][u]
            self.genotypes[u][j] =  (*cn_tup, 0,0)
        
    

    def clones(self):
        return set(self.tree.nodes)
    
    def is_leaf(self,node):
        return self.tree.out_degree[node] ==0
    
    def add_rho(self, rho):
        self.rho = self.rho | rho

    
    # def get_psi(self):
    #      self.update_mappings()
      
    #      return self.psi

      
    def get_psi(self):
        psi = {}

        for v in self.preorder():
            if v == self.root:
                continue
            parent = self.parent(v)
            g_par = self.genotypes[parent]
            for j, g in self.genotypes[v].items(): 
                if  self.mut_copies(g_par[j] ) ==0 and self.mut_copies(g) >0 and j not in psi:
                    psi[j] = v 
  
        for j in self.get_all_muts():
            if j not in psi:
                print(f"warning SNV {j} never introduced, appending to root")
                psi[j] = self.root
        self.psi = psi 
        return psi
            


    def get_muts(self,node):
        if self.mut_mapping is None:
            self.update_mappings()
        if node not in self.mut_mapping:
            return []
        else:
            return self.mut_mapping[node]
    

   

    def reindex_snvs(self,mut_lookup):
        new_genotypes = {}
        new_mut_to_seg = {}
        mapping_dict = {}
        for index, mut_label in mut_lookup.items():
            mapping_dict[mut_label] = index 
            new_mut_to_seg[index] = self.mut_to_seg[mut_label]
        
        for v in self.clones():
            new_genotypes[v] = {}
            for m in self.genotypes[v]:
                if m in mapping_dict:
                    new_genotypes[v][mapping_dict[m]] = self.genotypes[v][m]
        
        self.genotypes = new_genotypes
        self.mut_to_seg = new_mut_to_seg
        self.seg_to_muts = {}
        for m, seg in self.mut_to_seg.items():
            if seg in self.seg_to_muts:
                self.seg_to_muts[seg].append(m)
            else:
                self.seg_to_muts[seg] = [m]
        
        self.update_mappings()

    def get_snv_assignment(self, m):
        if self.mut_mapping is None:
            self.update_mappings()

        for node in self.mut_mapping:
            if m in self.mut_mapping[node]:
                return node  

    def get_clone_snvs(self,n):
        if self.mut_mapping is None:
            self.update_mappings()

        clone = nx.shortest_path(self.tree, source=self.root, target=n)
        snvs = []
        for c in clone:
            snvs += self.mut_mapping[c]
        return snvs
        
    def set_key(self, key):
        self.key = key
    
    def get_key(self):
        return self.key 


    def has_snvs(self, n):
        if self.mut_mapping is None:
            self.update_mappings()
        return len(self.mut_mapping[n]) > 0
    

    def collapse(self, phi, k, cell_threshold=10):
    

        while True:
            self.update_mappings()
            collapsed = False
            for u in nx.dfs_postorder_nodes(self.tree, self.root):
                if u == self.root:
                    continue
            

                #check if a candidate should be collapsed
                #must meet cell threshold, not be an SNV cluster and the outdegree of the parent mut be 1
                if len(phi.get_cells(u)) < cell_threshold and u >= k and self.tree.out_degree[u] ==1:
                    child = self.children(u)[0]
                    collapsed = True
                    phi.move_cells(u, child)
                    self.prune_tree(u)
                    phi.update_clones(self.clones())
                    break 
            if not collapsed:
                self.update_mappings()
                break 
        # self.update_mappings()

    #PRIVATE 
    def prune_tree(self, node_to_remove):
        
        parent = self.parent(node_to_remove) #next(self.tree.predecessors(node_to_remove))

        # Get the children of the node to be removed
        children = list(self.tree.successors(node_to_remove))

        # Remove the node
        self.tree.remove_node(node_to_remove)

        # Reattach the children to the parent
        for child in children:
            self.tree.add_edge(parent, child)

        # if node_to_remove in self.cell_mapping:
        #     del cellAssign.cell_mapping[node_to_remove]
        
        if node_to_remove in self.genotypes:
            if len(self.mut_mapping[node_to_remove]) > 0:
                for j in self.mut_mapping[node_to_remove]:
                    ell = self.mut_to_seg[j]
                    del self.mut_to_seg[j]
                    self.seg_to_muts[ell] = [j2 for j2 in self.seg_to_muts[ell] if j != j2]

            # self.genotypes[parent] = self.genotypes[node_to_remove]
            del self.genotypes[node_to_remove]

    def prune_leaves(self, phi, k):
        to_prune = [l for l in self.get_leaves() if len(phi.get_cells(l)) ==0 and l >= k]
        for l in to_prune:

            self.prune_tree(l)
        
            phi.update_clones(self.clones())
        self.update_mappings()



    ##### METHODS TO MODIFY TREE TOPOLOGY
    # def collapse(self, phi, k, cell_threshold=10):
    #     updated = False
    #     cell_counts = phi.get_cell_count()
    #     candidates = [n for n in cell_counts if cell_counts[n] < cell_threshold and 
    #                   n > k and n != self.root and self.tree.out_degree[n] == 1]
    #     chain_dict = {}
    #     for n in self.preorder():
    #         if len(candidates) == 0:
    #             break 
    #         if n in candidates:
    #             chain_dict[n] = []
    #             candidates.remove(n)
       
    #             for u in self.preorder(n):
    #                 if u == n:
    #                     continue
    #                 if u in candidates and self.tree.out_degree[u] ==1:
    #                     chain_dict[n].append(u)
    #                     candidates.remove(u)
    #                 else:
    #                     break

    #     for keep, remove in chain_dict.items():
    #         for u in remove:
    #             updated = True
    #             if cell_counts[u] > 0:
    #                 for i in phi.get_cells(u):
    #                     phi.move(i, keep)
    #             self.prune_tree(u, phi)
    #         if self.tree.out_degree[keep] == 1:
   
    #             if cell_counts[keep] > 0:
    #                 child = list(self.tree.successors(keep))[0]
    #                 for i in phi.get_cells(u):
    #                     phi.move(i, child)
    #             self.prune_tree(keep, phi)

    #     if updated:
    #         self.preorder_nodes = nx.dfs_preorder_nodes(self.tree, self.root)


    # def prune(self, cellAssign):
    #     updated = False
    #     counts = cellAssign.get_cell_count()
    #     self.update_mappings()
    #     for u in self.preorder():
    #         if u ==self.root:
    #             continue
    #         if len(self.children(u)) ==1:    
    #             if counts[u] ==0 and len(self.get_muts(u)) ==0 and len(self.mut_loss_mapping[u]) ==0:
    #                 updated = True
    #                 self.prune_tree(u,cellAssign)
   
    #     if updated:
    #         self.update_mappings()
    #         self.preorder_nodes = nx.dfs_preorder_nodes(self.tree, self.root)
        


             
     
    ######################################
  


    

    
    def update_mappings(self):
        self.psi = self.get_psi()
        self.mut_mapping = {n: [] for n in self.tree}
        self.mut_loss_mapping =  {n: [] for n in self.tree}
        for j,q in self.psi.items():
            self.mut_mapping[q].append(j)
        
        for q,snvs in self.mut_mapping.items():
            snvs = self.mut_mapping[q]
            for u in self.preorder(q):
                genos = self.genotypes[u]
                if q != u:
                    par = self.parent(u)
                    pgenos = self.genotypes[par]
                    for j in snvs:
                        if self.mut_copies(pgenos[j]) >0 and self.mut_copies(genos[j])==0:
                            self.mut_loss_mapping[u].append(j)
                      



        # #TODO: update psi
        # self.mut_mapping, self.mut_loss_mapping =self.get_mut_mapping()

        # for n, snvs in self.mut_mapping.items():
        #     for j in snvs:
        #         self.psi[j] = n

        
   
    # def get_desc_cells(self, node):
    #     '''
    #     node: int node in the tree

    #     returns a list of cell indices that are either assigned to that node or are descendents of given node
    #     '''

    #     cells = []
    #     path = nx.dfs_preorder_nodes(self.tree, node)
    #     for p in path:
    #         cells += self.get_cells(p)
    #     return cells 



    def get_ancestral_muts(self, node):
        if self.mut_mapping is None:
            self.update_mappings()
       
        path = list(nx.shortest_path(self.tree, self.root, node))
        present_muts =list(chain.from_iterable([self.mut_mapping[p] for p in path if p in self.mut_mapping]))
        # lost_muts = list(chain.from_iterable([self.mut_loss_mapping[p]
        #                 for p in path if p in self.mut_loss_mapping]))
        present_muts = set(present_muts) #- set(lost_muts))

        return present_muts
    
  
    
  
    
    def get_latent_genotypes(self,v):
        return {m: geno for m,geno in self.genotypes[v].items()} 
    
    def get_latent_vafs(self, latent_genos):
        return  {m: self.vaf(geno) for m, geno in latent_genos.items()}               
     
    
    def get_latent_cna_genos(self, latent_genos):
        return  {m: (geno[0], geno[1]) for m, geno in latent_genos.items()}    


  

    # def node_cna_cost(self, v, cells, data, latent_genos):
    #     #dictionary is snvs as keys and cna_geno types as values
    #     latent_cna_genos= self.get_latent_cna_genos(latent_genos)
    #     # snv_keys = list(latent_cna_genos.keys())
    #     # snv_to_seg = {s: data.snv_to_seg[s] for s in snv_keys}

    #     seg_geno = {}
    #     for key, seg in self.mut_to_seg.items():
    #         if seg not in seg_geno:
    #            seg_geno[seg] = latent_cna_genos[key]
    #     latent_x, latent_y = [], []
    #     for seg,cna_geno in seg_geno.items():
    #         latent_x.append(cna_geno[0])
    #         latent_y.append(cna_geno[1])
        
    #     latent_x = np.array(latent_x).reshape(1,-1)
    #     latent_y = np.array(latent_y).reshape(1,-1)


    #     segs = list(seg_geno.keys())
    #     obs_copy_x, obs_copy_y =  data.copy_profiles_by_seg(segs, cells)
    #     x_diff = np.abs(obs_copy_x - latent_x).sum(axis=1)
    #     y_diff = np.abs(obs_copy_y - latent_y).sum(axis=1)
    #     return x_diff + y_diff
 
    def node_cna_cost(self, v, cells, data, cna_genos, weights=None):
        # Get latent CNA genotypes for each segment
    

        # Construct latent_x and latent_y arrays
        latent_x = np.array([cna_genos[key][v][0] for key in cna_genos]).reshape(1,-1)
        latent_y = np.array([cna_genos[key][v][1] for key in cna_genos]).reshape(1,-1)

        # Get copy profiles for the segments
        segs = list(cna_genos.keys())
        obs_copy_x, obs_copy_y = data.copy_profiles_by_seg(segs, cells)

        weights = None # data.weights.reshape(1,-1)
        # Compute absolute differences
        x_diff = np.abs(obs_copy_x - latent_x)
        y_diff = np.abs(obs_copy_y - latent_y)
        diff = x_diff + y_diff

        if weights is not None:
            diff = diff*weights

    
        return diff.sum(axis=1)
    # def node_snv_cost(self,  v, cells, data, latent_genos):

    #     snvs = self.get_all_muts()
    #     latent_vafs =self.get_latent_vafs(latent_genos)
    #     lvafs = np.array(list(latent_vafs.values())).reshape(1,-1)
    #     snvs = list(latent_vafs.keys())
    #     #dims = cells by snvs
    #     obs_vafs = data.obs_vafs(cells, snvs)

    #     cell_scores = np.nansum(np.abs(obs_vafs - lvafs), axis=1)


    #     return cell_scores
    
    # def pooled_node_snv_cost(self, v, cells, data, latent_genos):

 
    #     latent_vafs = self.get_latent_vafs(latent_genos)
    #     lvafs = np.array(list(latent_vafs.values()))
    #     snvs = list(latent_vafs.keys())
    #     #dims = cells by snvs
    #     obs_vafs = data.compute_vafs(cells, snvs)

    #     comp_df =  pd.DataFrame({'obs_vafs': obs_vafs, 'latent_vafs': lvafs, 'snvs': snvs})
    #     comp_df["node"] = v
    #     comp_df["ncells"] = len(cells)
    #     self.all_dfs[v]=comp_df
    #     cell_scores = np.nansum(np.abs(obs_vafs - lvafs))


    #     return cell_scores




    # def compute_node_likelihoods(self, data, cells=None, lamb=0):
    #     if cells is None:
    #         cells = data.cells
    #     nodes = np.array(self.tree)
    #     node_cell_likes = []
    #     for u in nodes: 
    #         # vaf = self.get_latent_vafs(self.get_latent_genotypes(u))
    #         # if snvs is None:
    #         #     snvs = list(vaf.keys())
    #         # vafs = np.array([vaf[j] for j in snvs]).reshape(-1,1)
        
    #         node_cell_likes.append(self.node_likelihood(u, data, cells))
        
    #     cell_scores = np.vstack(node_cell_likes)
    #     cell_cna_scores = np.vstack([self.node_cna_cost(v, data.cells, data, self.get_latent_genotypes(v)) for v in nodes])
    #     cell_scores = cell_scores + lamb*cell_cna_scores

    #     return cell_scores, nodes
    
    @staticmethod
    def get_indices_map(vaf):
            snvs = np.array(list(vaf.keys()))
            vafs = np.array([vaf[j] for j in snvs])
            unique_values = list(set(vafs))
   
            # unique_values, indices = np.unique(vafs, return_inverse=True)

            # Create a dictionary to map each unique value to its indices
            indices_map = {}
           
            for vaf in unique_values:
                indices_map[vaf] = snvs[np.where(vafs==vaf)[0]]
            
            return indices_map
    
    def node_likelihood(self, data, cells, vaf):

  
            # snvs = np.array(list(vaf.keys()))
            # vafs = [vaf[j] for j in snvs]
            # unique_values = list(set(vafs))
   
            # # unique_values, indices = np.unique(vafs, return_inverse=True)

            # # Create a dictionary to map each unique value to its indices
            # indices_map = {}
            # vafs = np.array(vafs)
            # for vaf in unique_values:
               
            #     indices_map[vaf] = snvs[np.where(vafs==vaf)[0]]
            indices_map = self.get_indices_map(vaf)
            cell_scores = data.compute_cell_likelihoods(indices_map, cells=cells)
        
        
            return cell_scores.sum()
    
    # @utils.timeit_decorator
    # def test_likelihood(self,  data,  vaf):

  
    #         snvs = list(vaf.keys())
    #         vafs = np.array([vaf[j] for j in snvs]).reshape(-1,1)
        
    #         return data.binomial_likelihood(data.cells, snvs,vafs )
    

    def compute_node_likelihoods(self, data, cells=None, lamb=0):
        #use vectorization to compute all logpmfs at the same time
        #mask values where total=0, use sparse matrix?
        
        if cells is None:
            cells = data.cells
        nodes = np.array(self.tree)
        node_cell_likes = []
        node_cna_scores = []
        cna_genos = self.get_cna_genos()
    
        for u in nodes:
            latent_genos = self.get_latent_genotypes(u)
            latent_vaf = self.get_latent_vafs(latent_genos)
        
            indices_map = self.get_indices_map(latent_vaf)
            cell_scores = data.compute_cell_likelihoods(indices_map)
            node_cell_likes.append(cell_scores)
   

            # node_cell_likes.append(self.node_likelihood(u, data, cells,latent_vaf))
            # self.test_likelihood(data, latent_vaf)
            node_cna_scores.append(self.node_cna_cost(u, cells, data, cna_genos))
        cell_scores = np.vstack(node_cell_likes)
        cell_cna_scores = np.vstack(node_cna_scores)
        cell_scores += lamb * cell_cna_scores
        return cell_scores, nodes


    # def constrained_cell_assignment(self, data, lamb, dcfs):
    #         cell_scores, nodes = self.compute_node_likelihoods(data, lamb=lamb)
    #         obj, ca = ConstrainedCellAssign(self, cell_scores, nodes, dcfs_lb=dcfs).solve()
    #         return obj,ca 


    def assign_cells_by_likelihood(self, data, cells=None, lamb=1000, cellassign=True):
        cell_scores, nodes = self.compute_node_likelihoods(data, cells,lamb)
        node_assign = np.argmin(cell_scores, axis=0)
        obj =cell_scores.min(axis=0)
        obj = obj.sum()
      
        if cellassign:
            phi = {}
            for i, k in enumerate(node_assign):
                phi[i] = nodes[k]
            ca = CellAssign(phi, self.clones())
            self.cost = obj
            return obj, ca 
        else:
            return obj 
    
    
    def get_path_map(self):
            all_shortest_paths = dict(nx.all_pairs_shortest_path(self.tree))

      
            clones = self.clones()
            has_path = {(u, v): v in all_shortest_paths[u] for u in clones for v in clones}
            return has_path
    
    @staticmethod
    def get_vafs(q, snv_tree, clones, cna_genos, has_path):
        vafs = np.zeros(len(clones))
        root= [n for n in snv_tree if snv_tree.in_degree[n] ==0][0]


        for idx, u in enumerate(clones):
      
        
            if has_path[u,q] and u != q:
                continue
            if has_path[q,u]:
                # if u not in cna_genos:
                #     print(u)
                #     print(cna_genos.keys())
                    
                cn_state = cna_genos[u]
                for v in nx.dfs_preorder_nodes(snv_tree, root):

                    if (v[0], v[1]) == cn_state and v[2] + v[3] > 0:
                        vafs[idx] = (v[2] + v[3])/(v[0] + v[1])
                        break
        return vafs 
    
    def optimize(self, dat, lamb, rho={}, max_iterations=10, dcfs={}):

        has_path = self.get_path_map()
        
        self.psi = self.get_psi()
        if not rho and self.rho is not None:
            rho = self.rho 
        if not rho:
            raise ValueError("A mapping rho of SNV clusters to valid SNV trees for each segment must be provided.")
        
        opt_cost = np.inf
        best_phi = None
        cna_genos = self.get_cna_genos()
        for i in range(max_iterations):

            obj, ca = self.assign_cells_by_likelihood(dat, lamb=lamb)

            
            if obj < opt_cost:
                best_phi = ca.phi.copy()
                opt_cost  = obj 
                best_geno = {v: self.genotypes[v].copy() for v in self.genotypes}
            else:
                    break

            self.assign_genotypes(dat, ca, rho=rho, cna_genos=cna_genos, has_path= has_path)
    
       
        
        if best_phi is not None:
            best_ca = CellAssign(best_phi, self.clones())
            self.genotypes = best_geno
            self.cost = opt_cost 
            self.snv_cost = np.nan
            self.cna_cost = np.nan

            # _ = self.compute_likelihood(dat, best_ca, lamb)
      
        
            self.update_mappings()

        return opt_cost, best_ca 
    
    @staticmethod
    def snv_tree_has_loss(snv_tree):
        for u,v in snv_tree.edges:
            if u[2]+u[3] > 0 and v[2] + v[3]==0:
                return True 
        return False
    
    def compute_genotype_costs(self, dat, phi):
        self.update_mappings()
        has_path = self.get_path_map()
        cell_counts  = phi.get_cell_count()
        clones  = [u for u in self.clones() if cell_counts[u] >0]
   
   
        cna_genos_all = self.get_cna_genos()

        def snv_tree_to_string(snv_tree):
            return ''.join([f"{u}->{v}|" for u,v in snv_tree.edges])
        
        def snv_tree_has_loss(snv_tree):
            for u,v in snv_tree.edges:
                if u[2]+u[3] > 0 and v[2] + v[3]==0:
                    return True 
            return False

        #store the results in the dataframe for analysis and comparison
        #snv | cluster | tree_id | tree_edge_list | is_lost | cost 

        results = []
        for ell, snvs in self.seg_to_muts.items():
            snv_cluster_mapping = self.rho[ell]
            snv_clusters = list(snv_cluster_mapping.keys()
            )
            cna_genos = cna_genos_all[ell]
 

            #order the SNV clusters in preorder
            for q in snv_clusters:
                snv_tree_list = []
                vaf_tree_list =[]
                for id, snv_tree in enumerate(snv_cluster_mapping[q]):
                    vafs = set([(n[2] + n[3])/ (n[0] + n[1]) for n in snv_tree])
                    if vafs not in vaf_tree_list:
                        snv_tree_list.append(id)
                        vaf_tree_list.append(vafs)
                snv_trees = [snv_cluster_mapping[q][i] for i in snv_tree_list]  
                 
        
                for id, snv_tree in enumerate(snv_trees):
                    tree_str  = snv_tree_to_string(snv_tree)
                    has_loss = snv_tree_has_loss(snv_tree)
                    all_vafs =self.get_vafs(q,  snv_tree, clones, cna_genos, has_path)
          
                    cluster_costs = np.vstack([dat.compute_snv_likelihoods(all_vafs[i], snvs, phi.get_cells(u)) for i,u in enumerate(clones)])
                    snv_likes = cluster_costs.sum(axis=0)
                    for j, cost in zip(snvs, snv_likes):
                            assigned_tree = snv_tree_to_string(self.get_snv_tree(j))
                            results.append([j,q, id, tree_str, has_loss, cost, self.psi[j], assigned_tree])
        return pd.DataFrame(results, columns=["snv", "cluster", "tree_id", "tree_edges", "is_lost", "cost", "assigment", "assigned_tree"])

    def get_snv_tree(self, j):
        edge_list = []
        for u,v in self.tree.edges:
            u_geno = self.genotypes[u][j]   
            v_geno = self.genotypes[v][j]
            if u_geno  != v_geno:
                edge_list.append((u_geno, v_geno))
        return nx.DiGraph(edge_list)       

    def compute_snv_likelihoods(self, dat, phi):
 
        """
        For a given cell assignment phi and observed data dat=(C,A,D), find the optimal
        assignment of SNV to SNV cluster, i.e., node where SNV is first introduced,
        and corresponding genotypes for all nodes and SNVs that maximizes likelihood of 
        read counts of the read counts A,D.

        Data dat: a Pharming data object containing (C,A,D) for each SNV and cell
        CellAssign phi: a Pharming cell assignment object containing the mapping phi of cells to nodes
        dict: rho: a dictionary of dictionaries for each segment, each containing a 
                    mapping of SNV clusters to a list of valid SNV trees (networx DiGraph)
        
        returns void (genotypes are modified in place)
        """

        self.update_mappings()
        rho = self.rho 
        


        has_path = self.get_path_map()


        cell_counts  = phi.get_cell_count()
        clones  = [u for u in self.clones() if cell_counts[u] >0]
        cna_genos_all = self.get_cna_genos()

        results = []
      
        for ell, snvs in self.seg_to_muts.items():
            for j in snvs:
                snv_tree  = self.get_snv_tree(j)
                tree_str = ''.join([f"{u}->{v}|" for u,v in snv_tree.edges])
                vafs = self.get_vafs(self.psi[j], snv_tree, clones, cna_genos_all[ell], has_path)
                snv_node_costs= np.vstack([dat.compute_snv_likelihoods(vafs[i], [j], phi.get_cells(u)) for i,u in enumerate(clones)])
                snv_like = snv_node_costs.sum()
                results.append([ell, j, snv_like, tree_str, self.psi[j]])

        return pd.DataFrame(results, columns=["segment", "snv", "cost", "tree_edges", "cluster"])





    # @utils.timeit_decorator
    def assign_genotypes(self, dat, phi, rho={}, seg_to_snvs={}, states=None,
                          cna_genos={}, has_path={}, start_state=None):
        """
        For a given cell assignment phi and observed data dat=(C,A,D), find the optimal
        assignment of SNV to SNV cluster, i.e., node where SNV is first introduced,
        and corresponding genotypes for all nodes and SNVs that maximizes likelihood of 
        read counts of the read counts A,D.

        Data dat: a Pharming data object containing (C,A,D) for each SNV and cell
        CellAssign phi: a Pharming cell assignment object containing the mapping phi of cells to nodes
        dict: rho: a dictionary of dictionaries for each segment, each containing a 
                    mapping of SNV clusters to a list of valid SNV trees (networx DiGraph)
        
        returns void (genotypes are modified in place)
        """

        if not rho and self.rho is not None:
            rho = self.rho 
        
        if not rho:
            raise ValueError("A mapping rho of SNV clusters to valid SNV trees for each segment must be provided.")
        
     
        if not seg_to_snvs:
            seg_to_snvs = self.seg_to_muts

        if not has_path:
            has_path = self.get_path_map()


        cell_counts  = phi.get_cell_count()
        clones  = [u for u in self.clones() if cell_counts[u] >0]
   

        if cna_genos:
            cna_genos_all = cna_genos
        else:
            cna_genos_all = self.get_cna_genos()

        bfs_nodes= list(nx.bfs_tree(self.tree, self.root))
      
        for ell, snvs in seg_to_snvs.items():
            snv_cluster_mapping = rho[ell]
            snv_clusters = list(snv_cluster_mapping.keys())

            #reorder the nodes so that we priorize snv clusters closer to the root
            snv_clusters =[u for u in bfs_nodes if u in snv_clusters]

            
            if ell in cna_genos_all:
                cna_genos = cna_genos_all[ell]
            else:
                cna_genos = {v: states[ell] for v in self.tree}
                # cna_genos[self.root] = (1,1)

                if start_state is not None:
                    cna_genos[self.root] = start_state
            all_cluster_costs = []
            tree_assign = {}

            #order the SNV clusters in preorder
            for q in snv_clusters:
            
                all_tree_costs = []
                
                #trees should be ordered from those with no loss to those with loss in rho
                for snv_tree in snv_cluster_mapping[q]:
                    all_vafs =self.get_vafs(q,  snv_tree, clones, cna_genos, has_path)
             
                    cluster_costs = np.vstack([dat.compute_snv_likelihoods(all_vafs[i], snvs, phi.get_cells(u)) for i,u in enumerate(clones)])
                    all_tree_costs.append(cluster_costs.sum(axis=0))
               
                if len(all_tree_costs) == 0:
                    raise Exception(f"No valid SNV trees for cluster {q} in segment {ell}.")
                tree_costs= np.vstack(all_tree_costs)
                tree_assign[q] = tree_costs.argmin(axis=0)
                all_cluster_costs.append(tree_costs.min(axis=0))
                
        
            clust_costs = np.vstack(all_cluster_costs)

            snv_cluster_assign = clust_costs.argmin(axis=0)
            clust_assign = defaultdict(list) 
        
            for j,q, index in zip(snvs, snv_cluster_assign, range(len(snvs))): 
                opt_clust = snv_clusters[q]
                if j not in self.psi:
                    self.mut_to_seg[j] = ell 
                    if ell not in self.seg_to_muts:
                        self.seg_to_muts[ell] = [j]
                    else:
                        self.seg_to_muts[ell].append(j)
                if j not in self.psi or self.psi[j] != opt_clust: 
                    self.psi[j] = opt_clust
                    clust_assign[opt_clust, tree_assign[opt_clust][index]].append(j)
            
            for q, index in clust_assign:
                snvs = clust_assign[q, index]
                snv_tree = snv_cluster_mapping[q][index]
                self.batch_update(q, snvs, snv_tree, cna_genos)




        
               
             
    
                #     self.update_genotype(opt_clust, j,snv_cluster_tree[j,opt_clust], cna_geno=cna_genos)
                
                # elif  psi[j] != opt_clust:     
                #     self.update_genotype(opt_clust, j,snv_cluster_tree[j,opt_clust], cna_genos)

    # def assign_genotypes(self, dat, phi, rho=None, seg_to_snvs=None, states=None, cna_genos=None):
    #     """
    #     For a given cell assignment phi and observed data dat=(C,A,D), find the optimal
    #     assignment of SNV to SNV cluster, i.e., node where SNV is first introduced,
    #     and corresponding genotypes for all nodes and SNVs that maximizes likelihood of 
    #     read counts of the read counts A,D.

    #     Data dat: a Pharming data object containing (C,A,D) for each SNV and cell
    #     CellAssign phi: a Pharming cell assignment object containing the mapping phi of cells to nodes
    #     dict: rho: a dictionary of dictionaries for each segment, each containing a 
    #                 mapping of SNV clusters to a list of valid SNV trees (networx DiGraph)
        
    #     returns void (genotypes are modified in place)
    #     """

    #     if rho is None and self.rho is not None:
    #         rho = self.rho 
        
    #     if rho is None:
    #         raise ValueError("A mapping rho of SNV clusters to valid SNV trees for each segment must be provided.")
    #     snv_cluster_tree = {}
     
    #     if seg_to_snvs is None:
    #         seg_to_snvs = self.seg_to_muts
    #     has_path = self.get_path_map()
    #     cell_counts  = phi.get_cell_count()
    #     clones  = [u for u in self.clones() if cell_counts[u] >0]
    #     psi = self.get_psi_test()
    #     if cna_genos is not None:
    #         cna_genos_all = cna_genos
    #     else:
    #         cna_genos_all = self.get_cna_genos()

      
    #     for ell, snvs in seg_to_snvs.items():
    #         snv_cluster_mapping = rho[ell]
    #         snv_clusters = list(snv_cluster_mapping.keys())
            
       
    #         if ell in cna_genos_all:
    #             cna_genos = cna_genos_all[ell]
    #         else:
    #             cna_genos = {v: states[ell] for v in self.tree}

        
    #         all_cluster_costs = []
    #         for q in snv_clusters:
           
    #             all_tree_costs = []
    #             for snv_tree in snv_cluster_mapping[q]:
    #                 all_vafs =self.get_vafs(q,  snv_tree, clones, cna_genos, has_path)

    

    #                 cluster_costs = np.zeros(shape= len(snvs))
    #                 assert all_vafs.shape[0] == len(clones)

    #                 for i,u in enumerate(clones):
    #                     cells = phi.get_cells(u)
    #                     clone_cost = dat.compute_snv_likelihoods(all_vafs[i], snvs, cells)
    #                     cluster_costs += clone_cost
                
    #                 all_tree_costs.append(cluster_costs)
                    
    #             tree_costs= np.vstack(all_tree_costs)
    #             tree_assign = tree_costs.argmin(axis=0)
    #             clust_costs = tree_costs.min(axis=0)
    #             all_cluster_costs.append(clust_costs)

    #             for j, index  in zip(snvs, tree_assign):
    #                 snv_cluster_tree[j,q]  =rho[ell][q][index]
                
        
    #         clust_costs = np.vstack(all_cluster_costs)

    #         snv_cluster_assign = clust_costs.argmin(axis=0)
      
    #         for j,q in zip(snvs, snv_cluster_assign): 
        
    #             opt_clust = snv_clusters[q]
    #             if j not in psi:
    #                 self.mut_to_seg[j] = ell 
    #                 if ell not in self.seg_to_muts:
    #                     self.seg_to_muts[ell] = [j]
    #                 else:
    #                     self.seg_to_muts[ell].append(j)
    
    #                 self.update_genotype(opt_clust, j,snv_cluster_tree[j,opt_clust], cna_geno=cna_genos)
                
    #             elif  psi[j] != opt_clust:     
    #                 self.update_genotype(opt_clust, j,snv_cluster_tree[j,opt_clust])
            


            


    def compute_likelihood(self, data,  cellAssign, lamb=0):



        self.node_cost = {}
        self.snv_node_cost = {}
        self.cna_node_cost = {}
        self.cost = 0
        cna_genos = self.get_cna_genos()
        
        for v in self.tree:
                cells = cellAssign.get_cells(v)
                if len(cells) ==0:
                    continue      

                else:
                    lat_genos = self.get_latent_genotypes(v)
                    self.snv_node_cost[v] = self.node_likelihood(data, cells, self.get_latent_vafs(lat_genos))
                    self.cna_node_cost[v] = self.node_cna_cost(v, cells, data,cna_genos, weights=data.seg_weights).sum()
                    self.node_cost[v]=  self.snv_node_cost[v] + lamb * self.cna_node_cost[v]
             
             


        self.cost = sum([score for _, score in self.node_cost.items()])

        self.snv_cost = sum(self.snv_node_cost[u] for u in self.snv_node_cost)
        self.cna_cost = sum(self.cna_node_cost[u] for u in self.snv_node_cost)

        return self.cost 
    def ICL(self, phi, data, lamb):
        self.update_mappings()
        self.bic = self.BIC(phi, data, lamb)
        # self.entropy = self.ENT(phi, data, lamb)
        return self.bic #+ self.entropy

    def BIC(self, phi, data, lamb):
        snvs = self.get_all_muts()
        m = len(snvs)
  
        _ = self.compute_likelihood(data, phi, lamb)
        like= -1* self.snv_cost 


        k = len(set([self.psi[j] for j in snvs]))
  
        self.bic = np.log(m)*(k) - 2*like
        # print(self.bic)
        return self.bic 


    def ENT(self, phi, dat, lamb):
        
       
        rho = self.rho 
        
        seg_to_snvs = self.seg_to_muts
        has_path = self.get_path_map()
        cell_counts  = phi.get_cell_count()
        clones  = [u for u in self.clones() if cell_counts[u] >0]

        cna_genos_all = self.get_cna_genos()

        like_by_segment = []
        for ell, snvs in seg_to_snvs.items():
            snv_cluster_mapping = rho[ell]
            snv_clusters = list(snv_cluster_mapping.keys())
            cna_genos = cna_genos_all[ell]
            all_cluster_costs = []
            tree_assign = {}

            for q in snv_clusters:
            
                all_tree_costs = []
                for snv_tree in snv_cluster_mapping[q]:
                    all_vafs =self.get_vafs(q,  snv_tree, clones, cna_genos, has_path)
             
                    cluster_costs = np.vstack([dat.compute_snv_likelihoods(all_vafs[i], snvs, phi.get_cells(u)) for i,u in enumerate(clones)])
                    all_tree_costs.append(cluster_costs.sum(axis=0))
               
                if len(all_tree_costs) == 0:
                    raise Exception(f"No valid SNV trees for cluster {q} in segment {ell}.")
                tree_costs= np.vstack(all_tree_costs)
                tree_assign[q] = tree_costs.argmin(axis=0)
                all_cluster_costs.append(tree_costs.min(axis=0))
                
        
            clust_costs = np.vstack(all_cluster_costs)
            like_by_segment.append(clust_costs)
        #should be a k x m numpy array
        likes = np.hstack(like_by_segment)
        log_norm_likes = logsumexp(likes, axis=0)  # Shape (1, m)

        # Compute the posterior probabilities
        log_gamma = likes - log_norm_likes  # Broadcasting subtraction
        gamma = np.exp(log_gamma)  # Shape (k, m)

        # Compute the entropy term
        ent = -np.sum(gamma * log_gamma, axis=0)  # Summing over clusters (axis=0)

        # Total entropy term for all segments
        total_ent = np.sum(ent) 
        self.ent = total_ent
        return self.ent 





    # def assign_cells(self, data, lamb=0, lamb_vaf= 1, cna_only=False):

    #     nodes = np.array(self.tree)
    #     if cna_only:
    #         cell_scores =0
    #         lamb=1
    #     else:
    #         cell_scores = np.vstack([self.node_snv_cost(v, data.cells, data, self.get_latent_genotypes(v)) for v in nodes])

    #     if lamb > 0:
    #         cell_cna_scores = np.vstack([self.node_cna_cost(v, data.cells, data, self.get_latent_genotypes(v)) for v in nodes])
    #         cell_scores = lamb_vaf*cell_scores + lamb*cell_cna_scores
    #     assignments = np.argmin(cell_scores, axis=0)
    #     self.cost = np.nansum(np.nanmin(cell_scores,axis=0)).sum()
        
        
    #     phi = {i: v for i,v in zip(data.cells, nodes[assignments])}
  
    #     # self.phi_to_cell_mapping(self.phi)
    #     return CellAssign(phi, self.clones()), self.cost, cell_scores, nodes   
    #     # self.cell_mapping = defaultdict(list)
    #     # for i, v in self.phi.items():
    #     #     self.cell_mapping[v].append(i)
    #     # self.cell_mapping = dict(self.cell_mapping)

    def filter_snvs(self, snvs_to_keep):

        snvs_to_remove = set(self.get_all_muts()) - set(snvs_to_keep)
  
        for j in snvs_to_remove:
                self.del_snv(j)
        
        self.update_mappings()
        
        
    def del_snv(self, j):
        seg = self.mut_to_seg[j]
        del self.mut_to_seg[j]
        self.seg_to_muts[seg].remove(j)
        # node = self.psi[j]
        # if j in self.mut_mapping[node]:
        #     self.mut_mapping[node].remove(j)
        
        for v in self.tree:
            del self.genotypes[v][j]
    
    def batch_update(self, node, snvs, snv_tree, cna_geno):
        """
        assume all SNVS are in the same segment
        """
        # snv_tree_root = [n for n in snv_tree if snv_tree.in_degree[n]==0][0]
        pres_nodes = set()
    
        for u in self.preorder(node):
            cn_state = cna_geno[u]
            for v in snv_tree:
             
                if (v[0], v[1]) == cn_state and self.mut_copies(v) > 0:
                    for j in snvs:
                        self.genotypes[u][j] = v
                    pres_nodes.add(u)
                    break 
        
        for u in self.clones().difference(pres_nodes):
            for j in snvs:
                self.genotypes[u][j] = (*cna_geno[u], 0, 0)
        
        


    # def update_genotype(self, node, j, snv_tree, cna_geno=None):
    #     snv_tree_root = [n for n in snv_tree if snv_tree.in_degree[n]==0][0]
    #     if cna_geno is None:
    #         cna_geno = self.get_cna_genos()[self.mut_to_seg[j]]
    #     # if j == 68:
    #     #     self.draw("test/self..png", segments=[2])
    #     #     snv_tree.draw("test/snv_tree68.png", segments=[2])

    #     pres_nodes = []
    #     for u in sorted(nx.descendants(self.tree, node) | {node}):
    #         cn_state = cna_geno[u]
    #         # cn_state = cn_geno.to_tuple()
      
    #         added = False
    #         for v in nx.dfs_preorder_nodes(snv_tree, snv_tree_root):
             
    #             if (v[0], v[1]) == cn_state and self.mut_copies(v) > 0:
    #                 self.genotypes[u][j] = v
    #                 added = True
    #                 pres_nodes.append(u)
    #                 break 
    #         # if added:
    #         #     break  # Early termination if condition is met
   
    #     for u in self.clones().difference(pres_nodes):
    #             cn_state = cna_geno[u]
    #             self.genotypes[u][j] = (*cn_state, 0, 0)

    #     self.psi[j] = node
    
    # def update_genotype(self, node, j, snv_tree):
    #     cna_geno = self.get_cna_genos()[self.mut_to_seg[j]]

    #     pres_nodes = []
    #     for u in sorted(nx.descendants(self.tree, node) | {node}):
    #         cn_geno = cna_geno[u]
    #         cn_state= cn_geno.to_tuple()
    #         pres_nodes.append(u)
    #         added = False
    #         for v in snv_tree.preorder():
    #             geno = snv_tree.genotypes[v][j]
    #             if (geno.x, geno.y) == cn_state and geno.z > 0:
    #                 self.genotypes[u][j] = genotype(*geno.to_tuple())
    #                 added = True
    #                 break 
    #         if not added:
    #             # if j ==384:
    #             #     print(f"SNV {j} not added to node {node}")
    #             #     print(snv_tree)
    #             #     snv_tree.draw("test/snv_tree.png", segments=[10])
    #             self.genotypes[u][j] = genotype(*cn_state, 0,0)
        
    #     for u in self.clones().difference(pres_nodes):
    #         cn_state = cna_geno[u].to_tuple()
    #         self.genotypes[u][j] = genotype(*cn_state, 0,0)
        
    #     self.psi[j] = node 
        # self.mut_mapping[old_node].remove(j)
    
     

     

    # def compute_likelihood(self, data,  cellAssign, lamb=0, cna_only=False):
    #     if cna_only:
    #         lamb = 1
    #     # if  cellAssign is None:
    #     #     self.assign_cells(data, lamb)

    #     self.node_cost = {}
    #     self.cost = 0
    #     for v in self.tree:
    #             cells = cellAssign.get_cells(v)
    #             if len(cells) ==0:
    #                 continue
    #             latent_genos = self.get_latent_genotypes(v)
            
    #             if cna_only:
    #                 self.node_cost[v] =0
    #             else:
                    
    #                 self.node_cost[v]= self.snv_cost_func(v, cells, data, latent_genos).sum() 
    #             if lamb >0:
                    
    #                 self.node_cost[v]+=lamb*self.node_cna_cost(v, cells,data, latent_genos).sum()
             
        
    #     self.cost = sum([score for _, score in self.node_cost.items()])

    #     return self.cost 

    # def compute_costs(self, data,  cellAssign, lamb=0, cna_only=False):
    #     if cna_only:
    #         lamb = 1
    #     # if  cellAssign is None:
    #     #     self.assign_cells(data, lamb)

    #     self.node_cost = {}
    #     self.cost = 0
    #     for v in self.tree:
    #             cells = cellAssign.get_cells(v)
    #             if len(cells) ==0:
    #                 continue
    #             latent_genos = self.get_latent_genotypes(v)
            
    #             if cna_only:
    #                 self.node_cost[v] =0
    #             else:
                    
    #                 self.node_cost[v]= self.snv_cost_func(v, cells, data, latent_genos).sum() 
    #             if lamb >0:
                    
    #                 self.node_cost[v]+=lamb*self.node_cna_cost(v, cells,data, latent_genos).sum()
             
        
    #     self.cost = sum([score for _, score in self.node_cost.items()])

    #     return self.cost 
    
    # def compute_pooled_costs(self, data, cellAssign, lamb=0, cna_only=False):
    #     self.all_dfs = {}
    #     self.snv_cost_func = self.pooled_node_snv_cost
    #     cost = self.compute_costs(data,cellAssign, lamb, cna_only)
    #     self.snv_cost_func = self.node_snv_cost
    #     # self.comp_df = pd.concat([vals for v, vals in self.all_dfs.items()])
    #     # print(self.comp_df.head())

    #     return cost 

    




   #-------------------------- Save Methods ---------------------------------------#
    def draw(self, fname, cellAssign=None, mapping=None,segments=None, include_dcfs=False, cmap='Set3'):
        self.update_mappings()

        if segments is not None:
            cna_genos = self.get_cna_genos()
        mut_count = {n : len(self.mut_mapping[n]) for n in self.mut_mapping}

        if include_dcfs and cellAssign is not None:
            dcfs = self.compute_dcfs(cellAssign)
            include_dcfs = True 
        else:
            include_dcfs = False 

        if cellAssign is not None:
            cell_count = cellAssign.get_cell_count()
            if hasattr(cellAssign, 'n'):
                ncells = cellAssign.n 
            else:
                ncells = sum(cell_count[u] for u in cell_count)
        else:
            cell_count = {n: 0 for n in self.clones()}
            ncells = 0
        labels = {}
        # color_values, colormap = self.create_color_map(cmap)
        for n in self.tree:
                if mapping is not None:
                    labels[n] = str(mapping[n])
                else:
                    labels[n] = str(n)
                if include_dcfs and n in dcfs:
                    labels[n] += f"\ndelta:{round(dcfs[n],3)}"
                if n in cell_count:
                    if cell_count[n] > 0:
                        labels[n] += "\nCells:" + str(cell_count[n])
                # SNV
                if n in self.mut_mapping:
                    if mut_count[n] > 0:
                        labels[n] += "\n+SNVs:" + str(mut_count[n])
                    if len(self.mut_loss_mapping[n]) > 0:
                        labels[n] += "\n-SNVs:" + str(len(self.mut_loss_mapping[n]))
                if segments is not None:
                    labels[n] += "\n" + ",".join([f"{cna_genos[s][n][0]}|{cna_genos[s][n][1]}" for s in segments if s in self.seg_to_muts])

        like_label = f"Segment\n"
        tree = pgv.AGraph(strict=False, directed=False)
        tree.node_attr['style']='filled'
        segs = self.get_segments()
        segs = [str(ell) for ell in segs ]
        if self.cost is not None:
            total_like = np.round(self.cost)
            if hasattr(self, "snv_cost"):
                snv_cost = np.round(self.snv_cost)
                cna_cost = np.round(self.cna_cost)
            else:
                snv_cost = 0
                cna_cost = 0
        
            # tree.graph_attr["label"] = f"Objective: {total_like}\nSegments: {','.join(segs)}\nn={cellAssign.n} cells\nm={len(self.get_all_muts())} SNVs"
            tree.graph_attr["label"] = f"Objective: {total_like}\nSNV:{snv_cost}, CNA:{cna_cost}\nSegments: {','.join(segs)}\nn={ncells} cells\nm={len(self.get_all_muts())} SNVs"

 
        # colormap = cm.get_cmap(cmap)
        for n in self.tree:

            tree.add_node(n, label=labels[n])
      
            node_attr = tree.get_node(n)
            color_value = None
            # try:
            #     x,y = self.tree.nodes[n]["genotype"]
            #     color_value = x+y
            # except:
            #     color_value = None
        
            # if color_value is not None:
               
                
            #     # color = colormap(color_value)
            #     # hex_color = mcolors.rgb2hex(color)
            #     node_attr.attr['fillcolor'] =SET3_HEX[color_value]
                # node_attr['fillcolor'] = hex_color
    
        tree.add_edges_from(list(self.tree.edges))
        tree.layout("dot")
        tree.draw(fname)
  
    def write_psi(self, fname):
        psi = self.get_psi()
        df = pd.DataFrame(list(psi.items()), columns=['mutation', 'cluster'])
        df.to_csv(fname, index=False)


    def write_loss_mapping(self, fname):
        self.update_mappings()
        loss = {}
        for u, snvs in self.mut_loss_mapping.items():
            for j in snvs:
                loss[j] = u
        
        df = pd.DataFrame(list(loss.items()), columns=['mutation', 'cluster'])
        df.to_csv(fname, index=False)
    
    def write_likelihood(self, fname):
        with open(fname, "w+") as file:
            file.write("likelihood,snv,cna\n")
            file.write(f"{self.cost},{self.snv_cost},{self.cna_cost}\n")

    def write_genotypes(self, fname):
        colnames= ["node", "snv", "x", "y", "xbar", "ybar", "segment"]

        with open(fname, "w+") as file:
            file.write(",".join(colnames) + "\n")
            for n in self.tree:
                for ell, snvs in self.seg_to_muts.items():
                    for j in snvs:
                        g =self.genotypes[n][j]
                        file.write(f"{n},{j},{g[0]},{g[1]},{g[2]},{g[3]},{ell}\n")

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
    
    def cna_loss(self,  snvs:list):
        """Calculate loss precision, recall and accurracy by segments."""
        snvs = set(snvs)
        tp, tn, fn, fp = 0,0,0,0
        inf_lost_snvs = set(self.get_lost_snvs())
        for ell, muts in self.seg_to_muts.items():
            
            gt_seg_lost = len(set(muts).intersection(snvs)) > 0
            inf_seg_lost = len(set(muts).intersection(inf_lost_snvs)) > 0
            if gt_seg_lost and inf_seg_lost:
                tp += 1
            elif gt_seg_lost and not inf_seg_lost:
                fn += 1
            elif not gt_seg_lost and inf_seg_lost:
                fp += 1
            else:
                tn += 1   
        return tp, tn, fn, fp
    
    def generate_mut_dataframe( self,lookup):
    #convert mapping to series in order of mutations
        pred_df= self.mapping_to_dataframe(self.mut_mapping, "mutation_id")

        pred_df["mutation"] = lookup[pred_df['mutation_id']].values

        pred_df = pred_df.drop(['mutation_id'], axis=1)
        return pred_df
    
    # def generate_cell_dataframe( self,lookup):
    # #convert mapping to series in order of cells
    #     pred_df= self.mapping_to_dataframe(self.cell_mapping, "cell_id")

    #     pred_df["cell"] = lookup[pred_df['cell_id']].values

    #     pred_df = pred_df.drop(['cell_id'], axis=1)
    #     return pred_df 
    

    def generate_results(self, cell_lookup, mut_lookup):
        pcell = self.generate_cell_dataframe(cell_lookup)
        pmut = self.generate_mut_dataframe(mut_lookup)
        # ploss = generate_mut_dataframe(self.mut_loss_mapping, mut_lookup)


        return pcell, pmut

    #--------------------- End Save Functions ------------------------------------#
  

    def get_node_muts(self, t):
        if self.mut_mapping is None:
            self.update_mappings()
        if t in self.mut_mapping:
            return self.mut_mapping[t]
        return []
    
    
    # def get_node_cells(self, t):
    #     if t in self.cell_mapping:
    #         return self.cell_mapping[t]
    #     return []

    def get_leaves(self):
        leaves = [l for l in list(self.tree.nodes())
                  if self.tree.out_degree(l) == 0]
        return leaves

   

    # def get_cell_assignments(self):
    #     n_cell = len(self.get_all_cells())
    #     clusters = np.zeros(n_cell, dtype=int)
    #     for cluster, cells in self.cell_mapping.items():
    #         if len(cells) > 0:
    #             clusters[cells] = cluster
    #     return clusters

    # def get_all_muts(self):
    #     muts = list(chain.from_iterable([mut for n,mut in self.mut_mapping.items()]))
    #     muts.sort()
    #     return muts
    def get_all_muts(self):
       return  list(self.mut_to_seg.keys())
    

    def get_mut_assignments(self, n_mut=None):
        '''

        '''
        if self.mut_mapping is None:
            self.update_mappings()
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
   