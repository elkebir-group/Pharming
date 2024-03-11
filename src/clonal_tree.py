
import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter, defaultdict 
from itertools import product, chain, combinations
import pickle 
import pygraphviz as pgv
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from sklearn.metrics.cluster import adjusted_rand_score
from genotype import genotype 
from cell_mapping import CellAssign
from tree_to_json import convertToJson
from scipy.stats import binom
from scipy.stats import beta

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
    
    genotypes: dict (node ids as keys) of dict (snv index as keys) of genotype datastructs
        genotypes[v][s][m]
    
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

    def __init__(self, tree, genotypes: dict, seg_to_muts:dict=None,  key=0):
        self.tree: nx.DiGraph = tree
        self.genotypes = genotypes 
        if seg_to_muts is None:
            self.seg_to_muts = {}
        else:
            self.seg_to_muts = seg_to_muts
        
        self.mut_to_seg = {}
        for seg, muts in self.seg_to_muts.items():
            for m in muts:
                self.mut_to_seg[m] = seg
        
        self.root= self.find_root()

        #double check that all snvs have the same CNA copy number profile for each node
        # for v in self.tree:
        #     cna_genos =[self.genotypes[v][j].to_CNAgenotype().to_tuple() for j in self.genotypes[v]]
        #     assert all(t == cna_genos[0] for t in cna_genos)

        

        # if cell_mapping is None:
        #     self.cell_mapping = {}
            
        # else:    
        #     self.cell_mapping = cell_mapping


        self.mut_mapping, self.mut_loss_mapping = self.get_mut_mapping()
        self.psi = self.get_psi()

        self.key = key
        self.cost = np.Inf

        # self.snv_cost_func = self.compute_node_likelihoods

    
    def __str__(self) -> str:
        # all_cells = self.get_all_cells()
        all_snvs  = self.get_all_muts()
        mystr =f"\nClonal Tree {self.key}"
        mystr += "\n--------------------------------"

        mystr +=f"\n{len(list(self.tree.nodes))} Nodes and {len(list(self.tree.edges))} Edges"
        mystr += f"\n{len(all_snvs)} SNVs"
        # mystr+= f"\n{len(all_cells)} cells\n"


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
               
        mystr += "\n\nNode\tSNVs\tGenotype"
        mystr += "\n--------------------------------"
        for n in  node_order:
            # genotype = self.tree.nodes[n]['genotype']
            genotype = ""
            # if n in self.cell_mapping:
            #     ncells = len(self.cell_mapping[n])
            # else:
            #     ncells = 0
            if n in self.mut_mapping:
                nmuts = len(self.mut_mapping[n])
            else:
                nmuts = 0
            mystr += f"\n{n}\t{nmuts}\t{genotype}"
        mystr += "\n"
        mystr += f"\nJ(T,phi): {self.cost}\n"
        return mystr
    
    def preorder(self, source=None):
        if source is None:
            source = self.root 
    
        return list(nx.dfs_preorder_nodes(self.tree, source=source))

    def get_tree(self):
        return self.tree.copy()

    def get_genotypes(self):
        return self.genotypes.copy()

    def get_seg_to_muts(self):
        return self.seg_to_muts.copy()
    
    def get_segments(self):
        return set(self.seg_to_muts.keys())
    
    def compute_dcfs(self, ca):
        dcfs = {}
        counts = ca.get_cell_count()
        ncells = len(ca.get_all_cells())
        for n in self.tree:
            
            if len(self.mut_mapping[n]) > 0:
                snvs = self.mut_mapping[n]
                m = snvs[0]
                dcfs[n] = sum([counts[u] for u in self.preorder(source=n)])/ncells 
        return dcfs 
    


    def posterior_dcf(self, j:int, dcf:float,a: int,d:int ,cn_prop:dict ):

        # if j == 290:
        #     print(f"{j}: {dcf}")
        if d==0:
            posterior = EPSILON
        else:
            parent, u, p_geno, u_geno = self.get_split_nodes(j)
            m_star = u_geno.z 
            #compute fractional copy number as weighted sum 
            F = sum([ (cn[0]+ cn[1])*cn_prop[cn] for cn in cn_prop])
            desc_geno = self.get_desc_genotypes(u, j)
            # cn_prop = {x+y: 1.0 if x+y==cn else 0.0 for x,y,m in self.gamma}
            #v = self.dcf(dcf, cn)
            #copy code so that posterior 
            v = (dcf*m_star)/F 
            for geno in desc_geno:
                cn_state = geno.to_CNAgenotype().to_tuple()
                v += (1/F)* (geno.z-m_star)*cn_prop[cn_state]
                # + (1/F)*sum([(geno.w-m_star)*cn_prop[cn_state] for x,y,m in self.desc_genotypes])
    
            # posterior = max(beta.logpdf(v, a+1, d-a + 1),EPSILON)
            posterior = max(binom.logpmf(a, d,v), EPSILON)
            # if posterior > 0:
            #     print(f"{posterior}: {dcf} : {v} : {a/d}")
        return posterior



    # def compute_dcfs(self, ca, segment):
    #     counts = ca.get_cell_count()
    #     total_cells = sum(counts[u] for u in counts)
    #     cna_genos = self.get_cna_genos()[segment]
    #     # cna_states = set(cna_genos[segment][v].to_tuple() for v in self.tree)
    #     dcfs = {}
    #     cell_counts = {n: 0 for n in self.tree}
    #     res = []
    #     for n in self.tree:
    #         if len(self.mut_mapping[n]) > 0:
    #             # snvs = self.mut_mapping[n]
    #             # m = snvs[0]
  
    #             for u in self.preorder(source=n):
    #                         count = counts[u]
    #                         res.append([n, cna_genos[u].to_tuple(), count])
    #                         cell_counts[n] += count
        
    #     dcfs = {u: cell_counts[u]/total_cells for u in cell_counts}
    #     df = pd.DataFrame(res, columns=["cluster",  "state", "count"])
    #     dcf_by_state = df.groupby(["cluster", "state"]).sum("count")
    #     print(dcf_by_state)
    #     dcf_by_state = dcf_by_state/total_cells
    #     print(dcf_by_state)
    #     return dcfs, dcf_by_state
        


    def mutation_cluster_tree(self):
        cluster_nodes = [n for n in self.tree if len(self.mut_mapping[n]) > 0]

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
    
    def get_cna_genos(self):
        self.cna_genotypes = {}
        

        for seg, snvs in self.seg_to_muts.items():
            m = snvs[0]
            self.cna_genotypes[seg] = {}
            for v in self.tree:
                self.cna_genotypes[seg][v] = self.genotypes[v][m].to_CNAgenotype()
         
        return self.cna_genotypes
    


    def get_cna_tree(self, seg):
        self.get_cna_genos()
        S = nx.DiGraph()
        for u,v in self.tree.edges:
            u_cna = self.cna_genotypes[seg][u] 
            v_cna = self.cna_genotypes[seg][v]
            if u_cna != v_cna:
                S.add_edge(u_cna.to_tuple(), v_cna.to_tuple())
        if len(S) ==0:
            S.add_node(self.cna_genotypes[seg][self.root].to_tuple())
        return S
    
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

    def get_split_nodes(self, j):
        
        for u in self.preorder():
            if u != self.root:
                parent= self.parent(u)
                p_geno =self.genotypes[parent][j]
                u_geno = self.genotypes[u][j]
                if p_geno.cna_eq(u_geno) and u_geno.z > 0:
                    return parent, u, p_geno, u_geno

    def trim(self, cellAssign):
        '''
        remove nodes that only have cell assignments and pushees the cell assignments up the tree
        '''
        cell_mapping = cellAssign.cell_mapping
        to_del = []
        for k in cell_mapping:
            if len(self.mut_mapping[k])==0 and len(self.mut_loss_mapping[k])==0:
                if len(cell_mapping[k]) > 0 and self.is_leaf(k):
                    path =nx.shortest_path(self.tree, source=self.root, target=k)
                    path = path[::-1]
                    for p in path:
                        if len(self.mut_mapping[p]) >0 or len(self.mut_loss_mapping[p]) > 0:
                            cell_mapping[p] += cell_mapping[k] 
                            to_del.append(k)
                            self.tree.remove_node(k)
                            break
        for k in  to_del:
            del cell_mapping[k]
        cellAssign.update_phi()
        cellAssign.update_clones(self.clones())
    

    
            





    def toJson(self, fname, segment_csv=None, snv_csv=None):
        json = convertToJson(self, segment_csv, snv_csv)
        with open(fname, "w") as f:
            f.write(str(json))
    
    def filter_snvs(self,  snvs):
        snvs = set(snvs)
        # snvs = set(self.seg_to_muts[seg])
        keys = set(self.genotypes[self.root].keys())
        to_del  = keys - snvs 


        for v in self.genotypes:
            self.mut_mapping[v] = [j for j in self.mut_mapping[v] if j in snvs]

            for s in to_del:
                del self.genotypes[v][s]
        
        for s in to_del:
            seg = self.mut_to_seg[s]
            self.seg_to_muts[seg].remove(s)
            del self.mut_to_seg[s]
        
    
    
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
        self.mut_mapping[node].append(j)
        self.psi[j] = node 


        def match_snv_genotype(tree, cn_state):
            for u in tree:
                geno = genotype(*u)
                if geno.x_bar + geno.y_bar > 0 and geno.x == cn_state[0] and geno.y==cn_state[1]:
                    return geno
        
        def update_desc_genotype(u, snv_tree):
            descendants.append(u)
       
            cn_state = self.cna_genotypes[seg][u].to_tuple()
            self.genotypes[u][j] =match_snv_genotype(snv_tree, cn_state)

            for v in self.tree.successors(u):
                update_desc_genotype(v, snv_tree)
        
        descendants = []
        self.get_cna_genos()
        update_desc_genotype(node, snv_tree)

        not_present = set(self.clones())  - set(descendants)
        for u in not_present:
            self.genotypes[u][j] =  self.cna_genotypes[seg][u].to_genotype(0,0)
        




        
    

    def clones(self):
        return set(self.tree.nodes)
    
    def is_leaf(self,node):
        return self.tree.out_degree[node] ==0
    
    # def get_phi(self):
    #      self.phi = {i : k  for k, cells in self.cell_mapping.items() for i in cells}
    #      return self.phi
    
    def get_psi(self):
         self.psi = {m : k  for k, snvs in self.mut_mapping.items()  for m in snvs}
         return self.psi

    # def phi_to_cell_mapping(self, phi):
    #     self.phi = phi 
    #     self.cell_mapping = {v: [] for v in self.tree}
    #     for i, v in self.phi.items():
    #         self.cell_mapping[v].append(i)
      
    
    def get_muts(self,node):
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
        
        self.mut_mapping, self.mut_loss_mapping = self.get_mut_mapping()
        self.psi = self.get_psi()

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

    # def has_loss(self):
    #     return len(self.mut_loss_mapping) > 0

    def has_snvs(self, n):
        return len(self.mut_mapping[n]) > 0
    
    # def set_cell_mapping(self, cell_mapping):
    #     self.cell_mapping  = cell_mapping
    #     # self.psi = self.get_psi()
    #     self.phi = self.get_phi()
    #     # print("Warning, cost is not automatically updated, ensure cost is manually set.")

    # def clear_cell_mapping(self):
    #     self.phi = {}
    #     self.cell_mapping = {}
    #     self.cost = np.Inf

    def prune(self, cellAssign):
        counts = cellAssign.get_cell_count()
        self.update_mappings()
        for u in self.preorder():
            if u ==self.root:
                continue
            if len(self.children(u)) ==1:    
                if counts[u] ==0 and len(self.get_muts(u)) ==0 and len(self.mut_loss_mapping[u]) ==0:
                    self.prune_tree(u,cellAssign)
        

    def prune_tree(self, node_to_remove, cellAssign):
        
        parent = next(self.tree.predecessors(node_to_remove))

        # Get the children of the node to be removed
        children = list(self.tree.successors(node_to_remove))

        # Remove the node to be removed from the graph
        self.tree.remove_node(node_to_remove)

        # Reattach the children to the parent
        for child in children:
            self.tree.add_edge(parent, child)

        # if node_to_remove in self.cell_mapping:
        #     del cellAssign.cell_mapping[node_to_remove]
        
        if node_to_remove in self.genotypes:
            del self.genotypes[node_to_remove]
        
        if node_to_remove in self.mut_mapping:
            del self.mut_mapping[node_to_remove]

        if node_to_remove in self.mut_loss_mapping:
            del self.mut_loss_mapping[node_to_remove]
        
        cellAssign.update_clones(self.clones())

             
     
    
  

    def get_mut_mapping(self):
        if len(self.genotypes) ==0:
            return {v: [] for v in self.tree}, {v: [] for v in self.tree}
        gained= set()
        lost = set()
    
        muts= set(list(chain.from_iterable([self.seg_to_muts[ell] for ell in self.seg_to_muts])))
      
        mut_mapping = {v: [] for v in self.tree}
        mut_loss_mapping = {v: [] for v in self.tree}
        for u in self.preorder():
                for v in self.children(u):
                    for j in muts:
                        geno_v = self.genotypes[v][j]
                        geno_u = self.genotypes[u][j]
                        if geno_u.z ==0 and geno_v.z > 0:
                            mut_mapping[v].append(j)
                            gained.add(j)
                        if geno_u.z > 0 and geno_v.z == 0:
                            mut_loss_mapping[v].append(j)   
                            lost.add(j)


                    # for m, geno in self.genotypes[v].items():
                    # if m not in muts:
                    #     muts.append(m)
                    # # if m in [523, 792, 451,831]:
                    # #     print(m)
                    # if geno.z > 0 and m not in gained:
                    #     mut_mapping[v].append(m)
                    #     gained.append(m)
                    # if m in gained and geno.z ==0 and m not in lost:
                    #     mut_loss_mapping[v].append(m)
                    #     lost.append(m)
        missing= muts - (gained.union(lost))
        if len(missing) > 0:
            print(f"Warning: {len(missing)} SNVs never gained (w+z > 0) in any node, appending SNVs to root with 0 mutation state")
            # print(missing)
        for m in missing:
            mut_mapping[self.root].append(m)
        return mut_mapping, mut_loss_mapping 
    
    def update_mappings(self):

        #TODO: update psi
        self.mut_mapping, self.mut_loss_mapping =self.get_mut_mapping()

        for n, snvs in self.mut_mapping.items():
            for j in snvs:
                self.psi[j] = n

        
   
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
       
        path = list(nx.shortest_path(self.tree, self.root, node))
        present_muts =list(chain.from_iterable([self.mut_mapping[p] for p in path if p in self.mut_mapping]))
        # lost_muts = list(chain.from_iterable([self.mut_loss_mapping[p]
        #                 for p in path if p in self.mut_loss_mapping]))
        present_muts = set(present_muts) #- set(lost_muts))

        return present_muts
    
  
    
  
    
    def get_latent_genotypes(self,v):
        return {m: geno for m,geno in self.genotypes[v].items()} 
    
    def get_latent_vafs(self, latent_genos):
        return  {m: geno.vaf for m, geno in latent_genos.items()}               
     
    
    def get_latent_cna_genos(self, latent_genos):
        return  {m: geno.to_CNAgenotype() for m, geno in latent_genos.items()}    

  

    def node_cna_cost(self, v, cells, data, latent_genos):
        #dictionary is snvs as keys and cna_geno types as values
        latent_cna_genos= self.get_latent_cna_genos(latent_genos)
        # snv_keys = list(latent_cna_genos.keys())
        # snv_to_seg = {s: data.snv_to_seg[s] for s in snv_keys}

        seg_geno = {}
        for key, seg in self.mut_to_seg.items():
            if seg not in seg_geno:
               seg_geno[seg] = latent_cna_genos[key]
        latent_x, latent_y = [], []
        for seg,cna_geno in seg_geno.items():
            latent_x.append(cna_geno.x)
            latent_y.append(cna_geno.y)
        
        latent_x = np.array(latent_x).reshape(1,-1)
        latent_y = np.array(latent_y).reshape(1,-1)


        segs = list(seg_geno.keys())
        obs_copy_x, obs_copy_y =  data.copy_profiles_by_seg(segs, cells)
        x_diff = np.abs(obs_copy_x - latent_x).sum(axis=1)
        y_diff = np.abs(obs_copy_y - latent_y).sum(axis=1)
        return x_diff + y_diff
 
    def node_snv_cost(self,  v, cells, data, latent_genos):

        snvs = self.get_all_muts()
        latent_vafs =self.get_latent_vafs(latent_genos)
        lvafs = np.array(list(latent_vafs.values())).reshape(1,-1)
        snvs = list(latent_vafs.keys())
        #dims = cells by snvs
        obs_vafs = data.obs_vafs(cells, snvs)

        cell_scores = np.nansum(np.abs(obs_vafs - lvafs), axis=1)


        return cell_scores
    
    def pooled_node_snv_cost(self, v, cells, data, latent_genos):

 
        latent_vafs = self.get_latent_vafs(latent_genos)
        lvafs = np.array(list(latent_vafs.values()))
        snvs = list(latent_vafs.keys())
        #dims = cells by snvs
        obs_vafs = data.compute_vafs(cells, snvs)

        comp_df =  pd.DataFrame({'obs_vafs': obs_vafs, 'latent_vafs': lvafs, 'snvs': snvs})
        comp_df["node"] = v
        comp_df["ncells"] = len(cells)
        self.all_dfs[v]=comp_df
        cell_scores = np.nansum(np.abs(obs_vafs - lvafs))


        return cell_scores




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
    
    def node_likelihood(self, u, data, cells, vaf):

  
            snvs = list(vaf.keys())
            vafs = np.array([vaf[j] for j in snvs]).reshape(-1,1)
        
            return data.binomial_likelihood(cells, snvs,vafs )
    
    def compute_node_likelihoods(self, data, cells=None, lamb=0):
        if cells is None:
            cells = data.cells
        nodes = np.array(self.tree)
        node_cell_likes = []
        node_cna_scores = []
        for u in nodes:
            latent_genos = self.get_latent_genotypes(u)
            latent_vaf = self.get_latent_vafs(latent_genos)
            node_cell_likes.append(self.node_likelihood(u, data, cells,latent_vaf))
            node_cna_scores.append(self.node_cna_cost(u, cells, data, latent_genos))

        cell_scores = np.vstack(node_cell_likes)
        cell_cna_scores = np.vstack(node_cna_scores)
        cell_scores += lamb * cell_cna_scores
        return cell_scores, nodes


    def assign_cells_by_likelihood(self, data, cells=None, lamb=1000):
        cell_scores, nodes = self.compute_node_likelihoods(data, cells,lamb)
        node_assign = np.argmin(cell_scores, axis=0)
        obj =cell_scores.min(axis=0)
        obj = obj.sum()
        phi = {}
        for i, k in enumerate(node_assign):
            phi[i] = nodes[k]
        ca = CellAssign(phi, self.clones())
        self.cost = obj
        return obj, ca 

    def compute_likelihood(self, data,  cellAssign, lamb=0):



        self.node_cost = {}
        self.snv_node_cost = {}
        self.cna_node_cost = {}
        self.cost = 0
        for v in self.tree:
                cells = cellAssign.get_cells(v)
                if len(cells) ==0:
                    continue      

                else:
                    lat_genos = self.get_latent_genotypes(v)
                    self.snv_node_cost[v] = self.node_likelihood(v,data, cells, self.get_latent_vafs(lat_genos)).sum()
                    self.cna_node_cost[v] = self.node_cna_cost(v, cells, data,lat_genos ).sum()
                    self.node_cost[v]=  self.snv_node_cost[v] + lamb * self.cna_node_cost[v]
             
             


        self.cost = sum([score for _, score in self.node_cost.items()])

        return self.cost 


    def assign_cells(self, data, lamb=0, lamb_vaf= 1, cna_only=False):

        nodes = np.array(self.tree)
        if cna_only:
            cell_scores =0
            lamb=1
        else:
            cell_scores = np.vstack([self.node_snv_cost(v, data.cells, data, self.get_latent_genotypes(v)) for v in nodes])

        if lamb > 0:
            cell_cna_scores = np.vstack([self.node_cna_cost(v, data.cells, data, self.get_latent_genotypes(v)) for v in nodes])
            cell_scores = lamb_vaf*cell_scores + lamb*cell_cna_scores
        assignments = np.argmin(cell_scores, axis=0)
        self.cost = np.nansum(np.nanmin(cell_scores,axis=0)).sum()
        
        
        phi = {i: v for i,v in zip(data.cells, nodes[assignments])}
  
        # self.phi_to_cell_mapping(self.phi)
        return CellAssign(phi, self.clones()), self.cost, cell_scores, nodes   
        # self.cell_mapping = defaultdict(list)
        # for i, v in self.phi.items():
        #     self.cell_mapping[v].append(i)
        # self.cell_mapping = dict(self.cell_mapping)

    def filter_snvs(self, snvs_to_keep):
        snvs_to_remove = set(self.get_all_muts()) - set(snvs_to_keep)
        for s in snvs_to_remove:
            del self.psi[s]
            for v in self.tree:
                self.mut_mapping[v] = np.intersect1d(self.mut_mapping[v], snvs_to_keep).tolist()
                self.mut_loss_mapping[v] = np.intersect1d(self.mut_loss_mapping[v], snvs_to_keep).tolist()
                del self.genotypes[v][s]
        
            del self.mut_to_seg[s]

            self.seg_to_muts = {}
            for m, seg in self.mut_to_seg.items():
                if seg in self.seg_to_muts:
                    self.seg_to_muts[seg].append(m)
                else:
                    self.seg_to_muts[seg] = [m]

        
    def del_snv(self, j):
        seg = self.mut_to_seg[j]
        del self.mut_to_seg[j]
        self.seg_to_muts[seg].remove(j)
        node = self.psi[j]
        if j in self.mut_mapping[node]:
            self.mut_mapping[node].remove(j)
        
        for v in self.tree:
            del self.genotypes[v][j]
    def update_genotype(self, node, j, snv_tree):
        cna_geno = self.get_cna_genos()[self.mut_to_seg[j]]

        pres_nodes = []
        for u in sorted(nx.descendants(self.tree, node) | {node}):
            cn_geno = cna_geno[u]
            cn_state = cn_geno.to_tuple()
            pres_nodes.append(u)
            added = False
            for v in snv_tree.preorder():
                geno = snv_tree.genotypes[v][j]
                if (geno.x, geno.y) == cn_state and geno.z > 0:
                    self.genotypes[u][j] = genotype(*geno.to_tuple())
                    added = True
                    break 
            if added:
                break  # Early termination if condition is met
        else:
            for u in sorted(self.clones().difference(pres_nodes)):
                cn_state = cna_geno[u].to_tuple()
                self.genotypes[u][j] = genotype(*cn_state, 0, 0)

        self.psi[j] = node
    
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

    def compute_costs(self, data,  cellAssign, lamb=0, cna_only=False):
        if cna_only:
            lamb = 1
        # if  cellAssign is None:
        #     self.assign_cells(data, lamb)

        self.node_cost = {}
        self.cost = 0
        for v in self.tree:
                cells = cellAssign.get_cells(v)
                if len(cells) ==0:
                    continue
                latent_genos = self.get_latent_genotypes(v)
            
                if cna_only:
                    self.node_cost[v] =0
                else:
                    
                    self.node_cost[v]= self.snv_cost_func(v, cells, data, latent_genos).sum() 
                if lamb >0:
                    
                    self.node_cost[v]+=lamb*self.node_cna_cost(v, cells,data, latent_genos).sum()
             
        
        self.cost = sum([score for _, score in self.node_cost.items()])

        return self.cost 
    
    def compute_pooled_costs(self, data, cellAssign, lamb=0, cna_only=False):
        self.all_dfs = {}
        self.snv_cost_func = self.pooled_node_snv_cost
        cost = self.compute_costs(data,cellAssign, lamb, cna_only)
        self.snv_cost_func = self.node_snv_cost
        # self.comp_df = pd.concat([vals for v, vals in self.all_dfs.items()])
        # print(self.comp_df.head())

        return cost 

    




   #-------------------------- Save Methods ---------------------------------------#
    def draw(self, fname, cellAssign=None, mapping=None,segments=None, cmap='Set3'):
        if segments is not None:
            cna_genos = self.get_cna_genos()
        mut_count = {n : len(self.mut_mapping[n]) for n in self.mut_mapping}
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
                    labels[n] += "\n" + ",".join([str(cna_genos[s][n]) for s in segments if s in self.seg_to_muts])

        like_label = f"Segment {self.key}\n"
        tree = pgv.AGraph(strict=False, directed=False)
        tree.node_attr['style']='filled'
        segs = self.get_segments()
        segs = [str(ell) for ell in segs ]
        if self.cost is not None:
            total_like = np.round(self.cost)
            # tree.graph_attr["label"] = f"Objective: {total_like}\nSegments: {','.join(segs)}\nn={cellAssign.n} cells\nm={len(self.get_all_muts())} SNVs"
            tree.graph_attr["label"] = f"Objective: {total_like}\nSegments: {','.join(segs)}\nn={ncells} cells\nm={len(self.get_all_muts())} SNVs"

        tree
 
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
    
    
    # def get_node_cells(self, t):
    #     if t in self.cell_mapping:
    #         return self.cell_mapping[t]
    #     return []

    def get_leaves(self):
        leaves = [l for l in list(self.tree.nodes())
                  if self.tree.out_degree(l) == 0]
        return leaves

   

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
        gt_cell = pd.Series(self.phi).values
        pred_cell = pd.Series(obj.phi).values
        return self.ari(gt_cell, pred_cell)

    
    def compute_mut_ari(self,obj) -> float:
        #  gt_mut = self.get_mut_clusters()
         gt_mut =   pd.Series(self.get_psi())
         pred_mut = pd.Series(obj.get_psi())
         df = pd.concat([gt_mut, pred_mut], axis=1, keys=['gt', 'pred'])
        #  pred_mut = obj.get_mut_clusters()
      
         return self.ari(df["gt"].values, df["pred"].values)


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
                 'n_assigned': len(self.get_all_muts())
         
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

    def gt_latent_genos(self, u):
        lat_genos = self.get_latent_genotypes(u)
        latent_cna_genos = self.get_latent_cna_genos(lat_genos)
        seg_geno = {}
        for key, seg in self.mut_to_seg.items():
            if seg not in seg_geno:
               seg_geno[seg] = latent_cna_genos[key]
        latent_x, latent_y = [], []
        for seg,cna_geno in sorted(seg_geno.items()):
            latent_x.append(cna_geno.x)
            latent_y.append(cna_geno.y)
        
        latent_x = np.array(latent_x).reshape(-1,1)
        latent_y = np.array(latent_y).reshape(-1, 1)
        return latent_x, latent_y

    def cna_genotype_similarity(self, gt_phi, inf_ct, inf_phi):

        total_dist = 0
        # if len(gt_phi.get_all_cells() ^ inf_phi.get_all_cells()) ==0:

        
        all_u_dist = []
        for u in gt_phi.clones:
            gt_cells = gt_phi.get_cells(u)
            if len(gt_cells) ==0:
                continue
            gt_x, gt_y = self.gt_latent_genos(u)
            
            node_mapping = {}
            for i in gt_cells:
                node = inf_phi.phi[i]
                if node in node_mapping:
                    node_mapping[node].append(i)
                else:
                    node_mapping[node] = [i]
            for v, inf_cells in node_mapping.items():
                ncells = len(inf_cells)
                inf_x, inf_y = inf_ct.gt_latent_genos(v)
                x_dist =  np.abs(inf_x - gt_x)
                y_dist  = np.abs(inf_y - gt_y)
                geno_diff = (x_dist + y_dist).sum() 
                # if geno_diff > 0:
                #     print(f"gt node: {u} inf node: {v}")
                #     print(x_dist)
                #     print(y_dist)
    
                total_dist  += ncells*geno_diff
        
    
        mad = total_dist/ len(inf_phi.get_all_cells())
        return  mad
                    






        # self.mut2seg = {}
        # for v in self.genotypes:
        #     for s in self.genotypes[v]:
        #         for m in self.genotypes[v][s]:
        #             self.mut2seg[m]=s
        # self.mut2seg = {m: s  for v in self.genotypes for s, muts in self.genotypes[v].items() for m in muts}
        
        # self.seg2muts = {}
        # for m,s in self.mut2seg.items():
        #     if s in self.seg2muts:
        #         self.seg2muts[s].append(m)
        #     else:
        #         self.seg2muts[s] = [m]

    

    
