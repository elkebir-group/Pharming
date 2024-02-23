
import networkx as nx
import numpy as np
import pandas as pd 
from itertools import chain, combinations, permutations
from clonal_tree import ClonalTree
import pygraphviz as pgv
from genotype import genotype, CNAgenotype
import clonelib
from cell_mapping import CellAssign
from constrained_cell_assignment import ConstrainedCellAssign
from dataclasses import dataclass
import gurobipy as gp
from gurobipy import GRB
import pickle 
from copy import deepcopy
from enumerate import Enumerate


def draw(tree, fname):
    ptree = pgv.AGraph(strict=False, directed=True)
    ptree.add_edges_from(list(tree.edges))
    ptree.layout("dot")
    ptree.draw(fname)

def pickle_object(obj, file_path):
    """
    Pickle an object and save it to a file.

    Args:
    - obj: The object to pickle.
    - file_path: The file path where the pickled object will be saved.
    """
    with open(file_path, 'wb') as f:
        pickle.dump(obj, f)

def load_pickled_object(file_path):
    """
    Load a pickled object from a file.

    Args:
    - file_path: The file path from which to load the pickled object.

    Returns:
    - The unpickled object.
    """
    with open(file_path, 'rb') as f:
        obj = pickle.load(f)
    return obj


@dataclass
class STISol:
    cost: float
    ct: ClonalTree
    phi: CellAssign

    def png(self, fname, segments=None):
        self.ct.draw(fname,self.phi, segments=segments )


class STI:
    '''
    a class to solve the segment tree inference problem
    
    S: a nx:Digraph represnting the CNA tree
    k: int representing the number of SNV clusters
    r: int representing the number of cell clusters
    seed: int representing the random number seed 

    '''
    def __init__(self, S, T_m, delta, lamb1 =5, lamb2=1000) -> None:
    

        nodes = list(S)
        if len(nodes) == 0:
                 raise ValueError("CNA Tree S is not valid!")
        else:
            self.S = S 
    
        if len(delta) != len(T_m):
                 raise ValueError("Mutation cluster tree must match the number of DCFs")
        
        else:
            self.T_m = T_m
            self.delta = delta

        self.lamb1 = lamb1
        self.lamb2  =lamb2

 
        self.S_root = [n for n in self.S if S.in_degree[n]==0][0]
        self.cn_states = {}
        self.k = len(self.delta)
        for i, u in enumerate(self.S):
            self.cn_states[u] = self.k +i + 1 
        
        self.cn_states_inv = {val: key for key,val in self.cn_states.items()}



        if len(list(self.S)) > 0:
            snv_tree_edgelists = clonelib.get_genotype_trees(list(self.S.edges))
        else:
            snv_tree_edgelists = []
         
            root = nodes[0]
            for i in range(2):
                if i == 0:
                    snv_tree_edgelists.append( [((root[0], root[1], 0, 0), 
                                                 (root[0], root[1], 1, 0))])
                else:
                    snv_tree_edgelists.append( [((root[0], root[1], 0, 0),
                                                  (root[0], root[1], 0, 1))])

  
        self.T_SNVs = [nx.DiGraph(edges) for edges in snv_tree_edgelists]
        self.T_SNV_groups, self.group_desc = self.group_snv_trees()

    

    
    @staticmethod
    def to_inverse_dict(mydict):
        '''
        reverses an assignment/mapping dictionary
        '''
        rev_dict = {}
        for key, v in mydict.items():
            # for v in val:
                if v not in rev_dict:
                    rev_dict[v] = [key]
                else:
                    rev_dict[v].append(key) 
        return rev_dict
    
    @staticmethod
    def to_mapping_dict(mydict):
        '''
        reverses an inverse dictionary 
        '''
        rev_dict = {}
        for key, val in mydict.items():
            for v in val:
                rev_dict[v] = key

        return rev_dict
    
 

    
    def convert_to_clonal_tree(self, t, j):
        '''
        Converts a networkx digraph SNV tree t labeled by SNV genotypes (x,y,x_bar, y_bar)
        a ClonalTree for SNV j
        '''
        geno_dict = {}
        relabel = {}
        for i, v in enumerate(t):
            geno_dict[i] = {j : genotype(*v) }
            relabel[v] = i

        t_seg_to_muts = {self.ell: [j]}
        t_copy = nx.relabel_nodes(t, relabel)
        ct = ClonalTree(t_copy, geno_dict, t_seg_to_muts)
        return ct, {n: v for v,n in relabel.items()}
    

    @staticmethod
    def solve(model,vars, threads=1, timelimit=120):
        '''
        solves a gurobi MILP model to either optimality or time out (timelimit)
        returns the objective value and specified value of the decision variables vars
        '''
        model.Params.Threads = threads
        model.Params.TimeLimit = timelimit
        model.Params.OutputFlag = 0
        model.Params.LogToConsole = 0
     
        model.optimize()
        solutions = []
        if model.Status == GRB.OPTIMAL or (model.Status==GRB.TIME_LIMIT and model.SolCount>0):
            score = model.objVal
            for v in vars:
                solutions.append(model.getAttr('X', v))
 
     
            return score, solutions
        
        else:
             print("warning: model infeasible!")
             return np.Inf, solutions 


    def group_snv_trees(self):
        '''
        Given a set (self.T_SNVs) of SNV trees, break the trees into sets
        based on the split copy number state and the CNA states that SNV is present
        
        Returns a list of lists of GenotypeTree ids that are in the same set 
        '''
  

        def powerset(iterable):
            s = list(iterable)
            return list(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))

        clusters = {}
        index = 0
        for v in self.S:
            children = self.S.successors(v)
            for states in powerset(children):
                clusters[index] = {'sscn': (v[0], v[1]), 'children': set(state for state in  states)}
                index += 1
        
        tree_clusters = [[] for _ in clusters]
        for tree in self.T_SNVs:
            for u in tree:
                parent = list(tree.predecessors(u))
                if len(parent) > 0:
                    parent = genotype(*parent[0])
                    child = genotype(*u)
                
                    if parent.cna_eq(child):
                        sscn = child.to_CNAgenotype().to_tuple()
                        desc = list(tree.successors(u))
                        if len(desc) ==0:
                            children_states = set({})
                        else:
                            children_states = set((v[0], v[1]) for v in desc)
                        
                        break
            for index, clust in clusters.items():
                if clust['sscn'] == sscn and clust['children'] == children_states:
                    tree_clusters[index].append(tree)
                    break
            
        return tree_clusters, clusters

    
    
    def compute_tree_cost(self, j, dcf, tree):
 
          
        ct, mapping = self.convert_to_clonal_tree(tree,j)

        # ct.draw("test/snv_tree.png")

        # if j == 290 and dcf == 0.179:
   
        #     ct.draw("test/ct.png", segments=[self.ell])
        _, u, _, _ =  ct.get_split_nodes(j)
        # ca, obj, cell_scores,nodes = ct.assign_cells(self.data, self.lamb1)
        cell_scores, nodes = ct.compute_node_likelihoods(self.data,self.data.cells, lamb=self.lamb1)
        node_assign = np.argmin(cell_scores, axis=0)
        obj =cell_scores.min(axis=0)
        obj = obj.sum()
        phi = {}
        for i, k in enumerate(node_assign):
            phi[i] = nodes[k]
        ca = CellAssign(phi, ct.clones())
        # if j == 290 and dcf == 0.179:
        #     ct.draw("test/ct.png", ca,  segments=[self.ell])
        alt = self.data.var[:, j].sum()
        total = self.data.total[:,j].sum()
        posterior_dcf = ct.posterior_dcf(j, dcf, alt, total, self.cn_props)
        # dcfs = ct.compute_dcfs(ca)
        # obj += self.lamb2*np.abs(dcf - dcfs[u])
        
        return obj, ct, ca, posterior_dcf
    
    
    def compute_snv_cluster_tree_cost(self):
        all_costs = []
        cst = {}
        # delta_hat = {}
        tree_assign = {}
        # num_groups = len(self.T_SNV_groups)
        '''
        Enumerate costs for each SNV cluster and each tree/group
        '''
        for q, dcf in self.delta.items():
            cst[q] = {}
            for g, trees in enumerate(self.T_SNV_groups):
                cst[q][g] = {}
                tree_assign[g] = {}
              
                    # if j==1108 and q ==3 and g==2:
                    # if j==1108 and q ==4 and g==:
                    #     print("here")
                for j in self.snvs:
                    jg_cost = np.Inf
                    for p,t in enumerate(trees):
                     
                        cost, ct, ca,post_dcf = self.compute_tree_cost(j, dcf, t)
                        
                        all_costs.append([j,q,dcf,g,p,cost, post_dcf])
                        if cost < jg_cost:
                            cst[q][g][j]  = cost 
                            jg_cost = cost
                            # delta_hat[j,q,g] = np.abs(obs_dcf - dcf)
                            tree_assign[g][j] = (ct,ca)
        
        df = pd.DataFrame(all_costs, columns=["snv", "snv_clust", "dcf", "group", "tree", "cost", "posterior_dcf"])
        return cst, df, tree_assign
   


    def identify_snv_cluster_trees(self,T):
        '''
        Given a networkx tree labeled by (q, (x,y)),
        return a set of tuples (q, g) where the SNV clusters are introduced
        '''
  
        root = (-1, (self.S_root))
        sscn = {}
        children = {}
        for u in nx.dfs_preorder_nodes(T, source=root):
            u_q, u_cn = u
            for  v_q, v_cn in T.successors(u):
                  if u_q != v_q and u_cn == v_cn:
                      sscn[v_q] = u_cn
                      children[v_q] = []
                      for w_q, w_cn in nx.dfs_preorder_nodes(T, source = (v_q, v_cn)):
                          if w_cn != v_cn and w_cn not in children[v_q]:
                              children[v_q].append(w_cn)
        groups = {}
        for q in sscn:
            for g in self.group_desc:
                g_sscn = self.group_desc[g]['sscn']
                g_children = self.group_desc[g]['children']
                if sscn[q] == g_sscn and set(children[q]) == set(g_children):
                    if g in groups:
                  
                        groups[g].append(q)
                    else:
                        groups[g] = [q]
        return groups 
    
                          
    def cluster_snvs(self, df, cluster_groups, tree_assign):
        
        all_dfs = []
        for g in cluster_groups:
            temp = df[df['group']==g]
            temp = temp[temp['snv_clust'].isin(cluster_groups[g])]
            all_dfs.append(temp)
        df_filt= pd.concat(all_dfs)

        grouped = df_filt.groupby('snv')

        # Define a function to sort and select the first row of each group
        def sort_select_first(group):
            sorted_group = group.sort_values(by=['cost', 'posterior_dcf'], ascending=[True, False])
            return sorted_group.iloc[0]

        # Apply the function to each group and concatenate the results
        result = grouped.apply(sort_select_first)

        # Reset index if needed
        result = result.reset_index(drop=True)
        result['snv']  = result['snv'].astype(int)
        result['snv_clust']  = result['snv_clust'].astype(int)
        result['group']  = result['group'].astype(int)
        #     # Group by SNV and find the index of the row with minimum cost in each group
        # min_cost_indices = df_filt.groupby('snv')['cost'].idxmin()

        # # Use the indices to select the rows with minimum cost for each SNV
        # min_cost_rows = df_filt.loc[min_cost_indices]
        psi = dict(zip(result['snv'], result['snv_clust']))
        snv_group = dict(zip(result['snv'], result['group']))

        omega = {}
        for j,g in snv_group.items():
            omega[j] = tree_assign[g][j]
             


        return self.to_inverse_dict(psi), omega


    def construct_segment_tree(self, T, alpha_inv, omega):
        genotypes = {}
        seg_to_snvs = {self.ell : self.snvs}
        ct_tree = nx.DiGraph()
        r = self.cn_states[self.S_root]
        ct_tree.add_node(r)


        mapping = {(-1,self.S_root): r }
        for u in nx.dfs_preorder_nodes(T, source=(-1,self.S_root)):
            u_q, u_cn = u 
            parent = mapping[u]
            for v in T.successors(u):
                v_q, v_cn = v 
                if u_q != v_q and u_cn == v_cn:
                    ct_tree.add_edge(parent, v_q)
                    mapping[v] = v_q
                else:
                    ct_tree.add_edge(parent, self.cn_states[v_cn])
                    mapping[v] = self.cn_states[v_cn]
        # for n in ct_tree:
        #     if n != r:
          
        rev_mapping = {val: key for key,val in mapping.items()}

        for n in ct_tree:
            # if n ==10:
            #     print("here")
            cn_state = rev_mapping[n][1]
            cna_geno = CNAgenotype(*cn_state)
            genotypes[n] = {}
            anc = list(nx.ancestors(ct_tree, source=n))
            anc.append(n)
            present_clusters = [q for q in alpha_inv if q in anc]
            pres_snvs = set([j  for q in present_clusters for j in alpha_inv[q]])

            for j in pres_snvs:
                added = False
                ct = omega[j][0]
                cna_genos = ct.get_cna_genos()[self.ell]
                for v in ct.preorder():
                    geno = ct.genotypes[v][j]
                    if cna_genos[v] == cna_geno and geno.z > 0:
                        genotypes[n][j] = genotype(*geno.to_tuple())
                        added = True
                        break 
                if not added:
                    genotypes[n][j] =genotype(*cn_state, 0, 0)
            
            for j in set(self.snvs).difference(pres_snvs):
                genotypes[n][j] = genotype(*cn_state, 0,0)
        return ClonalTree(ct_tree, genotypes, seg_to_snvs)
                    



        # for q, snvs in alpha_inv.items():
        #     desc_nodes  = list(nx.dfs_preorder_nodes(ct_tree, source=q))
        #     for n in ct_tree:
        #         cn_state = rev_mapping[n][1]
        #         cna_geno = CNAgenotype(*cn_state)
        #         if n not in desc_nodes:
             
        #             for j in self.snvs:
        #                 genotypes[n][j] = genotype(*cn_state, 0,0)
        #         else:
        #             for j in self.snvs:
        #                 if j in snvs:
                   
        #                 else:
        #                     genotypes[n][j] = genotype(*cn_state, 0,0)


    @staticmethod
    def solve(model,vars, threads=1, timelimit=120):
        model.Params.Threads = threads
        model.Params.TimeLimit = timelimit
     
        model.optimize()
        solutions = []
        if model.Status == GRB.OPTIMAL or (model.Status==GRB.TIME_LIMIT and model.SolCount>0):
            score = model.objVal
            for v in vars:
                solutions.append(model.getAttr('X', v))
 
     
            return score, solutions
        
        else:
             print("warning: model infeasible!")
             return np.Inf, solutions 

    def assign_cell_clusters(self,ct, snv_clusters):

    
        delta =self.delta
        # _, _, cell_scores, nodes = ct.assign_cells(self.data, self.lamb1)
        cell_scores, nodes = ct.compute_node_likelihoods(self.data, self.data.cells, self.lamb1)
        nodes = nodes.tolist()
        valid_nodes = list(self.cn_states_inv.keys()) + snv_clusters
        clust_score = {}
        for i in range(cell_scores.shape[1]):
            for q,n in enumerate(nodes):
                if n in valid_nodes:
                    clust_score[i,n] = cell_scores[q,i]
  
            # print(f"{q}: {self.data.compute_vafs(cell_clust, list(ct.psi.keys()))}")
            # obs_vafs = self.data.obs_vafs(cell_clust, list(ct.psi.keys()))
            # print(obs_vafs)
            # clust = cell_scores[:, cell_clust].sum(axis=1) 

        #     clust_score.append(clust)
        #     for i, u in enumerate(nodes):
        #         clust_cost[q, u] = clust[i]
        
        # cs = np.vstack(clust_score)
        # best_clust = cs.argmin(axis=1)
        # objval = cs.min(axis=1).sum()
        phi ={}
        # for p, i in zip(beta_inv, best_clust):
        #     phi[p] = nodes[i]
        model = gp.Model("MIP")
        cells = np.arange(cell_scores.shape[1])
        cell_assign = [(i,u) for i in cells  for u in valid_nodes]
        x = model.addVars(cell_assign, vtype = GRB.BINARY )
        z = model.addVars(delta.keys(), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS )


        # model.setObjective(gp.quicksum(clust_score[i,u]*x[i,u] for i,u in cell_assign) + 
        #                    self.lamb2*gp.quicksum(z[u] for u in delta), gp.GRB.MINIMIZE)
        
        model.setObjective(gp.quicksum(clust_score[i,u]*x[i,u] for i,u in cell_assign))
        #every cell is assigned to exactly 1 node
        model.addConstrs(gp.quicksum(x[i,u] for u in valid_nodes)==1 for i in cells)



        #constraints for absolute values
        for u in delta:
   
            
                descendants = set(ct.preorder(u))
                desc = descendants.intersection(set(valid_nodes))
                # model.addConstr( (1/self.data.N)*gp.quicksum(x[i,u] for i in cells for u in desc ) - delta[u] <= z[u] )
                # model.addConstr( (1/self.data.N)*gp.quicksum(x[i,u] for i in cells for u in desc) - delta[u] >= -1*z[u] )
                model.addConstr( (1/self.data.N)*gp.quicksum(x[i,u] for i in cells for u in desc ) >= delta[u])

        
        
        objval, sol = self.solve(model, [x, z])
        if objval == np.Inf:
            return None, None
        x_sol = sol[0]
        z_sol = sol[1]
  
    
        for i,u in cell_assign:
    
                if x_sol[i,u] > 0.5:
                   phi[i] = u 
    
        ca = CellAssign(phi, ct.clones())
        # dcf = self.compute_dcfs(ct.get_all_muts()[0], ct, ca)
    
        return objval, ca

    def assign_genotypes(self, segtree, ca, tree_assign, snv_clusters):
        T=  deepcopy(segtree)
        cna_genos = segtree.get_cna_genos()[self.ell]
        cell_counts  = ca.get_cell_count()
        
        def get_group(q):    
            sscn = cna_genos[q].to_tuple()
            children_states = []
            for u in segtree.preorder(q):
                if u !=q:
                    cn_state = cna_genos[u].to_tuple()
                    if cn_state not in children_states and cn_state != sscn:
                        children_states.append(cn_state)
            for g in self.group_desc:
                if sscn == self.group_desc[g]['sscn']:
                    if set(children_states) == self.group_desc[g]['children']:
                        return g
  
                
        def get_vafs(q,j,snv_tree):
            vafs = {}
            
            for u in segtree.clones():
                if cell_counts[u] > 0:
                    if nx.has_path(segtree.tree, u, q) and u !=q:
                        vafs[u] =0
                        continue
                    if nx.has_path(segtree.tree, q,u):
                        cn_state = cna_genos[u].to_tuple()
                        added = False
                        for v in snv_tree.preorder():
                            geno = snv_tree.genotypes[v][j]
                            if (geno.x, geno.y) == cn_state and geno.z > 0:
                                vafs[u] = geno.vaf
                                added = True
                                break 
                        if not added:
                            vafs[u] =0 #loss occurred 
                    else:
                        vafs[u] =0
            return vafs 
        
        psi = segtree.get_psi()
        for j in self.snvs:
            best_cost = np.Inf
            vafs = {}
            for q in snv_clusters:
                #get group of q 
                g = get_group(q)
                #get the optimal snv tree for snv j in group of q
                snv_tree = tree_assign[g][j][0]
                vafs = get_vafs(q, j, snv_tree)
           
                #get latent vaf of each node in nonempty clones
                cost = 0
                for u in ca.clones:
                    if cell_counts[u] > 0:
                        cells = ca.get_cells(u)
                        cell_costs = self.data.binomial_likelihood(cells, [j], vafs[u])
                        cost += cell_costs.sum()
                if cost < best_cost:
                    new_node = q 
                    best_cost = cost 
            if psi[j] != new_node:
                # snv_tree.draw("test/snv_tree.png", segments=[self.ell])
                T.update_genotype(new_node, j, snv_tree)
        
        T.update_mappings()

        return T 
                

        
    def fit(self, data, segment):
        self.data = data
        self.ell = segment
        self.snvs = self.data.seg_to_snvs[segment]
        refinements = Enumerate( self.T_m, self.S,).solve()
        
        #TODO: fix to account for errors in observed cn states
        self.cn_props = self.data.cn_proportions(self.ell)

        cost, df, tree_assign = self.compute_snv_cluster_tree_cost()
        # df.to_csv("test/all_costs.csv", index=False)
        # pickle_object(cost, "test/costs.pkl")
        # pickle_object(tree_assign, "test/tree_assign.pkl")
        results = []
        opt_cost = np.Inf 
        print(f"Total refinements: {len(refinements)}")
        for f, T in enumerate(refinements):
            print(f"Starting refinement {f}...")
            cluster_groups = self.identify_snv_cluster_trees(T)
            alpha_inv, omega = self.cluster_snvs(df, cluster_groups, tree_assign)
            snv_clusters = list(alpha_inv.keys())
            
            segment_tree= self.construct_segment_tree(T, alpha_inv, omega)
            best_phi = None
            cost  = np.Inf

            for _ in range(10):
                ca_cost, ca = self.assign_cell_clusters(segment_tree, snv_clusters )
                if ca_cost is None:
                    break
                #TODO: come up with a bound to not explore viable solutions
                # segment_tree.draw("test/segtree.png", ca, segments=[self.ell])
                segment_tree = self.assign_genotypes(segment_tree, ca, tree_assign, snv_clusters)
                # segment_tree.draw("test/segtree.png", ca, segments=[self.ell])
                updated_cost = segment_tree.compute_likelihood(self.data, ca, lamb=self.lamb1)
                if updated_cost < cost:
                    best_segtree = deepcopy(segment_tree)
                    best_phi = deepcopy(ca)
                    cost  = updated_cost
                    if cost < opt_cost:
                        opt_cost = cost 
                else:
                    break
            print(f"Ending refinement {f}: cost: {cost} opt_cost: {opt_cost}")
            if best_phi is not None:
                results.append(STISol(cost, best_segtree, best_phi))

        return sorted(results, key= lambda x: x.cost)



        

    
   
    
    
    
        

  


    





        

    

        

    # def update_cell_clusters(self, ct, delta):
    #     cell_cost = {}
    #     cells = self.data.cells 
    #     N = len(cells)
    #     _, _, cell_scores, nodes = ct.assign_cells(self.data, self.lamb1)
    #     nodes = nodes.tolist()
    #     for j,u in enumerate(nodes):
    #         for i, c in enumerate(cells):
    #             cell_cost[c, u] = cell_scores[j,i]
        
    #     model = gp.Model("MIP")
    #     cell_assign = [(i,u) for i in cells for u in nodes]
    #     x = model.addVars(cell_assign, vtype = GRB.BINARY )
    #     z = model.addVars(delta.keys(), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS )
        
    #     #indicates if node u has any cells assigned 
    #     y = model.addVars(nodes, vtype=GRB.BINARY)



    #     model.setObjective(gp.quicksum(cell_cost[i,u]*x[i,u] for i,u in cell_assign) + 
    #                        self.lamb2*gp.quicksum(z[u] for u in delta), gp.GRB.MINIMIZE)
        
    #     #a node can only have assigned cells if it is used 
    #     for u in nodes:
    #         model.addConstrs(x[i,u] <= y[u] for i in cells)
        
    #     model.addConstrs((gp.quicksum(x[i,u] for u in nodes)==1) for i in cells )
        
    #     #if a node is used it must have assigned cells (non-empty cell cluster)
    #     for u in nodes:
    #         model.addConstr(gp.quicksum(x[i,u] for i in cells) >= y[u])

    #     #exactly r nodes must be used 
    #     model.addConstr(gp.quicksum(y[u] for u in nodes)==self.r)
    #     for u in delta:
    #         descendants = list(ct.preorder(u))
    #         model.addConstr( (1/N)*gp.quicksum(x[i,u] for i in cells for u in descendants ) - delta[u] <= z[u] )
    #         model.addConstr( (1/N)*gp.quicksum(x[q,u] for q in cells  for u in descendants ) - delta[u] >= -1*z[u] )
    #     objval, sol = self.solve(model, [x,y])
    #     x_sol, y_sol = sol[0], sol[1]
    #     cluster_id = 0

    #     phi_inv = {}
    #     for u in nodes:
    #         if y_sol[u] > 0.5:
    #             phi_inv[u] = cluster_id 
    #             cluster_id += 1
    #     beta_inv = {q: [] for q in range(self.r)}
     
    #     for i,u in cell_assign:
    #         if x_sol[i,u] > 0.5:
    #             beta_inv[phi_inv[u]].append(i)
        
    #     phi = {val: key for key,val in phi_inv.items()}
    #     return objval,  beta_inv, phi
            


        


    # def assign_cell_clusters(self,ct, beta_inv, delta,allow_multiple=True):
    #     # if 78 in ct.get_all_muts():
    #     #     ct.draw("test/snv_tree78.png")
    #     #     print("here")
    #     clust_cost = {}
    #     _, _, cell_scores, nodes = ct.assign_cells(self.data, self.lamb1, lamb_vaf=5)
    #     # cell_scores, nodes = ct.compute_node_likelihoods(self.data, lamb=100)
    #     nodes = nodes.tolist()
    #     clust_score = []
    #     for q, cell_clust in beta_inv.items():
  
    #         # print(f"{q}: {self.data.compute_vafs(cell_clust, list(ct.psi.keys()))}")
    #         # obs_vafs = self.data.obs_vafs(cell_clust, list(ct.psi.keys()))
    #         # print(obs_vafs)
    #         clust = cell_scores[:, cell_clust].sum(axis=1) /len(cell_clust)

    #         clust_score.append(clust)
    #         for i, u in enumerate(nodes):
    #             clust_cost[q, u] = clust[i]
        
    #     cs = np.vstack(clust_score)
    #     best_clust = cs.argmin(axis=1)
    #     objval = cs.min(axis=1).sum()
    #     phi ={}
    #     # for p, i in zip(beta_inv, best_clust):
    #     #     phi[p] = nodes[i]
    #     model = gp.Model("MIP")
    #     cell_assign = [(q,u) for q in beta_inv for u in nodes]
    #     x = model.addVars(cell_assign, vtype = GRB.BINARY )
    #     z = model.addVars(delta.keys(), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS )


    #     model.setObjective(gp.quicksum(clust_cost[q,u]*x[q,u] for q,u in cell_assign) + 
    #                        self.lamb2*gp.quicksum(z[u] for u in delta), gp.GRB.MINIMIZE)
        
        
    #     #every cluster is assigned to exactly 1 node
    #     model.addConstrs(gp.quicksum(x[q,u] for u in nodes)==1 for q in  beta_inv)



    #     #TODO: set up constraints for absolute values
    #     for u in delta:
    #         descendants = list(ct.preorder(u))
    #         model.addConstr( (1/self.data.N)*gp.quicksum(len(beta_inv[q])*x[q,u] for q in beta_inv for u in descendants ) - delta[u] <= z[u] )
    #         model.addConstr( (1/self.data.N)*gp.quicksum(len(beta_inv[q])*x[q,u] for q in beta_inv for u in descendants ) - delta[u] >= -1*z[u] )
        
        
    #     if not allow_multiple:
    #         #allow at most 1 cluster 
    #         for u in nodes:
    #             model.addConstr(gp.quicksum(x[q,u] for q in beta_inv) <=1)
        
    #     objval, sol = self.solve(model, [x, z])
    #     x_sol = sol[0]
    #     z_sol = sol[1]
  
    
    #     for q,u in cell_assign:
    
    #             if x_sol[q,u] > 0.5:
    #                phi[q] = u 
        
    #     cell_phi = {}
    #     for q in phi:
    #         for i in beta_inv[q]:
    #             cell_phi[i] = phi[q]
    #     ca = CellAssign(cell_phi, ct.clones())
    #     dcf = self.compute_dcfs(ct.get_all_muts()[0], ct, ca)
    
    #     return objval, phi, dcf


    def compute_snv_tree_cost(self, j, dcf, tree, beta_inv):
 
          
        ct, mapping = self.convert_to_clonal_tree(tree,j)

        # ct.draw("test/snv_tree.png")


        _, u, _, _ =  ct.get_split_nodes(j)
        obj, phi ,dcf = self.assign_cell_clusters(ct, beta_inv, delta={u: dcf})
        
        return obj, ct, phi, dcf

 

        
                
        # return v_to_dcf(vaf[0], F, cn_prop)
        # for state in states:
        #     cells= None
        #     # for s in cn_prop:
        #     #     if s == state:
        #     #         cn_prop[s] =1
        #     #     else:
        #     #         cn_prop[s] = 0
        #     # cells = cells_by_state[state]
        #    
        #     
        
        # return dcfs
        

    
    def update_SNVclusters(self, delta, beta_inv ):
        all_costs = []
        cst = {}
        delta_hat = {}
        tree_assign = {}
        num_groups = len(self.T_SNV_groups)
        '''
        Enumerate costs for each SNV cluster and each tree/group
        '''
        for j in self.snvs:
            tree_assign[j] = {}
            for q, dcf in delta.items():
                for g, trees in enumerate(self.T_SNV_groups):
                    g_cost = np.Inf
                    # if j==1108 and q ==3 and g==2:
                    # if j==1108 and q ==4 and g==:
                    #     print("here")
                    for p,t in enumerate(trees):
                        if j ==78 and q==4 and g==4:
                            print("here")
                        cost, ct, tree_phi, obs_dcf = self.compute_snv_tree_cost(j, dcf, t, beta_inv)
                        
                        # tphi =self.to_mapping_dict(tree_phi)
                        # tphi = {}
                        # for k, u in tree_phi.items():
                        #     for i in beta_inv[k]:
                        #         tphi[i] = u
                        # ca = CellAssign(tphi, ct.clones())
                        # dc, _ = ct.compute_dcfs(ca, self.ell)
                        # print(f"{q}: {dcf} vs {dc}")
                        all_costs.append([j,q,dcf,g,p,cost])
                        if cost < g_cost:
                            cst[(j,q,g)] = cost 
                            g_cost = cost
                            delta_hat[j,q,g] = np.abs(obs_dcf - dcf)
                            tree_assign[j][g] = (ct,self.to_inverse_dict(tree_phi))
        
        df = pd.DataFrame(all_costs, columns=["snv", "snv_clust", "dcf", "group", "tree", "cost"])
        df.to_csv("test/all_costs.csv", index=False)
        pickle_object(cst, "test/costs.pkl")
        pickle_object(tree_assign, "test/tree_assign.pkl")

        cst = load_pickled_object("test/costs.pkl")
        tree_assign = load_pickled_object("test/tree_assign.pkl")
        model = gp.Model("MIP")

        clust_group_assign = [(j,q,g) for j in self.snvs for q in range(self.k) for g in range(num_groups)]
        clust_assign = [(j,g) for j in self.snvs for g in range(num_groups)]
        group_assign =  [(q,g) for q in range(self.k) for g in range(num_groups)]
        snv_clust = [(j,q) for j in self.snvs for q in range(self.k)]

        # for j in [1108, 1309,1410, 1745]:
        #     for q, g in group_assign:
        #         print(f"{j}, {q}, {g}: {cst[j,q,g]}")
        x = model.addVars(clust_group_assign, vtype = GRB.BINARY )
        
        #assignment of cluster to group
        y =model.addVars(group_assign, vtype=GRB.BINARY)
        
        #assignment of SNV to group
        z = model.addVars(clust_assign, vtype = GRB.BINARY )
        
        #assignment of SNV to a cluster 
        w = model.addVars(snv_clust, vtype= GRB.BINARY)
       
  
         #every snv is assigned to exactly 1 cluster and 1 group
        model.addConstrs((gp.quicksum(x[j,q,g] for q in range(self.k) for g in range(num_groups)) == 1) for j in self.snvs)
        
        #every snv cluster q has at least 1 SNV assigned (clusters are non-empty)
        model.addConstrs((gp.quicksum(w[j,q] for j in self.snvs) >= 1) for q in range(self.k))

        #every snv is assinged to exactly 1 cluster
        model.addConstrs((gp.quicksum(w[j,q] for q in range(self.k)) ==1) for j in self.snvs)

        #every snv is assigned to one group
        model.addConstrs((gp.quicksum(z[j,g] for g in range(num_groups)) ==1) for j in self.snvs)

        #every snv cluster is assigned to one group
        model.addConstrs((gp.quicksum(y[q,g] for g in range(num_groups)) ==1) for q in range(self.k))

        model.addConstrs((gp.quicksum(w[j,q] for j in self.snvs)>=1) for q in range(self.k)) 


        #TODO: add pairwise compatible of selected groups 

        for j, q in snv_clust:
            for g in range(num_groups):
            # If cluster q is not assigned to group g, then SNV j in cluster q cannot be assigned to group g
                 model.addConstr(x[j, q, g] <= y[q, g])
        

        for j, g in clust_assign:
            for q in range(self.k):
            # If cluster q is not assigned to group g, then SNV j in cluster q cannot be assigned to group g
                 model.addConstr(x[j, q, g] <= z[j, g])
        
        for j, q in snv_clust:
            for g in range(num_groups):
            # If cluster q is not assigned to group g, then SNV j in cluster q cannot be assigned to group g
                 model.addConstr(x[j, q, g] <= w[j, q])
        # for j,g in clust_assign:
        #         model.addConstr(gp.quicksum(x[j,q,g] for q in range(self.k)) == z[j,g])
        # for j,q in snv_clust:
        #     model.addConstr(gp.quicksum(x[j,q,g] for g in range(num_groups))  == w[j,q])
        
        #every snv assigned to snv cluster q must all have the same group assignment 
   

        #minimize cost
        model.setObjective(gp.quicksum(cst[j,q,g]*x[j,q,g] for j,q,g in clust_group_assign), gp.GRB.MINIMIZE)

        objval, sols = self.solve(model, [x,y,w])

        x, y,  w = sols[0], sols[1], sols[2]
        alpha ={}
        omega = {}
        clust_to_group= {}
        for j,q,g in clust_group_assign:
            if x[j,q,g] > 0.5:
                alpha[j] = q 
                omega[j] = tree_assign[j][g]
        
        for q,g in group_assign:
            if y[q,g] > 0.5:
                clust_to_group[q]=g
        
        return self.to_mapping_dict(alpha), omega, clust_to_group




        

        
    @staticmethod
    def scalar_obj(dcf, snvs, omega, beta_inv, n):
        obj = 0
        for j in snvs:
            ct, tree_phi= omega[j]
           
            _, u, _, _ = ct.get_split_nodes(j)
            desc = set(ct.preorder(u))
            desc = desc.intersection(tree_phi.keys())
            desc_frac = sum(len(beta_inv[q])  for u in desc for q in tree_phi[u])/n
            obj += np.abs(desc_frac - dcf)
  
        return obj
            



    def update_tree(self, clonal_trees, beta_inv, delta):
        
     
        phi_list = []
        obj_vals = []
        for ct in clonal_trees:
            obj, phi = self.assign_cell_clusters(ct, beta_inv, delta, allow_multiple=False)
            obj_vals.append(obj)
            phi_list.append(phi) 
        
        index = obj_vals.index(min(obj_vals))
        return  clonal_trees[index], phi_list[index]
        

            
            




       




