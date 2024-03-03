
import networkx as nx
import numpy as np
import pandas as pd 
from itertools import chain, combinations, permutations
from clonal_tree import ClonalTree
import pygraphviz as pgv
from genotype import genotype, CNAgenotype
import clonelib
from cell_mapping import CellAssign
import gurobipy as gp
from gurobipy import GRB
import pickle 
from copy import deepcopy
from enumerate import Enumerate
from solution import Solution



def draw(tree, fname):
    ptree = pgv.AGraph(strict=False, directed=True)
    ptree.add_edges_from(list(tree.edges))
    ptree.layout("dot")
    ptree.draw(fname)




# @dataclass
# class Solution:
#     cost: float
#     ct: ClonalTree
#     phi: CellAssign

#     def png(self, fname, segments=None):
#         self.ct.draw(fname,self.phi, segments=segments )


class STI:
    '''
    a class to solve the segment tree inference problem
    
    S: a nx:Digraph represnting the CNA tree
    k: int representing the number of SNV clusters
    r: int representing the number of cell clusters
    seed: int representing the random number seed 

    '''
    def __init__(self, S, Tm_edges, delta, lamb1 =5, niter=10) -> None:
    
        T_m = nx.DiGraph(Tm_edges)
        nodes = list(S)
        if len(nodes) == 0:
                 raise ValueError("CNA Tree S is not valid!")
        else:
            self.S = S 
    
        if len(delta) != len(T_m):
                 raise ValueError("Mutation cluster tree must \
                                  match the number of DCFs")
        
        else:
            self.T_m = T_m
            self.delta = delta

        self.lamb1 = lamb1
     

        self.max_iterations = niter
        self.S_root = [n for n in self.S if S.in_degree[n]==0][0]
        self.cn_states = {}
        self.k = max(self.T_m)
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

        self.cn_delta = {}

    

    
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

        _, u, _, _ =  ct.get_split_nodes(j)

        obj, ca  = ct.assign_cells_by_likelihood(self.data,self.data.cells, lamb=self.lamb1)


        alt = self.data.var[:, j].sum()
        total = self.data.total[:,j].sum()
        posterior_dcf = ct.posterior_dcf(j, dcf, alt, total, self.cn_props)
    
        
        return obj, ct, ca, posterior_dcf
    
    
    def compute_snv_cluster_tree_cost(self):
        all_costs = []
        cst = {}

        tree_assign = {}

        '''
        Enumerate costs for each SNV cluster and each tree/group
        '''
        for q, dcf in self.delta.items():
            cst[q] = {}
            for g, trees in enumerate(self.T_SNV_groups):
                cst[q][g] = {}
                tree_assign[g] = {}
              
               
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
        
        df = pd.DataFrame(all_costs,
                           columns=["snv", "snv_clust", "dcf", "group", "tree", "cost", "posterior_dcf"])
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

          
        rev_mapping = {val: key for key,val in mapping.items()}

        for n in ct_tree:

            cn_state = rev_mapping[n][1]
            cna_geno = CNAgenotype(*cn_state)
            genotypes[n] = {}
            anc = list(nx.ancestors(ct_tree, source=n))
            anc.append(n)
            present_clusters = [q for q in alpha_inv if q in anc]
            pres_snvs = set([j  for q in present_clusters for j in alpha_inv[q]])
            not_added_muts = []
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
                    not_added_muts.append(j)
            
            for j in set(self.snvs).difference(pres_snvs):
                genotypes[n][j] = genotype(*cn_state, 0,0)
        ct = ClonalTree(ct_tree, genotypes, seg_to_snvs)
        if len(ct.mut_mapping[ct.root]) >0:
            print(ct.mut_mapping[ct.root])
            print(not_added_muts)

        return ct
                    





    @staticmethod
    def solve(model,vars, threads=1, timelimit=120):
        model.Params.Threads = threads
        model.Params.TimeLimit = timelimit
        model.setParam('LogFile', '')
     
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
  

        phi ={}

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

    #TODO: fix
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
            # if j ==384:
            #     print("here")
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
                    best_group = g
            if psi[j] != new_node:
                # snv_tree.draw("test/snv_tree.png", segments=[self.ell])
                
                T.update_genotype(new_node, j,tree_assign[best_group][j][0])
        
        T.update_mappings()
        if   len(T.mut_mapping[T.root]) > 0:
            T.draw("test/bad_ct.png", ca, segments = [self.ell])

        return T 
                
    def get_cn_dcfs(self):
        obs_copy_x, obs_copy_y =  self.data.copy_profiles_by_seg([self.ell], self.data.cells)
        def node_cna_cost(cn):
    

            latent_x = np.array([cn[0]]).reshape(1,-1)
            latent_y = np.array([cn[1]]).reshape(1,-1)
      
            x_diff = np.abs(obs_copy_x - latent_x).sum(axis=1)
            y_diff = np.abs(obs_copy_y - latent_y).sum(axis=1)
            return x_diff + y_diff
 

        cell_cna_scores = np.vstack([node_cna_cost(cn) for cn, _ in self.cn_states.items()])
        nodes = [v for _,v in self.cn_states.items()]
        cell_assign = np.argmin(cell_cna_scores, axis=0)
        cell_totals = {v: 0 for v in nodes}
        for i, v in zip(self.data.cells, cell_assign):
            cell_totals[nodes[v]] += 1
        
        cell_fractions = {v: cell_totals[v]/self.data.N for v in cell_totals}
        cn_dcfs = {v: cell_fractions[v] for v in nodes}
      
        for u in self.S:
            for v in nx.descendants(self.S,u):
                cn_dcfs[self.cn_states[u]] += cn_dcfs[self.cn_states[v]]

        return cn_dcfs

    @staticmethod
    def check_dcfs(T, merged_delta):
        for u in T:
            children_dcf = 0
            for v in T.successors(u):
                children_dcf += merged_delta[v]
            if merged_delta[u] < children_dcf or children_dcf > 1:
                return False
        return True
        
     
        
    def fit(self, data, segment):
        self.data = data
        self.ell = segment
        self.snvs = self.data.seg_to_snvs[segment]
        self.cn_dcfs = self.get_cn_dcfs()

        refinements = Enumerate( self.T_m, self.S,).solve()
        
        #TODO: fix to account for errors in observed cn states
        self.cn_props = self.data.cn_proportions(self.ell)

        cost, df, tree_assign = self.compute_snv_cluster_tree_cost()

        results = []
        opt_cost = np.Inf 
        print(f"Total refinements: {len(refinements)}")
        merged_dcfs = self.cn_dcfs | self.delta
        for f, T in enumerate(refinements):

     
            cluster_groups = self.identify_snv_cluster_trees(T)
            alpha_inv, omega = self.cluster_snvs(df, cluster_groups, tree_assign)
            snv_clusters = list(alpha_inv.keys())
            
            segment_tree= self.construct_segment_tree(T, alpha_inv, omega)
            
            if not self.check_dcfs(segment_tree.tree, merged_dcfs):
                continue

            print(f"Starting refinement {f}...")
            best_phi = None
            cost  = np.Inf
 
            for _ in range(self.max_iterations):
                ca_cost, ca = self.assign_cell_clusters(segment_tree, snv_clusters )
                if ca_cost is None:
                    break
                segment_tree = self.assign_genotypes(segment_tree, ca, tree_assign, snv_clusters)
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
                results.append(Solution(cost, best_segtree, best_phi))

        return sorted(results, key= lambda x: x.cost)



 