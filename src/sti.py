
import networkx as nx
import numpy as np
import pandas as pd 
from itertools import chain, combinations, permutations
from clonal_tree import ClonalTree
import seaborn as sns
import matplotlib.pyplot as plt
import pygraphviz as pgv
from genotype import genotype, CNAgenotype
from copy import deepcopy
import clonelib
from cell_mapping import CellAssign
from constrained_cell_assignment import ConstrainedCellAssign
from dataclasses import dataclass
import gurobipy as gp
from gurobipy import GRB
import pickle 
from scipy.optimize import minimize_scalar
import itertools


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
    phi: dict
    alpha: dict 
    beta: dict 
    delta: dict 


class STI:
    '''
    a class to solve the segment tree inference problem
    
    S: a nx:Digraph represnting the CNA tree
    k: int representing the number of SNV clusters
    r: int representing the number of cell clusters
    seed: int representing the random number seed 

    '''
    def __init__(self, S, k=4, r=3,seed=11, max_iter=3, nrestarts=20, tolerance=1.5, lamb1 =5, lamb2=100) -> None:
        self.k = k 
        self.r = r

        nodes = list(S)
        if len(nodes) == 0:
                 raise ValueError("CNA Tree S is not valid!")
        else:
            self.S = S 
        

        self.cn_states = {}
        for i, u in enumerate(self.S):
            self.cn_states[u] = self.k +i + 1 
        
        self.cn_states_inv = {val: key for key,val in self.cn_states.items()}
  
        self.max_iter = max_iter
        self.nrestarts = nrestarts
        self.tolerance = tolerance
        self.rng = np.random.default_rng(seed)
        self.lamb1 = 10
        self.lamb2  =lamb2
        self.lamb2 =100

        if len(list(self.S)) > 0:
            snv_tree_edgelists = clonelib.get_genotype_trees(list(self.S.edges))
        else:
            snv_tree_edgelists = []
         
            root = nodes[0]
            for i in range(2):
                if i == 0:
                    snv_tree_edgelists.append( [((root[0], root[1], 0, 0), (root[0], root[1], 1, 0))])
                else:
                    snv_tree_edgelists.append( [((root[0], root[1], 0, 0), (root[0], root[1], 0, 1))])

  
        self.T_SNVs = [nx.DiGraph(edges) for edges in snv_tree_edgelists]
        self.T_SNV_groups, self.group_desc = self.group_snv_trees()
        
        for g, trees in enumerate(self.T_SNV_groups):
            for i,t in enumerate(trees):
                draw(t,f"test/tree_g{g}_t{i}.png")

    
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
    
 
    def initialize_dcfs(self):
        '''
        randomly initialize DCFs for each of the k SNV clusters
        '''
        return {q: self.rng.random() for  q in range(self.k)}
    

    def initialize_cell_clusters(self, cells):
        C_x, C_y = self.data.copy_profiles_by_seg( [self.ell], cells)
        
        x = np.array([u[0] for u in self.S]).reshape(1,-1)
        y = np.array([u[1] for u in self.S]).reshape(1,-1)

        abs_diff = np.abs(C_x - x) + np.abs(C_y - y)
        init_cell_assign = np.argmin(abs_diff, axis=1)
        beta_inv = {}
        for q, _ in enumerate(self.S):
           cell_indices =  np.where(init_cell_assign==q)[0]
           beta_inv[q] = [cells[i] for i in cell_indices]
        
        '''
        randomly split cell clusters until we have r cell clusters
        '''
        while len(beta_inv) < self.r:
            clusters = list(beta_inv.keys())
            q = self.rng.choice(clusters)
            cells_to_split = np.array(beta_inv[q])
            random_indices = self.rng.choice([0, 1], size=len(cells_to_split))

            # Use boolean indexing to split the array based on random_indices
            beta_inv[q] = cells_to_split[random_indices == 0].tolist()
            beta_inv[max(clusters) + 1] = cells_to_split[random_indices == 1].tolist()
        

        return  self.to_mapping_dict(beta_inv), beta_inv
         



    def fit(self, data, segment, delta=None, beta_inv=None):
        self.data = data
        self.ell = segment
        self.snvs = self.data.seg_to_snvs[segment]

        sol_list = []
        for w in range(self.nrestarts):
            
            #a dictionary of len k with the current DCF of each SNV cluster q
            if delta is  None:
                delta = self.initialize_dcfs()
            print("initializing DCFs....")
            for q,dcf in delta.items():
                print(f"{q}: {dcf}")
            
            #a cell to cluster mapping alpha: [n] -> [r]
            if beta_inv is None:
                beta, beta_inv  = self.initialize_cell_clusters(self.data.cells)      
            
            #a mapping of cell cluster to node in clonal tree for the segment
            phi = {}
            prev_cost = np.Inf
            for i in range(self.max_iter):
           
                alpha_inv, omega, snv_clust_to_group = self.update_SNVclusters(delta, beta_inv)
                alpha = self.to_mapping_dict(alpha_inv)
                for q, val in alpha_inv.items():
                    print(f"k: {q} dcf: {delta[q]} nsnvs: {len(val)}")
                return alpha

                delta = self.update_dcfs(alpha_inv, beta_inv, omega)
                G, node_cn_states = self.construct_clonal_graph(delta, snv_clust_to_group)
                clonal_trees = self.enumerate_trees(G, alpha_inv, delta, omega, node_cn_states)
            
                ct, phi = self.update_tree(clonal_trees, beta_inv, delta)

                current_cost, beta_inv, phi = self.update_cell_clusters(ct, delta)
                print(f"iteration: {i} diff: {np.abs(current_cost-prev_cost)} current cost: {current_cost} previous cost: {prev_cost}")
                sol_list.append(STISol(current_cost, deepcopy(ct), phi, 
                                           self.to_mapping_dict(alpha_inv), 
                                           self.to_mapping_dict(beta_inv), delta))
                if np.abs(current_cost - prev_cost) <= self.tolerance:
                    break 
                else:
                    prev_cost = current_cost
        
        return sorted(sol_list, key= lambda x: x.cost)
    
   
    def group_to_graph(self, group, cluster_ids, delta): 
        sscn =self.group_desc[group]['sscn']
        children = self.group_desc[group]['children']
    
      
        centers = np.array([delta[q] for q in cluster_ids])
        G= nx.DiGraph()


        for u in self.cn_states:
            # G.add_node(self.cn_states[u])
            if u == sscn:
                root = self.cn_states[u]
                G.add_node(self.cn_states[u])

        if len(children) ==0:
            # for u,v in self.S.edges:
            #     G.add_edge(cn_states[u], cn_states[v])

            for q in cluster_ids:
                G.add_edge(root, q)

            for u in G:
                if u in delta:
                    for v in G:
                        if v in delta:
                            if delta[u] > delta[v]:
                                G.add_edge(u,v)
                        
        else:
                #sort the centers 
                sorted_clusters = np.argsort(centers)[::-1]
                parent = root
                for i,k in enumerate(sorted_clusters):
                    G.add_edge(parent, cluster_ids[k])
        
        
                    parent = cluster_ids[k] 
                
                # for u,v in self.S.edges:
                #     if u == sscn and v in children:
                #             G.add_edge(parent, cn_states[v])
                #     else:
                #         G.add_edge(cn_states[u], cn_states[v])
        
                
        # draw(G, f"test/graph{group}.png")
        return G


    def construct_clonal_graph(self, delta, snv_clust_to_group):

        group_to_snv_clust = self.to_inverse_dict(snv_clust_to_group)

        # cn_states = {}
        # for i, u in enumerate(self.S):
        #     cn_states[u] = self.k +i + 1 
        refinements ={}
        groups_by_sscn = {u: [] for u in self.S}
        for g in group_to_snv_clust:
            refinements[g] = self.group_to_graph(g, group_to_snv_clust[g], delta)
            groups_by_sscn[self.group_desc[g]['sscn']].append(g)

        merged_graphs = []
        for sscn in groups_by_sscn:
            if len(groups_by_sscn[sscn]) > 0:
                G= self.merge_groups(sscn,  refinements, groups_by_sscn[sscn], group_to_snv_clust,delta )
                merged_graphs.append(G)
        # clonal_graph= self.integrate_refinements(refinements,group_to_snv_clust)
        G= nx.compose_all(merged_graphs)

        node_cn_states = {val: key for key, val in self.cn_states.items()}
        for u in G:
            if u not in node_cn_states:
                grp = snv_clust_to_group[u]
                node_cn_states[u] = self.group_desc[grp]['sscn']
        # draw(G, "test/clonal_graph.png")

        for u in self.cn_states_inv:
            if u not in  G:
                parent = list(self.S.predecessors(self.cn_states_inv[u]))
                if len(parent) == 1:
                    G.add_edge(self.cn_states[parent[0]], u)

        if not nx.is_weakly_connected(G):
            components = list(nx.weakly_connected_components(G))
            components =sorted(components, key=lambda x: len(x), reverse=True)
            main_comp = components.pop(0)
            main_comp_cna = [u for u in main_comp if u in self.cn_states_inv]
     
            for comp in components:
                roots = [u for u in comp if u in self.cn_states_inv ]

                for u in main_comp_cna:
                    for v in roots:
                        if (node_cn_states[u], node_cn_states[v]) in list(self.S.edges):
                            G.add_edge(u,v)
                            break
            
        draw(G, "test/clonal_graph.png")


        #finally merge the merged graphs into a clonal graph
        return G, node_cn_states


    def merge_groups(self, sscn, refinements, groups, group_to_snv_clust, delta):
     
        #sort in reverse order from the most number of children to least
        sorted_groups = sorted([(g, self.group_desc[g]["children"]) for g in groups], key=lambda x: len(x[1]), reverse=True)
        
        leaf_group= None

        # sorted_groups = sorted_groups[::-1]

        g, children = sorted_groups.pop(0)
  

        G  = refinements[g]
        def find_last(graph, g):
            for q in group_to_snv_clust[g]:
                if graph.out_degree[q] ==0:
                        return q
             
        children_last = []
        new_last = (find_last(G,g))
        children_last.append((new_last, children))


        if len(children) == 1:
            child = children.pop()
            G.add_edge(new_last, self.cn_states[child])

        

       

        while len(sorted_groups) >0:
            g, children = sorted_groups.pop(0)
            G2 = refinements[g]
            # draw(G, "test/G.png")
            # draw(G2, "test/G2.png")
            added = False 
            for end, kids in children_last:
                if len(children.intersection(kids)) > 0:
                    new_root = list(G2.successors(self.cn_states[sscn]))[0]
                    new_last = find_last(G2, g)
                    G2.remove_node(self.cn_states[sscn])
                    G= nx.union(G,G2)
                    G.add_edge(end,new_root)
                    children_last.append((new_last, children))
                    added =True
                    break

                 
            if not added:
                #then its the leafset or on a distinct branch
                G = nx.compose(G, G2)
                new_last = find_last(G2, g)
            
             
            if len(children) == 1:
                child = list(children.copy())[0]
                G.add_edge(new_last, self.cn_states[child])
            
            elif len(children)==0:
                leaf_group = g 
                leafs_added = group_to_snv_clust[leaf_group]

            # draw(G, "test/merge_graph.png")
        for g in groups:
            if leaf_group is not None:
                if g != leaf_group:
                    for q in group_to_snv_clust[g]:
                        for u in leafs_added:
                            if delta[q] > delta[u]:
                                G.add_edge(q,u)
        return G




    
        

  


    





        

    
    def convert_to_clonal_tree(self, t, j):
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
        

    def update_cell_clusters(self, ct, delta):
        cell_cost = {}
        cells = self.data.cells 
        N = len(cells)
        _, _, cell_scores, nodes = ct.assign_cells(self.data, self.lamb1)
        nodes = nodes.tolist()
        for j,u in enumerate(nodes):
            for i, c in enumerate(cells):
                cell_cost[c, u] = cell_scores[j,i]
        
        model = gp.Model("MIP")
        cell_assign = [(i,u) for i in cells for u in nodes]
        x = model.addVars(cell_assign, vtype = GRB.BINARY )
        z = model.addVars(delta.keys(), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS )
        
        #indicates if node u has any cells assigned 
        y = model.addVars(nodes, vtype=GRB.BINARY)



        model.setObjective(gp.quicksum(cell_cost[i,u]*x[i,u] for i,u in cell_assign) + 
                           self.lamb2*gp.quicksum(z[u] for u in delta), gp.GRB.MINIMIZE)
        
        #a node can only have assigned cells if it is used 
        for u in nodes:
            model.addConstrs(x[i,u] <= y[u] for i in cells)
        
        model.addConstrs((gp.quicksum(x[i,u] for u in nodes)==1) for i in cells )
        
        #if a node is used it must have assigned cells (non-empty cell cluster)
        for u in nodes:
            model.addConstr(gp.quicksum(x[i,u] for i in cells) >= y[u])

        #exactly r nodes must be used 
        model.addConstr(gp.quicksum(y[u] for u in nodes)==self.r)
        for u in delta:
            descendants = list(ct.preorder(u))
            model.addConstr( (1/N)*gp.quicksum(x[i,u] for i in cells for u in descendants ) - delta[u] <= z[u] )
            model.addConstr( (1/N)*gp.quicksum(x[q,u] for q in cells  for u in descendants ) - delta[u] >= -1*z[u] )
        objval, sol = self.solve(model, [x,y])
        x_sol, y_sol = sol[0], sol[1]
        cluster_id = 0

        phi_inv = {}
        for u in nodes:
            if y_sol[u] > 0.5:
                phi_inv[u] = cluster_id 
                cluster_id += 1
        beta_inv = {q: [] for q in range(self.r)}
     
        for i,u in cell_assign:
            if x_sol[i,u] > 0.5:
                beta_inv[phi_inv[u]].append(i)
        
        phi = {val: key for key,val in phi_inv.items()}
        return objval,  beta_inv, phi
            









        # for q, cell_clust in beta_inv.items():
        #     clust = cell_scores[:, cell_clust].sum(axis=1)
        #     for i, u in enumerate(nodes):
        #         clust_cost[q, u] = clust[i]
        


    def assign_cell_clusters(self,ct, beta_inv, delta,allow_multiple=True):
        # if 78 in ct.get_all_muts():
        #     ct.draw("test/snv_tree78.png")
        #     print("here")
        clust_cost = {}
        _, _, cell_scores, nodes = ct.assign_cells(self.data, self.lamb1, lamb_vaf=5)
        # cell_scores, nodes = ct.compute_node_likelihoods(self.data, lamb=100)
        nodes = nodes.tolist()
        clust_score = []
        for q, cell_clust in beta_inv.items():
  
            # print(f"{q}: {self.data.compute_vafs(cell_clust, list(ct.psi.keys()))}")
            # obs_vafs = self.data.obs_vafs(cell_clust, list(ct.psi.keys()))
            # print(obs_vafs)
            clust = cell_scores[:, cell_clust].sum(axis=1) /len(cell_clust)

            clust_score.append(clust)
            for i, u in enumerate(nodes):
                clust_cost[q, u] = clust[i]
        
        cs = np.vstack(clust_score)
        best_clust = cs.argmin(axis=1)
        objval = cs.min(axis=1).sum()
        phi ={}
        # for p, i in zip(beta_inv, best_clust):
        #     phi[p] = nodes[i]
        model = gp.Model("MIP")
        clone_assign = [(q,u) for q in beta_inv for u in nodes]
        x = model.addVars(clone_assign, vtype = GRB.BINARY )
        z = model.addVars(delta.keys(), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS )


        model.setObjective(gp.quicksum(clust_cost[q,u]*x[q,u] for q,u in clone_assign) + 
                           self.lamb2*gp.quicksum(z[u] for u in delta), gp.GRB.MINIMIZE)
        
        
        #every cluster is assigned to exactly 1 node
        model.addConstrs(gp.quicksum(x[q,u] for u in nodes)==1 for q in  beta_inv)



        #TODO: set up constraints for absolute values
        for u in delta:
            descendants = list(ct.preorder(u))
            model.addConstr( (1/self.data.N)*gp.quicksum(len(beta_inv[q])*x[q,u] for q in beta_inv for u in descendants ) - delta[u] <= z[u] )
            model.addConstr( (1/self.data.N)*gp.quicksum(len(beta_inv[q])*x[q,u] for q in beta_inv for u in descendants ) - delta[u] >= -1*z[u] )
        
        
        if not allow_multiple:
            #allow at most 1 cluster 
            for u in nodes:
                model.addConstr(gp.quicksum(x[q,u] for q in beta_inv) <=1)
        
        objval, sol = self.solve(model, [x, z])
        x_sol = sol[0]
        z_sol = sol[1]
  
    
        for q,u in clone_assign:
    
                if x_sol[q,u] > 0.5:
                   phi[q] = u 
        
        cell_phi = {}
        for q in phi:
            for i in beta_inv[q]:
                cell_phi[i] = phi[q]
        ca = CellAssign(cell_phi, ct.clones())
        dcf = self.compute_dcfs(ct.get_all_muts()[0], ct, ca)
    
        return objval, phi, dcf


    def compute_snv_tree_cost(self, j, dcf, tree, beta_inv):
 
          
        ct, mapping = self.convert_to_clonal_tree(tree,j)

        # ct.draw("test/snv_tree.png")


        _, u, _, _ =  ct.get_split_nodes(j)
        obj, phi ,dcf = self.assign_cell_clusters(ct, beta_inv, delta={u: dcf})
        
        return obj, ct, phi, dcf

 

    def compute_dcfs(self, j, clonal_tree, cell_assign):



        # clonal_tree.draw("test/ct.png")
        # for u in clonal_tree.preorder():
        #     print(f"{u}: {clonal_tree.genotypes[u][j]}")
        # print(cell_assign.get_cell_count())

        states= list(self.S.nodes)
        parent, child, p_geno, c_geno = clonal_tree.get_split_nodes(j)
        desc_genotypes = clonal_tree.get_desc_genotypes(child, j)
        cn_star = p_geno.to_CNAgenotype().to_tuple()
        cna_genos= clonal_tree.get_cna_genos()[self.ell]
        ncells = len(cell_assign.get_all_cells())

        def compute_lambda(v, F, cn_prop):
            mu_x_y = cn_prop[cn_star]
            term1 = 1/(c_geno.z*mu_x_y)
            term2 = v*F
            term3 = 0
            for g in desc_genotypes:
                g_cna = g.to_CNAgenotype()
                if  g_cna != p_geno.to_CNAgenotype():
                    term3 += g.z*cn_prop[g_cna.to_tuple()]
            return term1*(term2 - term3)


        
        def v_to_dcf(v, F, cn_prop):
        
            # cn = genotype(*state, 0,0)
            #if split state
            # if state == cn_star:
            d = (1/c_geno.z)*(v*F -sum([(g.z-c_geno.z)*cn_prop[g.to_CNAgenotype().to_tuple()] for g in desc_genotypes]))
            # else:
            #     d = sum([cn_prop[g.to_CNAgenotype().to_tuple()] for g in desc_genotypes if g.cna_eq(cn)])
     
            if d > 1 and d <= 1.5:
                 return 1
            elif d > 1.5 or d < 0:
                return 1e8
            else:
                
                return d
            # if d < 0: return 0
            # elif d > 1: return 1
            # return d
        
        state_to_nodes = {state: [] for state in states}
        for n, state in cna_genos.items():
            state = state.to_tuple()
            state_to_nodes[state].append(n)
   
       
        
        cells_by_state = {state: [] for state in state_to_nodes}
    
        dcfs = {}
    
        for state, nodes in state_to_nodes.items():
            for n in nodes:
                cells_by_state[state] += cell_assign.get_cells(n)
        
        cn_prop = {state: len(cells_by_state[state])/ncells for state in cells_by_state}
        F = sum([(state[0]+ state[1])* cn_prop[state] for state in cn_prop])
        vaf =self.data.compute_vafs(None, snvs=[j])
        lamb = compute_lambda(vaf[0], F, cn_prop)
        if lamb >=0 and lamb <=1:
            dcf= (1-lamb)*cn_prop[cn_star] 
            desc_cn_states = []
            for g in desc_genotypes:
                if g.to_CNAgenotype() != p_geno.to_CNAgenotype():
                    desc_cn_states.append(g.to_CNAgenotype().to_tuple())
            for cn in set(desc_cn_states):
                dcf += cn_prop[cn]
            return dcf 
        else:
            return 1e8
        
                
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
        
        return self.to_inverse_dict(alpha), omega, clust_to_group




        

        
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
            

    def update_dcfs(self,alpha_inv, beta_inv, omega):
        new_dcfs = {}    
        for q, snvs in alpha_inv.items(): 
            new_dcfs[q]= minimize_scalar(self.scalar_obj,
                                 args=(snvs, omega, beta_inv, self.data.N),
                                   method='bounded', bounds=[0,1]).x

        return new_dcfs
    

    def enumerate_trees(self, G, alpha_inv, delta, omega, node_to_cn_state):
        # draw(G, "test/clonal_graph.png")
        root = [u  for u in G if G.in_degree[u] ==0][0]
        for u,v in G.edges:
            G[u][v]["weight"] = 1

        def check_sum_condition(tree):
            return True 
            for u, dcf in delta.items():
                kids = list(tree.successors(u))
                if len(kids) >0:
                    children_dcfs = sum([delta[v] for v in kids if v in delta])
                    if dcf < children_dcfs:
                        return False
            return True 
           
        def generate_genotypes(tree, root):
            genotypes ={}
            for u in nx.dfs_preorder_nodes(tree):
                genotypes[u] ={}
                x,y = node_to_cn_state[u]
                nodeCN = CNAgenotype(x,y)
                path = nx.shortest_path(tree, root, u)
                all_snvs = [alpha_inv[q] for q in path if q in alpha_inv]
                present_snvs = set(itertools.chain.from_iterable(all_snvs))

                not_present_snvs = set(self.snvs) - present_snvs
                for j in not_present_snvs:
                    genotypes[u][j] = genotype(x,y,0,0)
                for j in present_snvs:
                    ct, _= omega[j]
                    cna_genos = ct.get_cna_genos()[self.ell]
                    for v in ct.preorder():
                        geno = ct.genotypes[v][j]
                        if cna_genos[v] == nodeCN and geno.z > 0:
                          
                           genotypes[u][j] = genotype(*geno.to_tuple())
            return genotypes
                  

        clonal_trees = []
        trees = nx.algorithms.tree.branchings.ArborescenceIterator(G)
        for tree in trees:
    
            if check_sum_condition(tree):
        
                genotypes = generate_genotypes(tree,root)
                clonal_trees.append(ClonalTree(tree,genotypes,{self.ell : self.snvs}))
        
        return clonal_trees
                


    def update_tree(self, clonal_trees, beta_inv, delta):
        
     
        phi_list = []
        obj_vals = []
        for ct in clonal_trees:
            obj, phi = self.assign_cell_clusters(ct, beta_inv, delta, allow_multiple=False)
            obj_vals.append(obj)
            phi_list.append(phi) 
        
        index = obj_vals.index(min(obj_vals))
        return  clonal_trees[index], phi_list[index]
        

            
            




       




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