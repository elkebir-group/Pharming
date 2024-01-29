
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



def draw(tree, fname):
    ptree = pgv.AGraph(strict=False, directed=False)
    ptree.add_edges_from(list(tree.edges))
    ptree.layout("dot")
    ptree.draw(fname)

@dataclass
class STISol:
    cost: float
    ct: ClonalTree
    phi: dict
    psi: dict 
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
    def __init__(self, S, k=4, r=3,seed=1026, max_iter=1, nrestarts=1, tolerance=1.5, lamb1 =5, lamb2=5) -> None:
        self.k = k 
        self.r = r

        nodes = list(S)
        if len(nodes) == 0:
                 raise ValueError("CNA Tree S is not valid!")
        else:
            self.S = S 
  
        self.max_iter = max_iter
        self.nrestarts = nrestarts
        self.tolerance = tolerance
        self.rng = np.random.default_rng(seed)
        self.lamb1 = lamb1 
        self.lamb2  =lamb2

        if len(list(self.S)) > 0:
            snv_tree_edgelists = clonelib.get_genotype_trees(list(self.S))
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
        
     

    
    @staticmethod
    def to_inverse_dict(mydict):
        '''
        reverses an assignment/mapping dictionary
        '''
        rev_dict = {}
        for key, val in mydict.items():
            for v in val:
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
            rev_dict[val] = key

        return rev_dict
    
 
    def initialize_dcfs(self):
        return {q: self.rng.random() for  q in range(self.k)}
    

    def initialize_cell_clusters(self, cells):
        C_x, C_y = self.data.copy_profiles_by_seg(cells, self.ell)
        
        x = np.array([u[0] for u in self.S]).reshape(1,-1)
        y = np.array([u[1] for u in self.S]).reshape(1,-1)

        abs_diff = np.abs(C_x - x) + np.abs(C_y - y)
        init_cell_assign = np.argmin(abs_diff, axis=1)
        beta_inv = {}
        for q, state in enumerate(self.S):
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
         



    def fit(self, data, segment):
        self.data = data
        self.ell = segment
        self.snvs = self.data.seg_to_snvs(segment)

        sol_list = []
        for w in range(self.restarts):
            
            #a dictionary of len k with the current DCF of each SNV cluster q
            delta = self.initialize_dcfs()
            
            #a cell to cluster mapping alpha: [n] -> [r]
            beta, beta_inv  = self.initialize_cell_clusters(self.data.cells)      
            
            #a mapping of cell cluster to node in clonal tree for the segment
            phi = {}
            prev_cost = np.Inf
            for i in range(self.max_iter):
           
                alpha_inv, omega = self.update_SNVclusters(delta, beta_inv)
                delta = self.upate_dcfs(alpha_inv, beta_inv)
                ct, phi, psi = self.update_segment_tree(alpha_inv,  beta_inv, delta)
                beta_inv, current_cost = self.update_cell_clusters(ct, beta_inv, psi)

                if np.abs(current_cost - prev_cost) <= self.tolerance or i == self.max_iter-1:
                    sol_list.append(STISol(current_cost, ct.deepcopy(), phi, psi, 
                                           self.to_mapping_dict(alpha_inv), 
                                           self.to_mapping_dict(beta_inv), delta))
                    break 
                else:
                    prev_cost = current_cost


   
    
    
    def convert_to_clonal_tree(self, t, j):
        geno_dict = {}
        relabel = {}
        for i, v in enumerate(t):
            geno_dict[i] = {j : genotype(*v) }
            relabel[v] = i

        t_seg_to_muts = {self.seg: [j]}
        t_copy = nx.relabel_nodes(t, relabel)
        ct = ClonalTree(t_copy, geno_dict, t_seg_to_muts)
        return ct, {n: v for v,n in relabel.items()}
    

    @staticmethod
    def solve(model,vars =[x], threads=1, timelimit=60):
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
        
        
        


    


    def compute_snv_tree_cost(self, dcf, j, tree, beta_inv):
        clust_cost = {}
          
        ct, mapping = self.convert_to_clonal_tree(tree,j)

        #TODO: identify split node id in clonal tree and identify descednant
        split_node =  1
        descendants = 2

        _, _, cell_scores, nodes = ct.assign_cells(self.data, self.lamb1)
        nodes = nodes.tolist()
        for q, cell_clust in beta_inv.items():
            clust = cell_scores[:, cell_clust].sum(axis=1)
            for i, u in enumerate(nodes):
                clust_cost[q, u] = clust[i]
        
        model = gp.Model("MIP")
        clone_assign = [(i,u) for q in beta_inv for u in nodes]
        x = self.model.addVars(self.clone_assign, vtype = GRB.BINARY )
        x_pos = self.model.addVar( vtype = GRB.BINARY )
        x_neg =self.model.addVar(vtype = GRB.BINARY )

        model.setObjective(gp.quicksum(clust_cost[q,u]*self.x[q,u] for q,u in clone_assign) + 
                           self.lamb2*(x_pos) + self.lamb2*x_neg,
                             gp.GRB.MINIMIZE)
        
        #every cluster is assigned to exactly 1 node
        model.addConstrs(gp.quicksum(x[q,u] for q in beta_inv)==1 for u in nodes)

        #TODO: set up constraints for absolute values

        objval, sol = self.solve(model, x)

        tree_phi ={}
        for q,u in clone_assign:
    
                if sol[q,u] > 0.5:
                   tree_phi[q] = u 
        return objval, tree_phi



        

    
    def update_SNVclusters(self, delta, beta_inv ):
        
        cst = {}
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
                    for t in trees:
                        
                        cost, ct, tree_phi = self.compute_snv_tree_cost(j, dcf, t, beta_inv)
                        if cost < g_cost:
                            cst[(j,q,g)] = cost 
                            tree_assign[j][g] = (t, ct,tree_phi)
        
        model = gp.Model("MIP")
        clust_group_assign = [(j,q,g) for j in self.snvs for q in range(self.k) for g in range(num_groups)]
        clust_assign = [(j,g) for j in self.snvs for q in range(self.k)]
        group_assign =  [(q,g) for q in range(self.k) for g in range(num_groups)]
        snv_clust = [(j,q) for j in self.snvs for q in range(self.k)]
        x = self.model.addVars(clust_group_assign, vtype = GRB.BINARY )
        y =self.model.addVars(group_assign, vtype=GRB.BINARY)
        z = self.model.addVars(clust_assign, vtype = GRB.BINARY )
        w = self.model.addVars(snv_clust, vtype= GRB.BINARY)
       
        #minimize cost
        model.setObjective(gp.quicksum(cst[j,q,g]*x[j,q,g] for j,q,g in clust_group_assign), gp.GRB.MINIMIZE)
        
        #every snv cluster q has at least 1 SNV assigned (clusters are non-empty)
        model.addConstrs(gp.quicksum(z[j,q] for j in self.snvs) > 0 for q in range(self.k))

        #every snv is assinged to exactly 1 cluster
        model.addConstrs(gp.quicksum(z[j,q] for q in range(self.k)) ==1 for j in range(self.snvs))

        #every snv cluster is assigned to one group
        model.addConstrs(gp.quicksum(y[q,g] for g in range(num_groups)) ==1 for q in range(self.k))

        #TODO: add pairwise compatible of selected groups 

        for j,q,g in clust_group_assign:
            
            #linking variables
            model.addConstr(x[j,q,g] >= y[j,g])
            model.addConstr(x[j,q,g] >= z[q,g])
            model.addConstr(x[j,q,g] >= w[j,q])
            #every snv assigned to snv cluster q must all have the same group assignment z[j,g
            model.addConstr(w[j, q] == z[q, g])


        objval, sols = self.solve(model, [x,y,z,w])

        

        


            

    def update_dcfs(self,beta_inv, omega):
        pass 

    def update_tree(alpha_inv, beta_inv, delta):
        pass 

    def update_cell_clusters(ct, delta, psi):
        pass 


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