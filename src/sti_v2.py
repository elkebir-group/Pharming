
import networkx as nx
import numpy as np
import pandas as pd 
from itertools import chain, combinations
from clonal_tree import ClonalTree
import clonelib
from cell_mapping import CellAssign
from enumerate import Enumerate
from solution import Solution
from utils import load_pickled_object, pickle_object, timeit_decorator,draw
from genotype_tree import GenotypeTree








class STI:
    '''
    a class to solve the segment tree inference problem
    
    S: a nx:Digraph represnting the CNA tree
    k: int representing the number of SNV clusters
    r: int representing the number of cell clusters
    seed: int representing the random number seed 

    '''
    def __init__(self, ell, S, delta, lamb =5, niter=10, ilp=False) -> None:
    
        self.ell = ell 
        nodes = list(S)
        if len(nodes) == 0:
                 raise ValueError("CNA Tree S is not valid!")
        else:
            self.S = S 
        self.delta = delta
    
  

        self.lamb = lamb
        self.ilp = ilp 

        self.max_iterations = niter
        self.S_root = [n for n in self.S if S.in_degree[n]==0][0]
        self.cn_states = {}
        self.k = len(self.delta)
        for i, u in enumerate(nx.dfs_preorder_nodes(self.S, self.S_root)):
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
        self.cost1, self.cost2 = None, None

        # print(self.group_desc)
        # for g, trees in enumerate(self.T_SNV_groups):
        
        #     for idx, tree in enumerate(trees):
        #         draw(tree, f"test/group{g}_{idx}.png")



    
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
            geno_dict[i] = {j : v}
            relabel[v] = i

        t_seg_to_muts = {self.ell: [j]}
        t_copy = nx.relabel_nodes(t, relabel)
        ct = ClonalTree(t_copy, geno_dict, t_seg_to_muts)
        return ct, {n: v for v,n in relabel.items()}
    

  

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
                    parent = parent[0]
                    child = u
                
                    # if parent.cna_eq(child):
                    if (parent[0], parent[1]) == (child[0], child[1]):
                        sscn = (child[0], child[1])  #child.to_CNAgenotype().to_tuple()
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

    
    
    def compute_tree_cost(self, j, tree):
 
          
        ct, mapping = self.convert_to_clonal_tree(tree,j)

        # _, u, _, _ =  ct.get_split_nodes(j)

        obj = ct.assign_cells_by_likelihood(self.data,self.data.cells, lamb=self.lamb, cellassign=False)
        
        return obj, ct
    
    # #TODO: FIX THIS!!!!!
    # @timeit_decorator
    def precompute_costs_old(self, data):
        self.snvs = data.seg_to_snvs[self.ell]
        self.cn_props = data.cn_proportions(self.ell)
        if self.S_root not in self.cn_props:
            self.cn_props[self.S_root] = 0.0
        self.data  = data

        self.tree_assign = {}

        self.cost2 = {}
        alt = data.var[:, self.snvs].sum(axis=0)
        total =data.total[:, self.snvs].sum(axis=0)
        jg_list = []
        self.all_costs = []
        for g, trees in enumerate(self.T_SNV_groups):
     
            self.tree_assign[g] = {}
            # jg_cost = np.full(len(self.snvs), fill_value=np.Inf)
 
            for j in self.snvs:
                jg_cost = np.Inf
                for p, t in enumerate(trees):
                    # gt = GenotypeTree(t.edges)
        
                    cost, ct = self.compute_tree_cost(j, t)
                    if cost < jg_cost:
                        self.tree_assign[g][j] = ct
                        jg_cost = cost
                        # if cost < jg_cost[idx]:
                        #     self.tree_assign[g][j] = (ct, ca)
                        #     jg_cost[idx] = cost
         
            for j, a,d in zip(self.snvs, alt, total):
                for q, dcf in self.delta.items():
                    ct = self.tree_assign[g][j]
                    post_dcf = -1*ct.posterior_dcf(j, dcf, a,d, self.cn_props)
                    self.all_costs.append({"snv": j, "snv_clust": q,  "group": g, "cost": cost, "posterior_dcf": post_dcf})




            


  

        # return all_costs, tree_assign
    @timeit_decorator
    def precompute_costs(self, data):
        self.snvs = data.seg_to_snvs[self.ell]
        self.cn_props = data.cn_proportions(self.ell)
        if self.S_root not in self.cn_props:
            self.cn_props[self.S_root] = 0.0
        self.data  = data

        self.tree_assign = {}

        self.cost2 = {}
        alt = data.var[:, self.snvs].sum(axis=0)
        total =data.total[:, self.snvs].sum(axis=0)
        jg_list = []

        for g, trees in enumerate(self.T_SNV_groups):
     
            self.tree_assign[g] = {}
            jg_cost = np.full(len(self.snvs), fill_value=np.Inf)

            for p, t in enumerate(trees):
                gt = GenotypeTree(t.edges)
            
                for idx, j in enumerate(self.snvs):

    
                    cost, ct = self.compute_tree_cost(j, t)
                    if cost < jg_cost[idx]:
                        self.tree_assign[g][j] = ct
                    
                        jg_cost[idx] = cost
         
            
                for q, dcf in self.delta.items():
                    post_dcf = -1*gt.vectorized_posterior(dcf, alt, total, self.cn_props)
                    for j,post in zip( self.snvs, post_dcf):
                        self.cost2[g,q,p,j] = post 


            jg_list.append(jg_cost.reshape(1,-1))
        self.cost1 = np.vstack(jg_list)




    def identify_snv_cluster_trees(self,T):
        '''
        Given a networkx DiGraph tree labeled by (q, (x,y)),
        return a set of tuples (q, g) where the SNV clusters are introduced
        '''
        rho = {self.ell: {}}
        root = (-1, (self.S_root))
        sscn = {}
        children = {}
        groups = []
        for v in nx.dfs_preorder_nodes(T, source=root):
            if v != root:
                v_q, v_cn  = v
                u_q, u_cn = list(T.predecessors(v))[0]  #parent of v
                if u_q != v_q and u_cn == v_cn:
                    sscn[v_q] = v_cn   #split copy number state of cluster q
                    children[v_q] = []
                    for d in nx.descendants(T, v):
                        _, d_cn = d 
                        if d_cn != v_cn:
                            _, d_par = list(T.predecessors(d))[0]
                            if d_par == v_cn:
                                children[v_q].append(d_cn)
                    for g in self.group_desc:
                        g_sscn = self.group_desc[g]['sscn']
                        g_children = self.group_desc[g]['children']
                        if sscn[v_q] == g_sscn and set(children[v_q]) == set(g_children):
                            groups.append((g,v_q))

        for g,q in groups:
            rho[self.ell][q] = self.T_SNV_groups[g]
        
        # if len(groups) ==0:
        #     print(self.ell)
            # draw(T, "test/no_valid_clusters.png")
            # pickle_object(self, "test/no_valid_clusters.pkl")
        return groups, rho 

     

        # for u in nx.dfs_preorder_nodes(T, source=root):
        #     u_q, u_cn = u
        #     for  v_q, v_cn in T.successors(u):
        #           if u_q != v_q and u_cn == v_cn:
        #               sscn[v_q] = u_cn
        #               children[v_q] = []
        #               kids = list(T.successors((v_q, v_cn)))
        #               for w in kids:
        #                   next_cn_state = self.get_next_cn_child(T, w)
        #                   if len
        #                   children[v_q]
               
                          

        #               for w_q, w_cn in nx.dfs_preorder_nodes(T, source = (v_q, v_cn)):
        #                   if w_cn != v_cn and w_cn not in children[v_q]:
        #                       children[v_q].append(w_cn)
   


    def cluster_snvs(self, valid_groups_snvclusts):
        groups = list( set(g for g, _ in valid_groups_snvclusts))
        group_to_clust = {g: [q for g1, q in valid_groups_snvclusts if g1==g ] for g in groups}
        group_costs = self.cost1[groups, :]
        group_assign = np.argmin(group_costs, axis=0)
        omega = {}
        psi = {}
        for j, idx in zip(self.snvs, group_assign):
            g= groups[idx]
            bestval = np.Inf 
            bestq  = None
            omega[j] = self.tree_assign[g][j]
            for q in group_to_clust[g]:

                for p,t in enumerate(self.T_SNV_groups[g]):
                    if self.cost2[g,q,p,j] <= bestval:
                        bestq = q 
                        bestval = self.cost2[g,q,p,j] 
            
            if bestq is None:
                print(bestval)

            psi[j] = bestq
 
        return self.to_inverse_dict(psi), omega



    def cluster_snvs_old(self, valid_groups_snvclusts):
        filtered_costs = [cost for cost in self.all_costs if (cost["group"], cost["snv_clust"]) in valid_groups_snvclusts]
        sorted_costs = sorted(filtered_costs, key=lambda x: (x["cost"], x["posterior_dcf"]))
        
        first_rows = {}
        for cost in sorted_costs:
            if cost["snv"] not in first_rows:
                first_rows[cost["snv"]] = cost
        
        psi = {row["snv"]: row["snv_clust"] for row in first_rows.values()}
        omega = {row["snv"]: self.tree_assign[row["group"]][row["snv"]] for row in first_rows.values()}

        return self.to_inverse_dict(psi), omega
    
 


    def construct_segment_tree(self, T, alpha_inv, omega, rho):
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
            cna_geno = cn_state
            genotypes[n] = {}
            anc = list(nx.ancestors(ct_tree, source=n))
            anc.append(n)
            present_clusters = [q for q in alpha_inv if q in anc]
            pres_snvs = set([j  for q in present_clusters for j in alpha_inv[q]])
            not_added_muts = []
            for j in pres_snvs:
                added = False
                ct = omega[j]
                cna_genos = ct.get_cna_genos()[self.ell]
                for v in ct.preorder():
                    geno = ct.genotypes[v][j]
                
                    if cna_genos[v] == cna_geno and ct.mut_copies(geno) > 0:
                        genotypes[n][j] =geno
                        added = True
                        break 
                if not added:
                    genotypes[n][j] =(*cn_state, 0, 0)
                    not_added_muts.append(j)
            
            for j in set(self.snvs).difference(pres_snvs):
                genotypes[n][j] = (*cn_state, 0,0)
        ct = ClonalTree(ct_tree, genotypes, seg_to_snvs, rho=rho)
        # if len(ct.mut_mapping[ct.root]) >0:
        #     print(ct.mut_mapping[ct.root])
        #     print(not_added_muts)

        return ct
                    


                
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
        
     
    # @timeit_decorator
    def fit(self, Tm_edges, data, segment):

        # if self.cost1 is None or self.cost2 is None:
        # self.precompute_costs(data)
        # self.precompute_costs_old(data)
   
        T_m = nx.DiGraph(Tm_edges)
        if len(self.delta) != len(T_m):
                 raise ValueError("Mutation cluster tree must \
                                  match the number of DCFs")
        self.data = data
        self.ell = segment
        self.snvs = self.data.seg_to_snvs[segment]
        self.cn_dcfs = self.get_cn_dcfs()
        # pickle_object(self, "test/sti.pkl")
        refinements = Enumerate( T_m, self.S,).solve()
        
        #TODO: fix to account for errors in observed cn states

   
        results = []
        if not self.ilp:
            dcfs = {}
        else:
            dcfs = self.delta

        merged_dcfs = self.cn_dcfs | self.delta

        for  T in refinements:

            valid_group_snvclusts, rho = self.identify_snv_cluster_trees(T)
          
            alpha_inv, omega = self.cluster_snvs(valid_group_snvclusts)
  
            segment_tree= self.construct_segment_tree(T, alpha_inv, omega, rho)
            if not self.check_dcfs(segment_tree.tree, merged_dcfs):
                continue

            cost, best_ca = segment_tree.optimize(self.data, self.lamb, 
                                                  max_iterations= self.max_iterations, dcfs=dcfs)

            results.append(Solution(cost, segment_tree, best_ca))
        
        return sorted(results, key= lambda x: x.cost)
            # if segment_tree.has_loss():
            #     print(i)
            #     print(list(self.S.edges))
            #     pickle_object(T, "test/T.pkl")
            #     pickle_object(self, "test/sti.pkl")
            #     draw(T, "test/T.png")

        
        # if len(results) ==0:
        #     print("pause")







 

     # def compute_snv_cluster_tree_cost(self):
    #     all_costs = []
 

    #     tree_assign = {}

    #     '''
    #     Enumerate costs for each SNV cluster and each tree/group
    #     '''

    #     for g, trees in enumerate(self.T_SNV_groups):
    #         tree_assign[g] = {}
    #         for j in self.snvs:
    #             jg_cost = np.Inf
    #             for p,t, in enumerate(trees):
    #                 cost, ct, ca = self.compute_tree_cost(j, t)
    #                 if cost < jg_cost:
    #                     tree_assign[g][j] = (ct,ca)
    #                     jg_cost = cost 
    #                 for q,dcf in self.delta.items():
    #                     alt = self.data.var[:, j].sum()
    #                     total = self.data.total[:,j].sum()
    #                     post_dcf = ct.posterior_dcf(j, dcf, alt, total, self.cn_props)
    #                     all_costs.append([j,q,dcf,g,p,cost, post_dcf])
                   
            
    #     df = pd.DataFrame(all_costs,
    #                        columns=["snv", "snv_clust", "dcf", "group", "tree", "cost", "posterior_dcf"])

    #     df['snv'] = df['snv'].astype(int)
    #     df['snv_clust'] = df['snv_clust'].astype(int)
    #     df['group'] = df['group'].astype(int)
    #     return df, tree_assign         

        # for q, dcf in self.delta.items():
        #     cst[q] = {}
        #     for g, trees in enumerate(self.T_SNV_groups):
        #         cst[q][g] = {}
        #         tree_assign[g] = {}
              
               
        #         for j in self.snvs:
        #             jg_cost = np.Inf
        #             for p,t in enumerate(trees):
                     
        #                 cost, ct, ca,post_dcf = self.compute_tree_cost(j, dcf, t)
                        
        #                 all_costs.append([j,q,dcf,g,p,cost, post_dcf])
        #                 if cost < jg_cost:
        #                     cst[q][g][j]  = cost 
        #                     jg_cost = cost
        #                     # delta_hat[j,q,g] = np.abs(obs_dcf - dcf)
        #                     tree_assign[g][j] = (ct,ca)

    # @timeit_decorator
    # def compute_snv_cluster_tree_cost(self):
    #     all_costs = []
    #     cst = {}

    #     tree_assign = {}

    #     '''
    #     Enumerate costs for each SNV cluster and each tree/group
    #     '''
    #     for q, dcf in self.delta.items():
    #         cst[q] = {}
    #         for g, trees in enumerate(self.T_SNV_groups):
    #             cst[q][g] = {}
    #             tree_assign[g] = {}
              
               
    #             for j in self.snvs:
    #                 jg_cost = np.Inf
    #                 for p,t in enumerate(trees):
                     
    #                     cost, ct, ca,post_dcf = self.compute_tree_cost(j, dcf, t)
                        
    #                     all_costs.append([j,q,dcf,g,p,cost, post_dcf])
    #                     if cost < jg_cost:
    #                         cst[q][g][j]  = cost 
    #                         jg_cost = cost
    #                         # delta_hat[j,q,g] = np.abs(obs_dcf - dcf)
    #                         tree_assign[g][j] = (ct,ca)
        
    #     df = pd.DataFrame(all_costs,
    #                        columns=["snv", "snv_clust", "dcf", "group", "tree", "cost", "posterior_dcf"])
    #     return cst, df, tree_assign


        # def assign_cell_clusters(self,ct, snv_clusters ):

    
        # delta =self.delta
        # if not self.ilp:
        #     return ct.assign_cells_by_likelihood(self.data, self.data.cells, self.lamb1)
        # # _, _, cell_scores, nodes = ct.assign_cells(self.data, self.lamb1)
        # cell_scores, nodes = ct.compute_node_likelihoods(self.data, self.data.cells, self.lamb1)
        # nodes = nodes.tolist()
        # valid_nodes = list(self.cn_states_inv.keys()) + snv_clusters
        # clust_score = {}
        # for i in range(cell_scores.shape[1]):
        #     for q,n in enumerate(nodes):
        #         if n in valid_nodes:
        #             clust_score[i,n] = cell_scores[q,i]
  

        # phi ={}

        # model = gp.Model("MIP")
        # cells = np.arange(cell_scores.shape[1])
        # cell_assign = [(i,u) for i in cells  for u in valid_nodes]
        # x = model.addVars(cell_assign, vtype = GRB.BINARY )
        # z = model.addVars(delta.keys(), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS )


        # # model.setObjective(gp.quicksum(clust_score[i,u]*x[i,u] for i,u in cell_assign) + 
        # #                    self.lamb2*gp.quicksum(z[u] for u in delta), gp.GRB.MINIMIZE)
        
        # model.setObjective(gp.quicksum(clust_score[i,u]*x[i,u] for i,u in cell_assign))
        # #every cell is assigned to exactly 1 node
        # model.addConstrs(gp.quicksum(x[i,u] for u in valid_nodes)==1 for i in cells)



        # #constraints for absolute values
        # for u in delta:
   
            
        #         descendants = set(ct.preorder(u))
        #         desc = descendants.intersection(set(valid_nodes))
        #         # model.addConstr( (1/self.data.N)*gp.quicksum(x[i,u] for i in cells for u in desc ) - delta[u] <= z[u] )
        #         # model.addConstr( (1/self.data.N)*gp.quicksum(x[i,u] for i in cells for u in desc) - delta[u] >= -1*z[u] )
        #         model.addConstr( (1/self.data.N)*gp.quicksum(x[i,u] for i in cells for u in desc ) >= delta[u])

        
        
        # objval, sol = self.solve(model, [x, z])
        # if objval == np.Inf:
        #     return None, None
        # x_sol = sol[0]
        # z_sol = sol[1]
  
    
        # for i,u in cell_assign:
    
        #         if x_sol[i,u] > 0.5:
        #            phi[i] = u 
    
        # ca = CellAssign(phi, ct.clones())
 
    
        # return objval, ca
    

    #   @staticmethod
    # def solve(model,vars, threads=1, timelimit=120):
    #     '''
    #     solves a gurobi MILP model to either optimality or time out (timelimit)
    #     returns the objective value and specified value of the decision variables vars
    #     '''
    #     model.Params.Threads = threads
    #     model.Params.TimeLimit = timelimit
    #     model.Params.OutputFlag = 0
    #     model.Params.LogToConsole = 0
     
    #     model.optimize()
    #     solutions = []
    #     if model.Status == GRB.OPTIMAL or (model.Status==GRB.TIME_LIMIT and model.SolCount>0):
    #         score = model.objVal
    #         for v in vars:
    #             solutions.append(model.getAttr('X', v))
 
     
    #         return score, solutions
        
    #     else:
    #          print("warning: model infeasible!")
    #          return np.Inf, solutions 


    # def get_group(self, nxtree, cna_genos, q):
    #     sscn = cna_genos[q]
    #     children_states = {cna_genos[u] for u in nx.descendants(nxtree, q)} - set([sscn])

        
    #     for g, desc in self.group_desc.items():
    #         if sscn == desc['sscn'] and children_states == set(desc['children']):
    #             return g
    #     return None