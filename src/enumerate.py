#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 2021

@author: Palash Sashittal
"""

import gurobipy as gp
import numpy as np
import pandas as pd
import networkx as nx
import itertools



# minimum correction tree parsimonious clone reconciliation problem
class Enumerate:

    def __init__(self, T, S, threads = 1, timelimit = None, verbose = True):  
        # self.snv_mat = snv_mat
        # self.cna_mat = cna_mat
 
        self.cna_clones = {}
        for n in nx.dfs_preorder_nodes(T):
            if T.in_degree[n] ==0:
                T.add_edge(-1, n)
                break
        sorted_nodes = list(T.nodes)
        sorted_nodes.sort()
        snv_mapping= {n: i for i,n in enumerate(sorted_nodes)}
        T= nx.relabel_nodes(T, snv_mapping)
     
        
        self.snv_clones = {val: key for key,val in snv_mapping.items()}
        self.nsnv = len(list(T.nodes))#self.snv_mat.shape[0]
        snv_edges = list(T.edges)

        cna_mapping= {n: i for i,n in enumerate(S.nodes) }
        S= nx.relabel_nodes(S, cna_mapping)
        self.cna_clones = {val: key for key,val in cna_mapping.items()}
        self.ncna = (len(list(S.nodes)))
        cna_edges = list(S.edges)

        G = nx.DiGraph()
        G.add_edges_from(snv_edges)
        self.snv_dag = nx.algorithms.transitive_closure_dag(G)
        G.clear()
        G.add_edges_from(cna_edges)
        self.cna_dag = nx.algorithms.transitive_closure_dag(G)

        self.threads = threads
        self.timelimit = timelimit
        self.verbose = verbose

        self.nsamples  =  1 #self.snv_mat.shape[1]


        self.sol_clones = None
        self.sol_props = None

        # snv/cna parent dictionary
        self.snv_parent_dict = {}
        for edge in snv_edges:
            child = edge[1]
            parent = edge[0]
            if child not in self.snv_parent_dict.keys():
                self.snv_parent_dict[child] = [parent]
            else:
                self.snv_parent_dict[child].append(parent)
        self.cna_parent_dict = {}
        for edge in cna_edges:
            child = edge[1]
            parent = edge[0]
            if child not in self.cna_parent_dict.keys():
                self.cna_parent_dict[child] = [parent]
            else:
                self.cna_parent_dict[child].append(parent)

        self.snv_edges = snv_edges
        self.cna_edges = cna_edges

        for j in range(self.nsnv):
            if j not in self.snv_parent_dict.keys():
                self.snv_root = j
                break
        for k in range(self.ncna):
            if k not in self.cna_parent_dict.keys():
                self.cna_root = k
                break

        print(f'snv root is {self.snv_root} and cna root is {self.cna_root}')

    def solve(self, max_sol=1000):
        model = gp.Model('solveMCTPCR')
        model.Params.PoolSearchMode = 2
        model.Params.PoolSolutions = max_sol
        model.Params.PoolIgnore =  1
        # nsamples =  1 #self.snv_mat.shape[1]
        # assert nsamples == self.cna_mat.shape[1], 'SNV and CNA matrix sizes do not match up.'
        nsamples =self.nsamples
        # nsnv = self.snv_mat.shape[0]
        # ncna = self.cna_mat.shape[0]
        nsnv = self.nsnv
        ncna = self.ncna


        x = model.addVars(nsnv, ncna, vtype=gp.GRB.BINARY, name='x')
        w = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'w')
        y = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'y')
        z_snv = model.addVars(nsnv-1, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'z_snv')
        z_cna = model.addVars(nsnv, ncna-1, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'z_cna')
        d_snv = model.addVars(nsamples, nsnv, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_snv')
        d_cna = model.addVars(nsamples, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_cna')

        # encode product w[i,j,k] = y[i,j,k] * x[j,k]
        # for i in range(nsamples):
        #     for j in range(nsnv):
        #         for k in range(ncna):
        #             model.addConstr(w[i,j,k] <= y[i,j,k])
        #             model.addConstr(w[i,j,k] <= x[j,k])
        #             model.addConstr(w[i,j,k] >= x[j,k] + y[i,j,k] - 1)

        # # encode abundance constraint for snv with correction
        # for i in range(nsamples):
        #     for j in range(nsnv):
        #         sum = gp.LinExpr()
        #         for k in range(ncna):
        #             sum += w[i,j,k]
        #         #model.addConstr(sum == self.snv_mat[j,i])
        #         model.addConstr(self.snv_mat[j,i] - sum <= d_snv[i,j])
        #         model.addConstr(sum - self.snv_mat[j,i] <= d_snv[i,j])

        # # encode abundance constraint for cna
        # for i in range(nsamples):
        #     for k in range(ncna):
        #         sum = gp.LinExpr()
        #         for j in range(nsnv):
        #             sum += w[i,j,k]
        #         #model.addConstr(sum == self.cna_mat[k,i])
        #         model.addConstr(self.cna_mat[k,i] - sum <= d_cna[i,k])
        #         model.addConstr(sum - self.cna_mat[k,i] <= d_cna[i,k])

        # # encode total abundance constraint
        # for i in range(nsamples):
        #     sum = gp.LinExpr()
        #     for j in range(nsnv):
        #         for k in range(ncna):
        #             sum += w[i,j,k]
        #     model.addConstr(sum == 1)

        # encode z_snv[j,k] = x[parent(j), k] * x[j,k]
        # encode z_cna[j,k] = x[j, parent(k)] * x[j,k]
        for edge_idx, edge in enumerate(self.snv_edges):
            parent = edge[0]
            child = edge[1]
            for k in range(ncna):
                model.addConstr(z_snv[edge_idx, k] <= x[parent, k])
                model.addConstr(z_snv[edge_idx, k] <= x[child, k])
                model.addConstr(z_snv[edge_idx, k] >= x[parent, k] + x[child, k] - 1)

        for edge_idx, edge in enumerate(self.cna_edges):
            parent = edge[0]
            child = edge[1]
            for j in range(nsnv):
                model.addConstr(z_cna[j, edge_idx] <= x[j, parent])
                model.addConstr(z_cna[j, edge_idx] <= x[j, child])
                model.addConstr(z_cna[j, edge_idx] >= x[j, parent] + x[j, child] - 1)

        # encode sum_{k} z_snv[j, k] == 1
        # encode sum_{j} z_cna[j, k] == 1
        for edge_idx in range(nsnv-1):
            sum = gp.LinExpr()
            for k in range(ncna):
                sum += z_snv[edge_idx, k]
            #model.addConstr(sum <= 1)
            model.addConstr(sum == 1)
        for edge_idx in range(ncna - 1):
            sum = gp.LinExpr()
            for j in range(nsnv):
                sum += z_cna[j, edge_idx]
            #model.addConstr(sum <= 1)
            model.addConstr(sum == 1)

        # encode x[j,k] <= x[parent(j), k] + x[j, parent(k)]
        for j in range(nsnv):
            for k in range(ncna):
                if j in self.snv_parent_dict.keys() or k in self.cna_parent_dict.keys():
                    sum = gp.LinExpr()
                    if j in self.snv_parent_dict.keys():
                        sum += x[self.snv_parent_dict[j][0], k]
                    if k in self.cna_parent_dict.keys():
                        sum += x[j, self.cna_parent_dict[k][0]]
                    model.addConstr(x[j,k] <= sum)

#         # set objective function
#         obj_sum = gp.LinExpr()
# #        for j in range(nsnv):
# #            for k in range(ncna):
# #                obj_sum += x[j,k]
# #        model.setObjective(obj_sum, gp.GRB.MINIMIZE)
#         for i in range(nsamples):
#             for j in range(nsnv):
#                 obj_sum += d_snv[i,j]
#             for k in range(ncna):
#                 obj_sum += d_cna[i,k]
        model.setObjective(1, gp.GRB.MINIMIZE)

#        model.write('mctpcr.lp')


        model.setParam(gp.GRB.Param.Threads, self.threads)
        model.optimize()


        num_solutions = model.SolCount
        self.trees = []
        labels = [(i,j) for i in range(nsnv) for j in range(ncna)]
        for i in range(num_solutions):
            # model.setParam(gp.GRB.Param.SolutionNumber, i)
            model.setParam(gp.GRB.Param.SolutionNumber, i)
    # Retrieve the values of decision variables for the i-th solution
            x_variables = [var for var in model.getVars() if var.varName.startswith("x")]

# Retrieve the solution values for the filtered decision variables
            solx = model.getAttr('Xn', x_variables)
            self.sol_clones = [labels[i] for i, val in enumerate(solx) if val > 0.5]
            # self.sol_clones = [key for key, val in solx.items() if val >= 0.5]

            # x_variables = [var for var in model.getVars() if var.varName.startswith("x")]

# Retrieve the solution values for the filtered decision variables
# solution_values_for_x = model.getAttr('Xn', x_variables)
            # print(self.sol_clones)
            self.trees.append(self.getCloneTree())

        return self.trees 
        # if model.status == gp.GRB.OPTIMAL:
        #     solx = model.getAttr('x', x)

            # self.sol_props = model.getAttr('x', w)


    def writeCloneFile(self, clone_file, snv_clones = None, cna_clones = None):
        clone_data = []
        for clone in self.sol_clones:
            #for sample in range(self.nsamples):
            if snv_clones:
                snv_clone = snv_clones[clone[0]]
            else:
                snv_clone = clone[0]
            if cna_clones:
                cna_clone = cna_clones[clone[1]]
            else:
                cna_clone = clone[1]

            clone_data.append([clone, snv_clone, cna_clone] + [self.sol_props[sample, clone[0], clone[1]] for sample in range(self.nsamples)])
        df_clone = pd.DataFrame(clone_data, columns=['clone', 'snv_clone', 'cna_clone'] + [f'sample_{idx}' for idx in range(self.nsamples)])
        df_clone.to_csv(clone_file, sep='\t', index=False)

    def writeCloneTree(self,T, clone_tree_file):
    
        clone_edges = list(T.edges)

        with open(clone_tree_file, 'w') as output:
            for clone_edge in clone_edges:
                output.write(f'{clone_edge[0]}\t{clone_edge[1]}\n')
    
    def getCloneTree(self):
        clone_edges = []
        for clone1, clone2 in itertools.permutations(self.sol_clones, 2):
            snv_clone1 = clone1[0]
            snv_clone2 = clone2[0]
            cna_clone1 = clone1[1]
            cna_clone2 = clone2[1]
 
            # else:
     

   
            # else:


            if clone1[0] == clone2[0]:
                if clone1[1] in self.cna_parent_dict.keys():
                    if clone2[1] in self.cna_parent_dict[clone1[1]]:
         
                        snv_clone1 = self.snv_clones[clone1[0]]
                        snv_clone2 = self.snv_clones[clone2[0]]
                        cna_clone1 = self.cna_clones[clone1[1]]
                        cna_clone2 = self.cna_clones[clone2[1]]
                        clone_edges.append(((snv_clone2, cna_clone2), (snv_clone1, cna_clone1)))

            if clone1[1] == clone2[1]:
                if clone1[0] in self.snv_parent_dict.keys():
                    if clone2[0] in self.snv_parent_dict[clone1[0]]:
            
                        snv_clone1 = self.snv_clones[clone1[0]]
                        snv_clone2 = self.snv_clones[clone2[0]]
      
                        cna_clone1 = self.cna_clones[clone1[1]]
                        cna_clone2 = self.cna_clones[clone2[1]]
                        clone_edges.append(((snv_clone2, cna_clone2), (snv_clone1, cna_clone1)))
        

        return nx.DiGraph(clone_edges)

        # with open(clone_tree_file, 'w') as output:
        #     for clone_edge in clone_edges:
        #         output.write(f'{clone_edge[0]}\t{clone_edge[1]}\n')

    # def writeDOT(self, dot_file, snv_clones = None, cna_clones = None):

    #     with open(dot_file, 'w') as output:

    #         output.write(f'digraph N {{\n')
    #         output.write(f"\toverlap=\"false\"\n")
    #         output.write(f"\trankdir=\"TB\"\n")

    #         idx_dict = {}
    #         idx = 0
    #         for clone in self.sol_clones:
    #             if snv_clones:
    #                 snv_clone = snv_clones[clone[0]]
    #             else:
    #                 snv_clone = clone[0]
    #             if cna_clones:
    #                 cna_clone = cna_clones[clone[1]]
    #             else:
    #                 cna_clone = clone[1]

    #             idx_dict[clone] = idx
    #             output.write(f'\t{idx} [label=\"{snv_clone}, {cna_clone}\", style=\"bold\"];\n')

    #             idx += 1

    #         for clone1, clone2 in itertools.permutations(self.sol_clones, 2):
    #             if snv_clones:
    #                 snv_clone1 = snv_clones[clone1[0]]
    #                 snv_clone2 = snv_clones[clone2[0]]
    #             else:
    #                 snv_clone1 = clone1[0]
    #                 snv_clone2 = clone2[0]

    #             if cna_clones:
    #                 cna_clone1 = cna_clones[clone1[1]]
    #                 cna_clone2 = cna_clones[clone2[1]]
    #             else:
    #                 cna_clone1 = clone1[1]
    #                 cna_clone2 = clone2[1]

    #             if clone1[0] == clone2[0]:
    #                 if clone1[1] in self.cna_parent_dict.keys():
    #                     if clone2[1] in self.cna_parent_dict[clone1[1]]:
    #                         output.write(f"\t{idx_dict[clone2]} -> {idx_dict[clone1]} [style=\"bold\"];\n")

    #             if clone1[1] == clone2[1]:
    #                 if clone1[0] in self.snv_parent_dict.keys():
    #                     if clone2[0] in self.snv_parent_dict[clone1[0]]:
    #                         output.write(f"\t{idx_dict[clone2]} -> {idx_dict[clone1]} [style=\"bold\"];\n")

    #         output.write(f'}}')


# import pygraphviz as pgv
# def draw(tree, fname):
#         ptree = pgv.AGraph(strict=False, directed=False)
#         ptree.add_edges_from(list(tree.edges))
#         ptree.layout("dot")
#         ptree.draw(fname)
# # S_edge  = [(0,1), (0,2)]
# S = nx.DiGraph([((1,1), (1,3)), ((1,1), (2,0)) ])
# T = nx.DiGraph([ (0,1), (0,3), (1,2)])
# trees = Enumerate(T,S).solve()
# [draw(trees[i], f"test/refinement{i}.png") for i in range(len(trees))]
# cna_clones = {0: (1,1), 1: (1,3), 2: (2,0)}
# T_edges = [(0, 1), (1,2), (2,3), (1, 4) ]
# snv_clones = {0: -1, 1: 0, 2: 1, 3: 2, 4:3}

# trees = obj.solve(cna_clones=cna_clones, snv_clones=snv_clones)

# obj.writeCloneTree("test/clone_tree.txt")

