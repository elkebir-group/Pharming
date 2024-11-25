#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 2021

@author: Palash Sashittal
"""



import networkx as nx
import itertools
from pyomo.environ import *
from pyomo.opt import SolverFactory


# minimum correction tree parsimonious clone reconciliation problem
class Enumerate:

    def __init__(self, T, S, threads = 1, timelimit = None, verbose = False, same_root=False):  
   
 
        self.cna_clones = {}
        if not same_root:
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

        # print(f'snv root is {self.snv_root} and cna root is {self.cna_root}')

    # def solve_gurobi(self, max_sol=10000):
    #     import gurobipy as gp
    #     model = gp.Model('solveMCTPCR')
    #     model.Params.PoolSearchMode = 2
    #     model.Params.PoolSolutions = max_sol

    #     nsamples =self.nsamples
 
    #     nsnv = self.nsnv
    #     ncna = self.ncna


    #     x = model.addVars(nsnv, ncna, vtype=gp.GRB.BINARY, name='x')
    #     w = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'w')
    #     y = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'y')
    #     z_snv = model.addVars(nsnv-1, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'z_snv')
    #     z_cna = model.addVars(nsnv, ncna-1, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'z_cna')
    #     d_snv = model.addVars(nsamples, nsnv, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_snv')
    #     d_cna = model.addVars(nsamples, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_cna')


    #     for edge_idx, edge in enumerate(self.snv_edges):
    #         parent = edge[0]
    #         child = edge[1]
    #         for k in range(ncna):
    #             model.addConstr(z_snv[edge_idx, k] <= x[parent, k])
    #             model.addConstr(z_snv[edge_idx, k] <= x[child, k])
    #             model.addConstr(z_snv[edge_idx, k] >= x[parent, k] + x[child, k] - 1)

    #     for edge_idx, edge in enumerate(self.cna_edges):
    #         parent = edge[0]
    #         child = edge[1]
    #         for j in range(nsnv):
    #             model.addConstr(z_cna[j, edge_idx] <= x[j, parent])
    #             model.addConstr(z_cna[j, edge_idx] <= x[j, child])
    #             model.addConstr(z_cna[j, edge_idx] >= x[j, parent] + x[j, child] - 1)

    #     # encode sum_{k} z_snv[j, k] == 1
    #     # encode sum_{j} z_cna[j, k] == 1
    #     for edge_idx in range(nsnv-1):
    #         sum = gp.LinExpr()
    #         for k in range(ncna):
    #             sum += z_snv[edge_idx, k]
    #         #model.addConstr(sum <= 1)
    #         model.addConstr(sum == 1)
    #     for edge_idx in range(ncna - 1):
    #         sum = gp.LinExpr()
    #         for j in range(nsnv):
    #             sum += z_cna[j, edge_idx]
    #         #model.addConstr(sum <= 1)
    #         model.addConstr(sum == 1)

    #     # encode x[j,k] <= x[parent(j), k] + x[j, parent(k)]
    #     for j in range(nsnv):
    #         for k in range(ncna):
    #             if j in self.snv_parent_dict.keys() or k in self.cna_parent_dict.keys():
    #                 sum = gp.LinExpr()
    #                 if j in self.snv_parent_dict.keys():
    #                     sum += x[self.snv_parent_dict[j][0], k]
    #                 if k in self.cna_parent_dict.keys():
    #                     sum += x[j, self.cna_parent_dict[k][0]]
    #                 model.addConstr(x[j,k] <= sum)


    #     model.setObjective(1, gp.GRB.MINIMIZE)


    #     model.setParam(gp.GRB.Param.Threads, self.threads)
    #     model.optimize()


    #     num_solutions = model.SolCount
    #     self.trees = []
    #     labels = [(i,j) for i in range(nsnv) for j in range(ncna)]
    #     for i in range(num_solutions):
    #         # model.setParam(gp.GRB.Param.SolutionNumber, i)
    #         model.setParam(gp.GRB.Param.SolutionNumber, i)
    #     # Retrieve the values of decision variables for the i-th solution
    #         x_variables = [var for var in model.getVars() if var.varName.startswith("x")]

    #     # Retrieve the solution values for the filtered decision variables
    #         solx = model.getAttr('Xn', x_variables)
    #         self.sol_clones = [labels[i] for i, val in enumerate(solx) if val > 0.5]

    #         self.trees.append(self.getCloneTree())

    #     return self.trees 
 
    
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

 


    def solve(self, max_sol=10000):
        # Create Pyomo model
        model = ConcreteModel()

        nsamples = self.nsamples
        nsnv = self.nsnv
        ncna = self.ncna

        # Define sets
        model.SNV = RangeSet(0, nsnv - 1)
        model.CNA = RangeSet(0, ncna - 1)
        model.SAMPLES = RangeSet(0, nsamples - 1)
        model.SNV_EDGES = RangeSet(0, len(self.snv_edges) - 1)
        model.CNA_EDGES = RangeSet(0, len(self.cna_edges) - 1)

        # Define variables
        model.x = Var(model.SNV, model.CNA, within=Binary)
        model.w = Var(model.SAMPLES, model.SNV, model.CNA, bounds=(0, 1))
        model.y = Var(model.SAMPLES, model.SNV, model.CNA, bounds=(0, 1))
        model.z_snv = Var(model.SNV_EDGES, model.CNA, bounds=(0, 1))
        model.z_cna = Var(model.SNV, model.CNA_EDGES, bounds=(0, 1))
        model.delta_snv = Var(model.SAMPLES, model.SNV, bounds=(0, 1))
        model.delta_cna = Var(model.SAMPLES, model.CNA, bounds=(0, 1))

        # Add constraints for SNV edges
        def snv_edge_constraint_rule(model, edge_idx, k):
            parent, child = self.snv_edges[edge_idx]
            return [
                model.z_snv[edge_idx, k] <= model.x[parent, k],
                model.z_snv[edge_idx, k] <= model.x[child, k],
                model.z_snv[edge_idx, k] >= model.x[parent, k] + model.x[child, k] - 1
            ]
        model.snv_edge_constraints = ConstraintList()
        for edge_idx, edge in enumerate(self.snv_edges):
            for k in range(ncna):
                for constr in snv_edge_constraint_rule(model, edge_idx, k):
                    model.snv_edge_constraints.add(constr)

        # Add constraints for CNA edges
        def cna_edge_constraint_rule(model, j, edge_idx):
            parent, child = self.cna_edges[edge_idx]
            return [
                model.z_cna[j, edge_idx] <= model.x[j, parent],
                model.z_cna[j, edge_idx] <= model.x[j, child],
                model.z_cna[j, edge_idx] >= model.x[j, parent] + model.x[j, child] - 1
            ]
        model.cna_edge_constraints = ConstraintList()
        for edge_idx, edge in enumerate(self.cna_edges):
            for j in range(nsnv):
                for constr in cna_edge_constraint_rule(model, j, edge_idx):
                    model.cna_edge_constraints.add(constr)

        # Sum constraints for z_snv
        model.z_snv_sum_constraints = ConstraintList()
        for edge_idx in range(nsnv - 1):
            model.z_snv_sum_constraints.add(
                sum(model.z_snv[edge_idx, k] for k in range(ncna)) == 1
            )

        # Sum constraints for z_cna
        model.z_cna_sum_constraints = ConstraintList()
        for edge_idx in range(ncna - 1):
            model.z_cna_sum_constraints.add(
                sum(model.z_cna[j, edge_idx] for j in range(nsnv)) == 1
            )

        # Constraints for x
        model.x_constraints = ConstraintList()
        for j in range(nsnv):
            for k in range(ncna):
                if j in self.snv_parent_dict.keys() or k in self.cna_parent_dict.keys():
                    expr = 0
                    if j in self.snv_parent_dict.keys():
                        expr += model.x[self.snv_parent_dict[j][0], k]
                    if k in self.cna_parent_dict.keys():
                        expr += model.x[j, self.cna_parent_dict[k][0]]
                    model.x_constraints.add(model.x[j, k] <= expr)

        # Objective function
        model.objective = Objective(expr=1, sense=minimize)

        # Solve the model
        opt = SolverFactory('glpk')
        result = opt.solve(model)


        # Extract solutions
        self.trees = []

        for _ in range(max_sol):
            result = opt.solve(model)
            if result.solver.termination_condition != TerminationCondition.optimal:
                break  # No more solutions

            # Extract solution
            x_solution = {(i, j): model.x[i, j].value for i in range(nsnv) for j in range(ncna)}
            sol_clones = [(i, j) for (i, j), val in x_solution.items() if val > 0.5]
            self.sol_clones = sol_clones
            self.trees.append(self.getCloneTree())

            # Add no-good cut to exclude the current solution
            model.x_constraints.add(
                sum(1 - model.x[i, j] if val > 0.5 else model.x[i, j]
                    for (i, j), val in x_solution.items()) >= 1
            )

       
        return self.trees