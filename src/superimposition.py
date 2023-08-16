import gurobipy as gp
import numpy as np

import networkx as nx
import itertools


class Superimposition:
    def __init__(self, G:nx.digraph, T1:nx.digraph, T2:nx.digraph,  threads = 1, timelimit = None, verbose = True) -> None:
        self.G = G 
        self.T1 = T1
        self.T2 = T2 
        self.threads = threads 
        self.timelimit = timelimit
        self.verbose  = verbose 

        self.T1_clones = list(self.T1.nodes) 
        self.T2_clones = list(self.T2.nodes)
        self.K1 = len(self.T1)
        self.K2 = len(self.T2)
        self.edges = list(self.G.edges)

    def solve(self):
        model = gp.Model('mergeTrees')
        
        x = model.addVars(self.K1, self.K2, vtype=gp.GRB.BINARY, name='x')
        y = model.addVars(len(self.edges), vtype=gp.GRB.CONTINUOUS, name='y')
        z1 = model.addVars(len(list(self.T1.edges)), self.K2, vtype=gp.GRB.CONTINUOUS, name='z1')
        z2 = model.addVars(self.K1, len(list(self.T2.edges)), vtype=gp.GRB.CONTINUOUS, name='z2')

               # encode z1[j,k] = x[parent(j), k] * x[j,k]
                #encode z2[j,k] = x[j, parent(k)] * x[j,k]
        self.t1_edge_index = {}
        for edge_idx, edge in enumerate(self.T1.edges):
            parent = edge[0]
            child = edge[1]
            self.t1_edge_index[child] = edge_idx
            for k in range(self.K2):
             
                model.addConstr(z1[edge_idx, k] <= x[parent, k])
                model.addConstr(z1[edge_idx, k] <= x[child, k])
                model.addConstr(z1[edge_idx, k] >= x[parent, k] + x[child, k] - 1)
        
        self.t2_edge_index = {}
        for edge_idx, edge in enumerate(self.T2.edges):
            parent = edge[0]
            child = edge[1]
            self.t2_edge_index[child] = edge_idx
            for j in range(self.K1):
                model.addConstr(z2[j, edge_idx] <= x[j, parent])
                model.addConstr(z2[j, edge_idx] <= x[j, child])
                model.addConstr(z2[j, edge_idx] >= x[j, parent] + x[j, child] - 1)
        
        
        for i, i_clone in enumerate(self.T1_clones):
            for j, j_clone in enumerate(self.T2_clones):
                out = (0,i_clone)
                in_ = (1,j_clone)
                j_edge_index = self.t2_edge_index[j]
                if (out, in_) in self.G.edges:
                    model.addConstr(y[i,j]==z2[i, j_edge_index])
    
        
        for i, i_clone in enumerate(self.T1_clones):
            for j, j_clone in enumerate(self.T2_clones):
                out = (0,i_clone)
                in_ = (1,j_clone)
                i_edge_index = self.t1_edge_index[i]
                if (out, in_) in self.G.edges:
                    model.addConstr(y[j,i]==z1[i_edge_index, j])
        




        




