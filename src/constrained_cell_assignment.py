from clonal_tree import ClonalTree
import gurobipy as gp 
from gurobipy import GRB
import numpy as np 

class ConstrainedCellAssign:
    def __init__(self, tree, cell_scores, nodes, dcfs, segment) -> None:
              
        self.model = gp.Model("MIP")

    
        cna_genos = tree.get_cna_genos()[segment]
  
        self.clones = nodes 

        self.n = cell_scores.shape[0]
        self.cells = [i for i in range(self.n)]
        self.clone_assign = [(i,u) for i in self.cells for u in range(cell_scores.shape)]
        self.x = self.model.addVars(self.clone_assign, vtype = GRB.BINARY )

        self.model.setObjective(gp.quicksum(cell_scores[i,u]*self.x[i,u] for i,u in self.clone_assign), gp.GRB.MINIMIZE)
        
        for q in range(cell_scores.shape[0]):
            # desc = tree.get_descendants(u)
            if tree.is_leaf(nodes[q]):
                cna_state = cna_genos[nodes[q]].to_tuple()

    
                self.model.addConstr( gp.quicksum(self.x[i,q] for i in self.cells) >= self.n *dcfs[nodes[q]] )


    def solve(self, threads=1, timelimit=60):
        self.model.Params.Threads = threads
        self.model.Params.TimeLimit = timelimit
     
        self.model.optimize()
        phi = {}
        if self.model.Status == GRB.OPTIMAL or (self.model.Status==GRB.TIME_LIMIT and self.model.SolCount>0):
            solution = self.model.getAttr('X', self.x)
 
            score = self.model.objVal
        
        # elif self.model.Status == GRB.TIME_LIMIT:
        #     score = self.model.getBestObjVal()
        #     solution = self.model.getBestSol()

            
 
        else:
             print("warning: model infeasible!")
             return np.Inf, phi


        for i,q in self.clone_assign_assign:
    
                if solution[i,q] > 0.5:
                    phi[i] =self.clones[q]
    
        
        return score, phi 
