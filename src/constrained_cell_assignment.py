
import gurobipy as gp 
from gurobipy import GRB
import numpy as np 
from cell_mapping import CellAssign

class ConstrainedCellAssign:
    def __init__(self, tree, cell_scores, nodes, dcfs_lb, dcfs_ub=None) -> None:
              
        self.model = gp.Model("MIP")

    
        # cna_genos = tree.get_cna_genos()[segment]
  
        self.clones = nodes 

        self.n = cell_scores.shape[1]
        self.cells = [i for i in range(self.n)]
        self.clone_assign = [(i,u) for i in self.cells for u in range(cell_scores.shape[0])]
        self.x = self.model.addVars(self.clone_assign, vtype = GRB.BINARY )

        self.model.setObjective(gp.quicksum(cell_scores[u,i]*self.x[i,u] for i,u in self.clone_assign), gp.GRB.MINIMIZE)
        
        #every cells is assigned to exactly 1 node
        self.model.addConstrs(gp.quicksum(self.x[i,q] for q in range(cell_scores.shape[0]))==1 for i in self.cells)
        node_list = nodes.tolist()
        for q in range(cell_scores.shape[0]):
            if nodes[q] in dcfs_lb:
                desc =tree.preorder(nodes[q])
           
                desc_indices = [node_list.index(u) for u in desc]
            # if tree.is_leaf(nodes[q]) and nodes[q] in dcfs:
              
                
                # self.model.addConstr( gp.quicksum(self.x[i,q] for i in self.cells for q in desc_indices) <= self.n *dcfs_ub[nodes[q]] )

    
                self.model.addConstr( gp.quicksum(self.x[i,q] for i in self.cells for q in desc_indices) >= self.n *dcfs_lb[nodes[q]] )
                if dcfs_ub is not None:
                    self.model.addConstr( gp.quicksum(self.x[i,q] for i in self.cells for q in desc_indices) <= self.n *dcfs_ub[nodes[q]] )


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


        for i,q in self.clone_assign:
    
                if solution[i,q] > 0.5:
                    phi[i] =self.clones[q]
    
 
        return score,  CellAssign(phi, set(self.clones))
