import gurobipy as gp 
from gurobipy import GRB
import numpy as np

class CellAssignment:
    def __init__(self) -> None:
        self.model = gp.Model("MIQP")
    

    def build(self, A, D, t, mutations= None, init_phi=None, z_vals= None, y_vals=None):
        '''
        A and D are n x m count matrices of variant and total reads respectively
        t is a dictionary of dictionaries with clone as key and mutation as key
        '''
        self.nnodes = len(t)
        self.n, self.m = A.shape
        self.clones = list(t.keys())


        if mutations is None:
            self.muts = list(t[self.clones[0]])
        else:
            self.muts = mutations 
        target_dict = {}
        for k in t:
            for m in self.muts:
                target_dict[k,m] = t[k][m]
        
        indices, target = gp.multidict(target_dict)



    


        self.cluster_assign = [(i,k) for i in range(self.n) for k in self.clones]
        # self.y_indices = [(k,j)   for k in self.clones for j in self.muts]
        self.x = self.model.addVars(self.cluster_assign, vtype = GRB.BINARY )
        self.z = self.model.addVars(indices, lb=0.0, ub=1.0)
        self.y = self.model.addVars(indices, lb=0.0, ub=1.0)
        # self.q =self.model.addVars(self.clones,vtype=GRB.BINARY)

        for k,j in z_vals:
            self.z[k,j].start = z_vals[k,j]
            # self.y[k,j].start = y_vals[k,j]

        if init_phi is not None:
            for i,k in init_phi.items():
                for c in self.clones:
                    if c ==k:
                        self.x[i,c].start =1
                    else:
                        self.x[i,c].start =0

        self.model.setObjective(gp.quicksum(self.z[(k,j)]  for k,j in indices), gp.GRB.MINIMIZE)

        # self.model.addConstrs(self.q[k] >= self.x[i,k] for i,k in self.cluster_assign)
        #each cell is assigned to exactly one cluster 
        self.model.addConstrs(gp.quicksum(self.x[i,k] for k in self.clones) ==1 for i in range(self.n))
        
        # every cluster has at least x cells assigned 
        for k in self.clones:
            self.model.addConstr(gp.quicksum(self.x[i,k] for i in range(self.n)) >= 25)

        for k,j in indices:
            self.model.addConstr(gp.quicksum(D[i, j] * self.x[i, k] for i in range(self.n)) >= 1e-4)

            self.model.addConstr(self.y[k,j]* gp.quicksum(D[i,j]*self.x[i,k] for i in range(self.n))==1)
            self.model.addConstr(self.z[k,j]  >= (self.y[k,j]*gp.quicksum(A[i,j] * self.x[i, k] for i in range(self.n)))- target[k,j])
            self.model.addConstr(self.z[k, j] >= -1*((self.y[k,j]*gp.quicksum(A[i,j] * self.x[i, k] for i in range(self.n)))  -target[k,j]))

        # for k,j in indices:

        #     self.model.addConstr(self.y[k,j]* gp.quicksum(D[i,j]*self.x[i,k] for i in range(self.n))==self.q[k])
        #     self.model.addConstr(self.z[k,j]  >= (self.y[k,j]*gp.quicksum(A[i,j] * self.x[i, k] for i in range(self.n)))- self.q[k]*target[k,j])
        #     self.model.addConstr(self.z[k, j] >= -1*((self.y[k,j]*gp.quicksum(A[i,j] * self.x[i, k] for i in range(self.n)))  - self.q[k]*target[k,j]))


    def solve(self, threads=1, timelimit=500):
        self.model.Params.Threads = threads
        self.model.Params.NonConvex = 2
        self.model.Params.TimeLimit = timelimit
        self.model.Params.StartNodeLimi = 20000000
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


        for i,k in self.cluster_assign:
    
                if solution[i,k] > 0.5:
                    phi[i] =k
    
        
        return score, phi 

        

       



# #test cases 

# A = np.array([[0,1,1], [0,2,1], [1,3,3], [0,0,0], [0,0,0]])
# D = np.full((5,5), 3)

# target = np.array([[0,0,0,0],[0,0,0.5, 0.5], [0,1/3, 0,0], [1/3, 1/3,0,0]])

# ca = CellAssignment()
# ca.build(A,D, target)
# ca.solve()




        

