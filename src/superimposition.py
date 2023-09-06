import gurobipy as gp 
from gurobipy import GRB
import numpy as np
from clonal_tree_new import ClonalTreeNew
import networkx as nx



class Superimposition:
    def __init__(self,T1, T2, verbose = True) -> None:
        
        self.T1 = T1
        self.T2 = T2
        self.T1_root =self.T1.find_root()
        self.T2_root = self.T2.find_root()

        self.verbose  = verbose 

        self.T1_clones = self.T1.clones()
        self.T2_clones = self.T2.clones()
        self.T1_edges = self.T1.edges()
        self.T2_edges = self.T2.edges()

        self.T = nx.DiGraph()
        self.genotypes = {}
        self.cell_mapping = {}


        self.PI = [ (u,v) for u in self.T1_clones for v in self.T2_clones]



    def solve(self, data, threads=1, timelimit=500):

        # rng = np.random.default_rng(12345)
        model = gp.Model('mergeTrees')
        if not self.verbose:
            model.Params.LogToConsole = 0
        
        model.Params.Threads = threads 
        model.Params.TimeLimit = timelimit
        
        cells = data.cells 

        M = len(cells)
        
        x = model.addVars(self.T1_clones, self.T2_clones, vtype=GRB.BINARY, name='x')
        z1 = model.addVars(self.T1_edges, self.T2_clones, vtype=GRB.CONTINUOUS, name='z1',lb=0, ub=1)
        z2 = model.addVars(self.T1_clones, self.T2_edges, vtype=GRB.CONTINUOUS, name='z2', lb=0, ub=1)
        phi = model.addVars(cells, self.T1_clones, self.T2_clones, vtype = GRB.CONTINUOUS, name='phi', lb=0, ub=1)
        y = model.addVars(self.T1_clones, self.T2_clones, vtype=GRB.BINARY, name='y')

        
        for i in cells: 
            #every cell is assigned to exactly one clone
            model.addConstr(sum(phi[i,u,v] for u,v in self.PI)==1)
            
            #cells can only be assigned t o clones included in the refined tree
            model.addConstrs(phi[i,u,v] <= x[u,v] for u,v in self.PI)

        #tree refinement constraints
        for v in self.T2_clones:
            for u, u_prime in self.T1_edges:
            
                model.addConstr(z1[u, u_prime,v] <= x[u,v])
                model.addConstr(z1[u, u_prime,v] <= x[u_prime, v])
                model.addConstr(z1[u, u_prime,v] >= x[u_prime, v] + x[u, v] - 1)
        
        for u in self.T1_clones:
            for v, v_prime in self.T2_edges:
  
                model.addConstr(z2[u,v, v_prime] <= x[u,v])
                model.addConstr(z2[u,v, v_prime] <= x[u, v_prime])
                model.addConstr(z2[u,v, v_prime] >= x[u,v] + x[u, v_prime] - 1)


        for u, u_prime in self.T1_edges:
            model.addConstr(sum(z1[u, u_prime,v] for v in self.T2_clones) ==1)
        
        for v, v_prime in self.T2_edges:
            model.addConstr(sum(z2[u,v, v_prime] for u in self.T1_clones) ==1)

        # encode x[u,v] <= x[parent(u), v] + x[u, parent(v)] 
        for u,v in self.PI:
            if u != self.T1_root or v != self.T2_root:
                sm = gp.LinExpr()
                if u !=  self.T1_root:
                    u_parent = self.T1.parent(u)
                    sm += x[u_parent, v]
                if v != self.T2_root:
                    v_parent = self.T2.parent(v)
                    sm += x[u, v_parent]
                # if u != self.T1_root and v != self.T2_edges:
                #     sm += x[u_parent, v_parent]
                model.addConstr(x[u,v] <=  sm)

        
        cst ={}
        for u in self.T1_clones:
            u_costs = self.T1.node_snv_cost(u,cells, data)
            for v in self.T2_clones:
                v_costs = self.T2.node_snv_cost(v, cells, data)
                tot = u_costs + v_costs 
                for i,c in enumerate(cells):
                    cst[c,u,v]  = tot[i] 
        
        # for u,v in self.PI:
        #     model.addConstr(-1*sum(phi[i,u,v] for i in cells) - M*y[u,v] <= -40)
        #     model.addConstr(sum(phi[i,u,v] for i in cells) - M*(1-y[u,v]) <= 0)
           
        # for u,v in self.PI:
        #     if (u ==0 and v ==0) or v !=u: 
        #         continue
        
        #     model.addConstr(sum(phi[i,u,v] for i in cells) ==0)


        #DEBUGGING
        # myclones = [(0,0), (0,4), (4,4), (4,3), (3,3), (4,2), (2,2) ,(1,0), (1,1)]
        
        # model.addConstrs(x[u,v]==1 for u,v in myclones)

        # rest_clones = [(0,4),(4,3), (4,2), (1,0)]
        # #rest_clones = [(4,3), (4,2)]
        # model.addConstrs(sum(phi[i,u,v] for i in cells)==0 for u,v in rest_clones)

        model.setObjective( sum(cst[i,u,v]* phi[i,u,v] for i in cells for u,v in self.PI), GRB.MINIMIZE)
        # # + 1000*sum(x[u,v] for u,v in self.PI

        #solve ILP
        model.optimize()
        
        if model.Status == GRB.OPTIMAL:
                clones= model.getAttr('X', x)
                cell_assign = model.getAttr('X', phi)
                z1_vals = model.getAttr('X', z1)
                z2_vals = model.getAttr('X', z2)

                score = model.objVal
                # running_count = 0
                print(f"score: {score}")
                mapping = {}
                orig_clones = []
                index  =0 
                assigned_cells = 0
                for u,v in self.PI:  
                    if clones[u,v] > 0.5:
                        self.cell_mapping[index] = []
                        mapping[(u,v)] = index
                        orig_clones.append((u,v))
                        self.T.add_node(index)
                        #TODO: fix seg ids 
                        self.genotypes[index] = {self.T2.key: self.T2.genotypes[v][self.T2.key],self.T1.key: self.T1.genotypes[u][self.T1.key] }
                        for i in cells:
                            
                            if cell_assign[i,u,v] > 0.5:
                                assigned_cells += 1
                                # if u==3 and v ==0:
                                #     running_count+=1 
                                self.cell_mapping[index].append(i)
                        index+=1
          
                for u,v in orig_clones:
                    
                    if u == self.T1_root and v ==self.T2_root:
                        continue
                    u_prime = self.T1.parent(u)
                    v_prime = self.T2.parent(v)
                    if (u_prime, v) in orig_clones:
                        self.T.add_edge(mapping[u_prime,v], mapping[u,v])
                    elif (u, v_prime) in orig_clones:
                        self.T.add_edge(mapping[u, v_prime], mapping[u,v])
                    else:
                        self.T.add_edge(mapping[u_prime, v_prime], mapping[u,v])






        else:
            raise ValueError("Model is Infeasible!")

        self.CT = ClonalTreeNew(self.T, self.genotypes, self.cell_mapping, key=(self.T1.key, self.T2.key), cost=score)
        return score, self.CT, mapping 
    

            
   



                

       
         
      
        
   

                     

  
  


        






        




