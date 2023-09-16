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

        self.T1_leafs = self.T1.get_leaves()
        self.T2_leafs = self.T2.get_leaves()

        self.T = nx.DiGraph()
        self.mapping = {}



        self.PI = [ (u,v) for u in self.T1_clones for v in self.T2_clones]
      

    def solve_refinement(self, K):
        model = gp.Model('mergeTrees')
       
        
        model.Params.Threads = self.threads 
        model.Params.TimeLimit = self.timelimit
        if not self.verbose:
            model.Params.LogToConsole = 0
        

        x = model.addVars(self.T1_clones, self.T2_clones, vtype=GRB.BINARY, name='x')
        z1 = model.addVars(self.T1_edges, self.T2_clones, vtype=GRB.CONTINUOUS, name='z1',lb=0, ub=1)
        z2 = model.addVars(self.T1_clones, self.T2_edges, vtype=GRB.CONTINUOUS, name='z2', lb=0, ub=1)
        phi = model.addVars(self.cells, self.T1_clones, self.T2_clones, vtype = GRB.CONTINUOUS, name='phi', lb=0, ub=1)
        y = model.addVars(self.T1_clones, self.T2_clones, vtype=GRB.BINARY, name='y')

        
        for i in self.cells: 
            #every cell is assigned to exactly one clone
            model.addConstr(sum(phi[i,u,v] for u,v in self.PI)==1)
            
            #cells can only be assigned t o clones included in the refined tree
            model.addConstrs(phi[i,u,v] <= x[u,v] for u,v in self.PI)
            model.addConstrs(phi[i,u,v] <= y[u,v] for u,v in self.PI)
        
        model.addConstrs(y[u,v] <= x[u,v] for u,v in self.PI)

        model.addConstr(sum(y[u,v] for u,v in self.PI) ==K)

        # every clone in each of the trees must appear at least once 
        model.addConstrs(sum(phi[i,u,v] for i in self.cells for u in self.T1_clones) >=1 for v in self.T2_leafs)
        model.addConstrs(sum(phi[i,u,v] for i in self.cells for v in self.T2_clones)  >=1 for u in self.T1_leafs)


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

        model.setObjective( sum(self.cst[i,u,v]* phi[i,u,v] for i in self.cells for u,v in self.PI), GRB.MINIMIZE)

        model.optimize()
        # if True:
        if model.SolCount > 0:
        # model.Status in [GRB.OPTIMAL, GRB.TIME_LIMIT]:
            clones= model.getAttr('X', x)
            cell_assign = model.getAttr('X', phi)

            score = model.objVal
            
            self.mapping, cell_mapping, genotypes = {}, {}, {}
            orig_clones = []
            index  =0 
            assigned_cells = 0
            for u,v in self.PI:  
                if clones[u,v] > 0.5:
                    cell_mapping[(u,v)] = []
                    self.mapping[(u,v)] = index
                    orig_clones.append((u,v))
                    self.T.add_node((u,v))
                    #TODO: fix seg ids 
                    for i in self.cells:
                        # genotypes[(u,v)] = {self.T2.key: self.T2.genotypes[v][self.T2.key],self.T1.key: self.T1.genotypes[u][self.T1.key] }

                        if cell_assign[i,u,v] > 0.5:
                            assigned_cells += 1
                            cell_mapping[(u,v)].append(i)
                    index+=1
        
            for u,v in orig_clones:
                
                if u == self.T1_root and v ==self.T2_root:
                    continue
                u_prime = self.T1.parent(u)
                v_prime = self.T2.parent(v)
                if (u_prime, v) in orig_clones:
                    self.T.add_edge((u_prime,v), (u,v))
                elif (u, v_prime) in orig_clones:
                    self.T.add_edge((u, v_prime), (u,v))
                else:
                    # if (u_prime, v_prime) not in orig_clones:


                    self.T.add_edge(self.mapping[u_prime, v_prime], self.mapping[u,v])
            for u,v in self.T:
                genotypes[(u,v)] = {self.T2.key: self.T2.genotypes[v][self.T2.key],self.T1.key: self.T1.genotypes[u][self.T1.key] }

            print(f"clones: {len(self.T)} assigned_cells: {assigned_cells} score: {score}")
            
        else:
            print(f"Warning, model is infeasible for K={K}")
            return ClonalTreeNew(nx.DiGraph(), {}, cost=np.Inf)
            # raise ValueError("Model is Infeasible!")
        
        return ClonalTreeNew(self.T.copy(), genotypes, cell_mapping, cost=score)
           

    def solve_pruning(self, RT,lamb=100, threshold=10):
        leafs = [self.rev_mapping[l] for l in RT.get_leaves()]
        orig_clones = [self.rev_mapping[n] for n in  RT.clones()]

        pruning = gp.Model('pruneTrees')
        pruning.Params.Threads = self.threads 
        pruning.Params.TimeLimit = self.timelimit
        if not self.verbose:
            pruning.Params.LogToConsole = 0
        
        x = pruning.addVars(self.T1_clones, self.T2_clones, vtype=GRB.BINARY, name='x')
        phi = pruning.addVars(self.cells, self.T1_clones, self.T2_clones, vtype = GRB.CONTINUOUS, name='phi', lb=0, ub=1)

        for u,v in self.PI:
            if (u,v) not in orig_clones:
                pruning.addConstr(x[u,v]==0)
        
        for i in self.cells: 
            #every cell is assigned to exactly one clone
            pruning.addConstr(sum(phi[i,u,v] for u,v in self.PI)==1)
        
            #cells can only be assigned t o clones included in the refined tree
            pruning.addConstrs(phi[i,u,v] <= x[u,v] for u,v in self.PI)

        #leaf nodes must have cells 
        pruning.addConstrs(sum(phi[i,u,v] for i in self.cells) >= threshold for u,v in leafs)

        # every clone in each of the trees must appear at least once 
        pruning.addConstrs(sum(x[u,v] for u in self.T1_clones) >=1 for v in self.T2_clones)
        pruning.addConstrs(sum(x[u,v] for v in self.T2_clones)  >=1 for u in self.T1_clones)
        
        
        pruning.setObjective( sum(self.cst[i,u,v]* phi[i,u,v] for i in self.cells for u,v in self.PI) \
                              + lamb*sum(x[u,v] for u,v in self.PI), GRB.MINIMIZE)

        pruning.optimize()
        if True:
        # if pruning.Status == GRB.OPTIMAL:
            clones= pruning.getAttr('X', x)
            cell_assign = pruning.getAttr('X', phi)
            cell_mapping = {}
            score, assigned_cells = 0, 0
            for u,v in orig_clones:
                if clones[u,v] < 0.5:
                    RT.prune_tree(self.mapping[u,v])
                else:
                    cell_mapping[self.mapping[u,v]] = []
                    for i in self.cells:
                        if cell_assign[i,u,v] > 0.5:
                            assigned_cells += 1
                            score += self.cst[i,u,v]
                            cell_mapping[self.mapping[u,v]].append(i)
            print(f"clones: {len(self.T)} assigned_cells: {assigned_cells} score: {score}")


        else:
            raise ValueError("Model is Infeasible!")
        RT.set_cell_mapping(cell_mapping) 
        RT.set_cost(score)

    

    def solve(self, data, K, lamb=0, threshold=10, threads=1, timelimit=100):
        self.threads = threads 
        self.timelimit =  timelimit
   

        self.cells = data.cells 
        self.M = len(self.cells)

        self.cst ={}
        for u in self.T1_clones:
            u_costs = self.T1.node_snv_cost(u,self.cells, data)
            for v in self.T2_clones:
                v_costs = self.T2.node_snv_cost(v, self.cells, data)
                tot = u_costs + v_costs 
                for i,c in enumerate(self.cells):
                    self.cst[c,u,v]  = tot[i] 
     
        RT = self.solve_refinement(K)
    

        self.rev_mapping = {val: key for key,val in self.mapping.items()}
        # RT.draw("test/refined_tree.png", self.rev_mapping)

        if lamb > 0:

            self.solve_pruning(RT, lamb, threshold)
            # RT.draw("test/refined_tree.png", self.rev_mapping)

        return RT, self.rev_mapping

