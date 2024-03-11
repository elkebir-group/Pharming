
# Created by: L.L. Weber
# Created on: 2023-11-27 11:48:53

import numpy as np
from data import Data, load_from_pickle
from clonal_tree import ClonalTree
import pandas as pd

class Eval:
    def __init__(self, gt) -> None:
        self.gt_phi = gt.get_phi()
    def accuracy(self, inf_phi):
        ncells = len(self.gt_phi)
        total =0
        for i in self.gt_phi:
            total += self.gt_phi[i] == inf_phi[i]
        
        return total/ncells 
        # merged_phi = pd.concat([pd.Series(self.gt_phi), pd.Series(inf_phi)], axis=1)
        # merged_phi.columns = ["gt", "obj1"]

        # merged_phi["close"] =0
        # merged_phi["min3_val"] = False


        

class CellAssignEvol:
    def __init__(self, population_size=400, tol=0.01, parents=50, iterations=100, seed=102, mutation_rate=0.1) -> None:
        self.pop_size = population_size -1
        self.tol = tol
        self.iter = iterations
        self.rng = np.random.default_rng(seed)
        self.parents = parents
        self.mut_rate = mutation_rate
        self.overall_fitness = []
        self.overall_assign = []
    
    def create_individual(self):
        cell_assign = self.rng.choice(self.clones, self.data.N, replace=True)
        phi = {i: cell_assign[i] for i in self.data.cells}
        return phi

        # phi = {}
        # for i, assign in self.base_phi.items():
        #     if self.rng.random() < self.mut_rate:
        #         phi[i]  = self.rng.choice(self.clones)
        #     else:
        #         phi[i] = assign
        # return phi
    def initialize(self):
        self.populations = []
    
        self.populations.append(self.tree.get_phi())

        for p in range(self.pop_size):
            phi = self.create_individual()
        
            self.populations.append(phi)
    
    def add_new_parents(self):
        current_size = len(self.populations)
        for p in range(current_size, self.pop_size):
            phi= self.create_individual()
            self.populations.append(phi)


    def eval_fitness(self):
        self.fitness = []
        min_fitness = np.Inf
        for j,phi in enumerate(self.populations):
            self.tree.phi_to_cell_mapping(phi)
            cost = self.tree.compute_pooled_costs(self.data, self.lamb)
            if cost < min_fitness:
                min_fitness = cost 
                best_assign = phi
            self.fitness.append(cost)
        return min_fitness, best_assign


    @staticmethod 
    def get_indices_of_smallest_n(lst, n):
        # Enumerate the list to get (index, value) pairs
        indexed_list = list(enumerate(lst))
        
        # Sort the list of (index, value) pairs based on values
        sorted_list = sorted(indexed_list, key=lambda x: x[1])
        
        # Get the first n elements from the sorted list
        smallest_n = sorted_list[:n]
        
        # Extract the indices from the (index, value) pairs
        indices = [index for index, value in smallest_n]
        
        return indices    


    def randomly_pair_items(self, my_list):
        # Make a copy of the list to avoid modifying the original list
     
        
        # Shuffle the list randomly
        # self.rng.shuffle(my_list)
        
        # Iterate through the shuffled list two elements at a time
        paired_items = [(my_list[i], my_list[i + 1]) for i in range(0, len(my_list), 2)]
        
        return paired_items

    def crossover(self):
        new_pop = []
        most_fit = self.get_indices_of_smallest_n(self.fitness, self.parents)
        
        new_parents = self.randomly_pair_items(most_fit)
        for x, y in new_parents:
            new_phi = {}
            for i in self.data.cells:
                if self.rng.random() > 0.5:
                    new_phi[i] = self.populations[x][i]
                else:
                    new_phi[i] = self.populations[y][i]
            if self.rng.random() < 0.05:
                new_phi[i] = self.rng.choice(self.clones)
                
            new_pop.append(new_phi)
        self.populations= new_pop

        
   
    def fit(self, tree, data, lamb):
        self.tree = tree 
        self.clones =self.tree.clones()
        self.data = data 
        self.lamb = lamb
        self.tree.assign_cells(self.data)
        cost = self.tree.compute_pooled_costs(self.data, self.lamb)
        self.tree.comp_df.to_csv("test/pooled_vafs_greedy1.csv", index=False)

        print(f"Greedy cost: {cost}")
        for v, score in self.tree.node_cost.items():
            print(f"{v}: {score}")
        self.base_phi = self.tree.get_phi()
        self.initialize()
        cur_fitness = np.Inf 
        for j in range(self.iter):
            min_fitness, best_assign = self.eval_fitness()
            self.overall_fitness.append(min_fitness)
            self.overall_assign.append(best_assign)
            print(f"{j}: min fitness: {min_fitness} current fitness: {cur_fitness} gap: {np.abs(cur_fitness-min_fitness)}")
            # if np.abs(min_fitness - cur_fitness) <= self.tol:
        
            # else:
            cur_fitness = min_fitness
            
            self.crossover()
            self.add_new_parents()
        
        return self.overall_fitness, self.overall_assign

                
                

class HeuristPooled:
    def __init__(self, iter=100, tol=0.01, mut_rate=0.025, restarts=100, seed=1026) -> None:
        self.iter = iter
        self.tol = tol
        self.mut_rate = mut_rate
        self.restarts = restarts
        self.rng = np.random.default_rng(seed)
    
    def create_individual(self, base_phi):
        # cell_assign = self.rng.choice(self.clones, self.data.N, replace=True)
        # phi = {i: cell_assign[i] for i in self.data.cells}
        # return phi

        phi = {}
        for i, assign in base_phi.items():
            if self.rng.random() < self.mut_rate:
                phi[i]  = self.rng.choice(self.clones)
            else:
                phi[i] = assign
        return phi
    
    def initialize(self,phi):
        self.populations = []
    
        self.populations.append(phi)

        for p in range(self.restarts):
            phi = self.create_individual(phi)
        
            self.populations.append(phi)
    
    def compute_centroids(self):
    
        self.centroids = {}
        for v in self.tree.all_dfs:
               
               df = self.tree.all_dfs[v]
               df = df.set_index(df['snvs'])
               self.centroids[v] = df['obs_vafs']




    def reassign_cells(self, phi):
        df =pd.concat([ser for key,ser in self.centroids.items()], axis=1)
        cols = [v for v in self.centroids]
        centroids =df.values
        new_phi= {}
        for i in phi:
            vaf =self.data.compute_vafs([i], df.index.to_numpy())
            vaf = vaf.reshape(-1,1)
            total_dev = np.nansum(np.abs(vaf - centroids), axis=0)
            new_phi[i]=cols[total_dev.argmin()]
        return new_phi

    def fit(self, tree, data, lamb):
        self.tree = tree 
        ev = Eval(self.tree)
        self.clones =self.tree.clones()
        self.data = data 
        self.lamb = lamb
        self.tree.assign_cells(self.data)
        init_greedy= self.tree.compute_costs(self.data, self.lamb)
        phi = self.tree.get_phi()
        self.initialize(phi)
        all_res = []
        best_cost = np.Inf 
        best_acc = 0
        for k, phi in enumerate(self.populations):
            self.tree.phi_to_cell_mapping(phi)
            cost = self.tree.compute_pooled_costs(self.data, self.lamb)


            init_acc = ev.accuracy(phi)

            acc = []
            for j in range(self.iter):
                self.compute_centroids()
                phi = self.reassign_cells(phi)
                self.tree.phi_to_cell_mapping(phi)
                new_cost =self.tree.compute_pooled_costs(self.data,self.lamb)
                greedy_cost = self.tree.compute_costs(self.data, lamb)
                if new_cost < best_cost:
                    best_cost = new_cost
                    best_cost_phi = phi 
                
                acc.append(ev.accuracy(phi))
                if acc[j] > best_acc:
                    best_acc = acc[j]
                    best_acc_phi = phi
                # print(f"restart: {k}: iter:{j}: accuracy: {acc[j]} new cost: {new_cost} prev cost: {cost}")
                all_res.append([k, j, acc[j], new_cost, init_acc, greedy_cost])
                if np.abs(cost- new_cost) <= self.tol:
                    break 
                else:
                    cost = new_cost 
        df = pd.DataFrame(all_res, columns=["restart", "iteration", "acc", "cost", "init_acc", "greedy"])
        df["init_greedy"] =init_greedy
        return best_cost, best_cost_phi, best_acc, best_acc_phi, df








def main():
    pth = "/scratch/data/leah/pharming/simulation_study/input/s14_m10000_k25_l7/n1000_c0.1"
    tree = load_from_pickle(f"{pth}/gt.pickle")
    dat = load_from_pickle(f"{pth}/data.pickle")
    lamb = 0
    gt_cost = tree.compute_pooled_costs(dat, lamb)
    gt_greedy = tree.compute_costs(dat, lamb)
    tree.comp_df.to_csv("test/pooled_vafs_gt1.csv", index=False)
    for v, score in tree.node_cost.items():
        print(f"{v}: {score}")
    print(f"Ground truth pooled cost: {gt_cost}")
    hp = HeuristPooled(restarts=25)
    cost, cost_phi, acc, acc_phi, df = hp.fit(tree, dat, lamb)
    print(f"K-means cost: {cost}")
    df["gt_cost"] = gt_cost
    df["gt_greedy"] = gt_greedy
    df.to_csv("test/all_scores_0.25.csv", index=False)
    # cae = CellAssignEvol()

    # all_fit, all_assign = cae.fit(tree,dat,lamb)


if __name__ == '__main__':
    main()


#kmeans pooled assignment