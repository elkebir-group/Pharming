from dataclasses import dataclass
import pickle 
import pandas as pd 
import numpy as np
import argparse 
from sklearn.metrics.cluster import adjusted_rand_score

@dataclass
class CellAssign:
    phi: dict
    clones: set 


    def __post_init__(self):

        self.cell_mapping = self.to_mapping()
        self.n = len(self.phi)

    def update(self, phi):
        self.phi = phi
        self.cell_mapping = self.to_mapping()

    def update_clones(self, clones):
        self.clones= clones
    
    def update_phi(self):
        self.phi ={i : v  for v in self.clones for i in self.cell_mapping[v] if v in self.cell_mapping}

    def move_cells(self, u,v):
        if u in self.clones and v in self.clones:
            cells = self.cell_mapping[u]
            for i in cells:
                self.phi[i] = v
            del self.cell_mapping[u] 
            self.cell_mapping[v] = np.union1d(self.cell_mapping[v], cells)
            self.cell_mapping[v] = self.cell_mapping[v].astype(int)
        else:
            print("Warning: node does not exist, cells not moved!")
 

    def move(self, cell, node):
        if node in self.clones:
            cur_node = self.phi[cell]
            self.phi[cell] = node 
            arr = self.cell_mapping[cur_node]
            self.cell_mapping[cur_node] =arr[arr != cell]
            self.cell_mapping[node] = np.append(self.cell_mapping[node], cell)
        else:
            print("Warning: node does not exist. Cell not moved")
            
    def to_mapping(self):

        cell_mapping = {v: [] for v in self.clones}
        for i, v in self.phi.items():
            cell_mapping[v].append(i)
        cell_mapping = {v: np.array(val) for v, val in cell_mapping.items()}
        return cell_mapping
    
    def from_mapping(self, cell_mapping):
          self.cell_mapping = cell_mapping

          self.update_phi()
          return self.phi
    
    def get_cells(self,node):
        if node not in self.cell_mapping:
            return []
        else:
            return self.cell_mapping[node]
    
    # def get_all_cells(self):
    #     all_cells = []
    #     for n in self.cell_mapping:
    #         all_cells += self.cell_mapping[n]
    #     return set(all_cells)
    
    def get_cell_count(self):
         return {n : self.cell_mapping[n].shape[0] for n in self.clones if n in self.cell_mapping}
        

    def save(self, path):
        with open(path, "wb") as file:
              pickle.dump(self, file)

    def relabel(self, cell_lookup):
        new_phi = {}
        for index, cell_label in cell_lookup.items():
            new_phi[index]  = self.phi[cell_label]
        self.update(new_phi)

    
    def relabel_clones(self, mapping):
        phi = {}
        for i, u in self.phi.items():
            if u in mapping:
                phi[i] = mapping[u]
        self.update(phi)

    def to_dataframe(self):
         df = pd.DataFrame(list(self.phi.items()), columns=['cell', 'cluster'])
         return df
    
    def write_phi(self, fname):

        df = self.to_dataframe()
        df.to_csv(fname, index=False)

    def compute_ari(self,obj) -> float:
        #  gt_mut = self.get_mut_clusters()
         gt =   pd.Series(self.phi)
         pred = pd.Series(obj.phi)
         df = pd.concat([gt, pred], axis=1, keys=['gt', 'pred'])
        #  pred_mut = obj.get_mut_clusters()
 
         return adjusted_rand_score(df["gt"].values, df["pred"].values)


    def __str__(self):
        cell_count = self.get_cell_count()
        mystr = ""
        for n, count in cell_count.items():
            mystr += f"{n}: {count}\n"
        return mystr

# def load_from_pickle(fname):
#     return pd.read_pickle(fname)



# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-p", "--phi", required=True,
#                         help="input file for cell assignment")
#     parser.add_argument("-t", "--tree", required=True,
#                         help="pickled clonal tree")
#     parser.add_argument("-P", "--assign",
#                         help="input file for cell assignment")
        
    
#     args = parser.parse_args()

#     ct = load_from_pickle(args.tree)
    
#     cell_assignment = pd.read_csv(args.phi)
#     phi = {}

#     # print(cell_assignment.head())
#     for index, row in cell_assignment.iterrows():

#         i = row['Cell']
#         v = row['Cluster']

        
    
#         phi[index] = v


#     ca = CellAssign(phi, ct.clones())
#     if args.assign is not None:
#         ca.save(args.assign)

         