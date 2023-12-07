from dataclasses import dataclass

@dataclass
class CellAssign:
    phi: dict
    clones: set 


    def __post_init__(self):

        self.cell_mapping = self.to_mapping()

    def update(self, phi):
        self.phi = phi
        self.cell_mapping = self.to_mapping()

    
    def update_phi(self):
        self.phi ={i : v  for v in self.clones for i in self.cell_mapping[v] if v in self.cell_mapping}

    def to_mapping(self):

        cell_mapping = {v: [] for v in self.clones}
        for i, v in self.phi.items():
            cell_mapping[v].append(i)
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
    
    def get_all_cells(self):
        all_cells = []
        for n in self.cell_mapping:
            all_cells += self.cell_mapping[n]
        return all_cells
    
    def get_cell_count(self):
         return {n : len(self.cell_mapping[n]) for n in in self.clones if n self.cell_mapping}
    
    # def get_cell_assigment(self, i):
    #     return 
    #     for node in self.cell_mapping:
    #         if i in self.cell_mapping[node]:
    #             return node
         