
from dataclasses import dataclass
from clonal_tree import ClonalTree
from cell_mapping import CellAssign


@dataclass
class Solution:
    cost: float
    ct: ClonalTree
    phi: CellAssign

    def __post_init__(self):
        self.segments = self.ct.get_segments()
        self.m = len(self.ct.get_all_muts())

    def png(self, fname, segments=None):
        if segments is None:
            segments = self.segments

        self.ct.draw(fname,self.phi, segments=segments, include_dcfs=True)
    
    def toJson(self, fname, segment_csv=None, snv_csv=None):
        self.ct.toJson(fname, segment_csv, snv_csv)

    def get_tree(self):
        return self.ct 
    
    def optimize(self, data, lamb):
        self.cost, self.phi = self.ct.optimize(data, lamb)
    
    def collapse(self, k, cell_threshold=10):
        self.ct.collapse(self.phi, k, cell_threshold)
    
