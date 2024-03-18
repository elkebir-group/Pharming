
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

    def png(self, fname, segments=None):
        if segments is None:
            segments = self.segments

        self.ct.draw(fname,self.phi, segments=segments, include_dcfs=True)
    
    def toJson(self, fname, segment_csv=None, snv_csv=None):
        self.ct.toJson(fname, segment_csv, snv_csv)

    def get_tree(self):
        return self.ct 
    
    def prune(self):
        self.ct.prune(self.phi)


