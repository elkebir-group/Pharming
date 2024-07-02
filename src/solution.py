
from dataclasses import dataclass
from clonal_tree import ClonalTree
from cell_mapping import CellAssign
from draw_pretty_tree import DrawPrettyTree
from cmb import cmb_all

@dataclass
class Solution:
    cost: float
    ct: ClonalTree
    phi: CellAssign

    def __post_init__(self):
        self.update_segments()
        self.m = len(self.ct.get_all_muts())

    def png(self, fname, segments=None):
        if segments is None:
            segments = self.segments

        self.ct.draw(fname,self.phi, segments=segments, include_dcfs=True)
    
    def toJson(self, fname, segment_csv=None, snv_csv=None):
        self.ct.toJson(fname, segment_csv, snv_csv)

    def get_tree(self):
        return self.ct 
    
    def optimize(self, data, lamb, cell_threshold=0):
        # if self.segments == set([145, 178, 237]):
        #     self.png("dlp/test/pre.png")
        self.cost, self.phi = self.ct.optimize(data, lamb)
    
    def prune_leaves(self, k):
        self.ct.prune_leaves(self.phi, k)

    def post_process(self, dat, lamb, k, cell_threshold=0):
        self.optimize(dat, lamb)
        if cell_threshold > 0:
            self.collapse(k, cell_threshold)
        # self.prune_leaves()
        self.compute_likelihood(dat, lamb)


    def collapse(self, k, cell_threshold=10):
        self.ct.collapse(self.phi, k, cell_threshold)
    
    def compute_likelihood(self, data, lamb):
        cost= self.ct.compute_likelihood(data, self.phi, lamb)
        snv = self.ct.snv_cost
        cna =  self.ct.cna_cost
        return cost, snv, cna
    
    def update_segments(self):
        self.segments = self.ct.get_segments()
    
    def write_flat_files(self, pth, data=None, lamb=1000):
        self.phi.write_phi(f"{pth}/pred_cell.csv")
        self.ct.write_psi(f"{pth}/pred_mut.csv")
        self.ct.write_loss_mapping(f"{pth}/pred_mut_loss.csv")
        self.ct.write_genotypes(f"{pth}/pred_genos.csv")
        self.ct.save_text(f"{pth}/pred_tree.txt")
        if data is not None:
            cost, snv, cna = self.compute_likelihood(data, lamb)
            with open(f"{pth}/likelihood.csv", "w+") as file:
                file.write("cost,snv,cna\n")
                file.write(f"{cost},{snv},{cna}\n")

    def ICL(self, data, lamb):
       icl =  self.ct.ICL(self.phi, data, lamb)
       return icl, self.ct.bic #, #self.ct.entropy

    def drawPrettyTree(self, tname, lname=None, include_CNAs=True):
        dtp = DrawPrettyTree(self, include_CNAs=include_CNAs)
        dtp.draw(tname)
        if lname is not None:
            dtp.save_labels(lname)


    def computeCMB(self, dat, mincells=10, fname=None):
        cmb_df = cmb_all(self.ct, self.phi, dat, min_cells=mincells)
        if fname is not None:
             cmb_df.to_csv(fname, index=False)
     
        return cmb_df