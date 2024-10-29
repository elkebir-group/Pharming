
from dataclasses import dataclass
from clonal_tree import ClonalTree
from cell_mapping import CellAssign
from draw_pretty_tree import DrawPrettyTree
from cmb import cmb_all, cmb_loss_all, vaf_validation
import pandas as pd 
import numpy as np
from scipy.special import logsumexp
import networkx as nx 

@dataclass
class Solution:
    """
    Solution to the CTI problem.

    Attributes
    ----------
    cost : float
        the cost of the solution
    ct : ClonalTree
        the inferred clonal tree
    phi : CellAssign
        the inferred assignment of cells to clones in the clonal tree
    """
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
    
    def write_flat_files(self, pth, data=None, lamb=1000, stub="pred"):
        self.phi.write_phi(f"{pth}/{stub}_cell.csv")
        self.ct.write_psi(f"{pth}/{stub}_mut.csv")
        self.ct.write_loss_mapping(f"{pth}/{stub}_mut_loss.csv")
        self.ct.write_genotypes(f"{pth}/{stub}_genos.csv")
        self.ct.save_text(f"{pth}/{stub}_tree.txt")
        if data is not None:
            cost, snv, cna = self.compute_likelihood(data, lamb)
            with open(f"{pth}/{stub}_likelihood.csv", "w+") as file:
                file.write("cost,snv,cna\n")
                file.write(f"{cost},{snv},{cna}\n")

    def ICL(self, data, lamb):
       icl =  self.ct.ICL(self.phi, data, lamb)
       return icl, self.ct.bic #, #self.ct.entropy
    
    def model_selection(self, data, lamb, fname=None):
    
        bic = self.ct.BIC(self.phi, data, lamb)

        if fname is not None:
            with open(fname, "w+") as file:
                file.write("solution,SNV_likelihood,BIC,cost\n")
                file.write(f"0,{self.ct.snv_cost},{bic},{self.cost}\n")
           

    def drawPrettyTree(self, tname, lname=None, include_CNAs=True, 
                       min_cells=1, segment_dict=None, include_cell_counts=True,
                       horizontal=False):
        dtp = DrawPrettyTree(self, include_CNAs=include_CNAs, min_cells=1,
                              segment_dict=segment_dict, include_cell_counts=include_cell_counts)
        dtp.draw(tname, horizontal=horizontal)
        if lname is not None:
            dtp.save_labels(lname)


    def computeCMB(self, dat, mincells=10, fname=None):
        cmb_df = cmb_all(self.ct, self.phi, dat, min_cells=mincells)
        if fname is not None:
             cmb_df.to_csv(fname, index=False)
     
        return cmb_df
    
    def computelossCMB(self, dat, mincells=10, fname=None):
        cmb_df = cmb_loss_all(self.ct, self.phi, dat, min_cells=mincells)
        if fname is not None:
             cmb_df.to_csv(fname, index=False)
     
        return cmb_df
    
    def computeVAFs(self, dat, mincells=10, fname=None):
        vaf = vaf_validation(self.ct, self.phi, dat, min_cells=mincells)
        if fname is not None:
             vaf.to_csv(fname, index=False)
     
        return vaf
    
    def get_lost_snvs(self):
        return self.ct.get_lost_snvs()
    
    def get_snv_cluster_tree(self):
        return self.ct.mutation_cluster_tree()

    def assess_loss(self, snvs:list ):
        prec = self.ct.loss_precision(snvs)
        recall = self.ct.loss_recall(snvs)
        fps_leaves = self.ct.fp_leaves(snvs)
        return prec, recall, fps_leaves

    def cna_loss(self, snvs:list):
        tp, tn, fn, fp = self.ct.cna_loss(snvs)
        prec = tp / (tp + fp)
        recall = tp / (tp + fn)
        acc =   (tp + tn) / (tp + tn + fn + fp)

        return prec, recall, acc
    
    def cna_metrics(self, gt, fname=None):
        res=  self.ct.cna_metrics(gt)
        if fname is not None:
            
            df = pd.DataFrame([res])
            df.to_csv(fname, index=False)
        else:
            return res

    def snv_assignment_posterior(self, data, fname=None):
        df = self.all_snv_tree_costs(data)

        df_sorted = df.sort_values(['snv', 'cluster', 'cost', 'tree_id'])
        filtered_df = df_sorted.groupby(['snv', 'cluster'], as_index=False).first()
        filtered_df = filtered_df[['snv', 'cluster', 'cost', 'tree_id']]
        filtered_df = pd.merge( df,filtered_df, on=['snv', 'cluster', 'cost', 'tree_id'])

        filtered_df['cost'] = -1 * filtered_df['cost']


     

        # Step 3: Compute the log denominator using logsumexp for each group of 'folder' and 'snvs'
        log_denom_df = filtered_df.groupby(['snv'], as_index=False).agg(
            log_denom=('cost', lambda x: logsumexp(x))
        )

        # Step 4: Merge log_denom back into snv_tree_min DataFrame
        filtered_df = pd.merge(filtered_df, log_denom_df, on=['snv'])
  

        # Step 5: Calculate log posterior and posterior probabilities
        filtered_df['log_posterior'] = filtered_df['cost'] - filtered_df['log_denom']
        filtered_df['posterior'] = np.exp(filtered_df['log_posterior'])

        # Step 6: Avoid zero posterior by using the smallest float value if posterior is zero
        filtered_df['posterior'] =filtered_df['posterior'].replace(0, np.finfo(float).eps)

        # Step 7: Filter where 'cluster' equals 'assignment'
        snv_tree_post =filtered_df[filtered_df['cluster'] ==filtered_df['assigment']]
        snv_tree_post = snv_tree_post[[ 'snv','cluster', 'posterior']]

        if fname is not None:
            snv_tree_post.to_csv(fname, index=False)
        return snv_tree_post



    def all_snv_tree_costs(self, data):
        return self.ct.compute_genotype_costs(data, self.phi)
    
    def compute_snv_likelihoods(self, data):
        return self.ct.compute_snv_likelihoods(data, self.phi)

    # def refit_segment(self, segment, data,lamb, cna_tree=None, threshold=0.05):
    #     self.ct.filter_segments([segment])
    #     Tm = self.get_snv_cluster_tree()
    #     if not cna_tree:
    #         cn_prop = self.data.thresholded_cn_prop(segment, threshold,
    #                                                        start_state= (1,1), include_start_state=False)
    #         if len(cn_prop) == 1:


    def is_ancestral(self, i,j):
        
        u = self.phi.phi[i]
        v = self.phi.phi[j]
        if u==v: return True
        return u in nx.ancestors(self.ct.tree, v)

    def is_clustered(self, i,j):
        u = self.phi.phi[i]
        v = self.phi.phi[j]
        return u == v

    def is_incomparable(self, i,j):
        u = self.phi.phi[i]
        v = self.phi.phi[j]
        return self.ct.is_incomparable(u,v)

    