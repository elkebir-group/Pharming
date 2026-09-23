
# Created by: L.L. Weber
# Created on: 2024-02-29 18:40:59
import itertools
import networkx as nx 
import numpy as np
import multiprocessing
import sys
import time
from .sti_v2 import STI
from . import clonelib
from .tree_merging import ClonalTreeMerging
from .utils import concat_and_sort, timeit_decorator, get_top_n
from .dcf_clustering_v2 import DCF_Clustering


class Pharming:
    def __init__(self, 
                dcfs=None, 
                k=3,
                start_state=(1,1), 
                seed=102,
                verbose=False, 
                top_n=3, 
                collapse=False,
                order = "weighted-random",
                ninit_segs = None,
                ninit_Tm = None,
                cell_threshold=5,
                thresh_prop = 0,
                sum_condition = True
                  ) -> None:
        
        
        self.verbose =verbose
        self.verbose = True 
        self.rng = np.random.default_rng(seed)

        if dcfs is not None:
            self.k = len(dcfs)
            self.delta  = dcfs 
        else:
            self.k = k 
            self.delta = None
        
 

      
        self.cnatrees = {} 


        
        self.start_state = start_state
        print(f"Start state: {self.start_state}")
        self.ninit_segs= ninit_segs
        print(f"Max # of initial segments: {self.ninit_segs}")
        self.ninit_Tm = ninit_Tm
        print(f"Max # of full inference mutation cluster trees {self.ninit_Tm}")
        self.top_n = top_n
        print(f"Top n: {self.top_n}")

        self.collapse = collapse
        if cell_threshold is not None:
            self.cell_threshold = cell_threshold
        else:
            cell_threshold = 0

        print(f"Collapse CNA internal nodes: {self.collapse} with cell threshold: {self.cell_threshold}")
        self.order = order
        print(f"Integration ordering: {self.order}")
      
        self.thresh_prop = thresh_prop
        print(f"CN states threshold proportion of cells per segment: {self.thresh_prop}")
        self.sum_condition = sum_condition
        print(f"Filtering mutation cluster trees using the sum condition: {self.sum_condition}")

                
     
    @staticmethod
    def check_dcfs(T, delta):
   
        for u in T:
            desc_dcf = sum( delta[u] for u in  sorted(T.successors(u)))
            if delta[u] < desc_dcf or desc_dcf > 1:
                return False 
        return True 

    def enumerate_mutcluster_trees(self, delta):
        '''
        Enumerate the set of mutation cluster trees that respect 
        the sum condition of the DCFs
    
        '''

        G = nx.DiGraph()
        G.add_nodes_from([q for q in range(self.k)])
        for u,v in itertools.combinations(range(self.k),2):
            if delta[u] >= delta[v]:
                G.add_edge(u,v, weight=1)
            elif delta[u] < delta[v]:
                G.add_edge(v,u, weight=1)
  
            
        trees = nx.algorithms.tree.branchings.ArborescenceIterator(G)

        if self.sum_condition:
            return [tree for tree in trees if self.check_dcfs(tree,delta)]
        else:
            return [tree for tree in trees]
        
   
    def enumerate_cna_trees_python(self, cn_states):
        cn_states = [cn for cn in cn_states if cn != self.start_state]

        G = nx.DiGraph()
        G.add_nodes_from(cn_states)
        for u in G.nodes:
            for v in G.nodes:
                if u == v:
                    continue
                if (u[0] == 0 and not v[0] > 0) or (u[1] == 0 and not v[1] > 0):
                    continue
                G.add_edge(u,v, weight=1)

        G.add_node(self.start_state)
        for u in cn_states:
            G.add_edge(self.start_state, u,weight=1)

        cnatrees = nx.algorithms.tree.branchings.ArborescenceIterator(G)

        #check the CNA  tree to make sure we don't get a ressurrection of an allele
        def check_tree(T):
            for u,v in T.edges:
                if u[0] ==0 and v[0] > 0:
                    return False 
                if u[1] ==0 and v[1] > 1:
                    return False
            return True

        scriptS = [S for S in cnatrees if check_tree(S)]


     
        return scriptS
       


    def enumerate_cna_trees(self, cn_states):
   
        trees = clonelib.get_cna_trees(cn_states, *self.start_state )
    
        def convert_to_CNAtree(tree):
    
            S = nx.DiGraph()
            if len(tree) == 0:

                S.add_node(self.start_state)
            else:
                S.add_edges_from(tree)
       
            return S

        T_CNAS = [convert_to_CNAtree(tree) for tree in trees]
        

        return T_CNAS


   

    

        
    def fit_segment(self, ell,T_m, stis):
        """
        Fits a segment tree for the data in segment ell given mutation cluster tree T_m.
        
        This function takes a segment ID `ell` and a mutation cluster tree `T_m`, and 
        fits a segment tree for the data in the specified segment. 

        :param int ell: Segment ID.
        :param nx.DiGraph T_m: Mutation cluster tree with nodes labeled by cluster q in [k].
        :return: A list of the top n solutions.
        :rtype: list of Solutions
        """
        if self.verbose:
            print(f"Starting segment {ell}...")
        
        segtrees = []
        Tm_edges = list(T_m.edges)
        for st in stis:
            trees = st.fit(Tm_edges, self.data, ell)
            if len(trees) > 0:
                segtrees.append(trees)
            
        if len(segtrees) ==0:
            print(f"Segment {ell} failed for {Tm_edges}!")
        else:
            segtrees = concat_and_sort(segtrees)

            print(f" segment | cost | snv | cna")
        
            cost,snv, cna = segtrees[0].compute_likelihood(self.data, self.lamb)
            print(f"|{ell} | {cost} | {snv} | {cna} |")

        if self.verbose:
            print(f"Segment {ell} complete!")

        # best_tree = segtrees[0]
        # best_tree.png(f"test/s11/inf{ell}.png")
        return segtrees


    def segment_trees_inference(self, Tm, stis):
        """
        nx:DiGraph Tm: a mutation cluster tree
        iterable segments: an iterable of segments for which segtrees should be inferred, given Tm

        return list of lists of segment trees, one list per each segment
        """


        if self.cores  <= 1:
        # if True:
            segtrees = []
            for ell in stis:
                # try:
                    if self.data.num_cn_states(ell,self.thresh_prop) > 1:
                
                        segtrees.append(self.fit_segment(ell, Tm, stis[ell]))
                        
        else:
            arguments = [(ell, Tm, stis[ell]) for ell in stis 
                            if self.data.num_cn_states(ell, self.thresh_prop) > 1]
            with multiprocessing.Pool(processes=self.cores) as pool:                
                segtrees = pool.starmap(self.fit_segment, arguments, chunksize=1)
        
 
        return segtrees
   
    @timeit_decorator
    def integrate(self, Tm_edges, segtrees, restarts=1):
            
            print(f"Starting integration for {len(segtrees)} segments...")
            # ctm = ClonalTreeMerging(self.k, self.rng, top_n=self.top_n, order=self.order,
            #                         collapse=self.collapse, cell_threshold=self.cell_threshold)
            
            all_trees = []
            for i in range(restarts):
                Tm = nx.DiGraph(Tm_edges)

                #add the normal clone as the root
                root = [n for n in Tm if Tm.in_degree[n]==0][0]
                Tm.add_edge(max(Tm)+2, root) 

                
                ctm = ClonalTreeMerging(self.k, self.rng, top_n=self.top_n, order='in-place',
                                        collapse=self.collapse, cell_threshold=self.cell_threshold)
                
                top_trees = ctm.fit(segtrees, Tm, self.data, self.lamb, cores=self.cores)
                all_trees.extend(top_trees)
            
            return top_trees

    def order_segments(self, segs:set):
            
            segs = [ell for ell in segs]
            weights = [self.data.num_snvs(ell) for ell in segs]
            weights = weights/np.sum(weights)
           
            ordering = self.rng.choice(len(segs), size=len(segs), replace=False, p=weights)   
            seg_order = [segs[i] for i in ordering]
            return  seg_order

    def partition_segments(self, segments,  min_cn_states=2):
        if self.ninit_segs is None or self.ninit_segs > len(segments):
            init_segs = [ell for ell in segments 
                         if self.data.num_cn_states(ell, self.thresh_prop, include_start_state=False) >= min_cn_states and self.data.num_snvs(ell) > 0]

        else:
      
            init_segs = sorted([ell for ell in segments if 
                                self.data.num_cn_states(ell, self.thresh_prop, include_start_state=False) >= min_cn_states],
                                reverse=True, key= lambda x: self.data.num_snvs(x))
            if len(init_segs) > self.ninit_segs:
                init_segs = init_segs[:self.ninit_segs]

        
        init_segs = set(init_segs)
        

        
        remaining_segs = set(segments) -  init_segs
        infer_segs = set([ell for ell in remaining_segs if self.data.num_cn_states(ell, self.thresh_prop, include_start_state=False) > 1 and self.data.num_snvs(ell)> 0])
        place_segs = set([ell for ell in segments if self.data.num_cn_states(ell, self.thresh_prop, include_start_state=False)==1 and self.data.num_snvs(ell)> 0])

        no_snvs_segs = set([ell for ell in segments if  self.data.num_snvs(ell) ==0])

        print("Init Segs:")
        print(init_segs)

        print("Infer Segs:")
        print(infer_segs)

        print("Place Segs:")
        print(place_segs)

        print("No Snvs Segs:")
        print(no_snvs_segs)
        return init_segs, infer_segs, place_segs, no_snvs_segs
      

    def infer(self, Tm_list, stis, init_trees=None, init_order=None):
        """
        infer a clonal tree for a list of mutation cluster trees
        and a subset of segments.  

        list Tm_list: a list of networkx DiGraphs of mutation cluster trees
        iterable seg_list: a set/list of segments to infer clonal trees for each segment
        list init_trees: a list of lists of clonal trees previoyly inferred on disjoint segments from seg_list
        """
        costs = []
        all_trees = []

        #if given a list of initial trees, then initialize the integrated tree with the first list of tres in the list
        init_on_first = init_trees is not None

        if init_order is not None:
            order_dict = {ell: i for i,ell in enumerate(init_order)}

        for i,Tm in enumerate(Tm_list):
            print(f"Staring Tm {i}: {list(Tm.edges)}")
            # if set([(2, 0), (2, 1), (3, 2)]) != set(Tm.edges):
            #     continue

            segtrees = self.segment_trees_inference(Tm, stis)
      
            print("All specified segment tree inference complete!")
            if init_order is not None:
                segtrees = sorted(segtrees, key=lambda tl: order_dict[list(tl[0].ct.get_segments())[0]])
        
            if init_trees is None:
                tree_list = segtrees
            else:
                tree_list =  [init_trees[i]] + segtrees
            print("Starting integration...")
            top_trees = self.integrate(Tm, tree_list)
            print("Integration complete...")
            if len(top_trees) > 0:
                best_cost = top_trees[0].cost
                costs.append(best_cost)
                all_trees.append(top_trees)
            else:
                print(f"Integration failed for Tm_{i}")
                all_trees.append([])
                costs.append(np.Inf)
            print(f"Tm {i} complete!")

        return all_trees, costs 

    def place_snvs(self, solutions, segments):
        """
        Place SNVs that occur in segments with only a single copy number state
        in the clonal tree of each solution.

        list solutions: a list of solutions
        iterable segments: an iterable of segments that consist of only 1 copy number state
        """
        rho = {ell: {} for ell in segments}
        seg_to_snvs = {}
        states = {}
        for ell in segments:
     
            rho[ell] = {}
            cn_prop = self.data.thresholded_cn_prop(ell, self.thresh_prop,
                                                          self.start_state, include_start_state=False)
            assert len(cn_prop) == 1
            for k in cn_prop:
                state= k 
                
            states[ell] = state
            T_SNV = nx.DiGraph()
     

            T_SNV.add_edge((*state, 0,0), (*state, 1,0))
            if state != self.start_state:
                T_SNV.add_edge((*self.start_state, 0,0), (*state, 0,0))
       
            rho[ell] = {q: [T_SNV] for q in range(self.k)}
            seg_to_snvs[ell] = self.data.seg_to_snvs[ell]

        for sol in solutions:
            sol.ct.assign_genotypes(self.data, sol.phi, rho, seg_to_snvs, states, start_state=self.start_state)
            sol.ct.add_rho(rho)
            sol.update_segments()

    def place_cnas(self, solutions, segments):


        
        for ell in segments:
           states, counts = self.data.cn_states_by_seg(ell)
           setS = self.enumerate_cna_trees(states)



    def preprocess_helper(self, ell, delta):
        print(f"Segment {ell}: starting preprocessing..")
  

        stis = []
        cn_prop = self.data.thresholded_cn_prop(ell, self.thresh_prop, self.start_state)
        cn_states = set(cn_prop.keys())

        # try:
        
        if len(self.cnatrees) == 0 or ell not in self.cnatrees:
            cnatrees =  self.enumerate_cna_trees_python(cn_states)
            #  cnatrees =  self.enumerate_cna_trees(cn_states)
        else:
            cnatrees = [self.cnatrees[ell]]
        for S in cnatrees:
            st = STI(ell, S, delta, lamb=self.lamb, prop_thresh=self.thresh_prop)
            st.precompute_costs(self.data)
            stis.append(st)
    

        
        print(f"Segment {ell}: preprocessing complete")
        # except Exception as e:
            # print(f"Segment {ell} timed out enumerating CNA trees, skipping..")

 
        return ell, stis


    
    """Paste these methods inside the Pharming class in pharming.py."""

    def _export_dcf_bootstrap_analysis(self,matrix_norm_boot,normal_sorted_dcfs,boot_sorted_dcfs,output_dir="bootstrap_results",):
        """Save DCF tables, co-clustering heatmaps, and raw analysis arrays."""
        from pathlib import Path

        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import squareform

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        matrix_norm_boot = np.asarray(matrix_norm_boot)
        normal_sorted_dcfs = np.asarray(normal_sorted_dcfs, dtype=float)
        boot_sorted_dcfs = np.asarray(boot_sorted_dcfs, dtype=float)

        if matrix_norm_boot.ndim != 4 or matrix_norm_boot.shape[0] != 2:
            raise ValueError(
                "matrix_norm_boot must have shape (2, snv_count, snv_count, 2)."
            )
        if matrix_norm_boot.shape[-1] != 2:
            raise ValueError(
                "The last matrix dimension must contain [different, same] counts."
            )
        if normal_sorted_dcfs.ndim != 2 or boot_sorted_dcfs.ndim != 2:
            raise ValueError(
                "The sorted DCF arrays must have shape (replicates, clusters)."
            )
        if normal_sorted_dcfs.shape[1] != boot_sorted_dcfs.shape[1]:
            raise ValueError(
                "Normal and bootstrap runs must have the same number of clusters."
            )

        snv_count = matrix_norm_boot.shape[1]
        number_of_clusters = normal_sorted_dcfs.shape[1]
        dcf_columns = [
            f"dcf_rank_{rank}_descending"
            for rank in range(1, number_of_clusters + 1)
        ]

        normal_dcf_table = pd.DataFrame(
            normal_sorted_dcfs,
            columns=dcf_columns,
        )
        normal_dcf_table.insert(
            0,
            "replicate",
            np.arange(1, len(normal_dcf_table) + 1),
        )
        normal_dcf_table.to_csv(
            output_dir / "normal_sorted_dcfs.csv",
            index=False,
        )

        boot_dcf_table = pd.DataFrame(
            boot_sorted_dcfs,
            columns=dcf_columns,
        )
        boot_dcf_table.insert(
            0,
            "replicate",
            np.arange(1, len(boot_dcf_table) + 1),
        )
        boot_dcf_table.to_csv(
            output_dir / "bootstrap_sorted_dcfs.csv",
            index=False,
        )

        normal_counts = matrix_norm_boot[0]
        boot_counts = matrix_norm_boot[1]
        valid_normal = normal_counts.sum(axis=2)
        valid_boot = boot_counts.sum(axis=2)

        normal_rates = np.divide(
            normal_counts[:, :, 1],
            valid_normal,
            out=np.full((snv_count, snv_count), np.nan),
            where=valid_normal > 0,
        )
        boot_rates = np.divide(
            boot_counts[:, :, 1],
            valid_boot,
            out=np.full((snv_count, snv_count), np.nan),
            where=valid_boot > 0,
        )

        np.fill_diagonal(normal_rates, 1.0)
        np.fill_diagonal(boot_rates, 1.0)
        rate_difference = boot_rates - normal_rates

        # Determine one ordering from the normal runs, then use that same
        # ordering for all panels so the heatmaps remain directly comparable.
        if snv_count > 1:
            ordering_rates = np.nan_to_num(normal_rates, nan=0.0)
            ordering_rates = (ordering_rates + ordering_rates.T) / 2.0
            np.fill_diagonal(ordering_rates, 1.0)

            distance = np.clip(1.0 - ordering_rates, 0.0, 1.0)
            np.fill_diagonal(distance, 0.0)
            snv_order = leaves_list(
                linkage(
                    squareform(distance, checks=False),
                    method="average",
                )
            )
        else:
            snv_order = np.array([0], dtype=int)

        normal_ordered = normal_rates[np.ix_(snv_order, snv_order)]
        boot_ordered = boot_rates[np.ix_(snv_order, snv_order)]
        difference_ordered = rate_difference[np.ix_(snv_order, snv_order)]

        figure, axes = plt.subplots(
            1,
            3,
            figsize=(19, 6),
            constrained_layout=True,
        )

        normal_image = axes[0].imshow(
            normal_ordered,
            vmin=0.0,
            vmax=1.0,
            cmap="viridis",
            interpolation="nearest",
            aspect="auto",
        )
        axes[0].set_title("Normal co-clustering")
        axes[0].set_xlabel("SNVs")
        axes[0].set_ylabel("SNVs")

        axes[1].imshow(
            boot_ordered,
            vmin=0.0,
            vmax=1.0,
            cmap="viridis",
            interpolation="nearest",
            aspect="auto",
        )
        axes[1].set_title("Bootstrap co-clustering")
        axes[1].set_xlabel("SNVs")
        axes[1].set_ylabel("SNVs")

        finite_differences = difference_ordered[
            np.isfinite(difference_ordered)
        ]
        maximum_difference = (
            float(np.max(np.abs(finite_differences)))
            if finite_differences.size
            else 1.0
        )
        if maximum_difference == 0.0:
            maximum_difference = 1.0

        difference_image = axes[2].imshow(
            difference_ordered,
            vmin=-maximum_difference,
            vmax=maximum_difference,
            cmap="coolwarm",
            interpolation="nearest",
            aspect="auto",
        )
        axes[2].set_title("Bootstrap minus normal")
        axes[2].set_xlabel("SNVs")
        axes[2].set_ylabel("SNVs")

        # Thousands of SNV labels would be unreadable, so only display them
        # for small datasets. These labels are internal zero-based SNV IDs.
        if snv_count <= 40:
            positions = np.arange(snv_count)
            labels = [str(int(snv_id)) for snv_id in snv_order]
            for axis in axes:
                axis.set_xticks(positions)
                axis.set_yticks(positions)
                axis.set_xticklabels(labels, rotation=90, fontsize=6)
                axis.set_yticklabels(labels, fontsize=6)
        else:
            for axis in axes:
                axis.set_xticks([])
                axis.set_yticks([])

        figure.colorbar(
            normal_image,
            ax=[axes[0], axes[1]],
            label="Proportion of runs assigned together",
            shrink=0.8,
        )
        figure.colorbar(
            difference_image,
            ax=axes[2],
            label="Change in co-clustering rate",
            shrink=0.8,
        )
        figure.savefig(
            output_dir / "coclustering_heatmaps.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(figure)

        # Plot the median and interquartile range of each descending DCF rank.
        ranks = np.arange(1, number_of_clusters + 1)
        normal_median = np.nanmedian(normal_sorted_dcfs, axis=0)
        normal_lower = np.nanpercentile(normal_sorted_dcfs, 25, axis=0)
        normal_upper = np.nanpercentile(normal_sorted_dcfs, 75, axis=0)
        boot_median = np.nanmedian(boot_sorted_dcfs, axis=0)
        boot_lower = np.nanpercentile(boot_sorted_dcfs, 25, axis=0)
        boot_upper = np.nanpercentile(boot_sorted_dcfs, 75, axis=0)

        figure, axis = plt.subplots(
            figsize=(9, 6),
            constrained_layout=True,
        )
        axis.plot(
            ranks,
            normal_median,
            marker="o",
            linewidth=2,
            color="#2166ac",
            label="Normal median",
        )
        axis.fill_between(
            ranks,
            normal_lower,
            normal_upper,
            alpha=0.20,
            color="#2166ac",
            label="Normal interquartile range",
        )
        axis.plot(
            ranks,
            boot_median,
            marker="o",
            linewidth=2,
            color="#b2182b",
            label="Bootstrap median",
        )
        axis.fill_between(
            ranks,
            boot_lower,
            boot_upper,
            alpha=0.20,
            color="#b2182b",
            label="Bootstrap interquartile range",
        )
        axis.set_title("DCF distributions by descending rank")
        axis.set_xlabel("DCF rank, largest to smallest")
        axis.set_ylabel("DCF")
        axis.set_xticks(ranks)
        axis.set_ylim(0.0, 1.05)
        axis.grid(alpha=0.25)
        axis.legend()
        figure.savefig(
            output_dir / "dcf_rank_distributions.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(figure)

        np.savez_compressed(
            output_dir / "bootstrap_analysis_arrays.npz",
            matrix_norm_boot=matrix_norm_boot,
            normal_sorted_dcfs=normal_sorted_dcfs,
            boot_sorted_dcfs=boot_sorted_dcfs,
            normal_coclustering_rates=normal_rates,
            bootstrap_coclustering_rates=boot_rates,
            bootstrap_minus_normal=rate_difference,
            snv_order=snv_order,
        )

        print(f"Bootstrap analysis saved in: {output_dir.resolve()}")

        return {
            "normal_rates": normal_rates,
            "boot_rates": boot_rates,
            "rate_difference": rate_difference,
            "snv_order": snv_order,
        }


    def infer_dcfs(self):
        start_time = time.perf_counter()
        num_replicates = 10
        snv_count = self.data.m

        print("start dcf")
        dcf_clust = DCF_Clustering(
            rng=self.rng,
            nrestarts=18,
            cna_restriction=1,
        )

        normal_replicates = []
        normal_sorted_dcfs = []
        matrix_norm_boot = np.zeros(
            (2, snv_count, snv_count, 2),
            dtype=np.uint16,
        )
        # matrix_norm_boot = np.zeros(
        #     (2, snv_count, snv_count, 2, 2),
        #     dtype=np.uint16,
        # )

        # ------------------------------------------------------------
        # Full-data replicates
        # ------------------------------------------------------------
        for replicate_id in range(num_replicates):
            like, dcfs, _, clusters, _, _, _ = dcf_clust.run(
                self.data,
                k_vals=[self.k],
                cores=6,
            )

            dcfs = np.asarray(dcfs, dtype=float)
            dcf_order = np.argsort(-dcfs)
            dcfs_sorted = dcfs[dcf_order]
            cluster_to_rank = {
                int(cluster_id): rank
                for rank, cluster_id in enumerate(dcf_order, start=1)
            }

            normal_replicates.append(
                {
                    "like": like,
                    "dcfs": dcfs,
                    "dcfs_sorted": dcfs_sorted,
                    "cluster_to_rank": cluster_to_rank,
                    "clusters": clusters,
                }
            )
            normal_sorted_dcfs.append(dcfs_sorted)

            print(
                f"Normal replicate {replicate_id + 1}: "
                f"likelihood={like}, sorted DCFs={dcfs_sorted}"
            )

        # ------------------------------------------------------------
        # Cell-bootstrap replicates
        # ------------------------------------------------------------
        boot_replicates = []
        boot_sorted_dcfs = []

        for bootstrap_id in range(num_replicates):
            indices = self.rng.choice(
                self.data.N,
                size=self.data.N,
                replace=True,
            )
            bootstrap_data = self.data.resample_cells(indices)

            like, dcfs, _, clusters, _, _, _ = dcf_clust.run(
                bootstrap_data,
                k_vals=[self.k],
                cores=6,
            )

            dcfs = np.asarray(dcfs, dtype=float)
            dcf_order = np.argsort(-dcfs)
            dcfs_sorted = dcfs[dcf_order]
            cluster_to_rank = {
                int(cluster_id): rank
                for rank, cluster_id in enumerate(dcf_order, start=1)
            }
        # self.reality_matrix_norm_boot=self.reality_matrix_norm_boot
            boot_replicates.append(
                {
                    "like": like,
                    "dcfs": dcfs,
                    "dcfs_sorted": dcfs_sorted,
                    "cluster_to_rank": cluster_to_rank,
                    "clusters": clusters,
                }
            )
            boot_sorted_dcfs.append(dcfs_sorted)

            print(
                f"Bootstrap replicate {bootstrap_id + 1}: "
                f"likelihood={like}, sorted DCFs={dcfs_sorted}"
            )

        # ------------------------------------------------------------
        # Pairwise normal co-clustering counts
        # ------------------------------------------------------------
        for normal_rep in normal_replicates:
            for snv_i in range(snv_count):
                segment_i = self.data.snv_to_seg[snv_i]
                cluster_i = normal_rep["clusters"].get(
                    segment_i,
                    {},
                ).get(snv_i)
                # real_cluster_j = real_cluster_assign[snv_i]
                if cluster_i is None:
                    continue

                for snv_j in range(snv_i + 1, snv_count):
                    segment_j = self.data.snv_to_seg[snv_j]
                    cluster_j = normal_rep["clusters"].get(
                        segment_j,
                        {},
                    ).get(snv_j)
                    
                    # real_cluster_j = real_cluster_assign[snv_j]
                    if cluster_j is None:
                        continue

                    same_cluster = int(cluster_i == cluster_j)
                    # reality_cluster = int(real_cluster_i == real_cluster_j)
                    matrix_norm_boot[0, snv_i, snv_j, same_cluster] += 1
                    matrix_norm_boot[0, snv_j, snv_i, same_cluster] += 1
                    # reality_matrix_norm_boot[0, snv_i, snv_j, same_cluster, reality_cluster] += 1
                    # reality_matrix_norm_boot[0, snv_j, snv_i, same_cluster, reality_cluster] += 1
        # ------------------------------------------------------------
        # Pairwise bootstrap co-clustering counts
        # ------------------------------------------------------------
        for boot_rep in boot_replicates:
            for snv_i in range(snv_count):
                segment_i = self.data.snv_to_seg[snv_i]
                cluster_i = boot_rep["clusters"].get(
                    segment_i,
                    {},
                ).get(snv_i)

                # real_cluster_j = real_cluster_assign[snv_i]
                if cluster_i is None:
                    continue

                for snv_j in range(snv_i + 1, snv_count):
                    segment_j = self.data.snv_to_seg[snv_j]
                    cluster_j = boot_rep["clusters"].get(
                        segment_j,
                        {},
                    ).get(snv_j)
                    # real_cluster_j = real_cluster_assign[snv_j]
                    if cluster_j is None:
                        continue

                    same_cluster = int(cluster_i == cluster_j)
                    # reality_cluster = int(real_cluster_i == real_cluster_j)

                    matrix_norm_boot[1, snv_i, snv_j, same_cluster] += 1
                    matrix_norm_boot[1, snv_j, snv_i, same_cluster] += 1
                    # reality_matrix_norm_boot[1, snv_i, snv_j, same_cluster, reality_cluster] += 1
                    # reality_matrix_norm_boot[1, snv_j, snv_i, same_cluster, reality_cluster] += 1
        normal_sorted_dcfs = np.asarray(normal_sorted_dcfs, dtype=float)
        boot_sorted_dcfs = np.asarray(boot_sorted_dcfs, dtype=float)

        self._export_dcf_bootstrap_analysis(
            matrix_norm_boot=matrix_norm_boot,
            normal_sorted_dcfs=normal_sorted_dcfs,
            boot_sorted_dcfs=boot_sorted_dcfs,
            output_dir="bootstrap_results",
            # reality_matrix_norm_boot=reality_matrix_norm_boot
        )

        # Use the median of each descending bootstrap rank. Axis 0 aggregates
        # across bootstrap replicates while preserving the DCF ranks.
        consensus_dcfs = np.median(boot_sorted_dcfs, axis=0)
        self.delta = {
            rank: float(dcf)
            for rank, dcf in enumerate(consensus_dcfs)
        }

        # Retain the results on the object for later notebook analysis.
        self.matrix_norm_boot = matrix_norm_boot
        self.normal_sorted_dcfs = normal_sorted_dcfs
        self.boot_sorted_dcfs = boot_sorted_dcfs
        # self.reality_matrix_norm_boot=self.reality_matrix_norm_boot
        elapsed_time = time.perf_counter() - start_time
        print("Bootstrap consensus DCFs:", self.delta)
        print(f"Executed in {elapsed_time:.4f} seconds")
        sys.exit()
            
        def preprocess(self, seg_list, delta):

            if self.cores <= 1:
                stis = {}
                dictlist = []
                for ell in seg_list:
                    _, stis[ell]  = self.preprocess_helper(ell, delta)
            
            else:
                args = [(ell, delta) for ell in seg_list]
                with multiprocessing.Pool(self.cores) as pool:
                    dictlist = pool.starmap(self.preprocess_helper, args)
                
                stis = {ell: vals for ell, vals in dictlist}
            return stis
        
    

    @timeit_decorator
    def fit(self, data, lamb=1e3, segments= None, cores=1, Tm=None):
        '''
        @params Data data: the input data (C,A,D) to fit
        @params float lamb (float): a regularization parameter the cost function
        @params list segments: the list of segment ids to fit
        @params int cores: the number of processors to use 

        Fits a clonal tree T (with assosciated genotypes) and a mapping phi of cells to clones with minimum cost
        for the input data and specified segments. 

        returns a list of the top_n Solutions to the Clonal Tree Inference with Copy Number problem (CTICN)
        ''' 
        
        self.data = data


        self.lamb = lamb 
        self.cores = cores 
        


        if segments is None:
            segments = data.segments

        init_segs, infer_segs, place_segs, no_snvs_segs = self.partition_segments(segments, min_cn_states=2)
        print(f"\nSegment partition:\ninitial segments: {len(init_segs)}\ninference segments: {len(infer_segs)}\nplace segments: {len(place_segs)}\n")
        
        print("Plowing the field.... ")
        if self.delta is None:
            self.infer_dcfs()
        print("self.delta: ", self.delta)
        if Tm is not None:
            scriptTm = [Tm]
            for T_m in scriptTm:
                for n in T_m:
                    if n not in self.delta:
                        raise ValueError(f" Node {n} does not match a cluster id. \
                                         Each node label in the mutation cluster \
                                         tree must map to a unique value in [k] ")
        else:
            scriptTm = self.enumerate_mutcluster_trees(self.delta)


                    
    
  
       
        
        all_best_trees = []

        delta = self.delta.copy()
 
        # while loop < self.max_loops and len(scriptTm) > 0:
        print(f"DCFs delta: {delta}")
        print(f"Starting mutation cluster trees iteration with {len(scriptTm)} trees...")
        
        stis_init = self.preprocess(init_segs, delta)
    
        print("Planting the seeds.... ")
        init_order = self.order_segments(init_segs)
        print("Segment integration order:")
        print(init_order)
        init_trees, costs = self.infer(scriptTm, stis_init, init_order=init_order)

        self.clonal_trees = init_trees
    

        # return utils.concat_and_sort(init_trees)

        #identify the mutation cluster trees that yield minimum cost over the initial segments
        sorted_indices = sorted(range(len(costs)), key=lambda i: costs[i])
        if self.ninit_Tm is None or len(sorted_indices) <= self.ninit_Tm:
            smallest_indices = sorted_indices
        else:
            smallest_indices = sorted_indices[:self.ninit_Tm]
    
        
        print(f"Best mutation cluster trees:")
        for i in smallest_indices:
            print(f"{i}: {list(scriptTm[i].edges)}")
            # if i == self.ground_truth_tm:
            #     print("Including the ground truth mutation cluster tree!")
        
        best_tree_int = get_top_n(self.clonal_trees, self.top_n)
        print("Best trees after initial integration ")
        for i,sol in enumerate(best_tree_int):
            cost, snv, cna = sol.compute_likelihood(self.data, self.lamb)
            

    
        print("\nWatering the fields.... ")
        if len(infer_segs) > 0:
            init_Tm = [scriptTm[i] for i in smallest_indices]
            init_order = self.order_segments(infer_segs)
            stis_infer = self.preprocess(infer_segs, delta)
            self.clonal_trees, costs = self.infer(init_Tm, stis_infer, 
                                                    [init_trees[i] for i in smallest_indices], 
                                                    init_order = init_order )
        
        best_trees =  get_top_n(self.clonal_trees, self.top_n)

        # print(f" tree | cost | snv | cna")
        for i,b in enumerate(best_trees):
            cost, snv, cna = b.compute_likelihood(self.data, self.lamb)
            # print(f"|{i} | {cost} | {snv} | {cna} |")
            
        self.post_process(best_trees)

        self.place_snvs(best_trees, place_segs)

    
    
        print("\nHarvesting....")
        for sol in best_trees:
            sol.optimize(self.data, self.lamb)
        
        all_best_trees.append(best_trees)
        best_trees = sorted(best_trees, key=lambda x: x.cost)
        # print(f" tree | cost | snv | cna")
        for i,b in enumerate(best_trees):
            cost,snv, cna = b.compute_likelihood(self.data, self.lamb)
            # print(f"|{i} | {cost} | {snv} | {cna} |")
 
    
        best_trees = get_top_n(all_best_trees, self.top_n)
        for sol in best_trees:
            sol.prune_leaves(self.k)



        return  best_trees
        


  

    def post_process(self, sol_list):
        """
        Removes linear chain from each tree with no cells assigned and
        maps SNVs in segments with only 1 CN states 
        
        This function takes in a Solution list and post-processes the clonal tree of each solution
        such that linear chains with no cell assignments are removed from the tree. Each clonal
        tree object is modified in place. 

        :param list sol_list: a list of Solution objects
 

        """

        for sol in sol_list:
            sol.optimize(self.data, self.lamb)
            
       
        
        