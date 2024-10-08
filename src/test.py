import cProfile
import utils
import networkx as nx
import clonelib
import pandas as pd 
from scipy.stats import binomtest
import itertools
import numpy as np

sim = "s10_m5000_k25_l5_d2/n1000_c0.01_e0"
instance = f"dcf_clustk6/weighted-random/isegs10_tm10_top5_lamb1000/{sim}"

dat = pd.read_pickle(f"simulation_study/sims3/{sim}/data.pkl")
gt = pd.read_pickle(f"simulation_study/sims3/{sim}/gt.pkl")
sols = pd.read_pickle(f"simulation_study/pharming_test/{instance}/solutions.pkl")
sol = sols[0]

inf_lost = pd.read_csv(f"simulation_study/pharming_test/{instance}/pred_mut_loss.csv")
gt_lost = pd.read_csv(f"simulation_study/sims3/{sim}/mut_loss_clust_gt.csv")

print(inf_lost.head())
print(gt_lost.head())

inf_lost_dict = dict(zip(inf_lost["mutation"], inf_lost["cluster"]))
gt_lost_dict = dict(zip(gt_lost["mutation"], gt_lost["cluster"]))
results = []
for key, val in inf_lost_dict.items():
    status = "TP" if key in gt_lost_dict and gt_lost_dict[key] == val else "FP"

    desc_nodes = sorted(nx.descendants(sol.ct.tree, val) | {val})
    desc_cells = list(itertools.chain.from_iterable(sol.phi.get_cells(v) for v in desc_nodes))
    if len(desc_cells) == 0:
        print(val)
    var_total = dat.var[desc_cells,:][:,key].sum()
    total = dat.total[desc_cells,:][:,key].sum()
    vaf = var_total/total
    if total > 0:
        res = binomtest(var_total, total, 0.001, alternative='greater')
        pval = res.pvalue
    else:
        pval = np.nan
    
    results.append((key, val, status, var_total, total, vaf, pval, len(desc_cells)))

    print(f"Key: {key}, status: {status} Var: {var_total}, Total: {total} VAF: {vaf} p-value: {pval}")

pd.DataFrame(results, columns=["mutation", "cluster", "status", "var", "total", "vaf", "pval", "cells"]).to_csv(f"simulation_study/pharming_test/{instance}/lost_snv_prob.csv", index=False)
# result = binomtest(k, n, p)

print(len(sols))
print('done')



# gtpth = "dlp"
# ell = 342

# def enumerate_cna_trees_python(cn_states, start_state):
#     cn_states = [cn for cn in cn_states if cn != start_state]


#     G = nx.DiGraph()
#     G.add_nodes_from(cn_states)
#     for u in G.nodes:
#         for v in G.nodes:
#             if u == v:
#                 continue
#             if (u[0] == 0 and not v[0] > 0) or (u[1] == 0 and not v[1] > 0):
#                 continue
#             G.add_edge(u,v, weight=1)

#     G.add_node(start_state)
#     for u in cn_states:
#         G.add_edge(start_state, u,weight=1)

#     cnatrees = nx.algorithms.tree.branchings.ArborescenceIterator(G)
#     scriptS= list(cnatrees)
#     return scriptS

# dat = pd.read_pickle(f"dlp/input/data.pkl")
# cn_prop = dat.thresholded_cn_prop(ell, 0.05, (1,1), include_start_state=False)
# print(cn_prop)
# # clonelib.get_cna_trees(set(cn_prop.keys()), 1,1)

# S= [((3, 7), (1, 2)), ((3, 7), (2, 0)), ((3, 0), (3, 7)), ((1, 1), (3, 0))]
# # scriptS = enumerate_cna_trees_python(list(cn_prop.keys()), (1,1))
# utils.draw(nx.DiGraph(S), "test/S.png")
# snv_trees =clonelib.get_genotype_trees(S)
# for i,t in enumerate(snv_trees):
#     utils.draw(nx.DiGraph(t), f"test/T{i}.png")
# # print(len(scriptS))

#226,94, 341]]




# S= nx.DiGraph([((0,3), (0,1)), ((1,1), (0,3))])


# utils.draw(S, "test/S.png")
# T1 = nx.DiGraph([((1,1,0,0), (0,3,0,0)), ((0,3,0,0), (0,3,0,1)), ((0,3,0,1),(0,1,0,1))])
# T2 = nx.DiGraph([((1,1,0,0), (0,3,0,0)), ((0,3,0,0), (0,3,0,1)), ((0,3,0,1),(0,1,0,0))])
# utils.draw(T1, "missing_tree0.png")
# utils.draw(T2, "missing_tree1.png")


# trees = [nx.DiGraph(T) for T in snv_trees]
# for i,t in enumerate(trees):
#     utils.draw(t, f"test/T{i}.png")
# instance = "s11_m5000_k25_l5"
# folder = "n1000_c0.05_e0" 
# pth = f"simulation_study/input"

# gtpth = "test"
# profiler = cProfile.Profile()
# profiler.enable()

# f"{pth}/{instance}/{folder}/data.pkl"
# Tm_edges = [(3,4), (1,3), (1,0), (1,2)]
#         #  
#         #         pickle_object(self, "test/sti.pkl")
#         #         draw(T, "test/T.png")
# # data = utils.load_pickled_object("test/data19.pkl" )

# # scriptT = utils.load_pickled_object("test/scriptTm.pkl")
# obj = utils.load_pickled_object("test/sti.pkl")

# sol = obj.fit(Tm_edges, obj.data, 2)

# utils.pickle_object(sol, "test/solutions.pkl")

# # for i, Tm in enumerate(scriptT):
# #     print(i)
# #     Tm_edges = list(Tm.edges)
# #     sols = 

# profiler.disable()

# profiler.dump_stats("test/profile.prof")
# print("done")
# # sols[0].png("test/sol.png")

# #TM 4
# # Tm = scriptT[13:][4]
# # utils.draw(Tm, "test/mut_cluster_tree.png")
# #refinement 12

# # gt = cnmerge.genotypes1
# # cnmerge.verbose = True
# # utils.draw(cnmerge.T2, "test/T2.png")
# # utils.draw(cnmerge.T1, "test/T1.png")

# # trees = cnmerge.fit(data, 1e3, 3)



# # cnmerge.verbose = True

# # tree_list = cnmerge.fit(data, 1e3, 3)






