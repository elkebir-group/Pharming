# Created by: L.L. Weber
# Created on: 2023-11-14 11:23:09



from cell_assignment_problem import CellAssignment
from data import load_from_pickle, Data
from clonal_tree import ClonalTree
from copy import deepcopy
import numpy as np 
import pandas as pd 
import argparse




def main(args):

    # sample_size=100
    # timelimit = 1200
    # seed = 1026
    
    cna_only = args.cna_only
    # rng = np.random.default_rng(seed)

    dat = load_from_pickle(args.data)
    gt = load_from_pickle(args.tree)
   
    cost1_gt = gt.compute_costs(dat,cna_only=cna_only)
    cost2_gt= gt.compute_pooled_costs(dat, cna_only=cna_only)
    gt_phi = gt.get_phi()

    obj1 = deepcopy(gt)
    all_scores,nodes = obj1.assign_cells(dat, args.lamb, cna_only=cna_only)
    # print(all_scores.shape)
    obj1_phi = obj1.get_phi()
    cell_mapping = obj1.cell_mapping
    cost1 = obj1.compute_costs(dat, lamb=args.lamb, cna_only=cna_only)
    cost2_init = obj1.compute_pooled_costs(dat,args.lamb, cna_only=cna_only)

    merged_phi = pd.concat([pd.Series(gt_phi), pd.Series(obj1_phi)], axis=1)
    merged_phi.columns = ["gt", "obj1"]

    merged_phi["close"] =0
    merged_phi["min3_val"] = False


    node_index = {v: i for i,v in enumerate(nodes)}
    # for v in gt.tree:
    #     latent_genos = gt.get_latent_genotypes(v)
    #     cell_score_list.append(gt.node_snv_cost(v, dat.cells, dat, latent_genos))
    
    # all_scores = np.vstack(cell_score_list)



    for i,row in merged_phi.iterrows():

        if row["gt"] != row["obj1"]:
            scores = all_scores[:, i]
            partitioned_arr = np.partition(scores, 3)
            gt_ind =node_index[gt_phi[i]]

# Check if the value is one of the three smallest values
            merged_phi.at[i, "min3_val"] = scores[gt_ind] in partitioned_arr[:3]

            parent= gt.parent(gt_phi[i])
            children = gt.children(gt_phi[i])
            if parent == obj1_phi[i]:
                merged_phi.at[i, "close"] = 1
                continue
            for c in children:
                if c == obj1_phi[i]:
                    merged_phi.at[i, "close"] = 1
                    break
        else:
            merged_phi.at[i, "close"] = 1
            merged_phi.at[i, "min3_val"] = True
            

    # scores1 =gt.score_cells(obj1)
    scores1= {}
    scores1["perc_close"] = merged_phi["close"].sum() / merged_phi.shape[0]
    scores1["perc_exact"] = (merged_phi["gt"] == merged_phi["obj1"]).sum()/merged_phi.shape[0]
    scores1["perc_top3"] = merged_phi["min3_val"].sum() / merged_phi.shape[0]

    scores1["obj"] = 1
    scores1["cost1"] = cost1
    scores1["cost2"] = cost2_init 
    scores1["cost1_gt"] =  cost1_gt 
    scores1["cost2_gt"] = cost2_gt   
    df1 = pd.DataFrame(scores1, index=[0])
    if args.out is not None:
        df1.to_csv(args.out, index=False)

    # # print(df1.head())
    # target = {}
    # init_phi = {}
    # for k in cell_mapping:
    #    if len(cell_mapping[k])>= 25:

    #     lat_genos = gt.get_latent_genotypes(k)
    #     target[k] = gt.get_latent_vafs(lat_genos)
    #     for i in cell_mapping[k]:
    #         init_phi[i] =k
    #    else:
    #         for i in cell_mapping[k]:
    #             init_phi[i] =0
           
    
    # obj2 = deepcopy(gt)

   


    # ca = CellAssignment()


    # gt.phi_to_cell_mapping(init_phi)
    # t_vals = dat.total.sum(axis=0)
    # t_vals_prob = t_vals/t_vals.sum()
    # muts = dat.muts
    # sampled_muts = rng.choice(muts, size=sample_size, p=t_vals_prob, replace=False)
    # z_vals = {}
    # y_vals = {}
    # j_remove= []
    # obj2_pre = 0
    # for k in target:
    #     # if k in cell_mapping:
    #     #     if len(cell_mapping[k])> 10:
    #     for j in sampled_muts:
    #         var = dat.var[gt.cell_mapping[k],:]
    #         total = dat.total[gt.cell_mapping[k], :]
    #         d = total[:,j].sum(axis=0)
    #         v = var[:, j].sum(axis=0)
    #         if d > 0:
    #             y_vals[k,j] = 1/d
    #             z_vals[k,j]=np.abs(v/d - target[k][j])
    #             obj2_pre += z_vals[k,j]
    #         # else:
    #         #     if j not in j_remove:
    #         #         j_remove.append(j)


    # print(obj2_pre)

    # for k,j in z_vals:
    #     if z_vals[k,j] < 0 or z_vals[k,j] > 1:
    #         print(f"warnging cluster {k} and mut {j} outside of bounds")
    #     if y_vals[k,j] < 0 or y_vals[k,j] > 1:
    #         print(f"warning cluster {k} and mut {j} outside of bounds")

    # for k in target:
    #     if len(gt.cell_mapping[k]) < 5:
    #         print("warning too few cells in cluster k")       



    # sampled_muts = np.setdiff1d(sampled_muts, j_remove)  
    # ca.build(dat.var, dat.total, target, sampled_muts,init_phi, z_vals, y_vals)
    # obj2_val, obj2_phi = ca.solve(threads=8, timelimit=timelimit)

    # if obj2_val != np.Inf:
  
    #     obj2.phi_to_cell_mapping(obj2_phi)
    #     phi_post = obj2.get_phi()
    #     for i,k in phi_post.items():
    #         if obj1_phi[i] !=k:
    #             print(f"i:{i} phi1: {obj1_phi[i]} phi2: {k}")


    #     pooled_final = obj2.compute_pooled_costs(dat)
    #     greedy_final = obj2.compute_costs(dat)
    #     print(f"ilp obj: {obj2_val} init obj: {obj2_pre}")
    #     print(f"pre {cost2_init} post: {pooled_final} gt: {cost2_gt}")

    #     scores2 = gt.score_cells(obj2)
    #     scores2["obj"] = 2
    #     scores2["cost1"] = greedy_final
    #     scores2["cost2"] = pooled_final
    #     scores2["cost2_gt"] = cost2_gt 
    #     scores2["cost1_gt"] =  cost1_gt 
   
    # else:
    #     scores2={"obj": 2, "cost2": np.Inf,
    #             "cost2_gt": cost2_gt,
    #             "cost1_gt": cost1_gt ,
    #             "feature": "cell" 
    #             }
    

    
    # df1 = pd.DataFrame(scores1, index=[0])
    # df2 = pd.DataFrame(scores2, index=[0])
    # df_final = pd.concat([df1, df2], axis=0)
    # df_final.to_csv(f"test/results_comp_{instance}_{dpath}.csv", index=False)
    




  








if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=True,
                        help="pharming data object")
    parser.add_argument("-t", "--tree", type=str, required=True,
        help="simulated gt tree")
    parser.add_argument("-o", "--out", type=str, 
        help="outfile name")
    parser.add_argument("-l", "--lamb", type=float, 
        help="lambda value")
    parser.add_argument("-c", "--cna-only", action='store_true', 
        help="only consider cna portion of objective")
    


    
    args = parser.parse_args()
    # pth = "simulation_study/input/s10_m5000_k25_l7"
    # folder ="n500_c0.01_e0"
    # args = parser.parse_args([
    #     "-d", f"{pth}/{folder}/data.pickle",
    #     "-t", f"{pth}/{folder}/gt.pickle",
    #     "-l", "0",
    #     "--cna-only"

    # ])

    main(args)
    print("done")