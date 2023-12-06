from clonal_tree import load, ClonalTree
from data import load_from_pickle, Data
# from tree_merging import ClonalTreeMerge
from superimposition import Superimposition
from copy import deepcopy
import networkx as nx 
import pygraphviz as pgv
import clonelib
import numpy as np 
from genotype import genotype
import pandas as pd
import argparse




def draw(tree, fname):
    ptree = pgv.AGraph(strict=False, directed=False)
    ptree.add_edges_from(list(tree.edges))
    ptree.layout("dot")
    ptree.draw(fname)

def get_trees(segtree, m):
    cna_states = []
    cna_edge_list = []
    snv_edge_list = []

    for v in segtree.clones():

        snv_geno = segtree.genotypes[v][m]
        cna_geno = snv_geno.to_CNAgenotype()
        if cna_geno.to_tuple() not in cna_states:
            cna_states.append(cna_geno.to_tuple())
        if v != segtree.root:
            par = segtree.parent(v)
            par_snv_geno = segtree.genotypes[par][m]
            par_cna_geno = par_snv_geno.to_CNAgenotype()
            if snv_geno != par_snv_geno:
                snv_edge_list.append((par_snv_geno.to_tuple(), snv_geno.to_tuple()))
            if cna_geno != par_cna_geno:
                cna_edge_list.append((par_cna_geno.to_tuple(), cna_geno.to_tuple()))
    cna_tree, snv_tree = nx.DiGraph(), nx.DiGraph()
    if len(cna_states) > 1:
        cna_tree.add_edges_from(cna_edge_list)
    else:
        cna_tree.add_nodes_from(cna_states)
    snv_tree.add_edges_from(snv_edge_list)
    return cna_tree, snv_tree

def find_min_indices(lst):
    if not lst:
        return []

    min_value = min(lst)
    min_indices = [index for index, value in enumerate(lst) if value == min_value]
    
    return min_indices

def get_costs(seg, m, c, dat, lamb=0):

        t= nx.DiGraph()
        t.add_edges_from(c)
        geno_dict = {}
        relabel = {}
        for j, v in enumerate(t):
            geno_dict[j] = {m : genotype(*v)}
            relabel[v] = j

        t_seg_to_muts = {seg: [m]}
        t_copy = nx.relabel_nodes(t, relabel)
        ct = ClonalTree(t_copy, geno_dict, t_seg_to_muts)
        ct.assign_cells(dat, lamb=lamb)
        return ct.compute_costs(dat, lamb=lamb), t, ct
 
# print(gt.seg_to_muts)
def assign_snvs(gt, dat):
    res = []
    for lamb in [0, 0.1, 1]:
        # print(f"lambda {lamb}")
        for seg, snvs in dat.seg_to_snvs.items():
            segtree = deepcopy(gt)
            muts =segtree.seg_to_muts[seg]
            segtree.filter_snvs(muts)
            for m in snvs:
                cna_tree, snv_tree = get_trees(segtree, m)
                if len(cna_tree) ==1:
                    continue
                cand_trees = clonelib.get_genotype_trees(list(cna_tree.edges))
                vals = []
                for i,c in enumerate(cand_trees):
                    cost, t, ct = get_costs(seg, m, c, dat, lamb=lamb)
                    vals.append(cost)
                    if set(t.edges) == set(snv_tree.edges):
                        snv_gt  = i
                best_trees= find_min_indices(vals)
                res.append([lamb, seg, m, min(vals), snv_gt in best_trees ])

    df = pd.DataFrame(res, columns=["lamb", "segment", "snv","cost", "best_tree" ])

    return df 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=True,
                        help="pharming data object")
    parser.add_argument("-t", "--tree", type=str, required=True,
        help="simulated gt tree")
    parser.add_argument("-o", "--out", type=str, 
        help="outfile name")
    # parser.add_argument("-l", "--lamb", type=float, 
    #     help="lambda value")
    # parser.add_argument("-c", "--cna-only", action='store_true', 
    #     help="only consider cna portion of objective")
    args = parser.parse_args()
    # instance = "s10_m1000_k25_l7"
    # folder = "n1000_c1_e0.15" 
    # pth = f"/scratch/data/leah/pharming/simulation_study/input"

    # args = parser.parse_args([
    #     "-d", f"{pth}/{instance}/{folder}/data.pickle",
    #     "-t", f"{pth}/{instance}/{folder}/gt.pickle",
    #     "-o", f"test/snv_assignments.csv"

    # ])
    gt = load(args.tree)
    dat = load_from_pickle(args.data)
    df = assign_snvs(gt, dat)
    if args.out is not None:
        df.to_csv(args.out, index=False)
    # draw(cna_tree, "test/cna_tree.png")
    # draw(snv_tree, "test/snv_tree.png")

        # vals.append(cost)


        # draw(t, f"test/cand_tree{i}.png")

  




  



    # print(cna_tree)










print("Done")

    # SegTrees = [load(f"{pth}/segtrees/tree_g{g}.pickle") for g in segs]


# sp = Superimposition(SegTrees[1],SegTrees[0], verbose=True)
# RT = sp.solve(dat, threshold=10)
# costs = []
# for k in range(3,11):
# # k = 4
# # print(f"K: {k}\n")
#     Tree, mapping = sp.solve(dat, k)
#     # print(Tree.cost)
#     costs.append(Tree.cost)
#     if Tree.cost < np.Inf:
#         Tree.draw(f"test/tree_k{k}.png")


# Tree.draw("test/merged_tree.png", rev_mapping)






# paction = PactionSegments()
# tree, clones = paction.fit(SegTrees)


# merg = TreeMerge()
# merg.merge(SegTrees[1], SegTrees[9], dat)


# import networkx as nx
# import pygraphviz as pgv

# def draw_graph(G, fname, type_edges=False)-> None:
#         T_pg = pgv.AGraph(directed=True)
#         T_pg.add_edges_from(G.edges)
#         T_pg.layout(prog='dot')
#         edge_labels = nx.get_edge_attributes(G, 'weight')
#         edge_labels = {key: round(val,3) for key, val in edge_labels.items()}
#         if type_edges:
#             edge_types = nx.get_edge_attributes(G, 'etype')

#         def get_color_from_etype(etype):
#             if etype is None:
#                 return '#000000'
#             if etype == "ANCESTRAL":
#                 return '#FF0000'
#             else:
#                 return '#008000'
    
#         for edge, label in edge_labels.items():
        
#             T_pg.get_edge(*edge).attr['label'] = str(label)
#             if type_edges:
#                 etype = edge_types[edge]
#                 T_pg.get_edge(*edge).attr['color'] = get_color_from_etype(etype)

        
#         T_pg.layout("dot")


#         T_pg.draw(fname)


# import networkx as nx

# import gurobipy as gp
# from gurobipy import GRB

# def solve_ilp(graph, t1, t2):
#     # Create a Gurobi model
#     model = gp.Model("Minimum_Refined_Arborescence")

#     # Variable: Whether edge is selected in the arborescence
#     edge_vars = {}
#     for u, v in graph.edges:
#         edge_vars[(u, v)] = model.addVar(vtype=GRB.BINARY, name=f"edge_{u}_{v}")

#     # Objective: Minimize the total weight of selected edges
#     model.setObjective(gp.quicksum(graph[u][v]['weight'] * edge_vars[(u, v)] for u, v in graph.edges), GRB.MINIMIZE)

#     for u, v in graph.edges:
#         reachable_ancestors_t1 = set()
#         reachable_ancestors_t2 = set()

#         bfs_t1 = list(nx.bfs_edges(t1, u))
#         bfs_t2 = list(nx.bfs_edges(t2, u))

#         for ancestor in bfs_t1:
#             reachable_ancestors_t1.add(ancestor[1])

#         for ancestor in bfs_t2:
#             reachable_ancestors_t2.add(ancestor[1])

#         if v in reachable_ancestors_t1 or v in reachable_ancestors_t2:
#             model.addConstr(edge_vars[(u, v)] == 1, f"ancestor_{u}_{v}")

#     # Optimize the model
#     model.optimize()

#     if model.status == GRB.OPTIMAL:
#         arborescence = nx.DiGraph()
#         for u, v in graph.edges:
#             if edge_vars[(u, v)].X > 0.5:
#                 arborescence.add_edge(u, v, weight=graph[u][v]['weight'])
#         return arborescence
#     else:
#         return None
# # G = nx.DiGraph()
# # G.add_weighted_edges_from([(1, 2, 2), (1, 3, 4), (2, 4, 5), (2, 5, 3), (3, 6, 2), (3, 7, 6)])
# # T1 = nx.DiGraph([(1, 2), (2, 4), (3, 7)])
# # T2 = nx.DiGraph([(1, 3), (3, 6)])

# # combined_graph = nx.compose(G, T1)
# # combined_graph = nx.compose(combined_graph, T2)

# # root_node = 1
# # arborescence = chu_liu_edmonds(combined_graph, root_node)
# # print(arborescence.edges(data=True))


# # Example usage
# G = nx.DiGraph()

# T1 = nx.DiGraph()
# T2 = nx.DiGraph()
# T1.add_weighted_edges_from([(1, 2,2), (1,3,2)])
# T2.add_weighted_edges_from([(4, 5,2), (5,6,2)])


        
        

# G = nx.compose(T1, T2)
# # draw_graph(G, "test/combined_graph.png")

# for u in T1:
#     for v in T2:
#         G.add_edge(u,v, weight=1)

# for u in T2:
#     for v in T1:
#         G.add_edge(u,v,weight=1)
    
# # draw_graph(G, "test/combined_graph.png")
# root_node = 1


# # Refine the arborescence to include edges from T1 and T2
# refined_arborescence =solve_ilp(G, T1, T2)
# draw_graph(refined_arborescence, "test/refined_arbor.png")
# # print(arborescence.edges(data=True))


# # seg_gt_file = "/scratch/data/leah/pharming/sim_study/input/s13_n1000_m15000_c5_p0.25_l0/SegTrees/tree_g1.pickle"
# # seg_ph_file = "/scratch/data/leah/pharming/sim_study/pharming/s13_n1000_m15000_c5_p0.25_l0/SegTrees/tree_g1.pickle"


# # def load(fname):
# #     with open(fname, "rb") as file:
# #         ct = pickle.load(file)
# #     return ct 
# # segs = [1,9]
# # SegTrees = {g: load(f"{pth}/tree_g{g}.pickle") for g in segs}

# # def jaccard_similarity(set1, set2):
# #     # Calculate the size of the intersection
# #     intersection_size = len(set1.intersection(set2))
    
# #     # Calculate the size of the union
# #     union_size = len(set1.union(set2))
    
# #     # Calculate the Jaccard similarity
# #     jaccard_similarity = intersection_size / union_size
# #     return jaccard_similarity
    
# #     # Calculate the Jaccard distance (complement of similarity)
# #     # jaccard_distance = 1 - jaccard_similarity
    
# #     # return jaccard_distance
# # G = nx.Graph()
# # seg1_nodes = [f"g1_{n}" for n in SegTrees[1].cell_mapping]
# # seg9_nodes = [f"g9_{n}" for n in SegTrees[9].cell_mapping]
# # G.add_nodes_from(seg1_nodes, bipartite=0)  # Nodes in set1 have bipartite=0
# # G.add_nodes_from(seg9_nodes, bipartite=1)
# # for n1, cells1 in SegTrees[1].cell_mapping.items():
# #     for n9, cells9 in SegTrees[9].cell_mapping.items():
# #         if len(cells1) >0 or len(cells9) >0:
# #             jacc = jaccard_similarity(set(cells1), set(cells9))
# #             if jacc > 0.05:
# #                 g1_n = f"g1_{n1}"
# #                 g9_n = f"g9_{n9}"
                
# #                 G.add_edge(g1_n, g9_n, weight= jacc)




# # def draw_graph( G, fname)-> None:
# #     graph = pgv.AGraph(directed=False)

# # # Add nodes to the graph
# #     for node in G.nodes():
# #         graph.add_node(node)

# # # Add edges to the graph with labels (using "weight" attribute as edge labels)
# #     for edge in G.edges(data=True):
# #         source, target, data = edge
# #         weight = data["weight"]  # Get the "weight" attribute from the edge data
# #         graph.add_edge(source, target, label=str(round(weight,3)))

# # # Set layout and other attributes (optional)
# #     graph.layout(prog="dot")
# #     graph.node_attr["shape"] = "circle"

# #     graph.draw(fname, format='png', prog='dot')

# # # scores= gt.score(ph)

# # # print(scores)

# # draw_graph(G, "test/overlap.png")

