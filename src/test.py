import cProfile
import utils

instance = "s11_m5000_k25_l5"
folder = "n1000_c0.05_e0" 
pth = f"simulation_study/input"

gtpth = "test"
profiler = cProfile.Profile()
profiler.enable()

f"{pth}/{instance}/{folder}/data.pkl"
Tm_edges = [(3,4), (1,3), (1,0), (1,2)]
        #  
        #         pickle_object(self, "test/sti.pkl")
        #         draw(T, "test/T.png")
# data = utils.load_pickled_object("test/data19.pkl" )

# scriptT = utils.load_pickled_object("test/scriptTm.pkl")
obj = utils.load_pickled_object("test/sti.pkl")

sol = obj.fit(Tm_edges, obj.data, 2)

utils.pickle_object(sol, "test/solutions.pkl")

# for i, Tm in enumerate(scriptT):
#     print(i)
#     Tm_edges = list(Tm.edges)
#     sols = 

profiler.disable()

profiler.dump_stats("test/profile.prof")
print("done")
# sols[0].png("test/sol.png")

#TM 4
# Tm = scriptT[13:][4]
# utils.draw(Tm, "test/mut_cluster_tree.png")
#refinement 12

# gt = cnmerge.genotypes1
# cnmerge.verbose = True
# utils.draw(cnmerge.T2, "test/T2.png")
# utils.draw(cnmerge.T1, "test/T1.png")

# trees = cnmerge.fit(data, 1e3, 3)



# cnmerge.verbose = True

# tree_list = cnmerge.fit(data, 1e3, 3)






