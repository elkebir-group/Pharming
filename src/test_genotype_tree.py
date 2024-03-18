from genotype_tree import GenotypeTree
import clonelib

S_edges =[((1,1), (3,3))]
snv_trees = clonelib.get_genotype_trees(S_edges)

genotype_trees = [GenotypeTree(snv_tree) for snv_tree in snv_trees]


cn_prop = {(1,1): 0.004, (3,3): 0.996}

dcf = 0.781

print(f"dcf: {dcf}, cn proportions: {cn_prop}")
print(f"F: {genotype_trees[0].compute_F(cn_prop)}")
for i in range(len(genotype_trees)):

    T = genotype_trees[i]
    print(f"Genotype Tree {i}:")
    print(T)
    v_minus  = T.v_minus(cn_prop)
    v = genotype_trees[i].dcf_to_v(dcf, cn_prop)
    v_plus = T.v_plus(cn_prop, v_minus)
    isvalid = v >= v_minus and v <= v_plus
    print(f"v: {v} v-: {v_minus} v+: {v_plus}, IS VALID: {isvalid}\n")
    print(T.posterior_dcf(dcf, 5,10, cn_prop))
    
    




