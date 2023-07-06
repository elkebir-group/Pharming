import networkx as nx
from dataclasses import dataclass
class GenotypeTree:

    tree: nx.DiGraph
    node_mapping: dict
    id: int = 0

    # def __init__(self, tree:nx.DiGraph=None, node_mapping: dict=None, id=0) -> None:


    def __post_init__(self):


        #dynamically find the root node
        for n in self.tree:
            if self.tree.in_degree[n] ==0:
                self.root = n 
                break

        

    def __str__(self):
        mystr = ""
        
        for u,v in self.tree.edges:
        
            mystr += f"Node {u}: {self.tree.nodes[u]['genotype']} -> Node {v}: {self.tree.nodes[v]['genotype']}\n"
        return(mystr)
    
    def preorder(self):
        return list(nx.dfs_preorder_nodes(self.tree, source=self.root))
    
    def generate_CNA_tree(self):

        node_id = -1
        T_CNA = nx.DiGraph()
            
        node_mapping = {}
        for u in nx.dfs_preorder_nodes(self.tree, source=self.root):

            u_x, u_y, u_m  = self.tree.nodes[u]["genotype"]
            if (u_x, u_y) not in node_mapping:
                node_id +=1
                u_node = node_id
                T_CNA.add_node(u_node, genotype=(u_x, u_y,0))
                node_mapping[(u_x, u_y)] = u_node
         
            else:
                u_node = node_mapping[(u_x,u_y)]
            for v in self.tree.neighbors(u):
                v_x, v_y, v_m  = self.tree.nodes[v]["genotype"]
                if (v_x, v_y) not in node_mapping:
                    node_id +=1
                    v_node =node_id
                    T_CNA.add_node(v_node, genotype=(v_x, v_y,0))
                    T_CNA.add_edge(u_node, v_node)

                    node_mapping[(v_x, v_y)] = v_node
        return GenotypeTree(T_CNA, node_mapping=node_mapping)

    def find_split_pairs(self):
        split_cand =[]
        for u, geno in self.tree.nodes("genotype"):
            x, y, m = geno
            if m==0 and geno not in split_cand:
                split_cand.append((u,geno))
            else:
                for v, v_geno in split_cand:
                    v_x, v_y, v_m = v_geno
                    if x==v_x and y==v_y:
                        return (v,v_geno), (u, geno)
    
    def find_path(self, parent_cna_geno, child_cna_geno):
        p_x, p_y = parent_cna_geno
        c_x, c_y = child_cna_geno
        parent_node = self.node_mapping[(p_x, p_y,0)]
        for x,y,m in self.node_mapping:
            if x==c_x and y==c_y:
                child_node = self.node_mapping[(x,y,m)]
        shortest_path = nx.shortest_path(self.tree, source=parent_node, target=child_node)
        genotypes_in_shortest_path = [self.tree.nodes[node]["genotype"] for node in shortest_path]
        return shortest_path, genotypes_in_shortest_path
       






    
    


    