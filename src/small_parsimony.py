import networkx as nx
import numpy as np


# from polytomy_resolver import PolytomyResolver as pr
np.seterr(all="ignore")

class SmallParsimony:
    def __init__(self, T, root,  alphabet=None, cost=None):
        
        self.T =T
        self.root =root
        self.cost = cost
        if alphabet is None:
            self.alphabet = ["A", "C", "G", "T", "N", "-"]
        else:
            self.alphabet = alphabet

        if self.cost is None:
            self.cost= {}

            for i in self.alphabet:
                for j in self.alphabet:
                    if (i,j) not in self.cost:
                        if i == "-" and j == "-":
                            self.cost[i,j] =np.Inf
                        if i == "-" or j == "-":
                            self.cost[i,j] = 1
                        elif  i ==j:
                            self.cost[i,j] = 0
                        
                        else:
                            self.cost[i,j] = 1

        self.postorder = self.postorder_traversal()
        self.nodes = list(self.T.nodes())

      
        # self.att_name = attribute
        # if attribute == "SEQ":
        #     self.att=sequences
        # else:
        #     self.att = isotypes
         #nx.get_node_attributes(self.T, attribute)
        
        
    

    def postorder_traversal(self):
        return list(nx.dfs_postorder_nodes(self.T, source=self.root))
    

    def preorder_traversal(self):
        return list(nx.dfs_preorder_nodes(self.T, source=self.root))
    
    def is_leaf(self, node):
        return self.T.out_degree(node) ==0

 

    
    def sankoff_poly(self, isotypes, weights, states, ighd=None):
        total_nodes_added = 0
        self.att= isotypes 
        dp_matrix = {n : 0 for n in self.nodes}
        dp_bt = {}



    


        for n in self.postorder:
            
            if self.is_leaf(n):
                dp_matrix[n] =0
                dp_bt[n] = self.att[n]
            
            
            
            else:
                children = [c  for c in self.T.successors(n)]
                if n == self.root:
                    dp_bt[n] = 0
                else:
                    child_states = [dp_bt[c] for c in children]
                    if ighd is not None:
                 
                        if any([c ==ighd for c in child_states]):
                            if  all([c ==ighd for c in child_states]):
                                dp_bt[n] = ighd 
                            else:
                                dp_bt[n] = min(states)
                        else:
                            dp_bt[n] = min(child_states)
                    else:
                        dp_bt[n] = min(child_states)
             

                s = dp_bt[n]
                state_counts = {t:0 for t in states}
                state_lists = {t: [] for t in states}
                for c in children:
                    state_counts[dp_bt[c]] += 1
                    state_lists[dp_bt[c]].append(c)
            
                if len(children) > 2:
                    for t in state_counts:
                        
                        if state_counts[t] > 1 and t != s:
                            
                            #add in new node to partially resolve polytomy
                            poly_node = len(list(self.T.nodes)) + 1
                            while poly_node in self.T:
                                poly_node += 1
                            total_nodes_added += 1
                            dp_bt[poly_node] = t
                            dp_matrix[n] += weights[s,t]
                            for kid in state_lists[t]:
                                    self.T.remove_edge(n, kid)
                                
                                    self.T.add_edge(n, poly_node)
                                    self.T.add_edge(poly_node, kid)
                                    dp_matrix[n] += dp_matrix[kid] + weights[t,t] 
                        else:
                            for kid in state_lists[t]:
                                dp_matrix[n] += dp_matrix[kid] + weights[s,t]
                
                else:
                    
                    for c in children:
                        t = dp_bt[c]
                        dp_matrix[n] += dp_matrix[c] + weights[s,t]
        # print(total_nodes_added)
        return dp_matrix, dp_bt

       


               

        

        
    
    def polytomy_resolver(self, leaf_labels,  weights, states, ighd=None):
            dp_mat, dp_bt = self.sankoff_poly(leaf_labels, weights, states, ighd=ighd)

            # # dp_mat, dp_bt, dp_poly = self.sankoff_polytomy2(leaf_labels, weights, states)
            # score, labels, tree = self.polytomy_backtrace(dp_mat, dp_bt,dp_poly, start_state)

 
            return dp_mat[self.root], dp_bt, self.T




    def sankoff_dp(self, pos):
        dp_matrix = {n : {}  for n in self.nodes}
     
        for n in self.postorder:

            if self.is_leaf(n):
                for a in self.alphabet:

                    if self.att[n][pos] == a:
                        dp_matrix[n][a] = 0
                    else:
                        dp_matrix[n][a] = np.Inf
            
  

            elif n == self.root:
                base = self.att[n][pos]
                for a in self.alphabet:
                    if a == base:
                        dp_matrix[n][a] = np.sum([[self.min_cost(c,a, dp_matrix) for c in self.T.successors(n)]])
                    else:
                        dp_matrix[n][a] = np.Inf
            
            else:
                
          
                for a in self.alphabet:
                    score_array = np.array([self.min_cost(c,a, dp_matrix) for c in self.T.successors(n)])
                    
                
                    dp_matrix[n][a] = score_array.sum()

        opt_score = np.array([dp_matrix[self.root][a] for a in self.alphabet]).min()
        
        return opt_score, dp_matrix

    @staticmethod
    def argmin( my_dict):
        minval = np.Inf
        
        for key,val in my_dict.items():
            if val <= minval:
                minval = val
                minarg = key
            
        return minarg


    def traceback(self, dp_matrix):
        labels = {}
        labels[self.root] = self.argmin(dp_matrix[self.root])
        for c in self.T.successors(self.root):
            self.traceback_helper(labels[self.root], c, dp_matrix, labels)
        
        return labels
        
     
            
    def traceback_helper(self,state, root, dp_matrix, labels):
        self.find_state(state, root, dp_matrix, labels)
        children = list(self.T.successors(root))

        for c in children:
            self.traceback_helper(labels[root], c,dp_matrix, labels)
    


    def find_state(self, parent_state, node,  dp_matrix, labels):

        min_cost = np.Inf
        for b in self.alphabet:
           
        
            trans_cost = self.cost[parent_state,b] + dp_matrix[node][b]
            if trans_cost < min_cost:
                min_cost = trans_cost
                labels[node] = b

    @staticmethod 
    def concat_labels(all_labs, nodes):
        seq_assign = {}
        for k in nodes:
            seq_assign[k] = [all_labs[i][k] for i in all_labs if k in all_labs[i]]
          
            
        return seq_assign


    def sankoff(self, sequences):
        self.att = sequences
        for key, seq in self.att.items():
            seq_length = len(seq)
            break
        # if type(self.att[self.root]) not in [int, np.int64]:

    
        # else:
        #     seq_length = 1
        #     self.att = {key: [val] for key,val in self.att.items()}

        all_labels = {}
 
        min_scores =np.zeros(seq_length, dtype=float)
        for i in range(seq_length):
            min_scores[i], dp_mat = self.sankoff_dp(i)
            label_dict = self.traceback(dp_mat)
            all_labels[i] = label_dict
        
        node_labels= self.concat_labels(all_labels, self.nodes)
        if seq_length == 1:
            # for key in node_labels:
            #     if key in self.att:
            #         node_labels[key] = self.att[key]
            node_labels = {key: val[0] for key,val in node_labels.items()}
        # nx.set_node_attributes(self.T, node_labels, self.att_name)

        return min_scores.sum(),node_labels

            

    def min_cost(self, child, to_state, dp_matrix):
   

        scores = np.array([self.cost[to_state,a] + dp_matrix[child][a]for a in self.alphabet])
        
        return scores.min()
    

                
    
  
    

            
    


                        


   



    

        

