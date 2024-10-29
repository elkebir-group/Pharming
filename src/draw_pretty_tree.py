SET3 =  [
    "#8DD3C7",  # Light turquoise
    "#FFFFB3",  # Light yellow
    "#BEBADA",  # Light purple
    "#FB8072",  # Light coral
    "#80B1D3",  # Light blue
    "#FDB462",  # Light orange
    "#B3DE69",  # Light green
    "#FCCDE5",  # Light pink
    "#D9D9D9",  # Light grey
    "#BC80BD",  # Light lavender
    "#CCEBC5",  # Light mint green
    "#FFED6F"   # Light gold
]

SET1 = [
            "#E41A1C",
            "#377EB8",
            "#4DAF4A",
            "#984EA3",
            "#FF7F00",
            "#FFFF33",
            "#A65628",
            "#F781BF",
            "#999999"
    ]

DARK2 = [
    "#1B9E77",  # Dark green
    "#D95F02",  # Dark orange
    "#7570B3",  # Dark purple
    "#E7298A",  # Dark pink
    "#66A61E",  # Dark lime green
    "#E6AB02",  # Dark yellow
    "#A6761D",  # Dark brown
    "#666666"   # Dark grey
]

COMBINED = SET3 + DARK2
import pygraphviz as pgv
import networkx as nx 


import string

def get_alphabet_labels(n):
    alphabet = list(string.ascii_uppercase)
    labels = []
    
    for i in range(n):
        # Handle the wrap-around using divmod (like excel columns: A, B, ..., Z, AA, AB, ...)
        div, mod = divmod(i, 26)
        label = (alphabet[div-1] if div > 0 else '') + alphabet[mod]
        labels.append(label)
        
    return labels


class DrawPrettyTree:
    def __init__(self, sol, colors:list=COMBINED, min_cells=5, include_CNAs=True, segment_dict=None, include_cell_counts=True) -> None:
        self.ct = sol.ct
        self.phi = sol.phi    
        self.ct.update_mappings()
        self.mut_map= self.ct.mut_mapping
        self.loss = self.ct.mut_loss_mapping
        self.colors = colors 
        self.cell_count = self.phi.get_cell_count()
        self.cell_nodes = sorted([n for n,count in self.cell_count.items() if count >= min_cells ])
        non_cell_nodes = [n for n in self.ct.clones() if n not in self.cell_nodes]
        node_ids = self.cell_nodes + non_cell_nodes
        self.alphabet_labels = get_alphabet_labels(len(node_ids))
        # print(len(self.alphabet_labels))
        # print(self.alphabet_labels)
        self.idTolabel = {id: label for label, id in enumerate(node_ids)}
        # self.idToOrder = {id: label for label, id in enumerate(node_ids)}   
        self.labelToid = {label: id for id, label in self.idTolabel.items()}
        self.segment_dict = segment_dict

        self.T = pgv.AGraph(strict=False, directed=True)
        self.include_CNAs= include_CNAs
        self.include_cell_counts = include_cell_counts

    def add_node_to_tree(self, nodeID):
        lab = self.idTolabel[nodeID]
        letter = self.alphabet_labels[lab]
        nodeLab = f"{letter}\n"
  

        if nodeID in self.cell_nodes and lab < len(self.colors):
            col = self.colors[self.idTolabel[nodeID]]   
        else:
            col ="#FFFFFF"
        
        if self.cell_count[nodeID] > 0 and self.include_cell_counts:
                nodeLab += f"{self.cell_count[nodeID]} cells\n\n"
        if self.include_CNAs:
            nodeLab += self.new_cn_states(nodeID)
        self.T.add_node(nodeID, label=nodeLab, fillcolor=col, style="filled", fontsize=20) 
# G.add_node(3, label="Node 3", fillcolor="#3498DB", style="filled")
        if nodeID != self.ct.root:
            self.add_edge_to_tree(nodeID)
  

            
    def add_edge_to_tree(self,nodeID):
        lab = ''
        parID = self.ct.parent(nodeID)
        if nodeID in self.mut_map and len(self.mut_map[nodeID]) > 0:
            lab += f'+{len(self.mut_map[nodeID])}\n'
        if nodeID in self.loss and len(self.loss[nodeID]) > 0:
             lab += f'-{len(self.loss[nodeID])}'
        self.T.add_edge(parID, nodeID, label=lab, fontsize=20)

        
    def save_labels(self, fname):
        with open(fname, "w+") as file:
            file.write("id,label\n")
            for id, lab in self.idTolabel.items():
                file.write(f"{id},{self.alphabet_labels[lab]}\n")

    def new_cn_states_discrete(self, nodeID):
        mystr = ""
        if nodeID == self.ct.root:
            return mystr
        changes = {}
        par = self.ct.parent(nodeID)
        cna_genos = self.ct.get_cna_genos()
        for ell in cna_genos:
            parCN = cna_genos[ell][par]
            cn= cna_genos[ell][nodeID]
            if parCN != cn:
                if self.segment_dict is not None:   
                    label_ell = self.segment_dict[ell]
                else:
                    label_ell = ell
                changes[ell] =f"{label_ell}: {parCN[0]}|{parCN[1]} to {cn[0]}|{cn[1]}\n"

        changes = dict(sorted(changes.items()))
        for ell, cstr in changes.items():
            mystr += cstr
        return mystr
    
    def new_cn_states(self, nodeID):
        mystr = ""
        if nodeID == self.ct.root:
            return mystr
        
        changes = {}
        par = self.ct.parent(nodeID)
        cna_genos = self.ct.get_cna_genos()
        
        for ell in cna_genos:
            parCN = cna_genos[ell][par]
            cn = cna_genos[ell][nodeID]
            if parCN != cn:
                # Get the (chr, start, end) tuple from segment_dict if available, else default to ell
                if self.segment_dict is not None:
                    chr, start, end = self.segment_dict[ell]
                    label_ell = f"chr{chr} {start:.1f}-{end:.1f} Mb"
                else:
                    label_ell = ell
                changes[ell] = ((chr, start, end), f"{parCN[0]}|{parCN[1]} to {cn[0]}|{cn[1]}")
        
        # Sort changes by segment ID for consecutive grouping
        sorted_changes = dict(sorted(changes.items()))
        
        # Group consecutive segments with the same state change
        grouped_changes = []
        current_range = []
        current_change = None

        for ell, ((chr, start, end), change) in sorted_changes.items():
            if current_change is None:
                # Initialize the first segment range
                current_range = [(chr, start, end)]
                current_change = change
            elif change == current_change and chr == current_range[-1][0] and start == current_range[-1][2]:
                # Continue current range if the change and chromosome are the same
                current_range.append((chr, start, end))
            else:
                # Append current range to grouped_changes with combined label
                range_start = current_range[0][1]
                range_end = current_range[-1][2]
                grouped_changes.append(f"chr{chr} {range_start:.1f} Mb - {range_end:.1f} Mb : {current_change}\n")
                
                # Reset for the new range
                current_range = [(chr, start, end)]
                current_change = change

        # Add the last range
        if current_range:
            range_start = current_range[0][1]
            range_end = current_range[-1][2]
            grouped_changes.append(f"chr{chr} {range_start:.1f} Mb - {range_end:.1f} Mb : {current_change}\n")
        
        # Concatenate all grouped changes
        mystr = "".join(grouped_changes)
        return mystr


    def draw(self, fname, horizontal=False):
        for n in self.ct.preorder():
            self.add_node_to_tree(n) 
        if horizontal:
            self.T.layout(prog="dot", args="-Grankdir=LR")
        else:
            self.T.layout("dot")
        self.T.draw(fname)





