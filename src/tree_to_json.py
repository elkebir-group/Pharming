import csv
import json
import networkx as nx
from clonal_tree import ClonalTree

def convertToJson(clonal_tree, segment_csv, snv_csv):

    # This will become the output JSON format of the tree
    output = { "segments":  [],
               "snvs":      [],
               "tree": {
                   "nodes": [],
                   "edges": []
               }
              }
    
    # TODO: Get list of segments from segment dictionary (to be created), then add to "segments"
    # TODO: Get list of SNVs from SNV dictionary (to be created), then add to "snvs"

    # Start by reading in the segment CSV and extracting the segment info
    with open(segment_csv, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            segment_id, chromosome, start, stop = row

            # Look for the segment in the list, then add its chromosome, start, and stop to its info
            for segment in output["segments"]:
                if segment["segment_id"] == segment_id:
                    segment["chromosome"] = chromosome
                    segment["start"] = start
                    segment["stop"] = stop

            # NOTE: Should we still add the segment info if we don't find the segment in the segment dictionary?
            # NOTE: Will the segment dictionary only include segments that are in the tree?

    # Then read in the SNV CSV and extract the SNV info
    with open(snv_csv, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            snv_id, segment_id, pos, ref, alt = row

            # Look for the SNV in the list, then add its ids, pos, ref, and alt to its info
            for snv in output["snvs"]:
                if snv["snv_id"] == snv_id:
                    snv["segment_id"] = segment_id
                    snv["pos"] = pos
                    snv["ref"] = ref
                    snv["alt"] = alt

            # NOTE: same question as above!

    # Now get the nodes on the tree along with their metadata
    for node_id in clonal_tree.genotypes:
        segments = []
        snvs = []
        for snv_id in clonal_tree.genotypes[node_id]:

            # TODO: Add segment id corresponding to this SNV using the SNV dictionary
            segment_id = 123    # Placeholder id

            x, y, x_bar, y_bar = clonal_tree.genotypes[node_id][snv_id]
            segments.append({"segment_id": segment_id, "x": x, "y": y})
            snvs.append({"segment_id": segment_id, "snv_id": snv_id, "x_bar": x_bar, "y_bar": y_bar})

        output["tree"]["nodes"].append({"node_id": node_id, "segments": segments, "snvs": snvs})

    # Get the edges on the tree
    # TODO: Check if this method works
    edge_tuples = clonal_tree.tree.edges()
    edge_list = [list(edge) for edge in edge_tuples]
    output["tree"]["edges"] = edge_list

    # Convert final output dictionary to JSON
    json_output = json.dumps(output)