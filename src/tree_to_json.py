import csv
import json
# import networkx as nx
# import pickle
# from clonal_tree import ClonalTree

def map_mut_to_seg(seg_to_muts):
    mut_to_seg = {}
    mut_and_seg = []
    for seg in seg_to_muts:
        muts = seg_to_muts[seg]
        for mut in muts:
            mut_to_seg[mut] = seg
            combined = {
                "snv_id": mut,
                "segment_id": seg
            }
            mut_and_seg.append(combined)

    return mut_to_seg, mut_and_seg

def convertToJson(clonal_tree, segment_csv=None, snv_csv=None):

    # This will become the output JSON format of the tree
    output = { "segments":  [],
               "snvs":      [],
               "tree": {
                   "nodes": [],
                   "edges": []
               }
              }
    
    # output["segments"] = list(clonal_tree.seg_to_muts.keys())
    segments = list(clonal_tree.seg_to_muts.keys())
    for seg in segments:
        output["segments"].append(
            {
                "segment_id": seg
            }
        )

    mut_to_seg, mut_and_seg = map_mut_to_seg(clonal_tree.seg_to_muts)
    # print(mut_and_seg)
    output["snvs"] = mut_and_seg

    # Start by reading in the segment CSV and extracting the segment info
    if segment_csv != None:
        with open(segment_csv, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                segment_id, chromosome, start, stop = row
                if segment_id == "segment_id":
                    continue

                # Look for the segment in the list, then add its chromosome, start, and stop to its info
                for i in range(len(output["segments"])):
                    segment = output["segments"][i]
                    if segment["segment_id"] == int(segment_id):
                        segment["chromosome"] = int(chromosome)
                        segment["start"] = start
                        segment["stop"] = stop

    # Then read in the SNV CSV and extract the SNV info
    if snv_csv != None:
        with open(snv_csv, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                snv_id, segment_id, pos, ref, alt = row
                if snv_id == "snv_id":
                    continue

                # Look for the SNV in the list, then add its ids, pos, ref, and alt to its info
                for snv in output["snvs"]:
                    if snv["snv_id"] == int(snv_id) and snv["segment_id"] == int(segment_id):
                        snv["pos"] = int(pos)
                        snv["ref"] = ref
                        snv["alt"] = alt

    # Now get the nodes on the tree along with their metadata
    for node_id in clonal_tree.genotypes:
        segments = []
        snvs = []
        for snv_id in clonal_tree.genotypes[node_id]:

            # TODO: Add segment id corresponding to this SNV using the SNV dictionary
            segment_id = mut_to_seg[snv_id]

            x, y, x_bar, y_bar = clonal_tree.genotypes[node_id][snv_id].to_tuple()
            genotype = clonal_tree.genotypes[node_id][snv_id]
            # if genotype.x_bar == 1:
            #     print("found a 1")

            segments.append({"segment_id": int(segment_id), "x": int(x), "y": int(y)})
            snvs.append({"segment_id": int(segment_id), "snv_id": int(snv_id), "x_bar": int(x_bar), "y_bar": int(y_bar)})

        output["tree"]["nodes"].append({"node_id": node_id, "segments": segments, "snvs": snvs})

    # Get the edges on the tree
    # TODO: Check if this method works
    edge_list = list(clonal_tree.tree.edges())
    output["tree"]["edges"] = edge_list

    # Convert final output dictionary to JSON
    json_output = json.dumps(output)
    return json_output
    # with open("test.json", "w") as f:
    #     f.write(str(json_output))


# Testing on a tree
# with open('/Users/divya/Pharming/example/gt2.pickle', 'rb') as f:
#     data = pickle.load(f)
#     # print(type(data))
#     convertToJson(data, "/Users/divya/Pharming/example/segment_csv.csv", "/Users/divya/Pharming/example/snv_csv.csv")
    # print(data)