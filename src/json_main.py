import argparse
from clonal_tree import load


if __name__ == "__main__":

  
    parser = argparse.ArgumentParser()
    # parser.add_argument("-d", "--data", required=False,
    #                     help="input file of preprocessed data pickle")
    parser.add_argument("-t", "--tree", required=False,
                        help="input file of ground truth tree")
    parser.add_argument("-c", "--cell-assign", required=False,
                        help="input file of ground truth  cell assignment")
    parser.add_argument("-o", "--json", required=False,
                        help="output filename for json file")
    
    
    args = parser.parse_args()


    instance = "s14_m5000_k25_l7"
    # instance = "s12_m5000_k25_l7"
    folder = "n1000_c0.25_e0" 
    pth = f"simulation_study/input"

    args = parser.parse_args([
        "-t", f"{pth}/{instance}/gt.pkl",
        "-c", f"{pth}/{instance}/{folder}/cellAssign.pkl",
        "-o", f"/Users/leah/Documents/Research/projects/Pharming/demo_json/{instance}_{folder}.json"
    ])


    gt = load(args.tree)

    gt.toJson(args.json)


