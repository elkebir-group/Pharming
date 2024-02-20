from data import Data, load_from_pickle
import argparse



def write_input(dat, outfile):
    # outfile = "simulation_study/decifer_test/test.tsv"
    s_index = "0"
    s_label = "0"
    with open(outfile, "w+") as file:
        file.write('0 #characters\n')
        file.write('0 #samples\n')
        file.write('#sample_index\tsample_label\tcharacter_index\tcharacter_label\tref\tvar\n')

        for k, snvs in dat.seg_to_snvs.items():
            cn_props = dat.cn_proportions(k)
            prop_str = "\t".join([f"{s[0]}\t{s[1]}\t{val}" for s, val in cn_props.items()])
            # print(prop_str)
            for j in snvs:
                var  = dat.var[:, j].sum(axis=0)
                if var == 0:
                    continue
                total = dat.total[:, j].sum(axis=0)
                ref= total- var 

                file.write(f"{s_index}\t{s_label}\t{j}\t{j}\t{ref}\t{var}\t{prop_str}\n")



if __name__ == "__main__":

  
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=True,
                        help="input file of preprocessed data pickle")
    parser.add_argument("-o" ,"--out", required=True, type=str,
                        help="where the output file should be written")
    
    args = parser.parse_args()


    # instance = "s11_m5000_k25_l7"
    # # instance = "s12_m5000_k25_l7"
    # folder = "n1000_c0.25_e0" 
    # pth = f"simulation_study/input"

    # args = parser.parse_args([

    #     "-d", f"{pth}/{instance}/{folder}/data.pkl",
    #     "--out", f"test/dec.in",

    # ])

    dat = load_from_pickle(args.data)
    write_input(dat, args.out)