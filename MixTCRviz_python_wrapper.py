import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser
import rpy2.robjects.packages as rpackages
import rpy2.robjects as ro
import pickle



# the R object is composed of nested lists 

def convert_matrix_to_dataframe(matrix):
    array = np.array(matrix)
    rownames = ro.r['rownames'](matrix)
    colnames = ro.r['colnames'](matrix)
    
    if rownames.rclass[0] == 'NULL':
        rownames = None
    if colnames.rclass[0] == 'NULL':
        colnames = None
        
    df = pd.DataFrame(array, index=rownames, columns=colnames)
    return df

def convert_named_vector_to_series(named_vector):
    names = ro.r['names'](named_vector)
    if names.rclass[0] != 'NULL':
        return pd.Series(named_vector, index=names)
    else:
        return list(named_vector)

def parse_element(element):
    if isinstance(element, ro.vectors.ListVector):
        return {k: parse_element(element.rx2(k)) for k in element.names}
    elif isinstance(element, ro.vectors.Matrix):
        return convert_matrix_to_dataframe(element)
    elif isinstance(element, ro.vectors.Vector):
        return convert_named_vector_to_series(element)
    else:
        return element

def parse_r_object(r_object):
    return {k: parse_element(r_object.rx2(k)) for k in r_object.names}

def save_data_to_pickle(file_path, data):
    with open(file_path, 'wb') as file:
        pickle.dump(data, file)

## run MixTCRviz function from the MixTCRviz package
#run ?MixTCRviz.MixTCRviz

if __name__ == '__main__':

    parser = ArgumentParser()
    #model params
    parser.add_argument('-i', '--input1', default = None)
    parser.add_argument('-o', '--output_path', default =None)
    parser.add_argument('--input2', default ="")
    parser.add_argument('--baseline_file', default ="")
    parser.add_argument('--use_allele', default =0)
    parser.add_argument('--correct_gene_names', default =1)
    parser.add_argument('--use_mouse_strain', default =0)
    parser.add_argument('--check_cdr3_mode', default =1)
    parser.add_argument('--renormVJ', default =1)
    parser.add_argument('--N_min', default =10)
    parser.add_argument('--output_stat', default =1)
    parser.add_argument('--set_cdr3a_length', default ="")
    parser.add_argument('--set_cdr3b_length', default ="")
    parser.add_argument('--species_default', default ="HomoSapiens")
    parser.add_argument('--model_default', default ="Model_default")
    parser.add_argument('--verbose', default =1)
    parser.add_argument('--plot', default =1)
    parser.add_argument('--plot_cdr12_motif', default =0)
    parser.add_argument('--plot_oneline', default =0)
    parser.add_argument('--plot_logo_length', default =0)
    parser.add_argument('--plot_cdr3_norm', default =0)
    parser.add_argument('--chain_list_output', default ="AB")
    parser.add_argument('--input1_name', default ="Input")
    parser.add_argument('--input2_name', default ="")
    parser.add_argument('--output_format', default ="pdf")
    args = parser.parse_args()

    #if args.help:
    #    print("####################################################################################################")
    #    print("####################################################################################################")
    #    print("Usage: python MixTCRviz_python.py -i inp")
    #    sys.exit(0)

    # for a quick test
    #args.input1 = './test/test1.csv'
    #args.output_path = './test/out/test1'
    #args.plot_cdr12_motif= 1

    MixTCRviz = rpackages.importr('MixTCRviz')

    #Python wrapper of MixTCRviz funcion.
    MixTCRviz.MixTCRviz(input1 = args.input1, output_path = args.output_path, input2 = args.input2, baseline_file = args.baseline_file, 
                        use_allele= args.use_allele, correct_gene_names= args.correct_gene_names, use_mouse_strain= args.use_mouse_strain, 
                        check_cdr3_mode=args.check_cdr3_mode, renormVJ=args.renormVJ, N_min=args.N_min, output_stat= args.output_stat, 
                        set_cdr3a_length= args.set_cdr3a_length, set_cdr3b_length= args.set_cdr3b_length, species_default= args.species_default, 
                        model_default= args.model_default, verbose= args.verbose, plot= args.plot, plot_cdr12_motif= args.plot_cdr12_motif, 
                        plot_oneline= args.plot_oneline, plot_logo_length= args.plot_logo_length, plot_cdr3_norm= args.plot_cdr3_norm, 
                        chain_list_output= args.chain_list_output, input1_name= args.input1_name, input2_name= args.input2_name, output_format= args.output_format)



    ### now convert the rds objects from the stat folder into python pickle dictionaries

    # Define the path to your RDS file
    path_stats = '{0}/stats'.format(args.output_path)
    for rds_file in os.listdir(path_stats):
        if rds_file.split(".")[-1] != 'rds':
            continue
        # Load the RDS file
        data = ro.r['readRDS'](os.path.join(path_stats, rds_file))
        data_python = parse_r_object(data)


        name_pickle_file = os.path.join(path_stats, rds_file.replace('.rds', '.pkl'))
        save_data_to_pickle(name_pickle_file, data_python)               




