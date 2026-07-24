"""tealeaf single-cell clustering pipeline.

Reads alevin-fry equivalence-class output, runs EM transcript quantification,
aggregates barcodes into pseudobulk samples, and applies the same intron
clustering logic as the bulk pipeline.
"""

import numpy as np
import pandas as pd
import sys
import warnings
import scipy
import scipy.sparse
from scipy.sparse import csr_matrix, save_npz, load_npz, vstack, coo_array
from optparse import OptionParser
import random
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import os
import re

warnings.simplefilter(action='ignore', category=pd.errors.DtypeWarning)

from tealeaf.utils import timing_decorator, write_options_to_file
from tealeaf.sc import sc_utils
from tealeaf.shared_functions import build_init_cluster, process_clusters, compute_ratio

from joblib import Parallel, delayed


def barcode_group_print(barcode_group_dic, out_prefix = ''):
    """
    print a file with barcodes and its corresponding group
    """
    with open(f'{out_prefix}barcodes_to_pseudobulk.txt', 'w') as output:
        for group in barcode_group_dic:
            for barcode in barcode_group_dic[group]:
                print(f'{barcode},{group}', file = output)



def metacell_generation(type_dic, n, out_prefix = ''):
    """
    this is the metacell generate function
    type_dic: a dict that contain type/cluster:[barcodes]
    n: number of cell merge to a metacell
    out_prefix: the prefix of the output file
    """
    n = int(n)
    metacell_dic = {}
    with open(f'{out_prefix}meta_cells_record.txt', 'w') as output, open(f'{out_prefix}meta_group.txt', 'w') as group_file:
        print('type/cluster num', file=output)
        for key in type_dic:
            barcodes = type_dic[key]
            random.shuffle(barcodes) # randomize the barcodes
            num_meta_cell = len(barcodes)//n  #integer division
            # if avalible cell < n, then 0 metacell will be generated
    
            print(f'{key} {num_meta_cell}', file=output)
            for i in range(num_meta_cell):
                print(f'{key}_{i} {key}', file = group_file)
                metacell_dic[f'{key}_{i}'] = barcodes[n*i: n*(i+1)]
    
    return metacell_dic
    

def bootstrapping(type_dic, n, k, out_prefix = ''):
    """
    this is the bootstrapping generate function
    type_dic: a dict that contain type/cluster:[barcodes]
    n: number of cell to merge
    k: round of bootstrapping for each cell type or cluster
    out_prefix: the prefix of the output file
    """
    n = int(n)
    bootstrap_dic = {}
    k = int(k)
    with open(f'{out_prefix}bootstrapping_record.txt', 'w') as output, open(f'{out_prefix}bootstrapping_group.txt', 'w') as group_file:
        print('type/cluster num', file=output)
        for key in type_dic:
            barcodes = type_dic[key]
            
            
            if len(barcodes) < n: # if avalible cell < n, then skip this cell type
                print(f'{key} 0', file=output)
                
            else:
                
                for i in range(k):
                    print(f'{key}_{i} {key}', file = group_file)
                    bootstrap_dic[f'{key}_{i}'] = random.choices(barcodes , k=n)

    return bootstrap_dic
    
    
    



def pseudo_group_generation(barcodes_type_file, n, k = 30, group_method = 'metacells', out_prefix = ''):
    '''
    barcodes_type_file: a csv file that contain two columns barcodes cell_type/cluster
    group_method: can be metacells or bootstrapping, otherwise retrun error message
    n: number of barcodes to be merged into a single pseudobulk sample
    k: number of bootstrapping round when generating pseudobulk samples, unused when group_method == metacells
    Returns
    -------
    None.

    '''
    df = pd.read_csv(barcodes_type_file, index_col = 0, header = None, names = ['type/cluster'])
    df = df.sort_values(by='type/cluster')
    # ensure cluster number will be interpreted as str
    df['type/cluster'] = df['type/cluster'].astype(str) 
    # replace the possible space in the cluster or cell type name. Avoid problem in Leafcutter
    # some downstream analysis in leafcutter use space to sepearte values
    df['type/cluster'] = df['type/cluster'].str.replace(" ", "_") 
    

    types = df['type/cluster'].unique()
    type_dic = {}
    for name in types:
        type_dic[name] = list(df[df['type/cluster'] == name].index)




    if group_method == 'metacells':
        merged_dic = metacell_generation(type_dic, n, out_prefix)
        barcode_group_print(merged_dic, f'{out_prefix}meta_')
    else: 
        merged_dic = bootstrapping(type_dic, n, k, out_prefix)
        barcode_group_print(merged_dic, f'{out_prefix}bootstrapping_')






def pseudo_dic_generation(barcode_pseudo_df):
    """
    This function return a dict that map pseudo sample to a list of barcodes in that sample

    """
    
    samples = barcode_pseudo_df['sample'].unique()
    sample_dic = {}
    for name in samples:
        sample_dic[name] = list(barcode_pseudo_df[barcode_pseudo_df['sample'] == name].barcode)
    return sample_dic

def pseudo_index_dic_generation(sample_dic, barcodes):
    """
    This function generate a index dict that all barcodes in the experiment
    then find the indexes of barcodes in a sample among all barcodes

    """
    sample_index_dic = {}
    barcodes_index_dic = {}
    for index, barcode in enumerate(barcodes):
        barcodes_index_dic[barcode] = index
    
    for key in sample_dic:
        sample_index_dic[key] = [barcodes_index_dic.get(barcode) for barcode in sample_dic[key]]
    return sample_index_dic


def process_batch(sparse_matrix, batch_row_groups):
    """
    This function used to process sample groups in batch
    it take a sparse_matrix for all barcodes, and using indexes for each sample
    it generate a pseudobulk row for each sample and stack them to a matrix
    """
    
    
    result_matrix =csr_matrix((0, sparse_matrix.shape[1]))
    for group_indices in batch_row_groups:
        # Sum rows efficiently and ensure the result is a 2-D matrix
        summed_row = sparse_matrix[group_indices].sum(axis=0)
        if summed_row.ndim == 1:
            summed_row = summed_row.reshape(1, -1)  # Reshape to 2-D if necessary
        result_matrix = vstack([result_matrix, summed_row]) 

    return result_matrix


def check_barcodes_exsitent(barcode_pseudo_df, valid_barcodes, out_prefix):
    '''
    This function check whether the barcodes in barcode_pseudo_df is in the valid barcodes list 
    Will write a file that contain the final number of barcodes for each pseudo sample and a list 
    of eliminated samples
    will return a new df that only contain the valid barcodes

    '''
    valid_barcodes = set(valid_barcodes)
    """
    filtered_index_list = [x for x in barcode_pseudo_df.index if x in valid_barcodes]
    filter_out_index_list = [x for x in barcode_pseudo_df.index if x not in valid_barcodes]
    filtered_df = barcode_pseudo_df.loc[filtered_index_list]
    """

    filtered_df = barcode_pseudo_df[barcode_pseudo_df['barcode'].isin(valid_barcodes)]

    filter_out_df = barcode_pseudo_df[~ barcode_pseudo_df['barcode'].isin(valid_barcodes)]


    with open(f'{out_prefix}eliminated_barcodes.txt', 'w') as file:
        for barcode in list(filter_out_df['barcode']):
            file.write(f"{barcode}\n")
            
            
    pseudo_count = filtered_df['sample'].value_counts()
    # sorting by sample names
    split_columns = pseudo_count.index.to_series().str.extract(r'(.*)_(\d+)$')
    split_columns[1] = split_columns[1].astype(int)
    sorted_indices = split_columns.sort_values(by=[0, 1]).index
    pseudo_count = pseudo_count.reindex(sorted_indices)



    with open(f'{out_prefix}pseudo_barcodes_counts.txt', 'w') as file:
        for pseudo, count in pseudo_count.items():
            file.write(f"{pseudo} {count}\n")


    sys.stderr.write(f"There are {len(filter_out_df)} barcodes get eliminated \n")
    
    return  filtered_df



def alpha_to_tpm_count(alpha, total_count, w, normalize_mode, read_length, paired_end,
                       overhang, sizing_factor, bool_spliced_trans):
    """Convert an isoform probability vector into tealeaf TPM/count outputs."""
    TPM = alpha * 1000000

    if normalize_mode == 'global':
        count = (alpha * total_count) * sizing_factor
    elif normalize_mode == 'snRNA':
        count = np.where(bool_spliced_trans == 1, alpha, 0)
        count = count / w
        count_sum = count.sum()
        if count_sum > 0:
            count = (count / count_sum) * total_count

        if paired_end:
            effective_read_len = (read_length - overhang) * 2 * sizing_factor
        else:
            effective_read_len = (read_length - overhang) * sizing_factor

        count = count * w * effective_read_len
    else:
        count = alpha / w
        count_sum = count.sum()
        if count_sum > 0:
            count = (count / count_sum) * total_count

        if paired_end:
            effective_read_len = (read_length - overhang) * 2 * sizing_factor
        else:
            effective_read_len = (read_length - overhang) * sizing_factor

        count = count * w * effective_read_len

    return TPM, count


def process_pseudo_row(i, cell_ec_sparse_pseudo_filt, ec_transcript_input, w, normalize_mode, read_length, paired_end , overhang, \
                       sizing_factor,bool_spliced_trans, quant_method='em', nnls_max_iter=200, nnls_tol=1e-6):
    """
    This is the helper function for pseudo_eq_conversion, this function will be use to enable multiprocessing
    This function will be called by parallel_EM_processing

    """
    temp_sample = cell_ec_sparse_pseudo_filt.getrow(i).toarray().ravel()
    total_count = temp_sample.sum()
    temp_sample += 0.0000001  # add 1e-6 to avoid null value

    if quant_method == 'em':
        alpha = sc_utils.EM(temp_sample, ec_transcript_input, w)
    elif quant_method == 'nnls':
        alpha = sc_utils.NNLS(temp_sample, ec_transcript_input, w,
                              max_iter=nnls_max_iter, tol=nnls_tol)
    else:
        raise ValueError(f"Unsupported quant_method for row processing: {quant_method}")

    # we don't need the actual count, we just need a scale that represents the support level of the TPM value
    # TPM * total_count gives a count-like value that retains the TPM ratio of transcripts
    TPM, count = alpha_to_tpm_count(
        alpha, total_count, w, normalize_mode, read_length, paired_end,
        overhang, sizing_factor, bool_spliced_trans,
    )
    return TPM, count



def parallel_EM_processing(cell_ec_sparse_pseudo_filt, ec_transcript_input, w, thread,normalize_mode, read_length, paired_end , overhang, sizing_factor,\
                           bool_spliced_trans):
    return parallel_quant_processing(
        cell_ec_sparse_pseudo_filt, ec_transcript_input, w, thread,
        normalize_mode, read_length, paired_end, overhang, sizing_factor,
        bool_spliced_trans, quant_method='em',
    )


def parallel_quant_processing(cell_ec_sparse_pseudo_filt, ec_transcript_input, w, thread, normalize_mode,
                              read_length, paired_end, overhang, sizing_factor, bool_spliced_trans,
                              quant_method='em', nnls_max_iter=200, nnls_tol=1e-6,
                              nucnorm_lambda=0.01, nucnorm_max_iter=50,
                              nucnorm_tol=1e-4, nucnorm_rank=50,
                              nucnorm_max_dense_entries=100000000,
                              glm_device='auto', glm_batch_cells=4096,
                              glm_rank=64, admm_rho=1.0, admm_inner_iter=25,
                              admm_adaptive_rho=True,
                              admm_rho_update_interval=10,
                              admm_rho_balance=10.0, admm_rho_scale=2.0,
                              nucnorm_tau=None, fw_nonnegative_penalty=1.0,
                              regularization_target='phi',
                              glm_design_input=None, ec_design='legacy'):
    if quant_method == 'nnls_nucnorm':
        alphas = sc_utils.NNLS_nucnorm(
            cell_ec_sparse_pseudo_filt,
            ec_transcript_input,
            w,
            regularization=nucnorm_lambda,
            max_iter=nucnorm_max_iter,
            tol=nucnorm_tol,
            svd_rank=nucnorm_rank,
            max_dense_entries=nucnorm_max_dense_entries,
            regularization_target=regularization_target,
            design_matrix=glm_design_input,
        )
        tpm_rows = []
        count_rows = []
        totals = np.asarray(cell_ec_sparse_pseudo_filt.sum(axis=1)).ravel()

        for i, alpha in enumerate(alphas):
            tmp_tpm, tmp_count = alpha_to_tpm_count(
                alpha, totals[i], w, normalize_mode, read_length, paired_end,
                overhang, sizing_factor, bool_spliced_trans,
            )
            tpm_rows.append(csr_matrix(tmp_tpm.reshape(1, -1)))
            count_rows.append(csr_matrix(tmp_count.reshape(1, -1)))
        return vstack(tpm_rows), vstack(count_rows)

    if quant_method in {
        'admm', 'admm_factorized', 'frank_wolfe',
        'frank_wolfe_penalized', 'factorized',
    }:
        from tealeaf.sc import glm_solvers

        compatibility = glm_design_input
        if compatibility is None:
            compatibility = sc_utils.weighted_ec_transcript_matrix(
                ec_transcript_input, w, parameterization=regularization_target
            )
        result = glm_solvers.fit_glm(
            cell_ec_sparse_pseudo_filt,
            compatibility,
            quant_method,
            rank=glm_rank,
            regularization=nucnorm_lambda,
            max_iter=nucnorm_max_iter,
            tol=nucnorm_tol,
            rho=admm_rho,
            adaptive_rho=admm_adaptive_rho,
            rho_update_interval=admm_rho_update_interval,
            rho_balance=admm_rho_balance,
            rho_scale=admm_rho_scale,
            inner_iter=admm_inner_iter,
            max_dense_entries=nucnorm_max_dense_entries,
            tau=nucnorm_tau,
            nonnegative_penalty=fw_nonnegative_penalty,
            device=glm_device,
            batch_cells=glm_batch_cells,
        )
        result.diagnostics['regularization_target'] = regularization_target
        result.diagnostics['ec_design'] = ec_design
        alpha_matrix = glm_solvers.result_to_csr(
            result, 0, cell_ec_sparse_pseudo_filt.shape[0], threshold=0.0
        ).toarray()
        tpm_rows = []
        count_rows = []
        totals = np.asarray(cell_ec_sparse_pseudo_filt.sum(axis=1)).ravel()
        for alpha, total_count in zip(alpha_matrix, totals):
            tmp_tpm, tmp_count = alpha_to_tpm_count(
                alpha, total_count, w, normalize_mode, read_length, paired_end,
                overhang, sizing_factor, bool_spliced_trans,
            )
            tpm_rows.append(csr_matrix(tmp_tpm.reshape(1, -1)))
            count_rows.append(csr_matrix(tmp_count.reshape(1, -1)))
        return vstack(tpm_rows), vstack(count_rows)

    results = Parallel(n_jobs=thread)(delayed(process_pseudo_row)(i,cell_ec_sparse_pseudo_filt, ec_transcript_input, w, \
                                                                  normalize_mode, read_length, paired_end , overhang, sizing_factor,bool_spliced_trans,\
                                                                  quant_method, nnls_max_iter, nnls_tol) \
                                  for i in range(cell_ec_sparse_pseudo_filt.shape[0]))

    tpm_rows = [csr_matrix(value[0].reshape(1, -1)) for value in results]
    count_rows = [csr_matrix(value[1].reshape(1, -1)) for value in results]
    return vstack(tpm_rows), vstack(count_rows)






def pseudo_eq_conversion(alevin_dir, salmon_ref, barcode_pseudo_file, min_EC = 5, min_transcript = 0, out_prefix= '', threshold = 0.1, thread =8, \
                         normalize_mode = 'junction', read_length = 100, paired_end = True, overhang = 2, sizing_factor = 1,\
                         quant_method = 'em', nnls_max_iter = 200, nnls_tol = 1e-6, nucnorm_lambda = 0.01,\
                         nucnorm_max_iter = 50, nucnorm_tol = 1e-4, nucnorm_rank = 50,\
                         nucnorm_max_dense_entries = 100000000, glm_device='auto',\
                         glm_batch_cells=4096, glm_rank=64, admm_rho=1.0,\
                         admm_inner_iter=25, admm_adaptive_rho=True,\
                         admm_rho_update_interval=10, admm_rho_balance=10.0,\
                         admm_rho_scale=2.0, nucnorm_tau=None,\
                         fw_nonnegative_penalty=1.0,\
                         regularization_target='phi', ec_design='legacy',\
                         eq_probabilities=None, eq_weight_cache=None):
    """
    

    Parameters
    ----------
    salmon_dir : the dir for where the salmon barcdeos eq matrix and eq transcript matrix are
    salmon_ref: the reference .fa file that salmon used for aligenment

    Returns
    -------
    None.

    """

    threshold = float(threshold) # ensure the threshold is a float number
    min_EC = float(min_EC)
    min_transcript = float(min_transcript)
    thread = int(thread)
    nnls_max_iter = int(nnls_max_iter)
    nnls_tol = float(nnls_tol)
    nucnorm_lambda = float(nucnorm_lambda)
    nucnorm_max_iter = int(nucnorm_max_iter)
    nucnorm_tol = float(nucnorm_tol)
    nucnorm_rank = int(nucnorm_rank)
    nucnorm_max_dense_entries = int(nucnorm_max_dense_entries)
    if regularization_target not in {'phi', 'theta'}:
        sys.exit("Error: regularization_target must be phi or theta.\n")
    if ec_design not in {'legacy', 'binary', 'weighted'}:
        sys.exit("Error: pseudobulk ec_design must be legacy, binary, or weighted.\n")
    if quant_method not in {
        'em', 'nnls', 'nnls_nucnorm', 'admm', 'admm_factorized',
        'frank_wolfe', 'frank_wolfe_penalized', 'factorized',
    }:
        sys.exit("Error: invalid quantification method...\n")
    # step 1: data loading
    transcript_lengths_dic = sc_utils.get_transcript_lengths(Path(salmon_ref))    
    
    
    
    barcodes_pseudo = pd.read_csv(barcode_pseudo_file, header = None, sep = ',',   names = ['barcode','sample'])
    alevin_dir = Path(alevin_dir)
    
    map_cache = alevin_dir / "gene_eqclass.npz"
    if map_cache.is_file(): 
        ec_transcript_mat = scipy.sparse.load_npz(map_cache)
    else:
        num_genes, num_ec, ecs = sc_utils.read_alevin_ec(alevin_dir/ "gene_eqclass.txt.gz") 
        ecs_list = [ ecs[k] for k in range(len(ecs)) ] # this works because indexing is dense
        ec_transcript_mat = sc_utils.to_coo(ecs_list)
        scipy.sparse.save_npz(map_cache, ec_transcript_mat)


    # Load cells x EC counts
    npz_cache = alevin_dir  / "geqc_counts.npz"

    if npz_cache.is_file(): 
        cell_ec_sparse = scipy.sparse.load_npz(npz_cache)

    else:
        cell_ec_sparse = scipy.io.mmread(alevin_dir / "geqc_counts.mtx") 
        scipy.sparse.save_npz(npz_cache, cell_ec_sparse)  
        
    cell_ec_sparse = cell_ec_sparse.tocsr()  # convert coo to csr format, necessary for row operation
    
    
    with open(alevin_dir / "quants_mat_cols.txt", 'r') as file:
        features = np.array([line.strip() for line in file.readlines()])
    with open(alevin_dir / "quants_mat_rows.txt", 'r') as file:
        barcodes = np.array([line.strip() for line in file.readlines()])

    glm_phi_design = None
    regularized_methods = {
        'nnls_nucnorm', 'admm', 'admm_factorized', 'frank_wolfe',
        'frank_wolfe_penalized', 'factorized',
    }
    if quant_method in regularized_methods:
        _, full_w = sc_utils.get_feature_weights(features, transcript_lengths_dic)
        probability_file = eq_probabilities or (alevin_dir / "gene_eqclass_probs.tsv.gz")
        cache_file = eq_weight_cache or (alevin_dir / "gene_eqclass_fixed_weights.npz")
        glm_phi_design = sc_utils.glm_design_matrix(
            ec_transcript_mat,
            full_w,
            parameterization='phi',
            design=ec_design,
            probability_file=probability_file,
            cache_file=cache_file,
        )
    
    
    # step 2: first round of filtering, to facilitate computation speed
    #aggregate all barcodes to give a single pseudobulk sample
    pseudobulk = sc_utils.sparse_sum(cell_ec_sparse,0) 
    
    transcript_count = pseudobulk @ ec_transcript_mat
    ec_sizes = sc_utils.sparse_sum(ec_transcript_mat,1)
    
    
    # first round of filtering, this will speed up the pseudobulk matrix generation
    ECs_to_keep = pseudobulk >= min_EC
    ec_transcript_filt = ec_transcript_mat.tocsr()[ECs_to_keep,:]
    cell_ec_filt = cell_ec_sparse[:,ECs_to_keep]
    features_to_keep = transcript_count > min_transcript
    features_filt = features[features_to_keep]
    ec_transcript_ff = ec_transcript_filt[:,features_to_keep]
    if glm_phi_design is not None:
        glm_phi_design = glm_phi_design[ECs_to_keep, :][:, features_to_keep]
    
    
    sys.stderr.write("start barcodes checking\n")
    
    barcodes_pseudo = check_barcodes_exsitent(barcodes_pseudo,barcodes, out_prefix= f'{out_prefix}')
    pseudo_dict = pseudo_dic_generation(barcodes_pseudo)
    # get a dict that contain information about pseudo_sample: [barcode1,barcode2,.......]
    pseudo_index_dic = pseudo_index_dic_generation(pseudo_dict, barcodes)
    
    



    # step 3: generate the pseudo eq matrix by merging barcodes
    batch_size = 200 # process in batch to avoid crash
    row_group_values = list(pseudo_index_dic.values())
    row_names = list(pseudo_index_dic.keys())
    batches = [row_group_values[i:i + batch_size] for i in range(0, len(row_group_values), batch_size)]

    # Process each batch and incrementally build the new sparse matrix
    pseudo_ec_sparse = csr_matrix((0, cell_ec_filt.shape[1]))  # Initialize an empty sparse matrix
    for batch in batches:
        batch_result = process_batch(cell_ec_filt, batch)
        pseudo_ec_sparse = vstack([pseudo_ec_sparse, batch_result])




    pseudo_ec_sparse = pseudo_ec_sparse.tocsr() 
    
    # step 4: additional filtering to further reduce the matrix size
    aggragate_cell_ec = sc_utils.sparse_sum(pseudo_ec_sparse,0)
    transcript_count = aggragate_cell_ec @ ec_transcript_ff
    

    ECs_to_keep = aggragate_cell_ec >= min_EC
    ec_transcript_fff = ec_transcript_ff[ECs_to_keep,:]
    cell_ec_sparse_pseudo_filt = pseudo_ec_sparse.tocsr()[:,ECs_to_keep]
    
    features_to_keep = transcript_count > min_transcript
    features_ff = features_filt[features_to_keep]
    ec_transcript_ffff = ec_transcript_fff[:,features_to_keep]
    if glm_phi_design is not None:
        glm_phi_design = glm_phi_design[ECs_to_keep, :][:, features_to_keep]


    feature_lengths, w = sc_utils.get_feature_weights(features_ff, transcript_lengths_dic)
    glm_design_input = None
    if glm_phi_design is not None:
        glm_design_input = sc_utils.parameterize_glm_design(
            glm_phi_design,
            w,
            regularization_target,
            normalize_columns=ec_design in {'binary', 'weighted'},
        )
    
    
    bool_spliced_trans = np.array([0 if '-I' in value else 1 for value in features_ff]) # Knowing what isofrom is spliced and which is unspliced


    # step 5: compute pseudo transcript matrix:
    ec_transcript_input = ec_transcript_ffff.tocoo()
    
    
    pseudo_TPM, pseudo_count = parallel_quant_processing(
        cell_ec_sparse_pseudo_filt, ec_transcript_input, w, thread,
        normalize_mode, read_length, paired_end, overhang, sizing_factor,
        bool_spliced_trans, quant_method, nnls_max_iter, nnls_tol,
        nucnorm_lambda, nucnorm_max_iter, nucnorm_tol, nucnorm_rank,
        nucnorm_max_dense_entries, glm_device, glm_batch_cells, glm_rank,
        admm_rho, admm_inner_iter, admm_adaptive_rho,
        admm_rho_update_interval, admm_rho_balance, admm_rho_scale,
        nucnorm_tau, fw_nonnegative_penalty,
        regularization_target,
        glm_design_input, ec_design,
    )
       

    # this will ensure that the extremly small value won't be detected in the final result
    pseudo_TPM.data[pseudo_TPM.data < 0.1] = 0 
    pseudo_TPM.eliminate_zeros()
    
    pseudo_count.data[pseudo_count.data < 0.1] = 0 
    pseudo_count.eliminate_zeros()
   
    

    # step 6: store matrices and eliminate unspliced isoform as they doesn't excised intron    
    pseudo_names = row_names
    
    spliced_indices = [i for i, s in enumerate(features_ff) if "-I" not in s]    
    spliced_feature = [s for i,s in enumerate(features_ff) if "-I" not in s]
    unspliced_indices = [i for i, s in enumerate(features_ff) if "-I" in s]

    spliced_pseudo_TPM = pseudo_TPM.tocsr()[:,spliced_indices]
    spliced_pseudo_count = pseudo_count.tocsr()[:,spliced_indices]


    with open(f'{out_prefix}pseudo_cols.txt', 'w') as output:
        for trans in features_ff:
            print(trans, file=output)
    with open(f'{out_prefix}pseudo_rows.txt', 'w') as output:
        for pseudo in pseudo_names:
            print(pseudo, file=output)
    with open(f'{out_prefix}pseudo_spliced_cols.txt', 'w') as output:
        for trans in spliced_feature:
            print(trans, file=output)
            
    
    

    save_npz(f'{out_prefix}pseudo_TPM.npz', pseudo_TPM)
    save_npz(f'{out_prefix}pseudo_count.npz', pseudo_count)
    save_npz(f'{out_prefix}pseudo_spliced_TPM.npz', spliced_pseudo_TPM)
    save_npz(f'{out_prefix}pseudo_spliced_count.npz', spliced_pseudo_count)


def single_cell_glm_conversion(options):
    """Fit a genome-wide GLM directly to raw cells and write sparse cell chunks.

    The downstream intron and clustering stages expect a manageable pseudobulk
    matrix.  They are intentionally not run in ``single_cell`` mode.
    """
    from tealeaf.sc import glm_cv, glm_solvers

    if options.quant_method not in glm_solvers.SCALABLE_METHODS:
        raise ValueError(
            "--cell_mode single_cell requires a scalable GLM quantification method"
        )
    if options.primer_pairs:
        prepared = glm_cv.prepare_paired_primer_glm_data(
            options.alevin_dir,
            options.salmon_ref,
            options.primer_pairs,
            ec_design=options.ec_design,
            regularization_target=options.regularization_target,
            min_eq=options.min_eq,
            min_half_umis=options.min_half_umis,
            primer_sampling_model=options.primer_sampling_model,
            probability_file=options.eq_probabilities,
        )
    else:
        prepared = glm_cv.prepare_alevin_glm_data(
            options.alevin_dir,
            options.salmon_ref,
            ec_design=options.ec_design,
            regularization_target=options.regularization_target,
            min_eq=options.min_eq,
            probability_file=options.eq_probabilities,
            weight_cache=options.eq_weight_cache,
        )
    selected = glm_cv.sample_cells_by_count(
        prepared.counts,
        0,
        min_count=options.min_cell_umis,
        totals=prepared.cell_umi_totals,
    )
    if not len(selected):
        raise ValueError("no cells meet --min_cell_umis")
    counts = prepared.counts[selected].tocsr()
    barcodes = prepared.barcodes[selected]
    initial_factors = None
    if options.glm_initial_factors:
        factor_path = Path(options.glm_initial_factors)
        suffix = "glm_factors.npz"
        if not factor_path.name.endswith(suffix):
            raise ValueError("--glm_initial_factors must end in glm_factors.npz")
        prefix = factor_path.parent / factor_path.name[:-len(suffix)]
        initial_rows = np.loadtxt(f"{prefix}glm_rows.txt", dtype=str)
        initial_cols = np.loadtxt(f"{prefix}glm_cols.txt", dtype=str)
        if not np.array_equal(initial_rows, barcodes):
            raise ValueError("warm-start cell rows do not match prepared GLM rows")
        if not np.array_equal(initial_cols, prepared.features):
            raise ValueError("warm-start transcript columns do not match prepared GLM columns")
        with np.load(factor_path) as factors:
            initial_factors = (factors["left"], factors["right"])
    result = glm_solvers.fit_glm(
        counts,
        prepared.compatibility,
        options.quant_method,
        rank=options.glm_rank,
        regularization=options.nucnorm_lambda,
        rho=options.admm_rho,
        adaptive_rho=options.admm_adaptive_rho,
        rho_update_interval=options.admm_rho_update_interval,
        rho_balance=options.admm_rho_balance,
        rho_scale=options.admm_rho_scale,
        max_iter=options.nucnorm_max_iter,
        min_iter=options.nucnorm_min_iter,
        patience=options.nucnorm_patience,
        inner_iter=options.admm_inner_iter,
        tol=options.nucnorm_tol,
        max_dense_entries=options.nucnorm_max_dense_entries,
        tau=options.nucnorm_tau,
        max_atoms=options.fw_max_atoms,
        power_iter=options.fw_power_iter,
        nonnegative_penalty=options.fw_nonnegative_penalty,
        device=options.glm_device,
        batch_cells=options.glm_batch_cells,
        data_backend=options.glm_data_backend,
        polish_max_iter=options.glm_polish_max_iter,
        minibatch=options.glm_minibatch,
        exact_inner_steps=options.glm_exact_inner_steps,
        initial_factors=initial_factors,
    )
    result.diagnostics['regularization_target'] = options.regularization_target
    result.diagnostics['ec_design'] = options.ec_design
    if prepared.metadata:
        result.diagnostics.update({
            key: value for key, value in prepared.metadata.items()
            if key not in {"half_umi_totals", "source_rows"}
        })
    glm_solvers.write_chunked_result(
        result,
        options.outprefix,
        barcodes,
        prepared.features,
        batch_cells=options.glm_batch_cells,
        threshold=options.glm_output_threshold,
        write_chunks=options.glm_write_chunks,
    )
    sys.stderr.write("Genome-wide single-cell GLM outputs were written\n")


def extract_order(col):
    '''
    This is the helper function that ensure name in row will be sorted correctly

    '''
    match = re.match(r"(.*)_(\d+)$", col)
    if match:
        prefix = match.group(1)  # Everything before the last underscore
        number = int(match.group(2))  # The numeric suffix after the last underscore
        return (prefix, number)
    return (col, 0)  # Default return for non-matching patterns



def sc_intron_count(reference_directory, pseudo_matrix_file, pseudo_cols_file, pseudo_rows_file, in_prefix = '' , out_prefix = '', threshold = 10):
    """
    

    Parameters
    ----------
    reference_directory : Str
        The directory that contain the reference
    pseudo_matrix_file : Str
        
    pseudo_cols_file : TYPE
        DESCRIPTION.
    pseudo_rows_file : TYPE
        DESCRIPTION.
    in_prefix : TYPE, optional
        DESCRIPTION. The default is ''.
    out_prefix : TYPE, optional
        DESCRIPTION. The default is ''.

    Returns
    -------
    None.

    """
    
    #loading reference matrix
    intron_isoform_matrix = load_npz(f'{reference_directory}/{in_prefix}isoform_intron_matrix.npz').tocsr()
    exon_isoform_matrix = load_npz(f'{reference_directory}/{in_prefix}isoform_exon_matrix.npz').tocsr()
    
    isoform_rows = []
    with open(f'{reference_directory}/{in_prefix}isoform_rows.txt', 'r') as input_file:
        for line in input_file:
            isoform_rows.append(line.split('.')[0]) # get rid of the version number
            
            
    intron_cols = []
    with open(f'{reference_directory}/{in_prefix}intron_cols.txt', 'r') as input_file:
        for line in input_file:
            intron_cols.append(line[:-1]) # get rid of the \n


    exon_cols = []
    with open(f'{reference_directory}/{in_prefix}exon_cols.txt', 'r') as input_file:
        for line in input_file:
            exon_cols.append(line[:-1]) 



    pseudo_matrix = load_npz(pseudo_matrix_file)
    
    
    pseudo_cols = [] #isoforms
    with open(pseudo_cols_file, 'r') as input_file:
        for line in input_file:
            pseudo_cols.append(line.split('.')[0]) 

    
    


    pseudo_rows = [] #barcodes
    with open(pseudo_rows_file, 'r') as input_file:
        for line in input_file:
            pseudo_rows.append(line[:-1]) 





    common_isoforms = set(isoform_rows).intersection(pseudo_cols) # get the common isoforoms

    # get the index for the isoforms 
    reference_index = {isoform: idx for idx, isoform in enumerate(isoform_rows)}




    # Step 3: Rearrange Matrix 1
    # Filter and reorder columns to align with Matrix 2
    filtered_indices = [i for i, isoform in enumerate(pseudo_cols) if isoform in common_isoforms]
    filtered_psudo_cols = [isoform for i, isoform in enumerate(pseudo_cols) if isoform in common_isoforms]
    filtered_pseudo_matrix = pseudo_matrix[:, filtered_indices]


    # Create a new empty matrix with the same number of rows and the correct number of columns
    new_pseudo_matrix = scipy.sparse.lil_matrix((pseudo_matrix.shape[0], len(reference_index)))


    for idx, isoform in enumerate(filtered_psudo_cols):
        if isoform in reference_index:
            col_new_index = reference_index[isoform]
            new_pseudo_matrix[:, col_new_index] = filtered_pseudo_matrix[:, idx]



    intron_count_matrix = new_pseudo_matrix @ intron_isoform_matrix 
    exon_count_matrix = new_pseudo_matrix @ exon_isoform_matrix 



    filtered_columns = intron_count_matrix.sum(axis=0) > threshold  # Filter out lowerly expressed intron
    # Step 2: Filter the matrix to keep only non-empty columns
    filtered_intron_count_matrix = intron_count_matrix[:, filtered_columns.A.ravel()]
    # Step 3: Filter the column names list
    filtered_intron_cols = [name for name, keep in zip(intron_cols, filtered_columns.A.ravel()) if keep]



    filtered_columns = exon_count_matrix.sum(axis=0) > threshold  # Filter out lowerly expressed exon
    # Step 2: Filter the matrix to keep only non-empty columns
    filtered_exon_count_matrix = exon_count_matrix[:, filtered_columns.A.ravel()]
    # Step 3: Filter the column names list
    filtered_exon_cols = [name for name, keep in zip(exon_cols, filtered_columns.A.ravel()) if keep]



    dense_intron = filtered_intron_count_matrix.toarray()
    df_intron = pd.DataFrame(dense_intron).T
    df_intron.index = filtered_intron_cols
    df_intron.columns = pseudo_rows
    
    

    

    dense_exon = filtered_exon_count_matrix.toarray()
    df_exon = pd.DataFrame(dense_exon).T
    df_exon.index = filtered_exon_cols
    df_exon.columns = pseudo_rows
    #df_exon = df_exon.sort_index(axis=1) # sort the columns
    
    
    sorted_columns = sorted(df_intron.columns, key=extract_order)
    df_intron = df_intron[sorted_columns] # sort the columns
    df_exon = df_exon[sorted_columns]
    



    df_intron[['Chr', 'Start', 'End']] = df_intron.index.to_series().str.extract(r'(chr\w+):(\d+)-(\d+)')
    new_order = list(df_intron.columns)[-3:] + list(df_intron.columns)[:-3] 
    df_intron = df_intron[new_order]
    
    
    df_exon[['Chr', 'Start', 'End']] = df_exon.index.to_series().str.extract(r'(chr\w+):(\d+)-(\d+)')
    
    new_order = list(df_exon.columns)[-3:] + list(df_exon.columns)[:-3] 
    df_exon = df_exon[new_order]
    
    df_intron.index.name = 'Name'
    df_exon.index.name = 'Name'

    df_intron.sort_values(by=['Chr', 'Start', 'End'], inplace=True)
    df_exon.sort_values(by=['Chr', 'Start', 'End'], inplace=True)
    df_intron.to_csv(f'{out_prefix}count_intron', sep = ' ')
    df_exon.to_csv(f'{out_prefix}count_exon', sep = ' ')    



@timing_decorator
def tealeaf_sc(options):
    """
    This is the main function for tealeaf_sc

    """
    out_prefix = options.outprefix
    if options.cell_mode == 'single_cell':
        single_cell_glm_conversion(options)
        return
    if options.preprocessed == False:
    #  step 1: generate or read the barcodes to pseudobulk samples file
        if options.pseudobulk_samples == None:
    
            pseudo_group_generation(options.barcodes_clusters, options.num_cell, options.num_bootstrapping, \
                                options.pseudobulk_method, options.outprefix)
            if options.pseudobulk_method == 'metacells':
                pseudobulk_name = f'{out_prefix}meta_barcodes_to_pseudobulk.txt'
            else:
                pseudobulk_name = f'{out_prefix}bootstrapping_barcodes_to_pseudobulk.txt'
            
            sys.stderr.write("Barcodes to pseudobulk samples generated\n")

        else:
            pseudobulk_name = options.pseudobulk_samples 
            sys.stderr.write("Barcodes to pseudobulk samples read in\n")
    


        # step 2: compute the pseudobulk isoform matrix 
    
        paired_end = (not options.not_paired_end)
        min_transcript = 0
        
        pseudo_eq_conversion(options.alevin_dir, options.salmon_ref, pseudobulk_name, options.min_eq, min_transcript , \
                             out_prefix, options.samplecutoff,  options.thread, options.normalization_scale,
                             options.read_length, paired_end, options.overhang, options.sizing_factor,
                             options.quant_method, options.nnls_max_iter, options.nnls_tol,
                             options.nucnorm_lambda, options.nucnorm_max_iter,
                             options.nucnorm_tol, options.nucnorm_rank,
                             options.nucnorm_max_dense_entries, options.glm_device,
                             options.glm_batch_cells, options.glm_rank,
                             options.admm_rho, options.admm_inner_iter,
                             options.admm_adaptive_rho,
                             options.admm_rho_update_interval,
                             options.admm_rho_balance, options.admm_rho_scale,
                             options.nucnorm_tau, options.fw_nonnegative_penalty,
                             options.regularization_target,
                             options.ec_design, options.eq_probabilities,
                             options.eq_weight_cache)
    
    
    
    
    
        sys.stderr.write("pseudobulk eq matrix were computed\n")
    
    # step 3: counting intron and exon
    ref_prefix = f'{options.ref_dir}/{options.ref_prefix}'
        

    if options.use_TPM == False:
        pseudo_matrix_file = f'{out_prefix}pseudo_spliced_count.npz'
    else:
        pseudo_matrix_file = f'{out_prefix}pseudo_spliced_TPM.npz'


    pseudo_cols_file = f'{out_prefix}pseudo_spliced_cols.txt'
    pseudo_rows_file = f'{out_prefix}pseudo_rows.txt'

    sc_intron_count(options.ref_dir, pseudo_matrix_file,pseudo_cols_file, pseudo_rows_file,  in_prefix= options.ref_prefix ,\
                    out_prefix = out_prefix, threshold = options.introncutoff)
    sys.stderr.write("Intron and exon were quantified\n")
    # a TSV file compatible with the rest of the tealeaf pipeline is obtained

    # step 4: clustering, using shared functions from tealeaf

    if options.with_virtual == False:
        connect_file = f'{ref_prefix}intron_exon_connectivity.tsv'
    else:
        connect_file = f'{ref_prefix}intron_exon_connectivity_with_virtual.tsv'



    build_init_cluster(f'{out_prefix}count_intron', connect_file)
    
    sys.stderr.write("Finished Initial Clustering\n")



    process_clusters(f'{out_prefix}count_intron', f'{out_prefix}count_exon', \
                     connect_file, \
                     out_prefix = out_prefix,\
                     cutoff = options.introncutoff, \
                     percent_cutoff = options.mincluratio,\
                     min_cluster_val = options.minclucounts)

    
    
    sys.stderr.write("Finished Cluster refinement\n")
    
    compute_ratio(f'{out_prefix}refined_cluster', out_prefix)

    sys.stderr.write("Finished PSI calculation\n")


    sys.stderr.write('CLustering Fnished\n')
    







if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()


    
    # necessary parameter
    parser.add_option("--alevin_dir", dest="alevin_dir", default = None,
                  help="The directory for alevin results, should contain the eq matrix and other files (default: None)")
    
    
    parser.add_option("--salmon_ref", dest="salmon_ref", default = None,
                  help="The reference used for salmon index, The salmon reference,  maybe spliceu or splicei (default: None)")

    parser.add_option("--ref_dir", dest="ref_dir", default=None,
                      help="tealeaf reference directory containing isoform-to-intron and isoform-to-exon "
                           "matrices produced by tealeaf-map (default: None)")

    parser.add_option("--barcodes_clusters", dest="barcodes_clusters", default = None,
                  help="The file that records which barcodes belong to which cluster/cell type in the format 'barcode,cluster' \
                        this file will be used to generate pseudobulk samples")
                  
    parser.add_option('--pseudobulk_samples', dest="pseudobulk_samples", default = None,
                      help="a txt file with barcodes to pseudobulk sample are expected in format 'barcode pseudobulk_ample', if \
                      this option!= None, then it will overwrite the input to --barcodes_cluster, and use the file in this option \
                      for computation (default: None)")
                      

    # optional parameter

    parser.add_option("--ref_prefix", dest="ref_prefix", default='',
                      help="file prefix used when generating the isoform-to-intron map with "
                           "tealeaf-map (default: '')")            
               
    parser.add_option("-n",'--num_cell', dest="num_cell", default = 100, type="int",
                  help="the number of cell/barcode that would like to include in a pseudobulk sample, cluster/cell type that have fewer cell than\
                      this number will not included in the computation (default: 100)")
    
    parser.add_option("-k",'--num_bootstrapping', dest="num_bootstrapping", default = 30, type="int",
                  help="the number of bootstrapping samples generated for each cluster/cell type (default: 30)")
    parser.add_option("--min_eq", dest="min_eq", default = 5,
                  help="minimum count for each eq class for it to be included in the EM (default: 5)")


    parser.add_option('--pseudobulk_method', dest="pseudobulk_method", default = 'metacells',
                  help="the pseudobulk sample generate method, could be metacells or bootstrapping (default: metacells)")

    parser.add_option("--thread", dest="thread", default = 8, type="int",
                  help="the thread use for computation")            

    parser.add_option("--quant_method", dest="quant_method", default='em',
                  help="EC quantification method, including scalable frank_wolfe_penalized (default: em)")

    parser.add_option("--cell_mode", dest="cell_mode", default='pseudobulk',
                  help="fit pseudobulks or raw cells: pseudobulk or single_cell (default: pseudobulk)")

    parser.add_option("--min_cell_umis", dest="min_cell_umis", default=0, type="float",
                  help="minimum raw UMI count for single-cell GLM fitting (default: 0)")

    parser.add_option("--primer_pairs", dest="primer_pairs", default=None,
                  help="TSV pairing polydT and ranhex barcodes for a shared-cell GLM")

    parser.add_option("--min_half_umis", dest="min_half_umis", default=500, type="float",
                  help="minimum UMI count required in each paired primer half (default: 500)")

    parser.add_option("--primer_sampling_model", dest="primer_sampling_model",
                  default="effective_length",
                  help="paired-primer sampling: effective_length, oligodt_tpm, or all_tpm")

    parser.add_option("--glm_device", dest="glm_device", default='auto',
                  help="Torch device for scalable GLMs: auto, cuda, or cpu (default: auto)")

    parser.add_option("--glm_batch_cells", dest="glm_batch_cells", default=4096, type="int",
                  help="number of cells per Torch GLM optimization block (default: 4096)")

    parser.add_option("--glm_data_backend", dest="glm_data_backend", default='auto',
                  help="cell-block storage: auto, cuda, or pinned (default: auto)")

    parser.add_option("--glm_polish_max_iter", dest="glm_polish_max_iter", default=32, type="int",
                  help="maximum deterministic polishing iterations for direct factorization (default: 32)")

    parser.add_option("--glm_exact_inner_steps", dest="glm_exact_inner_steps", default=32, type="int",
                  help="FISTA steps per factor in each exact factorization epoch (default: 32)")

    parser.add_option("--glm_minibatch", dest="glm_minibatch", default=False,
                  action="store_true",
                  help="use stochastic cell-minibatch epochs before deterministic polishing")

    parser.add_option("--glm_no_minibatch", dest="glm_minibatch",
                  action="store_false",
                  help="use deterministic accelerated factorization epochs only (default)")

    parser.add_option("--glm_initial_factors", dest="glm_initial_factors", default=None,
                  help="existing glm_factors.npz to validate and continue")

    parser.add_option("--glm_rank", dest="glm_rank", default=64, type="int",
                  help="low-rank capacity for factorized and Frank-Wolfe GLMs (default: 64)")

    parser.add_option("--glm_output_threshold", dest="glm_output_threshold", default=1e-8, type="float",
                  help="drop smaller normalized values from single-cell sparse output (default: 1e-8)")

    parser.add_option("--glm_no_write_chunks", dest="glm_write_chunks", default=True,
                  action="store_false",
                  help="write compact factors without reconstructing cell-by-transcript chunks")

    parser.add_option("--nnls_max_iter", dest="nnls_max_iter", default=200, type="int",
                  help="maximum iterations for scipy bounded least-squares NNLS (default: 200)")

    parser.add_option("--nnls_tol", dest="nnls_tol", default=1e-6, type="float",
                  help="tolerance for scipy bounded least-squares NNLS (default: 1e-6)")

    parser.add_option("--nucnorm_lambda", dest="nucnorm_lambda", default=0.01, type="float",
                  help="nuclear-norm penalty for quant_method=nnls_nucnorm (default: 0.01)")

    parser.add_option("--regularization_target", dest="regularization_target", default='phi',
                  help="coefficient matrix receiving low-rank regularization: phi or theta (default: phi)")

    parser.add_option("--ec_design", dest="ec_design", default='legacy',
                  help="fixed EC design: legacy, binary, weighted, or paired positional (default: legacy)")

    parser.add_option("--eq_probabilities", dest="eq_probabilities", default=None,
                  help="alevin-fry per-UMI probability sidecar for --ec_design weighted")

    parser.add_option("--eq_weight_cache", dest="eq_weight_cache", default=None,
                  help="cache path for the fixed weighted EC design")

    parser.add_option("--nucnorm_max_iter", dest="nucnorm_max_iter", default=50, type="int",
                  help="maximum proximal-gradient iterations for quant_method=nnls_nucnorm (default: 50)")

    parser.add_option("--nucnorm_min_iter", dest="nucnorm_min_iter", default=10, type="int",
                  help="minimum iterations before scalable GLM convergence checks (default: 10)")

    parser.add_option("--nucnorm_patience", dest="nucnorm_patience", default=5, type="int",
                  help="consecutive convergence checks required for scalable GLMs (default: 5)")

    parser.add_option("--nucnorm_tol", dest="nucnorm_tol", default=1e-4, type="float",
                  help="relative convergence tolerance for quant_method=nnls_nucnorm (default: 1e-4)")

    parser.add_option("--nucnorm_rank", dest="nucnorm_rank", default=50, type="int",
                  help="number of singular vectors used by truncated SVT for quant_method=nnls_nucnorm (default: 50)")

    parser.add_option("--nucnorm_max_dense_entries", dest="nucnorm_max_dense_entries", default=100000000, type="int",
                  help="maximum dense Phi entries allowed for quant_method=nnls_nucnorm (default: 100000000)")

    parser.add_option("--nucnorm_tau", dest="nucnorm_tau", default=None, type="float",
                  help="nuclear-norm radius for frank_wolfe; default is sqrt(number of cells)")

    parser.add_option("--fw_max_atoms", dest="fw_max_atoms", default=None, type="int",
                  help="maximum stored Frank-Wolfe atoms; defaults to glm_rank")

    parser.add_option("--fw_power_iter", dest="fw_power_iter", default=3, type="int",
                  help="power iterations per Frank-Wolfe linear oracle (default: 3)")

    parser.add_option("--fw_nonnegative_penalty", dest="fw_nonnegative_penalty", default=1.0, type="float",
                  help="negative-mass penalty relative to squared design spectral norm for frank_wolfe_penalized (default: 1)")

    parser.add_option("--admm_rho", dest="admm_rho", default=1.0, type="float",
                  help="ADMM penalty for admm and admm_factorized (default: 1.0)")

    parser.add_option("--admm_fixed_rho", dest="admm_adaptive_rho", default=True,
                  action="store_false", help="disable residual-balanced ADMM rho updates")

    parser.add_option("--admm_rho_update_interval", dest="admm_rho_update_interval", default=10, type="int",
                  help="iterations between adaptive ADMM rho updates (default: 10)")

    parser.add_option("--admm_rho_balance", dest="admm_rho_balance", default=10.0, type="float",
                  help="primal/dual residual imbalance required to change rho (default: 10)")

    parser.add_option("--admm_rho_scale", dest="admm_rho_scale", default=2.0, type="float",
                  help="multiplicative adaptive ADMM rho update (default: 2)")

    parser.add_option("--admm_inner_iter", dest="admm_inner_iter", default=25, type="int",
                  help="projected-gradient inner iterations for dense admm (default: 25)")
                  
    parser.add_option("--use_TPM", dest="use_TPM", default = False, action="store_true",
                  help="whether to performance normalization, if not use TPM directly (default: True)")
    
    
    parser.add_option("--preprocessed", dest="preprocessed", default = False, action="store_true",
                  help="Whether pseudobulk generation and EM were done, if true, \
                      then the pipeline start from counting intron (default: False)")

    parser.add_option("-v","--with_virtual", dest="with_virtual", action="store_true", default = False,
                  help="Whether the map that contain virtual intron to capture AFE and ALE(default: False)")
    
    
    parser.add_option("--cluster_def", dest="cluster_def", default = 3, type="int",
                  help="three def available, 1: overlap, 2: overlap+share_intron_splice_site, \
                      3: overlap+share_intron_splice_site+shared_exon_splice_site")
    
    parser.add_option("-o", "--outprefix", dest="outprefix", default='tealeaf_',
                      help="output file prefix; include directory path if not the current directory "
                           "(default: tealeaf_)")    

    parser.add_option("--samplecutoff", dest="samplecutoff", default = 0.1, type="float",
                  help="minimum count for an isoform in a sample to count as exist(default: 0.1)")

    parser.add_option("--introncutoff", dest="introncutoff", default = 80,type="float",
                  help="minimum count for an intron to count as exist(default 5)")
    
    parser.add_option("-m", "--minclucounts", dest="minclucounts", default = 100,type="float",
                  help="minimum count in a cluster (default 30 normalized count)")
    
    parser.add_option("-r", "--mincluratio", dest="mincluratio", default = 0.01,type="float",
                  help="minimum fraction of reads in a cluster that support a junction (default 0.01)")


    parser.add_option("--normalization_scale", dest="normalization_scale", default = 'junction',
                  help="The mode use for normaliztion, whether the count/TPM scale is based on junction count simulation,\
                      snRNA (attribute unspliced to spliced), and global. Can only input junction, snRNA, or global (default: junction)")     
    
    parser.add_option("--read_len", dest="read_length", default = 100, type = "int",
                  help="The read length of sequencing data, use to simulate junction count, only work when normalization_scale = junction\
                  (default: 100)")
                  
    parser.add_option("--overhang", dest="overhang", default = 2, type = "int",
                      help="The oeverhand that would like to use, could be set to zero, use to simulate junction count, \
                      only work when normalization_scale = junction (default: 2)")

    parser.add_option("--sizing_factor", dest="sizing_factor", default = 1, type = "float",
                      help="The sizing factor for junction simulation normalization to better calibrate the p-values (default: 1)")
                  
    parser.add_option("--not_paired_end", dest="not_paired_end", default = False,action="store_true",
                  help="Whether the reads are not paired end use to simulate junction count, only work when normalization_scale = junction\
                      (default: False)")
    
    
    
    
    
    (options, args) = parser.parse_args()



    if options.salmon_ref is None:
        sys.exit("Error: no salmon reference is provided...\n")
    
    if options.ref_dir is None:
        sys.exit("Error: no tealeaf reference directory provided (use --ref_dir).\n")
    
    
    if options.cell_mode not in {'pseudobulk', 'single_cell'}:
        sys.exit("Error: --cell_mode must be pseudobulk or single_cell.\n")

    if options.cell_mode == 'pseudobulk' and options.barcodes_clusters == None and options.pseudobulk_samples == None and options.preprocessed == False:
        sys.exit("Error: no barcodes to cluster/cell_type or pseudobulk file is provided...\n")
        
        
    
    sys.stderr.write(f"Processing alevin-fry results in: {options.alevin_dir}\n")
    sys.stderr.write(f"Salmon reference: {options.salmon_ref}\n")
    sys.stderr.write(f"tealeaf reference directory: {options.ref_dir}\n")
    sys.stderr.write(f"tealeaf reference prefix: {options.ref_prefix}\n")
    sys.stderr.write(f"Output prefix: {options.outprefix}\n")
    
    if options.preprocessed != False:
        sys.stderr.write(f'reading preprocessed pseudobulk to EC matrix using prefix {options.outprefix}\n')
    elif options.pseudobulk_samples != None:
        sys.stderr.write(f'reading barcodes to pseudobulk sample from {options.pseudobulk_samples}\n')
    else: 
        sys.stderr.write(f'reading barcodes to pseudobulk sample from {options.barcodes_clusters}\n')
        sys.stderr.write(f'merging barcodes using {options.pseudobulk_method}\n')
    
    
    if options.normalization_scale != 'junction' and options.normalization_scale != 'global' and options.normalization_scale != 'snRNA':
        sys.exit("Error: invalid normalization scale...\n")

    if options.quant_method not in {
        'em', 'nnls', 'nnls_nucnorm', 'admm', 'admm_factorized',
        'frank_wolfe', 'frank_wolfe_penalized', 'factorized',
    }:
        sys.exit("Error: invalid quantification method...\n")

    if options.regularization_target not in {'phi', 'theta'}:
        sys.exit("Error: --regularization_target must be phi or theta.\n")

    if options.ec_design not in {'legacy', 'binary', 'weighted', 'positional'}:
        sys.exit(
            "Error: --ec_design must be legacy, binary, weighted, or positional.\n"
        )
    if options.ec_design == 'positional' and not options.primer_pairs:
        sys.exit("Error: --ec_design positional requires --primer_pairs.\n")
    if options.primer_sampling_model not in {
        'effective_length', 'oligodt_tpm', 'all_tpm',
    }:
        sys.exit(
            "Error: --primer_sampling_model must be effective_length, "
            "oligodt_tpm, or all_tpm.\n"
        )

    if options.glm_data_backend not in {'auto', 'cuda', 'pinned'}:
        sys.exit("Error: --glm_data_backend must be auto, cuda, or pinned.\n")

        
    record = f'{options.outprefix}clustering_parameters.txt'
    sys.stderr.write(f'Detailed parameters will be saved to {record}\n')
    write_options_to_file(options, record)

    tealeaf_sc(options)
