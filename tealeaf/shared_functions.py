"""
Shared clustering functions used by both the bulk and single-cell pipelines.

Intron representation used throughout this module:
    [start, end, total_count, name, exon_set]
where *exon_set* is the set of neighbouring exon identifiers (strings).
"""

import pandas as pd
import warnings

from tealeaf.utils import timing_decorator, write_options_to_file

warnings.simplefilter(action='ignore', category=pd.errors.DtypeWarning)



@timing_decorator
def build_init_cluster(intron_count_file, connect_file):
    """
    Parameters
    ----------
    intron_count_file : str
        Path to a CSV file containing information about introns and their counts.
        Sample row format: Chr Start End Gene sample_1 sample_2........

    Adds a 'Cluster' column to the CSV with cluster assignments for introns.
    """
    
    
    # Read in the file
    df = pd.read_csv(intron_count_file, sep=' ')
    
    
    connect_df = pd.read_csv(connect_file, sep = ' ')
    intron_strand_dic = dict(zip(connect_df['intron'], connect_df['strand']))
    df["Strand"] = df.Name.map(intron_strand_dic)
                            
    
    
    # Sort so that overlapping-intron detection works correctly (separate + and - strand)
    df.sort_values(by=['Chr', 'Strand', 'Start', 'End'], inplace=True)

    # State variables for the sequential cluster-assignment pass
    cluster_num = 1
    cluster_end = df.iloc[0]['End']
    cluster_chr = df.iloc[0]['Chr']
    cluster_strand = df.iloc[0]['Strand']

    # Function to apply to each row to determine cluster
    def check_new_cluster(row):
        nonlocal cluster_num
        nonlocal cluster_chr
        nonlocal cluster_end
        nonlocal cluster_strand

        # Start a new cluster when chromosome, strand, or genomic position breaks the overlap chain.
        # No need to check Start explicitly because rows are sorted by Start.
        if row['Chr'] != cluster_chr or row['Strand'] != cluster_strand or row['Start'] > cluster_end:
            cluster_num += 1
            cluster_end = row['End']
            cluster_chr = row['Chr']
            cluster_strand = row['Strand']
            
        else:
            cluster_end = max(cluster_end, row['End'])
        return f'cluster_{cluster_num}'

    df.insert(1, 'Cluster', df.apply(check_new_cluster, axis=1))
    df.drop(columns = ['Strand'], inplace = True)
    # Save the modified DataFrame back to a file
    df.to_csv(intron_count_file, sep=' ', index=False)


@timing_decorator
def process_clusters(init_clus_file, exon_count_file, intron_to_exon_file, out_prefix = '', mode = 3,
                     cutoff = 0.1, percent_cutoff = 0.01, min_cluster_val = 1):
    """
    
    Parameters
    ----------
    init_clus_file : Str, the file contain the init cluster result
    exon_count_file : Str, the file contain exon counts
    intron_to_exon_file : Str, the file contain information about          

    out_prefix : prefix of output file
    mode: mode that use to refine clusters, 1 for overlap, 2 for overlap+share_intron_splice_site, \
        3 for overlap+share_intron_splice_site+shared_exon_splice_site

    Returns
    -------
    None.

    """
    # test run code: process_clusters('count_intron', 'count_exon', 'intron_to_exon.tsv')
    
    # format of a cluster: {start: ; end: ; splicesites: ; exon_splicesites: ; introns: [[start,end,count,name, [exon_slicesites]]}}

    
    exon_count = pd.read_csv(exon_count_file, sep = ' ', index_col = 0 )
    initial_clus = pd.read_csv(init_clus_file, sep = ' ', index_col = 0)
    intron_to_exon = pd.read_csv(intron_to_exon_file, sep = ' ', index_col = 0)
    
    num_cluters = 0
    if "Gene" in initial_clus.columns.tolist(): # deal with the result from different version
        samples = initial_clus.columns.tolist()[5:]
        contained_gene = True
    else:
        samples = initial_clus.columns.tolist()[4:]
        contained_gene = False
    
    output = open(f'{out_prefix}refined_cluster', 'w')
    
    out_str = ''
    for sample in samples:
        out_str += f'{sample} '
    out_str = out_str[:-1]
    print(out_str, file = output)

    # Sum the values for the samples. It will create a columns at -1 position 
    initial_clus['Sum'] = initial_clus[samples].sum(axis=1)
    clu_start = 0
    num_rows = len(initial_clus)
    
    
    
    for idx in range(num_rows + 1): # +1 to deal with -1 row have a diff clu than -2 row
    
        if idx != num_rows and initial_clus.iloc[clu_start]['Cluster'] == initial_clus.iloc[idx]['Cluster']:
            
            continue
        
        else: #
            if idx == num_rows:
                cluster_df = initial_clus.iloc[clu_start:]
            
            else: # there are more rows left
                cluster_df = initial_clus.iloc[clu_start:idx]
        
            clu_start = idx # update idx
            
            if len(cluster_df) < 2: # skip the clu that only contain one intron
                continue
            
            introns = []
            
            # build intron from initial cluster df
            for index, row in cluster_df.iterrows():
                intron = build_intron(index, row, samples, exon_count, intron_to_exon)
                # in format [start,end,count,name, {exon_slicesites}]
                introns.append(intron)
            
            
            cluster_dic = {}
            
            # perform filtering and reclu
            if mode == 1: #overlap only
                process_clu(introns, cluster_dic, cutoff, percent_cutoff, mode)
            else: 
                new_clus = []
                
                if mode == 2: #shared splice site
                    new_clus = refine_links(introns, exon_connection = False)
                else: 
                    new_clus = refine_links(introns, exon_connection = True)
            
                for clu in new_clus:
                    process_clu(clu, cluster_dic, cutoff, percent_cutoff, mode)
            
            

            # we should have a complete cluster_dic at this point
            for cluster_num in cluster_dic:
                clu_val = 0
                out_str = ''
                
                
                for clu_intron in cluster_dic[cluster_num]:
                    
                    infor = initial_clus.loc[clu_intron[3]]
                    strand = str(intron_to_exon.loc[clu_intron[3]]['strand'])
                    
                    out_str += f'{infor.Chr}:{infor.Start}:{infor.End}:clu_{num_cluters}_{strand}'
                    
                    if contained_gene == True:
                        for value in list(infor)[5:-1]: # start of sample values, and not include the sum columns
                            out_str += f' {value}'
                            clu_val += float(value)
                    else:
                        for value in list(infor)[4:-1]: # start of sample values, and not include the sum columns
                            out_str += f' {value}'
                            clu_val += float(value)
                            
                            
                    out_str += '\n'
            
                if clu_val <= min_cluster_val:
                    continue
                out_str = out_str[:-1] # get rid of the last \n
                print(out_str, file = output)
                    
          
                num_cluters += 1

    output.close()










def build_intron(index, row, samples, exon_count, intron_to_exon):
    """Build the intron list representation used by the cluster-processing functions.

    Parameters
    ----------
    index : str
        Row index of the intron (its name, e.g. ``chr1:1000-2000``).
    row : pandas.Series
        Row from the initial-cluster DataFrame (Cluster Chr Start End [Gene] samples… Sum).
    samples : list[str]
        Sample column names (unused here, kept for API symmetry).
    exon_count : pandas.DataFrame
        Exon-level count table, indexed by exon name.
    intron_to_exon : pandas.DataFrame
        Intron-to-exon connectivity table from ``intron_exon_connectivity.tsv``.

    Returns
    -------
    list : [start, end, total_count, name, exon_set]
    """
    intron = [row['Start'], row['End'], row['Sum'], index]

    # Virtual introns (intron_0 / intron_-1) have no connected exons.
    if index in intron_to_exon.index and pd.notna(intron_to_exon.at[index, 'near_exons']):
        exons = intron_to_exon.at[index, 'near_exons'].split(',')
    else:
        exons = []

    intron.append(set(exons))
    return intron
    




def process_clu(clu, output_dic, cutoff, percent_cutoff, mode):
    """ 
    this is the helper function of process_clusters 
    clu: a list that contain intron in format [[start,end,count,name, {exon_slicesites}]]
    output_dic: a dic that contain final cluster information
    cutoff: float or int, cutoff value for an intron
    percent_cutoff: float, percent cutoff value
    mode: str, the processing mode of the cluster
    
    """
    
    introns, reclu = filter_introns(clu, cutoff, percent_cutoff)
    
    if reclu == True and len(introns) > 1:
        
        clus = cluster_intervals(introns)
        
        
        for cluster in clus:
            
            if mode == 'overlap':
                process_clu(cluster,output_dic, cutoff, percent_cutoff, mode)
            
            elif mode == 'exon_connection' : 
                temp_clus = refine_links(cluster, True)
            
                for i in temp_clus:
                    process_clu(i,output_dic, cutoff, percent_cutoff, mode)
            else:
                temp_clus = refine_links(cluster, False)
            
                for i in temp_clus:
                    process_clu(i,output_dic, cutoff, percent_cutoff, mode)
            
    else: 
        if len(introns) > 1: # skip cluster that is less than 2 intron
            num_clus = len(output_dic)
            output_dic[num_clus] = introns
    
    
    
    
    


def filter_introns(introns, cutoff, percent_cutoff):
    """ 
    this is the helper function of process_clusters, will return a new introns list
    and True if reclu needed, flase if not needed
    introns:a list in format [[start,end,count,name, {exon_slicesites}]]
    cutoff: int or float
    percent_cutoff: float
    """
    
    total_TPM = 0
    for intron in introns: #loop over to get total TPM
        total_TPM += intron[2]
    
    new_introns = []
    reclu = False
    
    for intron in introns: 
        
        count = intron[2]
        if (count/total_TPM >= percent_cutoff and count >= cutoff):
            new_introns.append(intron)
        else:
            reclu = True
    if len(cluster_intervals(introns)) >= 2: 
        # this is necessary after refined link, as refined link doesn't ensure overlap between introns
        reclu = True

    return new_introns, reclu



def cluster_intervals(introns):
    """ 
    this is the helper function of process_clusters, it will generate a list of introns clu
    """
    introns.sort()
    
    clusters = []
    


    tmp_clu = [introns[0]] # init with the first intron
    current_intron = introns[0]    
    
    tmp_start = current_intron[0]
    tmp_end = current_intron[1]
    
    for i in range(1,len(introns)):
        
        current_intron = introns[i]
        
        if overlaps([tmp_start, tmp_end], [current_intron[0], current_intron[1]]):
            tmp_clu.append(current_intron)
        else:
            clusters.append(tmp_clu)
            tmp_clu = [current_intron]
        tmp_start = current_intron[0]
        tmp_end = max(tmp_end, current_intron[1])
            

    if len(tmp_clu) > 0:
        clusters.append(tmp_clu)



    return clusters




def overlaps(A,B):
    '''
    Checks if range A overlap with range B
    range in format [start,end]
    '''

    if A[1] < B[0] or B[1] < A[0]:
        return False
    else: return True




def refine_links(introns, exon_connection = True):

    """ 
    this is the helper function of process_clusters, it will generate a list of Intron_cluster object
    input should be a list of introns in format [[start,end,count,name, {exon_slicesites}]]
    """
    
    unassigned = introns[1:]
    current = [introns[0]]
    

    splicesites = set([current[0][0],current[0][1]])
    exon_splicesites = set(current[0][4])
    newClusters = []
    
    
    while len(unassigned) > 0:
        finished = False

        while not finished:
            finished = True
            torm = []
            for intron in unassigned:
                count = intron[2]
                start = intron[0]
                end = intron[1]
                
                
                
                if start in splicesites or end in splicesites:
                    current.append(intron)
                    splicesites.add(start)
                    splicesites.add(end)
                    finished = False
                    torm.append(intron)
                    exon_splicesites = exon_splicesites | intron[4] # union of exon_sites
                    
                    
                    
                elif exon_connection == True: 
                    if len(intron[4] & exon_splicesites) > 0:
                        current.append(intron)
                        splicesites.add(start)
                        splicesites.add(end)
                        finished = False
                        torm.append(intron)
                        exon_splicesites = exon_splicesites | intron[4] # union of exon_sites
               
                    
               
            for intron in torm:
                unassigned.remove(intron)

                
                
        newClusters.append(current)
        current = []
        if len(unassigned) > 0:
            current = [unassigned[0]]
            splicesites = set([current[0][0],current[0][1]])
            exon_splicesites = set(current[0][4])
            unassigned = unassigned[1:]
    
    if len(current) > 0:
        newClusters.append(current)

    return newClusters


@timing_decorator
def compute_ratio(cluster_file, output_prefix=''):
    """Compute per-intron PSI (ratio) values from a refined-cluster count file.

    Divides each intron's count by the total count of its cluster for every
    sample, producing a ratio in [0, 1].  Results are written to
    ``{output_prefix}ratio_count`` in a format compatible with
    ``leafcutter_ds``.

    Parameters
    ----------
    cluster_file : str
        Path to the ``refined_cluster`` file produced by :func:`process_clusters`.
    output_prefix : str, optional
        Prefix prepended to the output file name.
    """
    df = pd.read_csv(cluster_file, sep=' ', index_col=0)
    samples = df.columns
    
    # Split the index to get the cluster information
    df['Cluster'] = df.index.to_series().apply(lambda x: x.split(':')[3])
    
    # Calculate the sum of each sample within each cluster
    cluster_sums = df.groupby('Cluster').sum()
    
    # Join the cluster sums back to the original dataframe
    df = df.join(cluster_sums, on='Cluster', rsuffix='_sum')
    
    
    # Calculate ratios in a vectorized way
    for sample in samples:  
        df[sample] = df[sample] / df[f'{sample}_sum']
    
    # Prepare the output dataframe, dropping sum columns and the Cluster column
    df.fillna(0, inplace=True)
    df = df.round(6)
    output_df = df[[col for col in df.columns if '_sum' not in col and col != 'Cluster']]
    # Exclude the last two columns which are 'Cluster' and sums
    
    # Write the result to a file
    output_df.to_csv(f'{output_prefix}ratio_count', sep=' ')

    # Strip the leading space on the header line so the file is compatible with leafcutter_ds.
    with open(f'{output_prefix}ratio_count', 'r') as file:
        lines = file.readlines()

    lines[0] = lines[0].lstrip()

    with open(f'{output_prefix}ratio_count', 'w') as file:
        file.writelines(lines)































