"""tealeaf bulk / pseudobulk clustering pipeline.

Reads Salmon quant.sf files, maps isoform abundances to introns, builds
initial overlap-based clusters and refines them.  The final outputs are
``{outprefix}refined_cluster`` and ``{outprefix}ratio_count``, which are
compatible with downstream LeafCutter tools (leafcutter_ds, leafviz).
"""

import numpy as np
import pandas as pd
import sys
import warnings

from tealeaf.utils import timing_decorator, write_options_to_file
from tealeaf.shared_functions import (
    build_init_cluster,
    process_clusters,
    compute_ratio,
)

warnings.simplefilter(action='ignore', category=pd.errors.DtypeWarning)


chr_dic = {'chr1': 0,'chr2': 1,'chr3': 2,'chr4': 3,'chr5': 4,'chr6': 5,'chr7': 6,'chr8': 7,'chr9': 8,'chr10': 9,
               'chr11': 10,'chr12': 11,'chr13': 12, 'chr14': 13,'chr15': 14,'chr16': 15,'chr17': 16,'chr18': 17,
               'chr19': 18,'chr20': 19,'chr21': 20, 'chr22': 21, 'chrX':22, 'chrY':23}




def transcript_to_intron_counts(row, trans_int_map, sample_name, introns_dic, exons_dic,TPM_count, unmapped_transcript):
    """
    
    Parameters
    ----------

    row: a df row from quant.sf
    trans_int_map : a dataframe that contain map for transcript to introns, name of transcript 
                    will be index
                    row in format: transcript_id chr gene_id support_introns support_exons
                    
    
    sample: str, sample name

    introns_dic : a dictionary in format {intron:[gene,{sample_name:value}]}
    
    exons_dic : a dictionary in format {exon:[gene,{sample_name:value}]}
    normalize_count: a boolean value, whether use normalize count or TPM
    unmapped_transcript: a dic use to record unmapped_transcript that appear 
    
    
    Returns
    -------
    None.
    The dic will be updated, no return is needed
    """

  
    transcript_name = row['Name']
    if TPM_count == False:
        value = row['normalized_count']  
    else: 
        value = row['TPM']

  
    if (transcript_name in trans_int_map.index) == False:
        unmapped_transcript.add(transcript_name)
        return None
    
    transcript = trans_int_map.loc[transcript_name]
    chromosome = transcript['Chr']
    gene = transcript['Gene']

    
  
    if transcript['support_introns'] != None:    # if there have no support intron, skip to exon section
        
        introns = transcript['support_introns'].split(',') # list of introns
        for intron in introns:
            if intron in introns_dic:
                if sample_name in introns_dic[intron][1]: #{sample_name:TPM}
                    introns_dic[intron][1][sample_name] += value
                else: 
                    introns_dic[intron][1][sample_name] = value
            else:
                introns_dic[intron] = [gene, {sample_name:value}]

        
    if transcript['support_exons'] != None:   # if there have no support intron, skip to next section
        exons = transcript['support_exons'].split(',') # list of exons

        for exon in exons:
            if exon in exons_dic:
                if sample_name in exons_dic[exon][1]: #{sample_name:TPM}
                    exons_dic[exon][1][sample_name] += value
                else: 
                    exons_dic[exon][1][sample_name] = value
            else:
                exons_dic[exon] = [gene, {sample_name:value}]
            
    return None
    


def dic_to_csv(dic, output_name, samples):
    """

    Parameters
    ----------
    dic : an intron or exon dic in format {intron:[gene,{sample_name:TPM}]}
    
    output_name : name of output 
    
    samples: name of samples 

    Returns
    -------
    None
    wirte to disk

    """

    output = open(output_name, 'w')
    
    out_str = 'Name Chr Start End Gene'
    for i in samples:
        out_str += f' {i}'
    
    
    
    print(out_str, file = output)
    
    for key in dic: # key = intron/exon e.g. chr1:10000-12000
    
        tmp = key.split(':')
        chromosome = tmp[0] 
        start, end = tmp[1].split('-')
        
        gene_name = dic[key][0]
        val_dic = dic[key][1]  # {sample:value} 

        tmp_str = f'{key} {chromosome} {start} {end} {gene_name}'
        
        for sample in samples:
            
            
            if sample in val_dic:
                tmp_str += f' {val_dic[sample]}'
            else:
                tmp_str += ' 0'
        
        print(tmp_str, file = output)
                

    output.close()






@timing_decorator
def count_introns(samples, trans_int_map_file, out_prefix = '', threshold = 0, TPM_count = False):
    """
    Parameters
    ----------
    samples : list, names of sample files
          
    trans_int_map_file : str, location of transcript_to_intron map file
        
    out_prefix : str, optional, prefix of output file
    
    threshold: the threshold use for filtering out introns, should decide based on TPM or normalize_count
        
    example: count_introns（['liver_quant.sf', 'heart_quant.sf'],'hg38_transcript_intron_map.tsv' ）
    Returns
    -------
    df of introns count infor

    """
    
    
    try:
        trans_map = pd.read_csv(trans_int_map_file, sep=' ')
    except Exception:
        sys.exit("%s does not exist... check your map file.\n" % trans_int_map_file)

    trans_map = trans_map.replace({np.nan: None})
    trans_map.index = trans_map.Transcript    # use transcript id as index to facilitate search
    introns_dic = {}
    exons_dic = {}                            
    unmapped_transcript = set()
    
    
    base_names = []
    for name in samples:
        try:
            df = pd.read_csv(name, sep = '\t')
            if TPM_count == False:
                df = df[df.normalized_count > threshold].reset_index(drop = True)   # drop trancript expressed too low
            else:
                df = df[df.TPM > threshold].reset_index(drop = True)
        except:
            sys.exit("%s does not exist... check your count files.\n"%name) 
          
        base_name = name.split("/")[-1]
        base_names.append(base_name)

        _ = df.apply(transcript_to_intron_counts, axis = 1, args=(trans_map,base_name, introns_dic, exons_dic, TPM_count, unmapped_transcript))

    dic_to_csv(introns_dic, f'{out_prefix}count_intron', base_names)  # save the dic in csv format
    dic_to_csv(exons_dic, f'{out_prefix}count_exon', base_names)
    not_mapped_transcript = open(f'{out_prefix}dropped_trancripts', 'w')
    for droped in unmapped_transcript:
        print(droped, file = not_mapped_transcript)
    not_mapped_transcript.close()
    
    
    
    if len(unmapped_transcript) == 0:
        print('All transcript mapped to intron')
    
    elif len(unmapped_transcript) <= 500:
        print(f'There are {len(unmapped_transcript)} transcript not found in Transcript to intron map, like due to transcript id difference or filtering step')
    else:
        print(f'There are {len(unmapped_transcript)} transcript not found in Transcript to intron map, please check the version of annotation or may due to filtering step')
    
    return pd.read_csv(f'{out_prefix}count_intron', sep=' ')

def input_file_processing(input_file):
    """Read a plain-text file of sample paths (one per line) and return them as a list."""
    result = []
    with open(input_file, 'r') as file:
        for line in file:
            result.append(line.rstrip('\n'))
    return result




def count_TPM_normalization_local(quant_file, gft_file, out_prefix= None):
    # Your function implementation
    df = pd.read_csv(quant_file, sep='\t')
    df = df[df['TPM'] > 0]
    
    transcript_to_gene_dic = extract_transcript_to_gene_map(gft_file)
    
    gene_dic = {}
    df['gene'] = df['Name'].map(transcript_to_gene_dic)
    
    for index, row in df.iterrows():
        TPM = row['TPM']
        count = row['NumReads']
        gene = row['gene']
        
        gene_dic.setdefault(gene, {'TPM': 0, 'count': 0})
        gene_dic[gene]['TPM'] += TPM
        gene_dic[gene]['count'] += count

    df['normalized_count'] = df.apply(
        lambda x: (x['TPM'] / gene_dic[x['gene']]['TPM']) * gene_dic[x['gene']]['count'], axis=1)
    
    # Use the out_prefix to determine the output file name
    output_file = f'{out_prefix}_normalized.sf' if out_prefix else 'normalized.sf'
    df.to_csv(output_file, sep='\t', index=False)

def count_TPM_normalization_global(quant_file, gft_file, out_prefix= None):
    # Your function implementation
    df = pd.read_csv(quant_file, sep='\t')
    df = df[df['TPM'] > 0]


    scale = df['NumReads'].sum() / df['TPM'].sum()
  

    df['normalized_count'] = scale * df['TPM'] 

    # Use the out_prefix to determine the output file name
    output_file = f'{out_prefix}_normalized.sf' if out_prefix else 'normalized.sf'
    df.to_csv(output_file, sep='\t', index=False)
    
    
    

def count_TPM_normalization_junction_simulation(quant_file, read_length, paired_end = True, overhang = 2, sizing_factor = 1, out_prefix = None):
    """
    This function take a quant.sf file from salmon output and normalize it based on
    expected_count * effective_lenth/read_length to simulate a count distribution for a junction 
    The read_length double when the read is pair_end
    
    """
    
    df = pd.read_csv(quant_file, sep='\t')
    df = df[df['TPM'] > 0]
    if paired_end == True: 
        eff_read_length = (read_length - overhang) * 2 * sizing_factor
    else:
        eff_read_length = (read_length - overhang) * sizing_factor
        
    df['normalized_count'] = df['NumReads'] * (eff_read_length/df['EffectiveLength'])
    
    output_file = f'{out_prefix}_normalized.sf' if out_prefix else 'normalized.sf'
    df.to_csv(output_file, sep='\t', index=False)
    



def extract_transcript_to_gene_map(gtf_file):
    """
    generate a map that map transcript_id to gene_id

    Parameters
    ----------
    gtf_file : a string that contain the name for the gtf file

    Returns
    -------
    map_transcript_to_gene : a dic that map transcript_id to gene_id

    """
    
    map_transcript_to_gene = {}

    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip header lines

            fields = line.strip().split('\t')
            
            if fields[2] == 'transcript':
                attributes = fields[8]
                attribute_parts = attributes.split(';')

                transcript_id = ''
                gene_id = ''

                for attribute in attribute_parts:
                    if 'transcript_id' in attribute:
                        transcript_id = attribute.split('"')[1]
                    elif 'gene_id' in attribute:
                        gene_id = attribute.split('"')[1]

                if transcript_id and gene_id:
                    map_transcript_to_gene[transcript_id] = gene_id

    return map_transcript_to_gene



@timing_decorator
def tealeaf_clustering(options):
    """
    This is the main function for tealeaf clustering
    """

    samples = input_file_processing(options.count_files)
    
    # perform normalization
    if options.use_TPM == False and options.preprocessed == False:
        new_samples = []
        for sample in samples:
            prefix = sample.split(".")[0] # get rid of the .sf

            if options.normalization_scale == 'junction':
                
                paired_end = (not options.not_paired_end)
                
                count_TPM_normalization_junction_simulation(sample, options.read_length, paired_end, \
                                                            options.overhang, options.sizing_factor, prefix)
            elif options.normalization_scale == 'global':
                count_TPM_normalization_global(sample, options.annot, prefix)
            elif options.normalization_scale == 'local':
                count_TPM_normalization_local(sample, options.annot, prefix)
            else: 
                sys.exit("Error: invalid normalization scale...\n")
              
            new_samples += [f'{prefix}_normalized.sf']
               
        samples = new_samples
        # saved the normalized sample for furture usage
        with open(f'{options.count_files.split(".")[0]}_normalized.txt', 'w') as file:
            for sample in samples:
                file.write(f"{sample}\n")
            
        sys.stderr.write("Finished Normalization\n")
        

    out_prefix = options.outprefix
    
    # count_introns
    count_introns(samples, \
                  options.map, \
                  out_prefix= out_prefix,\
                  threshold = options.samplecutoff,\
                  TPM_count = options.use_TPM)
    sys.stderr.write("Finished Introns Counting\n")
    
    
    # build initial cluster
    build_init_cluster(f'{out_prefix}count_intron', options.connect_file)
    sys.stderr.write("Finished Initial Clustering\n")


    process_clusters(f'{out_prefix}count_intron', f'{out_prefix}count_exon', \
                     options.connect_file, \
                     out_prefix = options.outprefix,\
                     cutoff = options.introncutoff, \
                     percent_cutoff = options.mincluratio,\
                     min_cluster_val = options.minclucounts)

    sys.stderr.write("Finished Cluster refinement\n")
    
    compute_ratio(f'{out_prefix}refined_cluster', out_prefix)

    sys.stderr.write("Finished PSI calculation\n")

    record = f'{options.outprefix}clustering_parameters.txt'
    sys.stderr.write(f'Saving parameters to {record}\n')
    write_options_to_file(options, record)

    sys.stderr.write('Clustering finished\n')


if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("--map", dest="map", default=None,
                      help="isoform-to-intron map file generated by tealeaf-map (default: None)")

    parser.add_option("--count_files", dest="count_files", default=None,
                      help="text file listing one Salmon quant.sf path per line (default: None)")

    parser.add_option("--connect_file", dest="connect_file", default=None,
                      help="intron-exon connectivity file from tealeaf-map (default: None)")

    parser.add_option("-a", "--annot", dest="annot", default=None,
                      help="transcriptome annotation GTF, required when not using --use_TPM (default: None)")

    parser.add_option("--cluster_def", dest="cluster_def", default=3, type="int",
                      help="cluster refinement mode: 1=overlap, 2=overlap+shared_intron_splice_site, "
                           "3=overlap+shared_intron_splice_site+shared_exon_splice_site (default: 3)")

    parser.add_option("-o", "--outprefix", dest="outprefix", default='tealeaf_',
                      help="output file prefix; include directory path if not the current directory "
                           "(default: tealeaf_)")

    parser.add_option("--normalization_scale", dest="normalization_scale", default='junction',
                      help="normalisation mode: 'junction' (simulate junction counts), "
                           "'local' (gene-level), or 'global' (default: junction)")

    parser.add_option("--use_TPM", dest="use_TPM", action="store_true", default=False,
                      help="use Salmon TPM values directly instead of normalised counts (default: False)")

    parser.add_option("--preprocessed", dest="preprocessed", default=False, action="store_true",
                      help="skip normalisation if count files are already normalised (default: False)")

    parser.add_option("--samplecutoff", dest="samplecutoff", default=0, type="float",
                      help="minimum count for an intron in a sample to be considered expressed (default: 0)")

    parser.add_option("--introncutoff", dest="introncutoff", default=5, type="float",
                      help="minimum summed count for an intron across samples (default: 5)")

    parser.add_option("-m", "--minclucounts", dest="minclucounts", default=30, type="float",
                      help="minimum total count to retain a cluster (default: 30)")

    parser.add_option("-r", "--mincluratio", dest="mincluratio", default=0.01, type="float",
                      help="minimum fraction of cluster reads supporting an intron (default: 0.01)")

    parser.add_option("--read_len", dest="read_length", default=100, type="int",
                      help="read length for junction-count simulation; used when "
                           "normalization_scale=junction (default: 100)")

    parser.add_option("--overhang", dest="overhang", default=2, type="int",
                      help="overhang bases for junction-count simulation (default: 2)")

    parser.add_option("--sizing_factor", dest="sizing_factor", default=1, type="float",
                      help="scaling factor for junction-count simulation (default: 1)")

    parser.add_option("--not_paired_end", dest="not_paired_end", default=False, action="store_true",
                      help="treat reads as single-end for junction-count simulation (default: False)")

    (options, args) = parser.parse_args()

    if options.count_files is None:
        sys.exit("Error: no --count_files provided.\n")

    if options.map is None:
        sys.exit("Error: no --map (isoform-to-intron map) provided.\n")

    if options.connect_file is None:
        sys.exit("Error: no --connect_file provided.\n")

    if not options.use_TPM and options.annot is None:
        sys.exit("Error: --annot is required when not using --use_TPM.\n")

    tealeaf_clustering(options)


