"""tealeaf map-generation pipeline.

Parses a GTF transcriptome annotation and produces the isoform→intron and
isoform→exon mapping files consumed by ``tealeaf-cluster`` and ``tealeaf-sc``.
Only needs to be run once per annotation file.
"""

import pyranges as pr
import numpy as np
import pandas as pd
import sys
from pathlib import Path
import warnings

from tealeaf.utils import timing_decorator, write_options_to_file
import scipy
import scipy.sparse
from scipy.sparse import csr_matrix, save_npz, load_npz
from optparse import OptionParser

warnings.simplefilter(action='ignore', category=pd.errors.DtypeWarning)


####### isoform to intron and exon map generation
def compute_transcript_intron_map(annot, out_prefix = '', annot_type = 'gencode', min_length = 50, max_length = 5000000, no_quality_control = False):
    """

    Parameters
    ----------
    annot : a gtf format file that contain the annotation for the genome
    out_prefix : the output prefix, example: "mm10_"
    annot_type : support gencode, Stringtie (in progress, will add support to ref_seq, ensembl)
    min_length: int, the minimun intron length for intron to be consider 
    max_length: int, the maximun intron length to be consider
    quality_control: whether to get rid of pseudogene, decay transcript

    Returns
    -------
    None.
    Result will be saved to disk

    """
    num_gene = 0
    num_transcript = 0
    
    if not Path(annot).is_file():
        raise FileNotFoundError(f"{annot} does not exist... check your annotation files.")
    
    gr = pr.read_gtf(annot)
    df = gr.df
    df = df.replace({np.nan: None})
    df = df[df['Chromosome'].str.contains('chr')].reset_index(drop=True) #eliminate gene and transcript not mapped into a chr
    df = df[~ df['gene_type'].str.contains('artifact', na = False)]
    df = df[df['Chromosome'] != ('chrM')].reset_index(drop=True) # eliminate mito genes as they do not undergo alternative splicing

    # deal with output from stringtie
    if annot_type == 'Stringtie':
        gene_dic = {}
        name_dic = {}

        filtered_df = df[df['ref_gene_id'].notna()]

        # Use `set` to get unique gene_ids and update dictionaries
        for gene_id, gene_name, ref_gene_id in zip(filtered_df['gene_id'], filtered_df['gene_name'], filtered_df['ref_gene_id']):
            if gene_id not in gene_dic or gene_id not in name_dic:
                name_dic[gene_id] = gene_name
                gene_dic[gene_id] = ref_gene_id

        df['ref_gene_id'] = df['gene_id'].map(gene_dic).fillna(df['ref_gene_id'])
        df['gene_name'] = df['gene_id'].map(name_dic).fillna(df['gene_name'])   
    
        with open(f'{out_prefix}stringtie_transcriptome_reference.tsv', 'w') as reference_list:
            print('Stringtie_id reference_id reference_name', file=reference_list)
            for key in gene_dic:
                print(f'{key} {gene_dic[key]} {name_dic[key]}', file=reference_list)

    elif no_quality_control == False:
        df = df[~ df['gene_type'].str.contains('pseudogene', na = False)]
        df = df[~ df['transcript_type'].str.contains('decay', na = False)]
        # df = df[~ df['transcript_type'].str.contains('retained_intron', na = False)]
        # maybe useful to account for retain_intron isoforms

    gene_start_index = 0
    num_rows = len(df) - 1

    output = open(f'{out_prefix}isoform_intron_map.tsv', 'w')

    print('Chr Gene Transcript support_introns support_exons Transcript_type', file = output)

    near_exon_dic = {}  # the dic that keep infor about the nearby exon of a intron
    # should in format {intron: {exon: counts}} counts not use 
    intron_str_dic = {}

    for idx in range(len(df)): # don't need to deal with the case that -1 row have a different gene than -2 row, since a gene must have at lease one transcript
        if df.iloc[gene_start_index]['gene_id'] == df.iloc[idx]['gene_id'] and idx != num_rows:
            continue
        else: #
            if idx == num_rows:
                gene_df = df.iloc[gene_start_index:]
            
            else: # there are more rows left
                gene_df = df.iloc[gene_start_index:idx]

            num_gene += 1
            gene = gene_df['gene_id'].unique()[0]
            chromosome = gene_df['Chromosome'].unique()[0]
            transcripts = list(gene_df['transcript_id'].unique()) #get all gene_id, will include None 

            for trans in transcripts:
                if trans is None: # skip the None value 
                    continue 
                num_transcript += 1
                exons = []
                trans_df = gene_df[(gene_df['transcript_id'] == trans) & (gene_df['Feature'] == 'exon')] # only care about exon for the transcript
                trans_type = trans_df['transcript_type'].unique()[0] if annot_type != 'Stringtie' else 'unknown'
        
                for index,row in trans_df.iterrows():
                    exon = [int(row['exon_number']), row['Chromosome'], row['Start'], row['End']] # exon number is the order of exon for this transcript
                    exons += [exon]
        
                # computer introns:
                exons = sorted(exons)  # guranttee that exon will follow the order
                introns = []
                for i in range(len(exons) - 1):  # -1 as only len(exons) - 1 introns exist 
                    current_exon = exons[i] 
                    next_exon = exons[i + 1]
                    current_exon_name = f'{current_exon[1]}:{current_exon[2]}-{current_exon[3]}'
                    next_exon_name = f'{next_exon[1]}:{next_exon[2]}-{next_exon[3]}'
                    if current_exon[3] < next_exon[2]: # current_exon end < next_exon start, foward direction
                        introns += [[i+1,current_exon[1], (current_exon[3]), next_exon[2] + 1]] # +1 to align with leafcutter intron
                        # the precise location could be problematic, but it is defined th
                        temp_name = f'{current_exon[1]}:{current_exon[3]}-{next_exon[2] + 1}'
                        add_near_exon_dic(near_exon_dic, temp_name, current_exon_name, next_exon_name)
                        intron_str_dic[temp_name] = '+'

                    else: # current_exon start > next_exon end, backward direction
                        introns += [[i+1,current_exon[1], (next_exon[3]), current_exon[2] + 1]] # +1 to align with leafcutter intron
                        temp_name = f'{current_exon[1]}:{next_exon[3]}-{current_exon[2] + 1}'
                        add_near_exon_dic(near_exon_dic, temp_name, current_exon_name, next_exon_name)
                        intron_str_dic[temp_name] = '-'

                # the order of intron doesn't really matter, but will keep in case have usage in future
                out_introns = ','.join(
                    f'{chromosome}:{i[2]}-{i[3]}' for i in introns
                    if min_length < abs(i[2] - i[3]) < max_length
                )
                out_exons = ','.join(f'{chromosome}:{i[2]}-{i[3]}' for i in exons)
                print(f'{chromosome} {gene} {trans} {out_introns} {out_exons} {trans_type}', file = output)
            gene_start_index = idx
        
    output.close()
    print_near_exon_dic(near_exon_dic, intron_str_dic, out_prefix)
    
    sys.stderr.write("\nFinished process %d genes and %d isofroms\n"% (num_gene, num_transcript))



def add_near_exon_dic(near_exon_dic, intron_name, current_exon_name, next_exon_name):
    exon_counts = near_exon_dic.setdefault(intron_name, {})
    exon_counts[current_exon_name] = exon_counts.get(current_exon_name, 0) + 1
    exon_counts[next_exon_name] = exon_counts.get(next_exon_name, 0) + 1
                        

def print_near_exon_dic(near_exon_dic, intron_str_dic, out_prefix):
    with open(f'{out_prefix}intron_exon_connectivity.tsv', 'w') as output:
        print('intron near_exons strand', file=output)
        for intron, exon_counts in near_exon_dic.items():
            exons = ','.join(exon_counts.keys())
            print(f'{intron} {exons} {intron_str_dic[intron]}', file=output)
    

def isoform_intron_exon_sparse_generation(isoform_intron_map_file, out_prefix = ''):
    """
    This function will generate two sparse matrices that map isoform to intron and exon     

    Returns
    -------
    None.

    """
    df = pd.read_csv(isoform_intron_map_file, sep = ' ')

    isoform_to_introns = {
        iso: [] if pd.isna(v) else v.split(',') for iso, v in zip(df['Transcript'], df['support_introns'])
    }
    isoform_to_exons = {
        iso: [] if pd.isna(v) else v.split(',') for iso, v in zip(df['Transcript'], df['support_exons'])
    }

    all_isoforms = list(isoform_to_introns.keys())
    all_introns = list({i for introns in isoform_to_introns.values() for i in introns})
    all_exons = list({e for exons in isoform_to_exons.values() for e in exons})

    # Step 2: Create index mappings
    isoform_to_index = {iso: i for i, iso in enumerate(all_isoforms)}
    intron_to_index = {intron: i for i, intron in enumerate(all_introns)}
    exon_to_index = {exon: i for i, exon in enumerate(all_exons)}

    # Step 3: Build sparse matrices
    def _build_matrix(mapping, item_to_index, num_cols):
        rows, cols = zip(*[
            (isoform_to_index[iso], item_to_index[item])
            for iso, items in mapping.items() for item in items
        ]) if any(mapping.values()) else ([], [])
        return scipy.sparse.coo_matrix(
            ([1] * len(rows), (list(rows), list(cols))),
            shape=(len(all_isoforms), num_cols)
        )

    isoform_intron_matrix = _build_matrix(isoform_to_introns, intron_to_index, len(all_introns))
    isoform_exon_matrix = _build_matrix(isoform_to_exons, exon_to_index, len(all_exons))

    save_npz(f'{out_prefix}isoform_intron_matrix.npz', isoform_intron_matrix)
    save_npz(f'{out_prefix}isoform_exon_matrix.npz', isoform_exon_matrix)

    for path, items in [
        (f'{out_prefix}isoform_rows.txt', all_isoforms),
        (f'{out_prefix}intron_cols.txt', all_introns),
        (f'{out_prefix}exon_cols.txt', all_exons),
    ]:
        with open(path, 'w') as f:
            f.write('\n'.join(items) + '\n')


##############################################################################################
# these function is some additional features, not fully tested yet
def add_virtual_first_last_introns(trancript_to_intron_map, intron_to_exon_map, out_prefix = '', include_exon = False):
    """
    This is the function that use to add virtual intron 0 and -1 to capture 
    AFE and ALE usage. Will use a location before or after all transcript to generate 
    a virtual exon and compute virtual intron between it and the first exon of the 
    transcript

    """
    df = pd.read_csv(trancript_to_intron_map, sep = ' ')
    df = df.where(pd.notna(df), None)
    gene_dic = {}
    intron_exon_df = pd.read_csv(intron_to_exon_map, sep = ' ')
    
    # this round to get the start and end for a gene
    for index, row in df.iterrows():
        current_gene = row['Gene']
        current_exons= row['support_exons'].split(',')
        
        # exon in format "chr1:3143475-3144545"
        first_exon = current_exons[0].split(':')[1].split('-')
        last_exon = current_exons[-1].split(':')[1].split('-')

        # tmp_start always smaller than tmp_end
        if int(first_exon[0]) <= int(last_exon[0]): # forward direction
            tmp_start = int(first_exon[0])
            tmp_end = int(last_exon[1]) 
            tmp_direction = '+'
        else: # reverse direction
            tmp_start = int(last_exon[0]) 
            tmp_end = int(first_exon[1])
            tmp_direction = '-'
        
        if current_gene not in gene_dic:
            gene_dic[current_gene] = {'start':tmp_start, 'end': tmp_end, 'direction': tmp_direction}
            
        else:
            if tmp_start < gene_dic[current_gene]['start']:
                gene_dic[current_gene]['start'] = tmp_start
            if tmp_end > gene_dic[current_gene]['end']:
                gene_dic[current_gene]['end'] = tmp_end
            
            if tmp_direction != gene_dic[current_gene]['direction'] and len(current_exons) >1:
                #correct the possible direction mistake when single exon appear
                gene_dic[current_gene]['direction'] = tmp_direction

    intron_exon_dic = {}
    intron_dic = {}

    for index, row in df.iterrows():
        current_chr = row['Chr']
        current_gene = row['Gene']
        current_exons= row['support_exons'].split(',')
        # exon in format "chr1:3143475-3144545"
        first_exon = current_exons[0].split(':')[1].split('-')
        last_exon = current_exons[-1].split(':')[1].split('-')
        
        direction = gene_dic[current_gene]['direction']
        
        if direction == '+':
            intron_0 = f'{current_chr}:{gene_dic[current_gene]["start"] - 500}-{int(first_exon[0])}'
            intron_minus1 = f'{current_chr}:{int(last_exon[1])}-{gene_dic[current_gene]["end"] + 500}'

        else:
            intron_0 = f'{current_chr}:{int(first_exon[1])}-{gene_dic[current_gene]["end"] + 500}'
            intron_minus1 = f'{current_chr}:{gene_dic[current_gene]["start"] - 500}-{int(last_exon[0])}'

        if include_exon:
            if intron_0 not in intron_exon_dic:
                intron_exon_dic[intron_0] = {'exons': {current_exons[0]}, 'strand': direction}
            else:
                intron_exon_dic[intron_0]['exons'].add(current_exons[0])
            if intron_minus1 not in intron_exon_dic:
                intron_exon_dic[intron_minus1] = {'exons': {current_exons[-1]}, 'strand': direction}
            else:
                intron_exon_dic[intron_minus1]['exons'].add(current_exons[-1])
        else:
            intron_exon_dic.setdefault(intron_0, {'exons': None, 'strand': direction})
            intron_exon_dic.setdefault(intron_minus1, {'exons': None, 'strand': direction})

        cur_introns = df.loc[index, 'support_introns']
        if cur_introns is None:
            df.loc[index, 'support_introns'] = f'{intron_0},{intron_minus1}'
        else:
            df.loc[index, 'support_introns'] = f'{intron_0},{cur_introns},{intron_minus1}'

        intron_dic.setdefault(intron_0, 'intron_0')
        intron_dic.setdefault(intron_minus1, 'intron_-1')

    intron_exon_list = [
        [intron,
         ','.join(d['exons']) if include_exon else d['exons'],
         d['strand']]
        for intron, d in intron_exon_dic.items()
    ]

    with open(f'{out_prefix}virtual.tsv', 'w') as list_virtual:
        print('intron type', file=list_virtual)
        for intron, itype in intron_dic.items():
            print(f'{intron} {itype}', file=list_virtual)
    df.to_csv(f'{out_prefix}isoform_intron_map_with_virtual.tsv', sep=' ', index=False)
    intron_exon_df = pd.concat(
        [intron_exon_df, pd.DataFrame(intron_exon_list, columns=['intron', 'near_exons', 'strand'])],
        ignore_index=True
    )
    intron_exon_df.to_csv(f'{out_prefix}intron_exon_connectivity_with_virtual.tsv', sep=' ', index=False)


def intron_source_generation(transcript_intron_map, out_prefix = ''):
    df = pd.read_csv(transcript_intron_map, sep = ' ')
    df = df[df['support_introns'].notna()]
    df['support_introns'] = df['support_introns'].str.split(',')
    df = df.explode('support_introns')
    df = df.rename(columns={'support_introns': 'intron', 'Transcript_type': 'source_type'})
    df = df[['intron', 'source_type']]
    # we now get a intron source type map, but there are duplucation
    df = df.groupby('intron')['source_type'].agg(lambda x: ','.join(set(x))).reset_index()
    # get rid of duplication of intron and merge all source_type that can mapped to an intron
    
    df.to_csv(f'{out_prefix}intron_source_map.tsv')



#######################################################################################################################################
@timing_decorator
def tealeaf_map_generation(options):
    out_prefix = options.outprefix
    compute_transcript_intron_map(
        options.annot,
        out_prefix=out_prefix,
        annot_type=options.annot_source,
        min_length=options.minintronlen,
        max_length=options.maxintronlen,
        no_quality_control=options.no_quality_control,
    )
    
    if options.single_cell:
        isoform_intron_exon_sparse_generation(f'{out_prefix}isoform_intron_map.tsv', out_prefix)
        
    if options.annot_source == 'gencode':
        intron_source_generation(f'{out_prefix}isoform_intron_map.tsv',out_prefix)
        
    if options.virtual_intron:
        add_virtual_first_last_introns(
            f'{out_prefix}isoform_intron_map.tsv',
            f'{out_prefix}intron_exon_connectivity.tsv',
            out_prefix,
        )


if __name__ == "__main__":

    parser = OptionParser()

    parser.add_option('-a', "--annot", dest="annot", default=None,
                      help="transcriptome annotation GTF file (required)")

    parser.add_option("--annot_source", dest="annot_source", default='gencode',
                      help="annotation source: 'gencode' or 'Stringtie' (default: gencode)")

    parser.add_option("-o", "--outprefix", dest="outprefix", default='tealeaf_',
                      help="output file prefix; include directory path if not the current directory "
                           "(default: tealeaf_)")

    parser.add_option("--maxintronlen", dest="maxintronlen", default=5000000, type="int",
                      help="maximum intron length in bp (default: 5,000,000)")

    parser.add_option("--minintronlen", dest="minintronlen", default=50, type="int",
                      help="minimum intron length in bp (default: 50)")

    parser.add_option("--no_quality_control", dest="no_quality_control", default=False,
                      action="store_true",
                      help="retain pseudogenes and decay transcripts (default: filter them out)")

    parser.add_option("-v", "--virtual_intron", dest="virtual_intron", action="store_true",
                      default=False,
                      help="add virtual introns to capture alternative first/last exon usage "
                           "(experimental; default: False)")

    parser.add_option("--single_cell", dest="single_cell", default=True,
                      help="build sparse isoform→intron/exon matrices required for tealeaf-sc "
                           "(default: True)")

    (options, args) = parser.parse_args()

    if options.annot is None:
        sys.exit("Error: no annotation file provided (use -a / --annot).\n")

    sys.stderr.write(f"Processing transcriptome annotation: {options.annot}\n")
    sys.stderr.write(f"Annotation source: {options.annot_source}\n")
    sys.stderr.write(f"Output prefix: {options.outprefix}\n")
    sys.stderr.write(f"Max intron length: {options.maxintronlen}\n")
    sys.stderr.write(f"Min intron length: {options.minintronlen}\n")
    sys.stderr.write(f"Skip quality control: {options.no_quality_control}\n")
    sys.stderr.write(f"Virtual introns: {options.virtual_intron}\n")

    record = f'{options.outprefix}map_parameters.txt'
    sys.stderr.write(f'Saving parameters to {record}\n')
    write_options_to_file(options, record)

    try:
        tealeaf_map_generation(options)
    except FileNotFoundError as e:
        sys.exit(f"Error: {e}\n")

    sys.stderr.write('Finished building isoform-to-intron map\n')



