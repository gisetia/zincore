from io import StringIO
import subprocess as sp
import os
import pandas as pd


def filter_peaks(peak_dir: str) -> pd.DataFrame:
    ''' Filter peaks to only include those significantly lost in SEPHS1KO 
        and are located in promoter regions of protein coding genes

    Parameters:
        peak_dir (str): Path to the directory containing the refined
                                peaks
    Returns:
        pd.DataFrame:         DataFrame containing with filtered peaks
    '''

    # Load the differential binding peaks and annotation data
    diff_peaks = pd.read_csv(f'{peak_dir}/db/diffbind_peaks.csv')
    annot_peaks = pd.read_csv(f'{peak_dir}/annotations/annot_peaks.csv')

    merged_peaks = pd.merge(diff_peaks, annot_peaks, how='outer',
                            on=['peak_id', 'Chr', 'Start', 'End'])

    # Filter peaks
    filt_peaks = merged_peaks.query('`SEPHS1_KO-vs-WT-fold` < -1 and '
                                    '`SEPHS1_KO-vs-WT-fdr` < 0.05 and '
                                    '`annotation` == "Promoter (<=1kb)" and '
                                    '`transcriptBiotype` == "protein_coding"')

    return filt_peaks[['Chr', 'Start', 'End', 'peak_id']]


def merge_promoters(promoters_gene: pd.DataFrame) -> pd.DataFrame:
    """
    Merges overlapping promoter regions for a given gene using bedtools.

    Parameters:
    -----------
    promoters_gene : pandas.DataFrame
        A DataFrame containing promoter regions for a single gene. 
        It is expected to have the following columns:
        - 'SYMBOL': Gene symbol.
        - Additional columns representing promoter region information.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the merged promoter regions. If there is only
        one promoter region for the gene, the original DataFrame is returned
        with an additional 'gene_SYMBOL' column. If multiple promoter regions
        exist, they are merged using bedtools, and the resulting DataFrame
        contains the following columns:
        - 'Chr': Chromosome.
        - 'Start': Start position of the merged region.
        - 'End': End position of the merged region.
        - 'gene_SYMBOL': Gene symbol.

        - 'strand': Strand information.
        - 'gene_id': Gene ID.
        - 'tx_id': Transcript ID.
    """

    # Extract the gene symbol from the first row
    gene = promoters_gene['SYMBOL'].iloc[0]

    if len(promoters_gene) > 1:  # Check if there are multiple promoter regions for the gene
        # Save the promoter regions to a temporary BED file
        promoters_gene.to_csv(f'tmp/promoters_{gene}.bed',  sep='\t',
                              index=False, header=False)
        # Use bedtools to merge overlapping regions and capture the output
        cmd = (f'bedtools merge -i tmp/promoters_{gene}.bed '
               '-c 4,5,6,7 -o distinct')
        # Execute the bedtools command
        merged = sp.check_output(cmd, shell=True)
        # Read the merged output into a DataFrame
        merged = pd.read_csv(StringIO(str(merged, 'utf-8')),
                             sep='\t', header=None)
        # Assign column names to the merged DataFrame
        merged.columns = ['Chr', 'Start', 'End', 'gene_SYMBOL', 'strand',
                          'gene_id', 'tx_id']
        # Remove the temporary BED file
        os.remove(f'tmp/promoters_{gene}.bed')
    else:
        merged = promoters_gene  # If only one region, use the original DataFrame
        merged['gene_SYMBOL'] = gene  # Add the gene symbol column

    return merged


def generate_control_regions(filtered_peaks: pd.DataFrame,
                             refined_peak_dir: str) -> pd.DataFrame:
    """
    Generate control sequences by processing filtered peaks and promoter regions.
    This function takes a DataFrame of filtered peaks, processes it to generate
    promoter regions based on transcription start sites (TSS), and applies
    filtering criteria to refine the annotations. The resulting promoter regions
    are merged, assigned unique region IDs, and saved to a BED file.
    Args:
        filtered_peaks (pd.DataFrame): A DataFrame containing filtered peaks with
            columns ['Chr', 'Start', 'End', 'peak_id'].
    Returns:
        pd.DataFrame: A DataFrame containing merged promoter regions with columns:
            ['Chr', 'Start', 'End', 'region_id', 'strand', 'gene_id', 'tx_id'].
    Notes:
        - The function assumes the presence of an annotation database file
          located at '{refined_peak_dir}/annotation_dbs/edb_full.csv'.
        - Promoter regions are calculated using padding values for upstream
          and downstream regions.
        - Overlapping promoter regions are merged, and unique region IDs are
          generated for each promoter region.
        - Intermediate results, such as filtered peaks and promoters, are saved
          to temporary BED files in the 'tmp/' directory.
    """

    # Ensure the 'tmp' directory exists
    if not os.path.exists('tmp'):
        os.mkdir('tmp')

    # Save filtered peaks to a BED file
    peaks_bed = 'tmp/filtered_peaks.bed'
    filtered_peaks[['Chr', 'Start', 'End', 'peak_id']].to_csv(
        peaks_bed, sep='\t', index=False, header=False)

    # Define padding values for upstream and downstream regions
    upstream_pad = 800
    downstream_pad = 200

    # Define filtering criteria for transcript annotations
    tx_biotype = 'protein_coding'
    tx_canonical = 1
    tx_support_level = 1

    # Load the annotation database
    annots_db = pd.read_csv(f'{refined_peak_dir}/annotation_dbs/edb_full.csv')

    # Filter annotations to include only valid chromosomes
    chroms = [str(x) for x in range(1, 23)] + ['X', 'Y']
    annots_db = annots_db[annots_db['Chr'].astype(str).isin(chroms)]

    # Apply additional filtering criteria for transcripts
    annots_db = annots_db.query(
        f'tx_biotype == "{tx_biotype}" '
        f' and (tx_canonical == {tx_canonical}'
        f' and tx_support_level <= {tx_support_level})')

    # Calculate the Transcription Start Site (TSS) based on strand orientation
    annots_db['TSS'] = annots_db.apply(lambda x: x['Start'] if x['strand']
                                       == 1 else x['End'], axis=1)

    # Extract promoter regions with relevant columns
    promoters = annots_db[['Chr', 'SYMBOL', 'gene_id', 'tx_id', 'TSS',
                           'strand']]

    # Calculate the start and end positions of promoter regions based on padding
    promoters['Start'] = promoters.apply(lambda x: x['TSS'] - upstream_pad if
                                         x['strand'] == 1 else x['TSS']
                                         - downstream_pad, axis=1)
    promoters['End'] = promoters.apply(lambda x: x['TSS'] + downstream_pad if
                                       x['strand'] == 1 else x['TSS']
                                       + upstream_pad,  axis=1)

    # Sort promoters by chromosome and start position, and remove duplicates
    promoters = promoters.sort_values(['Chr', 'Start']).drop_duplicates(
        subset=['Chr', 'Start', 'End'])

    # Fill missing SYMBOL values with the corresponding gene_id
    promoters['SYMBOL'] = promoters['SYMBOL'].fillna(promoters['gene_id'])

    # Select relevant columns for further processing
    promoters = promoters[['Chr', 'Start', 'End', 'SYMBOL', 'gene_id',
                           'tx_id', 'TSS', 'strand']]

    # Merge overlapping promoter regions for each gene
    merged = promoters.groupby('SYMBOL').apply(merge_promoters)

    # Generate unique region IDs for each promoter region
    merged['region_id'] = merged.groupby('gene_SYMBOL')['Chr'].transform(
        lambda x: x.reset_index().index + 1)
    merged['region_id'] = merged.apply(lambda x:
                                       f'{x["gene_SYMBOL"]}_{x["region_id"]}',
                                       axis=1)

    # Ensure chromosome column is of string type
    merged = merged.astype({'Chr': str})

    # Reorder columns and sort by chromosome and start position
    merged = merged.reindex(columns=['Chr', 'Start', 'End', 'region_id',
                                     'strand', 'gene_id', 'tx_id']).sort_values(
        ['Chr', 'Start'])

    # Save the merged promoters to a BED file
    promoters_bed = 'tmp/promoters.bed'
    merged.to_csv(promoters_bed, sep='\t', index=False, header=False)

    # Specify the reference genome file
    ref_genome = ('resources/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai')

    seed = 0
    # Construct the bedtools shuffle command to generate control regions
    cmd = (f'bedtools shuffle -i {peaks_bed} '
           f'-excl {peaks_bed} '
           f'-incl {promoters_bed} '
           f'-g {ref_genome} '
           f'-seed {seed} '
           '-noOverlapping '
           f'-chrom')

    # Execute the bedtools shuffle command
    shuffled = sp.check_output(cmd, shell=True)

    # Read the shuffled output into a DataFrame
    shuffled = pd.read_csv(StringIO(str(shuffled, 'utf-8')),
                           sep='\t', header=None)
    shuffled.columns = ['Chr', 'Start', 'End', 'peak_id']
    shuffled.to_csv('tmp/shuffled_peaks.bed', sep='\t',
                    index=False, header=False)

    return shuffled

def get_sequences(regions_dict) -> None:

    files_dict = dict()

    ref_genome = ('resources/Homo_sapiens.GRCh38.dna.primary_assembly.fa')
    for name, regions in regions_dict.items():

        # Create a temporary file for the peaks
        regions_bed = f'tmp/regions_{name}.bed'
        regions_fasta = f'tmp/regions_{name}.fasta'
        regions.to_csv(regions_bed, sep='\t', index=False, header=False)

        cmd = (f'bedtools getfasta -fi {ref_genome} -bed {regions_bed} '
        f'-fo {regions_fasta} -name')

        sp.run(cmd, shell=True)

        files_dict[name] = regions_fasta

    return files_dict


def compare_counts(counts_dict: dict) -> pd.DataFrame:
    """
    Compare the counts of k-mers in the filtered peaks and control regions.

    Parameters:
    -----------
    counts_dict : dict
        A dictionary containing DataFrames with k-mer counts for each set of
        sequences. The keys are the names of the sets (e.g., 'peaks', 'control').

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the k-mer counts for each set of sequences.
        The DataFrame is indexed by k-mer sequences and contains the following
        columns:
        - 'peaks': Count of k-mers in the filtered peaks.
        - 'control': Count of k-mers in the control regions.
        - 'total': Total count of k-mers across all sets.
        - 'peak_fraction': Fraction of k-mers in the filtered peaks relative
                          to the total count.
    """

    # Merge the dataframes in counts_dict into a single dataframe
    merged_counts = pd.concat(counts_dict.values(), keys=counts_dict.keys(),
                              names=['Source', 'Index'])
    merged_counts = merged_counts.reset_index(level='Index', drop=True).T

    # Calculate the total k-mer counts and peak fraction
    merged_counts['total'] = merged_counts.sum(axis=1)
    merged_counts['peak_fraction'] = (merged_counts['peaks'] /
                                      merged_counts['total'])

    return merged_counts