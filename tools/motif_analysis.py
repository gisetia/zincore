import numpy as np
import pandas as pd
from itertools import product
from Bio.Seq import Seq
from Bio.Seq import IUPACData

from .parallel import func_parallel
from typing import List
from typing import Dict, List



def extend_ambiguous_dna(seq: str) -> List[str]:
    """Return list of all possible sequences given an ambiguous DNA input."""
    d = IUPACData.ambiguous_dna_values
    return list(map("".join, product(*map(d.get, seq))))


def reverse_complement(motif_list: List[str]) -> Dict[str, str]:
    """
    Generate the reverse complement of a list of DNA motifs.
    This function takes a list of DNA motifs (strings) and computes their reverse complements.
    The reverse complement of a DNA sequence is formed by reversing the sequence and replacing
    each nucleotide with its complement (A <-> T, C <-> G).
    Args:
        motif_list (List[str]): A list of DNA motifs as strings.
    Returns:
        Dict[str, str]: A dictionary where the keys are the original motifs and the values are
                        their corresponding reverse complements.
    """

    rev_dict = dict()
    for motif in motif_list:
        rev_dict[motif] = str(Seq(motif).reverse_complement())

    return rev_dict


def get_rev_comp(motif: str) -> List[str]:
    '''Get sorted lists with given motifs and their reverse complement'''

    rev_comp = str(Seq(motif).reverse_complement())

    first_motif = sorted([motif, rev_comp])[0]
    second_motif = sorted([motif, rev_comp])[1]

    return first_motif, second_motif


def get_motif_counts(motif: str, seqs: pd.DataFrame, 
                     motif_revcomp: Dict[str, str]) -> List[str]:
    '''Get the number of sequences with at least one occurrence of a motif
    seqs: DataFrame with columns 'seq_id' and 'seq'
    motif: ATCG motif
    motif_revcomp: dictionary with motif as key and reverse complement as value

    Returns a list with the sequence ids of the sequences with at least one
    occurrence of the motif or its reverse complement
    '''

    motifs = [motif, motif_revcomp[motif]]

    # Get array with counts for both motif and reverse complement for each
    # query sequence
    raw_counts = sum([np.char.count(seqs['seq'].to_list(), motif)
                     for motif in motifs])
    seqs_with_counts = list(np.nonzero(raw_counts)[0])
    seqs_with_counts = [seqs['seq_id'].iloc[x] for x in seqs_with_counts]

    return seqs_with_counts


def read_fasta(fasta_file: str) -> pd.DataFrame:
    '''Read fasta file and return a dataframe with two columns: 
    sequence ids and sequences'''

    query_fasta = pd.read_csv(fasta_file, header=None)

    even_rows = [x for x in query_fasta.index if x % 2 == 0]
    odd_rows = [x for x in query_fasta.index if x % 2 == 1]
    seqs = pd.DataFrame(query_fasta.loc[even_rows]).reset_index(drop=True)
    seqs.columns = ['seq_id']
    seqs['seq_id'] = seqs['seq_id'].apply(
        lambda x: x.lstrip('>')).astype('<U10')
    seqs['seq'] = query_fasta.loc[odd_rows].reset_index(drop=True).astype(str)

    return seqs


def count_kmers_in_seqs(fasta_file: str, length: int, nucs: str = 'ACTG', 
                        cores: int = 1) -> pd.DataFrame:
    """
    Count k-mers of a specified length in sequences from a FASTA file.
    This function generates all possible k-mers of a given length using the specified
    nucleotide alphabet, identifies unique k-mers (considering reverse complements),
    and counts their occurrences in the sequences provided in the FASTA file.
    Args:
        fasta_file (str): Path to the input FASTA file containing sequences.
        length (int): Length of the k-mers to generate and count.
        nucs (str, optional): Nucleotide alphabet to use for generating k-mers. 
                                Defaults to 'ACTG'.
        cores (int, optional): Number of CPU cores to use for parallel processing. 
                                Defaults to 1.
    Returns:
        pd.DataFrame: A DataFrame containing the counts of unique k-mers for each sequence.
                        The DataFrame is indexed by sequence IDs from the FASTA file.
    """
                        
    # Print a message indicating the start of k-mer counting
    print(f'-- Counting all {length}-mers with {nucs}')

    # Generate all possible k-mers of the specified length using the nucleotide alphabet
    motifs = list(product(nucs, repeat=length))
    motifs = [''.join(x) for x in motifs]

    # Compute reverse complements for the generated motifs
    if cores == 1:
        motif_revcomp = reverse_complement(motifs)
    else:
        motif_revcomp_list = func_parallel(get_rev_comp, motifs, cores=cores, 
                                           return_dict=False, chunk_size=100000)
        motif_revcomp = {motif[0]: motif[1] for motif in motif_revcomp_list}
    
    # Identify unique motifs (considering reverse complements)
    unique_motifs = set(motif_revcomp.keys())

    # Read sequences from the input FASTA file
    seqs = read_fasta(fasta_file)

    # Count occurrences of unique motifs in the sequences using parallel processing
    counts = func_parallel(get_motif_counts, unique_motifs, seqs=seqs, 
                           motif_revcomp=motif_revcomp, cores=cores)

    # Create a DataFrame to store the counts of motifs for each sequence
    counts_df = pd.DataFrame([pd.Series(counts)], index=['seq_id']).map(len)

    return counts_df
