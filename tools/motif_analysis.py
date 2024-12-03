import numpy as np
import pandas as pd
from itertools import product
from Bio.Seq import Seq
from Bio.Seq import IUPACData

from .parallel import func_parallel


def extend_ambiguous_dna(seq):
    """return list of all possible sequences given an ambiguous DNA input"""
    d = IUPACData.ambiguous_dna_values
    return list(map("".join, product(*map(d.get, seq))))


def reverse_complement(motif_list):
    rev_dict = dict()
    for motif in motif_list:
        rev_dict[motif] = str(Seq(motif).reverse_complement())

    return rev_dict


def get_rev_comp(motif):
    '''Get sorted lists with given motifs and their reverse complement'''
    # if type(motifs) == str:
    #     motifs = [motifs]

    # print(len(motifs))

    # for motif in motifs:
    rev_comp = str(Seq(motif).reverse_complement())

    first_motif = sorted([motif, rev_comp])[0]
    second_motif = sorted([motif, rev_comp])[1]

    return first_motif, second_motif


def get_motif_counts(motif, seqs, motif_revcomp):
    '''Get the number of sequences with at least one ocurrence of a motif
    seqs: list of ATCG sequences
    motif: ATCG motif
    motif_revcomp: dictionary with motif as key and reverse complement as value
    '''

    print(motif)

    motifs = [motif, motif_revcomp[motif]]

    # Get array with counts for both motif and reverse complement for each
    # query sequence
    # raw_counts = sum([np.char.count(seqs, motif) for motif in motifs])
    raw_counts = sum([np.char.count(seqs['seq'].to_list(), motif)
                     for motif in motifs])
    seqs_with_counts = list(np.nonzero(raw_counts)[0])
    seqs_with_counts = [seqs['seq_id'].iloc[x] for x in seqs_with_counts]

    return seqs_with_counts


def read_fasta(fasta_file):
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


def count_kmers_in_seqs(fasta_file, length, nucs='ACTG', cores=1):

    print(f'-- Getting all possible {length}-mers with {nucs}')

    motifs = list(product(nucs, repeat=length))
    motifs = [''.join(x) for x in motifs]

    print('-- K-mers found:', len(motifs))

    print('-- Finding unique k-mers')
    motif_revcomp = func_parallel(get_rev_comp, motifs,
                                  cores=cores, return_dict=False,
                                  chunk_size=100000)
    motif_revcomp = {motif[0]: motif[1] for motif in motif_revcomp}
    unique_motifs = set(motif_revcomp.keys())

    print('-- Unique k-mers found:', len(unique_motifs))

    print('-- Reading fasta file')
    seqs = read_fasta(fasta_file)

    print('-- Counting motifs')
    counts = func_parallel(get_motif_counts, unique_motifs,
                        #    seqs=seqs['seq'].to_list(),
                           seqs=seqs,
                           motif_revcomp=motif_revcomp, cores=cores)

    print('-- Creating dataframe with counts')
    
    # Create dataframe with counts for each motif and each sequence
    # seq_ind = set(val for sublist in counts.values() for val in sublist)

    # counts_df = pd.DataFrame(index=seq_ind, columns=counts.keys())
    # for index, cols in counts.items():
    #     counts_df.loc[cols, index] = 1
    # counts_df = counts_df.fillna(0)

    # seq_ids = seqs['seq_id'].to_dict()
    # counts_df.index = [seq_ids[x] for x in counts_df.index]
    # counts_df.index.name = 'seq_id'
    # counts_df.columns.name = 'motif'

    counts_df = pd.DataFrame([pd.Series(counts)], index=['seq_id'])
    # counts_df.loc['count'] = counts_df.apply(lambda x: len(x.loc['seq_id']))

    # # Create dataframe with counts for each motif and each sequence
    # seq_ind = set(val for sublist in counts.values() for val in sublist)

    # max_cells = 5e8
    # # max_cells = 1e5
    # col_num = int(max_cells/len(seq_ind))


    # def break_into_sublists(lst, length):
    #     return [lst[i:i+length] for i in range(0, len(lst), length)]

    # motifs_list = break_into_sublists(list(counts.keys()), col_num)

    # counts_list = list()
    # for idx, motifs in enumerate(motifs_list):
    #     print('-- Processing batch:', idx+1, 'of', len(motifs_list))
    #     counts_df = pd.DataFrame(index=seq_ind, columns=motifs)
    #     for key, val in counts.items():
    #         if key in motifs:
    #             counts_df.loc[val, key] = 1
    #     counts_df = counts_df.fillna(0)
    #     counts_list.append(counts_df)

    # counts_df = pd.concat(counts_list, axis=1)

    # seq_ids = seqs['seq_id'].to_dict()
    # counts_df.index = [seq_ids[x] for x in counts_df.index]
    # counts_df.index.name = 'seq_id'
    # counts_df.columns.name = 'motif'

    return counts_df
    # return counts, seqs