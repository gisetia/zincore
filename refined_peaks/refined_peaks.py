from typing import Tuple, List
import os
import pysam
import pandas as pd
import numpy as np
import subprocess as sp

from typing import Optional


def call_peaks(wt_bam_dir: str, save_dir: str,
               input_file: Optional[bool] = False,
               ctrl_bam_dir: Optional[str] = None) -> None:
    ''' Call peaks using MACS3 for ChIP-seq data

    Parameters:
        wt_bam_dir (str):   Path to the directory containing the wild-type
                            ChIP-seq BAM files
        save_dir (str):     Path to the directory where the called peaks will 
                            be saved
        input_file (bool):  Whether to use an input control file for peak 
                            calling (default is False)
        ctrl_bam_dir (str): Path to the directory containing the control
                            ChIP-seq BAM files (default is None)   
    '''

    # If input_file is True, set the control BAM directory
    if input_file:
        ctrl_bam_dir = 'sample_data/ChIPseq/ctrl_BAM'

    # Define the path to the MACS3 virtual environment
    macs_env_path = '/Users/gis/Documents/venvs/macs3_env/bin'
    macs_env = os.environ.copy()
    macs_env['PATH'] += os.pathsep + macs_env_path

    # Get a sorted list of BAM files in the wild-type BAM directory
    wt_bam_files = sorted([x for x in os.listdir(wt_bam_dir) if
                           x.endswith('.bam')])
    
    # If input_file is True, get the first BAM file in the control BAM directory
    if input_file:
        ctrl_bam_file = [x for x in os.listdir(ctrl_bam_dir) if
                         x.endswith('.bam')][0]

    # Create the save directory if it does not exist
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    # Loop through each BAM file and call peaks using MACS3
    for bam in wt_bam_files:
        bam_name = bam.split('.')[0]

        print('Processing:', bam_name)

        # Construct the MACS3 command
        macs_cmd = (f'macs3 callpeak '
                    f'-t {wt_bam_dir}/{bam} '
                    f'-f BAMPE '
                    f'--outdir {save_dir} '
                    f'--name {bam_name} '
                    )

        # Add the control BAM file to the command if input_file is True
        if input_file:
            macs_cmd += f'-c {ctrl_bam_dir}/{ctrl_bam_file}'

        # Execute the MACS3 command
        sp.call(macs_cmd, shell=True, env=macs_env)

    print('Done: peak calling...')


def merge_narrow_peaks(narrow_peaks_dir: str,
                       summits_save_dir: str) -> pd.DataFrame:
    ''' Merge called narrow peaks from WT samples

    parameters: narrow_peaks_dir (str): Path to the directory containing the
                                        called narrow peaks
                summits_save_dir (str): Path to the directory where the merged 
                                        peaks will be saved
            '''

    # Create the summits save directory if it does not exist
    if not os.path.exists(summits_save_dir):
        os.mkdir(summits_save_dir)

    # Get a list of all narrowPeak files in the narrow peaks directory
    peak_files = [x for x in os.listdir(narrow_peaks_dir) if
                  x.endswith('.narrowPeak')]

    # Initialize a list to store called peaks
    called_peaks = list()
    for peak_file in peak_files:
        # Read each narrowPeak file and extract relevant columns
        peaks = pd.read_csv(f'{narrow_peaks_dir}/{peak_file}', sep='\t',
                            header=None)
        called_peaks.append(peaks[[0, 1, 2]])

    # Concatenate all called peaks into a single DataFrame
    called_peaks = pd.concat(called_peaks)

    # Sort the peaks by chromosome and start position, then save to a temporary file
    called_peaks.sort_values([0, 1]).to_csv(f'{summits_save_dir}/called_peaks.bed',
                                            sep='\t', index=False, header=None)

    # Use bedtools to merge overlapping peaks and save the result
    cmd = (f'bedtools merge -i {summits_save_dir}/called_peaks.bed > '
           f'{summits_save_dir}/WT_called_peaks_merged.bed')
    sp.run(cmd, shell=True)

    # Remove the temporary file
    os.remove(f'{summits_save_dir}/called_peaks.bed')

    # Read the merged peaks file into a DataFrame
    bed_file = pd.read_csv(f'{summits_save_dir}/WT_called_peaks_merged.bed',
                           sep='\t', header=None, names=['chrom', 'start',
                                                         'end']
                           ).sort_values(['chrom', 'start'])

    # Add a peak ID column to the DataFrame
    bed_file['peak_id'] = bed_file.index + 1

    # Save the final merged peaks file
    bed_file.to_csv(f'{summits_save_dir}/WT_called_peaks_merged.bed', sep='\t',
                    index=False)

    # Return the merged peaks DataFrame
    return bed_file


def find_max(arr: np.ndarray, N: int = 30, min_counts: int = 25,
             n: int = 30) -> Tuple[np.ndarray, np.ndarray]:
    '''
    Find local maxima in a 1D array using a two-step process:
    1. Smooth the array using a moving average with a window size of N  
    2. Identify local maxima by finding the negative peaks in the second 
    derivative of the smoothed array

    Parameters: 
        arr (np.ndarray): 1D numpy array of values
        N (int): Window size for smoothing the array (default is 30)
        min_counts (int): Minimum value for a peak (default is 25)
        n (int): Window size for smoothing the first derivative (default is 30)
    '''

    # Replace values in arr that are less than or equal to min_counts with 0
    arr = np.where(arr <= min_counts, 0, arr)

    # Pad the array with its edge values to handle boundary conditions
    padded_arr = np.pad(arr, (int(N/2), int(N/2)), 'constant',
                        constant_values=(arr[0], arr[-1]))

    # Calculate the cumulative sum of the padded array
    cumsum_vec = np.cumsum(padded_arr)

    # Smooth the array using a moving average with a window size of N
    smooth_arr = (cumsum_vec[N:] - cumsum_vec[:-N]) / N

    # Compute the first derivative of the smoothed array
    diff1 = np.diff(smooth_arr)

    # Pad the first derivative and smooth it using a moving average with a
    # window size of n
    padded_diff1 = np.pad(diff1, (int(n/2), int(n/2)), 'constant',
                          constant_values=(diff1[0], diff1[-1]))
    cumsum_diff1 = np.cumsum(padded_diff1)
    diff1 = (cumsum_diff1[n:] - cumsum_diff1[:-n]) / n

    # Calculate the sign of the smoothed first derivative
    sign_diff1 = np.sign(diff1)
    sign_diff1 = np.where(sign_diff1 == 0, 1, sign_diff1)

    # Compute the second derivative by taking the difference of the sign
    # of the first derivative
    diff = np.diff(sign_diff1)

    # Insert a zero at the beginning of the sign_diff1 array
    sign_diff1 = np.insert(sign_diff1, 0, 0)

    # Pad the second derivative
    diff = np.pad(diff, (1, 1), 'constant', constant_values=(0, 0))

    # Identify indices where the second derivative is negative, indicating local maxima
    maxima = np.where(diff < 0)[0]

    # Return the indices of the local maxima and the smoothed array
    return maxima, smooth_arr


def find_peak_summits(wt_bam_dir: str, merged_peaks: pd.DataFrame,
                      summits_save_dir: str) -> pd.DataFrame:
    ''' Find summits (ie, local maxima) within each called peak in the WT 
    samples

    Parameters:
        wt_bam_dir (str): Path to the directory containing the wild-type
                            ChIP-seq BAM files
    '''

    # Get a sorted list of BAM files in the wild-type BAM directory
    wt_bam_files = sorted([x for x in os.listdir(wt_bam_dir) if
                           x.endswith('.bam')])

    # Initialize a list to store coverage data for each sample
    coverage_sample = []
    for bam_file in wt_bam_files:
        bam_name = bam_file.split('.')[0]

        print('Processing:', bam_name)

        # Dictionary to store coverage regions for the current BAM file
        coverage_regions = {}

        # Open the BAM file for reading
        with pysam.AlignmentFile(f'{wt_bam_dir}/{bam_file}', 'rb') as bam:

            # Iterate through each peak in the merged peaks DataFrame
            for idx, row in merged_peaks.iterrows():
                if idx % 5000 == 0:
                    print(f'-- Processing peak coverage {idx+1} of '
                          f'{len(merged_peaks)}')
                    
                # Calculate coverage for the current peak region
                coverage = bam.count_coverage(str(row['chrom']), row['start'],
                                              row['end'])

                # Sum the coverage across all strands
                coverage = np.sum(coverage, axis=0)
                coverage_regions.update({row['peak_id']: coverage})

        # Convert coverage regions to a pandas Series
        coverage_regions = pd.Series(coverage_regions.values(),
                                     index=coverage_regions.keys(),
                                     name=bam_file.split('/')[-1])

        # Append the coverage data for the current sample
        coverage_sample.append(coverage_regions)

    # Combine coverage data across all samples and calculate total coverage
    coverage_sample = pd.concat(coverage_sample, axis=1).sum(
        axis=1).to_frame(name='WT_coverage')
    # Join the coverage data with the merged peaks DataFrame
    coverage_sample = merged_peaks.set_index('peak_id').join(coverage_sample)

    # Find local maxima and smoothed coverage for each peak
    coverage_sample['WT_max'] = coverage_sample['WT_coverage'].apply(
        find_max)
    coverage_sample['WT_smooth'] = [x[1] for x in coverage_sample['WT_max']]
    coverage_sample['WT_max'] = [x[0] for x in coverage_sample['WT_max']]
    # Extract the maximum value for each peak
    coverage_sample['WT_max_val'] = coverage_sample.apply(
        lambda x: x['WT_smooth'][x['WT_max']], axis=1)
    # Calculate the genomic position of the local maxima
    coverage_sample['WT_max'] = (coverage_sample['start']
                                 + coverage_sample['WT_max'])

    # Create a DataFrame for the summit positions
    max_bed = coverage_sample.set_index('chrom', append=True)[
        'WT_max'].explode().to_frame().assign(
        summit_id=lambda df: df.groupby(level=0).cumcount() + 1).rename(
        columns={'WT_max': 'start'}).reset_index()

    # Add the maximum value and end position for each summit
    max_bed['max_val'] = coverage_sample.set_index('chrom', append=True)[
        'WT_max_val'].explode().values
    max_bed['end'] = max_bed['start'] + 1
    # Create a unique summit ID
    max_bed['summit_id'] = (max_bed['peak_id'].astype(str) + '_'
                            + max_bed['summit_id'].astype(str))
    # Reorder and clean up the columns
    max_bed = max_bed[['chrom', 'start', 'end', 'summit_id', 'max_val',
                       'peak_id']].dropna()

    # Save the summit positions to a BED file
    max_bed.to_csv(f'{summits_save_dir}/WT_summits.bed', sep='\t', index=False,
                   header=False)

    # Return the DataFrame containing summit positions
    return max_bed


def expand_bed(bed: pd.DataFrame, save_dir: str, 
               width: Optional[int] = 100) -> pd.DataFrame:
    ''' Expand the regions in a BED file by a specified width

    Parameters:
        width (int):      The width by which to expand the peaks
        bed (pd.DataFrame):   DataFrame containing the peaks in BED format
        save_dir (str):   The directory where the expanded BED file 
                            will be saved
    '''

    # Expand the start and end positions of the peaks by half the specified width
    bed['start'] = bed['start'] - width / 2
    bed['end'] = bed['end'] + width / 2

    # Convert the start and end positions to integers
    bed['start'] = bed['start'].astype(int)
    bed['end'] = bed['end'].astype(int)

    # Save the expanded BED file to the specified directory
    bed.to_csv(f'{save_dir}/WT_refined_peaks.bed', sep='\t',
               index=False, header=False)

    return bed