import os
import pyBigWig
import pysam
import pandas as pd
import numpy as np
import subprocess as sp
from scipy.signal import argrelextrema


narrow_peaks_dir = '../narrow_peaks'

wt_narrow_peaks = pd.read_csv('wt_narrow_peaks.bed', sep='\t', header=None)


from typing import Tuple, List
import numpy as np

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