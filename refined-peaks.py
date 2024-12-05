# %%
import os
import subprocess as sp

import refined_peaks.refined_peaks as rp
# %%

chipseq_dir = 'sample_data/ChIPseq'

wt_bam_dir = f'{chipseq_dir}/WT_BAM'

narrow_peaks_save_dir = f'{chipseq_dir}/WT_called_peaks'
summits_save_dir = f'{chipseq_dir}/WT_summits'

# Call WT peaks using MACS3
rp.call_peaks(wt_bam_dir, narrow_peaks_save_dir)

# Aggregate called peaks from all WT samples
merged_peaks = rp.merge_narrow_peaks(narrow_peaks_save_dir, summits_save_dir)

# Find summits (ie, local maxima) within the merged peaks
summits_bed = rp.find_peak_summits(wt_bam_dir, merged_peaks, summits_save_dir)

# Expand the summits by a specified width in base pairs
refined_peaks = rp.expand_bed(summits_bed, summits_save_dir, width=100)

# Perfom differential binding analysis using DiffBind
# Requires samples to be defined in <chipseq_dir>/diffbind/sample_sheet.csv
# following the format required by the DiffBind package: 
# SampleID,Factor,Replicate,bamReads,Peaks,PeakCaller
sp.run(['Rscript', 'Rscripts/differential_binding.R'])

# %%
