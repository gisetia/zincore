#%%
import refined_peaks.refined_peaks as rp


wt_bam_dir = 'sample_data/ChIPseq/WT_BAM'

narrow_peaks_save_dir = 'sample_data/ChIPseq/WT_called_peaks'
summits_save_dir = 'sample_data/ChIPseq/WT_summits'

# Call WT peaks using MACS3
rp.call_peaks(wt_bam_dir, narrow_peaks_save_dir)

# Aggregate called peaks from all WT samples
merged_peaks = rp.merge_narrow_peaks(narrow_peaks_save_dir, summits_save_dir)

# Find summits (ie, local maxima) within the merged peaks
summits_bed = rp.find_peak_summits(wt_bam_dir, merged_peaks, summits_save_dir)

# Expand the summits by a specified width in base pairs
refined_peaks = rp.expand_bed(summits_bed, summits_save_dir, width=100)


# %%
