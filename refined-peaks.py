#%%
import refined_peaks.refined_peaks as rp


wt_bam_dir = 'sample_data/ChIPseq/WT_BAM'
save_dir = 'sample_data/ChIPseq/WT_called_peaks'

rp.call_peaks(wt_bam_dir, save_dir)
# %%

rp.