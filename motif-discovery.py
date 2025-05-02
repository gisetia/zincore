# %%
import os
from matplotlib import pyplot as plt
import motif_discovery.motif_discovery as md
import tools.motif_analysis as ma

# Define dirctory with the results of the refined peaks
# and differential binding analysis performed in refined_peaks.py
refined_peak_dir = 'sample_data/ChIPseq/diffbind'

# Filter refined peaks to only include those significantly lost in SEPHS1KO
# and are located in promoter regions of protein coding genes
filtered_peaks = md.filter_peaks(refined_peak_dir)

# Generate control regions for the filtered peaks
control_regions = md.generate_control_regions(filtered_peaks, refined_peak_dir)

# Get sequences for the filtered peaks and control regions
seq_files = md.get_sequences({'peaks': filtered_peaks,
                              'control': control_regions})

# Get the counts of k-mers in each set of sequences
k = 5
counts_dict = dict()
for key, seq_file in seq_files.items():
    counts = ma.count_kmers_in_seqs(seq_file, length=k, nucs='ACTG',
                                    cores=1)
    counts_dict[key] = counts

# Compare the counts of k-mers in the filtered peaks and control regions
merged_counts = md.compare_counts(counts_dict)

# Save the merged counts to a CSV file
merged_counts.to_csv(os.path.join(refined_peak_dir, 
                                  'kmer_counts.csv'), index=True)

# Plot fraction of k-mers in the filtered peaks relative to the total count
merged_counts.plot.scatter(x='total', y='peak_fraction')
plt.xscale('log')
# %%
