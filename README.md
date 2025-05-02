# Zincore: Refined ChIP-seq peak analysis and discovery of additional zinc-finger proteins

Code related to the publication: "**Zincore, an atypical coregulator, binds zinc finger transcription factors to control gene expression**".

This repository consists of two main scripts, namely `refined-peaks.py` and `motif-discovery.py`, which cover the following analysis sections:
1. **Refined peaks**: Standard peak callers for ChIP-seq datasets do not provide enough resolution to identify the small individual peaks found in Zincore ChIP-seq data. 
`refined-peaks.py` processes BAM files from WT samples to produce refined peaks consisting of 100bp intervals centered on summits (ie, local maxima) within MACS3 called peaks. Differential binding analysis between WT and KO samples is then performed based on the refined peaks.

2. **Motif discovery**: `motif-discovery.py` takes the refined Zincore ChIP-seq peak regions (produced by `refined-peaks.py`) as input. It analyzes the DNA sequences in these regions to identify k-mers that are overrepresented relative to a control set of sequences. The control set is constructed by randomly sampling regions from the promoter regions of canonical protein-coding genes (spanning from 800 bp upstream to 200 bp downstream of annotated transcription start sites). These control regions are selected to explicitly avoid overlap with the filtered Zincore peak regions, and they are matched to the Zincore peak set in both the number of sequences and their length distribution. For each identified k-mer, an enrichment score is computed as the ratio of its frequency in the Zincore peak sequences to its frequency in the aggregate of Zincore and control sequences.

## Requirements

* Python 3.9.6
* Python libraries in `python_requirements.txt` file
* R 4.4.2
    * DiffBind 3.16.0
    * ChIPseeker 1.42.0
    * ensembldb 2.30.0
    * AnnotationHub 3.14.0
    * rg.Hs.eg.db 3.20.0
* bedtools 2.29.2
* MACS3


## Contact
[Guizela Huelsz Prince](https://www.linkedin.com/in/g-huelsz-prince/)  
[Danielle Bianchi](https://www.nki.nl/research/find-a-researcher/researchers/danielle-bianchi/)  
[Lodewyk Wessels](https://www.nki.nl/research/research-groups/lodewyk-wessels/)  
[Thijn Brummelkamp](https://www.nki.nl/research/research-groups/thijn-brummelkamp/)  
The Netherlands Cancer Institute (NKI)

https://github.com/gisetia/zincore
