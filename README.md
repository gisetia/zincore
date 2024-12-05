# Zincore: Refined ChIP-seq peak analysis and discovery of additional zinc-finger proteins

Code related to the publication: "**Zincore, an atypical coactivator complex that controls zinc finger transcription factors**". (manuscript in preparation)

This repository consists of three main scripts, covering the following analysis sections:
1. **Refined peaks**: Standard peak callers for ChIP-seq datasets do not provide enough resolution to identify the small individual peaks found in Zincore ChIP-seq data. 
`refined-peaks.py` processes BAM files from WT samples to produce refined peaks consisting of 100bp intervals centered on summits (ie, local maxima) within MACS3 called peaks. Differential binding analysis between WT and KO samples is then performed based on the refined peaks.

2. **ChIP overlaps**: Tools for discovery of associated transcription factors through comparison with public human ChIP-seq datasets. <span style="color:red">(in progress)</span>
3. **Motif discovery**: Tools for finding enriched DNA binding motifs within Zincore ChIP-seq peaks. <span style="color:red">(in progress)</span>



## Requirements

* Python 3.9.6
* Python libraries in `python_requirements.txt` file
* R 4.4.2
    * DiffBind 3.16.0
    * <span style="color:red">ChIPseeker 1.32.1</span>
* <span style="color:red">meme 5.5.0 (??)</span>
* bedtools v2.29.2
* MACS3


## Contact
[Guizela Huelsz Prince](https://www.linkedin.com/in/g-huelsz-prince/)  
[Danielle Bianchi](https://www.nki.nl/research/find-a-researcher/researchers/danielle-bianchi/)  
[Lodewyk Wessels](https://www.nki.nl/research/research-groups/lodewyk-wessels/)  
[Thijn Brummelkamp](https://www.nki.nl/research/research-groups/thijn-brummelkamp/)  
The Netherlands Cancer Institute (NKI)
