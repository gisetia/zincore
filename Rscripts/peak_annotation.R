library(ensembldb)
library(AnnotationHub)
library(ChIPseeker)
library(org.Hs.eg.db)

diffbind_dir <- "sample_data/ChIPseq/diffbind"

annots_dir <- paste0(diffbind_dir, "/annotations")
annots_db_dir <- paste0(diffbind_dir, "/annotation_dbs")

dir.create(annots_dir)
dir.create(annots_db_dir)

peaks_df <- read.csv(paste(diffbind_dir, 
"/consensus/affinity_matrix.csv",
    sep = ""
))[, 1:4]

peaks <- makeGRangesFromDataFrame(peaks_df,
    keep.extra.columns = TRUE
)

ah <- AnnotationHub()
edb_full <- ah[["AH98047"]] # Ensembl version 105

annot_peaks <- annotatePeak(
    peaks,
    TxDb = edb_full,
    annoDb = "org.Hs.eg.db"
)

annot_peaks_df <- as.data.frame(annot_peaks)
print(head(annot_peaks_df))
colnames(annot_peaks_df)[1:3] <- c("Chr", "Start", "End")
annot_peaks_df <- annot_peaks_df %>%
    dplyr::relocate(peak_id, .before = Chr)
write.table(
    annot_peaks_df,
    file = paste(
        annots_dir, "/annot_peaks.csv",
        sep = ""
    ), row.names = FALSE, sep = ","
)

# Save databases
save_col_names <- c(
    "Chr", "Start", "End", "SYMBOL", "strand", "entrez_id", "gene_id", "tx_id",
    "tx_version", "tx_canonical", "gene_biotype",
    "tx_biotype", "tx_external_name", "tx_support_level"
)
bed_cols <- c(
    "Chr", "Start", "End", "SYMBOL", "entrez_id", "strand",
    "tx_version", "gene_biotype"
)

ens_cols <- c(
    "SEQNAME", "TXSEQSTART", "TXSEQEND", "SYMBOL", "SEQSTRAND",
    "ENTREZID", "GENEID", "TXID", "TXIDVERSION", "TXISCANONICAL",
    "GENEBIOTYPE", "TXBIOTYPE", "TXEXTERNALNAME", "TXSUPPORTLEVEL"
)

edb_full_data <- select(edb_full,
    columns = ens_cols, keytype = "TXID",
    keys = keys(edb_full, keytype = "TXID")
)

colnames(edb_full_data) <- save_col_names
edb_full_data <- edb_full_data[
    order(edb_full_data[, 1], edb_full_data[, 2]),
]
write.csv(
    as.data.frame(edb_full_data),
    paste(annots_db_dir, "/edb_full.csv",
        sep = ""
    ),
    row.names = FALSE
)
write.table(
    edb_full_data[bed_cols],
    file = paste(
        annots_db_dir, "/edb_full.bed",
        sep = ""
    ), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
)

# Get list of gene names
hsdb <- org.Hs.eg.db
hsdb_cols <- c("GENETYPE", "SYMBOL","ENSEMBL","ENTREZID")

hsdb_data <- select(hsdb,
    keys = keys(hsdb, keytype = "ENSEMBL"),
    columns = hsdb_cols, keytype = "ENSEMBL"
)

write.csv(
    as.data.frame(hsdb_data),
    paste(annots_db_dir, "/hsdb_names.csv",
        sep = ""
    ),
    row.names = FALSE
)