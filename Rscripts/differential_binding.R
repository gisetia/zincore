library(DiffBind)
library(dplyr)

diffbind_dir <- "sample_data/ChIPseq/diffbind"

dir.create(paste0(diffbind_dir, "/Robjects"))
dir.create(paste0(diffbind_dir, "/consensus"))
dir.create(paste0(diffbind_dir, "/db"))

samples <- read.csv(
        paste(diffbind_dir, "sample_sheet.csv", sep = "/")
    )

# Get consensus peaks
peaks <- dba(sampleSheet = samples, config = data.frame(
        AnalysisMethod = DBA_DESEQ2, doBlacklist = TRUE, doGreylist = FALSE,
        cores = 2
    ))
peaks <- dba.blacklist(peaks, blacklist = DBA_BLACKLIST_GRCH38)
dba.save(peaks,
        file = "peaks",
        dir = paste0(diffbind_dir, "/Robjects")
    )

# Get counts
counts <- dba.count(peaks, summits = FALSE)
dba.save(counts,
        file = "_counts",
        dir = paste0(diffbind_dir, "/Robjects")
    )

# Perform differential analysis
diff <- dba.normalize(counts)
diff <- dba.contrast(diff, 
        categories = DBA_FACTOR, minMembers=2
    )


diff <- dba.analyze(diff, method = DBA_DESEQ2)
dba.save(diff,
        file = "diff",
        dir = paste0(diffbind_dir, "/Robjects")
    )

# Get binding affinity matrix
aff_matrix <- dba.peakset(diff, bRetrieve = TRUE)
aff_matrix <- keepStandardChromosomes(
    aff_matrix,
    pruning.mode = "coarse"
)
aff_matrix_df <- as.data.frame(aff_matrix)
aff_matrix_df$peak_id <- seq.int(nrow(aff_matrix_df))
colnames(aff_matrix_df)[0:3] <- c("Chr", "Start", "End")
aff_matrix_df <- aff_matrix_df %>%
    dplyr::relocate(peak_id, .before = Chr)
write.table(
    aff_matrix_df,
    paste(diffbind_dir, "/consensus/affinity_matrix.csv", sep = ""),
    row.names = FALSE, sep = ",", quote = FALSE
)
write.table(
        aff_matrix_df[c("Chr", "Start", "End", "peak_id")],
        file = paste(
            diffbind_dir, "/consensus/affinity_peaks.bed",
            sep = ""
        ), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
    )

# Get called peaks
cons_peaks <- aff_matrix_df[, 1:4]
cons_peaks$idx <- rownames(cons_peaks)
called_peaks <- as.data.frame(diff$called)
called_peaks$idx <- rownames(called_peaks)
called_peaks <- merge(
    cons_peaks,
    called_peaks
)[, -c(1)]
called_peaks <- called_peaks[
    order(called_peaks$peak_id),
]
write.table(
    called_peaks,
    paste(diffbind_dir, "/consensus/called_peaks.csv", sep = ""),
    row.names = FALSE, sep = ",", quote = FALSE
)

# Report diffbind results
report <- dba.report(diff,
    method = DBA_DESEQ2, bDB = TRUE, th = 1
)

db_peaks_all <- aff_matrix_df[, 1:4]
# Merge diffbind results for all contrasts
for (i in 1:length(report$peaks)) {
    # Get diff info for contrast

    db_peaks <- report$peaks[[i]][
        c("Chr", "Start", "End", "Fold", "FDR")
    ]

    # Change column names to add contrast i
    contrast <- gsub(" ", "", gsub(
        "vs.", "-vs-", report$class[i]
    ))

    colnames(db_peaks)[4:5] <- c(
        paste(contrast, "-fold", sep = ""), paste(contrast, "-fdr", sep = "")
    )

    db_peaks_all <- merge(db_peaks_all,
        db_peaks,
        by = c("Chr", "Start", "End"), all.x = TRUE
    )
}


db_peaks_all <- db_peaks_all %>%
    dplyr::relocate(peak_id, .before = Chr)

write.table(
    db_peaks_all,
    file = paste(
        diffbind_dir, "/db/diffbind_peaks.csv",
        sep = ""
    ), row.names = FALSE, sep = ",", quote = FALSE
)