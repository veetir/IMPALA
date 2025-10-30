## ---------------------------------------------------------------------------
## ADD THE RPKM AND FILTER THE ASE GENES
## ---------------------------------------------------------------------------

suppressMessages(library(optparse))
suppressMessages(library(dplyr))

## ---------------------------------------------------------------------------
## OPTIONS
## ---------------------------------------------------------------------------

# Make help options
option_list = list(
  make_option(c("-b", "--mbased"), type="character", default=NULL,
              help="mbased rds file", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default = NULL,
              help="Sample name from the RPKM matrix (HTMCP written like e.g. HTMCP.03.06.02109)", metavar="character"),
  make_option(c("-r", "--rpkm"), type="character", default = NULL,
              help="RPKM matrix", metavar="character"),
  make_option(c("-m", "--min"), type="numeric", default = 1,
              help="Minimum RPKM value", metavar="numeric"),
  make_option(c("-t", "--maf_threshold"), type="numeric", default = 0.60,
              help="Threshold for MAF to consider as ASE", metavar="numeric"),
  make_option(c("-o", "--outdir"), type="character", default = "mBASED",
              help="Output directory name", metavar="character")
)

## ---------------------------------------------------------------------------
## VARIABLES
## ---------------------------------------------------------------------------

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

out <- opt$outdir
sample <- opt$sample
cat(sprintf("[info] mbased RDS: %s\n", opt$mbased))
cat(sprintf("[info] RPKM matrix: %s\n", opt$rpkm))
cat(sprintf("[info] sample (input): %s\n", sample))
cat(sprintf("[info] min RPKM threshold (--min): %s (reporting only; no filtering)\n", as.character(opt$min)))
cat(sprintf("[info] MAF threshold (--maf_threshold): %s\n", as.character(opt$maf_threshold)))
cat(sprintf("[info] outdir: %s\n", out))

rpkm <- read.delim(opt$rpkm, header = T, stringsAsFactors = F) 
cat(sprintf("[info] RPKM matrix dims: %d genes x %d columns\n", nrow(rpkm), ncol(rpkm)))
if (!"gene" %in% colnames(rpkm)) {
  cat("[warn] 'gene' column not found in RPKM matrix; subsequent matching will fail.\n")
}

results <- readRDS(opt$mbased)
if (!is.list(results) || !"geneOutput" %in% names(results)) {
  cat("[warn] Unexpected structure in mbased RDS; 'geneOutput' not found.\n")
} else {
  cat(sprintf("[info] genes in results$geneOutput (pre-expression merge): %d\n", nrow(results$geneOutput)))
}

min <- opt$min
maf_threshold <- opt$maf_threshold 

## ---------------------------------------------------------------------------
## Expression + Bins
## ---------------------------------------------------------------------------

print("Adding expression")

# fix sample name
sample_before <- sample
sample <- ifelse(length(grep("-", sample)) == 0, sample, gsub("-", ".", sample))
if (!identical(sample_before, sample)) {
  cat(sprintf("[info] normalized sample name: '%s' -> '%s'\n", sample_before, sample))
} else {
  cat(sprintf("[info] sample name unchanged after normalization: '%s'\n", sample))
}

# select the RPKM of this sample
if (!sample %in% colnames(rpkm)) {
  cat(sprintf("[warn] sample column '%s' not found in RPKM matrix columns.\n", sample))
}
rpkm_sample <- rpkm[,c("gene", sample)] 
cat(sprintf("[info] rpkm_sample rows: %d\n", nrow(rpkm_sample)))

# expressed genes in the sample
# (Log matching stats before assigning)
.match_idx <- match(results$geneOutput$gene, gsub(" ", "", rpkm_sample$gene, fixed = TRUE))
cat(sprintf("[info] matching genes between results$geneOutput and RPKM: matched=%d, unmatched=%d\n",
            sum(!is.na(.match_idx)), sum(is.na(.match_idx))))

results$geneOutput$RPKM <- rpkm_sample[.match_idx, 2]

# Summarize RPKM availability before thresholding
.na_before <- sum(is.na(results$geneOutput$RPKM))
.non_na_before <- sum(!is.na(results$geneOutput$RPKM))
cat(sprintf("[info] RPKM assigned: non-NA=%d, NA=%d\n", .non_na_before, .na_before))

# Distribution around the threshold (only for non-NA)
.gt_min_before <- sum(results$geneOutput$RPKM > min, na.rm = TRUE)
.le_min_before <- sum(results$geneOutput$RPKM <= min, na.rm = TRUE)
cat(sprintf("[info] RPKM > %s: %d; RPKM <= %s: %d (NAs excluded in counts)\n",
            as.character(min), .gt_min_before, as.character(min), .le_min_before))

# Tertile computation
rpkm_all <- results$geneOutput$RPKM
rpkm_pos <- rpkm_all[!is.na(rpkm_all) & rpkm_all > 0]

if (length(rpkm_pos) >= 3 && length(unique(rpkm_pos)) >= 3) {
  qs <- as.numeric(quantile(rpkm_pos, probs = c(1/3, 2/3), na.rm = TRUE, type = 8))
  q1 <- qs[1]; q2 <- qs[2]
} else {
  if (length(rpkm_pos) == 0) {
    q1 <- NA_real_; q2 <- NA_real_
  } else {
    pos_min <- min(rpkm_pos); pos_max <- max(rpkm_pos)
    q1 <- pos_min + (pos_max - pos_min)/3
    q2 <- pos_min + 2*(pos_max - pos_min)/3
  }
}
cat(sprintf("[info] expression tertiles (computed on RPKM>0): q1=%s, q2=%s\n",
            ifelse(is.na(q1), "NA", format(q1, digits=6)),
            ifelse(is.na(q2), "NA", format(q2, digits=6))))

expr_bin <- character(length(rpkm_all))
expr_bin[is.na(rpkm_all)] <- "Missing"
expr_bin[!is.na(rpkm_all) & rpkm_all == 0] <- "Zero"
expr_bin[!is.na(rpkm_all) & rpkm_all > 0 & !is.na(q1) & rpkm_all <= q1] <- "Low"
expr_bin[!is.na(rpkm_all) & rpkm_all > 0 & !is.na(q2) & rpkm_all > q1 & rpkm_all <= q2] <- "Medium"
expr_bin[!is.na(rpkm_all) & rpkm_all > 0 & !is.na(q2) & rpkm_all > q2] <- "High"
# If q1/q2 are NA (no/too few positives), mark all positives as "Low"
expr_bin[!is.na(rpkm_all) & rpkm_all > 0 & is.na(q1)] <- "Low"

results$geneOutput$expr_bin <- factor(expr_bin, levels = c("Missing","Zero","Low","Medium","High"))

# Keep all rows (no filtering by RPKM)
results_filt <- results$geneOutput
cat(sprintf("[info] rows kept without RPKM filtering: %d (RPKM zeros and NAs are retained with explicit bins)\n",
            nrow(results_filt)))

# MAF filter (labels only, no drops)
cat(sprintf("[info] labeling MAF using threshold %s (no filtering occurs here)\n", as.character(maf_threshold)))
results_filt$MAF <- as.factor(ifelse(results_filt$majorAlleleFrequency > maf_threshold, paste0("MAF > ", maf_threshold), paste0("MAF < ", maf_threshold)))
results_filt$aseResults <- as.factor(ifelse(results_filt$majorAlleleFrequency > maf_threshold & results_filt$padj < 0.05, "ASE", "BAE"))

# quick label counts
.ase_n <- sum(results_filt$aseResults == "ASE", na.rm = TRUE)
.bae_n <- sum(results_filt$aseResults == "BAE", na.rm = TRUE)
cat(sprintf("[info] label counts: ASE=%d, BAE=%d\n", .ase_n, .bae_n))

# rearrange columns to a logical order
results_filt <- results_filt[,c("gene", "geneBiotype", "RPKM", "expr_bin", "allele1IsMajor","majorAlleleFrequency", 
                                "pValueASE", "pValueHeterogeneity", "padj",
                                "significance", "MAF", "aseResults")]

# save the data frame as a table
.outfile <- paste0(out, "/MBASED_expr_gene_results.txt")
cat(sprintf("[info] writing output: %s (rows=%d, cols=%d)\n", .outfile, nrow(results_filt), ncol(results_filt)))
write.table(results_filt, .outfile, quote = F, col.names = T, row.names = F, sep = "\t")
cat("[info] done.\n")
