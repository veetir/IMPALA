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
cat(sprintf("[info] min RPKM threshold (--min): %s (strict '>')\n", as.character(opt$min)))
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
## VARIABLES
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

results_filt <- results$geneOutput[results$geneOutput$RPKM > min, ]
cat(sprintf("[info] rows kept after 'RPKM > %s': %d (this step drops both RPKM<=min and RPKM==NA)\n",
            as.character(min), nrow(results_filt)))

# filter for genes that have an RPKM calculated
.na_after_thresh <- sum(is.na(results_filt$RPKM))
results_filt <- results_filt[!is.na(results_filt$RPKM),] 
cat(sprintf("[info] rows removed by explicit '!is.na(RPKM)' step: %d; rows remaining: %d\n",
            .na_after_thresh, nrow(results_filt)))

# MAF filter (labels only, no drops)
cat(sprintf("[info] labeling MAF using threshold %s (no filtering occurs here)\n", as.character(maf_threshold)))
results_filt$MAF <- as.factor(ifelse(results_filt$majorAlleleFrequency > maf_threshold, paste0("MAF > ", maf_threshold), paste0("MAF < ", maf_threshold)))
results_filt$aseResults <- as.factor(ifelse(results_filt$majorAlleleFrequency > maf_threshold & results_filt$padj < 0.05, "ASE", "BAE"))

# quick label counts
.ase_n <- sum(results_filt$aseResults == "ASE", na.rm = TRUE)
.bae_n <- sum(results_filt$aseResults == "BAE", na.rm = TRUE)
cat(sprintf("[info] label counts: ASE=%d, BAE=%d\n", .ase_n, .bae_n))

# rearrange columns to a logical order
results_filt <- results_filt[,c("gene", "geneBiotype", "RPKM", "allele1IsMajor","majorAlleleFrequency", 
                                "pValueASE", "pValueHeterogeneity", "padj",
                                "significance", "MAF", "aseResults")]

# save the data frame as a table
.outfile <- paste0(out, "/MBASED_expr_gene_results.txt")
cat(sprintf("[info] writing output: %s (rows=%d, cols=%d)\n", .outfile, nrow(results_filt), ncol(results_filt)))
write.table(results_filt, .outfile, quote = F, col.names = T, row.names = F, sep = "\t")
cat("[info] done.\n")
