
## ---------------------------------------------------------------------------
## Allelic Imbalance in Expression using MBASED, part 1
## Vanessa Porter, Oct. 2021
## Modified in 2026 (github: veetir)
## ---------------------------------------------------------------------------

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
#suppressMessages(library(prob))
suppressMessages(library(tidyr))
suppressMessages(library(MBASED))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(BiocParallel))
suppressMessages(library(stats))
suppressMessages(library(tibble))

## ---------------------------------------------------------------------------
## USER FUNCTIONS
## ---------------------------------------------------------------------------

#' Initialize BiocParallel backend
#' @description Choose a BiocParallel backend based on thread count.
#' @param threads Integer number of threads.
#' @return A BiocParallelParam object.
init_bpparam <- function(threads) {
  message(sprintf("[init_bpparam] start threads=%s", threads))
  if (!is.null(threads) && threads > 1) {
    bpparam <- BiocParallel::MulticoreParam(workers = threads)
    message(sprintf("[init_bpparam] end backend=MulticoreParam workers=%s", threads))
    return(bpparam)
  }
  bpparam <- BiocParallel::SerialParam()
  message("[init_bpparam] end backend=SerialParam workers=1")
  return(bpparam)
}

#' Ensure output directory exists
#' @description Create output directory if needed and log the path.
#' @param outdir Character output directory path.
#' @return Normalized output directory path.
ensure_outdir <- function(outdir) {
  message(sprintf("[ensure_outdir] start outdir=%s", outdir))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outdir_norm <- normalizePath(outdir, mustWork = FALSE)
  message(sprintf("[ensure_outdir] end outdir=%s", outdir_norm))
  return(outdir_norm)
}

#' Read RNA TSV
#' @description Read RNA TSV and parse allele counts from AD.
#' @param path Character path to RNA TSV.
#' @return Data frame with parsed counts and variant ID.
read_rna_tsv <- function(path) {
  message(sprintf("[read_rna_tsv] start path=%s", path))
  rna_df <- read.delim(path, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
  if (ncol(rna_df) != 7) {
    stop("RNA TSV should have exactly 7 columns.")
  }
  colnames(rna_df) <- c("CHROM", "POS", "AD", "REF", "ALT", "gene", "gene_biotype")
  rna_df$variant <- paste0(rna_df$CHROM, ":", rna_df$POS)

  ad_split <- strsplit(rna_df$AD, ",", fixed = TRUE)
  rna_df$REF.COUNTS <- suppressWarnings(as.numeric(vapply(
    ad_split,
    function(x) if (length(x) >= 1) x[[1]] else NA_character_,
    FUN.VALUE = character(1)
  )))
  rna_df$ALT.COUNTS <- suppressWarnings(as.numeric(vapply(
    ad_split,
    function(x) if (length(x) >= 2) x[[2]] else NA_character_,
    FUN.VALUE = character(1)
  )))

  rows <- nrow(rna_df)
  parsed_ok <- sum(!is.na(rna_df$REF.COUNTS) & !is.na(rna_df$ALT.COUNTS))
  genes <- length(unique(rna_df$gene))
  message(sprintf("[read_rna_tsv] end rows=%d parsed_ad_ok=%d unique_genes=%d",
                  rows, parsed_ok, genes))
  return(rna_df)
}

#' Read phased VCF
#' @description Read WhatsHap phased VCF and extract GT and PS per FORMAT.
#' @param path Character VCF path (.vcf or .vcf.gz).
#' @param sample_col Integer sample column index (default 10).
#' @return Data frame with phased SNPs and phase blocks.
#' @details Filters to biallelic SNPs and phased heterozygous GT in {0|1,1|0}.
read_phased_vcf <- function(path, sample_col = 10) {
  message(sprintf("[read_phased_vcf] start path=%s sample_col=%s", path, sample_col))
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }
  on.exit(close(con))

  wh_raw <- read.delim(con, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  rows_read <- nrow(wh_raw)
  if (ncol(wh_raw) < sample_col) {
    stop("VCF has fewer than ", sample_col, " columns; expected at least one sample column.")
  }
  if (ncol(wh_raw) > sample_col) {
    message(sprintf("[read_phased_vcf] warn multiple sample columns ncol=%d using=%d",
                    ncol(wh_raw), sample_col))
  }

  wh <- wh_raw[, c(1, 2, 4, 5, 9, sample_col)]
  colnames(wh) <- c("CHROM", "POS", "REF", "ALT", "FORMAT", "SAMPLE")

  wh <- dplyr::filter(wh, nchar(REF) == 1 & nchar(ALT) == 1)
  snp_rows <- nrow(wh)

  fmt_list <- strsplit(wh$FORMAT, ":", fixed = TRUE)
  sample_list <- strsplit(wh$SAMPLE, ":", fixed = TRUE)
  get_format_value <- function(field_name) {
    idx <- lapply(fmt_list, function(x) match(field_name, x))
    return(mapply(function(v, i) {
      if (is.na(i) || i > length(v)) NA_character_ else v[[i]]
    }, sample_list, idx, USE.NAMES = FALSE))
  }

  wh$GT <- get_format_value("GT")
  wh$phaseBlock <- get_format_value("PS")

  phased_het <- wh$GT %in% c("0|1", "1|0")
  ps_valid <- !is.na(wh$phaseBlock) & wh$phaseBlock != "" & wh$phaseBlock != "."
  phased_het_count <- sum(phased_het, na.rm = TRUE)
  ps_valid_count <- sum(ps_valid, na.rm = TRUE)

  wh$variant <- paste0(wh$CHROM, ":", wh$POS)
  wh <- wh[phased_het & ps_valid, ]
  retained <- nrow(wh)
  message(sprintf("[read_phased_vcf] end rows=%d snp_rows=%d phased_het=%d ps_valid=%d retained=%d",
                  rows_read, snp_rows, phased_het_count, ps_valid_count, retained))

  if (retained == 0) {
    message("[read_phased_vcf] warn no usable phased loci after filtering GT/PS")
    stop("No usable phased loci after filtering GT and PS.")
  }

  return(wh)
}

#' Merge phase data into RNA
#' @description Map GT and phaseBlock from VCF into RNA by variant ID.
#' @param rna_df RNA data frame.
#' @param vcf_df VCF data frame.
#' @return Updated RNA data frame with GT and phaseBlock.
merge_phase_into_rna <- function(rna_df, vcf_df) {
  message(sprintf("[merge_phase_into_rna] start rna_rows=%d vcf_rows=%d",
                  nrow(rna_df), nrow(vcf_df)))
  match_idx <- match(rna_df$variant, vcf_df$variant)
  rna_df$GT <- vcf_df$GT[match_idx]
  rna_df$phaseBlock <- vcf_df$phaseBlock[match_idx]

  matched <- sum(!is.na(match_idx))
  unmatched <- nrow(rna_df) - matched
  na_gt <- sum(is.na(rna_df$GT))
  gt_01 <- sum(rna_df$GT == "0|1", na.rm = TRUE)
  gt_10 <- sum(rna_df$GT == "1|0", na.rm = TRUE)
  if (matched == 0 || matched / max(1, nrow(rna_df)) < 0.01) {
    ratio <- matched / max(1, nrow(rna_df))
    message(sprintf("[merge_phase_into_rna] warn low overlap ratio=%.4f; check chr naming (e.g. chr1 vs 1)",
                    ratio))
    message(sprintf("[merge_phase_into_rna] warn RNA example variants: %s",
                    paste(head(rna_df$variant, 2), collapse = ", ")))
    message(sprintf("[merge_phase_into_rna] warn VCF example variants: %s",
                    paste(head(vcf_df$variant, 2), collapse = ", ")))
  }

  message(sprintf("[merge_phase_into_rna] end matched=%d unmatched=%d gt_na=%d gt_0|1=%d gt_1|0=%d",
                  matched, unmatched, na_gt, gt_01, gt_10))
  return(rna_df)
}

#' Identify single unphased genes
#' @description Find genes with a single unphased variant.
#' @param rna_df RNA data frame.
#' @param phased_variants Character vector of phased variant IDs.
#' @return Data frame of single unphased loci.
identify_single_unphased_genes <- function(rna_df, phased_variants) {
  message(sprintf("[identify_single_unphased_genes] start rna_rows=%d", nrow(rna_df)))
  rna_phase <- dplyr::mutate(rna_df, phase = variant %in% phased_variants)
  gene_counts <- dplyr::summarise(
    dplyr::group_by(rna_df, gene),
    n = dplyr::n(),
    .groups = "drop"
  )
  single_unphased <- dplyr::filter(
    dplyr::left_join(rna_phase, gene_counts, by = "gene"),
    !phase & n == 1
  )

  genes_count <- length(unique(single_unphased$gene))
  loci_count <- nrow(single_unphased)
  message(sprintf("[identify_single_unphased_genes] end genes=%d loci=%d",
                  genes_count, loci_count))
  return(single_unphased)
}

#' Apply forced genotype for single unphased loci
#' @description Assign a forced GT to specified variants.
#' @param rna_df RNA data frame.
#' @param single_unphased_df Data frame of loci to update.
#' @param forced_gt Character genotype to assign.
#' @return Updated RNA data frame.
apply_single_unphased_gt <- function(rna_df, single_unphased_df, forced_gt = "1|0") {
  n_loci <- nrow(single_unphased_df)
  message(sprintf("[apply_single_unphased_gt] start forced_gt=%s loci=%d",
                  forced_gt, n_loci))
  idx <- rna_df$variant %in% single_unphased_df$variant
  rna_df$GT[idx] <- forced_gt
  message(sprintf("[apply_single_unphased_gt] end updated=%d", sum(idx)))
  return(rna_df)
}

#' Add haplotype alleles and counts
#' @description Create alleleA/alleleB and their counts based on GT.
#' @param rna_df RNA data frame.
#' @return Updated RNA data frame.
add_haplotype_alleles_and_counts <- function(rna_df) {
  message(sprintf("[add_haplotype_alleles_and_counts] start rna_rows=%d", nrow(rna_df)))
  rna_df$alleleA <- NA_character_
  rna_df$alleleB <- NA_character_
  rna_df$alleleA.counts <- NA_real_
  rna_df$alleleB.counts <- NA_real_

  idx_01 <- !is.na(rna_df$GT) & rna_df$GT == "0|1"
  idx_10 <- !is.na(rna_df$GT) & rna_df$GT == "1|0"


  rna_df$alleleA[idx_01] <- rna_df$REF[idx_01]
  rna_df$alleleB[idx_01] <- rna_df$ALT[idx_01]
  rna_df$alleleA.counts[idx_01] <- rna_df$REF.COUNTS[idx_01]
  rna_df$alleleB.counts[idx_01] <- rna_df$ALT.COUNTS[idx_01]

  rna_df$alleleA[idx_10] <- rna_df$ALT[idx_10]
  rna_df$alleleB[idx_10] <- rna_df$REF[idx_10]
  rna_df$alleleA.counts[idx_10] <- rna_df$ALT.COUNTS[idx_10]
  rna_df$alleleB.counts[idx_10] <- rna_df$REF.COUNTS[idx_10]

  n_01 <- sum(idx_01)
  n_10 <- sum(idx_10)
  n_total <- n_01 + n_10
  message(sprintf("[add_haplotype_alleles_and_counts] end gt_0|1=%d gt_1|0=%d haplotype_mapped=%d",
                  n_01, n_10, n_total))
  return(rna_df)
}

#' Filter phased loci with required fields
#' @description Retain phased loci with required non-missing fields.
#' @param rna_df RNA data frame.
#' @return Filtered data frame with phased loci.
filter_phased_complete <- function(rna_df) {
  message(sprintf("[filter_phased_complete] start rows=%d", nrow(rna_df)))
  required_cols <- c(
    "CHROM", "POS", "gene", "alleleA", "alleleB",
    "alleleA.counts", "alleleB.counts", "phaseBlock"
  )
  missing_cols <- setdiff(required_cols, colnames(rna_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  keep <- complete.cases(rna_df[, required_cols])
  phased_df <- rna_df[keep, ]
  message(sprintf("[filter_phased_complete] end retained=%d", nrow(phased_df)))
  return(phased_df)
}

#' Filter dominant phase block per gene
#' @description Retain loci from the dominant phase block per gene.
#' @param phased_df Phased loci data frame.
#' @return Filtered data frame.
#' @details Ties are resolved deterministically by phaseBlock ordering.
filter_dominant_phaseblock_per_gene <- function(phased_df) {
  message(sprintf("[filter_dominant_phaseblock_per_gene] start rows=%d genes=%d",
                  nrow(phased_df), length(unique(phased_df$gene))))
  phase_counts <- dplyr::count(phased_df, gene, phaseBlock, name = "n")
  phase_counts <- dplyr::arrange(phase_counts, gene, dplyr::desc(n), phaseBlock)

  dominant <- dplyr::slice(dplyr::group_by(phase_counts, gene), 1)
  dominant <- dplyr::ungroup(dominant)
  dominant <- dplyr::select(dominant, gene, phaseBlock)

  filtered <- dplyr::inner_join(phased_df, dominant, by = c("gene", "phaseBlock"))

  retained_loci <- nrow(filtered)
  retained_genes <- length(unique(filtered$gene))
  message(sprintf("[filter_dominant_phaseblock_per_gene] end retained_loci=%d retained_genes=%d",
                  retained_loci, retained_genes))

  if (nrow(filtered) == 0) {
    stop("No phased loci available for MBASED after filtering.")
  }

  return(filtered)
}

#' Label SNVs per gene
#' @description Add SNV labels and IDs per gene.
#' @param df Input data frame.
#' @return Data frame with SNV.ID column.
label_snvs <- function(df) {
  message(sprintf("[label_snvs] start rows=%d", nrow(df)))
  df <- dplyr::arrange(df, CHROM, POS)
  df <- dplyr::group_by(df, gene)
  df <- dplyr::mutate(df, label = paste0("SNV", dplyr::row_number()))
  df <- dplyr::ungroup(df)
  df$SNV.ID <- paste0(df$gene, ":", df$label)
  message(sprintf("[label_snvs] end loci=%d genes=%d",
                  nrow(df), length(unique(df$gene))))
  return(df)
}

#' Run MBASED
#' @description Build SummarizedExperiment and run MBASED.
#' @param df Input loci data frame.
#' @param phased Logical; TRUE for phased mode.
#' @param bpparam BiocParallelParam.
#' @param outdir Output directory path.
#' @return MBASED SummarizedExperiment object.
run_mbasesed <- function(df, phased, bpparam, outdir) {
  message(sprintf("[run_mbasesed] start phased=%s loci=%d", phased, nrow(df)))
  if (nrow(df) == 0) {
    stop("No loci available for MBASED.")
  }

  if (phased) {
    allele1 <- df$alleleA
    allele2 <- df$alleleB
    counts1 <- df$alleleA.counts
    counts2 <- df$alleleB.counts
    out_file <- file.path(outdir, "ASEresults_1s_haplotypesKnown.rds")
  } else {
    allele1 <- df$REF
    allele2 <- df$ALT
    counts1 <- df$REF.COUNTS
    counts2 <- df$ALT.COUNTS
    out_file <- file.path(outdir, "ASEresults_1s_haplotypesUnknown.rds")
  }

  mySNVs <- GenomicRanges::GRanges(
    seqnames = df$CHROM,
    ranges = IRanges::IRanges(start = df$POS, width = 1),
    aseID = df$gene,
    allele1 = allele1,
    allele2 = allele2
  )
  names(mySNVs) <- df$SNV.ID

  mySample <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      lociAllele1Counts = matrix(counts1, ncol = 1, dimnames = list(names(mySNVs), "mySample")),
      lociAllele2Counts = matrix(counts2, ncol = 1, dimnames = list(names(mySNVs), "mySample"))
    ),
    rowRanges = mySNVs
  )

  ASEresults <- MBASED::runMBASED(
    ASESummarizedExperiment = mySample,
    isPhased = phased,
    numSim = 0,
    BPPARAM = bpparam
  )

  saveRDS(ASEresults, file = out_file)
  message(sprintf("[run_mbasesed] end output=%s loci=%d", out_file, nrow(df)))
  return(ASEresults)
}

#' Postprocess MBASED results
#' @description Summarize and annotate MBASED results.
#' @param ASEresults MBASED SummarizedExperiment.
#' @param single_unphased_genes Data frame of single unphased loci.
#' @param rna_df RNA data frame.
#' @return List with gene and locus outputs.
postprocess_results <- function(ASEresults, single_unphased_genes, rna_df) {
  message("[postprocess_results] start")
  results <- summarizeASEResults_1s(ASEresults)

  results$geneOutput$padj <- stats::p.adjust(p = results$geneOutput$pValueASE, method = "BH")
  results$geneOutput$significance <- as.factor(
    ifelse(results$geneOutput$padj < 0.05, "padj < 0.05", "padj > 0.05")
  )
  results$geneOutput$gene <- rownames(results$geneOutput)

  if (nrow(single_unphased_genes) > 0) {
    results$geneOutput$allele1IsMajor[results$geneOutput$gene %in% single_unphased_genes$gene] <- NA
  }

  results$geneOutput$geneBiotype <- rna_df$gene_biotype[
    match(results$geneOutput$gene, rna_df$gene)
  ]

  n_genes <- nrow(results$geneOutput)
  n_sig <- sum(results$geneOutput$padj < 0.05, na.rm = TRUE)
  message(sprintf("[postprocess_results] end genes=%d padj_lt_0.05=%d",
                  n_genes, n_sig))
  return(results)
}

#' Summarize ASE results
#' @description Convert MBASED output into gene and locus summaries.
#' @param MBASEDOutput MBASED SummarizedExperiment object.
#' @return List with geneOutput and locusOutput.
summarizeASEResults_1s <- function(MBASEDOutput) {
  geneOutputDF <- data.frame(
    majorAlleleFrequency = SummarizedExperiment::assays(MBASEDOutput)$majorAlleleFrequency[, 1],
    pValueASE = SummarizedExperiment::assays(MBASEDOutput)$pValueASE[, 1],
    pValueHeterogeneity = SummarizedExperiment::assays(MBASEDOutput)$pValueHeterogeneity[, 1]
  )

  geneAllele <- as.data.frame(
    SummarizedExperiment::assays(S4Vectors::metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor
  )
  geneAllele <- tibble::rownames_to_column(geneAllele, var = "rowname")
  geneAllele <- dplyr::mutate(
    geneAllele,
    gene = unlist(lapply(strsplit(rowname, split = ":"), function(x) { x = x[1] }))
  )
  geneAllele <- dplyr::group_by(geneAllele, gene)
  geneAllele <- dplyr::summarise(geneAllele, allele1IsMajor = unique(mySample), .groups = "drop")

  geneOutputDF$allele1IsMajor <- geneAllele$allele1IsMajor[
    match(rownames(geneOutputDF), geneAllele$gene)
  ]

  lociOutputGR <- SummarizedExperiment::rowRanges(S4Vectors::metadata(MBASEDOutput)$locusSpecificResults)
  lociOutputGR$allele1IsMajor <- SummarizedExperiment::assays(
    S4Vectors::metadata(MBASEDOutput)$locusSpecificResults
  )$allele1IsMajor[, 1]
  lociOutputGR$MAF <- SummarizedExperiment::assays(
    S4Vectors::metadata(MBASEDOutput)$locusSpecificResults
  )$MAF[, 1]
  lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels = unique(lociOutputGR$aseID)))

  return(
    list(
      geneOutput = geneOutputDF,
      locusOutput = lociOutputList
    )
  )
}

## ---------------------------------------------------------------------------
## LOAD INPUT
## ---------------------------------------------------------------------------

# Make help options
option_list = list(
  optparse::make_option(c("-p", "--phase"), type = "character", default = NULL,
              help = "Phased VCF file (from WhatsHap)", metavar = "character"),
  optparse::make_option(c("-r", "--rna"), type = "character", default = NULL,
              help = "Tumour RNA vcf file (from Strelka2)", metavar = "character"),
  optparse::make_option(c("-o", "--outdir"), type = "character", default = "mBASED",
              help = "Output directory name", metavar = "character"),
  optparse::make_option(c("-t", "--threads"), type = "integer", default = 1,
              help = "Threads used for mbased", metavar = "integer")
)

# load in options
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

out <- ensure_outdir(opt$outdir)
bpparam <- init_bpparam(opt$threads)

## ---------------------------------------------------------------------------
## READ IN THE RNA SNV CALLS
## ---------------------------------------------------------------------------

rna_filt <- read_rna_tsv(opt$rna)

## ---------------------------------------------------------------------------
## MBASED
## ---------------------------------------------------------------------------

### WITH PHASING
# Phased mode: read WhatsHap VCF, map GT/PS into RNA, and build haplotype-aware
# alleles/counts. Filter to complete phased loci and dominant phase blocks per 
# gene, label SNVs, run MBASED (isPhased=TRUE), then summarize results.
if (!is.null(opt$phase)) {

  vcf_df <- read_phased_vcf(opt$phase, sample_col = 10)
  rna_filt <- merge_phase_into_rna(rna_filt, vcf_df)

  single_unphased <- identify_single_unphased_genes(rna_filt, vcf_df$variant)
  rna_filt <- apply_single_unphased_gt(rna_filt, single_unphased)
  rna_filt <- add_haplotype_alleles_and_counts(rna_filt)

  rna_phased <- filter_phased_complete(rna_filt)
  rna_phased <- filter_dominant_phaseblock_per_gene(rna_phased)
  rna_phased <- label_snvs(rna_phased)

  ASEresults_1s_haplotypesKnown <- run_mbasesed(
    rna_phased,
    phased = TRUE,
    bpparam = bpparam,
    outdir = out
  )

  results <- postprocess_results(ASEresults_1s_haplotypesKnown, single_unphased, rna_filt)

### WITHOUT PHASING
# Unphased mode: label SNVs, run MBASED (isPhased=FALSE), and summarize results.
# Skips phase-specific filtering and haplotype allele construction.
} else {

  rna_filt <- label_snvs(rna_filt)
  single_unphased <- data.frame(gene = character(0), variant = character(0), stringsAsFactors = FALSE)

  ASEresults_1s_haplotypesUnknown <- run_mbasesed(
    rna_filt,
    phased = FALSE,
    bpparam = bpparam,
    outdir = out
  )

  results <- postprocess_results(ASEresults_1s_haplotypesUnknown, single_unphased, rna_filt)
}

# save the results
saveRDS(results, file = file.path(out, "MBASEDresults.rds"))
message("[main] Finished MBASED")
