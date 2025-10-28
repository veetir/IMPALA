#!/usr/bin/env Rscript

## ---------------------------------------------------------------------------
## MBASED (phased only) — single-threaded driver for phased VCF + RNA TSV
## ---------------------------------------------------------------------------

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(MBASED))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(GenomicRanges))
suppressMessages(library(IRanges))
suppressMessages(library(BiocParallel))

option_list <- list(
  make_option(c("-p","--phase"),   type="character", help="Phased VCF(.gz), single-sample"),
  make_option(c("-r","--rna"),     type="character", help="RNA SNV TSV"),
  make_option(c("-o","--outdir"),  type="character", default="mBASED", help="Output dir"),
  make_option(c("-t","--threads"), type="integer",   default=1, help="Number of parallel workers for MBASED"),
  make_option(c("--num-sim"),      type="integer",   default=1e6, help="MBASED numSim (default 1e6)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$phase) || is.null(opt$rna)) {
  stop("Missing required --phase and/or --rna")
}
outdir <- opt$outdir
if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

threads <- suppressWarnings(as.integer(opt$threads))
if (is.na(threads) || threads < 1L) threads <- 1L
cat(sprintf("[info] threads requested=%d\n", threads))

numSim <- suppressWarnings(as.integer(opt$`num-sim`))
if (is.na(numSim)) numSim <- 1e6L
cat(sprintf("[info] numSim=%d\n", numSim))

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

bp <- if (threads > 1L) {
  BiocParallel::MulticoreParam(workers = threads)
} else {
  BiocParallel::SerialParam()
}

## helpers
list_n_item <- function(lst, n) sapply(lst, function(x) if (length(x) >= n) x[[n]] else NA)
to_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x, perl=TRUE), x, paste0("chr", x))
}

summarizeASEResults_1s <- function(MBASEDOutput) {
  geneOutputDF <- data.frame(
    majorAlleleFrequency = assays(MBASEDOutput)$majorAlleleFrequency[,1],
    pValueASE            = assays(MBASEDOutput)$pValueASE[,1],
    pValueHeterogeneity  = assays(MBASEDOutput)$pValueHeterogeneity[,1]
  )
  geneAllele <- assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="rowname") %>%
    mutate(gene = vapply(strsplit(rowname, ":"), `[`, "", 1L)) %>%
    group_by(gene) %>%
    summarise(allele1IsMajor = unique(mySample), .groups="drop")
  geneOutputDF$gene <- rownames(geneOutputDF)
  geneOutputDF$allele1IsMajor <- geneAllele$allele1IsMajor[match(geneOutputDF$gene, geneAllele$gene)]
  lociOutputGR <- rowRanges(metadata(MBASEDOutput)$locusSpecificResults)
  lociOutputGR$allele1IsMajor <- assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[,1]
  lociOutputGR$MAF            <- assays(metadata(MBASEDOutput)$locusSpecificResults)$MAF[,1]
  list(geneOutput = geneOutputDF,
       locusOutput = split(lociOutputGR, factor(lociOutputGR$aseID)))
}

## ---------------------------------------------------------------------------
## READ RNA TSV
## ---------------------------------------------------------------------------
rna <- tryCatch(
  read.delim(opt$rna, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE),
  error=function(e) stop("Failed to read RNA TSV: ", conditionMessage(e))
)
need_cols <- c("CHROM","POS","GEN[0].AD","REF","ALT","ANN[0].GENE","ANN[0].BIOTYPE")
if (!all(need_cols %in% colnames(rna))) {
  stop("RNA TSV missing columns: ", paste(setdiff(need_cols, colnames(rna)), collapse=", "))
}
colnames(rna)[match(need_cols, colnames(rna))] <- c("CHROM","POS","AD","REF","ALT","gene","gene_biotype")

rna$CHROM   <- to_chr(rna$CHROM)
rna$POS     <- suppressWarnings(as.integer(rna$POS))
rna <- rna[!is.na(rna$POS), , drop=FALSE]
rna$variant <- paste0(rna$CHROM, ":", rna$POS)

ad <- strsplit(rna$AD, ",", fixed=TRUE)
rna$REF.COUNTS <- suppressWarnings(as.numeric(list_n_item(ad, 1)))
rna$ALT.COUNTS <- suppressWarnings(as.numeric(list_n_item(ad, 2)))
rna$REF.COUNTS[is.na(rna$REF.COUNTS)] <- 0
rna$ALT.COUNTS[is.na(rna$ALT.COUNTS)] <- 0

# basic QC
rna <- rna[is.finite(rna$REF.COUNTS) & is.finite(rna$ALT.COUNTS) &
             rna$REF.COUNTS >= 0 & rna$ALT.COUNTS >= 0, , drop=FALSE]

cat(sprintf("[info] RNA rows after basic QC: %d\n", nrow(rna)))

## ---------------------------------------------------------------------------
## READ PHASED VCF (first 10 columns)
## ---------------------------------------------------------------------------
vcf <- tryCatch({
  con <- gzfile(opt$phase, open="rt")
  on.exit(close(con), add=TRUE)
  df <- read.delim(con, header=FALSE, comment.char="#",
                   stringsAsFactors=FALSE, check.names=FALSE)
  if (ncol(df) < 10) stop("Phased VCF has < 10 columns; got ", ncol(df))
  df <- df[, 1:10, drop=FALSE]
  colnames(df) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
  df
}, error=function(e) stop("Failed to read phased VCF: ", conditionMessage(e)))

vcf$CHROM <- to_chr(vcf$CHROM)
vcf$POS   <- suppressWarnings(as.integer(vcf$POS))
vcf <- vcf[!is.na(vcf$POS), , drop=FALSE]

# Extract GT by matching FORMAT field
fmt_list <- strsplit(vcf$FORMAT, ":", fixed=TRUE)
smp_list <- strsplit(vcf$SAMPLE, ":", fixed=TRUE)
gt_idx   <- sapply(fmt_list, function(x) match("GT", x))
vcf$GT   <- mapply(function(sf, ix) if (!is.na(ix) && ix <= length(sf)) sf[[ix]] else NA_character_,
                   smp_list, gt_idx, USE.NAMES=FALSE)

# Keep phased bi-allelic SNPs
vcf <- vcf[!is.na(vcf$GT) & grepl("\\|", vcf$GT) &
             nchar(vcf$REF)==1 & nchar(vcf$ALT)==1, , drop=FALSE]
vcf$variant <- paste0(vcf$CHROM, ":", vcf$POS)

cat(sprintf("[info] phased SNPs in VCF (with '|'): %d\n", nrow(vcf)))

## ---------------------------------------------------------------------------
## JOIN & ORIENT BY HAPLOTYPE
## ---------------------------------------------------------------------------
rna$GT <- vcf$GT[match(rna$variant, vcf$variant)]
overlap_n <- sum(!is.na(rna$GT))
cat(sprintf("[info] loci overlapping RNA & phased VCF: %d\n", overlap_n))

het <- rna$GT %in% c("0|1","1|0")
dat <- rna[het &
             nchar(rna$REF)==1 & nchar(rna$ALT)==1 &
             (rna$REF.COUNTS + rna$ALT.COUNTS) > 0, , drop=FALSE]

if (!"gene" %in% colnames(dat)) dat$gene <- NA_character_
dat$gene[is.na(dat$gene) | dat$gene==""] <- "NA_gene"

cat(sprintf("[info] heterozygous phased SNPs retained: %d\n", nrow(dat)))

if (nrow(dat) == 0L) {
  cat("[warn] No heterozygous phased SNPs after overlap/filter. Writing empty results.\n")
  saveRDS(list(), file=file.path(outdir, "ASEresults_1s_haplotypesKnown.rds"))
  saveRDS(list(geneOutput=data.frame(), locusOutput=list()), file=file.path(outdir, "MBASEDresults.rds"))
  quit(save="no", status=0)
}

dat$alleleA        <- ifelse(dat$GT == "1|0", dat$ALT, dat$REF)  # hap1
dat$alleleB        <- ifelse(dat$GT == "1|0", dat$REF, dat$ALT)  # hap2
dat$alleleA.counts <- ifelse(dat$GT == "1|0", dat$ALT.COUNTS, dat$REF.COUNTS)
dat$alleleB.counts <- ifelse(dat$GT == "1|0", dat$REF.COUNTS, dat$ALT.COUNTS)

dat <- dat %>%
  arrange(CHROM, POS) %>%
  group_by(gene) %>%
  mutate(label = paste0("SNV", dplyr::row_number())) %>%
  ungroup()
dat$SNV.ID <- paste0(dat$gene, ":", dat$label)

# Deduplicate SNV IDs just in case
dup <- duplicated(dat$SNV.ID)
if (any(dup)) {
  dat$SNV.ID[dup] <- paste0(dat$SNV.ID[dup], "_", seq_len(sum(dup)))
}

cat(sprintf("[info] genes covered (>=1 het SNP): %d\n", dplyr::n_distinct(dat$gene)))
cat(sprintf("[info] loci after QC: %d\n", nrow(dat)))

## ---------------------------------------------------------------------------
## BUILD SummarizedExperiment & RUN MBASED
## ---------------------------------------------------------------------------
mySNVs <- GRanges(
  seqnames = dat$CHROM,
  ranges   = IRanges(start = dat$POS, width = 1),
  aseID    = dat$gene,
  allele1  = dat$alleleA,
  allele2  = dat$alleleB
)
names(mySNVs) <- dat$SNV.ID

l1 <- matrix(as.integer(dat$alleleA.counts), ncol=1, dimnames=list(names(mySNVs),'mySample'))
l2 <- matrix(as.integer(dat$alleleB.counts), ncol=1, dimnames=list(names(mySNVs),'mySample'))

mySample <- SummarizedExperiment(
  assays    = list(lociAllele1Counts = l1, lociAllele2Counts = l2),
  rowRanges = mySNVs
)

cat(sprintf("[info] Running MBASED (isPhased=TRUE, %s)...\n",
            if (threads > 1L) paste0("parallel, workers=", threads) else "serial"))
ASE <- tryCatch(
  runMBASED(
    ASESummarizedExperiment = mySample,
    isPhased = TRUE,
    numSim   = numSim,
    BPPARAM  = bp
  ),
  error = function(e) {
    message("[error] runMBASED failed: ", conditionMessage(e))
    quit(save="no", status=1)
  }
)

saveRDS(ASE, file=file.path(outdir, "ASEresults_1s_haplotypesKnown.rds"))

res <- summarizeASEResults_1s(ASE)
res$geneOutput$padj <- p.adjust(res$geneOutput$pValueASE, method="BH")
res$geneOutput$significance <- ifelse(res$geneOutput$padj < 0.05, "padj<0.05", "padj>=0.05")
if (!"gene_biotype" %in% colnames(dat)) dat$gene_biotype <- NA_character_
res$geneOutput$geneBiotype <- dat$gene_biotype[match(res$geneOutput$gene, dat$gene)]

saveRDS(res, file=file.path(outdir, "MBASEDresults.rds"))
cat("[info] Finished MBASED (phased, serial).\n")
