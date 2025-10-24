#!/usr/bin/env Rscript

## ---------------------------------------------------------------------------
## MBASED (phased only) — minimal driver for SHAPEIT/VCF + RNA TSV
## Inputs:
##   --phase : phased VCF (single-sample), header lines start with "##"
##             columns: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <SAMPLE>
##             FORMAT must contain GT with "|" (e.g., 0|1 or 1|0)
##   --rna   : TSV with columns:
##             CHROM POS GEN[0].AD REF ALT ANN[0].GENE ANN[0].BIOTYPE
## Outputs:
##   <outdir>/ASEresults_1s_haplotypesKnown.rds
##   <outdir>/MBASEDresults.rds
## ---------------------------------------------------------------------------

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(MBASED))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(GenomicRanges))
suppressMessages(library(BiocParallel))

option_list <- list(
  make_option(c("-p","--phase"),   type="character", help="Phased VCF file"),
  make_option(c("-r","--rna"),     type="character", help="RNA SNV TSV"),
  make_option(c("-o","--outdir"),  type="character", default="mBASED", help="Output dir"),
  make_option(c("-t","--threads"), type="integer",   default=1, help="Threads")
)
opt <- parse_args(OptionParser(option_list=option_list))
stopifnot(!is.null(opt$phase), !is.null(opt$rna))
outdir <- opt$outdir; if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
threads <- if (is.null(opt$threads) || is.na(opt$threads) || opt$threads < 1) 1L else as.integer(opt$threads)

cat(sprintf("[info] threads=%d\n", threads))

## helpers
list_n_item <- function(lst, n) sapply(lst, `[`, n)
to_chr <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))

summarizeASEResults_1s <- function(MBASEDOutput) {
  geneOutputDF <- data.frame(
    majorAlleleFrequency = assays(MBASEDOutput)$majorAlleleFrequency[,1],
    pValueASE            = assays(MBASEDOutput)$pValueASE[,1],
    pValueHeterogeneity  = assays(MBASEDOutput)$pValueHeterogeneity[,1]
  )
  geneAllele <- as.data.frame(assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor) %>%
    tibble::rownames_to_column(var="rowname") %>%
    mutate(gene = vapply(strsplit(rowname, ":"), `[`, "", 1L)) %>%
    group_by(gene) %>% summarise(allele1IsMajor = unique(mySample), .groups="drop")
  geneOutputDF$gene <- rownames(geneOutputDF)
  geneOutputDF$allele1IsMajor <- geneAllele$allele1IsMajor[match(geneOutputDF$gene, geneAllele$gene)]
  lociOutputGR <- rowRanges(metadata(MBASEDOutput)$locusSpecificResults)
  lociOutputGR$allele1IsMajor <- assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[,1]
  lociOutputGR$MAF            <- assays(metadata(MBASEDOutput)$locusSpecificResults)$MAF[,1]
  list(geneOutput = geneOutputDF, locusOutput = split(lociOutputGR, factor(lociOutputGR$aseID)))
}

## ---------------------------------------------------------------------------
## READ RNA TSV
## ---------------------------------------------------------------------------
rna <- read.delim(opt$rna, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(all(c("CHROM","POS","GEN[0].AD","REF","ALT","ANN[0].GENE","ANN[0].BIOTYPE") %in% colnames(rna)))
colnames(rna) <- c("CHROM","POS","AD","REF","ALT","gene","gene_biotype")

rna$CHROM   <- to_chr(as.character(rna$CHROM))
rna$POS     <- as.integer(rna$POS)
rna$variant <- paste0(rna$CHROM, ":", rna$POS)

ad <- strsplit(rna$AD, ",", fixed=TRUE)
rna$REF.COUNTS <- suppressWarnings(as.numeric(list_n_item(ad, 1)))
rna$ALT.COUNTS <- suppressWarnings(as.numeric(list_n_item(ad, 2)))

cat(sprintf("[info] RNA rows: %d\n", nrow(rna)))

## ---------------------------------------------------------------------------
## READ PHASED VCF (gz-aware)
## ---------------------------------------------------------------------------
con <- gzfile(opt$phase, open = "rt")
vcf <- read.delim(con, header = FALSE, comment.char = "#",
                  stringsAsFactors = FALSE, check.names = FALSE)
close(con)

stopifnot(ncol(vcf) >= 10)
colnames(vcf)[1:10] <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
vcf$CHROM   <- to_chr(as.character(vcf$CHROM))
vcf$POS     <- as.integer(vcf$POS)
vcf$variant <- paste0(vcf$CHROM, ":", vcf$POS)


## parse GT by name from FORMAT/SAMPLE
fmt_list <- strsplit(vcf$FORMAT, ":", fixed=TRUE)
smp_list <- strsplit(vcf$SAMPLE, ":", fixed=TRUE)
gt_idx   <- sapply(fmt_list, function(x) match("GT", x))
vcf$GT   <- mapply(function(sf, ix) if (!is.na(ix) && ix <= length(sf)) sf[[ix]] else NA_character_, smp_list, gt_idx,
                   USE.NAMES=FALSE)

## keep phased genotypes with a pipe and SNPs only
vcf <- vcf[!is.na(vcf$GT) & grepl("\\|", vcf$GT) &
             nchar(vcf$REF)==1 & nchar(vcf$ALT)==1, , drop=FALSE]

cat(sprintf("[info] phased SNPs in VCF (with '|'): %d\n", nrow(vcf)))

## ---------------------------------------------------------------------------
## JOIN & ORIENT BY HAPLOTYPE
## ---------------------------------------------------------------------------
rna$GT <- vcf$GT[match(rna$variant, vcf$variant)]
overlap_n <- sum(!is.na(rna$GT))
cat(sprintf("[info] loci overlapping RNA & phased VCF: %d\n", overlap_n))

## heterozygous only (informative for ASE)
het_pat <- rna$GT %in% c("0|1","1|0")
dat <- rna[het_pat &
             !is.na(rna$REF.COUNTS) & !is.na(rna$ALT.COUNTS) &
             nchar(rna$REF)==1 & nchar(rna$ALT)==1, , drop=FALSE]

cat(sprintf("[info] heterozygous phased SNPs retained: %d\n", nrow(dat)))

if (nrow(dat) == 0L) {
  cat("[warn] No heterozygous phased SNPs after overlap/filter. Writing empty results.\n")
  saveRDS(list(), file=file.path(outdir, "ASEresults_1s_haplotypesKnown.rds"))
  saveRDS(list(geneOutput=data.frame(), locusOutput=list()), file=file.path(outdir, "MBASEDresults.rds"))
  quit(save="no", status=0)
}

## alleleA = haplotype 1; alleleB = haplotype 2 (based on GT)
dat$alleleA        <- ifelse(dat$GT == "1|0", dat$ALT, dat$REF)  # hap1
dat$alleleB        <- ifelse(dat$GT == "1|0", dat$REF, dat$ALT)  # hap2
dat$alleleA.counts <- ifelse(dat$GT == "1|0", dat$ALT.COUNTS, dat$REF.COUNTS)
dat$alleleB.counts <- ifelse(dat$GT == "1|0", dat$REF.COUNTS, dat$ALT.COUNTS)

## label SNVs per gene
dat <- dat %>%
  arrange(CHROM, POS) %>%
  group_by(gene) %>%
  mutate(label = paste0("SNV", dplyr::row_number())) %>%
  ungroup()
dat$SNV.ID <- paste0(dat$gene, ":", dat$label)

cat(sprintf("[info] genes covered (>=1 het SNP): %d\n", dplyr::n_distinct(dat$gene)))

## ---------------------------------------------------------------------------
## BUILD SE & RUN MBASED (isPhased=TRUE)
## allele1/allele2 are set to haplotype-consistent alleleA/alleleB
## ---------------------------------------------------------------------------
mySNVs <- GRanges(
  seqnames = dat$CHROM,
  ranges   = IRanges(start = dat$POS, width = 1),
  aseID    = dat$gene,
  allele1  = dat$alleleA,
  allele2  = dat$alleleB
)
names(mySNVs) <- dat$SNV.ID

mySample <- SummarizedExperiment(
  assays = list(
    lociAllele1Counts = matrix(dat$alleleA.counts, ncol=1,
                               dimnames=list(names(mySNVs),'mySample')),
    lociAllele2Counts = matrix(dat$alleleB.counts, ncol=1,
                               dimnames=list(names(mySNVs),'mySample'))
  ),
  rowRanges = mySNVs
)

cat("[info] Running MBASED (isPhased=TRUE)...\n")
ASE <- runMBASED(
  ASESummarizedExperiment = mySample,
  isPhased = TRUE,
  numSim   = 10^6,
  BPPARAM  = MulticoreParam(workers = threads)
)

saveRDS(ASE, file=file.path(outdir, "ASEresults_1s_haplotypesKnown.rds"))

res <- summarizeASEResults_1s(ASE)
res$geneOutput$padj <- p.adjust(res$geneOutput$pValueASE, method="BH")
res$geneOutput$significance <- ifelse(res$geneOutput$padj < 0.05, "padj<0.05", "padj>=0.05")
# add biotype if available
res$geneOutput$geneBiotype <- dat$gene_biotype[match(res$geneOutput$gene, dat$gene)]

saveRDS(res, file=file.path(outdir, "MBASEDresults.rds"))
cat("[info] Finished MBASED (phased).\n")
