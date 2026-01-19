
## ---------------------------------------------------------------------------
## Allelic Imbalance in Expression using MBASED, part 1
## Vanessa Porter, Oct. 2021
## ---------------------------------------------------------------------------

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(prob))
suppressMessages(library(tidyr))
suppressMessages(library(MBASED))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(BiocParallel))
suppressMessages(library(stats))
suppressMessages(library(tibble))

## ---------------------------------------------------------------------------
## LOAD INPUT 
## ---------------------------------------------------------------------------

# Make help options
option_list = list(
  make_option(c("-p", "--phase"), type="character", default=NULL,
              help="Phased VCF file (from WhatsHap)", metavar="character"),
  make_option(c("-r", "--rna"), type="character", default=NULL,
              help="Tumour RNA vcf file (from Strelka2)", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default = "mBASED",
              help="Output directory name", metavar="character"),
  make_option(c("-t", "--threads"), type="integer", default = 1,
              help="Threads used for mbased", metavar="integer")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir
threads <- opt$threads
dir.create(out, recursive=TRUE, showWarnings=FALSE)
message("Output directory: ", out)
if (!is.null(threads) && threads > 1) {
  message("Using MulticoreParam with ", threads, " workers")
  bpparam <- MulticoreParam(workers = threads)
} else {
  message("Using SerialParam")
  bpparam <- SerialParam()
}

## ---------------------------------------------------------------------------
## USER FUNCTIONS
## ---------------------------------------------------------------------------

# extract info from a list
list_n_item <- function(list, n){
  sapply(list, `[`, n)
}

# define function to print out the summary of ASE results
summarizeASEResults_1s <- function(MBASEDOutput) {
  
  geneOutputDF <- data.frame(
    majorAlleleFrequency = assays(MBASEDOutput)$majorAlleleFrequency[,1],
    pValueASE = assays(MBASEDOutput)$pValueASE[,1],
    pValueHeterogeneity = assays(MBASEDOutput)$pValueHeterogeneity[,1])
  
  geneAllele <- as.data.frame(assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor) %>%
    rownames_to_column(var = "rowname") %>%
    dplyr::mutate(gene = unlist(lapply(strsplit(rowname, split = ":"),function(x){x = x[1]}))) %>%
    dplyr::group_by(gene) %>%
    summarise(allele1IsMajor = unique(mySample))
  
  geneOutputDF$allele1IsMajor <- geneAllele$allele1IsMajor[match(rownames(geneOutputDF), geneAllele$gene)]
  
  lociOutputGR <- rowRanges(metadata(MBASEDOutput)$locusSpecificResults)
  lociOutputGR$allele1IsMajor <- assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[,1]
  lociOutputGR$MAF <- assays(metadata(MBASEDOutput)$locusSpecificResults)$MAF[,1]
  lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels=unique(lociOutputGR$aseID)))
  
  return(
    list(
      geneOutput=geneOutputDF,
      locusOutput=lociOutputList
    )
  )
}

## ---------------------------------------------------------------------------
## READ IN THE RNA SNV CALLS
## ---------------------------------------------------------------------------

# read in the RNA calls
rna_filt <- read.delim(opt$rna, header = T, comment.char = "#", stringsAsFactors = F)
colnames(rna_filt) <- c("CHROM", "POS", "AD","REF","ALT","gene", "gene_biotype") 
rna_filt$variant <- paste0(rna_filt$CHROM, ":", rna_filt$POS)
message("RNA TSV rows read: ", nrow(rna_filt))


## ---------------------------------------------------------------------------
## EXTRACT REF/ALT READ COUNTS
## ---------------------------------------------------------------------------

# Extract and add the read counts
expr <- strsplit(rna_filt$AD, ",")
rna_filt$REF.COUNTS <- as.numeric(list_n_item(expr, 1))
rna_filt$ALT.COUNTS <- as.numeric(list_n_item(expr, 2))
rna_count_ok <- sum(!is.na(rna_filt$REF.COUNTS) & !is.na(rna_filt$ALT.COUNTS))
message("RNA SNP rows with REF/ALT counts: ", rna_count_ok)
message("Unique genes in RNA: ", length(unique(rna_filt$gene)))

## ---------------------------------------------------------------------------
## MBASED WITH OR WITHOUT PHASING
## ---------------------------------------------------------------------------

### WITH PHASING
if (!is.null(opt$phase)){
  
  ### 
  ### PHASING
  ###
  
  # WhatsHap phased VCF from ONT sequencing pipeline
  wh_con <- if (grepl("\\.gz$", opt$phase, ignore.case=TRUE)) {
    gzfile(opt$phase, "rt")
  } else {
    file(opt$phase, "rt")
  }
  wh_raw <- read.delim(wh_con, header = F, comment.char = "#", stringsAsFactors = F)
  close(wh_con)
  if (ncol(wh_raw) < 10) {
    stop("VCF has fewer than 10 columns; expected at least one sample column.")
  }
  if (ncol(wh_raw) > 10) {
    message("VCF has multiple sample columns; using the first sample only.")
  }
  wh <- wh_raw[, c(1,2,4,5,9,10)]
  colnames(wh) <- c("CHROM", "POS", "REF", "ALT", "FORMAT", "SAMPLE")
  wh$variant <- paste0(wh$CHROM, ":", wh$POS)
  message("VCF rows read: ", nrow(wh))
  
  # remove indels
  wh <- wh %>% dplyr::filter(nchar(REF) == 1 & nchar(ALT) == 1)
  message("VCF SNP rows after removing indels: ", nrow(wh))
  
  # add genotype and phase block from the SAMPLE column using FORMAT keys
  fmt_list <- strsplit(wh$FORMAT, ":", fixed=TRUE)
  sample_list <- strsplit(wh$SAMPLE, ":", fixed=TRUE)
  get_format_value <- function(field_name) {
    idx <- lapply(fmt_list, function(x) match(field_name, x))
    mapply(function(v, i) {
      if (is.na(i) || i > length(v)) NA_character_ else v[[i]]
    }, sample_list, idx, USE.NAMES=FALSE)
  }
  wh$GT <- get_format_value("GT")
  wh$phaseBlock <- get_format_value("PS")
  n_phased <- sum(grepl("|", wh$GT, fixed=TRUE), na.rm=TRUE)
  n_ps <- sum(!is.na(wh$phaseBlock) & wh$phaseBlock != "")
  message("VCF phased SNPs (GT contains '|'): ", n_phased)
  message("VCF SNPs with non-empty PS: ", n_ps)
  if (n_ps == 0) {
    stop("No PS (phase set) values found in VCF; phased mode requires PS in FORMAT.")
  }
  
  # keep only phased genotypes
  wh <- wh[!is.na(wh$GT) & grepl("|", wh$GT, fixed=TRUE),]
  message("VCF phased SNPs after filtering: ", nrow(wh))
  
  # Add the genotype from WhatsHap 
  rna_filt$GT <- wh$GT[match(rna_filt$variant, wh$variant)]
  rna_filt$phaseBlock <- wh$phaseBlock[match(rna_filt$variant, wh$variant)]
  
  matched_variants <- sum(rna_filt$variant %in% wh$variant)
  message("RNA variants matched to phased VCF variants: ", matched_variants)
  if (matched_variants == 0 || matched_variants / max(1, nrow(rna_filt)) < 0.01) {
    message("Low overlap between RNA and VCF variants; check chromosome naming (e.g. chr1 vs 1).")
    message("RNA example variants: ", paste(head(rna_filt$variant, 2), collapse = ", "))
    message("VCF example variants: ", paste(head(wh$variant, 2), collapse = ", "))
  }
  
  # Find unphased genes with one variant (test)
  singleUnphased <- rna_filt %>%
    mutate(phase = variant %in% wh$variant) %>%
    left_join(rna_filt %>% group_by(gene) %>% summarize(n=n())) %>%
    dplyr::filter(!phase & n == 1)
  message("Single unphased genes with one variant: ",
          length(unique(singleUnphased$gene)), " genes, ",
          nrow(singleUnphased), " variants; allele1IsMajor set to NA later.")
  
  # Add genotype to unphased gene with one variant (test)
  rna_filt$GT[which(rna_filt$variant %in% singleUnphased$variant)] <- "1|0"
  
  # annotate the phased variants as alleleA and alleleB
  rna_filt$alleleA <- NA_character_
  rna_filt$alleleB <- NA_character_
  rna_filt$alleleA.counts <- NA_real_
  rna_filt$alleleB.counts <- NA_real_
  idx_10 <- rna_filt$GT == "1|0"
  idx_01 <- rna_filt$GT == "0|1"
  rna_filt$alleleA[idx_10] <- rna_filt$ALT[idx_10]
  rna_filt$alleleB[idx_10] <- rna_filt$REF[idx_10]
  rna_filt$alleleA[idx_01] <- rna_filt$REF[idx_01]
  rna_filt$alleleB[idx_01] <- rna_filt$ALT[idx_01]
  rna_filt$alleleA.counts[idx_10] <- rna_filt$ALT.COUNTS[idx_10]
  rna_filt$alleleB.counts[idx_10] <- rna_filt$REF.COUNTS[idx_10]
  rna_filt$alleleA.counts[idx_01] <- rna_filt$REF.COUNTS[idx_01]
  rna_filt$alleleB.counts[idx_01] <- rna_filt$ALT.COUNTS[idx_01]
  
  # phased only variants
  rna_phased <- rna_filt[complete.cases(rna_filt),]
  message("Phased loci retained after complete.cases filter: ", nrow(rna_phased))
  
  # Get phase block ID with most phased gene body variant 
  # (Used to account for phase blocks that occur in the middle of gene body)
  phaseBlockID <- rna_phased %>%
    group_by(gene, phaseBlock) %>%
    summarize(n=n()) %>%
    group_by(gene) %>%
    top_n(1, n) %>%
    pull(phaseBlock)
    
  rna_phased <- rna_phased %>%
    dplyr::filter(phaseBlock %in% phaseBlockID)
  message("Phased loci after dominant phaseBlock filter: ", nrow(rna_phased))
  message("Genes after dominant phaseBlock filter: ", length(unique(rna_phased$gene)))
  if (nrow(rna_phased) == 0) {
    stop("No phased loci available for MBASED after filtering.")
  }
  
  # make SNV IDs
  rna_phased <- rna_phased %>%
    arrange(CHROM, POS) %>%
    group_by(gene) %>%
    mutate(label = paste0("SNV",1:n()))
  rna_phased$SNV.ID <- paste0(rna_phased$gene, ":", rna_phased$label)
  
  ### 
  ### MBASED
  ###
  
  message("Beginning MBASED ...")
  
  # make the GRanges object of the loci
  mySNVs <- GRanges(seqnames=rna_phased$CHROM,
                     ranges=IRanges(start=rna_phased$POS, width=1),
                     aseID=rna_phased$gene,
                     allele1=rna_phased$alleleA,
                     allele2=rna_phased$alleleB)
  names(mySNVs) <- rna_phased$SNV.ID
  
  # create input RangedSummarizedExperiment object
  mySample <- SummarizedExperiment(
    assays=list(lociAllele1Counts=matrix(rna_phased$alleleA.counts,
                                         ncol=1,
                                         dimnames=list(names(mySNVs),'mySample')),
                lociAllele2Counts=matrix(rna_phased$alleleB.counts,
                                         ncol=1,
                                         dimnames=list(names(mySNVs),'mySample'))),
    rowRanges=mySNVs
  )
  
  # run MBASED
  ASEresults_1s_haplotypesKnown <- runMBASED(ASESummarizedExperiment=mySample,
                                             isPhased=TRUE,
                                             numSim=10^6,
                                             BPPARAM = bpparam)
  
  saveRDS(ASEresults_1s_haplotypesKnown, file=paste0(out, "/ASEresults_1s_haplotypesKnown.rds"))
  # extract results
  results <- summarizeASEResults_1s(ASEresults_1s_haplotypesKnown)
  
  # adjust the pvalue with BH correction
  results$geneOutput$padj <- p.adjust(p = results$geneOutput$pValueASE, method = "BH")
  results$geneOutput$significance <- as.factor(ifelse(results$geneOutput$padj < 0.05, "padj < 0.05", "padj > 0.05"))
  results$geneOutput$gene <- rownames(results$geneOutput)
  
  results$geneOutput$allele1IsMajor[results$geneOutput$gene %in% singleUnphased$gene] = NA
  
  # add the locus
  results$geneOutput$geneBiotype <- rna_filt$gene_biotype[match(results$geneOutput$gene, rna_filt$gene)]

### WITHOUT PHASING
} else {
  
  # make SNV labels
  rna_filt <- rna_filt %>%
    arrange(CHROM, POS) %>%
    group_by(gene) %>%
    mutate(label = paste0("SNV",1:n()))
  rna_filt$SNV.ID <- paste0(rna_filt$gene, ":", rna_filt$label)
  
  ### 
  ### MBASED
  ###
  
  message("Beginning MBASED ...")
  
  # make the GRanges object of the loci
  mySNVs <- GRanges(seqnames=rna_filt$CHROM,
                    ranges=IRanges(start=rna_filt$POS, width=1),
                    aseID=rna_filt$gene,
                    allele1=rna_filt$REF,
                    allele2=rna_filt$ALT)
  names(mySNVs) <- rna_filt$SNV.ID
  
  ## create input RangedSummarizedExperiment object
  mySample <- SummarizedExperiment(
    assays=list(lociAllele1Counts=matrix(rna_filt$REF.COUNTS,
                                         ncol=1,
                                         dimnames=list(names(mySNVs),'mySample')),
                lociAllele2Counts=matrix(rna_filt$ALT.COUNTS,
                                         ncol=1,
                                         dimnames=list(names(mySNVs),'mySample'))),
    rowRanges=mySNVs
  )
  
  # run MBASED
  if (nrow(rna_filt) == 0) {
    stop("No loci available for MBASED in unphased mode.")
  }
  ASEresults_1s_haplotypesUnknown <- runMBASED(ASESummarizedExperiment=mySample,
                                               isPhased=FALSE,
                                               numSim=10^6,
                                               BPPARAM = bpparam)
  saveRDS(ASEresults_1s_haplotypesUnknown, file=paste0(out, "/ASEresults_1s_haplotypesUnknown.rds"))
  
  # extract results
  results <- summarizeASEResults_1s(ASEresults_1s_haplotypesUnknown)
  
  # adjust the pvalue with BH correction
  results$geneOutput$padj <- p.adjust(p = results$geneOutput$pValueASE, method = "BH")
  results$geneOutput$significance <- as.factor(ifelse(results$geneOutput$padj < 0.05, "padj < 0.05", "padj > 0.05"))
  results$geneOutput$gene <- rownames(results$geneOutput)
  
  # add the locus
  results$geneOutput$geneBiotype <- rna_filt$gene_biotype[match(results$geneOutput$gene, rna_filt$gene)]
  
} 

# save the results 
saveRDS(results, file=paste0(out, "/MBASEDresults.rds"))
message("Finished MBASED")
