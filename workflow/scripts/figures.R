
## ---------------------------------------------------------------------------
## Allelic Imbalance in Expression using MBASED, part 2
## Vanessa Porter, Mar. 2022
## ---------------------------------------------------------------------------

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tibble))
suppressMessages(library(chromPlot))
suppressMessages(library(networkD3))
suppressMessages(library(htmlwidgets))
suppressMessages(library(reshape2))
suppressMessages(library(ggrepel))
suppressMessages(library(ggsci))

## ---------------------------------------------------------------------------
## FIGURES
## ---------------------------------------------------------------------------


# Make help options
option_list = list(
  make_option(c("-b", "--mbased"), type="character", default=NULL,
              help="mbased dataframe output", metavar="character"),
  make_option(c("-r", "--rpkm"), type="character", default=NULL,
              help="RPKM matrix", metavar="character"),
  make_option(c("-g", "--gene"), type="character", default = NULL,
              help="Ensembl gene annotation", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default = NULL,
              help="Sample name", metavar="character"),
  make_option(c("-m", "--min"), type="numeric", default = 1,
              help="Minimum RPKM value", metavar="numeric"),
  make_option(c("-t", "--maf_threshold"), type="numeric", default = 0.75,
              help="Minimum RPKM value", metavar="numeric"),
  make_option(c("-o", "--outdir"), type="character", default = NULL,
              help="Output directory", metavar="character")
)

### 
### SET UP THE DATAFRAME
###

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
df <- read.delim(opt$mbased, header = T, stringsAsFactors = F)
rpkm <- read.delim(opt$rpkm, header = T, stringsAsFactors = F)
all_genes <- read.delim(opt$gene, header = F, stringsAsFactors = F)
out <- opt$outdir
sample <- opt$sample
min <- opt$min
maf_threshold <- opt$maf_threshold

mafG_padjG <- paste0("MAF > ", maf_threshold, " & padj > 0.05")
mafL_padjG <- paste0("MAF < ", maf_threshold, " & padj > 0.05")
mafG_padjL <- paste0("MAF > ", maf_threshold, " & padj < 0.05")
mafL_padjL <- paste0("MAF < ", maf_threshold, " & padj < 0.05")


# fix sample name
sample <- ifelse(length(grep("-", sample)) == 0, sample, gsub("-", ".", sample))

# select sample
rpkm_sample <- rpkm[,c("gene", sample)] 
colnames(rpkm_sample) <- c("gene", "expr")

# make a colour filter
df$colour_filt <- ifelse(df$padj < 0.05 & df$majorAlleleFrequency > maf_threshold, mafG_padjL,
                         ifelse(df$padj > 0.05 & df$majorAlleleFrequency > maf_threshold, mafG_padjG,
                                ifelse(df$padj > 0.05 & df$majorAlleleFrequency < maf_threshold, mafL_padjG, mafL_padjL)))

# add the chromosome
df$chr <- all_genes$V1[match(df$gene, all_genes$V4)]

# set the factor levels
df$colour_filt <- factor(df$colour_filt, levels = c(mafG_padjL, mafL_padjL, 
                                                    mafG_padjG, mafL_padjG))
df$chr <- factor(df$chr, levels = c(paste0("chr", 1:22), "chrX"))
#df <- df[!is.na(df$chr),]

print("Beginning figures ...")

#### 
#### DOTPLOT
####

dotplot <- ggplot(df, aes(x = majorAlleleFrequency, y = padj, colour = colour_filt)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("#e74645", "black", "black", "grey")) +
  theme_bw() + 
  geom_hline(yintercept = 0.05, linetype = 2) +
  geom_vline(xintercept = maf_threshold, linetype = 2) +
  geom_text(aes(label = paste0(table(colour_filt)[mafG_padjL], " ASE genes"), x = 0.9, y = 0.75), size = 4.5, colour = "#e74645") +
  labs(x = "major allele frequency", y = "adjusted pvalue", colour = NULL) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 12, face = "bold", colour = "black"),
        axis.text = element_text(size = 10, colour = "black"))

ggsave(filename = paste0(out, "/aseGenesDot.pdf"), plot = dotplot, width = 5, height = 4, units = "in")

####
#### BARPLOT
main_chr <- c(paste0("chr", 1:22), "chrX")
normalize_chr <- function(chr_vals) {
  chr_vals <- as.character(chr_vals)
  chr_vals[chr_vals == "" | is.na(chr_vals)] <- NA_character_
  has_prefix <- grepl("^chr", chr_vals, ignore.case = TRUE)
  chr_vals[!is.na(chr_vals) & !has_prefix] <- paste0("chr", chr_vals[!is.na(chr_vals) & !has_prefix])
  chr_vals[has_prefix] <- paste0("chr", sub("^chr", "", chr_vals[has_prefix], ignore.case = TRUE))
  chr_vals[!chr_vals %in% main_chr] <- NA_character_
  factor(chr_vals, levels = main_chr)
}

cat(sprintf("[figures] df rows=%d, unique genes=%d\n", nrow(df), length(unique(df$gene))))
cat(sprintf("[figures] all_genes rows=%d, unique gene symbols (V4)=%d\n",
            nrow(all_genes), length(unique(all_genes$V4))))

# Map df$gene -> chromosome via all_genes (V4 = gene symbol, V1 = chrom)
df_chr <- df
midx <- match(df_chr$gene, all_genes$V4)
raw_chr <- all_genes$V1[midx]

cat(sprintf("[figures] gene->anno mapping: matched=%d (%.1f%%), unmapped=%d (%.1f%%)\n",
            sum(!is.na(midx)),
            100 * mean(!is.na(midx)),
            sum(is.na(midx)),
            100 * mean(is.na(midx))))

# Normalize chromosomes and drop non-main
df_chr$chr_raw  <- raw_chr
df_chr$chr      <- normalize_chr(raw_chr)
pre_n <- nrow(df_chr)
df_chr <- df_chr[!is.na(df_chr$chr), , drop = FALSE]
cat(sprintf("[figures] kept %d/%d rows after chr filter to %s\n",
            nrow(df_chr), pre_n, paste(main_chr, collapse=",")))

# Quick value counts
cat("[figures] colour_filt counts (post-filter):\n")
print(sort(table(df_chr$colour_filt), decreasing = TRUE))
cat("[figures] chr counts (post-filter):\n")
print(table(df_chr$chr))

# Write a tiny debug tsv
dbg_path <- file.path(out, "debug_barplot_input.tsv")
utils::write.table(
  df_chr[, c("gene", "colour_filt", "chr_raw", "chr")],
  file = dbg_path, sep = "\t", quote = FALSE, row.names = FALSE
)
cat(sprintf("[figures] wrote %s (n=%d)\n", dbg_path, nrow(df_chr)))

# Plot or placeholder
out_pdf <- file.path(out, "aseGenesBar.pdf")
if (nrow(df_chr) == 0L) {
  cat("[figures] no rows to plot for barplot; writing placeholder PDF.\n")
  placeholder <- ggplot() + theme_void() +
    annotate("text", x = 0, y = 0, hjust = 0,
      label = "No genes available for barplot after filtering.\nSee debug_barplot_input.tsv for details.")
  ggsave(filename = out_pdf, plot = placeholder, width = 12, height = 5, units = "in")
} else {
  barplot <- ggplot(df_chr, aes(x = chr, fill = colour_filt)) +
    geom_bar() +
    scale_fill_manual(values = rev(c("#e0f0ea", "#574f7d", "#95adbe", "#e74645"))) +
    theme_bw() +
    labs(x = "chromosome", y = "number of genes", fill = "ASE results") +
    theme(axis.title = element_text(size = 12, face = "bold", colour = "black"),
          axis.text.y = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(size = 10, colour = "black"),
          legend.text = element_text(size = 10, colour = "black"),
          legend.title = element_text(size = 12, face = "bold", colour = "black"))
  ggsave(filename = out_pdf, plot = barplot, width = 12, height = 5, units = "in")
}


####
#### SANKEY PLOT
####

# set filters on the RPKM matrix
rpkm_sample$gene_biotype <- all_genes$V7[match(rpkm_sample$gene, all_genes$V4)]
rpkm_sample_filt1 <- rpkm_sample[rpkm_sample$gene_biotype %in% c("lincRNA", "miRNA", "protein_coding"),]
rpkm_sample_filt2 <- rpkm_sample_filt1[rpkm_sample_filt1$expr > min,] 

# get the input values for the plot
a <- nrow(rpkm_sample_filt1)
b <- c(nrow(rpkm_sample_filt2),nrow(rpkm_sample_filt1[rpkm_sample_filt1$expr <= min,] ))
c <- c(nrow(df),nrow(rpkm_sample_filt2[!rpkm_sample_filt2$gene %in% df$gene,]))
d <- c(nrow(df[df$padj < 0.05 & df$majorAlleleFrequency > maf_threshold,]),
       sum(nrow(df[df$padj >= 0.05 & df$majorAlleleFrequency <= maf_threshold,]),
           nrow(df[df$padj >= 0.05 & df$majorAlleleFrequency > maf_threshold,]),
           nrow(df[df$padj < 0.05 & df$majorAlleleFrequency <= maf_threshold,])))

# create a connection data frame
links <- data.frame(
  source=c(rep(paste0("All Genes (n=", a, ")"), 2), rep(paste0("Expressed (n=", b[1], ")"), 2), rep(paste0("Phased Genes (n=", c[1], ")"), 2)),
  target=c(paste0("Expressed (n=", b[1], ")"), paste0("Not Expressed (n=", b[2], ")"), 
           paste0("Phased Genes (n=", c[1], ")"), paste0("Unphased Genes (n=", c[2], ")"), 
           paste0("ASE Genes (n=", d[1], ")"), paste0("Biallelic Genes (n=", d[2], ")")), 
  value=c(b, c, d)
)

# create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# Reformat the links
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
sankey <- sankeyNetwork(Links = links, Nodes = nodes,
                        Source = "IDsource", Target = "IDtarget",
                        Value = "value", NodeID = "name", 
                        sinksRight=FALSE, fontSize = 18)

saveWidget(sankey, file=paste0(out, "/sankeyPlot.html"), selfcontained = F)


print("Figures completed")
