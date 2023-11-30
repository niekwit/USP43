# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(ggplot2)
library(DESeq2)
library(RColorBrewer)
library(ggrepel)
library(cowplot)

# read in data
load(snakemake@input[[1]])

# log transform data
rld <- rlog(dds)

# select appropriate colour palette
if (length(unique(rld$treatment)) == 2) {
  palette <- "Paired"
} else {
  palette <- "Dark2"
}

# change sample labels for publication
for (i in seq(rld$sample)) {
  normoxia <-  " 21% O2"
  rld$sample[[i]] <- gsub("-normoxia", normoxia, rld$sample[[i]])
  hypoxia <- " 1% O2"
  rld$sample[[i]] <- gsub("-hypoxia", hypoxia, rld$sample[[i]])
  rld$sample[[i]] <- gsub("HIF1B", "HIF-1\U03B2", rld$sample[[i]])
}


#create PCA plot
pca <- plotPCA(rld, intgroup=c("genotype", "treatment")) +
  geom_label_repel(aes(label = rld$sample)) + 
  guides(colour = "none") +
  theme_cowplot(16) +
  scale_color_brewer(palette = palette) 

# save plot to file
cairo_pdf(snakemake@output[[1]],
          width = 10,
          height = 10)
print(pca)
dev.off()

# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")


