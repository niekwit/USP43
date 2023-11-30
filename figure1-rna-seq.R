library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(readr)

# set working directory (where git repo has been cloned)
setwd("/mnt/4TB_SSD/analyses/Projects/USP43")

# plotting parameter
lfc <- 0.5
pval <- 2
font.size <- 18

#### RNA-Seq volcano plots ####

# get DUB genes names from CRISPR sgRNA library fasta file
fasta <- "/mnt/4TB_SSD/analyses/Projects/USP43/data/CRISPR_screens/exp1/resources/DUBonly.fasta"
dubs <- system(paste0("grep -o -P '(?<=>)[^_]*' ", fasta, " | sort | uniq | sed '/Control/d'"), intern=TRUE)

# load data and filter for DUB only genes
data.wt <- read_csv("data/RNA-Seq/results/deseq2/Control_hypoxia_vs_Control_normoxia.csv") %>%
  filter(external_gene_name %in% dubs) %>%
  mutate(sample = "WT 1% O2 vs. 21% O2")

data.hif1b <- read_csv("data/RNA-Seq/results/deseq2/HIF1B_hypoxia_vs_Control_hypoxia.csv") %>%
  filter(external_gene_name %in% dubs) %>%
  mutate(sample = "HIF-1\U03B2 vs. WT 1% O2")

# create data frame for plotting
df <- rbind(data.wt, data.hif1b) %>%
  mutate(sample = factor(sample, levels = c("WT 1% O2 vs. 21% O2", "HIF-1\U03B2 vs. WT 1% O2")))

# create df for text labels in plot
df.label <- df %>%
  filter(log2FoldChange >= lfc | log2FoldChange <= -lfc,
         -log10(padj) > pval)

#df <- df %>%
#  mutate(colour = )

# create volcano plot
p <- ggplot(df, aes(x = log2FoldChange, 
                   y = -log10(padj),
                   colour = case_when(
                     log2FoldChange >= lfc & -log10(padj) > pval ~ "red",
                     log2FoldChange <= -lfc & -log10(padj) > pval ~ "navy",
                     TRUE ~ "black"))) + 
  theme_cowplot(font.size) +
  theme(plot.title = element_text(hjust = 0.50,),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(2, "lines"),
        strip.background = element_rect(colour="white", fill="white")) +
  geom_point(shape = 21,
             size = 7.5,
             stroke = 0.75) +
  facet_grid(. ~ sample) +
  xlab("log2(fold change)") +
  ylab("-log10(Adj. p value)") +
  xlim(-1.5, 1.5) +
  geom_hline(yintercept=pval, linetype="dashed", 
             color = "red", size = 0.5) +
  geom_vline(xintercept=c(lfc,-lfc), linetype="dashed", 
             color = "red", size = 0.5) +
  geom_label_repel(size=5,
                   aes(x = log2FoldChange,
                       y = -log10(padj),
                       label = external_gene_name), 
                   data = df.label,
                   nudge_x=-0.125,
                   nudge_y=0.05,
                   colour = "black") + 
  scale_colour_manual(values=c("black", "blue", "red"),
                      guide = "none")
  
# save plot (ggsave does not save unicode properly)
cairo_pdf("figures/figure1_rna-seq.pdf",
          width = 10,
          height = 5)
print(p)
dev.off()


# save output sessionInfo() to file
writeLines(capture.output(sessionInfo()), "sessionInfo/sessionInfo_figure1-rna-seq.txt")

