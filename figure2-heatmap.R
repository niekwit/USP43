library(dplyr)
library(readr)
library(forcats)
library(scales)
library(ggplot2)

font.size <- 16

# set working directory (where git repo has been cloned)
setwd("/mnt/4TB_SSD/analyses/Projects/USP43")

# genes for heatmap
genes <- c("HMOX1","PDGFB","ENO2","SERPINE1","PLIN2","LOX","HEY1","VEGFA",
           "EGLN3","P4HA2","CA9","ANKRD37","ALDOC","INSIG2","SLC2A1","SLC2A3", "NDRG1")
genes <- sort(genes,decreasing = FALSE)

# read DESeq2 data
df.wt <- read_csv("data/RNA-Seq/results/deseq2/Control_hypoxia_vs_Control_normoxia.csv") %>%
  filter(external_gene_name %in% genes) %>%
  dplyr::select(c("external_gene_name","log2FoldChange")) %>%
  mutate(sample = "WT")

df.hif <- read_csv("data/RNA-Seq/results/deseq2/HIF1B_hypoxia_vs_Control_hypoxia.csv") %>%
  filter(external_gene_name %in% genes) %>%
  dplyr::select(c("external_gene_name","log2FoldChange")) %>%
  mutate(sample = "HIF-1\U03B2 KO")

df.usp43 <- read_csv("data/RNA-Seq/results/deseq2/USP43_hypoxia_vs_Control_hypoxia.csv") %>%
  filter(external_gene_name %in% genes) %>%
  dplyr::select(c("external_gene_name","log2FoldChange")) %>%
  mutate(sample = "USP43 KO")

# create df for plotting
df <- rbind(df.wt, df.hif, df.usp43) %>%
  mutate(sample = factor(sample, levels = c("WT", "HIF-1\U03B2 KO", "USP43 KO"))) %>%
  dplyr::rename("log2(FC)" = log2FoldChange)

# create heatmap
p <- ggplot(df, aes(x = sample, 
                    y = forcats::fct_rev((external_gene_name)))) + 
  geom_tile(color = "black",
            aes(fill = `log2(FC)`)) + 
  scale_fill_gradientn(colours = c("navy", "white", "red"),
                       values = rescale(c(-4,0,7)), 
                       limits = c(-4,7),
                       breaks = seq(-4,7),
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black",
                                              barheight = 13.5)) +
  scale_x_discrete(expand = c(0, 0),
                   position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw(font.size) +
  theme(legend.title = element_text(size = font.size - 4),
        legend.text=element_text(size = font.size - 4),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, 
                                   hjust = -0.05,
                                   vjust = -0.25),
        axis.text.y = element_text(face = "italic")) +
  xlab(NULL) +
  ylab(NULL) +
  coord_equal()

# save plot (ggsave does not save unicode properly)
cairo_pdf("figures/figure2_heatmap.pdf",
          width = 3.5,
          height = 5)
print(p)
dev.off()

# save output sessionInfo() to file
writeLines(capture.output(sessionInfo()), "sessionInfo/sessionInfo_figure2-heatmap.txt")


