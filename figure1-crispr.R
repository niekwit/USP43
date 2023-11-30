library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(readxl)

# set working directory (where git repo has been cloned)
setwd("/mnt/4TB_SSD/analyses/Projects/USP43")

# plotting parameter
lfc <- 0.5
pval <- -log10(0.001)
font.size <- 18

# function to create dot plots from MAGeCK output data
dot.plot <- function(data){
  
  # load data
  data <- read.delim(data)
  
  # data for text labels in plot
  df.label <- data %>%
    filter(pos.lfc >= lfc,
           -log10(pos.p.value) > pval)
  
  # create dot plot
  p <- ggplot(data, 
              aes(x = pos.lfc, 
                  y = -log10(pos.p.value),
                  colour = pos.lfc >= lfc & -log10(pos.p.value) > pval)) + 
    geom_point(size = 7.5,
               stroke = 0.75,
               shape = 21) +
    scale_colour_manual(values=c("black", "red"),
                        guide = "none") +
    theme_cowplot(font.size) +
    labs(y = "-log10(p value)", 
         x = "log2(fold change)") +
    geom_vline(xintercept = c(-lfc,lfc), 
               linetype = "dashed", 
               colour = "red") +
    geom_hline(yintercept = pval, 
               linetype = "dashed", 
               colour = "red") +
    geom_text_repel(data=df.label,
                    aes(x=pos.lfc,
                        y=-log10(pos.p.value),
                        label=id),
                    nudge_y = 0.33,
                    colour = "black",
                    size = 5) +
    scale_y_continuous(limits=c(0,6)) 
  
  return(p)
  
}

#### Generate CRISPR screen figures ####
p1 <- dot.plot("data/CRISPR_screens/exp1/results/mageck/S8_vs_L8/S8_vs_L8.gene_summary.txt")
p2 <- dot.plot("data/CRISPR_screens/exp2/results/mageck/S7_vs_L7/S7_vs_L7.gene_summary.txt")

# combine plots
p <- plot_grid(p1, p2, labels = c('C', 'D'), label_size = font.size)
ggsave("figures/figure1_crispr.pdf", p, 
       width = 10, 
       height = 5)

# save output sessionInfo() to file
writeLines(capture.output(sessionInfo()), "sessionInfo/sessionInfo_figure1-crispr.txt")


