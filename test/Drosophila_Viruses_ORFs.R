#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly = TRUE)

file_name <- args[1]
#file_name <- "data/Drosophila_virus_ORFs.csv"

x <- read.csv(file_name, header = TRUE)


fig <- ggplot(data = x, aes(x = Viruses, y = Length, fill = Type)) +
  geom_boxplot() +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("plots/Drosophila_Viruses.ORFs.pdf",width=7,height=4)
print(fig)
dev.off()
