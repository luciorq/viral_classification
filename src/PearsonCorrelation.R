#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
# Lista de bibliotecas para ser carregadas, ja instala automaticamente caso não tenha instalado
list_of_packages <- c("corrplot")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(list_of_packages, require, character.only = TRUE)
library(corrplot)

#db_file <- "data/virus_ORF.db"
db_file <- args[1]
db <- read.table(db_file, header = TRUE)
db <- data.matrix(db[1:80,2:length(db)])

corr <- cor(db, use="complete.obs", method="pearson")

file_name <- strsplit(db_file, "/")[[1]]
fig_file = paste("plots/",file_name[2],"_corr.pdf",sep = "")

pdf(fig_file)
fig <- corrplot(corr, order = "hclust")
dev.off()


corr_di <- cor(db[1:16,], use="complete.obs", method="pearson")
pdf(paste("plots/",file_name[2],"_di_corr.pdf",sep = ""))
fig2 <- corrplot(corr_di, order = "hclust")
dev.off()

corr_tri <- cor(db[17:80,], use="complete.obs", method="pearson")
pdf(paste("plots/",file_name[2],"_tri_corr.pdf",sep = ""))
fig3 <- corrplot(corr_tri, order = "hclust")
dev.off()