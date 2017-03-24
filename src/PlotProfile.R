#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Lista de bibliotecas para ser carregadas, ja instala automaticamente caso não tenha instalado
list_of_packages <- c("ggplot2")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(list_of_packages, require, character.only = TRUE)

#db_file <- "data/Aedes.db"
db_file <- args[1]

#lib <- "Drosophila_viruses_SPO-81"
lib <- args[2]

db <- read.table(db_file, header = TRUE)

fig <- ggplot() + 
    geom_bar(data = db, position = "stack", stat = "identity", mapping = aes(x = db[[1]] , y = db[[length(db)]]), colour = "skyblue") + 
    theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + labs(x = "", y = "Score", title = lib)

pdf(paste("plots/",lib,".profile",".pdf",sep = ""),width=7,height=4)
print(fig)
dev.off()
