#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

list_of_packages <- c("rpart")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(list_of_packages, require, character.only = TRUE)

library(rattle)
library(rpart.plot)
library(RColorBrewer)
library(partykit)
library(kohonen)

db_file <- "data/db1_decision_tree.db"
#db_file <- args[1]

#new_data <- data.frame(new_data)
#new_data <- args[2]

db <- read.table(db_file, header = TRUE)

db <- subset(db, Type != "Unknown")

db <- subset(db, ATG != 0 )

Formula = Type ~ AA + AC + AG + AT + CA + CC + CG + CT + GA + GC + GG + GT + TA + TC + TG + TT + AAA + AAC + AAG + AAT + ACA + ACC + ACG + ACT + AGA + AGC + AGG + AGT + ATA + ATC + ATG + ATT + CAA + CAC + CAG + CAT + CCA + CCC + CCG + CCT + CGA + CGC + CGG + CGT + CTA + CTC + CTG + CTT + GAA + GAC + GAG + GAT + GCA + GCC + GCG + GCT + GGA + GGC + GGG + GGT + GTA + GTC + GTG + GTT + TAA + TAC + TAG + TAT + TCA + TCC + TCG + TCT + TGA + TGC + TGG + TGT + TTA + TTC + TTG + TTT
fit <- rpart(formula = Formula, method = "class", data = db, control = rpart.control("minsplit" = 2, cp = 0)) #cp se refere ao limite de corte da importancia de uma variavel, minsplit se refere ao minimo de elementos para criar uma divisao
fit <- rpart(formula = Formula, method = "class", data = db, control = rpart.control("minsplit" = 2))
plotcp(fit) # visualize cross-validation results 
printcp(fit) # Plotcp e printcp mostram um bom ponto de corte para o valor de cp (limiar de corte do valor de cp)
summary(fit) # detailed summary of splits

# plot tree 
plot(fit, uniform=TRUE, main="Classification Tree for Contigs")
text(fit, use.n=TRUE, all=TRUE, cex=.8)
fancyRpartPlot(fit)

fitp <- prune(fit, cp = 0.2)
plot(fitp, uniform=TRUE, main="Classification Tree for Contigs")
text(fit, use.n=TRUE, all=TRUE, cex=.8)


plot(as.party(fit), type="simple", gp = gpar(fontsize = 10),     # font size changed to 6
     inner_panel=node_inner,
     ip_args=list(
         abbreviate = FALSE, 
         id = FALSE))

new.fit <- prp(fit,snip=TRUE)$obj
fancyRpartPlot(new.fit)
######################################
#PREDICTING
Prediction <- predict(fit, db, type = "class" )
result <- data.frame(Contig = db$lib_name, Type = Prediction)
Virus <- subset(result, Type == "Virus")
head(result)
Virus
write.csv(result, file = "data/prediction_tree.csv", row.names = FALSE)

