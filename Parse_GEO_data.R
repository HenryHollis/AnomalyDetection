library(tidyverse)

library(GEOquery)
library(data.table)
setwd("~/Documents/R/Anomoly_detection/")
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*100)
myGSE = "GSE48113" # qin, demo 

gset <- getGEO(myGSE, GSEMatrix =TRUE)
if (length(gset) > 1) idx <- grep("GPL90", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset) #This is the expression matrix
dict = gset@featureData@data[, c('ID', 'GENE_SYMBOL')]
dict = drop_na(dict)
genes = as.data.frame(ex)

pat_dict = gset@phenoData@data[, c("geo_accession", "time point:ch1", "subject:ch1", "sleep timing:ch1") ]
#pat_dict = pat_dict[pat_dict$`sleep timing:ch1` == "In phase with respect to melatonin" & pat_dict$`time point:ch1` !=7, ]
pat_dict = pat_dict[pat_dict$`sleep timing:ch1` == "Out of phase with respect to melatonin" & pat_dict$`time point:ch1` !=7, ]

transpose = t(genes)
transpose = as.data.frame(transpose)
transpose$timepoints= pat_dict$`time point:ch1`[match(unlist(rownames(transpose)), pat_dict$geo_accession)]

transpose = transpose[!is.na(transpose$timepoint), ]
genes = as.data.frame(t(transpose))

# df %>% group_by(a, b) %>% summarise(n = n())
genes = select(genes, !c("GSM1168684", "GSM1168685", "GSM1168686", "GSM1168687", "GSM1168688"))

genes$symbols = dict$GENE_SYMBOL[match(unlist(as.numeric(rownames(genes))), dict$ID)]

genes = genes[complete.cases(genes), ]
genes = select(genes, symbols, everything())
GOI = c("PER1","PDK1","NPEPL1","EPHX2","GPCPD1","MS4A3","GNG2","MUM1","IL13RA1",
"IL1B","STIP1","ID3","DHRS13","CHSY1","MEGF6","NR1D1","AK5","TCN1",
"ZNF438","CYB561","NSUN3","NR1D2","SLPI","POLH","CD38","PARP2","SYT11",
"TIAM2","PGPEP1","SH2D1B","CD1C","C12orf75","REM2","LLGL2","FKBP4")

#genes_reduced = filter(genes, symbols %in% GOI)



#get every nth column:

samples_1_4 = genes[,seq(2, ncol(genes), 3)]
samples_2_5 = genes[,seq(3, ncol(genes), 3)]
samples_3_6 = genes[,seq(4, ncol(genes), 3)]

stack_cols = function(data){
  col_odd <- seq_len(ncol(data)) %% 2 
  bottom <- data[ , col_odd == 1]
  top <- data[ , col_odd == 0]
  names(bottom)<-names(top)
  return(rbind(bottom, top))
}
 temp1= stack_cols(samples_1_4)
 temp2 =stack_cols(samples_2_5)
 temp3 = stack_cols(samples_3_6)
 
 partial_cat = bind_cols(temp1, temp2)
 out_data = bind_cols(partial_cat, temp3)
 out_data = t(out_data)
 write.csv(out_data, "GSE48113_csv_data_Anomalous.csv", row.names = F )
