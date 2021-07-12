

library(GEOquery)
library(biomaRt)


#require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_hugene_1_0_st_v1",
    "gene_biotype",
    "external_gene_name"),
  filter = "affy_hugene_1_0_st_v1",
  values = rownames(merged_data), uniqueRows=TRUE,useCache = FALSE)
#Error in UseMethod("filter_") : 
#no applicable method for 'filter_' 
#applied to an object of class "c('tbl_SQLiteConnection', 
#'tbl_dbi', 'tbl_sql', 'tbl_lazy', 'tbl')" SOLVED with useCache argument added.

indicesLookup <- match(rownames(merged_data),annotLookup$affy_hugene_1_0_st_v1)
rownames(merged_data) <- paste(annotLookup[indicesLookup, "external_gene_name"], c(1:length(indicesLookup)), sep="_")

head(rownames(merged_data),20)
extGeneNames<-gsub("_[0-9]*$", "", rownames(merged_data))
gene.exp<-merged_data
dim(gene.exp)

rownames(gene.exp)<-extGeneNames
dim(gene.exp)

complete_record<-gene.exp[rownames(gene.exp) != "NA" & rownames(gene.exp) != "" , ]
dim(complete_record) 

write.csv(complete_record,"mergedData_externalGeneNames_mat.csv")

saveRDS(complete_record,"merged_data2.RDS")



