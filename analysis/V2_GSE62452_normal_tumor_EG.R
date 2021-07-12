#############################################################################
# Gene expression analysis of GSE155353
# 
# Organism	Homo sapiens
###
# 
# R code: Emine Guven
#############################################################################

library(gplots)
library(biomaRt)
library(matrixStats)
library(Biobase)
library(GEOquery)
library(rgr)

#load series and platform data from GEO
# load series and platform data from GEO

# gset <- getGEO("GSE62452", GSEMatrix =TRUE, getGPL=FALSE)
# if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
# gset_GSE62452 <- gset[[idx]]
# 
# saveRDS(gset_GSE62452,"gset_GSE62452.RDS")
gset_GSE62452<-readRDS("gset_GSE62452.RDS")
ex_GSE62452 <- exprs(gset_GSE62452)

ex_GSE62452_normal<-ex_GSE62452[,cbind("GSM1527195","GSM1527106", "GSM1527108","GSM1527110","GSM1527112", "GSM1527114",
                                    "GSM1527116", "GSM1527118", "GSM1527120", "GSM1527122", "GSM1527124",
                                    "GSM1527126", "GSM1527128", "GSM1527130", "GSM1527132", "GSM1527134",
                                    "GSM1527136","GSM1527138" ,"GSM1527140", "GSM1527142" ,"GSM1527144",
                                    "GSM1527146", "GSM1527148", "GSM1527150", "GSM1527152", "GSM1527154",
                                    "GSM1527156", "GSM1527158", "GSM1527160", "GSM1527162", "GSM1527164",
                                    "GSM1527166" ,"GSM1527168", "GSM1527170", "GSM1527172", "GSM1527174",
                                    "GSM1527176", "GSM1527178", "GSM1527180", "GSM1527182", "GSM1527184",
                                    "GSM1527186" ,"GSM1527188", "GSM1527190", "GSM1527192", "GSM1527194",
                                    "GSM1527197", "GSM1527201", "GSM1527203", "GSM1527206", "GSM1527208",
                                    "GSM1527211", "GSM1527214", "GSM1527217", "GSM1527221", "GSM1527222", 
                                    "GSM1527224","GSM1527226","GSM1527229", "GSM1527231", "GSM1527233")]
# log2 transform
qx <- as.numeric(quantile(ex_GSE62452_normal, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex_GSE62452_normal[which(ex_GSE62452_normal <= 0)] <- NaN
ex_GSE62452_normal <- log2(ex_GSE62452_normal)}

# box-and-whisker plot
#dev.new(width=3+ncol(gset)/6, height=5)
par(mar=c(7,4,2,1))
#title <- paste ("GSE62452_normal", "/", annotation(gset), sep ="")
boxplot(ex_GSE62452_normal, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2,cex.axis=0.6)
#dev.off()

# expression value distribution plot
par(mar=c(4,4,2,1))
#title <- paste ("GSE62452_normal", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex_GSE62452_normal, legend=F)

# mean-variance trend
ex_GSE62452 <- na.omit(ex_GSE62452_normal) # eliminate rows with NAs
plotSA(lmFit(ex_GSE62452_normal), main="Mean variance trend, GSE62452")

# UMAP plot (multi-dimensional scaling)
ex_GSE62452_normal <- ex_GSE62452_normal[!duplicated(ex_GSE62452_normal), ]  # remove duplicates
ump <- umap(t(ex_GSE62452_normal), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6,col="blue")

ex_GSE62452_normal<-ex_GSE62452_normal[1:length(ex_GSE78229_tumor[,1]),]
merged_data<-cbind(ex_GSE78229_tumor,ex_GSE62452_normal)
saveRDS(merged_data,"merged_data.RDS")
merged_data<-readRDS("merged_data.RDS")


par(mar=c(7,4,2,1))
#title <- paste ("GSE78229", "/", annotation(gset_tumor), sep ="")
boxplot(merged_data, boxwex=0.7,  notch=T, main=title, outline=FALSE, las=2,cex.axis=0.6)
# expression value distribution plot
par(mar=c(4,4,2,1))
#title <- paste ("GSE78229", "/", annotation(gset), " value distribution", sep ="")
#plotDensities(merged_data, legend=F)

merged_data <- merged_data[!duplicated(merged_data), ]  # remove duplicates
ump <- umap(t(merged_data), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

ex_GSE78229_tumor<-readRDS("ex_GSE78229_tumor.RDS")

ex_GSE62452_normal <- ex_GSE62452_normal[!duplicated(ex_GSE62452_normal), ]  # remove duplicates
ump <- umap(t(ex_GSE62452_normal), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="",yaxt="none",xaxt="none", pch=20, cex=1.5,col="darkgreen")
#library("maptools")  # point labels without overlaps
#pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6,col="blue")
par(new=TRUE)
ex_GSE78229_tumor <- ex_GSE78229_tumor[!duplicated(ex_GSE78229_tumor), ]  # remove duplicates
ump <- umap(t(ex_GSE78229_tumor), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="",pch=20, cex=1.5,col="purple")
#library("maptools")  # point labels without overlaps
#pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6,col="red")
legend("topleft", inset=.02, title="group",
       c("normal","tumor"), col=c("darkgreen", "purple"),pch=c(19,19), cex=0.8)
# # # set parameters and draw the plot




#raw_data=data.matrix(GSEdata,rownames.force = NA)
#rownames(raw_data)=geneIDs


##1) BEGUM e not: TEZ dosyasinda plots diye alt dosya olustur buradakileri pdf() fonksiyonu
### histogramGSE.pdf  adi altinda gen expresyonlarinin histogramini kaydedecek plots dosyasina
pdf("histogramGSEmerged_data2.pdf")
hist(merged_data2, col = "gray", main="GSE78229+GSE62452 - Histogram")
dev.off()



merged_data2<-readRDS("merged_data2.RDS")
##############You are here Jun 25, Friday########
######
pdf("boxPlot.pdf")
par(mar=c(7,4,2,1))
boxplot(merged_data2,main="tumor and normal tissues", las=2,cex.axis=0.6,
        col=c(rep("purple",50),rep("darkgreen",61)),outline=FALSE,border="black")
                             
legend("topright", inset=.02, legend=c("tumor","normal"),
                             fill=c("purple","darkgreen"),cex=0.6)
dev.off()

# Hierarchical clustering of the "samples" based on
# the correlation coefficients of the expression values
hc = hclust(as.dist(1-cor(merged_data2)))
pdf("HclustMerged_data2.pdf")
plot(hc, main="GSE78229+GSE62452 - Hierarchical Clustering",cex=0.6)
#dev.off()
#######################################
# Differential expression (DE) analysis
#######################################

######## Separate each of the conditions into two data2 frames
### Eda ve beg?m e not ?nceki ?al??ma 
### 2) ?imdi siz buray? bir tumor bir hasta seklinde ayr?mal?s?n?z

##ben tumor? duzelttim size normal hastalar? b?rak?yorum
tumor<-merged_data2[,1:50]
normal<-merged_data2[,51:111]
# Compute the means of the samples of each condition
tumor.mean = apply(tumor, 1, mean)
normal.mean = apply(normal, 1, mean)

head(tumor.mean)

head(normal.mean)

# Just get the maximum of all the means
limit = max(tumor.mean, normal.mean)

# Scatter plot
pdf("Scatter.pdf")
plot(normal.mean ~ tumor.mean, xlab = "tumor", ylab = "normal",
     main = "normal and tumor - Scatter", xlim = c(0, limit), ylim = c(0, limit))
# Diagonal line
abline(0, 1, col = "red")
dev.off()

# Compute fold-change (biological significance)
# Difference between the means of the conditions
fold = tumor.mean - normal.mean

# Histogram of the fold differences
pdf("hist.pdf")
hist(fold, col = "gray")
#dev.off()
# Compute statistical significance (using t-test)
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics

for(i in 1 : nrow(merged_data2)) { # For each gene : 
  x = tumor[i,] # tumor of gene number i
  y = normal[i,] # normal of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}

head(pvalue)

# Histogram of p-values (-log10)
#pdf("plots/pvalueHist.pdf")
hist(-log10(pvalue), col = "gray")
#dev.off()
# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
#pdf("plots/volcano.pdf")
plot(fold, -log10(pvalue), main = "normal vs tumor - Volcano")

fold_cutoff = 5
pvalue_cutoff = 0.05
abline(v = fold_cutoff, col = "blue", lwd = 3)
abline(v = -fold_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)
#dev.off()
# Screen for the genes that satisfy the filtering criteria

# Fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
fc_genes<-merged_data2[filter_by_fold,]

dim(fc_genes)

write.csv(fc_genes,"fc_genes.csv")


# P-value filter for "statistical" significance
filter_by_pvalue = pvalue <= pvalue_cutoff

pvalue_genes<-merged_data2[filter_by_pvalue,]

#pValue_genes<-apply(pvalue_genes,1,mean)
#pvalue_genes=data2.frame(rownames(pValue_genes),pValue_genes)

write.csv(pvalue_genes,"pvalue_genes.csv")


# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue

filtered = merged_data2[filter_combined,]
dim(filtered)

write.csv(filtered,"normalvstumor_DEGs_samples.csv")

head(filtered)

expres<-apply(filtered,1,mean)

DEGs=data.frame(rownames(filtered),expres)

write.csv(DEGs,"all_DEGs.csv")

pdf("DEGs.pdf")
plot(fold, -log10(pvalue), main = "normalvstumor - Volcano #2")
points (fold[filter_combined], -log10(pvalue[filter_combined]),
        pch = 16, col = "red")
dev.off()

# Highlighting up-regulated in red and down-regulated in blue
pdf("up&down.pdf")
plot(fold, -log10(pvalue), main = "GSE155353 - Volcano #3")
abline(v = fold_cutoff, col = "blue", lwd = 3)
abline(v = -fold_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)
points (fold[filter_combined & fold < 0],
        -log10(pvalue[filter_combined & fold < 0]),
        pch = 16, col = "red")
points (fold[filter_combined & fold > 0],
        -log10(pvalue[filter_combined & fold > 0]),
        pch = 16, col = "blue")
legend(-1.5,150,c("NO","down", "up"),cex=.8,
       col=c("black","red","blue"),pch=c(21,19,19),title="Genes")
dev.off()


upRegulated_DEGs = merged_data2[filter_combined & fold > 0,]
dim(upRegulated_DEGs)

expres<-apply(upRegulated_DEGs,1,mean)
upRegulated_DEGs=data.frame(rownames(upRegulated_DEGs),expres)
write.csv(upRegulated_DEGs,"upRegulated_DEGs.csv")

downRegulated_DEGs=  merged_data2[filter_combined & fold < 0,]
dim(downRegulated_DEGs)

# Cluster the rows (genes) & columns (samples) by correlation

expres<-apply(downRegulated_DEGs,1,mean)
downRegulated_DEGs=data.frame(rownames(downRegulated_DEGs),expres)
write.csv(downRegulated_DEGs,"downRegulated_DEGs.csv")

library(gplots)
pdf("heatmap_DEGs.pdf")
heatmap.2(filtered, scale = "none", col = bluered(10000), 
                    trace = "none", density.info = "none",
                    cexCol = 0.4,cexRow=0.4)
dev.off()

par(mar = c(7, 4, 2, 2) + 0.2) #add room for the rotated labels
end_point = 0.5 + nrow(top10) + nrow(top10) - 1 #this is the line which does the trick (together with barplot "space = 1" parameter)

pdf("barplot_pathways.pdf")
par(mar = c(7, 4, 2, 2) + 0.2) 
barplot(top10$Count,las=2,cex.axis=0.2,
        col=c(rep("purple",10),
        rep("darkgreen",10),rep("blue",10),rep("red",10)),ylim=c(0,200),horiz=TRUE,
        names.arg = top10$Term)
dev.off()

terms<-top10$Term
term2=gsub("GO:[0-9]{7}~","",terms)
top10$Term<-term2

ggsave("Goplots.pdf")#width=15,height = 10,units="cm")
ggplot(data=top10, aes(x=Term, y=Count, fill = Category)) +
  geom_bar(stat="identity")+
  scale_colour_gradient2()+
  coord_flip()+
  ylim(0, 180)+
  scale_x_discrete(limits = top10$Term)+
  theme_classic()
dev.off()
#text(seq(2, end_point, by = 2), par("usr")[3]-0.25, 
  #   srt = 60, adj = 1, xpd = TRUE,
   #  labels = paste(top10$Term), cex = 0.65)
#############################################################################

