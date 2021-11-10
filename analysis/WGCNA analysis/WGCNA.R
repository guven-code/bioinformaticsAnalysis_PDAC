


#Data<-read.csv("mergedData2.csv")
Data<-read.csv("all_DEGs.csv")
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator
library(WGCNA)  
library(dplyr)
#data <- tibble::rownames_to_column(testData, "geneID")
dim(Data);
names(Data);

datExpr0 = as.data.frame(t(Data[, -1]))
names(datExpr0) = Data$GeneId                      #Use the first column name after '$'
rownames(datExpr0) = names(Data)[-1]

# Checking data for excessive missing values and identification of outlier microarray samples
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

input_mat = datExpr0[,length(Data[1,])]
input_mat=df
#exp_mat=t(data[1:length(Data[1,]),])
names(input_mat)

allowWGCNAThreads() 

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;
#sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


softPower = 5;
adjacency = adjacency(input_mat, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#input_mat


temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = softPower,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

netwk$colors[netwk$blockGenes[[1]]]
table(netwk$colors)

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]

write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

module_df<-read.csv("gene_modules.csv")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

MEs1<-MEs0[c(2,85),]
# Add treatment names
#MEs0$treatment = row.names(MEs0)[c(2,62)]

# tidy & plot data
mME = MEs1 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() + 
  theme_bw() +
  geom_text(aes(label = round(value, 2)))+
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

modules_of_interest = c("blue","grey", "turquoise", "yellow","brown")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes


subexpr = Data[submod$gene_id,]
geneid<-module_df$gene_id
submod_df = data.frame(subexpr) %>%
  mutate(
    GeneId = row.names(.)
  ) %>%
  pivot_longer(-GeneId) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=GeneId)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")


bluemod<-module_df[module_df$colors=="blue",]
blue=length(bluemod[,1])
greymod<-module_df[module_df$colors=="grey",]
grey=length(greymod[,1])
yellowmod<-module_df[module_df$colors=="yellow",]
yellow=length(yellowmod[,1])
turqmod<-module_df[module_df$colors=="turquoise",]
turquoise=length(turqmod[,1])
brownmod<-module_df[module_df$colors=="brown",]
brown=length(brownmod[,1])

t<- c(80,76,52,93,76)

barplot(t, main="Distribution of Significant Modules",xlim=c(0,10),ylim= c(0,100),
      col = c("blue", "gray", "yellow", "turquoise","brown"), 
      xlab = "modules")


modules<-read.csv("gene_modules.csv")
down<-read.csv("downRegulated_DEGs.csv")
up<-read.csv("upregulated_DEGs.csv")
all<-read.csv("all_DEGs.csv")

turq<-modules[modules$blue=="turquoise",]

turq_expr<-all[turq$ACTB,]

up_expr<-up[turq$ACTB,]
down_expr<-down[turq$ACTB,]

down_expr<-na.omit(down_expr)
up_expr<-na.omit(up_expr)



write.csv(down_expr,'down_expr_module.csv')
write.csv(up_expr,'up_expr_module.csv')
write.csv(turq_expr,'Degs_expr_module.csv')


