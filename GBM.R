load("/data1/guanchenchen/2.cancer_subtype/4.GBM/TCGA-GBM/Rdata/rank.5.Rdata")

setwd("/data/wangyuting/8.glioma/GBM")
index <- extractFeatures(nmf.rank,"max")
sig.order <- unlist(index)
NMF.Exp.rank <- nmf.rank[sig.order,]
NMF.Exp.rank <- na.omit(NMF.Exp.rank)
group <- predict(nmf.rank)
pdf("consensus.pdf")
consensusmap(nmf.rank,labRow = NA,labCol = NA,annCol = data.frame("cluster"=group[colnames(NMF.Exp.rank)]))
dev.off()



label_data <- read.table("GBM.IDH.label.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(label_data) <- c("Sample", "IDH_Status")

cluster_df <- data.frame("cluster" = group[colnames(NMF.Exp.rank)])
label_df <- label_data  

common_samples <- intersect(rownames(cluster_df), label_df$Sample)
cluster_common <- cluster_df[common_samples, , drop = FALSE]  
label_common <- label_df[match(common_samples, label_df$Sample), ]
names <- label_common[,1]
label_common <- label_common[,2]
names(label_common) <- names

pdf("consensus.pdf")
consensusmap(nmf.rank,labRow = NA,labCol = NA,annCol = label_common)
dev.off()





setwd("/data/wangyuting/8.glioma/GBM")
library(ggplot2)
library(pheatmap)
load("/data1/guanchenchen/2.cancer_subtype/7.GBM.methylation/Rdata/DMP.Rdata")
load("/data1/guanchenchen/2.cancer_subtype/7.GBM.methylation/Rdata/GBM.methylation.group.Rdata")
load("/data1/guanchenchen/2.cancer_subtype/7.GBM.methylation/Rdata/Norm.Rdata")
load("/data1/guanchenchen/2.cancer_subtype/7.GBM.methylation/Rdata/DMP_DE.Rdata")
sp <-as.data.frame(c(rep("EIC",length(EIC)),rep("Rest",length(other))))
myNorm <- myNorm[,c(EIC,other)]
table(colnames(myNorm)==c(EIC,other))
rownames(sp) <- colnames(myNorm)
myNorm <- myNorm[rownames(DMP_DE),]
colnames(sp) <- 'Cluster'
bk <- seq(0,1,by = 0.01)

label_df <- read.table("GBM.IDH.label.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(label_df) <- c("Sample", "IDH_Status")
common_samples <- intersect(colnames(as.matrix(myNorm)), label_df$Sample)
label_common <- label_df[match(common_samples, label_df$Sample), ]
names <- label_common[,1]
label_common <- label_common[,2]
names(label_common) <- names

matched_labels <- label_common[colnames(as.matrix(myNorm))]

annotation_df <- data.frame(
  Cluster = sp$Cluster,
  IDH_Status = matched_labels,
  row.names = colnames(as.matrix(myNorm))
)

annotation_colors <- list(
  Cluster = c(EIC = "#CC2626", Rest = "#3232A2"),
  IDH_Status = c("IDH-mut" = "black", "IDH-wt" = "grey")
)

p5 <- pheatmap(as.matrix(myNorm),
               fontsize_col = 5,
               clustering_distance_rows = "euclidean",
               show_rownames = F,
               show_colnames = F,
               clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               cluster_cols = T,
               cluster_rows = T,
               annotation_legend = T,
               scale = "none",
               annotation_col = annotation_df,
               annotation_colors = annotation_colors,
               color = c(colorRampPalette(colors = c("#F4F6FD","#FFFFFF"))(length(bk)/2),
                        colorRampPalette(colors = c("#FFFFFF","#FF8080"))(length(bk)/2)))

pdf("m3_with_IDH.pdf")
print(p5)
dev.off()
