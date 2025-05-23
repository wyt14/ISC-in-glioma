# load("/data1/guanchenchen/2.cancer_subtype/0.LGG/Rdata/rank.Rdata")
load("/data1/guanchenchen/2.cancer_subtype/0.LGG/TCGA-LGG/Rdata/rank.5.Rdata")
setwd("/data/wangyuting/8.glioma/LGG")
index <- extractFeatures(nmf.rank,"max")
sig.order <- unlist(index)
NMF.Exp.rank <- nmf.rank[sig.order,]
NMF.Exp.rank <- na.omit(NMF.Exp.rank)
group <- predict(nmf.rank)
pdf("consensus.pdf")
consensusmap(nmf.rank,labRow = NA,labCol = NA,annCol = data.frame("cluster"=group[colnames(NMF.Exp.rank)]))
dev.off()



label_data <- read.table("LGG.IDH.label.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(label_data) <- c("Sample", "IDH_Status")

cluster_df <- data.frame("cluster" = group[colnames(NMF.Exp.rank)])
label_df <- label_data  # 包含Sample和IDH_Status列

clinical_data <- read.table("/data/wangyuting/8.glioma/LGG/clinical.project-tcga-lgg.2025-04-25/clinical.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
sex_data <- clinical_data[,c(10,22)]


cluster_df <- cbind(rownames(cluster_df),cluster_df[,1])
colnames(cluster_df) <- c("Sample","cluster") 
colnames(sex_data)[1] <- "Sample"  
sex_data_unique <- sex_data[!duplicated(sex_data$Sample), ]

merged_data <- merge(cluster_df, label_df, by = "Sample", all = FALSE)  # all=FALSE 只保留共有样本
final_data <- merge(merged_data, sex_data_unique, by = "Sample", all = TRUE)


ann_colors <- list(
  IDH_Status = c("IDH-mut" = "black", "IDH-wt" = "grey"),
  demographic.gender = c("male" = "blue", "female" = "pink")  # 连续变量用colorRampPalette
)

pdf("consensus.pdf")
consensusmap(nmf.rank,labRow = NA,labCol = NA,annCol = final_data,annColors = ann_colors )
dev.off()



setwd("/data/wangyuting/8.glioma/LGG")
library(ggplot2)
library(pheatmap)
load("/data1/guanchenchen/2.cancer_subtype/1.LGG.methylation/Rdata/DMP.Rdata")
load("/data1/guanchenchen/2.cancer_subtype/1.LGG.methylation/Rdata/LGG.methylation.group.Rdata")
load("/data1/guanchenchen/2.cancer_subtype/1.LGG.methylation/Rdata/Norm.Rdata")
load("/data1/guanchenchen/2.cancer_subtype/1.LGG.methylation/Rdata/DMP_DE.Rdata")
sp <-as.data.frame(c(rep("EIC",length(EIC)),rep("Rest",length(other))))
myNorm <- myNorm[,c(EIC,other)]
table(colnames(myNorm)==c(EIC,other))
rownames(sp) <- colnames(myNorm)
myNorm <- myNorm[rownames(DMP_DE),]
colnames(sp) <- 'Cluster'
bk <- seq(0,1,by = 0.01)
p5 <- pheatmap(as.matrix(myNorm),fontsize_col = 5,clustering_distance_rows = "euclidean",
         show_rownames = F,show_colnames = F,clustering_distance_cols = "euclidean",
         clustering_method = "complete",cluster_cols = T,cluster_rows = T,
         annotation_legend = T,scale = "none",annotation_col = sp,annotation_colors = list(Cluster= c(EIC="#CC2626",Rest="#3232A2")),
         color = c(colorRampPalette(colors = c("#3232A2","#FFFFFF"))(length(bk)/2),colorRampPalette(colors = c("#FFFFFF","#CC2626"))(length(bk)/2)))
pdf("m3.pdf")
p5
dev.off()


label_df <- read.table("LGG.IDH.label.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
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