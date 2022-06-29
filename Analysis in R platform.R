
# Analysis in R platform
# R v4.0.3
# Tximport import quantitative files of kallisto output
setwd("../files")
library(tximport)
library(tximeta)
base_dir = "../kallisto_out"
sample_id = dir(file.path(base_dir))
files = file.path(base_dir, sample_id, "abundance.h5")
names(files) = paste0(sample_id)
tx2gene = read.table(file = "tx2gene", header = TRUE)
# "tx2gene" : Documents corresponding to transcripts id and genes id
all_sample_quant.tsv = tximport(files, type="kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene=tx2gene, ignoreTxVersion = T)


# Data were normalized and VST corrected using DESeq2
library(DESeq2)
s2c = read.csv("../files/design_matrix.csv", header = TRUE)
s2c = dplyr::mutate(s2c, path = files)
s2c$sample = as.factor(s2c$sample)
s2c$condition = as.factor(s2c$condition)
s2c$reps = as.factor(s2c$reps)
s2c$path = as.factor(s2c$path)
s2c$type = as.factor(s2c$type)

dds1 = DESeqDataSetFromTximport(all_sample_quant.tsv, s2c, design= ~ condition)
dds1 = DESeq(dds1)
vsd1 = varianceStabilizingTransformation(dds1, blind=FALSE)

# Batch correction of data using limma
library(limma)
g=factor(s2c$condition)
g=relevel(g,'BMYJ')
design=model.matrix(~g)
Normolized_Counts_after_limma_del_batch = limma::removeBatchEffect(assay(vsd1), batch = s2c$type, design = design)

# T-SNE analysis and plot
info = read.csv(file = "../files/design_matrix.csv", header = T)
data = Normolized_Counts_after_limma_del_batch
colnames(data) = info$sample
library(Rtsne)
data = as.data.frame(t(data))
data$class = info$condition
data = as.data.frame(data)
data_unique = unique(data) # Remove duplicates
data_matrix = as.matrix(data_unique[,1:17190]) # Retention of gene expression data
# Set a seed if you want reproducible results
set.seed(10)
tsne_out = Rtsne(data_matrix, pca=T, dims=2, perplexity=15, theta=0.8) # Run TSNE

# Visualize T-SNE analysis results using ggplot2
library(ggplot2)
library(ggforce)
tsne_res = as.data.frame(tsne_out$Y)
colnames(tsne_res) = c("tSNE1","tSNE2")
tsne_res$sample = info$sample
tsne_res$condition = info$condition
head(tsne_res)
P = ggplot(tsne_res, aes(tSNE1, tSNE2))+ 
  geom_point(aes(shape = condition, colour = condition), size = 2)+
  geom_mark_ellipse(aes(color = condition), expand = unit(0,"mm"))+  # Mark minimum ellipse
  scale_fill_manual(values=c('#ff595e', '#ffca3a', '#8ac926', "#1982c4", "#6a4c93",
                             "#006633", '#990033', '#6699CC', '#99c1b9', "#8e7dbe"))+
  scale_color_manual(values=c('#ff595e', '#ffca3a', '#8ac926', "#1982c4", "#6a4c93",
                              "#006633", '#990033', '#6699CC', '#99c1b9', "#8e7dbe"))+
  scale_shape_manual(values = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))+
  theme_test()+
  coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10))+
  geom_hline(yintercept = 0, lty=2, col="grey", lwd=0.75)+ 
  geom_vline(xintercept = 0, lty=2, col="grey", lwd=0.75)+
  xlab("t-SNE 1")+
  ylab("t-SNE 2")+
  labs(title = "t-SNE plot")

ggsave(P, filename = "t-SNE plot.pdf", width = 7, height = 7)

# Spearman correlation heatmap (no clustering)
manual_order = c("BY-1-YWXJ-three", "BY-4-YWXJ-three", "BY-5-YWXJ-three", "BY-6-YWXJ-three",
                 "BY-1-BZCJ-one", "BY-2-BZCJ-two", "BY-3-BZCJ-one", "BY-4-BZCJ-two", "BY-5-BZCJ-one", "BY-6-BZCJ-one",
                 "BY-1-XJ-two", "BY-2-XJ-two", "BY-3-XJ-one","BY-4-XJ-one", "BY-5-XJ-two", "BY-6-XJ-four",
                 "BY-1-YDJ-three", "BY-2-YDJ-one", "BY-3-YDJ-three", "BY-4-YDJ-three", "BY-5-YDJ-one", "BY-6-YDJ-one",
                 "BY-1-GSTJ-two", "BY-2-GSTJ-two", "BY-3-GSTJ-three", "BY-4-GSTJ-three", "BY-5-GSTJ-three", "BY-6-GSTJ-two",
                 "BY-1-JGQJ-one", "BY-2-JGQJ-two", "BY-3-JGQJ-one", "BY-4-JGQJ-three", "BY-5-JGQJ-two", "BY-6-JGQJ-two",
                 "BY-1-FCJ-one", "BY-2-FCJ-one", "BY-3-FCJ-one", "BY-4-FCJ-two", "BY-5-FCJ-two", "BY-6-FCJ-one",
                 "BY-1-BMYJ-three", "BY-2-BMYJ-three", "BY-3-BMYJ-three", "BY-4-BMYJ-one", "BY-5-BMYJ-two", "BY-6-BMYJ-two",
                 "BY-1-ZCSJ-three", "BY-2-ZCSJ-three", "BY-3-ZCSJ-three", "BY-4-ZCSJ-two", "BY-5-ZCSJ-two", "BY-6-ZCSJ-three",
                 "BY-1-ZDQJ-three", "BY-2-ZDQJ-two", "BY-3-ZDQJ-three", "BY-4-ZDQJ-three", "BY-5-ZDQJ-three", "BY-6-ZDQJ-three")
# Calculation of Spearman correlation coefficient matrix
sampleCor_spearman_matrix = cor(Normolized_Counts_after_limma_del_batch, method = "spearman")
sampleCor_spearman_matrix = as.matrix(sampleCor_spearman_matrix)
rownames(sampleCor_spearman_matrix) = paste(vsd1$sample, vsd1$type, sep="-")
colnames(sampleCor_spearman_matrix) = paste(vsd1$sample, vsd1$type, sep="-")
sampleCor_spearman_matrix = sampleCor_spearman_matrix[c(manual_order),]
sampleCor_spearman_matrix = sampleCor_spearman_matrix[,c(manual_order)]

# Spearman correlation heatmap plot
library(pheatmap)
library(RColorBrewer)
pdf(file = "sampleCor_spearman_heatmap.pdf", width = 9.0, height = 7.5)
pheatmap(sampleCor_spearman_matrix,
         breaks=seq(0.9, 1, 0.001),
         cluster_cols = F,
         cluster_rows = F,
         border_color= NA,
         show_rownames = T,
         show_colnames = T,
         main = "Correlation Heatmap of all Samples",
         display_numbers=F,
         color=colorRampPalette(c("blue", "lightblue", "red"))(100))+theme_test()
dev.off()