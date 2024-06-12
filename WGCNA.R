################################################
################################################
### 作者：果子
### 更新时间：2020-05-22
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/


### 教程地址
### 文章链接
### https://mbio.asm.org/content/7/3/e00027-16.abstract
### 脚本：https://github.com/iscb-dc-rsg/2016-summer-workshop/blob/master/3B-Hughitt-RNASeq-Coex-Network-Analysis/tutorial/README.md
### 视频: https://www.youtube.com/watch?v=OdqDE5EJSlA
### PPT：http://khughitt.github.io/2016-iscb-dc-rsg-workshop-presentation/#1

rm(list = ls())
### 安装R包 
library('gplots')
library('ggplot2')
library('limma')
library('reshape2')
library('RColorBrewer')
library('WGCNA')
set.seed(1)
### 样本信息
samples <- data.table::fread('sample_metadata.csv',data.table = F)
### 表达量信息
raw_counts <- data.table::fread('count_table.csv', data.table = F)
rownames(raw_counts) <- raw_counts[,1]
raw_counts <- raw_counts[,-1]
head(raw_counts)
dim(raw_counts)
### 基因注释
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if(!require("Homo.sapiens")) BiocManager::install("Homo.sapiens",update = F,ask = F)
library('Homo.sapiens')
keytypes(Homo.sapiens)
columns(Homo.sapiens)
# Example query
gene_ids <- head(keys(Homo.sapiens, keytype='ENSEMBL'), 2)
select(Homo.sapiens, keytype='ENSEMBL', keys=gene_ids, columns=c('ALIAS', 'TXCHROM', 'TXSTART', 'TXEND'))
#####################################################################################
### 数据预处理
### 1.样本检测
# add a colorbar along the heatmap with sample condition
### 配色
num_conditions <- nlevels(as.factor(samples$condition))
pal <- colorRampPalette(brewer.pal(num_conditions, "Set1"))(num_conditions)
cond_colors <- pal[as.integer(as.factor(samples$condition))]

heatmap.2(cor(raw_counts), 
          RowSideColors=cond_colors,
          trace='none', 
          main='Sample correlations (raw)')
### 2.低counts数据过滤
# Remove all rows with less than n counts across all samples, where n=#samples
low_count_mask <- rowSums(raw_counts) < ncol(raw_counts)
raw_counts <- raw_counts[!low_count_mask,]
sprintf("Removing %d low-count genes (%d remaining).", sum(low_count_mask), sum(!low_count_mask))
### 3.log2转换
log_counts <- log2(raw_counts + 1)
### 数据变形
x = melt(as.matrix(log_counts))
colnames(x) = c('gene_id', 'sample', 'value')

ggplot(x, aes(x=value, color=sample)) + 
  geom_density()
### 未转化的
x1 = melt(as.matrix(raw_counts))
colnames(x1) = c('gene_id', 'sample', 'value')
ggplot(x1, aes(x=value, color=sample)) + 
  geom_density()

heatmap.2(cor(log_counts), 
          RowSideColors=cond_colors,
          trace='none', main='Sample correlations (log2-transformed)')

### 去除没有改变的基因
### 4.Remove non differentially-expressed genes
# first, let's remove any genes with _zero_ variance since these are not
# going to help us, and may cause problems with some of the models
log_counts <- log_counts[apply(log_counts, 1, var) > 0,]

# create design matrix for differential expression analysis;
# if you wanted to account for batch here, you could simply include a batch
# term in the linear model at this step, e.g.:
# mod <- model.matrix(~0+samples$condition+samples$batch)
mod <- model.matrix(~0+samples$condition)
# make model terms easier to work with
colnames(mod) <- levels(as.factor(samples$condition))
fit <- lmFit(log_counts, design=mod)
# generate a list of all possible pairwise contrasts
### 产生比较的组合
condition_pairs <- t(combn(levels(as.factor(samples$condition)), 2))                                                                                                                               
condition_pairs
comparisons <- list()                                                                                                                                          
for (i in 1:nrow(condition_pairs)) {                                                                                                                                     
  comparisons[[i]] <- as.character(condition_pairs[i,])                                                                                                      
}    
# vector to store de genes
sig_genes <- c()
# iterate over the contrasts, and perform a differential expression test for
# each pair
for (conds in comparisons) {
  # generate string contrast formula, "infLM24 - infLM4"
  contrast_formula <- paste(conds, collapse='-')
  
  contrast_mat <- makeContrasts(contrasts=contrast_formula, levels=mod)
  contrast_fit <- contrasts.fit(fit, contrast_mat)
  eb <- eBayes(contrast_fit)
  
  # Grab highly ranked genes; this is a pretty stringent p-value cutoff, but
  # it serves to limit the total number of genes we will use for this
  # tutorial
  sig_genes <- union(sig_genes, 
                     rownames(topTable(eb, number=Inf, p.value=0.005)))
}

# Filter out genes which were not differentially expressed for any contrast
log_counts <- log_counts[rownames(log_counts) %in% sig_genes,]

##########################################################################
### 共表达网络构建
cordist <- function(dat) {
  cor_matrix  <- cor(t(dat))
  
  dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))
  dist_matrix <- log1p(dist_matrix)
  dist_matrix <- 1 - (dist_matrix / max(dist_matrix))
  
  sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)
}
sim_matrix <- cordist(log_counts)

### 截取部分数据检测
heatmap_indices <- sample(nrow(sim_matrix), 500)

heatmap.2(t(sim_matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA, 
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Similarity matrix',
          density.info='none', revC=TRUE)

########################################################
### Construct adjacency matrix
# Construct adjacency matrix
library(WGCNA)
adj_matrix <- adjacency.fromSimilarity(sim_matrix, power=12, type='signed')

# Delete similarity matrix to free up memory
rm(sim_matrix)
gc()
# Convert to matrix
gene_ids <- rownames(adj_matrix)

adj_matrix <- matrix(adj_matrix, nrow=nrow(adj_matrix))
rownames(adj_matrix) <- gene_ids
colnames(adj_matrix) <- gene_ids
heatmap.2(t(adj_matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA, 
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Adjacency matrix',
          density.info='none', revC=TRUE)

##########################################################################
### 检测模块:Co-expression module detection
# Cluster gene expression profiles; the flashClust function from
# the authors of WGCNA is another options for larger datasets.
# For input, we use the reciprocal of the adjacency matrix; hierarchical
# clustering works by comparing the _distance_ between objects instead of the
# _similarity_.
gene_tree <- hclust(as.dist(1 - adj_matrix), method="average")

# we will use the cuttreeDynamicTree method to break apart the hc dendrogram
# into separate modules
module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=15,
                                   deepSplit=TRUE)

# assign a color to each module for easier visualization and referencing
module_colors <- labels2colors(module_labels)

########################################################################
sizeGrWindow(12, 9)

table(module_colors)
# Plot the dendrogram and the module colors underneath

plotDendroAndColors(gene_tree, 
                    module_colors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


module = "blue"
moduleGenes = module_colors==module
dd <- log_counts[moduleGenes,]

library(dplyr)
library(tibble)
library(tidyr)


samples$condition <- factor(samples$condition,levels = c("infLM4","infLM24","infLM48","infLM72"))
test <- dd %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  mutate(group = samples$condition) %>% 
  select(group,everything()) %>% 
  pivot_longer(cols = 3:ncol(.),
               names_to = "gene",
               values_to = "expression") 


library(ggplot2)
ggplot(test,aes(x = group,y = expression)) +
  geom_boxplot(aes(color=group)) +
  geom_jitter()
  labs(x = "Days after injury", y = "Expression Level", title = "black") + 
  theme(legend.position = "none")
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  theme(axis.line = element_line(colour = "grey")) + 
  theme(axis.ticks = element_blank()) 
