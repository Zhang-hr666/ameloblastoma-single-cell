#inferCNV
# setwd("GBM_Recur/")
rm(list = ls())
library(tidyverse)
library(infercnv)
library(Seurat)
scRNA <- qs::qread("tmp/03-3.AMall_and_subepi.seurat.qs")
###table(scRNA$celltype)
table(scRNA$subcelltype)
###Idents(scRNA) <- "celltype"
Idents(scRNA) <- "subcelltype"
scRNA <- subset(scRNA,idents = c("E00","E01","E02","E03","E04","E05","E06","T/NK"))
scRNA <- subset(scRNA, downsample = 200)

counts_matrix = Seurat::GetAssayData(scRNA, slot="counts")

cellanno = FetchData(scRNA, vars = "subcelltype" ) %>% tibble::rownames_to_column(var = "cellbarcode")
write.table(cellanno, "output/cnv_cellanno.txt", sep = "\t", col.names = F,row.names =FALSE, quote =F )
head(cellanno)

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,  # 可以直接提供矩阵对象
                                    annotations_file="./output/cnv_cellanno.txt",
                                    delim="\t",
                                    gene_order_file="./tmp/hg38_gencode_v27.txt",
                                    ref_group_names=c("T/NK"))  # 用于设置参考组，正常的细胞类型组别
# qs::qsave(infercnv_obj, file = "tmp/infercnv_obj.qs")

library(infercnv)
options(scipen = 100)
options(error = function() traceback(2))
# 
# #options(bitmapType="Xlib")
# 
# infercnv_obj <- qs::qread("tmp/infercnv_obj.qs")

infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir= "infercnv_output_EIP",  # dir is auto-created for storing outputs
                              cluster_by_groups=T ,   # cluster
                              hclust_method="ward.D2",
                              plot_steps=F,
                              denoise = TRUE,
                              output_format="pdf")
qs::qsave(infercnv_obj2, file = "tmp/infercnv_obj_10000.qs")




# infer_CNV_obj <- qs::qread("tmp/infercnv_EPI_ALL.qs")
# infer_CNV_obj<-readRDS('./infercnv_output_EPI_ALL/run.final.infercnv_obj')
# expr<-infer_CNV_obj@expr.data
# expr[1:4,1:4]
# data_cnv<-as.data.frame(expr)
# dim(expr)
# colnames(data_cnv)
# rownames(data_cnv)
# sce.all.int <- qs::qread("tmp/OSCC_annotated_STref.qs")
# sce.all.int
# sce.all.int=sce.all.int[,colnames(sce.all.int) %in% rownames(phe)]
# sce.all.int
# identical(colnames(sce.all.int), rownames(phe))
# sce.all.int$celltype = phe$celltype
# table(sce.all.int$celltype )
# meta <- sce.all.int@meta.data