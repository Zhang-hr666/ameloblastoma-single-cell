rm(list = ls());gc()
library(CellChat)
library(patchwork)
library(Seurat)
source("scrpit/00utils.R")
library(Seurat)
library(tidyverse)
library(patchwork)
## 导入Seurat对象
seu <- qs::qread("data/03-3.AMsubepi_fib.seurat.qs")
table(seu$eip_fibsubcell)
####  02单个cellchat对象的构建  ####
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis下面的是用全部，上面是某一个过程
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor", key = "annotation")
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact", key = "annotation")####只用细胞接触
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
#####示例样本的输入
# data.input = data_humanSkin$data # normalized data matrix
# meta = data_humanSkin$meta # a dataframe with rownames containing cell mata data
# cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data

# # Prepare input data for CelChat analysis
# data.input = data.input[, cell.use]
# meta = meta[cell.use, ]
# # meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
# unique(meta$labels) # check the cell labels
# cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

#####自己样本的输入
cellchat <- createCellChat(object = seu, 
                           group.by = "subcelltype", assay = "RNA")
# cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "subcelltype") # set "labels" as default cell identity
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
# set the used database in the object
cellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
# future::plan("multisession", workers = 15) # do parallel 可以不加，
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, type = "triMean")

#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2023-11-03 06:32:38.834754]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2023-11-03 06:33:16.367342]"

cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)


####  03可视化方案  ####
##  所有细胞交互作用的可视化
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
####热图
gg1 <- netVisual_heatmap(cellchat,color.heatmap = c("#2166ac","#b2182b"))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.heatmap = c("#2166ac","#b2182b"))
gg1 + gg2


##  单个细胞发送的信号可视化
mat <- cellchat@net$weight
par(mfrow = c(3,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


##  特定通路信号可视化
cellchat@netP$pathways

pathways.show <- c("WNT") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",arrow.size=0.5)
##画所有Circle plot
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
for (i in 1:length(pathways.show.all)) {
  gg <- netVisual_aggregate(cellchat, signaling = pathways.show.all[i], layout = "circle",arrow.size=0.5)
  pdf(paste0(pathways.show.all[i], "_Circle.pdf"), width = 8, height = 8)
  print(gg)
  dev.off()
}
#####保存为PNG
pathways.show.all <- cellchat@netP$pathways
for (i in 1:length(pathways.show.all)) {
  print(i)
  gg <- netVisual_aggregate(cellchat, signaling = pathways.show.all[i], layout = "circle",arrow.size=0.5)
  png(paste0(pathways.show.all[i], "_Circle.png"), width = 2300, height = 2300, units = "px", res = 300)
  print(gg)
  dev.off()
}
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

##  将细胞簇分组为不同的细胞类型
# Chord diagram
group.cellType <- c("E00","E01", "E02", "E03",
 "E04","E05","E06","E07" ,"E08" ,"E09","E10" ,"E11"  ,"E12" ,'Endothelial' ,"F00","F01"         
,"F02", "F03" , "F04" , "Macrophages" ,"T/NK") # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level

##  信号通路中具有贡献的配受体可视化
netAnalysis_contribution(cellchat, signaling = pathways.show)

##  单个配体-受体对介导的细胞间通讯可视化
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(4,10) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle",arrow.size=1)

netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")


##  自动保存所有推断网络的图，以便快速探索
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(4,10)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}


plotGeneExpression(cellchat, signaling = "WNT")


####  04计算作用强度，识别信号主要的细胞来源  ####
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 12, height = 4, font.size = 10)


##  主要贡献来源的散点图可视化
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("FGF"))####可以放几个通路
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width = 9,height = 30)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width = 9,height = 30)
ht1 + ht2


# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("WNT"))
ht
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", signaling = c("WNT"))
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", signaling = c("WNT"))
ht1 + ht2

save(cellchat,file = './output/Cell-Cell Contact-ReceptorCellChat.Rdata')
####  05两组间比较  ####

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

humanSkin_ser_list <- SplitObject(seu,split.by = 'group')
##  构建对象
cellchat.Primary <- Cellchat_proceeding(object =humanSkin_ser_list$Primary , 
                                   group.by = "eip_fibsubcell", assay = "RNA",
                                   workers = 10)
cellchat.Recurrent <- Cellchat_proceeding(object = humanSkin_ser_list$Recurrent, 
                                   group.by = "eip_fibsubcell", assay = "RNA",
                                   workers = 10)


# rm(humanSkin_ser,humanSkin_ser_list,data_humanSkin);gc()
object.list <- list(Primary = cellchat.Primary, Recurrent = cellchat.Recurrent)
save(cellchat.Primary,cellchat.Recurrent,object.list,file = './output/05CellChat_list.Rdata')
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
#> An object of class CellChat created from a merged object with multiple datasets 
#>  593 signaling genes.
#>  7563 cells. 
#> CellChat analysis of single cell RNA-seq data!


##  比较交互总数和交互强度
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


##  不同细胞群之间的相互作用次数或相互作用强度的差异
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2


weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2], 
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


table(object.list$Primary@meta$eip_fibsubcell)
group.cellType <- c(rep("Eipthelial", 13), rep("Endothelial", 1), rep("Fibroblast", 5), rep("Marcrophages", 1), rep("T/NK", 1))
group.cellType <- factor(group.cellType, levels = c("Eipthelial", "Endothelial", "Fibroblast","Marcrophages","T/NK"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

##  对比看一下两组的强弱
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], 
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)



##  比较每个信号通路的整体信息流
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2


##  通过比较通信概率来识别功能失调的信号
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)



##  分别看增强和减弱的配受体对
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  
                        comparison = c(1, 2), max.dataset = 2, 
                        title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  
                        comparison = c(1, 2), max.dataset = 1, 
                        title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2




##  使用差异表达分析识别功能失调的信号转导

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Recurrent"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS and cellchat@var.features$LS.info. 
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, 
                                       features.name = features.name, only.pos = FALSE, 
                                       thresh.pc = 0.1, thresh.fc = 0.05, thresh.p = 1) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "Recurrent",ligand.logFC = 0.1, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "Primary",ligand.logFC = -0.05, receptor.logFC = -0.05)


##  进一步进行反卷积以获得单个信号基因
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


##  使用气泡图或和弦图可视化上调和下调的信号配体-受体对
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2


cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
plotGeneExpression(cellchat, signaling = "WNT", split.by = "datasets", colors.ggplot = T)

