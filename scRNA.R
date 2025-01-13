library(tidyverse)
library(Seurat)
library(scales)
library(ggpubr)
library(harmony)
library(patchwork)
library(Matrix)
library(ggforce)

# estimation and removal of cell free mRNA contamination in scRNA-seq data
library(SoupX)
library(Seurat)
library(tidyverse)
library(scales)
library(ggpubr)

dirSam <- list.dirs(path = "./cellRanger", recursive = F, full.names = F)
samName <- gsub("_mm39", "", dirSam)
dir.create("soupX")
dirMtx <- paste0("./cellRanger/", dirSam, "/outs/")

for (i in 1:length(dirSam)) {
  filtX <- Read10X( paste0(dirMtx[i], "filtered_feature_bc_matrix"))
  rawX <- Read10X( paste0(dirMtx[i], "raw_feature_bc_matrix") )
  gfaps <- CreateSeuratObject(counts = filtX)
  gfaps <- SCTransform(gfaps, method = "glmGamPoi",  assay = "RNA", new.assay.name = "SCT", variable.features.n = 5000, do.scale = T, do.center = T, return.only.var.genes = F, verbose = T)  
  gfaps <- RunPCA(gfaps, assay = DefaultAssay(gfaps), features = VariableFeatures(gfaps), npcs = 50)
  gfaps <- RunUMAP(gfaps, dims = 1:50)
  gfaps <- FindNeighbors(gfaps, dims = 1:50)
  gfaps <- FindClusters(gfaps, resolution = 0.5)
  soupA <- SoupChannel(rawX, filtX, metaData = NULL, calcSoupProfile = F)
  soupA <- estimateSoup(soupA, soupRange = c(1, 100), keepDroplets = F)
  soupA <- setClusters(soupA, setNames(gfaps$seurat_clusters, rownames(gfaps@meta.data)))
  soupA <- setDR(soupA, gfaps@reductions$umap@cell.embeddings)
  soupA <- autoEstCont(soupA)
  counts <- adjustCounts(soupA, roundToInt = T)
  saveRDS(counts, paste0("./soupX/", samName[i], "_soupX.rds"))
}


# filter low quality and Doublet cell

library(DoubletFinder)
### remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
### devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')

sampleX <- 'E13Ctx'
counts <- readRDS("./soupX/", sampleX, "_soupX.rds") )
gfaps <- CreateSeuratObject(counts = counts, project = sampleX, assay = "RNA", meta.data = NULL, min.cells = 3)
gfaps[["percent.mt"]] <- PercentageFeatureSet(gfaps, pattern = "^mt-")
gfaps[["percent.rb"]] <- PercentageFeatureSet(gfaps, pattern = "^Rp[sl]") # 
VlnPlot(gfaps, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb" ), ncol = 4, cols = color$col7, group.by = "orig.ident") & geom_hline(yintercept = 8)
gfaps <- subset(gfaps, subset = nFeature_RNA > 1300 & nFeature_RNA < 6000 & nCount_RNA > 2000 & nCount_RNA < 30000 & percent.mt < 10) # 选择每个样本的过滤参数

gfaps <- NormalizeData(gfaps)
gfaps <- CellCycleScoring(gfaps, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = T, slot = "data")
gfaps <- SCTransform(gfaps, vst.flavor = "v2", assay = "RNA", new.assay.name = "SCT", 
                     variable.features.n = 4000,
                     do.scale = T, do.center = T, return.only.var.genes = T, verbose = T,
                     vars.to.regress = c( "S.Score", "G2M.Score") ) # 
gfaps <- RunPCA(gfaps, assay = DefaultAssay(gfaps), features = VariableFeatures(gfaps), npcs = 50, reduction.name = "pca") 
nPC <- 50
gfaps <- RunUMAP(gfaps, dims = 1:nPC, reduction = "pca", reduction.name = "umap")
gfaps <- FindNeighbors(gfaps, dims = 1:nPC,  assay = "SCT", reduction = "pca" )
gfaps <- FindClusters(gfaps, resolution = 0.5, algorithm = 1, n.iter = 30)
DimPlot(gfaps, reduction = "umap", label = T, label.size = 5, pt.size = 0.4, cols = color$col24)

###find pK
pNpK_List <- paramSweep_v3(gfaps, PCs = 1:40, sct = T )  # 耗时 用10000个细胞，进行pNpK参数扫描分析
pNpK_stats <- summarizeSweep(pNpK_List, GT = F )  # pN-pK生物学相关性
# pNpK_stats <- summarizeSweep(pNpK_List, GT = TRUE, GT.calls = classifications)  # GT ground-truth doublet classifications 
pK_BioCo <- find.pK(pNpK_stats) #  0.001-0.3每个pK的 mean BC, BC variance, and BCmvn scores for each pK value.
mpK <- as.numeric(as.vector(pK_BioCo$pK[which.max(pK_BioCo$BCmetric)])) # 生物学相关性BCmetric最大值 所对应的pK值
ggplot(data = pK_BioCo, mapping = aes(x = pK, y = BCmetric, group = 1)) + geom_point() + geom_line(linewidth =0.5, color = "#FD6C9E" ) + geom_text(data = pK_BioCo, aes(label = pK))

###find nExp
annotations <- gfaps@meta.data$seurat_clusters 
homoDoublet <- modelHomotypic(annotations)   # 同型双细胞比列 the proportion of homotypic doublets
doubletN <- round(ncol(gfaps)/1000 * 0.008 * ncol(gfaps)) # 根据比列计算双细胞数目， Assuming 7.6% doublet formation rate - tailor for your dataset
heterDoubletN <- round(doubletN*(1-homoDoublet))

###find doublet
gfaps <- doubletFinder_v3(gfaps, PCs = 1:33, pN = 0.25, pK = mpK, nExp = heterDoubletN, reuse.pANN = F, sct = T)
keep <- gfaps@meta.data[, paste("DF.classifications_0.25", mpK, heterDoubletN, sep = "_" )] == "Singlet"
gfaps <- gfaps[, keep]

#remove non neural cells
gfaps <- subset(gfaps, idents = c())

counts <- counts[, colnames(gfaps)]
saveRDS(counts, paste0("./doublet/", sampleX, "_RMdoublet.rds"))


# merge samples
library(Matrix)
E17a <- readRDS(paste0("./doublet/", sampleX, "_RMdoublet.rds"))
colnames(E17a) <- paste0("E17-a_", colnames(E17a))
colnames(E17b) <- paste0("E17-b_", colnames(E17b))
colnames(E18) <- paste0("E18_", colnames(E18))
````
counts <- merge(E17a, y = c(E17b, E18), by = "row.names", all = T, sort = F)  # 默认row.name = T 顺序会改变，会按字母排序
# counts <- merge(E17a, y = c(E17b, E18), by.x = "row.names", by.y = "row.names", all = T, sort = F, all.x = all, all.y = all) 
rownames(counts) <- counts$Row.names
counts <- counts[, -1]
counts <- as.matrix(counts)
counts <- Matrix(counts, sparse=TRUE)


# seurat analysis

dirCode <- "D:"
color <- readxl::read_excel( paste0(dirCode, "/R-code/color.xlsx"))
ccgenes <- read.csv(paste0(dirCode, "/R-code/ccgenes.csv") )
cc.genes$s.genes <- tolower(cc.genes$s.genes) %>% str_to_title()
cc.genes$g2m.genes <- tolower(cc.genes$g2m.genes) %>% str_to_title()

rmFeature <- list()
rmFeature[['mtGene']] <- c(rownames(counts)[grep("^mt-", rownames(counts))], rownames(counts)[grep("^Rp[sl]", rownames(counts))], ccgenes$gender[1:30], ccgenes$hbgenes[1:11], 'Prm1')
rmFeature[['ccGene']] <- c( ccgenes$cc.genes[1:240], ccgenes$Histone_mm39[1:36] )
# unwanted.feature <- c(rownames(gfaps)[grep("^mt-", rownames(gfaps))], rownames(gfaps)[grep("^Rp[sl]", rownames(gfaps))], ccgenes$gender[1:30], ccgenes$hbgenes[1:11], ccgenes$cc.genes[1:240], ccgenes$Histone_mm39[1:36] )


sampleX <- "hGFAP-FT"
dir.create
setwd("G:/scRNA_mm39/hGFAP-FT")

gfaps <- CreateSeuratObject(counts = counts, project = sampleX, assay = "RNA", meta.data = NULL, min.cells = 3)
gfaps[["percent.mt"]] <- PercentageFeatureSet(gfaps, pattern = "^mt-")
gfaps[["percent.rb"]] <- PercentageFeatureSet(gfaps, pattern = "^Rp[sl]") # 
VlnPlot(gfaps, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb" ), ncol = 4, cols = color$col7, group.by = "orig.ident") & geom_hline(yintercept = 8)
# gfaps <- subset(gfaps, subset = nFeature_RNA > 1300 & nFeature_RNA < 6000 & nCount_RNA > 2000 & nCount_RNA < 30000 & percent.mt < 10)

gfaps <- NormalizeData(gfaps)
gfaps <- CellCycleScoring(gfaps, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = T, slot = "data")
gfaps$ccDif  <- gfaps$S.Score - gfaps$G2M.Score

gfaps <- SCTransform(gfaps, vst.flavor = "v2", assay = "RNA", new.assay.name = "SCT", 
                     variable.features.n = 4000,
                     do.scale = T, do.center = T, return.only.var.genes = T, verbose = T,
                     vars.to.regress = c( "S.Score", "G2M.Score") ) # 

VariableFeatures(gfaps) <- setdiff(VariableFeatures(gfaps), c(rmFeature$mtGene, rmFeature$ccGene) )
gfaps <- RunPCA(gfaps, assay = DefaultAssay(gfaps), features = VariableFeatures(gfaps), npcs = 50 )
gfaps <- RunHarmony(gfaps, group.by.vars = "orig.ident", reduction.use = "pca", dim.use = 1:50, reduction.save = "harmony", max_iter = 20, verbose = T, plot_convergence = T, assay.use = "SCT" )
# Cluster the cells
nPC <- 50 # [gfaps@reductions$harmony@stdev > 1.5]
gfaps <- RunUMAP(gfaps, dims = 1:nPC, reduction = "harmony", reduction.name = "umap.har")
gfaps <- FindNeighbors(gfaps, dims = 1:nPC, reduction = "harmony")
gfaps <- FindClusters(gfaps, resolution = 0.8, algorithm = 2)
DimPlot(gfaps, reduction = "umap.har", label = T, label.size = 5, pt.size = 0.3, cols = color$col22) +labs(title = "har_res1.2") + guides( color = guide_legend(ncol = 1, bycol=T, override.aes = list(size = 4)))
