library(tidyverse)
library(ggpubr)
library(patchwork)
library(Signac)
library(Seurat)
set.seed(1234)
library(GenomeInfoDb)
library(ensembldb)  # library(EnsDb.Mmusculus.v79) # BiocManager::install("ensembldb")
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm39) # BiocManager::install('BSgenome.Mmusculus.UCSC.mm39')
# library(JASPAR2024)
library(TFBSTools)
library(motifmatchr)
library(LSD)

# monocle3 作图，去掉点的黑边
trace(plot_cells, edit = T) # 第301行 size = 1.5 * cell_size, 将1.5改成1，以去掉点的黑边

mm39E108 <- readRDS("/home/dell/app/genome/mm39/mm39E108_24388.rds")
dir1 <- "/media/dell/data/ATACsn/cellRanger/E18hGFAP_mm39/outs/"
sampleX <- 'E18Ctx-hGFAP'

counts <- Read10X_h5( paste0(dir1,"filtered_peak_bc_matrix.h5" ) )
# 过滤掉非standered chromosome, Y染色体，线粒体, 防止基因组片段命名引起的错误
aa <- rownames(counts) %>% as.data.frame() 
aa <- separate(aa, ., sep = ':', into = c('chr', 'start-end'), remove = F) %>% dplyr::filter(chr %in% standardChromosomes(BSgenome.Mmusculus.UCSC.mm39)[1:20] )
counts <- counts[aa$.,]


meta1 <- read.csv(file = paste0(dir1, "singlecell.csv" ), header = T, row.names = 1)
# frags1 <- CreateFragmentObject( path = paste0(dir1, "fragments.tsv.gz" ), cells = rownames(meta1))

brain_assay <- CreateChromatinAssay(counts = counts,  sep = c(":", "-"),  genome = "mm39", 
                                    fragments = paste0(dir1, 'fragments.tsv.gz' ),  min.cells = 10, min.features = 200)
brain <- CreateSeuratObject(counts = brain_assay,  assay = 'ATAC',  project = sampleX,  meta.data = meta1)
granges(brain)


Annotation(brain) <- mm39E108

brain <- NucleosomeSignal(brain)####Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to fragments < 147 bp (nucleosome-free)
brain$nucleosome_group <- ifelse(brain$nucleosome_signal > 4, 'NS > 4', 'NS < 4') 
FragmentHistogram(object = brain, group.by = 'nucleosome_group', region = 'chr1-1-10000000') 

brain <- TSSEnrichment(brain, fast = F)
brain$high.tss <- ifelse(brain$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(brain, group.by = 'high.tss') + NoLegend()

brain$'pct.peakFragment' <- brain$peak_region_fragments / brain$passed_filters * 100
VlnPlot(object = brain,  features = c('pct.peakFragment',  'peak_region_fragments', "nCount_ATAC", 'TSS.enrichment', 'nucleosome_signal'),
        pt.size = 0.1,  ncol = 3) & geom_hline(yintercept = c(2, 3,3000,5000))

# 过滤，参数谨慎设置
brain <- subset(x = brain, subset = pct.peakFragment > 45 & peak_region_fragments > 3000 & peak_region_fragments < 100000
                & TSS.enrichment > 2 & nucleosome_signal < 3)
brain

VlnPlot(object = brain,  features = c('pct.peakFragment', 'peak_region_fragments', "nCount_ATAC", "nFeature_ATAC", 'TSS.enrichment', 'nucleosome_signal'), pt.size = 0.1,  ncol = 6)
ggsave(paste0(sName, "_QC.png" ),plot = last_plot(), path = getwd(), scale = 1.5, width = 30, height = 10, units = "cm", dpi = 300, limitsize = TRUE)

brain <- RunTFIDF(brain, assay = 'ATAC', scale.factor = 10000, method = 1, idf = NULL)
brain <- FindTopFeatures(brain, min.cutoff = 'q0')
brain <- RunSVD(object = brain, n =50, reduction.key = "LSI_", reduction.name = "lsi", assay = 'ATAC' )



DepthCor(brain) ###LSI
# The first LSI component often captures sequencing depth
brain <- RunUMAP(object = brain,  reduction = 'lsi',  dims = 2:30)
brain <- FindNeighbors(object = brain,  reduction = 'lsi',  dims = 2:30)
brain <- FindClusters(object = brain, algorithm = 1, resolution = 0.2, verbose = T)
DimPlot(object = brain, label = TRUE, label.size = 6, pt.size = 0.6 ) + guides(color=guide_legend(ncol = 1, bycol=T, override.aes = list(size = 4)))

ggsave(paste0(sName, "_res1.2.png" ),plot = last_plot(), path = getwd(), scale = 1.2, width = 20, height = 18, units = "cm", dpi = 300, limitsize = TRUE)
CoveragePlot( object = brain, region = "Dlx1",  annotation = T,  peaks = F, extend.upstream = 20000, extend.downstream = 20000, heights = c(20,1)) & scale_fill_manual(values = color$col20)

gliaCell <- colnames(glia)
SplitFragments(glia)
# bash
sampleN='glia'	
sort -k1,1V -k2,2n -k3,3n ${sampleN}.bed > ${sampleN}_sort.bed
size='/home/dell/app/genome/mm39/mm39m31_size.txt'
bedtools genomecov -i ${sampleN}_sort.bed -g ${size} -bg > ${sampleN}.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n ${sampleN}.bedgraph > ${sampleN}_sort.bedgraph
bedGraphToBigWig ${sampleN}_sort.bedgraph ${size} ${sampleN}.bw
lanceotron callPeaks ${sampleN}.bw -t 4 -w 300 -c 0.001 -f ./peak_lance --format Web

Sample <- 'glia'
gliaCell <- read.csv('gliaCell.csv')[,2]
peak <- read.delim("/media/dell/data/ATACsn/analysis/E18hGFAP/glia/peak_lance/glia_L-tron.bed")
peak <- makeGRangesFromDataFrame(peak[1:3])
peak <- keepStandardChromosomes(peak, species = 'Mus_musculus', pruning.mode="coarse") # 查询支持的物种名
peak <- peak[width(peak) > 200]
write.table(peak,'glia_peak.bed', col.names = F, row.names = F, quote = F, sep = '\t')
# bash

mm39E108 <- readRDS("/home/dell/app/genome/mm39/mm39E108_24388.rds")
dir1 <- "/media/dell/data/ATACsn/cellRanger/E18hGFAP_mm39/outs/"

meta <- read.csv(file = paste0(dir1, "singlecell.csv" ), header = T, row.names = 1)[gliaCell,]
frags <- CreateFragmentObject( path = paste0(dir1, 'fragments.tsv.gz' ), cells = gliaCell )
counts <- FeatureMatrix(fragments = frags, features = peak, cells = gliaCell)

assay  <- CreateChromatinAssay(counts, fragments = frags, genome = "mm39")
brain <- CreateSeuratObject(counts = assay, assay = 'ATAC', project = Sample,  meta.data = meta )

ids <- c('RG', 'IGP', 'As', 'OPC', 'IN-OB')
names(ids) <- levels(brain)
brain <- RenameIdents(brain, ids)
brain$CellType <- brain@active.ident
saveRDS(brain, 'E18Ctx_hGFAP_glia.rds')


library(monocle3)
library(Signac)
library(Seurat)
library(tidyverse)
library(SeuratWrappers)
gfapx <- as.cell_data_set(gfaps) # as.CellDataSet()
gfapx <- cluster_cells(cds = gfapx, reduction_method = "UMAP")
cluX <- gfapx@colData$ident
names(cluX) <- rownames(gfapx@colData)
gfapx@clusters$UMAP$clusters <- cluX
gfapx@clusters$UMAP$clusters <- factor(clusters(gfapx), levels = c('RG', 'IGP', 'As', 'OPC', 'IN-OB'))

gfapx <- learn_graph(gfapx, close_loop = T, verbose = T, # trajectory_graph_color = 'steelblue',
                     learn_graph_control = list(minimal_branch_len= 10,  nn.k=12, prune_graph=T, # ncenter=400,
                                                nn.cores = 12, geodesic_distance_ratio = 0.1, maxiter = 50, nn.method ='annoy', orthogonal_proj_tip = F) ) 

p1 <- plot_cells(gfapx, reduction_method = "UMAP", color_cells_by = "cluster", label_cell_groups = T, group_label_size = 8,
                 label_branch_points = F, label_roots = T, cell_size = 0, cell_stroke = 3, alpha = 1, label_leaves = F,
                 trajectory_graph_segment_size = 1.2, graph_label_size = 4) & theme_void() &
  scale_color_manual(values =  c("#FF7185", "#fcaf17", "#1AA260", "#7ACAF2", "#DDA0DD"))


'#d71345', '#bb505d', '#CD853F', '#FF8C00','#FFD700'

gfapx <- order_cells(gfapx, reduction_method = "UMAP")
p2 <- plot_cells(gfapx, reduction_method = "UMAP", color_cells_by = "pseudotime", cell_size = 0, cell_stroke = 3, alpha = 1,
                 label_branch_points = F, label_leaves = F, label_roots = T,  
                 trajectory_graph_segment_size = 1.2, graph_label_size = 6) & theme_void() &
  scale_color_gradient2( low = '#FF7185', mid = '#FFD700', high ='#30B21A', midpoint = 5)

ggsave( paste0('E18hGFAP_snATAC_glia', "_trajectory.pdf"), plot = p1|p2 , path = getwd(), scale = 1.2, width = 20, height =  20, units = "cm", dpi = 300, limitsize = T)


library(universalmotif) # BiocManager::install("universalmotif")
library('chromVAR') # BiocManager::install('chromVAR')
library(motifmatchr) # BiocManager::install('motifmatchr')
library('BSgenome.Mmusculus.UCSC.mm39') 
# options(BioC_mirror = 'https://mirrors.tuna.tsinghua.edu.cn/bioconductor')
# BiocManager::install('BSgenome.Mmusculus.UCSC.mm39')
library(TFBSTools)

motifX <- convert_motifs( motifs = motifATAC$motif, class = 'TFBSTools-PFMatrix')
motifX <- do.call(TFBSTools::PFMatrixList, motifX) # 
brain <- AddMotifs(object = brain, genome = BSgenome.Mmusculus.UCSC.mm39, pfm = motifX )
# saveRDS(brain, 'E18Ctx_hGFAP_glia_motif.rds')

library(BiocParallel)
register(SerialParam())

brain <- RunChromVAR( object = brain, genome = BSgenome.Mmusculus.UCSC.mm39)
# 此过程中motif名字中的特殊字符 | _ 会被转换成- 会导致ATAC assay和chromvar assay中的motif名字不一致
DefaultAssay(brain) <- 'chromvar'

# 过滤NA  需要谨慎，会去除ATAC assay, 只保留chroVAR assay
# 为什么会产生NA？？？
FeatureX <- setdiff(rownames(brain), rownames(brain@assays$chromvar@data)[is.na(brain@assays$chromvar@data[,1]  )])
brain <- brain[FeatureX,]

# 寻找cluster的marker motif
DefaultAssay(brain) <- 'chromvar'
TFacti <- FindAllMarkers(brain, only.pos = T, mean.fxn = rowMeans,  fc.name = "avg_diff")
TFacti$'mlogQ' <- -log((TFacti$p_val_adj + 1.028776e-323), base = 10) 
TFacti <- TFacti %>% dplyr::group_by(cluster) %>% dplyr::filter(mlogQ > 10)

for (i in levels(brain)) {
  MotifPlot(brain, motifs = TFacti[TFacti$cluster==i,]$gene[1:100], assay = 'ATAC')  + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_blank()) + 
  scale_fill_manual(values = c("#1AA260", "#FF7185", "#fcaf17",  "#7ACAF2"))
  ggsave( paste0('E18hGFAP_snATAC_glia_motif_', i, ".pdf"), plot = last_plot(), path = getwd(), 
          scale = 1.2, width = 40, height =  20, units = "cm", dpi = 300, limitsize = T)
  
}




theme1 <- theme( plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'), text = element_text(size = 16), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), panel.background = element_rect(fill = 'transparent'), plot.background = element_rect(fill = 'transparent'))

# look at the activity of Lhx2
geneX <- motifATAC$meta$Name[motifATAC$meta$Symbol=='Gli3'][1]
p1 <- FeaturePlot( object = brain, features = geneX, min.cutoff = 'q10', max.cutoff = 'q95', pt.size = 1, ncol = 1, label = T, order = T) & 
  theme1 & theme(legend.position = c(0.9,0.3) ) & scale_y_reverse()
p2 <- VlnPlot(brain, features = geneX, pt.size = 0.1)  & theme(axis.title = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1, size = 14), plot.title = element_blank()) & NoLegend()
p3 <- MotifPlot(brain, motifs = geneX, assay = 'ATAC') & theme1 
ggarrange(p1, ggarrange(p3,p2, nrow = 1, widths = c(1,1)),ncol = 1, heights = c(1, 0.5))

ggsave( paste0('E18hGFAP_snATAC_glia_', geneX, ".pdf"), plot = last_plot(), path = getwd(), scale = 1.2, width = 10, height =  20, units = "cm", dpi = 300, limitsize = T)

MotifPlot( brain,  motifs = TFacti, assay = 'ATAC' )

motifData <- brain@assays$chromvar@data %>% t()
mon3 <- readRDS("/media/dell/xiaosu41/ATACsn/analysis/E18hGFAP/E18hGFAP-glia-snATAC-monocle3.rds")
mon3 <- readRDS()

geneA <- motifATAC$meta$Name[motifATAC$meta$Symbol=='Hes1'] 
geneA
geneX <- geneA[1]
mon3@colData[[geneX]] <- as.data.frame(motifData) [[geneX]]
p1 <- plot_cells(mon3, color_cells_by = geneX, cell_size = 3, label_cell_groups = , cell_stroke = 1, # min_expr = 2, scale_to_range = T,
           trajectory_graph_color = 'black', trajectory_graph_segment_size = 1,
           label_roots = F, label_leaves = F, label_branch_points = F) + theme_void() + theme1  + scale_y_reverse() + NoLegend()+ 
  scale_color_gradient2(low = 'transparent', mid = 'gray90', high = 'blue', midpoint = 0) 

p2 <- VlnPlot(brain, features = geneX, pt.size = 0.1, cols = c("#FF7185", "#fcaf17", "#1AA260", "#7ACAF2", "#DDA0DD"))  & 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1, size = 14), plot.title = element_blank()) & NoLegend()
p3 <- MotifPlot(brain, motifs = geneX, assay = 'ATAC') & theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_blank())
# ggarrange(p1, ggarrange(p3,p2, nrow = 1, widths = c(1,1)),ncol = 1, heights = c(1, 0.5))
ggarrange(p1, p3,p2,ncol = 1, heights = c(1, 0.3, 0.4))
ggsave( paste0('E18hGFAP_snATAC_glia_', geneX, ".pdf"), plot = last_plot(), path = getwd(), scale = 1, width = 10, height =  20, units = "cm", dpi = 300, limitsize = T)

for (i in colnames(motifData)[1220:2858]) {
  geneX <- i
  geneA <- gsub("[:|\\?|\\|]", "-", geneX)
  mon3@colData[[geneX]] <- as.data.frame(motifData) [[geneX]]
  p1 <- plot_cells(mon3, color_cells_by = geneX, cell_size = 1.5, label_cell_groups = , 
                   trajectory_graph_color = 'black', trajectory_graph_segment_size = 1,
                   label_roots = F, label_leaves = F, label_branch_points = F) + 
    theme_void() + theme1 + ggtitle(label = geneX) + theme(plot.title = element_text(size = 24)) + scale_y_reverse() + NoLegend()+ 
    scale_color_gradient2(low = 'gray85', mid = '#e4e4ff', high = 'blue', midpoint = 0) 
  
  p2 <- VlnPlot(brain, features = geneX, pt.size = 0.1, cols = c("#FF7185", "#fcaf17", "#1AA260", "#7ACAF2", "#DDA0DD"))  & theme_void() &
    theme(axis.line = element_line(linewidth = 0.5), axis.text = element_text(size = 10),  axis.text.x = element_text(angle = 45, hjust = 0.5, size = 14), plot.title = element_blank()) & NoLegend()
  p3 <- MotifPlot(brain, motifs = geneX, assay = 'ATAC') + theme_void() + theme(text = element_blank()) + scale_color_manual(values = c("#1AA260", "#FF7185", "#fcaf17",  "#7ACAF2") )
  # ggarrange(p1, ggarrange(p3,p2, nrow = 1, widths = c(1,1)),ncol = 1, heights = c(1, 0.5))
  ggarrange(p1, p3,p2,ncol = 1, heights = c(1, 0.3, 0.4))
  ggsave( paste0('E18hGFAP_snATAC_glia_', geneA, ".png"), plot = last_plot(), path = getwd(), scale = 1, width = 10, height =  24, units = "cm", dpi = 300, limitsize = T)
  
}





p1 <- plot_cells(mon3, reduction_method = "UMAP", color_cells_by = "cluster", label_cell_groups = T, group_label_size = 8,
                 label_branch_points = F, label_roots = T, cell_size = 0, cell_stroke = 3, alpha = 1, label_leaves = F,
                 trajectory_graph_segment_size = 1.2, graph_label_size = 4) & theme_void() & scale_y_reverse() &
  scale_color_manual(values =  c("#FF7185", "#fcaf17", "#1AA260", "#7ACAF2", "#DDA0DD")) 

 
p2 <- plot_cells(mon3, reduction_method = "UMAP", color_cells_by = "pseudotime", cell_size = 0, cell_stroke = 3, alpha = 1,
                 label_branch_points = F, label_leaves = F, label_roots = T,  
                 trajectory_graph_segment_size = 1.2, graph_label_size = 6) & theme_void() & scale_y_reverse() &
  scale_color_gradient2( low = '#FF7185', mid = '#FFD700', high ='#30B21A', midpoint = 5)


theme1 <- theme( plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'), text = element_text(size = 16), 
                 axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), 
                 axis.line.x = element_blank(), axis.line.y = element_blank(),
                 panel.background = element_rect(fill = 'transparent'), plot.background = element_rect(fill = 'transparent'))

# 线图
plotLine1 <- function( fitValue = pdTimeFit, gene1 = 'Olig2', linewidth = 2, textSize = 12, alpha = 0.8, legendSize = 4 ) {
  dataA <- data.frame( )
  for (i in names(pdTimeFit)){
  dataB <- pdTimeFit[[i]]$fitValue[, c('Time', gene1)]
  dataB[, 'lineage'] <- i
  dataA <- rbind(dataA, dataB)    
  }
  
  theme2 <- theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 24), 
        panel.background = element_rect(fill = 'white'), 
        axis.line = element_line(linewidth = 0.5),  
        text = element_text(size = textSize))
  
  ggplot( data = dataA ) + geom_line(aes_string(x = 'Time', y = gene1, color = 'lineage' ), linewidth = linewidth, alpha = alpha) + 
    ggtitle(label = gene1) + xlab(label = 'pseudotime') + ylab(label = 'expression levels') +
    theme2 + guides(color=guide_legend(ncol = 1, bycol=T, override.aes = list(linewidth = legendSize)))

}

gene1 <- 'Olig2'
pA <- plot_cells(gfapx, genes = gene1, cell_size = 0.5, scale_to_range = F, trajectory_graph_segment_size = 0.3,
           label_cell_groups = F, label_branch_points = F, label_roots = F, label_leaves = F)   +  
   theme1 + theme(text  = element_text(size = 24, face = 'bold') )  + NoLegend() + 
   scale_color_gradient(low = 'gray90',  high = 'blue') 

pB <- plotLine1( fitValue = pdTimeFit, gene1 = gene1, linewidth = 1.5, alpha = 0.8, legendSize = 4, textSize = 12) + theme(plot.title = element_blank())

#ggarrange(pA, pB, ncol = 1, heights = c(1,0.4) )

geneA <- colnames(motifData)[grep(gene1, colnames(motifData), ignore.case = T)]     
geneA
MotifPlot(brain, motifs = geneA, assay = 'ATAC') + theme_void() + theme(text = element_blank()) + 
  scale_fill_manual(values = c("#1AA260", "#FF7185", "#fcaf17",  "#7ACAF2") )
geneX <- geneA[1]
mon3@colData[[geneX]] <- as.data.frame(motifData)[[geneX]]


p1 <- plot_cells(mon3, color_cells_by = geneX, cell_size = 1.5, label_cell_groups = , 
                 trajectory_graph_color = 'black', trajectory_graph_segment_size = 0.25,
                 label_roots = F, label_leaves = F, label_branch_points = F) + 
  theme1 + ggtitle(label = geneX) + theme(plot.title = element_text(size = 20) ) + scale_y_reverse() + NoLegend()+ 
  scale_color_gradient2(low = 'gray85', mid = '#e4e4ff', high = 'blue', midpoint = 0) 

p2 <- VlnPlot(brain, features = geneX, pt.size = 0.1, cols = c("#FF7185", "#fcaf17", "#1AA260", "#7ACAF2", "#DDA0DD"))  & theme_void() &
  theme(axis.line = element_line(linewidth = 0.5), axis.text = element_text(size = 10),  axis.text.x = element_text(angle = 45, hjust = 0.5, size = 12), plot.title = element_blank()) & NoLegend()
p3 <- MotifPlot(brain, motifs = geneX, assay = 'ATAC') + theme_void() + theme(text = element_blank()) + scale_fill_manual(values = c("#1AA260", "#FF7185", "#fcaf17",  "#7ACAF2") )
# ggarrange(p1, ggarrange(p3,p2, nrow = 1, widths = c(1,1)),ncol = 1, heights = c(1, 0.5))
#ggarrange(p1, p3,p2,ncol = 1, heights = c(1, 0.3, 0.4))


ggarrange(ggarrange(pA, pB, ncol = 1, heights = c(1,0.4) ), ggarrange(p1, p3,p2,ncol = 1, heights = c(1, 0.3, 0.4)), nrow = 1 )

ggsave( paste0('E18hGFAP_snATAC_glia_', gene1, ".pdf"), plot = last_plot(), path = getwd(), scale = 1, width = 20, height =  20, units = "cm", dpi = 300, limitsize = T)





CoveragePlot(brain, region = gene1, assay = 'ATAC', extend.upstream = 1000, extend.downstream = 1000, peaks = F) &
  theme_void() & NoLegend() & theme(text = element_text(size = 14)) &
  scale_fill_manual(values = c("#FF7185", "#fcaf17", "#1AA260", "#7ACAF2", "#DDA0DD"))

ggsave( paste0('E18hGFAP_snATAC_glia_peak_', gene1, ".pdf"), plot = last_plot(), path = getwd(), scale = 1.2, width = 12, height = 8, units = "cm", dpi = 300, limitsize = T)







  
  
  
  
