4#Paper Figures
setwd("/data/pearce/group/Stanczak/Project_1713/")

#Libraries and functions####
#load libraries
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(clustree)
  library(ggplot2)
  library(plyr)
  library(scales)
  library(ggrepel)
  library(biomaRt)
  library(tidyverse)
  library(patchwork)
  library(future)
  library(reticulate)
  library(SeuratWrappers)
  library(destiny)
  library(gam)
  library(slingshot)
  library(TSP)
  library(pheatmap)
  library(reshape2)
  library(goseq)
  library(RColorBrewer)
  library(ggcorrplot)
  library(GO.db)
  library(KEGGREST)
})
plan(strategy = 'multiprocess', workers = 4)

#Functions
color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#In terminal####
# ##change the snakefile to the right path
# #login
# ssh maximus
# module load slurm snakemake cellranger
# cd /data/pearce/group/Stanczak/Project_1713
# 
# ##rum command
# snakemake -c "SlurmEasy -t {threads} --no_log" -s Snakefile_mouse_Michal -j 6

#Data####
files <- list.files(path = "/data/pearce/group/Stanczak/Project_1713/CellRanger", 
                    full.names = F, pattern = 'mouse_done_')
files <- gsub('mouse_done_','',files)

target <- files[12:16]

d <- list()
p <- list()
path <- "/data/pearce/group/Stanczak/Project_1713/CellRanger/mouse_"
for(i in target){
  #Read data
  t <- Read10X(data.dir = paste0(path,i,'/outs/filtered_feature_bc_matrix'))
  
  t <- CreateSeuratObject(counts = t, 
                          min.cells = 3,
                          min.features = 200)
  t$condition <- i
  
  #Get mitochondrial percentage
  t[['percent.mito']] <- PercentageFeatureSet(t, 
                                              pattern = '^mt-')
  p1 <- FeatureScatter(object = t,
                       feature1 = "nCount_RNA",
                       feature2 = "percent.mito",
                       pt.size = 0.2) + 
    NoLegend() +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE), 
                       n.breaks = 4)
  
  p2 <- FeatureScatter(object = t,
                       feature1 = "nCount_RNA",
                       feature2 = "nFeature_RNA",
                       pt.size = 0.2) + 
    NoLegend() +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE), 
                       n.breaks = 4)
  
  p[[i]] <- p1+p2+plot_layout(guides = 'collect') & 
    theme(plot.title = element_blank())
  
  d[[i]] <- t
}

ggsave(filename = "QC_mito_UMI.png",
       plot = wrap_plots(p, nrow = 2), 
       width = 400, height = 160,
       units = 'mm', dpi = 300)

dir.create('Individual')

p <- list()
for(i in target){
  p[[i]] <- FeatureScatter(object = t,
                       feature1 = "nCount_RNA",
                       feature2 = "percent.mito",
                       pt.size = 0.2) + 
    NoLegend() +
    geom_hline(yintercept = 10, linetype = 'dashed') +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE), 
                       n.breaks = 4)
}

wrap_plots(p, nrow = 2)

#Filtering, QC, scaling, normalizing and tissue integration####
max.f = 20000
for(i in target){
  #Filter
  d[[i]] <- subset(x = d[[i]],
                   subset = percent.mito <= 10)
  
  #Scale and normalize
  d[[i]] <- SCTransform(object = d[[i]], assay = "RNA")
  
  #Get max genes
  if(nrow(d[[i]]) < max.f) {max.f <- nrow(d[[i]])}
  
  #save files
  saveRDS(d[[i]], file = paste0('Individual/',i,'.rds'))
}

#Select integration features
features <- SelectIntegrationFeatures(object.list = d, 
                                      nfeatures = max.f)

#Run PCA on individual datasets
for(i in target){
  d[[i]] <- RunPCA(d[[i]], verbose = FALSE,
                   npcs = 50,
                   features = features)
}

#Prep objects for integration
max.size <- 310000*max.f #To correct the globals
options(future.globals.maxSize= max.size)
d <- PrepSCTIntegration(object.list = d, 
                        anchor.features = features)

#Integrate using rpca
anchors <- FindIntegrationAnchors(object.list = d,
                                  anchor.features = features,
                                  normalization.method = 'SCT',
                                  reduction = "rpca", 
                                  dims = 1:50)

d <- IntegrateData(anchorset = anchors, 
                   normalization.method = "SCT",
                   dims = 1:50)

#Perform dimentionality reduction and cluster####
d <- RunPCA(d, verbose = FALSE)
d <- RunUMAP(d, dims = 1:50)
d <- FindNeighbors(object = d, dims = 1:50)

res <- list(seq(0, 1.2, by=0.1))
for(i in 1:length(res)) {
  d <- FindClusters(object = d, resolution = res[[i]])
}

colnames(d@meta.data) <- gsub('([A-Za-z]+)(_)(snn_res.)([0-9.]+)',
                              "res_\\4",
                              colnames(d@meta.data))
clust <- d@meta.data[,grep('res_',colnames(d@meta.data))]

ggsave('Clustree.png',
       plot = clustree(clust, prefix = 'res_'),
       dpi = 300, height = 29, width = 21,
       units = 'cm')

#Select a resolution####
res <- 'res_1'
Idents(object = d) <- d[[res]]

#Make dim map and clustree
p1 <- DimPlot(object = d, reduction = "umap",
              pt.size = 0.1, label.size = 5, 
              label = T) + 
  guides(color = guide_legend(title = "Cluster"))
p2 <- DimPlot(object = d, reduction = "umap",
              pt.size = 0.1, group.by = 'condition') + 
  guides(color = guide_legend(title = "Treatment"))

p1 <- p1 + labs(title = paste0("Resolution ",res))

fig <- p1+p2

ggsave("DimMap_clustree_rpca.png",
       plot = fig,
       dpi = 300, height = 12,
       width = 28, units = 'cm')

#Save file####
saveRDS(file = "ICI_MS.rds", object = d)

s <- subset(d, cells = sample(colnames(d), 5000))
s <- RunTSNE(s, dims = 1:50)
saveRDS(file = "/data/pearceed/repository/Database/scRNA/ICI_MS", 
        object = s, version = 2)

#Find markers####
markers <- FindAllMarkers(object = d, only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25,
                          return.thresh = 0.01)

#Filter ribosomal genes
markers <- markers[-grep('Rpl',markers$gene),]
markers <- markers[-grep('Rps',markers$gene),]

write.table(x = markers, quote = F, sep = "\t",
            file = "Markers.txt")

#will print just the cluster IDS instead of every cell name
top3 <- markers %>% group_by(cluster) %>% top_n(3, avg_logFC)

s <- subset(d, cells = sample(colnames(d), 1000))

p1 <- DoHeatmap(object = s, features = top3$gene) + NoLegend()

ggsave(filename = "Top3_markers.png",
       plot = p1, units = 'mm', dpi = 300, 
       width = 360, height = 360)

#Cell count####
out <- table(Idents(d), d$condition)
write.table(out,"Proportions_condtion.txt",
            sep="\t",quote=F,col.names=NA)

out <- t(t(out)/colSums(out)*100)
out <- data.frame(out)
colnames(out) <- c('Cluster','Condition','Percentage')

p1 <- DimPlot(object = d, reduction = "umap",
              pt.size = 0.01, split.by = 'condition',
              combine = F, ncol = 2)

p1 <- lapply(p1, function(x){
  x+theme(text = element_text(face = 'plain', size = 12),
          strip.text = element_text(face = 'plain', size = 12))
})

p2 <- ggplot() +
  geom_bar(aes(y = Percentage, x = Condition, fill = Cluster), 
           color ='black', data = out, stat = "identity") +
  scale_y_continuous(labels = dollar_format(suffix = "%", prefix = "")) +
  labs(x = "Condition", y = "Percentage") +
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank(),
        legend.position = 'none',
        text=element_text(family="Arial"),
        axis.text.y = element_text(colour="black", size = 10),
        axis.text.x = element_text(colour="black", size = 10, 
                                   angle = 30, vjust = 1, 
                                   hjust = 1),
        axis.title =element_text(colour="black", size = 12))

ggsave(filename = 'Distribution_condition.png', 
       plot = (p1[[1]]|p2)+plot_layout(widths = c(3,1)),
       device = 'png', width = 210, height = 290,
       units = 'mm', dpi = 300)

#Cell cycle genes####
s.genes <- readLines(con = '/home/sanin/Code/scRNA_seq/Cell_cycle_genes/s.genes.mouse.txt')
g2m.genes <- readLines(con = '/home/sanin/Code/scRNA_seq/Cell_cycle_genes/g2m.genes.mouse.txt')

d <- CellCycleScoring(d, s.features = s.genes, 
                      g2m.features = g2m.genes, 
                      set.ident = FALSE)

p2 <- VlnPlot(object = d, 
              features = c('S.Score','G2M.Score'), 
              pt.size = 0.01, 
              combine = F)

p2 <- lapply(p2, function(x){
  x + theme(text = element_text(face = 'plain', size = 12),
            plot.title = element_text(face = 'plain', size = 12), 
            axis.text.x = element_text(face = 'plain', size = 12, 
                                       angle = 0, hjust = 0.5),
            axis.title.x = element_blank()) +
    NoLegend()
})

ggsave(filename = 'Cell_cycle_scoring.png', 
       plot = wrap_plots(p2, nrow = 2),
       device = 'png', width = 210, height = 160,
       units = 'mm', dpi = 300)

#Macrophage clusters####
d <- AddModuleScore(object = d, name = 'Mac',
                    features = list('Mac' = c('H2-Ab1','Lyz2',
                                              'Csf1r','Adgre1',
                                              'Mertk','Cd164',
                                              'Cd68','Itgam')))

p1 <- DimPlot(object = d, reduction = "umap",
              pt.size = 0.1, label.size = 5,
              label = T) +
  NoLegend()
p2 <- FeaturePlot(d, features = 'Mac1',
                  min.cutoff = 'q20', pt.size = 0.1,
                  max.cutoff = 'q80') + NoLegend()

d <- AddModuleScore(object = d, name = 'NK',
                    features = list('NK' = c('Itgam','Cd27','Il2rb',
                                             'Il7ra','Klrb1','Itga2',
                                             'Klra1','Klrk1','Ncr1')))

p3 <- FeaturePlot(d, features = 'NK1',
                  min.cutoff = 'q20', pt.size = 0.1,
                  max.cutoff = 'q80') + NoLegend()

d <- AddModuleScore(object = d, name = 'Th',
                    features = list('Th' = c('Cd3e','Cd4','Cd8a','Il2ra')))

p4 <- FeaturePlot(d, features = 'Th1',
                  min.cutoff = 'q20', pt.size = 0.1,
                  max.cutoff = 'q80') + NoLegend()

p <- list(p1,p2,p3,p4)

p <- lapply(p, function(x){
  x + theme(text = element_text(face = 'plain', size = 12),
            plot.title = element_text(face = 'plain', size = 12), 
            axis.text = element_text(face = 'plain', size = 12, 
                                       angle = 0, hjust = 0.5)) +
    NoLegend() +
    xlab('UMAP 1') +
    ylab('UMAP 1')
})

ggsave(filename = 'Cell_Type_scoring.png', 
       plot = wrap_plots(p, nrow = 2),
       device = 'png', width = 210, height = 160,
       units = 'mm', dpi = 300)

#Subsetting based on cell populations####
tc <- list('Th' = c(21,9,24,8,20,12,26,2,25,0,10,7,17),
           'NK' = c(11,13,15,16,5,4),
           'Mac' = c(31,18,3,19,28,29,6,27,14))

p <- list()
ds <- list()
z = 1
for(i in tc){
  t <- subset(d, idents = i)
  t@meta.data <- t@meta.data[,1:4]
  t <-SplitObject(object = t, split.by = 'condition')
  n = 1
  max.f = 100000
  for(j in t){
    t[[n]] <- SCTransform(object = j, assay = "RNA")
    if(nrow(t[[n]]) < max.f){
      max.f <- nrow(t[[n]])
    }
    n = n+1
  }
  features <- SelectIntegrationFeatures(object.list = t, nfeatures = max.f)
  max.size = 350000*max.f #To correct the globals
  options(future.globals.maxSize= max.size)
  t <- PrepSCTIntegration(object.list = t, anchor.features = features)
  
  #Run PCA in each object with integration features
  t <- lapply(X = t, FUN = function(x) {
    x <- RunPCA(x, verbose = FALSE, npcs = 50, features = features)
  })
  
  #Find anchors
  anchors <- FindIntegrationAnchors(object.list = t, 
                                    normalization.method = "SCT",
                                    anchor.features = features,
                                    reduction = "rpca", 
                                    dims = 1:50)
  
  #integrate dataset
  t <- IntegrateData(anchorset = anchors, 
                       new.assay.name = 'integrated',
                       normalization.method = "SCT",
                       dims = 1:50)
  
  t <- RunPCA(t, verbose = FALSE)
  t <- RunUMAP(t, dims = 1:50)
  t <- FindNeighbors(object = t, dims = 1:50)
  res <- list(seq(0, 1.2, by=0.1))
  for(j in 1:length(res)) {
    t <- FindClusters(object = t, resolution = res[[j]])
  }
  colnames(t@meta.data) <- gsub('([A-Za-z]+)(_)(snn_res.)([0-9.]+)',
                                     "res_\\4",
                                     colnames(t@meta.data))
  clust <- t@meta.data[,grep('res_',colnames(t@meta.data))]
  p[[z]] <- clustree(clust, prefix = 'res_')
  ds[[names(tc)[z]]] <- t
  z = z+1
}

trees <- p
res <- list('Th' = 'res_0.6',
            'NK' = 'res_0.3',
            'Mac' = 'res_1')

s.genes <- readLines(con = '/home/sanin/Code/scRNA_seq/Cell_cycle_genes/s.genes.mouse.txt')
g2m.genes <- readLines(con = '/home/sanin/Code/scRNA_seq/Cell_cycle_genes/g2m.genes.mouse.txt')
goi <- list()
acc <- c("mmu00010","mmu00190")
for(i in acc){
  input <- keggGet(i)
  term <- input[[1]]$PATHWAY_MAP
  genes <- gsub('([A-Za-z0-9]+)(;)(.+)','\\1',input[[1]]$GENE)[c(F,T)]
  goi[[term]] <- genes
}

term <- 'Macrophage angiogenesis'
genes <- c('Fgf1','Ang2','Plg','Ace','Fgf2','Ccl2',
           'Cxcl5','Csf2','Csf3','Igf1',
           'Igf2','Il1a','Il1b','Il10','Il19',
           'Il20','Il6','Mmp2','Mmp9','Pdgf',
           'Tymp','Tsp1','Tsp2','Tgfa','Tgfb',
           'Tnf','Plau','Vegfa','Vegfb','Vegfc',
           'Vegfd')
goi[[term]] <- genes

term <- 'M2'
genes <- c("Retnla","Clec10a","Mgl2",'Mrc1','Cd163',
           'Il10','Il4ra','Lipa','Apoe')
goi[[term]] <- genes

for(i in names(ds)){
  Idents(object = ds[[i]]) <- ds[[i]][[res[[i]]]]
  saveRDS(ds[[i]], 
          file = paste0("Subsets/",i,"_subset.rds"))
  
  markers <- FindAllMarkers(object = ds[[i]], only.pos = TRUE,
                            min.pct = 0.25,
                            logfc.threshold = 0.25,
                            return.thresh = 0.01)
  
  #Filter ribosomal genes
  if(sum(str_detect(rownames(markers),'Rpl')) > 0){
    markers <- markers[-grep('Rpl',markers$gene),]
  }
  if(sum(str_detect(rownames(markers),'Rps')) > 0){
    markers <- markers[-grep('Rps',markers$gene),]
  }
  
  write.table(x = markers, quote = F, sep = "\t",
              file = paste0(i,"_markers.txt"))
  
  #will print just the cluster IDS instead of every cell name
  top8 <- markers %>% group_by(cluster) %>% top_n(8, avg_logFC)
  p1 <- DoHeatmap(object = ds[[i]], features = top8$gene) + NoLegend()
  
  ggsave(filename = paste0(i,"_Top8_markers.png"),
         plot = p1, units = 'mm', dpi = 300, 
         width = 360, height = 360)
  
  #Cell count####
  out <- table(Idents(ds[[i]]), ds[[i]]$condition)
  write.table(out,paste0(i,"_proportions_condtion.txt"),
              sep="\t",quote=F,col.names=NA)
  
  out <- t(t(out)/colSums(out)*100)
  out <- data.frame(out)
  colnames(out) <- c('Cluster','Condition','Percentage')
  
  p1 <- DimPlot(object = ds[[i]], reduction = "umap",
                pt.size = 0.01, split.by = 'condition',
                combine = F, ncol = 2)
  
  p1 <- lapply(p1, function(x){
    x+theme(text = element_text(face = 'plain', size = 12),
            strip.text = element_text(face = 'plain', size = 12))
  })
  
  p2 <- ggplot() +
    geom_bar(aes(y = Percentage, x = Condition, fill = Cluster), 
             color ='black', data = out, stat = "identity") +
    scale_y_continuous(labels = dollar_format(suffix = "%", prefix = "")) +
    labs(x = "Condition", y = "Percentage") +
    theme(axis.line = element_line(size=0.5, colour = "black"),
          panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_blank(),
          legend.position = 'none',
          text=element_text(family="Arial"),
          axis.text.y = element_text(colour="black", size = 10),
          axis.text.x = element_text(colour="black", size = 10, 
                                     angle = 30, vjust = 1, 
                                     hjust = 1),
          axis.title =element_text(colour="black", size = 12))
  
  ggsave(filename = paste0(i,'_distribution_condition.png'), 
         plot = (p1[[1]]|p2)+plot_layout(widths = c(3,1)),
         device = 'png', width = 210, height = 290,
         units = 'mm', dpi = 300)
  
  ds[[i]] <- CellCycleScoring(ds[[i]], s.features = s.genes, 
                              g2m.features = g2m.genes, 
                              set.ident = FALSE)
  
  ds[[i]] <- AddModuleScore(object = ds[[i]],
                            features = goi,
                            assay = 'SCT',
                            name = 'Signatures')
}

for(i in names(ds)){
  colnames(ds[[i]]@meta.data)[
    grep('Signature',colnames(ds[[i]]@meta.data))] <- c('Glycolysis',
                                                        'OCPHOS',
                                                        'Angiogenesis',
                                                        'M2')
  saveRDS(ds[[i]],
          file = paste0("Subsets/",i,"_subset.rds"))
  ds[[i]] <- RunTSNE(ds[[i]], dims = 1:50)
  saveRDS(file = paste0("/data/pearceed/repository/Database/scRNA/ICI_MS_",i), 
          object = ds[[i]], version = 2)
}

#Graphs for paper####
ds <- c(readRDS(file = "ICI_MS.rds"),
        readRDS(file = "Subsets/Mac_subset.rds"),
        readRDS(file = "Subsets/NK_subset.rds"),
        readRDS(file = "Subsets/Th_subset.rds"))

for(i in ds){
  p <- list()
  p[[1]] <- DimPlot(object = i, reduction = "umap",
                pt.size = 0.01, label = T)
  
  for(j in unique(i$condition)){
    cells <- rownames(i@meta.data[i$condition == j,])
    p[[j]] <- DimPlot(object = i, reduction = "umap",
                       pt.size = 0.01, label = F, 
                       combine = FALSE,
                       cells = cells)[[1]] + 
      labs(title = paste(j))
  }
  
  p <- lapply(p, function(x){
    x+theme(text = element_text(face = 'plain', size = 12),
            strip.text = element_text(face = 'plain', size = 12),
            axis.text = element_text(face = 'plain', size = 12),
            axis.title = element_text(face = 'plain', size = 12),
            plot.title = element_text(face = 'plain', size = 12)) +
      labs(x = 'UMAP 1', y = "UMAP 2") + 
      NoLegend()
    })
  
  ggsave(paste0("paper/UMAPs_",ncol(i),".pdf"),
         plot = wrap_plots(p, nrow = 2), 
         height = 160, width = 240, units = 'mm',
         dpi = 300)
}

goi <- list()

term <- 'Macrophage angiogenesis'
genes <- c('Fgf1','Ang2','Plg','Ace','Fgf2','Ccl2',
           'Cxcl5','Csf2','Csf3','Igf1',
           'Igf2','Il1a','Il1b','Il10','Il19',
           'Il20','Il6','Mmp2','Mmp9','Pdgf',
           'Tymp','Tsp1','Tsp2','Tgfa','Tgfb',
           'Tnf','Plau','Vegfa','Vegfb','Vegfc',
           'Vegfd')
goi[[term]] <- genes

term <- 'M2'
genes <- c("Retnla","Clec10a","Mgl2",'Mrc1','Cd163',
           'Il10','Il4ra','Lipa','Apoe','Arg1')
goi[[term]] <- genes

genes <- c("Il1b", "Il6", "Il10", "Il12a", "Il12b", "Tgfb1",
           "Tnf", "Nos2", "Cd33", "Irf4", "Irf8", "Junb", 
           "Zeb1", "Bhlhe41", "Nfkbia", "Nfkbiz", "Nfkb1", 
           "Cd274","Pdcd1lg1", "Cd80", "Cd86", "H2-Ab1", 
           "H2-Q10", "Lyz2", "Csf1r", "Csf2r", "Csf3r", 
           "Adgre1", "Mertk", "Cd164", "Cd68", "Itgam", 
           "Arg2", "Csf1", "Csf2", "Csf3", "Ccl1", "Ccl2", 
           "Ccl3", "Ccl5", "Cxcl2", "Cxcl3", "Cxcl9", "Cxcl10",
           "Cxcr2", "Ccr1", "Ccr2", "Ccr3", "Ccr5", "Slc2a3",
           "Slc7a11", "Chil3", "Mki67", "Vegfa", "Mmp9", "S100a8",
           "S100a9", "Gzmm", "Tlr1", "Tlr2", "Tlr3", "Tlr4", "Tlr5",
           "Tlr6","Tlr7", "Tlr8", "Tlr9", "Il1r1", "Il1r2")
term <- 'Random'
goi[[term]] <- genes

markers <- read.table('Mac_markers.txt', sep = '\t') %>%
  dplyr::filter(p_val_adj < 0.01, avg_logFC > 0.3) %>%
  distinct(gene) %>%
  dplyr::select('gene')

p <- list()
for(i in names(goi)){
  genes <- goi[[i]][goi[[i]] %in% markers$gene]
  print(length(genes))
  x <- DotPlot(object = ds[[2]], features = genes)$data %>%
    dplyr::select('avg.exp', 'id', 'gene' = 'features.plot') %>%
    reshape2::acast(value.var = 'avg.exp', formula = gene ~ id)
  write.table(x,paste0('paper/dotplot.data.',i,'.txt'), sep = '\t')
  genes <- rownames(x)[hclust(dist(x))$order]
  p[[i]] <- DotPlot(object = ds[[2]], features = genes, dot.scale = 4,
                    scale.min = 0.01, cols = c('#e52165','#0d1137'))
}

p <- lapply(p, function(x){
  x+theme(text = element_text(face = 'plain', size = 12),
          strip.text = element_text(face = 'plain', size = 12),
          axis.text.x = element_text(face = 'plain', size = 6, 
                                     angle = 90, vjust = 0.5,
                                     hjust = 1),
          axis.text.y = element_text(face = 'italic', size = 6),
          axis.title = element_blank(),
          plot.title = element_text(face = 'plain', size = 12)) +
    coord_flip()
})

ggsave("paper/Dot_plots.pdf",
       plot = wrap_plots(p, nrow = 3, 
                         guides = 'collect',
                         heights = c(12,9,62)), 
       height = 280, width = 120, units = 'mm',
       dpi = 300)

genes <- c("Siglece", "Siglecf", "Siglecg")
p <- VlnPlot(object = ds[[2]], features = genes, 
             pt.size = 0.01, combine = F)

p[[3]] <- p[[3]]+ 
  theme(axis.text.x = element_text(face = 'plain', size = 8, 
                                   angle = 0, hjust = 0.5),
        axis.text.y = element_text(face = 'plain', size = 8),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = 'italic', size = 12),
        plot.title = element_blank()) +
  NoLegend() +
  ylab(genes[[3]])

p[1:2] <- lapply(p[1:2], function(x){
  x+ 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(face = 'plain', size = 8),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = 'italic', size = 12),
          plot.title = element_blank()) +
    NoLegend()
})

p[[1]] <- p[[1]]+ylab(genes[[1]])
p[[2]] <- p[[2]]+ylab(genes[[2]])

ggsave("paper/Vlnplots.pdf",
       plot = wrap_plots(p, nrow = 3), 
       height = 90, width = 90, units = 'mm',
       dpi = 300)

p1 <- FeaturePlot(object = ds[[1]], features = 'Siglece', 
                  cols = colorRampPalette(c('#191970'	,'#008B8B',	'#00FF00'))(255), 
                  pt.size = 0.01, min.cutoff = 'q40', 
                  max.cutoff = 'q80', combine = F, order = T)[[1]] +
  theme(text = element_text(face = 'plain', size = 12),
        strip.text = element_text(face = 'plain', size = 12),
        axis.text = element_text(face = 'plain', size = 12),
        axis.title = element_text(face = 'plain', size = 12),
        plot.title = element_text(face = 'italic', size = 12)) +
  labs(x = 'UMAP 1', y = "UMAP 2") + 
  NoLegend()

p2 <- FeaturePlot(object = ds[[2]], features = 'Siglece', 
                  cols = colorRampPalette(c('#191970'	,'#008B8B',	'#00FF00'))(255), 
                  pt.size = 0.01, min.cutoff = 'q40', 
                  max.cutoff = 'q80', combine = F, order = T)[[1]] +
  theme(text = element_text(face = 'plain', size = 12),
        strip.text = element_text(face = 'plain', size = 12),
        axis.text = element_text(face = 'plain', size = 12),
        axis.title = element_text(face = 'plain', size = 12),
        plot.title = element_text(face = 'italic', size = 12)) +
  labs(x = 'UMAP 1', y = "UMAP 2") + 
  NoLegend()

ggsave("paper/Siglece.pdf",
       plot = p1+p2, 
       height = 80, width = 160, units = 'mm',
       dpi = 300)  

#smaller dots
p <- list()
for(i in names(goi)){
  genes <- goi[[i]][goi[[i]] %in% markers$gene]
  print(length(genes))
  x <- DotPlot(object = ds[[2]], features = genes)$data %>%
    dplyr::select('avg.exp', 'id', 'gene' = 'features.plot') %>%
    reshape2::acast(value.var = 'avg.exp', formula = gene ~ id)
  #write.table(x,paste0('paper/dotplot.data.',i,'.txt'), sep = '\t')
  genes <- rownames(x)[hclust(dist(x))$order]
  p[[i]] <- DotPlot(object = ds[[2]], features = genes, dot.scale = 3,
                    scale.min = 0.01, cols = c('#e52165','#0d1137'))
}

p <- lapply(p, function(x){
  x+theme(text = element_text(face = 'plain', size = 12),
          strip.text = element_text(face = 'plain', size = 12),
          axis.text.x = element_text(face = 'plain', size = 6, 
                                     angle = 90, vjust = 0.5,
                                     hjust = 1),
          axis.text.y = element_text(face = 'italic', size = 6),
          axis.title = element_blank(),
          plot.title = element_text(face = 'plain', size = 12)) +
    coord_flip()
})

ggsave("paper/Dot_plots_small.pdf",
       plot = wrap_plots(p, nrow = 3, 
                         guides = 'collect',
                         heights = c(12,9,62)), 
       height = 230, width = 120, units = 'mm',
       dpi = 300)

#Reviewer info####
d <- readRDS(file = 'ICI_MS.rds')
df <- d@meta.data %>%
  group_by(condition) %>%
  dplyr::summarise(
    mean_gene = mean(nFeature_RNA),
    mean_molec = mean(nCount_RNA),
    count = n()
  )

markers <- read.table(file = 'Mac_markers.txt', 
                      sep = '\t', header = 1, row.names = 1)

genes <- c("Il1b", "Il6", "Il10", "Il12a", "Il12b", "Tgfb1",
           "Tnf", "Nos2", "Cd33", "Irf4", "Irf8", "Junb", 
           "Zeb1", "Bhlhe41", "Nfkbia", "Nfkbiz", "Nfkb1", 
           "Cd274","Pdcd1lg1", "Cd80", "Cd86", "H2-Ab1", 
           "H2-Q10", "Lyz2", "Csf1r", "Csf2r", "Csf3r", 
           "Adgre1", "Mertk", "Cd164", "Cd68", "Itgam", 
           "Arg2", "Csf1", "Csf2", "Csf3", "Ccl1", "Ccl2", 
           "Ccl3", "Ccl5", "Cxcl2", "Cxcl3", "Cxcl9", "Cxcl10",
           "Cxcr2", "Ccr1", "Ccr2", "Ccr3", "Ccr5", "Slc2a3",
           "Slc7a11", "Chil3", "Mki67", "Vegfa", "Mmp9", "S100a8",
           "S100a9", "Gzmm", "Tlr1", "Tlr2", "Tlr3", "Tlr4", "Tlr5",
           "Tlr6","Tlr7", "Tlr8", "Tlr9", "Il1r1", "Il1r2")

mark <- markers %>%
  dplyr::filter(cluster == 2| cluster == 13,
                p_val_adj <= 0.05,
                pct.1 >= 0.25,
                avg_logFC >= 0.25,
                gene %in% genes)

mark <- unique(mark$gene)

dm <- readRDS(file = 'Subsets/Mac_subset.rds')
table(Idents(dm),dm$condition) %>% 
  data.frame() %>%
  dplyr::filter(
    Var1 == 2 | Var1 == 13
  )


#Done####
