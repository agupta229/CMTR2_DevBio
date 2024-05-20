#######################################################################
# Final cleaned up code for Yermavolich 2024 CMTR2 publication
# Authors: Akansha Gupta (akansha@broadinstitute.org) and Bruno Giotti (brunogiotti@gmail.com)
# Contact: Matthew Meyerson, matthew_meyerson@dfci.harvard.edu
# Last updated: May 20, 2024

# License: GNU GPL3, Copyright (C) 2024 Dana-Farber Cancer Institute
#######################################################################

library (Seurat)
library (ggplot2)
library (ggpubr)
library (patchwork)
library (tidyr)
library (plyr)
library (tidyverse)
library (rstatix)
library (circlize)
library (paletteer)
library (readxl)
library (fgsea)
library (ComplexHeatmap)
library (ggrepel)
library (cowplot)
library (reshape2)
library (S4Vectors)
library (SingleCellExperiment)
library (pheatmap)
library (apeglm)
library (RColorBrewer)
library (data.table)
setwd("/meyersonlab/akansha/Alena")

# Global variables
res = c(0.2, 0.8, 2, 3, 5) # denovo cluster resolutions
# Set variables for UMAP
reductionSave = 'pca'
reductionGraphKnn = 'RNA_knn'
reductionGraphSnn = 'RNA_snn' 
reductionName = 'umap'
projdir <- "/meyersonlab/akansha/Alena/results/plots/" # Location to store plots / results
date <- "04.17.2024" # To store results
# dir.create("/meyersonlab/akansha/Alena/results/plots", recursive = T) # Only once

# NOTE: PART 1 only needs to be completed ONCE, then skip directly to lower parts as needed

#######################################################################
# PART 1: snRNA-seq Analysis: QC, Normalization, etc
# Only need to run through these steps once
#######################################################################

# (1a) SETUP
#############################################
# Specify folder, metadata project name and other variables
meta = read.csv("/meyersonlab/akansha/Alena/Shared/metadata.csv", row.names=NULL)
proj_name = 'Alena_snRNA'

samples_path = sapply(meta$sampleID, function(x) 
  paste0("/meyersonlab/akansha/Alena/Shared/Data/", x, "_raw_feature_bc_matrix"))
meta$sample_path = samples_path

samples_path = samples_path
nFeat = 400 # Number of features per cells. default 400

# Import count matrices
scrna_mats = sapply (meta$sample_path, function (dir) 
{
  if (grepl('.h5',dir)) 
  {
    Seurat::Read10X_h5 (dir, use.names = TRUE, unique.features = TRUE)    
  } else {
    Seurat::Read10X (data.dir = dir)
  }
})

# Create seurat object
srt = sapply (seq_along (scrna_mats), function(i) CreateSeuratObject (counts = scrna_mats[[i]], min.cells = 0, min.features = nFeat, project = proj_name))  
names (srt) = meta$sampleID
#############################################

# (1b) BEGIN ANALYSIS AND DO QC
#############################################
meta = read.csv("/meyersonlab/akansha/Alena/Shared/metadata.csv", row.names=NULL)
samples_path = sapply(meta$sampleID, function(x) 
  paste0("/meyersonlab/akansha/Alena/Shared/Data/", x, "_raw_feature_bc_matrix"))
meta$sample_path = samples_path


# QC AND CELL SELECTION
# We first went to select cells for downstream analysis
# Common metrics used include: 
# - Number of unique genes detected in cells (low quality often have few genes; multiplets have high)
# - Total number of molecules detected within cell (correlates strongly with unique genes)
# - Percentage of reads that map to mitochondrial genome (# low quality cells have high p.mito)

# Filters we will use
nCounts = 800 # Number of UMI per cell. Default 800
pchM = 25 # Percent mitochondrial genes. Default 25 
variablefeatures = 'seurat'
nfeat = 2000 # number of variable genes to consider for dimentionality reduction
sigPCs = 15
ccRegress = FALSE # Regress cell cycle gene expression 
metaGroupNames = c('Genotype','sampleID','Litter.ID','Sex')
org = 'mouse'
nFeat = 400
proj_name = 'Alena_snRNA'

# Calculate mitochondrial QC metrics with PercentageFeatureSet()
# which calculates percent of counts originating from set of features
if (org == 'mouse') for (i in seq_along (srt)) srt[[i]]$percent.mt = PercentageFeatureSet (srt[[i]], pattern = "^mt-")

# Filter data
cell_idx = sapply (lapply (srt, function(x) x$nFeature_RNA > nFeat & x$percent.mt < pchM & x$nCount_RNA > nCounts), sum)
srt = lapply (srt[cell_idx > 0], function (x) x[,x$nFeature_RNA > nFeat
                                                & x$percent.mt < pchM & x$nCount_RNA > nCounts])
# End of QC
#############################################

# (1c) NORMALIZE DATA AND IDENTIFY HIGHLY VARIABLE FEATURES
#############################################
# After removing unwwanted cells, normalize data

# Employ global-scaling normalization method "LogNormalize" that normalizes feature expression
# measurements for each cell by total expression, multiplies by scale factor (10000 = default),
# and log transformst the result

# Normalized values stored in srt[["RNA]]$data

# Note: relies on assumption that each cell originally contains same number of RNA molecules

# Merge samples into one seurat object and map metadata
set.seed (1234)
message ('Merge samples in one seurat object and map metadata')
srtM = merge (srt[[1]], y = srt[2:length(srt)],
              add.cell.ids = names (srt), project = proj_name)        
srtM$sampleID = rep (meta$sampleID, sapply (srt, ncol)) # attach sampleIDs to merged object
srtM$orig.ident = rep (meta$sampleID, sapply (srt, ncol)) # attach sampleIDs to merged object
srt = srtM

# Add additional metadata to merged obj
if (exists('meta')) {
  if (ncol (meta) > 1) {
    for (x in colnames (meta[,colnames(meta) != 'sampleID', drop=F])) { 
      metaGroup = meta[,x]
      names(metaGroup) = meta$sampleID
      srt@meta.data[,x] = metaGroup[srt$sampleID]
    }
  }
} 

# Normalize data
srt = NormalizeData (object = srt, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features: exhibit high cell-to-cell variation in dataset
# By default, returns 2000 features per dataset
if (variablefeatures == 'seurat') srt = FindVariableFeatures (srt, selection.method = "vst", nfeat = nfeat)

#############################################

# (1d) SCALE DATA AND PERFORM LINEAR DIMENSIONAL REDUCTION
#############################################
# Apply linear transformation (scaling) that is standard prior to dimensional 
# reduction techniques like PCA

# Calculate cell cycle phase scores based on canonical markers
# Cell cycle scoring inputs -- marker sets should be anticorrelated in expression level
# Cells expressing neither are likely not cycling and in G1 phase
if (org == 'mouse') { 
  s.genes = readRDS ('/meyersonlab/akansha/Alena/Shared/cellcycle_mouse.s.genes.rds') # S phase markers
  g2m.genes = readRDS ('/meyersonlab/akansha/Alena/Shared/cellcycle_mouse.g2m.genes.rds') # G2/M phase markers
}

srt <- JoinLayers(srt)

srt = CellCycleScoring (object = srt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
srt$cc = srt$S.Score + srt$G2M.Score

# Apply linear transformation (scaling); no unwanted sources of variation -> no vars_to_regress
srt = ScaleData (srt, features = VariableFeatures (object=srt)) 

# Perform PCA on scaled data
srt = RunPCA (srt, features = VariableFeatures (object = srt), npcs = ifelse(ncol(srt) <= 30,ncol(srt)-1,30), 
              ndims.print = 1:5, nfeat.print = 5, verbose = FALSE)
#############################################

# (1e) CREATE QC PLOTS
#############################################
# PLOT 1: Plot distribution of nFeature, nCounts, p.mito (Supplemental Figure)

ccomp_df = as.data.frame (table(srt$sampleID))
# Replace sample ID labels with WT vs KO # labels
# Order basically same as in meta table --> manually updating
ccomp_df$simple_id <- c('KO1', 'WT1', 'WT2', 'KO2', 'KO3', 'KO4', 'WT3', 'KO5', 'WT4')
ccomp_df$simple_id <- factor(ccomp_df$simple_id, levels = c('WT1', 'WT2', 'WT3', 'WT4', 'KO1', 'KO2',
                                                            'KO3', 'KO4', 'KO5'))

cc_p2 = ggplot (ccomp_df, aes (x= simple_id, y= Freq)) +
  geom_bar (position="stack", stat="identity") +
  labs (title = paste('Total cells',ncol(srt)), y = "Frequency") + 
  scale_y_continuous(limits = c(0, 5000), expand = c(0, 0)) +
  theme_classic() + 
  theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1, size = 14, color = "black"), 
         axis.title.x = element_blank(),
         axis.text.y = element_text(size = 14, color = "black"),
         axis.title.y = element_text(size = 16, color = "black"),
         plot.title = element_text(size = 16, color = "black", vjust = 0.5, hjust = 0.5)) 
show(cc_p2)

srt$nFeature_RNAL = log10 (srt$nFeature_RNA)
srt$nCount_RNAL = log10 (srt$nCount_RNA)

vln_p = VlnPlot (srt, features = c("nFeature_RNAL", "nCount_RNAL", "percent.mt"), combine=T,group.by = 'sampleID',pt.size = 0, ncol = 3)

vln_p1 <- vln_p[[1]]$data
vln_p1 <- left_join(vln_p1, ccomp_df, by = c("ident" = "Var1"))
vln_p1 <- ggplot(vln_p1, aes(x = simple_id, y = nFeature_RNAL, fill = simple_id)) + 
  geom_violin( adjust =1,trim=TRUE, scale = "width") +
  labs (title = "Number of Features", y = expression("log"[10]*"(# of Features)")) + 
  scale_y_continuous(limits = c(2, 5), expand = c(0, 0)) +
  theme_classic() + 
  theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1, size = 14, color = "black"), 
         axis.title.x = element_blank(),
         axis.text.y = element_text(size = 14, color = "black"),
         axis.title.y = element_text(size = 16, color = "black"),
         plot.title = element_text(size = 16, color = "black", vjust = 0.5, hjust=0.5),
         legend.position = "none")
show(vln_p1)

vln_p2 <- vln_p[[2]]$data
vln_p2 <- left_join(vln_p2, ccomp_df, by = c("ident" = "Var1"))
vln_p2 <- ggplot(vln_p2, aes(x = simple_id, y = nCount_RNAL, fill = simple_id)) + 
  geom_violin( adjust =1,trim=TRUE, scale = "width") +
  labs (title = "Counts", y = expression("log"[10]*"(nCounts)")) + 
  scale_y_continuous(limits = c(2, 5), expand = c(0, 0)) +
  theme_classic() + 
  theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1, size = 14, color = "black"), 
         axis.title.x = element_blank(),
         axis.text.y = element_text(size = 14, color = "black"),
         axis.title.y = element_text(size = 16, color = "black"),
         plot.title = element_text(size = 16, color = "black", vjust = 0.5, hjust=0.5),
         legend.position = "none")
show(vln_p2)

vln_p3 <- vln_p[[3]]$data
vln_p3 <- left_join(vln_p3, ccomp_df, by = c("ident" = "Var1"))
vln_p3 <- ggplot(vln_p3, aes(x = simple_id, y = percent.mt, fill = simple_id)) + 
  geom_violin( adjust =1,trim=TRUE, scale = "width") +
  labs (title = "Mitochondrial Percentages", y = "Percentage") + 
  scale_y_continuous(limits = c(0, 25), expand = c(0, 0)) +
  theme_classic() + 
  theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1, size = 14, color = "black"), 
         axis.title.x = element_blank(),
         axis.text.y = element_text(size = 14, color = "black"),
         axis.title.y = element_text(size = 16, color = "black"),
         plot.title = element_text(size = 16, color = "black", vjust = 0.5, hjust=0.5),
         legend.position = "none")
show(vln_p3)

ggsave(paste0("/meyersonlab/akansha/Alena/results/plots/QC_nFeat_nCount_m.percent_vlnPlot_", date, ".pdf"), (cc_p2 | vln_p1 | vln_p2 | vln_p3) + plot_layout (widths=c(2,2,2,2)), device='pdf', width = 15, height = 5)

# PLOTS 2 AND 3 (broken down below): Plot distribution of nFeature, nCounts, p.mito via UMAP
srt = RunUMAP (object = srt, reduction = reductionSave, dims = 1:sigPCs)

# PLOT 2: UMAP of non corrected and harmony corrected clusters
humap_p = lapply (metaGroupNames, function(y) DimPlot (object = srt, reduction = reductionName, pt.size = .01, group.by = y) + theme_classic())
png (paste0(projdir, paste(metaGroupNames, collapse='_'),'_umap.png'), width = 2500, height = 1200, pointsize=10, res = 300, type="cairo")
print (wrap_plots (humap_p), ncol=3)
dev.off()

# PLOT 3: Color UMAP by number of detected genes and color from 5% to 95% of counts
qc_metrics = c('nCount_RNA', 'nFeature_RNA', 'percent.mt')
umap_df = data.frame (srt[[reductionName]]@cell.embeddings, nCount_RNA = log10(srt$nCount_RNA+1), nFeature_RNA = log10(srt$nFeature_RNA+1), percent.mt = srt$percent.mt)
umap_p = lapply (qc_metrics, function(x) ggplot(data = umap_df) + 
                   geom_point (mapping = aes_string (x = colnames(umap_df)[1], y=colnames(umap_df)[2], color = x), size = .01) + 
                   scale_colour_gradientn (colours = rainbow(7)) +
                   theme_classic() + 
                   theme(
                     plot.background = element_blank()
                   ))
png (paste0(projdir,'QC_umap.png'), width = 3500, height = 1200, pointsize=10, res = 300, type="cairo")
print (wrap_plots (umap_p), ncol=3)
dev.off()
#############################################

# (1f) CLUSTER CELLS
#############################################
# Graph based clustering approach: embed cells in graph structure and partition graph into 
# highly interconnected 'quasi-cliques' or 'communities'

# FindNeighbors() constructrs KNN graph based on euclidean distance in PCA space
# and refines edge weights between cells based on overlab in local neighborhoods
# Run denovo clustering on non-adjusted reductions
srt = FindNeighbors (object = srt, reduction = reductionSave, dims = 1:sigPCs, k.param = 30,
                     verbose = TRUE, force.recalc = T, graph.name=c(reductionGraphKnn,reductionGraphSnn))

# FindClusters() applies modularity optimization to iteratively group cells together 
for (i in seq_along(res)) srt = FindClusters (srt, resolution = res[i], verbose = T, n.start = 100, graph.name=reductionGraphSnn)
#############################################

# (1g) VISUALIZE CLUSTERS -- SETUP
#############################################
# Visualize clustered cells using UMAP
clust_p1 = list()
for (i in seq_along(res)) clust_p1[[i]] = DimPlot (srt, pt.size = .01, label = T, group.by= paste0(reductionGraphSnn,'_res.',res[i]), reduction = reductionName) + NoLegend()
png (paste0(projdir,'denovo_clusters_',reductionSave,'_umaps.png'), 3000, 3000, pointsize=10, res = 300, type="cairo")
print (wrap_plots (clust_p1, ncol=2))
dev.off() 

# Use JShendure mouse dev data to annotate cell types
# Set additional column metadata
srt$Genotype2 = ifelse (srt$Genotype == 'Cmtr2 WT','WT','KO')
ccomp_df = as.data.frame (table(srt$sampleID))
ccomp_df$sampleID2 = srt$Genotype2 [match (ccomp_df$Var1, srt$sampleID)]
ccomp_df = ccomp_df[order (ccomp_df$sampleID2),]
ccomp_df$sampleID2 = paste0(ccomp_df$sampleID2, c(1,2,3,4,5,1,2,3,4))
srt$sampleID2 = ccomp_df$sampleID2[factor(srt$sampleID, levels = unique(ccomp_df$Var1))]

# Color palette
genotype_pal = c(WT = '#BDBDBE', KO = '#733bbe')
celltype_pal = DiscretePalette(25, palette = "polychrome", shuffle = FALSE)

# Annotate cell clusters using label transfer from JShendure data and Marioni after subsampling
srt_js = readRDS ('/meyersonlab/akansha/Alena/Shared/srt_JShendure_mouse_dev_E9.5.rds')
srt_js$Main_cell_type[srt_js$Main_cell_type == 'Chondroctye progenitors'] = 'Chondrocyte progenitors'
srt_refs = list (srt_js)
names (srt_refs) = c('JShendure')
metaGroupNames = c('Main_cell_type')
subsample = Inf
set.seed (123)
small_ann = 10
for (i in seq_along (srt_refs)) {
  ref = srt_refs[[i]]
  metaGroup = ref@meta.data[,metaGroupNames[i]]
  metaGroup[is.na(metaGroup)] = 'not_assigned' # rename NA  
  sub_cells = lapply (unique(metaGroup), function(x) {
    sub_celltype = colnames(ref)[metaGroup == x]
    if (length(sub_celltype) < subsample) {
      sub_celltype } else {
        sub_celltype = sample (sub_celltype, subsample)
      }})
  sub_cells = unlist (sub_cells)
  ref = ref[,sub_cells[sub_cells != 'not_assigned']] # remove NA
  
  metaGroup = ref@meta.data[,metaGroupNames[i]]
  query = srt
  DefaultAssay (query) = 'RNA'
  anchors = FindTransferAnchors (reference = ref, query = query, dims = 1:30)
  predictions = TransferData (anchorset = anchors, refdata = metaGroup,
                              dims = 1:30)
  small_ann2 = table (predictions$predicted.id) < small_ann
  predictions$predicted.id[predictions$predicted.id %in% names(which(small_ann2))] = 'Others'
  write.csv (predictions, paste0('/meyersonlab/akansha/Alena/results/prediction_scores_',names (srt_refs)[i],'_subsampled_',subsample,'.csv'))
  srt@meta.data [,paste0('predicted_',names (srt_refs)[i], '_',subsample)] = predictions$predicted.id
}

saveRDS(srt, file = "/meyersonlab/akansha/Alena/results/labeled_seurat_object.rds") # This is a very large object -> may take time

#############################################

#######################################################################
# PART 2: snRNA-seq Analysis: UMAP Cluster Labeling, Cell Composition, etc
#######################################################################

# (2a) VISUALIZE CLUSTERS -- LOADING
# BEGIN HERE WHEN YOU COMPLETE RUNNING 
# THROUGH 'VISUALIZE CLUSTERS -- SETUP'
#############################################
# Load in Seurat object (avoids running through all of above every time)
srt <- readRDS("/meyersonlab/akansha/Alena/results/labeled_seurat_object.rds")

# Color palette
genotype_pal = DiscretePalette(2, palette = "glasbey", shuffle = FALSE)
celltype_pal = DiscretePalette(25, palette = "polychrome", shuffle = FALSE)

# Annotate cell clusters using label transfer from JShendure data and Marioni after subsampling
srt_js = readRDS ('/meyersonlab/akansha/Alena/Shared/srt_JShendure_mouse_dev_E9.5.rds')
srt_js$Main_cell_type[srt_js$Main_cell_type == 'Chondroctye progenitors'] = 'Chondrocyte progenitors'
srt_refs = list (srt_js)
names (srt_refs) = c('JShendure')
metaGroupNames = c('Main_cell_type')
subsample = Inf
set.seed (123)
small_ann = 10

# Generate barplots or boxplots of meta groups proportions specified
# meta_groups: vector or meta group names:
# 1) first meta_group is the group on which proportions are calculated
# 2) second meta_group split first meta_group on x axes  
# 3) third meta_group will group barplots separately
# if splits include only one value runs barplot instead
cellComp = function (
    seurat_obj = NULL, 
    metaGroups = NULL, # vector of at least 3 metaGroups e.g. c('orig.ident','celltypes','celltypes'),
    plot_as = 'box', # box or bar 
    pal = NULL,
    prop = TRUE,
    ptable_factor = 1, # specify which column of the data.frame or seurat object metadata should be used to compute proportions
    facet_ncol = 20,
    facet_scales = 'free',
    subset_prop = NULL, # subset prop table by any group in any column
    removeNA = TRUE,
    returnDF = FALSE
) {
  if (is.data.frame (seurat_obj)) {
    meta_groups_df = seurat_obj[,metaGroups]  
  } else {
    meta_groups_df = seurat_obj@meta.data[,metaGroups]
  }
  # Refactor to remove 0 groups
  if(is.null(pal)) pal = rainbow (length(unique(meta_groups_df[,2])))
  if (prop) {
    ccomp_df = as.data.frame (prop.table (table (meta_groups_df),ptable_factor))
    ccomp_df = na.omit (ccomp_df) # this is to remove NaN somehow produced from the line above 
  } else {
    ccomp_df = as.data.frame (table (meta_groups_df)) 
  }
  
  if(removeNA) ccomp_df = ccomp_df[ccomp_df$Freq != 0, ] # remove 0s from proportions
  if (!is.null (subset_prop)) {
    subset_col = unlist(sapply (seq(ncol(ccomp_df)), function(x) if(any(ccomp_df[,x] %in% subset_prop)) colnames(ccomp_df)[x]))
    ccomp_df = ccomp_df[ccomp_df[,subset_col] %in% subset_prop,]
  }
  if (plot_as == 'box') {
    p = ggplot (ccomp_df, aes_string (x= metaGroups[2], y= 'Freq')) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position = "bottom") +
      scale_fill_manual (breaks = c("WT", "KO"), values= pal) + 
      xlab (metaGroups[2]) + 
      ylab (ifelse (prop, 'proportion','counts'))
    if (length(metaGroups > 2)) p = p + geom_boxplot(aes_string (fill= metaGroups[3]), outlier.size=.2, alpha = 0.7, lwd=.2) 
    else p = p + geom_boxplot(aes_string (fill= metaGroups[2]), outlier.size=.2, alpha = 0.7, lwd=.2)   
    if (length(metaGroups) > 3) p = p + facet_wrap (as.formula(paste("~", metaGroups[4])), scales=facet_scales, ncol=facet_ncol)
  } 
  if (plot_as == 'bar') {
    p = ggplot (ccomp_df, aes_string (x= metaGroups[1], y= 'Freq')) +
      geom_bar(position="stack", stat="identity", aes_string(fill= metaGroups[2])) +
      theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      scale_fill_manual (values= pal) + xlab (metaGroups[2]) + ylab (ifelse (prop, 'proportion','counts'))
    if (length(metaGroups) == 3) p = p + facet_wrap (as.formula(paste("~", metaGroups[3])), scales=facet_scales, ncol=facet_ncol)
  }
  if (returnDF) return(ccomp_df) else 
    return (p)
}
#############################################

# (2b) VISUALIZE CLUSTERS -- ANALYSIS
#############################################
# Label clusters on UMAP (any aesthetic updates after this made in Illustrator)
for (i in seq_along (srt_refs)) {
  predictions = read.csv (paste0 ('/meyersonlab/akansha/Alena/results/prediction_scores_',names (srt_refs)[i],'_subsampled_',subsample,'.csv'))
  umap_df = data.frame (srt[[reductionName]]@cell.embeddings, prediction = srt@meta.data[,paste0('predicted_',names (srt_refs)[i], '_',subsample)])
  direction_col = umap_df %>% group_by(prediction) %>% summarise (mean (umap_1))
  umap_df$prediction = factor (umap_df$prediction, levels = direction_col$prediction[order(unlist(-direction_col[,2]))])
  labeltransfer_p1 = ggplot(data = umap_df) + 
    geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], color = 'prediction'), size = .1) + 
    ggtitle (names(srt_refs)[1]) + 
    theme_void() + 
    scale_color_manual (values = setNames(celltype_pal, levels(umap_df$prediction))) +
    theme(legend.key.size = unit(0.2, "cm")) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme(legend.position="bottom")
  ggsave(paste0("/meyersonlab/akansha/Alena/results/plots/labeltransfer_p1_", date, ".pdf"), labeltransfer_p1, device='pdf', width = 12, height = 8)
  
  umap_df = data.frame (srt[[reductionName]]@cell.embeddings, Genotype = srt@meta.data[,'Genotype2'])
  labeltransfer_p2 = ggplot(data = umap_df) + 
    geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], 
                                      fill = 'Genotype'), size = 1.25, shape = 21, stroke = NA) + 
    ggtitle (names(srt_refs)[1]) + 
    theme_classic() +
    theme(legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12)) +
    scale_fill_manual (breaks = c("WT", "KO"), values= genotype_pal)
  ggsave(paste0("/meyersonlab/akansha/Alena/results/plots/labeltransfer_p2_", date, ".pdf"), labeltransfer_p2, device='pdf', width = 13, height = 8.5)
  
  pred_mat = data.frame (predicted.id = predictions$predicted.id, score = predictions$prediction.score.max)
  predictions_p = ggplot (pred_mat, aes (x= predicted.id, y = score)) +
    geom_boxplot () +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  png (paste0(projdir,'label_tranfer_',names(srt_refs)[i],'_subsampled_',subsample,'UMAP_', date, '.png'), width=7500, height=3000,res=500)
  print (wrap_plots (labeltransfer_p2))
  dev.off()
}

# Plot both label transfers UMAPs: JShendure only
metaGroupNames_Jonly = c('predicted_JShendure_Inf')
labeltransfer_p1_Jonly = lapply (metaGroupNames_Jonly, function(x) DimPlot (srt, pt.size = .05, 
                                                                            label = T, group.by= x, repel = T,
                                                                            reduction = reductionName) 
                                 + ggplot2::theme(legend.position = "bottom") +
                                   ggplot2::scale_color_manual (values = setNames(celltype_pal, levels(umap_df$prediction))))
ggsave(paste0("/meyersonlab/akansha/Alena/results/plots/cell_annotations_p1_", date, ".pdf"), wrap_plots(labeltransfer_p1_Jonly, ncol = 1), device='pdf', width = 13, height = 8.5)

#############################################

# (2c) CELL COMPOSITION ANALYSIS
#############################################
metaGroupName = 'predicted_JShendure_Inf' 
metaGroupName2 = 'Genotype2'
cc_df = srt@meta.data[srt@meta.data[,metaGroupName] != 'Others',]
cc_df[,metaGroupName] = factor (cc_df[,metaGroupName], levels =names(table (cc_df[,metaGroupName])[order (-table (cc_df[,metaGroupName]))]))
cc_df$Genotype2 <- fct_rev(cc_df$Genotype2) # Added 02/20/2024 to reverse WT and KO order
cc_box1 = cellComp (
  seurat_obj = cc_df, 
  metaGroups = c('sampleID', metaGroupName, 'Genotype2'),
  plot_as = 'box',
  ptable_factor = c(1,3),
  pal = genotype_pal,
) + labs(y = "Proportion") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"), legend.title = element_blank())

ggsave(paste0("/meyersonlab/akansha/Alena/results/plots/cell_composition_JShendure_", date, ".pdf"), cc_box1, device='pdf', width = 8, height = 5)

# use speckle package to test for cell composition differences
library (speckle)
library (limma)

cell_info <- speckle_example_data()
head(cell_info)

metaGroupNames = c('sampleID','predicted_JShendure_Inf','Genotype2')

ccomp_df = srt@meta.data [,metaGroupNames]
cc_res = propeller (clusters = ccomp_df[,metaGroupNames[2]], sample = ccomp_df[,metaGroupNames[1]], group = ccomp_df[,metaGroupNames[3]])
write.csv (cc_res, 'cell_composition_analysis.csv')
#############################################

# (2d) CELL CYCLE COMPOSITION ANALYSIS
#############################################
# Box plot for cell cycle score between WT and KO 
metaGroupName1 = 'predicted_JShendure_Inf'
metaGroupName2 = 'sampleID2'
metaGroupName3 = 'Genotype2'

# Examine distributions to set manual threshold for analysis
pdf (paste0(projdir, 'cycling_distribution.pdf'))
hist (srt$S.Score)
hist (srt$G2M.Score)
dev.off()
srt$S.Score2 = ifelse (srt$S.Score > .2, 'cycling','G1')
srt$G2M.Score2 = ifelse (srt$G2M.Score > .2, 'cycling','G1')
srt$Phase2 = ifelse (srt$S.Score2 != 'G1' & srt$G2M.Score2 != 'G1', 'cycling','G1')
metaGroupNames = c(metaGroupName2, 'Phase2', metaGroupName3, metaGroupName1)
# Compare WT and KO cycling proportions across cell types 
ccc_box1 = cellComp (
  seurat_obj = srt, 
  metaGroups = metaGroupNames,
  plot_as = 'box',
  prop = TRUE,
  ptable_factor = c(1,4),
  subset_prop = 'cycling',
  facet_ncol = 6,
  facet_scales = 'free'
) + theme_classic()
show(ccc_box1)

# Look at relative cycling vs G1 proportions in each sample based on cell type
metaGroupNames = c(metaGroupName2, 'Phase2', metaGroupName1)
ccc_bar1 = cellComp (
  seurat_obj = srt, 
  metaGroups = metaGroupNames,
  plot_as = 'bar',
  pal = c(cycling = 'blue', G1= 'grey'),
  prop = TRUE,
  ptable_factor = c(1,3),
  facet_ncol = 6,
  facet_scales = 'free'
) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+ NoLegend() + ggtitle (i)
show(ccc_bar1)

# Look at cycling WT vs KO across all cell types --> figure 4C 
metaGroupNames = c(metaGroupName2, 'Phase2', metaGroupName3)
metaGroups = metaGroupNames
meta_groups_df = srt@meta.data[,metaGroupNames]
pal = genotype_pal
ccomp_df = as.data.frame (prop.table (table (meta_groups_df),1))
ccomp_df = na.omit (ccomp_df) # this is to remove NaN somehow produced from the line above 
ccomp_df = ccomp_df[ccomp_df$Freq != 0, ] # remove 0s from proportions
ccomp_df = ccomp_df[ccomp_df$Phase2 == 'cycling',]
ccomp_df$Phase2 = as.factor (as.character(ccomp_df$Phase2))
ccomp_df$Genotype2 <- factor(ccomp_df$Genotype2, levels = c("WT", "KO"))
ccc_box2 = ggplot (ccomp_df, aes (x= Genotype2, y= Freq)) +
  geom_boxplot(aes(fill = Genotype2), outlier.size=.2)  + 
  theme_classic() + 
  theme (axis.text.x = element_text(vjust = 0.5, hjust=0.5),
         legend.position = "none", axis.title.x = element_blank()) +
  xlab (metaGroups[2]) + 
  ylab ('Cell Cycling Score') + 
  geom_point(aes_string(fill = 'Genotype2'), size = 1.5, shape = 21, position = position_jitterdodge()) +
  ylim (c(0, 0.008)) +
  scale_fill_manual(breaks = c("WT", "KO"), values= genotype_pal)

stat.test2 = ccomp_df %>%
  t_test(Freq ~ Genotype2) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test2 <- stat.test2 %>% add_xy_position(step.increase = 0.0001, dodge = 1)

ccc_box2 <- ccc_box2 + stat_pvalue_manual(stat.test2, label = "p.adj", remove.bracket = F,
                                          bracket.nudge.y = 0.0011, hide.ns = F)
ggsave(paste0(projdir, 'cc_genotype_composition2_', date, '.pdf'), ccc_box2, device='pdf', width = 700, height = 850, units = "px")

#############################################

#######################################################################
# PART 3: Differential Expression Analysis
#######################################################################

# (3a) Run differential expression analysis using Muscat
#############################################
# Load in Seurat object
srt <- readRDS("/meyersonlab/akansha/Alena/results/labeled_seurat_object.rds")
srt$Genotype2 = ifelse (srt$Genotype == 'Cmtr2 WT','WT','KO')

# Convert seurat obj to singleCellExperiment obj
counts <- srt@assays$RNA@layers$counts
metaGroupName1 = 'sampleID'
metaGroupName2 = 'predicted_JShendure_Inf' # specify here the metadata of predicted celltypes from JShendure
metaGroupName3 = 'Genotype2' # Genotype2 is = Genotype but with labels with special characters removed

ei = data.frame (
  sample_id = as.factor(as.character(srt@meta.data[,metaGroupName1])),
  cluster_id = as.factor(as.character(srt@meta.data[,metaGroupName2])),
  group_id = as.factor(as.character(srt@meta.data[,metaGroupName3])))
rownames(ei) = colnames (srt)
levels(ei$group_id) <- c("WT", "KO")

genes <- data.frame(gene = rownames(srt))
rownames(genes) <- genes$gene

sce = SingleCellExperiment (list(counts = counts),
                            colData=ei,
                            # rowData=rownames(srt),
                            rowData=genes,
                            metadata=list(experiment_info = ei))

ds_pb <- readRDS("/meyersonlab/akansha/Alena/results/ds_pb.rds")
frq <- readRDS("/meyersonlab/akansha/Alena/results/frq.rds")

gids = levels (sce$sample_id)
frq10 = vapply (as.list(assays(frq)), 
                function(u) apply(u[, gids] > 0.01, 1, any), 
                logical(nrow(sce)))

# Filter low % expressed genes
tbl_fil1 = lapply(names(tbl), function(k)
  dplyr::filter(tbl[[k]], 
                gene %in% names(which(frq10[, k]))))
ds_pb[['pch_filtered']] = tbl_fil1

tbl_df = do.call (rbind, ds_pb[['table']][[1]])
colnames (tbl_df)[colnames (tbl_df) == 'logFC'] = 'avg_log2FC'
colnames (tbl_df)[colnames (tbl_df) == 'p_adj.loc'] = 'p_val_adj'
colnames (tbl_df)[colnames (tbl_df) == 'cluster_id'] = 'cluster'

write.table(tbl_df, file = "master_muscat_table.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

sig_muscat <- filter(tbl_df, p_val_adj < 0.05 & abs(avg_log2FC) > 1)

write.table(sig_muscat, file = "sig_muscat_table.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
#############################################

# (3b) Number of differentially expressed genes by cell type
#############################################
sig_muscat <- read.table("sig_muscat_table.tsv", sep = "\t", header = TRUE)
sig_muscat <- sig_muscat %>% dplyr::select(cluster, gene, avg_log2FC, p_val_adj)
colnames(sig_muscat)[3:4] <- c("log2FoldChange", "padj")

genotype_pal <- c('#EF8A62', '#67A9CF')

num_deg_per_cluster_plot <- function(sig_genes, method = "method") {
  pos_genes <- as.data.frame(sig_genes %>% group_by(cluster) %>% summarise(count_pos = sum(log2FoldChange > 1)) %>% arrange(desc(count_pos)))
  neg_genes <- as.data.frame(sig_genes %>% group_by(cluster) %>% summarise(count_neg = -1 * sum(log2FoldChange < 1)))
  
  joint_table <- full_join(pos_genes, neg_genes, by = "cluster")
  joint_table$cluster <- factor(joint_table$cluster, levels = joint_table$cluster)
  
  deg_bar = ggplot(joint_table, aes(x = cluster)) +
    geom_bar(aes(y = count_pos), fill = genotype_pal[1],color='black', stat = "identity") +
    geom_bar(aes(y = count_neg), fill = genotype_pal[2],color='black', stat = "identity", position = "identity") +
    ylim (c(-max(abs(na.omit(c(joint_table$count_neg, joint_table$count_pos)))),
            max(na.omit(abs(c(joint_table$count_neg, joint_table$count_pos)))))) +
    labs(title = "Number of DEGs Per Cell Type", y = "Count") +
    theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = "black"),
                            axis.text.y = element_text(color = "black"),
                            axis.title.x = element_blank())
  
  ggsave(paste0("/meyersonlab/akansha/Alena/results/plots/number_deg_cluster_sig_", method, "_", date, ".pdf"), deg_bar, device='pdf', width = 8, height = 5)
  
  return(joint_table)
}

muscat_table <- num_deg_per_cluster_plot(sig_muscat, method = "muscat")
#############################################

# (3c) Heatmap for top differentially expressed genes
#############################################
all_muscat <- read.table("master_muscat_table.tsv", sep = "\t", header = TRUE)
all_muscat <- all_muscat %>% dplyr::select(cluster, gene, avg_log2FC, p_val_adj)
colnames(all_muscat)[3:4] <- c("log2FoldChange", "padj")

# Make heatmap of the expression difference of between two conditions (e.g. KO vs WT) aggregated by cluster and top DEG across the clusters
# def_genes_df must contain following columns: "p_val_adj", "avg_log2FC", gene", "cluster"
# addGene specifies whether to add gene of interest in heatmap
# 05/20/24: Update so only includes genes up through Gsdme:
# Eda2r, Ano3, Ccdc7a, Ptprt, Tmtc1, Zmat3, Pvt1, Samd12, Map3k20, Serpine2,
# Cdkn1a, Svop, Lrrc69, Muc15, Armh2, Phlda3, Sned1, Kif27, Dscaml1, Arap2,
# Ccng1, Plin5, Ahnak2, Angpt2, Tmem132cos, Zfp746, Gsdme
diffClustHeat = function (deg_genes_df = NULL, topGenes = 20, pvalAdjThreshold = 0.05, sort_by = 'padj',
                          only_pos = F, addGene = NULL, plotcol = NULL, col_limit = NULL,
                          return_mat = F, ...) {
  deg_genes_df <- filter(deg_genes_df, !is.na(padj))
  deg_genes_df1 <- deg_genes_df
  sig_genes <- unique (deg_genes_df1$gene[deg_genes_df1$padj < pvalAdjThreshold])
  deg_genes_df <- deg_genes_df1[deg_genes_df1$gene %in% sig_genes, ]
  
  geneNames <- c('Eda2r', 'Ano3', 'Ccdc7a', 'Ptprt', 'Tmtc1', 'Zmat3', 'Pvt1',
                 'Samd12', 'Map3k20', 'Serpine2', 'Cdkn1a', 'Svop', 'Lrrc69',
                 'Muc15', 'Armh2', 'Phlda3', 'Sned1', 'Kif27', 'Dscaml1', 'Arap2',
                 'Ccng1', 'Plin5', 'Ahnak2', 'Angpt2', 'Tmem132cos', 'Zfp746', 'Gsdme')
  
  # if (sort_by == 'padj') {
  #   geneNames = unique (unlist (lapply (split (deg_genes_df, deg_genes_df$cluster), function(x) 
  #   {
  #     x = x[order (x[,sort_by]),]	
  #     head(x, topGenes)$gene
  #   })))
  # }
  # if (sort_by == 'log2FoldChange') {
  #   geneNames = unique (unlist (lapply (split (deg_genes_df, deg_genes_df$cluster), function(x) {
  #     x = x[order (-x[,sort_by]),]
  #     if (only_pos) {
  #       head(x, topGenes)$gene
  #     } else {
  #       c(head(x, topGenes)$gene, tail(x, 2)$gene)				
  #     }	
  #   })))
  #   print(geneNames)
  # }
  
  # Add gene specified in 'addGene' argument if needed
  if (!is.null(addGene)) geneNames = unique(append (geneNames, addGene))
  # Make mat of pvalues
  pval_mat = split (deg_genes_df1, deg_genes_df1$cluster)
  pval_mat = as.data.frame(sapply (pval_mat, function(x) x$padj[match (geneNames, x$gene)]))
  if (ncol(pval_mat) == 1) pval_mat = as.data.frame(t(pval_mat))
  rownames(pval_mat) = geneNames
  pval_mat[is.na(pval_mat)] = 1
  
  lfc_mat = split (deg_genes_df1, deg_genes_df1$cluster)
  lfc_mat = as.data.frame (sapply (lfc_mat, function(x) x$log2FoldChange[match (geneNames, x$gene)]))
  if (ncol(lfc_mat) == 1) lfc_mat = as.data.frame(t(lfc_mat))
  rownames (lfc_mat) = geneNames
  lfc_mat[is.na(lfc_mat)] = 0
  
  if (is.null(col_limit)) col_limit = max (abs(lfc_mat)) 
  if (is.null(plotcol)) {	
    plotcol = colorRamp2(c(-max (abs(lfc_mat)), 0, max (abs(lfc_mat))), c("green", "white", "red"))
  } else {
    plotcol = colorRamp2(c(-col_limit,0, col_limit), plotcol)
  }
  labels_annotation = HeatmapAnnotation(
    text = anno_text(rownames(lfc_mat), rot = 45, location = unit(1, "npc"), just = "right",gp =gpar(fontsize = 6)),
    annotation_height = max_text_width(rownames(lfc_mat))
  )
  ht_opt$HEATMAP_LEGEND_PADDING <- unit(0.6, "cm")
  ht = Heatmap (t(lfc_mat), 
                clustering_distance_columns = 'euclidean',
                clustering_distance_rows = 'euclidean',
                column_title = 'Top Enriched Genes CMTR2 KO vs WT',
                col = plotcol,
                heatmap_legend_param = list(title = "LFC", direction = "horizontal",
                                            title_position = "lefttop"), 
                bottom_annotation = labels_annotation,
                column_names_gp = gpar(fontsize = 0),
                row_names_gp = gpar(fontsize = 6),
                cell_fun = function (j, i, x, y, width, height, fill)
                {
                  if(t(pval_mat)[i, j] < pvalAdjThreshold) {
                    grid.text("*", x, y, just='center', vjust=.8,
                              gp = gpar(fontsize = 6, col='black'))
                  }
                }, ...
  )
  if (!return_mat) return (ht) 
  else return (list (lfc = lfc_mat,pval = pval_mat))
}

# Plot diff heatmap of significant genes
method <- "muscat"
tbl_df <- all_muscat
  
# Remove genes starting with "Gm-" or ending with "Rik" -- are annotated genes
# that don't have a canonical name (yet)
tbl_df = tbl_df[!grepl ('Gm', tbl_df$gene),]
tbl_df = tbl_df[!grepl ('Rik', tbl_df$gene),]
  
dhm = diffClustHeat (
  deg_genes_df = tbl_df,
  sort_by = 'log2FoldChange',
  topGenes = 5, # topGenes, 
  pvalAdjThreshold = 0.05,
  col_limit = NULL,
  plotcol = rev(RColorBrewer::brewer.pal(3,'RdBu')),
  cluster_columns=T,
  cluster_rows=T,
  border=T,
  only_pos = T)
  
pdf (paste0("/meyersonlab/akansha/Alena/results/plots/deg_heatmap_top_5_pos_", method, "_", date, ".pdf"), width = 8, height= 5)
draw (dhm, heatmap_legend_side = "bottom")
dev.off()	
#############################################

# (3d) Volcano plot of top differentially expressed genes
#############################################
all_muscat <- read.table("/meyersonlab/akansha/Alena/results/DESeq2/master_muscat_table.tsv", sep = "\t", header = TRUE)
endothelial <- all_muscat %>% filter(cluster == "Endothelial cells")

qval_most_sig <- c(endothelial %>% filter(abs(avg_log2FC) >= 1, p_val_adj < .05) %>%
                     arrange(p_val_adj) %>% .$gene %>% .[1:10],
                   endothelial %>% filter(abs(avg_log2FC) >= 1, p_val_adj < .05) %>%
                     arrange(desc(abs(avg_log2FC))) %>% .$gene %>% .[1:10]) %>% unique()

volcano_plot <- ggplot(data = endothelial %>% select(gene, avg_log2FC, p_val_adj) %>% na.omit,
                       aes(x = avg_log2FC, y = -log10(p_val_adj))) + 
  geom_point(aes(colour=ifelse(gene %in% qval_most_sig, "red", "black")), size = 3) +
  scale_color_identity() +
  geom_vline(xintercept=c(-1, 1), col="red", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "dashed") +
  geom_text_repel(aes(label=ifelse(gene %in% qval_most_sig, gene, '')), hjust=0, vjust=0, size = 5, fontface = "italic") +
  labs(title = "Endothelial Cells") + #CHANGE title
  theme_classic() + 
  xlab(expression(log[2]*FoldChange)) +
  ylab(expression(-log[10]*(p[adj]))) +
  theme(
    plot.title = element_text(size=16, color=1, face="bold", hjust=0.5, margin = margin(t=0, r=5, b=15, l=5, unit="pt")),
    axis.text.x=element_text(size=14, color=1, hjust=0.5, margin = margin(t=15, r=0, b=0, l=0, unit="pt")),
    axis.title.x=element_text(size=16, color=1, margin=margin(t=5, r=0, b=5, l=0, unit="pt")),
    axis.text.y=element_text(size=14, color=1, margin = margin(t=0, r=15, b=0, l=0, unit="pt")),
    axis.title.y=element_text(size=16, color=1, margin = margin(t=0, r=5, b=0, l=5, unit="pt")),
    legend.position="None",
    strip.text=element_text(size=14, color=1),
    strip.background = element_blank(), 
    strip.placement = "outside",
  )

show(volcano_plot)

ggsave(paste0("/meyersonlab/akansha/Alena/results/plots/endothelial_cells_volcano_plot_", date, ".pdf"), volcano_plot, device='pdf', width = 8, height = 5)

#######################################################################
# PART 4: Apoptosis Analysis
#######################################################################

# (4a) Bulk RNA seq apoptosis analysis 
#############################################
bulk_orig <- read_excel("data/Cmtr2KOvsCmtr2WT_DEG.xlsx")

# Apoptotic signaling
apop_gmt_file <- 'CMTR2_Scripts/data/h.all.v7.1.symbol.gmt'
pathways <- gmtPathways(apop_gmt_file)
apop_genes <- pathways$HALLMARK_APOPTOSIS # 156 genes

bulk_apop <- bulk_orig %>% filter(gene_name %in% apop_genes) # 151 of 156 genes in list
bulk_apop_mat <- bulk_apop[, c(13, 2:9)]
bulk_apop_mat$WT_avg <- rowMeans(bulk_apop_mat[, 6:9], )
bulk_apop_mat$KO_avg <- rowMeans(bulk_apop_mat[, 2:5], )

bulk_apop_lfc <- bulk_apop[, c(13, 10, 12)]
bulk_apop_lfc$padjust <- as.numeric(bulk_apop_lfc$padjust)
bulk_apop_lfc$sig <- ifelse(bulk_apop_lfc$padjust < 0.05 & abs(bulk_apop_lfc$log2foldchange) > 1, TRUE, FALSE)

bulk_apop_sig_lfc <- filter(bulk_apop_lfc, sig == TRUE)

write.table(bulk_apop_lfc, file = "results/apoptosis/bulk_apop_lfc.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(bulk_apop_sig_lfc, file = "results/apoptosis/bulk_apop_lfc_sig.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
#############################################

# (4b) snRNA seq apoptosis analysis 
#############################################
muscat_orig <- read.table("results/DESeq2/master_muscat_table.tsv", sep = "\t", header = TRUE)
sc_apop <- muscat_orig %>% filter(gene %in% apop_genes) # 151 of 156 genes in list

sc_apop_lfc <- as.data.frame(sc_apop %>% dplyr::select(gene, cluster, avg_log2FC, p_val_adj))
sc_apop_lfc <- sc_apop_lfc[with(sc_apop_lfc, order(gene, cluster)), ]
sc_apop_lfc$sig <- ifelse(sc_apop_lfc$p_val_adj < 0.05 & abs(sc_apop_lfc$avg_log2FC) > 1, TRUE, FALSE)

sc_apop_sig_lfc <- filter(sc_apop_lfc, sig == TRUE)

write.table(sc_apop_lfc, file = "results/apoptosis/sc_apop_lfc.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(sc_apop_sig_lfc, file = "results/apoptosis/sc_apop_lfc_sig.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
#############################################


