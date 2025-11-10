replicates\_&\_integration
================
Shounak Kadam
2025-11-08

## Overview 🧬

This sessioon focuses on handling multiple single-cell RNA-seq
replicates — a critical step for producing biologically meaningful
insights in large-scale experiments. In practice, single-sample analyses
are rarely sufficient; integrating data from several replicates or
experimental conditions ensures that observed patterns reflect true
biology rather than technical variation. The workflow presented here
mirrors real-world single-cell analysis pipelines used in research and
industry settings. It demonstrates how to structure a clear,
reproducible, and version-controlled integration process in Seurat (v5),
with full annotation and interpretability.

## Learning & Analysis Goals

• Efficiently load and manage multiple Seurat objects (biological
replicates)

• Assess pre-integration variation using unintegrated UMAPs

• Perform data integration using Seurat’s anchor-based approach to
correct batch effects

• Re-cluster and re-annotate integrated data to confirm alignment of
shared cell populations

## Part 1 - Loading Multiple Samples

### Data overview

In this experiment, we analyse human peripheral blood mononuclear cells
(PBMCs)

Unlike earlier sample analyses, this dataset includes four independant
scRNA-seq samples representing two biological conditions:

Control - 2 replicates - Untreated PBMC samples

IFNG treated - 2 replicates - PBMCs stimulated with interferon gamma

### Results Images
Please Note: All the images for this particular workflow are stored in the figures_generated folder inside this directiory

## Set the working directory to the 1k folder and load each sample individually, giving each a unique name

``` r
con_r1=readRDS("1K/con_r1.rds")
con_r2=readRDS("1K/con_r2.rds")
treat_r1=readRDS("1K/treat_r1.rds")
treat_r2=readRDS("1K/treat_r2.rds")
```

## Task 1.2 - Adding group information to metadata

Before merging or integrating, it’s essential to include group-level
annotations (e.g., Control vs IFNG treated, or con vs alz). Doing this
now ensures the information is preserved throughout all downstream
analyses and prevents errors or confusion later.

``` r
con_r1 = AddMetaData(con_r1,"con",col.name="sample_group")
con_r2 = AddMetaData(con_r2,"con",col.name="sample_group")
treat_r1 = AddMetaData(treat_r1,"ifng",col.name="sample_group")
treat_r2 = AddMetaData(treat_r2,"ifng",col.name="sample_group")
```

## Task 1.3 - Combining All samples into one seurat object

Now that each sample has its sample_group metadata column, the next step
is to combine all four Seurat objects into a single object that holds
every cell from every replicate and condition. This merged object will
serve as the foundation for integration and downstream comparative
analysis.

``` r
sc.data = list("con_r1" = con_r1, "con_r2" = con_r2, "treat_r1" = treat_r1, "treat_r2" = treat_r2)

#Merging the first sample to a list of samples
sc.data = merge(x = sc.data$con_r1, y = sc.data[-1])
```

### The merged Seurat object (sc.data) contains 4,000 PBMC cells from four replicates, correctly annotated with both sample and condition information. It includes two RNA layers (counts, data), confirming successful consolidation and readiness for integration.

## Task 1.4 - As we no longer need the per sample Seurat objects, we remove them to save space. for this we use the rm() function and enter the objects to get rid off. Then we run the garbage collection -gc() to free the space immediately

``` r
rm(con_r1,con_r2,treat_r1,treat_r2)
gc()
```

    ##            used  (Mb) gc trigger  (Mb) limit (Mb) max used  (Mb)
    ## Ncells  8233767 439.8   14793905 790.1         NA 10069355 537.8
    ## Vcells 23041392 175.8   47192004 360.1      16384 38819489 296.2

## Task 1.5 - Quality Control (QC) Across All Samples

Now that all four samples are merged into a single Seurat object
(sc.data), we’ll perform standard QC checks to evaluate mitochondrial
content and overall cell quality across replicates. Because the merged
object retains metadata such as orig.ident, each plot will automatically
separate the four samples — allowing visual comparison of QC metrics per
sample.

``` r
# get% mitochondria
sc.data = PercentageFeatureSet(sc.data,"^MT", col.name="percent_mito")

# plot to determine filtering cut-offs
png("QC_feature_violin.png", height=1000,width=1000)
VlnPlot(sc.data,features=c("nCount_RNA","nFeature_RNA","percent_mito"),pt.size = 0)+NoLegend()+ggtitle("Violin plot for RNA counts, features and mitochondrial percentage")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# Scatter QC plots
png("QC_feature_scatter.png", height=750,width=750)
FeatureScatter(sc.data,feature1="nCount_RNA",feature2="percent_mito")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

# Task 1.6- filter RNA and set cuttoffs by analysing the violin plot and the scatter plot

``` r
sc.data = subset(sc.data, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 & percent_mito < 8)
sc.data
```

    ## An object of class Seurat 
    ## 14053 features across 3972 samples within 1 assay 
    ## Active assay: RNA (14053 features, 0 variable features)
    ##  2 layers present: counts, data

\#Cells are filtered using thresholds of nFeature_RNA \> 200,
nFeature_RNA \< 1500, and percent_mito \< 8, which effectively removes
low-quality, doublet, or high-mitochondrial cells. These conservative
limits ensure a clean dataset for downstream integration and
differential analysis.

\#After QC filtering, the merged Seurat object (sc.data) contains 14,053
detected genes across 3,972 high-quality cells in one RNA assay. Two
layers are available — counts (raw data) and data (normalized) — and
variable feature selection will be performed in the next step before
integration.

## Task 2.1 - 🧬 Task 2 – Unintegrated UMAP

\##Normalisation using SCTransform

Before generating an unintegrated UMAP, we need to normalize the merged
dataset so that gene expression values are comparable across cells.
Although Seurat supports several normalization methods
(e.g. NormalizeData() with log-normalization), here we use SCTransform,
a model-based normalization approach that performs variance
stabilization and removes sequencing-depth effects. This method is
preferred because it provides: Better control of technical noise,
Improved detection of variable features, Compatibility with Seurat’s
integration workflow

``` r
# Perform SCTransform normalization on the merged datset
sc.data = SCTransform(sc.data, vars.to.regress = "percent_mito", verbose = TRUE)
```

## Task 2.2 - After normalization, the next step is dimensionality reduction.

We use Principal Component Analysis (PCA) to capture the major sources
of variation in the SCT-normalized data. However, because we will later
perform another PCA after integration, we explicitly name this reduction
pca_unintegrated to avoid overwriting the default “PCA” slot.

``` r
#Perform PCA on the SCTnormalised data
sc.data = RunPCA(object = sc.data, assay = "SCT", reduction.name = "pca_unintegrated")
#pca_unintegrated will be saved in the reduction slot in the sc.data
```

## Task 2.3 - Make all plots and decide the optimal dimensions to carry forwards

``` r
#save PCA loadings
png("pca_unintegrated_loadings.png", height = 800, width = 1000)
VizDimLoadings(sc.data, reduction = "pca_unintegrated", dims = 1:2)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
#The VizDimLoadings plot shows that the top genes contributing to PC1 include ISG15, CXCL10, CCL2, and other interferon-stimulated genes, reflecting a strong IFNG-treatment response. PC2 is dominated by GZMB, NKG7, CD14, S100A8/A9, which mark cytotoxic and myeloid populations, indicating cell-type heterogeneity.Together, these components confirm that the PCA is capturing both biological signal (treatment and immune identity) rather than technical noise — validating the quality of the unintegrated dataset before UMAP visualization.

#save ElbowPlot
png("pca_unintegrated_elbow.png", height = 800, width = 1000)
ElbowPlot(sc.data, reduction = "pca_unintegrated", ndims  = 50)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
#The ElbowPlot shows a clear inflection around PC20, beyond which the explained variance plateaus.Therefore, 20 principal components were selected for downstream clustering and UMAP visualization to capture both treatment-driven and cell-type variation while avoiding overfitting.

#save PCA scatter plot 
png("pca_unintegrated_scatter.png")
DimPlot(sc.data, reduction = "pca_unintegrated")+
  ggtitle("Unintegrated PCA: Cells coloured by conddition")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# The unintegrated PCA shows clear separation between IFNG-treated and control PBMC samples along PC1, confirming that treatment drives major transcriptional variance.Within each condition, mild replicate-specific clustering is visible, suggesting technical batch effects that will be corrected in the subsequent integration step.Overall, the PCA confirms that the dataset captures strong biological signal while maintaining balanced replicate quality.

#save PCA heatmap
png("pca_unintegrated_heatmap.png", height = 1000 , width = 1200)
DimHeatmap(sc.data, reduction = "pca_unintegrated", dims = 1:15)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
#The PCA heatmap for the first 15 components highlights interferon-response genes (e.g., CXCL10, ISG15, IFIT1) and immune lineage markers (GZMB, CD14, S100A8/A9) driving the main axes of variation.
#Clear structured colour gradients across cells indicate strong biological signal, while higher PCs show diffuse patterns consistent with noise.
#This confirms the unintegrated PCA successfully captures both treatment-induced and cell-type variation, validating the choice of ~30 PCs for downstream UMAP and clustering.
#

#decide dims to use, set as a parameter
dims_to_use = 25
```

## Task 2.4 - Now we cluster and UMAP. Again, we need to explicitly state that we are using the reduction pca_unintegrated. For RunUMAP() we also need to deliberatley name the UMAP that is created

``` r
sc.data = FindNeighbors(object=sc.data, reduction = "pca_unintegrated", dims = 1:dims_to_use)
sc.data = FindClusters(object=sc.data, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3972
    ## Number of edges: 142752
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9014
    ## Number of communities: 11
    ## Elapsed time: 0 seconds

``` r
sc.data = RunUMAP(sc.data, dims = 1:dims_to_use, reduction = "pca_unintegrated", reduction.name = "umap_unintegrated", min.dist = 0.2, spread = 0.5)
```

## Task 2.5 - Plotting the UMAPSs. One is grouped by the orig id’s of the 4 samples. The other one is grouped by the number of clusters

``` r
png("UMAP_unintegrated.png", height = 1000, width = 2000)
p1 = DimPlot(sc.data, reduction = "umap_unintegrated", group.by = "orig.ident")
p2 = DimPlot(sc.data, reduction = "umap_unintegrated", group.by = "seurat_clusters")
p1 + p2
dev.off()
```

    ## quartz_off_screen 
    ##                 2

\##The unintegrated UMAP reveals that although the overall structure
reflects the IFNG treatment effect and immune cell diversity,
substantial per-sample batch effects remain. Replicates form parallel,
non-overlapping layers, and identical cell types appear as separate
clusters across samples.Integration is therefore required to align these
datasets, remove technical variance, and recover unified biological
structure.

## Task 3 - Integration

## Task 3.1 - To integrate we take the PCA and effectively “correct” by sample. By far the best and quickest tool for this is Harmony. We choose a reduction (e.g. pca_unintegrated) and an item in the meta.data to harmonise on, in this case orig.ident. We select a harmonise strength – theta, where the default is 3 and \< 3 means stronger “correction”. Then we choose to save the harmonised reduction as e.g. pca_integrated

``` r
sc.data = RunHarmony(sc.data, assay.use = "SCT", group.by.vars = c("orig.ident"), theta = c(3), reduction = "pca_unintegrated", reduction.save = "pca_integrated")
```

## Task 3.2 - Rerunning clustering and UMAP for integration

``` r
sc.data = FindNeighbors(object=sc.data, reduction = "pca_integrated", dims = 1:dims_to_use)
sc.data = FindClusters(object=sc.data, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3972
    ## Number of edges: 150711
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8972
    ## Number of communities: 12
    ## Elapsed time: 0 seconds

``` r
sc.data = RunUMAP(sc.data, dims = 1:dims_to_use, reduction = "pca_integrated", reduction.name = "umap_integrated", min.dist = 0.2, spread = 0.5, labels = TRUE)
```

## Task 3.3 - Once again replotting the UMAPSs. One is grouped by the orig id’s of the 4 samples. The other one is grouped by the number of clusters. This time we save it as “umap_integrated”

``` r
png("UMAP_integrated.png", height = 1000, width = 2000)
p1 = DimPlot(sc.data, reduction = "umap_integrated", group.by = "orig.ident", label = TRUE)
p2 = DimPlot(sc.data, reduction = "umap_integrated", group.by = "seurat_clusters", label = TRUE)
p1 + p2
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Task 4 - Identification of cell types

## Task 4.1 - Marker Optimization and Preparation for Cell Type Identification. Once the UMAP embedding and clustering look satisfactory, the next step is to assign biological cell identities to the detected clusters. This dataset contains human peripheral blood mononuclear cells (PBMCs), which typically include T cells, B cells, NK cells, and monocytes.Before identifying marker genes, we first need to optimise the SCT assay for differential expression testing using PrepSCTFindMarkers()

``` r
#prep markers - needed to combine the 4 samples 
sc.data = PrepSCTFindMarkers(sc.data)

#get the marker genes
sc.data.markers = FindAllMarkers(sc.data, only.pos = TRUE)

#get the top 15 markers for each cluster
sc.data.markers%>%group_by(cluster)%>%dplyr::filter(avg_log2FC>1)%>%slice_head(n=15)%>%ungroup()->topn_markers

#heatmap
png("Markers_heatmap.png",height = 3000, width = 3000)
DoHeatmap(sc.data, features = topn_markers$gene)+NoLegend()
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Task 4.2 - Reference based auto-annotation with SinleR

``` r
BiocManager::install("SingleR")
BiocManager::install("celldex")
library(SingleR)
library(celldex)

#Download available human primary cell atlas markers
hpca.se=celldex::HumanPrimaryCellAtlasData()

#To run single R we need to extract the normalised expression data from Seurat. Please note scaled data won't work, as it needs to know relative expression of genes. SingleR compares relative expression; it needs log-normalized data 

exp.ma=sc.data@assays$SCT@data

#Now to get the predictions

predictions = SingleR(test = exp.ma, ref = hpca.se, labels = hpca.se$label.fine)

sc.data = AddMetaData(sc.data, predictions$labels, col.name = "predictions.celltype")
```

### plot the UMAP

``` r
png("UMAP_SingleR.png",height=1000,width=2000)
p1 = DimPlot(sc.data,group.by="predictions.celltype",reduction="umap_integrated",label=TRUE, repel = TRUE, label.size=5)+NoLegend()
p2 = DimPlot(sc.data,group.by="seurat_clusters",reduction="umap_integrated",label=TRUE,label.size=5)+NoLegend()
p1 + p2
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### Combining all information, update the meta-data with your final cell types

``` r
new.cluster.ids = c("Monocyte/monocyte-derived_DC", "CD8+_T_cell", "Macrophage(M-CSF/IFNa)","CD4+_T_cell", "B_cell", "CD8+_T_cell(effector/memory)", "CD16+_monocyte", "NK/ᵞ𝛅T_cell","platelets","erythrocytes","erythrocytes","platelets")

length(new.cluster.ids)
```

    ## [1] 12

``` r
length(levels(sc.data))
```

    ## [1] 12

``` r
names(new.cluster.ids)=levels(sc.data)
sc.data = RenameIdents(sc.data, new.cluster.ids)
sc.data = AddMetaData(sc.data, sc.data@active.ident,col.name="cell.types")
```

### plot UMAP

``` r
png("Cells_UMAP.png",height=1500,width=1500)
DimPlot(sc.data,group.by="cell.types",reduction="umap_integrated",label=TRUE,pt.size=1,label.size = 10)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

### frequencies

``` r
table(sc.data@active.ident,sc.data@meta.data$orig.ident)
```

    ##                               
    ##                                con_r1 con_r2 treat_r1 treat_r2
    ##   Monocyte/monocyte-derived_DC    345    320        9       13
    ##   CD8+_T_cell                     152    138      167      162
    ##   Macrophage(M-CSF/IFNa)            2      1      251      260
    ##   CD4+_T_cell                      78     81      154      187
    ##   B_cell                           85     92      109       89
    ##   CD8+_T_cell(effector/memory)     66     73       80       79
    ##   CD16+_monocyte                   86     77       63       66
    ##   NK/ᵞ𝛅T_cell                      63     78       81       60
    ##   platelets                        71     88       52       51
    ##   erythrocytes                     43     43       32       25

### plot frequencies

``` r
png("Cells_Frequency.png",height=1000,width=1000)
dittoBarPlot(object=sc.data,var="cell.types",group.by="orig.ident")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Interpretation

\###A proportional cell-type analysis using dittoBarPlot revealed clear
shifts in immune cell composition following IFNG treatment.Control
replicates (con_r1, con_r2) were enriched in monocytes and CD4⁺ T cells,
while IFNG-treated samples (treat_r1, treat_r2) showed increased
proportions of CD8⁺ T cells and macrophage populations. These findings
indicate an IFNG-driven activation response favouring cytotoxic and
pro-inflammatory cell states, consistent with established
interferon-mediated immune modulation.

## Task 5 - saving the seurat object to disk, so that it doesn’t have to be reloaded again

``` r
saveRDS(sc.data,file="s4_final.rds")
```
