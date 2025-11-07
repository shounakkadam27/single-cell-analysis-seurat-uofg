Identifying_cell_types
================
Shounak Kadam
2025-11-4

# Load processed object from Tutorial 2

``` r
# Load required libraries
library(Seurat)
```

    ## Loading required package: SeuratObject

    ## Loading required package: sp

    ## 'SeuratObject' was built under R 4.4.1 but the current version is
    ## 4.4.3; it is recomended that you reinstall 'SeuratObject' as the ABI
    ## for R may have changed

    ## 
    ## Attaching package: 'SeuratObject'

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, t

``` r
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(SeuratData)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
```

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following object is masked from 'package:SeuratObject':
    ## 
    ##     intersect

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:sp':
    ## 
    ##     %over%

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 

``` r
library(clusterProfiler)
```

    ## 

    ## clusterProfiler v4.14.6 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.
    ## clusterProfiler: an R package for comparing biological themes among
    ## gene clusters. OMICS: A Journal of Integrative Biology. 2012,
    ## 16(5):284-287

    ## 
    ## Attaching package: 'clusterProfiler'

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(msigdbr)
library(enrichplot)
```

    ## enrichplot v1.26.6 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## G Yu. Thirteen years of clusterProfiler. The Innovation. 2024,
    ## 5(6):100722

``` r
library(EnhancedVolcano)
```

    ## Loading required package: ggrepel

``` r
library(dittoSeq)
library(scDblFinder)
```

    ## Loading required package: SingleCellExperiment

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     rowMedians

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'SummarizedExperiment'

    ## The following object is masked from 'package:Seurat':
    ## 
    ##     Assays

    ## The following object is masked from 'package:SeuratObject':
    ## 
    ##     Assays

``` r
library(GOSemSim)
```

    ## GOSemSim v2.32.0 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## Guangchuang Yu, Fei Li, Yide Qin, Xiaochen Bo, Yibo Wu and Shengqi
    ## Wang. GOSemSim: an R package for measuring semantic similarity among GO
    ## terms and gene products. Bioinformatics. 2010, 26(7):976-978

``` r
library(dittoSeq)


# Load processed Seurat object from Tutorial 2
sc.data <- readRDS("../Tut2/output/seurat/02_sct_pca_umap_clusters.rds")
sc.data
```

    ## An object of class Seurat 
    ## 45257 features across 2638 samples within 2 assays 
    ## Active assay: SCT (12519 features, 3000 variable features)
    ##  3 layers present: counts, data, scale.data
    ##  1 other assay present: RNA
    ##  2 dimensional reductions calculated: pca, umap

## Overview 🧬

In this session we identify per-cluster biomarkers with
FindAllMarkers(), visualise canonical markers (DotPlot, FeaturePlot),
and annotate clusters with biologically meaningful cell types.

## Task 1.1 - Cluster Biomarkers

Objective Our goal is to identify genes that are uniquely upregulated in
each cluster. These serve as biomarkers that help determine the
biological identity of each cluster. We begin by testing one cluster
(cluster 2) against all other cells, followed by additional clusters.
Approach We use FindMarkers() from Seurat, which performs a differential
expression test comparing cells in one cluster (ident.1) to all other
clusters.

``` r
cluster1.markers = FindMarkers(sc.data, ident.1 = 2)
```

## ✅ Interpretation

This marker signature very strongly indicates: Cluster 2 = Classical
Monocytes / Inflammatory Monocytes

🧬 Biological reasoning

> S100A8/S100A9 → hallmark of early inflammatory myeloid cells

> CD14 → classical monocyte receptor

> FCN1 → monocyte–macrophage lineage

> LGALS2, TYROBP → innate immune activation genes

The top biomarkers for Cluster 2 included S100A8, S100A9, LGALS2, FCN1,
CD14 and TYROBP. These genes are strongly associated with the classical
monocyte lineage, particularly CD14⁺ inflammatory monocytes. Based on
this expression profile, we assigned Cluster 2 as Classical Monocytes.

## running code for other clusters

``` r
cluster0.markers=FindMarkers(sc.data, ident.1 = 1)
cluster2.markers=FindMarkers(sc.data, ident.1 = 3)
cluster3.markers=FindMarkers(sc.data, ident.1 = 4)
cluster4.markers=FindMarkers(sc.data, ident.1 = 5)
cluster5.markers=FindMarkers(sc.data, ident.1 = 6)
cluster6.markers=FindMarkers(sc.data, ident.1 = 7)
cluster7.markers=FindMarkers(sc.data, ident.1 = 8)
```

## Interpretation of cluster 1, 3, 4 and 5

### Cluster 0

Cluster 1 showed high expression of RPS27, RPL32, S100A4, RPS6, RPS12.
likely low-signal / mixed / QC Ribosomal high-cluster

### Cluster 1

Cluster 3 expressed CD79A, MS4A1, TCL1A, CD79B, LINC00926 characteristic
of classic B-cell lineage. TCL1A = transistional/naive B-cell

### Cluster 3

Cluster 4 displayed CCL5, NKG7, GZMK, CST7, GZMA. Cytotoxic T cell
signature (GZMK-high effector CD8 T cells)

### Cluster 4

Cluster 5 displayed GZMB, FGFBP2, SPON2, PRF1, GNLY. NK Cells signature

# Extract top 5 gene names for cluster 2 and plot using FeaturePlot()

``` r
top5_cluster2 = row.names(head(cluster2.markers, 5))
top5_cluster2
```

    ## [1] "CD79A"     "MS4A1"     "TCL1A"     "CD79B"     "LINC00926"

``` r
FeaturePlot(sc.data, features = top5_cluster2, reduction="umap", ncol = 3)
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Task 1.2 - We experimented with different colour scales and point aesthetics in FeaturePlot. Changing colour gradients, point size, and transparency can improve visual clarity when highlighting gene expression patterns across clusters. The modified plot confirms strong and distinct expression of the selected monocyte marker in the cluster of interest.

``` r
S100A8 = top5_cluster2[1]
S100A8
```

    ## [1] "CD79A"

``` r
FeaturePlot(sc.data, features = S100A8, cols = c("white", "red", "purple"), pt.size = 1, alpha = 0.5)
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
FCN1 = top5_cluster2[3]
FCN1
```

    ## [1] "TCL1A"

``` r
FeaturePlot(sc.data, features = FCN1, cols = c("white", "red", "purple"), pt.size = 1, alpha = 0.5)
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

## Task 1.3 - Violin plots for top 3 biomarkers of cluster 2

``` r
top3_cluster2=row.names(cluster2.markers)[1:3]
VlnPlot(sc.data, features = c(top3_cluster2), pt.size = 0)+NoLegend()
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Task 1.4 - We used FindAllMarkers() to compute biomarkers across all clusters simultaneously.

Setting only.pos = TRUE returns genes that are specifically upregulated
in each cluster relative to all others. The resulting table includes: •
cluster — which cluster the gene marks • gene — the marker gene •
p_val_adj — FDR‐adjusted significance • avg_log2FC — effect size • pct.1
/ pct.2 — % of cells in vs out of the cluster expressing the gene The
table shows clear sets of enriched markers per cluster, confirming
distinct transcriptional programs and supporting our manual cell‐type
assignments.

``` r
sc.data.markers = FindAllMarkers(sc.data, only.pos = TRUE)
```

## Task 1.5 - We can simplify this table using some complicated code from the Seurat website. Instead of the full list of biomarkers only the top 5-50 or so are useful. The topn_markers will have 15 biomarkers per cluster

``` r
sc.data.markers %>%group_by(cluster)%>%dplyr::filter(avg_log2FC>1)%>%slice_head(n=15)%>%ungroup()->topn_markers
```

## Task 1.6 - Plotting the most useful plot in cell type identification - a heatmap of the top 15 biomarkers by cell cluster

``` r
png("/Users/shounakkadam/Documents/Module C - Single cell and spatial transcriptomics/Tut3/Session_3_heatmap.png", height = 3000, width = 1500)
top15_genes=topn_markers$gene
DoHeatmap(sc.data, features = top15_genes)+NoLegend()
```

## Task 1.7 - Interpreting the heatmap vs UMAP

The heatmap displays the top markers across clusters, with each row
representing a gene and each coloumn representing a cell, grouped by
cluster. Inspecting the heatmap alongside the UMAP from Session 2 shows
string agreement between the two visualisations. Clusters that appear
spatially close on the UMAP also share highly similar transcriptional
profiles in the heamtap, whereas well-separated clusters show distinct
patterns of marker expression

The top-marker heatmap showed clear lineage-specific transcriptional
patterns and closely matched cluster separations seen in the UMAP.
Clusters that form a single broad “blob” on the UMAP (e.g., clusters 0
and 1) also exhibited similar heatmap signatures, suggesting that they
represent sub-states of the same major cell type, likely naive and
transitional B cells. Similarly, clusters 4 and 5 show overlapping
profiles with cytotoxic effector signatures, consistent with NK cell
subsets. Clusters 2 and 6 demonstrate strong S100A8/A9 and CD14
expression, confirming classical monocytes with potential
activation/state differences. Cluster 3 forms an isolated island and
shows strong GZMK/NKG7 expression, representing cytotoxic T cells. Tiny
clusters (7 and 8) likely reflect rare subpopulations or state
transitions. Overall, the heatmap and UMAP are highly consistent,
indicating ~4 major immune cell types with biologically meaningful sub
clusters.

Please Note: The UMAP of the clusters and the heatmap of the top 15
genes for each cluster have been saved in the “output_images” folder

## Task 2 - Identifying cell types

Our next step is to identify cell types from the biomarkers. Though
there are “automated” methods that can assisit you with this, by far the
best method is to manually curate cluster by comparing the bimarkers
from Task 1 and Task 2 to our own expert knowledge. This step required
careful thought, reviews and literature discussion with your collegues.
It is often iterative

## Task 2.1 - A simple table created manually for for identifying the cell types from the cluster

### Part 2 — Cluster to Cell Type Mapping

| Cluster | Assigned Cell Type                                  |
|:-------:|:----------------------------------------------------|
|    0    | Cytotoxic low quality T cells / Low-quality B-cells |
|    1    | T cells - CD4+, CD8+                                |
|    2    | T cells - CD4+, CD8+                                |
|    3    | B cells - Naive, immature, plasma,                  |
|    4    | NK cells                                            |
|    5    | NK cells                                            |
|    6    | Classical Monocytes                                 |
|    7    | Rare / Transitional subset                          |
|    8    | Rare / Low-quality / Unclassified                   |

## Task 2.2 - Reference based auto-annotation with SinleR

Initially used, the automated method SingleR to try to auto-assign cell
types to clusters. Still please note, it always needs to be manually
curated, and is dependant on having a suitable reference database.

``` r
#install and load library
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

sc.data = AddMetaData(sc.data, predictions$labels, col.name = "predictions")
```

## Task 2.3 - Inspect the object to see what SingleR produced

``` r
head(predictions$labels)
```

    ## [1] "T_cell:CD4+_central_memory"  "B_cell"                     
    ## [3] "T_cell:CD4+_central_memory"  "Monocyte:CD16-"             
    ## [5] "NK_cell"                     "T_cell:CD4+_effector_memory"

``` r
table(sc.data$predictions)
```

    ## 
    ##                         B_cell                B_cell:immature 
    ##                             19                            225 
    ##                  B_cell:Memory                   B_cell:Naive 
    ##                             17                             75 
    ##             B_cell:Plasma_cell                            CMP 
    ##                              3                              5 
    ##            DC:monocyte-derived                            GMP 
    ##                              1                              1 
    ##                       Monocyte                 Monocyte:CD16- 
    ##                              5                            380 
    ##                 Monocyte:CD16+                        NK_cell 
    ##                            279                            124 
    ##           NK_cell:CD56hiCD62L+                    NK_cell:IL2 
    ##                             31                              7 
    ##                      Platelets               Pre-B_cell_CD34- 
    ##                              3                             35 
    ##                    T_cell:CD4+     T_cell:CD4+_central_memory 
    ##                             27                            611 
    ##    T_cell:CD4+_effector_memory              T_cell:CD4+_Naive 
    ##                            390                            120 
    ##                    T_cell:CD8+    T_cell:CD8+_effector_memory 
    ##                            229                             41 
    ## T_cell:CD8+_effector_memory_RA              T_cell:CD8+_naive 
    ##                              4                              4 
    ##             T_cell:gamma-delta 
    ##                              2

## Task 2.3.1 - Visual checl on UMAP, if SingleR labels form coherent islands on UMAP, the’re liekly sensible

``` r
old_id = Idents(sc.data)
Idents(sc.data)= factor(sc.data$predictions)
DimPlot(sc.data, reduction = "umap", label = TRUE) +
ggtitle("SingleR (HPCA) - per cell labels")
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
Idents(sc.data) = old_id 
ggsave("output-images/03_umap_singleR_predictions.png", width = 10, height = 5, dpi = 300)
```

We can change the active.indent slot in the sc.data so that we can plot
the predictions instead of the clusters

``` r
sc.data=SetIdent(sc.data, value=sc.data$predictions)
```

We can change it manually to be anything in the meta.data. To change the
active.indent to be seurat_clusters again:

``` r
sc.data=SetIdent(sc.data, value=sc.data$seurat_clusters)
```

## Task 2.4 - As well as being able to change the active.ident manually we can specify the meta.data slot to use when making a DimPLot. We look at both the original clusters and the predictions side by side.

``` r
p1 = DimPlot(sc.data, reduction = "umap", group.by = "predictions", label = TRUE, pt.size = 3, raster = TRUE)+NoLegend()
p2 = DimPlot(sc.data,group.by="seurat_clusters",label=TRUE,label.size=3)+NoLegend()
p1 + p2
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

## Updated the Task 2 - cluster to cell type mapping table after closely matching the predicted and the seurat_cluster UMAPs

## Task 2.5 - Pathway analysis of cluster markers (using ClusterProfiler)

In this section we explore whether the marker genes for each cluster are
enriched for particular biological processes This serves as an
additional validation for the SingleR and manual annotations

Load the required packages

``` r
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "msigdbr"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
```

Pick the cluster to be selected. In this case 1

``` r
markers = rownames(subset(cluster1.markers, avg_log2FC > 1))
markers = bitr(markers,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
```

Run pathway analysis and plot

``` r
enriched = enrichGO(gene = markers$ENTREZID, OrgDb = org.Hs.eg.db, ont ="BP")
barplot(enriched,showCategory=10)
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

The enrichment analysis shows strong over-representation of:

• T cell differentiation Lymphocyte differentiation

• Lymphocyte differentiation

Both pathways are hallmarks of adaptive-immune T-cell biology. These
results support the annotation of cluster 1 as a T-cell population,
confirming that the marker-basedand SingleR based annotations are
biologically consistent

Alternative method using MSigDB C8 enrichment check

Get markers

``` r
cluster.markers.1 = rownames(subset(cluster1.markers, avg_log2FC > 1)) #cluster1
```

Load the databse

``` r
msigdbr_df=msigdbr(species="human",category="C8")
msigdbr_t2g=msigdbr_df%>%dplyr::distinct(gs_name,gene_symbol)%>%as.data.frame()
```

Run enricher and plot

``` r
enriched_1 = enricher(gene=cluster.markers.1,TERM2GENE = msigdbr_t2g)
barplot(enriched_1,showCategory = 10)
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

## Task 2.5.1 Performing same steps for all clusters

Get markers

``` r
cluster.markers.2 = rownames(subset(cluster2.markers, avg_log2FC > 1))
cluster.markers.3 = rownames(subset(cluster3.markers, avg_log2FC > 1))
cluster.markers.4 = rownames(subset(cluster4.markers, avg_log2FC > 1))
cluster.markers.5 = rownames(subset(cluster5.markers, avg_log2FC > 1))
cluster.markers.6 = rownames(subset(cluster6.markers, avg_log2FC > 1))
cluster.markers.7 = rownames(subset(cluster7.markers, avg_log2FC > 1))
```

Run enricher and plot

``` r
enriched_2 = enricher(gene=cluster.markers.2,TERM2GENE = msigdbr_t2g)
barplot(enriched_2,showCategory = 10) #cluster 2
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
enriched_3 = enricher(gene=cluster.markers.3,TERM2GENE = msigdbr_t2g)
barplot(enriched_3,showCategory = 10)#cluster 3
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->

``` r
enriched_4 = enricher(gene=cluster.markers.4,TERM2GENE = msigdbr_t2g)
barplot(enriched_4,showCategory = 10)#cluster 4
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-23-3.png)<!-- -->

``` r
enriched_5 = enricher(gene=cluster.markers.5,TERM2GENE = msigdbr_t2g)
barplot(enriched_5,showCategory = 10)#cluster 5
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-23-4.png)<!-- -->

``` r
enriched_6 = enricher(gene=cluster.markers.6,TERM2GENE = msigdbr_t2g)
barplot(enriched_6,showCategory = 10)#cluster 6
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-23-5.png)<!-- -->

``` r
enriched_7 = enricher(gene=cluster.markers.7,TERM2GENE = msigdbr_t2g)
barplot(enriched_7,showCategory = 10)#cluster 7
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-23-6.png)<!-- -->

## Task 2.6 - manual annotation

After performing unsupervised clustering and pathway enrichment, we
manually annotated each cluster by comparing its top marker genes
(identified using FindAllMarkers()) with canonical PBMC markers (Table
X). Manual curation is considered the gold standard for single-cell
annotation, as it relies on biological context and literature-supported
marker knowledge rather than pre-defined references.

✅ Validation of Manual Annotation The manually annotated identities
showed strong concordance with both: The unsupervised seurat_clusters
structure, indicating that transcriptomic clustering accurately
reflected cell-type distinctions. The reference-based predictions from
SingleR, which independently assigned nearly identical labels (Naive
CD4⁺ T, CD14⁺ Monocytes, B cells, CD8⁺ T, NK, FCGR3A⁺ Monocytes, and
Platelets). This agreement between data-driven clustering,
reference-based annotation, and manual biological curation reinforces
the robustness of the analysis and confirms that each cluster represents
a distinct immune cell population typically observed in PBMC datasets.

## Final Manual Annotation Table

| Cluster | Marker Evidence | Assigned cell type |
|:-------:|:----------------|:-------------------|
|    0    | IL7R, CCR7      | Naive CD4+ T cells |
|    1    | IL7R, CCR7      | Naive CD4+ T cells |
|    2    | CD14, LYZ       | CD14+ Monocytes    |
|    3    | MS4A1           | B cells            |
|    4    | CD8A            | CD8+ T cells       |
|    5    | GNLY, NKG7      | NK Cells           |
|    6    | FCGR3A, MS4A7   | FCGR3A+ Monocytes  |
|    7    | PPBP            | Platelets          |

## Task 2.7 - Finlaising Cell-Type Labels and Updating Metadata

### Goal: Use all the evidence gathered from the analysis (marker expression, enrichment results, manual annotation, and SingleR predictions) to finalise the cluster-to-cell-type mapping and store it permanently in the Seurat object’s metadata

### This ensures that a;; downstream analyses (differential expression, pathway analyis and visualisation) utilise the correct biological identities

``` r
new.cluster.ids = c("Naive_CD4_T_cells", "Naive_CD4_T_cells", "CD14_Monocytes", "B_cells", "CD8_T_cells","NK_cells", "FCGR3A_Monocytes", "Platelets", "Dendritic_cells")
names(new.cluster.ids)=levels(sc.data)
sc.data=RenameIdents(sc.data,new.cluster.ids)
sc.data=AddMetaData(sc.data, sc.data@active.ident, col.name="cell.types")
```

## Visual confirmation of cluster type correspondence

``` r
p1 = DimPlot(sc.data, group.by = "seurat_clusters", label = TRUE, label.size = 3)+NoLegend()
p2 = DimPlot(sc.data, group.by = "cell.types", label = TRUE, label.size = 3)+NoLegend()
p1 + p2
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

### Interpretaion: The two UMAPs have the same cluster structure. The left plot shows numeric cluster IDs, and the right plot shows your biological cell-type labels.

### Every cluster label aligns with the same UMAP region, and proves annotation is fully validated.

## Task 2.8 - Visualising cell-type frequencies

### To quantify and visualise the frequency of each annotated cell type across samples or conditions. This provides an overview of cellular compositions within the dataset and helps iidentify any sample-specific biases in cell representation.

``` r
dittoBarPlot(object=sc.data,var="cell.types",group.by="orig.ident")+ggtitle("Cell-type Frequency distribution among samples")
```

![](Identifying_cell_types_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

### The cell-type frequency distribution indicates that Naive CD4⁺ T cells and CD14⁺ Monocytes constitute the majority of cells in the Control_R1 sample, with smaller fractions of CD8⁺ T cells, NK cells, B cells, FCGR3A⁺ Monocytes, Dendritic cells, and Platelets.This composition closely matches expected PBMC cellular profiles, supporting the accuracy of the annotation and sample integrity.

### For raw counts

``` r
table(sc.data@active.ident,sc.data@meta.data$orig.ident)
```

    ##                    
    ##                     Control_R1
    ##   Naive_CD4_T_cells       1107
    ##   CD14_Monocytes           486
    ##   B_cells                  347
    ##   CD8_T_cells              347
    ##   NK_cells                 153
    ##   FCGR3A_Monocytes         152
    ##   Platelets                 34
    ##   Dendritic_cells           12

## Save the disk so that there is no need to run the entire sript once again

``` r
dir.create("output/seurat", recursive = TRUE, showWarnings = FALSE)
saveRDS(sc.data,file="output/seurat/03-identifying_cell_types.rds")
```
