# Task 1 - Part 1 - Loading all the tables in the environment
em = read.table("em.csv", header=TRUE,row.names=1,sep="\t")
anno = read.table("annotations.csv", header=TRUE,row.names=1,sep="\t")
de = read.table("de_duct_vs_gut.csv", header=TRUE,row.names=1,sep="\t")
ss = read.table("sample_sheet.csv", header=TRUE,row.names=1,sep="\t")

#Task 1 - Part 2 - Making a list to store items 
bulk_rna = list("em" = em, "anno" = anno, "de" = de, "ss" = ss)

#Task 1 - Part 4 - accessing individual items in the list
bulk_rna$em
bulk_rna$anno
bulk_rna$de

# Task 1 - Part 5 - adding new items to the bulk_rna list

assay = "Bulk RNA"
design = "RNA-seq profiling of dendritic cells from mouse gut, lymph duct and node. N=3"
bulk_rna$assay = assay
bulk_rna$design = design 

#Task 2 - Loading single cell datasets and handling Seurat objects

#Task 2 - Part 1
#The dataset we will use for the first 3 sesions is the 10x example dataset, of 
#2700 human peripheral blood mononuclear cells (PBMC), taken from a healthy patient. 
#It has the advantage of being a small but complex dataset, which will speed things up

#Task 2 - Part 2 - loading the Seurat and ggplot2 libraries 
install.packages("Seurat")
library(Seurat)
library(ggplot2)

#Task 2 - Part 3 - Reading the data from the tsv files and h5 files 
sc.data=Read10X(data.dir="/Users/shounakkadam/Documents/Module C - Single cell and spatial transcriptomics/sc_data/filtered_feature_bc_matrix")
install.packages("hdf5r")
sc.data=Read10X_h5("/Users/shounakkadam/Documents/Module C - Single cell and spatial transcriptomics/sc_data/filtered_feature_bc_matrix.h5")

#Task 2 - Part 5
sc.data = CreateSeuratObject(sc.data)
sc.data$orig.ident <- "Control_R1"

#Task 2 - Part 6 - Pulling out key inforation from the object
meta.table = data.frame(sc.data@meta.data)
counts.table = data.frame(sc.data@assays$RNA$counts)

#Reason why there are so many 0s in the counts is because it's the per cell quality 
#of scRNA is always very poor. We only sequenced 40k reads per cell, or a mean of roughly one read per gene

#Task 2 - Part 7 - Making plots from the meta-data using ggplot2
ggplot(meta.table, aes(x=nCount_RNA))+geom_histogram() #for number of reads per cell
ggplot(meta.table, aes(x=nFeature_RNA))+geom_histogram()#for number of genes with at least one read 

ggplot(counts.table, aes(x=AAACATACAACCAC.1))+geom_histogram()#number of reads per gene for the first individual cell in the counts table

#Task 3 - Quality Control - Part 1 and 2
VlnPlot(sc.data,features=c("nCount_RNA"),pt.size=0)+NoLegend()
VlnPlot(sc.data,features=c("nCount_RNA","nFeature_RNA"),pt.size=0)+NoLegend()

FeatureScatter(sc.data,feature1="nCount_RNA",feature2="nFeature_RNA")

#Task 3 - Part 4 - need to identify dead or dying cells.
# these cells can be detected based on having abnormally high% of reads at mitochondrial genes compared to all other genes 

sc.data=PercentageFeatureSet(sc.data,pattern="^MT-",col.name="percent_mito")

#Task 3 - Part 5 
sc.data=PercentageFeatureSet(sc.data,pattern="^ZNF",col.name="percent_znf")

my_genes=c("MALAT1", "ACTB")
sc.data=PercentageFeatureSet(sc.data,features= my_genes,col.name="percent_my_genes")

#Task 3 - Part 6
#To assess the distribution of dead or dying cells, so we can define suitable cut-offs
VlnPlot(sc.data,features=c("nCount_RNA","nFeature_RNA","percent_mito"),pt.size=0)+NoLegend()

#Task 3 - Part 7 
#highlighting the relationship between cell death and reads captured
#high% mt expression have low features and counts overall

FeatureScatter(sc.data, feature1 = "nCount_RNA", feature2 = "percent_mito")
FeatureScatter(sc.data, feature1 = "nFeature_RNA", feature2 = "percent_mito")

#Task 3 - Part 8
#removing unwanted cells 
n_before <- ncol(sc.data)
sc.data = subset(sc.data,subset=nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mito < 5)
n_after <- ncol(sc.data)
c(before = n_before, after = n_after, removed = n_before - n_after)

VlnPlot(sc.data, features = c("nCount_RNA","nFeature_RNA","percent_mito"), ncol = 3, pt.size = 0.1)
FeatureScatter(sc.data, feature1 = "nCount_RNA",   feature2 = "percent_mito")
FeatureScatter(sc.data, feature1 = "nFeature_RNA", feature2 = "percent_mito")

