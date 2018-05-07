### Version of analysis8.R that uses the top n most variable genes
### instead of the ones with fewest 0s

library(Seurat)
library(dplyr)
library(stringr)
library(liger)
library(GSA)
library(pheatmap)
library(parallel)
library(matrixStats)

## Mapping of barcode suffix number to experiment name
exp.names <- c("r06_5", "r18_5", "s02", "s03", "r14_5", "r06_15",
              "r16_15", "r14_15", "r16_5", "caov3", "r18_15", "parent")
names(exp.names) <- 1:12

## Uses the above list to map a barcode to its experiment
get.pop.name <- function(bc) {
  exp.names[as.numeric(strsplit(bc, "-")[[1]][2])]
}

## Load 10X data w/Seurat
scdata <- Read10X(data.dir="subset_500/")
sc <- CreateSeuratObject(raw.data=scdata, min.cells=3, min.genes=3, project="10X_CPTRES")

## Normalize and scale 10X data w/Seurat
sc <- NormalizeData(object=sc, normalization.method="LogNormalize", scale.factor=10000)
sc <- ScaleData(object=sc)

## Build a dataframe of barcodes and experiment name
barcode.src <- data.frame(sc@raw.data@Dimnames[[2]])
names(barcode.src) <- c("src")
barcode.src$src <- lapply(as.character(barcode.src$src), get.pop.name)
rownames(barcode.src) <- sc@raw.data@Dimnames[[2]]

## Convert sc@data to a basic matrix
expr.mat <- as.matrix(sc@data)

## Build a dataframe of genes and standard deviations across cells
genes.stdev <- data.frame(gene = attr(expr.mat, "dimnames")[[1]],
                          stdev = rowSds(expr.mat))
genes.stdev <- genes.stdev[order(-genes.stdev$stdev),]

## Get the top n variable genes
top.var.genes <- head(as.character(genes.stdev$gene), 200)

## Subset expr.mat to only contain top variable genes
expr.mat.var <- expr.mat[top.var.genes,]

## Mapping of original 14 clusters to consolidated version
cluster.names <- c("Large", "Small-1", "Medium-1", "Large", "Small-3",
                   "Large", "Medium-1", "Medium-2", "Small-2", "Medium-2",
                   "Small-4", "Medium-1", "Small-2", "Medium-2")
names(cluster.names) <- 1:14

## Function to map cluster number to combined group
get.cluster.name <- function(num) {
  cluster.names[num]
}

## Load original 14 clusters
full.14.clusters <- read.csv("full_28k/orig_graphclust_14/clusters.csv")

## Get the subset of cells from the clusters that are in the sample
short.14.clusters <- full.14.clusters[which(full.14.clusters$Barcode %in% colnames(expr.mat.var)),]

## Combine the cells, clusters, and experiment IDs into one dataframe
clusts.exps <- merge(barcode.src, short.14.clusters, by.x="row.names", by.y="Barcode")

## Sort the combined dataframe by expr.mat.var order
clusts.exps <- clusts.exps[match(colnames(expr.mat.var), clusts.exps$Row.names),]

## Build the pheatmap annotation column dataframe
anno.col <- data.frame(
  Experiment=factor(unlist(clusts.exps$src)),
  Orig.Cluster=factor(unlist(lapply(clusts.exps$Cluster, get.cluster.name)))
)
rownames(anno.col) <- colnames(expr.mat.var)
pheatmap(expr.mat.var, show_rownames=F, show_colnames=F, 
         cluster_rows=F, annotation_col=anno.col)
