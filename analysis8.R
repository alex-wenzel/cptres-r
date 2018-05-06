### This analysis loads the 28k cell data into seurat as before,
### selects a random subset of 5000 cells, and does a pheatmap
### clustering without regard for previous graphclust results.
###
### This analysis should also introduce better GSEA code (add 
### ES/NES, run in bulk, etc)

library(Seurat)
library(dplyr)
library(stringr)
library(liger)
library(GSA)
library(pheatmap)
library(parallel)

scdata <- Read10X(data.dir="subset_500/")
sc <- CreateSeuratObject(raw.data=scdata, min.cells=3, min.genes=3, project="10X_CPTRES")

exp.names <- c("r06_5", "r18_5", "s02", "s03", "r14_5", "r06_15",
               "r16_15", "r14_15", "r16_5", "caov3", "r18_15", "parent")
names(exp.names) <- 1:12

get.pop.name <- function(bc) {
  exp.names[as.numeric(strsplit(bc, "-")[[1]][2])]
}

barcode.src <- data.frame(sc@raw.data@Dimnames[[2]])
names(barcode.src) <- c("src")
barcode.src$src <- lapply(as.character(barcode.src$src), get.pop.name)
rownames(barcode.src) <- sc@raw.data@Dimnames[[2]]

sc <- NormalizeData(object=sc, normalization.method="LogNormalize", scale.factor=10000)

## This needs a lot of memory
sc <- ScaleData(object=sc)

## Convert to matrix preserving rows with < n% zeros
n.zero.lim = 0.001*sc@data@Dim[1]
prevalent.genes <- as.matrix(sc@data[rowSums(sc@data==0)<=n.zero.lim,])

## Cell Clusters
cluster.names <- c("Large", "Small-1", "Medium-1", "Large", "Small-3",
                   "Large", "Medium-1", "Medium-2", "Small-2", "Medium-2",
                   "Small-4", "Medium-1", "Small-2", "Medium-2")
names(cluster.names) <- 1:14

get.cluster.name <- function(num) {
  cluster.names[num]
}

full.14.clusters <- read.csv("full_28k/orig_graphclust_14/clusters.csv")
short.14.clusters <- full.14.clusters[which(full.14.clusters$Barcode %in% colnames(prevalent.genes)),]

clusts.exps <- merge(barcode.src, short.14.clusters, by.x="row.names", by.y="Barcode")
clusts.exps <- clusts.exps[match(colnames(prevalent.genes), clusts.exps$Row.names),]

anno.col <- data.frame(
  Experiment = factor(unlist(clusts.exps$src)),
  Orig.Cluster = factor(unlist(lapply(clusts.exps$Cluster, get.cluster.name)))
)
rownames(anno.col) <- colnames(prevalent.genes)
pheatmap(prevalent.genes, show_rownames=F, show_colnames=F, 
         cluster_rows=T, annotation_col=anno.col)
