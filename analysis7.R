### Analysis 7
###
### This script starts with the k=14 clustering and builds
### a heatmap and dendrogram based on Seurat's FindAllMarkers

library(Seurat)
library(dplyr)
library(stringr)
library(liger)
library(GSA)
library(pheatmap)

scdata <- Read10X(data.dir="full_28k/raw28k/")
sc <- CreateSeuratObject(raw.data=scdata, min.cells=3, min.genes=3, project="10X_CPTRES")

get.pop.name <- function(bc) {
  strsplit(bc, "-")[[1]][2]
}

barcode.src <- data.frame(sc@raw.data@Dimnames[[2]])
names(barcode.src) <- c("src")
barcode.src$src <- as.numeric(lapply(as.character(barcode.src$src), get.pop.name))
rownames(barcode.src) <- sc@raw.data@Dimnames[[2]]
sc <- AddMetaData(object=sc, metadata=barcode.src, col.name="barcode.src")

sc <- NormalizeData(object=sc, normalization.method="LogNormalize", scale.factor=10000)

## Swaps out sc@ident to the external clustering
new.ident <- as.factor(read.csv("full_28k/orig_graphclust_14/clusters.csv")$Cluster)
attributes(new.ident)$names <- attributes(sc@ident)$names
sc@ident <- new.ident

## Filter the sc@data matrix to only contain at most
## n.zero.lim fraction of 0s

n.zero.lim = 0.10*sc@data@Dim[1]
prevalent.genes <- as.matrix(sc@data[rowSums(sc@data==0)<=n.zero.lim,])

## Build a matrix where rows are genes and columns are
## clusters by averaging cell level expression in prevelant.genes.mat
## by cluster

cluster.ids <- data.frame(as.numeric(sc@ident))
cluster.ids$cell <- sc@data@Dimnames[[2]]
names(cluster.ids) <- c("cluster", "cell")
cluster.expr <- matrix(nrow=nrow(prevalent.genes), ncol=max(cluster.ids$cluster))
rownames(cluster.expr) <- attr(prevalent.genes, "dimnames")[[1]]
colnames(cluster.expr) <- 1:14
for (i in 1:nrow(prevalent.genes)) {
  for (j in 1:max(cluster.ids$cluster)) {
    cluster.expr[i,j] <- mean(prevalent.genes[i,cluster.ids[cluster.ids$cluster==j,]$cell])
  }
}

pheatmap(cluster.expr, show_rownames=F, show_colnames=T)
