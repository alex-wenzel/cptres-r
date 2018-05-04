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
library(parallel)

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

## This needs a lot of memory
sc <- ScaleData(object=sc)

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

## Now replace the clusters with Olivier's combination
## Make combined cluster names using the following mappings
##
##  1: 146
##  3: 1237
##  4: 146
##  6: 146
##  7: 1237
##  8: 10148
##  9: 913
## 10: 10148
## 12: 1237
## 13: 913
## 14: 10148

combined.ident <- c(length(sc@ident))
for (i in 1:length(sc@ident)) {
  if (sc@ident[i]==1 || sc@ident[i]==4 || sc@ident[i]==6) {
    combined.ident[i] <- 146
  } else if (sc@ident[i]==12 || sc@ident[i]==3 || sc@ident[i]==7) {
    combined.ident[i] <- 1237
  } else if (sc@ident[i]==8 || sc@ident[i]==10 || sc@ident[i]==14) {
    combined.ident[i] <- 10148
  } else if (sc@ident[i]==9 || sc@ident[i]==13) {
    combined.ident[i] <-913 
  } else {
    combined.ident[i] <- sc@ident[i]
  }
}
combined.ident <- as.factor(combined.ident)
attributes(combined.ident)$names <- attributes(sc@ident)$names
sc@ident <- combined.ident

## DEG
sc.markers <- FindAllMarkers(object=sc, min.pct=0.25, thresh.use=0.25)
top10 <- sc.markers %>% group_by(cluster) %>% top_n(10,avg_logFC)
DoHeatmap(object=sc, genes.use=top10$gene, slim.col.label=TRUE, remove.key=TRUE)

## GSEA
n.clusters <- max(as.numeric(sc.markers$cluster))
n.cores <- detectCores()

hallmarks.genes <- GSA.read.gmt("msigdb/h.all.v6.1.symbols.gmt")
hallmarks.res <- read.csv(text="cluster,geneset.name,n.genes.set,pvalue")

cluster.names <- as.numeric(levels(combined.ident))

for (i in 1:n.clusters) {
  cluster.genes <- sc.markers[sc.markers$cluster==cluster.names[i],]$avg_logFC
  names(cluster.genes) <- sc.markers[sc.markers$cluster==cluster.names[i],]$gene
  for (j in 1:length(hallmarks.genes$genesets)) {
    gs.name <- hallmarks.genes$geneset.names[j]
    gs <- hallmarks.genes$genesets[[j]]
    n.genes.set <- length(intersect(attributes(cluster.genes)$names, hallmarks.genes$genesets[[j]]))
    hallmarks.res[nrow(hallmarks.res)+1,] = c(cluster.names[i], gs.name, n.genes.set, 
                                              gsea(cluster.genes, gs, plot=FALSE, mc.cores=n.cores))
  }
}
hallmarks.res <- hallmarks.res[hallmarks.res$n.genes.set>2,]
hallmarks.res <- hallmarks.res[order(hallmarks.res$cluster, as.numeric(hallmarks.res$pvalue)),]
write.table(hallmarks.res, "gsea_res/graphclust_comb_hallmarks_gsea_a7.tsv", quote=FALSE, row.names=FALSE, sep='\t')

reactome.genes <- GSA.read.gmt("msigdb/c2.cp.reactome.v6.1.symbols.gmt")
reactome.res <- read.csv(text="cluster,geneset.name,n.genes.set,pvalue")

for (i in 1:n.clusters) {
  cluster.genes <- sc.markers[sc.markers$cluster==cluster.names[i],]$avg_logFC
  names(cluster.genes) <- sc.markers[sc.markers$cluster==cluster.names[i],]$gene
  for (j in 1:length(reactome.genes$genesets)) {
    gs.name <- reactome.genes$geneset.names[j]
    gs <- reactome.genes$genesets[[j]]
    n.genes.set <- length(intersect(attributes(cluster.genes)$names, reactome.genes$genesets[[j]]))
    reactome.res[nrow(reactome.res)+1,] = c(cluster.names[i], gs.name, n.genes.set, 
                                            gsea(cluster.genes, gs, plot=FALSE, mc.cores=n.cores))
  }
}
reactome.res <- reactome.res[reactome.res$n.genes.set>2,]
reactome.res <- reactome.res[order(reactome.res$cluster, as.numeric(reactome.res$pvalue)),]
write.table(reactome.res, "gsea_res/graphclust_comb_reactome_gsea_a7.tsv", quote=FALSE, row.names=FALSE, sep='\t')
