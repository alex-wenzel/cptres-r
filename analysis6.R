### Same as analysis5 but with k=7

### This analysis uses the original clustering from 10X K=7 at
### /mnt/oncogxA/Projects/CPTRES/RNAExpression/10x/aggregates/10x_5_15/outs/analysis/clustering/kmeans_9_clusters

library(Seurat)
library(dplyr)
library(stringr)
library(liger)
library(GSA)

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

VlnPlot(object=sc, features.plot=c("nGene", "nUMI"), nCol=2, point.size.use=0.5)

GenePlot(object=sc, gene1="nUMI", gene2="nGene")

sc <- NormalizeData(object=sc, normalization.method="LogNormalize", scale.factor=10000)

sc <- FindVariableGenes(object=sc, mean.function=ExpMean, dispersion.function=LogVMR,
                        x.low.cutoff=0.0125, x.high.cutoff=3, y.cutoff=0.5)

## This needs a lot of memory
sc <- ScaleData(object=sc)

## Swaps out sc@ident to the external clustering
new.ident <- as.factor(read.csv("full_28k/orig_k7/clusters.csv")$Cluster)
attributes(new.ident)$names <- attributes(sc@ident)$names
sc@ident <- new.ident

## PCA
sc <- RunPCA(object=sc, pc.genes=sc@var.genes)
PCAPlot(object=sc, dim.1=1, dim.2=2)

## tSNE
sc <- RunTSNE(object=sc, dims.use=1:10, do.fast=TRUE)
TSNEPlot(object=sc)
TSNEPlot(object=sc, group.by='src')

## DEG
sc.markers <- FindAllMarkers(object=sc, min.pct=0.25, thresh.use=0.25)
top10 <- sc.markers %>% group_by(cluster) %>% top_n(10,avg_logFC)
DoHeatmap(object=sc, genes.use=top10$gene, slim.col.label=TRUE, remove.key=TRUE)

## GSEA
top50 <- sc.markers %>% group_by(cluster) %>% top_n(50,avg_logFC)
n.clusters <- max(as.numeric(top50$cluster))

hallmarks.genes <- GSA.read.gmt("msigdb/h.all.v6.1.symbols.gmt")
hallmarks.res <- read.csv(text="cluster,geneset.name,n.genes.set,pvalue")

for (i in 1:n.clusters) {
  cluster.genes <- top50[top50$cluster==i,]$avg_logFC
  names(cluster.genes) <- top50[top50$cluster==i,]$gene
  for (j in 1:length(hallmarks.genes$genesets)) {
    gs.name <- hallmarks.genes$geneset.names[j]
    gs <- hallmarks.genes$genesets[[j]]
    n.genes.set <- length(intersect(attributes(cluster.genes)$names, hallmarks.genes$genesets[[j]]))
    hallmarks.res[nrow(hallmarks.res)+1,] = c(i, gs.name, n.genes.set, gsea(cluster.genes, gs, plot=FALSE))
  }
}
hallmarks.res <- hallmarks.res[hallmarks.res$n.genes.set>2,]
hallmarks.res <- hallmarks.res[order(hallmarks.res$cluster, as.numeric(hallmarks.res$pvalue)),]
write.table(hallmarks.res, "k7_hallmarks_gsea.tsv", quote=FALSE, row.names=FALSE, sep='\t')

reactome.genes <- GSA.read.gmt("msigdb/c2.cp.reactome.v6.1.symbols.gmt")
reactome.res <- read.csv(text="cluster,geneset.name,n.genes.set,pvalue")

for (i in 1:n.clusters) {
  cluster.genes <- top50[top50$cluster==i,]$avg_logFC
  names(cluster.genes) <- top50[top50$cluster==i,]$gene
  for (j in 1:length(reactome.genes$genesets)) {
    gs.name <- reactome.genes$geneset.names[j]
    gs <- reactome.genes$genesets[[j]]
    n.genes.set <- length(intersect(attributes(cluster.genes)$names, reactome.genes$genesets[[j]]))
    reactome.res[nrow(reactome.res)+1,] = c(i, gs.name, n.genes.set, gsea(cluster.genes, gs, plot=FALSE))
  }
}
reactome.res <- reactome.res[reactome.res$n.genes.set>2,]
reactome.res <- reactome.res[order(reactome.res$cluster, as.numeric(reactome.res$pvalue)),]
write.table(reactome.res, "k7_reactome_gsea.tsv", quote=FALSE, row.names=FALSE, sep='\t')