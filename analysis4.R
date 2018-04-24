### This script inspired by http://satijalab.org/seurat/pbmc3k_tutorial.html

library(Seurat)
library(dplyr)
library(stringr)

## This function implements cell-cycle scoring
## should be run after scaling/variable genes

scdata <- Read10X(data.dir="subset_5000/")
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

## cell cycle regression
## inspired by http://satijalab.org/seurat/cell_cycle_vignette.html

cc.genes <- readLines(con="satija_cell_cycle/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]  #S phase genes
g2m.genes <- cc.genes[44:97]

## linear dimensional reduction
sc <- RunPCA(object=sc, pc.genes=sc@var.genes)
PCAPlot(object=sc, dim.1=1, dim.2=2)
## commented heatmap because it takes forever
#PCHeatmap(object=sc, pc.use=1, do.balanced=TRUE, label.columns=FALSE)

## Clustering
sc <- FindClusters(object=sc, reduction.type="pca", dims.use=1:10,
                   resolution=0.6, print.output=0, save.SNN=TRUE)

sc <- RunTSNE(object=sc, dims.use=1:10, do.fast=TRUE)
TSNEPlot(object=sc)
TSNEPlot(object=sc, group.by='src')

## Cluster markers (DEG)
sc.markers <- FindAllMarkers(object=sc, min.pct=0.25, thresh.use=0.25)
top10 <- sc.markers %>% group_by(cluster) %>% top_n(10,avg_logFC)
DoHeatmap(object=sc, genes.use=top10$gene, slim.col.label=TRUE, remove.key=TRUE)