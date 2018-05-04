### This script inspired by http://satijalab.org/seurat/pbmc3k_tutorial.html

library(Seurat)
library(dplyr)

scdata <- Read10X(data.dir="subset_test/")
sc <- CreateSeuratObject(raw.data=scdata, min.cells=3, min.genes=3, project="10X_CPTRES")

mito.genes <- grep(pattern="^MT-", x=rownames(x=sc@data), value=TRUE)
percent.mito <- Matrix::colSums(sc@raw.data[mito.genes, ])/Matrix::colSums(sc@raw.data)

sc <- AddMetaData(object=sc, metadata=percent.mito, col.name="percent.mito")
VlnPlot(object=sc, features.plot=c("nGene", "nUMI", "percent.mito"), nCol=3, point.size.use=0.5)

par(mfrow=c(1,2))
GenePlot(object=sc, gene1="nUMI", gene2="percent.mito")
GenePlot(object=sc, gene1="nUMI", gene2="nGene")

sc <- NormalizeData(object=sc, normalization.method="LogNormalize", scale.factor=10000)

par(mfrow=c(1,1))
sc <- FindVariableGenes(object=sc, mean.function=ExpMean, dispersion.function=LogVMR,
                        x.low.cutoff=0.0125, x.high.cutoff=3, y.cutoff=0.5)

## This needs a lot of memory
sc <- ScaleData(object=sc, vars.to.regress=c("nUMI", "percent.mito"))

## linear dimensional reduction
sc <- RunPCA(object=sc, pc.genes=sc@var.genes)
PCAPlot(object=sc, dim.1=1, dim.2=2)
PCHeatmap(object=sc, pc.use=1, do.balanced=TRUE, label.columns=FALSE)

## Clustering
sc <- FindClusters(object=sc, reduction.type="pca", dims.use=1:10,
                   resolution=0.6, print.output=0, save.SNN=TRUE)

sc <- RunTSNE(object=sc, dims.use=1:10, do.fast=TRUE)
TSNEPlot(object=sc)
