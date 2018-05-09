### Same as Analysis11 but use top genes w/o zeros instead of Seurat variable

library(Seurat)
library(dplyr)
library(stringr)
library(liger)
library(GSA)
library(pheatmap)
library(parallel)
library(matrixStats)
library(RColorBrewer)

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
#sc <- ScaleData(object=sc)

## Load cell cycle data
cc.genes <- readLines(con="satija_cell_cycle/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

## Assign cell cycle scores
sc <- CellCycleScoring(object=sc, s.genes=s.genes, g2m.genes=g2m.genes, set.ident=T)

## Regress out cell cycle
sc <- ScaleData(object=sc, vars.to.regress=c("S.Score", "G2M.Score"))

## Build a dataframe of barcodes and experiment name
barcode.src <- data.frame(sc@raw.data@Dimnames[[2]])
names(barcode.src) <- c("src")
barcode.src$src <- lapply(as.character(barcode.src$src), get.pop.name)
rownames(barcode.src) <- sc@raw.data@Dimnames[[2]]

## Convert to matrix preserving rows with < n% zeros
n.zero.lim = 10
prevalent.genes <- as.matrix(sc@scale.data[rowSums(sc@data==0)<=n.zero.lim,])

## Find Seurat Clusters
sc <- RunPCA(object=sc, pc.genes=rownames(prevalent.genes), do.print=F)
sc <- FindClusters(object=sc, reduction.type="pca", dims.use=1:10,
                   resolution=0.6, print.output=0, save.SNN=T)

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
short.14.clusters <- full.14.clusters[which(full.14.clusters$Barcode %in% colnames(prevalent.genes)),]

## Combine the cells, clusters, and experiment IDs into one dataframe
clusts.exps <- merge(barcode.src, short.14.clusters, by.x="row.names", by.y="Barcode")

## Sort the combined dataframe by expr.mat.var order
clusts.exps <- clusts.exps[match(colnames(prevalent.genes), clusts.exps$Row.names),]

## Build the pheatmap annotation column dataframe
anno.col <- data.frame(
  Experiment=factor(unlist(clusts.exps$src)),
  Orig.Cluster=factor(unlist(lapply(clusts.exps$Cluster, get.cluster.name))),
  Cell.Cycle=sc@meta.data$Phase,
  Seurat.Cluster=sc@ident
)
rownames(anno.col) <- colnames(prevalent.genes)
all.colors <- c(
  "#c8d24d",
  "#8a49c7",
  "#6dcb60",
  "#c64d91",
  "#87cab7",
  "#d1503b",
  "#7388c5",
  "#c69a56",
  "#442f59",
  "#4b643e",
  "#cc9fb0",
  "#723a33")
hm.anno.colors <- list(
  Experiment=all.colors,
  Orig.Cluster=all.colors[1:7],
  Cell.Cycle=all.colors[1:3],
  Seurat.Cluster=all.colors[1:4]
)
names(hm.anno.colors$Experiment) <- attr(anno.col$Experiment, "levels")
names(hm.anno.colors$Orig.Cluster) <- attr(anno.col$Orig.Cluster, "levels")
names(hm.anno.colors$Cell.Cycle) <- attr(anno.col$Cell.Cycle, "levels")
names(hm.anno.colors$Seurat.Cluster) <- attr(anno.col$Seurat.Cluster, "levels")
pheatmap(prevalent.genes, show_rownames=F, show_colnames=F, 
         cluster_rows=F, annotation_col=anno.col,
         annotation_colors=hm.anno.colors)
