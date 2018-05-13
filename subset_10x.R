### This script subsets the 10X output 
### NOTE: Because I'm still learning R, after running this script, you need to
### open the output matrix and add these lines:
###
### %%MatrixMarket matrix coordinate integer general
### %

## Input paths for matrix, genes, and barcodes files, and the number of cells to sample
input.mtx.path <- "10x_c14agg/matrix.mtx"
input.barcodes.path <- "10x_c14agg/barcodes.tsv"
input.genes.path <- "10x_c14agg/genes.tsv"
output.mtx.path <- "10x_c14agg_500/matrix.mtx"
output.barcodes.path <- "10x_c14agg_500/barcodes.tsv"
output.genes.path <- "10x_c14agg_500/genes.tsv"
num.cells <- 500

## Load the input matrix, extract number of cells in the matrix total
input.mtx <- read.csv(input.mtx.path, skip=2, sep=' ')
total.cells <- as.numeric(substring(colnames(input.mtx)[2], 2))

## Rename the input matrix columns for convenience
names(input.mtx) <- c("row", "col", "val")

## Load the input barcodes
input.barcodes <- readLines(input.barcodes.path)

## Load the input genes
input.genes <- read.csv(input.genes.path, sep='\t', header=FALSE)

## Randomly choose num.cells cells to sample from teh matrix
rand.col.indx <- sort(sample(total.cells, num.cells), decreasing=FALSE)

## Subset the matrix, barcodes, and genes
## barcodes are chosen by their index from the random sampling
## genes are chosen by the maximum row number of original genes with non-zero
## values remaining after sampling
## Due to Seurat requirements, cells are renamed to be 1:num.cells
mtx.subset <- input.mtx[input.mtx$col %in% rand.col.indx,]
mtx.subset$col <- match(mtx.subset$col, rand.col.indx)
barcodes.subset <- input.barcodes[rand.col.indx]

## Compute the number of genes and non-zero entries in the new matrix
num.genes <- max(mtx.subset$row)
num.entries <- length(mtx.subset[,1])
genes.subset <- input.genes[1:num.genes,]

## Add the genes, number of cells (columns) and total non-zero entries
mtx.subset <- rbind(c(num.genes, num.cells, num.entries), mtx.subset)

## Write the resulting lists/dataframes
write.table(mtx.subset, output.mtx.path, sep=' ', quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(barcodes.subset, output.barcodes.path, quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(genes.subset, output.genes.path, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
