library(Seurat)
library(tidyverse)
library(magrittr)
library(monocle)

combined <- readRDS('data/Demo_CombinedSeurat_SCT_Preprocess_FilterLQCells.rds')
Idents(combined) <- 'anno'
DefaultAssay(combined) <- 'RNA'
combined <- NormalizeData(combined)

cds <- as.CellDataSet(combined)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~anno")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
