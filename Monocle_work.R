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
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree',auto_param_selection = F)
cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "anno")
plot_cell_trajectory(cds, color_by = "seurat_clusters")
plot_cell_trajectory(cds, color_by = "Pseudotime")
#orderCells(HSMM_myo, root_state = c(a,b,c))
#plot_cell_trajectory(cds, color_by = "Pseudotime")

BEAM_res <- BEAM(cds, branch_point = 1, cores = 1)

plot_genes_branched_heatmap(lung[row.names(subset(BEAM_res,qval < 1e-4)),],branch_point = 1,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
#plot_genes_branched_heatmap(cds[c('Prdx1','Prdx4',"Stmn1",'Top2a','Mki67','Cd19','Jchain',"Ighd",'Ms4a1','Cd22',"Ighg1",'Sell','Il7r','Tcf7','Igkc','Ighm','Hmgb2'),],branch_point = 1,num_clusters = 4,cores = 1,use_gene_short_name = T,show_rownames = T)
