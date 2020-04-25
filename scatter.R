ARGS<- commandArgs(TRUE)

library("scater")
ms=as.matrix(read.table(ARGS[1],header=T,row.names=1))
ms=t(ms)
m=ms[3:nrow(ms)-1,]

sce=SingleCellExperiment(assays = list(counts = m))
per.cell <- perCellQCMetrics(sce)
colData(sce) <- cbind(colData(sce), per.cell)
sce<-logNormCounts(sce)

pdf("RNA_sum_vs_detected.pdf")
plotColData(sce, x = "sum", y="detected")
dev.off()
pdf("RNA_highest_count.pdf")
plotHighestExprs(sce)
dev.off()
pdf("RNA_count_20_genes.pdf")
plotExpression(sce, rownames(sce)[1:20])
dev.off()
pdf("RNA_count_scatter_plot.pdf")
plotExpression(sce, rownames(sce)[1:20],x = rownames(sce)[20])
dev.off()
