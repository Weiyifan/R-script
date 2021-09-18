args<-commandArgs(TRUE)

res=read.table(args[1],header=F,row.names=1,sep="\t")

library("edgeR")
group=c(1,1,2,2)#no replicates(1,2,3,4,5...)
y=DGEList(counts=res, group=group)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=TRUE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)# with replicates
#y <- exactTest(y, dispersion=bcv^2)# with no replicates
ye<- exactTest(y)
#ye<- exactTest(y)
#ou=topTags(ye,n=nrow(keep))
out=ye$table


library("ggplot2")

out$sig=ifelse(out$PValue < 0.05 & abs(out$logFC) > 1, "significant","Non-sig")
outsig=out[which(out$sig=="significant"),]
write.table(outsig,file=args[2],sep="\t",quote=FALSE)
p<-ggplot(out,aes(logFC,-log10(PValue),colour=sig))+geom_point(size=0.65,alpha=0.5)
p+xlim(-5,5)+ylim(0,5)+scale_colour_manual(values=c("grey","red"))+geom_hline(yintercept=0.05,linetype='dashed',size=0.3)+ geom_vline(xintercept=1,linetype='dashed',size=0.3)+ geom_vline(xintercept=-1,linetype='dashed',size=0.3)+theme_classic()

name=paste(args[2],"pdf",sep=".")
ggsave(name,width=8,height=7)
#res$signi=ifelse(res$pvalue <0.05 & abs(res$correlation) > 0.2,"significant","Non-sig")
#p <- ggplot(res,aes(correlation,-log10(pvalue),colour=signi))+ geom_point(size=0.5,alpha=0.5)

#p+ylim(0,10)+scale_colour_manual(values=c("grey","red"))+theme_classic()
#name=paste(args[1],'pdf',sep=".")
#ggsave(name,width=8,height=7)
