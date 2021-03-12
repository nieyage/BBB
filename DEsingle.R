# 
# library(DESeq2)
# #Young<- BBB_EC.integrated@meta.data[BBB_EC.integrated@meta.data$orig.ident=="a_Young",]
# #Old<- BBB_EC.integrated@meta.data[BBB_EC.integrated@meta.data$orig.ident=="c_Old",]
# Young <-BBB_EC.integrated[,which(BBB_EC.integrated@meta.data$orig.ident=="a_Young")]
# Young_count <- data.frame(Young@assays$RNA@counts)
# write.table(Young_count,"Young_count.tsv")
# 
# Old <-BBB_EC.integrated[,which(BBB_EC.integrated@meta.data$orig.ident=="c_Old")]
# Old_count <- data.frame(Old@assays$RNA@counts)
# write.table(Old_count,"Old_count.tsv")


library("DEsingle")
a<-read.table("/public/home/yangjw28/projects/BBB/BBB_sc/YMO/Young_count.tsv",head=T)
b<-read.table("/public/home/yangjw28/projects/BBB/BBB_sc/YMO/Old_count.tsv",head=T)
counts=cbind(a[,2:length(a[1,])],b[,2:length(b[1,])])
group=factor(c(rep(1,length(a[1,])-1),rep(2,length(b[1,])-1)))
resutls<-DEsingle(counts=counts, group=group)
pdf("DE.pdf")
plot(log2(resutls$foldChange),-log10(resutls$FDR_LR2))
results.classified <- DEtype(results = resutls, threshold = 0.05)
dev.off()


