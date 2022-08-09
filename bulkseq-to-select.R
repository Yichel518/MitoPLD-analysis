BiocManager::install("GEOquery")
BiocManager::install("umap")
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE57046", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17795", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "00000111111100001111111XXXXXXXXXXXXXXXXXXXXXXXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("dndMO14","WT14"))
levels(gs) <- groups

gset$
        design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
cont.matrix
# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
View(fit2)
tT <- topTable(fit2, adjust="fdr", sort.by="B",number = 1000000)
?topTable
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC","SPOT_ID"))
tT <- subset(tT, adj.P.Val < 0.05 & (logFC > 0 |logFC < 0))
View(tT)
write.table(tT, file="/home/tuyx/Paper/Gonad_bulkseq/PGC_trunk_dndvswt14.txt", row.names=F, sep="\t")
deg_dnd_up <- subset(tT, logFC > 0)
deg_wt_up <- subset(tT, logFC < 0)
BP_dnd_up <- enrichGO(gene=deg_dnd_up$ID,keyType = "SYMBOL",OrgDb= org.Dr.eg.db, ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.05,readable=F)
BP_dnd_up <- simplify(BP_dnd_up)
BP_dnd_up_table=as.data.frame(BP_dnd_up)
View(BP_dnd_up_table)
BP_wt_up <- enrichGO(gene=deg_wt_up$ID,keyType = "SYMBOL",OrgDb= org.Dr.eg.db, ont="BP",
                     pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1,readable=F)
BP_wt_up <- simplify(BP_wt_up)
BP_wt_up_table=as.data.frame(BP_wt_up)
View(BP_wt_up_table)
##########################22#################
# load series and platform data from GEO

gset <- getGEO("GSE57046", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17795", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "XXXXXXXXXXXXXXXXXXXXXXX000011111111000000111111"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("dndMO22","WT22"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT2 <- topTable(fit2, adjust="fdr", sort.by="B",number = 1000000)
?topTable
tT2 <- subset(tT2, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC","SPOT_ID"))
tT2 <- subset(tT2, adj.P.Val < 0.05 & (logFC > 0 |logFC < 0))
View(tT2)
write.table(tT2, file="/home/tuyx/Paper/Gonad_bulkseq/PGC_trunk_dndvswt22.txt", row.names=F, sep="\t")
deg22_dnd_up <- subset(tT2, logFC > 0)
deg22_wt_up <- subset(tT2, logFC < 0)
BP22_dnd_up <- enrichGO(gene=deg22_dnd_up$ID,keyType = "SYMBOL",OrgDb= org.Dr.eg.db, ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.05,readable=F)
BP22_dnd_up <- simplify(BP22_dnd_up)
BP22_dnd_up_table=as.data.frame(BP22_dnd_up)
View(BP_dnd_up_table)
BP22_wt_up <- enrichGO(gene=deg22_wt_up$ID,keyType = "SYMBOL",OrgDb= org.Dr.eg.db, ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.05,readable=F)
BP22_wt_up <- simplify(BP22_wt_up)
BP22_wt_up_table=as.data.frame(BP22_wt_up)
View(BP22_wt_up_table)




#######################################ALL
gset <- getGEO("GSE57046", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL17795", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "00000111111100001111111222233333333222222333333"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
View(ex)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

expL <- exprs(gset)
colnames(expL) <- c("dndMO14_1","dndMO14_2","dndMO14_3","dndMO14_4","dndMO14_5",
                    "WT14_1","WT14_2","WT14_3","WT14_4","WT14_5","WT14_5","WT14_7",
                    "dndMO14_6","dndMO14_7","dndMO14_8","dndMO14_9",
                    "WT14_8","WT14_9","WT14_10","WT14_11","WT14_12","WT14_13","WT14_14",
                    "dndMO22_1","dndMO22_2","dndMO22_3","dndMO22_4",
                    "WT22_1","WT22_2","WT22_3","WT22_4","WT22_5","WT22_5","WT22_7","WT22_8",
                    "dndMO22_5","dndMO22_6","dndMO22_7","dndMO22_8","dndMO22_9","dndMO22_10",
                    "WT22_9","WT22_10","WT22_11","WT22_12","WT22_13","WT22_14")
name <- c(rep("MO_14dpf",5),rep("WT_14dpf",7),rep("MO_14dpf",4),rep("WT_14dpf",7),
          rep("MO_22dpf",4),rep("WT_22dpf",8),rep("MO_22dpf",6),rep("WT_22dpf",6))
dim(expL)
m.vars=apply(expL,1,var)
expL<- expL[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4]),]
meanL=as.data.frame(t(apply(expL,1,function(a){
        tapply(a,name,mean)})))
n <- t(scale(t(meanL)))
n[n>2]=2
n[n<-2]=-2
head(n)
dim(n)
write.table(n,"/home/tuyx/Paper/Gonad_bulkseq/public_14_22_scale.txt",col.names=T,row.names=T,quote=F,sep="\t")

set.seed(200)
p=pheatmap::pheatmap(n,kmeans_k = 6,border_color = "white",fontsize = 18,cluster_cols = F,
                     clustering_method = "ward.D2",
                     color = colorRampPalette(c("#336699", "white", "#CC3333"))(50) )
saveRDS(p,"kmeans6_pgc_plot.rds")
a=as.data.frame(p$kmeans$cluster)
a$gene <- rownames(a)
View(a)
write.table(a,"cluster_pgc_table.txt",col.names=T,row.names=T,quote=F,sep="\t")
cluster3.1=a[a$`p$kmeans$cluster`=="3",]
BP_cluster3.1 <- enrichGO(gene=rownames(cluster3.1),keyType = "SYMBOL",OrgDb= org.Dr.eg.db, ont="BP",
                        pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.05,readable=F)
BP_cluster3.1 <- clusterProfiler::simplify(BP_cluster3.1)
BP_cluster3.1_table=as.data.frame(BP_cluster3.1)
View(BP_cluster3.1_table)
ggplot(arrange(BP_cluster3.1_table,desc(Count)),
       aes(Count / as.numeric(sub("/\\d+", "", BgRatio)),reorder(Description,Count))) + 
        geom_point(aes(color=qvalue,size=Count))+
        scale_colour_gradient(low="red",high= "gold")+
        labs(size="Count",x="RichFactor",y="GO Term",title="Biological processes for 
WT differentiation")+
        theme_classic()+
        theme(title =element_text(size = 24),
              axis.text =element_text(size = 24, color = 'black'))+
        theme(axis.text.x = element_text(size = 24))+
        theme(axis.ticks = element_blank(),legend.text = element_text(size=16),legend.title = element_text(size=18))


cluster5.1=a[a$`p$kmeans$cluster`=="5",]
BP_cluster5.1 <- enrichGO(gene=rownames(cluster5.1),keyType = "SYMBOL",OrgDb= org.Dr.eg.db, 
                          ont="BP",pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1,readable=F)
BP_cluster5.1 <- clusterProfiler::simplify(BP_cluster5.1)
BP_cluster5.1_table=as.data.frame(BP_cluster5.1)%>%subset(.,BP_cluster5.1$pvalue<0.05)
View(BP_cluster5.1_table)

ggplot(arrange(BP_cluster5.1_table,desc(pvalue))[1:20,],
       aes(Count / as.numeric(sub("/\\d+", "", BgRatio)),reorder(Description,Count))) + 
        geom_point(aes(color=pvalue,size=Count))+
        scale_colour_gradient(low="red",high= "gold")+
        labs(size="Count",x="RichFactor",y="GO Term",title="Biological processes for 
MO differentiation")+
        theme_classic()+
        theme(title =element_text(size = 24),
              axis.text =element_text(size = 24, color = 'black'))+
        theme(axis.text.x = element_text(size = 24))+
        theme(axis.ticks = element_blank(),legend.text = element_text(size=16),
              legend.title = element_text(size=18))

write.table(BP_cluster5.1_table,"BP_dnd_C5_table.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(BP_cluster3.1_table,"BP_wt_C3_table.txt",col.names=T,row.names=T,quote=F,sep="\t")
