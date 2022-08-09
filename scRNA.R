################################ovary####################################################
setwd("/home/tuyx/scRNA_pld6/")
BiocManager::install("Seurat")
data1 <- read_csv("/home/tuyx/scRNA_pld6/obj_ssportal_Matrix.csv")
data1[1:5,1:5]
rownames(data1) <- data1$GENE

# csv矩阵转换成数据框
data1 = data.frame(data1)
rownames(data1) <- data1$GENE
data1$GENE <- NULL
# 数据框转换成稀疏矩阵matrix
dim(data1)
obj <- readRDS("obj.rds")
obj@reductions$umap@cell.embeddings %>%class()
dataan <- as(as.matrix(data1), "dgCMatrix")
dim(dataan)

obj <- CreateSeuratObject(counts = dataan, project = "test", min.cells = 3, min.features = 200)
cluster <- read.delim("zx124_40com_ssportal_cluster.txt",row.names = 1)
cluster <- cluster[-1,]
head(cluster)
cluster$X <- as.numeric(cluster$X)
cluster$Y <- as.numeric(cluster$Y)
meta <- read.delim("~/scRNA_pld6/zx124_40com_ssportal_meta.txt",row.names = 1)
meta <- meta[-1,]
head(meta)
View(obj@meta.data)
obj@meta.data$gene <- rownames(obj@meta.data)
meta$gene <- rownames(meta)
mmeta <- merge(obj@meta.data,meta[,4:5],all.x=T,by="gene")
obj@meta.data <- mmeta
dim(mmeta)
mmeta$gene
data2 <- data1[,colnames(data1)%in%mmeta$gene]
dataan <- as(as.matrix(data2), "dgCMatrix")
dim(dataan)
write.table(mmeta,"mmeta.txt",col.names=T,row.names=T,quote=F,sep="\t")
obj <- CreateSeuratObject(counts = dataan, project = "test", min.cells = 3, min.features = 200)

#rownames(mmeta) <- mmeta$gene
#mmeta$gene <- NULL
#obj@meta.data <- mmeta
obj <- NormalizeData(obj)
obj <- ScaleData(obj, display.progress = F)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

obj <- RunPCA(obj, features = VariableFeatures(object = obj))

obj <- JackStraw(obj, num.replicate = 100,dims = 50)
obj <- ScoreJackStraw(obj, dims = 1:50)
JackStrawPlot(obj, dims = 1:50)

obj <- FindNeighbors(obj, dims = 1:22)
obj <- FindClusters(obj, resolution = 0.5)
#obj@reductions$umap <- NULL
obj <- RunUMAP(obj, dims = 1:22, min.dist = 0.3)
DimPlot(obj, label = T,pt.size = 1.2,label.size = 5)
View(obj@meta.data)
mmeta
#obj$seurat_clusters <- factor(obj$active.ident)
View(obj@meta.data)
#obj@meta.data$seurat_clusters <- obj@meta.data$RNA_snn_res.0.5
View(obj@meta.data)

obj@reductions$umap
?DimPlot
gene=as.data.frame(obj@assays$RNA@counts)%>%rownames(.)
DimPlot(obj, label = T, pt.size = 1,
        cols=c("#2E8B57","#88ada6","#D8BFD8","#800000","#ffc773","#a1afc9",
        "#4682B4","#808080","#556B2F",
        "#5F9EA0","#DAA520","#BDB76B","#F4A460","#A0522D",
        "#C0C0C0","#FFB6C1","#4169E1","#008B8B","#FF6347"))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+ NoLegend()
  
DotPlot(object =obj, features = c("nanos2",'pld6',"ddx4","dmc1","sycp1","sycp2","sycp3","zp3","zar1",
                                  "mfn1b","mfn2","opa1",
                                  "dnm1l",
                                  "tfam","ppargc1b","pink1","map1lc3b","sqstm1"),
        cluster.idents = F)+scale_y_discrete(limits=c("GSC+GC_Pro","Early_Meio_OO","Meio","Late_Meio_OO",
                                                      "Early_OO_1","Early_OO_2","Early_OO_3",
                                                      "Stromal_cxcl12a","Stromal_scara3","Stromal_fgf24","Stromal_tbx2b",
                                                      "Follicle_2","Follicle_1","Follicle_lhx9","Theca","Macrophage","NK-like",
                                                      "Vasculature","Neutrophils"))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(size = 20,angle = 90, vjust = 0.5,hjust = 1))
View(obj@meta.data)

tt <- data.frame(cell=rownames(obj@meta.data),cluster=obj@meta.data[,5])
head(tt)
mmeta$cell <- rownames(mmeta)
mm <- mmeta[,4:5]
dim(mm)
dim(tt)
cluster_name <- merge(tt,mm,all.x=T,by="cell")
View(cluster_name)
new.cluster.ids <- c("Stromal_cxcl12a","GSC+GC_Pro","Follicle_2","GSC+GC_Pro","Follicle_1",
                     "Early_Meio_OO","Late_Meio_OO","NK-like","Early_OO_3","Stromal_scara3",
                     "Vasculature","Theca","Early_OO_2","Stromal_scara3","Early_OO_1",
                     "Macrophage","Follicle_lhx9","Meio","Stromal_fgf24","Stromal_tbx2b",
                     "Neutrophils","NK-like","Neutrophils")
names(x=new.cluster.ids)=levels(x = obj)
obj<- RenameIdents(object = obj, new.cluster.ids)

df1=as.data.frame(obj@active.ident)
obj@meta.data$celltype=df1$`obj@active.ident`
unique(obj$celltype)

FeaturePlot(object =obj, c("nanos2",'pld6',"dmc1","ddx4","zp3","zar1","sycp1","sycp2","sycp3"), cols = c("grey90", "blue"), pt.size = 0.5)
FeaturePlot(object =obj, c("gk5"), cols = c("grey90", "blue"), pt.size = 0.5)

saveRDS(obj,"obj.rds")
obj_new.markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(obj_new.markers,"obj_new.markers.txt",col.names=T,row.names=T,quote=F,sep="\t")
top50_obj.new=obj_new.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) #将ovar.markers传给group_by，group_by按cluster 排序，再将结果传给top_n，top_n按avg_logFC排序，显示每个类中的前两个
write.table(top50_obj.new,"top50_obj.txt",col.names=T,row.names=T,quote=F,sep="\t")
dim(obj_new.markers)
dim(top50_obj.new)
obj <- readRDS("obj.rds")

FeaturePlot(object =obj, c('pld6'), cols = c("grey90", "red"), pt.size = 0.5)+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(size = 20))+theme(
    plot.title = element_text(family = "Arial",color="black", size=24, face="italic"))+theme(
      plot.title = element_text(family = "Arial",color="black", size=24, face="italic"))
FeaturePlot(object =obj, c('gk5'), cols = c("grey90", "red"), pt.size = 0.5)+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(size = 20))+theme(
    plot.title = element_text(family = "Arial",color="black", size=24, face="italic"))+theme(
      plot.title = element_text(family = "Arial",color="black", size=24, face="italic"))

avgData_top50 <- as.data.frame(t(apply(obj@assays$RNA@data[top50_obj.new$gene,],1,function(x){
    tapply(x, obj$celltype, mean) # ExpMean
  })))
head(avgData_top50)
avgData <- as.data.frame(t(apply(obj@assays$RNA@data[obj_new.markers$gene,] ,1,function(x){
    tapply(x, obj$celltype, mean) # ExpMean
  })))
write.table(avgData_top50,"avgData_top50.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(avgData,"avgData.txt",col.names=T,row.names=T,quote=F,sep="\t")
DotPlot(object =obj, features = l,
        cluster.idents = F)+scale_y_discrete(limits=c("GSC+GC_Pro","Early_Meio_OO","Meio","Late_Meio_OO",
                                                      "Early_OO_1","Early_OO_2","Early_OO_3",
                                                      "Stromal_cxcl12a","Stromal_scara3","Stromal_fgf24","Stromal_tbx2b",
                                                      "Follicle_2","Follicle_1","Follicle_lhx9","Theca","Macrophage","NK-like",
                                                      "Vasculature","Neutrophils"))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(size = 20,angle = 90, vjust = 0.5,hjust = 1))+theme(
    plot.title = element_text(family = "Times",color="black", size=14, face="italic"))


germ_cell <-subset(x = obj,idents=c("GSC+GC_Pro","Early_Meio_OO","Meio","Late_Meio_OO",
                                    "Early_OO_1","Early_OO_2","Early_OO_3"))
DotPlot(object =germ_cell, features =  c("nanos2","nanos3",'pld6',"dmc1","ddx4","sycp1","sycp2","sycp3","zp3","zar1"),
        cluster.idents = F)+scale_y_discrete(limits=c("GSC+GC_Pro","Early_Meio_OO","Meio","Late_Meio_OO",
                                                      "Early_OO_1","Early_OO_2","Early_OO_3"))+theme(axis.title =element_text(size = 20),
                                                                                                     axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(size = 20,angle = 90, vjust = 0.5,hjust = 1))
obj <- readRDS("/home/tuyx/scRNA_pld6/obj.rds")
diff_20_big <- read.delim("~/scRNA_pld6/diff_20_big.txt", row.names=1)
diff_20 <- read.delim("~/scRNA_pld6/diff_d20.txt", row.names=1)
diff_20_small <- subset(diff_20,log2FoldChange > 0)
sym <- AnnotationDbi::select(org.Dr.eg.db,keys=rownames(diff_20_small),columns=c("ENSEMBL","SYMBOL"),keytype="ENSEMBL")

DotPlot(object =obj, features =diff_20_big$Gene.name,
        cluster.idents = F)+scale_y_discrete(limits=c("GSC+GC_Pro","Early_Meio_OO","Meio","Late_Meio_OO",
                                                      "Early_OO_1","Early_OO_2","Early_OO_3",
                                                      "Stromal_cxcl12a","Stromal_scara3","Stromal_fgf24","Stromal_tbx2b",
                                                      "Follicle_2","Follicle_1","Follicle_lhx9","Theca","Macrophage","NK-like",
                                                      "Vasculature","Neutrophils"))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(size = 20,angle = 90, vjust = 0.5,hjust = 1))+theme(
    plot.title = element_text(family = "Times",color="black", size=14, face="italic"))

DotPlot(object =obj, features =sym$SYMBOL,
        cluster.idents = F)+scale_y_discrete(limits=c("GSC+GC_Pro","Early_Meio_OO","Meio","Late_Meio_OO",
                                                      "Early_OO_1","Early_OO_2","Early_OO_3",
                                                      "Stromal_cxcl12a","Stromal_scara3","Stromal_fgf24","Stromal_tbx2b",
                                                      "Follicle_2","Follicle_1","Follicle_lhx9","Theca","Macrophage","NK-like",
                                                      "Vasculature","Neutrophils"))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(size = 20,angle = 90, vjust = 0.5,hjust = 1))+theme(
    plot.title = element_text(family = "Times",color="black", size=14, face="italic"))

############################################ testis ###################################################
# 1.Setup the Seurat Object（创建Seurat对象）
library(dplyr)
library(Seurat)
rm(list = ls())
tes.data <- Read10X(data.dir = "~/scRNA_pld6/t1/")
testis <- CreateSeuratObject(counts = tes.data, project = "testis", min.cells = 3, 
                             min.features =200) # gene:13714,cells：2700

testis[["percent.mito"]] <- PercentageFeatureSet(object = testis, pattern = "^mt-") #线粒体基因以MT开头的选择出来
testis.mt=testis[["percent.mito"]]
head(x = testis@meta.data, 5)
VlnPlot(object = testis, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3) 
testis <- subset(x = testis, subset = nFeature_RNA >200& nFeature_RNA < 2000& percent.mito < 10)
testis <- NormalizeData(object = testis, normalization.method = "LogNormalize", scale.factor = 1e4)
testis <- FindVariableFeatures(object = testis,selection.method = 'vst', nfeatures = 2000)
all.genes <- rownames(x = testis)
testis <- ScaleData(object = testis, vars.to.regress =  c("nCount_RNA", "percent.mito"))
testis <- RunPCA(object = testis, features = VariableFeatures(object = testis))
testis <- JackStraw(object = testis, num.replicate = 100)
testis <- ScoreJackStraw(object = testis, dims = 1:20)
ElbowPlot(object = testis) #基于每一个解释的方差百分比（“ElbowPlot”函数）对主成分进行排序

testis <- FindNeighbors(object = testis, dims = 1:18)
testis <- FindClusters(object = testis, resolution =0.6)
testis <- RunUMAP(object = testis, dims = 1:18)
DimPlot(object = testis, reduction = 'umap',label = TRUE,pt.size = 1)+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'),legend.text = element_text(size=14))


# TSNE
testis <- RunTSNE(object = testis, dims = 1:18)
DimPlot(object = testis, reduction = 'tsne',label = TRUE )
DotPlot(object = testis, features = c("star","gsdf","dazl","pcna","nanos3","pld6",
                                      "sycp2","sycp3","tekt1"))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =90,hjust = 0.5,vjust = 0.5),legend.text = element_text(size=16))

new.cluster.ids <- c("spermatocyte-late","spermatocyte","spermatocyte-early",
                     "spermatids","spermatogonia-late","spermatogonia-late", "spermatids",
                     "spermatogonia","sertoli","leydig")
names(x=new.cluster.ids)=levels(x = testis)
testis <- RenameIdents(object = testis, new.cluster.ids)
FeaturePlot(object = testis, features = c("pld6"))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0),legend.text = element_text(size=16))

FeaturePlot(object = testis, features = c("sycp3"))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=16))

saveRDS(testis, file = "~/scRNA_pld6/testis1.rds")
testis1 <- readRDS("~/scRNA_pld6/testis1.rds")
FeaturePlot(object =testis1, c('pld6'), cols = c("grey90", "red"), pt.size = 0.5)+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(size = 20))+theme(
    plot.title = element_text(family = "Arial",color="black", size=24, face="italic"))+theme(
      plot.title = element_text(family = "Arial",color="black", size=24, face="italic"))
# 1.Setup the Seurat Object（创建Seurat对象）
library(dplyr)
library(Seurat)
rm(list = ls())
gc()
tes.data <- Read10X(data.dir = "~/scRNA_pld6/t2/")
testis <- CreateSeuratObject(counts = tes.data, project = "testis", min.cells = 3, 
                             min.features =200) # gene:13714,cells：2700

testis[["percent.mito"]] <- PercentageFeatureSet(object = testis, pattern = "^mt-") #线粒体基因以MT开头的选择出来
testis.mt=testis[["percent.mito"]]
head(x = testis@meta.data, 5)
VlnPlot(object = testis, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3) 
testis <- subset(x = testis, subset = nFeature_RNA >200& nFeature_RNA < 2000& percent.mito < 10)
testis <- NormalizeData(object = testis, normalization.method = "LogNormalize", scale.factor = 1e4)
testis <- FindVariableFeatures(object = testis,selection.method = 'vst', nfeatures = 2000)
all.genes <- rownames(x = testis)
testis <- ScaleData(object = testis, vars.to.regress =  c("nCount_RNA", "percent.mito"))
testis <- RunPCA(object = testis, features = VariableFeatures(object = testis))
testis <- JackStraw(object = testis, num.replicate = 100)
testis <- ScoreJackStraw(object = testis, dims = 1:20)
ElbowPlot(object = testis) #基于每一个解释的方差百分比（“ElbowPlot”函数）对主成分进行排序
testis <- FindNeighbors(object = testis, dims = 1:18)
testis <- FindClusters(object = testis, resolution =0.6)
testis <- RunUMAP(object = testis, dims = 1:18)
DimPlot(object = testis1, reduction = 'umap',label = TRUE,pt.size = 1)+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'),legend.text = element_text(size=14))+guides(color=FALSE)
# TSNE
testis <- RunTSNE(object = testis, dims = 1:18)
DimPlot(object = testis, reduction = 'tsne',label = TRUE )
DotPlot(object = testis, features = c("gsdf","star","dazl","pcna","ddx4","piwil1","pld6",
                                      "sycp2","sycp3","tekt1","mfn1b","mfn2","opa1",
                                      "dnm1l",
                                      "tfam","ppargc1b","pink1","map1lc3b","sqstm1"))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =90,hjust = 1,vjust=0.5),legend.text = element_text(size=16))+
  scale_y_discrete(limits=c("sertoli","leydig","spermatogonia","spermatogonia-late","spermatocyte-early",
                            "spermatocyte","spermatocyte-late","spermatids"))
new.cluster.ids <- c("spermatocyte-late","spermatocyte","spermatocyte-early",
                     "spermatids","spermatogonia-late","spermatogonia-late", "spermatids",
                     "spermatogonia","sertoli","leydig")
names(x=new.cluster.ids)=levels(x = testis)
testis <- RenameIdents(object = testis, new.cluster.ids)
FeaturePlot(object = testis, features = c("pld6"))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0),legend.text = element_text(size=16))
FeaturePlot(object = testis, features = c("sycp3"))+
  theme(axis.title =element_text(size = 20),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=16))
saveRDS(testis, file = "~/scRNA_pld6/testis.rds")


