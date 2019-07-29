rm(list = ls())

library(data.table)
library(Seurat)

library(plyr) 

RC = fread("count-gene.txt",header = T)


RC$V2[duplicated(RC$V2)]


RC.S = ddply(RC, "V2", numcolwise(sum)) 


rownames(RC.S) = RC.S$V2

RC.S = RC.S[,-1]

RC.S = as.matrix(RC.S)

culture.condition = ifelse(substr(colnames(RC.S),1,1) == "C", "Co-culture" ,"Only_embryo" ) 

cell.type = ifelse(substr(colnames(RC.S),1,2) == "Ch", "Endomaterial" ,"Embryonic_cells" ) 

library(stringr)

embryo.stage = 
  substr( colnames(RC.S), (str_locate(colnames(RC.S), fixed("."))[,1])-1 , (str_locate(colnames(RC.S), fixed("."))[,1])+1) 

embryo.stage[is.na(embryo.stage)] = "Endomaterial"


embryo.no = colnames(RC.S)
embryo.no[1:69] = substr(paste("x",colnames(RC.S)[1:69],sep=""),1,5)
embryo.no = substr(embryo.no,1,5)

meta.data = data.frame(culture.condition,cell.type,embryo.stage,embryo.no)

###################################
# readcount数据预处理完成，开始构建 Seurat object
###################################

mr <- CreateSeuratObject(raw.data = RC.S, min.cells = 30, min.genes = 500, 
                         project = "Human coculture")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = mr@data), value = TRUE)

percent.mito <- Matrix::colSums(mr@raw.data[mito.genes, ])/Matrix::colSums(mr@raw.data)

meta.data = cbind(meta.data,percent.mito)

mr <- AddMetaData(object = mr, metadata = meta.data, col.name = c("culture.condition","cell.type","embryo.stage","percent.mito") )

VlnPlot(object = mr, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


#VlnPlot(object = mr, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


par(mfrow = c(1, 2))
GenePlot(object = mr, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = mr, gene1 = "nUMI", gene2 = "nGene")
par(mfrow = c(1, 1))

##################################

mr <- FilterCells(object = mr, subset.names = c("nGene", "percent.mito","nUMI"), 
                  low.thresholds = c(7000, -Inf,0), high.thresholds = c(Inf, 0.3, 2000000))

VlnPlot(object = mr, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


mr <- NormalizeData(object = mr, normalization.method = "LogNormalize", 
                    scale.factor = 10000)


par(mfrow = c(1, 1))
mr <- FindVariableGenes(object = mr, mean.function = ExpMean, dispersion.function = LogVMR)

length(x = mr@var.genes)

date()
mr <- ScaleData(object = mr, vars.to.regress = c("nUMI", "percent.mito"))
date()

###############################

# 整体的Seurat object 构建完成

#############################

mr <- RunPCA(object = mr, pc.genes = mr@var.genes, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5,pcs.compute = 40)


PCElbowPlot(object = mr)

mr <- JackStraw(object = mr, num.replicate = 100, do.print = FALSE,num.pc = 40)

JackStrawPlot(object = mr, PCs = 1:40) # 这里，选择PC12

#mr <- RunTSNE(object = mr, dims.use = 1:12, do.fast = TRUE)

#TSNEPlot(object = mr,group.by = "embryo.stage")

## 全部细胞的 PCA

PCAPlot(object = mr, dim.1 = 1, dim.2 = 2,group.by = "embryo.stage")

table(mr@meta.data$embryo.stage)

table(mr@meta.data[ mr@meta.data$culture.condition =="Only_embryo" ,6:7])

table(mr@meta.data[ mr@meta.data$culture.condition =="Co-culture" ,6:7])

PCAPlot(object = mr, dim.1 = 1, dim.2 = 2,group.by = "culture.condition")

which(mr@dr$pca@cell.embeddings[,1]>50)


length(mr@var.genes)

dim(mr@scale.data)
dim(mr@raw.data)
################################################

# 只选择embryo cells进行分析

#mr@dr$pca@cell.embeddings

mr.emb = SubsetData(object = mr, cells.use = colnames(mr@scale.data)[mr@dr$pca@cell.embeddings[,1] < 25])

PCAPlot(object = mr.emb, dim.1 = 1, dim.2 = 2) # 这里面只有embryo cells

dim(mr.emb@scale.data)

################################################
# 根据 cell paper 的 300个 marker 对所有embryo cell 重新PCA

cell.marker = read.csv("cell marker.csv",header = F,stringsAsFactors = F)$V1

mr.emb <- RunPCA(object = mr.emb, pc.genes = cell.marker, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5,pcs.compute = 40)

PCElbowPlot(object = mr.emb)

mr.emb <- JackStraw(object = mr.emb, num.replicate = 100, do.print = FALSE,num.pc = 40)

JackStrawPlot(object = mr.emb, PCs = 1:40) # 这里，选择PC10

PCAPlot(object = mr.emb, dim.1 = 1, dim.2 = 2,group.by = "embryo.stage") # 这里面只有embryo cells


##################################################

### 显示具有代表性的lineage marker的表达，尝试标注ICM/TE


FeaturePlot(object = mr.emb, features.plot = c("CDX2","GATA2","GATA3","KRT19","CLDN4","CLDN10","PTGES","PDGFA","DAB2"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca")

FeaturePlot(object = mr.emb, features.plot = c("NANOG","POU5F1","PDGFRA","GDF3","LAMA4","DPPA5"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca")

FeaturePlot(object = mr.emb, features.plot = c("GATA2","GATA3","KRT19","CLDN4"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca")

###################################################

## 尝试使用 cellcyclescoring标注细胞类型

CellLineageScoring <- function (object, TE.genes, PE.genes, EPI.genes, set.ident = FALSE) 
{
  enrich.name <- "Cell lineage"
  genes.list <- list(TE.Score = TE.genes, PE.Score = PE.genes, EPI.Score = EPI.genes)
  object.cc <- AddModuleScore(object = object, genes.list = genes.list, 
                              enrich.name = enrich.name, ctrl.size = min(vapply(X = genes.list, 
                                                                                FUN = length, FUN.VALUE = numeric(1))))
  cc.columns <- grep(pattern = enrich.name, x = colnames(x = object.cc@meta.data))
  cc.scores <- object.cc@meta.data[, cc.columns]
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores, 
                                                                 first = "TE", second = "PE", third = "EPI",null = "Untypical") {
    if (all(scores < 0)) {
      return(null)
    }
    else {
      return(c(first, second, third)[which(x = scores == max(scores))])
    }
  })
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), 
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "TE.Score", "PE.Score", "EPI.Score",
                               "Phase")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("TE.Score", "PE.Score", "EPI.Score","Phase")]
  object <- AddMetaData(object = object, metadata = cc.scores)
  if (set.ident) {
    object <- StashIdent(object = object, save.name = "old.ident")
    object <- SetAllIdent(object = object, id = "Phase")
  }
  return(object)
}

cell.marker.table = read.csv("cell marker.csv",header = F,stringsAsFactors = F)

mr.emb <- CellLineageScoring(mr.emb,TE.genes = cell.marker.table$V1[cell.marker.table$V2 == "TE"], # G2M 对应 TE
                    PE.genes =  cell.marker.table$V1[cell.marker.table$V2 == "PE"],
                    EPI.genes =  cell.marker.table$V1[cell.marker.table$V2 == "EPI"],set.ident = T)

PCAPlot(object = mr.emb, dim.1 = 1, dim.2 = 2,group.by = "Phase")
PCAPlot(object = mr.emb, dim.1 = 1, dim.2 = 2)
PCAPlot(object = mr.emb, dim.1 = 1, dim.2 = 2,group.by = "embryo.stage")

#table(mr.emb@meta.data)
table(mr.emb@meta.data[ mr.emb@meta.data$culture.condition =="Only_embryo" ,c(6,12)])

table(mr.emb@meta.data[ mr.emb@meta.data$culture.condition =="Co-culture" ,c(6,12)])


###################################################


# marker 在各lineage的表达


DotPlot(mr.emb,genes.plot=rev(c("NANOG","POU5F1","GATA6","GATA3","CDX2","CGB","CGB5","CGB8","KRT23")) ,
        group.by = "Phase",
        dot.scale=10,x.lab.rot=T,plot.legend =T)


VlnPlot(object = mr.emb, features.plot = 
          rev(c("NANOG","POU5F1","GATA6","GATA3","CDX2","CGB","CGB5","CGB8","KRT23")), 
        x.lab.rot = TRUE,point.size.use=0)

#################################################

mr.emb <- FindVariableGenes(object = mr.emb, mean.function = ExpMean, dispersion.function = LogVMR)

mr.emb <- RunPCA(object = mr.emb, pc.genes = mr.emb@var.genes, do.print = TRUE, pcs.print = 1:5, 
                 genes.print = 5,pcs.compute = 40)

PCElbowPlot(object = mr.emb)

mr.emb <- JackStraw(object = mr.emb, num.replicate = 100, do.print = FALSE,num.pc = 40)

JackStrawPlot(object = mr.emb, PCs = 1:40) # 这里，选择PC10

PCAPlot(object = mr.emb, dim.1 = 1, dim.2 = 2,group.by = "embryo.stage") # 这里面只有embryo cells

PCAPlot(object = mr.emb, dim.1 = 1, dim.2 = 2,group.by = "Phase") # 这里面只有embryo cells

FeaturePlot(object = mr.emb, features.plot = c("GATA2","GATA3","KRT19","CLDN4"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca")

FeaturePlot(object = mr.emb, features.plot = c("CGB5"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca")

###############
# 找一下unknown 组的marker gene


Untypical.group.markers = FindMarkers(mr.emb,ident.1 = "Untypical",logfc.threshold=0.5849625)

Untypical.group.markers$avg_logFC>0.5849625 & Untypical.group.markers$p_val_adj < 0.05

Untypical.group.markers.gene.name = rownames(Untypical.group.markers)[Untypical.group.markers$avg_logFC>0.5849625 & Untypical.group.markers$p_val_adj < 0.05]

VlnPlot(object = mr.emb, features.plot = c(Untypical.group.markers.gene.name))

##############################################################
# 对所有TE和untypical的细胞，拟合一个lineage 时间点，看发育情况

PCAPlot(object = mr.emb, dim.1 = 1, dim.2 = 2,group.by = "embryo.stage") # 这里面只有embryo cells
PCAPlot(object = mr.emb, dim.1 = 1, dim.2 = 2,group.by = "culture.condition")

mr.emb.TE.UN = SubsetData(object = mr.emb, 
                          cells.use = colnames(mr.emb@scale.data)[    (mr.emb@ident == "TE" | mr.emb@ident == "Untypical")   ])


mr.emb.TE.UN <- FindVariableGenes(object = mr.emb.TE.UN, mean.function = ExpMean, dispersion.function = LogVMR)

mr.emb.TE.UN <- RunPCA(object = mr.emb.TE.UN, pc.genes = mr.emb.TE.UN@var.genes, do.print = TRUE, pcs.print = 1:5, 
                 genes.print = 5,pcs.compute = 40)

mr.emb.TE.UN <- JackStraw(object = mr.emb.TE.UN, num.replicate = 100, do.print = FALSE,num.pc = 40)

JackStrawPlot(object = mr.emb.TE.UN, PCs = 1:40) # 这里，选择PC15

PCAPlot(object = mr.emb.TE.UN, dim.1 = 1, dim.2 = 2,group.by = "embryo.stage") 
PCAPlot(object = mr.emb.TE.UN, dim.1 = 1, dim.2 = 2,group.by = "culture.condition")
PCAPlot(object = mr.emb.TE.UN, dim.1 = 1, dim.2 = 2,group.by = "Phase") 



###################

# 拟合pseudotime

# library(monocle)
# 
# Moc.TE.UN = importCDS( mr.emb.TE.UN )
# 
# Moc.TE.UN <- setOrderingFilter(Moc.TE.UN, mr.emb.TE.UN@var.genes)
# 
# Moc.TE.UN <- reduceDimension(Moc.TE.UN, max_components = 2,
#                             method = 'DDRTree')
# 

library(princurve)

#ddd = mr.emb.TE.UN@scale.data[match(mr.emb.TE.UN@var.genes,rownames(mr.emb.TE.UN@scale.data)),]

ddd = mr.emb.TE.UN@dr$pca@cell.embeddings[,c(1,2)]

pcurve = principal.curve(ddd, smoother="smooth.spline")

plot(pcurve)

points(pcurve)

lines(pcurve)

###################

## 画PCA和pseudotime

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color_hue(3)

### PCA 和 pseudotime 的图
plot(ddd,col=c(gg_color_hue(3))[mr.emb.TE.UN@meta.data$embryo.stage],pch = 16, cex = 1)
lines(pcurve)
points(pcurve)


plot(ddd,col=c(gg_color_hue(4))[as.factor(paste(mr.emb.TE.UN@meta.data$embryo.stage,
                                                mr.emb.TE.UN@meta.data$culture.condition))],pch = 16, cex = 1)
lines(pcurve)
points(pcurve)


#################################################
# 画 pseudotime的boxplot

# 构建画图用的数据框

df.boxplot = cbind(mr.emb.TE.UN@meta.data,lambda = pcurve$lambda)

boxplot(df.boxplot$lambda[ df.boxplot$embryo.stage == "6.5"],
        df.boxplot$lambda[ df.boxplot$embryo.stage == "7.5" & df.boxplot$culture.condition == "Only_embryo"],
        df.boxplot$lambda[ df.boxplot$embryo.stage == "7.5" & df.boxplot$culture.condition == "Co-culture"],
        df.boxplot$lambda[ df.boxplot$embryo.stage == "8.5" & df.boxplot$culture.condition == "Co-culture"],
        col=gg_color_hue(4),outline=F)

boxplot(df.boxplot$lambda[ df.boxplot$embryo.stage == "8.5" & df.boxplot$culture.condition == "Co-culture"],
        df.boxplot$lambda[ df.boxplot$embryo.stage == "7.5" & df.boxplot$culture.condition == "Co-culture"],
        df.boxplot$lambda[ df.boxplot$embryo.stage == "7.5" & df.boxplot$culture.condition == "Only_embryo"],
        df.boxplot$lambda[ df.boxplot$embryo.stage == "6.5"],
        col=gg_color_hue(4),horizontal=F)

wilcox.test(df.boxplot$lambda[ df.boxplot$embryo.stage == "7.5" & df.boxplot$culture.condition == "Only_embryo"],
df.boxplot$lambda[ df.boxplot$embryo.stage == "7.5" & df.boxplot$culture.condition == "Co-culture"])

####################################################

# 这个是day 8.5的TE subtype，typical和untypical
 
mr.emb.TE.UN.8.5 = SubsetData(object = mr.emb.TE.UN, 
                              cells.use = colnames(mr.emb.TE.UN@scale.data)[ mr.emb.TE.UN@meta.data$embryo.stage ==   "8.5"  ])

mr.emb.TE.UN.8.5 <- FindVariableGenes(object = mr.emb.TE.UN.8.5, mean.function = ExpMean, dispersion.function = LogVMR)

mr.emb.TE.UN.8.5 <- RunPCA(object = mr.emb.TE.UN.8.5, pc.genes = mr.emb.TE.UN.8.5@var.genes, do.print = TRUE, pcs.print = 1:5, 
                       genes.print = 5,pcs.compute = 40)

mr.emb.TE.UN.8.5 <- JackStraw(object = mr.emb.TE.UN.8.5, num.replicate = 100, do.print = FALSE,num.pc = 40)

JackStrawPlot(object = mr.emb.TE.UN.8.5, PCs = 1:40) # 这里，选择PC8

PCAPlot(object = mr.emb.TE.UN.8.5, dim.1 = 1, dim.2 = 2) 

mr.emb.TE.UN.8.5 <- RunTSNE(mr.emb.TE.UN.8.5,dims.use = 1:8)

TSNEPlot(mr.emb.TE.UN.8.5,group.by = "embryo.no")

mr.emb.TE.UN.8.5 <- FindClusters(object = mr.emb.TE.UN.8.5, reduction.type = "pca", dims.use = 1:8, 
             resolution = 0.6, print.output = 0, save.SNN = TRUE)

TSNEPlot(mr.emb.TE.UN.8.5)

mr.emb.TE.UN.8.5.markers <- FindAllMarkers(object = mr.emb.TE.UN.8.5, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)

mr.emb.TE.UN.8.5.markers.filtered = mr.emb.TE.UN.8.5.markers[mr.emb.TE.UN.8.5.markers$avg_logFC > log2(1.5),]

write.csv(mr.emb.TE.UN.8.5.markers, file = "mr.emb.TE.UN.8.5.markers.csv",quote=F)

FeaturePlot(object = mr.emb.TE.UN.8.5, features.plot = c("CGB","CGB5","KRT23","PGF","GCM1"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

VlnPlot(object = mr.emb.TE.UN.8.5, features.plot = c("E2F8","KRT23","PGF","GCM1","PRDX3","LRP2"), use.raw = TRUE, y.log = TRUE,
        nCol = 2)

VlnPlot(object = mr.emb.TE.UN.8.5, features.plot = c("GATA3","NANOG"), use.raw = TRUE, y.log = TRUE,
        nCol = 2)


VlnPlot(object = mr.emb.TE.UN.8.5, features.plot = c("CGB","CGB5","CGB8"), use.raw = TRUE, y.log = TRUE,
        nCol = 2)

match(c("E2F8","KRT23","PGF","GCM1"),mr.emb.TE.UN.8.5.markers$gene)



complete.cluster.8.5.nocombat = hclust(as.dist(1-abs(cor(mr.emb.TE.UN.8.5@scale.data[match(mr.emb.TE.UN.8.5@var.genes,rownames(mr.emb.TE.UN.8.5@scale.data)),],method="spearman"))), 
                              method="ward.D")

plot(complete.cluster.8.5.nocombat, hang = -1)


############################################

# 用combat对day 8.5进行更精细的分类

library(sva)

modcombat = model.matrix(~1, data=mr.emb.TE.UN.8.5@meta.data)

combat_edata = ComBat(dat=(mr.emb.TE.UN.8.5@scale.data), batch=  as.numeric(mr.emb.TE.UN.8.5@meta.data$embryo.no)  ,
                      mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

mr.emb.TE.UN.8.5.combat = mr.emb.TE.UN.8.5
mr.emb.TE.UN.8.5.combat@scale.data = combat_edata


mr.emb.TE.UN.8.5.combat <- FindVariableGenes(object = mr.emb.TE.UN.8.5.combat, mean.function = ExpMean, dispersion.function = LogVMR)

mr.emb.TE.UN.8.5.combat <- RunPCA(object = mr.emb.TE.UN.8.5.combat, do.print = TRUE, pcs.print = 1:5,  
                           genes.print = 5,pcs.compute = 40) 
mr.emb.TE.UN.8.5.combat <- RunTSNE(mr.emb.TE.UN.8.5.combat)

PCAPlot(mr.emb.TE.UN.8.5.combat,group.by = "embryo.no")

TSNEPlot(object = mr.emb.TE.UN.8.5.combat,group.by = "embryo.no")


##################################################

# hclust 聚类画图

complete.cluster.8.5 = hclust(as.dist(1-abs(cor(mr.emb.TE.UN.8.5.combat@scale.data[match(mr.emb.TE.UN.8.5.combat@var.genes,rownames(mr.emb.TE.UN.8.5.combat@scale.data)),],method="spearman"))), 
                              method="ward.D")

plot(complete.cluster.8.5, hang = -1)

mr.TE.8.5.combat.leaf = cutree(complete.cluster.8.5,k = 5)

rect.hclust(complete.cluster.8.5,k = 5)

table(mr.TE.8.5.combat.leaf)


day8.5.exp.boxplot <- function(the.gene){
  require(ggplot2)
  dat = data.frame(exp = mr.emb.TE.UN.8.5.combat@scale.data[match(the.gene,rownames(mr.emb.TE.UN.8.5.combat@scale.data)),],
                   mr.TE.8.5.combat.leaf)
  ggplot(dat, aes(x = factor(mr.TE.8.5.combat.leaf), y = exp, fill = factor(mr.TE.8.5.combat.leaf))) +  
    geom_boxplot() + ggtitle(the.gene)
}

##########################
# 对5类进行画图，看某一个基因在5类中的表达

day8.5.exp.boxplot("CDX2")
day8.5.exp.boxplot("GATA3")

day8.5.exp.boxplot("GATA6")
day8.5.exp.boxplot("POU5F1")

day8.5.exp.boxplot("PARD6A")
day8.5.exp.boxplot("KRT7")
day8.5.exp.boxplot("CGB5")
day8.5.exp.boxplot("KRT23")
day8.5.exp.boxplot("CGB8")
day8.5.exp.boxplot("CGB")

day8.5.exp.boxplot("TUBA1C")
day8.5.exp.boxplot("NR6A1")
day8.5.exp.boxplot("PRDX3")
day8.5.exp.boxplot("STMN1")
day8.5.exp.boxplot("MYH10")
day8.5.exp.boxplot("NR6A1")
day8.5.exp.boxplot("SHISA5")
day8.5.exp.boxplot("FABP5")
day8.5.exp.boxplot("MYH9")

day8.5.exp.boxplot("CCR7")
day8.5.exp.boxplot("GCM1")



##########################################

mr.emb.TE.UN.8.5.combat.groups = mr.emb.TE.UN.8.5.combat

mr.emb.TE.UN.8.5.combat.groups@ident = factor(mr.TE.8.5.combat.leaf)

VlnPlot(mr.emb.TE.UN.8.5.combat.groups,features.plot = c("CGB5"))
VlnPlot(mr.emb.TE.UN.8.5.combat.groups,features.plot = c("PRDX3"))

# 包括高表达低表达两种marker
mr.TE.8.5.markers <- FindAllMarkers(mr.emb.TE.UN.8.5.combat.groups, min.pct = 0.25, thresh.use = 0.05)
mr.TE.8.5.markers.sig = mr.TE.8.5.markers[mr.TE.8.5.markers$p_val_adj<0.05,]
write.csv(mr.TE.8.5.markers.sig,file = "mr.TE.8.5.markers.sig.positive and negative.csv",quote = F)


# 仅包括高表达marker，用来画图
mr.TE.8.5.markers <- FindAllMarkers(mr.emb.TE.UN.8.5.combat.groups, only.pos = T ,min.pct = 0.25, thresh.use = 0.05)
mr.TE.8.5.markers.sig = mr.TE.8.5.markers[mr.TE.8.5.markers$p_val_adj<0.05 & abs(mr.TE.8.5.markers$avg_logFC) > log2(1.5),]
png(file = "b4.jpg",width = 2000,height = 1800)
DoHeatmap(mr.emb.TE.UN.8.5.combat.groups, genes.use = mr.TE.8.5.markers.sig$gene, 
          slim.col.label = TRUE, remove.key = TRUE)
dev.off()
write.csv(mr.TE.8.5.markers.sig,file = "mr.TE.8.5.markers.sig.positive only.csv",quote = F)



mr.emb.TE.UN.8.5 <- SetIdent(mr.emb.TE.UN.8.5, ident.use = mr.TE.8.5.combat.leaf)

VlnPlot(mr.emb.TE.UN.8.5,features.plot = c("CGB"))
VlnPlot(mr.emb.TE.UN.8.5,features.plot = c("CGB8"))
VlnPlot(mr.emb.TE.UN.8.5,features.plot = c("KRT23"))
VlnPlot(mr.emb.TE.UN.8.5,features.plot = c("CCR7"))

VlnPlot(mr.emb.TE.UN.8.5.combat.groups,features.plot = c("GATA6"))
VlnPlot(mr.emb.TE.UN.8.5.combat.groups,features.plot = c("GATA3"))
VlnPlot(mr.emb.TE.UN.8.5.combat.groups,features.plot = c("KRT7"))
VlnPlot(mr.emb.TE.UN.8.5.combat.groups,features.plot = c("CGB"))
VlnPlot(mr.emb.TE.UN.8.5.combat.groups,features.plot = c("CDX2"))

VlnPlot(mr.emb.TE.UN.8.5.combat.groups,features.plot = c("GATA3"))
################################################################
# 看不同时期TE之间的比较

# mr.emb.TE.UN

# 修改一下ident，方便比较gene

mr.emb.TE.UN <- SetIdent(mr.emb.TE.UN, ident.use = paste(mr.emb.TE.UN@meta.data$embryo.stage,mr.emb.TE.UN@meta.data$culture.condition))

TSNEPlot(mr.emb.TE.UN)

# 1 
Seven.emb.to.Seven.co <- FindMarkers(mr.emb.TE.UN, ident.1 = "7.5 Only_embryo", ident.2 = "7.5 Co-culture", 
                                 test.use = "wilcox")
write.csv(Seven.emb.to.Seven.co,file = "Day7.5 embryo only to day 7.5 co-culture.csv",quote=F)


# 2
Six.emb.to.Seven.emb <- FindMarkers(mr.emb.TE.UN, ident.1 = "6.5 Only_embryo", ident.2 = "7.5 Only_embryo", 
                                     test.use = "wilcox")
write.csv(Six.emb.to.Seven.emb,file = "Day6.5 embryo only to day 7.5 embryo only.csv",quote=F)

#3
Six.emb.to.Seven.co <- FindMarkers(mr.emb.TE.UN, ident.1 = "6.5 Only_embryo", ident.2 = "7.5 Co-culture", 
                                    test.use = "wilcox")
write.csv(Six.emb.to.Seven.co,file = "Day6.5 embryo only to day 7.5 co-culture.csv",quote=F)

#4
Seven.emb.to.Eight.co <- FindMarkers(mr.emb.TE.UN, ident.1 = "7.5 Only_embryo", ident.2 = "8.5 Co-culture", 
                                   test.use = "wilcox")
write.csv(Seven.emb.to.Eight.co,file = "Day7.5 embryo only to day 8.5 co-culture.csv",quote=F)

#5
Seven.co.to.Eight.co <- FindMarkers(mr.emb.TE.UN, ident.1 = "7.5 Co-culture", ident.2 = "8.5 Co-culture", 
                                     test.use = "wilcox")
write.csv(Seven.co.to.Eight.co,file = "Day7.5 Co-culture to day 8.5 co-culture.csv",quote=F)


#############################

# polar/mural TE 区分。day 8.5的分出来了好像，现在来份day7.5，cultured

mr.emb.TE.7.5 = SubsetData(object = mr.emb.TE.UN, 
                              cells.use = colnames(mr.emb.TE.UN@scale.data)[ mr.emb.TE.UN@meta.data$embryo.stage ==   "7.5"  ])

mr.emb.TE.7.5.CoCul = SubsetData(object = mr.emb.TE.7.5, 
                                 cells.use = colnames(mr.emb.TE.7.5@scale.data)[ mr.emb.TE.7.5@meta.data$culture.condition ==  
                                                                                   "Co-culture"  ])
mr.emb.TE.7.5.Emb = SubsetData(object = mr.emb.TE.7.5, 
                                 cells.use = colnames(mr.emb.TE.7.5@scale.data)[ mr.emb.TE.7.5@meta.data$culture.condition ==  
                                                                                   "Only_embryo"  ])

VlnPlot(mr.emb.TE.7.5,features.plot = c("CCR7"))
VlnPlot(mr.emb.TE.7.5,features.plot = c("CGB5"))
VlnPlot(mr.emb.TE.7.5,features.plot = c("CGB7"))
VlnPlot(mr.emb.TE.7.5,features.plot = c("CGB8"))

########################################

# 现在来份day7.5，ccultured

mr.emb.TE.7.5.CoCul
mr.emb.TE.7.5.CoCul <- FindVariableGenes(object = mr.emb.TE.7.5.CoCul, mean.function = ExpMean, dispersion.function = LogVMR)

# 使用combat 修正embryo effect

library(sva)

modcombat = model.matrix(~1, data=mr.emb.TE.7.5.CoCul@meta.data)

combat_edata = ComBat(dat=(mr.emb.TE.7.5.CoCul@scale.data), batch=  as.numeric(mr.emb.TE.7.5.CoCul@meta.data$embryo.no)  ,
                      mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

mr.emb.TE.7.5.CoCul.combat = mr.emb.TE.7.5.CoCul
mr.emb.TE.7.5.CoCul.combat@scale.data = combat_edata


mr.emb.TE.7.5.CoCul.combat <- FindVariableGenes(object = mr.emb.TE.7.5.CoCul.combat, mean.function = ExpMean, dispersion.function = LogVMR)

mr.emb.TE.7.5.CoCul.combat <- RunPCA(object = mr.emb.TE.7.5.CoCul.combat, do.print = TRUE, pcs.print = 1:5,  
                                  genes.print = 5,pcs.compute = 40) 

PCAPlot(mr.emb.TE.7.5.CoCul.combat,group.by = "embryo.no")

complete.cluster.7.5.cocul.combat = hclust(as.dist(1-abs(cor(mr.emb.TE.7.5.CoCul.combat@scale.data[match(mr.emb.TE.7.5.CoCul.combat@var.genes,rownames(mr.emb.TE.7.5.CoCul.combat@scale.data)),],method="spearman"))), 
                              method="ward.D")

plot(complete.cluster.7.5.cocul.combat, hang = -1)

complete.cluster.7.5.cocul.combat.leaf = cutree(complete.cluster.7.5.cocul.combat,k = 4)
rect.hclust(complete.cluster.7.5.cocul.combat,k = 4)

mr.emb.TE.7.5.CoCul.combat <- SetIdent(mr.emb.TE.7.5.CoCul.combat,ident.use = complete.cluster.7.5.cocul.combat.leaf)

VlnPlot(mr.emb.TE.7.5.CoCul.combat,features.plot = c("CGB8"))
VlnPlot(mr.emb.TE.7.5.CoCul.combat,features.plot = c("CCR7"))

PCAPlot(mr.emb.TE.7.5.CoCul)
PCAPlot(mr.emb.TE.7.5.CoCul,group.by = "embryo.no")

###################################
# C7.5 combat修正过的 四组，每一组的marker


mr.emb.TE.7.5.CoCul.combat.markers <- FindAllMarkers(mr.emb.TE.7.5.CoCul.combat, only.pos = T ,min.pct = 0.25, thresh.use = 0.05)
mr.emb.TE.7.5.CoCul.combat.markers.sig = mr.emb.TE.7.5.CoCul.combat.markers[mr.emb.TE.7.5.CoCul.combat.markers$p_val_adj<0.05 & abs(mr.emb.TE.7.5.CoCul.combat.markers$avg_logFC) > log2(1.5),]
png(file = "b6.jpg",width = 2000,height = 1800)
DoHeatmap(mr.emb.TE.7.5.CoCul.combat, genes.use = mr.emb.TE.7.5.CoCul.combat.markers.sig$gene, 
          slim.col.label = TRUE, remove.key = TRUE)
dev.off()
write.csv(mr.emb.TE.7.5.CoCul.markers.sig,file = "mr.TE.7.5.combat.markers.sig.positive only.csv",quote = F)
# 这个效果很差

###################################
# C7.5 combat修正过的 四组，第三组 polar TE 和其他组之间差异表达的基因

mr.emb.TE.7.5.CoCul.combat.3.markers = FindMarkers(mr.emb.TE.7.5.CoCul.combat,ident.1 = "3")

write.csv(mr.emb.TE.7.5.CoCul.combat.3.markers,file = "mr.emb.TE.7.5.CoCul.combat.3.markers.csv",quote=F)

VlnPlot(mr.emb.TE.7.5.CoCul.combat,features.plot = c("TNFAIP3"))
VlnPlot(mr.emb.TE.7.5.CoCul.combat,features.plot = c("GREM2"))

VlnPlot(mr.emb.TE.7.5.CoCul.combat,features.plot = c("CGB"))

######################################
# 不同时期的mural/polar TE之间的比较：C7.5和C8.5

mr.emb.TE.7.5.CoCul.8.5.combat <- MergeSeurat(object1 = mr.emb.TE.7.5.CoCul.combat, object2 = mr.emb.TE.UN.8.5.combat.groups, 
                             add.cell.id1 = "7.5", 
                             add.cell.id2 = "8.5", project = "Lvbo")

mr.emb.TE.7.5.CoCul.8.5.combat <- FindVariableGenes(object = mr.emb.TE.7.5.CoCul.8.5.combat, 
                                                    
                                                    mean.function = ExpMean, dispersion.function = LogVMR)
mr.emb.TE.7.5.CoCul.8.5.combat <- ScaleData(object = mr.emb.TE.7.5.CoCul.8.5.combat, vars.to.regress = c("nUMI", "percent.mito"))

mr.emb.TE.7.5.CoCul.8.5.combat <- RunPCA(object = mr.emb.TE.7.5.CoCul.8.5.combat, do.print = TRUE, pcs.print = 1:5,  
                                     genes.print = 5,pcs.compute = 40) 

PCAPlot(mr.emb.TE.7.5.CoCul.8.5.combat)

FeaturePlot(object = mr.emb.TE.7.5.CoCul.8.5.combat, features.plot = c("CGB","CGB5","CCR7"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca")

mr.emb.TE.7.5.CoCul.8.5.combat <- SetIdent(mr.emb.TE.7.5.CoCul.8.5.combat, ident.use = c( 
  paste("7.5",mr.emb.TE.7.5.CoCul.combat@ident,sep="_"),
  paste("8.5",mr.emb.TE.UN.8.5.combat.groups@ident,sep="_")
  ))

PCAPlot(mr.emb.TE.7.5.CoCul.8.5.combat) # lineage segregation analysis

mr.emb.TE.7.5.CoCul.8.5.combat <- AddMetaData(object = mr.emb.TE.7.5.CoCul.8.5.combat, 
                                              metadata = mr.emb.TE.7.5.CoCul.8.5.combat@ident, 
                                              col.name = c("cell.stages") )

################################################
# ST cell dominate 差异gene。 去掉8.5_5, 看得到的差异gene 和monocle 是不是好一点

names(mr.emb.TE.7.5.CoCul.8.5.combat@ident)[( mr.emb.TE.7.5.CoCul.8.5.combat@ident != "8.5_5")]

mr.emb.TE.7.5.CoCul.8.5.noST = SubsetData(mr.emb.TE.7.5.CoCul.8.5.combat, 
                                          cells.use = names(mr.emb.TE.7.5.CoCul.8.5.combat@ident)[( mr.emb.TE.7.5.CoCul.8.5.combat@ident != "8.5_5")])

mr.emb.TE.7.5.CoCul.8.5.noST <- FindVariableGenes(object = mr.emb.TE.7.5.CoCul.8.5.noST, 
                                                    mean.function = ExpMean, dispersion.function = LogVMR,
                                                  x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 2)

length(mr.emb.TE.7.5.CoCul.8.5.noST@var.genes)

mr.emb.TE.7.5.CoCul.8.5.noST <- ScaleData(object = mr.emb.TE.7.5.CoCul.8.5.noST, vars.to.regress = c("nUMI", "percent.mito"))


mr.emb.TE.7.5.CoCul.8.5.noST <- RunPCA(object = mr.emb.TE.7.5.CoCul.8.5.noST, do.print = TRUE, pcs.print = 1:5,  
                                         genes.print = 5,pcs.compute = 40) 


PCAPlot(mr.emb.TE.7.5.CoCul.8.5.noST)


mr.emb.TE.7.5.CoCul.8.5.combat.newPCA = RunPCA(object = mr.emb.TE.7.5.CoCul.8.5.combat, 
                                               pc.genes =  mr.emb.TE.7.5.CoCul.8.5.noST@var.genes,
                                               do.print = TRUE, pcs.print = 1:5,  
                                               genes.print = 5,pcs.compute = 40) 

PCAPlot(mr.emb.TE.7.5.CoCul.8.5.combat.newPCA)

################################################
# use monocole to analyze cell lineage

cutree(d7.5.hclust,k=2)

inside.1 = names(cutree(d7.5.hclust,k=2)) %in% substr(names(mr.emb.TE.7.5.CoCul.8.5.combat@ident  ),5,20)

polar.tag = cutree(d7.5.hclust,k=2)[inside.1]

old.tag = as.character(mr.emb.TE.7.5.CoCul.8.5.combat@ident)

old.tag[match(names(polar.tag), substr(names(mr.emb.TE.7.5.CoCul.8.5.combat@ident  ),5,20) )] = polar.tag

mr.emb.TE.7.5.CoCul.8.5.combat.polarLabeled = SetIdent(mr.emb.TE.7.5.CoCul.8.5.combat,ident.use =  old.tag)

PCAPlot(mr.emb.TE.7.5.CoCul.8.5.combat.polarLabeled)

mr.emb.TE.7.5.CoCul.8.5.combat.polarLabeled@meta.data$polar.ident =  mr.emb.TE.7.5.CoCul.8.5.combat.polarLabeled@ident


###############################################
library(monocle)

Mono.emb.TE.7.5.CoCul.8.5.combat <- importCDS(mr.emb.TE.7.5.CoCul.8.5.combat.polarLabeled, import_all = TRUE)

Mono.emb.TE.7.5.CoCul.8.5.combat <- estimateSizeFactors(Mono.emb.TE.7.5.CoCul.8.5.combat)

Mono.emb.TE.7.5.CoCul.8.5.combat <- 
  setOrderingFilter(Mono.emb.TE.7.5.CoCul.8.5.combat, 
                    mr.emb.TE.7.5.CoCul.8.5.combat.polarLabeled@var.genes)


Mono.emb.TE.7.5.CoCul.8.5.combat <- 
  reduceDimension(Mono.emb.TE.7.5.CoCul.8.5.combat, max_components = 2,num_dim = 2,
                            method = 'DDRTree')

Mono.emb.TE.7.5.CoCul.8.5.combat <- orderCells(Mono.emb.TE.7.5.CoCul.8.5.combat)

plot_cell_trajectory(Mono.emb.TE.7.5.CoCul.8.5.combat,color_by="polar.ident")

plot_cell_trajectory(Mono.emb.TE.7.5.CoCul.8.5.combat)

Mono.emb.TE.7.5.CoCul.8.5.combat <- orderCells(Mono.emb.TE.7.5.CoCul.8.5.combat,root_state = 1)

plot_cell_trajectory(Mono.emb.TE.7.5.CoCul.8.5.combat,color_by="cell.stages")

plot_cell_trajectory(Mono.emb.TE.7.5.CoCul.8.5.combat,markers  = c("CGB5","GCM1","CCR7") )

plot_cell_trajectory(Mono.emb.TE.7.5.CoCul.8.5.combat,color_by = "Pseudotime")


blast_genes <- row.names(subset(fData(Mono.emb.TE.7.5.CoCul.8.5.combat),
                                gene_short_name %in% c("CGB5","GCM1","CCR7","CTNNB1","ITGA6","LRP5","TP63","FGFR2","FZD5",
                                                       "FN1","MMP2","ITGA5","CD9")))

plot_genes_jitter(Mono.emb.TE.7.5.CoCul.8.5.combat[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1,
                  color_by = "State",
                  relative_expr = F)


plot_genes_branched_pseudotime(Mono.emb.TE.7.5.CoCul.8.5.combat[c("CGB5","GCM1","CCR7","CTNNB1","ITGA6","LRP5","TP63","FGFR2","FZD5",
                                                                  "FN1","MMP2","ITGA5","CD9"),],
                               branch_point = 2,
                               color_by = "cell.stages",
                               ncol = 1)

plot_genes_branched_pseudotime(Mono.emb.TE.7.5.CoCul.8.5.combat[c("FN1","MMP2","ITGA5","CD9","ITGA1"),],
                               branch_point = 2,
                               color_by = "cell.stages",
                               ncol = 1)


plot_genes_branched_pseudotime(Mono.emb.TE.7.5.CoCul.8.5.combat[c("CGB5","WNT7A"),],
                               branch_point = 2,
                               color_by = "cell.stages",
                               ncol = 1)


to_be_tested <- row.names(subset(fData(Mono.emb.TE.7.5.CoCul.8.5.combat),
                                 gene_short_name %in% mr.emb.TE.7.5.CoCul.8.5.combat@var.genes ))

cds_subset <- Mono.emb.TE.7.5.CoCul.8.5.combat[to_be_tested,]

#plot_genes_in_pseudotime(cds_subset, color_by = "cell.stages")


BEAM_res <- BEAM(Mono.emb.TE.7.5.CoCul.8.5.combat, branch_point = 2, cores = 4)

rownames(BEAM_res[(BEAM_res$use_for_ordering == T & BEAM_res$qval<1e-150),])

#BEAM_res <- BEAM_res[order(BEAM_res$qval),]

#BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

head(BEAM_res,50)

pdf(width = 6,height = 50,file = "a.pdf")
plot_genes_branched_heatmap(Mono.emb.TE.7.5.CoCul.8.5.combat[rownames(BEAM_res[(BEAM_res$use_for_ordering == T & BEAM_res$qval<1e-150),]),],
                            branch_point = 2,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()
######################################

# day 6.5 分类
mr.emb.TE.6.5 = SubsetData(object = mr.emb.TE.UN, 
                           cells.use = colnames(mr.emb.TE.UN@scale.data)[ mr.emb.TE.UN@meta.data$embryo.stage ==   "6.5"  ])


library(sva)

modcombat = model.matrix(~1, data=mr.emb.TE.6.5@meta.data)

combat_edata = ComBat(dat=(mr.emb.TE.6.5@scale.data), batch=  as.numeric(mr.emb.TE.6.5@meta.data$embryo.no)  ,
                      mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

mr.emb.TE.6.5.combat = mr.emb.TE.6.5
mr.emb.TE.6.5.combat@scale.data = combat_edata


mr.emb.TE.6.5.combat <- FindVariableGenes(object = mr.emb.TE.6.5.combat, mean.function = ExpMean, dispersion.function = LogVMR)

mr.emb.TE.6.5.combat <- RunPCA(object = mr.emb.TE.6.5.combat, do.print = TRUE, pcs.print = 1:5,  
                                     genes.print = 5,pcs.compute = 10) 

PCAPlot(mr.emb.TE.6.5.combat,group.by = "embryo.no")

complete.cluster.7.5.cocul.combat = hclust(as.dist(1-abs(cor(mr.emb.TE.7.5.CoCul.combat@scale.data[match(mr.emb.TE.7.5.CoCul.combat@var.genes,rownames(mr.emb.TE.7.5.CoCul.combat@scale.data)),],method="spearman"))), 
                                           method="ward.D")

plot(complete.cluster.7.5.cocul.combat, hang = -1)

complete.cluster.7.5.cocul.combat.leaf = cutree(complete.cluster.7.5.cocul.combat,k = 4)
rect.hclust(complete.cluster.7.5.cocul.combat,k = 4)

mr.emb.TE.7.5.CoCul.combat <- SetIdent(mr.emb.TE.7.5.CoCul.combat,ident.use = complete.cluster.7.5.cocul.combat.leaf)


VlnPlot(mr.emb.TE.7.5.CoCul.combat,features.plot = c("CCR7"))

########
# 看 mr.emb.TE.7.5.CoCul.combat ， 四类细胞都是啥

mr.emb.TE.7.5.CoCul.4group = SetIdent(mr.emb.TE.7.5.CoCul,ident.use = mr.emb.TE.7.5.CoCul.combat@ident)
  
marker.mr.emb.TE.7.5.CoCul.combat = FindAllMarkers(mr.emb.TE.7.5.CoCul.4group)
marker.mr.emb.TE.7.5.CoCul.combat.sig = marker.mr.emb.TE.7.5.CoCul.combat[marker.mr.emb.TE.7.5.CoCul.combat$p_val_adj<0.1 & abs(marker.mr.emb.TE.7.5.CoCul.combat$avg_logFC)>log2(1.5),]

DoHeatmap(mr.emb.TE.7.5.CoCul.combat,genes.use = marker.mr.emb.TE.7.5.CoCul.combat.sig$gene)



############################################

# tangfuchou 的 polar/Mural

library(readxl)

tang.data = read_excel("nsmb.2660-S2-mural-polar-TE.xlsx",sheet = "t1")
tang.data2 = read_excel("nsmb.2660-S2-mural-polar-TE.xlsx",sheet = "t1")
rownames(tang.data) = tang.data$Gene_ID
tang.data = tang.data[,-1]
dim(tang.data)
tang.data = as.matrix(tang.data)

rownames(tang.data) = tang.data2$Gene_ID
colnames(tang.data)
head(tang.data)

tang.data.log = log10(tang.data+1)

p.value = apply(tang.data.log,MARGIN = 1,FUN = function(x){t.test(x[which(str_sub(colnames(tang.data),1,1) == "P")],
                                                                  x[which(str_sub(colnames(tang.data),1,1) == "M")]
                                                                  )$p.value})


fold = rowMeans(tang.data[,which(str_sub(colnames(tang.data),1,1) == "P")])/rowMeans(tang.data[,which(str_sub(colnames(tang.data),1,1) == "M")])

polar.markers = intersect(names(which((p.value) < 0.01)),names(which(fold > 2) ))
polar.markers

mural.markers = intersect(names(which((p.value) < 0.05)),names(which(fold < 0.5) ))
mural.markers

mr.emb.TE.UN



CellLineageScoring.polar.mural <- function (object, Polar.genes,  set.ident = FALSE) 
{
  enrich.name <- "Cell lineage"
  genes.list <- list(Polar.Score = Polar.genes)
  object.cc <- AddModuleScore(object = object, genes.list = genes.list, 
                              enrich.name = enrich.name, ctrl.size = min(vapply(X = genes.list, 
                                                                                FUN = length, FUN.VALUE = numeric(1))))
  cc.columns <- grep(pattern = enrich.name, x = colnames(x = object.cc@meta.data))
  cc.scores <- object.cc@meta.data[, cc.columns]
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores, 
                                                                 first = "Polar", null = "Untypical") {
    if (all(scores < 0)) {
      return(null)
    }
    else {
      return(c(first, second, third)[which(x = scores == max(scores))])
    }
  })
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), 
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "Polar.Score", 
                               "Phase")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("Polar.Score",  "Phase")]
  object <- AddMetaData(object = object, metadata = cc.scores)
  if (set.ident) {
    object <- StashIdent(object = object, save.name = "old.ident")
    object <- SetAllIdent(object = object, id = "Phase")
  }
  return(object)
}

mr.emb.TE.UN.polar.mural = CellLineageScoring.polar.mural(mr.emb.TE.UN,Polar.genes  = polar$polar,set.ident = T)

PCAPlot(mr.emb.TE.UN.polar.mural)



polar = read_excel("Polar.xlsx",sheet = "Sheet1")
polar$polar



te6.5 = mr.emb.TE.6.5@scale.data[na.omit(match(c(polar$polar),rownames(mr.emb.TE.6.5@scale.data))),]
library(gplots)
palette.breaks <- seq(-1.5, 1.5, 0.1)
color.palette = colorRampPalette(c("dodgerblue4","dodgerblue1","white","firebrick1","firebrick3"), 
                                 space="Lab")(length(palette.breaks) - 1)

heatmap.2(as.matrix( te6.5 ),
          trace="none",density="none",
          Colv = as.dendrogram(hclust(as.dist(1-abs(cor( te6.5,method="spearman"))), 
                                      method="ward.D")),
          Rowv = F,
          key=T,scale="row",
          dendrogram="both",
          margins=c(10.5,6),
          labRow=NA,
          col=color.palette,
          main="Day 6.5",
          
)



te7.5 = mr.emb.TE.7.5@scale.data[na.omit(match(c(polar$polar),rownames(mr.emb.TE.7.5@scale.data))),]
library(gplots)
palette.breaks <- seq(-1.5, 1.5, 0.1)
color.palette = colorRampPalette(c("dodgerblue4","dodgerblue1","white","firebrick1","firebrick3"), 
                                 space="Lab")(length(palette.breaks) - 1)

d7.5.hclust = hclust(as.dist(1-abs(cor( te7.5,method="spearman"))), 
                     method="ward.D") 

d7.5.hclust.tree.group = cutree(d7.5.hclust,k = 2)

d7.5.hclust.condition = substr(names(d7.5.hclust.tree.group),1,1)

table(data.frame(d7.5.hclust.tree.group,d7.5.hclust.condition))


heatmap.2(as.matrix( te7.5 ),
          trace="none",density="none",
          Colv = as.dendrogram(hclust(as.dist(1-abs(cor( te7.5,method="spearman"))), 
                                      method="ward.D")),
          Rowv = F,
          key=T,scale="row",
          dendrogram="both",
          margins=c(10.5,6),
          labRow=NA,
          breaks = seq(-2.5, 2.5, 0.1),
          col=colorRampPalette(c("#FF00FF","#000000","#FFFF00"), space="Lab"),
          main="Day 7.5",
          ColSideColors = c("orange","royalblue")[as.factor(d7.5.hclust.condition)]
          
)

te.clust.tree = (hclust(as.dist(1-abs(cor( te7.5,method="spearman"))), 
                        method="ward.D"))

plot(te.clust.tree,hang = -1)

library(RColorBrewer)
n <- 7
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
dendextend::colored_bars( col_vector[as.numeric(factor(substr(te.clust.tree$labels,1,5)))][te.clust.tree$order]   )


###########################################

VlnPlot(object = mr.emb.TE.7.5, features.plot = 
          rev(c("CCR7","CYP19A1","DLX5","ERVFRD-1","GCM1","GREM2","MUC15","OVOL1")), 
        x.lab.rot = TRUE,point.size.use=0,cols.use=c("royalblue","orange"))


##########################################

te8.5 = mr.emb.TE.UN.8.5@scale.data[na.omit(match(c(polar$polar),rownames(mr.emb.TE.UN.8.5@scale.data))),]
library(gplots)
palette.breaks <- seq(-1.5, 1.5, 0.1)
color.palette = colorRampPalette(c("dodgerblue4","dodgerblue1","white","firebrick1","firebrick3"), 
                                 space="Lab")(length(palette.breaks) - 1)

heatmap.2(as.matrix( te8.5 ),
          trace="none",density="none",
          Colv = as.dendrogram(hclust(as.dist(1-abs(cor( te8.5,method="spearman"))), 
                                      method="ward.D")),
          Rowv = F,
          key=T,scale="row",
          dendrogram="both",
          margins=c(10.5,6),
          labRow=NA,
          col=color.palette,
          main="Day 8.5",
          
)


############################

# day 6.5 and day 7.5 WGCNA



mr.emb.TE.6.5
mr.emb.TE.7.5

mr.emb.TE.6.5and7.5 <- MergeSeurat(object1 = mr.emb.TE.6.5, object2 = mr.emb.TE.7.5, add.cell.id1 = "6.5", 
                             add.cell.id2 = "7.5", project = "WGCNA")

RPM.emb.TE.6.5and7.5 = 
  sweep(mr.emb.TE.6.5and7.5@raw.data,2,colSums(mr.emb.TE.6.5and7.5@raw.data),"/" ) * 1000000

sum(rowMeans(RPM.emb.TE.6.5and7.5) > 2)

RPM.emb.TE.6.5and7.5.fil = RPM.emb.TE.6.5and7.5[rowMeans(RPM.emb.TE.6.5and7.5) > 2,]

#RPM.emb.TE.6.5and7.5.fil = RPM.emb.TE.6.5and7.5[apply(RPM.emb.TE.6.5and7.5,1,function(x) sum(x>1))>=40,]

library(WGCNA)
library("flashClust")
options(stringsAsFactors = FALSE);
#enableWGCNAThreads()
disableWGCNAThreads()


mydata=log2(RPM.emb.TE.6.5and7.5.fil+1)
dim(mydata)
gene.names = rownames(mydata)

mydata.trans=t(mydata);

datExpr=mydata.trans
SubGeneNames=gene.names

######
powers = c(1:30)
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,
                      powerVector = powers,corFnc = cor,
                      corOptions = list(use = 'p'),
                      networkType = "signed")


par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main = paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##################

softPower = 6;

adj= adjacency(datExpr,type = "signed", power = softPower);

#Turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr,
                          networkType = "signed", 
                          TOMType = "signed", power = softPower);

colnames(TOM) =rownames(TOM) = SubGeneNames
dissTOM = 1-TOM

######

#Hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");

#Plot the resulting clustering tree (dendrogram)
par(mfrow = c(1,1))
plot(geneTree, xlab="", sub="",cex=0.3);

##############

minModuleSize = 7;

# Module identification using dynamic tree cut, you can also choose the hybrid method

dynamicMods = cutreeDynamic(dendro = geneTree,  
                            method="tree", 
                            minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, 
#                            distM = dissTOM, method="hybrid", deepSplit = 3, 
#                            pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#Get the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

## dynamicMods
##   0   1   2 
## 253 159  88

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

## dynamicColors
##      blue      grey turquoise 
##        88       253       159

plotDendroAndColors((geneTree),  dynamicColors,
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

datTrait

###################################################

table(mr.emb.TE.6.5and7.5@meta.data$embryo.stage)

m = c(rep("6.5",28),paste(d7.5.hclust.tree.group,
                      mr.emb.TE.7.5@meta.data$culture.condition))

# 2: polar

paste(d7.5.hclust.tree.group,
mr.emb.TE.6.5and7.5@meta.data$culture.condition)

mr.emb.TE.6.5and7.5.ident <- SetIdent(mr.emb.TE.6.5and7.5, 
                                      ident.use = m)

mr.emb.TE.6.5and7.5.ident <- NormalizeData(object = mr.emb.TE.6.5and7.5.ident, normalization.method = "LogNormalize", 
                    scale.factor = 10000)

mr.emb.TE.6.5and7.5.ident <- ScaleData(object = mr.emb.TE.6.5and7.5.ident, vars.to.regress = c("nUMI", "percent.mito"))


U6toU7.5polar = FindMarkers(mr.emb.TE.6.5and7.5.ident, ident.1 = "6.5", 
                           ident.2 = "2 Only_embryo")

U6toC7.5polar = FindMarkers(mr.emb.TE.6.5and7.5.ident, ident.1 = "6.5", 
                         ident.2 = "2 Co-culture")

U6toU7.5polar.U6.5high = U6toU7.5polar[(U6toU7.5polar$avg_logFC>log2(1.5) & U6toU7.5polar$p_val_adj<0.05),]

U6toU7.5polar.U6.5low = U6toU7.5polar[(U6toU7.5polar$avg_logFC<log2(1/1.5) & U6toU7.5polar$p_val_adj<0.05),]

###

U6toC7.5polar.U6.5high = U6toC7.5polar[(U6toC7.5polar$avg_logFC>log2(1.5) & U6toC7.5polar$p_val_adj<0.05),]

U6toC7.5polar.U6.5Low = U6toC7.5polar[(U6toC7.5polar$avg_logFC<log2(1/1.5) & U6toC7.5polar$p_val_adj<0.05),]

##

write.csv(U6toU7.5polar.U6.5high,file = "U6toU7.5polar.U6.5high.csv",quote=F)
write.csv(U6toU7.5polar.U6.5low,file = "U6toU7.5polar.U6.5low.csv",quote=F)
write.csv(U6toC7.5polar.U6.5high,file = "U6toC7.5polar.U6.5high.csv",quote=F)
write.csv(U6toC7.5polar.U6.5Low,file = "U6toC7.5polar.U6.5Low.csv",quote=F)


################################################


rownames(U6toC7.5polar.U6.5Low)
rownames(U6toU7.5polar.U6.5low)

intersect(rownames(U6toC7.5polar.U6.5Low),
          rownames(U6toU7.5polar.U6.5low))

###############################################



U6toUC7.5.diff.genes = (c(rownames(U6toU7.5polar.U6.5high), rownames(U6toU7.5polar.U6.5low),
                                rownames(U6toC7.5polar.U6.5high), rownames(U6toC7.5polar.U6.5Low)))

"CCR7" %in% U6toUC7.5.diff.genes 

"CCR7" %in% rownames(U6toU7.5polar.U6.5low)

"CCR7" %in% rownames(U6toC7.5polar.U6.5Low)

rownames(U6toU7.5polar.U6.5low)
rownames(U6toC7.5polar.U6.5Low)

a1= intersect(rownames(U6toU7.5polar.U6.5low),
          rownames(U6toC7.5polar.U6.5Low))

a2 = setdiff(rownames(U6toU7.5polar.U6.5low),
             rownames(U6toC7.5polar.U6.5Low))
  
a3 = setdiff(rownames(U6toC7.5polar.U6.5Low),rownames(U6toU7.5polar.U6.5low)
             )

U6toUC7.5.diff.genes.data = mr.emb.TE.6.5and7.5.ident@scale.data[match(U6toUC7.5.diff.genes,rownames(mr.emb.TE.6.5and7.5.ident@scale.data)),]

#U6toUC7.5.diff.genes.data = mr.emb.TE.6.5and7.5.ident@scale.data[match(c(a1,a2,a3),rownames(mr.emb.TE.6.5and7.5.ident@scale.data)),]





U6toUC7.5.diff.genes.data = U6toUC7.5.diff.genes.data[,c(which(mr.emb.TE.6.5and7.5.ident@ident == "6.5"), which(mr.emb.TE.6.5and7.5.ident@ident == "2 Only_embryo"),
                                                         which(mr.emb.TE.6.5and7.5.ident@ident == "2 Co-culture"))]

ccc = rep(c("orange","blue","green"),  times = c(length(which(mr.emb.TE.6.5and7.5.ident@ident == "6.5")),
                                           length(which(mr.emb.TE.6.5and7.5.ident@ident == "2 Only_embryo")),
                                           length(which(mr.emb.TE.6.5and7.5.ident@ident == "2 Co-culture"))))

heatmap.2(as.matrix( U6toUC7.5.diff.genes.data ),
          trace="none",density="none",
          Rowv = F,
          Colv = F,
          key=T,scale="row",
          dendrogram="both",
          margins=c(10.5,6),
          labRow=NA,
          col=color.palette,
          main="Day6.5 to 7.5 polar",
          ColSideColors = ccc

          
)

#############################################


U6toU7.5mural = FindMarkers(mr.emb.TE.6.5and7.5.ident, ident.1 = "6.5", 
                            ident.2 = "1 Only_embryo")

U6toC7.5mural = FindMarkers(mr.emb.TE.6.5and7.5.ident, ident.1 = "6.5", 
                            ident.2 = "1 Co-culture")
##
U6toU7.5mural.U6.high = U6toU7.5mural[(U6toU7.5mural$avg_logFC>log2(1.5) & U6toU7.5mural$p_val_adj<0.05),]

U6toU7.5mural.U6.Low = U6toU7.5mural[(U6toU7.5mural$avg_logFC<log2(1/1.5) & U6toU7.5mural$p_val_adj<0.05),]

###

U6toC7.5mural.U6.high = U6toC7.5mural[(U6toC7.5mural$avg_logFC>log2(1.5) & U6toC7.5mural$p_val_adj<0.05),]

U6toC7.5mural.U6.Low = U6toC7.5mural[(U6toC7.5mural$avg_logFC<log2(1/1.5) & U6toC7.5mural$p_val_adj<0.05),]

##

write.csv(U6toU7.5mural.U6.high,file = "U6toU7.5mural.U6.high.csv",quote=F)
write.csv(U6toU7.5mural.U6.Low,file = "U6toU7.5mural.U6.Low.csv",quote=F)
write.csv(U6toC7.5mural.U6.high,file = "U6toC7.5mural.U6.high.csv",quote=F)
write.csv(U6toC7.5mural.U6.Low,file = "U6toC7.5mural.U6.Low.csv",quote=F)




U6toUC7.5.mural.diff.genes = (c(rownames(U6toU7.5mural.U6.high), rownames(U6toU7.5mural.U6.Low),
                          rownames(U6toC7.5mural.U6.high), rownames(U6toC7.5mural.U6.Low)))


U6toUC7.5.mural.diff.genes.data = mr.emb.TE.6.5and7.5.ident@scale.data[match(U6toUC7.5.mural.diff.genes,rownames(mr.emb.TE.6.5and7.5.ident@scale.data)),]

#U6toUC7.5.diff.genes.data = mr.emb.TE.6.5and7.5.ident@scale.data[match(c(a1,a2,a3),rownames(mr.emb.TE.6.5and7.5.ident@scale.data)),]


U6toUC7.5.mural.diff.genes.data = U6toUC7.5.mural.diff.genes.data[,c(which(mr.emb.TE.6.5and7.5.ident@ident == "6.5"), which(mr.emb.TE.6.5and7.5.ident@ident == "2 Only_embryo"),
                                                         which(mr.emb.TE.6.5and7.5.ident@ident == "2 Co-culture"))]

ccc = rep(c("orange","blue","green"),  times = c(length(which(mr.emb.TE.6.5and7.5.ident@ident == "6.5")),
                                                 length(which(mr.emb.TE.6.5and7.5.ident@ident == "2 Only_embryo")),
                                                 length(which(mr.emb.TE.6.5and7.5.ident@ident == "2 Co-culture"))))

heatmap.2(as.matrix( U6toUC7.5.mural.diff.genes.data ),
          trace="none",density="none",
          Rowv = F,
          Colv = F,
          key=T,scale="row",
          dendrogram="both",
          margins=c(10.5,6),
          labRow=NA,
          col=color.palette,
          main="Day6.5 to 7.5 mural",
          ColSideColors = ccc
          
          
)

#######################################

U7.5.C7.5.polar = FindMarkers(mr.emb.TE.6.5and7.5.ident, ident.1 = "2 Co-culture", 
            ident.2 = "2 Only_embryo")

U7.5.C7.5.mural = FindMarkers(mr.emb.TE.6.5and7.5.ident, ident.1 = "1 Co-culture", 
                              ident.2 = "1 Only_embryo")


write.csv(U7.5.C7.5.polar,file = "U7.5.C7.5.polar.csv",quote=F)

write.csv(U7.5.C7.5.mural,file = "U7.5.C7.5.mural.csv",quote=F)

##################################


mr.emb.TE.UN.8.5.combat.groups@ident

te.8.5.cluseter.using.polar.genes = hclust(as.dist(1-abs(cor( te8.5,method="spearman"))), 
                                           method="ward.D")

heatmap.2(as.matrix( te8.5 ),
          trace="none",density="none",
          Colv = as.dendrogram(hclust(as.dist(1-abs(cor( te8.5,method="spearman"))), 
                                      method="ward.D")),
          Rowv = F,
          key=T,scale="row",
          dendrogram="both",
          margins=c(10.5,6),
          labRow=NA,
          col=color.palette,
          main="Day 8.5",
          ColSideColors = c("red","green","black","orange","blue")[mr.emb.TE.UN.8.5.combat.groups@ident]
          
)

##############################


Merged.TE.7.8 = MergeSeurat(mr.emb.TE.7.5,mr.emb.TE.UN.8.5)

Merged.TE.7.8 = ScaleData(object = Merged.TE.7.8, vars.to.regress = c("nUMI", "percent.mito"))

c(colnames(te7.5),colnames(te8.5))

te7.8 = Merged.TE.7.8@scale.data[na.omit(match(c(polar$polar),rownames(Merged.TE.7.8@scale.data))),]



heatmap.2(as.matrix( te7.8 ),
          trace="none",density="none",
          Colv = as.dendrogram(hclust(as.dist(1-abs(cor( te7.8,method="spearman"))), 
                                      method="ward.D")),
          Rowv = F,
          key=T,scale="row",
          dendrogram="both",
          margins=c(10.5,6),
          labRow=NA,
          col=color.palette,
          main="Day 7.5 and 8.5",
          ColSideColors = c("red","green","black","orange","blue","pink","yellow")[c(d7.5.hclust.tree.group+5,mr.emb.TE.UN.8.5.combat.groups@ident)]
          
)

###############

a1 = colnames(mr.emb@scale.data)[which(mr.emb@meta.data$Phase=="Untypical")]

a1

a1[18:36]

a2 = names(mr.emb.TE.UN.8.5.combat.groups@ident)[mr.emb.TE.UN.8.5.combat.groups@ident %in% c(5,4)]

intersect(a1[18:36],a2)

length(a1)
length(a2)

###############################



DoHeatmap(mr.emb.TE.UN.8.5.combat.groups, genes.use = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                                        "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                                        "CGA","CGB","PSG1","CSH1","HSD3B1","CYP19A1","SDC1","INHA"))
ttt = as.numeric(mr.emb.TE.UN.8.5.combat.groups@ident)
ttt[ttt==2] = 1
ttt[ttt==5] = 4
names(ttt) = names(mr.emb.TE.UN.8.5.combat.groups@ident)

mr.emb.TE.UN.8.5.combat.groups.3group = AddMetaData(mr.emb.TE.UN.8.5.combat.groups, metadata = ttt,
                                                    col.name = "TE.3.biogroups")

mr.emb.TE.UN.8.5.combat.groups.3group = SetAllIdent(mr.emb.TE.UN.8.5.combat.groups.3group,id = "TE.3.biogroups")
mr.emb.TE.UN.8.5.combat.groups.3group@ident

mr.emb.TE.UN.8.5.combat.groups.3group.markers = FindAllMarkers(mr.emb.TE.UN.8.5.combat.groups.3group,
                                                               only.pos=T )

mr.emb.TE.UN.8.5.combat.groups.3group.markers.sig = mr.emb.TE.UN.8.5.combat.groups.3group.markers[
  abs(mr.emb.TE.UN.8.5.combat.groups.3group.markers$avg_logFC) > log2(1.5) ,]

DoHeatmap(mr.emb.TE.UN.8.5.combat.groups.3group,genes.use = 
            mr.emb.TE.UN.8.5.combat.groups.3group.markers$gene)

write.csv(mr.emb.TE.UN.8.5.combat.groups.3group.markers,
            file = "mr.emb.TE.UN.8.5.combat.groups.3group.markers.csv",
            quote=F)


ST.markers = mr.emb.TE.UN.8.5.combat.groups.3group.markers.sig$gene[mr.emb.TE.UN.8.5.combat.groups.3group.markers.sig$cluster==4]


write.table(ST.markers,file = "ST.markers.txt",quote=F,row.names = F,col.names = F)
#############################################

# 用8.5天数据做WGCNA, 找ST module

#RPM.emb.TE.6.5and7.5.fil = RPM.emb.TE.6.5and7.5[apply(RPM.emb.TE.6.5and7.5,1,function(x) sum(x>1))>=40,]

mr.emb.TE.UN.8.5.combat.groups.3group = FindVariableGenes(mr.emb.TE.UN.8.5.combat.groups.3group)

PCAPlot(mr.emb.TE.UN.8.5.combat.groups.3group)

library(WGCNA)
library("flashClust")
options(stringsAsFactors = FALSE);
#enableWGCNAThreads()
disableWGCNAThreads()

#mydata=(mr.emb.TE.UN.8.5.combat.groups.3group@scale.data[match(mr.emb.TE.UN.8.5.combat.groups.3group.markers$gene,rownames(mr.emb.TE.UN.8.5.combat.groups.3group@scale.data)),])
mydata=(mr.emb.TE.UN.8.5.combat.groups.3group@scale.data[match(mr.emb.TE.UN.8.5.combat.groups.3group@var.genes,rownames(mr.emb.TE.UN.8.5.combat.groups.3group@scale.data)),])

dim(mydata)
gene.names = rownames(mydata)

mydata.trans=t(mydata);

datExpr=mydata.trans
SubGeneNames=gene.names

######
powers = c(1:30)
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,
                      powerVector = powers,corFnc = cor,
                      corOptions = list(use = 'p'),
                      networkType = "signed")


par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main = paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##################

softPower = 6;

adj= adjacency(datExpr,type = "signed", power = softPower);

#Turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr,
                          networkType = "signed", 
                          TOMType = "signed", power = softPower);

colnames(TOM) =rownames(TOM) = SubGeneNames
dissTOM = 1-TOM

######

#Hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");

#Plot the resulting clustering tree (dendrogram)
par(mfrow = c(1,1))
plot(geneTree, xlab="", sub="",cex=0.3);

##############

minModuleSize = 7;

# Module identification using dynamic tree cut, you can also choose the hybrid method

dynamicMods = cutreeDynamic(dendro = geneTree,  
                            method="tree", 
                            minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, 
#                            distM = dissTOM, method="hybrid", deepSplit = 3, 
#                            pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#Get the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

## dynamicMods
##   0   1   2 
## 253 159  88

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

## dynamicColors
##      blue      grey turquoise 
##        88       253       159

plotDendroAndColors((geneTree),  dynamicColors,
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")


Three.group.dummy = model.matrix(~  mr.emb.TE.UN.8.5.combat.groups.3group@ident  -1)


trait.exp.cor.mat = cor(Three.group.dummy,datExpr,use = "p")

library(fifer)
trait.exp.cor.mat.color.1 = number.to.colors(trait.exp.cor.mat[1,], colors = c("green", "red"), num = 100)
trait.exp.cor.mat.color.2 = number.to.colors(trait.exp.cor.mat[2,], colors = c("green", "red"), num = 100)
trait.exp.cor.mat.color.3 = number.to.colors(trait.exp.cor.mat[3,], colors = c("green", "red"), num = 100)



plotDendroAndColors(geneTree, cbind(dynamicColors,
                                    trait.exp.cor.mat.color.1,
                                    trait.exp.cor.mat.color.2,
                                    trait.exp.cor.mat.color.3
),
c("Dynamic Tree Cut", "1","2","3","4","5"
  ,"6"
),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

intersect(colnames(datExpr)[dynamicColors == "turquoise"],
          mr.emb.TE.UN.8.5.combat.groups.3group.markers$gene[mr.emb.TE.UN.8.5.combat.groups.3group.markers$cluster=="4"])

WGCNA.turquoise.gene = mr.emb.TE.UN.8.5.combat.groups.3group.markers$gene[mr.emb.TE.UN.8.5.combat.groups.3group.markers$cluster=="4"]

write.table(WGCNA.turquoise.gene,file = "WGCNA.turquoise.gene.txt",row.names = F,col.names = F,quote=F)


# 找这个 injury induced module 的 hub gene

ADJ1=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, dynamicColors)
head(Alldegrees1)

GS1=as.numeric(cor(Three.group.dummy[,3],datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, dynamicColors, mean, na.rm=T)

par(mar = c(3,7,3,1))
plotModuleSignificance(GeneSignificance,dynamicColors,las=2)



max(ModuleSignificance)


datME=moduleEigengenes(datExpr,dynamicColors)$eigengenes


colorlevels=unique(dynamicColors)
sizeGrWindow(9,6)
plotModuleSignificance(GeneSignificance,dynamicColors)
#par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))

for (i in c(1:length(colorlevels))){
  whichmodule=colorlevels[[i]];
  restrict1 = (dynamicColors==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=dynamicColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}

datKME=signedKME(datExpr, datME, outputColumnName="MM.")
# Display the first few rows of the data frame
head(datKME)

FilterGenes= abs(GS1)> 0.6 & abs(datKME$MM.turquoise)> 0.6

table(FilterGenes)

hubs.MM.turquoise = dimnames(data.frame(datExpr))[[2]][FilterGenes]

hubs.MM.turquoise

################

source("/Users/anqin/Box Sync/Projects/DRG_project/2018_01_16_merged_rnaseq/tutorialFunctions.R")

visantPrepOverall(dynamicColors, moduleColor = "turquoise", 
                  (datExpr), colnames(datExpr), 500, softPower, TRUE)


# 
# #################################################
# # day 6.5 分类
# 
# mr.day6.5 = SubsetData(object = mr, cells.use = colnames(mr@scale.data)[mr@meta.data$embryo.stage == "6.5"])
# 
# #mr.day6.5 = SubsetData(object = mr, )
# 
# mr.day6.5 <- FindVariableGenes(object = mr.day6.5, 
#                            mean.function = ExpMean, 
#                            dispersion.function = LogVMR)
# 
# 
# mr.day6.5 <- RunPCA(object = mr.day6.5,
#                 pc.genes = cell.marker,
#                 do.print = TRUE, pcs.print = 1:5, 
#                 genes.print = 5,pcs.compute = 30)
# 
# PCElbowPlot(object = mr.day6.5, num.pc = 40) + ylim(0,60) +
#   geom_hline(yintercept = 2) 
# 
# PCAPlot(object = mr.day6.5, dim.1 = 1, dim.2 = 2)
# 
# 
# FeaturePlot(object = mr.day6.5, features.plot = c("CDX2","GATA2","GATA3","KRT19","CLDN4","CLDN10","PTGES","PDGFA","DAB2"), 
#             cols.use = c("grey", "blue"), 
#             reduction.use = "pca")
# 
# FeaturePlot(object = mr.day6.5, features.plot = c("NANOG","POU5F1","PDGFRA","GDF3","LAMA4","DPPA5"), 
#             cols.use = c("grey", "blue"), 
#             reduction.use = "pca")
# 
# FeaturePlot(object = mr.day7.5, features.plot = c("CGB5","CGB","CGB8"), 
#             cols.use = c("grey", "blue"), 
#             reduction.use = "pca")
# 
# ###########################################
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #####################################################
# 
# # day 6.5 分类
# 
# mr.day7.5 = SubsetData(object = mr, cells.use = colnames(mr@scale.data)[mr@meta.data$embryo.stage == "7.5"])
# 
# #mr.day6.5 = SubsetData(object = mr, )
# 
# mr.day7.5 <- FindVariableGenes(object = mr.day7.5, 
#                                mean.function = ExpMean, 
#                                dispersion.function = LogVMR)
# 
# 
# mr.day7.5 <- RunPCA(object = mr.day7.5,
#                     pc.genes = cell.marker,
#                     do.print = TRUE, pcs.print = 1:5, 
#                     genes.print = 5,pcs.compute = 30)
# 
# PCElbowPlot(object = mr.day7.5, num.pc = 40) + ylim(0,60) +
#   geom_hline(yintercept = 2) 
# 
# PCAPlot(object = mr.day7.5, dim.1 = 1, dim.2 = 2)
# 
# mr.day7.5 <- RunTSNE(object = mr.day7.5, dims.use = 1:5, do.fast = TRUE)
# 
# 
# FeaturePlot(object = mr.day7.5, features.plot = c("CDX2","GATA2","GATA3","KRT19","CLDN4","CLDN10","PTGES","PDGFA","DAB2"), 
#             cols.use = c("grey", "blue"), 
#             reduction.use = "pca")
# 
# FeaturePlot(object = mr.day7.5, features.plot = c("NANOG","POU5F1","PDGFRA","GDF3","LAMA4","DPPA5"), 
#             cols.use = c("grey", "blue"), 
#             reduction.use = "pca")
# 
# 
# ##########################
# 
# 
# mr <- RunPCA(object = mr, pc.genes = cell.marker, do.print = TRUE, pcs.print = 1:5, 
#              genes.print = 5,pcs.compute = 40)
# PCAPlot(object = mr, dim.1 = 1, dim.2 = 2,group.by = "embryo.stage")
# 
# PCAPlot(object = mr, dim.1 = 1, dim.2 = 2,group.by = "culture.condition")
# 
# mr <- RunTSNE(mr)
# TSNEPlot(mr,group.by = "embryo.stage")
# TSNEPlot(mr,group.by = "culture.condition")
# 
# 
# FeaturePlot(object = mr, features.plot = c("NANOG","POU5F1","PDGFRA","GDF3","LAMA4","DPPA5"), 
#             cols.use = c("grey", "blue"), 
#             reduction.use = "tsne")
# 
# FeaturePlot(object = mr, features.plot = c("CDX2","GATA2","GATA3","KRT19","CLDN4","CLDN10","PTGES","PDGFA","DAB2"), 
#             cols.use = c("grey", "blue"), 
#             reduction.use = "tsne")
# 
# mr <- FindClusters(object = mr, reduction.type = "pca", dims.use = 1:5, 
#                       resolution = 0.6, print.output = 0,force.recalc = T,prune.SNN = 0.3)
# 
# TSNEPlot(object = mr,do.label = T)
# 
# 
# ###############################################
# 
# x = read.csv("cell marker.csv",header = F,stringsAsFactors = F)
# TE = x$V1[x$V2 == "TE"]
# 
# 
# 
# ####
# mr.TE = SubsetData(object = mr, ident.remove = "5")
# mr.TE = SubsetData(object = mr.TE, cells.use = colnames(mr.TE@data)[mr.TE@meta.data$cell.type!= "Endomaterial" ] )
# ####
# 
# TSNEPlot(object = mr.TE,do.label = T)
# 
# PCAPlot(object = mr.TE,do.label = T)
# 
# PCAPlot(mr.TE,group.by = "embryo.stage")
# PCAPlot(mr.TE,group.by = "culture.condition")
# 
# DimHeatmap(object = mr.TE, reduction.type = "pca", 
#            dim.use = 1:9, do.balanced = TRUE)
# 
# 
# mr.TE <- FindVariableGenes(object = mr.TE, mean.function = ExpMean, dispersion.function = LogVMR, 
#                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
# 
# mr.TE <- RunPCA(object = mr.TE, do.print = TRUE, pcs.print = 1:5,  #### 这里用了整体基因，之后可以看TE的时间上表达的差异
#              genes.print = 5,pcs.compute = 40) ### 之后fit 一个 principle curve
# 
# mr.TE <- RunTSNE(mr.TE)
# 
# PCAPlot(mr.TE,group.by = "embryo.stage")
# 
# TSNEPlot(object = mr.TE,group.by = "embryo.stage")
# 
# TSNEPlot(object = mr.TE,group.by = "culture.condition")
# 
# FeaturePlot(object = mr.TE, features.plot = c("CCR7","DLX5","GCM1","GREM2","MUC15","OVOL1"), 
#             cols.use = c("grey", "blue"), 
#             reduction.use = "tsne")
# #########################################################
# 
# library(readxl)
# 
# TESUB = read_excel("TESUB.xlsx",col_names = F)
# 
# library(sva)
# 
# modcombat = model.matrix(~1, data=mr.TE@meta.data)
# 
# combat_edata = ComBat(dat=(mr.TE@scale.data), batch=  as.numeric(mr.TE@meta.data$embryo.no)  ,
#                       mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
# 
# mr.TE.combat = mr.TE
# mr.TE.combat@scale.data = combat_edata
# 
# 
# mr.TE.combat <- FindVariableGenes(object = mr.TE.combat, mean.function = ExpMean, dispersion.function = LogVMR, 
#                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
# 
# mr.TE.combat <- RunPCA(object = mr.TE.combat, do.print = TRUE, pcs.print = 1:5,  #### 这里用了整体基因，之后可以看TE的时间上表达的差异
#                 genes.print = 5,pcs.compute = 40) ### 之后fit 一个 principle curve
# 
# mr.TE.combat <- RunTSNE(mr.TE.combat)
# 
# PCAPlot(mr.TE.combat,group.by = "embryo.stage")
# 
# 
# TSNEPlot(object = mr.TE.combat,group.by = "embryo.stage")
# 
# ####################################
# 
# mr.TE.8.5 = SubsetData(object = mr.TE, cells.use = colnames(mr.TE@data)[mr.TE@meta.data$embryo.stage == "8.5" ] )
# 
# FeaturePlot(object = mr.TE.8.5, features.plot = c("NANOG","POU5F1","PDGFRA","GDF3","LAMA4","DPPA5"), 
#             cols.use = c("grey", "blue"), 
#             reduction.use = "tsne")
# 
# FeaturePlot(object = mr.TE.8.5, features.plot = c("CDX2","GATA2","GATA3","KRT19","CLDN4","CLDN10","PTGES","PDGFA","DAB2"), 
#             cols.use = c("grey", "blue"), 
#             reduction.use = "tsne")
# 
# PCAPlot(mr.TE.8.5,group.by = "embryo.stage")
# 
# modcombat = model.matrix(~1, data=mr.TE.8.5@meta.data)
# 
# combat_edata = ComBat(dat=(mr.TE.8.5@scale.data), batch=  as.numeric(mr.TE.8.5@meta.data$embryo.no)  ,
#                       mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
# 
# mr.TE.8.5.combat = mr.TE.8.5
# mr.TE.8.5.combat@scale.data = combat_edata
# 
# 
# mr.TE.8.5.combat <- FindVariableGenes(object = mr.TE.8.5.combat, mean.function = ExpMean, dispersion.function = LogVMR, 
#                                   x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 1)
# 
# mr.TE.8.5.combat <- RunPCA(object = mr.TE.8.5.combat, do.print = TRUE, pcs.print = 1:5,  
#                        genes.print = 5,pcs.compute = 40) 
# mr.TE.8.5.combat <- RunTSNE(mr.TE.8.5.combat)
# 
# PCAPlot(mr.TE.8.5.combat,group.by = "embryo.no")
# 
# 
# TSNEPlot(object = mr.TE.8.5.combat,group.by = "embryo.no")
# 
# DimHeatmap(object = mr.TE.8.5.combat, reduction.type = "pca", 
#            dim.use = 1:9, do.balanced = TRUE)
# 
# 
# complete.cluster.8.5 = hclust(as.dist(1-abs(cor(mr.TE.8.5.combat@scale.data[match(mr.TE.8.5.combat@var.genes,rownames(mr.TE.8.5.combat@scale.data)),],method="spearman"))), 
#                           method="ward.D")
# 
# 
# plot(complete.cluster.8.5, hang = -1)
# 
# mr.TE.8.5.combat.leaf = cutree(complete.cluster.8.5,k = 5)
# 
# rect.hclust(complete.cluster.8.5,k = 5)
# 
# table(mr.TE.8.5.combat.leaf)
# 
# day8.5.exp.boxplot <- function(the.gene){
#   require(ggplot2)
# dat = data.frame(exp = mr.TE.8.5.combat@scale.data[match(the.gene,rownames(mr.TE.8.5.combat@scale.data)),],
#                  mr.TE.8.5.combat.leaf)
# 
# ggplot(dat, aes(x = factor(mr.TE.8.5.combat.leaf), y = exp)) +  geom_boxplot() + ggtitle(the.gene)
# }
# 
# 
# day8.5.exp.boxplot("CDX2")
# day8.5.exp.boxplot("GATA3")
# 
# 
# day8.5.exp.boxplot("GATA6")
# day8.5.exp.boxplot("POU5F1")
# 
# day8.5.exp.boxplot("PARD6A")
# day8.5.exp.boxplot("KRT7")
# day8.5.exp.boxplot("CGB5")
# 
# 
# 
# 
# ##############################################
# 
# 
# 
# mr.TE.7.5 = SubsetData(object = mr.TE, cells.use = colnames(mr.TE@data)[mr.TE@meta.data$embryo.stage == "7.5" ] )
# 
# 
# 
# modcombat = model.matrix(~1, data=mr.TE.7.5@meta.data)
# 
# combat_edata = ComBat(dat=(mr.TE.7.5@scale.data), batch=  as.numeric(mr.TE.7.5@meta.data$embryo.no)  ,
#                       mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
# 
# mr.TE.7.5.combat = mr.TE.7.5
# mr.TE.7.5.combat@scale.data = combat_edata
# 
# 
# mr.TE.7.5.combat <- FindVariableGenes(object = mr.TE.7.5.combat, mean.function = ExpMean, dispersion.function = LogVMR, 
#                                       x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 1)
# 
# mr.TE.7.5.combat <- RunPCA(object = mr.TE.7.5.combat, do.print = TRUE, pcs.print = 1:5,  
#                            genes.print = 5,pcs.compute = 40) 
# mr.TE.7.5.combat <- RunTSNE(mr.TE.7.5.combat)
# 
# 
# PCAPlot(mr.TE.7.5.combat,group.by = "embryo.no")
# 
# TSNEPlot(object = mr.TE.7.5.combat,group.by = "embryo.no")
# 
# 
# ###############################
# intersect(mr.TE.7.5.combat@var.genes,mr.TE.8.5.combat@var.genes)
# 
# length(mr.TE.7.5.combat@var.genes)
# length(mr.TE.8.5.combat@var.genes)
# 
# complete.cluster.7.5 = hclust(as.dist(1-abs(cor(mr.TE.7.5.combat@scale.data[match(mr.TE.7.5.combat@var.genes,rownames(mr.TE.7.5.combat@scale.data)),],method="spearman"))), 
#                           method="ward.D")
# 
# plot(complete.cluster.7.5, hang = -1)
# 
# mr.TE.7.5.combat.leaf = cutree(complete.cluster.7.5,k = 5)
# rect.hclust(complete.cluster.7.5,k = 5)
# 
# 
# day7.5.exp.boxplot <- function(the.gene){
#   require(ggplot2)
#   dat = data.frame(exp = mr.TE.7.5.combat@scale.data[match(the.gene,rownames(mr.TE.7.5.combat@scale.data)),],
#                    mr.TE.7.5.combat.leaf)
#   
#   ggplot(dat, aes(x = factor(mr.TE.7.5.combat.leaf), y = exp)) +  geom_boxplot() + ggtitle(the.gene)
# }
# 
# day7.5.exp.boxplot("KRT7")
# day7.5.exp.boxplot("CGB5")
# 
# day7.5.exp.boxplot("CDX2")
# 
# ########################################
# 
# var.genes = union(mr.TE.7.5.combat@var.genes,mr.TE.8.5.combat@var.genes)
# pos = match(var.genes,rownames(mr.TE.8.5.combat@scale.data))
# 
# plot.dat = mr.TE.8.5.combat@scale.data[pos,]
# 
# library(gplots)
# png(file = "b3.jpg",width = 2000,height = 1800)
# heatmap(as.matrix(  plot.dat   ),
#         trace = "none",density = "none",
#         col=color.palette,
#         Colv = as.dendrogram(complete.cluster.8.5),
#         scale = c("row"),
#         dendrogram = "both")
# dev.off()
# 
# GenePlot(object = mr.TE.8.5, gene1 = "KRT7", gene2 = "CGB5")
# 
# 
# 
# ##############################################
# 
# # 8.5 天的marker
# mr.TE.8.5.combat.groups = mr.TE.8.5.combat
# 
# mr.TE.8.5.combat.groups@ident = factor(mr.TE.8.5.combat.leaf)
# 
# mr.TE.8.5.markers <- FindAllMarkers(mr.TE.8.5.combat.groups, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.05)
# 
# mr.TE.8.5.markers.sig = mr.TE.8.5.markers[mr.TE.8.5.markers$p_val_adj<0.05,]
# 
# png(file = "b4.jpg",width = 2000,height = 1800)
# DoHeatmap(mr.TE.8.5.combat.groups, genes.use = mr.TE.8.5.markers.sig$gene, 
#           slim.col.label = TRUE, remove.key = TRUE)
# dev.off()
# 
# write.csv(mr.TE.8.5.markers.sig,file = "mr.TE.8.5.markers.sig.csv",quote = F)
# ############################################
# 
# 
# 
# 
# 
# 
# ###############################################
# 
# 
# 
# # 画图 RPM
# 
# grep(pattern = "^MT-", x = rownames(x = mr@data), value = F)
# 
# RC.S.noMT = RC.S[grep(pattern = "^MT-", x = rownames(RC.S), value = F),]
# 
# #mr@scale.data
# 
# palette.breaks <- seq(-1.5, 1.5, 0.1)
# color.palette = colorRampPalette(c("dodgerblue4","dodgerblue1","white","firebrick1","firebrick3"), 
#                                  space="Lab")(length(palette.breaks) - 1)
# 
# var.genes = union(mr.TE.7.5.combat@var.genes,mr.TE.8.5.combat@var.genes)
# 
# pos = match(var.genes,rownames(mr@scale.data))
# 
# plot.dat = mr@scale.data[pos,]
# 
# library(gplots)
# png(file = "b.jpg",width = 1000,height = 800)
# heatmap(as.matrix(  plot.dat   ),
#           trace = "none",density = "none",
#           col=color.palette,
#           breaks = palette.breaks,
#         ColSideColors = c("red","brown","green","cyan","blue","pink")[as.numeric(mr@ident)],
#           scale = c("row"),
#           dendrogram = "both")
# dev.off()
# 
# ###############################################
# 
# mr.TE = SubsetData(object = mr, ident.remove = "5")
# mr.TE = SubsetData(object = mr.TE, cells.use = colnames(mr.TE@data)[mr.TE@meta.data$cell.type!= "Endomaterial" ] )
# 
# modcombat = model.matrix(~1, data=mr.TE@meta.data)
# 
# combat_edata = ComBat(dat=(mr.TE@scale.data), batch=  as.numeric(mr.TE@meta.data$embryo.no)  ,
#                       mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
# 
# mr.TE.combat = mr.TE
# mr.TE.combat@scale.data = combat_edata
# 
# PCAPlot(mr.TE.combat,group.by = "embryo.no")
# 
# pos = match(var.genes,rownames(mr.TE@scale.data))
# 
# plot.dat = mr.TE@scale.data[pos,]
# 
# library(gplots)
# png(file = "b2.jpg",width = 1000,height = 800)
# heatmap(as.matrix(  plot.dat   ),
#         trace = "none",density = "none",
#         col=color.palette,
#         breaks = palette.breaks,
#         ColSideColors = c("red","brown","green","cyan","blue","pink")[as.numeric(mr.TE@ident)],
#         scale = c("row"),
#         dendrogram = "both")
# dev.off()
# ############################################
# 
# mr.TE.7.5.combat.diff1 = mr.TE.7.5.combat
# 
# mr.TE.7.5.combat.diff1@ident = factor(ifelse( substr( names(mr.TE.7.5.combat.diff1@ident),1,2)=="7.","Co","X"  ))
# 
# names(mr.TE.7.5.combat.diff1@ident) = names(mr.TE.7.5.combat@ident)
# 
# X.C.7.5.diff.genes <- FindMarkers(mr.TE.7.5.combat.diff1,ident.1 = "Co","X")
# 
# X.C.7.5.diff.genes.sig = X.C.7.5.diff.genes[X.C.7.5.diff.genes$p_val_adj < 0.05,]
# 
# png(file = "x.jpg",width = 1000,height = 800)
# DoHeatmap(mr.TE.7.5.combat.diff1, genes.use = X.C.7.5.diff.genes.sig$gene, 
#           slim.col.label = TRUE, remove.key = TRUE)
# dev.off()
# 
# 
# 
# write.csv(X.C.7.5.diff.genes.sig,file = "X.C.7.5.diff.genes.sig.csv",quote = F)
# 
# 
# 






