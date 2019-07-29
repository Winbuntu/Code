rm(list = ls())

library(data.table)
library(Seurat)

library(plyr) 

RC = as.data.frame(fread("READ_COUNT_table_lvbo_TE_gencode.txt",header = T))

#RC[which.max(RC$`6A-13_R1_val_1.fq.gz_STARAligned.out.bam`),]

rownames(RC) = RC$Geneid

#########
# 去除mt-rna, ERCC

MT.index = which(substr(RC$Geneid,1,3) == "MT-")

ERCC.index = which(substr(RC$Geneid,1,5) == "ERCC-")

RC.clean = RC[-c(MT.index,ERCC.index),-c(1:6)]

#rownames(RC.clean) = rownames(RC)[-c(MT.index,ERCC.index)]

#############

# 构建 annotation information

culture.condition = ifelse(substr(colnames(RC.clean),1,1) == "C", "Coculture" ,"OnlyEmbryo" ) 
culture.condition[which(substr(colnames(RC.clean),1,4)=="ChEC")] = "EndoC"

cell.type = ifelse(substr(colnames(RC.clean),1,4) == "ChEC", "EndoCell" ,"EmbryoCell" ) 

#write.csv(colnames(RC.clean), file = "sample.name.csv",quote=F)

embryo.day = read.csv("sample.name.csv",stringsAsFactors = F)$EmbryoDay

meta.data = data.frame(culture.condition,cell.type,embryo.day)

rownames(meta.data) = colnames(RC.clean)
#############################

# 开始构建Seurat object

mr <- CreateSeuratObject(raw.data = RC.clean, min.cells = 30, min.genes = 500, 
                         project = "Human coculture")


mr <- AddMetaData(object = mr, metadata = meta.data, col.name = c("culture.condition","cell.type","embryoday") )

VlnPlot(object = mr, features.plot = c("nGene", "nUMI"), nCol = 3)


mr <- FilterCells(object = mr, subset.names = c("nGene","nUMI"), 
                  low.thresholds = c(5000,0), high.thresholds = c(Inf,  2000000))


mr <- NormalizeData(object = mr, normalization.method = "LogNormalize", 
                    scale.factor = 10000)

mr <- FindVariableGenes(object = mr, mean.function = ExpMean, dispersion.function = LogVMR,do.plot = F)

length(x = mr@var.genes)

date()
mr <- ScaleData(object = mr, vars.to.regress = c("nUMI"))
date()


###############################

# 整体的Seurat object 构建完成

#############################

mr <- RunPCA(object = mr, pc.genes = mr@var.genes, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5,pcs.compute = 40)


PCElbowPlot(object = mr)

mr <- JackStraw(object = mr, num.replicate = 100,num.pc = 20)

JackStrawPlot(object = mr, PCs = 1:20) # 这里，选择PC12

# 4.54 * 3.58

PCAPlot(object = mr, dim.1 = 1, dim.2 = 2,group.by = "embryo.day", pt.size=1)

PCAPlot(object = mr, dim.1 = 1, dim.2 = 2,group.by = "cell.type") # 分离Endocell 的PCA

PCAPlot(object = mr, dim.1 = 1, dim.2 = 2,group.by = "culture.condition") # 分离Endocell 的PCA


pca.cell.type = ifelse(mr@dr$pca@cell.embeddings[,1] < 15, "EmbryoCell" ,"EndoCell" ) 

mr <- AddMetaData(object = mr, metadata = pca.cell.type, col.name = c("pca.cell.type") )

PCAPlot(object = mr, dim.1 = 1, dim.2 = 2,group.by = "pca.cell.type") # 分离Endocell 的PCA, 使用pca results

#################################
# 只选择embryo cell进行分析
# 去掉Day 14的细胞。没有卵用。


mr.Embryo = SubsetData(mr,cells.use = rownames(mr@meta.data)[which(mr@meta.data$pca.cell.type == "EmbryoCell")],
                       do.center = T,do.scale = T)

mr.Embryo = SubsetData(mr.Embryo,cells.use = rownames(mr.Embryo@meta.data)[which(mr.Embryo@meta.data$embryo.day != 14)],
                       do.center = T,do.scale = T)

mr.Embryo = FindVariableGenes(object = mr.Embryo, mean.function = ExpMean, dispersion.function = LogVMR,do.plot = F)

mr.Embryo <- RunPCA(object = mr.Embryo, pc.genes = mr.Embryo@var.genes, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5,pcs.compute = 40)

PCAPlot(object = mr.Embryo, dim.1 = 1, dim.2 = 2,group.by = "embryo.day") # 所有 embryo cell 的PCA。这里面包含TE/ICM

VlnPlot(mr.Embryo)

##############################

# 只保留TE细胞


cell.marker = read.csv("../cell marker.csv",header = F,stringsAsFactors = F)$V1

mr.Embryo.300markerPCA <- RunPCA(object = mr.Embryo, pc.genes = cell.marker, do.print = TRUE, pcs.print = 1:5, 
                 genes.print = 5,pcs.compute = 40)

#PCElbowPlot(object = mr.emb)

#mr.emb <- JackStraw(object = mr.emb, num.replicate = 100, do.print = FALSE,num.pc = 40)

#JackStrawPlot(object = mr.emb, PCs = 1:40) # 这里，选择PC10

PCAPlot(object = mr.Embryo.300markerPCA, dim.1 = 1, dim.2 = 2,group.by = "embryo.day") # 这里面只有embryo cells


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

cell.marker.table = read.csv("../cell marker.csv",header = F,stringsAsFactors = F)

mr.Embryo.300markerPCA <- CellLineageScoring(mr.Embryo.300markerPCA,TE.genes = cell.marker.table$V1[cell.marker.table$V2 == "TE"], # G2M 对应 TE
                             PE.genes =  cell.marker.table$V1[cell.marker.table$V2 == "PE"],
                             EPI.genes =  cell.marker.table$V1[cell.marker.table$V2 == "EPI"],set.ident = T)

PCAPlot(object = mr.Embryo.300markerPCA, dim.1 = 1, dim.2 = 2,group.by = "Phase")
PCAPlot(object = mr.Embryo.300markerPCA, dim.1 = 1, dim.2 = 2)

PCAPlot(object = mr.Embryo.300markerPCA, dim.1 = 1, dim.2 = 2,group.by = "embryo.day")

table(mr.Embryo.300markerPCA@meta.data$Phase)

mr.TE = SubsetData(mr.Embryo.300markerPCA, cells.use = rownames(mr.Embryo.300markerPCA@meta.data)[which(mr.Embryo.300markerPCA@meta.data$Phase != "EPI" & 
                                                                                                          mr.Embryo.300markerPCA@meta.data$Phase != "PE")],
                   do.center = T,do.scale = T)


mr.TE = FindVariableGenes(object = mr.TE, mean.function = ExpMean, dispersion.function = LogVMR,do.plot = F)

mr.TE <- RunPCA(object = mr.TE, pc.genes = mr.TE@var.genes, do.print = TRUE, pcs.print = 1:5, 
                    genes.print = 5,pcs.compute = 40)

PCAPlot(object = mr.TE, dim.1 = 1, dim.2 = 2,group.by = "embryo.day") # 所有 embryo cell 的PCA。这里面包含TE/ICM


TSNEPlot(object = mr.TE, dim.1 = 1, dim.2 = 2,group.by = "embryo.day") # 所有 embryo cell 的PCA。这里面包含TE/ICM

PCAPlot(object = mr.TE, dim.1 = 1, dim.2 = 2)
PCAPlot(object = mr.TE, dim.1 = 1, dim.2 = 2,group.by = "Phase") # 所有 embryo cell 的PCA。这里面包含TE/ICM

mr.TE.1 = SubsetData(mr.TE,ident.use = 1,do.center = T,do.scale = T)

mr.TE.1 = FindVariableGenes(mr.TE.1,do.plot=F)

mr.TE.1@var.genes

mr.TE <- FindClusters(object = mr.TE, reduction.type = "pca", dims.use = 1:10, 
                                 resolution = 0.6, print.output = 0, save.SNN = TRUE)

PCAPlot(object = mr.TE, dim.1 = 1, dim.2 = 2) 


mr.TE <- JackStraw(object = mr.TE, num.replicate = 100, num.pc = 40)

JackStrawPlot(object = mr.TE, PCs = 1:40) # 这里，选择PC8


mr.TE <- RunTSNE(mr.TE,dims.use = 1:40)

TSNEPlot(object = mr.TE, dim.1 = 1, dim.2 = 2,group.by = "embryo.day") # 所有 embryo cell 的PCA。这里面包含TE/ICM
TSNEPlot(object = mr.TE, dim.1 = 1, dim.2 = 2,group.by = "culture.condition") # 所有 embryo cell 的PCA。这里面包含TE/ICM

FeaturePlot(object = mr.TE, features.plot = c("LRP2","TPM1","DPPA3","FOS","ERVW-1","TET1","TET3","TET2","DPPA3","TBX3","DNMT3B",
                                              "DNMT1","DNMT3A"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca")

VlnPlot(object = mr.day9.10, features.plot = c("LRP2","TPM1","DPPA3","FOS","ERVW-1","TET1","TET3","TET2","DPPA3","TBX3","DNMT3B",
                                              "DNMT1","DNMT3A","ZFP42"))

VlnPlot(object = mr.day9.10, features.plot = c("ACTB","TET2","TET1","TET3"))

PCAPlot(object = mr.TE, dim.1 = 1, dim.2 = 2) # 所有 embryo cell 的PCA。这里面包含TE/ICM
TSNEPlot(object = mr.TE, dim.1 = 1, dim.2 = 2) # 所有 embryo cell 的PCA。这里面包含TE/ICM


mr.TE.marker = FindAllMarkers(mr.TE,only.pos = T)

mr.TE.marker.filter = mr.TE.marker[mr.TE.marker$avg_logFC > 1,]


dim(mr.TE.marker.filter)
  
mr.TE.marker[mr.TE.marker$cluster==5,]

mr.TE.2.compare = FindMarkers(mr.TE,ident.1 = 5,ident.2 = 4)

mr.TE.2.compare[mr.TE.2.compare$avg_logFC < -1.5,]

mr.TE.2.compare[mr.TE.2.compare$avg_logFC > 1.5,]

mr.TE.cluster.average.exp = AverageExpression(mr.TE)

library(corrplot)
M =cor(mr.TE.cluster.average.exp[rownames(mr.TE.cluster.average.exp) %in% mr.TE@var.genes,],method = "spearman")
corrplot(M, method = "circle")

corrplot(M, order = "hclust", addrect = 4)

###################################################################

# endogenous virus

RC.virus = as.data.frame(fread("Lvbo_RC_repbase.tsv",header = T))

rownames(RC.virus) = RC.virus$Geneid

RC.virus.clean = RC.virus[,-c(1:6)]

colnames(RC.virus.clean) =  substr(colnames(RC.virus.clean),1,DescTools::StrPos(colnames(RC.virus.clean),pattern = "_")-1)

mr.TE.mRNA.RC = mr.TE@raw.data[names(mr.TE@ident)]

RC.virus.clean.keep = RC.virus.clean[,match(substr(names(mr.TE@ident),1,DescTools::StrPos(names(mr.TE@ident),pattern = "_")-1),colnames(RC.virus.clean))]


RPM.virus.clean.keep = sweep(  RC.virus.clean.keep,2,  colSums( mr.TE.mRNA.RC),"/" ) * 1000000

RPM.virus.clean.keep.filtered = RPM.virus.clean.keep[(rowMeans(RPM.virus.clean.keep) > 5),]


color.palette = colorRampPalette(c("green","black","orangered"), space="Lab")

library(gplots)
pdf(file = "a.pdf",width = 30,height = 10)
heatmap.2(as.matrix(  log2(RPM.virus.clean.keep.filtered+1)),
          trace="none",density="none",
          Colv = T,
          Rowv = T,
          key=T,scale="row",
          dendrogram="none",
          col=color.palette,
          margins=c(10.5,6),
          labRow=NA,
          main="Tissue specific genes"
          
)
dev.off()

ggplot2.dat = data.frame(mr.TE@dr$pca@cell.embeddings[,c(1,2)],t(RPM.virus.clean.keep.filtered))

sort(rowMeans(RPM.virus.clean.keep.filtered),decreasing = T)

ggplot(ggplot2.dat, aes(x=PC1, y=PC2,colour = L1HS.L1.Homo_sapiens)) + geom_point()


#######################

# 找出9-10天两组，然后看差异


CT.samples = rownames(mr.TE@meta.data)[which((mr.TE@meta.data$embryo.day == 9 | mr.TE@meta.data$embryo.day == 10 ) & mr.TE@dr$pca@cell.embeddings[,2] > 0)]

ST.samples = rownames(mr.TE@meta.data)[which((mr.TE@meta.data$embryo.day == 9 | mr.TE@meta.data$embryo.day == 10 ) & mr.TE@dr$pca@cell.embeddings[,2] < 0)]


mr.day9.10 = SubsetData(mr.TE,cells.use = c(CT.samples,ST.samples),do.center = T,do.scale = T)

iii = names(mr.day9.10@ident)

iii[iii %in% CT.samples] = "CT"
iii[iii %in% ST.samples] = "ST"

names(iii) = names(mr.day9.10@ident)
  
mr.day9.10 = AddMetaData(mr.day9.10, metadata = iii,col.name = "CTST")
mr.day9.10 = SetAllIdent(mr.day9.10,"CTST")

ST.CT.DEG = FindMarkers(mr.day9.10,ident.1 = "ST",ident.2 = "CT")

ST.CT.DEG[ST.CT.DEG$avg_logFC > 1.5,]

ST.CT.DEG[ST.CT.DEG$avg_logFC < -1.5,]



#################

#  看HERV 在两组间的差异

CT.sample.name = substr(names(mr.day9.10@ident)[which(mr.day9.10@ident == "CT")],1,DescTools::StrPos(names(mr.day9.10@ident)[which(mr.day9.10@ident == "CT")],pattern = "_")-1)

ST.sample.name = substr(names(mr.day9.10@ident)[which(mr.day9.10@ident == "ST")],1,DescTools::StrPos(names(mr.day9.10@ident)[which(mr.day9.10@ident == "ST")],pattern = "_")-1)

write.table(CT.samples, file = "CT.bams.txt",quote = F, row.names = F, col.names = F)
write.table(ST.samples, file = "ST.bams.txt",quote = F, row.names = F, col.names = F)


match(CT.sample.name,colnames(RPM.virus.clean.keep.filtered))
match(ST.sample.name,colnames(RPM.virus.clean.keep.filtered))

p.val = c()
ave.CT = c()
ave.ST =c()

for(i in c(1:(dim(RPM.virus.clean.keep.filtered)[1]))){
  
  print(i)
  p.val[i] = wilcox.test(as.numeric(RPM.virus.clean.keep.filtered[i,match(CT.sample.name,colnames(RPM.virus.clean.keep.filtered))]),
            as.numeric(RPM.virus.clean.keep.filtered[i,match(ST.sample.name,colnames(RPM.virus.clean.keep.filtered))]) )$p.value
  
  ave.CT[i] = mean(as.numeric(RPM.virus.clean.keep.filtered[i,match(CT.sample.name,colnames(RPM.virus.clean.keep.filtered))]))
  ave.ST[i]= mean(as.numeric(RPM.virus.clean.keep.filtered[i,match(ST.sample.name,colnames(RPM.virus.clean.keep.filtered))]))
  
}

test.herv.df = data.frame(p.val,ave.CT,ave.ST,p.adjust(p.val))

rownames(test.herv.df)=rownames(RPM.virus.clean.keep.filtered)

test.herv.df[test.herv.df$p.adjust.p.val.<0.01,]

#####################################################

mr.TE.cc =  CellCycleScoring(object = mr.TE, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes, 
                             set.ident = TRUE)


#mr.TE.cc = ScaleData(object = mr.TE.cc, vars.to.regress = c("S.Score", "G2M.Score"), 
#          display.progress = T)


#mr.TE.cc <- RunPCA(object = mr.TE.cc, pc.genes = mr.TE.cc@var.genes, genes.print = 10)

PCAPlot(object = mr.TE.cc, dim.1 = 1, dim.2 = 2,group.by = "Phase") # 所有 embryo cell 的PCA。这里面包含TE/ICM

#########################

PCAPlot(object = mr.TE.cc, dim.1 = 1, dim.2 = 2,group.by = "embryo.day") # 所有 embryo cell 的PCA。这里面包含TE/ICM



FeaturePlot(object = mr.TE, features.plot = c("LRP2","TPM1","DPPA3","FOS","ERVW-1","TET1","TET3","TET2","DPPA3","TBX3","DNMT3B",
                                              "DNMT1","DNMT3A","CGB8","CGB3","ZFP42","MKI67","OVOL1","CDH1","GCM1","PITX1","SP6"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca",pt.size = 2)


##############################################

write.table(rownames(mr.TE.2.compare[mr.TE.2.compare$avg_logFC < -1,]), file =  "CT.txt", quote=F,row.names = F)
write.table(rownames(mr.TE.2.compare[mr.TE.2.compare$avg_logFC > 1,]), file =  "ST.txt", quote=F,row.names = F)

mr.TE.2.compare[mr.TE.2.compare$avg_logFC > 1,]


###########################


mr.TE.day8 = SubsetData(mr.TE,cells.use = rownames(mr.TE@meta.data)[which(mr.TE@meta.data$embryo.day == 8)],
                       do.center = T,do.scale = T)

mr.TE.day8 =  CellCycleScoring(object = mr.TE.day8, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes, 
                             set.ident = TRUE)


mr.TE.day8.cc = ScaleData(object = mr.TE.day8, vars.to.regress = c("S.Score", "G2M.Score"), 
          display.progress = T)

mr.TE.day8.cc = FindVariableGenes(mr.TE.day8.cc,do.plot = F,y.cutoff =3)

mr.TE.day8.cc@var.genes

intersect(rownames(mr.TE.2.compare[mr.TE.2.compare$avg_logFC < -1.5,]), mr.TE.day8.cc@var.genes)

intersect(rownames(mr.TE.2.compare[mr.TE.2.compare$avg_logFC > 1.5,]), mr.TE.day8.cc@var.genes)

HTF = read.table("hTF.txt",stringsAsFactors = F)$V1

intersect(HTF, intersect(rownames(mr.TE.2.compare[mr.TE.2.compare$avg_logFC < -1.5,]), mr.TE.day8.cc@var.genes))
intersect(HTF, intersect(rownames(mr.TE.2.compare[mr.TE.2.compare$avg_logFC > 1.5,]), mr.TE.day8.cc@var.genes))

intersect(HTF, (rownames(mr.TE.2.compare[mr.TE.2.compare$avg_logFC > 1.5,])))

intersect(HTF, (rownames(mr.TE.2.compare[mr.TE.2.compare$avg_logFC < -1,])))

intersect(rownames(mr.TE.2.compare[mr.TE.2.compare$avg_logFC > 1.5,]), mr.TE.day8.cc@var.genes)

GenePlot(mr.TE, gene1 = "TBX3",gene2 = "DNMT3B")

as.matrix(mr.TE@data)[match("DNMT3B", rownames(as.matrix(mr.TE@data))),]


################################################

mr.TE@raw.data[match("TBX3",row.names(mr.TE@raw.data)),]

#################################################################

library(ggpubr)
data("ToothGrowth")
bxp <- ggboxplot(ToothGrowth, x = "dose", y = "len",
                 color = "dose", palette = "jco")
bxp
# Dot plot (dp)
dp <- ggdotplot(ToothGrowth, x = "dose", y = "len",
                color = "dose", palette = "jco", binwidth = 1)
dp
ggarrange( plotlist = list(bxp, dp), 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)


#########################



gene.length = RC$Length[match(rownames(mr.TE@raw.data),RC$Geneid)]

mr.RPM = sweep(x = mr.TE@raw.data,MARGIN = 2, STATS = colSums(mr.TE@raw.data),FUN = "/" ) * 1000000

mr.RPKM = sweep(x = mr.RPM,MARGIN = 1, STATS = gene.length,FUN = "/" ) * 1000


gene.list = c("LRP2","TPM1","DPPA3","FOS","ERVW-1","TET1","TET3","TET2","DPPA3","TBX3","DNMT3B",
              "DNMT1","DNMT3A","CGB8","CGB3","ZFP42","MKI67","OVOL1","CDH1","GCM1","PITX1")

plot.list=list()
counter = 1


for(the.gene in gene.list){
  
  
  pca.plot = data.frame(RPKM = as.numeric(mr.RPKM[match(the.gene, rownames(mr.RPKM)), match(names(mr.TE@ident), colnames(mr.RPKM)) ]),
                        PC1 = mr.TE@dr$pca@cell.embeddings[,1],
                        PC2 = mr.TE@dr$pca@cell.embeddings[,2]
  )
  
  p = ggplot(pca.plot, aes(x=PC1, y=PC2,colour=log10(RPKM+1) )) +
    geom_point()   + scale_color_gradient(low = "grey", high = "red") + ggtitle(the.gene)

  plot.list[[counter]] = p
    
    
  counter = counter +1
}


merged.plot = ggarrange(plotlist = plot.list,ncol = 3,nrow = length(gene.list)/3+1)
print(merged.plot)

ggsave(merged.plot,filename = "a.pdf",width = 15,height = 20)

############################################3

gene.length = RC$Length[match(rownames(mr.TE@raw.data),RC$Geneid)]

mr.RPM = sweep(x = mr.TE@raw.data,MARGIN = 2, STATS = colSums(mr.TE@raw.data),FUN = "/" ) * 1000000

mr.RPKM = sweep(x = mr.RPM,MARGIN = 1, STATS = gene.length,FUN = "/" ) * 1000


gene.list = c("LRP2","TPM1","DPPA3","FOS","ERVW-1","TET1","TET3","TET2","DPPA3","TBX3","DNMT3B",
              "DNMT1","DNMT3A","CGB8","CGB3","ZFP42","MKI67","OVOL1","CDH1","GCM1","PITX1")

plot.list=list()
counter = 1


for(the.gene in gene.list){
  
  
  pca.plot = data.frame(RPKM = as.numeric(mr.RPKM[match(the.gene, rownames(mr.RPKM)), match(names(mr.day9.10@ident), colnames(mr.RPKM)) ]),
                        type = mr.day9.10@ident
  )
  
  p = ggplot(pca.plot, aes(x=type, y=log10(RPKM+1), fill = type)) + 
    geom_violin() + ggtitle(the.gene) + geom_jitter()
  
  plot.list[[counter]] = p
  
  
  counter = counter +1
}


merged.plot = ggarrange(plotlist = plot.list,ncol = 3,nrow = length(gene.list)/3+1)
print(merged.plot)

ggsave(merged.plot,filename = "b.pdf",width = 15,height = 20)




PCAPlot(mr.TE)

PCAPlot(mr.TE,group.by="embryo.day")

########################################################
mr.TE.fullPCA = RunPCA(object = mr.TE, pcs.compute = 40,pc.genes = rownames(mr.TE@scale.data), genes.print = 10)

mr.TE.fullPCA = RunTSNE(object = mr.TE.fullPCA, dims.use = 1:40, do.fast = TRUE)


mr.TE.fullPCA = FindClusters(mr.TE.fullPCA)

PCAPlot(mr.TE)
FeaturePlot(object = mr.TE, features.plot = c("TBX3"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca",pt.size=2)


PCAPlot(mr.day9.10,group.by="embryo.day")


mr.day9.10 =  FindClusters(mr.day9.10)

PCAPlot(mr.day9.10)

TSNEPlot(mr.TE.fullPCA,group.by="embryo.day")


FeaturePlot(object = mr.TE.fullPCA, features.plot = c("LRP2","TPM1","DPPA3","FOS","ERVW-1","TET1","TET3","TET2","DPPA3","TBX3","DNMT3B",
                                              "DNMT1","DNMT3A","CGB8","CGB3","ZFP42","MKI67","OVOL1","CDH1","GCM1","PITX1","SP6"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca",pt.size = 2)




FeaturePlot(object = mr.TE, features.plot = c("HLA-G","CDH1"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca",pt.size = 2)


PCAPlot(mr.TE.cc, group.by = "Phase",pt.size = 2)

VlnPlot(mr.TE.cc, features.plot = c("HLA-G","CDH1"),group.by = "Phase")


GenePlot(mr.TE.fullPCA.day9.10,gene1 = "HLA-G",gene2 = "CDH1")



mr.TE.fullPCA.day9.10 = SubsetData(mr.TE.fullPCA,cells.use = c(CT.samples,ST.samples),do.center = T,do.scale = T)


mr.TE.fullPCA.day8 = SubsetData(mr.TE.fullPCA, cells.use =rownames( mr.TE.fullPCA@meta.data)[ mr.TE.fullPCA@meta.data$embryo.day==8]
                                ,do.center = T,do.scale = T)

mr.TE.fullPCA.day8 = FindVariableGenes(mr.TE.fullPCA.day8,do.plot = F,y.cutoff =2)

mr.TE.fullPCA.day7 = SubsetData(mr.TE.fullPCA, cells.use =rownames( mr.TE.fullPCA@meta.data)[ mr.TE.fullPCA@meta.data$embryo.day==7]
                                ,do.center = T,do.scale = T)

mr.TE.fullPCA.day7 = FindVariableGenes(mr.TE.fullPCA.day7,y.cutoff =2,do.plot = F)

p=intersect(mr.TE.fullPCA.day7@var.genes, rownames(mr.TE.2.compare[ mr.TE.2.compare$avg_logFC > 1.5,]) )
p
intersect(p,HTF)

intersect(intersect(mr.TE.fullPCA.day8@var.genes, rownames(mr.TE.2.compare[ abs(mr.TE.2.compare$avg_logFC) > 1.5 & mr.TE.2.compare$p_val_adj < 0.05,])),HTF)

intersect(intersect(mr.TE.fullPCA.day7@var.genes, rownames(mr.TE.2.compare[ abs(mr.TE.2.compare$avg_logFC) > 1.5 & mr.TE.2.compare$p_val_adj < 0.05,])),HTF)

rownames(mr.TE.2.compare[ abs(mr.TE.2.compare$avg_logFC) > 1.5 & mr.TE.2.compare$p_val_adj < 0.05,])

(rownames(mr.TE.2.compare[ (mr.TE.2.compare$avg_logFC)< -1.5  & mr.TE.2.compare$p_val_adj < 0.05,]))

(rownames(mr.TE.2.compare[ (mr.TE.2.compare$avg_logFC) > 1.5 & mr.TE.2.compare$p_val_adj < 0.05,]))


FeaturePlot(mr.TE, features.plot = (rownames(mr.TE.2.compare[ (mr.TE.2.compare$avg_logFC)< -1.5  & mr.TE.2.compare$p_val_adj < 0.05,])), 
            reduction.use = "pca",pt.size=2, cols.use = c("grey","blue"))


FeaturePlot(mr.TE, features.plot = (rownames(mr.TE.2.compare[ (mr.TE.2.compare$avg_logFC) > 2  & mr.TE.2.compare$p_val_adj < 0.05,])), 
            reduction.use = "pca",pt.size=2, cols.use = c("grey","blue"))



mr.TE.fullPCA.day9.10
GenePlot(mr.TE.fullPCA.day9.10, gene1 = "TBX3",gene2 = "TET2")

GenePlot(mr.TE.fullPCA.day9.10, gene1 = "DNMT3B",gene2 = "MKI67")

#######################################

mr.TE.fullPCA.noday9.10 = SubsetData(mr.TE.fullPCA, 
                                     cells.use =rownames( mr.TE.fullPCA@meta.data)[ !(mr.TE.fullPCA@meta.data$embryo.day    %in% c(9,10) )]
                                     ,do.center = T,do.scale = T)

mr.TE.fullPCA.noday9.10 = FindVariableGenes(mr.TE.fullPCA.noday9.10, do.plot = F)


mr.TE.fullPCA.using.early.genes = mr.TE.fullPCA

mr.TE.fullPCA.using.early.genes = RunPCA(mr.TE.fullPCA.using.early.genes, pc.genes =mr.TE.fullPCA.noday9.10@var.genes )


PCAPlot(mr.TE.fullPCA.using.early.genes, group.by = "embryo.day")


VlnPlot(mr.day9.10, features.plot = c("TBX3","ACTB","TET2","TET1","TET3","DNMT3A","DNMT3B","DNMT1","DPPA3"), group.by = "CTST")

VlnPlot(mr.day9.10, features.plot = c("GCM1"), group.by = "CTST")


##################


mr.day9.10.rePCA = FindVariableGenes(mr.day9.10, y.cutoff  =2, do.plot = F )
length(mr.day9.10.rePCA@var.genes)

mr.day9.10.rePCA = RunPCA(mr.day9.10.rePCA, pc.genes =mrs )

PCAPlot(mr.day9.10.rePCA)

FeaturePlot(object = mr.day9.10.rePCA, features.plot = c("HLA-G","CDH1","ERVW-1"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca",pt.size = 2)

EVT_CT = names(mr.day9.10.rePCA@dr$pca@cell.embeddings[,1] > -3)[mr.day9.10.rePCA@dr$pca@cell.embeddings[,1] > -3]

mrs = read.csv("../EVT_CT_ST_markers.csv",header=F,stringsAsFactors = F)$V1

mr.day9.10.rePCA.EVT_CT = SubsetData(mr.day9.10.rePCA,cells.use = EVT_CT)

PCAPlot(mr.day9.10.rePCA.EVT_CT)

mr.day9.10.rePCA.EVT_CT = FindVariableGenes(mr.day9.10.rePCA.EVT_CT, do.plot = F )

mr.day9.10.rePCA.EVT_CT@var.genes

mr.day9.10.rePCA.EVT_CT = RunPCA(mr.day9.10.rePCA.EVT_CT )

PCAPlot(mr.day9.10.rePCA.EVT_CT)

FeaturePlot(object = mr.day9.10.rePCA.EVT_CT, features.plot = c("HLA-G","CDH1","PLAC8"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca",pt.size = 2)


FeaturePlot(object = mr.TE, features.plot = c("XIST","DDX3Y"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca",pt.size = 2)

####################################################

PCAPlot(mr.Embryo.300markerPCA)

mr.Embryo.300markerPCA.cell.label = as.character(mr.Embryo.300markerPCA@meta.data$Phase)
             
mr.Embryo.300markerPCA.cell.label[mr.Embryo.300markerPCA.cell.label=="Untypical"] = "TE"

mr.Embryo.300markerPCA.cell.label.df = data.frame(lineage = mr.Embryo.300markerPCA.cell.label)

rownames(mr.Embryo.300markerPCA.cell.label.df) = rownames(mr.Embryo.300markerPCA@meta.data)
  
mr.Embryo.300markerPCA = AddMetaData(mr.Embryo.300markerPCA,metadata = mr.Embryo.300markerPCA.cell.label.df, col.name = "lineage")


PCAPlot(mr.Embryo.300markerPCA, group.by = "lineage")

PCAPlot(mr.Embryo.300markerPCA, group.by = "embryo.day")



#####################################################

# figure 1C

mr.Embryo.300markerPCA.TE = SubsetData(mr.Embryo.300markerPCA,subset.name = "lineage", accept.value = "TE")

PCAPlot(mr.Embryo.300markerPCA.TE, group.by = "embryo.day")

PCAPlot(mr.Embryo.300markerPCA.TE, group.by = "embryo.day")

FeaturePlot(object = mr.Embryo.300markerPCA.TE, features.plot = c("GATA2","GATA3"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca",pt.size = 2, nCol = 1, no.legend = F)


# not maintained markers

FeaturePlot(object = mr.Embryo.300markerPCA.TE, features.plot = c("DAB2","PTGES","TGFBR3","PDGFA"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca",pt.size = 1, nCol = 4, no.legend = F)

VlnPlot(mr.Embryo.300markerPCA.TE,features.plot = c("DAB2","PTGES","TGFBR3","PDGFA"), group.by = "embryo.day",
        nCol = 4)


# maintained_TE markers

VlnPlot(mr.TE, features.plot = c("KRT19","KRT18","KRT8","TFRC"), group.by = "embryo.day",nCol=4)

FeaturePlot(object = mr.Embryo.300markerPCA.TE, features.plot = c("KRT19","KRT18","KRT8","TFRC"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca",pt.size = 1, nCol = 4, no.legend = F)


###################################################



# Figure1 

VlnPlot(mr.TE, features.plot = c("nGene"), group.by = "embryo.day")

nGene.density.df = mr.TE@meta.data[,c("nGene","embryo.day")]

ggplot(nGene.density.df, aes(x=nGene,color = embryo.day ,fill=embryo.day))+
  geom_density(alpha = 0.3) + theme_bw()

################

mr.cell.type = SetAllIdent(mr, id = "pca.cell.type")

mr.cell.type.marker = FindMarkers(mr.cell.type,ident.1 ="EmbryoCell", ident.2 = "EndoCell"  )

mr.cell.type.marker.sig = mr.cell.type.marker[mr.cell.type.marker$p_val_adj<0.05 & abs(mr.cell.type.marker$avg_logFC)>1,]

library(gplots)
palette.breaks <- seq(-1.5, 1.5, 0.1)
color.palette = colorRampPalette(c("dodgerblue4","dodgerblue1","white","firebrick1","firebrick3"), 
                                 space="Lab")(length(palette.breaks) - 1)

heatmap.2(as.matrix( mr.cell.type@scale.data[rownames(mr.cell.type.marker.sig),] ),trace = "none",density = "none",
          Colv = NA,
          Rowv  = NA,
          #ColSideColors = c("grey","red","yellow","green")[ as.numeric(Xy.cc.noD26@ident) ] ,
          col=color.palette,
          breaks = palette.breaks,
          scale = c("row"),
          
          )

##################################

mr.embryo.dt = mr.Embryo.300markerPCA@raw.data[,colnames(mr.Embryo.300markerPCA@data)]

mr.embryo.dt.rpm = sweep(mr.embryo.dt,MARGIN = 2,STATS = colSums(mr.embryo.dt),FUN = "/") * 1000000

FeaturePlot(mr.Embryo.300markerPCA, features.plot = c("PDGFRA","DPPA5","GATA3"),reduction.use = "pca", nCol = 3, 
            cols.use = c("grey","blue"), pt.size = 2)


FeaturePlot(mr.Embryo.300markerPCA, features.plot = c( "PDGFRA","GATA2","GATA3", "NANOG", "POU5F1","GATA6"),reduction.use = "pca", nCol = 3, 
            cols.use = c("grey","blue"), pt.size = 2)

######################## lineage marker expression

FeaturePlot(mr.Embryo.300markerPCA, features.plot = c( "TDGF1","PDGFRA","GATA2","GATA3"),reduction.use = "pca", 
            nCol = 4, 
            cols.use = c("grey","blue"), pt.size = 2)

FeaturePlot(mr.Embryo.300markerPCA, features.plot = c( "TDGF1","PDGFRA","GATA2","GATA3"),reduction.use = "pca", 
            nCol = 4, 
            cols.use = c("grey","blue"), pt.size = 2)

"PTGES, EMP2, TGFBR3, and PDGFA"

FeaturePlot(mr.Embryo.300markerPCA, features.plot = c("GATA2","GATA3","DAB2"),reduction.use = "pca", 
            nCol = 4, 
            cols.use = c("grey","blue"), pt.size = 2)

######################################
# identify maintained TE markers 

mr.TE.id.day = SetAllIdent(mr.TE,id = "embryo.day")

TE.id.day.exp = AverageExpression(mr.TE.id.day)

TE.id.day.exp[c("GATA3","GATA2"),]


TE.cell.markers.100 = head(readxl::read_excel("mmc2-2_cell_ICMTEmarker.xlsx",sheet = "TE_rank",col_names = F)$X__1,100)

which(apply(TE.id.day.exp[TE.cell.markers.100,],1,min) > 5)

FeaturePlot(mr.Embryo.300markerPCA, features.plot = c("KRT19","GATA3","KRT18","GATA2","KRT8","TFRC"),reduction.use = "pca", 
            nCol = 4, 
            cols.use = c("grey","blue"), pt.size = 2)

VlnPlot(mr.TE, features.plot = c("KRT19","GATA3","KRT18","GATA2","KRT8","TFRC"), group.by = "embryo.day")

VlnPlot(mr.TE, features.plot = c("KRT19","KRT18","KRT8","GATA2","GATA3"), group.by = "embryo.day", nCol = 3)

####################



the.gene.vioplot = "GATA3"

the.rpkm = as.numeric(mr.RPKM[the.gene.vioplot, match(rownames(mr.TE@meta.data), colnames(mr.RPKM))  ])

vio.df = data.frame(exp = log10(the.rpkm+1), day = mr.TE@meta.data$embryo.day)

ggplot(vio.df, aes(x=day, y=exp, fill = day)) + 
  geom_violin(trim=F)+ geom_boxplot(width=0.2, fill = "white") + ggtitle(the.gene.vioplot)


mr.RPKM[the.gene.vioplot, match(rownames(mr.TE@meta.data), colnames(mr.RPKM))  ]



mr.RPKM.average.per.day = data.frame(d6 = rowMeans(mr.RPKM[,match(rownames(mr.TE@meta.data)[mr.TE@meta.data$embryo.day == 6], colnames(mr.RPKM))  ]),
                                     d7 = rowMeans(mr.RPKM[,match(rownames(mr.TE@meta.data)[mr.TE@meta.data$embryo.day == 7], colnames(mr.RPKM))  ]),
                                     d8 = rowMeans(mr.RPKM[,match(rownames(mr.TE@meta.data)[mr.TE@meta.data$embryo.day == 8], colnames(mr.RPKM))  ]),
                                     d9 = rowMeans(mr.RPKM[,match(rownames(mr.TE@meta.data)[mr.TE@meta.data$embryo.day == 9], colnames(mr.RPKM))  ]),
                                     d10 = rowMeans(mr.RPKM[,match(rownames(mr.TE@meta.data)[mr.TE@meta.data$embryo.day == 10], colnames(mr.RPKM))  ])
                                     )
  

which(apply(mr.RPKM.average.per.day[TE.cell.markers.100,],1,min) > 100)
# maintained TE markers





####################################

FeaturePlot(mr.Embryo.300markerPCA, features.plot = c("PTGES", "EMP2", "TGFBR3",  "PDGFA"),reduction.use = "pca", 
            nCol = 4, 
            cols.use = c("grey","blue"), pt.size = 2)


VlnPlot(mr.Embryo.300markerPCA,features.plot = c("PDGFRA","HNF1B","BMP2","DPPA5","POU5F1","TDGF1"), group.by="Phase")

lineage.marker.df = (as.data.frame(t(as.matrix(mr.Embryo.300markerPCA@data[c("PDGFRA","HNF1B","BMP2","DPPA5","POU5F1","TDGF1",
                                                                             "GATA2","GATA3","DAB2"),]))))

#lineage.marker.df = (as.data.frame(t(as.matrix(mr.embryo.dt.rpm[c("PDGFRA","HNF1B","BMP2","DPPA5","POU5F1","TDGF1"),]))))


lineage.marker.df$phase = as.character(mr.Embryo.300markerPCA@meta.data$Phase)
lineage.marker.df$phase[lineage.marker.df$phase == "Untypical"] = "TE"

ggplot(lineage.marker.df, aes(x=phase, y=DPPA5, fill = phase)) + 
  geom_boxplot()

library(ggpubr)

gen.boxplot <- function(gene.name){
  ggplot(lineage.marker.df, aes_string(x="phase", y=gene.name, fill = "phase")) + 
    geom_boxplot() + theme_bw()+ guides(fill=FALSE) + theme(text = element_text(size=20))+ ggtitle(paste(gene.name)) + 
    ylab("log10(CPM+1)") + xlab("")
}

gen.boxplot("HNF1B")

xp = ggarrange(plotlist = lapply(c("PDGFRA","HNF1B","BMP2","DPPA5","POU5F1","TDGF1"), FUN = gen.boxplot),nrow = 2, ncol = 3)


ggarrange(plotlist = lapply(c("GATA3","GATA2"), FUN = gen.boxplot),nrow = 2, ncol = 3)

##############################################

PCAPlot(mr.TE, group.by="embryo.day")

length(mr.TE@var.genes)

?clValid::sota

mr.TE.embryo.day = SetAllIdent(mr.TE, id = "embryo.day")


mr.TE.embryo.day

mr.TE.embryo.day.marker.6.7 = FindMarkers(mr.TE.embryo.day, logfc.threshold = log(1.5), ident.1 = "6", ident.2 = "7")
dim(mr.TE.embryo.day.marker.6.7)

mr.TE.embryo.day.marker.7.8 = FindMarkers(mr.TE.embryo.day, logfc.threshold = log(1.5), ident.1 = "7", ident.2 = "8")
dim(mr.TE.embryo.day.marker.7.8)

mr.TE.embryo.day.marker.8.9 = FindMarkers(mr.TE.embryo.day, logfc.threshold = log(1.5), ident.1 = "8", ident.2 = "9")
dim(mr.TE.embryo.day.marker.8.9)

mr.TE.embryo.day.marker.9.10 = FindMarkers(mr.TE.embryo.day, logfc.threshold = log(1.5), ident.1 = "9", ident.2 = "10")
dim(mr.TE.embryo.day.marker.9.10)


mr.TE.embryo.day.average.expression = AverageExpression(mr.TE.embryo.day)


mr.TE.embryo.day.average.expression.pair.gene = mr.TE.embryo.day.average.expression[
  unique(c( rownames(mr.TE.embryo.day.marker.6.7),  rownames(mr.TE.embryo.day.marker.7.8),
                                                           rownames(mr.TE.embryo.day.marker.8.9), rownames(mr.TE.embryo.day.marker.9.10))),]

mr.TE.embryo.day.average.expression.pair.gene.log = log10(mr.TE.embryo.day.average.expression.pair.gene+1)

cp = clValid::sota(mr.TE.embryo.day.average.expression.pair.gene.log, 4,distance = "correlation")

#### use WGCNA
#######################################
library(WGCNA)
library("flashClust")
options(stringsAsFactors = FALSE);
#enableWGCNAThreads()
disableWGCNAThreads()



# apply(mr.RPKM, MARGIN = 1, function(x){sum(x>1)>20}  )

expressed.genes.mr.TE = rownames(mr.RPKM)[(apply(mr.RPKM[,colnames(mr.TE.embryo.day@data)], MARGIN = 1, function(x){sum(x>1)>30}  ))]


mydata.day=log10(mr.TE.embryo.day@data+1)

dim(mydata.day)

datExpr.day = t(as.matrix(mydata.day))
datExpr.day.used = datExpr.day[,unique(c( rownames(mr.TE.embryo.day.marker.6.7),  rownames(mr.TE.embryo.day.marker.7.8),
                                          rownames(mr.TE.embryo.day.marker.8.9), rownames(mr.TE.embryo.day.marker.9.10)))]

all.g = unique(c((c( rownames(mr.TE.embryo.day.marker.6.7),  rownames(mr.TE.embryo.day.marker.7.8),
                    rownames(mr.TE.embryo.day.marker.8.9), rownames(mr.TE.embryo.day.marker.9.10))),
          mr.TE.embryo.day@var.genes))

datExpr.day.used = datExpr.day[,mr.TE.embryo.day@var.genes]

datExpr.day.used = datExpr.day[,all.g] # used this one

#datExpr.day.used = datExpr.day[,expressed.genes.mr.TE] 




#######

powers = c(1:20)
sft=pickSoftThreshold(datExpr.day.used,dataIsExpr = TRUE,
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

softPower = 8; # used this one
#softPower = 6;

adj= adjacency(datExpr.day.used,type = "signed", power = softPower);

#Turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr.day.used,
                          networkType = "signed", 
                          TOMType = "signed", power = softPower);

colnames(TOM) =rownames(TOM) = colnames(datExpr.day.used)
dissTOM = 1-TOM

######

#Hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");

#Plot the resulting clustering tree (dendrogram)
par(mfrow = c(1,1))
plot(geneTree, xlab="", sub="",cex=0.3);

##########


minModuleSize = 20;

# Module identification using dynamic tree cut, you can also choose the hybrid method

dynamicMods = cutreeDynamic(dendro = geneTree,  
                            method="tree", 
                            minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, 
#                            distM = dissTOM, method="hybrid", deepSplit = 3, 
#                            pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#Get the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

## dynamicModsƒnotch
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

#########

MEList = moduleEigengenes(datExpr.day.used, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "", hang=-1)
MEDissThres = 0.15
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
######################
# merge tree

merge = mergeCloseModules(datExpr.day.used, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


######################



###################

embryo.day.group.dummy = model.matrix(~  as.factor(mr.TE@meta.data$embryo.day)  -1)

trait.exp.cor.mat = cor(embryo.day.group.dummy,datExpr.day.used,use = "p")

library(fifer)
trait.exp.cor.mat.color.1 = number.to.colors(trait.exp.cor.mat[1,], colors = c("green", "red"), num = 100)
trait.exp.cor.mat.color.2 = number.to.colors(trait.exp.cor.mat[2,], colors = c("green", "red"), num = 100)
trait.exp.cor.mat.color.3 = number.to.colors(trait.exp.cor.mat[3,], colors = c("green", "red"), num = 100)
trait.exp.cor.mat.color.4 = number.to.colors(trait.exp.cor.mat[4,], colors = c("green", "red"), num = 100)
trait.exp.cor.mat.color.5 = number.to.colors(trait.exp.cor.mat[5,], colors = c("green", "red"), num = 100)
#trait.exp.cor.mat.color.6 = number.to.colors(trait.exp.cor.mat[6,], colors = c("green", "red"), num = 100)


plotDendroAndColors(geneTree, cbind(dynamicColors,
                                    trait.exp.cor.mat.color.1,
                                    trait.exp.cor.mat.color.2,
                                    trait.exp.cor.mat.color.3,
                                    trait.exp.cor.mat.color.4,
                                    trait.exp.cor.mat.color.5
                                    
),
c("Gene Module", "Day 6","Day 7","Day 8","Day 9","Day 10"

),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)


###########################################
# output modules

WGCNA.modules.full.gene.list = data.frame(module = dynamicColors,
  gene = colnames(datExpr.day.used) )

write.csv(WGCNA.modules.full.gene.list, file = "WGCNA.modules.full.gene.list.csv",quote=F)
##########
#intersect
##################
nGenes = ncol(datExpr.day.used);
nSamples = nrow(datExpr.day.used);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr.day.used, dynamicColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, embryo.day.group.dummy, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = 6:10,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#################################

######    output hub genes 
chooseTopHubInEachModule( datExpr.day.used, dynamicColors, power = softPower, type= "signed")

modules.to.plot = setdiff( unique(dynamicColors) ,"grey")

hub.gene.list = c()
hub.gene.assignment = c()
hub.gene.connectivity = c()

for(hub.this.module in modules.to.plot){
  
  nTopHubs = 30 
  
  hub.this.kIN = softConnectivity(datExpr.day.used[, dynamicColors == hub.this.module ])
  
  selectHubs = (rank (-hub.this.kIN) <= nTopHubs) 
  
  
  
  hub.gene.list = c(hub.gene.list, colnames(datExpr.day.used)[dynamicColors == hub.this.module][selectHubs])
  hub.gene.assignment = c(hub.gene.assignment, rep(hub.this.module,nTopHubs ))
  hub.gene.connectivity = c(hub.gene.connectivity,  hub.this.kIN[selectHubs]  )
  
  
}


hub.data.frame = data.frame(hub.gene.list, hub.gene.assignment, hub.gene.connectivity)

write.csv(hub.data.frame, file = "hub.data.frame__WGCNA hub.gene.list.csv",quote=F)

#########################
# 这些是画line plot
modules.to.plot = setdiff(unique(dynamicColors),"grey")
total.module.matrix.to.heatmap = NULL
total.module.matrix.to.line.plot = NULL
side.colour = c()
  
for(the.module in modules.to.plot){

  scaled.data.to.use = mr.TE@scale.data[colnames(datExpr.day.used)[dynamicColors == the.module],]
  
  this.module.mean.by.day = NULL
  
  for(l in  levels(mr.TE@meta.data$embryo.day) ){
    if(is.null(this.module.mean.by.day)){
      this.module.mean.by.day = rowMeans(scaled.data.to.use[,mr.TE@meta.data$embryo.day == l])
    }else{
      this.module.mean.by.day = cbind(this.module.mean.by.day, rowMeans(scaled.data.to.use[,mr.TE@meta.data$embryo.day == l]))
    }
    
    
  }
  # add heatmap.plot
  if(is.null(total.module.matrix.to.heatmap)){
    total.module.matrix.to.heatmap = this.module.mean.by.day
  }else{
    total.module.matrix.to.heatmap = rbind(total.module.matrix.to.heatmap,this.module.mean.by.day)
  }
  
  # add line.plot
  if(is.null(total.module.matrix.to.line.plot)){
    total.module.matrix.to.line.plot = colMeans(this.module.mean.by.day)
  }else{
    total.module.matrix.to.line.plot = rbind(total.module.matrix.to.line.plot,colMeans(this.module.mean.by.day))
  }
  
  
  side.colour = c(side.colour, rep(the.module, dim(scaled.data.to.use)[1]))
  
}
  
  
  
library(gplots)
palette.breaks <- seq(-1.5, 1.5, 0.1)
color.palette = colorRampPalette(c("dodgerblue4","dodgerblue1","white","firebrick1","firebrick3"), 
                                 space="Lab")(length(palette.breaks) - 1)

#color.palette = colorRampPalette(c("green","black","red"), 
#                                 space="Lab")(length(palette.breaks) - 1)

pdf(file = "add.pdf")
heatmap.2(as.matrix( total.module.matrix.to.heatmap ),trace = "none",density = "none",
          Colv = NA,
          Rowv  = NA,
          #ColSideColors = c("grey","red","yellow","green")[ as.numeric(Xy.cc.noD26@ident) ] ,
          RowSideColors = side.colour,
          col=color.palette,
          #breaks = palette.breaks,
          scale = c("row")
          
)
dev.off()

plot(total.module.matrix.to.line.plot[1,],type="l")
plot(total.module.matrix.to.line.plot[2,],type="l")
plot(total.module.matrix.to.line.plot[3,],type="l")
plot(total.module.matrix.to.line.plot[4,],type="l")

matplot(t(total.module.matrix.to.line.plot),type="l")


make.day.lineplot <- function(series, name){
  
  df = data.frame(exp = as.numeric(total.module.matrix.to.line.plot[series,]),
                  day = (6:10)  )
  ggplot(data=df, aes(x=day, y=exp, group=1)) +
    geom_line(color = "red")+
    geom_point(color = "red") + ggtitle(name) + ylab("") + xlab("")#+theme_classic()
}


make.day.lineplot(8,'a')

lineplot.list = list(make.day.lineplot(1,'yellow'),
                     make.day.lineplot(2,'blue'),
                     make.day.lineplot(3,'green'),
                     make.day.lineplot(4,'turquoise'),
                     
                     make.day.lineplot(5,'red'),
                     make.day.lineplot(6,'brown'),
                     make.day.lineplot(7,'pink'),
                     make.day.lineplot(8,'black')
                     )

library(ggpubr)

ggarrange(plotlist = lineplot.list, ncol = 1, nrow = 8)

############
# 第一种方法
modules.choosed = c("turquoise","black","red","pink")

genes.choosed = (colnames(datExpr.day.used))[dynamicColors %in% modules.choosed]

mr.RPKM.average.per.day.gene.choosed = mr.RPKM.average.per.day[genes.choosed,]

mr.RPKM.average.per.day.gene.choosed = mr.RPKM.average.per.day.gene.choosed[osj,]

mr.RPKM.average.per.day.gene.choosed.reshaped = data.frame(exp=c(mr.RPKM.average.per.day.gene.choosed[,1],
                                                                 mr.RPKM.average.per.day.gene.choosed[,2],
                                                                 mr.RPKM.average.per.day.gene.choosed[,3],
                                                                 mr.RPKM.average.per.day.gene.choosed[,4],
                                                                 mr.RPKM.average.per.day.gene.choosed[,5]) ,
                                                           day = rep(c(6:10), each = dim(mr.RPKM.average.per.day.gene.choosed)[1])
                                                           )


library(ggplot2)


ggplot(mr.RPKM.average.per.day.gene.choosed.reshaped, aes(x= factor(day) , y=  log10(exp+1)  )) + 
  geom_boxplot() + theme_bw()

######################
# 第二种方法

#modules.choosed = c("turquoise","black","red","pink")
#modules.choosed = c("yellow","blue")
modules.choosed = c("brown","green")

genes.choosed = (colnames(datExpr.day.used))[dynamicColors %in% modules.choosed]

mr.RPKM.clean.module.average = data.frame(module.exp = colMeans(mr.RPKM.clean[genes.choosed,  ],na.rm = T) )

mr.RPKM.clean.module.average.day.df = data.frame(module.exp = mr.RPKM.clean.module.average$module.exp, day = mr.TE@meta.data$embryo.day)

library(ggplot2)


ggplot(mr.RPKM.clean.module.average.day.df, aes(x= (day) , y=  (module.exp) , fill = day )) + 
  geom_boxplot(outlier.shape = NA) + theme_bw()+ guides(fill=FALSE) + ylab("Ave FPKM")  + xlab("") +theme(text = element_text(size=20)) + ylab("")


####################################
# 输出module list

gene.nodle.df = data.frame( gene = (colnames(datExpr.day.used)), module = dynamicColors   )

write.csv(gene.nodle.df,file = "gene.nodle.df.csv", quote=F)



##################################
# WGCNA module 的Go 分析图

library(clusterProfiler)

library(org.Hs.eg.db)

GO.plot <- function(gexntez,title){
  
  hs <- org.Hs.eg.db
  my.symbols <- gexntez
  ENTREZID.id = select(hs, 
                       keys = my.symbols,
                       columns = c("ENTREZID", "SYMBOL"),
                       keytype = "SYMBOL")$ENTREZID
  
  ego <- enrichGO(gene          = ENTREZID.id,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  dotplot(ego)+ggtitle(title)
}


#modules.choosed = c("turquoise","black","red","pink")
#modules.choosed = c("yellow","blue")
# modules.choosed = c("brown","green")


GO.plot((colnames(datExpr.day.used))[dynamicColors %in% c("turquoise","black","red","pink")]  , "Late")

GO.plot((colnames(datExpr.day.used))[dynamicColors %in% c("yellow","blue")]  , "Early")

GO.plot((colnames(datExpr.day.used))[dynamicColors %in% c("brown","green")]  , "Middle")

###################################
# 找hub gene，输出hub gene network


################################################

# figure 9 10 classification

PCAPlot(mr.TE, group.by = "embryo.day")


###############################
# day 9, 10 , cell type classification

mr.TE.fullPCA.day9.10

complete.cluster.day9.10 = stats::hclust(as.dist(1-abs(cor((   as.matrix(mr.TE.fullPCA.day9.10@data[mr.TE.subtypeA.day9.10@var.genes,])     ),method="spearman"))), 
                          method="ward.D")

#complete.cluster.day9.10 = stats::hclust(as.dist(1-abs(cor((   as.matrix(mr.TE.fullPCA.day9.10@data)   ),method="spearman"))), 
#                                         method="ward.D")


plot(complete.cluster.day9.10, hang = -1,labels=FALSE)

dddt = data.frame(col = dendextend::cutree(complete.cluster.day9.10,k=3,order_clusters_as_data = FALSE))

dendextend::colored_bars(colors = dddt, dend = complete.cluster.day9.10)


hclust.ident = cutree(complete.cluster.day9.10,k=3)

hclust.ident = data.frame(hclust.ident)
rownames( hclust.ident) = names(cutree(complete.cluster.day9.10,k=3))
mr.TE.fullPCA.day9.10 = AddMetaData(mr.TE.fullPCA.day9.10, metadata = hclust.ident, col.name = "hclust.ident")

#########################
# day 9 day 10 hclust with colour bar

par( oma = c(0,1,1,1), mgp = c(1,0.5,0), mar = c(10,2,2,2) )
library( squash )
library( dendextend )


hclust.ident = cutree(complete.cluster.day9.10,k=3)

myColors <- matrix(NA,nrow = length(hclust.ident), ncol=2)
myColors[,1] = gg_color_hue(3)[hclust.ident]

myColors[,2] =  c("#F0E442","#F0E442")[ as.numeric(as.character(mr.TE.fullPCA.day9.10@meta.data$embryo.day))-8]


complete.cluster.day9.10.2 = complete.cluster.day9.10
complete.cluster.day9.10.2$labels = ""

myDend <-  as.dendrogram(complete.cluster.day9.10.2) 
plot(myDend)
colored_bars(colors = myColors, dend = myDend, rowLabels = c("Cluster","Day"),y_scale = 7,cex.rowLabels=1)


#########################################

# distribution of pair-wise correlation


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

boxplot(cor(  log10(mr.RPKM.clean[,mr.TE@meta.data$embryo.day == 6]+1)  ,method = "spearman" )[upper.tri(cor(mr.RPKM.clean[,mr.TE@meta.data$embryo.day == 6]))],
        cor(  log10(mr.RPKM.clean[,mr.TE@meta.data$embryo.day == 7]+1)  ,method = "spearman" )[upper.tri(cor(mr.RPKM.clean[,mr.TE@meta.data$embryo.day == 7]))],
        cor(  log10(mr.RPKM.clean[,mr.TE@meta.data$embryo.day == 8] +1) ,method = "spearman" )[upper.tri(cor(mr.RPKM.clean[,mr.TE@meta.data$embryo.day == 8]))],
        cor(  log10(mr.RPKM.clean[,mr.TE@meta.data$embryo.day == 9]+1)  ,method = "spearman" )[upper.tri(cor(mr.RPKM.clean[,mr.TE@meta.data$embryo.day == 9]))],
        cor(  log10(mr.RPKM.clean[,mr.TE@meta.data$embryo.day == 10] +1) ,method = "spearman" )[upper.tri(cor(mr.RPKM.clean[,mr.TE@meta.data$embryo.day == 10]))],
        
        outline=FALSE, names = c(6:10), col = gg_color_hue(5))

PCAPlot(mr.TE.fullPCA.day9.10)



mr.TE.subtypeA.day9.10 <- NormalizeData(object = mr.TE.fullPCA.day9.10, normalization.method = "LogNormalize", 
                    scale.factor = 10000)


#########


mr.TE.subtypeA.day9.10 <- FindVariableGenes(object = mr.TE.subtypeA.day9.10, mean.function = ExpMean, dispersion.function = LogVMR,do.plot = F)

length(x = mr.TE.subtypeA.day9.10@var.genes)

mr.TE.subtypeA.day9.10 <- ScaleData(object = mr.TE.subtypeA.day9.10, vars.to.regress = c("nUMI"))


mr.TE.subtypeA.day9.10 <- RunPCA(object = mr.TE.subtypeA.day9.10, pc.genes = mr.TE.subtypeA.day9.10@var.genes, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5,pcs.compute = 40)

PCAPlot(mr.TE.subtypeA.day9.10)

FeaturePlot(mr.TE.subtypeA.day9.10, features.plot = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","HLA-G","ERVW-1","CD9","ITGA1"), reduction.use = "pca")



mr.TE.subtypeA.day9.10 = RunTSNE(mr.TE.subtypeA.day9.10, genes.use = mr.TE.subtypeA.day9.10@var.genes)

FeaturePlot(mr.TE.subtypeA.day9.10, features.plot = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","HLA-G","ERVW-1","CD9","ITGA1"), reduction.use = "tsne")

mr.TE.subtypeA.day9.10.merkers = FindAllMarkers(mr.TE.subtypeA.day9.10, only.pos = T )

mr.TE.subtypeA.day9.10.merkers.sig = mr.TE.subtypeA.day9.10.merkers[mr.TE.subtypeA.day9.10.merkers$p_val_adj < 0.05 & 
                                                                      mr.TE.subtypeA.day9.10.merkers$avg_logFC > 0.5,]



mr.TE.subtypeA.day9.10 = AddMetaData(mr.TE.subtypeA.day9.10, metadata = hclust.ident, col.name = "hclust.ident")



TSNEPlot(mr.TE, group.by = "embryo.day")


PCAPlot(mr.TE.subtypeA.day9.10, group.by = "hclust.ident")


PCAPlot(mr.TE.fullPCA.day9.10, group.by = "hclust.ident")

PCAPlot(mr.TE.subtypeA.day9.10, group.by = "hclust.ident")


mr.TE.subtypeA.day9.10.hclust.type = SetAllIdent(mr.TE.subtypeA.day9.10, id = "hclust.ident")

mr.TE.subtypeA.day9.10.hclust.type.markers = FindAllMarkers(mr.TE.subtypeA.day9.10.hclust.type, only.pos = T )


mr.TE.subtypeA.day9.10.hclust.type.markers.sig = mr.TE.subtypeA.day9.10.hclust.type.markers[mr.TE.subtypeA.day9.10.hclust.type.markers$p_val_adj < 0.05 & 
                                                                      mr.TE.subtypeA.day9.10.hclust.type.markers$avg_logFC >0.5,]


DoHeatmap(mr.TE.subtypeA.day9.10.hclust.type, genes.use = (mr.TE.subtypeA.day9.10.hclust.type.markers.sig$gene) )


mr.TE.subtypeA.day9.10.hclust.type.markers.sig[c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                               "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                               
                                               "PSG1","CSH1","HSD3B1","CYP19A1","SDC1","INHA",
                                               "ERVW-1","ERVV-1","CGA","CGB5"),]


mr.TE.subtypeA.day9.10.hclust.type = 
  CellCycleScoring(mr.TE.subtypeA.day9.10.hclust.type, g2m.genes = cc.genes$g2m.genes, s.genes = cc.genes$g2m.genes, set.ident = T)

PCAPlot(mr.TE.subtypeA.day9.10.hclust.type, group.by = "Phase")

mr.TE.subtypeA.day9.10.hclust.type <- ScaleData(object = mr.TE.subtypeA.day9.10.hclust.type, vars.to.regress = c("nUMI","S.Score","G2M.Score"))

mr.TE.subtypeA.day9.10.hclust.type <- RunPCA(object = mr.TE.subtypeA.day9.10.hclust.type, 
                                             pc.genes = mr.TE.subtypeA.day9.10.hclust.type@var.genes, do.print = TRUE, pcs.print = 1:5, 
                                 genes.print = 5,pcs.compute = 40)

PCAPlot(mr.TE.subtypeA.day9.10.hclust.type)
## cell cycle correction 不管用

## public data

library()
ITG6, TP63，CTNNB1,LRP5, TP63，TEAD4,ELF5,FGFR2,FZD5
EVT marker: ITG5, HLA-G,FN1,MMP3,CD9,ITGA1

VlnPlot(mr.TE.subtypeA.day9.10.hclust.type, features.plot = c("ERVW-1","CGA","CGB5", "ITGA6","TP63","LRP5", "ITGA5","HLA-G","ITGA1"),nCol = 3)

library(readxl)


mr.TE.subtypeA.day9.10.hclust.type.markers.sig[mr.TE.subtypeA.day9.10.hclust.type.markers.sig$cluster == 1,]
mr.TE.subtypeA.day9.10.hclust.type.markers.sig[mr.TE.subtypeA.day9.10.hclust.type.markers.sig$cluster == 2,]
mr.TE.subtypeA.day9.10.hclust.type.markers.sig[mr.TE.subtypeA.day9.10.hclust.type.markers.sig$cluster == 3,]

mr.TE.subtypeA.day9.10.hclust.type.markers.sig$gene[mr.TE.subtypeA.day9.10.hclust.type.markers.sig$cluster == 3]

c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5") %in% mr.TE.subtypeA.day9.10.hclust.type.markers.sig$gene  # cluster 3

c("FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1") %in% mr.TE.subtypeA.day9.10.hclust.type.markers.sig$gene # cluster 1 and 3


FeaturePlot(mr.TE.subtypeA.day9.10, features.plot = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                                      "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                                      
                                                      "ERVW-1","ERVV-1","CGA","CGB5"), reduction.use = "tsne")

DoHeatmap(mr.TE.subtypeA.day9.10.hclust.type, genes.use = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                                "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                                
                                                "PSG1","CSH1","HSD3B1","CYP19A1","SDC1","INHA",
                                                "ERVW-1","ERVV-1","CGA","CGB5"))

DoHeatmap(mr.TE, 
          disp.min = -2, disp.max = 2.5,genes.use = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                                            "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                                            
                                                            "PSG1","CSH1","HSD3B1","CYP19A1","SDC1","INHA",
                                                            "ERVW-1","ERVV-1","CGA","CGB5"))


mr.TE.subtypeA.day9.10.hclust.type.markers.sig[mr.TE.subtypeA.day9.10.hclust.type.markers.sig$gene == "ITGA5",]



##################


mr.TE.day8 <- NormalizeData(object = mr.TE.day8, normalization.method = "LogNormalize", 
                                        scale.factor = 10000)


mr.TE.day8 <- FindVariableGenes(object = mr.TE.day8, mean.function = ExpMean, dispersion.function = LogVMR,do.plot = F)

length(x = mr.TE.day8@var.genes)

mr.TE.day8 <- ScaleData(object = mr.TE.day8, vars.to.regress = c("nUMI"))



mr.TE.day8 = SetAllIdent(mr.TE.day8, "old.ident")

DoHeatmap(mr.TE.day8, genes.use = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                    "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                    "PSG1","CSH1","HSD3B1","CYP19A1","SDC1","INHA",
                                    "ERVW-1","ERVV-1","CGA","CGB5"))



complete.cluster.day8 = stats::hclust(as.dist(1-abs(cor((   as.matrix(mr.TE.day8@data[mr.TE.day8@var.genes,])     ),method="spearman"))), 
                                         method="ward.D")

complete.cluster.day8.2 = complete.cluster.day8
complete.cluster.day8.2$labels = ""


hclust.ident.day8 = cutree(complete.cluster.day8,k=2)
hclust.ident.day8 = data.frame(hclust.ident.day8)
rownames( hclust.ident.day8) = names(cutree(complete.cluster.day8,k=2))
mr.TE.day8 = AddMetaData(mr.TE.day8, metadata = hclust.ident.day8, col.name = "hclust.ident.day8")
mr.TE.day8 = SetAllIdent(mr.TE.day8, "hclust.ident.day8")

DoHeatmap(mr.TE.day8, genes.use = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                    "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                    
                                    "PSG1","CSH1","HSD3B1","CYP19A1","SDC1","INHA",
                                    "ERVW-1","ERVV-1","CGA","CGB5") )




# day 8 hclust, colour bar

par( oma = c(0,1,1,1), mgp = c(1,0.5,0), mar = c(10,2,2,2) )
library( squash )
library( dendextend )

hclust.ident.day8 = cutree(complete.cluster.day8.2,k=2)

myColors.day8 <- matrix(NA,nrow = length(hclust.ident.day8), ncol=2)

myColors.day8[,1] = gg_color_hue(2)[hclust.ident.day8]

myColors.day8[,2] =  c("red","green","blue")[ as.numeric(as.factor(substr(colnames(mr.TE.fullPCA.day8@data),1,3)))]


myDend <-  as.dendrogram(complete.cluster.day8.2) 
plot(myDend)
colored_bars(colors = myColors.day8, dend = myDend, rowLabels = c("Cluster","Embryo"),y_scale = 2,cex.rowLabels=1)



DoHeatmap(mr.TE.day8, genes.use = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                    "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                    
                                    "PSG1","CSH1","HSD3B1","CYP19A1","SDC1","INHA",
                                    "ERVW-1","ERVV-1","CGA","CGB5"), cells.use = colnames(mr.TE.day8@data)[complete.cluster.day8.2$order])



palette.breaks <- seq(-2.5, 2.5, 0.1)
color.palette = colorRampPalette(c("purple","black","yellow"), 
                                 space="Lab")(length(palette.breaks) - 1)

gplots::heatmap.2(as.matrix(  mr.TE.day8@scale.data[match(c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                                            "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                                            "PSG1","HSD3B1","CYP19A1","SDC1","INHA",
                                                            "ERVW-1","ERVV-1","CGA","CGB5"), rownames(mr.TE.day8@scale.data))
                                                    ,  complete.cluster.day8.2$order ] ),
          trace="none",density="none",
          Colv = F,
          Rowv = F,
          key=F,scale="row",
          dendrogram="none",
          breaks = palette.breaks,
          col = colorRampPalette(c("#FF00FF","#000000","#FFFF00"), space="Lab"),
          margins=c(10.5,6),
          ColSideColors = myColors.day8[,1][complete.cluster.day8.2$order],
          labRow=NA,
          main="Tissue specific genes"
          
)









################################


plot(mr.TE.day8@dr$pca@gene.loadings[,1], mr.TE.day8@dr$pca@gene.loadings[,2], 
     #col= alpha(col.vector, alpha.vector) ,
     pch=16, ylab = "PC2 loading", xlab = "PC1 loading" ,
     cex.main=1.25, cex.lab=1.5, cex.axis=1.25, main = "Gene loading")

VizPCA(mr.TE.day8,pcs.use = 1, num.genes = 50,do.balanced=F)




complete.cluster.day8 = stats::hclust(as.dist(1-abs(cor((   mr.TE.day8@data[names(c(head(sort(mr.TE.day8@dr$pca@gene.loadings[,1]),500), tail(sort(mr.TE.day8@dr$pca@gene.loadings[,1]),500))),]     ),method="spearman"))), 
                                      method="ward.D")


plot(complete.cluster.day8, hang = -1)


hclust.ident.day8 = cutree(complete.cluster.day8,k=2)
hclust.ident.day8 = data.frame(hclust.ident.day8)

rownames( hclust.ident.day8) = names(cutree(complete.cluster.day8,k=2))

mr.TE.day8 = AddMetaData(mr.TE.day8, metadata = hclust.ident.day8, col.name = "hclust.ident.day8")

mr.TE.day8 = SetAllIdent(mr.TE.day8, "hclust.ident.day8")

DoHeatmap(mr.TE.day8, genes.use = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                    "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                    "PSG1","CSH1","HSD3B1","CYP19A1","SDC1","INHA",
                                    "ERVW-1","ERVV-1","CGA","CGB5"))



mr.TE.day8.2.type.markers = FindMarkers(mr.TE.day8, ident.1 = 1, ident.2 = 2)

write.csv(mr.TE.day8.2.type.markers, file = "Dec_9_2018_day8.2.types.marker.csv", quote=F)

mr.TE.day8.2.type.markers.sig = mr.TE.day8.2.type.markers[mr.TE.day8.2.type.markers$p_val_adj<0.05 & abs(mr.TE.day8.2.type.markers$avg_logFC) > 0.5,]


DoHeatmap(mr.TE.day8, genes.use = c(rownames(mr.TE.day8.2.type.markers.sig[mr.TE.day8.2.type.markers.sig$avg_logFC > 0,]),
                                    rownames(mr.TE.day8.2.type.markers.sig[mr.TE.day8.2.type.markers.sig$avg_logFC < 0,])) )

##################


c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5") %in% mr.day9.10.rePCA.EVT_CT@var.genes

c("FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1") %in% mr.day9.10.rePCA.EVT_CT@var.genes

c("ERVV-1","ERVW-1","CGA","CGB5") %in% mr.day9.10.rePCA.EVT_CT@var.genes

PCAPlot(mr.day9.10.rePCA.EVT_CT)


#################


pp =FindAllMarkers(mr.TE.fullPCA.day9.10, only.pos = T)

DoHeatmap(mr.TE.fullPCA.day9.10, genes.use = (pp[pp$p_val_adj<0.05 & pp$avg_logFC > 0.5,]$gene) )

#############


mr.TE.markers = FindAllMarkers(mr.TE, only.pos = T)

pp = mr.TE.markers[mr.TE.markers$p_val_adj<0.05 & mr.TE.markers$avg_logFC> 0.5,]

DoHeatmap(mr.TE, genes.use = (pp$gene) )

pp[pp$cluster == 4,]

########################


csc.pub.data = read.csv("1-s2.0-S1934590917304563-mmc2.csv",header = T)

col.used = c(2:4)

get.marker = function(col.used, p.val = 0.05){
  
  col.used = c(5:7)
  p.val=0.05
  
  this.data = csc.pub.data[,col.used]
  compared.to.data = csc.pub.data[, setdiff(2:10, col.used) ]
  
  flag1 = rowMeans(this.data) > log2(10+1)
  
  flag2 = rowMeans(this.data) - rowMeans(compared.to.data) > 1
  
  t.p.value = NULL
  
  for(i in 1:24103){
    
    if(var(as.numeric(this.data[i,])) == 0 | var(as.numeric(compared.to.data[i,])) == 0 ){
      t.p.value[i] = 1
      next
    }
    t.p.value[i] = t.test(as.numeric(this.data[i,]),  as.numeric(compared.to.data[i,]))$p.value
    
  }
  
  flag3 = (t.p.value) < 0.01
  
  sum(flag1 & flag2 & flag3)
  
  csc.pub.data$Gene_symbol[flag1 & flag2 & flag3]
  
  
}

csc.CT.markers = get.marker(2:4)
csc.EVT.markers = get.marker(col.used= 5:7,p.val = 0.1)
csc.ST.markers = get.marker(8:10)


mr.TE.subtypeA.day9.10.hclust.type.marker.score  = CellCycleScoring(mr.TE.subtypeA.day9.10.hclust.type,
                                                            g2m.genes = csc.CT.markers,
                                                          s.genes = m)


GenePlot(mr.TE.subtypeA.day9.10.hclust.type.marker.score, gene1 = "S.Score",
         gene2 = "G2M.Score")


DoHeatmap(mr.TE.subtypeA.day9.10.hclust.type, genes.use = c(csc.CT.markers,m, csc.ST.markers))

########################
get.marker.2 = function(col.used,col.used2, p.val = 0.05){
  
  col.used = c(2:4)
  col.used2 = c(5:7)
  
  this.data = csc.pub.data[,col.used]
  compared.to.data = csc.pub.data[, col.used2 ]
  
  flag1 = rowMeans(this.data) > log2(10+1)
  
  flag2 = rowMeans(this.data) - rowMeans(compared.to.data) > 1
  
  t.p.value = NULL
  
  for(i in 1:24103){
    
    if(var(as.numeric(this.data[i,])) == 0 | var(as.numeric(compared.to.data[i,])) == 0 ){
      t.p.value[i] = 1
      next
    }
    t.p.value[i] = t.test(as.numeric(this.data[i,]),  as.numeric(compared.to.data[i,]))$p.value
    
  }
  
  flag3 = p.adjust(t.p.value) < 0.05
  
  sum(flag1 & flag2 & flag3)
  
  csc.pub.data$Gene_symbol[flag1 & flag2 & flag3]
  
  
}

m = csc.pub.data$Gene_symbol[flag1 & flag2 & flag3]

####################

hdsf = FindMarkers(mr.TE.subtypeA.day9.10.hclust.type,ident.1 = 1,ident.2 = 3 )

PCAPlot(mr.TE.subtypeA.day9.10.hclust.type)

#FeaturePlot(mr.TE.subtypeA.day9.10.hclust.type, features.plot = c("ERVW-1","CGA"), reduction.use = "pca")

DoHeatmap(mr.TE.subtypeA.day9.10.hclust.type, genes.use = c(rownames(hdsf[hdsf$p_val_adj<0.05 & hdsf$avg_logFC> 1,]),
                                                            rownames(hdsf[hdsf$p_val_adj<0.05 & hdsf$avg_logFC< -1,])))


hdsf["DPPA3",]

ds = FindAllMarkers(mr.TE.subtypeA.day9.10.hclust.type,only.pos = T)



DoHeatmap(mr.TE.subtypeA.day9.10.hclust.type, genes.use = ds$gene[ds$avg_logFC> 0.5 & ds$p_val_adj<0.05])

write.csv(ds, file = "Dec_9_2018_day9_10.3.type.markers.csv", quote=F)

c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5") %in%  ds$gene[ds$avg_logFC> 0.5 & ds$p_val_adj<0.05]

c("CD9","ITGA1","HLA-G","FN1","MMP2","ITGA5") %in%  ds$gene[ds$avg_logFC> 0.5 & ds$p_val_adj<0.05]

ds$gene[ds$avg_logFC> 0.5 & ds$p_val_adj<0.05 & ds$cluster==1]

VlnPlot(mr.TE.subtypeA.day9.10.hclust.type, features.plot = c("TET2","CD9","ITGA1","HLA-G","FN1","MMP2","ITGA5","ITGA6","TP63","CGB5"))

VlnPlot(mr.TE.subtypeA.day9.10.hclust.type, features.plot = c("CD9","ITGA1","HLA-G","FN1","MMP2","ITGA5"))
VlnPlot(mr.TE.subtypeA.day9.10.hclust.type, features.plot = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5"))

mr.TE.subtypeA.day9.10.hclust.type.CT.EVT


mr.TE.subtypeA.day9.10.cluster.type = FindClusters(mr.TE.subtypeA.day9.10.hclust.type)

PCAPlot(mr.TE.subtypeA.day9.10.hclust.type)
PCAPlot(mr.TE.subtypeA.day9.10.cluster.type)



PCHeatmap(object = mr.TE.subtypeA.day9.10.hclust.type, pc.use = 1:4, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)


##############################################################
# PCA curve

mr.TE.fullPCA.noday9.10

mr.TE.fullPCA.noday9.10 <- NormalizeData(object = mr.TE.fullPCA.noday9.10, normalization.method = "LogNormalize", 
                                        scale.factor = 10000)



mr.TE.fullPCA.noday9.10 <- FindVariableGenes(object = mr.TE.fullPCA.noday9.10, mean.function = ExpMean, dispersion.function = LogVMR,do.plot = F)

length(x = mr.TE.fullPCA.noday9.10@var.genes)

mr.TE.fullPCA.noday9.10 <- ScaleData(object = mr.TE.fullPCA.noday9.10, vars.to.regress = c("nUMI"))


mr.TE.fullPCA.noday9.10 <- RunPCA(object = mr.TE.fullPCA.noday9.10, pc.genes = rownames(mr.TE.fullPCA.noday9.10@data), do.print = TRUE, pcs.print = 1:5, 
                                 genes.print = 5,pcs.compute = 40)

PCAPlot(mr.TE.fullPCA.noday9.10, group.by = "culture.condition")
PCAPlot(mr.TE.fullPCA.noday9.10, group.by = "embryo.day")


mr.TE.fullPCA.noday9.10 <- JackStraw(object = mr.TE.fullPCA.noday9.10, num.replicate = 100, display.progress = T)

JackStrawPlot(object = mr.TE.fullPCA.noday9.10, PCs = 1:20)

mr.TE.fullPCA.noday9.10 = RunTSNE(mr.TE.fullPCA.noday9.10, genes.use = mr.TE.fullPCA.noday9.10@var.genes,dims.use = 1:10, dim.embed = 1)

TSNEPlot(mr.TE.fullPCA.noday9.10,group.by = "embryo.day")
TSNEPlot(mr.TE.fullPCA.noday9.10,group.by = "culture.condition")

#####

mr.TE.fullPCA.noday9.10@dr$tsne@cell.embeddings[,1]

PCAPlot(mr.TE.fullPCA.noday9.10,group.by = "embryo.day")


wilcox.test(pseudotime.dev.plot$pos[pseudotime.dev.plot$lab =="OnlyEmbryo_7" ], pseudotime.dev.plot$pos[pseudotime.dev.plot$lab =="Coculture_7" ])



#############################


VizPCA(mr.TE, pcs.use = 1:2)

PCHeatmap(object = mr.TE, pc.use = 1:2, do.balanced = F, label.columns = FALSE)


############################


day.and.condition.label = paste(as.character(mr.TE.fullPCA.noday9.10@meta.data$culture.condition),
                                as.character(mr.TE.fullPCA.noday9.10@meta.data$embryo.day), sep = "_")

pseudotime.dev.plot  = data.frame(day.and.condition.label, pos = mr.TE@dr$pca@cell.embeddings[match(rownames(mr.TE.fullPCA.noday9.10@meta.data), rownames(mr.TE@dr$pca@cell.embeddings)),2] )

pseudotime.dev.plot$lab = factor(day.and.condition.label, 
                                 levels = c("OnlyEmbryo_6","OnlyEmbryo_7","Coculture_7","Coculture_8"))

library(ggpubr)
ggplot(pseudotime.dev.plot, aes(x=  lab  , 
                                y= pos, fill = lab )) + 
  geom_violin(trim=T) + geom_boxplot(width=0.2,outlier.shape = NA, fill = "white") + theme_bw()



wilcox.test(pseudotime.dev.plot$pos[pseudotime.dev.plot$lab =="OnlyEmbryo_7" ], pseudotime.dev.plot$pos[pseudotime.dev.plot$lab =="Coculture_7" ])

###############################

day.and.condition.label.df.add = data.frame(day.and.condition.label)

rownames(day.and.condition.label.df.add) = rownames(mr.TE.fullPCA.noday9.10@meta.data)

mr.TE.fullPCA.noday9.10 = AddMetaData(mr.TE.fullPCA.noday9.10, metadata = day.and.condition.label.df.add, col.name = "day.and.condition.label")

mr.TE.fullPCA.noday9.10 = SetAllIdent(mr.TE.fullPCA.noday9.10,id = "day.and.condition.label")

######

day7.coculture.onlyembryo.marker = FindMarkers(mr.TE.fullPCA.noday9.10, ident.1 = "OnlyEmbryo_7", ident.2 = "Coculture_7")

day7.coculture.onlyembryo.marker.sig = day7.coculture.onlyembryo.marker[day7.coculture.onlyembryo.marker$p_val_adj<0.05 & 
                                                                          abs(day7.coculture.onlyembryo.marker$avg_logFC) > 0.5,]



dim(day7.coculture.onlyembryo.marker.sig)

HTF

intersect(HTF, rownames(day7.coculture.onlyembryo.marker.sig))

rownames(day7.coculture.onlyembryo.marker.sig)[order(day7.coculture.onlyembryo.marker.sig$avg_logFC)]


DoHeatmap(mr.TE, genes.use = rownames(day7.coculture.onlyembryo.marker.sig)[order(day7.coculture.onlyembryo.marker.sig$avg_logFC)] )





########################################

PCAPlot(mr.TE.fullPCA.noday9.10)


DoHeatmap(mr.TE.fullPCA.noday9.10, genes.use = 
            rownames(day7.coculture.onlyembryo.marker.sig)[order(day7.coculture.onlyembryo.marker.sig$avg_logFC)],
          group.order = c("OnlyEmbryo_6", "OnlyEmbryo_7", "Coculture_7", "Coculture_8" ))


mr.TE.day7 = SubsetData(mr.TE.fullPCA.noday9.10, cells.use = rownames(mr.TE.fullPCA.noday9.10@meta.data)[mr.TE.fullPCA.noday9.10@meta.data$embryo.day == "7"],
                          do.center = T, do.scale = T,do.clean = T)

FindVariableGenes()


VlnPlot(mr.TE.day7, features.plot = c("CDH5","CCR7","MMP2","DPPA3","SALL4","DAB2"))

intersect(  rownames(day7.coculture.onlyembryo.marker.sig)   , polar$polar)

#######################################

mr.TE@meta.data$culture.condition

mr.TE@meta.data$embryo.day

######################################


polar = readxl::read_excel("../Polar.xlsx",sheet = "Sheet1")
polar$polar

DoHeatmap(mr.TE.fullPCA.noday9.10, genes.use = polar$polar,
          group.order = c("OnlyEmbryo_6", "OnlyEmbryo_7", "Coculture_7", "Coculture_8" ))

mr.TE.fullPCA.noday9.10.group.mean = AverageExpression(mr.TE.fullPCA.noday9.10, show.progress = T , use.scale = T)

boxplot(mr.TE.fullPCA.noday9.10.group.mean[polar$polar,3],
        mr.TE.fullPCA.noday9.10.group.mean[polar$polar,4],
        mr.TE.fullPCA.noday9.10.group.mean[polar$polar,1],
        mr.TE.fullPCA.noday9.10.group.mean[polar$polar,2], ylim = c(-2,1), outline = F, col=gg_color_hue(4))


sdsd.df = data.frame(exp = c(mr.TE.fullPCA.noday9.10.group.mean[polar$polar,3], mr.TE.fullPCA.noday9.10.group.mean[polar$polar,4],
                       mr.TE.fullPCA.noday9.10.group.mean[polar$polar,1],mr.TE.fullPCA.noday9.10.group.mean[polar$polar,2]),
                     col.lab  = rep(  1:4, each = 129 )   )

ggplot(sdsd.df, aes(x= factor(col.lab), y=exp, fill = factor(col.lab))) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + ylim(c(-2,1)) + guides(fill=FALSE)

##################


gplots::heatmap.2(as.matrix(  mr.TE.day8@scale.data[match(c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                                            "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                                            "PSG1","HSD3B1","CYP19A1","SDC1","INHA",
                                                            "ERVW-1","ERVV-1","CGA","CGB5"), rownames(mr.TE.day8@scale.data))
                                                    ,  complete.cluster.day8.2$order ] ),
                  trace="none",density="none",
                  Colv = F,
                  Rowv = F,
                  key=F,scale="row",
                  dendrogram="none",
                  breaks = palette.breaks,
                  col = colorRampPalette(c("#FF00FF","#000000","#FFFF00"), space="Lab"),
                  margins=c(10.5,6),
                  ColSideColors = myColors.day8[,1][complete.cluster.day8.2$order],
                  labRow=NA,
                  main="Tissue specific genes"
                  
)


rgesrgrseb = FindAllMarkers(mr.TE.day8,only.pos = T)



huge.clust.data = mr.TE@scale.data[c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                              "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                              "PSG1","HSD3B1","CYP19A1","SDC1","INHA",
                              "ERVW-1","ERVV-1","CGA","CGB5"),]

huge.cluster = stats::hclust(as.dist(1-abs(cor((   as.matrix(huge.clust.data)     ),method="spearman"))), 
                                      method="ward.D")

plot(huge.cluster,hang=-1)


gplots::heatmap.2(as.matrix(  mr.TE@scale.data[match(c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                                            "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                                            "PSG1","HSD3B1","CYP19A1","SDC1","INHA",
                                                            "ERVW-1","ERVV-1","CGA","CGB5"), rownames(mr.TE@scale.data))
                                                    ,   ] ),
                  trace="none",density="none",
                  Colv = as.dendrogram(huge.cluster),
                  Rowv = F,
                  key=F,scale="row",
                  #dendrogram=as.dendrogram(huge.cluster),
                  breaks = palette.breaks,
                  col = colorRampPalette(c("#FF00FF","#000000","#FFFF00"), space="Lab"),
                  margins=c(10.5,6),
                  labRow=NA,
                  main="Tissue specific genes",
                  
                  
)


##########################################

table(mr.TE.fullPCA.day7@ident)

table(mr.TE.day8.cc@meta.data$old.ident)

table(mr.TE.subtypeA.day9.10@ident)


DoHeatmap(mr.TE.day8, genes.use = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                                "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                                
                                                "PSG1","CSH1","HSD3B1","CYP19A1","SDC1","INHA",
                                                "ERVW-1","ERVV-1","CGA","CGB5"), group.by ="old.ident" )


DoHeatmap(mr.TE.fullPCA.day9.10, genes.use = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                                            "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                                            
                                                            "PSG1","CSH1","HSD3B1","CYP19A1","SDC1","INHA",
                                                            "ERVW-1","ERVV-1","CGA","CGB5"))


DoHeatmap(mr.TE, genes.use = c("CTNNB1","ITGA6","LRP5","TP63","TEAD4","ELF5","FGFR2","FZD5",
                                                "FN1","HLA-G","MMP2","ITGA5","CD9","ITGA1",
                                                
                                                "PSG1","CSH1","HSD3B1","CYP19A1","SDC1","INHA",
                                                "ERVW-1","ERVV-1","CGA","CGB5"),use.scaled =T)


new.day9.10.markers = FindAllMarkers(mr.TE.subtypeA.day9.10, only.pos = T)

1+

sdcsd = new.day9.10.markers[new.day9.10.markers$p_val_adj<0.05 & new.day9.10.markers$avg_logFC>0.5 ,]

DoHeatmap(mr.TE.subtypeA.day9.10, genes.use = sdcsd$gene)

sdcsd[sdcsd$cluster==2,]

table(data.frame(group = mr.TE@ident, day = (mr.TE@meta.data$embryo.day)))


DoHeatmap(mr.TE, genes.use = polar$polar)

mr.TE.polar.gene.average.exp = AverageExpression(mr.TE, use.scale = T, 
                                                 genes.use = intersect(polar$polar, rownames(mr.TE@scale.data)))





boxplot(mr.TE.polar.gene.average.exp)

mr.TE.polar.gene.average.exp.long.df = data.frame(exp = c(mr.TE.polar.gene.average.exp[,1],
                                                          mr.TE.polar.gene.average.exp[,2],
                                                          mr.TE.polar.gene.average.exp[,3],
                                                          mr.TE.polar.gene.average.exp[,4],
                                                          mr.TE.polar.gene.average.exp[,5],
                                                          mr.TE.polar.gene.average.exp[,6]),
                                                  Cluster = factor(rep(c(0:5), each=122))
                                                  )

ggplot(mr.TE.polar.gene.average.exp.long.df, aes(x=Cluster, y=exp , fill = Cluster )) + 
  geom_boxplot(outlier.shape=NA)+ theme(legend.position = "none")  + ylab("Scaled Expression") + theme_bw()


wilcox.test(mr.TE.polar.gene.average.exp[,6], rowMeans(mr.TE.polar.gene.average.exp[,-6]))
wilcox.test(mr.TE.polar.gene.average.exp[,3], rowMeans(mr.TE.polar.gene.average.exp[,-3]))


##########################
# day , cell type distribution stack plot


PCAPlot(mr.TE)

trophoblast.type.day = ((data.frame(day = as.character(mr.TE@meta.data$embryo.day),  type = as.character(mr.TE@ident)  , stringsAsFactors = F)))

trophoblast.type.day.stackbar.df = data.frame(table(trophoblast.type.day)) 
trophoblast.type.day.stackbar.df$percent = trophoblast.type.day.stackbar.df$Freq / c(93,28,125,168,62)

trophoblast.type.day.stackbar.df$day = factor(trophoblast.type.day.stackbar.df$day, levels = c("6","7","8","9","10"))

ggplot(trophoblast.type.day.stackbar.df, aes(x = day, y = percent, fill = type)) + 
  geom_bar(stat = "identity") +  scale_y_continuous(labels = percent_format()) 





trophoblast.type.day.merged = trophoblast.type.day
trophoblast.type.day.merged$type[!(trophoblast.type.day.merged$type %in% c("4","2","5") )] = "Undiff"

trophoblast.type.day.merged$type[trophoblast.type.day.merged$type == "4"] = "CT"
trophoblast.type.day.merged$type[trophoblast.type.day.merged$type == "2"] = "EVT"
trophoblast.type.day.merged$type[trophoblast.type.day.merged$type == "5"] = "ST"

trophoblast.type.day.stack.df = data.frame(table(trophoblast.type.day.merged))

library(scales)

ggplot(trophoblast.type.day.stack.df, aes(x = day, y = Freq, fill = type)) + 
  geom_bar(stat = "identity") 

trophoblast.type.day.stack.df$percent = trophoblast.type.day.stack.df$Freq / c(93,28,125,168,62)

trophoblast.type.day.stack.df$day = factor(trophoblast.type.day.stack.df$day, levels = c("6","7","8","9","10"))

ggplot(trophoblast.type.day.stack.df, aes(x = day, y = percent, fill = type)) + 
  geom_bar(stat = "identity") +  scale_y_continuous(labels = percent_format()) 


######################################

# u-day7 co-day7 hclust





mr.TE.day7 <- NormalizeData(object = mr.TE.day7)
mr.TE.day7 <- FindVariableFeatures(object = mr.TE.day7)
mr.TE.day7 <- ScaleData(object = mr.TE.day7)


day7.average.exp.colcondition = AverageExpression(mr.TE.day7,group.by = "culture.condition" )




head(day7.coculture.onlyembryo.marker.sig)

day7.average.exp.colcondition$high.co =  rownames(day7.average.exp.colcondition) %in% rownames(day7.coculture.onlyembryo.marker.sig)[day7.coculture.onlyembryo.marker.sig$avg_logFC<0]
day7.average.exp.colcondition$high.U =  rownames(day7.average.exp.colcondition) %in% rownames(day7.coculture.onlyembryo.marker.sig)[day7.coculture.onlyembryo.marker.sig$avg_logFC>0]


ggplot(day7.average.exp.colcondition, aes(x= log10(OnlyEmbryo_7+1), y= log10(Coculture_7+1) )) +
  geom_point(color = "grey") + theme_bw() + geom_point(mapping = aes(x = day7.average.exp.colcondition$OnlyEmbryo_7 ,
                                                           y = day7.average.exp.colcondition$Coculture_7 , color = "red", colour) )


plot(y = log10(day7.average.exp.colcondition$Coculture_7+1),
     x = log10(day7.average.exp.colcondition$OnlyEmbryo_7+1) , pch = 16, col = "grey")

points(y = log10(day7.average.exp.colcondition$Coculture_7[day7.average.exp.colcondition$high.co]+1),
       x = log10(day7.average.exp.colcondition$OnlyEmbryo_7[day7.average.exp.colcondition$high.co]+1), pch = 16, col = "red")

points(y = log10(day7.average.exp.colcondition$Coculture_7[day7.average.exp.colcondition$high.U]+1),
       x = log10(day7.average.exp.colcondition$OnlyEmbryo_7[day7.average.exp.colcondition$high.U]+1), pch = 16, col = "blue")



stats::hclust(as.dist(1-abs(cor((   as.matrix(mr.TE.day7@scale.data)     ),method="spearman"))), 
                                         method="ward.D")

mr.TE.day7.unbiased.clust = stats::hclust(as.dist(1-abs(cor((   as.matrix(mr.TE.day7@scale.data[])     ),method="spearman"))), 
                     method="ward.D")
plot(mr.TE.day7.unbiased.clust, hang = -1)







##########################


"WEE1" %in% rownames(day7.coculture.onlyembryo.marker.sig)

DoHeatmap(mr.TE.fullPCA.noday9.10, genes.use = rownames(day7.coculture.onlyembryo.marker.sig)[order(day7.coculture.onlyembryo.marker.sig$avg_logFC)] )


##########################



GO.plot(rownames(day7.coculture.onlyembryo.marker.sig), title = "a")



PCAPlot(mr.TE, group.by = "culture.condition")
PCAPlot(mr.TE, group.by = "embryo.day")
PCAPlot(mr.TE)


############################

FeaturePlot(object = mr.TE, features.plot = c("TBX3","DPPA3","ERVW-1","CGA","CGB5",
                                              "TET1", "TET2","TET3",
                                              "DNMT1",
                                              "DNMT3A","DNMT3B"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca",no.legend = F, pt.size = 1.5)


VlnPlot(object = mr.day9.10, features.plot = c("TBX3","DPPA3","ERVW-1","CGA","CGB5",
                                               "TET1", "TET2","TET3",
                                               "DNMT1",
                                               "DNMT3A","DNMT3B"))


###########################


PCAPlot(mr.TE, group.by = "culture.condition")
PCAPlot(mr.TE, group.by = "embryo.day")
PCAPlot(mr.TE)

mr.TE.day6.day7 = SubsetData(mr.TE.fullPCA.noday9.10, cells.use = rownames(mr.TE.fullPCA.noday9.10@meta.data)[ as.character(mr.TE.fullPCA.noday9.10@meta.data$embryo.day) %in% c("6","7") ])

mr.TE.fullPCA.noday9.10@meta.data$res.0.6

table(data.frame(mr.TE.day6.day7@meta.data$embryo.day,mr.TE.day6.day7@meta.data$culture.condition))

stackplot.day.condition.distribution = data.frame(table(data.frame(   mr.TE.day6.day7@ident,  mr.TE.day6.day7@meta.data$res.0.6  )))[1:6,]

stackplot.day.condition.distribution$percent = stackplot.day.condition.distribution$Freq/c(93,93,93,59,59,59)

stackplot.day.condition.distribution$mr.TE.day6.day7.ident = factor(stackplot.day.condition.distribution$mr.TE.day6.day7.ident, 
                                                                    levels = c("OnlyEmbryo_6","OnlyEmbryo_7","Coculture_7"))

library(scales)
ggplot(stackplot.day.condition.distribution, aes(x = mr.TE.day6.day7.meta.data.res.0.6, y = percent, fill = mr.TE.day6.day7.ident)) + 
  geom_bar(stat = "identity") +  scale_y_continuous(labels = percent_format()) 


two.clusters.pseudotime = data.frame(pseudotime = as.numeric(mr.TE@dr$pca@cell.embeddings[rownames(mr.TE.day6.day7@meta.data),2]),
                                     id = mr.TE.day6.day7@meta.data$res.0.6
                                     )

ggplot(two.clusters.pseudotime[two.clusters.pseudotime$id != "5",], aes(x=id, y=pseudotime, fill = id)) + 
  geom_boxplot()


###########################################

DoHeatmap(mr.TE, genes.use = polar$polar)


FeaturePlot(mr.TE,features.plot = c("TP63","TEAD4","SALL4"),reduction.use = "pca")

########################################


TSC.marker = read.table("/Users/anqin/Box Sync/Projects/Lvbo/JGA_RNAseq/cvs.high.sig")

DoHeatmap(mr.TE,genes.use = rownames(TSC.marker$x))

#######################


FeaturePlot(object = mr.TE, features.plot = c("XIST","SRY"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca",no.legend = F, pt.size = 1.5)




chrY.genes = RC$Geneid[(substr(RC$Chr,1,4)) == "chrY"]



chrY.genes.exp = mr.TE@scale.data[  intersect(chrY.genes,rownames(mr.TE@scale.data) )  ,]

hist(colMeans(chrY.genes.exp))


mr.TE@meta.data$

rowMeans(chrY.genes.exp)

ggplot(as.data.frame(mr.TE@dr$pca@cell.embeddings), aes(x=PC1, y=PC2, color =colMeans(chrY.genes.exp) )) +
  geom_point()+scale_color_gradient2(low="blue", mid="white",
                                     high="red", space ="Lab" )


embryo.gender.by.Y = data.frame(exp = 
colMeans(chrY.genes.exp),

embryo.no = c(
substr(colnames(mr.TE@data)[1:64],1,2),

substr(colnames(mr.TE@data)[65:149],1,4),

substr(colnames(mr.TE@data)[150:447],1,3),

substr(colnames(mr.TE@data)[448:455],1,4),

substr(colnames(mr.TE@data)[456:476],1,3)

),

XIST = mr.TE@scale.data["XIST",]

)





ggplot(embryo.gender.by.Y, aes(x=embryo.no, y=exp)) + 
  geom_boxplot() + xlab("Embryo Number") + ylab("Y-chromosome linked gene exp") +
  geom_abline(intercept =0,slope = 0, color = "red")

ggplot(embryo.gender.by.Y, aes(x=embryo.no, y=XIST)) + 
  geom_boxplot() + xlab("Embryo Number") + ylab("XIST exp") 

###########################

write.table(rownames(mr@meta.data)[which(mr@meta.data$pca.cell.type == "EndoCell")], file = "endo.cell.list.txt",quote=F)

















