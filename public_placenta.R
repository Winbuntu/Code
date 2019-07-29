
library(data.table)
library(Seurat)

library(plyr) 

pub.placenta.RC = as.data.frame(fread("../public_placenta_single_cellRNAseq/READ_COUNT_table_public_placenta_gencode.txt",header = T,
                                      stringsAsFactors = F))


pub.placenta.RC.clean = pub.placenta.RC[,-c(1:6)]

colnames(pub.placenta.RC.clean) = substr(colnames(pub.placenta.RC.clean),1,10  )


# 构建 annotation information
sra1 = read.csv("../public_placenta_single_cellRNAseq/sra_result.csv",header = T,stringsAsFactors = F)

sra2 = read.csv("../public_placenta_single_cellRNAseq/SraRunInfo.csv",header = T,stringsAsFactors = F)

placenta.annot = merge(x = sra1, y = sra2, by.x = "Experiment.Accession", by.y = "Experiment")

write.csv(placenta.annot, "../public_placenta_single_cellRNAseq/placenta.annot.csv",quote=F)


placenta.annot.clean = placenta.annot[match(colnames(pub.placenta.RC.clean) , placenta.annot$Run),]

placenta.sample.term = unlist(lapply(strsplit(  placenta.sample.name  ,split = "_",fixed = T), function(x) {x[[1]]}))

placenta.sample.cell_type = substr(unlist(lapply(strsplit(  placenta.sample.name  ,split = "_",fixed = T), function(x) {x[[2]]})),1,3)

placenta.sample.cluster = paste(placenta.sample.term, placenta.sample.cell_type, sep= "_")

placenta.metadata = data.frame(term = placenta.sample.term, cell_type = placenta.sample.cell_type)
rownames(placenta.metadata) = placenta.annot.clean$Run

############################

pub.placenta.RC2 = as.data.frame(fread("../public_placenta_single_cellRNAseq/READ_COUNT_table_public_placenta_pooled_samples_gencode.txt",
                                       header = T,stringsAsFactors = F))

rownames(pub.placenta.RC2) = pub.placenta.RC2$Geneid

pub.placenta.RC2.clean = pub.placenta.RC2[,-c(1:6)]



placenta2.sample.name = unlist(lapply(strsplit(  colnames(pub.placenta.RC2.clean)  ,split = "_R1",fixed = T), function(x) {x[[1]]}))

placenta2.sample.term = paste(substr(unlist(lapply(strsplit(  colnames(pub.placenta.RC2.clean)  ,split = "_",fixed = T), function(x) {x[[1]]})),1,3),"W",sep="")

placenta2.sample.cell_type = substr(unlist(lapply(strsplit(  placenta2.sample.name  ,split = "_",fixed = T), function(x) {x[[2]]})),1,3)

placenta2.metadata = data.frame(term = placenta2.sample.term, cell_type = placenta2.sample.cell_type)
rownames(placenta2.metadata) = colnames(pub.placenta.RC2.clean)



###################################################
placenta.metadata.merged = rbind(placenta.metadata, placenta2.metadata)

placenta.RC.merged = cbind(pub.placenta.RC.clean,pub.placenta.RC2.clean )



#pub.placenta.RC = cbind(pub.placenta.RC,pub.placenta.RC2)


#pub.placenta.RC = pub.placenta.RC[,-c(1:6)]


#########
# 去除mt-rna, ERCC


MT.index.pub.placenta = which(substr(placenta.RC.merged$Geneid,1,3) == "MT-")

ERCC.index.pub.placenta = which(substr(placenta.RC.merged$Geneid,1,5) == "ERCC-")

placenta.RC.merged.clean = placenta.RC.merged[-c(MT.index,ERCC.index,ERCC.index.pub.placenta),]

#rownames(RC.clean) = rownames(RC)[-c(MT.index,ERCC.index)]


#############

##############################


Seu.pub <- CreateSeuratObject(raw.data = placenta.RC.merged.clean, min.cells = 30, min.genes = 500, 
                         project = "public placenta")

Seu.pub <- AddMetaData(object = Seu.pub, metadata = placenta.metadata.merged )



VlnPlot(object = Seu.pub, features.plot = c("nGene", "nUMI"), nCol = 3)

Seu.pub <- FilterCells(object = Seu.pub, subset.names = c("nGene","nUMI"), 
                  low.thresholds = c(3000,0), high.thresholds = c(Inf,  2000000))


Seu.pub <- NormalizeData(object = Seu.pub, normalization.method = "LogNormalize", 
                    scale.factor = 100000)

Seu.pub <- FindVariableGenes(object = Seu.pub, mean.function = ExpMean, dispersion.function = LogVMR,do.plot = F)

length(x = Seu.pub@var.genes)

date()
Seu.pub <- ScaleData(object = Seu.pub, vars.to.regress = c("nUMI"))
date()

###############################

# 整体的Seurat object 构建完成

#############################

Seu.pub <- RunPCA(object = Seu.pub, pc.genes = Seu.pub@var.genes, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5,pcs.compute = 40)

PCAPlot(object = Seu.pub, dim.1 = 1, dim.2 = 2,group.by = "term",do.label=T)

PCAPlot(object = Seu.pub, dim.1 = 1, dim.2 = 2,group.by = "cell_type")

PCElbowPlot(object = Seu.pub)

Seu.pub <- JackStraw(object = Seu.pub, num.replicate = 100,num.pc = 20)

JackStrawPlot(object = Seu.pub, PCs = 1:20) # 这里，选择PC12

Seu.pub <- RunTSNE(Seu.pub, dims.use = 1:9)

TSNEPlot(object = Seu.pub, dim.1 = 1, dim.2 = 2,group.by = "cell_type",do.label=T)
TSNEPlot(object = Seu.pub, dim.1 = 1, dim.2 = 2,group.by = "term",do.label=T)


Combined.RC = cbind(RC.clean,pub.placenta.RC.clean)


###############################################################


Seu.Merge.two.study = MergeSeurat(object1 = Seu.pub, object2 = mr.TE,min.cells = 20, do.normalize = T, scale.factor = 1000000)

Seu.Merge.two.study <- FindVariableGenes(object = Seu.Merge.two.study, mean.function = ExpMean, dispersion.function = LogVMR,do.plot = F)

length(x = Seu.Merge.two.study@var.genes)

date()
Seu.Merge.two.study <- ScaleData(object = Seu.Merge.two.study, vars.to.regress = c("nUMI"))
date()

###############################

# 整体的Seurat object 构建完成

#############################

hk = read.csv("../public_placenta_single_cellRNAseq/human_housekeeping.csv",header = F, stringsAsFactors = F)

Seu.Merge.two.study <- RunPCA(object = Seu.Merge.two.study, do.print = TRUE, pcs.print = 1:5, 
                  genes.print = 5,pcs.compute = 40, pc.genes = hk$V1 )

Seu.Merge.two.study =  AddMetaData(Seu.Merge.two.study, Seu.Merge.two.study@dr$pca@cell.embeddings[,1],"PC1_regress")


Seu.Merge.two.study <- ScaleData(object = Seu.Merge.two.study, vars.to.regress = c("nUMI","PC1_regress"))

Seu.Merge.two.study <- RunPCA(object = Seu.Merge.two.study, do.print = TRUE, pcs.print = 1:5, 
                              genes.print = 5,pcs.compute = 40)


PCAPlot(Seu.Merge.two.study, group.by = "orig.ident", dim.1 = 1,dim.2 = 2)

FeaturePlot(Seu.Merge.two.study, features.plot = c("CGB5","ERVW-1","HLA-G","ITGA6","ITGA5","TBX3","DPPA3","GAPDH","JUN","TBX3","DNMT3A"), 
            dim.1 = 2,dim.2 = 3, reduction.use = "pca")





PrintPCA(object = Seu.Merge, pcs.print = 1, genes.print = 50, use.full = FALSE)

VlnPlot(Seu.Merge, features.plot = "CHMP2A",group.by = "orig.ident")


#####################################


g.1 <- head(rownames(mr.TE@hvg.info), 1000)
g.2 <- head(rownames(Seu.Merge@hvg.info), 1000)

genes.use <- unique(c(g.1, g.2))

Seu.Merge.two.study <- RunCCA(mr.TE, Seu.pub, genes.use = genes.use, num.cc = 30)

DimPlot(object = Seu.Merge.two.study, reduction.use = "cca", group.by = "orig.ident", 
        pt.size = 0.5, do.return = TRUE)

VlnPlot(object = Seu.Merge.two.study, features.plot = "CC1", group.by = "orig.ident", 
        do.return = TRUE)

MetageneBicorPlot(Seu.Merge.two.study, grouping.var = "orig.ident", dims.eval = 1:30, 
                  display.progress = FALSE)


Seu.Merge.two.study <- AlignSubspace(Seu.Merge.two.study, reduction.type = "cca", grouping.var = "orig.ident", 
                                 dims.align = 1:20)

VlnPlot(object = Seu.Merge.two.study, features.plot = "ACC1", group.by = "orig.ident", 
              do.return = TRUE)

Seu.Merge.two.study <- RunTSNE(Seu.Merge.two.study, reduction.use = "cca.aligned", dims.use = 1:9, 
                           do.fast = T)

TSNEPlot(Seu.Merge.two.study, do.return = T, pt.size = 0.5, group.by = "orig.ident")

FeaturePlot(Seu.Merge.two.study, features.plot = c("CGB5","ERVW-1","HLA-G","ITGA6","ITGA5","TBX3","DPPA3","GAPDH","JUN"), 
            dim.1 = 1,dim.2 = 2, reduction.use = "cca.aligned")


DimPlot(object = Seu.Merge.two.study, reduction.use = "cca.aligned", group.by = "orig.ident", 
        pt.size = 0.5, do.return = TRUE)
