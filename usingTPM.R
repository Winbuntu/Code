TPM = as.data.frame(fread("../public_placenta_single_cellRNAseq/GSE89497_Human_Placenta_TMP_V2.txt",header = T))

Seu.pub <- CreateSeuratObject(raw.data = TPM, min.cells = 30, min.genes = 500, 
                              project = "public placenta")

Seu.pub <- AddMetaData(object = Seu.pub, metadata = placenta.metadata )

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

PCAPlot(object = Seu.pub, dim.1 = 1, dim.2 = 2,group.by = "placenta.sample.cluster",do.label=T)
