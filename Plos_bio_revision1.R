PCAPlot(mr.TE)


VlnPlot(mr.Embryo.300markerPCA.TE,features.plot = c("DAB2","PTGES","TGFBR3","PDGFA"), group.by = "embryo.day",
        nCol = 4)


wilcox.test( mr.Embryo.300markerPCA.TE@data["DAB2", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 6 ],
             mr.Embryo.300markerPCA.TE@data["DAB2", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 7 ]
             )
10^(mean(mr.Embryo.300markerPCA.TE@data["DAB2", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 6 ])/
             mean(mr.Embryo.300markerPCA.TE@data["DAB2", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 7 ]))



wilcox.test( mr.Embryo.300markerPCA.TE@data["PTGES", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 6 ],
             mr.Embryo.300markerPCA.TE@data["PTGES", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 7 ]
)

10^(mean(mr.Embryo.300markerPCA.TE@data["PTGES", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 6 ])/
      mean(mr.Embryo.300markerPCA.TE@data["PTGES", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 7 ]))



wilcox.test( mr.Embryo.300markerPCA.TE@data["TGFBR3", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 6 ],
             mr.Embryo.300markerPCA.TE@data["TGFBR3", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 7 ]
)
10^(mean(mr.Embryo.300markerPCA.TE@data["TGFBR3", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 6 ])/
      mean(mr.Embryo.300markerPCA.TE@data["TGFBR3", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 7 ]))



wilcox.test( mr.Embryo.300markerPCA.TE@data["PDGFA", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 6 ],
             mr.Embryo.300markerPCA.TE@data["PDGFA", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 7 ]
)
10^(mean(mr.Embryo.300markerPCA.TE@data["PDGFA", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 6 ])/
      mean(mr.Embryo.300markerPCA.TE@data["PDGFA", mr.Embryo.300markerPCA.TE@meta.data$embryo.day == 7 ]))


###################################


PCAPlot(mr.Embryo.300markerPCA, group.by = "Phase")

mr.Embryo.allgene.PCA = RunPCA(mr.Embryo.300markerPCA, pc.genes = rownames(mr.Embryo.300markerPCA@data ))

PCAPlot(mr.Embryo.allgene.PCA, group.by = "Phase")

mr.Embryo.allgene.PCA@meta.data$Phase[mr.Embryo.allgene.PCA@meta.data$Phase == "Untypical"] = "TE"

PCAPlot(mr.Embryo.allgene.PCA, group.by = "Phase")
PCAPlot(mr.Embryo.allgene.PCA, group.by = "embryo.day")

#####################################

# heatmap
# cell.marker.table

DoHeatmap(mr.Embryo.allgene.PCA, genes.use = cell.marker.table$V1, group.by = "Phase")

mr.Embryo.allgene.PCA

dim(mr.Embryo.allgene.PCA@data[cell.marker.table$V1,])

mr.Embryo.allgene.PCA = SetAllIdent(mr.Embryo.allgene.PCA, id = "lineage")

mr.Embryo.allgene.PCA.ave.expression = AverageExpression(mr.Embryo.allgene.PCA, use.scale = T,
                  genes.use = intersect(cell.marker.table$V1, rownames(mr.Embryo.allgene.PCA@scale.data)) )


library(gplots)

palette.breaks <- seq(-3, 3, 0.1)
color.palette = colorRampPalette(c("dodgerblue4","dodgerblue1","white","firebrick1","firebrick3"), 
                                 space="Lab")(length(palette.breaks) - 1)

heatmap.2(as.matrix( ( mr.Embryo.allgene.PCA.ave.expression ) ),trace = "none",density = "none",
          Colv = NA,
          Rowv = NA,
          #ColSideColors = c("grey","red")[ factor(type[complete.cluster$order] ) ] ,
          #ColSideColors =  c("grey","red")[factor(type)],
          col=color.palette,
          #breaks = palette.breaks,
          scale = c("row"),
          dendrogram = "both")
#dev.off()





##############################

FeaturePlot(mr.Embryo.300markerPCA, features.plot = c("GATA2","GATA3"), reduction.use = "pca", 
            cols.use = c("grey","blue"), no.legend = F, pt.size = 2)

#############################
library(dendextend)

library(RColorBrewer)
n <- 7
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

plot(mr.TE.day7.unbiased.clust, hang= -1)

colored_bars( col_vector[as.numeric(factor(substr(mr.TE.day7.unbiased.clust$labels,1,3)))][mr.TE.day7.unbiased.clust$order]   )




stats::hclust(as.dist(1-abs(cor((   as.matrix(mr.TE.day7@scale.data)     ),method="spearman"))), 
              method="ward.D")

mr.TE.day7.polar.clust = stats::hclust(as.dist(1-abs(cor((   as.matrix(mr.TE.day7@scale.data[intersect(rownames(mr.TE.day7@scale.data), polar$polar),])     ),method="spearman"))), 
                                                                                 method="ward.D")
plot(mr.TE.day7.polar.clust, hang = -1)

colored_bars( col_vector[as.numeric(factor(substr(mr.TE.day7.polar.clust$labels,1,3)))][mr.TE.day7.polar.clust$order]   )


###########################################################

pca.mr.TE = mr.TE@dr$pca
pca.mr.eigValues = (pca.mr.TE@sdev)^2  ## EigenValues
pca.mr.varExplained = pca.mr.eigValues / sum(pca.mr.eigValues)


pca.mr.df <- data.frame(pca.mr.varExplained,
                 PC = factor(c(1:40))
                 )

ggplot(data=pca.mr.df[1:10,], aes(x=PC, y=pca.mr.varExplained*100)) +
  geom_bar(stat="identity", fill = "grey") + xlab("Dimensions") + ylab("% of variance explained")+
  geom_line(data=pca.mr.df[1:10,], aes(x=c(1:10), y=pca.mr.varExplained*100), colour="blue")+
  geom_point(colour="blue")

##########################################################


mr.TE@data

PCAPlot(mr.Embryo.300markerPCA)

library(FactoMineR)
library(ggplot2)
library(dplyr)
library(ggrepel)


temp.aaa = t(   as.data.frame(as.matrix(mr.Embryo.300markerPCA@data)))


pca.res = PCA(  temp.aaa[,  na.omit(match(cell.marker.table$V1, colnames(temp.aaa)))  ]   ,graph = F)

PC1 <- as.numeric(pca.res$ind$coord[,1])
PC2 <- as.numeric(pca.res$ind$coord[,2])




PCs <- data.frame(PC1,PC2, lineage = mr.Embryo.300markerPCA@meta.data$lineage )

#### day, using only single cell samples
P<-ggplot(PCs, aes(PC1,PC2) , color = lineage  )
P +geom_point( size = 3) + 
  xlab(paste("PC1",as.character(round(pca.res$eig[,2][1],2)),"%")) + 
  ylab(paste("PC2",as.character(round(pca.res$eig[,2][2],2)),"%")) +
  ggtitle(   "PCA on scRNAseq"    ) + theme_bw() 


ggplot(PCs, aes(x=PC1, y=PC2, color=lineage)) + geom_point()

####################################

FeaturePlot(object = mr.TE, features.plot = c("TBX3"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "pca")

VlnPlot(object = mr.TE, features.plot = c("TBX3"))

PCAPlot(object = mr.TE)

PCAPlot(object = mr.TE, group.by = "embryo.day")

VlnPlot(object = mr.TE, features.plot = c("TBX3"), group.by = "embryo.day")

###################################


VlnPlot(object = mr.Embryo.allgene.PCA, features.plot = c("TBX3"), group.by = "Phase")

###################################

write.csv(file = "STCTEVT_day_table.csv", table(mr.TE@meta.data[,c("embryo.day","res.0.6")]))


################################


PCAPlot(object = mr.TE)

co.u.day7.lab = paste(mr.TE@meta.data$culture.condition, mr.TE@meta.data$embryo.day) 

co.u.day7.lab[!(co.u.day7.lab %in% c("Coculture 7", "OnlyEmbryo 7"))] = "Other"

names(co.u.day7.lab) = rownames(mr.TE@meta.data)

mr.TE = AddMetaData(mr.TE, metadata = co.u.day7.lab, col.name = "co.u.day7.lab")

PCAPlot(mr.TE, group.by = "co.u.day7.lab")


##################################################

# add GO plot for u- co-day7 cells

day7.coculture.onlyembryo.marker

day7.coculture.onlyembryo.marker.sig

GO.plot(rownames(day7.coculture.onlyembryo.marker.sig)[day7.coculture.onlyembryo.marker.sig$avg_logFC > 0], title = "a")

sum(day7.coculture.onlyembryo.marker$p_val_adj< 0.05 & day7.coculture.onlyembryo.marker$avg_logFC <0 )

GO.plot(rownames(day7.coculture.onlyembryo.marker)[day7.coculture.onlyembryo.marker$p_val_adj< 0.05 & day7.coculture.onlyembryo.marker$avg_logFC <0], title = "Co-day7 high genes")

GO.plot(rownames(day7.coculture.onlyembryo.marker.sig)[day7.coculture.onlyembryo.marker.sig$p_val_adj< 0.05 & day7.coculture.onlyembryo.marker.sig$avg_logFC > 0], title = "Co-day7 low genes")
