library(mclust)

data(diabetes)
class <- diabetes$class
table(class)
## class
## Chemical   Normal    Overt 
##       36       76       33
X <- diabetes[,-1]
head(X)
##   glucose insulin sspg
## 1      80     356  124
## 2      97     289  117
## 3     105     319  143
## 4      90     356  199
## 5      90     323  240
## 6      86     381  157
clPairs(X, class)


BIC <- mclustBIC(X)

plot(BIC)

mod1 <- Mclust(X, x = BIC)
summary(mod1, parameters = TRUE)
#
## sspg    -2486.208 -10393.00   2217.533

plot(mod1, what = "classification")

#######

irisBIC <- mclustBIC(mr.TE.fullPCA.day9.10@dr$pca@cell.embeddings[,c(1,2)], G = 1:2)
irisBIC
plot(irisBIC)
mod1iris <- Mclust(mr.TE.fullPCA.day9.10@dr$pca@cell.embeddings[,c(1,2)], x = irisBIC)
summary(mod1iris, parameters = TRUE)
plot(mod1iris, what = "classification")
plot(mod1iris, what = "density")
############



irisBIC <- mclustBIC(mr.TE.fullPCA.day8@dr$pca@cell.embeddings[,c(1,2)], G = 1:2)
irisBIC
plot(irisBIC)
mod1iris <- Mclust(mr.TE.fullPCA.day8@dr$pca@cell.embeddings[,c(1,2)], x = irisBIC)
summary(mod1iris, parameters = TRUE)
plot(mod1iris, what = "classification")
plot(mod1iris, what = "density")


irisBIC <- mclustBIC(mr.TE.fullPCA.day7@dr$pca@cell.embeddings[,c(1,2)], G = 1:2)
irisBIC
plot(irisBIC)
mod1iris <- Mclust(mr.TE.fullPCA.day7@dr$pca@cell.embeddings[,c(1,2)], x = irisBIC)
summary(mod1iris, parameters = TRUE)
plot(mod1iris, what = "classification")
plot(mod1iris, what = "density")

##########
library(cluster)
library(factoextra)

fviz_nbclust(mr.TE.fullPCA.day9.10@dr$pca@cell.embeddings[,c(1,2)], kmeans, method = "silhouette")
fviz_nbclust(mr.TE.fullPCA.day8@dr$pca@cell.embeddings[,c(1,2)], kmeans, method = "silhouette")
fviz_nbclust(mr.TE.fullPCA.day7@dr$pca@cell.embeddings[,c(1,2)], kmeans, method = "silhouette")

set.seed(123)
gap_stat <- clusGap(mr.TE.fullPCA.day8@dr$pca@cell.embeddings[,c(1,3)], FUN = kmeans, nstart = 25,
                    K.max = 10, B = 500)
fviz_gap_stat(gap_stat)

library("scatterplot3d") # load

scatterplot3d(mr.TE.fullPCA@dr$pca@cell.embeddings[,c(1:3)], pch = 16, color=mr.TE.fullPCA@meta.data$embryo.day)

library(rgl)

plot3d(mr.TE.fullPCA@dr$pca@cell.embeddings[,c(1:3)], col=mr.TE.fullPCA@meta.data$embryo.day, size=3)



library(destiny)
dm <- DiffusionMap( t(mr.TE.fullPCA@scale.data) )
plot(dm, col = c("red","green","black","blue","orange")[(mr.TE.fullPCA@meta.data$embryo.day)] )


head(mr.TE.fullPCA@raw.data)

mr.TE.fullPCA = RunPCA(mr.TE.fullPCA, pc.genes = rownames(mr.TE.fullPCA@scale.data))

PCAPlot(mr.TE.fullPCA, group.by = "embryo.day")

Six.group.markers = FindAllMarkers(mr.TE)
Six.group.markers = Six.group.markers[Six.group.markers$p_val_adj<0.05,]

ICAPlot(mr.TE.fullPCA)

RunPCA()

?RunICA

###############################

apply(mr.RPKM, MARGIN = 1, function(x){sum(x>1)>20}  )
sum(apply(mr.RPKM, MARGIN = 1, function(x){sum(x>1)>20}  ))

dim(mr.RPKM)

mr.RPKM.clean = mr.RPKM[rowMeans(mr.RPKM) > 1, match(    colnames(mr.TE@scale.data),   colnames(mr.RPKM)     )]

rpkm.pc.res = FactoMineR::PCA(t(   log10(mr.RPKM.clean+1)   ),graph = F)


plot(rpkm.pc.res$ind$coord[,1], rpkm.pc.res$ind$coord[,2], col = mr.TE@meta.data$embryo.day)



library(rgl)

plot3d(rpkm.pc.res$ind$coord[,1:3], col=mr.TE.fullPCA@meta.data$embryo.day, size=3)

write.csv(mr.RPKM.clean,file = "mr.RPKM.clean.csv",  quote=F)


Gene