require(cluster)
require(factoextra)
require(rlist)
require(princurve)


#require(FactoMineR)

opt = list()

opt$pc.comp = 5
opt$pc.clust = c(1:3)

###

dat = list()
#dat$exp_matrix = mr.TE.fullPCA@scale.data
dat$exp_matrix = mr.TE.fullPCA@data

dat$time = mr.TE.fullPCA@meta.data$embryo.day # levels should be ranked by time, increasingly


write.csv(mr.TE.fullPCA@scale.data, file = "mr.TE.fullPCA.scale.data.csv",quote = F)

write.csv(mr.TE.fullPCA@meta.data, file = "embryo.dat.csv",quote = F)


##### Gap statistic method

# exp_matrix = 
# exp_matrix = mr.RPKM
# perform PCA

pca.res = FactoMineR::PCA(t(  log2(  dat$exp_matrix +1)  ),graph = F)

# pca.res = FactoMineR::PCA(t(  (  dat$exp_matrix )  ),graph = F)

PC1 <- as.numeric(pca.res$ind$coord[,1])
PC2 <- as.numeric(pca.res$ind$coord[,2])
plot(PC1,PC2, col=mr.TE.fullPCA@meta.data$embryo.day)

library(rgl)

plot3d(pca.res$ind$coord[,1:3], col=mr.TE.fullPCA@meta.data$embryo.day, size=3)

set.seed(19930426)
gap_stat <- clusGap(pca.res$ind$coord[,1:3], FUN = kmeans, nstart = 25,
                    K.max = 10, B = 500)
fviz_gap_stat(gap_stat)

####################################
# for each time point, cut into clusters 

# tp = 10

tree.df = data.frame(level = 0, uniq.id = 0, number.of.cells = 0) # common ancestor

sample_name_list = list()

set.seed(0)

id.counter = 1
level.counter = 1


for(tp in levels(dat$time)){
  
  # tp = 8
  
  pc_score_tp = pca.res$ind$coord[dat$time == tp,opt$pc.clust]
  
  
  gap_stat <- clusGap(pc_score_tp, FUN = kmeans, nstart = 25,
                      K.max = 10, B = 500)
  
  # fviz_gap_stat(gap_stat)
  
  
  # first.min.cluster.num = min(which(  gap_stat$Tab[-(length(gap_stat$Tab[,"gap"])),"gap"] - gap_stat$Tab[-1,"gap"] <=0     ) - 1)
  
  first.min.cluster.num = maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstSEmax")
  
  #findpeaks(gap_stat$Tab)
  
  print(first.min.cluster.num)
  
  kmeans.res = kmeans(pc_score_tp, first.min.cluster.num, nstart = 25)
  
  print(tp)
  print(unique(kmeans.res$cluster))
  
  
  for(k in unique(kmeans.res$cluster)){
    
    # add node data
    
    node.samples = names(which(kmeans.res$cluster == k))
    
    # print(node.samples)
    
    temp = data.frame(level = level.counter, uniq.id = id.counter, number.of.cells = length(node.samples) )
    
    tree.df = rbind(tree.df, temp)
    
    sample_name_list[[id.counter]] = node.samples
    
    id.counter = id.counter + 1
    
  }
  
  level.counter = level.counter + 1
}


tree.df = tree.df[-1,]

########################
# construct mean expression matrix 

mean_exp_matrix = matrix(NA, ncol = length(sample_name_list), nrow = nrow(dat$exp_matrix) )

colnames(mean_exp_matrix) = c(1:(length(sample_name_list)))

rownames(mean_exp_matrix) = rownames(dat$exp_matrix)

for(node.id in c(1:(length(sample_name_list)))){
  
  #node.id = 1
  mean_exp_matrix[,node.id] = rowMeans(dat$exp_matrix[,match(sample_name_list[[node.id]], colnames(dat$exp_matrix))])
  
  
  
}


# construct path string

tree.df$path.strings = rep("1", nrow(tree.df))

#path.strings[tree.df$uniq.id[tree.df$level==1]] = "1"

# cor_mean_exp_matrix = cor(mean_exp_matrix,method = "spearman")

cor_mean_exp_matrix = lsa::cosine(mean_exp_matrix)

for(l in  (sort(unique(tree.df$level))[-1])  ){ # 每层循环
  
  
  
  for(s in tree.df$uniq.id[tree.df$level==l]){ # 对一个node，选择上层node中和他correlation最高的那一个
    
    
    parent.nodes = tree.df$uniq.id[tree.df$level==(l-1)]
    
    biggest.cor.so.far = -2
    biggest.p.node = NA
    for(p in parent.nodes){
      
      if(cor_mean_exp_matrix[p,s] > biggest.cor.so.far){
        biggest.cor.so.far = cor_mean_exp_matrix[p,s]
        biggest.p.node = p
      }
      
    }
    
    print(s)
    print(biggest.p.node)
    print(biggest.cor.so.far)
    
    # tree.df[s,'path.strings'] = paste(tree.df[s,'path.strings'], biggest.p.node, sep = "/")
    tree.df[s,'path.strings'] = biggest.p.node
    
  }
  
  

  
}

tree.df$path.strings = as.numeric(tree.df$path.strings)

build.path.string <-  function(x){
  

    
  path.string = ""
    
  pointer = x
  
  if(pointer ==1 ){
    return("1")
  }
  
  while(T){
  
    #print(tree.df[pointer,2])
    path.string = paste(tree.df[pointer,2], path.string,sep="/")
    pointer = tree.df[pointer,4]
  
  if(pointer == 1) {break} 
  
  }
  
  
  path.string = gsub(".$", "", path.string)
  path.string = gsub("^", "1/", path.string)
  return(path.string)
}
      

 build.path.string(1)      

tree.df$pathString = sapply(tree.df$uniq.id, build.path.string)


############ 转换成tree

tree.df.Tree <- as.Node(tree.df)

print(tree.df.Tree, "number.of.cells",limit = 20)

plot(tree.df.Tree)


#############################################

# select a branch, generate pc-curve, and look into gene expression 


get.branch.data <- function(branch.node.lists){
#branch.node.lists = c(1,2,3,4,7)

samples.in.branch = c()

for(n in branch.node.lists){
  
  samples.in.branch = c(samples.in.branch, sample_name_list[[n]])
}


branch.pcspace.to.curve = pca.res$ind$coord[  match(samples.in.branch,rownames(pca.res$ind$coord) ) ,opt$pc.clust]

branch.pc.curve  = principal.curve(branch.pcspace.to.curve, smoother="smooth.spline")

branch.exp = dat$exp_matrix[, match(  samples.in.branch  ,colnames(dat$exp_matrix) )]

return(list(branch.pc.curve, branch.exp, samples.in.branch, branch.pcspace.to.curve))

}


##################

b1 = get.branch.data(c(1,2,3,4,7))
b2 = get.branch.data(c(1,2,3,5,6))

##################

plot(b1[[4]][,1:2],pch=20)
points(b2[[4]][,1:2],col="red",pch=20)

lines(b1[[1]])
lines(b2[[1]])

plot(x=b1[[1]]$lambda,y=  b1[[2]][match( "TBX3",rownames(b1[[2]]) ),]   )
points(x=b2[[1]]$lambda,y=  b2[[2]][match( "TBX3",rownames(b2[[2]]) ),]  ,col="red" )


######################################
# for each node, determine variance of a gene

copmute.gene.node.var <- function(gene.name){
  
  var.vect = c()
  
  for(node.id.temp in tree.df$uniq.id){
  
    #node.id.temp = 4
    
    mr.te.temp = SubsetData(mr.TE, cells.use = sample_name_list[[node.id.temp]]  )  
  # node.gene.exp = dat$exp_matrix[  match(gene.name,rownames(dat$exp_matrix))    , match(  sample_name_list[[node.id.temp]]  ,colnames(dat$exp_matrix) )]
    mr.te.temp = FindVariableGenes(mr.te.temp, do.plot = F)
    
    dispersion = mr.te.temp@hvg.info[match( gene.name, rownames(mr.te.temp@hvg.info)     ),2]
  
    dispersion
    
    var.vect = c(var.vect, dispersion)
  
  }
  return(var.vect)

}

copmute.gene.node.var("TBX3")

copmute.gene.node.var("ERVW-1")

copmute.gene.node.var("CGB5")

copmute.gene.node.var("TET1")

copmute.gene.node.var("ACTB")

PCAPlot(mr.TE)

PCAPlot(mr.TE,group.by = "embryo.day")

plot.tree.df =  tree.df
plot.tree.df$dispersion = copmute.gene.node.var("TBX3")

plot(x=c(1:5),y=plot.tree.df$dispersion[c(1,2,3,4,7)],type="l",ylim=c(-3,3))
lines(x=c(1:5),y=plot.tree.df$dispersion[c(1,2,3,5,6)],col="red")



copmute.gene.node.exp <- function(gene.name){
  
  var.vect = c()
  
  for(node.id.temp in tree.df$uniq.id){
    
    #node.id.temp = 4
    
    # mr.te.temp = SubsetData(mr.TE, cells.use = sample_name_list[[node.id.temp]]  )  
    # node.gene.exp = dat$exp_matrix[  match(gene.name,rownames(dat$exp_matrix))    , match(  sample_name_list[[node.id.temp]]  ,colnames(dat$exp_matrix) )]
    #mr.te.temp = FindVariableGenes(mr.te.temp, do.plot = F)
    
    #dispersion = mr.te.temp@hvg.info[match( gene.name, rownames(mr.te.temp@hvg.info)     ),2]
    
    #dispersion
    
    exp.dat = dat$exp_matrix[, colnames(dat$exp_matrix) %in% sample_name_list[[node.id.temp]] ]
    
    
    
    var.vect = c(var.vect, mean(exp.dat[match(gene.name, rownames(exp.dat)),]))
    
  }
  return(var.vect)
  
}


#############

plot.tree.df$EXP = copmute.gene.node.exp("TBX3")
plot(x=c(1:5),y=plot.tree.df$EXP[c(1,2,3,4,7)],type="l",ylim=c(-1,3))
lines(x=c(1:5),y=plot.tree.df$EXP[c(1,2,3,5,6)], col="red")

plot.tree.df$dispersion = copmute.gene.node.var("TBX3")
plot(x=c(1:5),y=plot.tree.df$dispersion[c(1,2,3,4,7)],type="l",ylim=c(-3,3))
lines(x=c(1:5),y=plot.tree.df$dispersion[c(1,2,3,5,6)],col="red")

###############

plot.tree.df$EXP = copmute.gene.node.exp("ITGA5")
plot(x=c(1:5),y=plot.tree.df$EXP[c(1,2,3,4,7)],type="l", ,ylim=c(0,2))
lines(x=c(1:5),y=plot.tree.df$EXP[c(1,2,3,5,6)], col="red")

plot.tree.df$EXP = copmute.gene.node.exp("LRP5")
plot(x=c(1:5),y=plot.tree.df$EXP[c(1,2,3,4,7)],type="l",ylim=c(0,0.2))
lines(x=c(1:5),y=plot.tree.df$EXP[c(1,2,3,5,6)], col="red")


plot.tree.df$EXP = copmute.gene.node.exp("GATA3")
plot(x=c(1:5),y=plot.tree.df$EXP[c(1,2,3,4,7)],type="l")
lines(x=c(1:5),y=plot.tree.df$EXP[c(1,2,3,5,6)], col="red")


plot.tree.df$dispersion = copmute.gene.node.var("GATA3")
plot(x=c(1:5),y=plot.tree.df$dispersion[c(1,2,3,4,7)],type="l")
lines(x=c(1:5),y=plot.tree.df$dispersion[c(1,2,3,5,6)],col="red")


###########



plot.tree.df$EXP = copmute.gene.node.exp("TFAP2C")
plot(x=c(1:5),y=plot.tree.df$EXP[c(1,2,3,4,7)],type="l",ylim=c(-1,3))
lines(x=c(1:5),y=plot.tree.df$EXP[c(1,2,3,5,6)], col="red")


plot.tree.df$dispersion = copmute.gene.node.var("TFAP2C")
plot(x=c(1:5),y=plot.tree.df$dispersion[c(1,2,3,4,7)],type="l",ylim=c(-3,3))
lines(x=c(1:5),y=plot.tree.df$dispersion[c(1,2,3,5,6)],col="red")


# 
# # print(gap_stat, method = "firstmax")
# 
# 
# 
# 
# kmeans.res$centers
# 
# kmeans.res$cluster
# 
# plot(pc_score_tp[,1:2], col=kmeans.res$cluster,pch=20)
# 
# xx = list(list(c(1:2)),list(list(1:4,4:5)))
# library(data.tree)
# 
# yy = as.Node(xx)
# yy
# plot(yy)
# 
# print(yy, "level")
# 
# yy$Get('level', traversal = "level")



###################################
# 4.64, 5.01

t1x = rnorm(n = 100,mean = 0,sd = 1)
t1y =rnorm(n = 100,mean = 0,sd = 1)

t2x = rnorm(n = 100,mean = 0,sd = 2)
t2y = rnorm(n = 100,mean = 0,sd = 2)

t3x = c(rnorm(n = 100,mean = -3,sd = 1.5) , rnorm(n = 100,mean = 5,sd = 1))
t3y = c(rnorm(n = 100,mean = -3,sd = 1.5) , rnorm(n = 100,mean = 5,sd = 1))

plot(t1x, t1y, ylim=c(-10,10), xlim=c(-10,10),pch=19, col = "red")
plot(t2x, t2y, ylim=c(-10,10), xlim=c(-10,10),pch=19, col = "orange")

plot(    t3x , 
         
         t3y , pch=19, col=rep(c("gold","green"), each= 100))
