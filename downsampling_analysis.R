
library(mclust)


# PCAPlot(mr.TE)

accuracy.list = c()
subsamp.portion.list = rep(seq(from=0.9, to=0.3, by = -0.1), each=100)

mr.TE.forSampling = mr.TE


for(subsamp.portion in subsamp.portion.list){
  
  
  set.seed(as.numeric(format(Sys.time(), "%OS3"))*1000)
  
  sampled.index = sample(1:476, size = round(476*subsamp.portion))
  
  #print(sampled.index)
  
  mr.TE.forSampling.sampled = SubsetData(mr.TE.forSampling,
                                                               cells.use = rownames(mr.TE.forSampling@meta.data)[sampled.index])
  
  mr.TE.forSampling.sampled = FindClusters(mr.TE.forSampling.sampled,  reduction.type = "pca", dims.use = 1:10, 
                                           resolution = 0.6, print.output = 0,  force.recalc = T)
  

  
  # PCAPlot(Xy.cc.noD26.RPC.no.rePCA.forsubSampling.sampled, group.by = "res.0.6")
  
  ###########
  # compute clustering accuracy
  
  t.m.orig = mr.TE.forSampling@meta.data[sampled.index,]
  t.m.orig = t.m.orig[order(rownames(t.m.orig)),]
  origional.cluster = as.numeric(t.m.orig$res.0.6)
  
  t.m.subsampled = mr.TE.forSampling.sampled@meta.data
  t.m.subsampled = t.m.subsampled[order(rownames(t.m.subsampled)),]
  subsampled.cluster = as.numeric(t.m.subsampled$res.0.6)
  
  accuracy.list = c(accuracy.list , adjustedRandIndex(origional.cluster,subsampled.cluster))
  
  print(accuracy.list)
}

gc()

subsamp.portion.list.to.plot = data.frame(subsamp.portion.list ,
                                          accuracy.list)


new.row = data.frame(subsamp.portion.list = 1, accuracy.list = 1)

subsamp.portion.list.to.plot2 = rbind(subsamp.portion.list.to.plot, new.row)

ggplot(subsamp.portion.list.to.plot, aes(x=factor(subsamp.portion.list), y=(subsamp.portion.list.to.plot$accuracy.list*100)^0.5 * 10)) + 
  geom_boxplot() + ylim(c(0,1))


ggplot(subsamp.portion.list.to.plot2, aes(x=factor(subsamp.portion.list), y=accuracy.list, 
                                          fill = factor(subsamp.portion.list))) + 
  geom_boxplot(outlier.shape = NA)  + ylab("Adj Rank Index") + theme(legend.position = "none")  + xlab("Subsampling portion")


#########################################




