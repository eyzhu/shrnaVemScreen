library("edgeR")
library("GenomicRanges")
library("Rsamtools")
library("ggplot2")

setwd("/home/eliot/Dropbox/Projects/shRNA_vem")

##load("raw_counts_p3.Rdata")
load("raw_counts_p3_w_dups.Rdata")

run1F = "barcodes_p3.txt"

##run1F = 
run2F = "shRNA_run2_UP.txt"

run1 = read.csv(run1F, sep="\t", header=F)
run2 = read.csv(run2F, sep="\t", header=F)

run1[,1] = gsub("_.*", "", run1[,1])
run2[,1] = gsub("_.*", "", run2[,1])

run1 = run1[,-2]
run2 = run2[,-2]

finalCounts = counts

#rmvInd = c(which(colnames(finalCounts)=="23"), which(colnames(finalCounts)=="39"))
#finalCounts = finalCounts[,-rmvInd]

sampNames = colnames(finalCounts)
hpNames = rownames(finalCounts)

sampInfo = run1[match(sampNames, run1[,1]),]
sampInfo[,2] = sub("shRNA","pool",sampInfo[,2])

shRNAInfo = read.csv("shRNA_info.txt", header=T, sep="\t")
shRNAInfo = shRNAInfo[,1:3]
shRNAInfo[,1] = as.character(shRNAInfo[,1])
shRNAInfo = unique(shRNAInfo)
##matchedSHRNAinfo = shRNAInfo[match(hpNames, shRNAInfo[,1]),] 

## make edgeR object
dataObjs = list()
mypools = as.character(unique(sampInfo[,2]))
mypools = sub("shRNA","pool",mypools)

# qq = c()
# for (pool in mypools){
#         pool = mypools[1] 
#         sel = which(sampInfo[,2]==pool)
#         qq=cbind(qq, as.character(currSampleInfo[,3]))
# }

WXthresh = 0.05

# masterCountMat = matrix(0, ncol=length(sampList), nrow=length(hpList))
# rownames(masterCountMat) = hpList
# colnames(masterCountMat) = sampList
# 
# hpList= c()
# sampList = c()

myGenes = c()

#goodInds = match(c(1:16), colnames(finalCounts))

#goodCounts = finalCounts[ ,goodInds]

tropical = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=19)

## mm = log2(currER$counts[,13]+1) - log2(currER$counts[,14]+1)
## aa = log2(currER$counts[,13]+1) + log2(currER$counts[,14]+1)

for (pool in mypools){
  
  if (pool=="pool6"){next}
  ## select columns 
  ## if(!(pool %in% c("pool9","pool6","pool11"))){next}
  ## pool = mypools[1]
  ## print(pool)
  ##pool = mypools[10]
  ##pool = "pool8"
  
  ##pool="shRNA3"
  pool="pool3_1"
  
  sel = which(sampInfo[,2]==pool)
  currSampleInfo = sampInfo[sel,]
  
  currCount = finalCounts[,sel]
  ## select for theorhetical hairpins
  totsel = shRNAInfo[which(shRNAInfo[,3]==pool), 1]
  sampSel = totsel[totsel%in%rownames(currCount)]
  currCount = currCount[sampSel, ]
  
  ##currCount = goodCounts
  currER = new("DGEList") 
  currER$counts = currCount
  
  ## remove < 0.5 cpm
  sel = rowSums(cpm(currER$counts)>0.5)>=3
  currER = currER[sel,]
  
  ## the code below is only used to create sapce for the master matrix, but needs
  ## to commented out when actually filling the matrix
  # hpList = c(hpList, rownames(currER$counts))
  # sampList = c(sampList, colnames(currER$counts))
  
  ## fill out matrix
  cnames = colnames(currER$counts)
  rnames = rownames(currER$counts)
  # masterCountMat[rnames, cnames] = masterCountMat[rnames, cnames] + currER$counts 
  
  #######
  currER$samples = data.frame(SampleID=colnames(currER$counts),
                              group=as.factor(currSampleInfo[,3]),
                              lib.size=colSums(currER$counts))                                    
  currER$samples$norm.factors=1
  currNames = rownames(currER$counts)        
  currER$genes = shRNAInfo[match(currNames, shRNAInfo[,1]), ]
  ## barplots 
  #barplot(colSums(currER$counts), las=2)
  # tiff(filename=paste("./figures/hairpin_bar_plots/","hp_bp_", pool, ".tiff", sep=""),
  #      width = 6, height = 6, units = 'in', res = 250)
  # barplot(rowSums(currER$counts), main=pool, xaxt='n')
  # dev.off()
  ##
  print(pool)
  currER$counts = cpm(currER$counts)
  
  ##currER = calcNormFactors(currER, method="upperquartile")  
  # tiff(filename=paste("./figures/MDS/","MDS_", pool, ".tiff", seIp=""),
  #      width = 6, height = 6, units = 'in', res = 250)
  # MDS = plotMDS(currER, labels=currER$samples$group, 
  #         main=pool, 
  #         col=c(rep(1,sum(currSampleInfo[,3]==unique(currSampleInfo[,3])[1] )),
  #               rep(2,sum(currSampleInfo[,3]==unique(currSampleInfo[,3])[2])))
  #         )
  # dev.off()
  
  ## heatmaps
  myCountsMat = currER$counts
  colnames(myCountsMat) = currSampleInfo[,3]
  VemMat = myCountsMat[,currSampleInfo[,3]=="VEM"]
  vemSums = rowSums(VemMat)
  myCountsMatSorted = myCountsMat[order(vemSums),] 
  
  # tiff(filename=paste("./figures/heatmaps/","heatmap_", pool, ".tiff", sep=""),
  #      width = 6, height = 6, units = 'in', res = 250)
  # 
  # heatmap(myCountsMatSorted, Colv=NA, Rowv = NA, labRow = NA, main=pool)
  # 
  # dev.off()
  # 
  # ## sample bar plot
  # tiff(filename=paste("./figures/Sample_bar_plots/","sample_bp_", pool, ".tiff", sep=""),
  #      width = 6, height = 6, units = 'in', res = 250)
  # barplot(colSums(myCountsMat), main=pool)
  # dev.off()
  # 
  #colnames(myCountsMat) = currSampleInfo[,3]
  
  ### make heatmap
  VemSel = which(currSampleInfo[,3]=="VEM")
  VehSel = which(currSampleInfo[,3]=="VEH")
  
  # llVem = length(VemSel)
  # llVeh = length(VehSel)
  # 
  # vemNorm = llVem*(llVem+1)/2
  # vehNorm = llVeh*(llVeh+1)/2
  # 
  # rankedCounts = apply(myCountsMat, 1, rank)
  # ranksums = apply(rankedCounts[VehSel,],2,sum)-vehNorm
  # ranksums_vem = apply(rankedCounts[VemSel,],2,sum)-vemNorm
  # qq = sapply(ranksums_vem, function (x) sum(ranksums>x)/length(ranksums))
  
  # hist(ranksums)
  
  wilcoxResult = 
    apply(myCountsMat, 1, 
          function(x) wilcox.test( x[VemSel], x[VehSel], paired=F, alternative="greater") )
  
  pvals = sapply(wilcoxResult, '[[', 3)
  
  hist(pvals, nclass=20)
  
  pvals = pvals[order(pvals)]
  pvals = pvals[pvals<WXthresh]
  
  
  ##topPvals = head(pvals, 200)
  topPvals = data.frame(pvals)
  topPvalsNamed = cbind(topPvals, shRNAInfo[match(rownames(topPvals), shRNAInfo[,1]),2:3])
  
  myGenes = rbind(myGenes, topPvalsNamed)
  
}
  ##as.character(tppGenes[,1])%in%as.character(topPvalsNamed)
myGenes[,2] = as.character(myGenes[,2])
geness = unique(myGenes[,2])
  
finalOut = c()

for (gene in geness){
  ind = which(myGenes[,2]==gene)
  poolNums = data.frame(table(myGenes[ind,3]))
  outLine = c(gene, poolNums[,2])
  finalOut = rbind(finalOut, outLine)
}

finalOutt = rbind(c("gene", levels(poolNums[,1])), finalOut  )

outMat = matrix(NA, nrow=length(geness), ncol=length(unique(myGenes[,3])) )
rownames(outMat) = geness
colnames(outMat) = unique(myGenes[,3])

for (ii in 1:nrow(myGenes)){
  cGene = as.character(myGenes[ii, 2])
  cPool = as.character(myGenes[ii, 3])
  
  if( is.na(outMat[cGene, cPool]) ){ 
    outMat[cGene, cPool] = myGenes[ii, 1]  
  } else {
    outMat[cGene, cPool] = paste(outMat[cGene, cPool], myGenes[ii, 1], sep=";")  
  }
    
}

finalOutMat = cbind(rownames(outMat), outMat)
finalOutMat = rbind(c("gene", colnames(outMat)), finalOutMat)
write.table(finalOutMat, file="gene_hp_pvals_0.05_14AUG17.txt", sep="\t"
            ,col.names=F, row.names=F )


######################### print 0.05 - 0.1
myGenes[,2] = as.character(myGenes[,2])
myGenes = myGenes[myGenes[,1]>0.05,]
geness = unique(myGenes[,2])

finalOut = c()

for (gene in geness){
  ind = which(myGenes[,2]==gene)
  poolNums = data.frame(table(myGenes[ind,3]))
  outLine = c(gene, poolNums[,2])
  finalOut = rbind(finalOut, outLine)
}

finalOutt = rbind(c("gene", levels(poolNums[,1])), finalOut  )

outMat = matrix(NA, nrow=length(geness), ncol=length(unique(myGenes[,3])) )
rownames(outMat) = geness
colnames(outMat) = unique(myGenes[,3])

for (ii in 1:nrow(myGenes)){
  cGene = as.character(myGenes[ii, 2])
  cPool = as.character(myGenes[ii, 3])
  
  if( is.na(outMat[cGene, cPool]) ){ 
    outMat[cGene, cPool] = myGenes[ii, 1]  
  } else {
    outMat[cGene, cPool] = paste(outMat[cGene, cPool], myGenes[ii, 1], sep=";")  
  }
  
}

finalOutMat = cbind(rownames(outMat), outMat)
finalOutMat = rbind(c("gene", colnames(outMat)), finalOutMat)
write.table(finalOutMat, file="gene_hp_pvals_0.05to0.1_14AUG17.txt", sep="\t"
            ,col.names=F, row.names=F )

## MDS
  # tiff(filename=paste("./figures/MDS/","MDS_", pool, ".tiff", sep=""), 
  #      width = 6, height = 6, units = 'in', res = 250)
  # plotMDS(currER, labels=currER$samples$group, 
  #         main=pool, 
  #         col=c(rep(1,sum(currSampleInfo[,3]==unique(currSampleInfo[,3])[1] )),
  #               rep(2,sum(currSampleInfo[,3]==unique(currSampleInfo[,3])[2])))
  # )
  # dev.off() 
  
  # des = model.matrix(~currER$samples$group)
  # des = des[,-2]
  # xglm = estimateDisp(currER, des)
  # 
  # plotBCV(xglm, cex=0.5)        
  # ##
  # 
  # ##http://www.nonlinear.com/support/progenesis/comet/faq/v2.0/pq-values.aspx--explain q value
  # fit = glmFit(xglm, des)
  # lrt = glmLRT(fit,2)
  # plotSmear(lrt)
  # thresh = 0.5
  # lfc = 1
  # top15 = topTags(lrt, n=15, sort.by="logFC")
  # topp = topTags(lrt, n=Inf, sort.by="logFC")
  # topGenes = unique(as.character(topp$table[topp$table$FDR<thresh & topp$table$logFC>lfc,2]))        
  # 
  # myTPgenes = c(myTPgenes, topGenes)
  
  ## smear plot 
  
  ## CPM vs LFC 
# myTPgenes = as.data.frame(table(myTPgenes))
# tppGenes = myTPgenes[which(myTPgenes[,2]>1),]

