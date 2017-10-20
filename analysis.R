library("edgeR")
library("GenomicRanges")
library("Rsamtools")
library("ggplot2")

setwd("/home/eliot/Dropbox/Projects/shRNA_vem")

load("raw_counts_13Sep17.Rdata")
##load("raw_counts_7Aug17.Rdata")

# ## compare data w/ duplicates w/ data w/o duplicates
# finalCounts = oCount
# sel = rowSums(cpm(finalCounts)>0.5)>=3
# oldCounts = finalCounts[sel,] 
# pp = colSums(oldCounts)
# 
# finalCounts = nCount
# sel = rowSums(cpm(finalCounts)>0.5)>=3
# newCounts = finalCounts[sel,]
# qq = colSums(newCounts)
# 
# oldNames = names(pp)
# newNames = names
run1F = "shRNA_run1_sample_info_UP.csv"
run2F = "shRNA_run2_UP.txt"

run1 = read.csv(run1F, sep="\t", header=T)
run2 = read.csv(run2F, sep="\t", header=F)

run1[,1] = gsub("_.*", "", run1[,1])
run2[,1] = gsub("_.*", "", run2[,1])

run1 = run1[,-2]
run2 = run2[,-2]

cNames = colnames(counts)
finalCounts = counts[ ,cNames %in% run1[,1]]

## 23 and 39 (vem outlier) removed, 23 removed b/c 1 control w/ 40 tags w/ >50000 counts
## 164 is outlier
# rmvInd = c(which(colnames(finalCounts)=="23"), which(colnames(finalCounts)=="39"), which(colnames(finalCounts)=="164") )
rmvInd = c(which(colnames(finalCounts)=="39"), which(colnames(finalCounts)=="164") )

finalCounts = finalCounts[,-rmvInd]

### sum the counts from same pellet
pelletNums = read.csv("shRNA Barcodes Nov 18 2016.csv", header=F, sep="\t")
pelletNums = pelletNums[-grep("Pilot", pelletNums[,3]),]

## trbshoot code
# pelletNums[which(!(pelletNums[,1] %in% colnames(finalCountsBk))), ]
pelletNums2 = pelletNums[(pelletNums[,1] %in% colnames(finalCounts)), ]
pelletNums2[,3] = as.character(pelletNums2[,3])

plateIDL = unlist(strsplit(pelletNums2[,3], " "))
plateIDL = plateIDL[-grep("NEW", plateIDL)]
plateIDs = matrix(plateIDL, nrow=nrow(pelletNums2), byrow = T)

pelletNums3 = cbind(pelletNums2[,1], plateIDs)
pelletNums3 = pelletNums3[, -c(2,4)]

pelletNums3[,3] = sub( "\\..*", "", pelletNums3[,3] )
pelletNums3[,4] = sub(",", "", pelletNums3[,4])

# ## compare plate nums and sequencing nums
# pelletComp = pelletNums3[,c(1:2)]
# pelletComp[,2] = sub("#", "shRNA", pelletComp[,2])
# 
# sampInfoComp = run1[match(pelletComp[,1], run1[,1]),]
# pelletComp[which(pelletComp[,2]!=sampInfoComp[,2]),]

pelletNums4 = cbind(pelletNums3[,1], apply(pelletNums3[,c(2:4)], 1, paste, collapse = "_"  ))
## need to make correction
pelletNums4[which(pelletNums4[,1]=="63"),2] = "#3_3_VEM"

unqPlates = unique(pelletNums4[,2])

finalCountsBk = finalCounts

finalCounts = c()

for(unqPlate in unqPlates){
  ##unqPlate = unqPlates[45]
  exps = pelletNums4[pelletNums4[,2]==unqPlate, 1]
  
  if( length(exps)==1 ){
    finalCounts = cbind(finalCounts, finalCountsBk[,exps]*2 )
    print(exps)
  } else {
    finalCounts = cbind(finalCounts, rowSums(finalCountsBk[,exps]) )
  }
  
  colnames(finalCounts)[ncol(finalCounts)] = exps[1]
  
}

## need to adjust counts for orphan plates

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

WXthresh = 0.5

# masterCountMat = matrix(0, ncol=length(sampList), nrow=length(hpList))
# rownames(masterCountMat) = hpList
# colnames(masterCountMat) = sampList
# 
# hpList= c()
# sampList = c()

myGenes = c()
myRawCounts = c()

for (pool in mypools){
  
  if (pool=="pool6"){next}
  ## select columns 
  ## if(!(pool %in% c("pool9","pool6","pool11"))){next}
  ## pool = mypools[1]
  ## print(pool)
  ##pool = mypools[10]
  pool = "pool3"
  # pool="pool3"  
  sel = which(sampInfo[,2]==pool)
  currSampleInfo = sampInfo[sel,]
  
  currCount = finalCounts[,sel]
  ## select for theorhetical hairpins
  totsel = shRNAInfo[which(shRNAInfo[,3]==pool), 1]
  sampSel = totsel[totsel%in%rownames(currCount)]
  currCount = currCount[sampSel, ]
  
  currER = new("DGEList") 
  currER$counts = cpm(currCount)
  ##currER$counts = cpm(currER$counts)
  ## remove < 0.5 cpm
  sel = rowSums(currER$counts>0.5)>=2
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
  
  ##currER = calcNormFactors(currER, method="upperquartile")  
  tiff(filename=paste("./figures/MDS/","MDS_", pool, ".tiff", seIp=""),
       width = 6, height = 6, units = 'in', res = 250)
  MDS = plotMDS(currER, labels=currER$samples$group,
          main=pool,
          col=c(rep(1,sum(currSampleInfo[,3]==unique(currSampleInfo[,3])[1] )),
                rep(2,sum(currSampleInfo[,3]==unique(currSampleInfo[,3])[2])))
          )
  dev.off()
  
  vehInds = which(currSampleInfo[,3]=="VEH");
  vemInds = which(currSampleInfo[,3]=="VEM");
  
  plotCountMat = log2(currER$counts+1)[,c(vehInds, vemInds)]
  
  boxplot( plotCountMat, main=pool, col=c(rep(2,6),rep(1,6)), ylab="log2(counts per million)", xlab="Sample Number")
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
  rawCounts = myCountsMat[,c(VehSel,VemSel)]
  
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
  IDnames = rownames(topPvalsNamed)
  fcTab = myCountsMat[IDnames,]
  if( is.null(dim(fcTab)) ){
    fc = log2(mean(fcTab[VemSel])/(mean(fcTab[VehSel])+1))
  }else{
    fc = log2(rowMeans(fcTab[,VemSel])/(rowMeans(fcTab[,VehSel])+1))
  }
  
  pvalsF = cbind(fc,topPvalsNamed)
  pvalsF2 = pvalsF[order(-pvalsF$fc), ]
  
  myRawCounts = rbind(myRawCounts, rawCounts[rownames(pvalsF2),])
  
  myGenes = rbind(myGenes, pvalsF2)
  ##myGenes = rbind(myGenes, topPvalsNamed)
  
}

myGenes[,3] = as.character(myGenes[,3])
geneFreq = data.frame(table(myGenes[,3]))
dupGenes = geneFreq[geneFreq[,2]>1,1]

dupInds = myGenes[,3] %in% dupGenes
unqInds = !(myGenes[,3] %in% dupGenes) 

dupOut = myGenes[dupInds,]
dupOut2 = dupOut[order(dupOut[,"gene"]),]

unqOut = myGenes[unqInds,]
unqOut2 = unqOut[order(unqOut[,"pvals"]), ]

finalOut = rbind(dupOut2, unqOut2)
finalOut = cbind(tag=rownames(finalOut), finalOut[, c(3,1,2,4)], myRawCounts[rownames(finalOut),])

fcThresh = 10000

fcOut = finalOut[finalOut[,"fc"]>=8, ]
fcOutFiltered = fcOut[fcOut[,9]>=fcThresh | fcOut[,10]>=fcThresh | fcOut[,11]>=fcThresh ,]

valThresh = 50000

colBool = rowSums(cbind(fcOut[,9]>valThresh, fcOut[,10]>valThresh, fcOut[,11]>valThresh))
##validationOut = fcOut[colBool>=2, ]
validationOut = finalOut[(finalOut[,9]>valThresh & finalOut[,10]>valThresh & finalOut[,11]>valThresh), ]

valGenes = unique(validationOut[,"gene"])

## overlap w/ a375 vem study
a375study = read.csv("./supp_data_1.csv", sep="\t",header=T)
underExp = a375study[as.numeric(a375study[,4]) < -2,]
genes = as.character(underExp[,3])
valGenes[valGenes %in% genes]

##fcOut = finalOut[finalOut[,"fc"]>=5, ]

### String analysis
library(STRINGdb)

string_db = STRINGdb$new()
# mapped_genes = unique(string_db$map(validationOut, "gene"))
mapped_genes = unique(string_db$map(fcOut, "gene"))
subsetMapped = unique(mapped_genes[,c("gene","STRING_id")])
hits = subsetMapped[,2]

enrichMats = c()

for (i in 1:length(hits) ) {
  
  print(i)
  hit = hits[i]
  
  result = try( string_db$get_neighbors(hit) )
  
  if (class(result) == "try-error"){
    print(subsetMapped[i,1])
    next
  }
  
  interactors = string_db$get_neighbors(hit)
  interactions = string_db$get_interactions(c(hit, interactors))
  
  strongInteractions = interactions[ interactions[,"combined_score"]>=300, ] 
  primInteractions = 
    strongInteractions[(strongInteractions[,1]==hit | strongInteractions[,2]==hit), ]
  
  ## sort genes and take top 10
  primInteractionsSorted = primInteractions[order(-primInteractions[,"combined_score"]), ]
  
  if (nrow(primInteractionsSorted)>=10){
    primInteractionsSorted = primInteractionsSorted[1:10,]
  }
  
  primGenesUnq = unique(c(primInteractionsSorted[,1], primInteractionsSorted[,2]))
  
  # secondInteractions = 
  #   strongInteractions[( strongInteractions[,1] %in% primGenesUnq | strongInteractions[,2] %in% primGenesUnq), ]
  
  ##strongInteractorsList = unique(c(strongInteractions[,1], strongInteractions[,2]))
  ##strongInteractorsList = unique(c(second))
  
  # tiff(filename=paste("./figures/string_db_analysis/","STRING_", subsetMapped[i,1], ".tiff", sep=""),
  #      width = 15, height = 15, units = 'in', res = 250)
  # 
  # string_db$plot_network(primGenesUnq)
  # dev.off()
  
  if(length(primGenesUnq)>0){
  
    enrichProcess = string_db$get_enrichment(primGenesUnq, category = "Process", methodMT = "fdr", iea=T)
    enrichComponent = string_db$get_enrichment(primGenesUnq, category = "Component", methodMT = "fdr", iea=T)
    enrichFunction = string_db$get_enrichment(primGenesUnq, category = "Function", methodMT = "fdr", iea=T)
    
    enrichProcess = enrichProcess[1:5, c("pvalue_fdr", "term_description")]
    enrichComponent = enrichComponent[1:5, c("pvalue_fdr", "term_description")]
    enrichFunction = enrichFunction[1:5, c("pvalue_fdr", "term_description")]
    
    enrichMat = rbind(enrichProcess, enrichFunction, enrichComponent)
    enrichMat[,1] = -log10(enrichMat[,1])
    enrichMat = cbind(enrichMat, myType=factor(c( rep("Process", 5), rep("Function", 5), rep("Component", 5))), myGene = subsetMapped[i,1])
  
  }
  
  enrichMats = rbind(enrichMats, enrichMat)
  
  # 
  # nameOrder = enrichMat$term_description[order(enrichMat$myType, enrichMat$pvalue_fdr)]
  # 
  # enrichMat$term_description = factor(enrichMat$term_description, levels=nameOrder)
  
  # tiff(filename=paste("./figures/string_db_analysis/","GO_", subsetMapped[i,1], ".tiff", sep=""), 
  #      width = 8, height = 4, units = 'in', res = 250)
  
  # print(ggplot(enrichMat, aes(x=pvalue_fdr, y=term_description )) +
  #   geom_segment(aes(yend=term_description), xend=0, colour="grey50") + 
  #   geom_point(size=3, aes(colour=myType)) + 
  #   scale_colour_brewer(palette="Set1", limits=c("Process", "Function", "Component"), guide=FALSE)  +
  #   theme_bw()+ theme( panel.grid.major.y = element_blank() ) + 
  #   labs(x="-log10(pvalue_fdr)", y="", title=subsetMapped[i,1]) + 
  #   facet_grid(myType ~ ., scales="free_y", space="free_y"))
  # 
  # dev.off()
  
}



save(enrichMats, file="enrichMats.Rdata")
  ##string_db$plot_ppi_enrichment(hits, quiet=T)

library(ggplot2)
library(reshape2)

enrichMatsFilt = enrichMats[enrichMats[,"myGene"] %in% fcOutFiltered[,"gene"], ]
enrichMatsFilt = enrichMatsFilt[!is.na(enrichMatsFilt[,1]),]
enrichMatsFilt = enrichMatsFilt[enrichMatsFilt[,1]>3, ]

enrichMatsFiltProcess = enrichMatsFilt[which(enrichMatsFilt[,"myType"]=="Process"),]
enrichMatsFiltFunction = enrichMatsFilt[which(enrichMatsFilt[,"myType"]=="Function"),]
enrichMatsFiltComponent = enrichMatsFilt[which(enrichMatsFilt[,"myType"]=="Component"),]

processCounts = table(enrichMatsFiltProcess[,"term_description"])
processCounts = processCounts[order(-processCounts)]
pCounts = processCounts[1:10]

pGenes = unique(enrichMatsFiltProcess[enrichMatsFiltProcess[,"term_description"] %in% names(pCounts), "myGene"])

pMat = matrix(0,length(pCounts),length(pGenes))
colnames(pMat) = pGenes
rownames(pMat) = names(pCounts)
## build counts matrix
for (i in 1:nrow(pMat)){
  currGenes = as.character(enrichMatsFiltProcess[enrichMatsFiltProcess[,"term_description"]==names(pCounts[i]), "myGene"])
  
  pMat[i,currGenes] = enrichMatsFiltProcess[enrichMatsFiltProcess[,"term_description"]==names(pCounts[i]), "pvalue_fdr"]
  
}

pPlotDf = melt(pMat)

tiff(filename=paste("./figures/process", ".tiff", sep=""), width = 10, height = 10, units = 'in', res = 250)
     
print(ggplot(pPlotDf, aes(Var2, Var1)) + geom_tile(aes(fill=value), color="black") + scale_fill_gradient(low="white", high="steelblue") + 
  ylab("Biological Process") + xlab("Genes") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) + coord_equal())

dev.off()

functionCounts = table(enrichMatsFiltFunction[,"term_description"])
functionCounts = functionCounts[order(-functionCounts)]
fCounts = functionCounts[1:10]

fGenes = unique(enrichMatsFiltFunction[enrichMatsFiltFunction[,"term_description"] %in% names(fCounts), "myGene"])

fMat = matrix(0,length(fCounts),length(fGenes))
colnames(fMat) = fGenes
rownames(fMat) = names(fCounts)
## build counts matrix
for (i in 1:nrow(fMat)){
  currGenes = as.character(enrichMatsFiltFunction[enrichMatsFiltFunction[,"term_description"]==names(fCounts[i]), "myGene"])
  fMat[i,currGenes] = enrichMatsFiltFunction[enrichMatsFiltFunction[,"term_description"]==names(fCounts[i]), "pvalue_fdr"]
}

fPlotDf = melt(fMat)

tiff(filename=paste("./figures/function", ".tiff", sep=""), width = 10, height = 10, units = 'in', res = 250)

print(ggplot(fPlotDf, aes(Var2, Var1)) + geom_tile(aes(fill=value), color="black") + scale_fill_gradient(low="white", high="steelblue") + 
        ylab("Biological Function") + xlab("Genes") +
        theme(axis.text.x = element_text(angle=45, hjust=1)) + coord_equal())

dev.off()

componentCounts = table(enrichMatsFiltComponent[,"term_description"])
componentCounts = componentCounts[order(-componentCounts)]
cCounts = componentCounts[1:10]

cGenes = unique(enrichMatsFiltComponent[enrichMatsFiltComponent[,"term_description"] %in% names(cCounts), "myGene"])

cMat = matrix(0,length(cCounts),length(cGenes))
colnames(cMat) = cGenes
rownames(cMat) = names(cCounts)
## build counts matrix
for (i in 1:nrow(cMat)){
  currGenes = as.character( enrichMatsFiltComponent[enrichMatsFiltComponent[,"term_description"]==names(cCounts[i]), "myGene"])
  cMat[i, currGenes] = enrichMatsFiltComponent[enrichMatsFiltComponent[,"term_description"]==names(cCounts[i]), "pvalue_fdr"]
}

cPlotDf = melt(cMat)

tiff(filename=paste("./figures/component", ".tiff", sep=""), width = 10, height = 10, units = 'in', res = 250)

print(ggplot(cPlotDf, aes(Var2, Var1)) + geom_tile(aes(fill=value), color="black") + scale_fill_gradient(low="white", high="steelblue") + 
        ylab("Biological Component") + xlab("Genes") +
        theme(axis.text.x = element_text(angle=45, hjust=1)) + coord_equal())

dev.off()


## plot most enriched pathways
enrichMatsFiltOrder = enrichMatsFilt[order(-enrichMatsFilt[,1]),]

enrichMatMax = c()
unqGenes = unique(enrichMatsFiltOrder[,"myGene"])

for (gene in unqGenes){
  
  currMat = enrichMatsFiltOrder[enrichMatsFiltOrder[,"myGene"]==gene, ]
  
  maxPro = currMat[currMat[,"myType"]=="Process",]
  maxPro = maxPro[which.max(maxPro[,1]),]
  
  maxFun = currMat[currMat[,"myType"]=="Function",]
  maxFun = maxFun[which.max(maxFun[,1]),]
  
  maxCom = currMat[currMat[,"myType"]=="Component",]
  maxCom = maxCom[which.max(maxCom[,1]),]
  
  enrichMatMax = rbind(enrichMatMax, maxPro, maxFun, maxCom)
    
}

nameOrder = enrichMatMax$term_description[order(enrichMatMax$myType, enrichMatMax$pvalue_fdr)]
enrichMatMax$term_description = factor(enrichMatMax$term_description, levels=nameOrder)

print(ggplot(enrichMatMax[enrichMatMax[,"myType"]=="Process",], aes(x=pvalue_fdr, y=term_description )) +
  geom_segment(aes(yend=term_description), xend=0, colour="grey50") +
  geom_point(size=3, aes(colour=myType)) +
  scale_colour_brewer(palette="Set1", limits=c("Process", "Function", "Component"), guide=FALSE)  +
  theme_bw()+ theme( panel.grid.major.y = element_blank() ) +
  labs(x="-log10(pvalue_fdr)", y="") +
  facet_grid(myType ~ ., scales="free_y", space="free_y"))

dev.off()



##processCountsNum = as.numeric(processCounts)
##processInds = which(processCountsNum>2)



## plot pathways
# library("GOplot")
# 
# colNamesenrichMatsFilt
# chord = chord_dat()
# 
# GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
#   
## check expression of genes
library(biomaRt)

ensembl = useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl")

expression1 = read.table("./entrez_expression/ENCFF391EQZ.tsv", sep="\t", header = T)
expression2 = read.table("./entrez_expression/ENCFF742WED.tsv", sep="\t", header = T)

## get only rows w/ ensemble gene IDs
expression1 = expression1[grep("ENSG", expression1[,1]),]
expression2 = expression2[grep("ENSG", expression2[,1]),]

## remove version from geneIDs
expression1[,1] = sub("\\..*","",expression1[,1])
expression2[,1] = sub("\\..*","",expression2[,1])

## convert ensemble gene symbols to IDs
ensmblGenes = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = validationOut[,"gene"], mart = ensembl)
goodGenes = ensmblGenes[[1]]

expin1 = cbind(ensmblGenes[[2]], expression1[match(goodGenes, expression1[,1]), c("gene_id", "FPKM")])
expin2 = cbind(ensmblGenes[[2]], expression2[match(goodGenes, expression2[,1]), c("gene_id", "FPKM")])

##string_db$plot_network(hits)

##finalOut2 = cbind 

##finalOut = finalOut[-which(finalOut[,"fc"]<0), ]

# write.table(finalOut, file="raw_gene_list_19Sep17.txt", sep="\t"
#             ,col.names=T, row.names=F, quote=F )
# 
# finalOutFilt = finalOut[finalOut[,"pvals"]<=.1,]
# 
# # for (i in 1:nrow(finalOutFilt)){
# #   
# # }
# 
# write.table(finalOutFilt, file="filt_gene_list_19Sep17.txt", sep="\t"
#             ,col.names=T, row.names=F, quote=F )

#######################################################################3
# rawCountsOut = finalCounts[,as.character(c(pool3samps, pool2samps))]
# 
# finalOut2 = cbind(finalOut, rawCountsOut[rownames(finalOut),])
# 
# wantPools = c("pool2","pool3")
# wantshRNAInfo= shRNAInfo[shRNAInfo[,3]%in%wantPools,]

##rawCounts = finalCounts[rownames(finalOut),]
# for (gene in geness){
#   ind = which(myGenes[,3]==gene)
#   poolNums = data.frame(table(myGenes[ind,3]))
#   outLine = c(gene, poolNums[,2])
#   finalOut = rbind(finalOut, outLine)
# }
# 
# finalOutt = rbind(c("gene", levels(poolNums[,1])), finalOut  )
# 
# outMat = matrix(NA, nrow=length(geness), ncol=length(unique(myGenes[,3])) )
# rownames(outMat) = geness
# colnames(outMat) = unique(myGenes[,3])
# 
# for (ii in 1:nrow(myGenes)){
#   cGene = as.character(myGenes[ii, 2])
#   cPool = as.character(myGenes[ii, 3])
#   
#   if( is.na(outMat[cGene, cPool]) ){ 
#     outMat[cGene, cPool] = myGenes[ii, 1]  
#   } else {
#     outMat[cGene, cPool] = paste(outMat[cGene, cPool], myGenes[ii, 1], sep=";")  
#   }
#     
# }
# 
# finalOutMat = cbind(rownames(outMat), outMat)
# finalOutMat = rbind(c("gene", colnames(outMat)), finalOutMat)




# for (gene in geness){
#   ind = which(myGenes[,2]==gene)
#   poolNums = data.frame(table(myGenes[ind,3]))
#   outLine = c(gene, poolNums[,2])
#   finalOut = rbind(finalOut, outLine)
# }
# 
# finalOutt = rbind(c("gene", levels(poolNums[,1])), finalOut  )
# 
# outMat = matrix(NA, nrow=length(geness), ncol=length(unique(myGenes[,3])) )
# rownames(outMat) = geness
# colnames(outMat) = unique(myGenes[,3])
# 
# for (ii in 1:nrow(myGenes)){
#   cGene = as.character(myGenes[ii, 2])
#   cPool = as.character(myGenes[ii, 3])
#   
#   if( is.na(outMat[cGene, cPool]) ){ 
#     outMat[cGene, cPool] = myGenes[ii, 1]  
#   } else {
#     outMat[cGene, cPool] = paste(outMat[cGene, cPool], myGenes[ii, 1], sep=";")  
#   }
#     
# }
# 
# finalOutMat = cbind(rownames(outMat), outMat)
# finalOutMat = rbind(c("gene", colnames(outMat)), finalOutMat)
# write.table(finalOutMat, file="gene_hp_pvals_0.05_15Sep17.txt", sep="\t"
#             ,col.names=F, row.names=F )
# 
# 
# ######################### print 0.05 - 0.1
# myGenes[,2] = as.character(myGenes[,2])
# myGenes = myGenes[myGenes[,1]>0.05,]
# geness = unique(myGenes[,2])
# 
# finalOut = c()
# 
# for (gene in geness){
#   ind = which(myGenes[,2]==gene)
#   poolNums = data.frame(table(myGenes[ind,3]))
#   outLine = c(gene, poolNums[,2])
#   finalOut = rbind(finalOut, outLine)
# }
# 
# finalOutt = rbind(c("gene", levels(poolNums[,1])), finalOut  )
# 
# outMat = matrix(NA, nrow=length(geness), ncol=length(unique(myGenes[,3])) )
# rownames(outMat) = geness
# colnames(outMat) = unique(myGenes[,3])
# 
# for (ii in 1:nrow(myGenes)){
#   cGene = as.character(myGenes[ii, 2])
#   cPool = as.character(myGenes[ii, 3])
#   
#   if( is.na(outMat[cGene, cPool]) ){ 
#     outMat[cGene, cPool] = myGenes[ii, 1]  
#   } else {
#     outMat[cGene, cPool] = paste(outMat[cGene, cPool], myGenes[ii, 1], sep=";")  
#   }
#   
# }
# 
# finalOutMat = cbind(rownames(outMat), outMat)
# finalOutMat = rbind(c("gene", colnames(outMat)), finalOutMat)
# write.table(finalOutMat, file="gene_hp_pvals_0.05to0.1_14AUG17.txt", sep="\t"
#             ,col.names=F, row.names=F )

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

