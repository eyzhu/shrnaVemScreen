library(edgeR)

###################################################
### code chunk number 14: zuber1
###################################################
# Read in the table of counts
dat = read.table("zuber_screen.txt", sep="\t", header=TRUE, as.is=TRUE)
dat[1,]

# Make DGE list containing hairpin counts
x4 = new("DGEList")
x4$counts = as.matrix(dat[,9:12])

# Remove hairpins with zero counts in all samples 
selnonzero = rowSums(x4$counts)!=0 
x4$counts = x4$counts[selnonzero,]

# Add sample annotation data
x4$samples = data.frame("SampleID"=colnames(x4$counts), 
                        "group"=as.factor(rep(c("Day0", "Day14"), times=2)), 
                        "lib.size"=colSums(x4$counts))
x4$samples$norm.factors = 1
x4$genes = dat[selnonzero,1:5]
rownames(x4$counts) = dat[selnonzero,1]
dim(x4)

# Make an MDS plot to visualise relationships between replicate samples
par(mfrow=c(1,3))
plotMDS(x4, labels=gsub("Reads_","",colnames(x4)), 
        col=c(1,2,1,2), main="Zuber: MDS Plot")
legend("topright", legend=c("Day 0", "Day 14"), col=1:2, pch=15)

# Assess differential representation between Day 14 and Day 0 samples
# using GLM in edgeR
# Set up design matrix for GLM
des = model.matrix(~x4$samples$group)
colnames(des)[2] = "Day14"
des

# Estimate dispersions
xglm = estimateDisp(x4, des)

# Plot BCVs versus abundance
plotBCV(xglm, main="Zuber: BCV Plot")

# Fit negative bionomial GLM
fit = glmFit(xglm, des)
# Carry out Likelihood ratio test
lrt = glmLRT(fit, 2)

# Show top ranked hairpins
topTags(lrt, n=15)

# Select hairpins with FDR < 0.0001 and logFC < -1 to highlight on plot
thresh = 0.0001
lfc = -1
top = topTags(lrt, n=Inf, sort.by="logFC")

sum(top$table[,9]<thresh)
sum(top$table[,9]<thresh & top$table[,6]<lfc)

topids = as.character(top$table[top$table$FDR<thresh & top$table$logFC<lfc,1])

# Make a plot of logFC versus logCPM
plotSmear(lrt, de.tags=topids, pch=20, cex=0.6, 
          main="Zuber: logFC vs logCPM")


###################################################
### code chunk number 15: zuber2
###################################################
# Carry out roast gene-set analysis
# Begin with hairpins targeting Brd4
genesymbols = x4$genes[,2]
brd4 = genesymbols=="Brd4"
set.seed(6012014)
roast(xglm, index=brd4, des, contrast=2, nrot=9999)

# Make a barcode plot for Brd4
par(mfrow=c(1,1))
barcodeplot(lrt$table$logFC, index=brd4, 
            main="Barcodeplot for Brd4 (Day14 versus Day0)",
            labels=c("Positive logFC", "Negative logFC"))

# Repeat analysis for all genes using mroast
genesymbollist = list()
for(i in unique(genesymbols))
  genesymbollist[[i]] = which(genesymbols==i)

roast.res = mroast(xglm, index=genesymbollist, des, contrast=2, nrot=9999)
roast.res[1,]

# Display results for top ranked genes
roast.res[1:20,1:6]

