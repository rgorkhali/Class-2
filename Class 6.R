#Volcano Plot
#load package and data 

library(SpikeInSubset)
data(spikein95)

#Run functions for normalization and, output into object,  pulled out expression values from the set 

rma.eset <- rma(spikein95)
rma.e <- exprs(rma.eset)

#took rowmeans of 3 columns and substracted with ...
d <- rowMeans(rma.e[,1:3]) - rowMeans(rma.e[,4:6])
a <- rowMeans(rma.e)
source("https://bioconductor.org/biocLite.R")
biocLite("genefilter")
# volcano plot
library("genefilter")
#pData(rma.eset) <- pData(mas5.eset)
tt <- rowMeans(rma.e)
lod <- tt

#cex= dot size; yaxt= "n-none"-> don't use default y-axis,we will provide it
plot(d, lod, cex = 0.25, main = "Volcano plot for MA", xlim = c(-2, 2), xlab = "M", ylab = "A", yaxt = "n")
#at-> sequence starting from 0 to 3 and increament by 1
#labels-> 
# 2-> left side for labeling axis
axis(2, at = seq(0, 3, by = 1), labels = 10^(-seq(0, 3, by = 1)))
axis(2, at = seq(0, 6, by = 1), labels = 10^(-seq(0, 6, by = 1)))
points(d[spikedIndex], lod[spikedIndex], pch = 19, col = "red")
abline(h = 2, v = c(-1, 1))

library(affy)

#set working directory and "" is for giving path
#for path got to folder, properties and copy the path and paste
#for windows forward slash
setwd("C:/Users/raksh/Dropbox/1 Rakshya Personal/bioinformatics Class Shrikant/estrogen")

#reading the txt file and putting that into the object
targetsFile <- "estrogen.txt"
#read the file with function read.AnnotatedDataFrame, header is first line, separator..there  is none, eg for csv file give ",", 
#even if you don't put object name its fine
#pData> creates samples as rows, variables as columns
pd <- read.AnnotatedDataFrame(targetsFile, header=TRUE, sep="", row.names=1)
#check data
pData(pd)


##making lineat fit with the two matrices

# pData dataframe and $ is playing with the column (asking if estrogen is absent or present
ER <-pData(pd)$estrogen
#factor-> ?? pulling unique time frame
Time <- factor(pData(pd)$time.h)
#designing matrix binary data to your dataset
design <-model.matrix(~ER+Time)
design
#matix with multiplication
design2 <-model.matrix(~ER*Time)
design2

# Read .cel files in the estrogen folder
#ReadAffy detects the probe ID....in our case human cell line  : #annotation=hgu95av2
#Affy chip!
raw <- ReadAffy(celfile.path="C:/Users/raksh/Dropbox/1 Rakshya Personal/bioinformatics Class Shrikant/estrogen", 
filenames=rownames(pData(pd)), phenoData = pd)
raw

x11()
#boxplot 
boxplot(raw)
x11()
boxplot(raw, col = "red")

#write on linear regression/linear fit??
#putty download for windows
#ssh gsu ID  rgorkhali1@snowball.cs.gsu.edu


