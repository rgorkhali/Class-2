if(!requireNamespace("BiocManager", quietly=TRUE))
   install.packages("BiocManager")# Install the package Affy

if(!requireNamespace("BiocManager", quietly=TRUE))
   install.packages("BiocManager")
BiocManager::install("affy",version = "devel")

# OR

source("https://bioconductor.org/bioLite.R")
biocLite("affy")

source("https://bioconductor.org/biocLite.R")
biocLite("SpikeInSubset")

# Loading this package
library(SpikeInSubset)
library(affy)

# Load data called spikein95
data(spikein95)

# check the chips
image(spikein95)

# collect the gene ids from this dataset and put it into an object
ids<-geneNames(spikein95)
ids[1:10]

# perform Mas 5.0 micro arrary normalization
mas_log  Mas <- mas5(spikein95)

# Box plot for raw data,and boxplot for normalized data on log fold
boxplot(spikein95)
mas_log <- log2(exprs(Mas))
x11()ids<-geneNames(spikein95)
boxplot(mas_log)
boxplot(mas_log,col = 2:5)

# Density plot
density <- density(mas_log[,1])
plot(density,main="Mas normalization")

density2 <- density(mas_log[,2])
lines(density2,col="red")

density3<- density(mas_log[,3])
lines(density3,col="blue")


# Making MA plots
# M: difference in average log intensities
# A: average log intensities
#used 1,2,3 as treatment, and 4,5,6 as control

d <- rowMeans(mas_log[,1:3] - rowMeans (mas_log[,4:6]))
a <- rowMeans(mas_log)

# Plotting the data

plot(a,d, ylim = c(-5,5), main="Mas 5.0 normalized MA plot",
xlab = "A", ylab = "M", pch = ".")

abline(h = c(-1,1))
abline(h = c(-1.5,1.5), col="red")

#notmalization

rma.eset <- rma(spikein95)
rma.e <- exprs(rma.eset)

#checking the valued in object rma.e and mas_log
#rma.e function is in log scale 

head(rma.e)
head(mas_log)

#Plot two types of normalized plots
x11()
boxplot(mas_log,col = 2:5, main="Mas 5.0 Norm")
x11()
boxplot(rma.e,col = 2:5, main= "RMA Norm")

#Finding specific genes/probesets in MA plot
spikedn <- colnames(pData(spikein95))
spikedIndex <- match(spikedn, featureNames(rma.eset))
points(a[spikedIndex], d[spikedIndex], pch=19, col="red")

plot(a,d, ylim = c(-5,5), main="Mas 5.0 normalized MA plot",
xlab = "A", ylab = "M", pch = ".")

abline(h = c(-1,1))
abline(h = c(-1.5,1.5), col="red")


