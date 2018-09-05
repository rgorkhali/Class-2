library(SpikeInSubset)
data(spikein95)
image(spikein95)
ids <- geneNames(spikein95)
ids[1:10]
mas5.eset <- mas5(spikein95)
mas5.e <- log2(exprs(mas5.eset))
boxplot(spikein95)
x11()
boxplot(mas5.e, col = 2:5)
density1 <- density(mas5.e[, 1])
plot(density1, main = "MAS5 expression measure distributions")
density2 <- density(mas5.e[, 2])
lines(density2, col = "red")
density3 <- density(mas5.e[, 3])
lines(density3, col = "blue")

# Making MA plots 
# M: difference in average log intensities
# A: average log intensities

d <- rowMeans(mas5.e[,1:3]) - rowMeans(mas5.e[,4:6])
a <- rowMeans(mas5.e)
plot(a, d, ylim = c(-5, 5), main = "MAS 5.0 MA plot", xlab = "A", ylab = "M", pch = ".")
abline(h = c(-1, 1))

# How to find specific genes in the plot
#double braket is function
# colnames= column names which are probe ID
#match-matches


spikedin <- colnames(pData(spikein95))
spikedIndex <- match(spikedin, featureNames(mas5.eset))
points(a[spikedIndex], d[spikedIndex], pch=19, col="red")

#google affymetrix ID for the gene of interest of specific species to get the ID

ESR1 <- "211233_x_at_1"
spikedIndex <-match(ESR1, featureNames(mas5.eset))
points(a[spikedIndex], d[spikedIndex], pch=19, col="blue")

#c is for strings of interested 
#pch= point size
genes_of_interest <- c("100_g_at", "1000_at", "1002_f_at")
spikedIndex <- match(genes_of_interest, featureNames(mas5.eset))
points(a[spikedIndex], d[spikedIndex], pch=19, col="green")

#Plot 5 genes of interest as orange
genes_of_interest <- c("669_s_at", "670_s_at", "671_at", "672_at", "673_at")
spikedIndex <- match(genes_of_interest, featureNames(mas5.eset))
points(a[spikedIndex], d[spikedIndex], pch=19, col="orange")

#weizmann geneAnnot -> put in probe set by using affy IDs and this gives official gene names

#comparing RMA normalization and Mas 5 normalization
#RMA is more stringent
rma.eset <- rma(spikein95)
rma.e <-exprs(rma.eset)

x11()
d <- rowMeans(rma.e[,1:3]) - rowMeans(rma.e[,4:6])
a <- rowMeans(rma.e)
plot(a, d, ylim = c(-5, 5), main = "MAS 5.0 MA plot", xlab = "A", ylab = "M", pch = ".")
abline(h = c(-1, 1))

#draw another threshold line
 abline(h = c(-0.5, 0.5), col="red")

abline(h = c(-0.2, 0.2), col="orange")


#MA VOLCANO PLOT 
# library gene filter in chrome
#bioconductor
#download package
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("genefilter")

library("genefilter")

tt <- rowMeans(rma.e)
lod <- tt
plot(d, lod, cex = 0.25, main="Volcano plot", xlim = c(-2,2), yaxt = "n")
axis(2, at = seq(0,3, by=1), labels = 10^(-seq(0,3, by = 1)))
points(d[spikedIndex], lod[spikedIndex], pch=19, col="red")
abline(h=2, v = c(-1,1)) 

