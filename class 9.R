#In class assignment: 09/24/2018, Name: "Rakshya"

### try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
biocLite("limma")

#set your working directory to desktop (otherwise says you don't have permit)
setwd("C:/Users/raksh/Desktop")

#load library GEOquery, limma and affy

library(GEOquery)
library(limma)
library(affy)

#Use following FTP link to download sample microarray data from GEO database
#use the function getGEO and assign the data to object named gse.

url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1000/matrix/GSE1000_series_matrix.txt.gz"
filenm <- "GSE1000_series_matrix.txt.gz"
if(!file.exists(filenm)) download.file(url, destfile=filenm)
gse <- getGEO(filename=filenm)

#File stored at: 
#C:\Users\raksh\AppData\Local\Temp\RtmpiqzI7f/GPL96.soft

#check what you have in the gse object
Object can be a class:  data frame, matrix, or a list or a vector....
#to find 
class(gse)

#information of object
head(gse)

#for expression level of object
head(exprs(gse))

#make 2 objects called treatment and control, 
and put exprs of first 5 columns of gse in treatment and next 5 column in control...

treatment <- exprs(gse[,c("GSM15785", "GSM15786",  "GSM15787",  "GSM15788",  "GSM15789",  "GSM15790",  "GSM15791")])

#OR better
treatment <- exprs(gse[,1:5]) ## Better way because we don't need the column labels/names

control <- exprs(gse[,6:10])

#Assignment
#take row means of object treatment and control
#hint: use "for loop" to calculate rowmeans in a dataframe for each row
rowMeans(treatment)
rowMeans(control)

#
treatment_means <- rowMeans(treatment)
control_means <- rowMeans(control)

#fold change of treatment means/control means and put in object called fold
fold <- treatment_means/control_means

#use if else loop, if the fold change is >2 put all genes into new object called up regulation
 and if fold change < -2 put the genes into new object called down regulation

up_regulation =list()
down_regulation =list()
if (fold > 2) {up_regulation <- fold} else if (fold < -2) {down_regulation <- fold}

##***did not work

#subsetting in r
up_regulation <- fold[which(fold > 2)]
down_regulation <- fold[which(fold < -2)]







