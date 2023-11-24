#DAEDALUS P2 run example
#install.packages("R.matlab")
install.packages("pracma")
library(R.matlab)
library(pracma)

#Single run inputs:
country <- "singapore"
folder <- "C:/Users/David/Downloads/DAEDALUS-main"

################################################################################

#Convert .mat file to R list
filename <- paste(country,".mat", sep = "", collapse = NULL)
f = readMat(filename, fixNames=TRUE)
data = f$data

fieldNames <- readMat("dataFieldnames.mat", fixNames=TRUE)
fieldNames <- unlist(fieldNames$dataFieldnames)
names(data) <- fieldNames

#Save as .Rdata
newFilename <- paste(country,".Rdata", sep = "", collapse = NULL)
save(data, file=newFilename)

################################################################################

setwd(folder)
source("R_DAEDALUS_functions.R")
load("C:/Users/David/Downloads/DAEDALUS-main/singapore.Rdata")