#DAEDALUS P2 run example

#Install all of these packages if you haven't already:
#install.packages("R.matlab")
install.packages("pracma")
install.packages("deSolve")
install.packages("ggplot2")

#Access packages:
library(R.matlab)
library(pracma)
library(deSolve)
library(ggplot2)

################################################################################
#Single run inputs - literally everything oyu need toi change is in this block!

country <- "singapore"
folder <- "C:/Users/David/Downloads/DAEDALUS-main"

################################################################################
#Only use this block if you need to import any MATLAB objects

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
#Source/call DAEDALUS

setwd(folder)
source("R_DAEDALUS_diseases.R")
source("R_DAEDALUS_functions.R")
load("C:/Users/David/Downloads/DAEDALUS-main/singapore.Rdata")