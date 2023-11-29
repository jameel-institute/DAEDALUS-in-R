#DAEDALUS P2 run example

#Install all of these packages if you haven't already:
#install.packages("R.matlab")
#install.packages("pracma")
#install.packages("deSolve")
#install.packages("ggplot2")

#Access packages:
library(R.matlab)
library(pracma)
library(deSolve)
library(ggplot2)

################################################################################
#Single run inputs - literally everything you need to change is in this block!

country <- "singapore"
folder <- "C:/Users/David/Downloads/DAEDALUS-main"
setwd(folder)

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

source("R_DAEDALUS_diseases.R")
source("R_DAEDALUS_functions.R")
load("C:/Users/David/Downloads/DAEDALUS-main/singapore.Rdata")

listOut <- p2params(data, 'Covid Wildtype')
data <- listOut$data
dis <- listOut$dis
p2 <- listOut$p2

simOut <- p2Run(data, dis, Xit, p2)
  
