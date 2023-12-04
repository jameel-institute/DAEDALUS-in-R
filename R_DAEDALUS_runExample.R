#DAEDALUS P2 run example

################################################################################

#David's notes:
#search for "xx"
#Comment out additional definitions at the start of functions, used for line-by-line running
#odes - output g for admissions

#NAs:
#p2$startp2 <- p2$startp1+1; p2$startp3 <- p2$startp2+1; p2$startp4 <- p2$startp3+1; p2$startp5 <- p2$startp4+1; p2$end <- p2$startp5+1

################################################################################

#Install all of these packages if you haven't already:
#install.packages("R.matlab")
#install.packages("pracma")
#install.packages("deSolve")
#install.packages("ggplot2")
#install.packages("fastmatrix")

#Access packages:
library(R.matlab)
library(pracma)
library(deSolve)
library(ggplot2)
#library(fastmatrix)

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

lx        <- length(data$B)
data$int  <- 5
data$tvec <- c(-75, seq(data$Tres, 365*3+1, length.out = data$int)) #xx Produces warning

listOut <- p2params(data, 'Covid Wildtype')
data <- listOut$data
dis <- listOut$dis
p2 <- listOut$p2

#xx Missing info
p2$startp2 <- data$tvec[2]#p2$startp1+1; 
p2$startp3 <- data$tvec[3]#p2$startp2+1; 
p2$startp4 <- data$tvec[4]#p2$startp3+1; 
p2$startp5 <- data$tvec[5]#p2$startp4+1; 
p2$end <- data$tvec[6]#p2$startp5+1

#xoptim     <- cbind(matrix(1, lx, 2), matrix(0.5, lx, 1), matrix(1, lx, 2))
xoptim <- c(rep(1, lx), rep(1, lx), rep(.5, lx), rep(1, lx), rep(1, lx))
data$hw    <- matrix(0, data$int, lx) #No testing in period 5 or after vaccine rollout (whatever comes first)
data$imand <- Inf #No social distancing in period 1, but continues until end
#data.imand is the index of the configuration with a stay at home (SAH) order


[data,f,g] = p2Run(data, dis, xoptim, p2)

#Still to translate:

#[cost,~]   = p2Cost(data,dis,p2,g);
#sec(1)     = sum(cost([3,6,7:10],:),'all');
#sec(2)     = sum(cost(3,:),'all');
#sec(3)     = sum(cost(6,:),'all');
#sec(4)     = sum(cost(7:10,:),'all');

#p2Plot(data,f,p2,g,cost);

#simOut <- p2Run(data, dis, Xit, p2)
  
