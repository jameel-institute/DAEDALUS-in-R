#DAEDALUS P2 run example

################################################################################

#David's notes:
#search for "xx"
#odes - output g for admissions

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
library(fastmatrix)

################################################################################
#Single run inputs - literally everything you need to change is in this block!

#Pick a country:
country <- "singapore"
#
#Copy the file path to your working directory here:
folder <- "C:/Users/dhaw/Documents/DAEDALUS-main"
#This directory should contain all files downloaded from GitHub
#
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
load("C:/Users/dhaw/Documents/DAEDALUS-main/singapore.Rdata")

lx        <- length(data$B)
data$int  <- 5
data$tvec <- c(-75, seq(from=as.numeric(data$Tres), to=365*3+1, length.out=as.numeric(data$int)))

#Create all objects that can be made pre-simulation:
listOut <- p2params(data, 'Covid Wildtype')
data <- listOut$data
dis <- listOut$dis
p2 <- listOut$p2

#Define economic configuration for all time periods:
xoptim <- c(rep(1, lx), rep(1, lx), rep(.5, lx), rep(1, lx), rep(1, lx))
data$hw    <- matrix(0, data$int, lx) #No testing in period 5 or after vaccine rollout (whatever comes first)
data$imand <- Inf #No social distancing in period 1, but continues until end
#data.imand is the index of the configuration with a stay at home (SAH) order

#Run the model:
listOut <- p2Run(data, dis, xoptim, p2)
data <- listOut$data
f <- listOut$f
g <- listOut$g

#Simple plot:
plot(f[,1],f[,3],ylim=c(0,10000),type="l")

colnames(f) <- c("Time", "Infections", "Hospital occ.", "Deaths (cumulative)", "Vacc. coverage",
                 "Transmission modifier",  "V1", "V2", "V3", "V4", "D1", "D2", "D3", "D4")

#Calculate costs:
listCost   <- p2Cost(data,dis,p2,g)
cost <- listCost$cost
ccost_t <- listCost$ccost_t
sec <- rep(0, 4)
sec[1]    <- sum(cost[c(3,6,7:10), ])
sec[2]    <- sum(cost[3, ])
sec[3]    <- sum(cost[6, ])
sec[4]    <- sum(cost[c(7:10), ])

#Plot stuff:
#p2Plot(data,f,p2,g,cost)