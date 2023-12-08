#DAEDALUS P2 run example
#install.packages("R.matlab")
install.packages("pracma")
library(R.matlab)
library(pracma)

#Single run inputs:
country <- "singapore"

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

R_p2Run <- function() {
  
}

################################################################################

R_p2params <- function(data, inp2) {
  lnn <- length(nn)
  
  #Population by age:
  nn <- data$Npop
  nn <- c(nn[1:16], sum(nn[17:lnn]))
  nntot <- c(nn[1], sum(nn[2:4]), sum(nn[5:13]), sum(nn[14:lnn]))
  ranges <- c(1, 3, 9, 4)
  nntot <- rep.int(nntot, times=ranges)
  nnprop <- nn/nntot
  subs <- c(1:4)
  subs <- rep.int(subs, times=ranges)
  
  #Population by sector:
  adInd <- 3
  lx <- length(data$obj)
  ntot <- length(data$NNs)#Vector - ML code specifies dim=1
  data$NNs[which(data$NNs==0)] <- 1
  data$alp <- 1
  
  #Contact matrix:
  listOut <- R_p2MakeDs(data, data$NNs, rep(1, lx), rep(0,lx));#xx
  Dout <- listOut$Dout
  data <- listOut$data
  
  ## INITIAL DISEASE PARAMETERS:
  
  if (inp2=='Influenza 2009'){
    dis <- p2Params_Flu2009#xx access from somewhere
  } elseif(inp2=='Influenza 1957'){
    dis <-  p2Params_Flu1957
  } elseif(inp2,'Influenza 1918'){
    dis <-  p2Params_Flu1918
  }elseif(inp2,'Covid Wildtype'){
    dis <-  p2Params_CovidWT
  }elseif(inp2,'Covid Omicron'){
    dis <-  p2Params_CovidOM;
  }elseif(inp2,'Covid Delta'){
    dis <-  p2Params_CovidDE 
  }elseif(inp2,'SARS'){
    dis <-  p2Params_SARS
  }else{
    stop('Unknown Disease!')
  }
  
  #Probabilities
  phgs    <- dis$ihr./dis$ps    #4*1
  pdgh    <- dis$ifr./dis$ihr
  phgs    <- accumarray(subs, phgs*nnprop)
  dis$ph  <- c(rep(phgs(adInd),lx), phgs)
  nnh     <- nn*dis$ihr
  nnhtot  <- c(nnh[1], sum(nnh[2:4]), sum(nnh[5:13]), sum(nnh[14:length(nnh)])) #1*17
  nnhtot  <- rep.int(nnhtot, times=ranges)
  nnhprop <- nnh/nnhtot
  pdgh    <- accumarray(subs, pdgh*nnhprop)
  dis$pd  <- c(rep(pdgh(adInd),lx), pdgh)
  
  #Durations (for each sector)
  dis$Ts = ((1-dis$ph)*dis$Tsr)   + (dis$ph*dis$Tsh)
  dis$Th = ((1-dis$pd)*dis$Threc) + (dis$pd*dis$Thd)
  
  #Rates
  dis$sig1 <- (1-dis$ps)/dis$Tlat  #0.088
  dis$sig2 <- dis$ps/dis$Tlat      # 0.1293
  dis$g1   <- 1/dis$Tay            # 0.4762
  dis$g2   <- (1-dis.ph)/dis$Ts    #[49 x 1]
  dis$g3   <- (1-dis$pd)/dis$Th    #[49 x 1]
  dis$h    <- dis$ph./dis$Ts
  dis$mu   <- dis$pd/dis$Th
  dis$nu   <- 1/dis$Ti
  
  #Transmission
  Deff  <- Dout*matrix(rep(data$NNs,ntot), ntot, ntot)/t(matrix(rep(data$NNs, ntot), ntot, ntot)) # [49 x 49]
  onesn <- rep(1, ntot)
  F     <- matrix(0, 3*ntot, 3*ntot); 
  F[1:ntot,ntot+1:dim(F)[2]] <- cbind(dis$red*Deff, Deff);   # [147 x 147] #xx dis$red scalar
  
  vvec <- c((dis$sig1+dis$sig2)*onesn, dis$g1*onesn, (dis$g2+dis$h)*onesn) #g2 and h are vectors
  V    <- diag(vvec)
  V[ntot+1:2*ntot, 1:ntot]   <- diag(-dis$sig1*onesn)
  V[2*ntot+1:3*ntot, 1:ntot] <- diag(-dis$sig2*onesn)
  
  GD <- F%*%inv(V) # [147 x 147]
  ev <- eigen(GD)#largest in magnitude (+/-)  % 7.29
  d <- ev$values
  R0a <- max(d) # 7.29
  dis$beta <- dis$R0/R0a#beta scales actual R0 to fitted R0
  
  #Vaccination
  dis$hrv1 <- 1/28                       #time to develop v-acquired immunity
  dis$scv1 <- 0.60                       #infection-blocking efficacy
  dis$heff <- 0.87                       #severe-disease-blocking efficacy
  dis$hv1  <- 1-((1-dis.heff)/(1-dis.scv1))  # 0.6750 probability of reduction in infection/overall VE combining severe and symptomatic infection
  dis$trv1 <- 0.52                       #transmission-blocking efficacy %0.52
  dis$nuv1 <- 1/365                      #duration of v-acquired immunity 0.0027
  
  dis$Ts_v1 <- ((1-(1-dis$hv1)*dis$ph)*dis$Tsr)  +((1-dis$hv1)*dis$ph*dis$Tsh) # 4 days  [49x1]
  dis$g2_v1 <- (1-(1-dis$hv1)*dis$ph)/dis$Ts_v1 # [49 x 1]
  dis$h_v1  <- (1-dis$hv1)*dis$ph/dis$Ts_v1 # hospitalisation rate for vaccinated individuals [49x1]
  
  ## PREPAREDNESS PARAMETERS:
    
  p2 = struct;
  
  p2$t_tit <- data$t_tit                      #Test-Isolate-Trace Time
  p2$trate <- data$trate                      #Test-Isolate-Trace Rate
  p2$sdl   <- data$sdl                        #Social Distancing Asymptote
  p2$sdb   <- data$sdb                        #Social Distancing Steepness
  p2$Hmax  <- data$Hmax*sum(data$Npop)/10^5   #Hospital Capacity
  t_vax    <- data$t_vax                      #Vaccine Administration Time
  arate    <- data$arate*sum(data$Npop/10^5)  #Vaccine Administration Rate
  puptake  <- data$puptake                    #Vaccine Uptake
  
  #Response Time
  J                                   <- matrix(0, 7*ntot, 7*ntot)
  J[1:ntot, 2*ntot+1:3*ntot]          <- -dis$beta*dis$red*Dout
  J[1:ntot, 3*ntot+1:4*ntot]          <- -dis$beta*Dout
  J[1:ntot, 5*ntot+1:6*ntot]          <- diag(onesn*dis$nu)
  J[ntot+1:2*ntot, 1*ntot+1:2*ntot]   <- diag(onesn*(-dis$sig1-dis$sig2))
  J[ntot+1:2*ntot, 2*ntot+1:3*ntot]   <- dis$beta*dis$red*Dout
  J[ntot+1:2*ntot, 3*ntot+1:4*ntot]   <- dis$beta*Dout
  J[2*ntot+1:3*ntot, 1*ntot+1:2*ntot] <- diag(onesn*dis$sig1)
  J[2*ntot+1:3*ntot, 2*ntot+1:3*ntot] <- diag(onesn*-dis$g1)
  J[3*ntot+1:4*ntot, 1*ntot+1:2*ntot] <- diag(onesn*dis$sig2)
  J[3*ntot+1:4*ntot, 3*ntot+1:4*ntot] <- diag(onesn*(-dis$g2-dis$h))
  J[4*ntot+1:5*ntot, 3*ntot+1:4*ntot] <- diag(onesn*dis$h)
  J[4*ntot+1:5*ntot, 4*ntot+1:5*ntot] <- diag(onesn*(-dis$g3-dis$mu))
  J[5*ntot+1:6*ntot, 2*ntot+1:3*ntot] <- diag(onesn*dis$g1)
  J[5*ntot+1:6*ntot, 3*ntot+1:4*ntot] <- diag(onesn*dis$g2)
  J[5*ntot+1:6*ntot, 4*ntot+1:5*ntot] <- diag(onesn*dis$g3)
  J[5*ntot+1:6*ntot, 5*ntot+1:6*ntot] <- diag(onesn*-dis$nu)
  J[6*ntot+1:7*ntot, 4*ntot+1:5*ntot] <- diag(onesn*dis$mu)
  
  ev <- eigen(J)
  r       <- max(real(ev$values))
  Td      <- log(2)/r
  if ~exists(data$Td_CWT){
    data$Td_CWT <- Td
  }
  
  #Test-Isolate-Trace
  p2$dur    <- 1
  p2$qg1    <- 1/(dis$Tay-p2$dur)
  p2$qg2    <- (1-dis$ph)/(dis$Ts-p2$dur)
  p2$qg2_v1 <- (1-(1-dis$hv1)*dis$ph)/(dis$Ts_v1-p2$dur)
  p2$qh     <- dis$ph/(dis$Ts-p2$dur)
  p2$qh_v1  <- (1-dis$hv1)*dis$ph/(dis$Ts_v1-p2$dur)
  
  #Hospital Capacity
  p2$thl   <- max(1, 0.25*p2$Hmax)#lower threshold can't be less than 1 occupant
  p2$Hmax  <- max(4*p2$thl, p2$Hmax)
  p2$SHmax <- 2*p2$Hmax
  
  #Vaccine Uptake
  Npop    <- data$Npop
  NNage   <- c(Npop[1],sum(Npop[2:4]),sum(Npop[5:13]),sum(Npop[14:length(Npop)]))
  puptake <- min(0.99*(1-NNage[1]/sum(NNage)),puptake)#population uptake cannot be greater than full coverage in non-pre-school age groups
  up3fun  <- func(u3) puptake*sum(NNage) - u3*(NNage[2]/2 + NNage[3]) - min(1.5*u3, 1)*NNage[4]
  if (up3fun(0)*up3fun(1)<=0){
    u3  <- uniroot(up3fun, c(0, 1))#xx check same
  } else{
    u3  <- fminbnd(up3fun, 0, 1)
  }
  u4      <- min(1.5*u3, 1)
  u1      <- 0
  up2fun  <- func(u2) u2*NNage[2] + u3*NNage[3] + u4*NNage[4] - puptake*sum(NNage)
  u2      <- uniroot(up2fun, c(0, 1));#xx check same
  uptake  <- c(u1,u2,u3,u4);
  
  #Vaccine Administration Rate
  t_ages     = min((uptake.*NNage)/arate,Inf)#arate may be 0
  if (inp2=='Influenza 1918'){
    t_ages     <- c(t_ages[3], t_ages[4], t_ages[2], t_ages[1])
    p2$aratep1 <- c(0, 0, arate, 0)#Period 1 - working-age#to be split across all economic sectors in heSimCovid19vax.m
    p2$aratep2 <- c(0, 0, 0, arate)#Period 2 - retired-age
    p2$aratep3 <- c(0, arate, 0, 0)#Period 3 - school-age
    p2$aratep4 <- c(0, 0, 0, 0)    #Period 4 - pre-school-age
  } else{
    t_ages     <- c(t_ages[4], t_ages[3], t_ages[2], t_ages[1])
    p2$aratep1 <- c(0, 0, 0, arate)#Period 1 - retired-age
    p2$aratep2 <- c(0, 0, arate, 0)#Period 2 - working-age%to be split across all economic sectors in heSimCovid19vax.m
    p2$aratep3 <- c(0, arate, 0, 0)#Period 3 - school-age
    p2$aratep4 <- c(0, 0, 0, 0)    #Period 4 - pre-school-age
  }    
  tpoints    <- cumsum(c(t_vax, t_ages))
  p2$startp1 <- tpoints[1]
  p2$startp2 <- tpoints[2]
  p2$startp3 <- tpoints[3]
  p2$startp4 <- tpoints[4]
  p2$end     <- tpoints[5]#End of Rollout
  
  ## COST PARAMETERS:
    
    
  n    <- c(data$Npop[1:16], sum(data$Npop[17:end]))#length is 17 to match ifr
  la   <- c(data$la[1:16], dot(data$la[17:end], c(data$Npop(17), sum(data$Npop[18:length(Npop)])))/sum(data$Npop[17:end]))
  napd <- na*dis$ifr
  lg   <- c(dot(la[1], napd[1])/sum(napd[1]), 
            dot(la(2:4), napd(2:4))/sum(napd(2:4)), 
            dot(la[5:13], napd[5:13])/sum(napd[5:13]), 
            dot(la[14:length(la)], napd[14:length(napd)])/sum(napd[4:length(napd)]))

for (k in 1:length(lg)){
    lgh[k] <- sum(1/((1+0.03)^(1:lg(k))))
}  
data$lgh   <- c(rep(lgh[3],45), lgh)
  
}  

################################################################################
#Prepare objects for simulation

R_p2MakeDs <- function(data, NN, x, hw) {
  C <- data$CM
  Npop <- data$Npop
  Npop[16] <- sum(Npop[16:length(Npop)])
  Npop <- Npop[1:16]
  C <- cbind(C[,1],rowSums(C[,2:4]),rowSums(C[,5:13]),rowSums(C[,14:16]))#Column sums
  C <- rbind(C[1,],
             Npop[2:4]%*%C[2:4,]/sum(Npop[2:4]),
             Npop[5:13]%*%C[5:13,]/sum(Npop[5:13]),
             Npop[2:4]%*%C[2:4,]/sum(Npop[2:4]))
  Cav <- (c(Npop[1],sum(Npop[2:4]),sum(Npop[5:13]),sum(Npop[14:16]))/sum(Npop))%*%rowSums(C)#weighted average of rows
  C <- data[["comm"]][1,1]*(C/Cav[1,1])#xx better way than as.numeric
  
  ##
  
  adInd <- 3#Adult index
  CworkRow <- C[adInd,]
  lx <- nrow(x)#Number of sectors xx rows are sectors
  ln <- length(NN)
  
  NNrep <- t(matrix(rep(NN/sum(NN),ln),ln,ln))#Total population proportion matrix
  NNrel <- NN[c(1:lx,lx+addInd)]/sum(NN[c(1:lx,lx+addInd)])#Adult proportion proportion vector
  NNrea <- t(matrix(rep(NN[1:lx]/sum(NN[1:lx]),lx),lx,lx))#Workforce population proportion matrix
  
  #Make A:
  matA <- matrix(0,ln,ln)
  matA[lx+1:ln,lx+1:ln] <- C
  matA[1:lx,lx+1:ln] <- t(matrix(rep(CworkRow,lx),ln-lx,lx))
  matA[,c(1:lx,lx+adInd)] <- matrix(rep(matA[,lx+adInd],lx+1),lx+1,lx+1)*t(matrix(rep(NNrel,ln),ln,ln))

  ##
  
  data$EdInd <- 41#Education sector index
  data$hospInd <- c(32,43,44)#Hospitality sector indices
  w <- x^(1/data$alpha)
  
  if (lx==45){
    #Education:
    matA[lx+1,lx+1] <- matA[lx+1,lx+1] + w[data$EdInd]^2*data$schoolA1#Mixing within age groups only
    matA[lx+2,lx+2] <- matA[lx+2,lx+2] + w[data$EdInd]^2*data$schoolA2
    
    #Hospitality:
    psub <- data$NNs[data$HospInd]
    psub <- sum(psub*x[data$HospInd])/sum(psub)#Constant from 0-1, weighted measure of how much sectors are open
    matA[c(1:lx,lx+adInd),] <- matA[c(1:lx,lx+adInd),] + psub^2*data$hospA3*NNrep[c(1:lx,lx+adInd),]
    matA[lx+2,] <- matA[lx+2,] + psub^2*data$hospA2*NNrep[lx+2,]
    matA[ln,] <- matA[ln,] + psub^2*data$hospA4*NNrep[ln,]
  } else{
    stop("Unknown economic configuration!")
  }

  #Transport:
  matA[1:lx,1:lx] <- matA[1:lx,1:lx] + t(matrix(rep(w,lx),lx,lx))*data$traveA3[1]*NNrea*t(matrix(rep(1-hw,lx),lx,lx))#Home working has compound effect
  
  #Worker-worker and community-worker matrices:
  
  #Make B and C
  valB <- data[["B"]]
  valB <- valB*(1-hw)*(1-hw)
  valC <- data$C
  valC <- valC*(1-hw)
  valB[lx+1:ln] <- rep(0,ln-lx)
  valC[lx+1:ln] <- rep(0,ln-lx)
  x[lx+1:ln] <- rep(0,ln-lx)
  w[lx+1:ln] <- rep(0,ln-lx)
  matB <- diag(w*valB)
  matC <- t(matrix(rep(x*valC,ln),ln,ln))*NNrep
  
  if (!exists(data,"wnorm")){
    data$wnorm <- dot(rowSums(matB + matC),NN)/sum(NN[c(1:lx,lx+adInd)])
  }
  
  D <- matA + (data$workp/data$wnorm)%*%(matB + matC)
  
  listOut <- list("Dout" = D, "data" = data)
  return(listOut)
}
