# Translation of p2Run from MATLAB to R


library(deSolve)
library(ggplot2)

p2Run <- function(data, dis, Xit, p2) {
  
  adInd <- 3
  lx <- length(data$obj) # 45 sectors
  ln <- length(data$NNs) # 49 (45 sectors + 4 age groups)
  int <- data$int # 5 time periods
  NNbar <- data$NNs # Total workforce in each sector if economy is fully open
  XitMat <- reshape(Xit, lx, int)
  WitMat <- XitMat ^ (1 / data$alp)
  WitMat[data$EdInd, ] <- XitMat[data$EdInd, ]
  NNvec <- rep(NNbar[1:lx], times = int) * WitMat
  NNworkSum <- sum(NNvec, 1)
  NNvec[lx + 1:ln, ] <- rep(NNbar[lx + 1:ln], times = int)
  NNvec[lx + adInd, ] <- sum(NNbar[1:lx, lx + adInd]) - NNworkSum
  data$NNvec <- NNvec
  
  Dvec <- array(0, dim = c(ln, ln, int))
  for (i in 1:int) {
    Dtemp <- p2MakeDs(data, NNvec[, i], XitMat[, i], data$hw[i, ])
    Dvec[,, i] <- Dtemp
  }
  data$Dvec <- Dvec
  
  data <- p2SimVax(data, NNvec, Dvec, dis, NNvec[, 1], WitMat, p2)
  
  return(data)
}



p2SimVax <- function(data, NNvec, Dvec, dis, S0, WitMat, p2) {
  
  ntot <- length(data$NNs)
  adInd <- 3
  lx <- ntot - 4
  NNbar <- NNvec[, 1]
  sumWorkingAge <- sum(NNbar[1:lx, lx + 3])
  nc <- 20
  zn <- rep(0, ntot)
  y0 <- c(S0, rep(zn, 6), NNbar - S0, rep(zn, nc - 9), S0)
  Tout <- data$tvec[1]
  Iout <- zn
  Isaout <- zn
  Issout <- zn
  Insout <- zn
  Hout <- zn
  Dout <- zn
  Wout <- matrix(0, nrow = 0, ncol = lx)
  hwout <- matrix(0, nrow = 0, ncol = lx)
  poutout <- 0
  betamodout <- 1
  Vout <- zn
  rout <- 0

# Loop
  for (i in 1:data$int) {
    
    t0 <- data$tvec[i]
    tend <- data$tvec[i + 1]
    Wit <- WitMat[, i]
    NNfeed <- NNvec[, i]
    NNfeed[NNfeed == 0] <- 1
    D <- Dvec[,, i]
    
  # Vaccination roll out by sector
    
    NNnext <- NNvec[, i]  # total population in a given period i
    NNnext[lx + c(1, 2)] <- 1
    NNnext[1:lx, lx + 3] <- NNnext[1:lx, lx + 3] / sumWorkingAge
    NNnext[ntot] <- 1 # Retired age population
    
    p2$ratep1 <- NNnext * c(rep(p2$aratep1[3], lx, 1), p2$aratep1)
    p2$ratep2 <- NNnext * c(rep(p2$aratep2[3], lx, 1), p2$aratep2)
    p2$ratep3 <- NNnext * c(rep(p2$aratep3[3], lx, 1), p2$aratep3)
    p2$ratep4 <- NNnext * c(rep(p2$aratep4[3], lx, 1), p2$aratep4)
    
    
  # Solve the ODE system
    res <- ode(y = y0,
                        t = c(t0, tend),
                        func = p2Model,
                        parms = p2)
    
  # Extract the results
    Tout <- c(Tout, res$t[2:end])
    Iout <- c(Iout, res$y[2:end, 1])
    Isaout <- c(Isaout, res$y[2:end, 2])
    Issout <- c(Issout, res$y[2:end, 3])
    Insout <- c(Insout, res$y[2:end, 4])
    Hout <- c(Hout, res$y[2:end, 5])
    Dout <- c(Dout, res$y[2:end, 6])
  
  # Calculate the workforce
    W <- Wit * matrix(1, ncol = lx)
    Wout <- rbind(Wout, W[1:(nrow(W) - 1), ])
    
  # worker hours
    hw <- data$hw[i, ] * matrix(1, ncol = lx)
    hwout <- rbind(hwout, hw[1:(nrow(W) - 1), ]) 
    
  # vaccination coverage
    poutout <- c(poutout, res$y[2:end, 19])
    
  # transmission modifier
    betamodout <- c(betamodout, res$y[2:end, 20])
    
  # vaccination uptake
    Vout <- c(Vout, res$y[2:end, 17])
    
    
  # Update the workforce and worker hours for the next time period
    
    if (i < data$int) {
      
      y0 <- reshape(y0, ntot, nc)
      
      # Calculate the number of extra workers needed next period
      Xh2w <- NNvec[1:lx, i + 1] - NNvec[1:lx, i]
  
      
      # Set the negative elements of Xh2w and Xw2h to zero
      Xh2w[Xh2w < 0] <- 0
      Xw2h[Xw2h < 0] <- 0
      
      # Calculate the proportion of extra workers needed next period by sector
      Xh2w <- Xh2w / NNvec[lx + adInd, i]
      
      # Set the elements of Xh2w to zero for sectors that have no workers
      Xh2w[NNvec[lx + adInd, i] == 0] <- 0
      
      # Calculate the proportion of workers not needed next period by sector
      Xw2h <- Xw2h / NNvec[1:lx, i]
      
      # Set the elements of Xw2h to zero for sectors that have no workers
      Xw2h[NNvec[1:lx, i] == 0] <- 0
      
      # Calculate the number of people to be put at home
      y0w2h <- y0[1:lx, ] * matrix(Xw2h, nrow = 1, ncol = nc)
      
      # Add the people to be put at home to the non-working population and subtract them from the workforce
      y0w2h <- c(-y0w2h, sum(y0w2h, 1))
      
      # Calculate the number of people to be put at work
      y0h2w <- kron(y0[lx + adInd, ], Xh2w)
      
      # Add the people to be put at work to the workforce and subtract them from the non-working population
      y0h2w <- c(y0h2w, -sum(y0h2w, 1))
      
      # Update the y0 vector
      y0[1:lx, ] <- y0[1:lx, ] + y0w2h
      y0[lx + adInd, ] <- y0[lx + adInd, ] + y0h2w

      y0 = reshape(y0,ntot*nc,1)
      
    }
  
}

  # Outputs
  
  Wout <- rbind(Wout, Wout[nrow(Wout), ])
  
  hwout <- rbind(hwout, hwout[nrow(hwout), ])
  
  # Create a new matrix `g` that stores results of the ODE 
  g <- cbind(Tout, Wout, hwout, Isaout, Issout, Insout, Hout, Dout, Vout, betamodout)
  
  f <- cbind(Tout,
             sum(Iout, 2),
             sum(Hout, 2),
             sum(Dout, 2),
             poutout,
             betamodout,
             sum(Vout[,(lx + 1)], 2),
             sum(Vout[,(lx + 2)], 2),
             sum(Vout[, c(1:lx, lx + 3)], 2),
             sum(Vout[,(lx + 4)], 2),
             sum(Dout[,(lx + 1)], 2),
             sum(Dout[,(lx + 2)], 2),
             sum(Dout[, c(1:lx, lx + 3)], 2),
             sum(Dout[,(lx + 4)], 2))
  
  return(list(f = f, g = g))
  
}


# ODEs

integr8 <- function(data, NN0, D, i, t0, tend, dis, y0, p2) {
  ntot <- length(data$NNs)
  params <- list(data = data, NN0 = NN0, D = D, i = i, dis = dis, p2 = p2, b0 = 2.197, b1 = 0.1838, b2 = -1.024)
  
  # Define the ODE function
    odes <- function(y, t, params) {
    ode(y, t, params)
    }
  
  # Solve the ODEs
  yout <- ode(y = y0, times = seq(t0, tend, by = 0.1), func = odes, parms = params)
  
  tout <- yout[, "time"]
  y0new <- tail(yout, n = 1)[, -1] 
  
  # Output
  Iclass <- rowSums(yout[, (2*ntot + 1):(3*ntot)] + yout[, (3*ntot + 1):(4*ntot)] + ...)
  Isaclass = rowSums(yout[, (3*ntot+ 1): (4*ntot)] + yout[, (12*ntot + 1): (13*ntot) ] )
  Issclass = rowSums( yout[, (5*ntot+1): (6*ntot)] + yout[, (14*ntot+1) : (15*ntot) ] )
  Insclass = rowSums ( yout[, (4*ntot+1): (5*ntot)] + yout[, (13*ntot+1) : (14*ntot) ] ) 
  Hclass   = rowSums (yout [,  (6*ntot+1): (7*ntot)] + yout[, (15*ntot+1) :(16*ntot) ] )
  Dclass   = rowSums (yout [, (17*ntot+1) : (18*ntot) ] )
  Vclass   =  rowSums (yout [, (18*ntot+1) : (19*ntot) ] )
  
  
# Time - dependent parameters
  
  occ <- pmax(1, rowSums(Hclass))
  Hmax <- p2$Hmax
  SHmax <- p2$SHmax
  th0 <- pmax(1, 1 + 1.87 * ((occ - Hmax) / (SHmax - Hmax)))
  
  pd <- pmin(th0 * dis$pd, 1)
  Th <- (1 - pd) * dis$Threc + pd * dis$Thd
  mu <- pd / Th
  ddk <- 10^5 * rowSums(mu * Hclass) / sum(NN0)
  
  sd_fun <- function(l, b, x) {
    (l - b) + (1 - l + b) * (1 + ((l - 1) / (1 - l + b)))^(x / 10)
  }
  
  # betamod 
  if (i == 1) {
    betamod <- rep(1, length(occ))
  } else if (i %in% data$imand) {
    betamod <- pmin(pmax(p2$sdl, sd_fun(p2$sdl, p2$sdb, ddk)), pmax(p2$sdl, sd_fun(p2$sdl, p2$sdb, 2)))
  } else {
    betamod <- pmax(p2$sdl, sd_fun(p2$sdl, p2$sdb, ddk))
  }

  S <- yout[, 1:ntot]
  Sn <- yout[, (19*ntot + 1):(20*ntot)]
  Ins <- yout[, (4*ntot + 1):(5*ntot)]
  Iss <- yout[,(5*ntot+1): (6*ntot)]
  Insv1 <- yout[,(13*ntot+1):(14*ntot)]
  Issv1 <- yout[,(14*ntot+1):(15*ntot)]
  H     <- yout[,(6*ntot+1):(7*ntot)]
  Hv1   <- yout[, (15*ntot+1):(16*ntot)]
  
  
  amp <- (Sn + (1 - dis$heff) * (S - Sn)) / S
  ph <- amp * dis$ph
  Ts <- (1 - ph) * dis$Tsr + ph * dis$Tsh
  g3 <- (1 - pd) / Th
  h <- ph / Ts
  h_v1 <- dis$h_v1
  dur <- p2$dur
  qh <- ph / (Ts - dur)
  qh_v1 <- p2$qh_v1  
  
  # Hospitalisation
  Hdot <- (h * Ins) + (qh * Iss) - ((g3 + mu) * Hclass)
  Hv1dot <- h_v1 * Insv1 + qh_v1 * Issv1 - (g3 + mu) * Hv1
  occdot <- rowSums(Hdot + Hv1dot)
  r <- occdot/occ
  
  Ina <- yout[, (2*ntot+1):(3*ntot)]
  Isa <- yout[, (3*ntot+1):(4*ntot)]
  Inav1 <- yout[, (11*ntot+1):(12*ntot)]
  Isav1 <- yout[, (12*ntot+1):(13*ntot)]
  
  Ip  <-  10^5*sum(Ina+Ins+Isa+Iss+Inav1+Insv1+Isav1+Issv1,2)/sum(NN0)
  trate <- p2$trate
  
  pout <- ifelse(i != 5, 
                 (Ip < trate) * (1 / (1 + exp(b0 + b1 * Ip + b2 * log10(trate)))) / p2$dur + 
                   (Ip >= trate) * pmin(1 / (1 + exp(b0 + b1 * Ip + b2 * log10(trate))), trate / 10^5) / p2$dur * 
                   (tout > p2$t_tit) * (tout < p2$end), 
                 rep(0, length(tout)))
  
  results <- list(
    tout = tout,
    Iclass = Iclass,
    Isaclass = Isaclass,
    Issclass = Issclass,
    Insclass = Insclass,
    Hclass = Hclass,
    Dclass = Dclass,
    Vclass = Vclass,
    y0new = y0new,
    betamod = betamod,
    pout = pout
  )
  
  return(results)
}



odes <- function(y, t, params) {
  
  params <- list(data = data, NN0 = NN0, D = D, i = i, dis = dis, p2 = p2, ntot = length(NNs), 
                 b0 = 2.197, b1 = 0.1838, b2 = -1.024, phi = 1, betamod = 1 )
  
  S <- y[1:ntot]
  E <- y[(ntot + 1):(2 * ntot)]
  Ina <- y[(2 * ntot + 1):(3 * ntot)]
  Isa <- y[(3*ntot+1) : (4*ntot) ]
  Ins <- y[(4*ntot+1) : (5*ntot)]
  Iss <- y[(5*ntot+1): (6*ntot)]
  H   <- y[(6*ntot+1): (7*ntot)]
  R  <-  y[(7*ntot+1): (8*ntot)]
  Shv1 <- y[(8*ntot+1): (9*ntot)]
  Sv1 <- y[(9*ntot+1):(10*ntot)]
  Ev1 <- y[(10*ntot+1):(11*ntot)]
  Inav1 <- y[(11*ntot+1): (12*ntot)]
  Isav1 <- y[(12*ntot+1):(13*ntot)]
  Insv1 <- y[(13*ntot+1): (14*ntot)]
  Issv1 <- y[(14*ntot+1): (15*ntot)]
  Hv1   <- y[(15*ntot+1): (16*ntot)]
  Rv1   <- y[(16*ntot+1): (17*ntot)]
  DE    <- y[(17*ntot+1): (18*ntot)]
  V     <- y[(18*ntot+1):(19*ntot)]
  Sn    <- y[(19*ntot+1):(20*ntot)]
  
  # Hospital Occupancy
    
  occ   <- max(1,sum(H+Hv1))
  Hmax  <- p2$Hmax
  SHmax <- p2$SHmax
  
  # Time-dependent disease parameters 
    
  # Amplitudes
  
  amp <- (Sn+(1-dis$heff)*(S-Sn))/S
  th0 <- max(1, (1+1.87)*((occ-Hmax)/(SHmax-Hmax)))
  
  # Probabilities
  ph <- amp*(dis$ph)
  pd <- min(th0*dis$pd,1)
  
  # Calculations
  Ts <- ((1-ph)* dis$Tsr) + (ph*(dis$Tsh))
  Th <- ((1-pd)* dis$Threc)+ (pd*(dis$Thd))
  
  sig1 <- dis$sig1
  sig2 <- dis$sig2
  g1   <- dis$g1
  g2   <- (1-ph)/Ts
  g3   <- (1-pd)/Th
  h    <- ph/Ts
  mu   <- pd/Th
  nu   <- dis$nu
  
  # Transmission
  red  <- dis$red
  beta <- dis$beta
  
  # Vaccination
  hrv1  <- dis$hrv1    
  scv1  <- dis$scv1 
  g2_v1 <- dis$g2_v1
  h_v1  <- dis$h_v1
  trv1  <- dis$trv1
  nuv1  <- dis$nuv1
  
  # Preparedness
  
  dur   <- p2$dur
  qg1   <- p2$qg1
  qg2   <- (1-ph)/(Ts-dur)
  qg2_v1 <- p2$qg2_v1
  qh     <- ph/(Ts-dur)
  qh_v1  <- p2$qh_v1
  
  startp1 <- p2$startp1
  startp2 <- p2$startp2
  startp3 <- p2$startp3
  startp4 <- p2$startp4
  pend    <- p2$end
  
  ratep1 <- p2$ratep1
  ratep2 <- p2.ratep2
  ratep3 <- p2.ratep3
  ratep4 <- p2.ratep4
  
  # Force of Infection; phi = 1 in params; betamod = 1 in params
  
  ddk    = 10^5*sum(mu*(H+Hv1))/sum(NN0)
  
  sd_fun <- function(l, b, x) {
    (l - b) + (1 - l + b) * (1 + ((l - 1) / (1 - l + b)))^(x / 10)
  }
  

  I <- (red * Ina + Ins) + (1 - trv1) * (red * Inav1 + Insv1)
  foi <- phi * beta * betamod * (D * (I / NN0))
  
  seedvec <- rep(10^-15 * sum(data$Npop), ntot)
  seed <- phi * beta * betamod * (D * (seedvec / NN0))
  
  # Self-Isolation
  
  if (t < p2$t_tit) {
    p3 <- rep(0, ntot)
    p4 <- rep(0, ntot)
  } else if (t < p2$end && i != 5) {
  
    Ip  <- 10^5* sum(Ina+Ins+Isa+Iss+Inav1+Insv1+Isav1+Issv1)/sum(NN0)
    trate <- p2$trate
    #b0    = 2.197
    #b1    = 0.1838
    #b2    = -1.024
    p3  <- (Ip<trate) *   (1/(1+exp(b0+ (b1*Ip) + (b2*log10(trate)))))/dur + ...
    (Ip>=trate)*min(1 / (1+exp(b0+ (b1*Ip) + (b2*log10(trate)))),trate/10^5)/dur
    p4    <- p3
  
  } else {
    p3 <- rep(0, ntot)
    p4 <- rep(0, ntot)
  }
  
  
# Vaccination
  
  nonVax = NN0 - V
  
  if (t >= pend) {
    
    v1rates <- rep(0, ntot)
    v1rater <- rep(0, ntot)
    Vdot <- rep(0, ntot)
    
  } else if (t>= startp4) {
    
    v1rates <- ratep4*S/nonVax
    v1rater <- ratep4*R/nonVax
    Vdot    <-  ratep4
    
  } else if (t >= startp3) {
    
    v1rates <- ratep3*S/nonVax
    v1rater <- ratep3*R/nonVax
    Vdot   <- ratep3
    
  } else if (t>= startp2) {
  
    v1rates <- ratep2*S/nonVax
    v1rater <- ratep2*R/nonVax
    Vdot    <-  ratep2
    
  } else if (t>=startp1) {
    
   v1rates <- ratep1*S/nonVax
   v1rater <- ratep1*R/nonVax
   Vdot    <- ratep1
   
  }  else {
   
    v1rates <- rep(0, ntot) 
    v1rater <- rep(0, ntot) 
    Vdot    <- rep(0, ntot) 
  
  }


# Equations
  
  Sndot       <-     -Sn*(foi+seed) - ( (v1rates*Sn)/S )
  Sdot        <-     -S*(foi+seed)  + (nu*R)  - v1rates +  ( nuv1*Sv1 )
  Shv1dot     <-      v1rates  -  ( hrv1*Shv1 )  - ( Shv1*foi )
  Sv1dot      <-      ( hrv1*Shv1 ) - (Sv1*(1-scv1)*foi)  - (nuv1*Sv1)  
  Edot        <-      S*(foi+seed)  + (Shv1*foi)  - ( (sig1+sig2)*E )
  Ev1dot      <-      ( Sv1 * (1-scv1)* foi )  - ( (sig1+sig2)*Ev1 )
  Inadot      <-      ( sig1 * E )  - (g1 * Ina )   - (p3*Ina)
  Insdot      <-      ( sig2 * E )  - ( (g2+h)*Ins )  - ( p4*Ins)
  Inav1dot    <-      (sig1 * Ev1 )  - (g1*Inav1 )  - (p3*Inav1)
  Insv1dot    <-      ( sig2 * Ev1 )  - ( g2_v1 + h_v1)*Insv1 - (p4*Insv1)
  Isadot      <-      (p3*Ina )  - (qg1 * Isa)
  Issdot      <-      (p4*Ins )  - ( (qg2 + qh)*Iss )
  Isav1dot    <-      (p3*Inav1) - (qg1*Isav1)
  Issv1dot    <-      (p4*Insv1) - ( (qg2_v1+qh_v1)*Issv1 )
  Hdot        <-      (h*Ins)   + (qh*Iss)  - (g3+mu)*H
  Hv1dot      <-      (h_v1*Insv1) + (qh_v1*Issv1) - (g3+mu)*Hv1
  Rdot        <-      (g1*Ina)  + (qg1*Isa) + (g2*Ins)  + (qg2*Iss)   + (g3*H ) - (nu*R)  -v1rater
  Rv1dot      <-      (g1*Inav1) + (qg1*Isav1)  + (g2_v1 * Insv1 )  + (qg2_v1.*Issv1)  + (g3*Hv1)  + v1rater  
  DEdot       <-      (mu*H )    + (mu*Hv1)
  
  
  f <- c(Sdot, Edot, Inadot, Isadot, Insdot, Issdot, Hdot, Rdot, Shv1dot, Sv1dot, 
         Ev1dot, Inav1dot, Isav1dot, Insv1dot, Issv1dot, Hv1dot, Rv1dot, DEdot, Vdot, Sndot)
  
  # Keeping only the positives
  f <- pmax(f, 0) 
  
  # for g
  g <- h * (Ins + Iss) + h_v1 * (Insv1 + Issv1)
  
  return(list(f, g))
}
  