install.packages("pracma")
library(pracma)

p2Cost <- function(data, dis, p2, g) {
  
  t <- g[, 1]
  lx <- length(data$obj)
  ln <- lx + 4
  
  # Initialising the variables cost and ccost_t that p2Cost returns
  
  cost <- matrix(0, nrow = 10, ncol = ln)
  #ccost_t <- matrix(0, nrow = nrow(g),  ) 
  
  # VLYL
  
  deaths <-  g[ nrow(g), ( 1 + 2*lx + 4*ln + 1):(1 + 2*lx + 5*ln)]
  cost[1,] <- deaths
  
  lyl <- deaths*data$lgh
  cost[2,] <- lyl
  
  vlyl <- lyl*data$vly
  cost[3,] <- vlyl
  
  deaths <- g[, (1 + 2*lx + 4*ln +1):(1 +2*lx+5*ln)]
  ccost_t[,1:ln] <- deaths * data$lgh * data$vly
  
  # VSYL
  
  Stu <- lx + 2
  students <- data$NNs[Stu]
  cost[4, lx + 1:2] <- students
  isoasy <- g[ ,1 + 2 * lx + 0* ln+ Stu] * 14 / (dis$Tay - p2$dur)
  isosym <- g[ ,1 + 2 * lx + 1* ln+Stu]
  isorec <- g[ ,1 + 2 * lx + 1* ln+Stu]*(14-dis$Ts[Stu] + p2$dur)/(dis$Ts[Stu] - p2$dur)
  nissym <- g[ ,1 + 2 * lx + 2* ln+Stu]
  hospts <- g[ ,1 + 2 * lx + 3* ln+Stu]
  deaths <- g[ ,1 + 2 * lx + 4* ln+Stu]
  abs   <- isoasy + isosym + isorec + nissym + hospts + deaths
  absint <- trapz(t, abs)/365
  cost[5, lx + 1] <- absint
  vsyl_sts <- absint * data$vsy
  cost[6, lx + 1] <- vsyl_sts
  
  # Student Demand
   
  pres <- student - abs
  presl <- pres*(1 - g[, 1 + data$EdInd])
  preslint <- trapz(t, presl) / 365
  cost[5, lx + 2] <- preslint
  vsyl_std <- preslint * data$vsy
  cost[6, lx + 2] <- vsyl_std
  
  ccost_t[, ln + lx + 1] <- cumtrapz(t, abs, 1) / 365 * data$vsy
  ccost_t[, ln + lx + 2] <- cumtrapz(t, presl, 1) / 365 * data$vsy
  
  # SGDPL
  
  notEd <- c(1:(data$EdInd - 1), (data$EdInd + 1):lx)
  
  # Labour Supply
  
  hw <- g[, (1 + 1*lx + notEd)]
  isoasy <- g[, (1 + 2 * lx + 0* ln + notEd)] * (1 - hw) * 14 / (dis$Tay - p2$dur)  # 14-day isolation period
  isosym <- g[, (1 + 2 * lx + 0* ln + notEd)]
  isorec <- g[, (1 + 2 * lx + ln + notEd)] * (1 - hw) * (14 - dis$Ts[notEd] + p2$dur) / (dis$Ts[notEd] - p2$dur)
  nissym <- g[, (1 + 2 * lx + 2 * ln + notEd)]
  hospts <- g[, (1 + 2 * lx + 3 * ln + notEd)]
  deaths <- g[, (1 + 2 * lx + 4 * ln + notEd)] # number of workers absent
  abspc <- pmax((isoasy + isosym + isorec + nissym + hospts + deaths) / data$NNs[notEd], 0)  # Percentage of workers absent
  prespc <- 1 - abspc  # Percentage of workers present
  presx <- prespc^data$alp  # Percentage of GDP output
  absx <- 1 - presx  # Percentage of GDP lost
  absxint <- trapz(t, absx)   
  gdpl_lbs <- absxint * data$obj[notEd]
  cost[7, notEd] <- gdpl_lbs
  
  # Labour Demand
  
  w <- g[, notEd]
  x <- w^data$alp
  xint <- diff(t) * (1 - head(x, -1))
  gdpl_lbd <- xint * data$obj[notEd]
  cost[8, notEd] <- gdpl_lbd
  cost[9, notEd] <- 0
  
  # Medium-term
  
  cost[10, notEd] <- 0
  
  ccost_t[, (2 * ln + notEd)] <- cumtrapz(t, absx) / length(t) * data$obj[notEd]
  ccost_t[, (3 * ln + notEd)] <- cumtrapz(t, 1 - x) / length(t) * data$obj[notEd]
  
  return(cost, ccost_t)
  
}