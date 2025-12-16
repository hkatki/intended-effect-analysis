# Conceal-Reveal
# HKatki 5 May 2022

setwd("/Users/katkih/hkCurrent/Dropbox/Liquid-Biopsy/2021-MCED-trial-proposal/Conceal-Reveal/R")

packages <- c("xlsx","scales","lmtest","dplyr","ggplot2","survival","gmodels","coxph.risk","geepack",
              "MESS","psych","Hmisc","glmnet","boot","zoo", "scales")
lapply(packages, require, c = T)


conceal.reveal <- function(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02,
                           RR=0.9, RR_neg=1, RR_pos=0.86667,
                           frac.g0.plus=1, frac.g0.minus=1,
                           frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,
                           censor.screen.Dplus = 0, censor.control.Dplus = 0,
                           censor.screen.Dminus = 0, censor.control.Dminus = 0,
                           move.screen.Mplus.Dplus = 0, move.screen.Mminus.Dplus = 0,
                           sim=FALSE)
{
  # RR_neg = P(D+|M-,G=1) / P(D+|M-,G=0)
  if (RR>RR_neg) stop("unrealistic for RR>RR_neg")
  # RR_pos = P(D+|M+,G=1) / P(D+|M+,G=0)
  if (RR<RR_pos) stop("unrealistic for RR<RR_pos")
  
  # Fix P(D+|G=0) as population disease mortality and expected mortality RR=P(D+|G=1)/P(D+|G=0)
  # Calculate P(D+|G=1) = RR*P(D+|G=0), and true RD=P(D+|G=0)-P(D+|G=1)
  pD_g1 <- RR*pD_g0
  RD <- pD_g0 - pD_g1

  # P(D+|M-,G=0) = (RR_pos*P(D+|G=0) - P(D+|G=1)) / [P(M-)*(RR_pos-RR_neg)]
  pD_g0.neg <- ifelse(RR_pos-RR_neg==0,
                      pD_g0, # under the null, P(D+|M-,G=0)=P(D+|G=0)
                      (RR_pos*pD_g0 - pD_g1) / ( (1-p_Mplus)*(RR_pos-RR_neg) ) )
  
  # P(D+|M-,G=1) = RR_neg * P(D+|M-,G=0)
  pD_g1.neg <- RR_neg * pD_g0.neg
  
  # P(D+|M+,G=0) = (RR_neg*P(D+|G=0) - P(D+|G=1)) / [P(M+)(RR_neg-RR_pos)]
  pD_g0.pos <- ifelse(RR_neg-RR_pos==0,
                      pD_g0, # under the null, P(D+|M-,G=0)=P(D+|G=0)
                      (RR_neg*pD_g0 - pD_g1) / ( p_Mplus*(RR_neg-RR_pos) ) )
  
  if (pD_g0.pos>1) stop(paste("pD_g0.pos=",pD_g0.pos,"is >1 because (RR_neg*pD_g0 - pD_g1) > p_Mplus*(RR_neg-RR_pos)"))
  
  # P(D+|M+,G=1) = RR_pos * P(D+|M+,G=0)
  pD_g1.pos <- RR_pos * pD_g0.pos
  
  # P(D+|M-,G=1) = RR_neg * P(D+|M-,G=0) 
  pD_g1.neg <- RR_neg * pD_g0.neg
  
  # RD_neg = P(D+|M-,G=0) - P(D+|M-,G=1)
  RD_neg <- pD_g0.neg - pD_g1.neg
    
  # Given RD and RD_neg=P(D+|M-,G=0)-P(D+|M-,G=1), solve for RD_pos=P(D+|M+,G=0)-P(D+|M+,G=1)
  # RD_pos <- (RD-RD_neg*(1-p_Mplus)) / p_Mplus
  
  # RD_pos = P(D+|M+,G=0) - P(D+|M+,G=1)
  RD_pos <- pD_g0.pos - pD_g1.pos
  
  
  # Deal with constraints
  # stop if P(D+|M-,G) < 0, that is, we hit the boundary of no outcomes among never-positives
  if (pD_g0.neg < 0 | pD_g1.neg < 0) stop("P(D+|M-,G) <0, meaning negative outcomes in never-positives. Check RR and RRpos.")
  # stop if P(D+|M+,G) < 0, that is, we hit the boundary of no outcomes among ever-positives
  if (pD_g0.pos < 0 | pD_g1.pos < 0) stop("P(D+|M+,G) <0, meaning negative outcomes in ever-positives. Check RR and RRpos.")
  
  # Specify PPV=P(D+|M+,G=0)
  # Check for compatibilty with marginals by checking if sens=P(M+|D+,G=0)=PPV*P(M+)/P(D+|G=0) > 1
  # I use 1.001 because sometimes there are numerical issues and sens is slightly > 1
  # pD_g0.pos <- 0.3
  if (pD_g0.pos*p_Mplus/pD_g0>1.001) stop("incompatible PPV conditional because sens>1")
  
  
  # if (pD_g0.pos*p_Mplus/pD_g0<RD_pos) stop("PPV too small given large RD_pos")

  # # Now calculate the rest of the conditionals
  # # P(D+|M+,G=1) = P(D+|M+,G=0) - RD_pos
  # pD_g1.pos <- pD_g0.pos - RD_pos
  # paste("conceal-reveal RR:", RR_pos <- pD_g1.pos/pD_g0.pos)
  # # back out P(D+|M-,G=0) from P(D+|G=0) = P(D+|M-,G=0)P(M-) + P(D+|M+,G=0)P(M+)
  # pD_g0.neg <- (pD_g0 - pD_g0.pos*p_Mplus) / (1-p_Mplus)
  # # P(D+|M-,G=1) = P(D+|M-,G=0) - RD_pos
  # pD_g1.neg <- pD_g0.neg - RD_neg
  # RR_neg <- pD_g1.neg/pD_g0.neg # other conceal-reveal RR

  # Theoretical ratio of Z-statistics: RD_pos/RD * P(M+) * sqrt{P(M+)/(P(M+|D+)P(M+|D-)}, where:
  # P(M+|D+) = P(D+|M+)P(M+)/P(D+)
  # P(M+|D-) = P(D-|M+)(M+)/P(D-)
  # P(D+|M+) = P(D+|M+,G=0)P(G=0) + P(D+|M+,G=1)P(G=1)
  pD_pos <- pD_g0.pos*p_g0 + pD_g1.pos*(1-p_g0)
  
  # Make sure this expression agrees with the above, valid only for P(G=1)=0.5: equal sample size both arms
  # P(D+|M+) = {0.5*(1+RR_pos)*P(D+|G=0)/P(M+)} * (RR_neg-RR)/(RR_neg-RR_pos)
  pD_pos_check <- (0.5*(1+RR_pos)*pD_g0/p_Mplus) * (RR_neg-RR)/(RR_neg-RR_pos)
  
  # P(D+|M-) = P(D+|M-,G=0)P(G=0) + P(D+|M-,G=1)P(G=1)
  pD_neg <- pD_g0.neg*p_g0 + pD_g1.neg*(1-p_g0)
  # P(D+) = P(D+|M+)P(M+) + P(D+|M-)P(M-)
  pD <- pD_pos*p_Mplus + pD_neg*(1-p_Mplus)
  # P(M+|D+) = P(D+|M+)P(M+)/P(D+)
  pMplus_Dplus <- pD_pos*p_Mplus/pD
  # P(M+|D-) = P(D-|M+)P(M+)/P(D-)
  pMplus_Dminus <- (1-pD_pos)*p_Mplus/(1-pD)

  # theoretical Z-ratio
  # Z_ratio <- RD_pos/RD * p_Mplus * sqrt(p_Mplus/(pMplus_Dminus*pMplus_Dplus))
  Z_ratio <- (1 - (RD_neg/RD)*(1-p_Mplus)) * sqrt(p_Mplus/(pMplus_Dminus*pMplus_Dplus))
  
  
  ###
  # Print out expected value of 2x2 tables
  ###
  # M+ table
  # P(D+,G=1|M+) * #(M+)
  a1 <- (pD_g1.pos*(1-p_g0) * n*p_Mplus) * (1-censor.screen.Dplus) 
  move.a1.to.b1 <- (pD_g1.pos*(1-p_g0) * n*p_Mplus) * (censor.screen.Dplus)
  # P(D+,G=0|M+) * #(M+)
  a0 <- (pD_g0.pos*p_g0 * n*p_Mplus) * (1-censor.control.Dplus)
  move.a0.to.b0 <- (pD_g0.pos*p_g0 * n*p_Mplus) * (censor.control.Dplus)
  # P(D-,G=1|M+) * #(M+)
  c1 <- ((1-pD_g1.pos)*(1-p_g0) * n*p_Mplus) * (1-censor.screen.Dminus) 
  move.c1.to.d1 <- ((1-pD_g1.pos)*(1-p_g0) * n*p_Mplus) * (censor.screen.Dminus)
  # P(D-,G=0|M+) * #(M+)
  c0 <- ((1-pD_g0.pos)*p_g0 * n*p_Mplus) * (1-censor.control.Dminus)
  move.c0.to.d0 <- ((1-pD_g0.pos)*p_g0 * n*p_Mplus) * (censor.control.Dminus)

  # M- table
  # P(D+,G=1,M-) * #(M-)
  b1 <- pD_g1.neg*(1-p_g0) * n*(1-p_Mplus) + move.a1.to.b1
  # P(D+,G=0,M-) * #(M-)
  b0 <- pD_g0.neg*p_g0 * n*(1-p_Mplus) + move.a0.to.b0
  # P(D-,G=1,M-) * #(M-)
  d1 <- (1-pD_g1.neg)*(1-p_g0) * n*(1-p_Mplus) + move.c1.to.d1
  # P(D-,G=0,M-) * #(M-)
  d0 <- (1-pD_g0.neg)*p_g0 * n*(1-p_Mplus) + move.c0.to.d0
  
  # Now account for those flipping from D+ to D- in screen arm due to non-adherence
  # Do among ever-positives
  move.c1.to.a1 <- c1 * move.screen.Mplus.Dplus
  a1 <- a1 + move.c1.to.a1
  c1 <- c1 - move.c1.to.a1
  # Do among never-positives
  move.d1.to.b1 <- d1 * move.screen.Mminus.Dplus
  b1 <- b1 + move.d1.to.b1
  d1 <- d1 - move.d1.to.b1
  
  # compute margins and final tables
  n1pos <- a1 + c1 
  n0pos <- a0 + c0 
  Mpos <- matrix( c(a1,a0,c1,c0),nrow=2,byrow=TRUE,
                  dimnames=list(c("D+","D-"),c("screen","control")))
  n1neg <- b1 + d1 #n*(1-p_g0)*(1-p_Mplus)
  n0neg <- b0 + d0 #n*p_g0*(1-p_Mplus)
  # d1 <- n1neg - b1
  # d0 <- n0neg - b0
  Mneg <- matrix( c(b1,b0,d1,d0),nrow=2,byrow=TRUE,
                  dimnames=list(c("D+","D-"),c("screen","control")))
  
  # Now include censoring of ever-test-positivity in the test-positive table
  # These people will be misclassified into the test-negative table
  print("Ever-positive table")
  print(Mpos) ; Mpos.test <- prop.test(t(Mpos), correct=F)
  print(paste("RR_pos:", Mpos.test$estimate[1]/Mpos.test$estimate[2],"RD_pos:", diff(Mpos.test$estimate), 
              "p-value:",Mpos.test$p.val,
              "per-arm N for 90% power:",round((n/2)*((qnorm(1-0.9)-1.96)/qnorm(Mpos.test$p.val/2))^2)))
  
  print("Never-positive table")
  print(Mneg); Mneg.test <- prop.test(t(Mneg), correct=F)
  print(paste("RR_neg:", Mneg.test$estimate[1]/Mneg.test$estimate[2],"RD_neg:", diff(Mneg.test$estimate), 
              "p-value:",Mneg.test$p.val))
  
  print("Total table")
  print(all<-Mneg+Mpos); all.test <- prop.test(t(all), correct=F)
  print(paste("RR_all:", all.test$estimate[1]/all.test$estimate[2],"RD_all:", diff(all.test$estimate), 
              "p-value:",all.test$p.val,
              "per-arm N for 90% power:",round((n/2)*((qnorm(1-0.9)-1.96)/qnorm(all.test$p.val/2))^2)))
  
  # Do Peter Sasieni's "targeted IE"
  targeted <- matrix( c(a1,a0,d1+c1+b1,d0+c0+b0),nrow=2,byrow=TRUE,
                      dimnames=list(c("D+ and M+","otherwise"),c("screen","control")))
  print(targeted) ; targeted.test <- prop.test(t(targeted), correct=F)
  print(paste("Targeted RR_pos:", targeted.test$estimate[1]/targeted.test$estimate[2],
              "RD_pos:", diff(targeted.test$estimate), 
              "p-value:",targeted.test$p.val,
              "per-arm N for 90% power:",round((n/2)*((qnorm(1-0.9)-1.96)/qnorm(targeted.test$p.val/2))^2)))
  
  # # Check if RRpos = RR * P(M+|D+,G=1)/P(M+|D+,G=0)
  # pMplus_Dplus_g1 <- a1/(a1+b1)
  # pMplus_Dplus_g0 <- a0/(a0+b0)
  # RR_pos.check <- RR * pMplus_Dplus_g1/pMplus_Dplus_g0
  # print(paste("RR_pos.check=",RR_pos.check,"RR_pos=",RR_pos))
  # 
  # # Check if RRneg = RR * P(M-|D+,G=1)/P(M-|D+,G=0)
  # pMminus_Dplus_g1 <- 1-pMplus_Dplus_g1
  # pMminus_Dplus_g0 <- 1-pMplus_Dplus_g0
  # RR_neg.check <- RR * pMminus_Dplus_g1/pMminus_Dplus_g0
  # print(paste("RR_neg.check=",RR_neg.check,"RR_neg=",RR_neg))
  # # check if RRneg = (P(D+│G=0)⋅P(M-│D+,G=1)⋅RR)/(P(M-│G=1)-P(M-│D-,G=0)⋅P(D-|G=0))
  # pMminus_Dminus_g0 <- d0/(d0+c0)
  # print(paste("specificity in control arm=",pMminus_Dminus_g0))
  # RR_neg.check1 <- (pD_g0*pMminus_Dplus_g1*RR) / ((1-p_Mplus)-pMminus_Dminus_g0*(1-pD_g0))
  # print(paste("RR_neg.check1=",RR_neg.check1,"RR_neg=",RR_neg))
  # pMminus_Dminus_g1 <- d1/(d1+c1)
  # print(paste("specificity in screen arm=",pMminus_Dminus_g1))
  # RR_neg.spec <- (pD_g0*pMminus_Dplus_g1*RR) / ((1-p_Mplus)-pMminus_Dminus_g1*(1-pD_g0))
  # print(paste("RR_neg.spec=",RR_neg.spec,"RR_neg=",RR_neg))
  # 
  # # check if RRpos = (P(D+│G=0)⋅P(M+│D+,G=1)⋅RR)/(P(M+│G=1)-P(M+│D-,G=0)⋅P(D-|G=0))
  # RR_pos.check1 <- (pD_g0*(1-pMminus_Dplus_g1)*RR) / (p_Mplus-(1-pMminus_Dminus_g0)*(1-pD_g0))
  # print(paste("RR_pos.check1=",RR_pos.check1,"RR_pos=",RR_pos))
  # # Inputting 1-spec from the screen arm (1-pMminus_Dminus_g1), is the approximation good?
  # RR_pos.spec <- (pD_g0*(1-pMminus_Dplus_g1)*RR) / (p_Mplus-(1-pMminus_Dminus_g1)*(1-pD_g0))
  # print(paste("RR_pos.spec=",RR_pos.spec,"RR_pos=",RR_pos))
  
  ##
  # Output power for RDpos given fixed power for RD
  # This is pnorm(Zratio * mean_Z - 1.96)
  # mean Z is that from Goodman 1992 for replication probability
  # i.e. power=0.5 is mean Z stat = 1.96, power=0.8 is mean Z stat = 2.802, power=0.9 is mean Z stat = 3.24
  ##
  if (sim==FALSE) { 
    # RD_power <- 0.9
    # return(pnorm(Z_ratio*(qnorm(1-0.05/2)+qnorm(RD_power))-1.96)) 
    # return(Z_ratio)
    tables <- NULL
    tables$all <- all ; tables$Mpos <- Mpos ; tables$Mneg <- Mneg 
    return(tables)
  }
    # return(data.frame(
    #   RR=RR, RR_pos=RR_pos, RR_neg=RR_neg,
    #   pD_g1.pos=pD_g1.pos,
    #   pD_g0.pos=pD_g0.pos,
    #   pD_g1.neg=pD_g1.neg,
    #   pD_g0.neg=pD_g0.neg,
    #   pMplus_Dminus=pMplus_Dminus,
    #   pMplus_Dplus=pMplus_Dplus,
    #   pD_pos=pD_pos,
    #   pD_neg=pD_neg,
    #   pD=pD,
    #   RD=RD, RD_neg=RD_neg, RD_pos=RD_pos, Z_ratio=Z_ratio, 
    #   power.pos=pnorm(Z_ratio*2.802-1.96)
    # ))
  
  ##
  # initialize simulation
 
  # initialize simulation
  ##
  Z_neg <- Z_pos <- Z_all <- p.val_neg <- p.val_pos <- p.val_all <- RD_neg.sim <- RD_pos.sim <- RD_all.sim <- RR_neg.sim <- RR_pos.sim <- RR_all.sim <- numeric(nsim)
  Z_tar <- p.val_tar  <- RD_tar.sim  <- RR_tar.sim  <- RR_pos.noncomp.sim <- RR_tar.noncomp.sim <- RR_neg.noncomp.sim <- numeric(nsim)
  p.val_neg.subsamp <- p.val_pos.subsamp <- RD_neg.subsamp.sim <- RD_pos.subsamp.sim <- p.val_pos.noncomp <- p.val_tar.noncomp <- p.val_neg.noncomp <- numeric(nsim)
  RD_pos.loss.ss.sim <- RR_pos.loss.ss.sim <- RD_neg.loss.sim <- RD_pos.loss.sim <- RR_neg.loss.sim <- RR_pos.loss.sim <- p_Mplus.g0.hat <- p_Mplus.g1.hat <- check.RD_all.hat <- numeric(nsim)
  RD_pos.retest <- RD_pos.BF <- RD_neg.BF <- RD_neg.BF.retest <- RD_pos.BF.retest <- RR_neg.BF.retest <- RR_pos.BF.retest <- numeric(nsim)
  RD_pos_IE.sim <- RD_IE.sim <- RR_pos_IE.sim <- p.val.RD_pos_IE.sim <- numeric(nsim)
  
  # Simulate 2x2 tables
  set.seed(406551735)
  for(i in 1:nsim) {

    ##
    # Generate M+ and M- sample sizes, fixing total sample sizes n*p_g0 and n*(1-p_g0) in each group
    ##
    
    # #(M+,G=1) = n_1+
    g1.pos <- rbinom(1, n*(1-p_g0), p_Mplus) 
    
    # #(M+,G=0) = n_0+
    g0.pos <- rbinom(1, n*p_g0, p_Mplus) 
    
    # #(M-,G=1) = n_1-
    g1.neg <- n*(1-p_g0)-g1.pos
    
    # #(M-,G=0) = n_0-
    g0.neg <- n*p_g0-g0.pos
    
    ##
    # generate the M- and M+ tables
    ##
    
    # #(D+,G=1,M-) = b_1
    Dplus_g1.neg <- rbinom(1, g1.neg, pD_g1.neg)
    # #(D+,G=0,M-) = b_0
    Dplus_g0.neg <- rbinom(1, g0.neg, pD_g0.neg)
    # #(D+,G=1,M+) = a_1
    Dplus_g1.pos <- rbinom(1, g1.pos, pD_g1.pos)
    # #(D+,G=0,M+) = a_0
    Dplus_g0.pos <- rbinom(1, g0.pos, pD_g0.pos)

    # Calculate the D- cells of the tables
    Dminus_g1.pos <- g1.pos - Dplus_g1.pos # need to calculate true cell c_1
    Dminus_g0.pos <- g0.pos - Dplus_g0.pos # need to calculate true cell c_0
    Dminus_g1.neg <- g1.neg - Dplus_g1.neg # need to calculate true cell d_1
    Dminus_g0.neg <- g0.neg - Dplus_g0.neg # need to calculate true cell d_0
    

    
    ##
    # Simulate non-adherence to blood draws in both group: dropout
    # Now simulate dropout among the ever-positives in both screen and control arms
    # They get redistributed to the never-positives in both screen and control arms
    ##

    # First non-adherence among ever-positives in screen (G=1) and control (G=0) arms
    # Calculate the non-censoring fractions, separately for D+ and D-, in each arm (although these are unobservable)
    # The "move" variables are those who are unknown positivity but truly ever-positive - unobservable
    # #(D+,G=1,M+) = a_1
    move.a1.to.b1 <- rbinom(1, Dplus_g1.pos, censor.screen.Dplus)
    Dplus_g1.pos <-  Dplus_g1.pos - move.a1.to.b1
    # #(D+,G=0,M+) = a_0
    move.a0.to.b0 <- rbinom(1, Dplus_g0.pos, censor.control.Dplus)
    Dplus_g0.pos <- Dplus_g0.pos - move.a0.to.b0
    # #(D-,G=1,M+) = c_1
    move.c1.to.d1 <- rbinom(1, Dminus_g1.pos, censor.screen.Dminus)
    Dminus_g1.pos <- Dminus_g1.pos - move.c1.to.d1
    # #(D-,G=0,M+) = c_0
    move.c0.to.d0 <- rbinom(1, Dminus_g0.pos, censor.control.Dminus)
    Dminus_g0.pos <- Dminus_g0.pos - move.c0.to.d0
    
    # Second: non-adherence among never-positives in screen (G=1) and control (G=0) arms
    # rename the "move" variables from above as ".unknown" to clarify ("move" is only if unknown is treated as negative)
    # ".neg.unknown" is truly never-positive but observed as unknown (these are unobservable)
    # #(D+,G=1,M-) = b_1
    Dplus_g1.neg.unknown <- rbinom(1, Dplus_g1.neg, censor.screen.Dplus) 
    Dplus_g1.neg <-  Dplus_g1.neg - Dplus_g1.neg.unknown
    # #(D+,G=0,M-) = b_0
    Dplus_g0.neg.unknown <- rbinom(1, Dplus_g0.neg, censor.control.Dplus)
    Dplus_g0.neg <- Dplus_g0.neg - Dplus_g0.neg.unknown
    # #(D-,G=1,M-) = d_1
    Dminus_g1.neg.unknown <- rbinom(1, Dminus_g1.neg, censor.screen.Dminus)
    Dminus_g1.neg <- Dminus_g1.neg - Dminus_g1.neg.unknown
    # #(D-,G=0,M-) = d_0
    Dminus_g0.neg.unknown <- rbinom(1, Dminus_g0.neg, censor.control.Dminus)
    Dminus_g0.neg <- Dminus_g0.neg - Dminus_g0.neg.unknown
    
    # Third: Calculate the 2x2 table for group with unknown positivity ".unknown"
    # Key point: each term in the sums below are unobservable, but the sum is observable
    # First term is the 
    Dplus_g1.unknown <- move.a1.to.b1 + Dplus_g1.neg.unknown 
    Dplus_g0.unknown <- move.a0.to.b0 + Dplus_g0.neg.unknown
    Dminus_g1.unknown <- move.c1.to.d1 + Dminus_g1.neg.unknown
    Dminus_g0.unknown <- move.c0.to.d0 + Dminus_g0.neg.unknown
    
    # # NEED TO COMMENT THIS CODE OUT IF I WANT TO USE THE NON-COMPLIANCE CORRECTIONS TO CORRECT RR_NEG
    # # For analysis that treats unknown as "negative", move the dropped out in M+ into the never-positive
    # # #(D+,G=1,M-) = b_1
    # Dplus_g1.neg <- Dplus_g1.neg + move.a1.to.b1
    # # #(D+,G=0,M-) = b_0
    # Dplus_g0.neg <- Dplus_g0.neg + move.a0.to.b0
    # # #(D-,G=1,M-) = d_1
    # Dminus_g1.neg <- Dminus_g1.neg + move.c1.to.d1
    # # #(D-,G=0,M-) = d_0
    # Dminus_g0.neg <- Dminus_g0.neg + move.c0.to.d0
    
    
    # # I HAVEN'T HANDLED NON-COMPLIANCE CORRECTIONS THAT ACCOUNT FOR CHANGING OUTCOMES IN THE SCREEN ARM
    # # IN THE INPUTS, SET THIS TO 0 FOR NOW
    # # In screen arm, some of those who move have their outcomes change
    # # Redistribute some non-outcomes to outcomes in the screen arm
    # # Do separately for observed ever-positives and observed never-positives
    # # By "observed" I mean after non-adherence has already flipped some ever-positives to never-positive
    # 
    # # Observed ever-positives: move #(D-,G=1,M+) = c_1 to #(D+,G=1,M+) = a_1
    # move.c1.to.a1 <- rbinom(1, Dminus_g1.pos, move.screen.Mplus.Dplus)
    # Dplus_g1.pos <-  Dplus_g1.pos + move.c1.to.a1
    # Dminus_g1.pos <- Dminus_g1.pos - move.c1.to.a1
    # 
    # # Observed ever-positives: move #(D-,G=1,M-) = d_1 to #(D+,G=1,M-) = b_1
    # move.d1.to.b1 <- rbinom(1, Dminus_g1.neg, move.screen.Mminus.Dplus)
    # Dplus_g1.neg <- Dplus_g1.neg + move.d1.to.b1
    # Dminus_g1.neg <- Dminus_g1.neg - move.d1.to.b1
    # 
    # # recalculate margins, ASSUMING UNKNOWNS ARE IN NEVER-POSITIVE VARIABLES .NEG
    # g1.pos <- Dplus_g1.pos + Dminus_g1.pos
    # g0.pos <- Dplus_g0.pos + Dminus_g0.pos
    # g1.neg <- Dplus_g1.neg + Dminus_g1.neg
    # g0.neg <- Dplus_g0.neg + Dminus_g0.neg
    
    # recalculate margins, where the unknown are separated from the never-positive
    g1.pos <- Dplus_g1.pos + Dminus_g1.pos 
    g0.pos <- Dplus_g0.pos + Dminus_g0.pos 
    g1.neg <- Dplus_g1.neg + Dminus_g1.neg 
    g0.neg <- Dplus_g0.neg + Dminus_g0.neg 
    g1.unknown <- Dplus_g1.unknown + Dminus_g1.unknown
    g0.unknown <- Dplus_g0.unknown + Dminus_g0.unknown
    
    # Define the 4 cells for the all table and its margins
    Dplus_g1.all <- Dplus_g1.neg+Dplus_g1.pos+Dplus_g1.unknown
    Dplus_g0.all <- Dplus_g0.neg+Dplus_g0.pos+Dplus_g0.unknown
    Dminus_g1.all <- Dminus_g1.neg+Dminus_g1.pos+Dminus_g1.unknown
    Dminus_g0.all <- Dminus_g0.neg+Dminus_g0.pos+Dminus_g0.unknown
    g1.all <- g1.neg+g1.pos+g1.unknown
    g0.all <- g0.neg+g0.pos+g0.unknown
    
    # Uncomment these to print out example tables
    # print("ever-positive table with dropouts\n")
    # print(matrix( c(Dplus_g1.pos,Dplus_g0.pos,Dminus_g1.pos,Dminus_g0.pos),nrow=2,byrow=TRUE,
    #                 dimnames=list(c("D+","D-"),c("screen","control"))) )
    # print("never-positive table with dropouts\n")
    # print(matrix( c(Dplus_g1.neg,Dplus_g0.neg,Dminus_g1.neg,Dminus_g0.neg),nrow=2,byrow=TRUE,
    #                 dimnames=list(c("D+","D-"),c("screen","control"))) )
    # print("unknown-positive table with dropouts\n")
    # print(matrix( c(Dplus_g1.unknown,Dplus_g0.unknown,Dminus_g1.unknown,Dminus_g0.unknown),nrow=2,byrow=TRUE,
    #               dimnames=list(c("D+","D-"),c("screen","control"))) )
    # print("Targeted table with dropouts\n")
    # print(matrix( c(Dplus_g1.pos, Dplus_g0.pos, g1.all-Dplus_g1.pos, g0.all-Dplus_g0.pos),nrow=2,byrow=TRUE,
    #               dimnames=list(c("D+","D-"),c("screen","control"))) )
    # print(paste("#(M+,G=1) =",g1.pos))
    # print(paste("#(M+,G=0) =",g0.pos))
    # print(paste("#(M-,G=1) =",g1.neg))
    # print(paste("#(M-,G=0) =",g0.neg))
    # print(paste("#(M unknown,G=1) =",g1.unknown))
    # print(paste("#(M unknown,G=0) =",g0.unknown))

      
    ##
    # Analysis 
    ##
    
    # calculate statistics for M-, M+, and marginal tables
    neg.out <- prop.test(c(Dplus_g1.neg, Dplus_g0.neg), c(g1.neg, g0.neg), correct=F)
    pos.out <- prop.test(c(Dplus_g1.pos, Dplus_g0.pos), c(g1.pos, g0.pos), correct=F)
    all.out <- prop.test(c(Dplus_g1.all, Dplus_g0.all), c(g1.all, g0.all), correct=F)
    tar.out <- prop.test(c(Dplus_g1.pos, Dplus_g0.pos), c(g1.all, g0.all), correct=F)
    
    # Correction for ever-positives (RRpos and RRtar) using observable data from the unknown-positive table
    # The correction factor is the ratio of the noncensored fractions of cases in the screen vs control arms
    Dplus_g0.pos.corrected <- Dplus_g0.pos * (1-Dplus_g1.unknown/Dplus_g1.all) / 
                                             (1-Dplus_g0.unknown/Dplus_g0.all)
    Dminus_g0.pos.corrected <- Dminus_g0.pos * (1-Dminus_g1.unknown/Dminus_g1.all) / 
                                               (1-Dminus_g0.unknown/Dminus_g0.all)
    g0.pos.corrected <- Dplus_g0.pos.corrected + Dminus_g0.pos.corrected
    
    # Correction for never-positives (RRneg) using observable data from the unknown-positive table
    # The correction factor is the ratio of the noncensored fractions of cases in the screen vs control arms
    Dplus_g0.neg.corrected <- Dplus_g0.neg * (1-Dplus_g1.unknown/Dplus_g1.all) / 
      (1-Dplus_g0.unknown/Dplus_g0.all)
    Dminus_g0.neg.corrected <- Dminus_g0.neg * (1-Dminus_g1.unknown/Dminus_g1.all) / 
      (1-Dminus_g0.unknown/Dminus_g0.all)
    g0.neg.corrected <- Dplus_g0.neg.corrected + Dminus_g0.neg.corrected
    
    # Uncomment these to print out example tables
    # print("ever-positive table corrected for dropouts\n")
    # print(matrix( c(Dplus_g1.pos, Dplus_g0.pos.corrected, Dminus_g1.pos, Dminus_g0.pos.corrected),nrow=2,byrow=TRUE,
    #               dimnames=list(c("D+","D-"),c("screen","control"))) )
    # print("targeted table corrected for dropouts\n")
    # print(matrix( c(Dplus_g1.pos, Dplus_g0.pos.corrected, g1.all-Dplus_g1.pos, g0.all-Dplus_g0.pos.corrected),nrow=2,byrow=TRUE,
    #               dimnames=list(c("D+","D-"),c("screen","control"))) )


    # Calculate RRpos and RRneg with control-arm non-compliance corrected to equal that of the screen-arm
    pos.noncomp.out <- prop.test(c(Dplus_g1.pos, Dplus_g0.pos.corrected), c(g1.pos, g0.pos.corrected), correct=F)
    neg.noncomp.out <- prop.test(c(Dplus_g1.neg, Dplus_g0.neg.corrected), c(g1.neg, g0.neg.corrected), correct=F)
    # Calculate RRtar with control-arm non-compliance corrected to equal that of the screen-arm
    tar.noncomp.out <- prop.test(c(Dplus_g1.pos, Dplus_g0.pos.corrected), c(g1.all, g0.all), correct=F)
    
    
    
    # RDs and RRs
    RD_neg.sim[i] <- diff(neg.out$estimate)
    RD_pos.sim[i] <- diff(pos.out$estimate)
    RD_all.sim[i] <- diff(all.out$estimate)
    RD_tar.sim[i] <- diff(tar.out$estimate)
    RR_neg.sim[i] <- neg.out$estimate[1]/neg.out$estimate[2]
    RR_pos.sim[i] <- pos.out$estimate[1]/pos.out$estimate[2]
    RR_all.sim[i] <- all.out$estimate[1]/all.out$estimate[2]
    RR_tar.sim[i] <- tar.out$estimate[1]/tar.out$estimate[2]
    RR_pos.noncomp.sim[i] <- pos.noncomp.out$estimate[1]/pos.noncomp.out$estimate[2]
    RR_neg.noncomp.sim[i] <- neg.noncomp.out$estimate[1]/neg.noncomp.out$estimate[2]
    RR_tar.noncomp.sim[i] <- tar.noncomp.out$estimate[1]/tar.noncomp.out$estimate[2]
    
    # p-values for the RDs
    p.val_neg[i] <- neg.out$p.val
    p.val_pos[i] <- pos.out$p.val
    p.val_all[i] <- all.out$p.val
    p.val_tar[i] <- tar.out$p.val
    p.val_pos.noncomp[i] <- pos.noncomp.out$p.val
    p.val_neg.noncomp[i] <- neg.noncomp.out$p.val
    p.val_tar.noncomp[i] <- tar.noncomp.out$p.val
    
    # Z-stats: need to get the sign right because power<50% implies that mean Z-stat is negative
    Z_neg[i]  <- sign(RD_neg.sim[i])*sqrt(neg.out$stat)
    Z_pos[i]  <- sign(RD_pos.sim[i])*sqrt(pos.out$stat)
    Z_all[i]  <- sign(RD_all.sim[i])*sqrt(all.out$stat)
    Z_tar[i]  <- sign(RD_tar.sim[i])*sqrt(tar.out$stat)
    
    
    ##
    # Calculate RD_pos and its var under IE assumptions 
    # for a trial with no control arm test results
    ##
    
    # Usual IE assumption RDneg=0, but this can be varied to see how the variance changes
    # I set this to whatever is set for RD_neg for the simulation
    RD_neg_IE <- RD_neg

    # Observed ever-positive rate P(M+) in screen-arm (truth is p_Mplus), 
    # by IE assumed to be the same in the control-arm
    p_plus = g1.pos/(g1.neg+g1.pos) 
    
    # The RD calculated under IE assumptions: RD_IE
    # This requires RD_pos to be observed, which is not true if no control arm test results
    RD_IE.sim[i] = RD_pos.sim[i]*p_plus + RD_neg_IE*(1-RD_pos.sim[i])
    
    # Calculate RDpos_IE = (RD - RDneg_IE*P(M-)) / P(M+)
    # This formulation allows for non-zero RDneg_IE
    RD_pos_IE.sim[i] <-  (RD_all.sim[i] - RD_neg_IE*(1-p_plus) ) / p_plus
    
    # RRpos_IE = p_splus / (p_splus + RD_pos_IE), where usual RRpos = p_splus/p_cplus
    p_splus <- Dplus_g1.pos / g1.pos # outcome incidence in screen-arm: observed
    p_cplus <- Dplus_g0.pos / g0.pos # outcome incidence in control-arm: NOT observed
    RR_pos_IE.sim[i] <- p_splus / (p_splus + RD_pos_IE.sim[i])
    
    ##
    # Uncomment to print tables
    ##
    # print(paste("\n Simulation:",i))
    # 
    # # Usual analysis table
    # print(all <- matrix( c(Dplus_g1.neg+Dplus_g1.pos, Dplus_g0.neg+Dplus_g0.pos,
    #                        (g1.neg+g1.pos) - (Dplus_g1.neg+Dplus_g1.pos),
    #                        (g0.neg+g0.pos) - (Dplus_g0.neg+Dplus_g0.pos) ),
    #                      nrow=2,byrow=TRUE, dimnames=list(c("D+","D-"),c("screen","control"))))
    # print(paste("Observed RD_all:", diff(all.out$estimate), "p-value:",all.out$p.val))
    # 
    # # True Ever-Positive table
    # print(pos <- matrix( c(Dplus_g1.pos,Dplus_g0.pos,g1.pos-Dplus_g1.pos,g0.pos-Dplus_g0.pos),
    #                      nrow=2,byrow=TRUE, dimnames=list(c("D+","D-"),c("screen","control"))))
    # print(paste("Observed RD_pos:", diff(pos.out$estimate), "p-value:",pos.out$p.val))
    # 
    # # True never-positive table
    # print(neg <- matrix( c(Dplus_g1.neg,Dplus_g0.neg,g1.neg-Dplus_g1.neg,g0.neg-Dplus_g0.neg),
    #                      nrow=2,byrow=TRUE, dimnames=list(c("D+","D-"),c("screen","control"))))
    # print(paste("Observed RD_neg:", diff(neg.out$estimate), "p-value:",neg.out$p.val))
    # print("\n")
    # 
    # # IE Never-Positive table that observed screen arm but has to impute control arm
    # # Number of control-arm negatives under IE: #(G=1,M- |IE) = #control-arm * P(M+)
    g0.neg_IE <- (g0.neg+g0.pos)*(1-p_plus)
    #  Number of control-arm negatives with outcome under IE: #(G=1,M- |IE)  * P(D+|G=1,M-)
    Dplus_g0.neg_IE <- g0.neg_IE * (Dplus_g1.neg / g1.neg)
    # never-positives 2x2 matrix
    # print(neg_IE <- matrix( c(Dplus_g1.neg, Dplus_g0.neg_IE,
    #                           g1.neg-Dplus_g1.neg, g0.neg_IE-Dplus_g0.neg_IE),
    #                         nrow=2,byrow=TRUE, dimnames=list(c("D+","D-"),c("screen","control"))))
    # neg_IE.out <- prop.test(t(neg_IE), correct=F)
    # print(paste("IE RD_neg:", diff(neg_IE.out$estimate), "p-value:",neg_IE.out$p.val))

    # IE Ever-Positive table that observed screen arm but has to impute control arm
    #  Number of control-arm positives with outcome under IE: #control-arm-D+ - Dplus_g0.neg_IE
    g0.pos_IE <- (g0.neg+g0.pos)*p_plus
    Dplus_g0.pos_IE <- (Dplus_g0.pos+Dplus_g0.neg) - Dplus_g0.neg_IE
    pos_IE <- matrix( c(Dplus_g1.pos, Dplus_g0.pos_IE,
                              g1.pos-Dplus_g1.pos, g0.pos_IE-Dplus_g0.pos_IE),
                         nrow=2,byrow=TRUE, dimnames=list(c("D+","D-"),c("screen","control")))
    # print(pos_IE)
    if ( all(pos_IE>0) ) { 
      # only analyze IE ever-positive table if cell counts are positive
      pos_IE.out <- prop.test(t(pos_IE), correct=F)
      p.val.RD_pos_IE.sim[i] <- pos_IE.out$p.val
    }
    else {
      p.val.RD_pos_IE.sim[i] <- NA # Later on, fix problems with negative cell counts
    }
    # print(paste("IE RD_pos:", diff(pos_IE.out$estimate), "p-value:",pos_IE.out$p.val))
    

    
    ##
    # Now we only test a subsamples of G=0 for the marker M
    # frac.g0.plus:  The fraction of #G=0,D+ that we will observe
    # frac.g0.minus: The fraction of #G=0,D- that we will observe
    # We will estimate these fractions and weight by them
    # Note: we cannot estimate 4 sampling fractions (the above 2 by M+ and M-)
    #       because we cannot know the total M+ or M- in the G=0 arm
    ##
    
    # subsample among a_0=#(D+,M+,G=0)
    Dplus_g0.pos.subsamp <- rbinom(1,Dplus_g0.pos,frac.g0.plus)
    
    # subsample among c_0=#(D-,M+,G=0)
    Dminus_g0.pos.subsamp <- rbinom(1, g0.pos - Dplus_g0.pos,frac.g0.minus)

    # subsample among b_0=#(D+,M-,G=0)
    Dplus_g0.neg.subsamp <- rbinom(1,Dplus_g0.neg,frac.g0.plus)

    # subsample among d_0=#(D-,M-,G=0)
    Dminus_g0.neg.subsamp <- rbinom(1, g0.neg - Dplus_g0.neg,frac.g0.minus)

    # Estimate sampling fractions for D+ and for D- within arm G=0
    # These are P(subsampled|D+,G=0) and P(subsampled|D-,G=0)
    # Estimating sampling fractions is important for D+ because it is the rare event
    # Not important for D- because there are so many of them so true fraction is close 
    # to the estimated fraction
    frac.g0.plus.hat <- (Dplus_g0.pos.subsamp + Dplus_g0.neg.subsamp) /
                        (Dplus_g0.pos + Dplus_g0.neg)
    frac.g0.minus.hat<- (Dminus_g0.pos.subsamp + Dminus_g0.neg.subsamp) / 
                        ((g0.pos - Dplus_g0.pos) + (g0.neg - Dplus_g0.neg))

    ##
    # calculate stats for M-,M+ with a subsampled table in control arm, 
    # weight up each total
    ##
    
    # Do for M- table: RD_neg = p_control - p_screen
    # weight #(D+,M-,G=0,subsampled) * 1/P(subsampled|D+,G=0)
    Dplus_g0.neg.subsamp.weighted <- Dplus_g0.neg.subsamp*1/frac.g0.plus.hat
    # also weight #(D-,M-,G=0,subsampled) * 1/P(subsampled|D-,G=0)
    g0.neg.weighted <- Dplus_g0.neg.subsamp.weighted + 
                       Dminus_g0.neg.subsamp*1/frac.g0.minus.hat
    p_screen.neg  <- Dplus_g1.neg / g1.neg  # P(D+|M-,G=1)
    p_control.neg <- Dplus_g0.neg.subsamp.weighted / g0.neg.weighted  # P(D+|M-,G=0)
    RD_neg.subsamp.sim[i] <- p_control.neg - p_screen.neg
    
    # Do for M+ table: RD_pos = p_control - p_screen
    Dplus_g0.pos.subsamp.weighted <- Dplus_g0.pos.subsamp*1/frac.g0.plus.hat
    g0.pos.weighted <- Dplus_g0.pos.subsamp.weighted + 
                       Dminus_g0.pos.subsamp*1/frac.g0.minus.hat
    p_screen.pos  <- Dplus_g1.pos / g1.pos  # P(D+|M+,G=1)
    p_control.pos <- Dplus_g0.pos.subsamp.weighted / g0.pos.weighted  # P(D+|M+,G=0)
    RD_pos.subsamp.sim[i] <- p_control.pos - p_screen.pos
    
    
    ##
    # Now simulate loss of MCED signal in the control group, no signal strength available
    # Start with the subsampled group
    ##
    
    # Subsample among a_0=#(D+,M+,G=0): the fraction of these that test true positive M+
    Dplus_g0.pos.loss <- rbinom(1,Dplus_g0.pos.subsamp,frac.g0.plus.pos)

    # The a_0=#(D+,M+,G=0) that test false-negative, put them into b_0=#(D+,M-,G=0)
    Dplus_g0.neg.loss <- Dplus_g0.neg.subsamp + (Dplus_g0.pos.subsamp - Dplus_g0.pos.loss)

    # Subsample among c_0=#(D-,M+,G=0): the fraction of these that test true positive M+
    Dminus_g0.pos.loss <- rbinom(1, Dminus_g0.pos.subsamp ,frac.g0.minus.pos)

    # The c_0=#(D-,M+,G=0) that test false-negative, put them into d_0=#(D-,M-,G=0)
    Dminus_g0.neg.loss <- Dminus_g0.neg.subsamp + (Dminus_g0.pos.subsamp - Dminus_g0.pos.loss)
        
    # Estimate sampling fractions for D+ and for D- within arm G=0
    frac.g0.plus.hat <- (Dplus_g0.pos.loss + Dplus_g0.neg.loss) /
                        (Dplus_g0.pos + Dplus_g0.neg)
    frac.g0.minus.hat<- (Dminus_g0.pos.loss + Dminus_g0.neg.loss) /
                        ((g0.pos - Dplus_g0.pos) + (g0.neg - Dplus_g0.neg))

    # Do for M- table: RD_neg = p_control - p_screen
    Dplus_g0.neg.loss.weighted <- Dplus_g0.neg.loss*1/frac.g0.plus.hat
    g0.neg.weighted <- Dplus_g0.neg.loss.weighted + 
                       Dminus_g0.neg.loss*1/frac.g0.minus.hat
    p_screen.neg  <- Dplus_g1.neg / g1.neg  # P(D+|M-,G=1)
    p_control.neg <- Dplus_g0.neg.loss.weighted / g0.neg.weighted  # P(D+|M-,G=0)
    RD_neg.loss.sim[i] <- p_control.neg - p_screen.neg
    RR_neg.loss.sim[i] <- p_screen.neg / p_control.neg
    
    # Do for M+ table: RD_pos = p_control - p_screen
    Dplus_g0.pos.loss.weighted <- Dplus_g0.pos.loss*1/frac.g0.plus.hat
    g0.pos.weighted <- Dplus_g0.pos.loss.weighted + 
                       Dminus_g0.pos.loss*1/frac.g0.minus.hat
    p_screen.pos  <- Dplus_g1.pos / g1.pos  # P(D+|M+,G=1)
    p_control.pos <- Dplus_g0.pos.loss.weighted / g0.pos.weighted  # P(D+|M+,G=0)
    RD_pos.loss.sim[i] <- p_control.pos - p_screen.pos
    RR_pos.loss.sim[i] <- p_screen.pos / p_control.pos
    
    
    # Check if estimated P(M+|G=0) is close to the true p_Mplus we input
    # This will not be true if there is any dropout
    p_Mplus.g0.hat[i] <- g0.pos.weighted/(g0.pos+g0.neg)
    p_Mplus.g1.hat[i] <- g1.pos/(g1.pos+g1.neg)
    
    # Check if RD_pos and RD_neg are consistent with RD_all in the simulation
    # when using p_Mplus.g0.hat as the (possibly biased) estimate for P(M+)
    check.RD_all.hat[i] <- RD_pos.loss.sim[i]*p_Mplus.g0.hat[i] + 
                           RD_neg.loss.sim[i]*(1-p_Mplus.g0.hat[i])
    
    ##
    # Calculate corrected RD_neg using prior odds and assuming randomization worked
    # so can substitute P(M+|G=1) for P(M+|G=0)
    # But this requires no loss of signal, so it cannot generalize
    ##
    
    # prior odds = P(D+|G=0)/P(D-|G=0) = #(D+,G=0)/#(D-,G=0) 
    # = {#(D+,M-,G=0)+#(D+,M+,G=0)} / {#(M-,G=0)+#(M+,G=0) - [#(D+,M-,G=0)+#(D+,M+,G=0)]}
    Dplus_g0 <- Dplus_g0.neg + Dplus_g0.pos
    Dminus_g0 <- g0.neg + g0.pos - Dplus_g0
    odds.prior <- Dplus_g0 / Dminus_g0
    
    # Calculate P(M+|D+,G=0) = #(D+,M+,G=0)/#(D+,G=0), need to weight up numerator
    p_Mplus_Dplus.g0.weighted <- Dplus_g0.pos.loss.weighted / Dplus_g0
    
    # Calculate P(M+|D-,G=0) = #(D-,M+,G=0)/#(D-,G=0)
    # Dminus_g0.pos <- g0.pos - Dplus_g0.pos
    p_Mplus_Dminus.g0.weighted <- (Dminus_g0.pos.loss*1/frac.g0.minus.hat) / Dminus_g0

    # # For equal P(loss|M+,D+,G=0)=P(loss|M+,D-,G=0), the correction factor 
    # # is P(M+|G=1)/P(M+|G=0) ???
    # correction <- (g1.pos/(g1.pos+g1.neg)) / p_Mplus.g0.hat[i]
    # 
    # # BF = {1-P(M+|D+,G=0)} / {1-P(M+|D-,G=0)} 
    # BF_neg <- (1 - correction * p_Mplus_Dplus.g0.weighted) / (1 - correction * p_Mplus_Dminus.g0.weighted)
    # 
    # # Calculate corrected P(D+|M-,G=0)
    # pDplus_g0.neg <- 1/(1 + (odds.prior * BF_neg)^-1)
    
    # Assume randomization worked, so can substitute P(M+|G=1) for P(M+|G=0)
    # 1st calculate P(D+|M+,G=0) = P(M+|D+,G=0)P(D+|G=0)/P(M+|G=1)
    g0 <- g0.pos + g0.neg # #(G=0)
    p_Mplus.g1 <- g1.pos/(g1.pos+g1.neg)
    pDplus_g0.pos = p_Mplus_Dplus.g0.weighted * (Dplus_g0/g0) / p_Mplus.g1
    
    # 2nd calculate P(D+|M-,G=0) = {1-P(M+|D+,G=0)}P(D+|G=0)/{1-P(M+|G=1)}
    pDplus_g0.neg = (1-p_Mplus_Dplus.g0.weighted) * (Dplus_g0/g0) / (1-p_Mplus.g1)
    
    # Finally calculate RD_pos and RD_neg valid if randomization worked
    RD_pos.BF[i] <- pDplus_g0.pos - Dplus_g1.pos / g1.pos
    RD_neg.BF[i] <- pDplus_g0.neg - Dplus_g1.neg / g1.neg
    
    
    ##
    # Simulate validation study in screened arm (G=1) that retests all M+ 
    # for loss of signal
    ##
    
    # Retest in G=1 all M+, and loss of signal depends on D+ vs. D- 
    Dplus_g1.pos.retest <- rbinom(1,Dplus_g1.pos, frac.g0.plus.pos)
    Dminus_g1.pos.retest <- rbinom(1, g1.pos - Dplus_g1.pos, frac.g0.minus.pos)
    
    # Estimate retention of signal (i.e. 1-loss): the 2 retest-positive fractions
    # P(~M+|M+,D+,G=1)
    frac.g1.plus.pos.hat <- Dplus_g1.pos.retest / Dplus_g1.pos
    # P(~M+|M+,D-,G=1)
    frac.g1.minus.pos.hat <- Dminus_g1.pos.retest / (g1.pos - Dplus_g1.pos)
    
    # Now recalculate BF_neg, BF_pos, P(D+|M+/-,G=0), and RD_pos and RD_neg
    BF_neg <- (1 - p_Mplus_Dplus.g0.weighted*1/frac.g1.plus.pos.hat) / 
              (1 - p_Mplus_Dminus.g0.weighted*1/frac.g1.minus.pos.hat)
    BF_pos <- (p_Mplus_Dplus.g0.weighted*1/frac.g1.plus.pos.hat) / 
              (p_Mplus_Dminus.g0.weighted*1/frac.g1.minus.pos.hat)
    pDplus_g0.neg <- 1/(1 + (odds.prior * BF_neg)^-1)
    pDplus_g0.pos <- 1/(1 + (odds.prior * BF_pos)^-1)
    RD_neg.BF.retest[i] <- pDplus_g0.neg - Dplus_g1.neg / g1.neg
    RD_pos.BF.retest[i] <- pDplus_g0.pos - Dplus_g1.pos / g1.pos
    RR_neg.BF.retest[i] <- (Dplus_g1.neg / g1.neg) / pDplus_g0.neg
    RR_pos.BF.retest[i] <- (Dplus_g1.pos / g1.pos) / pDplus_g0.pos
    
    # Calculate estimated P(M+|G)
    p_Mplus.g1.hat[i] <- g1.pos/(g1.pos+g1.neg) # G=1 is easy because no loss of signal
    # G=0 has to account for loss of signal:
    # P(M+|G=0) = P(M+|D+,G=0)P(D+|G=0) + P(M+|D-,G=0)P(D-|G=0)
    p_Mplus.g0.hat[i] <- ( (p_Mplus_Dplus.g0.weighted*1/frac.g1.plus.pos.hat)*Dplus_g0 + 
                           (p_Mplus_Dminus.g0.weighted*1/frac.g1.minus.pos.hat)*Dminus_g0 ) /
                         (Dplus_g0+Dminus_g0)
    
    # Use retest-positive fractions directly to correct P(D+|~M+,G=0) = p_control.pos
    # so it estimates P(D+|M+,G=0). See 03_Intended-Effect-Stat.docx
    # First calculate P(~M+|M+,G=1) directly
    frac.g1.pos.hat <- (Dplus_g1.pos.retest + Dminus_g1.pos.retest) / g1.pos
    # Now calculate P(D+|M+,G=0) = P(D+|~M+,G=0) * {P(~M+|M+,G=1)/P(~M+|M+,D+,G=1)}
    # and plug into RDpos
    RD_pos.retest[i] <- p_control.pos*(frac.g1.pos.hat/frac.g1.plus.pos.hat) - Dplus_g1.pos / g1.pos
    
    
    
    ##
    # Now simulate loss of MCED signal in the control group: signal strength available
    # Start with the subsampled group
    ##
    
    # fraction of signal that is strong vs. weak
    frac.ss <- 0.25
    frac.ws <- 1-frac.ss

    ### First do loss-of-signal among D+
    
    # Subsample among a_0=#(D+,M+,G=0): the fraction of these that test true positive M+
    # Split this into numbers that are true-strong-signal vs true-weak-signal
    Dplus_g0.pos.ss <- rbinom(1,Dplus_g0.pos.subsamp, frac.ss)
    Dplus_g0.pos.ws <- Dplus_g0.pos.subsamp - Dplus_g0.pos.ss
    
    # Now among these do loss-of-signal
    # Strong-signal can degrade to weak or no signal; Weak-signal can degrade to no signal
    frac.ss.to.ws <- 0 #0.2
    frac.ss.to.neg <- 0.1
    Dplus_g0.pos.ss.to.ws <- rbinom(1,Dplus_g0.pos.ss, frac.ss.to.ws)
    Dplus_g0.pos.ss.to.neg <- rbinom(1,Dplus_g0.pos.ss, frac.ss.to.neg)
    frac.ws.to.neg <- 0.3
    Dplus_g0.pos.ws.to.neg <- rbinom(1,Dplus_g0.pos.ws, frac.ws.to.neg)
    
    # Total observed a_0=#(D+,M+,G=0) for strong signals and for weak signals
    Dplus_g0.pos.ss.loss <- Dplus_g0.pos.ss - (Dplus_g0.pos.ss.to.ws + Dplus_g0.pos.ss.to.neg)
    Dplus_g0.pos.ws.loss <- Dplus_g0.pos.ws + Dplus_g0.pos.ss.to.ws - Dplus_g0.pos.ws.to.neg 
    Dplus_g0.pos.loss <- Dplus_g0.pos.ss.loss + Dplus_g0.pos.ws.loss
    
    # The a_0=#(D+,M+,G=0) that test false-negative, put them into b_0=#(D+,M-,G=0)
    Dplus_g0.neg.loss <- Dplus_g0.neg.subsamp + (Dplus_g0.pos.ss.to.neg - Dplus_g0.pos.ws.to.neg)
    
    
    ### Second do loss-of-signal among D-
    
    # Subsample among c_0=#(D-,M+,G=0): the fraction of these that test true positive M+
    # Split this into numbers that are true-strong-signal vs true-weak-signal
    Dminus_g0.pos.ss <- rbinom(1,Dminus_g0.pos.subsamp, frac.ss)
    Dminus_g0.pos.ws <- Dminus_g0.pos.subsamp - Dminus_g0.pos.ss
    
    # Now among these do loss-of-signal
    # Strong-signal can degrade to weak or no signal; Weak-signal can degrade to no signal
    frac.ss.to.ws <- 0 #0.2
    frac.ss.to.neg <- 0.1
    Dminus_g0.pos.ss.to.ws <- rbinom(1,Dminus_g0.pos.ss, frac.ss.to.ws)
    Dminus_g0.pos.ss.to.neg <- rbinom(1,Dminus_g0.pos.ss, frac.ss.to.neg)
    frac.ws.to.neg <- 0.3
    Dminus_g0.pos.ws.to.neg <- rbinom(1,Dminus_g0.pos.ws, frac.ws.to.neg)
    
    # Total observed a_0=#(D+,M+,G=0) for strong signals and for weak signals
    Dminus_g0.pos.ss.loss <- Dminus_g0.pos.ss - (Dminus_g0.pos.ss.to.ws + Dminus_g0.pos.ss.to.neg)
    Dminus_g0.pos.ws.loss <- Dminus_g0.pos.ws + Dminus_g0.pos.ss.to.ws - Dminus_g0.pos.ws.to.neg 
    Dminus_g0.pos.loss <- Dminus_g0.pos.ss.loss + Dminus_g0.pos.ws.loss
    
    # The c_0=#(D-,M+,G=0) that test false-negative, put them into d_0=#(D-,M-,G=0)
    Dminus_g0.neg.loss <- Dminus_g0.neg.subsamp + (Dminus_g0.pos.ss.to.neg - Dminus_g0.pos.ws.to.neg)
    
    
    ### Third do the analysis, without correcting for loss-of-signal 
    
    # Estimate sampling fractions for D+ and for D- within arm G=0
    frac.g0.plus.hat <- (Dplus_g0.pos.loss + Dplus_g0.neg.loss) /
      (Dplus_g0.pos + Dplus_g0.neg)
    frac.g0.minus.hat<- (Dminus_g0.pos.loss + Dminus_g0.neg.loss) /
      ((g0.pos - Dplus_g0.pos) + (g0.neg - Dplus_g0.neg))
    
    # Do for M- table: RD_neg = p_control - p_screen, account for subsampling in control arm
    Dplus_g0.neg.loss.weighted <- Dplus_g0.neg.loss*1/frac.g0.plus.hat
    g0.neg.weighted <- Dplus_g0.neg.loss.weighted + 
      Dminus_g0.neg.loss*1/frac.g0.minus.hat
    p_screen.neg  <- Dplus_g1.neg / g1.neg  # P(D+|M-,G=1)
    p_control.neg <- Dplus_g0.neg.loss.weighted / g0.neg.weighted  # P(D+|M-,G=0)
    RD_neg.loss.sim[i] <- p_control.neg - p_screen.neg
    RR_neg.loss.sim[i] <- p_screen.neg / p_control.neg
    
    # Do for M+ table: RD_pos = p_control - p_screen, account for subsampling in control arm
    Dplus_g0.pos.loss.weighted <- Dplus_g0.pos.loss*1/frac.g0.plus.hat
    g0.pos.weighted <- Dplus_g0.pos.loss.weighted + 
      Dminus_g0.pos.loss*1/frac.g0.minus.hat
    p_screen.pos  <- Dplus_g1.pos / g1.pos  # P(D+|M+,G=1)
    p_control.pos <- Dplus_g0.pos.loss.weighted / g0.pos.weighted  # P(D+|M+,G=0)
    RD_pos.loss.ss.sim[i] <- p_control.pos - p_screen.pos
    RR_pos.loss.ss.sim[i] <- p_screen.pos / p_control.pos
    
    
  } # END SIMULATION

  
  
  # Pack up output for return
  out <- data.frame(RR=RR, RR_pos=RR_pos, RR_neg=RR_neg, RD_neg=RD_neg, RD_pos=RD_pos, RD=RD,
         RD_neg.sim=mean(RD_neg.sim), RD_pos.sim=mean(RD_pos.sim), RD.sim=mean(RD_all.sim),
         Z_ratio=Z_ratio,
         Z_ratio.sim = mean(Z_pos)/mean(Z_all),
         # power.neg=mean(p.val_neg <= .05), # power for Fisher exact test not Wald test
         # power.pos=mean(p.val_pos <= .05),
         # power.all=mean(p.val_all <= .05),
         power.RDneg=mean(2*pnorm(-abs(RD_neg.sim)/sd(RD_neg.sim)) <= 0.05),
         power.RDpos=mean(2*pnorm(-abs(RD_pos.sim)/sd(RD_pos.sim)) <= 0.05),
         power.RDall=mean(2*pnorm(-abs(RD_all.sim)/sd(RD_all.sim)) <= 0.05),
         power.RRneg=mean(2*pnorm(-abs(log(RR_neg.sim))/sd(log(RR_neg.sim))) <= 0.05),
         power.RRpos=mean(2*pnorm(-abs(log(RR_pos.sim))/sd(log(RR_pos.sim))) <= 0.05),
         power.RRall=mean(2*pnorm(-abs(log(RR_all.sim))/sd(log(RR_all.sim))) <= 0.05),
         # p.neg=median(p.val_neg),
         # p.pos= median(p.val_pos),
         # p.all= median(p.val_all),
         # p.ratio=median(p.val_all)/median(p.val_pos),
         p.pos.lt.p.all=sum(p.val_pos<p.val_all)/nsim,
         RD_neg.subsamp=mean(RD_neg.subsamp.sim), 
         RD_pos.subsamp=mean(RD_pos.subsamp.sim), 
         # power.neg.subsamp=mean(p.val_neg.subsamp <= .05),
         # power.pos.subsamp=mean(p.val_pos.subsamp <= .05),
         # p.neg.subsamp=median(p.val_neg.subsamp),
         # p.pos.subsamp=median(p.val_pos.subsamp),
         # p.ratio.subsamp=median(p.val_all)/median(p.val_pos.subsamp),
         # p.pos.gt.p.all.subsamp=sum(p.val_pos.subsamp>p.val_all)/nsim,
         sd.RD_neg = sd(RD_neg.sim),
         sd.RD_pos = sd(RD_pos.sim),
         sd.RD_pos.subsamp = sd(RD_pos.subsamp.sim),
         sd.ratio = sd(RD_pos.subsamp.sim)/sd(RD_pos.sim),
         sd.RD_neg.subsamp = sd(RD_neg.subsamp.sim),
         # p.sim=median(2*pnorm(-abs(RD_all.sim)/sd(RD_all.sim))), # RD_all p-val agrees with simulation
         # p.pos.sim=median(2*pnorm(-abs(RD_pos.sim)/sd(RD_pos.sim))),# RD_pos p-val agrees with simulation
         # p.pos.subsamp.sim = median(2*pnorm(-abs(RD_pos.subsamp.sim)/sd(RD_pos.subsamp.sim))),
         power.pos.subsamp = mean(2*pnorm(-abs(RD_pos.subsamp.sim)/sd(RD_pos.subsamp.sim)) <= 0.05),
         # p.neg.subsamp.sim = median(2*pnorm(-abs(RD_neg.subsamp.sim)/sd(RD_neg.subsamp.sim))),
         power.neg.subsamp = mean(2*pnorm(-abs(RD_neg.subsamp.sim)/sd(RD_neg.subsamp.sim)) <= 0.05),
         RD_pos.loss = mean(RD_pos.loss.sim),
         RD_neg.loss = mean(RD_neg.loss.sim),
         sd.RD_pos.loss = sd(RD_pos.loss.sim),
         power.RD_pos.loss = mean(2*pnorm(-abs(RD_pos.loss.sim)/sd(RD_pos.loss.sim)) <= 0.05),
         power.RD_neg.loss = mean(2*pnorm(-abs(RD_neg.loss.sim)/sd(RD_neg.loss.sim)) <= 0.05),
         RR_pos.loss = mean(RR_pos.loss.sim),
         RR_neg.loss = mean(RR_neg.loss.sim),
         sd.RR_pos.loss = sd(RR_pos.loss.sim),
         power.RR_pos.loss = mean(2*pnorm(-abs(log(RR_pos.loss.sim))/sd(log(RR_pos.loss.sim))) <= 0.05),
         power.RR_neg.loss = mean(2*pnorm(-abs(log(RR_neg.loss.sim))/sd(log(RR_neg.loss.sim))) <= 0.05),
         p_Mplus.g0.hat = mean(p_Mplus.g0.hat),
         p_Mplus.g1.hat = mean(p_Mplus.g1.hat),
         check.RD_all.hat = mean(check.RD_all.hat),
         # RD_neg.BF = mean(RD_neg.BF),         
         # sd.RD_neg.BF = sd(RD_neg.BF),
         # power.neg.BF = mean(2*pnorm(-abs(RD_neg.BF)/sd(RD_neg.BF)) <= 0.05),
         # RD_pos.BF = mean(RD_pos.BF),         
         # sd.RD_pos.BF = sd(RD_pos.BF),
         # power.pos.BF = mean(2*pnorm(-abs(RD_pos.BF)/sd(RD_pos.BF)) <= 0.05),
         RD_neg.BF.retest=mean(RD_neg.BF.retest),
         RD_pos.BF.retest=mean(RD_pos.BF.retest),
         sd.RD_neg.BF.retest = sd(RD_neg.BF.retest),
         sd.RD_pos.BF.retest = sd(RD_pos.BF.retest),
         power.RD_neg.BF.retest=mean(2*pnorm(-abs(RD_neg.BF.retest)/sd(RD_neg.BF.retest)) <= 0.05),
         power.RD_pos.BF.retest=mean(2*pnorm(-abs(RD_pos.BF.retest)/sd(RD_pos.BF.retest)) <= 0.05),
         RR_neg.BF.retest=mean(RR_neg.BF.retest),
         RR_pos.BF.retest=mean(RR_pos.BF.retest),
         sd.RR_neg.BF.retest = sd(RR_neg.BF.retest),
         sd.RR_pos.BF.retest = sd(RR_pos.BF.retest),
         power.RR_neg.BF.retest=mean(2*pnorm(-abs(log(RR_neg.BF.retest))/sd(log(RR_neg.BF.retest))) <= 0.05),
         power.RR_pos.BF.retest=mean(2*pnorm(-abs(log(RR_pos.BF.retest))/sd(log(RR_pos.BF.retest))) <= 0.05),
         RD_pos_IE = mean(RD_pos_IE.sim),
         sd.RD_pos_IE = sd(RD_pos_IE.sim),
         power.RD_pos_IE = mean(2*pnorm(-abs(RD_pos_IE.sim)/sd(RD_pos_IE.sim)) <= 0.05),
         RD_IE = mean(RD_IE.sim),
         sd.RD_IE = sd(RD_IE.sim),
         sd.RD_all = sd(RD_all.sim),
         # RR_pos_IE = mean(RR_pos_IE.sim),
         # sd.RR_pos_IE = sd(RR_pos_IE.sim),
         # power.RR_pos_IE = mean(2*pnorm(-abs(log(RR_pos_IE.sim))/sd(log(RR_pos_IE.sim))) <= 0.05),
         medianp.RD_pos = median(2*pnorm(-abs(RD_pos.sim)/sd(RD_pos.sim))),
         medianp.RD_pos_IE = median(p.val.RD_pos_IE.sim),
         medianp.RD_all = median(2*pnorm(-abs(RD_all.sim)/sd(RD_all.sim))),
         p.RD_pos_IE_lt_p.RD_all = sum(p.val.RD_pos_IE.sim < p.val_all)/nsim,
         p.RD_pos_IE_lt_p.RD_pos = sum(p.val.RD_pos_IE.sim < p.val_pos)/nsim,
         RR_all = mean(RR_all.sim),
         RD_pos.retest = mean(RD_pos.retest),
         power.RD_pos.retest = mean(2*pnorm(-abs(RD_pos.retest)/sd(RD_pos.retest)) <= 0.05),
         RR_pos.sim = mean(RR_pos.sim),
         RD_pos.loss.ss = mean(RD_pos.loss.ss.sim),
         # RD_neg.loss = mean(RD_neg.loss.sim),
         # sd.RD_pos.loss = sd(RD_pos.loss.sim),
         power.RD_pos.loss.ss = mean(2*pnorm(-abs(RD_pos.loss.ss.sim)/sd(RD_pos.loss.ss.sim)) <= 0.05),
         # power.RD_neg.loss = mean(2*pnorm(-abs(RD_neg.loss.sim)/sd(RD_neg.loss.sim)) <= 0.05),
         RR_pos.loss.ss = mean(RR_pos.loss.ss.sim),
         # RR_neg.loss = mean(RR_neg.loss.sim),
         # sd.RR_pos.loss = sd(RR_pos.loss.sim),
         power.RR_pos.loss.ss = mean(2*pnorm(-abs(log(RR_pos.loss.ss.sim))/sd(log(RR_pos.loss.ss.sim))) <= 0.05),
         # power.RR_neg.loss = mean(2*pnorm(-abs(log(RR_neg.loss.sim))/sd(log(RR_neg.loss.sim))) <= 0.05),
         RR.tar=mean(RR_tar.sim),
         power.RRtar=mean(2*pnorm(-abs(log(RR_tar.sim))/sd(log(RR_tar.sim))) <= 0.05),
         RR_pos.noncomp = mean(RR_pos.noncomp.sim),
         RR_tar.noncomp = mean(RR_tar.noncomp.sim),
         power.RRpos.noncomp = mean(2*pnorm(-abs(log(RR_pos.noncomp.sim))/sd(log(RR_pos.noncomp.sim))) <= 0.05),
         power.RRtar.noncomp = mean(2*pnorm(-abs(log(RR_tar.noncomp.sim))/sd(log(RR_tar.noncomp.sim))) <= 0.05),
         RR_neg.noncomp = mean(RR_neg.noncomp.sim),
         power.RRneg.noncomp = mean(2*pnorm(-abs(log(RR_neg.noncomp.sim))/sd(log(RR_neg.noncomp.sim))) <= 0.05),
                pD_g1.pos=pD_g1.pos,
                pD_g0.pos=pD_g0.pos,
                pD_g1.neg=pD_g1.neg,
                pD_g0.neg=pD_g0.neg,
                pMplus_Dminus=pMplus_Dminus,
                pMplus_Dplus=pMplus_Dplus,
                pD_pos=pD_pos,
         pD_pos_check,
                pD_neg=pD_neg,
                pD=pD
  )

  # print(cbind(RD_pos_IE.sim,RD_pos.sim, p.val_pos, p.val.RD_pos_IE.sim, p.val_all))

  return(out)

  # return(list(RD_neg.sim=RD_neg.sim,RD_pos.sim=RD_pos.sim,RD_all.sim=RD_all.sim,
  #              p.val_neg=p.val_neg,p.val_pos=p.val_pos,p.val_all=p.val_all,
  #              Z_neg=Z_neg,Z_pos=Z_pos,Z_all=Z_all)
  #        )
}

##
# Compare IE to targeted IE
##

# JNCI example, 20% vs 25% differential non-compliance in the ppt file Tiger-Team-6Sep2022.ppt
conceal.reveal(nsim=1, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos=0.86667,
               frac.g0.plus = 1, frac.g0.minus = 1,
               censor.screen.Dplus = 0.2,
               censor.control.Dplus = 0.2,
               censor.screen.Dminus = 0.2,
               censor.control.Dminus = 0.2,
               move.screen.Mplus.Dplus = 0,
               move.screen.Mminus.Dplus = 0,
               sim=FALSE)

# JNCI example, strongly differential non-compliance in the ppt file Tiger-Team-6Sep2022.ppt
conceal.reveal(nsim=60e3, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos=0.86667,
               frac.g0.plus = 1, frac.g0.minus = 1,
               censor.screen.Dplus = 0.4,
               censor.control.Dplus = 0.8,
               censor.screen.Dminus = 0.8,
               censor.control.Dminus = 0.4,
               move.screen.Mplus.Dplus = 0,
               move.screen.Mminus.Dplus = 0,
               sim=TRUE)

# Redo JNCI example with IE vs Targeted
conceal.reveal(nsim=10e3, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos=0.86667,
               frac.g0.plus = 1,
               frac.g0.minus = 1,
               sim=TRUE)



# JNCI example, strongly differential non-compliance *but null both overall and in ever/never-pos*
conceal.reveal(nsim=30e3, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=1,RR_neg=1, RR_pos=1,
               frac.g0.plus = 1, frac.g0.minus = 1,
               censor.screen.Dplus = 0.4,
               censor.control.Dplus = 0.8,
               censor.screen.Dminus = 0.8,
               censor.control.Dminus = 0.4,
               move.screen.Mplus.Dplus = 0,
               move.screen.Mminus.Dplus = 0,
               sim=TRUE)

# Increase P(M+) to 50%, have slightly differential non-compliance
conceal.reveal(nsim=10e3, n=100e3, p_Mplus=0.5, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos=0.86667,
               frac.g0.plus = 1, frac.g0.minus = 1,
               censor.screen.Dplus = 0.1,
               censor.control.Dplus = 0.13,
               censor.screen.Dminus = 0.12,
               censor.control.Dminus = 0.15,
               move.screen.Mplus.Dplus = 0,
               move.screen.Mminus.Dplus = 0,
               sim=TRUE)

# % reduction in test costs for IE vs targeted, 3 screens, $3000/person and $300/test
1-(53e3*(2*2000+ 2*3*200))/(73e3*(2*2000 + 1*3*200))
# % reduction in test costs for IE vs standard, 3 screens, $3000/person and $300/test
1-(53e3*(2*2000+ 2*3*200))/(98e3*(2*2000 + 1*3*200))
# Find the ratio of person/test cost so that IE and targeted are equal cost
# Finds that the ratio is 742.5/300=2.475
uniroot(function(x) 1-(53e3*(2*x+ 2*3*200))/(73e3*(2*x + 1*3*200)), c(0,3000))
# Above but say we only need to test half the control arm
uniroot(function(x) 1-(53e3*(2*x+ 1.5*3*200))/(73e3*(2*x + 1*3*200)), c(0,3000))

# Make sure I can reproduce Paul's figure 2 of IE sample sizes for 90% power in the JNCI paper
# Run this line, take the p-value for the RRpos, in this case p=1.55680431351154e-05
# Plug it into the line below to get the sample size for 90% power for the IE
# Ex: for P(M+)=5%, RRpos=1-(1-RR)/(1-P_Evpos)=0.875, N=55165 per arm, which looks like it matches to the figure
conceal.reveal(nsim=1, n=2*98e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos= 1-0.1/0.8, sim=FALSE)

conceal.reveal(nsim=1, n=2*50e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos= 0.86667, sim=FALSE)



##
# Loss of signal
# Fix n=100k, but now vary differential loss of signal from 0.5 to 1 and 0.4 to 0.9
##
basecase<-conceal.reveal(nsim=1e4,sim=TRUE)
# fracs.plus <-c(seq(0.5,1,by=0.1),1)
# fracs.minus <- c(seq(0.5,1,by=0.1),1.1) - 0.1 # allow for a base-case of 1.0 fraction
fracs.plus <-c(1,0.45,0.9) #c(seq(0.5,1,by=0.1),1)
fracs.minus <- c(1,0.45,0.45) # allow for a base-case of 1.0 fraction
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (i in 1:length(fracs.plus)) {
  # out <- rbind(out,
  #              as.matrix(conceal.reveal(
  #                nsim=10000, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
  #                RR=0.9,RR_neg=1, RR_pos= 0.86667,
  #                frac.g0.plus = 0.95, frac.g0.minus = 0.5, 
  #                frac.g0.plus.pos = fracs.plus[i], frac.g0.minus.pos = fracs.minus[i], 
  #                sim=TRUE))
  # )
  out <- rbind(out,
               as.matrix(conceal.reveal(
                 nsim=10000, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
                 RR=0.9,RR_neg=1, RR_pos= 0.86667,
                 frac.g0.plus = 1, frac.g0.minus = 1, 
                 frac.g0.plus.pos = fracs.plus[i], frac.g0.minus.pos = fracs.minus[i], 
                 sim=TRUE))
  )
}
print(out[-1,])
print(basecase <- as.data.frame(cbind(fracs.plus,fracs.minus,
                                      out[-1,4+c(24,25,27,28,29,30,32,33)])), digits=2, row.names=FALSE)
print(basecase <- as.data.frame(cbind(fracs.plus,fracs.minus,
                                      out[-1,4+c(36,37,38,40,41,42,43,46,47,65-4,66-4)])), digits=2, row.names=FALSE)




##
# Table 3: effect of non-adherence to blood draws (dropout)
# Compare IE to Targeted
##
# Simulate and save the data
basecase <- conceal.reveal(nsim=1,sim=TRUE)
# n <- c(25e3,50e3,75e3,100e3,150e3,200e3) #seq(120e3,150e3,by=10e3)
# RR_poss <- c(0.86667,0.8)
# RR_negs <- c(1,1.05,0.95)
# RR <- 0.9
n <- 100e3 #c(10e3, 15e3, 20e3, 25e3,30e3,40e3,50e3) #seq(120e3,150e3,by=10e3)
RR_poss <-0.86667 #c(0.8,0.86667)
RR_negs <- 1 #c(1,1.05,0.95)
RR <- 0.9
censor <- matrix(c(0.00,	0.00,	0.00,	0.00, 0.00, 0.00,
                   0.30,	0.30,	0.30,	0.30, 0.00, 0.00,
                   0.20,	0.20,	0.30,	0.30, 0.00, 0.00,
                   0.30,	0.10,	0.30,	0.10, 0.00, 0.00,
                   0.30,	0.15,	0.15,	0.30, 0.00, 0.00,
                   0.15,	0.30,	0.30,	0.15, 0.00, 0.00),ncol=6,byrow=TRUE)
colnames(censor)<- c("screen D+","screen D-","control D+","control D-", "screen M+ to D-", "screen M- to D-")
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (l in 1:nrow(censor)) {
  for (k in 1:length(RR_negs)) {
    for (j in 1:length(RR_poss)) {
      for (i in 1:length(n)) {
        out <- rbind(out,
                     as.matrix(conceal.reveal(nsim=2e4, n=n[i], p_Mplus=0.05, p_g0=0.5, pD_g0=0.02,
                                              RR=RR,RR_neg=RR_negs[k], RR_pos=RR_poss[j],
                                              frac.g0.plus = 1, 
                                              frac.g0.minus = 1,
                                              frac.g0.plus.pos = 1, 
                                              frac.g0.minus.pos = 1,
                                              censor.screen.Dplus = censor[l,1],
                                              censor.control.Dplus = censor[l,3],
                                              censor.screen.Dminus = censor[l,2],
                                              censor.control.Dminus = censor[l,4],
                                              move.screen.Mplus.Dplus = censor[l,5],
                                              move.screen.Mminus.Dplus = censor[l,6],
                                              sim=TRUE)) )
      }}}}
# basecase <- as.data.frame(cbind(censor, out[-1,4+c(38,63,43,60,48,34,35)]))
# names(basecase)<-c(colnames(censor),"RDpos","RRpos","RRneg","RR","power RRpos","P(M+|G=0)","P(M+|G=1)")
basecase <- as.data.frame(cbind(censor, out[-1,4+c(63,43,60,48,68,69,70,71,74,72,73,75,34,35)]))
names(basecase)<-c(colnames(censor),"RRpos","RRneg","RR","power RRpos","RRtar","power RRtar","RRposFix",
                   "RRtarFix","RRnegFix","power RRposFix","power RRtarFix","power RRnegFix","P(M+|G=0)","P(M+|G=1)")
print(t(basecase), digits=3, row.names=FALSE)


# censor <- matrix(c(0.00,	0.00,	0.00,	0.00, 0.00, 0.00,
#                    0.30,	0.30,	0.30,	0.30, 0.00, 0.00,
#                    0.30,	0.30,	0.30,	0.30, 1e-3, 1e-3,
#                    0.20,	0.20,	0.30,	0.30, 0.00, 0.00,
#                    0.20,	0.20,	0.30,	0.30, 1e-3, 1e-3,
#                    0.30,	0.10,	0.30,	0.10, 0.00, 0.00,
#                    0.30,	0.10,	0.30,	0.10, 1e-3, 1e-3,
#                    0.10,	0.05,	0.05,	0.10, 0.00, 0.00,
#                    0.10,	0.05,	0.05,	0.10, 1e-3, 1e-3,
#                    0.05,	0.10,	0.10,	0.05, 0.00, 0.00,
#                    0.05,	0.10,	0.10,	0.05, 1e-3, 1e-3),ncol=6,byrow=TRUE)
# censor <- matrix(c(0.30,	0.30,	0.30,	0.30, 0.00, 0.00,
#                    0.30,	0.10,	0.30,	0.10, 0.00, 0.00),ncol=6,byrow=TRUE)

censor <- matrix(c(0.00,	0.00,	0.00,	0.00, 0.00, 0.00,
                   0.10,	0.10,	0.10,	0.10, 0.00, 0.00,
                   0.10,	0.10,	0.12,	0.12, 0.00, 0.00,
                   0.10,	0.10,	0.13,	0.13, 0.00, 0.00,
                   0.20,	0.20,	0.20,	0.20, 0.00, 0.00,
                   0.20,	0.20,	0.21,	0.21, 0.00, 0.00,
                   0.20,	0.20,	0.22,	0.22, 0.00, 0.00,
                   0.20,	0.20,	0.23,	0.23, 0.00, 0.00),ncol=6,byrow=TRUE)


##
# Table 2: effect of loss-of-signal in stored bloods in control arm
##
# Simulate and save the data
basecase <- conceal.reveal(nsim=1,sim=TRUE)
# n <- c(25e3,50e3,75e3,100e3,150e3,200e3) #seq(120e3,150e3,by=10e3)
# RR_poss <- c(0.86667,0.8)
# RR_negs <- c(1,1.05,0.95)
# RR <- 0.9
n <- 100e3 #c(10e3, 15e3, 20e3, 25e3,30e3,40e3,50e3) #seq(120e3,150e3,by=10e3)
RR_poss <-0.86667 #c(0.8,0.86667)
RR_negs <- 1 #c(1,1.05,0.95)
RR <- 0.9
censor <- matrix(c(0.50,	0.20,	0.00,	0.00,
                   0.50,  0.20, 0.50, 0.20),ncol=4,byrow=TRUE)
colnames(censor)<- c("screen D+","screen D-","control D+","control D-")
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (l in 1:nrow(censor)) {
  for (k in 1:length(RR_negs)) {
    for (j in 1:length(RR_poss)) {
      for (i in 1:length(n)) {
        out <- rbind(out,
                     as.matrix(conceal.reveal(nsim=1e4, n=n[i], p_Mplus=0.05, p_g0=0.5, pD_g0=0.02,
                                              RR=RR,RR_neg=RR_negs[k], RR_pos=RR_poss[j],
                                              # frac.g0.plus = 0.95, 
                                              # frac.g0.minus = 0.5,
                                              # frac.g0.plus.pos = 0.9, 
                                              # frac.g0.minus.pos = 0.8,
                                              frac.g0.plus = 1, 
                                              frac.g0.minus = 1,
                                              frac.g0.plus.pos = 0.9, 
                                              frac.g0.minus.pos = 0.8,
                                              censor.screen.Dplus = 0,
                                              censor.control.Dplus = 0,
                                              censor.screen.Dminus = 0,
                                              censor.control.Dminus = 0,
                                              move.screen.Mplus.Dplus = 0,
                                              move.screen.Mminus.Dplus = 0,
                                              sim=TRUE)) )
      }}}}
basecase <- as.data.frame(cbind(censor, out[-1,4+c(38,44,43,48,34,35)]))
names(basecase)<-c(colnames(censor),"RDpos","RRpos","RRneg","power RRpos","P(M+|G=0)","P(M+|G=1)")
print(basecase, digits=3, row.names=FALSE)



# Print 1 simulation of 2x2 tables with dropout
conceal.reveal(nsim=1, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02,
               RR=0.9,RR_neg=1, RR_pos=0.86667,
               frac.g0.plus = 0.95,
               frac.g0.minus = 0.5,
               frac.g0.plus.pos = 0.9,
               frac.g0.minus.pos = 0.8,
               censor.screen.plus = 0.5,
               censor.control.plus = 0,
               censor.screen.minus = 0.2,
               censor.control.minus = 0,
               sim=TRUE)





##
# Figure 2, Stat Paper, power vs RRpos:
#
# Calculate Z_ratios (ratio increase in the noncentrality param) varying RRpos and P(M+)
# fix RRneg=1, RR, P(D+|G=0), and P(G=0)=0.5
basecase <- conceal.reveal(nsim=1,sim=TRUE)
# Note that P(M+)=0.025, you get a non-monotonic Z-ratio at RRpos=RR
# RR <- 0.8
# n <- 50e3 #c(10e3, 15e3, 20e3, 25e3,30e3,40e3,50e3)
# RR_poss <- c(0.8,0.7,0.6,0.3,0)
# p_Mpluss <- c(0.025,0.05,0.10, 0.25, 0.50, 0.75, 0.9)
RR <- 0.9
n <- 50e3 #c(10e3, 15e3, 20e3, 25e3,30e3,40e3,50e3)
RR_poss <- seq(0.9,0.55,by=-0.05) #c(0.9,0.8,0.7,0.5)
p_Mpluss <- c(0.025,0.05,0.50,0.90)
RR_negs <- c(1,1.05,0.95)
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (l in 1:length(RR_negs)) {
  for (k in 1:length(p_Mpluss)) {
    for (j in 1:length(RR_poss)) {
      for (i in 1:length(n)) {
        out <- rbind(out,
                     as.matrix(conceal.reveal(nsim=10000, n=n[i], p_Mplus=p_Mpluss[k], p_g0=0.5, pD_g0=0.02, 
                                              RR=RR,RR_neg=RR_negs[l], RR_pos=RR_poss[j],
                                              sim=TRUE))
        )
      }
    }
  }
}
out.fig1 <- as.data.frame(cbind(n,out[-1,c(2,10:11,13:14,73)],
                                rep(p_Mpluss,each=length(n)*length(RR_poss)),
                                rep(RR_negs,each=length(n)*length(RR_poss)*length(p_Mpluss))
))
names(out.fig1)<-c("n","RR_pos","Zratio","Zratio_sim","IEpower","Stdpower","Targetedpower","everpos","RR_neg")
# Replace everpos=0.9 with Targeted since they are essentially identical
out.fig1[out.fig1$everpos==0.9,"IEpower"] <- out.fig1[out.fig1$everpos==0.9,"Targetedpower"]
out.fig1[out.fig1$everpos==0.9,"everpos"] <- "Targeted"
print(out.fig1, row.names=FALSE)

# Helper function to make plots
plot.power.fig <- function(out.fig1) {
  
  # Don't think I need this
  # plot_data <- reshape2::melt(out.fig1,id.var=c("everpos","RR_pos","RR_neg"))
  # colnames(plot_data)[4:5] <- c("analysis","power")
  # plot_data$power[plot_data$power==1] <- 0.9999 # cannot plot a power of 1 on probit scale
  # print(plot_data)
  
  plot_data <- out.fig1
  
  # The color-blind palette with black as #3, move yellow to end:
  cbbPalette <- c("#0072B2","#D55E00","#000000","#009E73","#56B4E9","#E69F00","#CC79A7", 
                  "#F0E442")
  
  plotpowers <- 
    ggplot(plot_data, aes(x=RR_pos,y=IEpower,color=as.factor(everpos))) +
    scale_color_manual(values=cbbPalette) +
    #scale_colour_manual(values=c("red","black","blue")) +
    geom_point(size=3)+
    geom_line()+#aes(lty=as.factor(everpos))) +
    # scale_linetype_manual(values=c("dotted", "solid"))+ # rev so solid is the usual analysis
    # scale_color_discrete(na.translate = F) +
    scale_x_continuous(breaks=unique(out.fig1$RR_pos)) +
    # scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1)) +
    scale_y_continuous(trans = "probit", limits = c(0.1,0.9999), 
                       breaks=c(seq(0.2,0.9,0.1),0.95,0.98,0.99,0.999)) +
    geom_hline(yintercept = c(0.8,0.9),lty=3) +
    geom_hline(yintercept = c(0.36),lty=1) + # std analysis has power 36%+
    labs(x="RR_pos",y="power",color="ever positivity") +
    #theme with white background
    theme_bw() +
    #eliminates background, gridlines, #and chart border
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      # panel.border = element_blank()
    ) +
    NULL
}

# plot the data
plotpowers.rrneg.1 <- plot.power.fig(out.fig1[out.fig1$RR_neg==1,])
plotpowers.rrneg.2 <- plot.power.fig(out.fig1[out.fig1$RR_neg==1.05,])
plotpowers.rrneg.3 <- plot.power.fig(out.fig1[out.fig1$RR_neg==0.95,])

# Extract legend from one of the plots - they're all the same
legend_b <- cowplot::get_legend(
  plotpowers.rrneg.1 + 
    guides(color = guide_legend(nrow = 1, byrow=TRUE)) +
    theme(legend.position = "bottom")
)

# Plot single analysis and subgroup analysis
prow <- cowplot::plot_grid(plotpowers.rrneg.1 + theme(legend.position="none"), 
                                 plotpowers.rrneg.2 + theme(legend.position="none"),
                                 plotpowers.rrneg.3 + theme(legend.position="none"),
                                 nrow=1, ncol=3, rel_heights = c(1, 1,1),
                                 labels = c("RR_neg=1 (RR=0.9)","RR_neg=1.05 (RR=0.9)","RR_neg=0.95 (RR=0.9)"))
# Include legen
final_plot <- cowplot::plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))

# show the plot 
final_plot

# Save the plot
ggsave(filename="Fig-2-stat-paper.pdf", plot=final_plot, device="pdf",
       units="in", width=15, height=8,  dpi=500)




# Try examples, calculate Zratio
conceal.reveal(nsim=1, n=50e3, p_Mplus=0.025, p_g0=0.5, pD_g0=0.02,
               RR=0.9,RR_neg=0.95, RR_pos= 0.85,
               frac.g0.plus = 1, frac.g0.minus = 1,
               frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,sim=TRUE)





# Try examples where control arm gets MCD, so that RR and RRpos are higher
conceal.reveal(nsim=1, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos=0.86667,
               frac.g0.plus = 1,frac.g0.minus = 1,sim=FALSE) # std JNCI example

conceal.reveal(nsim=1, n=100e3, p_Mplus=0.03, p_g0=0.5, pD_g0=0.02, RR=0.95,RR_neg=1, RR_pos=0.5,
               frac.g0.plus = 1,frac.g0.minus = 1,sim=FALSE) #



##
# Calculate 2x2 tables for tiger team presentation
##
temp <- conceal.reveal(nsim=1e5, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
               RR=0.9,RR_neg=1, RR_pos= 0.86667,
               frac.g0.plus = 1, frac.g0.minus = 1, 
               frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,sim=TRUE)
print(temp)
2*(1-pnorm(1*(qnorm(1-0.05/2)+qnorm(temp$power.RDall))))
2*(1-pnorm(1*(qnorm(1-0.05/2)+qnorm(temp$power.RDpos))))
# binomial SD for RD_all
sd_RD_all <- 0.002/(qnorm(1-0.05/2)+qnorm(temp$power.RDall))
# binomial SD for RD_IE
sd_RD_all*(1/temp$Z_ratio)

# For RR_pos=0.86667 and 100k per arm
Mpos <- matrix( c(1300,1500,3700,3500),nrow=2,byrow=TRUE,
                dimnames=list(c("D+","D-"),c("screen","control")))
Mneg <- matrix( c(500,500,94500,94500),nrow=2,byrow=TRUE,
                dimnames=list(c("D+","D-"),c("screen","control")))


conceal.reveal(nsim=1, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
               RR=0.9,RR_neg=1, RR_pos= 0.867,
               frac.g0.plus = 1, frac.g0.minus = 1, 
               frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,sim=TRUE)
# # For RR_pos=0.8
# Mpos <- matrix( c(800,1000,4200,4000),nrow=2,byrow=TRUE,
#                 dimnames=list(c("D+","D-"),c("screen","control")))
# Mneg <- matrix( c(1000,1000,94000,94000),nrow=2,byrow=TRUE,
#                 dimnames=list(c("D+","D-"),c("screen","control")))

# Example with false-reassurance RR_neg=1.21
conceal.reveal(nsim=1, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02,
               RR=0.9,RR_neg=1.21, RR_pos= 0.866667,
               frac.g0.plus = 1, frac.g0.minus = 1,
               frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,sim=TRUE)
Mpos <- matrix( c(1565,1806,3435,3194),nrow=2,byrow=TRUE,
                dimnames=list(c("D+","D-"),c("screen","control")))
Mneg <- matrix( c(235,194,94765,94806),nrow=2,byrow=TRUE,
                dimnames=list(c("D+","D-"),c("screen","control")))

# Do 1/4 the sample size for the example
# Mpos <- Mpos/4
# Mneg <- Mneg/4
print(Mpos)
print(Mneg)
print(all<-Mneg+Mpos)
print(prop.test(t(all),correct=F))
print(fisher.test(t(all)))
prop.test(t(Mpos),correct=F)
fisher.test(t(Mpos))
prop.test(t(Mneg),correct=F)
fisher.test(t(Mneg))




##
# Paul's comment that we can use IE to reduce variance in the usual trial OR
# Not so because IE says OR is not homogeneous groups: ever-pos vs. never-pos
# The interaction means that we can't stratify and combine over the 2 groups
# to estimate a more accurate OR.  We have to estimate separate ORs.
## 
# Stratified (by ever/never-positive) adjusted estimate of common OR
# But OR is not common between those 2 groups: there is interaction
mantelhaen.test(simplify2array(list(Mpos,Mneg)), correct=FALSE)
# Do GLM including interaction
longdata <- data.frame(expand.grid(arm = c("control","screen"),test= c("ever-pos","never-pos")),
                       D = c(1500,1300,500,500), N=c(5000,5000,95000,95000)
)
longdata
# Usual trial compares arms, but this ignores test results
summary(glm(cbind(D,N-D)~arm,family=binomial,data=longdata))
# One can also compare ever-pos vs never-pos, but this ignores arm which determines treatment
summary(glm(cbind(D,N-D)~test,family=binomial,data=longdata))
# We are interested in the arm*test interaction
summary(glm(cbind(D,N-D)~arm*test,family=binomial,data=longdata))
# Within test-positives only
summary(glm(cbind(D,N-D)~arm,family=binomial,data=longdata,subset=(test=="ever-pos")))


# False-reassurance RRneg=1.05, 20% dropout in control arm, 10% dropout in screen arm
# RRpos=0.867 is unbiased (loses power), but RRneg=0.762 is badly biased (should be 1.05)
tables <- conceal.reveal(nsim=1, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02,
               RR=0.9,RR_neg=1, RR_pos= 0.866667,
               frac.g0.plus = 1, frac.g0.minus = 1,
               frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,
               censor.screen.plus = 0.1, censor.control.plus = 0.2,
               censor.screen.minus = 0.1, censor.control.minus = 0.2,
               sim=FALSE)
print(prop.test(t(tables$all),correct=F))
print(fisher.test(t(tables$all)))
prop.test(t(tables$Mpos),correct=F)
fisher.test(t(tables$Mpos))
prop.test(t(tables$Mneg),correct=F)
fisher.test(t(tables$Mneg))

# Example with no false-reassurance, but 20% dropout in control arm D+ and 10% in control D-
tables <- conceal.reveal(nsim=1, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02,
                         RR=0.9,RR_neg=1, RR_pos= 0.866667,
                         frac.g0.plus = 1, frac.g0.minus = 1,
                         frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,
                         censor.screen.plus = 0.15, censor.control.plus = 0.35,
                         censor.screen.minus = 0.1, censor.control.minus = 0.3,
                         sim=FALSE)
print(prop.test(t(tables$all),correct=F))
print(fisher.test(t(tables$all)))
prop.test(t(tables$Mpos),correct=F)
fisher.test(t(tables$Mpos))
prop.test(t(tables$Mneg),correct=F)
fisher.test(t(tables$Mneg))

##
# Example with larger P(M+)
##
tables <- conceal.reveal(nsim=10000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
               RR=0.9,RR_neg=1, RR_pos=0.866667,
               frac.g0.plus = 1, frac.g0.minus = 1, 
               frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,sim=FALSE)
# conceal.reveal(nsim=10000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
#                RR=0.9,RR_neg=1, RR_pos=0.866667,
#                frac.g0.plus = 1, frac.g0.minus = 1, 
#                frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,sim=TRUE)
print(overall<-prop.test(t(tables$all),correct=F))
sd.overall <- diff(overall$conf.int)/(2*1.96)
print(fisher.test(t(tables$all)))
print(Mpos <- prop.test(t(tables$Mpos),correct=F))
sd.Mpos <- diff(Mpos$conf.int)/(2*1.96)
# c(sd.overall, sd.Mpos, 1-sd.Mpos/sd.overall)
c(Z.overall<-diff(overall$estimate)/sd.overall, Z.Mpos<-diff(Mpos$estimate)/sd.Mpos, Z.Mpos/Z.overall)

# Test based on difference of Poisson counts, no rates
# % reduction in the stdev and pure increase in Z-stat
Z.overall.Poisson <- diff(tables$all[1,])/sqrt(sum(tables$all[1,]))
Z.Mpos.Poisson <- diff(tables$Mpos[1,])/sqrt(sum(tables$Mpos[1,]))
c(Z.overall.Poisson, Z.Mpos.Poisson, 1-Z.overall.Poisson/Z.Mpos.Poisson, Z.Mpos.Poisson/Z.overall.Poisson) 

# Test based on Poisson rates
overall.Poisson.rates <- poisson.test(tables$all[1,],apply(tables$all,2,sum))
Mpos.Poisson.rates <- poisson.test(round(tables$Mpos[1,]),apply(round(tables$Mpos),2,sum))
sd.overall.rates <- diff(log(overall.Poisson.rates$conf.int))/(2*1.96)
sd.Mpos.rates <- diff(log(Mpos.Poisson.rates$conf.int))/(2*1.96)
c(Z.overall.rates<-abs(log(overall.Poisson.rates$estimate)/sd.overall.rates), 
  Z.Mpos.rates<-abs(log(Mpos.Poisson.rates$estimate)/sd.Mpos.rates), 
  1-Z.overall.rates/Z.Mpos.rates, Z.Mpos.rates/Z.overall.rates)



# fisher.test(t(tables$Mpos))
# prop.test(t(tables$Mneg),correct=F)
# fisher.test(t(tables$Mneg))




##
# Fix base case parameters, vary RRneg , RRpos, and sample size
# For Tiger Team
##

# Simulate and save the data
basecase <- conceal.reveal(nsim=1,sim=TRUE)
# # Simulation params for RR=0.9
RR <- 0.9
n <- c(25e3,50e3,75e3,100e3,150e3,200e3)
RR_poss <- c(0.86667,0.8)
RR_negs <- c(1,1.05,0.95)
# # Simulation params for RR=0.8
# RR <- 0.8
# n <- c(10e3, 15e3, 20e3, 25e3,30e3,40e3,50e3)
# RR_poss <- c(0.7,0.76667)
# RR_negs <- c(1,1.05,0.95)
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (k in 1:length(RR_negs)) {
  for (j in 1:length(RR_poss)) {
    for (i in 1:length(n)) {
      out <- rbind(out,
         as.matrix(conceal.reveal(nsim=10000, n=n[i], p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
                                  RR=RR,RR_neg=RR_negs[k], RR_pos=RR_poss[j],
                                  frac.g0.plus = 0.95, frac.g0.minus = 0.5, 
                                  frac.g0.plus.pos = 0.9, frac.g0.minus.pos = 0.8, 
                                  censor.screen.plus = 0.05, censor.control.plus = 0.1,
                                  censor.screen.minus = 0.05, censor.control.minus = 0.1,
                                  sim=TRUE))
      )
    }
  }
}
print(out.fig1 <- as.data.frame(cbind(n,out[-1,c(2,3,13:14)])), row.names=FALSE)
names(out.fig1)<-c("n","RR_pos","RR_neg","IE","usual")
save(out.fig1,file="out.fig1.RR0.9.RData")
# save(out.fig1,file="out.fig1.RR0.8.RData")


# Helper function to make plots
plot.power.fig <- function(out.fig1) {
  plot_data <- reshape2::melt(out.fig1,id.var=c("n","RR_pos","RR_neg"))
  colnames(plot_data)[4:5] <- c("analysis","power")
  plot_data$power[plot_data$power==1] <- 0.9999 # cannot plot a power of 1 on probit scale
  
  print(plot_data)
  
  # remove excess rows for usual trial RD_all (keep only the first such rows)
  plot_data <- rbind(plot_data[plot_data$analysis=="IE",],
                     plot_data[plot_data$analysis=="usual",][1:length(unique(plot_data$n)),])
  # set RR_neg to -1 for the usual trial analysis to give it's own color
  plot_data$RR_neg[plot_data$analysis=="usual"] <- "Not applicable"
  plot_data$RR_pos[plot_data$analysis=="usual"] <- "Not applicable"
  
  # Calculate median expected 2-sided p-value
  plot_data$medianp <- 2*pnorm(-abs(qnorm(plot_data$power)+qnorm(1-0.05/2)))
  
  # The color-blind palette with black as #3, move yellow to end:
  cbbPalette <- c("#0072B2","#D55E00","#000000","#009E73","#56B4E9","#E69F00","#CC79A7", 
                  "#F0E442")
  
  print(plot_data)
  
  plotpowers <- 
    ggplot(plot_data, aes(x=n,y=power,color=as.factor(RR_neg),shape=as.factor(RR_pos),
                          type=analysis)) +
    scale_color_manual(values=cbbPalette) +
    #scale_colour_manual(values=c("red","black","blue")) +
    geom_point(size=3)+
    geom_line(aes(lty=analysis)) +
    scale_linetype_manual(values=c("dotted", "solid"))+ # rev so solid is the usual analysis
    # scale_color_discrete(na.translate = F) +
    scale_x_continuous(breaks=unique(out.fig1[,1])) +
    # scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1)) +
    scale_y_continuous(trans = "probit", limits = c(0.2,0.9999), 
                       breaks=c(seq(0.2,0.9,0.1),0.95,0.98,0.99,0.999)) +
    geom_hline(yintercept = c(0.8,0.9),lty=3) +
    labs(x="total trial sample size (half in each arm)",y="power",shape="RR_pos",
         color="RR_neg",lty="Analysis") +
    # theme(legend.position = "none")  + # Hide legend
    NULL
  
  # No one intuits median expected p-value
  # plotmedianps <- 
  #   ggplot(plot_data, aes(x=n,y=medianp,color=as.factor(RR_neg),shape=as.factor(RR_pos),
  #                         type=analysis)) +
  #   scale_color_manual(values=cbbPalette) +
  #   geom_point(size=3)+
  #   geom_line(aes(lty=analysis)) +
  #   scale_linetype_manual(values=c("dotted", "solid"))+ # rev so solid is the usual analysis
  #   scale_x_continuous(breaks=unique(out.fig1[,1])) +
  #   # scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1)) +
  #   scale_y_log10(breaks=10^seq(-8,0,by=1),limits=c(1e-8,1)) + 
  #   geom_hline(yintercept = c(0.05,0.005,0.001),lty=3) +
  #   labs(x="total trial sample size (half in each arm)",
  #        y="median expected p-value over trial replications",
  #        shape="RR_pos",color="RR_neg",lty="Analysis") +
  #   NULL
}

# Load and plot the data
# Don't both plotting RRneg=0.95 to remove clutter
load("out.fig1.RR0.9.RData")
out.fig1 <- out.fig1[out.fig1$RR_neg!=0.95,]
plotpowers.RR0.9 <- plot.power.fig(out.fig1)
load("out.fig1.RR0.8.RData")
out.fig1 <- out.fig1[out.fig1$RR_neg!=0.95,]
plotpowers.RR0.8 <- plot.power.fig(out.fig1)

# Plot single analysis and subgroup analysis
final_plot <- cowplot::plot_grid(plotpowers.RR0.9, plotpowers.RR0.8,
                                 nrow=1, ncol=2, rel_heights = c(1, 1),
                                 labels = c("power (trial RR=0.9)",
                                            "power (trial RR=0.8)"))
# show the plot 
final_plot

# Save the plot
ggsave(filename="power-p-n.RR.pdf", plot=final_plot, device="pdf", 
       units="in", width=15, height=8,  dpi=500)




##
# Figure 4: Fix base case parameters, vary RRpos, P(M+), and sample size
##

# Simulate and save the data
basecase <- conceal.reveal(nsim=1,sim=TRUE)
# Simulation params for RR=0.9
# RR <- 0.9
# n <- c(25e3,50e3,75e3,100e3,150e3,200e3)
# RR_poss <- c(0.86667,0.8)
# p_Mpluss <- c(0.03,0.05,0.10,0.25,0.50)
# # Simulation params for RR=0.8
RR <- 0.8
n <- c(10e3, 15e3, 20e3, 25e3,30e3,40e3,50e3)
RR_poss <- c(0.76667,0.7)
p_Mpluss <- c(0.03,0.05,0.10, 0.25, 0.50)
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (k in 1:length(p_Mpluss)) {
  for (j in 1:length(RR_poss)) {
    for (i in 1:length(n)) {
      out <- rbind(out,
                   as.matrix(conceal.reveal(nsim=100000, n=n[i], p_Mplus=p_Mpluss[k], p_g0=0.5, pD_g0=0.02, 
                                            RR=RR,RR_neg=1, RR_pos=RR_poss[j],
                                            sim=TRUE))
      )
    }
  }
}
print(out.fig1 <- as.data.frame(cbind(n,out[-1,c(2,13:14)],rep(p_Mpluss,each=length(n)*length(RR_poss)))), row.names=FALSE)
names(out.fig1)<-c("n","RR_pos","IE","Standard Analysis","everpos")
# save(out.fig1,file="out.fig1.everpos.RR0.9.RData")
save(out.fig1,file="out.fig1.everpos.RR0.8.RData")





# Helper function to make plots
plot.power.fig <- function(out.fig1) {
  plot_data <- reshape2::melt(out.fig1,id.var=c("n","RR_pos","everpos"))
  colnames(plot_data)[4:5] <- c("analysis","power")
  plot_data$power[plot_data$power==1] <- 0.9999 # cannot plot a power of 1 on probit scale
  
  print(plot_data)
  
  # remove excess rows for usual trial RD_all (keep only the first such rows)
  plot_data <- rbind(plot_data[plot_data$analysis=="IE",],
                     plot_data[plot_data$analysis=="Standard Analysis",][1:length(unique(plot_data$n)),])
  # set everpos to -1 for the usual trial analysis to give it's own color
  plot_data$everpos[plot_data$analysis=="Standard Analysis"] <- "Standard Analysis"
  plot_data$RR_pos[plot_data$analysis=="Standard Analysis"] <- "Standard Analysis"
  
  # Calculate median expected 2-sided p-value
  plot_data$medianp <- 2*pnorm(-abs(qnorm(plot_data$power)+qnorm(1-0.05/2)))
  
  # The color-blind palette with black as #5, move yellow to end:
  cbbPalette <- c("#0072B2","#D55E00","#009E73","#56B4E9","#000000","#E69F00","#CC79A7", "#F0E442")
  
  print(plot_data)
  
  # Plot power
  # plotpowers <-
  #   # ggplot(plot_data, aes(x=n,y=power,color=as.factor(everpos),shape=as.factor(RR_pos),
  #   #                       type=analysis)) +
  #   ggplot(plot_data, aes(x=n,y=power,color=as.factor(everpos),type=analysis)) +
  #   scale_color_manual(values=cbbPalette) +
  #   #scale_colour_manual(values=c("red","black","blue")) +
  #   geom_point(size=3)+
  #   geom_line(aes(lty=analysis)) +
  #   scale_linetype_manual(values=c("dotted", "solid"))+ # rev so solid is the usual analysis
  #   # scale_color_discrete(na.translate = F) +
  #   scale_x_continuous(breaks=unique(out.fig1[,1])) +
  #   # scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1)) +
  #   scale_y_continuous(trans = "probit", limits = c(0.2,0.9999),
  #                      breaks=c(seq(0.2,0.9,0.1),0.95,0.98,0.99,0.999)) +
  #   geom_hline(yintercept = c(0.8,0.9),lty=3) +
  #   # labs(x="total trial sample size (half in each arm)",y="power",shape="RR_pos",
  #   #      color="everpos",lty="Analysis") +
  #   labs(x="total trial sample size (half in each arm)",y="power",
  #        color="Ever-positivity",lty="Analysis") +
  #   # theme(legend.position = "none")  + # Hide legend
  #   NULL
  
  # Plot median expected p-value
  plotmedianps <-
    ggplot(plot_data, aes(x=n/2,y=medianp,color=as.factor(everpos), type=analysis)) +
    scale_color_manual(values=cbbPalette) +
    geom_point(size=3)+
    geom_line(aes(lty=analysis)) +
    scale_linetype_manual(values=c("dotted", "solid"))+ # rev so solid is the usual analysis
    scale_x_continuous(breaks=unique(out.fig1[,1]/2)) +
    # scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1)) +
    scale_y_log10(breaks=10^seq(-8,0,by=1),limits=c(1e-8,1)) +
    geom_hline(yintercept = c(0.05,0.005,0.001),lty=3) +
    labs(x="trial per-arm sample size",
         # y="median expected p-value over trial replications",
         y="p-value",
         shape="RR_pos",color="everpos",lty="Analysis") +
    #theme with white background
    theme_bw() +
    #eliminates background, gridlines, #and chart border
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      # panel.border = element_blank()
    ) +
    NULL
}

# Load and plot the data
# Don't bother plotting RRneg=0.95 or everpos=0.25 to remove clutter
load("out.fig1.everpos.RR0.9.RData")
# out.fig1 <- out.fig1[out.fig1$RR_neg!=0.95,]
out.fig1 <- out.fig1[out.fig1$everpos!=0.25,]
plotpowers.RR0.9.RRpos0.867 <- plot.power.fig(out.fig1[out.fig1$RR_pos==0.86667,])
plotpowers.RR0.9.RRpos0.8 <- plot.power.fig(out.fig1[out.fig1$RR_pos==0.8,])
load("out.fig1.everpos.RR0.8.RData")
# out.fig1 <- out.fig1[out.fig1$RR_neg!=0.95,]
out.fig1 <- out.fig1[out.fig1$everpos!=0.25,]
plotpowers.RR0.8.RRpos0.767 <- plot.power.fig(out.fig1[out.fig1$RR_pos==0.76667,])
plotpowers.RR0.8.RRpos0.7 <- plot.power.fig(out.fig1[out.fig1$RR_pos==0.7,])

# Extract legend from one of the plots - they're all the same
legend_b <- cowplot::get_legend(
  plotpowers.RR0.9.RRpos0.867 + 
    guides(color = guide_legend(nrow = 1, byrow=TRUE)) +
    theme(legend.position = "bottom")
)

# Plot single analysis and subgroup analysis
# final_plot <- cowplot::plot_grid(plotpowers.RR0.9, plotpowers.RR0.8,
#                                  nrow=1, ncol=2, rel_heights = c(1, 1),
#                                  labels = c("power (trial RR=0.9)",
#                                             "power (trial RR=0.8)"))
prow <- cowplot::plot_grid(plotpowers.RR0.9.RRpos0.867 + theme(legend.position="none"),
                           plotpowers.RR0.9.RRpos0.8 + theme(legend.position="none"), 
                           plotpowers.RR0.8.RRpos0.767 + theme(legend.position="none"),
                           plotpowers.RR0.8.RRpos0.7 + theme(legend.position="none"),
                                 nrow=1, ncol=4, rel_heights = c(1, 1, 1, 1),
                                 labels = c("RR=0.9, RRpos=0.867\n   P_EV-pos=75%",
                                            "RR=0.9, RRpos=0.8\n   P_EV-pos=50%",
                                            "RR=0.8, RRpos=0.767\n   P_EV-pos=86%",
                                            "RR=0.8, RRpos=0.7\n    P_EV-pos=67%"))
# Include legen
final_plot <- cowplot::plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))

# show the plot 
final_plot

# Save the plot
# ggsave(filename="power-RRpos.pdf", plot=final_plot, device="pdf", 
#        units="in", width=16, height=8,  dpi=500)
ggsave(filename="median-p-RRpos.pdf", plot=final_plot, device="pdf", 
       units="in", width=16, height=8,  dpi=500)














##
# Corner Case, Strange: Example where stronger RRpos has worse p-value
##
Mpos_0.8 <- matrix( c(200,250,112.5,62.5),nrow=2,byrow=TRUE,
                    dimnames=list(c("D+","D-"),c("screen","control")))
Mpos_0.7 <- matrix( c(116.6667,166.6667,195.8333,145.8333),nrow=2,byrow=TRUE,
                    dimnames=list(c("D+","D-"),c("screen","control")))
Mall <- matrix( c(200,250,12300,12250),nrow=2,byrow=TRUE,
                dimnames=list(c("D+","D-"),c("screen","control")))

print(Mpos)
str(prop.test(t(Mpos_0.8),correct=F))
str(fisher.test(t(Mpos_0.8)))
str(prop.test(t(Mpos_0.7),correct=F))
str(fisher.test(t(Mpos_0.7)))
str(prop.test(t(Mall),correct=F))
fisher.test(t(Mall))




##
# Fix base params, vary sample size for a large simple 1 screen trial, 3-year result
##
basecase <- conceal.reveal(nsim=1,sim=TRUE)
n <- seq(150e3,350e3,by=200e3) 
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (i in 1:length(n)) {
  out <- rbind(out,
     as.matrix(conceal.reveal(nsim=5000, n=n[i], p_Mplus=0.01, p_g0=0.5, pD_g0=3*(0.02/5), 
                              RR=0.96,RR_neg=1, RR_pos=  0.95,
                              frac.g0.plus = 0.95, frac.g0.minus = 0.5, 
                              frac.g0.plus.pos = 0.9, frac.g0.minus.pos = 0.8,sim=TRUE))
  )
}
print(cbind(n,out[-1,]))
print(basecase <- as.data.frame(cbind(n,out[-1,c(12:17)])), digits=2, row.names=FALSE)

# Vectorize conceal.reveal() so we can input pairs of RR_neg,RR_pos to get 
# theoretical power for RR_pos based on the computed Z_ratio
# Remember that Z_ratio doesn't account for sampling or loss of signal
# And presumes that you've set the parameters so that power for the RR is 90%
conceal.reveals <- Vectorize(function(RR_neg,RR_pos,nsim, n, p_Mplus, 
                                      p_g0, pD_g0, RR, sim) 
  tryCatch(conceal.reveal(RR_neg=RR_neg, RR_pos=RR_pos, nsim=nsim,n=n, p_Mplus=p_Mplus, p_g0=p_g0, pD_g0=pD_g0, 
                          RR=RR),
           error=function(e){return(NA)})
)

# Examine the range of RRneg and RRpos to see power for 1 baseline test in control arm
RR_negs <- seq(0.97,1.03,by=0.01)
RR_poss <- seq(0.905,0.965, by=0.01)
# parameters set so power=14% for usual RR (extreme case)
out <- outer(RR_negs, RR_poss, conceal.reveals, nsim=1, n=150e3, p_Mplus=0.01, p_g0=0.5, 
             pD_g0=3*(0.02/5), RR=0.96, sim=FALSE)
rownames(out) <- signif(RR_negs,digits=3) ; colnames(out) <- signif(RR_poss,digits=3)
signif(out,digits=2)






# Examine the range of RRneg and RRpos to see power for 1 baseline test in control arm
RR_negs <- seq(0.90,0.93,by=0.005)
RR_poss <- seq(0.5,0.9,by=1/30)
# parameters set so power 80% for usual RR???
out <- outer(RR_negs, RR_poss, conceal.reveals, nsim=1, n=150e3, p_Mplus=0.01, p_g0=0.5, 
             pD_g0=0.02, RR=0.9, sim=FALSE)
rownames(out) <- signif(RR_negs,digits=3) ; colnames(out) <- signif(RR_poss,digits=3)
signif(out,digits=2)



# Examine only a baseline test in control arm only
print(basecase<-conceal.reveal(nsim=10000,n=150e3, p_Mplus=0.01, p_g0=0.5, pD_g0=0.02, 
                               RR=0.9, RR_neg=0.925, RR_pos=0.86667,sim=TRUE))
print(basecase<-conceal.reveal(nsim=10000,n=150e3, p_Mplus=0.01, p_g0=0.5, pD_g0=0.02, 
                               RR=0.9, RR_neg=0.925, RR_pos=0.8,sim=TRUE))
print(conceal.reveal(nsim=10000,n=200e3, p_Mplus=0.01, p_g0=0.5, pD_g0=0.02, 
                     RR=0.9, RR_neg=0.915, RR_pos=0.8, sim=TRUE))




##
# Fix base case parameters, but vary sample size
##
basecase <- conceal.reveal(nsim=1,sim=TRUE)
n <- c(50e3,75e3,100e3,150e3,200e3) #seq(120e3,150e3,by=10e3)
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (i in 1:length(n)) {
  out <- rbind(out,
     as.matrix(conceal.reveal(nsim=5000, n=n[i], p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
                              RR=0.9,RR_neg=1, RR_pos=  0.86667,
                              frac.g0.plus = 0.95, frac.g0.minus = 0.5, 
                              frac.g0.plus.pos = 0.9, frac.g0.minus.pos = 0.8, sim=TRUE))
  )
}
print(cbind(n,out[-1,]))
print(basecase <- as.data.frame(cbind(n,out[-1,c(12:14)])), digits=2, row.names=FALSE)




##
# Examine this as the case of estimating test*arm interaction
##
print(trial <- data.frame(outcome=rep(c(0,1),4),
                   arm=rep(c("control","control","screen","screen"),2),
                   test=c(rep(0,4),rep(1,4)),
                   Freq=c(as.vector(Mneg[2:1,2:1]),as.vector(Mpos[2:1,2:1]))
))
# reproduce log(OR)=log(0.86667)=-0.1987 among test-positives, p=8.5e-06
summary(glm(outcome~ arm, weights=Freq, family=binomial(link="logit"), data=trial,
            subset=(test==1)))
# calculate the sd of the test+-only interaction manually, sd=0.04463
with(trial[trial$test==1,],sqrt(sum(1/Freq)))

# Interaction has same log(OR)=-0.1987 but p=0.01, showing power of test+-only test
summary(glm(outcome~ test*arm, weights=Freq, family=binomial(link="logit"), data=trial))
# calculate the sd of the interaction manually, sd=0.07754
with(trial,sqrt(sum(1/Freq)))


# # case-only analysis: OR=0.8667 exactly (which is the RR) p=0.052; unsure what this means
# # only works if RR_neg=1 exactly, makes no sense otherwise
# summary(glm(test~ arm, weights=Freq, family=binomial(link="logit"), data=trial,
#             subset=(outcome==1)))


##
# Do subsampling
# Fix n=100k, but now vary frac.g0.minus from 0 to 1
##
basecase<-conceal.reveal(nsim=1, sim=TRUE)
fracs <-seq(0.1,1,by=0.1)
RR_pos <- c(0.8667,0.8,0.7)
n <- c(100e3,75e3,50e3)
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (j in 1:length(RR_pos)) {
  for (i in 1:length(fracs)) {
    out <- rbind(out,
                 as.matrix(conceal.reveal(nsim=10000, n=n[j], p_Mplus=0.05, 
                                          p_g0=0.5, pD_g0=0.02, 
                                          RR=0.9,RR_neg=1, RR_pos= RR_pos[j],
                                          frac.g0.plus = 1, frac.g0.minus = fracs[i], 
                                          frac.g0.plus.pos = 1, frac.g0.minus.pos = 1, 
                                          sim=TRUE))
    )
  }
}
print(out[-1,])
print(basecase <- as.data.frame(cbind(fracs,out[-1,3+c(23,24)])), digits=2, row.names=FALSE)


print(out.fig1 <- as.data.frame(cbind(n,out[-1,c(2,3,13:14)])), row.names=FALSE)
names(out.fig1)<-c("n","RR_pos","RR_neg","ICE","usual")
# save(out.fig1,file="out.fig1.RData")
save(out.fig1,file="out.fig1.RR0.8.RData")




# print(basecase <- cbind(fracs,out[-1,c(7:9,12:36)]), row.names=FALSE )

##
# Fix n=100k, but now vary non-differential loss of signal from 0.5 to 1
##
basecase<-conceal.reveal(nsim=1,sim=TRUE)
fracs <-seq(0.5,1,by=0.1)
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (i in 1:length(fracs)) {
  out <- rbind(out,
               as.matrix(conceal.reveal(
                 nsim=10000, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
                 RR=0.9,RR_neg=1, RR_pos= 0.86667,
                 frac.g0.plus = 0.95, frac.g0.minus = 0.5, 
                 frac.g0.plus.pos = fracs[i], frac.g0.minus.pos = fracs[i], sim=TRUE))
  )
}
print(out[-1,])
print(basecase <- as.data.frame(cbind(fracs,out[-1,4+c(24,25,27,28,29,30,32,33)])), digits=2, row.names=FALSE)

##
# Fix n=100k, but now vary differential loss of signal from 0.5 to 1 and 0.4 to 0.9
##
basecase<-conceal.reveal(nsim=1,sim=TRUE)
fracs.plus <-c(seq(0.5,1,by=0.1),1)
fracs.minus <- c(seq(0.5,1,by=0.1),1.1) - 0.1 # allow for a base-case of 1.0 fraction
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (i in 1:length(fracs.plus)) {
  out <- rbind(out,
               as.matrix(conceal.reveal(
                 nsim=10000, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
                 RR=0.9,RR_neg=1, RR_pos= 0.86667,
                 frac.g0.plus = 0.95, frac.g0.minus = 0.5, 
                 frac.g0.plus.pos = fracs.plus[i], frac.g0.minus.pos = fracs.minus[i], 
                 sim=TRUE))
  )
}
print(out[-1,])
print(basecase <- as.data.frame(cbind(fracs.plus,fracs.minus,
                                      out[-1,4+c(24,25,27,28,29,30,32,33)])), digits=2, row.names=FALSE)
print(basecase <- as.data.frame(cbind(fracs.plus,fracs.minus,
                                      out[-1,4+c(36,37,40,41,42,43,46,47)])), digits=2, row.names=FALSE)


















# OLD
############ 

##
# Power for sample sizes for 20k to 200k and 10% sampling of D- control arm: base case
##
basecase<-conceal.reveal(nsim=50, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                         RR_neg=1, RR_pos= 0.86667,
                         frac.g0.plus = 1, frac.g0.minus = 0.1)
ns <- seq(20e3,200e3,by=20e3)
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (i in 1:length(ns)) {
  out <- rbind(out,
               as.matrix(conceal.reveal(nsim=5000, n=ns[i], p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                             RR_neg=1, RR_pos= 0.86667,
                             frac.g0.plus = 0.95, frac.g0.minus = 0.25))
  )
}
print(basecase <- cbind(ns,out[-1,12:26]), row.names=FALSE )

##
# Power for sample sizes for 20k to 200k: base case except RR_neg=1.1 false-reassurance
##
basecase<-conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                               RR_neg=1, RR_pos= 0.86667)
ns <- seq(20e3,200e3,by=20e3)
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (i in 1:length(ns)) {
  out <- rbind(out,
               as.matrix(conceal.reveal(nsim=1000, n=ns[i], p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                                        RR_neg=1.1, RR_pos= 0.86667))
  )
}
print(basecase.RR_neg <- cbind(ns,out[-1,12:19]), row.names=FALSE )

##
# Power for sample sizes for 20k to 200k: base case except RR_neg=0.91 "unreassuurance"
##
basecase<-conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                         RR_neg=1, RR_pos= 0.86667)
ns <- seq(20e3,200e3,by=20e3)
out <- matrix(NA,nrow=1,ncol=ncol(basecase))
for (i in 1:length(ns)) {
  out <- rbind(out,
               as.matrix(conceal.reveal(nsim=1000, n=ns[i], p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                                        RR_neg=0.92, RR_pos= 0.86667))
  )
}
print(basecase.RR_neg <- cbind(ns,out[-1,12:19]), row.names=FALSE )

  



##
# Old simulation code below
##


# base-case, cut sample size 200k, 150k, 100k
# 100k still has 90% power for RD_pos but only 64% power for RD
# Note that Z-ratio does not change - it is independent of sample size
# The ratio reduction in p-value goes from 155 to 33 to 10

# base case scenario with 90% power for standard trial analysis
print( out<-conceal.reveal(nsim=10000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                           RR_neg=1, RR_pos= 0.86667), row.names=FALSE )

#  RR RR_pos RR_neg RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim
# 0.9 0.8667      1      0   0.04 0.002  6.693e-06    0.04025 0.002012    1.36        1.36
# power_neg power_pos power_all median.p_neg median.p_pos median.p_all p.ratio p_pos.gt.p_all
#     0.052     0.991     0.911        0.493    7.129e-06      0.00111   155.7             43
# pD_g1.pos pD_g0.pos pD_g1.neg pD_g0.neg pMplus_Dminus pMplus_Dplus pD_pos   pD_neg    pD
#      0.26       0.3  0.005263  0.005263        0.0367       0.7369   0.28 0.005263 0.019

# cut sample size to 75k each arm
print( out1<-conceal.reveal(nsim=1000, n=150e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
                            RR=0.9, RR_neg=1, RR_pos= 0.86667), row.names=FALSE )

# cut sample size to 50k each arm
print( out2<-conceal.reveal(nsim=1000, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                            RR_neg=1, RR_pos= 0.86667), row.names=FALSE )

rbind(out,out1,out2)
# power.neg power.pos power.all  p.neg     p.pos    p.all p.ratio p.pos.gt.p.all
#     0.052     0.991     0.911 0.4930 7.129e-06 0.001110  155.72          0.043
#     0.050     0.971     0.794 0.5094 1.384e-04 0.004605   33.28          0.094
#     0.041     0.884     0.637 0.4931 1.784e-03 0.018601   10.43          0.143


##
# base-case then increase P(M+).  
# The Z-ratio is cut more than half by increasing P(M+) from 5% to 10% to 50%
# Going just from 5% to 10% cuts p-ratio from 155 to 21, sensitive because power is so high
# However power remains 97% at P(M+)=50%, so insensitive in that sense
##

out<-conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                           RR_neg=1, RR_pos= 0.86667)

out1<-conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.10, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                           RR_neg=1, RR_pos= 0.86667)

out2<-conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.50, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                            RR_neg=1, RR_pos= 0.86667)

rbind(out,out1,out2)

# power.neg power.pos power.all  p.neg     p.pos    p.all p.ratio p.pos.gt.p.all
#     0.052     0.991     0.911 0.4930 7.129e-06 0.001110 155.723          0.043
#     0.058     0.986     0.889 0.4785 4.482e-05 0.000967  21.583          0.112
#     0.052     0.969     0.905 0.5098 1.338e-04 0.001087   8.125          0.153



##
# base-case, reduce to 4 and 3 screens with P(D+)=1.5% and 1% respectively
# power for usual analysis drops to 80% and 66%, while for screen positives 96% to 86%
##

out<-conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                     RR_neg=1, RR_pos= 0.86667)

out1<-conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.04, p_g0=0.5, pD_g0=0.015, RR=0.9, 
                     RR_neg=1, RR_pos= 0.86667)

out2<-conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.03, p_g0=0.5, pD_g0=0.01, RR=0.9, 
                     RR_neg=1, RR_pos= 0.86667)

rbind(out,out1,out2)

# power.neg power.pos power.all  p.neg     p.pos    p.all p.ratio p.pos.gt.p.all
#     0.052     0.991     0.911 0.4930 7.129e-06 0.001110 155.723          0.043
#     0.047     0.967     0.805 0.4911 1.438e-04 0.004507  31.347          0.101
#     0.055     0.864     0.663 0.4924 2.304e-03 0.016397   7.116          0.164



##
# base-case, but now false reassurance RR_neg from 1.1 to 1.2
# Power for screen-positives increases, p-value sharply decreases
# No power to see there is false-reassurance in screen-positives: 20% to 50%
# Note that RR_neg=1.2, but RR_pos=0.8667 means that there is more harm to screen-negatives
# (1/1.2=0.83) than there is benefit to screen positives, which is perverse but possible
##

 out<-conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                     RR_neg=1, RR_pos= 0.86667)

out1<-conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                     RR_neg=1.1, RR_pos= 0.86667)

out2<-conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                     RR_neg=1.2, RR_pos= 0.86667)

rbind(out,out1,out2)
# power.neg power.pos power.all   p.neg     p.pos     p.all p.ratio p.pos.gt.p.all
#     0.052     0.991     0.911 0.49296 7.129e-06 0.0011102   155.7          0.043
#     0.219     0.999     0.905 0.22988 9.461e-07 0.0009181   970.3          0.005
#     0.502     1.000     0.895 0.04874 4.120e-07 0.0011737  2848.9          0.003




##
# Differential uptake of SOC screening
##

# Change base case so pD_g0.pos=0.2, otherwise RR_neg blows up sometimes
print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0,
                      pD_g0.pos=0.2), row.names=FALSE )
#  RR RR_pos RR_neg RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim
# 0.9    0.8      1      0   0.04 0.002  9.477e-06    0.04048 0.002035   1.589        1.58
# power_neg power_pos power_all median.p_neg median.p_pos median.p_all
#      0.05         1     0.907       0.4774    1.651e-07     0.000936

# More SOC uptake in control arm: RD_neg=-0.001 but also 
# decrease P(D|G=0) from 2% to 1.5%, and decrease P(D+|M+,G=0) from 20% to 15%, RR from 0.9 to 0.95
# RD_pos decreases from 0.04 to 0.034 (RR_pos from 0.8 to 0.83), so effect size decreases
# Power for RD_pos decreases, but power for RD collapses
print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.015, RR=0.95, RD_neg=-0.001,
                      pD_g0.pos=0.15), row.names=FALSE )
#   RR RR_pos RR_neg RD_neg RD_pos      RD RD_neg.sim RD_pos.sim    RD.sim Z_ratio Z_ratio.sim
# 0.95 0.7733  1.127 -0.001  0.034 0.00075  -0.001003    0.03397 0.0007467   3.584       3.405
# power_neg power_pos power_all median.p_neg median.p_pos median.p_all
#    0.6692     0.999    0.2854      0.01706    5.845e-07       0.1574

# More SOC uptake in screen arm: RD_neg=0.001 but also RR decreases from 0.9 to 0.85
# RD_pos increases from 0.04 to 0.041 but now RD_reg=0.001.  Effects nearly cancel on the Z-stat
# Power for RD_pos is same as base-case, but power for RD increases a lot
print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.85, RD_neg=0.001,
                      pD_g0.pos=0.2), row.names=FALSE )
#   RR RR_pos RR_neg RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim
# 0.85  0.795  0.905  0.001  0.041 0.003   0.001007     0.0411 0.003015   1.073        1.07
# power_neg power_pos power_all median.p_neg median.p_pos median.p_all
#    0.5963    0.9997    0.9988      0.02728    9.436e-08    5.516e-07












###
# OLD CODE
#####

#####
conceal.reveal.old <- function(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0,
                               pD_g0.pos=0.3)
{
  
  # P(M+) and P(G=0)
  # p_Mplus <- 0.05 ; p_g0 <- 0.5 # P(G=0) control arm
  
  # Fix P(D+|G=0) as population disease mortality and expected mortality RR=P(D+|G=1)/P(D+|G=0)
  # pD_g0 <- 0.02
  # RR <- 0.90
  # Calculate P(D+|G=1) = RR*P(D+|G=0), and true RD=P(D+|G=0)-P(D+|G=1)
  pD_g1 <- RR*pD_g0
  RD <- pD_g0 - pD_g1
  
  # Given RD and RD_neg=P(D+|M-,G=0)-P(D+|M-,G=1), solve for RD_pos=P(D+|M+,G=0)-P(D+|M+,G=1)
  # RD_neg <- 0
  RD_pos <- (RD-RD_neg*(1-p_Mplus)) / p_Mplus
  
  # Specify PPV=P(D+|M+,G=0)
  # Check for compatibilty with marginals by checking if sens=P(M+|D+,G=0)=PPV*P(M+)/P(D+|G=0) > 1
  # pD_g0.pos <- 0.3
  if (pD_g0.pos*p_Mplus/pD_g0>1) print("incompatible PPV conditional because sens>1")
  if (pD_g0.pos*p_Mplus/pD_g0<RD_pos) print("PPV too small given large RD_pos")
  
  # Now calculate the rest of the conditionals
  # P(D+|M+,G=1) = P(D+|M+,G=0) - RD_pos
  pD_g1.pos <- pD_g0.pos - RD_pos
  paste("conceal-reveal RR:", RR_pos <- pD_g1.pos/pD_g0.pos)
  # back out P(D+|M-,G=0) from P(D+|G=0) = P(D+|M-,G=0)P(M-) + P(D+|M+,G=0)P(M+)
  pD_g0.neg <- (pD_g0 - pD_g0.pos*p_Mplus) / (1-p_Mplus)
  # P(D+|M-,G=1) = P(D+|M-,G=0) - RD_pos
  pD_g1.neg <- pD_g0.neg - RD_neg
  RR_neg <- pD_g1.neg/pD_g0.neg # other conceal-reveal RR
  
  # Theoretical ratio of Z-statistics: RD_pos/RD * P(M+) * sqrt{P(M+)/(P(M+|D+)P(M+|D-)}, where:
  # P(M+|D+) = P(D+|M+)P(M+)/P(D+)
  # P(M+|D-) = P(D-|M+)(M+)/P(D-)
  # P(D+|M+) = P(D+|M+,G=0)P(G=0) + P(D+|M+,G=1)P(G=1)
  pD_pos <- pD_g0.pos*p_g0 + pD_g1.pos*(1-p_g0)
  # P(D+|M-) = P(D+|M-,G=0)P(G=0) + P(D+|M-,G=1)P(G=1)
  pD_neg <- pD_g0.neg*p_g0 + pD_g1.neg*(1-p_g0)
  # P(D+) = P(D+|M+)P(M+) + P(D+|M-)P(M-)
  pD <- pD_pos*p_Mplus + pD_neg*(1-p_Mplus)
  # P(M+|D+) = P(D+|M+)P(M+)/P(D+)
  pMplus_Dplus <- pD_pos*p_Mplus/pD
  # P(M+|D-) = P(D-|M+)P(M+)/P(D-)
  pMplus_Dminus <- (1-pD_pos)*p_Mplus/(1-pD)
  
  paste("theoretical Z-ratio: ",
        Z_ratio <- RD_pos/RD * p_Mplus * sqrt(p_Mplus/(pMplus_Dminus*pMplus_Dplus)) )
  
  # Simulate 2x2 tables
  set.seed(406551735)
  # nsim <- 1000
  Z_neg <- Z_pos <- Z_all <- p.val_neg <- p.val_pos <- p.val_all <- RD_neg.sim <- RD_pos.sim <- RD_all.sim <- numeric(nsim)
  # n <- 100e3
  for(i in 1:nsim) {
    
    # Generate M+ and M- sample sizes, fixing total sample size n in each group
    g1.pos <- rbinom(1, n*(1-p_g0), p_Mplus) ; g1.neg <- n*(1-p_g0)-g1.pos
    g0.pos <- rbinom(1, n*p_g0, p_Mplus) ; g0.neg <- n*p_g0-g0.pos
    
    # generate the M- and M+ tables
    Dplus_g1.neg <- rbinom(1, g1.neg, pD_g1.neg);  Dplus_g0.neg <- rbinom(1, g0.neg, pD_g0.neg)
    Dplus_g1.pos <- rbinom(1, g1.pos, pD_g1.pos);  Dplus_g0.pos <- rbinom(1, g0.pos, pD_g0.pos)
    
    # calculate statistics for M-, M+, and marginal tables
    neg.out <- prop.test(c(Dplus_g1.neg,Dplus_g0.neg), c(g1.neg,g0.neg), correct=F)
    pos.out <- prop.test(c(Dplus_g1.pos,Dplus_g0.pos), c(g1.pos,g0.pos), correct=F)
    all.out <- prop.test(c(Dplus_g1.neg+Dplus_g1.pos, Dplus_g0.neg+Dplus_g0.pos),
                         c(g1.neg+g1.pos, g0.neg+g0.pos), correct=F)
    
    # RDs
    RD_neg.sim[i] <- diff(neg.out$estimate)
    RD_pos.sim[i] <- diff(pos.out$estimate)
    RD_all.sim[i] <- diff(all.out$estimate)
    
    # p-values for the RDs
    p.val_neg[i] <- neg.out$p.val
    p.val_pos[i] <- pos.out$p.val
    p.val_all[i] <- all.out$p.val
    
    # Z-stats
    Z_neg[i]  <- sqrt(neg.out$stat)
    Z_pos[i]  <- sqrt(pos.out$stat)
    Z_all[i]  <- sqrt(all.out$stat)
  }
  
  # Pack up output for return
  out <- data.frame(RR=RR, RR_pos=RR_pos, RR_neg=RR_neg, RD_neg=RD_neg, RD_pos=RD_pos, RD=RD,
                    RD_neg.sim=mean(RD_neg.sim), RD_pos.sim=mean(RD_pos.sim), RD.sim=mean(RD_all.sim),
                    Z_ratio=Z_ratio,
                    Z_ratio.sim = mean(Z_pos)/mean(Z_all),
                    power_neg=mean(p.val_neg <= .05),
                    power_pos=mean(p.val_pos <= .05),
                    power_all=mean(p.val_all <= .05),
                    median.p_neg=median(p.val_neg),
                    median.p_pos=median(p.val_pos),
                    median.p_all=median(p.val_all),
                    p_pos.gt.p_all=sum(p.val_pos>p.val_all),
                    pD_g1.pos=pD_g1.pos,
                    pD_g0.pos=pD_g0.pos,
                    pD_g1.neg=pD_g1.neg,
                    pD_g0.neg=pD_g0.neg,
                    pMplus_Dminus=pMplus_Dminus,                    
                    pMplus_Dplus=pMplus_Dplus
                    
  )
  #               pD=pD, p_Mplus=p_Mplus, p_g0=p_g0, pD_pos=pD_pos, pD_neg=pD_neg,
  #               pMplus_Dplus=pMplus_Dplus, pMplus_Dminus=pMplus_Dminus,
  # )
  
  return(out)
  
  # return(list(RD_neg.sim=RD_neg.sim,RD_pos.sim=RD_pos.sim,RD_all.sim=RD_all.sim,
  #              p.val_neg=p.val_neg,p.val_pos=p.val_pos,p.val_all=p.val_all,
  #              Z_neg=Z_neg,Z_pos=Z_pos,Z_all=Z_all)
  #        )
}

# base case scenario with 90% power for standard trial analysis
print( out<-conceal.reveal.old(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0,
                               pD_g0.pos=0.3), row.names=FALSE )
print( out<-conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, 
                           RR_neg=1, RR_pos= 0.86667), row.names=FALSE )

#  RR RR_pos RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim power_neg
# 0.9 0.8667      0   0.04 0.002  7.992e-06    0.04019 0.002006    1.36       1.362      0.05
# power_pos power_all median.p_neg median.p_pos median.p_all
#     0.991     0.908       0.4941    7.237e-06     0.001147

# cut sample size to 75k each arm
print( conceal.reveal(nsim=1000, n=150e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0,
                      pD_g0.pos=0.3), row.names=FALSE )
#  RR RR_pos RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim power_neg
# 0.9 0.8667      0   0.04 0.002  3.274e-06    0.03993 0.002005    1.36       1.353     0.052
# power_pos power_all median.p_neg median.p_pos median.p_all
#     0.975     0.802       0.5104     0.000122     0.004245

# cut sample size to 50k each arm
print( conceal.reveal(nsim=1000, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0,
                      pD_g0.pos=0.3), row.names=FALSE )
#  RR RR_pos RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim power_neg
# 0.9 0.8667      0   0.04 0.002  -3.16e-06    0.04014 0.002014    1.36       1.347     0.042
# power_pos power_all median.p_neg median.p_pos median.p_all
#     0.885     0.639       0.4957     0.001675      0.01813


##
# base-case, reduce PPV=P(D+|M+,G=0): usual trial RD p-val is insensitive
##

print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0,
                      pD_g0.pos=0.3), row.names=FALSE )
#  RR RR_pos RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim power_neg
# 0.9 0.8667      0   0.04 0.002  7.992e-06    0.04019 0.002006    1.36       1.362      0.05
# power_pos power_all median.p_neg median.p_pos median.p_all
#   0.991     0.908       0.4941    7.237e-06     0.001147

print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0,
                      pD_g0.pos=0.2), row.names=FALSE )
#  RR RR_pos RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim power_neg
# 0.9    0.8      0   0.04 0.002  9.477e-06    0.04048 0.002035   1.589        1.58      0.05
# power_pos power_all median.p_neg median.p_pos median.p_all
#        1     0.907       0.4774    1.651e-07     0.000936

print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0,
                      pD_g0.pos=0.1), row.names=FALSE )
#  RR RR_pos RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim power_neg
# 0.9    0.6      0   0.04 0.002 -4.184e-08    0.03985 0.001996   2.251       2.244     0.055
# power_pos power_all median.p_neg median.p_pos median.p_all
#       1     0.891       0.4795    2.199e-13     0.001045


##
# base-case [reduce PPV to 0.05 and cut sample to 25k to see p-value reduction], increase P(M+).  
# The Z-ratio is cut more than half by increasing P(M+) from 5% to 20%, much less p-value reduction
##

print( conceal.reveal(nsim=1000, n=50e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0,
                      pD_g0.pos=0.05), row.names=FALSE )
#  RR RR_pos RD_neg RD_pos    RD RD_neg.sim RD_pos.sim  RD.sim Z_ratio Z_ratio.sim power_neg
# 0.9    0.2      0   0.04 0.002 -3.735e-05    0.03993 0.00196   3.579       3.516     0.052
# power_pos power_all median.p_neg median.p_pos median.p_all
#        1     0.366       0.4972      5.2e-09      0.09986

print( conceal.reveal(nsim=1000, n=50e3, p_Mplus=0.2, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0,
                      pD_g0.pos=0.05), row.names=FALSE )
#  RR RR_pos RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim power_neg
# 0.9    0.8      0   0.01 0.002  4.352e-06    0.00987 0.001975   1.473       1.449     0.059
# power_pos power_all median.p_neg median.p_pos median.p_all
#     0.654     0.351       0.4963      0.02106       0.1125


##
# base-case, reduce P(D+) to 1% and P(M+) to 3% as an "early stop after 3 screens" example
# RD_pos has 80% power, RD has only 64% power
##

print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0,
                      pD_g0.pos=0.3), row.names=FALSE )
#  RR RR_pos RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim power_neg
# 0.9 0.8667      0   0.04 0.002  7.992e-06    0.04019 0.002006    1.36       1.362      0.05
# power_pos power_all median.p_neg median.p_pos median.p_all
#    0.991     0.908       0.4941    7.237e-06     0.001147

print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.03, p_g0=0.5, pD_g0=0.01, RR=0.9, RD_neg=0,
                      pD_g0.pos=0.3), row.names=FALSE )
#  RR RR_pos RD_neg  RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim power_neg
# 0.9 0.8889      0 0.03333 0.001  4.635e-06     0.0333 0.001004   1.243       1.233      0.05
# power_pos power_all median.p_neg median.p_pos median.p_all
#    0.806     0.635       0.4797     0.004323      0.02045

# If you also reduce PPV=P(D+|M+,G=0) to 0.1, power for RD_pos shoots up to >0.999
print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.03, p_g0=0.5, pD_g0=0.01, RR=0.9, RD_neg=0,
                      pD_g0.pos=0.1), row.names=FALSE )
#  RR RR_pos RD_neg  RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim power_neg
# 0.9 0.6667      0 0.03333 0.001  9.852e-06    0.03338 0.001009   2.026       2.005     0.045
# power_pos power_all median.p_neg median.p_pos median.p_all
#        1     0.654       0.4859    2.889e-06      0.02005



##
# base-case, but now false reassurance RD_neg=-0.001.  This fixes RD, so the RD_pos gets stronger.
# The median p for RD_pos explodes down from 1e-6 to 1e-11.
# If RD_neg=+0.001, power for RD_pos is worse than for the RD.
#     * But RR_neg=0.81 is better than RR_pos=0.93, which is really unreasonable
# RD_neg=0.00055 implies Z_ratio=1, thus same power for RD_pos and RD
##

print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0,
                      pD_g0.pos=0.3), row.names=FALSE )
#  RR RR_pos RR_neg RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim
# 0.9 0.8667      1      0   0.04 0.002  7.992e-06    0.04019 0.002006    1.36       1.362
# power_neg power_pos power_all median.p_neg median.p_pos median.p_all
#      0.05     0.991     0.908       0.4941    7.237e-06     0.001147

print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=-0.001,
                      pD_g0.pos=0.3), row.names=FALSE )
#  RR RR_pos RR_neg RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim
# 0.9 0.8033   1.19 -0.001  0.059 0.002  -0.001008    0.05893 0.001986   2.027       2.039
# power_neg power_pos power_all median.p_neg median.p_pos median.p_all
#      0.82         1     0.908     0.003649    3.538e-11     0.001144

# Example of RD_pos>0, in particular 0.001, but this is unreasonable because RR_pos>RR_neg  
print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0.001,
                      pD_g0.pos=0.3), row.names=FALSE )
#  RR RR_pos RR_neg RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim
# 0.9   0.93   0.81  0.001  0.021 0.002   0.001003    0.02114 0.002009  0.7068      0.7108
# power_neg power_pos power_all median.p_neg median.p_pos median.p_all
#     0.891     0.632     0.913     0.001578      0.01937     0.001002

# Example of RD_pos>0, in particular 0.000055 has Z_ratio=1 so same power for RD_pos and RD
print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0.00055,
                      pD_g0.pos=0.3), row.names=FALSE )
#  RR RR_pos RR_neg  RD_neg  RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim
# 0.9 0.9015 0.8955 0.00055 0.02955 0.002  0.0005572    0.02957 0.002003   0.999      0.9977
# power_neg power_pos power_all median.p_neg median.p_pos median.p_all
#       0.4     0.907     0.911      0.08916      0.00108     0.001079


##
# Differential uptake of SOC screening
##

# Change base case so pD_g0.pos=0.2, otherwise RR_neg blows up sometimes
print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9, RD_neg=0,
                      pD_g0.pos=0.2), row.names=FALSE )
#  RR RR_pos RR_neg RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim
# 0.9    0.8      1      0   0.04 0.002  9.477e-06    0.04048 0.002035   1.589        1.58
# power_neg power_pos power_all median.p_neg median.p_pos median.p_all
#      0.05         1     0.907       0.4774    1.651e-07     0.000936

# More SOC uptake in control arm: RD_neg=-0.001 but also 
# decrease P(D|G=0) from 2% to 1.5%, and decrease P(D+|M+,G=0) from 20% to 15%, RR from 0.9 to 0.95
# RD_pos decreases from 0.04 to 0.034 (RR_pos from 0.8 to 0.83), so effect size decreases
# Power for RD_pos decreases, but power for RD collapses
print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.015, RR=0.95, RD_neg=-0.001,
                      pD_g0.pos=0.15), row.names=FALSE )
#   RR RR_pos RR_neg RD_neg RD_pos      RD RD_neg.sim RD_pos.sim    RD.sim Z_ratio Z_ratio.sim
# 0.95 0.7733  1.127 -0.001  0.034 0.00075  -0.001003    0.03397 0.0007467   3.584       3.405
# power_neg power_pos power_all median.p_neg median.p_pos median.p_all
#    0.6692     0.999    0.2854      0.01706    5.845e-07       0.1574

# More SOC uptake in screen arm: RD_neg=0.001 but also RR decreases from 0.9 to 0.85
# RD_pos increases from 0.04 to 0.041 but now RD_reg=0.001.  Effects nearly cancel on the Z-stat
# Power for RD_pos is same as base-case, but power for RD increases a lot
print( conceal.reveal(nsim=1000, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.85, RD_neg=0.001,
                      pD_g0.pos=0.2), row.names=FALSE )
#   RR RR_pos RR_neg RD_neg RD_pos    RD RD_neg.sim RD_pos.sim   RD.sim Z_ratio Z_ratio.sim
# 0.85  0.795  0.905  0.001  0.041 0.003   0.001007     0.0411 0.003015   1.073        1.07
# power_neg power_pos power_all median.p_neg median.p_pos median.p_all
#    0.5963    0.9997    0.9988      0.02728    9.436e-08    5.516e-07





