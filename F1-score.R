## Computation of F1 score and its confidence interval ##

## Takahashi K, Yamamoto K, Kuchiba A, Koyama T. Confidence interval for micro-averaged F1 and macro-averaged F1 scores.  Appl Intell (Dordr). 2022 Mar;52(5):4961-4972. doi: 10.1007/s10489-021-02635-5. Epub 2021 Jul 31.

## PMID: 35317080
## PMCID: PMC8936911
## https://doi.org/10.1007%2Fs10489-021-02635-5

f1scores <- function(mat, conf.level=0.95){
  ## This function computes point estimates and (conf.level*100%) confidence intervals
  ## for microF1, macroF1, and macroF1* scores.
  
  ## mat is an r by r matrix (confusion matrix).
  ## Rows indicate the predicted (fitted) conditions,
  ## and columns indicate the truth.
  ## miF1 is micro F1
  ## maF1 is macro F1
  ## maF2 is macro F1* (Sokolova and Lapalme)
  
  ## ###### ##
  ## Set up ##
  ## ###### ##
  r <- 1
  n <- sum(mat) ## Total sample size
  p <- mat/n ## probabilities
  pii <- diag(p)[2]
  pi. <- rowSums(p)[2]
  p.i <- colSums(p)[2]
  
  ## ############### ##
  ## Point estimates ##
  ## ############### ##
  miP <- miR <- sum(pii) ## MICRO precision, recall
  miF1 <- miP ## MICRO F1
  F1i <- 2*pii/(pi.+p.i)
  maF1 <- sum(F1i)/r ## MACRO F1
  maR <- pii/rowSums(p)[2] ## MACRO Precision  # This part has been changed due to structure of confusion matrix tables (it's not correct mathematically)
  maP <- pii/colSums(p)[2] ## MACRO Recall     # This part has been changed due to structure of confusion matrix tables (it's not correct mathematically)
  maF2 <- 2*(maP*maR)/(maP+maR) ## MACRO F1*
  
  ## ################## ##
  ## Variance estimates ##
  ## ################## ##
  
  ## ----------------- ##
  ## MICRO F1 Variance ##
  ## ----------------- ##
  miF1.v <- sum(pii)*(1-sum(pii))/n
  miF1.s <- sqrt(miF1.v)
  
  ## ----------------- ##
  ## MACRO F1 Variance ##
  ## ----------------- ##
  a <- sum(F1i*(pi.+p.i-2*pii)/(pi.+p.i)^2 * ((pi.+p.i-2*pii)/(pi.+p.i)+F1i/2))
  b <- 0
  
  for(i in 1:r){
    jj <- (1:r)[-i]
    for(j in jj){
      b <- b+ p[i,j]*F1i[i]*F1i[j]/((pi.[i]+p.i[i])*(pi.[j]+p.i[j]))
    }}
  maF1.v <- 2*(a+b)/(n*r^2)
  maF1.s <- sqrt(maF1.v)
  
  ## ------------------ ##
  ## MACRO F1* Variance ##
  ## ------------------ ##
  varmap <- pii*(pi.-pii)/pi.^3 / r^2 / n
  varmar <- (pii*(p.i-pii)/p.i^3) / r^2 / n
  covmpr1 <- ( ((pi.-pii) * pii * (p.i-pii)) / (pi.^2 * p.i^2) )
  covmpr2 <- 0
  for(i in 1:r){
    covmpr2 <- covmpr2 + sum(pii[i] * p[i,-i] * pii[-i] / pi.[i]^2 / p.i[-i]^2)
  }
  covmpr <- (covmpr1+covmpr2) / r^2 / n
  maF2.v <- 4 * (maR^4*varmap + 2*maP^2*maR^2*covmpr + maP^4*varmar) / (maP+maR)^4
  maF2.s <- sqrt(maF2.v)
  
  ## #################### ##
  ## Confidence intervals ##
  ## #################### ##
  z <- qnorm(1-(1-conf.level)/2)
  miF1.ci <- miF1 + c(-1,1)*z*miF1.s
  maF1.ci <- maF1 + c(-1,1)*z*maF1.s
  maF2.ci <- maF2 + c(-1,1)*z*maF2.s
  
  ## ################# ##
  ## Formattnig output ##
  ## ################# ##
  pr <- data.frame(microPrecision=miP, microRecall=miR, macroPrecision=maP, macroRecall=maR)
  fss <- data.frame(
    rbind(miF1=c(miF1, miF1.s, miF1.ci),
          maF1=c(maF1, maF1.s, maF1.ci),
          maF1.star=c(maF2, maF2.s, maF2.ci)))
  names(fss) <- c('PointEst','Sd','Lower','Upper')
  out <- list(pr, fss)
  names(out) <- c('Precision.and.Recall','Confidence.Interval')
  out
}
