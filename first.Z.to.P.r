first.Z.to.P <- function(D) {
 with(D, {
  Z <- model.matrix(~ 0 + Init_grp)
  attr(Z,"assign") <- NULL
  attr(Z,"contrasts") <- NULL 
  colnames(Z) <- NULL
  # Set up store structures for sufficient statistics.
  # Values given not important.
  dstat <- vector("list", length(dlink))
  ldstat <- vector("list", length(ldlink))
  cstat <- vector("list", length(cdep))
  cstat2 <- vector("list", length(cdep))
  cvar <- vector("list", length(cdep))
  cpstat <- vector("list", length(cdep))
  ccov <- vector("list", length(cdep))
  ###
  lcstat <- lcstat2 <- lcvar <- lcpstat <- lccov <-
    replicate(length(lcdep), list())
  ###
  MVMV <- list()
  for (i in seq_along(cdep) ) {
    MVMV[[i]] <- list()
    for (j in seq_len(q)) {
      MVMV[[i]][[j]] <- diag(length(cdep[[i]]))
    }}
  LMV <- list()
  for (i in seq_along(lcdep) ) {
    LMV[[i]] <- list()
    for (j in seq_len(q)) {
      LMV[[i]][[j]] <- diag(length(lcdep[[i]]) - 1)
    }}
  results <- matrix(0,nrow=nIt+1, ncol=q + 2)	
  ppi <- colSums(Z)/n
  W <- Z%*%diag(1/{n*ppi})

    for (i in seq_along(dlink)) dstat[[i]] <- crossprod(W, dvals[[i]])
    for (i in seq_along(lcdisc)) ldstat[[i]] <- crossprod(W, ldvals[[i]])
    ostat <- crossprod(W, ovals)
    #ostat2 <- crossprod(W, ovals2)
    #ovar <- ostat2 - ostat^2
    ovar <- ostat
    for (j in seq_len(q)) {
      ovar[j,] <- diag(varw(ovals, W[,j]))
    }
    # ns
    pistat <- ppi
    for (i in seq_along(cdep) ) {
      cstat[[i]] <- crossprod(W, cvals[[i]])
      # cstat2[[i]] <- crossprod(W, cvals2[[i]])
      # cpstat[[i]] <- crossprod(W, cprods[[i]])
    }
    for (i in seq_along(lcdep)) {
      for (lv in seq_len(ldlevs[i])){
        group <- ldvals[[i]][,lv]==1
        gtot <- colSums(W[group,])   
        ### or maybe pmin(colsums(W[group,]), minpstar)
        WW <- W[group,]%*%diag(1/gtot)
        lcstat[[i]][[lv]] <- crossprod(WW,lcvals[[i]][group,])
        #lcstat2[[i]][[lv]] <- crossprod(WW,lcvals2[[i]][group,])
        #lcpstat[[i]][[lv]] <- crossprod(WW,lcprods[[i]][group,])
      }
    }
  
#  cvar <- cstat
#  ccov <- cpstat
#  for (i in seq_along(cdep) ) {
#     lcdi <- length(cdep[[i]])
#    nxp <- lcdi*(lcdi - 1)/2  
#    cvar[[i]] <-  cstat2[[i]] - {cstat[[i]]}^2
#    for (j in 1:q) {
#      for (k in  seq_along(cdep[[i]])) {
#        MVMV[[i]][[j]][k, k] = cvar[[i]][j,k]
#      }
#    }	
#    for (ii in seq_len(nxp)) {
#      ccov[[i]][ ,ii] <- (cpstat[[i]][ ,ii] - 
#                            cstat[[i]][ ,left(ii)]*cstat[[i]][ ,right(ii)])
#    }
#    for (j in 1:q) {
#      for (ii in seq_len(nxp)) {
#        MVMV[[i]][[j]][left(ii), right(ii)] <-
#          MVMV[[i]][[j]][right(ii), left(ii)] <- ccov[[i]][j, ii]
#      }
#    }
#  }
#            
    for (j in 1:q) {
      for (k in  seq_along(cdep[[i]])) {
        MVMV[[i]][[j]] <- varw(cvals[[i]], W[,j])
      }
    }	
# ns    
    
  lcvar <- lcstat
  lccov <- lcpstat
  for (i in seq_along(lcdep) ) {
    lcdi <- length(lcdep[[i]]) - 1
    for (j in 1:q) {
      Temp <- 0*diag(lcdi)
      for (lv in seq_len(ldlevs[i])) {
        group <- ldvals[[i]][,lv]==1
        gtot <- colSums(W[group,])   
        ### or maybe pmin(colsums(W[group,]), minpstar)
        WW <- W[group,]%*%diag(1/gtot)
        Temp <- Temp + varw(lcvals[[i]][group,], WW[,lv])*ldstat[j,lv]
      } #lv  
      LMV[[i]][[j]] <- Temp
    } #j 
  } #i
  P <- list(dstat = dstat,
            ldstat = ldstat,
            ostat = ostat,
            ostat2 = ostat2,
            ovar = ovar,
            pistat = pistat,
            cstat = cstat,
            cstat2 = cstat2,
            cvar = cvar,
            cpstat = cpstat,
            lcstat = lcstat,
            lcstat2 = lcstat2,
            lcpstat = lcpstat,
            MVMV = MVMV,
            LMV = LMV,
			W = W)
  return(P)
  })		
 #  return(P)
}