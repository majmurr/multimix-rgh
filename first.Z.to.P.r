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
  ###
  lcstat <- replicate(length(lcdep), list())
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
    ovar <- ostat
    for (j in seq_len(q)) {
      ovar[j,] <- diag(varw(ovals, W[,j]))
    }
    # ns
    pistat <- ppi
    for (i in seq_along(cdep) ) {
      cstat[[i]] <- crossprod(W, cvals[[i]])
     for (j in 1:q) {
        MVMV[[i]][[j]] <- varw(cvals[[i]], W[,j])
     }	
    }
        
  
  for (i in seq_along(lcdep) ) {
    lcdi <- length(lcdep[[i]]) - 1
    for (j in 1:q) {
      Temp <- 0*diag(lcdi)
      for (lv in seq_len(ldlevs[i])) {
        group <- ldvals[[i]][,lv]==1
        gtot <- colSums(W[group,])   
        ### or maybe pmin(colsums(W[group,]), minpstar)
        WW <- W[group,]%*%diag(1/gtot)
        lcstat[[i]][[lv]] <- crossprod(WW,
                        lcvals[[i]][group, , drop=FALSE])
        Temp <- Temp + varw(lcvals[[i]][group, , drop=FALSE], WW[,j])*
                       ldstat[[i]][j,lv]
      } #lv  
      LMV[[i]][[j]] <- Temp
    } #j 
  } #i
P <- list(dstat = dstat,
            ldstat = ldstat,
            ostat = ostat,
            #ostat2 = ostat2,
            #ovar = ovar,
            pistat = pistat,
            cstat = cstat,
            #cstat2 = cstat2,
            #cvar = cvar,
            #cpstat = cpstat,
            lcstat = lcstat,
            #lcstat2 = lcstat2,
            #lcpstat = lcpstat,
            MVMV = MVMV,
            LMV = LMV,
	     W = W)
  return(P)
  })		 
}
