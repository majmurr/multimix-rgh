Z.to.P <- function(Z, D, P) {
  with( c(D, P), {
    ppi <- colSums(Z)/n
    W <- Z%*%diag(1/{n*ppi}) # Z scaled to have columns sum to 1 
    for (i in seq_along(dlink)) dstat[[i]] <- crossprod(W, dvals[[i]])
    for (i in seq_along(lcdisc)) ldstat[[i]] <- crossprod(W, ldvals[[i]])
    ostat <- crossprod(W, ovals)
    for (j in seq_len(q)) {
      ovar[j,] <- diag(varw(ovals, W[,j]))
    }
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
        lcqstat[[i]][[j]][lv, ] <- lcstat[[i]][[lv]][j, ] 						
        Temp <- Temp + varw(lcvals[[i]][group, , drop=FALSE], WW[,j])*
                       ldstat[[i]][j,lv]
      } #lv  
      LMV[[i]][[j]] <- Temp
    } #j 
  } #i
P <- list(dstat = dstat,
            ldstat = ldstat,
            ostat = ostat,
			ovar = ovar,
            pistat = pistat,
            cstat = cstat,
            lcstat = lcstat,
            MVMV = MVMV,
            LMV = LMV,
	     W = W)
  return(P)
  })		 
}



