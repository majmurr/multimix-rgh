P.to.Z <- function(P, D) {
    with(c(P, D), {
        ollq <- matrix(0, nrow=n, ncol=q)
        for(j in 1:q){
          ldens <- dnorm(as.vector(ovals),rep(ostat[j,],rep(n,op)), 
                    rep(sqrt(ovar)[j,],rep(n,op)), log = TRUE)
          lDENS <- matrix(ldens, nrow = n)
          ollq[,j] <- rowSums(lDENS) 	
        }	  

        cll <- matrix(0, nrow=n, ncol=q)
        for(j in 1:q){
          cno <- 0
          while(cno < length(cdep)) {
            cno <- cno + 1
            cll[,j] <- cll[,j] +
            dmvnorm(cvals[[cno]], mean = cstat[[cno]][j,], 
            sigma = MVMV[[cno]][[j]], log = TRUE)
          }
        }

        dllq <- matrix(0, nrow=n, ncol=q)  
        for(j in 1:q){
          for(v in seq_along(dstat)){
            for(k in seq_len(dlevs[v])){
              dvk <- dvals[[v]][,k]  
              dllq[,j] <- dllq[,j] +
                dvk*log(dstat[[v]][j,k] + !dvk)
            }
          }
        }	 

        ldllq <- matrix(0, nrow=n, ncol=q)  
        for(j in 1:q){
          for(v in seq_along(ldstat)){
            for(k in seq_len(ldlevs[v])){
			  ldvk <- ldvals[[v]][,k]
        	  ldllq[,j] <- ldllq[,j] +
        	  ldvk*log(ldstat[[v]][j,k] + !ldvk)
              }
        	 }
        }	 
        
        for (cno in seq_along(lcdep)) {
            for (j in 1:q) {
                lcqexp[[cno]][[j]] <- lcqstat[[cno]][[j]][ldvnum[[cno]], ]
            }
        }

        lcll <- matrix(0, nrow=n, ncol=q)
        for (j in 1:q) {
            cno <- 0
            w <- W[,j]
            while(cno < length(lcdep)) {
                cno <- cno + 1
                lcqexp[[cno]][[j]] <- lcqstat[[cno]][[j]][ldvnum[[cno]], ]
                lcll[,j] <- lcll[,j] +
                    dmvnorm(lcvals[[cno]] - lcqexp[[cno]][[j]],
                            mean = rep(0,dim(lcvals[[cno]])[2]), sigma = LMV[[cno]][[j]], log = TRUE)
            }
        }

        llx <- ollq + cll + dllq + ldllq + lcll
        rmx <- apply(llx, 1, max)
        expld <- exp(llx - rmx) 
        pf_ <- t(expld)*pistat
        pstar <- colSums(pf_)
        pf_[ ,pstar < minpstar ] <- 1/q
        pstar[pstar < minpstar ] <- 1
        Z <- t(pf_)/pstar
        llik <- sum( log(pstar) + rmx )	
        zll <- list(Z = Z, llik = llik)
        return(zll)	
	})	
}





