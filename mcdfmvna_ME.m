#library(gtools)
library(pbivnorm)

#sys_randper = 1;


mpdfmvna_ME <- function(a,r,s)
{
  K = length(a); # dimension of CDF.
  
  a = 	a * (a < 5.7) + 5.7*(a >= 5.7); # truncate a at 5.7.
  w =  order(a,decreasing=TRUE)
  
  Zj = t(a[w]);
  Rij = r[w,w];
  
  # prepare for derivatives
  npar = K+K*(K-1)/2
  dRij = array(0,dim = c(K,K,npar))
  dZj = matrix(0,K,npar)
  dlp = matrix(0,1,npar)
  dtZ = dZj;
  dtRij = dRij;
  dsig = dZj;
  dajjm1 = matrix(0,1,npar)
  
  fm = lower.tri(diag(K),diag=FALSE)
  rRowCol = which(fm ==TRUE, arr.in = TRUE)
  for (p in 1:K){
    dZj[p,p]= 1
  }
  for(p in (K+1):npar){
    dRij[rRowCol[p-K,1],rRowCol[p-K,2],p]=1
    dRij[rRowCol[p-K,2],rRowCol[p-K,1],p]=1
  }
  lp = log(pnorm(Zj[1])); # perform all calcs in logs.
  for (p in 1:npar) {
    dlp[p] = dnorm(Zj[1])/pnorm(Zj[1])*dZj[1,p];
  }
  for (j in 1:(K-1)){
    # update a_{j|j-1}
    ajjm1 = dnorm(Zj[j])/pnorm(Zj[j]) # a_{j|j-1}
    # deriv 
    for (p in 1:npar) {
      dajjm1[p] = -dZj[j,p] * (ajjm1*Zj[j]+ajjm1^2)
    }
    # update Zj 
    tZ = Zj + ajjm1*Rij[,j]        
    # deriv
    for (p in 1:npar) {
      dtZ[,p] = dZj[,p] + dajjm1[p]*Rij[,j] + ajjm1*dRij[,j,p]
    }

    # update Rij
    tRij = Rij - (Rij[,j] %o% Rij[j,])*(ajjm1+Zj[j]) *ajjm1
    for (p in 1:npar) {
      dtRij[,,p] = dRij[,,p] - (ajjm1+Zj[j]) *ajjm1 * (dRij[,j,p] %o% Rij[j,] + Rij[,j] %o% dRij[j,,p]) - (Rij[,j] %o% Rij[j,]) * (2*ajjm1 * dajjm1[p]+ Zj[j]*dajjm1[p] + ajjm1 *dZj[j,p])  
    }
    # adjust variances
    sig = sqrt(diag(tRij)) # variances
    for (p in 1:npar) {
      dsig[,p] = 0.5 * diag(dtRij[,,p])/sig
    }

    # normalize via std.
    Zj = tZ/sig
    # deriv 
    for (p in 1:npar) {
      dZj[,p] = dtZ[,p]/sig - tZ/sig^2*dsig[,p]
    }

    # normalize Rij. 
    Rij = t(tRij/sig)/sig
    diag(Rij)=1    

    # deriv
    for (p in 1:npar) {
      dRij[,,p] = t(dtRij[,,p]/sig)/sig - t(tRij/sig^2*dsig[,p]) - t(tRij/sig)/sig^2*dsig[,p]
    }
    
    lp = lp + log(pnorm(Zj[j+1]))
    # deriv
    for (p in 1:npar) {
      dlp[p] = dlp[p] + dnorm(Zj[j+1])/pnorm(Zj[j+1])*dZj[j+1,p]
    }
  }
  
  # resort in original order
  aaa = cbind(w,(1:K));
  aaa = aaa[sort.list(aaa[,1]),2];
  dlp1 = dlp[1:K]
  dlp1 = t(dlp1[aaa]);
  
  res = matrix(nrow=0,ncol=2);
  for(ir in (1:(K-1)))
  {
    jr = ir+1;
    temp = matrix((rep(w[ir],K-ir)),ncol=1);     
    temp = cbind(temp,matrix(w[jr:K],nrow=K-jr+1,ncol=1));          
    res = rbind(res,temp);                              
  }
  res2 = cbind(apply(res,1,min),apply(res,1,max));
  
  res2 = matrix(res2[,1] *(10^(floor (log(res2[,2])/log(10))+1)) + res2[,2],ncol=1);  
  
  aaa1 = cbind(res2,(1:nrow(res2)));
  aaa1 = matrix(aaa1[sort.list(aaa1[,1]),2],ncol=1);
  dlp2 = dlp[(K+1):npar]
  dlp2 = t(dlp2[aaa1])

  # return values.
  return(cbind(lp,dlp1,dlp2,2343));
}





