########## R function: twoLevBern ##########

# For illustration of Fortranisation of R code
# for multilevel model fitting via MFVB.

# Last changed: 09 FEB 2014
############################################

twoLevBern <- function(y,XG,XR,ZG,ZR=NULL,reBlockInds,ncZG,
                       doStreamlined=TRUE,
                       maxIter=35)
{
   # Obtain constant matrices and dimension variables:
  
   X <- cbind(XR,XG) ; yAdj <- y - 0.5
   CG <- cbind(X,ZG) ; ncX <- ncol(X) ; ncXR <- ncol(XR)
   ncCG <- ncol(CG)  
   m <- length(reBlockInds)
   nVec <- unlist(lapply(reBlockInds,length))
   numObs <- sum(nVec); L <- length(ncZG)
      
   # Compute constant matrices:
   
   if (!doStreamlined)
   {
      C <- cbind(CG,ZR) ; CTC <- crossprod(C)
      CTy <- crossprod(C,yAdj) ; ncC <- ncol(C)
   }
       
   if (doStreamlined)
   {
      CGTCG <- crossprod(CG) ; CGTy <- crossprod(CG,yAdj)
   }
   
   # Determine indsStt and indsEnd arrays:

   firstEntry <- function(x) return(x[1])
   lastEntry <- function(x) return(x[length(x)])
   
   indsStt <- lapply(reBlockInds,firstEntry)
   indsEnd <- lapply(reBlockInds,lastEntry)
   
   # Determine indsStt and indsEnd arrays for the
   # spline components:

   indsSplStt <- rep(0,L); indsSplEnd <- rep(0,L)
   indsSplStt[1] <- ncol(X) + 1
   for (ell in 1:L)
   {
      indsSplEnd[ell] <- indsSplStt[ell] + ncZG[ell] - 1
      indsSplStt[L] <- indsSplEnd[1] + 1                                              
   }

   if (!doStreamlined)
   {
      wtVec <- rep(0,numObs)
      xiVec <- rep(1,numObs)
      muqrecssqu <- rep(1,L)
      MqinvSigR <- diag(ncXR)
      MqinvSig <- matrix(0,ncC,ncC)
      Sigqbetau <- matrix(0,ncC,ncC)
      det <- rep(0,2)
      ipvtBig <- rep(0,ncC)
      workBig <- matrix(0,ncC,ncC)
      muqbetau <- rep(0,ncC)
      xiSqd <- rep(0,numObs)
      EsqMat <- matrix(0,numObs,ncC)
      BqaR <- rep(0,ncXR)
      muqrecaR <- rep(0,ncXR)
      BqSigmaR <- matrix(0,ncXR,ncXR)
      ipvt <- rep(0,ncXR)
      work <- matrix(0,ncXR,ncXR)
      Bqau <- rep(0,L)
      muqrecau <- rep(0,L)
      Bqsigsqu <- rep(0,L)
      Aqsigsqu <- rep(0,L)
      logMLgrid <- rep(0,maxIter)
      lgamAval <- rep(0,ncXR)
      lgamAmval <- rep(0,ncXR)
      lgamZG <- lgamma(0.5*(ncZG+1))
      lgamXR <- lgamma(0.5*(nuVal+ncXR))
      lgamObs <- lgamma(0.5*(numObs+1))
      
      ans <- .Fortran("twlvbfNA",as.integer(numObs),as.integer(ncXR),
                      as.integer(ncX),as.integer(L),as.integer(ncZG),as.integer(m),
                      as.integer(ncCG),as.integer(ncC),as.integer(indsStt),
                      as.integer(indsEnd),as.double(det),as.double(ipvtBig),
                      as.double(workBig),as.double(ipvt),as.double(work),
                      as.double(CTy),wtVec=as.double(wtVec),xiVec=as.double(xiVec),
                      EsqMat=as.double(EsqMat),as.double(xiSqd),
                      as.double(muqrecssqu), MqinvSigR=as.double(MqinvSigR),
                      MqinvSig=as.double(MqinvSig),C=as.double(C),
                      Sigqbetau=as.double(Sigqbetau),muqbetau=as.double(muqbetau),
                      BqaR=as.double(BqaR),muqrecaR=as.double(muqrecaR),
                      BqSigmaR=as.double(BqSigmaR),Bqau=as.double(Bqau),
                      muqrecau=as.double(muqrecau),Bqsigsqu=as.double(Bqsigsqu),
                      Aqsigsqu=as.double(Aqsigsqu),maxIter=as.integer(maxIter),
                      logMLgrid=as.double(logMLgrid),lgamAval=as.double(lgamAval),
                      lgamAmval=as.double(lgamAval),lgamZG=as.double(lgamZG),
                      lgamXR=as.double(lgamXR),lgamObs=as.double(lgamObs))

      return(list(mu.q.betauG=ans$muqbetau,Sigma.q.betauG=matrix(ans$Sigqbetau,ncC,ncC),
                  logMLgrid=ans$logMLgrid))             
   }
   
   if (doStreamlined)
   {
      # Create matrix in which to store ZRTy values:
   
      ZRTy <- matrix(0,m,ncXR)
   
      # Create G and H matrix arrays required for MFVB,
      # as well as working arrays required by LINPACK:

      nVecNew <- c(0,cumsum(nVec))[1:m]
      xiVec <- rep(1,numObs)
      wtVec <- rep(1,numObs)
      G <- array(0,c(ncCG,ncXR,m))
      H <- array(0,c(ncXR,ncXR,m))
      Hwork <- matrix(0,ncXR,ncXR)
      det <- rep(0,2)
      ipvt <- rep(0,ncXR)
      work <- matrix(0,ncXR,ncXR)
      ipvtBig <- rep(0,ncCG)
      workBig <- matrix(0,ncCG,ncCG)
      sVec <- rep(0,ncCG)
      Smat <- matrix(0,ncCG,ncCG)
      SigqbetauG <- matrix(0,ncCG,ncCG)
      ridgeVec <- rep(0,ncCG)
      muqbetauG <- rep(0,ncCG)
      HtimestG <- array(0,c(ncXR,ncCG,m))
      SigquR <- array(0,c(ncXR,ncXR,m))
      zetaVec <- array(0,c(m,ncXR))
      SigqsVec <- rep(0,ncCG)
      kappaVec <- array(0,c(m,ncXR))
      SigqCGTy <- rep(0,ncCG)
      HtimesKappa <- array(0,c(m,ncXR))
      muquR <- array(0,c(m,ncXR))
      MqinvSigR <- diag(ncXR)
      BqaR <- rep(0,ncXR)
      muqrecaR <- rep(1,ncXR)
      BqSigmaR <- matrix(0,ncXR,ncXR)
      muqrecssqu <- rep(1,ncXR)
      Bqau <- rep(0,L)
      muqrecau <- rep(0,L)
      Bqsigsqu <- rep(0,L)
      Aqsigsqu <- rep(0,L)      
      diagMidTermOne <- matrix(0,numObs,ncCG)
      diagTermOne <- rep(0,numObs)      
      diagMidTermTwo <- array(0,c(ncCG,ncXR,m))
      diagMidTwoSqd <- matrix(0,numObs,ncXR)
      diagTermTwo <- rep(0,numObs)      
      diagMidTermThree <- array(0,c(ncXR,ncXR,m)) 
      diagMidThreeSqd <- matrix(0,numObs,ncXR)
      diagTermThree <- rep(0,numObs)      
      xiSqd <- rep(0,numObs)
      logMLgrid <- rep(0,maxIter)
      lgamZG <- lgamma(0.5*(ncZG+1))
      lgamXR <- lgamma(0.5*(nuVal+ncXR))
      lgamObs <- lgamma(0.5*(numObs+1))
      
      ans <- .Fortran("twlvbfSL",as.double(yAdj),as.double(XR),
                      as.double(CG),as.double(CGTCG),as.double(CGTy),
                      as.integer(numObs),as.integer(m),as.integer(L),as.integer(nVec),
                      as.integer(ncXR),as.integer(ncX),as.integer(ncCG),as.integer(nVecNew),
                      as.integer(indsStt),as.integer(indsEnd),as.integer(indsSplStt),
                      as.integer(indsSplEnd),as.integer(ncZG),ZRTy=as.double(ZRTy),
                      xiVec=as.double(xiVec),wtVec=as.double(wtVec),
                      G=as.double(G),H=as.double(H),Hwork=as.double(Hwork),as.integer(ipvt),
                      as.integer(ipvtBig),as.double(det),as.double(work),
                      as.double(workBig),sVec=as.double(sVec),Smat=as.double(Smat),
                      SigqbetauG=as.double(SigqbetauG),muqbetauG=as.double(muqbetauG),
                      SigquR=as.double(SigquR),zetaVec=as.double(zetaVec),
                      kappaVec=as.double(kappaVec),muquR=as.double(muquR),
                      diagTermOne=as.double(diagTermOne),diagTermTwo=as.double(diagTermTwo),
                      diagTermThree=as.double(diagTermThree),
                      xiSqd=as.double(xiSqd),MqinvSigR=as.double(MqinvSigR),BqaR=as.double(BqaR),
                      muqrecaR=as.double(muqrecaR),BqSigmaR=as.double(BqSigmaR),
                      muqrecssqu=as.double(muqrecssqu),Bqau=as.double(Bqau),
                      muqrecau=as.double(muqrecau),Bqsigsqu=as.double(Bqsigsqu),
                      Aqsigsqu=as.double(Aqsigsqu),maxIter=as.integer(maxIter),
                      logMLgrid=as.double(logMLgrid),lgamZG=as.double(lgamZG),
                      lgamXR=as.double(lgamXR),lgamObs=as.double(lgamObs))

      return(list(mu.q.betauG=ans$muqbetauG,Sigma.q.betauG=matrix(ans$Sigqbetau,ncCG,ncCG),
                  logMLgrid=ans$logMLgrid))
   }
}   

######### End of twoLevBern ##########


