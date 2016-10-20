########## R function: twoLevGauss ##########

# For illustration of Fortranisation of R code
# for multilevel model fitting via MFVB.

# Last changed: 17 FEB 2014
#############################################

twoLevGauss <- function(y,XG,XR,ZG,ZR=NULL,reBlockInds,ncZG,
                        doStreamlined=TRUE,
                        maxIter=35)
{
   # Obtain constant matrices and dimension variables:

   nuVal <- 2
   X <- cbind(XR,XG) ; CG <- cbind(X,ZG) 
   ncX <- ncol(X) ; ncXR <- ncol(XR)
   ncCG <- ncol(CG)
   m <- length(reBlockInds)
   nVec <- unlist(lapply(reBlockInds,length))
   numObs <- sum(nVec)
   L <- length(ncZG)
   
   # Compute constant matrices:

   if (!doStreamlined)
   {
      C <- cbind(CG,ZR) ; CTC <- crossprod(C)
      CTy <- crossprod(C,y) ; ncC <- ncol(C)
   }
       
   if (doStreamlined)
   {
      CGTCG <- crossprod(CG) ; CGTy <- crossprod(CG,y)
   }
   
   # Determine indsStt and indsEnd arrays:

   firstEntry <- function(x) return(x[1])
   lastEntry <- function(x) return(x[length(x)])
   
   indsStt <- lapply(reBlockInds,firstEntry)
   indsEnd <- lapply(reBlockInds,lastEntry)

   # Determine indsStt and indsEnd arrays for the
   # spline components:
   
   if (L == 1)
   {
      indsSplStt <- ncol(X) + 1
      indsSplEnd <- indsSplStt + ncZG - 1
   }

   if (L > 1)
   {
      indsSplStt <- rep(0,L); indsSplEnd <- rep(0,L)
      indsSplStt[1] <- ncol(X) + 1
      for (ell in 1:L)
      {
         indsSplEnd[ell] <- indsSplStt[ell] + ncZG[ell] - 1
         indsSplStt[L] <- indsSplEnd[1] + 1
      }
   }
   
   if (!doStreamlined)
   {
      muqrecssqu <- rep(1,L)
      MqinvSigR <- diag(ncXR)
      MqinvSig <- matrix(0,ncC,ncC)
      muqrecsseps <- 1
      det <- rep(0,2)
      ipvtBig <- rep(0,ncC)
      workBig <- matrix(0,ncC,ncC)
      Sigqbetau <- matrix(0,ncC,ncC)
      muqbetau <- rep(0,ncC)
      residSS <- 0
      trTerm <- 0
      muqrecaeps <- 1
      Aqsigsqeps <- 0
      Bqsigsqeps <- 0
      Bqaeps <- 0
      BqaR <- rep(0,ncXR)
      muqrecaR <- rep(0,ncXR)
      BqSigmaR <- matrix(0,ncXR,ncXR)
      ipvt <- rep(0,ncXR)
      work <- matrix(0,ncXR,ncXR)
      Bqau <- rep(0,L)
      muqrecau <- rep(0,L)
      muqbetauG <- rep(0,ncCG)
      SigqbetauG <- matrix(0,ncCG,ncCG)
      Bqsigsqu <- rep(0,L)
      Aqsigsqu <- rep(0,L)
      logMLgrid <- rep(0,maxIter)
      lgamAval <- rep(0,ncXR)
      lgamAmVal <- rep(0,ncXR)
      lgamZG <- lgamma(0.5*(ncZG+1))
      lgamXR <- lgamma(0.5*(nuVal+ncXR))
      lgamObs <- lgamma(0.5*(numObs+1))

      ans <- .Fortran("twlvgfNA",as.integer(numObs),as.integer(ncXR),
                      as.integer(ncX),as.integer(L),as.integer(ncZG),as.integer(m),
                      as.integer(nVec),as.integer(ncCG),as.integer(ncC),as.double(CTC),
                      as.double(CTy),as.double(C),as.double(y),as.integer(indsStt),
                      as.integer(indsEnd),as.double(indsSplStt),as.double(indsSplEnd),
                      muqrecssqu=as.double(muqrecssqu),MqinvSigR=as.double(MqinvSigR),
                      MqinvSig=as.double(MqinvSig),as.double(det),as.integer(ipvt),
                      as.double(work),as.integer(ipvtBig),as.double(workBig),
                      muqrecsseps=as.double(muqrecsseps),Sigqbetau=as.double(Sigqbetau),
                      muqbetau=as.double(muqbetau),residSS=as.double(residSS),
                      trTerm=as.double(trTerm),muqrecaeps=as.double(muqrecaeps),
                      Aqsigsqeps=as.double(Aqsigsqeps),Bqsigsqeps=as.double(Bqsigsqeps),
                      Bqaeps=as.double(Bqaeps),BqaR=as.double(BqaR),
                      muqrecaR=as.double(muqrecaR),BqSigmaR=as.double(BqSigmaR),
                      Bqau=as.double(Bqau),muqrecau=as.double(muqrecau),
                      Bqsigsqu=as.double(Bqsigsqu),Aqsigsqu=as.double(Aqsigsqu),
                      maxIter=as.integer(maxIter),logMLgrid=as.double(logMLgrid),
                      lgamAval=as.double(lgamAval),lgamAmVal=as.double(lgamAmVal),
                      lgamZG=as.double(lgamZG),lgamXR=as.double(lgamXR),
                      lgamObs=as.double(lgamObs))

       return(list(mu.q.betau=ans$muqbetau,
                   Sigma.q.betau=matrix(ans$Sigqbetau,ncC,ncC),
                   A.q.sigsq.eps=ans$Aqsigsqeps,B.q.sigsq.eps=ans$Bqsigsqeps,
                   logMLgrid=ans$logMLgrid))
    }
   
   if (doStreamlined)
   {
      # Create matrix in which to store ZRTy values:
   
      ZRTy <- matrix(0,m,ncXR)
   
      # Create G and H matrix arrays required for MFVB,
      # as well as working arrays required by LINPACK:

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
      muqbetauG <- rep(0,ncCG)
      SigquR <- array(0,c(ncXR,ncXR,m))
      zetaVec <- array(0,c(m,ncXR))
      SigqsVec <- rep(0,ncCG)
      kappaVec <- array(0,c(m,ncXR))
      SigqCGTy <- rep(0,ncCG)
      HtimesKappa <- array(0,c(m,ncXR))
      muquR <- array(0,c(m,ncXR))
      Asigsqeps <- 0
      Bqsigsqeps <- 0
      muqrecsseps <- 1
      Bqaeps <- 0
      muqrecaeps <- 1
      MqinvSigR <- diag(ncXR)
      BqaR <- rep(0,ncXR)
      muqrecaR <- rep(1,ncXR)
      BqSigmaR <- matrix(0,ncXR,ncXR)
      muqrecssqu <- rep(1,ncXR)      
      Bqau <- rep(0,L)
      muqrecau <- rep(0,L)
      Bqsigsqu <- rep(0,L)
      Aqsigsqu <- rep(0,L)
      logMLgrid <- rep(0,maxIter)
      lgamZG <- lgamma(0.5*(ncZG+1))
      lgamXR <- lgamma(0.5*(nuVal+ncXR))
      lgamObs <- lgamma(0.5*(numObs+1))
      
      ans <- .Fortran("twlvgfSL",as.double(y),as.double(XR),
                   as.double(CG),as.double(CGTCG),as.double(CGTy),
                   as.integer(numObs),as.integer(m),as.integer(L),
                   as.integer(ncXR),as.integer(ncX),as.integer(ncCG),
                   as.integer(indsStt),as.integer(indsEnd),as.integer(indsSplStt),
                   as.integer(indsSplEnd),as.integer(ncZG),ZRTy=as.double(ZRTy),
                   G=as.double(G),H=as.double(H),Hwork=as.double(Hwork),as.integer(ipvt),
                   as.integer(ipvtBig),as.double(det),as.double(work),
                   as.double(workBig),sVec=as.double(sVec),
                   Smat=as.double(Smat),SigqbetauG=as.double(SigqbetauG),
                   muqbetauG=as.double(muqbetauG),SigquR=as.double(SigquR),
                   zetaVec=as.double(zetaVec),kappaVec=as.double(kappaVec),
                   muquR=as.double(muquR),Aqsigsqeps=as.double(Asigsqeps),
                   Bqsigsqeps=as.double(Bqsigsqeps),muqrecsseps=as.double(muqrecsseps),
                   Bqaeps=as.double(Bqaeps),muqrecaeps=as.double(muqrecaeps),
                   MqinvSigR=as.double(MqinvSigR),BqaR=as.double(BqaR),
                   muqrecaR=as.double(muqrecaR),BqSigmaR=as.double(BqSigmaR),
                   muqrecssqu=as.double(muqrecssqu),Bqau=as.double(Bqau),
                   muqrecau=as.double(muqrecau),Bqsigsqu=as.double(Bqsigsqu),
                   Aqsigsqu=as.double(Aqsigsqu),maxIter=as.integer(maxIter),
                   logMLgrid=as.double(logMLgrid),lgamZG=as.double(lgamZG),
                   lgamXR=as.double(lgamXR),lgamObs=as.double(lgamObs))

      return(list(mu.q.betauG=ans$muqbetauG,Sigma.q.betauG=matrix(ans$SigqbetauG,ncCG,ncCG),
                  A.q.sigsq.eps=ans$Aqsigsqeps,B.q.sigsq.eps=ans$Bqsigsqeps,logMLgrid=ans$logMLgrid))
   }
}   

######### End of twoLevGauss ##########


