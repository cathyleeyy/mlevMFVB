cccccccccc FORTRAN subroutine twlvbfNA.f cccccccccc

c For performing variational Bayes inference for a
c Bernoulli random effects model (naive approach).

c The following abbreviations are used with repect
c to the parent R code:
c
c "ssq" for "sigsq", "rec" for "recip".

c Last changed: 17 FEb 2014
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine twlvbfNA(numObs,ncXR,ncX,L,ncZG,m,ncCG,
     +                    ncC,indsStt,indsEnd,det,ipvtBig,workBig,
     +                    ipvt,work,CTy,wtVec,xiVec,EsqMat,xiSqd,
     +                    muqrecssqu,MqinvSigR,MqinvSig,C,
     +                    Sigqbetau,muqbetau,BqaR,muqrecaR,
     +                    BqSigmaR,Bqau,muqrecau,Bqsigsqu,Aqsigsqu,
     +                    maxIter,logMLgrid,lgamAval,lgamAmval,
     +                    lgamZG,lgamXR,lgamObs)
      double precision sigsqBeta,nuVal,AR,Au,
     +                 muqrecssqu(L),MqinvSigR(ncXR,ncXR),
     +                 MqinvSig(ncC,ncC),wtVec(numObs),xiVec(numObs),
     +                 C(numObs,ncC),Sigqbetau(ncC,ncC),det(2),
     +                 workBig(ncC,ncC),CTy(ncC),muqbetau(ncC),
     +                 EsqMat(numObs,ncC),xiSqd(numObs),BqaR(ncXR),
     +                 muqrecaR(ncXR),BqSigmaR(ncXR,ncXR),
     +                 work(ncXR,ncXR),Bqau(L),muqrecau(L),Bqsigsqu(L),
     +                 Aqsigsqu(L),VBtoler,logMLprev,logMLcurr,ASigma,
     +                 logMLgrid(maxIter),relerr,lgamAval(ncXR),
     +                 lgamAmval(ncXR),logCqRA,logCqRAm,phiVec(numObs),
     +                 ans,sumXi,detSig,detBqSigmaR,lgamZG(L),
     +                 lgamXR,lgamObs
      integer i,ii,jj,ll,kk,k,itnum,numObs,ncXR,ncX,L,ncZG(L),m,
     +        ncCG,ncC,indsStt(m),indsEnd(m),iStt,iEnd,ipvtBig(ncC),
     +        info,ipvt(ncXR),maxIter,KF,indcvgd
  
c     Set hyperparameters:

      sigsqBeta = 1000.0
      AR = 1000.0
      Au = 1000.0

c     Set scalar starting values:

      nuVal = 2.0
      info = 0.0
      KF = 0
      VBtoler = 0.0000001

c     Perform iterations:

      logMLcurr = -1.0e20
      indcvgd = 0
      itnum = 0

1     itnum = itnum + 1

c        Set current log(ML) value to previous one:
  
         logMLprev = logMLcurr

         do 2 ii = 1,ncC
            do 3 jj = 1,ncC
               Sigqbetau(ii,jj) = 0.0
3           continue
2        continue

c        Compute wtVec:

         do 5 ii = 1,numObs
            wtVec(ii) = tanh(xiVec(ii)/2)/(4*xiVec(ii))
            phiVec(ii) = (xiVec(ii)/2) - log(1+exp(xiVec(ii))) 
     +                   + (0.25*xiVec(ii)*tanh(xiVec(ii)/2))
5        continue

c        Set M.q.inv.Sigma:

         do 10 ii = 1,ncX
            MqinvSig(ii,ii) = 1.0/sigsqBeta
10       continue

         iStt = ncX + 1
         do 20 ll = 1,L
            do 30 k = iStt,(iStt+ncZG(ll)-1)
               MqinvSig(k,k) = muqrecssqu(ll)
30          continue
            iStt = iStt + ncZG(ll)
20       continue

         iStt = ncCG + 1
         do 40 ii = 1,m
           iEnd = iStt + ncXR - 1
            do 50 k = iStt,iEnd
               do 60 kk = iStt,iEnd
                 MqinvSig(k,kk) = MqinvSigR(k-(iStt-1),kk-(iStt-1))
60             continue
50          continue
            iStt = iStt + ncXR
40       continue

c        Performa updates for q*(beta,u) parameters:
        
         do 70 i = 1,m
            do 80 ii = 1,ncC
               do 90 jj = 1,ncC
                  do 100 k = indsStt(i),indsEnd(i)
                     Sigqbetau(ii,jj) = Sigqbetau(ii,jj)
     +                        + C(k,ii)*wtVec(k)*C(k,jj)
100               continue
90             continue
80          continue
70       continue

         do 110 ii = 1,ncC
            do 120 jj = 1,ncC
               Sigqbetau(ii,jj) = 2*Sigqbetau(ii,jj)
     +                          + MqinvSig(ii,jj)
120         continue
110      continue

c        Compute determinant of Sigma.q.betau:

         call dgefa(Sigqbetau,ncC,ncC,ipvtBig,info)
         call dgedi(Sigqbetau,ncC,ncC,ipvtBig,det,workBig,10)
         detSig = -(log(abs(det(1))) + det(2)*log(10.0))

c        Invert Sigma.q.betau:

         call dgedi(Sigqbetau,ncC,ncC,ipvtBig,det,workBig,01)

         do 130 ii = 1,ncC
            muqbetau(ii) = 0.0
            do 140 jj = 1,ncC
               muqbetau(ii) = muqbetau(ii) 
     +                      + Sigqbetau(ii,jj)*CTy(jj)
140         continue
130       continue

c        Update xiVec:

         do 150 i = 1,numObs
            xiSqd(i) = 0.0
            do 160 ii = 1,ncC
               EsqMat(i,ii) = 0.0
               do 170 jj = 1,ncC
                  EsqMat(i,ii) = EsqMat(i,ii) 
     +                           + C(i,jj)*(Sigqbetau(jj,ii) 
     +                           + muqbetau(ii)*muqbetau(jj))
170            continue
               xiSqd(i) = xiSqd(i) 
     +                  + EsqMat(i,ii)*C(i,ii)
160         continue
150      continue

         do 180 i = 1,numObs
            xiVec(i) = sqrt(xiSqd(i))
180      continue
c        call dblepr("ans",3,xiVec,numObs)

c        Perform updates for q*(a.R) parameters:

         do 190 ii = 1,ncXR  
            BqaR(ii) = nuVal*MqinvSigR(ii,ii) + (1/AR**2)
            muqrecaR(ii) = 0.5*(nuVal+ncXR)/BqaR(ii)

            do 200 jj = 1,ncXR
               BqSigmaR(ii,jj) = 0.0
200         continue
190      continue

         iStt = ncCG + 1
         do 210 i = 1,m
            iEnd = iStt + ncXR - 1
            do 220 k = iStt,iEnd
               do 230 kk = iStt,iEnd
                  BqSigmaR((k-(iStt-1)),(kk-(iStt-1))) =
     +            BqSigmaR((k-(iStt-1)),(kk-(iStt-1))) + Sigqbetau(k,kk)
     +            + muqbetau(k)*muqbetau(kk)
230            continue
220         continue
            iStt = iStt + ncXR
210      continue
 
         do 240 ii = 1,ncXR
            BqSigmaR(ii,ii) = BqSigmaR(ii,ii) 
     +                      + 2*nuVal*muqrecaR(ii)
240      continue

c        Compute determinat of B.q.SigmaR:

         call dgefa(BqSigmaR,ncXR,ncXR,ipvt,info)
         call dgedi(BqSigmaR,ncXR,ncXR,ipvt,det,work,10)
         detBqSigmaR = log(abs(det(1))) + det(2)*log(10.0)

c        Invert B.q.SigmaR:

         call dgedi(BqSigmaR,ncXR,ncXR,ipvt,det,work,01)

c        Perform updates for q*(SigmaR) parameters:

         do 250 ii = 1,ncXR
            do 260 jj = 1,ncXR
               MqinvSigR(ii,jj) = (nuVal + m + ncXR - 1)*BqSigmaR(ii,jj)
260          continue
250      continue

c        Perform updates for q*(sigsq.u) and q*(a.u) parameters:

         iStt = ncX + 1
         do 270 ii = 1,L
            Bqau(ii) = muqrecssqu(ii) + (1/Au**2)
            muqrecau(ii) = 1/Bqau(ii)
            Bqsigsqu(ii) = 0.0
            iEnd = iStt + ncZG(ii) - 1
            do 280 k = iStt,iEnd
               Bqsigsqu(ii) = Bqsigsqu(ii) + muqbetau(k)**2 
     +                      + Sigqbetau(k,k)
280         continue
            iStt = iStt + ncZG(ii)
            Bqsigsqu(ii) = Bqsigsqu(ii) + 2*muqrecau(ii)
            Aqsigsqu(ii) = ncZG(ii) + 1
            muqrecssqu(ii) = Aqsigsqu(ii)/Bqsigsqu(ii)
270      continue

c     Obtain current value of log(ML):

      logMLcurr = 0.0
      ASigma = nuVal + ncXR - 1.0

c     Obtain log gamma values:

      do 290 ii = 1,ncXR
         call LGAMA(KF,0.5*(ASigma + 1.0 - ii),lgamAval(ii))
         call LGAMA(KF,0.5*(ASigma + m + 1.0 - ii),lgamAmVal(ii))
290   continue

      logCqRA = 0.5*ASigma*ncXR*log(2.0) 
     +          + 0.25*ncXR*(ncXR-1)*log(4.0*atan(1.0))
     +          + sum(lgamAval)

      logCqRAm = 0.5*(ASigma+m)*ncXR*log(2.0)
     +           + 0.25*ncXR*(ncXR-1)*log(4.0*atan(1.0)) 
     +           + sum(lgamAmVal)

      logMLcurr = (0.5*ncXR*ASigma*log(2.0*nuVal)
     +          - 0.5*numObs*log(8.0*atan(1.0)) 
     +          - (0.5*ncXR+L+1.0)*log(4.0*atan(1.0)) 
     +          - 0.5*ncX*log(sigsqBeta))

      ans = 0.0
      do 300 ii = 1,ncC
         ans = ans + CTy(ii)*muqbetau(ii)
300   continue

      sumXi = 0.0
      do 310 ii = 1,numObs
         sumXi = sumXi + wtVec(ii)*xiVec(ii)*xiVec(ii) + phiVec(ii)
310   continue

      logMLcurr = logMLcurr + ans + sumXi

      do 320 ii = 1,ncX
         logMLcurr = logMLcurr - 0.5*(1/sigsqBeta)*(
     +               muqbetau(ii)**2 + Sigqbetau(ii,ii))
320   continue

      logMLcurr = logMLcurr + 0.5*detSig + 0.5*(sum(ncZG)+m+ncX)
      logMLcurr = logMLcurr -logCqRA + logCqRAm 
     +          - 0.5*(ASigma + m)*detBqSigmaR
      logMLcurr = logMLcurr + sum(lgamZG) 
     +          - 0.5*sum((ncZG+1)*log(Bqsigsqu))
      logMLcurr = logMLcurr - ncXR*log(AR) + ncXR*lgamXR
 
      do 330 ii = 1,ncXR
         logMLcurr = logMLcurr + nuVal*MqinvSigR(ii,ii)*muqrecaR(ii)
330   continue

      logMLcurr = logMLcurr - 0.5*(nuVal + ncXR)*sum(log(BqaR))
      logMLcurr = logMLcurr - L*log(Au) - sum(log(Bqau))
      logMLcurr = logMLcurr +  sum(muqrecau*muqrecssqu)
c      call dblepr("ans",3,logMLcurr,1)

      logMLgrid(itnum) = logMLcurr

c     Assess convergence:

      relerr = abs((logMLcurr/logMLprev) - 1.0)

      if (itnum.ge.maxIter.or.relerr.lt.VBtoler) then
         indcvgd = 1
      endif    

      if (indcvgd.eq.0) then
         goto 1
      endif

      return
      end

cccccccccccc End of twlvbfNA.f cccccccccccc
