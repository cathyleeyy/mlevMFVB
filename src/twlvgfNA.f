cccccccccc FORTRAN subroutine twlvgfNA.f cccccccccc

c For performing variational Bayes inference for a
c Gaussian random effects model (naive approach).

c The following abbreviations are used with repect
c to the parent R code:
c
c "ssq" for "sigsq", "rec" for "recip".

c Last changed: 17 FEb 2014
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine twlvgfNA(numObs,ncXR,ncX,L,ncZG,m,nVec,ncCG,
     +                    ncC,CTC,CTy,C,y,indsStt,indsEnd,indsSplStt,
     +                    indsSplEnd,muqrecssqu,MqinvSigR,MqinvSig,
     +                    det,ipvt,work,ipvtBig,workBig,muqrecsseps,
     +                    Sigqbetau,muqbetau,residSS,trTerm,muqrecaeps,
     +                    Aqsigsqeps,Bqsigsqeps,Bqaeps,BqaR,muqrecaR,
     +                    BqSigmaR,Bqau,muqrecau,Bqsigsqu,Aqsigsqu,
     +                    maxIter,logMLgrid,lgamAval,lgamAmVal,lgamZG,
     +                    lgamXR,lgamObs)
      double precision muqrecssqu(L),MqinvSigR(ncXR,ncXR),Aeps,nuVal,
     +                 AR,Au,CTC(ncC,ncC),CTy(ncC),C(numObs,ncC),
     +                 y(numObs),det(2),work(ncXR,ncXR),
     +                 workBig(ncC,ncC),MqinvSig(ncC,ncC),sigsqBeta,
     +                 muqrecsseps,Sigqbetau(ncC,ncC),muqbetau(ncC),
     +                 residSS,trTerm,muqrecaeps,Aqsigsqeps,Bqsigsqeps,
     +                 Bqaeps,BqaR(ncXR),muqrecaR(ncXR),Bqau(L),
     +                 BqSigmaR(ncXR,ncXR),muqrecau(L),Bqsigsqu(L),
     +                 Aqsigsqu(L),logMLcurr,logMLprev,ASigma,
     +                 logCqRA,logCqRAm,logMLgrid(maxIter),relerr,
     +                 VBtoler,detSig,detBqSigmaR,lgamAval(ncXR),
     +                 lgamAmVal(ncXR),lgamZG(L),lgamXR,lgamObs
      integer i,ii,jj,ll,kk,k,itnum,numObs,ncXR,ncX,L,ncZG(L),m,
     +        indsStt(m),indsEnd(m),indsSplStt(L),indsSplEnd(L),
     +        nVec(m),ncCG,ncC,iStt,iEnd,info,ipvt(ncXR),ipvtBig(ncC),
     +        maxIter,indcvgd,KF

c      Set hyperparameters:

       sigsqBeta = 1000.0
       Aeps = 1000.0
       AR = 1000.0
       Au = 1000.0

c      Set scalar starting values:

       nuVal = 2.0
       info = 0.0
       KF = 0
       VBtoler = 0.0000001

c      Perform iterations:

       logMLcurr = -1.0e20
       indcvgd = 0
       itnum = 0

1      itnum = itnum + 1

c      Set current log(ML) value to previous one:

       logMLprev = logMLcurr

c      Set M.q.inv.Sigma:

       do 10 ii = 1,ncX
          MqinvSig(ii,ii) = 1.0/sigsqBeta
10     continue

       iStt = ncX + 1
       do 20 ll = 1,L
          do 30 k = iStt,(iStt+ncZG(ll)-1)
             MqinvSig(k,k) = muqrecssqu(ll)
30        continue
          iStt = iStt + ncZG(ll)
20     continue

       iStt = ncCG + 1
       do 40 ii = 1,m
          iEnd = iStt + ncXR - 1
          do 50 k = iStt,iEnd
             do 60 kk = iStt,iEnd
                MqinvSig(k,kk) = MqinvSigR(k-(iStt-1),kk-(iStt-1))
60           continue
50        continue
          iStt = iStt + ncXR
40     continue

c      Perform updates for q*(beta,u) parameters:
        
       do 70 ii = 1,ncC
          do 80 jj = 1,ncC
             Sigqbetau(ii,jj) = muqrecsseps*CTC(ii,jj) 
     +                        + MqinvSig(ii,jj)
80        continue
70     continue

c      Compute determinant of Sigma.q.betau:

       call dgefa(Sigqbetau,ncC,ncC,ipvtBig,info)
       call dgedi(Sigqbetau,ncC,ncC,ipvtBig,det,workBig,10)
       detSig = -(log(abs(det(1))) + det(2)*log(10.0))
    
c      Invert Sigma.q.betau:

       call dgedi(Sigqbetau,ncC,ncC,ipvtBig,det,workBig,01)

       do 90 ii = 1,ncC
          muqbetau(ii) = 0.0
          do 100 jj = 1,ncC
             muqbetau(ii) = muqbetau(ii) 
     +                    + muqrecsseps*Sigqbetau(ii,jj)*CTy(jj)
100       continue
90     continue

c      Compute residual sum of square:

       residSS = 0.0
       do 110 i = 1,m
          do 120 k = indsStt(i),indsEnd(i)
             bCurr = y(k)
             do 130 jj = 1,ncC
                bCurr = bCurr - C(k,jj)*muqbetau(jj)
130          continue 
             residSS = residSS + bCurr**2
120       continue
110    continue

c      Compute trace term for B.q.sigsq.eps:

       trTerm = 0.0
       do 140 ii = 1,ncC
          do 150 jj = 1,ncC
             trTerm = trTerm + CTC(ii,jj)*Sigqbetau(jj,ii)
150       continue
140    continue

c     Perform updates for q*(sigsq.eps) parameters:

      Aqsigsqeps = 0.5*(numObs+1)
      Bqsigsqeps = muqrecaeps + 0.5*(residSS + trTerm)
      muqrecsseps = Aqsigsqeps/Bqsigsqeps

c     Perform updates for q*(a.eps) parameters:

      Bqaeps = muqrecsseps + (1/Aeps**2)
      muqrecaeps = 1/Bqaeps

c     Perform updates for q*(a.R) parameters:

      do 160 ii = 1,ncXR  
         BqaR(ii) = nuVal*MqinvSigR(ii,ii) + (1/AR**2)
         muqrecaR(ii) = 0.5*(nuVal+ncXR)/BqaR(ii)

         do 170 jj = 1,ncXR
            BqSigmaR(ii,jj) = 0.0
170      continue
160   continue
 
      iStt = ncCG + 1
      do 180 i = 1,m
         iEnd = iStt + ncXR - 1
         do 190 k = iStt,iEnd
            do 200 kk = iStt,iEnd
               BqSigmaR((k-(iStt-1)),(kk-(iStt-1))) =
     +         BqSigmaR((k-(iStt-1)),(kk-(iStt-1))) + Sigqbetau(k,kk)
     +         + muqbetau(k)*muqbetau(kk)
200         continue
190      continue
         iStt = iStt + ncXR
180   continue
 
      do 210 ii = 1,ncXR
         BqSigmaR(ii,ii) = BqSigmaR(ii,ii) 
     +                   + 2*nuVal*muqrecaR(ii)
210   continue

c     Compute determinant of B.q.SigmaR:

      call dgefa(BqSigmaR,ncXR,ncXR,ipvt,info)
      call dgedi(BqSigmaR,ncXR,ncXR,ipvt,det,work,10)
      detBqSigmaR = log(abs(det(1))) + det(2)*log(10.0)

c     Invert B.q.SigmaR:

      call dgedi(BqSigmaR,ncXR,ncXR,ipvt,det,work,01)

c     Perform updates for q*(SigmaR) parameters:

      do 220 ii = 1,ncXR
         do 230 jj = 1,ncXR
            MqinvSigR(ii,jj) = (nuVal + m + ncXR - 1)*BqSigmaR(ii,jj)
230      continue
220   continue

c     Perform updates for q*(sigsq.u) and q*(a.u) parameters:

      iStt = ncX + 1
      do 240 ii = 1,L
         Bqau(ii) = muqrecssqu(ii) + (1/Au**2)
         muqrecau(ii) = 1/Bqau(ii)
         Bqsigsqu(ii) = 0.0
         iEnd = iStt + ncZG(ii) - 1
         do 250 k = iStt,iEnd
            Bqsigsqu(ii) = Bqsigsqu(ii) + muqbetau(k)**2 
     +                   + Sigqbetau(k,k)
250      continue
         iStt = iStt + ncZG(ii)
         Bqsigsqu(ii) = Bqsigsqu(ii) + 2*muqrecau(ii)
         Aqsigsqu(ii) = ncZG(ii) + 1
         muqrecssqu(ii) = Aqsigsqu(ii)/Bqsigsqu(ii)
240   continue

c     Obtain current value of log(ML):

      logMLcurr = 0.0
      ASigma = nuVal + ncXR - 1.0   

c     Obtain log gamma values:

      do 260 ii = 1,ncXR
         call LGAMA(KF,0.5*(ASigma + 1.0 - ii),lgamAval(ii))
         call LGAMA(KF,0.5*(ASigma + m + 1.0 - ii),lgamAmVal(ii))
260   continue

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

      do 280 ii = 1,ncX
         logMLcurr = logMLcurr - 0.5*(1/sigsqBeta)*(
     +               muqbetau(ii)**2 + Sigqbetau(ii,ii))
280   continue

      logMLcurr = logMLcurr + 0.5*detSig + 0.5*(sum(ncZG)+m+ncX)
      logMLcurr = logMLcurr -logCqRA + logCqRAm 
     +          - 0.5*(ASigma + m)*detBqSigmaR
      logMLcurr = logMLcurr + sum(lgamZG) 
     +          - 0.5*sum((ncZG+1)*log(Bqsigsqu))
      logMLcurr = logMLcurr - ncXR*log(AR) + ncXR*lgamXR
 
      do 290 ii = 1,ncXR
         logMLcurr = logMLcurr + nuVal*MqinvSigR(ii,ii)*muqrecaR(ii)
290   continue

      logMLcurr = logMLcurr - 0.5*(nuVal + ncXR)*sum(log(BqaR))
      logMLcurr = logMLcurr - L*log(Au) - sum(log(Bqau))
      logMLcurr = logMLcurr +  sum(muqrecau*muqrecssqu)
      logMLcurr = logMLcurr - 0.5*(numObs + 1.0)*log(Bqsigsqeps) 
     +          + lgamObs
      logMLcurr = logMLcurr -log(Aeps) - log(Bqaeps) 
     +          + muqrecaeps*muqrecsseps

      logMLgrid(itnum) = logMLcurr

c     Assess convergence:

      relerr = abs((logMLcurr/logMLprev) - 1.0)
c     call intpr("itnum",5,itnum,1)
c     call intpr("maxIter",7,maxIter,1)
c     call dblepr("relerr",6,relerr,1)
c     call dblepr("VBtoler",7,VBtoler,1)

c     call intpr("indcvgd",7,indcvgd,1)
      if (itnum.ge.maxIter.or.relerr.lt.VBtoler) then
         indcvgd = 1
      endif    
c     call intpr("indcvgd",7,indcvgd,1)
  
      if (indcvgd.eq.0) then
        goto 1
      endif

      return
      end

cccccccccccc End of twlvgfNA.f cccccccccccc
