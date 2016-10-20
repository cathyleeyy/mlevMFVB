cccccccccc FORTRAN subroutine twlvgfSL.f cccccccccc

c For performing variational Bayes inference for a
c Gaussian random effects model (streamlined approach).

c The following abbreviations are used with repect
c to the parent R code:
c
c "ssq" for "sigsq", "rec" for "recip".

c Last changed: 17 FEB 2014
cccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine twlvgfSL(y,XR,CG,CGTCG,CGTy,numObs,m,L,ncXR,ncX,
     +                    ncCG,indsStt,indsEnd,indsSplStt,
     +                    indsSplEnd,ncZG,ZRTy,G,H,Hwork,ipvt,ipvtBig,
     +                    det,work,workBig,sVec,Smat,
     +                    SigqbetauG,muqbetauG,SigquR,zetaVec,kappaVec,
     +                    muquR,Aqsigsqeps,Bqsigsqeps,muqrecsseps,
     +                    Bqaeps,muqrecaeps,MqinvSigR,BqaR,muqrecaR,
     +                    BqSigmaR,muqrecssqu,Bqau,muqrecau,Bqsigsqu,
     +                    Aqsigsqu,maxIter,logMLgrid,lgamZG,lgamXR,
     +                    lgamObs)
      double precision y(numObs),XR(numObs,ncXR),CG(numObs,ncCG),
     +                 CGTCG(ncCG,ncCG),CGTy(ncCG),ZRTy(m,ncXR),
     +                 G(ncCG,ncXR,m),H(ncXR,ncXR,m),Hwork(ncXR,ncXR),
     +                 work(ncXR,ncXR),workBig(ncCG,ncCG),det(2),
     +                 sVec(ncCG),Smat(ncCG,ncCG),ridgeVec(ncCG),
     +                 SigqbetauG(ncCG,ncCG),muqbetauG(ncCG),sigsqBeta,
     +                 Aeps,AR,Au,nuVal,HtimestG(ncXR,ncCG,m),
     +                 SigquR(ncXR,ncXR,m),zetaVec(m,ncXR),
     +                 SigqsVec(ncCG),kappaVec(m,ncXR),SigqCGTy(ncCG),
     +                 HtimesKappa(m,ncXR),muquR(m,ncXR),residSS,
     +                 trTermOne,trTermTwo,trTermThree,Aqsigsqeps,
     +                 Bqsigsqeps,muqrecsseps,Bqaeps,muqrecaeps,
     +                 MqinvSigR(ncXR,ncXR),BqaR(ncXR),muqrecaR(ncXR),
     +                 BqSigmaR(ncXR,ncXR),muqrecssqu(L),
     +                 Bqau(L),muqrecau(L),Bqsigsqu(L),Aqsigsqu(L),
     +                 logMLgrid(maxIter),logMLcurr,logMLprev,relerr,
     +                 VBtoler,ASigma,lgamAval(ncXR),lgamAmVal(ncXR),
     +                 logCqRA,logCqRAm,detA,detAVec(m),
     +                 detAMat(ncCG,ncCG),detBMat(ncCG,ncCG),
     +                 detB,detSig,detBqSigmaR,
     +                 lgamZG(L),lgamXR,lgamObs
      integer i,ii,jj,k,kk,ll,L,m,itnum,numObs,ncXR,ncX,ncCG,
     +        indsStt(m),indsEnd(m),indsSplStt(L),indsSplEnd(L),ncZG(L),
     +        info,iStt,ipvt(ncXR),ipvtBig(ncCG),indcvgd,KF,maxIter

c     Set hyperparameters:

      sigsqBeta = 1000.0
      Aeps = 1000.0
      AR = 1000.0
      Au = 1000.0
      detB = 0.0
      detSig = 0.0
      detBqSigmaR = 0.0

c     Set scalar starting values:

      nuVal = 2.0
      info = 0.0
      KF = 0
      logCqRA = 0.0
      logCqRAm = 0.0

c     Compute ZRTy:

      do 1 i = 1,m
         detAVec(i) = 0.0
         do 2 jj = 1,ncXR
            ZRTy(i,jj) = 0.0
            lgamAval(jj) = 0.0
            lgamAmval(jj) = 0.0
            do 3 k = indsStt(i),indsEnd(i)     
               ZRTy(i,jj) = ZRTy(i,jj) + XR(k,jj)*y(k) 
3          continue
2        continue
1     continue

c     Perform iterations:

      logMLcurr = -1.0e20
      VBtoler = 0.0000001
      indcvgd = 0
      itnum = 0
5     itnum = itnum + 1

         detA = 0.0

c        Set current log)ML) value to previous one:

         logMLprev = logMLcurr

c        Initalise vectors and matrices:

         do 6 ii = 1,ncCG
            sVec(ii) = 0.0
            do 7 jj = 1,ncCG
               Smat(ii,jj) = 0.0
7           continue
6        continue

c        Compute G and H matrices:
             
         iStt = 1
         do 10 i = 1,m
            do 20 jj = 1,ncXR
               do 30 ii = 1,ncCG
                  G(ii,jj,i) = 0.0
                  do 40 k = indsStt(i),indsEnd(i)     
                     G(ii,jj,i) = G(ii,jj,i) + CG(k,ii)*XR(k,jj) 
40                continue
                  G(ii,jj,i) = muqrecsseps*G(ii,jj,i)
30             continue  
               do 50 ii = 1,ncXR
                  Hwork(ii,jj) = 0.0
                  do 60 k = indsStt(i),indsEnd(i)     
                     Hwork(ii,jj) = Hwork(ii,jj) + XR(k,ii)*XR(k,jj) 
60                continue
                  Hwork(ii,jj) = muqrecsseps*Hwork(ii,jj) 
     +                         + MqinvSigR(ii,jj)
50             continue  
20          continue 
      
c           Compute determinant of Hwork:

            call dgefa(Hwork,ncXR,ncXR,ipvt,info)
            call dgedi(Hwork,ncXR,ncXR,ipvt,det,work,10)
 
            detAVec(i) = log(abs(det(1))) + det(2)*log(10.0)

c           Compute the determinant A:

            detA = detA + detAVec(i)

c           Invert Hwork matrices for each i:   
   
            call dgedi(Hwork,ncXR,ncXR,ipvt,det,work,01)
 
c           Store inverted Hwork into current H array:

            do 70 ii = 1,ncXR
               do 80 jj = 1,ncXR 
                  H(ii,jj,i) = Hwork(ii,jj)
80             continue 
 
c              Compute H%*%t(G):

               do 90 jj = 1,ncCG
                  HtimestG(ii,jj,i) = 0.0
                  do 100 kk = 1,ncXR
                     HtimestG(ii,jj,i) = HtimestG(ii,jj,i) 
     +                                 + H(ii,kk,i)*G(jj,kk,i)
100              continue              
90            continue         
70          continue 
       
c           Update sVec and Smat: 
       
            do 110 ii = 1,ncCG
               do 120 kk = 1,ncXR
                  sVec(ii) = sVec(ii) 
     +                     + HtimestG(kk,ii,i)*ZRTy(i,kk)
120            continue     
               do 130 jj = 1,ncCG
                  do 140  kk = 1,ncXR
                     Smat(ii,jj) = Smat(ii,jj) 
     +                            + G(ii,kk,i)*HtimestG(kk,jj,i) 
140               continue 
130            continue
110         continue  
 
10       continue    
c        End of first i loop:
   
c        Create ridge vector for Sigma.q.betauG:

         do 150 i = 1,ncX
            ridgeVec(i) = 1.0/sigsqBeta
150      continue

         iStt = ncX + 1
         do 160 ll = 1,L
            do 170 k = iStt,(iStt+ncZG(ll)-1) 
               ridgeVec(k) = muqrecssqu(ll)
170         continue        
         iStt = iStt + ncZG(ll)
160      continue  

c        Perform updates for Sigma.q.betauG parameters:

         do 180 ii = 1,ncCG
            do 190 jj = 1,ncCG
               detBMat(ii,jj) = 0.0
               SigqbetauG(ii,jj) = muqrecsseps*CGTCG(ii,jj) 
     +                           - Smat(ii,jj) 
               detAMat(ii,jj) = muqrecsseps*CGTCG(ii,jj) 
               if (ii.eq.jj) then
                  SigqbetauG(ii,jj) = SigqbetauG(ii,jj) + ridgeVec(ii)
                  detAMat(ii,jj) = detAMat(ii,jj) + ridgeVec(ii)
               endif   
190         continue  
180      continue  
c        call dblepr("detAMat",7,sum(detAMat),1)

c        Invert Sigma.q.betauG:

         call dgefa(SigqbetauG,ncCG,ncCG,ipvtBig,info)
         call dgedi(SigqbetauG,ncCG,ncCG,ipvtBig,det,workBig,01)

c        Compute determinant B:
         
         do 191 i = 1,m
            do 192 k = 1,ncXR
               do 193 kk = 1,ncXR
                  do 194 ii = 1,ncCG
                     do 195 jj = 1,ncCG
                        detBMat(ii,jj) = detBMat(ii,jj)
     +                  + G(ii,k,i)*H(k,kk,i)*G(jj,kk,i)             
195                 continue
194               continue
193            continue
192         continue
191      continue
c         call dblepr("detBMat1",8,sum(detBMat),1)

         do 196 ii = 1,ncCG
            do 197 jj = 1,ncCG
               detBMat(ii,jj) = detAMat(ii,jj)-detBMat(ii,jj)
197        continue
196     continue
c        call dblepr("detBMat2",8,sum(detBMat),1)

c        Compute determinant of Sigma.q.uR:

         call dgefa(detBMat,ncCG,ncCG,ipvtBig,info)
         call dgedi(detBMat,ncCG,ncCG,ipvtBig,det,workBig,10)

         detB = log(abs(det(1))) + det(2)*log(10.0)
c         call dblepr("detB",4,detB,1)

         detSig = - (detA + detB)
c         call dblepr("detSig",6,detSig,1)

c        Perform updates for mu.q.betauG:

         do 200 ii = 1,ncCG
            muqbetauG(ii) = 0.0
            do 210 jj = 1,ncCG
               muqbetauG(ii) = muqbetauG(ii) 
     +                       + muqrecsseps*SigqbetauG(ii,jj)*(CGTy(jj)
     +                       -sVec(jj))
210         continue
200      continue   

c        Perform updates for Sigma.q.uR:

         do 220 i = 1,m
            do 230 ii = 1,ncXR
               zetaVec(i,ii) = 0.0
               do 240 jj = 1,ncXR
                  SigquR(ii,jj,i) = 0.0

c                 Compute zector vector:

                  zetaVec(i,ii) = zetaVec(i,ii) 
     +                          + H(ii,jj,i)*ZRTy(i,jj)

                do 250 kk = 1,ncCG
                   do 260 ll = 1,ncCG
                      SigquR(ii,jj,i) = SigquR(ii,jj,i) 
     +                            + HtimestG(ii,kk,i)*SigqbetauG(kk,ll)*
     +                              HtimestG(jj,ll,i)
260                continue
250              continue
                SigquR(ii,jj,i) = SigquR(ii,jj,i) + H(ii,jj,i)
240            continue
230         continue
220      continue 
c        End of second i loop:

c        Compute SigqsVecU and SigqCGTy:

         do 290 jj = 1,ncCG
            SigqsVec(jj) = 0.0
            SigqCGTy(jj) = 0.0
            do 300 kk = 1,ncCG
               SigqsVec(jj) = SigqsVec(jj) 
     +                       + SigqbetauG(jj,kk)*sVec(kk)

               SigqCGTy(jj) = SigqCGTy(jj)
     +                      + SigqbetauG(jj,kk)*CGTy(kk)   
300         continue
290      continue
    
c        Compute kappaVec and H%*%kappaVec:
         
         residSS = 0.0
         do 310 i = 1,m
            do 320 ii = 1,ncXR
               kappaVec(i,ii) = 0.0 
               do 330 jj = 1,ncCG
                  kappaVec(i,ii) = kappaVec(i,ii) 
     +                           + G(jj,ii,i)*SigqsVec(jj)
330            continue            
320         continue

            do 340 ii = 1,ncXR
               HtimesKappa(i,ii) = 0.0
               do 350 jj = 1,ncXR
                  HtimesKappa(i,ii) = HtimesKappa(i,ii)
     +                              + H(ii,jj,i)*kappaVec(i,jj)
350            continue  
340         continue 

c           Perform updates for mu.q.uR parameters:

            do 360 ii = 1,ncXR
               muquR(i,ii) = 0.0 
               do 370 jj = 1,ncCG
                  muquR(i,ii) = muquR(i,ii) + HtimestG(ii,jj,i)
     +                        *SigqCGTy(jj)
370            continue 
                  muquR(i,ii) = - muquR(i,ii) + zetaVec(i,ii)  
     +                        + HtimesKappa(i,ii)
                  muquR(i,ii) = muqrecsseps*muquR(i,ii)
360         continue 
 
c           Compute residual sum of square:

            do 380 k = indsStt(i),indsEnd(i)
               bCurr = y(k)
               do 390 jj = 1,ncXR
                  bCurr = bCurr - XR(k,jj)*muquR(i,jj)
390            continue 
               do 400 kk = 1,ncCG
                  bCurr = bCurr - CG(k,kk)*muqbetauG(kk)
400            continue
               residSS = residSS + bCurr**2
380         continue
310      continue 
c        End of i loop 3:

c        Compute trace terms for B.q.sigsq.eps:
     
         trTermOne = 0.0
         do 410 ii = 1,ncCG
            do 420 jj = 1,ncCG
               trTermOne = trTermOne
     +                   + CGTCG(ii,jj)*SigqbetauG(jj,ii)
420         continue
410      continue

         trTermTwo = trTermOne
         trTermThree = 0.0
         do 430 i = 1,m
            do 440 ii = 1,ncXR
               do 450 jj = 1,ncXR
                  do 460 k = indsStt(i),indsEnd(i)
                     trTermTwo = trTermTwo 
     +                         + XR(k,ii)*XR(k,jj)*SigquR(ii,jj,i)
460               continue
450            continue
440         continue

            do 470 ii = 1,nCCG
               do 480 jj = 1,ncXR
                  do 490 kk = 1,ncCG
                     trTermThree = trTermThree + 
     +                  + G(ii,jj,i)*HtimestG(jj,kk,i)*SigqbetauG(kk,ii)

490               continue
480            continue
470         continue
430      continue

c        Perform updates for q*(sigsq.eps) parameters:

         Aqsigsqeps = 0.5*(numObs+1)
         Bqsigsqeps = muqrecaeps + 0.5*(residSS + trTermTwo 
     +              - 2*(1/muqrecsseps)*trTermThree)
         muqrecsseps = Aqsigsqeps/Bqsigsqeps

c        Perform updates for q*(a.eps) parameters:

         Bqaeps = muqrecsseps + (1/Aeps**2)
         muqrecaeps = 1/Bqaeps

c        Perform updates fro q*(a.R) parameters:

         do 500 ii = 1,ncXR  
            BqaR(ii) = nuVal*MqinvSigR(ii,ii) + (1/AR**2)
            muqrecaR(ii) = 0.5*(nuVal+ncXR)/BqaR(ii)

            do 510 jj = 1,ncXR
               BqSigmaR(ii,jj) = 0.0
510         continue
500      continue

         do 520 i = 1,m
            do 530 ii = 1,ncXR
               do 540 jj = 1,ncXR
                  BqSigmaR(ii,jj) = BqSigmaR(ii,jj) 
     +                  + muquR(i,ii)*muquR(i,jj) + SigquR(ii,jj,i)
540            continue          
530         continue
520      continue
 
         do 550 ii = 1,ncXR
            BqSigmaR(ii,ii) = BqSigmaR(ii,ii) 
     +                      + 2*nuVal*muqrecaR(ii)
550      continue

c        Compute determinant of B.q.SigmaR:

         call dgefa(BqSigmaR,ncXR,ncXR,ipvt,info)
         call dgedi(BqSigmaR,ncXR,ncXR,ipvt,det,work,10)
         detBqSigmaR = log(abs(det(1))) + det(2)*log(10.0)

c        Invert B.q.SigmaR:

         call dgedi(BqSigmaR,ncXR,ncXR,ipvt,det,work,01)

c        Perform updates for q*(SigmaR) parameters:

         do 560 ii = 1,ncXR
            do 570 jj = 1,ncXR
               MqinvSigR(ii,jj) = (nuVal + m + ncXR - 1)*BqSigmaR(ii,jj)
570          continue
560      continue

c        Perform updates for q*(sigsq.u) and q*(q.u) parameters:

         do 580 ii = 1,L
            Bqau(ii) = muqrecssqu(ii) + (1/Au**2)
            muqrecau(ii) = 1/Bqau(ii)
            Bqsigsqu(ii) = 0.0
            do 590 k = indsSplStt(ii),indsSplEnd(ii)
               Bqsigsqu(ii) = Bqsigsqu(ii) + muqbetauG(k)**2 
     +                      + SigqbetauG(k,k) 
590         continue
            Bqsigsqu(ii) = Bqsigsqu(ii) + 2*muqrecau(ii)
            Aqsigsqu(ii) = ncZG(ii) + 1
            muqrecssqu(ii) = Aqsigsqu(ii)/Bqsigsqu(ii)
580      continue

c        Obtain current value of log(ML):

         logMLcurr = 0.0
         ASigma = nuVal + ncXR - 1.0

c        Obtain log gamma values:

         do 600 ii = 1,ncXR
            call LGAMA(KF,0.5*(ASigma + 1.0 - ii),lgamAval(ii))
            call LGAMA(KF,0.5*(ASigma + m + 1.0 - ii),lgamAmVal(ii))
600      continue

         logCqRA = 0.5*ASigma*ncXR*log(2.0) 
     +           + 0.25*ncXR*(ncXR-1)*log(4.0*atan(1.0))
     +           + sum(lgamAval)

         logCqRAm = 0.5*(ASigma+m)*ncXR*log(2.0)
     +            + 0.25*ncXR*(ncXR-1)*log(4.0*atan(1.0)) 
     +            + sum(lgamAmVal)

         logMLcurr = (0.5*ncXR*ASigma*log(2.0*nuVal)
     +             - 0.5*numObs*log(8.0*atan(1.0)) 
     +             - (0.5*ncXR+L+1.0)*log(4.0*atan(1.0)) 
     +             - 0.5*ncX*log(sigsqBeta))

      do 610 ii = 1,ncX
         logMLcurr = logMLcurr - 0.5*(1/sigsqBeta)*(
     +               muqbetauG(ii)**2 + SigqbetauG(ii,ii))
610   continue

      logMLcurr = logMLcurr + 0.5*detSig + 0.5*(sum(ncZG)+m+ncX)
c     call dblepr("ans3",4,logMLcurr,1)

      logMLcurr = logMLcurr -logCqRA + logCqRAm 
     +          - 0.5*(ASigma + m)*detBqSigmaR
      logMLcurr = logMLcurr + sum(lgamZG) 
     +          - 0.5*sum((ncZG+1)*log(Bqsigsqu))
      logMLcurr = logMLcurr - ncXR*log(AR) + ncXR*lgamXR

      do 620 ii = 1,ncXR
         logMLcurr = logMLcurr + nuVal*MqinvSigR(ii,ii)*muqrecaR(ii)
620   continue

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
        goto 5
      endif

      return
      end

cccccccccccc End of twlvgfSL.f cccccccccccc
