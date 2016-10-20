cccccccccc FORTRAN subroutine twlvbfSL.f cccccccccc

c For performing variational Bayes inference for a
c Bernoulli random effects model (streamlined approach).

c The following abbreviations are used with repect
c to the parent R code:
c
c "ssq" for "sigsq", "rec" for "recip".

c Last changed: 17 FEB 2014
cccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine twlvbfSL(yAdj,XR,CG,CGTCG,CGTy,numObs,m,L,nVec,ncXR,
     +                  ncX,ncCG,nVecNew,indsStt,indsEnd,indsSplStt,
     +                  indsSplEnd,ncZG,ZRTy,xiVec,wtVec,G,H,Hwork,ipvt,
     +                  ipvtBig,det,work,workBig,sVec,Smat,SigqbetauG,
     +                  muqbetauG,SigquR,zetaVec,kappaVec,muquR,
     +                  diagTermOne,diagTermTwo,diagTermThree,xiSqd,
     +                  MqinvSigR,BqaR,muqrecaR,
     +                  BqSigmaR,muqrecssqu,Bqau,muqrecau,Bqsigsqu,
     +                  Aqsigsqu,maxIter,logMLgrid,lgamZG,
     +                  lgamXR,lgamObs)
      double precision yAdj(numObs),XR(numObs,ncXR),CG(numObs,ncCG),
     +                 CGTCG(ncCG,ncCG),CGTy(ncCG),ZRTy(m,ncXR),
     +                 xiVec(numObs),wtVec(numObs),G(ncCG,ncXR,m),
     +                 H(ncXR,ncXR,m),Hwork(ncXR,ncXR),work(ncXR,ncXR),
     +                 workBig(ncCG,ncCG),det(2),sVec(ncCG),
     +                 Smat(ncCG,ncCG),ridgeVec(ncCG),
     +                 SigqbetauG(ncCG,ncCG),muqbetauG(ncCG),sigsqBeta,
     +                 Aeps,AR,Au,nuVal,HtimestG(ncXR,ncCG,m),
     +                 SigquR(ncXR,ncXR,m),zetaVec(m,ncXR),
     +                 SigqsVec(ncCG),kappaVec(m,ncXR),SigqCGTy(ncCG),
     +                 HtimesKappa(m,ncXR),muquR(m,ncXR),
     +                 diagMidTermOne(numObs,ncCG),diagTermOne(numObs),
     +                 diagMidTermTwo(ncCG,ncXR,m),
     +                 diagMidTwoSqd(numObs,ncXR),
     +                 diagTermTwo(numObs),
     +                 diagMidTermThree(ncXR,ncXR,m),
     +                 diagMidThreeSqd(numObs,ncXR),
     +                 diagTermThree(numObs),xiSqd(numObs),
     +                 MqinvSigR(ncXR,ncXR),BqaR(ncXR),muqrecaR(ncXR),
     +                 BqSigmaR(ncXR,ncXR),muqrecssqu(L),Bqau(L),
     +                 muqrecau(L),Bqsigsqu(L),Aqsigsqu(L),detB,
     +                 detSig,logCqRA,logCqRAm,detAVec(m),
     +                 lgamAval(ncXR),lgamAmval(ncXR),logMLcurr,
     +                 VBtoler,detA,logMLprev,phiVec(numObs),
     +                 detAMat(ncCG,ncCG),
     +                 detBMat(ncCG,ncCG),detBqSigmaR,ASigma,
     +                 logMLgrid(maxIter),relerr,fitVec(numObs),
     +                 ans,sumXi,lgamZG(L),lgamXR,lgamObs
      integer i,ii,jj,k,kk,ll,L,m,itnum,numObs,ncXR,ncX,ncCG,nVec(m),
     +        indsStt(m),indsEnd(m),indsSplStt(L),indsSplEnd(L),ncZG(L),
     +        info,iStt,ipvt(ncXR),ipvtBig(ncCG),KF,indcvgd,maxIter
     +        nVecNew(m),nVecc(m)

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

      do 2 i = 1,m
         detAVec(i) = 0.0
         nVecc(i) = 0
         do 3 jj = 1,ncXR
            ZRTy(i,jj) = 0.0
            lgamAval(jj) = 0.0
            lgamAmval(jj) = 0.0
            do 4 k = indsStt(i),indsEnd(i)     
               ZRTy(i,jj) = ZRTy(i,jj) + XR(k,jj)*yAdj(k) 
4          continue
3        continue
2      continue

       do 999 i = 1,m
         if (i.eq.1) then
            nVecc(i) = 0
         endif
         if (i.gt.1) then
            nVecc(i) = nVecc(i-1) + nVec(i-1)
         endif
999    continue

c      Perform iterations:

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
               detAMat(ii,jj) = 0.0
               Smat(ii,jj) = 0.0
               SigqbetauG(ii,jj) = 0.0
7           continue
6        continue

c        Compute lambda function:

         do 8 ii = 1,numObs
            fitVec(ii) = 0.0
            wtVec(ii) = tanh(xiVec(ii)/2)/(4*xiVec(ii))
            phiVec(ii) = (xiVec(ii)/2) - log(1+exp(xiVec(ii))) 
     +                   + (0.25*xiVec(ii)*tanh(xiVec(ii)/2))
 8       continue

c        Compute G and H matrices:
             
         iStt = 1
         do 10 i = 1,m
            do 20 jj = 1,ncXR
               do 30 ii = 1,ncCG
                  G(ii,jj,i) = 0.0
                  do 40 k = indsStt(i),indsEnd(i)     
                     G(ii,jj,i) = G(ii,jj,i) 
     +                          + CG(k,ii)*wtVec(k)*XR(k,jj) 
40                continue
                  G(ii,jj,i) = 2*G(ii,jj,i)
30             continue  
               do 50 ii = 1,ncXR
                  Hwork(ii,jj) = 0.0
                  do 60 k = indsStt(i),indsEnd(i)     
                     Hwork(ii,jj) = Hwork(ii,jj) 
     +                            + XR(k,ii)*wtVec(k)*XR(k,jj) 
60                continue
                  Hwork(ii,jj) = 2*Hwork(ii,jj) + MqinvSigR(ii,jj)
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
 
c           Update sVecBeta and Smat: 
       
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

         do 175 i = 1,m
            do 176 ii = 1,ncCG
               do 180 jj = 1,ncCG
                  do 190 k = indsStt(i),indsEnd(i)   
                     SigqbetauG(ii,jj) = SigqbetauG(ii,jj) 
     +                      + 2*CG(k,ii)*wtVec(k)*CG(k,jj)   
                     detAMat(ii,jj) = detAMat(ii,jj)  
     +                      + 2*CG(k,ii)*wtVec(k)*CG(k,jj)
190               continue  
180            continue  
176         continue
175     continue

        do 191 ii = 1,ncCG
           do 192 jj = 1,ncCG
              detBMat(ii,jj) = 0.0
              SigqbetauG(ii,jj) = SigqbetauG(ii,jj) 
     +                          - Smat(ii,jj)
              if (ii.eq.jj) then
                 SigqbetauG(ii,jj) = SigqbetauG(ii,jj) 
     +                             + ridgeVec(ii)
                 detAMat(ii,jj) = detAMat(ii,jj) + ridgeVec(ii)
              endif
192       continue
191    continue

c       Invert Sigma.q.betauG:

         call dgefa(SigqbetauG,ncCG,ncCG,ipvtBig,info)
         call dgedi(SigqbetauG,ncCG,ncCG,ipvtBig,det,workBig,01)

c        Compute determinant B:
         
         iStt = 1
         do 193 i = 1,m
            do 194 k = iStt,(iStt+ncXR-1)
               do 195 kk = iStt,(iStt+ncXR-1)
                  do 196 ii = 1,ncCG
                     do 197 jj = 1,ncCG
                        detBMat(ii,jj) = detBMat(ii,jj)
     +                  + G(ii,k-2*(i-1),i)*H(k-2*(i-1),kk-2*(i-1),i)
     +                  *G(jj,kk-2*(i-1),i)             
197                 continue
196             continue
195         continue
194     continue
            iStt = iStt + ncXR
193     continue

         do 198 ii = 1,ncCG
            do 199 jj = 1,ncCG
               detBMat(ii,jj) = detAMat(ii,jj)-detBMat(ii,jj)
199        continue
198      continue

c        Compute determinant of Sigma.q.uR:

         call dgefa(detBMat,ncCG,ncCG,ipvtBig,info)
         call dgedi(detBMat,ncCG,ncCG,ipvtBig,det,workBig,10)

         detB = log(abs(det(1))) + det(2)*log(10.0)
         detSig = - (detA + detB)

c        Perform updates for mu.q.betauG:

         do 200 ii = 1,ncCG
            muqbetauG(ii) = 0.0
            do 210 jj = 1,ncCG
               muqbetauG(ii) = muqbetauG(ii) 
     +                       + SigqbetauG(ii,jj)*(CGTy(jj) - sVec(jj))
210         continue
200      continue   

c        Perform updates for Sigma.q.uR:

         do 220 i = 1,m
            do 230 ii = 1,ncXR
               zetaVec(i,ii) = 0.0
               do 240 jj = 1,ncXR
                  SigquR(ii,jj,i) = 0.0

c              Compute zector vector:

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
     +                          *SigqCGTy(jj)
370            continue 
                  muquR(i,ii) = - muquR(i,ii) + zetaVec(i,ii)  
     +                        + HtimesKappa(i,ii)
360         continue 

            do 371 k = indsStt(i),indsEnd(i)
               do 372 kk = 1,ncCG
                  fitVec(k) = fitVec(k) 
     +                        + CG(k,kk)*muqbetauG(kk)
372            continue
               do 373 jj = 1,ncXR
                  fitVec(k) = fitVec(k) + XR(k,jj)*muquR(i,jj)
373            continue 
371        continue
310      continue 
c        End of i loop 3:

c        Update xiVec:

         do 380 i = 1,numObs
            diagTermOne(i) = 0.0
            do 390 ii = 1,ncCG
               diagMidTermOne(i,ii) = 0.0
               do 400 jj = 1,ncCG
                  diagMidTermOne(i,ii) = diagMidTermOne(i,ii) 
     +                                 + CG(i,jj)*(SigqbetauG(jj,ii) 
     +                                 + muqbetauG(ii)*muqbetauG(jj))
400            continue
               diagTermOne(i) = diagTermOne(i) 
     +                             + diagMidTermOne(i,ii)*CG(i,ii)
390         continue
380      continue

         do 410 i = 1,m
            do 420 ii = 1,ncCG
               do 430 jj = 1,ncXR
                  diagMidTermTwo(ii,jj,i) = 0.0
                  do 440 kk = 1,ncCG
                     do 450 ll = 1,ncXR
                       diagMidTermTwo(ii,jj,i) = diagMidTermTwo(ii,jj,i)
     +                         - SigqbetauG(ii,kk)*G(kk,ll,i)*H(ll,jj,i)
450                  continue
440               continue 
                  diagMidTermTwo(ii,jj,i) =  diagMidTermTwo(ii,jj,i)
     +                                   + muqbetauG(ii)*muquR(i,jj)   
430            continue             
420         continue
410      continue

         do 460 i = 1,m
            do 470 k = 1,nVec(i)
               diagTermTwo(k+nVecc(i)) = 0.0
               do 480 ii = 1,ncXR
                  diagMidTwoSqd(k+nVecc(i),ii) = 0.0
                  do 490 jj = 1,ncCG
                     ll = k + nVecc(i)
                     diagMidTwoSqd(ll,ii) = diagMidTwoSqd(ll,ii)
     +                       + 2*CG(ll,jj)*diagMidTermTwo(jj,ii,i)
490               continue
                  diagTermTwo(k+nVecc(i)) = 
     +                             diagTermTwo(k+nVecc(i)) 
     +                      + diagMidTwoSqd(k+nVecc(i),ii)
     +                                  *XR(k+nVecc(i),ii)
480            continue
470         continue
460      continue

         do 491 i = 1,m
            do 492 ii = 1,ncXR
               do 493 jj = 1,ncXR
                  diagMidTermThree(ii,jj,i) = 0.0
                  diagMidTermThree(ii,jj,i) = diagMidTermThree(ii,jj,i) 
     +                      + SigquR(ii,jj,i) + muquR(i,ii)*muquR(i,jj)
493            continue
492         continue
491      continue

         do 494 i = 1,m
            do 495 k = 1,nVec(i)
               diagTermThree(k+nVecc(i)) = 0.0
               do 496 ii = 1,ncXR
                  diagMidThreeSqd(k+nVecc(i),ii) = 0.0
                  do 497 jj = 1,ncXR
                     ll = k + nVecc(i)
                     diagMidThreeSqd(ll,ii) = diagMidThreeSqd(ll,ii) 
     +               + XR(ll,jj)*diagMidTermThree(jj,ii,i)
497               continue
                  diagTermThree(k+nVecc(i)) = 
     +                           diagTermThree(k+nVecc(i)) 
     +                    + diagMidThreeSqd(k+nVecc(i),ii)
     +                                  *XR(k+nVecc(i),ii)
496            continue
495         continue
494      continue

         do 498 i = 1,numObs
            xiSqd(i) = diagTermOne(i)
498      continue

         do 499 i = 1,numObs
            xiSqd(i) = xiSqd(i) + diagTermTwo(i) 
     +                 + diagTermThree(i)
499      continue

         do 500 i = 1,numObs
            xiVec(i) = sqrt(xiSqd(i))
500      continue
c        call dblepr("ans",3,xiVec,numObs)

c        Perform updates for q*(a.R) parameters:

         do 510 ii = 1,ncXR  
            BqaR(ii) = nuVal*MqinvSigR(ii,ii) + (1/AR**2)
            muqrecaR(ii) = 0.5*(nuVal+ncXR)/BqaR(ii)

            do 520 jj = 1,ncXR
               BqSigmaR(ii,jj) = 0.0
520         continue
510      continue

         do 530 i = 1,m
            do 540 ii = 1,ncXR
               do 550 jj = 1,ncXR
                  BqSigmaR(ii,jj) = BqSigmaR(ii,jj) 
     +                  + muquR(i,ii)*muquR(i,jj) + SigquR(ii,jj,i)
550            continue          
540         continue
530      continue
 
         do 560 ii = 1,ncXR
            BqSigmaR(ii,ii) = BqSigmaR(ii,ii) 
     +                      + 2*nuVal*muqrecaR(ii)
560      continue

c        Compute determinant of B.q.SigmaR:

         call dgefa(BqSigmaR,ncXR,ncXR,ipvt,info)
         call dgedi(BqSigmaR,ncXR,ncXR,ipvt,det,work,10)
         detBqSigmaR = log(abs(det(1))) + det(2)*log(10.0)

c        Invert B.q.SigmaR:

         call dgedi(BqSigmaR,ncXR,ncXR,ipvt,det,workBig,01)

c        Perform updates for q*(SigmaR) parameters:

         do 570 ii = 1,ncXR
            do 580 jj = 1,ncXR
               MqinvSigR(ii,jj) = (nuVal + m + ncXR - 1)*BqSigmaR(ii,jj)
580          continue
570      continue

c        Perform updates for q*(sigsq.u) and q*(q.u) parameters:

         do 590 ii = 1,L
            Bqau(ii) = muqrecssqu(ii) + (1/Au**2)
            muqrecau(ii) = 1/Bqau(ii)
            Bqsigsqu(ii) = 0.0
            do 600 k = indsSplStt(ii),indsSplEnd(ii)
               Bqsigsqu(ii) = Bqsigsqu(ii) + muqbetauG(k)**2 + 
     +                      + SigqbetauG(k,k) 
600         continue
            Bqsigsqu(ii) = Bqsigsqu(ii) + 2*muqrecau(ii)
            Aqsigsqu(ii) = ncZG(ii) + 1
            muqrecssqu(ii) = Aqsigsqu(ii)/Bqsigsqu(ii)
590      continue

c        Obtain current value of log(ML):

         logMLcurr = 0.0
         ASigma = nuVal + ncXR - 1.0

c        Obtain log gamma values:

         do 610 ii = 1,ncXR
            call LGAMA(KF,0.5*(ASigma + 1.0 - ii),lgamAval(ii))
            call LGAMA(KF,0.5*(ASigma + m + 1.0 - ii),lgamAmVal(ii))
610      continue

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

         do 620 ii = 1,ncX
            logMLcurr = logMLcurr - 0.5*(1/sigsqBeta)*(
     +                  muqbetauG(ii)**2 + SigqbetauG(ii,ii))
620      continue

         ans = 0.0
         do 630 ii = 1,numObs
            ans = ans + yAdj(ii)*fitVec(ii)
630      continue

         sumXi = 0.0
         do 640 ii = 1,numObs
            sumXi = sumXi + wtVec(ii)*xiVec(ii)*xiVec(ii) + phiVec(ii)
640      continue
         logMLcurr = logMLcurr + ans + sumXi

         logMLcurr = logMLcurr + 0.5*detSig + 0.5*(sum(ncZG)+m+ncX)
         logMLcurr = logMLcurr -logCqRA + logCqRAm 
     +             - 0.5*(ASigma + m)*detBqSigmaR
         logMLcurr = logMLcurr + sum(lgamZG) 
     +             - 0.5*sum((ncZG+1)*log(Bqsigsqu))
         logMLcurr = logMLcurr - ncXR*log(AR) + ncXR*lgamXR

         do 650 ii = 1,ncXR
            logMLcurr = logMLcurr + nuVal*MqinvSigR(ii,ii)*muqrecaR(ii)
650      continue

         logMLcurr = logMLcurr - 0.5*(nuVal + ncXR)*sum(log(BqaR))
         logMLcurr = logMLcurr - L*log(Au) - sum(log(Bqau))
         logMLcurr = logMLcurr +  sum(muqrecau*muqrecssqu)
c         call dblepr("ans",3,logMLcurr,1)

         logMLgrid(itnum) = logMLcurr

c        Assess convergence:

         relerr = abs((logMLcurr/logMLprev) - 1.0)
c        call intpr("itnum",5,itnum,1)
c        call intpr("maxIter",7,maxIter,1)
c        call dblepr("relerr",6,relerr,1)
c        call dblepr("VBtoler",7,VBtoler,1)

c        call intpr("indcvgd",7,indcvgd,1)
         if (itnum.ge.maxIter.or.relerr.lt.VBtoler) then
           indcvgd = 1
         endif    
c        call intpr("indcvgd",7,indcvgd,1)

         if (indcvgd.eq.0) then
           goto 5
         endif

      return
      end

cccccccccccc End of twlvbfSL.f cccccccccccc
