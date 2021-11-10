C     ================================================================
C     PROGRAM skyrme_rpa
C     FULLY SELF-CONSISTENT HF PLUS RPA FOR STANDARD SKYRME FORCES
C
C     Written by Gianluca Colo', Luigi Capelli, Nguyen Van Giai and
C     Ligang Cao
C     ================================================================

      implicit integer(i) 
      implicit double precision(d)

      common/C_TODO/iRPA

      open(unit=1,status='old',file='skyrme_rpa.in')
      open(unit=2,status='unknown',file='skyrme_rpa.out')

      call UNITS
      call Reader
      call SHF
      if(iRPA.ne.0)then
       call RPA
      end if
      close(2)
      close(1)

      stop
      end

C     ================================================================
C     Definition of useful constants
C     ================================================================
      subroutine UNITS()

      implicit double precision(d)

      common/C_UNITS/dpPI,dpHBARC,dpAMC2,dpHBDM,dpPMass,dpNMass

      dpPI=4.d0*datan(1.d0)
      dpHBARC=197.327053d0
      dpPMass=938.27231d0
      dpNMass=939.56563d0
      dpAMC2=(dpPMass+dpNMass)/2.d0
      dpHBDM=dpHBARC**2/2.d0/dpAMC2

      return
      end


C     ================================================================
C     Get input quantities from skyrme_rpa.in and set up main functio-
C     nal flags
C     ================================================================
      subroutine Reader()

      implicit integer(i)
      implicit double precision(d)

      character sSkyrme*5,sPurpose*80
      character sRange*80,sPotential*80

      common/C_TODO/iRPA
      common/C_ITER/dpDeltaE,iNumIter
      common/C_CTR/dpQBess,iISPIN,iORB,iPAR,iISFLIP,iRAD
      common/C_MAX/iRadialMax
      common/C_MESHR/dpStep
      common/C_OCCUP/iNOS,iNPO
      common/C_POTENTIAL/iCoul_D,iCoul_E,iISO
      common/C_CTRL/dpECFMin,dpECFMax
      common/C_EP/dpEPMax
      common/C_ISOTOPE/dpAA,dpZZ
      common/C_RC/dpEnCentrMin,dpEnCentrMax
      common/C_GAMMA/dpExIni,dpExFin,dpEStep,dpGamma

      read(1,err=1000,end=1000,FMT='(A5)')sSkyrme
      call GetForceParameters(sSkyrme)
      read(1,*,err=1000,end=1000)iRadialMax,dpStep
      read(1,*,err=1000,end=1000)dpAA,dpZZ
      write(2,140)
      write(2,100)int(dpAA),int(dpZZ)
      write(2,110)trim(sSkyrme)
      if(iRadialMax.le.0.and.dpStep.le.0.d0)then
       dpStep=0.1
       iRadialMax=int(3*1.2d0*dpAA**(1/3)/dpStep)
      end if
      read(1,*,err=1000,end=1000)dpDeltaE,iNumIter
      if(iNumIter.le.0)iNumIter=500
      iCoul_D=1
      iCoul_E=1
      iISO=1
      read(1,err=1000,end=1000,fmt='(A80)')sPotential
      if(index(sPotential,'NO COULOMB').ne.0) then
       iCoul_D=0
       iCoul_E=0
      end if
      if(index(sPotential,'NO COULOMB EXCHANGE').ne.0) then
       iCoul_D=1
       iCoul_E=0
      end if
      read(1,*,err=1000,end=1100)sPurpose
      iRPA=0
      if(index(sPurpose,'RPA').ne.0)iRPA=2
      read(1,*,err=1000,end=1000)iISPIN,iPAR
      iORB=iISPIN
      iISFLIP=0
      if(iPAR.ge.0)then
       write(2,120)iISPIN
      else
       write(2,130)iISPIN
      end if
      write(2,140)
      read(1,*,err=1000,end=1000)iRAD,dpQBess
      read(1,*,err=1000,end=1000)dpEPMax
      dpECFMin=0.d0
      dpECFMax=1000.d0
      read(1,*,err=1000,end=1000)dpExIni,dpExFin,dpEStep,dpGamma
      read(1,err=1000,end=1100,fmt='(A80)')sRange
      if(index(sRange,'CENTROID ENERGY RANGE').ne.0)then
       read(1,*,err=1000,end=1000)dpEnCentrMin,dpEnCentrMax
      end if
      if(iRAD.lt.0.or.iRAD.ge.5)then
       write(2,*)'CHANGE INPUT PARAMETER IRAD'
       stop
      end if
      if(iORB.eq.0.and.iRAD.eq.1)then
       write(2,*)'THE MONOPOLE CASE J=0 DOES NOT ADMIT IRAD=1'
       write(2,*)'CHANGE INPUT PARAMETER IRAD'
       stop
      end if
      if(iORB.eq.1.and.iRAD.eq.1)then
       write(2,*)'THE DIPOLE CASE J=1 DOES NOT ADMIT IRAD=1'
       write(2,*)'CHANGE INPUT PARAMETER IRAD'
       stop
      end if
      if(iORB.ne.1.and.iRAD.eq.3)then
       write(2,*)'IRAD=3 IS SPECIFIC FOR THE DIPOLE CASE'
       write(2,*)'CHANGE INPUT PARAMETER IRAD'
       stop
      end if
      if(iRAD.eq.4.and.dpQBess.le.1E-20)then
       write(2,*)'IRAD=4 NEEDS A NON ZERO VALUE FOR Q'
       write(2,*)'CHANGE INPUT PARAMETER DPQBESS'
       stop
      end if
      write(2,140)
 1100 continue

  100 format('HARTREE-FOCK FOR THE NUCLEUS A=',I4,1x,'AND Z=',I3)
  110 format('SKYRME PARAMETER SET: ',A7)
  120 format('RPA FOR J AND PARITY: ',I2,' +1')
  130 format('RPA FOR J AND PARITY: ',I2,' -1')
  140 format(75('*'))

      return
 1000 continue
      write(2,*)'>>> INVALID VALUES IN FILE skyrme_rpa.in <<<'
      stop '>>> PROBLEM READING FILE skyrme_rpa.in'
      end


C     ================================================================
C     Simulate BASIC MID function
C     ================================================================
      function sf_LCMid(sValue,iFrom,iCar)

      implicit character*50(s)
      implicit integer(i)

      iLen=50
      if(iFrom.le.0) iFrom=1
      if(iFrom.gt.iLen) iFrom=iLen
      if(iCar.le.0) iCar=iLen
      if(iCar.gt.iLen-iFrom+1) iCar=iLen-iFrom+1
      sf_LCMid=sValue(iFrom:iFrom+iCar-1)

      return
      end


C     ================================================================
C     Simulate BASIC LEFT function
C     ================================================================
      function sf_LCLeft(sValue,iTo)

      implicit character*50(s)
      implicit integer(i)

      iLen=50
      if(iTo.le.0) then
       iTo=0
       sTmp=''
      else
       if(iTo.gt.iLen) iTo=iLen
       sTmp=sValue(1:iTo)
      end if
      sf_LCLeft=sTmp

      return
      end


C     ================================================================
C     Get Skyrme parameters from skyrme.in file
C     ================================================================
      function bf_GetSkyrmeValue(iFile,sSkyrme)

      implicit integer(i)
      implicit double precision(d)
      implicit logical(b)
      implicit character*50(s)

      character sSkyrme*5

      common/C_SKYRME/dpT0,dpT1,dpT2,dpT3,dpX0,dpX1,dpX2,dpX3,dpW0_0
     &      ,dpW2P_0,dpAlpha,iIncludeJ2

      bFound=.false.
      bExit=.false.
      do while(bFound.eqv..false.)
       sStr=''
       sRifKey=trim(sSkyrme)
       read(iFile,*,err=100,end=100) sStr
       if(sf_LCLeft(sStr,2).ne.''.and.sf_LCLeft(sStr,2).ne.'\\')then
        if(sf_LCLeft(sStr,1).eq.'[') then
         sTmpStr=sf_LCMid(sStr,2,len(trim(sStr))-2)
         if(trim(sTmpStr).eq.trim(sRifKey))then
          read(iFile,*,err=100,end=100)dpT0
          read(iFile,*,err=100,end=100)dpT1
          read(iFile,*,err=100,end=100)dpT2
          read(iFile,*,err=100,end=100)dpT3
          read(iFile,*,err=100,end=100)dpX0
          read(iFile,*,err=100,end=100)dpX1
          read(iFile,*,err=100,end=100)dpX2
          read(iFile,*,err=100,end=100)dpX3
          read(iFile,*,err=100,end=100)dpW0_0
          read(iFile,*,err=100,end=100)dpW2P_0
          read(iFile,*,err=100,end=100)dpAlpha
          read(iFile,*,err=100,end=100)iIncludeJ2
          bFound=.true.
          exit
         end if
        end if
       end if
      end do
100   continue
      bf_GetSkyrmeValue=bFound

      return
      end


C     ================================================================
C     Search Skyrme parameterization among the known ones
C     ================================================================
      subroutine GetForceParameters(sSkyrme)

      implicit integer(i)
      implicit double precision(d)
      implicit logical(b)

      character sSkyrme*5

      open(unit=3,status='old',file='skyrme.in')

      bFound=bf_GetSkyrmeValue(3,sSkyrme)
      if(bFound.eqv..false.)then
       write(2,*)'>>> SKYRME PARAMETERS NOT FOUND: ',sSkyrme
       stop '>>> SKYRME PARAMETERS NOT FOUND'
      end if

      return
      end


C     ================================================================
C     Calculate Hartree-Fock mean field and single-particle states.
C     See Sec. 2.1 for more details
C     ================================================================
      subroutine SHF()

      implicit integer(i)
      implicit double precision(d)
      implicit logical(b)

      include 'param.rpa'

      common/C_MESHR/dpStep
      common/C_DER/dpaFCT(IP_NNP),dpaDF(IP_NNP)
      common/C_MAX/iRadialMax
      common/C_DENSITY/dpaDT(IP_NNP),dpaDN(IP_NNP),dpaDP(IP_NNP)
      common/C_SKYRME/dpT0,dpT1,dpT2,dpT3,dpX0,dpX1,dpX2,dpX3,dpW0_0
     &      ,dpW2P_0,dpAlpha,iIncludeJ2
      common/C_OCCUP/iNOS,iNPO
      common/C_POTENTIAL/iCoul_D,iCoul_E,iISO
      common/C_ITER/dpDeltaE,iNumIter
      common/C_RADII/dpaR2L(12)
      common/C_ISOTOPE/dpAA,dpZZ
      common/C_DIPO/dpaDipoKAP(7),dpaIPGL(7),dpaINBes(7)
      common/C_CTR/dpQBess,iISPIN,iORB,iPAR,iISFLIP,iRAD
      common/C_UNITS/dpPI,dpHBARC,dpAMC2,dpHBDM,dpPMass,dpNMass
      common/C_BESSL/dpaSJ(0:10),dpaDJ(0:10)

      dimension dpaU(IP_NNP),dpaGL(8,IP_NNP),dpaBes(8,IP_NNP)
      dimension dpaVSON(IP_NNP),dpaVSOP(IP_NNP),dpaVQRE(IP_NNP)
      dimension dpaRM2(IP_NNP),dpaRM3(IP_NNP)
      dimension dpaHMEN(IP_NNP),dpaHMEP(IP_NNP)
      dimension dpaDHMEN(IP_NNP),dpaDHMEP(IP_NNP)
      dimension dpaD2HMEN(IP_NNP),dpaD2HMEP(IP_NNP)
      dimension dpaAVP(IP_NNP,IP_NEOC),dpaDUNL(IP_NNP,IP_NEOC)
      dimension dpaUNL(IP_NNP,IP_NEOC),dpaDC(IP_NNP)
      dimension dpaVN(IP_NNP),dpaVP(IP_NNP),dpaVC(IP_NNP)
      dimension dpaDU(IP_NNP),dpaRM(IP_NNP),dpaVS(IP_NNP)
      dimension dpaDN1(IP_NNP),dpaDP1(IP_NNP),dpaDT1(IP_NNP)
      dimension dpaDN2(IP_NNP),dpaDP2(IP_NNP),dpaDT2(IP_NNP)
      dimension dpaDJN(IP_NNP),dpaDJP(IP_NNP),dpaDJT(IP_NNP)
      dimension dpaDJN1(IP_NNP),dpaDJP1(IP_NNP),dpaDJT1(IP_NNP)
      dimension dpaTAUN(IP_NNP),dpaTAUP(IP_NNP),dpaTAUT(IP_NNP)
      dimension dpaVCD(IP_NNP),dpaDVC(IP_NNP)
      dimension dpaMESMP(IP_NNP),dpaMESMN(IP_NNP)
      dimension dpaRMSM(IP_NNP),dpaEHF_Prev(IP_NEOC),DPAFO(IP_NEOC)
      dimension iaNN(IP_NEOC),iaLL(IP_NEOC),iaLJ(IP_NEOC),iaLT(IP_NEOC)
      dimension dpaDEG(IP_NEOC),dpaEHF(IP_NEOC),dpaEVSO(IP_NEOC)
      dimension iaNX(IP_ORB),iaLX(IP_ORB),iaJX(IP_ORB)
      dimension dpaVQ(IP_NNP),dpaUQ(IP_NNP)

      data dpUNS4PI/0.79577471d-01/,dpQP/12.566371d0/
      data iaNX/1,1,1,1,2,1,1,2,1,2,1,2,1,3,2,1,1,2,2,3,3,1,2,3
     &         ,1,2,4,3,1/
      data iaLX/0,1,1,2,0,2,3,1,3,1,4,2,4,0,2,5,5,3,3,1,1,6,4,2
     &         ,6,4,0,2,7/
      data iaJX/1,3,1,5,1,3,7,3,5,1,9,5,7,1,3,11,9,7,5,3,1,13,9
     &         ,5,11,7,1,3,15/

      open(unit=10,status='unknown',file='temp2.dat')
      open(unit=11,status='unknown',file='density.out')
      open(unit=12,status='unknown',file='temp3.dat')

      dpAZ=dpAA-dpZZ
      write(2,10)
      write(2,20)iRadialMax,dpStep
      write(2,30)
      dpW2=dpW0_0
      dpW2P=dpW2P_0
      if(dabs(dpAlpha).lt.1.e-04)dpAlpha=1.d0
      do iJ=1,iRadialMax
       dpaU(iJ)=0.d0
      end do
      do iJ=1,iRadialMax
       dpaRM(iJ)=dpStep*dfloat(iJ)
       dpaRM2(iJ)=dpaRM(iJ)*dpaRM(iJ)
       dpaRM3(iJ)=dpaRM(iJ)*dpaRM2(iJ)
      end do
      dpCRO=0.25*(dpT1*(1.d0+0.5*dpX1)+dpT2*(1.d0+0.5*dpX2))
      dpCROQ=0.125*(dpT2*(1.d0+2.d0*dpX2)-dpT1*(1.d0+2.d0*dpX1))
      iNTP=0
      do iCP=1,IP_ORB
       dpaDEG(iCP)=0.d0
       iaNN(iCP)=iaNX(iCP)
       iaLL(iCP)=iaLX(iCP)
       iaLJ(iCP)=iaJX(iCP)
       iaLT(iCP)=1
       iNTP=iNTP+iaLJ(iCP)+1
       dpaDEG(iCP)=iaLJ(iCP)+1
       if(iNTP.ge.int(dpZZ))exit
      end do
      if(iNTP.lt.int(dpZZ))then
       stop '>>> THE REQUESTED NUCLEUS HAS TOO MANY PROTON LEVELS <<<'
      end if
      iNPO=iCP
      if(iNTP.gt.int(dpZZ))dpaDEG(ICP)=dpaDEG(ICP)-(dfloat(iNTP)-dpZZ)
      iNTN=0
      do iCN=1,IP_ORB
       iCN1=iCN+iNPO
       dpaDEG(iCN1)=0.d0
       iaNN(iCN1)=iaNX(iCN)
       iaLL(iCN1)=iaLX(iCN)
       iaLJ(iCN1)=iaJX(iCN)
       iaLT(iCN1)=0
       iNTN=iNTN+iaLJ(iCN1)+1
       dpaDEG(iCN1)=iaLJ(iCN1)+1
       if(iNTN.ge.int(dpAZ))exit
      end do
      if(iNTN.lt.int(dpAZ))then
       stop '>>> THE REQUESTED NUCLEUS HAS TOO MANY NEUTRON LEVELS <<<'
      end if
      iNOS=iCN1
      if(iNTN.gt.int(dpAZ))dpaDEG(ICN1)=dpaDEG(ICN1)-(dfloat(iNTN)-dpAZ)
      dpZ_VER=0.d0
      do i=1,iNPO
       dpZ_VER=dpZ_VER+dpaDEG(i)
      end do
      dpAMZ_VER=0.d0
      do i=iNPO+1,iNOS
       dpAMZ_VER=dpAMZ_VER+dpaDEG(i)
      end do
      if(dabs(dpZ_VER-dpZZ).gt.1e-3
     &   .or.dabs((dpZ_VER+dpAMZ_VER)-dpAA).gt.1e-3)then
       write(2,*)
       write(2,*)'>>> WRONG COUNT OF PARTICLES <<<'
       write(2,*)dpZZ,dpZ_VER,dabs(dpZ_VER-dpZZ)
       write(2,*)dpAA,dpZ_VER+dpAMZ_VER,dabs((dpZ_VER+dpAMZ_VER)-dpAA)
       write(2,*)
       stop
      end if
      dpDMSHB=0.04823*dpAA/(dpAA-1)
      dpHB=1.d0/dpDMSHB
C     Starting Woods-Saxon potential
      dpR0=1.51
      dpA1=65.5
      dpA2=33.2
      dpA3=0.36
      dpA5=0.68
      dpCWPI=1.41
      dpCW2=dpCWPI*dpCWPI
      dpR=dpR0*dpAA**(1.d0/3.d0)
C     Initializations
      iIter=0
      do iJ=1,iRadialMax
       dpaHMEN(iJ)=dpHB
       dpaHMEP(iJ)=dpHB
       dpaDHMEN(iJ)=0.d0
       dpaDHMEP(iJ)=0.d0
       dpaD2HMEN(iJ)=0.d0
       dpaD2HMEP(iJ)=0.d0
       dpaMESMN(iJ)=1.d0
       dpaMESMP(iJ)=1.d0
      endDO
      dpE=0.d0
      dpE22=1.43986d0*dpStep
      dpF=0.d0
      dpRC=dpR*1.09d0/dpR0
      dpSYM=(dpAA-2*dpZZ)/dpAA
      dpZ=0.d0
      do iK=1,iRadialMax
       dpaDP(iK)=1.d0/(1.d0+dexp((dpaRM(iK)-dpRC)/0.55d0))
       dpZ=dpZ+dpaRM(iK)*dpaRM(iK)*dpaDP(iK)
      end do
      dpY=dpZZ/(dpQP*dpZ*dpStep)
      do iK=1,iRadialMax
       dpaDP(iK)=dpY*dpaDP(iK)
      end do
      do iK=1,iRadialMax
       dpY=dpQP*dpaRM(iK)*dpaRM(iK)
       dpE=dpE+dpY*dpaDP(iK)
       dpF=dpF+dpY*dpaDP(iK)/dpaRM(iK)
       dpaVC(iK)=dpE/dpaRM(iK)-dpF
      end do
      do iJ=1,iRadialMax
       dpaVC(iJ)=dpE22*(dpaVC(iJ)+dpF)
       dpPE=dexp((dpaRM(iJ)-dpR)/dpA5)
       dpFPE=-1.d0/(1.d0+dpPE)
       dpaVN(iJ)=dpFPE
       dpaVS(iJ)=-dpPE*dpFPE*dpFPE/dpA5
       do i=1,iNOS
        iL=iaLL(i)
        dpFL=(iaLJ(i)*(iaLJ(i)+2)-4*iL*(iL+1)-3)/8.d0
        iSig=2*iaLT(i)-1
        dpV=dpA1+iSig*dpA2*dpSYM
        dpaAVP(iJ,i)=dpV*dpaVN(iJ)
     &              +dpA3*dpCW2*dpaVS(iJ)*2.d0*dpFL*dfloat(iISO)
        if(i.le.iNPO)dpaAVP(iJ,i)=dpaAVP(iJ,i)+dfloat(iCoul_D)*dpaVC(iJ)
       end do
      end do
      do i=1,iNOS
       dpaEHF(i)=0.3*dpaAVP(1,i)
       if(i.le.iNPO)dpaEHF(i)=dpaEHF(i)+0.3*dfloat(iCoul_D)*dpaVC(1)
      end do
C     Iterative HF procedure
      dpXMU=0.15
      bConverge=.false.
      do while(.true.)
       dpEnHFMax=0 
       iIter=iIter+1
       if(iIter.eq.5)dpXMU=0.95
       if(iIter.eq.10)dpXMU=0.93
       if(iIter.eq.15)dpXMU=0.90
       if(iIter.eq.20)dpXMU=0.85
       if(iIter.eq.30)dpXMU=0.80
       if(iIter.eq.40)dpXMU=0.70
       if(iIter.eq.50)dpXMU=0.60
       if(iIter.eq.70)dpXMU=0.50
       if(iIter.eq.90)dpXMU=0.40
       if(iIter.eq.100)dpXMU=0.30
       if(iIter.eq.120)dpXMU=0.20
       if(iIter.eq.150)dpXMU=0.10
C     Solution of the radial Schroedinger equation
       dpEPote=0.d0
       do i=1,iNOS
        iN=iaNN(i)
        iLL1=iaLL(i)*(iaLL(i)+1)
        dpEI=dpaEHF(i)
        dpaEHF_Prev(i)=dpaEHF(i)
        iIndK=0
        do iJ=1,iRadialMax
         if(iaLT(i).eq.0)then
          dpaRMSM(iJ)=(1.d0-dpaMESMN(iJ))
          dpaVQRE(iJ)=dpaMESMN(iJ)*(dpaAVP(iJ,i)-0.25*dpaDHMEN(iJ)
     &               *dpaDHMEN(iJ)/dpaHMEN(iJ)+0.5*dpaD2HMEN(iJ))
         end if
         if(iaLT(i).eq.1)then
          dpaRMSM(iJ)=(1.d0-dpaMESMP(iJ))
          dpaVQRE(iJ)=dpaMESMP(iJ)*(dpaAVP(iJ,i)-0.25*dpaDHMEP(iJ)
     &               *dpaDHMEP(iJ)/dpaHMEP(iJ)+0.5*dpaD2HMEP(iJ))
         end if
         dpaVQRE(iJ)=dpaVQRE(iJ)/dpHB+iLL1/(dpaRM(iJ)*dpaRM(iJ))
        end do
        bCycle=.true.
        bFirstTime=.true.
        do while(bCycle)
         if(bFirstTime.eqv..false.)then
          iIndK=iIndK+2
          if(iIndK-iRadialMax.gt.0)exit
          dpEI=-0.5*dfloat(iIndK)
         else
          bFirstTime=.false.
         end if
         dpE=dpEI*dpDMSHB
         iNO=iN-1
         call NUM2B(iRadialMax-1,dpStep,dpE,dpaVQRE,dpaRMSM,dpaU
     &              ,iNO,dpHB)
         if(iNO.le.0)cycle
        end do
        dpaU(iRadialMax)=0.d0
        do iJ=1,iRadialMax
         dpaFCT(iJ)=dpaU(iJ)
        end do
        call DERIV
        dpS=0.d0
        do iJ=1,iRadialMax
         dpaDU(iJ)=dpaDF(iJ)
         dpS=dpS+dpaU(iJ)*dpaU(iJ)
        end do
        dpS=dsqrt(dpS*dpStep)
        dpSig=dpaU(1)/dabs(dpaU(1))
        dpSig=1.d0
        do iJ=1,iRadialMax
         dpaDU(iJ)=dpSig*dpaDU(iJ)/dpS
         dpaU(iJ)=dpSig*dpaU(iJ)/dpS
        end do
        dpE=dpE/dpDMSHB
        dpaEHF(i)=dpE
        dpEnHF=dabs(dpaEHF(i)-dpaEHF_Prev(i))
        if(dpEnHF.gt.dpEnHFMax)dpEnHFMax=dpEnHF
        do iJ=1,iRadialMax
         dpaDUNL(iJ,i)=dpaDU(iJ)
         dpaUNL(iJ,i)=dpaU(iJ)
        end do
        dpEPote=dpEPote+dpE*dpaDEG(i)
       end do
       dpZ_VER=0.d0
       do i=1,iNPO
        dpZ_VER=dpZ_VER+dpaDEG(i)
       end do
       dpAMZ_VER=0.d0
       do i=iNPO+1,iNOS
        dpAMZ_VER=dpAMZ_VER+dpaDEG(i)
       end do
       if(dabs(dpZ_VER-dpZZ).gt.1e-3
     &    .or.dabs((dpZ_VER+dpAMZ_VER)-dpAA).gt.1e-3)then
        write(2,*)
        write(2,*)'>>> WRONG COUNT OF PARTICLES <<<'
        write(2,*)dpZZ,dpZ_VER,dabs(dpZ_VER-dpZZ)
        write(2,*)dpAA,dpZ_VER+dpAMZ_VER,dabs((dpZ_VER+dpAMZ_VER)-dpAA)
        write(2,*)
        stop
       end if
       do i=1,iNOS
        dpaFO(i)=dpaDEG(i)/(iaLJ(i)+1)
       end do
       do iJ=1,iRadialMax
        dpE=0.d0
        dpF=0.d0
        dpE2=0.d0
        dpF2=0.d0
        dpTE=0.d0
        dpTF=0.d0
        do i=1,iNOS
         iLL1=iaLL(i)*(iaLL(i)+1)
         dpFL=(iaLJ(i)*(iaLJ(i)+2)-4*iaLL(i)*(iaLL(i)+1)-3.d0)/8.d0
         dpaEVSO(i)=dpFL
         dpY=dpaUNL(iJ,i)*dpaUNL(iJ,i)*dpaDEG(i)
         dpE=dpE+dpY*(1-iaLT(i))
         dpF=dpF+dpY*iaLT(i)
         dpEE=2*dpFL*dpY
         dpE2=dpE2+dpEE*(1-iaLT(i))
         dpF2=dpF2+dpEE*iaLT(i)
         dpY=(-dpaUNL(iJ,i)/dpaRM(iJ)+dpaDUNL(iJ,i))/dpaRM(iJ)
         dpY1=dpaUNL(iJ,i)/dpaRM2(iJ)
         dpY2=dpaDEG(i)*(dpY*dpY+iLL1*dpY1*dpY1)
         dpTE=dpTE+dpY2*(1-iaLT(i))
         dpTF=dpTF+dpY2*iaLT(i)
        end do
        dpaDP(iJ)=dpF*dpUNS4PI/dpaRM2(iJ)
        dpaDN(iJ)=dpE*dpUNS4PI/dpaRM2(iJ)
        dpaDT(iJ)=dpaDN(iJ)+dpaDP(iJ)
        dpaDJN(iJ)=dpE2*dpUNS4PI/dpaRM3(iJ)
        dpaDJP(iJ)=dpF2*dpUNS4PI/dpaRM3(iJ)
        dpaDJT(iJ)=dpaDJN(iJ)+dpaDJP(iJ)
        dpaTAUN(iJ)=dpTE*dpUNS4PI
        dpaTAUP(iJ)=dpTF*dpUNS4PI
        dpaTAUT(iJ)=dpaTAUN(iJ)+dpaTAUP(iJ)
       end do
C     Derivatives of the densities
       do iJ=1,iRadialMax
        dpaFCT(iJ)=dpaDN(iJ)
       end do
       call DERIV
       do iJ=1,iRadialMax
        dpaDN1(iJ)=dpaDF(iJ)
        dpaFCT(iJ)=dpaDP(iJ)
       end do
       call DERIV
       do iJ=1,iRadialMax
        dpaDP1(iJ)=dpaDF(iJ)
        dpaDT1(iJ)=dpaDN1(iJ)+dpaDP1(iJ)
        dpaFCT(iJ)=dpaDN1(iJ)
       end do
       call DERIV
       do iJ=1,iRadialMax
        dpaDN2(iJ)=dpaDF(iJ)
        dpaFCT(iJ)=dpaDP1(iJ)
       end do
       call DERIV
       do iJ=1,iRadialMax
        dpaDP2(iJ)=dpaDF(iJ)
        dpaDT2(iJ)=dpaDP2(iJ)+dpaDN2(iJ)
        dpaFCT(iJ)=dpaDJN(iJ)
       end do
       call DERIV
       do iJ=1,iRadialMax
        dpaDJN1(iJ)=dpaDF(iJ)
        dpaFCT(iJ)=dpaDJP(iJ)
       end do
       call DERIV
       do iJ=1,iRadialMax
        dpaDJP1(iJ)=dpaDF(iJ)
        dpaDJT1(iJ)=dpaDJP1(iJ)+dpaDJN1(iJ)
       end do
C     Calculation of the average fields
       dpE=0.d0
       dpF=0.d0
       do iJ=1,iRadialMax
        dpY=dpQP*dpaRM(iJ)*dpaRM(iJ)
        dpE=dpE+dpY*dpaDP(iJ)
        dpF=dpF+dpY*dpaDP(iJ)/dpaRM(iJ)
        dpaVCD(iJ)=dpE/dpaRM(iJ)-dpF
       end do
       dpE222=((12.d0/dpQP)**0.333333)*dpE22/dpStep
       do iJ=1,iRadialMax
        dpaVCD(iJ)=dpE22*(dpaVCD(iJ)+dpF)
        dpaVC(iJ)=dpaVCD(iJ)-dfloat(iCoul_E)
     &           *dpE222*(dpaDP(iJ)**0.333333)
       end do
       do iJ=1,iRadialMax
        dpaFCT(iJ)=dpaVC(iJ)
       end do
       call DERIV
       do iJ=1,iRadialMax
        dpaDVC(iJ)=dpaDF(iJ)
       end do
       do iJ=1,iRadialMax
        dpaHMEN(iJ)=(dpHB+dpCRO*dpaDT(iJ)+dpCROQ*dpaDN(iJ))
     &             *(1.d0-dpXMU)+dpaHMEN(iJ)*dpXMU
        dpaHMEP(iJ)=(dpHB+dpCRO*dpaDT(iJ)+dpCROQ*dpaDP(iJ))
     &             *(1.d0-dpXMU)+dpaHMEP(iJ)*dpXMU
        dpaMESMN(iJ)=dpHB/dpaHMEN(iJ)
        dpaMESMP(iJ)=dpHB/dpaHMEP(iJ)
        dpaDHMEN(iJ)=(dpCRO*dpaDT1(iJ)+dpCROQ*dpaDN1(iJ))
     &              *(1.d0-dpXMU)+dpaDHMEN(iJ)*dpXMU
        dpaDHMEP(iJ)=(dpCRO*dpaDT1(iJ)+dpCROQ*dpaDP1(iJ))
     &              *(1.d0-dpXMU)+dpaDHMEP(iJ)*dpXMU
        dpaD2HMEN(iJ)=(dpCRO*dpaDT2(iJ)+dpCROQ*dpaDN2(iJ))
     &               *(1.d0-dpXMU)+dpaD2HMEN(iJ)*dpXMU
        dpaD2HMEP(iJ)=(dpCRO*dpaDT2(iJ)+dpCROQ*dpaDP2(iJ))
     &               *(1.d0-dpXMU)+dpaD2HMEP(iJ)*dpXMU
       end do
       do iJ=1,iRadialMax
        dpVNT33=(1.d0-dpX3)*(dpAlpha+2.d0)*dpaDT(iJ)*dpaDT(iJ)
        dpVNT3=(dpVNT33+(2.d0+4.d0*dpX3)*dpaDP(iJ)*(dpaDT(iJ)
     &        +dpAlpha*dpaDN(iJ)))/24.d0
        dpVPT3=(dpVNT33+(2.d0+4.d0*dpX3)*dpaDN(iJ)*(dpaDT(iJ)
     &        +dpAlpha*dpaDP(iJ)))/24.d0
        if(dabs(dpAlpha-1.d0).ge.1.e-06)then
         if(dpaDT(iJ).eq.0) then
          dpVNT=0.d0
         else
          dpVNT=dpaDT(iJ)**(dpAlpha-1.d0)
         end if
         dpVNT3=dpVNT3*dpVNT
         dpVPT3=dpVPT3*dpVNT
        end if
        dpVNT3=-dpVNT3
        dpVPT3=-dpVPT3
        dpaVN(iJ)=dpT0*(1.d0+0.5*dpX0)*dpaDT(iJ)
     &           +0.5*(dpT2*(1.d0+0.5*dpX2)
     &           -dpT1*(1.d0+0.5*dpX1))*dpaDT1(iJ)/dpaRM(iJ)
     &           +0.125*(dpT2*(1.d0+0.5*dpX2)
     &           -3.d0*dpT1*(1.d0+0.5*dpX1))
     &           *dpaDT2(iJ)
     &           +dpCRO*dpaTAUT(iJ)
     &           -dpW2*(dpaDJT(iJ)/dpaRM(iJ)+dpaDJT1(iJ)/2.d0)
        dpaVSON(iJ)=2.d0*dpW0_0*dpaDT1(iJ)/(2.d0*dpaRM(iJ))
     &             -2.d0*dfloat(iIncludeJ2)
     &             *0.125*(dpT1*dpX1+dpT2*dpX2)*dpaDJT(iJ)/dpaRM(iJ)
        dpaVP(iJ)=dpaVN(iJ)-dpT0*(dpX0+0.5)*dpaDP(iJ)
     &           +0.25*(dpT1*(1.d0+2.d0*dpX1)
     &           +dpT2*(1.d0+2.d0*dpX2))*dpaDP1(iJ)/dpaRM(iJ)
     &           +0.0625*(3.d0*dpT1*(1.d0+2.d0*dpX1)
     &           +dpT2*(1.d0+2.d0*dpX2))*dpaDP2(iJ)
     &           +dpCROQ*dpaTAUP(iJ)
     &           -dpW2P*(dpaDJP(iJ)/dpaRM(iJ)+0.5*dpaDJP1(iJ))
     &           -0.25*dpT3*4.d0*dpVPT3
        dpaVSOP(iJ)=dpaVSON(iJ)+2.d0*dpW2P_0
     &             *dpaDP1(iJ)/(2.d0*dpaRM(iJ))
     &             +2.d0*dfloat(iIncludeJ2)
     &             *0.125*(dpT1-dpT2)*dpaDJP(iJ)/dpaRM(iJ)
        do i=1,iNPO
         dpaAVP(iJ,i)=(dpaVP(iJ)+dfloat(iISO)*dpaVSOP(iJ)*dpaEVSO(i)
     &               +dfloat(iCoul_D)*dpaVC(iJ))
     &               *(1.d0-dpXMU)+dpXMU*dpaAVP(iJ,i)
        end do
        dpaVN(iJ)=dpaVN(iJ)-dpT0*(dpX0+0.5)*dpaDN(iJ)
     &           +0.25*(dpT1*(1.d0+2.d0*dpX1)
     &           +dpT2*(1.d0+2.d0*dpX2))*dpaDN1(iJ)/dpaRM(iJ)
     &           +0.0625*(3.d0*dpT1*(1.d0+2.d0*dpX1)
     &           +dpT2*(1.d0+2.d0*dpX2))*dpaDN2(iJ)
     &           +dpCROQ*dpaTAUN(iJ)
     &           -dpW2P*(dpaDJN(iJ)/dpaRM(iJ)+0.5*dpaDJN1(iJ))
     &           -0.25*dpT3*4.d0*dpVNT3
        dpaVSON(iJ)=dpaVSON(iJ)+2.d0*dpW2P_0
     &             *dpaDN1(iJ)/(2.d0*dpaRM(iJ))+2.d0
     &             *dfloat(iIncludeJ2)*0.125*(dpT1-dpT2)
     &             *dpaDJN(iJ)/dpaRM(iJ)
        do i2=1,iNOS-iNPO
         i=i2+iNPO
         dpaAVP(iJ,i)=(dpaVN(iJ)+dfloat(iISO)*dpaVSON(iJ)*dpaEVSO(i))
     &               *(1.d0-dpXMU)+dpXMU*dpaAVP(iJ,i)
        end do
       end do
       if(dpEnHFMax.le.dpDeltaE) then
        bConverge=.true.
        exit
       end if
       if(iIter.ge.iNumIter)exit
      end do
C     End of iteration
      if(bConverge.eqv..true.)then
       write(2,40)iIter,dpEnHFMax
      else
       write(2,50)iNumIter,dpEnHFMax
      end if
      write(2,30)
      write(2,60)
      write(2,70)
      do i=1,iNPO
       write(2,80)i,iaNN(i),iaLL(i),iaLJ(i),dpaEHF(i)
      end do
      write(2,90)
      write(2,70)
      do i=iNPO+1,iNOS
       write(2,80)i,iaNN(i),iaLL(i),iaLJ(i),dpaEHF(i)
      end do
      write(2,30)
      dpFAC_EM_SO=0.022075d0
      open(unit=9,status='unknown',file='temp1.dat')
      do i=1,iNOS
       dpFAC_MU=2.39d0
       if(i.gt.iNPO)dpFAC_MU=-1.91d0
       dpFAC_L=-iaLL(i)-1
       if(iaLJ(i).gt.2*iaLL(i))dpFAC_L=iaLL(i)
       dpAInt=0.d0
       do iPRN=1,iRadialMax-2
        dpR=dfloat(iPRN)*dpStep
        dpAInt=dpAInt+dpStep*dpaUNL(iPRN,i)**2/dpR*dpaDVC(iPRN)
       end do
       dpAInt=dpFAC_EM_SO*dpFAC_MU*dpAInt*dpFAC_L
       if(dpaFO(i).ge.0.00001d0)then
        write(9,*)iaNN(i),iaLL(i),iaLJ(i),iaLT(i)
       end if
       if(dpaFO(i).ge.0.00001d0)then
        write(12,*)dpaEHF(i),iaNN(i),iaLL(i),iaLJ(i),iaLT(i)
        write(12,*)(dpaUNL(iJ,i),iJ=1,iRadialMax-2)
        write(12,*)(dpaDUNL(iJ,i),iJ=1,iRadialMax-2)
       end if
      end do
      close(9)
      call CREATE_LIST_UNOCC
      dpEMAXP=-1000.d0
      do i=1,iNPO
       if(dpaEHF(i).gt.dpEMAXP)dpEMAXP=dpaEHF(i)
      end do
      dpEMAXN=-1000.d0
      do i=iNPO+1,iNOS
       if(dpaEHF(i).gt.dpEMAXN)dpEMAXN=dpaEHF(i)
      end do
      dpECBD=0.d0
      dpECBE=0.d0
      do iJ=1,iRadialMax
       dpECBD=dpECBD+dpaVCD(iJ)*dpQP*dpaRM2(iJ)*dpaDP(iJ)*dpStep*0.5
       dpECBE=dpECBE+dpQP*dpaRM2(iJ)*(dpaDP(iJ)**(4.d0/3.d0))*dpStep
      end do
      dpCostAlpha=dpE22*((12.d0/dpQP)**(1.d0/3.d0))/dpStep
      dpECBE=dpECBE*dpCostAlpha
      dpECBE=-dpECBE*0.75
      if(iCoul_D.eq.0)dpECBD=0.d0
      if(iCoul_E.eq.0)dpECBE=0.d0
      write(2,140)dpECBD,dpECBE
      open(unit=9,status='old',file='temp1.dat')
      do i=1,iNOS
       read(9,*)i_DUM1,i_DUM2,i_DUM3,i_DUM4
      end do
      bToDo=.true.
      do while(.true.)
       if(bToDO.eqv..true.)then
        read(9,*,end=500)iL,iNJ,iNNT,iNND,iNNF
        iNN=iNNF-iNND+1
        iN=iNND-1
        iJK=0
       end if
       iN=iN+1
       iJK=iJK+1
       if(iN.le.iNNF)then
        iLL1=iL*(iL+1)
        dpESO=(iNJ*(iNJ+2)-4*iL*(iL+1)-3.d0)/8.d0
        do iKK=1,iRadialMax+1
         if(iKK-iRadialMax.le.0)then
          if(iNNT.eq.0)then
           dpaVQ(iKK)=dpaMESMN(iKK)*(dpaVN(iKK)
     &               +dfloat(iISO)*dpESO*dpaVSON(iKK)-0.25*dpaDHMEN(iKK)
     &               *dpaDHMEN(iKK)/dpaHMEN(iKK)+0.5*dpaD2HMEN(iKK))
          end if
          if(iNNT.eq.1)then
           dpaVQ(iKK)=dpaMESMP(iKK)*(dpaVP(iKK)
     &               +dfloat(iCoul_D)*dpaVC(iKK)
     &               +dfloat(iISO)*dpESO*dpaVSOP(iKK)-0.25*dpaDHMEP(iKK)
     &               *dpaDHMEP(iKK)/dpaHMEP(iKK)+0.5*dpaD2HMEP(iKK))
          end if
          dpaVQ(iKK)=dpaVQ(iKK)/dpHB+iLL1/(dpaRM(iKK)*dpaRM(iKK))
          dpaRMSM(iKK)=(1.d0-dpaMESMN(iKK))*(1.d0-iNNT)
     &                +(1.d0-dpaMESMP(iKK))*iNNT
          cycle
         end if
         dpaVQ(iKK)=iLL1/(dpStep*dpStep*dfloat(iKK*iKK))
         dpaRMSM(iKK)=0.d0
        end do
        iNO=iN-1
        dpEMIN=-500.d0
        iNPrim=0
        dpEPrim=0.d0
        if(iNPrim.eq.(iN-1))dpEMIN=dpEPrim
        iNPrim=iN
        dpEPrim=dpE
        call NUM2B(iRadialMax,dpStep,dpE,dpaVQ,dpaRMSM,dpaUQ,iNO,dpHB)
        if(iNO.lt.0)then
         write(2,150)iN,iL,iNJ,iNNT
         write(2,151)
         stop
        end if
        dpE=dpE/dpDMSHB
        do iKX=1,iRadialMax
         dpaVQ(iKX)=dpaVQ(iKX)/dpDMSHB
        end do
        iNOL=iNO+1
        write(12,*)dpE,iNOL,iL,iNJ,iNNT
        write(12,160)(dpaUQ(IK),IK=1,iRadialMax-2)
        do iPRN=1,iRadialMax-2
         dpR=dfloat(iPRN)*dpStep
        end do
        do iK=1,iRadialMax
         dpaFCT(iK)=dpaUQ(iK)
        end do
        call DERIV
        write(12,160)(dpaDF(iK),iK=1,iRadialMax-2)
        bToDo=.false.
        cycle
       end if
       bToDo=.true.
      end do
  500 continue
      write(10,*)iRadialMax-2
      write(10,170)dpT0,dpT1,dpT2,dpT3,dpX0
      write(10,170)dpAlpha,dpX3,dpX1,dpX2,dpW2,dpW2P
      write(10,170)dpAA,dpZZ,dpStep
      iNXXX=1
      iND=1
      iNX1=2
      iNX2=2
      iNX3=2
      iNNN=iRadialMax-2
      write(10,180)iNXXX,iND,iNX1,iNX2,iNX3
      write(10,190)(dpaDT(iJ),iJ=1,iNNN),(dpaTAUT(iJ),iJ=1,iNNN)
      write(10,190)(dpaDP(iJ),iJ=1,iNNN),(dpaDN(iJ),iJ=1,iNNN)
      write(10,190)(dpaHMEN(iJ),iJ=1,iNNN),(dpaHMEP(iJ),iJ=1,iNNN)
      write(10,190)(dpaVN(iJ),iJ=1,iNNN),(dpaVP(iJ),iJ=1,iNNN)
      write(10,190)(dpaDHMEN(iJ),iJ=1,iNNN),(dpaDHMEP(iJ),iJ=1,iNNN)
      write(10,190)(dpaVSON(iJ),iJ=1,iNNN),(dpaVSOP(iJ),iJ=1,iNNN)
      write(10,190)(dpaVC(iJ),iJ=1,iNNN)
      write(10,190)(dpaDT1(iJ),iJ=1,iNNN),(dpaDT2(iJ),iJ=1,iNNN)
      write(11,*)'    Step        P density     N density',
     &           '     P+N density  Charge density'
      do iJ=1,iRadialMax
C     Calculation of the charge density
       dpaDC(iJ)=dpf_DSC(dpaRM(iJ))
       write(11,200)dpStep*iJ,dpaDP(iJ),dpaDN(iJ),dpaDP(iJ)+dpaDN(iJ)
     &              ,dpaDC(iJ)
      end do
      write(11,*)
      write(11,*)
      write(11,*)'   SPURIOUS STATES'
      write(11,*)'   ---------------'
      write(11,*)
      write(11,*)'     Step       P density     N density',
     &           '     P+N density'
      do iJ=1,iRadialMax
       write(11,210)dpStep*iJ,dpaDP1(iJ),dpaDN1(iJ)
     &             ,dpaDP1(iJ)+dpaDN1(iJ)
      end do
      do i=1,iNOS
       write(12,*)dpaFO(i)
      end do
      dpEspo=0.d0
      do iJ=1,iRadialMax
       dpEspo=dpEspo+dpW2*dpStep*dpaDT(iJ)*(2.d0*dpaRM(iJ)*dpaDJT(iJ)
     &       +dpaRM2(iJ)*dpaDJT1(iJ))
     &       +dpW2P*dpStep*dpaDP(iJ)*(2.d0*dpaRM(iJ)*dpaDJP(iJ)
     &       +dpaRM2(iJ)*dpaDJP1(iJ))
     &       +dpW2P*dpStep*dpaDN(iJ)*(2.d0*dpaRM(iJ)*dpaDJN(iJ)
     &       +dpaRM2(iJ)*dpaDJN1(iJ))
      end do
      dpEspo=dpEspo*dpQP*(-0.5d0)
      write(2,100)dpEspo
C     Contributions to the total energy
      dpEREA=0.d0
      dpE_DI_T0=0.d0
      dpE_DI_T12=0.d0
      dpE_DI_T3=0.d0
      dpE_D_COUL=0.d0
      do iJ=1,iRadialMax
       dpEREA=dpEREA-dpaRM2(iJ)*(0.5*(1.d0-dpX3)*dpaDT(iJ)*dpaDT(iJ)
     &       +(2.d0*dpX3+1.d0)
     &       *dpaDN(iJ)*dpaDP(iJ))*(dpaDT(iJ)**dpAlpha)
     &       *dpT3*dpAlpha/24.d0
       dpE_DI_T0=dpE_DI_T0+0.25d0*dpT0*((2.d0+dpX0)*dpaDT(iJ)**2
     &          -(2.d0*dpX0+1.d0)*(dpaDP(iJ)**2
     &          +dpaDN(iJ)**2))*dpaRM2(iJ)
       dpE_DI_T3=dpE_DI_T3+(1.d0/24.d0)*dpT3*dpaDT(iJ)**dpAlpha
     &          *((2.d0+dpX3)*dpaDT(iJ)**2-(2.d0*dpX3+1.d0)
     &          *(dpaDP(iJ)**2+dpaDN(iJ)**2))*dpaRM2(iJ)
       dpCTR=(1.d0/8.d0)*(dpT1*(2.d0+dpX1)
     &      +dpT2*(2.d0+dpX2))*dpaTAUT(iJ)*dpaDT(iJ)
     &      +(1.d0/8.d0)*(dpT2*(2.d0*dpX2+1.d0)-dpT1*(2.d0*dpX1+1.d0))
     &      *(dpaTAUP(iJ)*dpaDP(iJ)+dpaTAUN(iJ)*dpaDN(iJ))
     &      +(1.d0/32.d0)*(3.d0*dpT1*(2.d0+dpX1)
     &      -dpT2*(2.d0+dpX2))*dpaDT1(iJ)**2
     &      -(1.d0/32.d0)*(3.d0*dpT1*(2.d0*dpX1+1.d0)
     &      +dpT2*(2.d0*dpX2+1.d0))
     &      *(dpaDP1(iJ)**2+dpaDN1(iJ)**2)
       dpE_DI_T12=dpE_DI_T12+dpCTR*dpaRM2(iJ)
       dpDINT_CDE=0.d0
       do iJJ=1,iRadialMax
        if(iJJ.le.iJ)dpFAC=dpaRM2(iJJ)/dpaRM(iJ)
        if(iJJ.gt.iJ)dpFAC=dpaRM(iJJ)
        dpDINT_CDE=dpDINT_CDE+(dpaDN(iJJ)-dpaDP(iJJ))*dpFAC*dpStep
       end do
       dpE_D_COUL=dpE_D_COUL+dpDINT_CDE*dpaDP(iJ)*dpaRM2(iJ)*dpStep
      end do
      dpE_DI_T0=dpE_DI_T0*dpStep*dpQP
      dpE_DI_T3=dpE_DI_T3*dpStep*dpQP
      dpE_DI_T12=dpE_DI_T12*dpStep*dpQP
      dpE_D_COUL=dpE_D_COUL*dpE22/dpStep*dpQP**2/(dpAA-2*dpZZ)
      dpERea=dpERea*dpStep*dpQP
      dpERea=dpERea+dpECBE/3.d0
      dpERea=dpERea/dpAA
      dpEcTot=0.d0
      do iJ=1,iRadialMax
       dpEcTot=dpEcTot+dpaTAUT(iJ)*dpaRM2(iJ)
      end do
      dpEcTot=dpEcTot*dpQP*dpStep*dpHB
      dpEPote=dpEPote/dpAA
      dpEcTot=dpEcTot/dpAA
      dpH0=0.5*(dpEPote+dpEcTot)
      dpETot=dpH0+dpERea
      write(2,220)dpEcTot,dpH0,dpERea,dpETot
C     Calculation of radii
      dpRN=0.d0
      dpRT=0.d0
      dpRP=0.d0
      dpRC=0.d0
      do iJ=1,iRadialMax
       dpX4=dpaRM2(iJ)*dpaRM2(iJ)
       dpRN=dpRN+dpaDN(iJ)*dpX4
       dpRT=dpRT+dpaDT(iJ)*dpX4
       dpRP=dpRP+dpaDP(iJ)*dpX4
       dpRC=dpRC+dpaDC(iJ)*dpX4
      end do
      dpRN=dsqrt(dpQP*(dpStep*dpRN/(dpAZ)))
      dpRT=dsqrt(dpQP*(dpStep*dpRT/dpAA))
      if(int(dpZZ).ne.0)dpRP=dsqrt(dpQP*(dpStep*dpRP/dpZZ))
      if(int(dpZZ).ne.0)dpRC=dsqrt(dpQP*(dpStep*dpRC/dpZZ))
      write(2,230)dpRN,dpRP,dpRC
      iLJ1=1
      write(2,240)
      do while(iLJ1.le.7)
       dpRP=0.d0
       do iJ=1,iRadialMax
        dpX4=dpaRM2(iJ)**iLJ1
        dpX4PN=dpaRM(iJ)**(2+iLJ1) 
        dpRP=dpRP+dpaDT(iJ)*dpX4
       end do
       dpRP=dpQP*dpStep*dpRP/dpAA
       write(2,250)iLJ1,dpRP
       dpaR2L(iLJ1)=dpRP
       do iJ=1,iRadialMax
        if(iLJ1.eq.1)dpaGL(1,iJ)=4.d0*dpaRM2(iJ)
        if(iLJ1.eq.2)dpaGL(2,iJ)=3.d0
        if(iLJ1.eq.3)dpaGL(3,iJ)=10.d0*dpaRM2(iJ)
        if(iLJ1.eq.4)dpaGL(4,iJ)=21.d0*dpaRM2(iJ)*dpaRM2(iJ)
        if(iLJ1.eq.5)dpaGL(5,iJ)=36.d0*dpaRM2(iJ)**3
        if(iLJ1.eq.6)dpaGL(6,iJ)=55.d0*dpaRM2(iJ)**4
        if(iLJ1.eq.7)dpaGL(7,iJ)=78.d0*dpaRM2(iJ)**5
        dpRDist=dfloat(iJ)*dpStep
        call SPHJ(dpQBess*dpRDist)
        dpaBes(ILJ1,IJ)=(dpaDJ(ILJ1-1)*dpQBess)**2+dfloat(iLJ1-1)
     &                 *(dfloat(iLJ1-1)+1.d0)*(dpaSJ(iLJ1-1)/dpRDist)**2
       end do
       dpaDipoKAP(iLJ1)=0.d0
       dpaIPGL(iLJ1)=0.d0
       dpaINBes(iLJ1)=0.d0
       dpSkyTerm=(dpT1*(1.d0+(dpX1/2.d0)))+(dpT2*(1.d0+(dpX2/2.d0)))
       dpDenom=(1.d0/0.04823)
       dpAAZZ=4.d0*dpZZ*(dpAA-dpZZ)/dpAA/dpAA
       dpNN=(dpAA-dpZZ)
       dpGLPU=0.d0
       dpGLPD=0.d0
       dpGLPS=0.d0
       dpGLPNP=0.d0
       do iJ=1,iRadialMax
        dpGLPU=dpGLPU+dpQP*dpaGL(iLJ1,iJ)*dpaRM2(iJ)*dpaDP(iJ)*dpaDN(iJ)
        dpGLPD=dpGLPD+dpQP*dpaGL(iLJ1,iJ)*dpaRM2(iJ)*(dpaDP(iJ)
     &        +dpaDN(iJ))
        dpGLPS=dpGLPS+dpQP*dpaBes(iLJ1,iJ)*dpaRM2(iJ)*(dpaDP(iJ)
     &        +dpaDN(iJ))
       end do
       dpaDipoKAP(iLJ1)=dpSkyTerm/dpDenom*(dpGLPU/dpGLPD)/dpAAZZ
       dpaIPGL(iLJ1)=dpGLPD*dpStep
       dpaINBes(iLJ1)=dpGLPS*dpStep
       iLJ1=iLJ1+1
      end do

      close(12)
      close(11)
      close(10)

   10 format('FULL HF POTENTIAL INCLUDED')
   20 format('NUMBER OF POINTS AND MESH: ',I4,2X,F6.4)
   30 format(75('*'))
   40 format('CONVERGENCE REACHED AFTER ',I4,' ITERATIONS'/
     &       'AT THE DESIRED ACCURACY ',E12.5)
   50 format('CONVERGENCE NOT REACHED AFTER ',I4,' ITERATIONS'/
     &       'WITH MAX. DIFF. BETWEEN HF STATES ',E12.5)
   60 format('PROTON STATES')
   70 format('  i    n    l     j       Energy')
   80 format(I3,2X,I3,2X,I3,2X,I3,'/2',3x,E12.5)
   90 format('NEUTRON STATES')
  100 format(/'      E(S.O.)=',E12.5)
  120 format(  //' ITERATION ', I3,'      XMU=',F5.3//)
  140 format('  ECB(DIRECT)=',E12.5,5X,/'ECB(EXCHANGE)=',E12.5)
  150 format(/2X,'SEARCH OF THE STATE N=',I4,' L=',I4,' J=',I4,
     &       /2X'CHARGE='I4,3X,'HAS NOT CONVERGED')
  151 format(2X,'CHANGE VALUES OF THE BOX RADIUS OR CUTOFF ENERGY')
  160 format(6(E12.6,1X))
  170 format(6E12.6)
  180 format(15I3)
  190 format(6E12.6)
  200 format(5(2X,E12.6))
  210 format(4(2X,E12.6))
  220 format(/'  EC/A=',E15.8,/' EHF/A=',E15.8,/'ERea/A=',E15.8,
     &       /'ETot/A=',E15.8)
  230 format(/'RN=',E13.5,/'RP=',E13.5,/'RC=',E13.5//)
  240 format('GROUND STATES EXPECTATION VALUES',/'   L   <r^(2L-2)>')
  250 format(2X,I2,2X,E12.5,2X,E12.5)

      return
      end


C     ================================================================
C     Determine quantum numbers of unoccupied states
C     ================================================================
      subroutine CREATE_LIST_UNOCC()

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'
      parameter(IP_LMX=15)

      common/C_HF/iNMax,iNOcc,iNUnocc,iNOrb
      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP)
     &      ,iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF)
     &      ,iaIPP(IP_NCF),iaINN(IP_NCF)
      common/C_CTR/dpQBess,iISPIN,iORB,iPAR,iISFLIP,iRAD
      common/C_OCCUP/iNOS,iNPO

      dimension iaChosen(0:IP_LMX,2,0:1),iaLastN(0:IP_LMX,2,0:1)

      open(unit=9,status='old',file='temp1.dat')

      iPR=0
      do i=1,1000
       read(9,*,end=1000)iaNN(i),iaLL(i),iaLJ(i),iaIQ(i)
       if(iaIQ(i).gt.0)iPR=i
       if((iaLJ(i)+2*iISPIN).gt.(2*IP_LMX+1))then
        write(2,*)'>>> INCREASE PARAMETER IP_LMX IN THE',
     &            ' SUBROUTINE CREATE_LIST_UNOCC <<<'
        stop
       end if
      end do
 1000 continue
      iNOcc=i-1
      iNOS=iNOcc
      iNPO=iPR
      do iL_UN=0,IP_LMX
       do iS=1,2
        do iQQ=0,1
         iaChosen(iL_UN,iS,iQQ)=0
         iaLastN(iL_UN,iS,iQQ)=0
        end do
       end do
      end do
      do i=1,iNOcc
       iJPMin=iabs(iaLJ(i)-2*iISPIN)
       iJPMax=iaLJ(i)+2*iISPIN
       iPAR_Occ=1-2*mod(iaLL(i),2)
       if(iaLJ(i).lt.(2*iaLL(i)))iSS1=1
       if(iaLJ(i).gt.(2*iaLL(i)))iSS1=2
       if(iaNN(i).gt.iaLastN(iaLL(i),iSS1,iaIQ(i)))
     &    iaLastN(iaLL(i),iSS1,iaIQ(i))=iaNN(i)
       do iL_UN=0,IP_LMX
        iPAR_Unocc=1-2*mod(iL_UN,2)
        iTest=iPAR_Occ*iPAR_Unocc*iPAR
        if(iTest.le.0)cycle
        do iS=1,2
         if(iL_UN.eq.0.and.iS.eq.1)cycle
         iSS=-1+2*(iS-1)
         iJ_UN=2*iL_UN+iSS
         if(iJ_UN.lt.iJPMin)cycle
         if(iJ_UN.gt.iJPMax)cycle
         iaChosen(iL_UN,iS,iaIQ(i))=1
        end do
       end do
      end do
      do iQQ=1,0,-1
       do iL_UN=0,IP_LMX
        do iS=1,2
         if(iaChosen(iL_UN,iS,iQQ).eq.1)then
          iJJJ=2*iL_UN-1+2*(iS-1)
          write(9,*)iL_UN,iJJJ,iQQ,iaLastN(iL_UN,iS,iQQ)+1,
     &              iaLastN(iL_UN,iS,iQQ)+IP_MOREN
         end if
        end do
       end do
      end do
      close(9)

      return
      end


C     ================================================================
C     Gaussian proton form factor EXP(-R^2/MU^2)
C     for the calculation of the charge density,
C     where MU=0.65 FM
C     ================================================================
      function dpf_DSC(dpZ1)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      dimension dpaYK(10),dpaPI(10),dpaSIG(2)

      data(dpaYK(i),i=1,10)/0.24534071d0,0.73747373d0
     &    ,1.2340762d0,1.7385377d0,2.2549740d0,2.7888061d0
     &    ,3.3478546d0,3.9447640d0,4.6036824d0,5.3874809d0/
      data(dpaPI(i),i=1,10)/0.49092150d0,0.49384339d0
     &    ,0.49992087d0,0.50967903d0,0.52408035d0,0.54485174d0
     &    ,0.57526244d0,0.62227870d0,0.70433296d0,0.89859196d0/
      data dpaSIG/1.d0,-1.d0/

      dpXMU=0.65
      dpZW=1.d0/(0.4431125*dpXMU**3)
      dpSG=0.d0
      dpXSM=dpZ1/dpXMU
      do i=1,10
       do iK=1,2
        dpY=dpZ1+dpXMU*dpaYK(i)*dpaSIG(iK)
        if(dpY.gt.0.d0)then
         dpF=dpZW*dpf_GUGUS(dpY)*dpY*dpY
         dpYSM=dpY/dpXMU
         dpSG=dpSG+dpXMU*dpaPI(i)*dpF*dpf_VKG(0,dpXSM,dpYSM)
        else
         exit
        end if
       end do
      end do
      dpf_DSC=dpSG

      return
      end


C     ================================================================
      function dpf_GUGUS(dpX)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_MESHR/dpStep
      common/C_MAX/iRadialMax
      common/C_DENSITY/dpaDT(IP_NNP),dpaDN(IP_NNP),dpaDP(IP_NNP)

      dpX2=2.d0*dpStep
      if(dpX.le.dpX2)then
       dpGUGUS=((4.d0*dpaDP(1)-dpaDP(2))
     &        +(dpaDP(2)-dpaDP(1))*dpX*dpX)/3.d0
      else
       i=dpX/dpStep
       if(i.lt.iRadialMax)then
        dpY=(dpX-i*dpStep)/dpStep
        dpS=dpaDP(i)+0.5*dpY*(dpaDP(i+1)-dpaDP(i-1))
        dpGUGUS=dpS+0.5*dpY*dpY*(dpaDP(i+1)+dpaDP(i-1)-2.d0*dpaDP(i))
       else
        dpGUGUS=0.d0
       end if
      end if
      dpf_GUGUS=dpGUGUS

      return
      end


C     ================================================================
      function dpf_VKG(iL,dpX,dpY)

      implicit integer(i)
      implicit double precision(d)
      implicit logical(b)

      dimension dpaFJ(50)

      dpDelta=dpY-dpX
      dpA=2.d0*dpX*dpY
      if(dpA-15.d0.le.0.d0)then
       bJump=.false.
       if(dpA-0.01.le.0.d0)then
        dpArg=dexp(-(dpX*dpX+dpY*dpY))
        dpaFJ0=dpArg*(1.d0+dpA*dpA/6.d0)
        if(iL-1.le.0)then
         dpaFJ(1)=dpaFJ0
         dpaFJ(2)=dpArg*dpA*(10.d0+dpA*dpA)/30.d0
         dpf_VKG=dfloat(iL+iL+1)*dpaFJ(iL+1)
         return
        else
         bJump=.true.
        end if
       end if
       if(bJump.eqv..false.)then
        dpU=dpDelta*dpDelta
        dpV=(dpX+dpY)*(dpX+dpY)
        dpU=dexp(-dpU)
        dpV=dexp(-dpV)
        dpaFJ0=(dpU-dpV)/(2.d0*dpA)
       end if
       iL2=iL+5
       iL2=iL+10
       dpaFJ(iL2)=1.e-10*dfloat(2*iL2+1)/dpA
       dpaFJ(iL2+1)=1.e-10
       iL3=iL2-1
       do iLL=1,iL3
        iL1=iL2-iLL
        dpaFJ(iL1)=dfloat(2*iL1+1)*dpaFJ(iL1+1)/dpA+dpaFJ(iL1+2)
        if(dpaFJ(iL1)-1.e+30.le.0.d0)cycle
        do iL4=iL1,iL2
         dpaFJ(iL4)=1.e-10*dpaFJ(iL4)
        end do
       end do
       dpZZ=dpaFJ0/dpaFJ(1)
       iL2=iL2-9
       do iL1=1,iL2
        dpaFJ(iL1)=dpZZ*dpaFJ(iL1)
       end do
       dpf_VKG=dfloat(iL+iL+1)*dpaFJ(iL+1)
       return
      end if
      dpU=dpDelta*dpDelta
      dpV=(dpX+dpY)*(dpX+dpY)
      dpU=dexp(-dpU)
      dpV=dexp(-dpV)
      dpaFJ(1)=(dpU-dpV)/(2.d0*dpA)
      dpaFJ(2)=((dpA-1.d0)*dpU+(dpA+1.d0)*dpV)/(2.d0*dpA*dpA)
      if(iL-1.gt.0)then
       do iL1=2,iL
        dpaFJ(iL1+1)=-dfloat(2*iL1-1)*dpaFJ(iL1)/dpA+dpaFJ(iL1-1)
       end do
      end if
      dpf_VKG=dfloat(iL+iL+1)*dpaFJ(iL+1)

      return
      end


C     ================================================================
C     Five point differentiation formula
C     ================================================================
      subroutine DERIV()

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_MESHR/dpStep
      common/C_DER/dpaFCT(IP_NNP),dpaDF(IP_NNP)
      common/C_MAX/iRadialMax

      dimension dpaA(5,5)
      data dpaA(1,1),dpaA(1,2),dpaA(1,3),dpaA(1,4),dpaA(1,5)
     &     /-50.d0,96.d0,-72.d0,32.d0,-6.d0/
      data dpaA(2,1),dpaA(2,2),dpaA(2,3),dpaA(2,4),dpaA(2,5)
     &     /-6.d0,-20.d0,36.d0,-12.d0,2.d0/
      data dpaA(3,1),dpaA(3,2),dpaA(3,3),dpaA(3,4),dpaA(3,5)
     &     /2.d0,-16.d0,0.d0,16.d0,-2.d0/
      data dpaA(4,1),dpaA(4,2),dpaA(4,3),dpaA(4,4),dpaA(4,5)
     &     /-2.d0,12.d0,-36.d0,20.d0,6.d0/
      data dpaA(5,1),dpaA(5,2),dpaA(5,3),dpaA(5,4),dpaA(5,5)
     &     /6.d0,-32.d0,72.d0,-96.d0,50.d0/
      data dpEMFact/24.d0/

      do iJ=1,iRadialMax
       iK=3
       if(iJ.lt.3)iK=iJ
       if(iJ.gt.iRadialMax-2)iK=iJ-iRadialMax+5
       dpSUM=0.d0
       do i=1,5
        iJJ=iJ+i-iK
        dpSUM=dpSUM+dpaA(iK,i)*dpaFCT(iJJ)
       end do
       dpaDF(iJ)=dpSUM/(dpStep*dpEMFact)
      end do

      return
      end


C     ================================================================
C     Numerov algorithm for the solution of the Schroedinger equation 
C     ================================================================
      subroutine NUM2B(iN,dpStep,dpE,dpaU,dpaRMS,dpaS,iNO,dpHB)

      implicit integer(i)
      implicit double precision(d)
      implicit logical(b)

      include 'param.rpa'

      dimension dpaU(IP_NNP),dpaRMS(IP_NNP),dpaS(IP_NNP)
      dimension dpaP(IP_NNP),dpaT(IP_NNP)

      iN1=iN+1
      dpH12=dpStep*dpStep/12.d0
      dpE0=-2.8938d0
      dpE1=0.d0
      dpDE=0.24115d0
      dpE0=dpE0+dpDE
      do i=1,iN1
       dpaP(i)=dpaU(i)+dpaRMS(i)*dpE0
      end do    
      dpUMin=dpaP(1)
      do i=2,iN1
       if(dpaP(i).lt.dpUMin)dpUMin=dpaP(i)
      end do
      do iE=1,2
       iKK=0
       if(iE.eq.1)then
        dpE0=dpUMin-dpDE 
       else if(iE.eq.2)then
        dpDE=0.24115d0
        dpE0=dpE1-dpDE
       end if
       do while(.true.)
        dpE0=dpE0+dpDE
        iKK=iKK+1
        if(iKK.gt.500)then
         iNO=-1
         return
        end if
        do i=1,iN1
         dpaP(i)=dpaU(i)+dpaRMS(i)*dpE0
        end do
        dpaS(1)=1.e-10
        dpB0=0.d0
        dpAA=dpH12*(dpaP(1)-dpE0)
        dpB1=dpaS(1)*(1.d0-dpAA)
        do i=2,iN1
         dpB2=12.d0*dpaS(i-1)-10.d0*dpB1-dpB0
         dpAA=dpH12*(dpaP(i)-dpE0)
         dpaS(i)=dpB2/(1.d0-dpAA)
         dpB0=dpB1
         dpB1=dpB2
        end do
        iNEL=0
        do i=8,iN
         if(dpaS(i-1)*dpaS(i).lt.0.d0)then
          iNEL=iNEL+2
          cycle
         else if(dpaS(i-1)*dpaS(i).eq.0)then
          iNEL=iNEL+1
         end if
        end do
        iNEL=iNEL/2
        if(iE.eq.1)then
         if(iNEL-iNO.eq.0)then
          dpE1=dpE0
          dpSE1=dpaS(iN)
          iNEL1=iNEL
          exit
         else if(iNEL-iNO.gt.0)then
          dpE0=dpE0-dpDE
          dpDE=0.5*dpDE
         end if
        else if(iE.eq.2)then
         if(iNEL-iNO-1.eq.0)then
          dpE2=dpE0
          dpSE2=dpaS(iN)
          iNEL2=iNEL
          exit
         else if(iNEL-iNO-1.gt.0)then
          dpE0=dpE0-dpDE
          dpDE=0.5*dpDE
         end if
        else
         exit
        end if
       end do
      end do
      iKK=0
      do i=1,iN1
       dpaT(i)=dpaS(i)
      end do
      do while(.true.)
       bDoThis=.true.
       if(iNEL2-iNEL1-1.ne.0)then
        dpET1=dpE1*dpHB
        dpET2=dpE2*dpHB
        write(2,100)dpET1,dpSE1,iNEL1,dpET2,dpSE2,iNEL2
        iNO=-1
        return
       end if
       if(dpSE1*dpSE2.gt.0.d0)then
        dpET1=dpE1*dpHB
        dpET2=dpE2*dpHB
        write(2,100)dpET1,dpSE1,iNEL1,dpET2,dpSE2,iNEL2
        iNO=-1
        return
       else if(dpSE1*dpSE2.eq.0.d0)then
        if(dpSE1.eq.0.d0)then
         dpET1=dpE1*dpHB
         dpET2=dpE2*dpHB
         write(2,100)dpET1,dpSE1,iNEL1,dpET2,dpSE2,iNEL2
         iNO=-1
         return
        end if
        dpE=dpE2
        iNEL=iNEL2
       else
        dpDE=0.5d0*(dpE2-dpE1)
        dpDEB=dpDE*dpHB
        if(dpDEB.le.1.e-12)then
         dpE=dpE2
         iNEL=iNEL2
         bDoThis=.false.
        end if
        if(bDoThis.eqv..true.)then
         dpE3=dpE1+dpDE
         iKK=iKK+1
         if(iKK.gt.50)then
          dpET1=dpE1*dpHB
          dpET2=dpE2*dpHB
          write(2,100)dpET1,dpSE1,iNEL1,dpET2,dpSE2,iNEL2
          iNO=-1
          return
         end if
         do i=1,iN1
          dpaP(i)=dpaU(i)+dpaRMS(i)*dpE3
         end do
         dpaS(1)=1.e-10
         dpB0=0.d0
         dpAA=dpH12*(dpaP(1)-dpE3)
         dpB1=dpaS(1)*(1.d0-dpAA)
         do i=2,iN1
          dpB2=12.d0*dpaS(i-1)-10.d0*dpB1-dpB0
          dpAA=dpH12*(dpaP(i)-dpE3)
          dpaS(i)=dpB2/(1.d0-dpAA)
          dpB0=dpB1
          dpB1=dpB2
         end do
         dpSE3=dpaS(iN)
         iNEL3=0
         do i=8,iN
          if(dpaS(i-1)*dpaS(i).lt.0.d0)then
           iNEL3=iNEL3+2
           cycle
          else if(dpaS(i-1)*dpaS(i).eq.0.d0)then
           iNEL3=iNEL3+1
          end if
         end do
         iNEL3=iNEL3/2
         if(dpSE1*dpSE3.lt.0.d0)then
          iNEL2=iNEL3
          dpSE2=dpSE3
          dpE2=dpE3
          do i=1,iN1
           dpaT(i)=dpaS(i)
          end do
          cycle
         else if(dpSE1*dpSE3.gt.0.d0)then
          iNEL1=iNEL3
          dpSE1=dpSE3
          dpE1=dpE3
          cycle
         else
          dpE=dpE3
          iNEL=iNEL3
          do i=1,iN1
           dpaT(i)=dpaS(i)
          end do
         end if
        end if
       end if
       exit
      end do
      iNO=iNEL-1
      do i=iN,4,-1
       if(dpaT(i-1)*dpaT(i).lt.0.d0)then
        iKK=i-1
        exit
       else if(dpaT(i-1)*dpaT(i).eq.0.d0)then
        iKK=i-2
        exit
       end if
      end do
      if(dpE.lt.0.d0)then
       do i=iKK,4,-1
        if(dabs(dpaT(i-1))-dabs(dpaT(i)).le.0.d0)exit
       end do
       iKL=(i+iKK)/2
       iK1=iN-1
       iK2=iN-2
       dpaS(iN)=0.d0
       dpaS(iK1)=1.e-10
       dpAA=dpH12*(dpaP(iK1)-dpE)
       dpB0=0.d0
       dpB1=dpaS(iK1)*(1.d0-dpAA)
       iK4=iK1-iKL
       do iK3=1,iK4
        i=iK2-iK3+1
        dpB2=12.d0*dpaS(i+1)-10.d0*dpB1-dpB0
        dpAA=dpH12*(dpaP(i)-dpE)
        dpaS(i)=dpB2/(1.d0-dpAA)
        dpB0=dpB1
        dpB1=dpB2
       end do
       dpFAC=dpaT(iKL)/dpaS(iKL)
       do i=iKL,iN
        dpaT(i)=dpaS(i)*dpFAC
       end do
      else
       do i=iKK,iN
        dpaT(i)=dpaT(i)-dfloat(i-iKK)*dpaT(iN)/dfloat(iN-iKK)
       end do
      end if
      dpSOM=0.d0
      do i=1,iN
       dpaT(i)=dpaT(i)*dsqrt(1.d0-dpaRMS(i))
       dpSom=dpSom+dpaT(i)*dpaT(i)
      end do
      dpSom=dsqrt(dpSom*dpStep)
      dpSig=dpaT(10)/dabs(dpaT(10))
      do i=1,iN
       dpaS(i)=dpSig*dpaT(i)/dpSom
      end do

  100 format(2X,2E12.5,I5,2E12.5,I5)

      return
      end


C     ================================================================
C     Main RPA procedure. It sets the RPA equations.
C     See Cap 2.2 for more details
C     ================================================================
      subroutine RPA()

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_DIM/iNV,iNV1,iProt,iNeutr

      dimension dpaA(IP_NCF,IP_NCF),dpaB(IP_NCF,IP_NCF)
      dimension dpaFreq(IP_NCF),dpaBel_IS(IP_NCF),dpaBel_IV(IP_NCF)
      dimension dpaFractS0_IS(IP_NCF),dpaFractS1_IS(IP_NCF)
      dimension dpaFractS0_IV(IP_NCF),dpaFractS1_IV(IP_NCF)

      call READ_RPA
      call LABEL
      call SKYRME
      call PHMATRIXELEMENTS
      call MATRIX(dpaA,dpaB)
      call SRPA(dpaA,dpaB,dpaFreq,iFlag)
      call EIGSRT(dpaFreq,dpaA,dpaB,iNV1,IP_NCF)
      if(iFlag.ne.1)then
       call ELECTRO(dpaA,dpaB,dpaFreq,dpaBel_IS,dpaBel_IV,
     &              dpaFractS0_IS,dpaFractS0_IV,
     &              dpaFractS1_IS,dpaFractS1_IV)
       if(dpaFreq(1).lt.1e-5.and.dpaFreq(1).gt.0.d0)stop
      end if

      return
      end


C     ================================================================
C     Read HF data 
C     ================================================================
      subroutine READ_RPA()

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_DIM/iNV,iNV1,iProt,iNeutr
      common/C_EP/dpEPMax
      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)
      common/C_ASI/dpaFOcc(IP_NSP),dpaSPE(IP_NSP)
      common/C_HF/iNMax,iNOcc,iNUnocc,iNOrb
      common/C_DENS/dpaDT(IP_NNP),dpaDN(IP_NNP),dpaDP(IP_NNP)
      common/C_WF/dpaWF(IP_NNP,IP_NSP),dpaDWF(IP_NNP,IP_NSP)
      common/C_SKYRME/dpT0,dpT1,dpT2,dpT3,dpX0,dpX1,dpX2,dpX3,dpW0_0,
     &       dpW2P_0,dpAlpha,iIncludeJ2
      common/C_NUC/dpANucl,dpZNucl
      common/C_RES_SO/dpW2,dpW2P
      common/C_MESHR/dpStep
      common/C_ENERG/dpaSPEP(IP_NSP),dpaSPEN(IP_NSP)

      dimension dpaTAUT(IP_NNP)

      open(unit=10,status='old',file='temp2.dat')
      open(unit=12,status='old',file='temp3.dat')
      read(10,*)iNMax
      read(10,200)dpT0,dpT1,dpT2,dpT3,dpX0
      read(10,200)dpAlpha,dpX3,dpX1,dpX2,dpW2,dpW2P
      read(10,200)dpANucl,dpZNucl,dpStep
      iProt=int(dpZNucl)
      iNeutr=int(dpANucl-dpZNucl)
      read(10,210)i1,i2,i3,i4,i5
      read(10,200)(dpaDT(i),i=1,iNMax),(dpaTAUT(i),i=1,iNMax)
      read(10,200)(dpaDP(i),i=1,iNMax),(dpaDN(i),i=1,iNMax)
      dpEHMax_P=-1000.0d0
      iHMax_P=-10
      dpEPMin_P=1000.0d0
      dpEHMax_N=-1000.0d0
      iHMax_N=-10
      dpEPMin_N=1000.0d0
      open(unit=9,status='old',file='temp1.dat')
      rewind(9)
      do i=1,iNOcc
       read(9,*)iaNN(i),iaLL(i),iaLJ(i),iaIQ(i)
       iaLev(i)=i
       read(12,*)dpaSPE(i),iaNN(i),iaLL(i),iaLJ(i),iaIQ(i)
       if(iaIQ(i).eq.1)dpaSPEP(i)=dpaSPE(i)
       if(iaIQ(i).eq.0)dpaSPEN(i)=dpaSPE(i)
       read(12,*)(dpaWF(iJ,i),iJ=1,iNMax)
       read(12,*)(dpaDWF(iJ,i),iJ=1,iNMax)
       if(iaIQ(i).eq.1)then
        if(dpaSPE(i).gt.dpEHMax_P)then
         dpEHMax_P=dpaSPE(i)
         iHMAX_P=i
        end if
       end if
       if(iaIQ(i).eq.0)then
        if(dpaSPE(i).gt.dpEHMax_N)then
         dpEHMax_N=dpaSPE(i)
         iHMAX_N=i
        end if
       end if
      end do
      i=iNOcc
      do while(.true.)
       read(9,*,end=300)iLP,iJP,iLTP,iN1,iN2
       do iI=iN1,iN2
        i=i+1
        if(i.gt.IP_NSP)then
         write(*,100)i,IP_NSP
         stop
        end if
        iaLev(i)=i
        read(12,*)dpaSPE(i),iaNN(i),iaLL(i),iaLJ(i),iaIQ(i)
        if(iI.eq.iN2.and.dpaSPE(i).lt.dpEPMax)then
         write(2,*)
     &    '>>> WARNING: THE CHOSEN VALUE OF PARTICLE ENERGY',
     &    ' CUTOFF'
         write(2,*)
     &    'IS SO LARGE THAT NOT ALL PARTICLE STATES ARE CALCULATED'
         write(2,*)
     &    '>>> INCREASE PARAMETER IP_MOREN IN FILE param.rpa'
        end if
        if(iaIQ(i).eq.1)dpaSPEP(i)=dpaSPE(i)
        if(iaIQ(i).eq.0)dpaSPEN(i)=dpaSPE(i)
        read(12,*)(dpaWF(iJ,i),iJ=1,iNMax)
        read(12,*)(dpaDWF(iJ,i),iJ=1,iNMax)
        if(iaIQ(i).eq.1)then
         if(dpaSPE(i).lt.dpEPMin_P)then
          dpEPMin_P=dpaSPE(i)
         end if
        end if
        if(iaIQ(i).eq.0)then
         if(dpaSPE(i).lt.dpEPMin_N)then
          dpEPMin_N=dpaSPE(i)
         end if
        end if
       end do
      end do
  300 continue
      do iN=1,iNOcc
       read(12,*)dpaFOcc(iN)
      end do
      if(dabs(dpaFOcc(iHMax_P)-1.0d0).gt.0.01d0)then
       dpEPMin_P=dpEHMax_P
      end if
      if(dabs(dpaFOcc(iHMax_N)-1.0d0).gt.0.01d0)then
       dpEPMin_N=dpEHMax_N
      end if
      write(2,*)
      iNOrb=i
      iNUnocc=iNOrb-iNOcc
      dpFP=(dpEPMin_P+dpEHMax_P)/2.0d0
      dpFN=(dpEPMin_N+dpEHMax_N)/2.0d0
      dpFP1=dpFP
      dpFN1=dpFN
      if(iNOrb.gt.IP_NSP)then
       write(2,110)
       stop
      else if(iNMax.gt.IP_NNP)then
       write(2,120)
       stop
      end if
      write(2,130)
      iFLAG=0
      do i=1,iNOrb
       if(iaIQ(i).eq.1)then
        if(dpaSPEP(i).gt.dpFP1.and.iFLAG.eq.0)then
         write(2,140)
         iFLAG=1
        endIF
        write(2,150)i,iaIQ(i),iaLL(i),iaLJ(i),dpaFOcc(i),dpaSPEP(i)
       end if
      end do
      iFLAG=0
      write(2,160)
      do i=1,iNOrb
       if(iaIQ(i).eq.0)then
        if(dpaSPEN(i).gt.dpFN1.and.iFLAG.eq.0)then
         write(2,140)
         iFLAG=1
        end if
        write(2,150)i,iaIQ(i),iaLL(i),iaLJ(i),dpaFOcc(i),dpaSPEN(i)
       end if
      end do
      close(12)
      close(10)

  100 format(2X,'NUMBER OF S.P. STATES I=',I3,' IS LARGER THAN',
     &       ' RESERVED dimension IP_NSP=',I3)
  110 format(/,10X,' WARNING : parameter IP_NSP ( = MAX NO',
     &       ' OF SINGLE PARTICLE STATES) IS TOO SMALL ',/)
  120 format(/,10X,' WARNING : parameter IP_NNP ( = MAX NO',
     &       ' OF POINTS FOR RADIAL QUANTITIES) IS TOO SMALL',/)
  130 format(/,'SINGLE-PARTICLE STATES:',//,
     &       5X,'- PROTON STATES -',/,
     &       5X,'    N   IQ    L     J       FOCC     SPE',/,
     &       5X,'   ---  --   --    ---      ----    -----',/)
  140 format(5X,'--FERMI SURFACE--',15('-'))
  150 format(5X,I5,3I5,'/2',2F10.3)
  160 format(/,5X,'- NEUTRON STATES -',/,
     &       5X,'    N   IQ    L     J       FOCC     SPE',/,
     &       5X,'   ---  --   --    ---      ----    -----',/)
  200 format(6E12.5)
  210 format(15I3)

      return
      end


C     ================================================================
C     Determination of the RPA p-h configurations
C     ================================================================
      subroutine LABEL()

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_DIM/iNV,iNV1,iProt,iNeutr
      common/C_CTR/dpQBess,iISPIN,iORB,iPAR,iISFLIP,iRAD
      common/C_CTRL/dpECFMin,dpECFMax
      common/C_EP/dpEPMax
      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)
      common/C_ASI/dpaFOcc(IP_NSP),dpaSPE(IP_NSP)
      common/C_HF/iNMAX,iNOcc,iNUnocc,iNOrb
      common/C_IQI/iaQI(IP_NCF)
      common/C_ENERG/dpaSPEP(IP_NSP),dpaSPEN(IP_NSP)

      write(2,100)iISPIN,iPAR
      iNV=0
      dpECF_F_Min=1000.d0
      dpECF_F_Max=0.d0
      do i1=1,iNOrb
       do i2=1,iNOrb
        iQTot=iaIQ(i1)+iaIQ(i2)
        iJMIN=iabs(iaLJ(i1)-iaLJ(i2))/2
        iJMAX=(iaLJ(i1)+iaLJ(i2))/2
        do iJT=iJMin,iJMax
         if(iQTot.eq.2)dpaECF=-dpaSPEP(i1)+dpaSPEP(i2)
         if(iQTot.eq.0)dpaECF=-dpaSPEN(i1)+dpaSPEN(i2)
         if(iQTot.eq.2)dpeTest=dpaSPEP(i2)
         if(iQTot.eq.0)dpeTest=dpaSPEN(i2)
         if(dpeTest.gt.dpEPMax)cycle
         if(dpaECF.lt.dpECFMin)cycle
         if(dpaECF.gt.dpECFMax)cycle
         if(iaIQ(i1).eq.1.and.iaIQ(i2).eq.1)then
          if(dpaFOcc(i1).eq.1.d0.and.dpaFOcc(i2).eq.1.d0)cycle
          if(dpaFOcc(i1).eq.0.d0.and.dpaFOcc(i2).eq.0.d0)cycle
         end if
         if(iaIQ(i1).eq.0.and.iaIQ(i2).eq.0)then
          if(dpaFOcc(i1).eq.1.d0.and.dpaFOcc(i2).eq.1.d0)cycle
          if(dpaFOcc(i1).eq.0.d0.and.dpaFOcc(i2).eq.0.d0)cycle
         end if
         iQTOT=iaIQ(i1)+iaIQ(i2)
         if(iQTOT.eq.1)cycle
         if(iaLev(i2).lt.iaLev(i1))cycle
         if(iJT.ne.iISPIN)cycle
         if(iPAR.eq.0)then
          iPRTY=0
         else if(iPAR.ne.0)then
          iPRTY=(-1)**(iaLL(i1)+iaLL(i2))
          if(iPRTY.ne.iPAR)cycle
         end if
         iNV=iNV+1
         if(iNV.gt.IP_NCF)then
          write(2,110)iNV
          stop
         end if
         iaIPP(iNV)=i2
         iaINN(iNV)=i1
         iaJJ(iNV)=iJT
         iaQI(iNV)=-2*iaIQ(i1)+1
         write(2,120)iNV,iaIPP(iNV),iaIQ(i2),iaLL(i2),iaLJ(i2),
     &               iaINN(iNV),iaIQ(i1),iaLL(i1),iaLJ(i1),dpaECF
         if(dpaECF.lt.dpECF_F_Min)dpECF_F_Min=dpaECF
         if(dpaECF.gt.dpECF_F_Max)dpECF_F_Max=dpaECF
        end do
       end do
      end do

  100 format(/,' PARTICLE HOLE CONFIGURATIONS',
     &       //,10X,' IISPIN =',I2,' IPAR =',I2,//,10X,
     &       '         I1 IQ1 L1 J1     I2 IQ2 L2 J2',
     &       '       ECONF',/,10X,
     &       '        --- --- -- ---   --- --- -- ---',
     &       '     --------'/)
  110 format(//,10X,' WARNING : PARAMETER IP_NCF(= MAX NO.',
     &       ' OF TWO BODY STATES) IS TOO SMALL ',/I7)
  120 format(10X,I5,3X,2(4I3,'/2',3X),F9.3)

      return
      end


C     ================================================================
C     Set up the residual force to calculate isoscalar or isovector
C     Skyrme particle-hole matrix elements in SUBPHME routine
C     ================================================================
      subroutine SKYRME()

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_HF/iNMax,iNOcc,iNUnocc,iNOrb
      common/C_SKYRME/dpT0,dpT1,dpT2,dpT3,dpX0,dpX1,dpX2,dpX3,
     &       dpW0_0,dpW2P_0,dpAlpha,iIncludeJ2
      common/C_ETALF/dpaCALF(7),dpaCBET(7),dpaCALFS(7),dpaCBETS(7),
     &       dpaCALFV(7),dpaCBETV(7)
      common/C_DENS/dpaDT(IP_NNP),dpaDN(IP_NNP),dpaDP(IP_NNP)
      common/C_FG/dpaF0(IP_NNP,3),dpaG0(IP_NNP,3),
     &       dpaCorrRes(IP_NNP),dpaF1(3),dpaG1(3)

      open(unit=14,status='unknown',file='temp4.dat')
      dpaF1(1)=0.75d0*dpT0
      dpaG1(1)=-0.25d0*dpT0*(1.d0-2.d0*dpX0)
      do iJ=1,iNMax
       dpaF0(iJ,1)=dpaF1(1)+(dpT3/48.d0)*(dpaDT(iJ)**dpAlpha)
     &            *(3.d0*(dpAlpha+1.d0)*(dpAlpha+2.d0)
     &            +dpAlpha*(1.d0-dpAlpha)*(1.d0+2.d0*dpX3)
     &            *(((dpaDN(iJ)-dpaDP(iJ))/dpaDT(iJ))**2))
       dpaG0(iJ,1)=dpaG1(1)-(dpT3/24.d0)
     &            *(1.d0-2.d0*dpX3)*(dpaDT(iJ)**dpAlpha)
       dpaCorrRes(iJ)=-(dpT3/24.d0)*(2.d0*dpX3+1.d0)*dpAlpha
     &               *(dpaDT(iJ)**(dpAlpha-1.d0))*(dpaDN(iJ)
     &               -dpaDP(iJ))
      end do
      dpaF1(2)=-0.25d0*dpT0*(1.d0+2.d0*dpX0)
      dpaG1(2)=-0.25d0*dpT0
      do iJ=1,iNMax
       dpaF0(iJ,2)=dpaF1(2)-(dpT3/24.d0)
     &            *(1.d0+2.d0*dpX3)*(dpaDT(iJ)**dpAlpha)
       dpaG0(iJ,2)=dpaG1(2)-(dpT3/24.d0)*(dpaDT(iJ)**dpAlpha)
      end do
      dpC0T=(-4.d0*dpT1+8.d0*dpT1*dpX1+4.d0*dpT2+8.d0*dpT2*dpX2)/64.d0
      dpC1T=(-4.d0*dpT1+4.d0*dpT2)/64.d0
      dpaCALFS(1)=-3.d0*dpT1/32.d0
      dpaCALFS(2)=6.d0*dpT1/32.d0
      dpaCALFS(3)=6.d0*dpT1/32.d0
      dpaCBETS(1)=-0.125d0*dpT1*(dpX1/2.d0-0.25d0)
      dpaCBETS(2)=0.25d0*dpT1*(dpX1/2.d0-0.25d0)
     &           -dfloat(1-iIncludeJ2)*dpC0T/2.d0
      dpaCBETS(3)=0.25d0*dpT1*(dpX1/2.d0-0.25d0)
     &           -dfloat(1-iIncludeJ2)*dpC0T/2.d0
      dpaCALFV(1)=(-dpT1/8.d0)*(-dpX1/2.d0-0.25d0)
      dpaCALFV(2)=-2.d0*(-dpT1/8.d0)*(-dpX1/2.d0-0.25d0)
      dpaCALFV(3)=-2.d0*(-dpT1/8.d0)*(-dpX1/2.d0-0.25d0)
      dpaCBETV(1)=-0.125d0*dpT1*(-0.25d0)
      dpaCBETV(2)=0.25d0*dpT1*(-0.25d0)-dfloat(1-iIncludeJ2)*dpC1T/2.d0
      dpaCBETV(3)=0.25d0*dpT1*(-0.25d0)-dfloat(1-iIncludeJ2)*dpC1T/2.d0
      do iY=4,7
       iYY=1-2*(iY/6)
       write(14,100)iY,(2.d0)**(-iYY),(iYY+1)/2,(2.d0)**((iYY+1)/2)
       dpaCALFS(iY)=-0.25d0*iYY*dpT2*(dpX2+5.d0/4.d0)
       dpaCBETS(iY)=-0.25d0*iYY*dpT2*(dpX2/2.d0+0.25d0)
     &             +dfloat(1-iIncludeJ2)*iYY*dpC0T/(2.d0)**((iYY+1)/2)
       dpaCALFV(iY)=-0.25d0*iYY*dpT2*(dpX2/2.d0+0.25d0)
       dpaCBETV(iY)=-0.25d0*iYY*dpT2*0.25d0
     &             +dfloat(1-iIncludeJ2)*iYY*dpC1T/(2.d0)**((iYY+1)/2)
      end do
      close(14)

  100 format(I2,5X,E15.8,/)

      return
      end


C     ================================================================
C     Calculate particle-hole matrix elements
C     of the Skyrme effective nuclear interaction,
C     coupled to total J
C     ================================================================
      subroutine PHMATRIXELEMENTS()

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_SKYRME/dpT0,dpT1,dpT2,dpT3,dpX0,dpX1,dpX2,dpX3,
     &       dpW0_0,dpW2P_0,dpAlpha,iIncludeJ2
      common/C_DIM/iNV,iNV1,iProt,iNeutr
      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)
      common/C_HF/iNMax,iNOcc,iNUnocc,iNOrb
      common/C_MESHR/dpStep
      common/C_MATEL/dpaPHME(IP_NCF,IP_NCF),dpaPHMEEX(IP_NCF,IP_NCF)
      common/C_ETALF/dpaCALF(7),dpaCBET(7),dpaCALFS(7),dpaCBETS(7),
     &       dpaCALFV(7),dpaCBETV(7)
      common/C_FG/dpaF0(IP_NNP,3),dpaG0(IP_NNP,3),
     &       dpaCorrRes(IP_NNP),dpaF1(3),dpaG1(3)
      common/C_IQI/iaQI(IP_NCF)
      common/C_POTENTIAL/iCoul_D,iCoul_E,iISO

      do i1=1,iNV
       iPRTY1=(-1)**(iaLL(iaIPP(i1))+iaLL(iaINN(i1)))
       do i2=i1,iNV
        iPRTY2=(-1)**(iaLL(iaIPP(i2))+iaLL(iaINN(i2)))
        if (iaJJ(i1).ne.iaJJ(i2))cycle
        if (iPRTY1.ne.iPRTY2)cycle
        dpaPHME(i1,i2)=0.d0
        dpaPHME(i2,i1)=0.d0
        dpaPHMEEX(i1,i2)=0.d0
        dpaPHMEEX(i2,i1)=0.d0
        dpRes=0.d0
        dpResEx=0.d0
        iPhase=iaQI(i1)*iaQI(i2)
        dpaF1(3)=dpaF1(1)+iPhase*dpaF1(2)
        dpaG1(3)=dpaG1(1)+iPhase*dpaG1(2)
        do iJ=1,iNMax
         dpaF0(iJ,3)=dpaF0(iJ,1)+iPhase*dpaF0(iJ,2)
     &              -dfloat((1+iPhase)/2)
     &              *dfloat(iaQI(i1))*2.d0*dpaCorrRes(iJ)
         dpaG0(iJ,3)=dpaG0(iJ,1)+iPhase*dpaG0(iJ,2)
        end do
        do i=1,7
         dpaCALF(i)=dpaCALFS(i)+iPhase*dpaCALFV(i)
         dpaCBET(i)=dpaCBETS(i)+iPhase*dpaCBETV(i)
        end do
        call SUBPHME(1,iaIPP(i1),iaINN(i2),iaIPP(i2),
     &               iaINN(i1),iaJJ(i1),iNMax,dpStep,dpRes)
        dpRes_Coul=0.d0
        dpRes_Coul_Ex=0.d0
        dpRes_Coul2=0.d0
        dpRes_Coul_Ex2=0.d0
        if(.not.(iaQI(i1).ne.-1.or.iaQI(i2).ne.-1))then
         if(iCoul_D.eq.1)then
          call SUBPHME_COULOMB(iaIPP(i1),iaINN(i2),iaIPP(i2),
     &                         iaINN(i1),iaJJ(i1),iNMax,
     &                         dpStep,dpRes_Coul)
          if(iCoul_E.eq.1)then
           call SUBPHME_COUL_EXC(iaIPP(i1),iaINN(i2),
     &                           iaIPP(i2),iaINN(i1),iaJJ(i1),
     &                           iNMax,dpStep,dpRes_Coul_Ex)
          end if
         end if
        end if
        dpRes_SO=0.d0
        dpRes_SO2=0.d0
        if(iISO.eq.1)then
         call SUBPHME_SO_N(iaJJ(i1),iaIPP(i1),iaINN(i2),
     &                     iaIPP(i2),iaINN(i1),dpRes_SO)
        end if
        dpRes=dpRes+dpRes_Coul+dpRes_Coul_Ex+dpRes_SO
        call SUBPHME(1,iaIPP(i1),iaIPP(i2),iaINN(i2),iaINN(i1),
     &               iaJJ(i1),iNMax,dpStep,dpResEx)
        if(iISO.eq.1)then
         call SUBPHME_SO_N(iaJJ(i1),iaIPP(i1),iaIPP(i2),
     &                     iaINN(i2),iaINN(i1),dpRes_SO2)
        end if
        if(iCoul_D.eq.1)then
         call SUBPHME_COULOMB(iaIPP(i1),iaIPP(i2),iaINN(i2),
     &                        iaINN(i1),iaJJ(i1),iNMax,
     &                        dpStep,dpRes_Coul2)
         if(iCoul_E.eq.1) then
          call SUBPHME_COUL_EXC(iaIPP(i1),iaIPP(i2),
     &                          iaINN(i2),iaINN(i1),iaJJ(i1),
     &                          iNMax,dpStep,dpRes_Coul_Ex2)
         end if
        end if
        dpResEx=dpResEx+dpRes_Coul2+dpRes_Coul_Ex2+dpRes_SO2
        dpRNorm1=1.d0
        dpRNorm2=1.d0
        if(iaIPP(i1).eq.iaINN(i1))dpRNorm1=1.d0/dsqrt(2.d0)
        if(iaIPP(i2).eq.iaINN(i2))dpRNorm2=1.d0/dsqrt(2.d0)
        dpRes=dpRNorm1*dpRNorm2*dpRes
        dpResEx=dpRNorm1*dpRNorm2*dpResEx
        dpaPHME(i1,i2)=dpRes
        dpaPHME(i2,i1)=dpRes
        dpaPHMEEX(i1,i2)=dpResEx
        dpaPHMEEX(i2,i1)=dpResEx
       end do
      end do

      return
      end


C     ================================================================
C     Contribution of the central part of the Skyrme force 
C     to p-h matrix elements
C     ================================================================
      subroutine SUBPHME(i_T3,iA,iB,iC,iD,iJT,iNR,dpDR,dpRes)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_ETALF/dpaCALF(7),dpaCBET(7),dpaCALFS(7),dpaCBETS(7),
     &       dpaCALFV(7),dpaCBETV(7)
      common/C_RINT/dpRI0,dpRIS,dpRX1,dpaRX2(2,2)
      common/C_FG/dpaF0(IP_NNP,3),dpaG0(IP_NNP,3),
     &       dpaCorrRes(IP_NNP),dpaF1(3),dpaG1(3)
      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)

      dimension dpaHH(9),dpaVV(7,3)

      dpRes=0.d0
      do i=1,9
       dpaHH(i)=0.d0
      end do
      iQA=1-2*iaIQ(iA)
      iQB=1-2*iaIQ(iB)
      iQC=1-2*iaIQ(iC)
      iQD=1-2*iaIQ(iD)
      if((iQA.eq.iQD).and.(iQB.eq.iQC))then
       dpELIAS=1.d0
       dpELIAV=dfloat(iQA*iQB)
      else
       dpELIAS=0.d0
       dpELIAV=0.d0
      end if
      if(((iQA+iQD).eq.0).and.((iQB+iQC).eq.0))then
       dpELIAV=dpELIAV+1.d0+dfloat(iQA*iQC)
      end if
      iELIAV0=(iQA-iQB)**2+(iQB-iQC)**2+(iQC-iQD)**2
      if(iELIAV0.eq.0)then
       iELIAV1=2*iQA
      else
       iELIAV1=0
      end if
      do iJ=1,IP_NNP
       dpaF0(iJ,3)=dpaF0(iJ,1)*dpELIAS+dpELIAV*dpaF0(iJ,2)
     &            +float(iELIAV1)*dpaCorrRes(iJ)
       dpaG0(iJ,3)=dpaG0(iJ,1)*dpELIAS+dpELIAV*dpaG0(iJ,2)
      end do
      call RDINT1(i_T3,iA,iB,iC,iD,iNR,dpDR)
      dpDJP11=2.d0*iJT+1.d0
      dpSLA=dfloat(iaLL(iA))
      dpSLB=dfloat(iaLL(iB))
      dpSLC=dfloat(iaLL(iC))
      dpSLD=dfloat(iaLL(iD))
      dpSJA=dfloat(iaLJ(iA))/2.d0
      dpSJB=dfloat(iaLJ(iB))/2.d0
      dpSJC=dfloat(iaLJ(iC))/2.d0
      dpSJD=dfloat(iaLJ(iD))/2.d0
      dpSJT=dfloat(iJT)
      dpaHH(8)=dpRI0*dpf_YL(iaLL(iA),iaLJ(iA),iaLL(iD),iaLJ(iD),iJT)
     &        *dpf_YL(iaLL(iC),iaLJ(iC),iaLL(iB),iaLJ(iB),iJT)/dpDJP11
      iLADM=iabs(iaLL(iA)-iaLL(iD))
      iLADP=iaLL(iA)+iaLL(iD)
      iLBCM=iabs(iaLL(iB)-iaLL(iC))
      iLBCP=iaLL(iB)+iaLL(iC)
      iL1=max0(iLADM,iLBCM,iabs(iJT-1))
      iL2=min0(iLADP,iLBCP,iJT+1)
      dpX9=0.d0
      if(iL1.le.iL2)then
       do iLL=iL1,iL2
        dpX9=dpX9+dpf_TJL(iaLL(iA),iaLJ(iA),iaLL(iD),iaLJ(iD),iJT,iLL)
     &      *dpf_TJL(iaLL(iC),iaLJ(iC),iaLL(iB),iaLJ(iB),iJT,iLL)
       end do
      end if
      dpaHH(9)=dpRIS*dpX9/dpDJP11
      iL1=max0(iLADM,iLBCM)
      iL2=min0(iLADP,iLBCP)
      if(iL1.le.iL2)then
       do i1=1,7
        do i2=1,3
         dpaVV(i1,i2)=0.d0
        end do
       end do
       dpCOFS=(iaLJ(iA)+1.d0)*(iaLJ(iB)+1.d0)
     &       *(iaLJ(iC)+1.d0)*(iaLJ(iD)+1.d0)
       dpCOFS=dsqrt(dpCOFS)*6.d0
       iS=iaLL(iA)+iaLL(iC)+(iaLJ(iB)+iaLJ(iD))/2+1
       dpSI=1.d0-2.d0*mod(iS,2)
       call SIXJ(dpSJA,dpSLA,0.5d0,dpSLD,dpSJD,dpSJT,dpC1)
       call SIXJ(dpSJC,dpSLC,0.5d0,dpSLB,dpSJB,dpSJT,dpC2)
       dpC1=dpC1*dpC2
       dpCOF0=dpSI*dpCOFS*dpC1/6.d0
       iLPT=1
       iLGD=3
       do iL=iLPT,iLGD
        iLL=iJT-2+iL
        if(iLL.lt.0)cycle
        if((iLL-iL1)*(iLL-iL2).gt.0)cycle
        dpaVV(1,iL)=dpRX1*dpf_REDY(iaLL(iA),iLL,iaLL(iD))
     &             *dpf_REDY(iaLL(iC),iLL,iaLL(iB))/(2.d0*iLL+1.d0)
       end do
       do iUP=2,5
        if(iUP.eq.1.or.iUP.eq.2)then
         iKA=iA
         iKB=iB
         iKC=iC
         iKD=iD
        else if(iUP.eq.3)then
         iKA=iD
         iKB=iC
         iKC=iB
         iKD=iA
        else if(iUP.eq.4)then
         iKA=iA
         iKB=iC
         iKC=iB
         iKD=iD
        else if(iUP.eq.5)then
         iKA=iB
         iKB=iD
         iKC=iA
         iKD=iC
        end if
        dpSLKA=dfloat(iaLL(iKA))
        dpSLKB=dfloat(iaLL(iKB))
        dpSLKC=dfloat(iaLL(iKC))
        dpSLKD=dfloat(iaLL(iKD))
        call RDINT2(iKA,iKB,iKC,iKD,iNR,dpDR)
        dpXLKAB=(2.d0*iaLL(iKA)+1.d0)*(2.d0*iaLL(iKB)+1.d0)
        dpXLKAB=dsqrt(dpXLKAB)
        do iL=iLPT,iLGD
         iLL=iJT-2+iL
         dpSL=dfloat(iLL)
         if(iLL.lt.0)cycle
         if((iLL-iL1)*(iLL-iL2).gt.0)cycle
         dpSGN1=1.d0-2.d0*mod(iLL,2)
         if(iUP.eq.4.or.iUP.eq.5)dpSGN1=1.d0
         do iK1=1,2
          iKK1=iaLL(iKA)-3+2*iK1
          if(iKK1.lt.0)cycle
          dpSLK1=dfloat(iKK1)
          call TROISJ(dpSLK1,dpSLKA,1.d0,0.d0,0.d0,0.d0,dpC1)
          dpCOFK1=(1.d0-2.d0*mod(iKK1,2))
     &           *dsqrt(2.d0*iKK1+1.d0)*dpC1
          do iK2=1,2
           ikK2=iaLL(iKB)-3+2*iK2
           if(ikK2.lt.0)cycle
           dpSLK2=dfloat(ikK2)
           call TROISJ(dpSLK2,dpSLKB,1.d0,0.d0,0.d0,0.d0,dpC2)
           dpCOFK2=(1.d0-2.d0*mod(ikK2,2))
     &            *dsqrt(2.d0*ikK2+1.d0)*dpC2
           iLLM=max0(iabs(iLL-1),iabs(iaLL(iKD)-iKK1),
     &               iabs(iaLL(iKC)-ikK2))
           iLLP=min0(iLL+1,iaLL(iKD)+iKK1,iaLL(iKC)+ikK2)
           if(iLLM.gt.iLLP)cycle
           do iL1L=iLLM,iLLP
            dpSLL=dfloat(iL1L)
            call SIXJ(dpSLL,1.d0,dpSL,dpSLKA,dpSLKD,dpSLK1,dpC1)
            call SIXJ(dpSLL,1.d0,dpSL,dpSLKB,dpSLKC,dpSLK2,dpC2)
            dpC3=dpC1*dpC2
            dpaVV(iUP,iL)=dpaVV(iUP,iL)+dpSGN1*dpXLKAB*dpCOFK1
     &                   *dpCOFK2*dpC3*dpf_REDY(iaLL(iKD),iL1L,iKK1)
     &                   *dpf_REDY(iaLL(iKC),iL1L,ikK2)
     &                   *dpaRX2(iK1,iK2)
           end do
          end do
         end do
        end do
       end do
       do iUP=6,7
        iKA=iA
        iKB=iB
        iKC=iC
        iKD=iD
        if(iUP.ne.6)then
         iKA=iB
         iKB=iA
         iKC=iD
         iKD=iC
        end if
        dpSLKA=dfloat(iaLL(iKA))
        dpSLKD=dfloat(iaLL(iKD))
        call RDINT2(iKA,iKD,iKB,iKC,iNR,dpDR)
        dpXLKAD=(2.d0*iaLL(iKA)+1.d0)*(2.d0*iaLL(iKD)+1.d0)
        dpXLKAD=dsqrt(dpXLKAD)
        do iL=iLPT,iLGD
         iLL=iJT-2+iL
         dpSL=dfloat(iLL)
         if(iLL.lt.0)cycle
         if((iLL-iL1)*(iLL-iL2).gt.0)cycle
         dpCOFL=dpf_REDY(iaLL(iKB),iLL,iaLL(iKC))/(2.d0*iLL+1.d0)
         do iK1=1,2
          iKK1=iaLL(iKA)-3+2*iK1
          if(iKK1.lt.0)cycle
          dpSLK1=dfloat(iKK1)
          call TROISJ(dpSLK1,dpSLKA,1.d0,0.d0,0.d0,0.d0,dpC1)
          dpCOFK1=dsqrt(2.d0*iKK1+1.d0)*dpC1
          do iK2=1,2
           ikK2=iaLL(iKD)-3+2*iK2
           if(ikK2.lt.0)cycle
           dpSLK2=dfloat(ikK2)
           call TROISJ(dpSLK2,dpSLKD,1.d0,0.d0,0.d0,0.d0,dpC2)
           dpCOFK2=dsqrt(2.d0*ikK2+1.d0)*dpC2
           call SIXJ(dpSLK2,dpSLK1,dpSL,dpSLKA,dpSLKD,1.d0,dpC3)
           dpCOFK12=dpf_REDY(iKK1,iLL,ikK2)*dpC3
           dpaVV(iUP,iL)=dpaVV(iUP,iL)
     &                  +dpXLKAD*dpCOFL*dpCOFK1*dpCOFK2*dpCOFK12
     &                  *dpaRX2(iK1,iK2)
          end do
         end do
        end do
       end do
       do iL=iLPT,iLGD
        iLL=iJT-2+iL
        if(iLL.lt.0)cycle
        dpSL=dfloat(iLL)
        dpDLP1=2.d0*iLL+1.d0
        call NEUFJ(dpSJA,dpSLA,0.5d0,dpSJD
     &            ,dpSLD,0.5d0,dpSJT,dpSL,1.d0,dpC3)
        call NEUFJ(dpSJC,dpSLC,0.5d0,dpSJB
     &            ,dpSLB,0.5d0,dpSJT,dpSL,1.d0,dpC4)
        dpC4=dpC3*dpC4
        do iUP=1,7
         dpaHH(iUP)=dpaHH(iUP)+dpCOFS*dpDLP1*dpC4*dpaVV(iUP,iL)
        end do
       end do
       do iUP=1,7
        dpaCBET(iUP)=dpELIAS*dpaCBETS(iUP)+dpELIAV*dpaCBETV(iUP)
        dpaCALF(iUP)=dpELIAS*dpaCALFS(iUP)+dpELIAV*dpaCALFV(iUP)
        dpaHH(iUP)=dpaCBET(iUP)*dpaHH(iUP)
     &            +dpaCALF(iUP)*dpCOF0*dpaVV(iUP,2)
       end do
      end if
      do iUP=1,9
       dpRes=dpRes+dpaHH(iUP)
      end do
      iExp=(-iaLL(iA)-iaLL(iB)+iaLL(iC)+iaLL(iD))/2
      dpRes=if_POT(iExp)*dpRes

      return
      end


C     ================================================================
C     Three-j symbols
C     ================================================================
      subroutine TROISJ(dpXJ1,dpXJ2,dpXJ3,dpXM1,dpXM2,dpXM3,dpC3J)

      implicit integer(i)
      implicit double precision(d)
      logical bExecute

      dimension dpaFLog(301)

      data(dpaFLog(i),i=2,31)/0.d0,.69314718d0,1.7917595d0,3.1780538d0,
     A4.7874917d0,6.5792511d0,8.5251613d0,10.604603d0,12.801827d0,
     B15.104413d0,17.502307d0,19.987214d0,22.552163d0,25.191221d0,
     C27.899271d0,30.671860d0,33.505072d0,36.395445d0,39.339884d0,
     D42.335616d0,45.380139d0,48.471180d0,51.606674d0,54.784729d0,
     E58.003604d0,61.261702d0,64.557537d0,67.889743d0,71.257038d0,
     F74.658235d0/
      data(dpaFLog(i),i=32,61)/78.092223d0,81.557959d0,85.054466d0,
     A88.580827d0,92.136175d0,95.719694d0,99.330612d0,102.96820d0,
     B106.63176d0,110.32064d0,114.03421d0,117.77188d0,121.53308d0,
     C125.31727d0,129.12393d0,132.95257d0,136.80272d0,140.67392d0,
     D144.56574d0,148.47776d0,152.40959d0,156.36083d0,160.33112d0,
     E164.32011d0,168.32744d0,172.35279d0,176.39584d0,180.45629d0,
     F184.53383d0,188.62817d0/
      data(dpaFLog(i),i=62,91)/192.73904d0,196.86618d0,201.00931d0,
     A205.16820d0,209.34258d0,213.53224d0,217.73693d0,221.95644d0,
     B226.19054d0,230.43904d0,234.70172d0,238.97839d0,243.26885d0,
     C247.57291d0,251.89040d0,256.22113d0,260.56494d0,264.92164d0,
     D269.29110d0,273.67312d0,278.06757d0,282.47429d0,286.89313d0,
     E291.32394d0,295.76659d0,300.22094d0,304.68685d0,309.16419d0,
     F313.65283d0,318.15264d0/
      data(dpaFLog(i),i=92,121)/322.66349d0,327.18529d0,331.71788d0,
     A336.26118d0,340.81505d0,345.37940d0,349.95411d0,354.53908d0,
     B359.13420d0,363.73937d0,368.35449d0,372.97946d0,377.61419d0,
     C382.25859d0,386.91255d0,391.57598d0,396.24881d0,400.93094d0,
     D405.62230d0,410.32277d0,415.03230d0,419.75080d0,424.47819d0,
     E429.21439d0,433.95932d0,438.71291d0,443.47508d0,448.24576d0,
     F453.02489d0,457.81238d0/
      data(dpaFLog(i),i=122,151)/462.60817d0,467.41220d0,472.22438d0,
     A477.04466d0,481.87298d0,486.70926d0,491.55345d0,496.40547d0,
     B501.26529d0,506.13282d0,511.00802d0,515.89082d0,520.78117d0,
     C525.67901d0,530.58428d0,535.49694d0,540.41692d0,545.34417d0,
     D550.27865d0,555.22029d0,560.16905d0,565.12488d0,570.08772d0,
     E575.05753d0,580.03427d0,585.01787d0,590.00830d0,595.00552d0,
     F600.00946d0,605.02010d0/
      data(dpaFLog(i),i=152,181)/610.03738d0,615.06126d0,620.09170d0,
     A625.12866d0,630.17208d0,635.22193d0,640.27818d0,645.34077d0,
     B650.40968d0,655.48486d0,660.56626d0,665.65385d0,670.74760d0,
     C675.84747d0,680.95341d0,686.06541d0,691.18340d0,696.30735d0,
     D3701.43726d0,706.57306d0,711.71472d0,716.86221d0,722.01551d0,
     E727.17456d0,732.33934d0,737.50983d0,742.68598d0,747.86776d0,
     F753.05516d0,758.24811d0/
      data(dpaFLog(i),i=182,211)/763.44661d0,768.65061d0,773.86010d0,
     A779.07503d0,784.29539d0,789.52114d0,794.75224d0,799.98869d0,
     B805.23044d0,810.47747d0,815.72973d0,820.98722d0,826.24991d0,
     C831.51778d0,836.79078d0,842.06890d0,847.35209d0,852.64036d0,
     D857.93366d0,863.23199d0,868.53529d0,873.84356d0,879.15676d0,
     E884.47488d0,889.79789d0,895.12577d0,900.45848d0,905.79603d0,
     F911.13836d0,916.48547d0/
      data(dpaFLog(i),i=212,241)/921.83732d0,927.19391d0,932.55521d0,
     A937.92118d0,943.29181d0,948.66710d0,954.04699d0,959.43148d0,
     B964.82056d0,970.21419d0,975.61235d0,981.01503d0,986.42220d0,
     C991.83385d0,997.24995d0,1002.6705d0,1008.0954d0,1013.5248d0,
     D1018.9585d0,1024.3966d0,1029.8389d0,1035.2857d0,1040.7367d0,
     E1046.1920d0,1051.6516d0,1057.1155d0,1062.5836d0,1068.0558d0,
     F1073.5323d0,1079.0129d0/
      data(dpaFLog(i),i=242,271)/1084.4977d0,1089.9866d0,1095.4797d0,
     A1100.9768d0,1106.4781d0,1111.9834d0,1117.4928d0,1123.0063d0,
     B1128.5237d0,1134.0452d0,1139.5706d0,1145.1001d0,1150.6335d0,
     C1156.1708d0,1161.7120d0,1167.2573d0,1172.8063d0,1178.3593d0,
     D1183.9161d0,1189.4768d0,1195.0413d0,1200.6097d0,1206.1818d0,
     E1211.7577d0,1217.3375d0,1222.9209d0,1228.5082d0,1234.0992d0,
     F1239.6939d0,1245.2924d0/
      data(dpaFLog(i),i=272,301)/1250.8944d0,1256.5003d0,1262.1097d0,
     A1267.7228d0,1273.3396d0,1278.9600d0,1284.5840d0,1290.2117d0,
     B1295.8429d0,1301.4777d0,1307.1160d0,1312.7580d0,1318.4034d0,
     C1324.0524d0,1329.7048d0,1335.3609d0,1341.0203d0,1346.6833d0,
     D1352.3497d0,1358.0196d0,1363.6929d0,1369.3697d0,1375.0499d0,
     E1380.7334d0,1386.4204d0,1392.1107d0,1397.8045d0,1403.5016d0,
     F1409.2020d0,1414.9058d0/
      data dpEPS1,dpEPS2/0.1d0,-0.2d0/
 
      dpXN=dpXJ2-dpXM2+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN1=int(dpXN)
      dpXIII=dabs(dpXJ1)+dabs(dpXJ2)+dabs(dpXJ3)
      if(dpXIII.ge.200)then
       write(9,*)'XJ1 XJ2 XJ3=',dpXJ1,dpXJ2,dpXJ3
       stop
      end if
      dpXN=dpXJ3+dpXM3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN2=int(dpXN)
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      dpXN=dpXJ3-dpXM3+dpEPS1
      iN3=int(dpXN)
      dpXN=dpXJ1+dpXM1+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN4=int(dpXN)
      dpXN=dpXJ2+dpXJ3-dpXJ1+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN5=int(dpXN)
      dpXN=dpXJ1+dpXJ3-dpXJ2+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN6=int(dpXN)
      dpXN=dpXJ1+dpXJ2-dpXJ3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN7=int(dpXN)
      dpXN=dpXJ1-dpXM1+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN8=int(dpXN)
      dpXN=dpXJ2+dpXM2+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN9=int(dpXN)
      dpXN=dpXJ3-dpXJ1-dpXM2+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN10=int(dpXN)
      dpXN=dpXJ3-dpXJ2+dpXM1+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN11=int(dpXN)
      dpXN=dpXJ1+dpXJ2+dpXJ3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN12=int(dpXN)+1
      iK=iN4*iN8
      if(iK.lt.0)then
       dpC3J=0.d0
       return
      end if
      iK=iN1*iN9
      if(iK.lt.0)then
       dpC3J=0.d0
       return
      end if
      iK=iN2*iN3
      if(iK.lt.0.or.iN5.lt.0.or.iN6.lt.0.or.iN7.lt.0)then
       dpC3J=0.d0
       return
      end if
      iL=iN1-iN2+iN3-iN4+iN8-iN9
      if(iL.ne.0)then
       dpC3J=0.d0
       return
      end if
      iK=iN12-1
      if(iK.le.0)then
       dpC3J=1.d0
       return
      end if
      iK=0
      iL=-iN10
      if(iL.gt.iK)iK=iL
      iL=-iN11
      if(iL.gt.iK)iK=iL
      iL=iN7
      if(iN8.gt.iL)iL=iN8
      if(iN9.gt.iL)iL=iN9
      dpF=1.d0
      dpS=1.d0
      i=iK+1
      do while(i.le.iL)
       iM1=i-1
       iNN=(iN7-iM1)*(iN8-iM1)*(iN9-iM1)
       iNND=i*(iN10+i)*(iN11+i)
       dpF=-dpF*dfloat(iNN)/dfloat(iNND)
       dpS=dpS+dpF
       i=i+1
      end do
      dpC2N=dpaFLog(iN1+1)+dpaFLog(iN2+1)+dpaFLog(iN3+1)
     &     +dpaFLog(iN4+1)+dpaFLog(iN5+1)+dpaFLog(iN6+1)
     &     +dpaFLog(iN7+1)+dpaFLog(iN8+1)+dpaFLog(iN9+1)
      dpC2N=0.5d0*dpC2N
      iKM1=iK-1
      iKP1=iK+1
      dpC2D=dpaFLog(iKP1)+dpaFLog(iN7-iKM1)+dpaFLog(iN8-iKM1)
     &     +dpaFLog(iN9-iKM1)+dpaFLog(iN10+iKP1)
     &     +dpaFLog(iN11+iKP1)+0.5d0*dpaFLog(iN12+1)
      dpF=dpC2D-dpC2N
      bExecute=.true.
      if(dpF.le.80.d0)then
       dpF=dpC2N/dpC2D
       if(.not.((dpF.lt.1.01d0).and.(dpF.gt.0.98d0)))then
        dpC3J=dpS*dexp(dpC2N-dpC2D)
        bExecute=.false.
       end if
      end if
      if(bExecute.eqv..true.)then
       if(dpS.eq.0)then
        dpC3J=0.d0
        return
       else if(dpS.lt.0.d0)then
        dpS=dlog(-dpS)
        dpC3J=-dexp(dpS+dpC2N-dpC2D)
       else
        dpS=dlog(dpS)
        dpC3J=dexp(dpS+dpC2N-dpC2D)
       end if
      end if
      dpXN=dpXJ1-dpXJ2-dpXM3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iL=int(dpXN)
      iL=iL+iK
      iK=iL/2
      iK=2*iK
      if(iL.ne.iK)dpC3J=-dpC3J

      return
      end


C     ================================================================
      subroutine RDINT1(i_T3,iA,iB,iC,iD,iNR,dpDR)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_RINT/dpRI0,dpRIS,dpRX1,dpaRX2(2,2)
      common/C_WF/dpaWF(IP_NNP,IP_NSP),dpaDWF(IP_NNP,IP_NSP)
      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)
      common/C_FG/dpaF0(IP_NNP,3),dpaG0(IP_NNP,3),
     &       dpaCorrRes(IP_NNP),dpaF1(3),dpaG1(3)

      dpRI0=0.d0
      dpRIS=0.d0
      dpRX1=0.d0
      dpXL=iaLL(iA)*(iaLL(iA)+1)+iaLL(iB)*(iaLL(iB)+1)
     &    +iaLL(iC)*(iaLL(iC)+1)+iaLL(iD)*(iaLL(iD)+1)
      do iJ=1,iNR
         dpXX1=iJ*dpDR
         dpXX2=dpXX1*dpXX1
         dpC=dpaWF(iJ,iA)*dpaWF(iJ,iB)*dpaWF(iJ,iC)*dpaWF(iJ,iD)
     &      *dpDR/dpXX2
         dpRI0=dpRI0+dfloat(i_T3)*dpaF0(iJ,3)
     &        *dpC+dfloat(1-i_T3)*dpaF1(3)*dpC
         dpRIS=dpRIS+dfloat(i_T3)*dpaG0(iJ,3)
     &        *dpC+dfloat(1-i_T3)*dpaG1(3)*dpC
         dpC1=-dpXL*dpC/dpXX2+2.d0*dpDR*((dpaDWF(iJ,iA)*dpaWF(iJ,iB)
     &       +dpaWF(iJ,iA)*dpaDWF(iJ,iB))
     &       *dpaWF(iJ,iC)*dpaWF(iJ,iD)+(dpaDWF(iJ,iC)*dpaWF(iJ,iD)
     &       +dpaWF(iJ,iC)*dpaDWF(iJ,iD))
     &       *dpaWF(iJ,iA)*dpaWF(iJ,iB))/(dpXX1*dpXX2)
         dpC1=dpC1-2.d0*dpDR*(dpaDWF(iJ,iA)*(dpaDWF(iJ,iB)
     &       *dpaWF(iJ,iC)*dpaWF(iJ,iD)
     &       +dpaWF(iJ,iB)*(dpaDWF(iJ,iC)*dpaWF(iJ,iD)
     &       +dpaWF(iJ,iC)*dpaDWF(iJ,iD)))
     &       +dpaWF(iJ,iA)*(dpaDWF(iJ,iB)*(dpaDWF(iJ,iC)*dpaWF(iJ,iD)
     &       +dpaWF(iJ,iC)*dpaDWF(iJ,iD))
     &       +dpaWF(iJ,iB)*dpaDWF(iJ,iC)*dpaDWF(iJ,iD)))/dpXX2
         dpRX1=dpRX1+dpC1
      end do

      return
      end


C     ================================================================
      subroutine RDINT2(iA,iB,iC,iD,iNR,dpDR)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_RINT/dpRI0,dpRIS,dpRX1,dpaRX2(2,2)
      common/C_WF/dpaWF(IP_NNP,IP_NSP),dpaDWF(IP_NNP,IP_NSP)
      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)

      dimension dpaUA(2),dpaUB(2)

      do iK1=1,2
       do iK2=1,2
        dpaRX2(iK1,iK2)=0.d0
       end do
      end do
      do iJ=1,iNR
       dpXX1=iJ*dpDR
       dpXX2=dpXX1*dpXX1
       do i=1,2
        dpSK=3.d0-2.d0*i
        dpaUA(i)=dpaDWF(iJ,iA)+dpSK*(iaLL(iA)+i-1)*dpaWF(iJ,iA)/dpXX1
        dpaUB(i)=dpaDWF(iJ,iB)+dpSK*(iaLL(iB)+i-1)*dpaWF(iJ,iB)/dpXX1
       end do
       do iK1=1,2
        do iK2=1,2
         dpC23=dpaUA(iK1)*dpaUB(iK2)*dpaWF(iJ,iC)*dpaWF(iJ,iD)
         dpaRX2(iK1,iK2)=dpaRX2(iK1,iK2)+dpC23*dpDR/dpXX2
        end do
       end do
      end do

      return
      end


C     ================================================================
C     Contribution of the direct Coulomb interaction
C     to p-h matrix elements
C     ================================================================
      subroutine SUBPHME_COULOMB(iA,iB,iC,iD,iJT,iNR,dpDR,dpRes)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)
      common/C_WF/dpaWF(IP_NNP,IP_NSP),dpaDWF(IP_NNP,IP_NSP)

      data dpE2/1.43986d0/

      dpRes=0.d0
      if(.not.(iaIQ(iA).ne.1.or.iaIQ(iB).ne.1.or.iaIQ(iC).ne.1
     &   .or.iaIQ(iD).ne.1))then
       dpDblInt=0.d0
       do i1=1,iNR
        dpR1=dpDR*dfloat(i1)
        do i2=1,iNR
         dpR2=dpDR*dfloat(i2)
         dpRMax=max(dpR1,dpR2)
         dpRMin=min(dpR1,dpR2)
         dpPS4=dpaWF(i1,iA)*dpaWF(i2,iB)*dpaWF(i2,iC)*dpaWF(i1,iD)
         dpDblInt=dpDblInt+dpPS4*dpRMin**iJT/dpRMax**(iJT+1)
        end do
       end do
       dpDblInt=dpDblInt*dpDR*dpDR
       dpRes=dpE2*16.d0*datan(1.d0)*dpDblInt
     &      *dpf_YL(iaLL(iA),iaLJ(iA),iaLL(iD),iaLJ(iD),iJT)
     &      *dpf_YL(iaLL(iC),iaLJ(iC),iaLL(iB),iaLJ(iB),iJT)
     &      /dfloat(2*iJT+1)**2
       iExp=(-iaLL(iA)-iaLL(iB)+iaLL(iC)+iaLL(iD))/2
       dpRes=if_POT(iExp)*dpRes
      end if

      return
      end


C     ================================================================
C     Contribution of the exchange Coulomb interaction
C     to p-h matrix elements
C     ================================================================
      subroutine SUBPHME_COUL_EXC(iA,iB,iC,iD,iJT,iNR,dpDR,dpRes)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_UNITS/dpPI,dpHBARC,dpAMC2,dpHBDM,dpPMass,dpNMass
      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)
      common/C_WF/dpaWF(IP_NNP,IP_NSP),dpaDWF(IP_NNP,IP_NSP)
      common/C_DENS/dpaDT(IP_NNP),dpaDN(IP_NNP),dpaDP(IP_NNP)

      data dpE2/1.43986d0/

      dpRes=0.d0
      if(.not.(iaIQ(iA).ne.1.or.iaIQ(iB).ne.1
     &   .or.iaIQ(iC).ne.1.or.iaIQ(iD).ne.1))then
       dpC=-(dpE2/3.d0)*(3.d0/dpPI)**(1.d0/3.d0)
       dpDJP1=2.d0*iJT+1.d0
       dpRI=0.d0
       do i=1,iNR
        dpR2=dpDR*dfloat(i)
        dpR2=dpR2*dpR2
        dpPS4=dpaWF(i,iA)*dpaWF(i,iB)*dpaWF(i,iC)*dpaWF(i,iD)
        dpRI=dpRI+dpDR*dpPS4/dpR2*dpaDP(i)**(-2.d0/3.d0)
       end do
       dpRI=dpRI*dpC
       dpRes=dpRI*dpf_YL(iaLL(iA),iaLJ(iA),iaLL(iD),iaLJ(iD),iJT)
     &      *dpf_YL(iaLL(iC),iaLJ(iC),iaLL(iB),iaLJ(iB),iJT)/dpDJP1
       iExp=(-iaLL(iA)-iaLL(iB)+iaLL(iC)+iaLL(iD))/2
       dpRes=if_POT(iExp)*dpRes
      end if

      return
      end


C     ================================================================
C     Contribution of the spin-orbit term of the Skyrme force
C     to p-h matrix elements
C     ================================================================
      subroutine SUBPHME_SO_N(iJ,iA,iB,iC,iD,dpSOPHME)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)
      common/C_RES_SO/dpW2,dpW2P

      iTA=1-2*iaIQ(iA)
      iTB=1-2*iaIQ(iB)
      iTC=1-2*iaIQ(iC)
      iTD=1-2*iaIQ(iD)
      dpSJ=dfloat(iJ)
      dpSLA=dfloat(iaLL(iA))
      dpSLB=dfloat(iaLL(iB))
      dpSJA=dfloat(iaLJ(iA))*0.5d0
      dpSJB=dfloat(iaLJ(iB))*0.5d0
      dpSLC=dfloat(iaLL(iC))
      dpSLD=dfloat(iaLL(iD))
      dpSJC=dfloat(iaLJ(iC))*0.5d0
      dpSJD=dfloat(iaLJ(iD))*0.5d0
      dpVPHJ=0.d0
      dpVPHJ_F1=0.0d0
      call SUBPHME_SO_F1(iJ,iA,iB,iC,iD,dpVPHJ_F1)
      dpC6J=0.0d0
      call SIXJ(dpSJC,dpSJB,dpSJ,dpSLB,dpSLC,0.5d0,dpC6J)
      dpVPHJ_F1=dpVPHJ_F1*if_POT((iaLJ(iB)-1)/2
     &         -iaLL(iB))*dsqrt(3.d0)*dpC6J
     &         *dsqrt(2.d0)*dsqrt(6.d0)*dsqrt(iaLJ(iA)+1.0d0)
     &         *dsqrt(iaLJ(iB)+1.0d0)*dsqrt(iaLJ(iC)+1.0d0)
     &         *dsqrt(iaLJ(iD)+1.0d0)/4.d0
      dpVPHJ_F8=0.0d0
      call SUBPHME_SO_F1(iJ,iB,iA,iD,iC,dpVPHJ_F8)
      dpC6J=0.0d0
      call SIXJ(dpSJD,dpSJA,dpSJ,dpSLA,dpSLD,0.5d0,dpC6J)
      dpVPHJ_F8=dpVPHJ_F8*if_POT((iaLJ(iB)+iaLJ(iC)+iaLJ(iD)+1)/2
     &         -iaLL(iA))*dsqrt(3.d0)*dpC6J
     &         *dsqrt(6.d0)*dsqrt(2.d0)*dsqrt(iaLJ(iA)+1.0d0)
     &         *dsqrt(iaLJ(iB)+1.0d0)*dsqrt(iaLJ(iC)+1.0d0)
     &         *dsqrt(iaLJ(iD)+1.0d0)/4.d0
      dpVPHJ_F3=0.0d0
      call SUBPHME_SO_F3(iJ,iA,iB,iC,iD,dpVPHJ_F3)
      dpC6J=0.0d0
      call SIXJ(dpSLA,dpSLD,dpSJ,dpSJD,dpSJA,0.5d0,dpC6J)
      dpVPHJ_F3=dpVPHJ_F3*if_POT((-iaLJ(iD)+1)/2
     &         -iaLL(ID)+iJ)*dsqrt(3.d0)*dpC6J
     &         *dsqrt(6.d0)*dsqrt(2.d0)*dsqrt(iaLJ(iA)+1.0d0)
     &         *dsqrt(iaLJ(iB)+1.0d0)*dsqrt(iaLJ(iC)+1.0d0)
     &         *dsqrt(iaLJ(iD)+1.0d0)/4.d0
      dpVPHJ_F6=0.0d0
      call SUBPHME_SO_F3(iJ,iB,iA,iD,iC,dpVPHJ_F6)
      dpC6J=0.0d0
      call SIXJ(dpSLB,dpSLC,dpSJ,dpSJC,dpSJB,0.5d0,dpC6J)
      dpVPHJ_F6=dpVPHJ_F6*if_POT((iaLJ(iA)+iaLJ(iD)-iaLJ(iB)-1)/2
     &         +iaLL(iC)+iJ)*dsqrt(3.d0)*dpC6J
     &         *dsqrt(6.d0)*dsqrt(2.d0)*dsqrt(iaLJ(iA)+1.0d0)
     &         *dsqrt(iaLJ(iB)+1.0d0)*dsqrt(iaLJ(iC)+1.0d0)
     &         *dsqrt(iaLJ(iD)+1.0d0)/4.d0
      dpVPHJ_F2=0.0d0
      call SUBPHME_SO_F2(iJ,iA,iB,iC,iD,dpVPHJ_F2)
      dpC6J=0.0d0
      call SIXJ(dpSLB,dpSLC,dpSJ,dpSJC,dpSJB,0.5d0,dpC6J)
      dpVPHJ_F2=dpVPHJ_F2*if_POT((-iaLJ(iB)+1)/2
     &         -iaLL(iC)+iJ)*dpC6J*dsqrt(3.d0)
     &         *dsqrt(6.d0)*dsqrt(2.d0)*dsqrt(iaLJ(iA)+1.0d0)
     &         *dsqrt(iaLJ(iB)+1.0d0)*dsqrt(iaLJ(iC)+1.0d0)
     &         *dsqrt(iaLJ(iD)+1.0d0)/4.d0
      dpVPHJ_F7=0.0d0
      call SUBPHME_SO_F2(iJ,iB,iA,iD,iC,dpVPHJ_F7)
      dpC6J=0.0d0
      call SIXJ(dpSLD,dpSLA,dpSJ,dpSJA,dpSJD,0.5d0,dpC6J)
      dpVPHJ_F7=dpVPHJ_F7*if_POT((iaLJ(iC)+iaLJ(iB)-iaLJ(iD)+1)/2
     &         +iaLL(iD)+iJ+1)*dpC6J*dsqrt(3.d0)
     &         *dsqrt(6.d0)*dsqrt(2.d0)*dsqrt(iaLJ(iA)+1.0d0)
     &         *dsqrt(iaLJ(iB)+1.0d0)*dsqrt(iaLJ(iC)+1.0d0)
     &         *dsqrt(iaLJ(iD)+1.0d0)/4.d0
      dpVPHJ_F4=0.0d0
      call SUBPHME_SO_F4(iJ,iA,iB,iC,iD,dpVPHJ_F4)
      dpC6J=0.0d0
      call SIXJ(dpSLA,dpSLD,dpSJ,dpSJD,dpSJA,0.5d0,dpC6J)
      dpVPHJ_F4=dpVPHJ_F4*if_POT((iaLJ(iC)+iaLJ(iB)-iaLJ(iD)+1)/2
     &         +iaLL(iD))*dpC6J*dsqrt(3.d0)
     &         *dsqrt(6.d0)*dsqrt(2.d0)*dsqrt(iaLJ(iA)+1.0d0)
     &         *dsqrt(iaLJ(iB)+1.0d0)*dsqrt(iaLJ(iC)+1.0d0)
     &         *dsqrt(iaLJ(iD)+1.0d0)/4.d0
      dpVPHJ_F5=0.0d0
      call SUBPHME_SO_F4(iJ,iB,iA,iD,iC,dpVPHJ_F5)
      dpC6J=0.0d0
      call SIXJ(dpSLB,dpSLC,dpSJ,dpSJC,dpSJB,0.5d0,dpC6J)
      dpVPHJ_F5=dpVPHJ_F5*if_POT((-iaLJ(iB)+1)/2
     &         -iaLL(iC)+1)*dpC6J*dsqrt(3.d0)
     &         *dsqrt(6.d0)*dsqrt(2.d0)*dsqrt(iaLJ(iA)+1.0d0)
     &         *dsqrt(iaLJ(iB)+1.0d0)*dsqrt(iaLJ(iC)+1.0d0)
     &         *dsqrt(iaLJ(iD)+1.0d0)/4.d0
      dpVPHJ=dpVPHJ_F1+dpVPHJ_F2+dpVPHJ_F3+dpVPHJ_F4
     &      +dpVPHJ_F5+dpVPHJ_F6+dpVPHJ_F7+dpVPHJ_F8
      if((iTA.eq.iTD).and.(iTB.eq.iTC))then
       dpCIS=1.d0
       dpCIV=dfloat(iTA*iTB)
      else
       dpCIS=0.d0
       dpCIV=0.d0
      end if
      if(((iTA+iTD).eq.0).and.((iTB+iTC).eq.0))then
       dpCIV=dpCIV+1.d0+dfloat(iTA*iTC)
      end if
      dpSOPHME=0.5d0*(2.d0*dpW2+dpW2P)*dpCIS*dpVPHJ
     &        +0.5d0*dpW2P*dpCIV*dpVPHJ
      iExp=(-iaLL(iA)-iaLL(iB)+iaLL(iC)+iaLL(iD))/2
      dpSOPHME=if_POT(iExp)*dpSOPHME

      return
      end


C     ================================================================
      subroutine SUBPHME_SO_F1(iJF,iA,iB,iC,iD,dpVPHLA1)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)

      dpSJF=dfloat(iJF)
      dpSLA=dfloat(iaLL(iA))
      dpSLB=dfloat(iaLL(iB))
      dpSJA=dfloat(iaLJ(iA))*0.5d0
      dpSJB=dfloat(iaLJ(iB))*0.5d0
      dpSLC=dfloat(iaLL(iC))
      dpSLD=dfloat(iaLL(iD))
      dpSJC=dfloat(iaLJ(iC))*0.5d0
      dpSJD=dfloat(iaLJ(iD))*0.5d0
      dpVPHLA1=0.0d0
      do iK=-1,1,2
       do iI=-1,1,2
        dpSKA=dpSLA+dfloat(iK)
        dpSID=dpSLD+dfloat(iI)
        if(iaLL(iA)+iK.lt.0.or.iaLL(iD)+iI.lt.0)cycle
        dpVAlpha=0.d0
        dpSLALP=dfloat(iJF)
        dpY_AD=dpf_REDY(iaLL(iA)+iK,iJF,iaLL(iD)+iI)
        dpY_BC=dpf_REDY(iaLL(iB),iJF,iaLL(iC))
        iJTPrimeMin=max(iabs(iaLL(iD)-iaLL(iA)),iabs(iJF-1))
        iJTPrimeMax=min(iaLL(iD)+iaLL(iA),iJF+1)
        do iLAMBP=iJTPrimeMin,iJTPrimeMax
         dpSLAMBP=dfloat(iLAMBP)
         dpC9J2=0.0d0
         call NEUFJ(dpSLA,0.5d0,dpSJA,dpSLD,0.5d0,dpSJD,
     &              dpSLAMBP,1.0d0,dpSLALP,dpC9J2)
         dpC9J1=0.0d0
         call NEUFJ(dpSLALP,dpSLAMBP,1.0d0,dpSID,dpSLD,1.0d0,
     &              dpSKA,dpSLA,1.0d0,dpC9J1)
         dpVAlpha=dpVAlpha+if_POT(iLAMBP)*(2.d0*dpSLAMBP+1.d0)
     &           *dpC9J1*dpC9J2*dpY_AD*dpY_BC/(2.d0*dpSLALP+1.d0)
        endDO
        dpS_AD=0.d0
        call SUBINTRAD_F1(iK,iI,iA,iB,iC,iD,iaLL(iA),iaLL(iD),dpS_AD)
        dpVPHLA1=dpVPHLA1+if_POT((iK+iI)/2-iK)*dpS_AD*dpVAlpha
       end do
      end do

      return
      end


C     ================================================================
      subroutine SUBPHME_SO_F2(iJF,iA,iB,iC,iD,dpVPHLA2)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP)
     &      ,iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF)
     &      ,iaIPP(IP_NCF),iaINN(IP_NCF)

      dpSJF=dfloat(iJF)
      dpSLA=dfloat(iaLL(iA))
      dpSLB=dfloat(iaLL(iB))
      dpSJA=dfloat(iaLJ(iA))*0.5d0
      dpSJB=dfloat(iaLJ(iB))*0.5d0
      dpSLC=dfloat(iaLL(iC))
      dpSLD=dfloat(iaLL(iD))
      dpSJC=dfloat(iaLJ(iC))*0.5d0
      dpSJD=dfloat(iaLJ(iD))*0.5d0
      dpVPHLA2=0.0d0
      do iK=-1,1,2
       do iI=-1,1,2
        dpSKB=dpSLB+dfloat(iK)
        dpSID=dpSLD+dfloat(iI)
        if(iaLL(iB)+iK.lt.0.or.iaLL(iD)+iI.lt.0)cycle
        iJTPrimeMin=max(dabs(dpSLA-dpSLD-iI),dabs(dpSLB+iK-dpSLC),
     &                  dabs(dpSJF-1))
        iJTPrimeMax=min(dpSJF+1,dpSLA+dpSLD+iI,dpSLB+iK+dpSLC)
        dpVALP=0.d0
        do iLALP=iJTPrimeMin,iJTPrimeMax
         dpSLALP=dfloat(iLALP)
         dpY_AD=dpf_REDY(iaLL(iA),iLALP,iaLL(iD)+iI)
         dpY_BC=dpf_REDY(iaLL(iB)+iK,iLALP,iaLL(iC))
         dpC6J=0.0d0
         call SIXJ(dpSLALP,1.0d0,dpSJF,dpSLB,dpSLC,dpSKB,dpC6J)
         iJTPrimeMinD=max(dabs(dpSLALP-1),dabs(dpSLA-dpSLD),
     &                    dabs(dpSJF-1))
         iJTPrimeMaxD=min(dpSJF+1,dpSLA+dpSLD,dpSLALP+1)
         do iLAM=iJTPrimeMinD,iJTPrimeMaxD
          dpSLAM=dfloat(iLAM)
          dpC6J1=0.0d0
          call SIXJ(dpSJF,dpSLALP,1.0d0,1.0d0,1.0d0,dpSLAM,dpC6J1)
          dpC6J2=0.0d0
          call SIXJ(1.0d0,dpSLALP,dpSLAM,dpSLA,dpSLD,dpSID,dpC6J2)
          dpC9J1=0.0d0
          call NEUFJ(dpSLA,dpSLAM,dpSLD,dpSJA,dpSJF,dpSJD,0.5d0,
     &               1.0d0,0.5d0,dpC9J1)
          dpVALP=dpVALP+if_POT(-iLAM)*(2.d0*dpSLAM+1.d0)
     &          *dpC6J1*dpC6J2*dpC9J1*dpY_AD*dpY_BC*dpC6J
         end do
        end do
        dpS_BD=0.d0
        call SUBINTRAD_F1(iK,iI,iB,iA,iC,iD,iaLL(iB),iaLL(iD),dpS_BD)
        dpVPHLA2=dpVPHLA2+if_POT((iK+iI)/2-iK)*dpVALP*dpS_BD
       end do
      end do

      return
      end


C     ================================================================
      subroutine SUBPHME_SO_F3(iJF,iA,iB,iC,iD,dpVPHLA3)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)

      dpSJF=dfloat(iJF)
      dpSLA=dfloat(iaLL(iA))
      dpSLB=dfloat(iaLL(iB))
      dpSJA=dfloat(iaLJ(iA))*0.5d0
      dpSJB=dfloat(iaLJ(iB))*0.5d0
      dpSLC=dfloat(iaLL(iC))
      dpSLD=dfloat(iaLL(iD))
      dpSJC=dfloat(iaLJ(iC))*0.5d0
      dpSJD=dfloat(iaLJ(iD))*0.5d0
      dpVPHLA3=0.0d0
      do iK=-1,1,2
       do iI=-1,1,2
        dpSKA=dpSLA+dfloat(iK)
        dpSID=dpSLD+dfloat(iI)
        if(iaLL(iA)+iK.lt.0.or.iaLL(iD)+iI.lt.0)cycle
        iLAMBP=iJF
        dpSLAMBP=dfloat(iLAMBP)
        iJTPrimeMin=max(iabs(iaLL(iA)+iK-iaLL(iD)-iI),
     &                  iabs(iaLL(iB)-iaLL(iC)),iabs(iJF-1))
        iJTPrimeMax=min(iaLL(iA)+iK+iaLL(iD)+iI,iaLL(iB)+iaLL(iC),iJF+1)
        dpVALP=0.d0
        do iLALP=iJTPrimeMin,iJTPrimeMax
         dpSLALP=dfloat(iLALP)
         dpSLAMB=dpSLALP
         dpY_AD=dpf_REDY(iaLL(iA)+iK,iLALP,iaLL(iD)+iI)
         dpY_BC=dpf_REDY(iaLL(iB),iLALP,iaLL(iC))
         dpC9J2=0.0d0
         call NEUFJ(dpSLALP,dpSLAMBP,1.0d0,dpSID,dpSLD,1.0d0,
     &              dpSKA,dpSLA,1.0d0,dpC9J2)
         dpC9J1=0.0d0
         call NEUFJ(dpSLC,0.5d0,dpSJC,dpSLB,0.5d0,dpSJB,
     &              dpSLAMB,1.0d0,dpSJF,dpC9J1)
         dpVALP=dpVALP+dpY_AD*dpY_BC*dpC9J1*dpC9J2
        end do
        dpS_AD=0.d0
        call SUBINTRAD_F1(iK,iI,iA,iB,iC,iD,iaLL(iA),iaLL(iD),dpS_AD)
        dpVPHLA3=dpVPHLA3+if_POT((iK+iI)/2-iI)*dpS_AD*dpVALP
       end do
      end do

      return
      end


C     ================================================================
      subroutine SUBPHME_SO_F4(iJF,iA,iB,iC,iD,dpVPHLA4)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)

      dpSJF=dfloat(iJF)
      dpSLA=dfloat(iaLL(iA))
      dpSLB=dfloat(iaLL(iB))
      dpSJA=dfloat(iaLJ(iA))*0.5d0
      dpSJB=dfloat(iaLJ(iB))*0.5d0
      dpSLC=dfloat(iaLL(iC))
      dpSLD=dfloat(iaLL(iD))
      dpSJC=dfloat(iaLJ(iC))*0.5d0
      dpSJD=dfloat(iaLJ(iD))*0.5d0
      dpVPHLA4=0.0d0
      do iK=-1,1,2
       do iI=-1,1,2
        dpSKB=dpSLB+dfloat(iK)
        dpSID=dpSLD+dfloat(iI)
        if(iaLL(iB)+iK.lt.0.or.iaLL(iD)+iI.lt.0)cycle
        iJTPrimeMin=max(dabs(dpSLA-dpSLD-iI),dabs(dpSLB+iK-dpSLC),
     &                  dabs(dpSJF-1))
        iJTPrimeMax=min(dpSJF+1,dpSLA+dpSLD+iI,dpSLB+iK+dpSLC)
        dpVALP=0.d0
        do iLALP=iJTPrimeMin,iJTPrimeMax
         dpSLALP=dfloat(iLALP)
         dpY_AD=dpf_REDY(iaLL(iA),iLALP,iaLL(iD)+iI)
         dpY_BC=dpf_REDY(iaLL(iB)+iK,iLALP,iaLL(iC))
         dpC6J2=0.0d0
         call SIXJ(1.0d0,dpSLALP,dpSJF,dpSLA,dpSLD,dpSID,dpC6J2)
         iJTPrimeMin1=max(dabs(dpSLB-dpSLC),dabs(dpSLALP-1),
     &                    dabs(dpSJF-1))
         iJTPrimeMax1=min(dpSJF+1,dpSLB+dpSLC,dpSLALP+1)
         do iLAMB=iJTPrimeMin1,iJTPrimeMax1
          dpSLAMB=dfloat(iLAMB)
          dpC6J1=0.0d0
          call SIXJ(dpSLAMB,dpSLALP,1.0d0,1.0d0,1.0d0,dpSJF,dpC6J1)
          dpC6J3=0.0d0
          call SIXJ(1.0d0,dpSLALP,dpSLAMB,dpSLC,dpSLB,dpSKB,dpC6J3)
          dpC9J=0.0d0
          call NEUFJ(dpSLB,0.5d0,dpSJB,dpSLC,0.5d0,dpSJC,
     &               dpSLAMB,1.0d0,dpSJF,dpC9J)
          dpVALP=dpVALP+(2.d0*dpSLAMB+1.d0)*dpC6J1*dpC6J3*dpC9J
     &          *dpY_AD*dpY_BC*dpC6J2
         end do
        end do
        dpS_BD=0.d0
        call SUBINTRAD_F1(iK,iI,iB,iA,iC,iD,iaLL(iB),iaLL(iD),dpS_BD)
        dpVPHLA4=dpVPHLA4+if_POT((iK+iI)/2+iI)*dpVALP*dpS_BD
       end do
      end do

      return
      end


C     ================================================================
      function dpf_REDY(iL1,iL2,iL3)

      implicit integer(i)
      implicit double precision(d)

      dpREDY=0.d0
      iL13M=iabs(iL1-iL3)
      iL13P=iL1+iL3
      if((iL2-iL13M)*(iL2-iL13P).le.0)then
       iL123=iL13P+iL2
       iL123=mod(iL123,2)
       if(iL123.eq.0)then
        dpQP=4.d0*3.14159265d0
        dpXX4=(2.d0*iL1+1.d0)*(2.d0*iL2+1.d0)*(2.d0*iL3+1.d0)
        dpXX4=dpXX4/dpQP
        dpSG=1.d0-2.d0*mod(iL1,2)
        dpXX1=dfloat(iL1)
        dpXX2=dfloat(iL2)
        dpXX3=dfloat(iL3)
        dpXM=0.d0
        call TROISJ(dpXX1,dpXX2,dpXX3,dpXM,dpXM,dpXM,dpC3J)
        dpREDY=dpSG*dsqrt(dpXX4)*dpC3J
       end if
      end if
      dpf_REDY=dpREDY

      return
      end


C     ================================================================
C     Six-j symbols
C     ================================================================
      subroutine SIXJ(dpXJ1,dpXJ2,dpXJ3,dpXL1,dpXL2,dpXL3,dpC6J)

      implicit integer(i)
      implicit double precision(d)

      dimension dpaFLog(301)

      data(dpaFLog(i),i=2,31)/0.d0,.69314718d0,1.7917595d0,3.1780538d0,
     A4.7874917d0,6.5792511d0,8.5251613d0,10.604603d0,12.801827d0,
     B15.104413d0,17.502307d0,19.987214d0,22.552163d0,25.191221d0,
     C27.899271d0,30.671860d0,33.505072d0,36.395445d0,39.339884d0,
     D42.335616d0,45.380139d0,48.471180d0,51.606674d0,54.784729d0,
     E58.003604d0,61.261702d0,64.557537d0,67.889743d0,71.257038d0,
     F74.658235d0/
      data(dpaFLog(i),i=32,61)/78.092223d0,81.557959d0,85.054466d0,
     A88.580827d0,92.136175d0,95.719694d0,99.330612d0,102.96820d0,
     B106.63176d0,110.32064d0,114.03421d0,117.77188d0,121.53308d0,
     C125.31727d0,129.12393d0,132.95257d0,136.80272d0,140.67392d0,
     D144.56574d0,148.47776d0,152.40959d0,156.36083d0,160.33112d0,
     E164.32011d0,168.32744d0,172.35279d0,176.39584d0,180.45629d0,
     F184.53383d0,188.62817d0/
      data(dpaFLog(i),i=62,91)/192.73904d0,196.86618d0,201.00931d0,
     A205.16820d0,209.34258d0,213.53224d0,217.73693d0,221.95644d0,
     B226.19054d0,230.43904d0,234.70172d0,238.97839d0,243.26885d0,
     C247.57291d0,251.89040d0,256.22113d0,260.56494d0,264.92164d0,
     D269.29110d0,273.67312d0,278.06757d0,282.47429d0,286.89313d0,
     E291.32394d0,295.76659d0,300.22094d0,304.68685d0,309.16419d0,
     F313.65283d0,318.15264d0/
      data(dpaFLog(i),i=92,121)/322.66349d0,327.18529d0,331.71788d0,
     A336.26118d0,340.81505d0,345.37940d0,349.95411d0,354.53908d0,
     B359.13420d0,363.73937d0,368.35449d0,372.97946d0,377.61419d0,
     C382.25859d0,386.91255d0,391.57598d0,396.24881d0,400.93094d0,
     D405.62230d0,410.32277d0,415.03230d0,419.75080d0,424.47819d0,
     E429.21439d0,433.95932d0,438.71291d0,443.47508d0,448.24576d0,
     F453.02489d0,457.81238d0/
      data(dpaFLog(i),i=122,151)/462.60817d0,467.41220d0,472.22438d0,
     A477.04466d0,481.87298d0,486.70926d0,491.55345d0,496.40547d0,
     B501.26529d0,506.13282d0,511.00802d0,515.89082d0,520.78117d0,
     C525.67901d0,530.58428d0,535.49694d0,540.41692d0,545.34417d0,
     D550.27865d0,555.22029d0,560.16905d0,565.12488d0,570.08772d0,
     E575.05753d0,580.03427d0,585.01787d0,590.00830d0,595.00552d0,
     F600.00946d0,605.02010d0/
      data(dpaFLog(i),i=152,181)/610.03738d0,615.06126d0,620.09170d0,
     A625.12866d0,630.17208d0,635.22193d0,640.27818d0,645.34077d0,
     B650.40968d0,655.48486d0,660.56626d0,665.65385d0,670.74760d0,
     C675.84747d0,680.95341d0,686.06541d0,691.18340d0,696.30735d0,
     D701.43726d0,706.57306d0,711.71472d0,716.86221d0,722.01551d0,
     E727.17456d0,732.33934d0,737.50983d0,742.68598d0,747.86776d0,
     F753.05516d0,758.24811d0/
      data(dpaFLog(i),i=182,211)/763.44661d0,768.65061d0,773.86010d0,
     A779.07503d0,784.29539d0,789.52114d0,794.75224d0,799.98869d0,
     B805.23044d0,810.47747d0,815.72973d0,820.98722d0,826.24991d0,
     C831.51778d0,836.79078d0,842.06890d0,847.35209d0,852.64036d0,
     D857.93366d0,863.23199d0,868.53529d0,873.84356d0,879.15676d0,
     E884.47488d0,889.79789d0,895.12577d0,900.45848d0,905.79603d0,
     F911.13836d0,916.48547d0/
      data(dpaFLog(i),i=212,241)/921.83732d0,927.19391d0,932.55521d0,
     A937.92118d0,943.29181d0,948.66710d0,954.04699d0,959.43148d0,
     B964.82056d0,970.21419d0,975.61235d0,981.01503d0,986.42220d0,
     C991.83385d0,997.24995d0,1002.6705d0,1008.0954d0,1013.5248d0,
     D1018.9585d0,1024.3966d0,1029.8389d0,1035.2857d0,1040.7367d0,
     E1046.1920d0,1051.6516d0,1057.1155d0,1062.5836d0,1068.0558d0,
     F1073.5323d0,1079.0129d0/
      data(dpaFLog(i),i=242,271)/1084.4977d0,1089.9866d0,1095.4797d0,
     A1100.9768d0,1106.4781d0,1111.9834d0,1117.4928d0,1123.0063d0,
     B1128.5237d0,1134.0452d0,1139.5706d0,1145.1001d0,1150.6335d0,
     C1156.1708d0,1161.7120d0,1167.2573d0,1172.8063d0,1178.3593d0,
     D1183.9161d0,1189.4768d0,1195.0413d0,1200.6097d0,1206.1818d0,
     E1211.7577d0,1217.3375d0,1222.9209d0,1228.5082d0,1234.0992d0,
     F1239.6939d0,1245.2924d0/
      data(dpaFLog(i),i=272,301)/1250.8944d0,1256.5003d0,1262.1097d0,
     A1267.7228d0,1273.3396d0,1278.9600d0,1284.5840d0,1290.2117d0,
     B1295.8429d0,1301.4777d0,1307.1160d0,1312.7580d0,1318.4034d0,
     C1324.0524d0,1329.7048d0,1335.3609d0,1341.0203d0,1346.6833d0,
     D1352.3497d0,1358.0196d0,1363.6929d0,1369.3697d0,1375.0499d0,
     E1380.7334d0,1386.4204d0,1392.1107d0,1397.8045d0,1403.5016d0,
     F1409.2020d0,1414.9058d0/
      data dpEPS1,dpEPS2/0.1d0,-0.2d0/

      dpXN=-dpXJ1+dpXJ2+dpXJ3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN1=int(dpXN)
      dpXN=-dpXL1+dpXL2+dpXJ3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN2=int(dpXN)
      dpXN=-dpXL1+dpXJ2+dpXL3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN3=int(dpXN)
      dpXN=-dpXJ1+dpXL2+dpXL3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN4=int(dpXN)
      dpXN=dpXJ1-dpXJ2+dpXJ3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN5=int(dpXN)
      dpXN=dpXL1-dpXL2+dpXJ3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN6=int(dpXN)
      dpXN=dpXL1-dpXJ2+dpXL3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN7=int(dpXN)
      dpXN=dpXJ1-dpXL2+dpXL3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN8=int(dpXN)
      dpXN=dpXJ1+dpXJ2-dpXJ3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN9=int(dpXN)
      dpXN=dpXL1+dpXL2-dpXJ3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN10=int(dpXN)
      dpXN=dpXL1+dpXJ2-dpXL3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN11=int(dpXN)
      dpXN=dpXJ1+dpXL2-dpXL3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN12=int(dpXN)
      dpXN=-dpXJ1-dpXL1+dpXJ3+dpXL3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN13=int(dpXN)
      dpXN=-dpXJ2-dpXL2+dpXJ3+dpXL3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN14=int(dpXN)
      dpXN=dpXJ1+dpXL1+dpXJ2+dpXL2+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN15=int(dpXN)+1
      dpXN=dpXJ1+dpXJ2+dpXJ3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN16=int(dpXN)+1
      dpXN=dpXL1+dpXL2+dpXJ3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN17=int(dpXN)+1
      dpXN=dpXL1+dpXJ2+dpXL3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN18=int(dpXN)+1
      dpXN=dpXJ1+dpXL2+dpXL3+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN19=int(dpXN)+1
      if(iN9.lt.0.or.iN5.lt.0.or.iN1.lt.0.or.iN10.lt.0.or.iN6.lt.0
     &   .or.iN2.lt.0.or.iN11.lt.0.or.iN7.lt.0.or.iN3.lt.0
     &   .or.iN12.lt.0.or.iN8.lt.0.or.iN4.lt.0)then
       dpC6J=0.d0
       return
      end if
      iK=iN17+iN18+iN19-3
      if(iK.le.0)then
       dpC6J=1.d0
       return
      end if
      iK=0
      iL=-iN13
      if(iL.gt.iK)iK=iL
      iL=-iN14
      if(iL.gt.iK)iK=iL
      iL=iN9
      if(iN10.lt.iL)iL=iN10
      if(iN11.lt.iL)iL=iN11
      if(iN12.lt.iL)iL=iN12
      if(iN15.lt.iL)iL=iN15
      dpF=1.d0
      dpS=1.d0
      i=iK+1
      do while(i.le.iL)
       iM1=i-1
       iNN=(iN9-iM1)*(iN10-iM1)*(iN11-iM1)*(iN12-iM1)
       iNND=i*(iN13+i)*(iN14+i)*(iN15-iM1)
       dpF=-dpF*dfloat(iNN)/dfloat(iNND)
       dpS=dpS+dpF
       i=i+1
      end do
      dpC2N=dpaFLog(iN1+1)+dpaFLog(iN2+1)+dpaFLog(iN3+1)
     &     +dpaFLog(iN4+1)+dpaFLog(iN5+1)+dpaFLog(iN6+1)
     &     +dpaFLog(iN7+1)+dpaFLog(iN8+1)+dpaFLog(iN9+1)
     &     +dpaFLog(iN10+1)+dpaFLog(iN11+1)+dpaFLog(iN12+1)
      dpC2N=0.5d0*dpC2N
      dpC2D=dpaFLog(iN16+1)+dpaFLog(iN17+1)
     &     +dpaFLog(iN18+1)+dpaFLog(iN19+1)
      dpC2D=0.5d0*dpC2D
      iKM1=iK-1
      iKP1=iK+1
      dpC2N=dpC2N+dpaFLog(iN15-iKM1)
      dpC2D=dpC2D+dpaFLog(iKP1)+dpaFLog(iN13+iKP1)+dpaFLog(iN14+iKP1)
     &     +dpaFLog(iN9-iKM1)+dpaFLog(iN10-iKM1)+dpaFLog(iN11-iKM1)
     &     +dpaFLog(iN12-iKM1)
      dpF=dpC2D-dpC2N
      if(dpF.le.80.d0)then
       dpF=dpC2N/dpC2D
       if(.not.((dpF.lt.1.01d0).and.(dpF.gt.0.98d0)))then
        dpC6J=dpS*dexp(dpC2N-dpC2D)
        iL=iN15+iKM1
        iK=iL/2
        iK=2*iK
        if(iL.ne.iK)dpC6J=-dpC6J
        return
       end if
      end if
      if(dpS.lt.0.d0)then
       dpS=dlog(-dpS)
       dpC6J=-dexp(dpS+dpC2N-dpC2D)
      else if(dpS.eq.0)then
       dpC6J=0.d0
       return
      else
       dpS=dlog(dpS)
       dpC6J=dexp(dpS+dpC2N-dpC2D)
      end if
      iL=iN15+iKM1
      iK=iL/2
      iK=2*iK
      if(iL.ne.iK)dpC6J=-dpC6J

      return
      end


C     ================================================================
C     Nine-j symbols
C     ================================================================
      subroutine NEUFJ(dpXJ11,dpXJ12,dpXJ13,dpXJ21,dpXJ22,dpXJ23,
     &                 dpXJ31,dpXJ32,dpXJ33,dpC9J)

      implicit integer(i)
      implicit double precision(d)

      dimension dpaFLog(301)

      data(dpaFLog(i),i=2,31)/0.d0,.69314718d0,1.7917595d0,3.1780538d0,
     A4.7874917d0,6.5792511d0,8.5251613d0,10.604603d0,12.801827d0,
     B15.104413d0,17.502307d0,19.987214d0,22.552163d0,25.191221d0,
     C27.899271d0,30.671860d0,33.505072d0,36.395445d0,39.339884d0,
     D42.335616d0,45.380139d0,48.471180d0,51.606674d0,54.784729d0,
     E58.003604d0,61.261702d0,64.557537d0,67.889743d0,71.257038d0,
     F74.658235d0/
      data(dpaFLog(i),i=32,61)/78.092223d0,81.557959d0,85.054466d0,
     A88.580827d0,92.136175d0,95.719694d0,99.330612d0,102.96820d0,
     B106.63176d0,110.32064d0,114.03421d0,117.77188d0,121.53308d0,
     C125.31727d0,129.12393d0,132.95257d0,136.80272d0,140.67392d0,
     D144.56574d0,148.47776d0,152.40959d0,156.36083d0,160.33112d0,
     E164.32011d0,168.32744d0,172.35279d0,176.39584d0,180.45629d0,
     F184.53383d0,188.62817d0/
      data(dpaFLog(i),i=62,91)/192.73904d0,196.86618d0,201.00931d0,
     A205.16820d0,209.34258d0,213.53224d0,217.73693d0,221.95644d0,
     B226.19054d0,230.43904d0,234.70172d0,238.97839d0,243.26885d0,
     C247.57291d0,251.89040d0,256.22113d0,260.56494d0,264.92164d0,
     D269.29110d0,273.67312d0,278.06757d0,282.47429d0,286.89313d0,
     E291.32394d0,295.76659d0,300.22094d0,304.68685d0,309.16419d0,
     F313.65283d0,318.15264d0/
      data(dpaFLog(i),i=92,121)/322.66349d0,327.18529d0,331.71788d0,
     A336.26118d0,340.81505d0,345.37940d0,349.95411d0,354.53908d0,
     B359.13420d0,363.73937d0,368.35449d0,372.97946d0,377.61419d0,
     C382.25859d0,386.91255d0,391.57598d0,396.24881d0,400.93094d0,
     D405.62230d0,410.32277d0,415.03230d0,419.75080d0,424.47819d0,
     E429.21439d0,433.95932d0,438.71291d0,443.47508d0,448.24576d0,
     F453.02489d0,457.81238d0/
      data(dpaFLog(i),i=122,151)/462.60817d0,467.41220d0,472.22438d0,
     A477.04466d0,481.87298d0,486.70926d0,491.55345d0,496.40547d0,
     B501.26529d0,506.13282d0,511.00802d0,515.89082d0,520.78117d0,
     C525.67901d0,530.58428d0,535.49694d0,540.41692d0,545.34417d0,
     D550.27865d0,555.22029d0,560.16905d0,565.12488d0,570.08772d0,
     E575.05753d0,580.03427d0,585.01787d0,590.00830d0,595.00552d0,
     F600.00946d0,605.02010d0/
      data(dpaFLog(i),i=152,181)/610.03738d0,615.06126d0,620.09170d0,
     A625.12866d0,630.17208d0,635.22193d0,640.27818d0,645.34077d0,
     B650.40968d0,655.48486d0,660.56626d0,665.65385d0,670.74760d0,
     C675.84747d0,680.95341d0,686.06541d0,691.18340d0,696.30735d0,
     D701.43726d0,706.57306d0,711.71472d0,716.86221d0,722.01551d0,
     E727.17456d0,732.33934d0,737.50983d0,742.68598d0,747.86776d0,
     F753.05516d0,758.24811d0/
      data(dpaFLog(i),i=182,211)/763.44661d0,768.65061d0,773.86010d0,
     A779.07503d0,784.29539d0,789.52114d0,794.75224d0,799.98869d0,
     B805.23044d0,810.47747d0,815.72973d0,820.98722d0,826.24991d0,
     C831.51778d0,836.79078d0,842.06890d0,847.35209d0,852.64036d0,
     D857.93366d0,863.23199d0,868.53529d0,873.84356d0,879.15676d0,
     E884.47488d0,889.79789d0,895.12577d0,900.45848d0,905.79603d0,
     F911.13836d0,916.48547d0/
      data(dpaFLog(i),i=212,241)/921.83732d0,927.19391d0,932.55521d0,
     A937.92118d0,943.29181d0,948.66710d0,954.04699d0,959.43148d0,
     B964.82056d0,970.21419d0,975.61235d0,981.01503d0,986.42220d0,
     C991.83385d0,997.24995d0,1002.6705d0,1008.0954d0,1013.5248d0,
     D1018.9585d0,1024.3966d0,1029.8389d0,1035.2857d0,1040.7367d0,
     E1046.1920d0,1051.6516d0,1057.1155d0,1062.5836d0,1068.0558d0,
     F1073.5323d0,1079.0129d0/
      data(dpaFLog(i),i=242,271)/1084.4977d0,1089.9866d0,1095.4797d0,
     A1100.9768d0,1106.4781d0,1111.9834d0,1117.4928d0,1123.0063d0,
     B1128.5237d0,1134.0452d0,1139.5706d0,1145.1001d0,1150.6335d0,
     C1156.1708d0,1161.7120d0,1167.2573d0,1172.8063d0,1178.3593d0,
     D1183.9161d0,1189.4768d0,1195.0413d0,1200.6097d0,1206.1818d0,
     E1211.7577d0,1217.3375d0,1222.9209d0,1228.5082d0,1234.0992d0,
     F1239.6939d0,1245.2924d0/
      data(dpaFLog(i),i=272,301)/1250.8944d0,1256.5003d0,1262.1097d0,
     A1267.7228d0,1273.3396d0,1278.9600d0,1284.5840d0,1290.2117d0,
     B1295.8429d0,1301.4777d0,1307.1160d0,1312.7580d0,1318.4034d0,
     C1324.0524d0,1329.7048d0,1335.3609d0,1341.0203d0,1346.6833d0,
     D1352.3497d0,1358.0196d0,1363.6929d0,1369.3697d0,1375.0499d0,
     E1380.7334d0,1386.4204d0,1392.1107d0,1397.8045d0,1403.5016d0,
     F1409.2020d0,1414.9058d0/
      data dpEPS1,dpEPS2/0.1d0,-0.2d0/

      dpXN=-dpXJ11+dpXJ21+dpXJ31+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN1=int(dpXN)
      dpXN=-dpXJ32+dpXJ33+dpXJ31+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN2=int(dpXN)
      dpXN=dpXJ11-dpXJ21+dpXJ31+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN5=int(dpXN)
      dpXN=dpXJ32-dpXJ33+dpXJ31+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN6=int(dpXN)
      dpXN=dpXJ11+dpXJ21-dpXJ31+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN9=int(dpXN)
      dpXN=dpXJ32+dpXJ33-dpXJ31+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN10=int(dpXN)
      dpXN=dpXJ11+dpXJ32+dpXJ21+dpXJ33+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN15=int(dpXN)+1
      dpXN=dpXJ11+dpXJ21+dpXJ31+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN16=int(dpXN)+1
      dpXN=dpXJ32+dpXJ33+dpXJ31+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN17=int(dpXN)+1
      dpXN=-dpXJ12+dpXJ22+dpXJ32+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN21=int(dpXN)
      dpXN=-dpXJ21+dpXJ22+dpXJ23+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN23=int(dpXN)
      dpXN=dpXJ12-dpXJ22+dpXJ32+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN25=int(dpXN)
      dpXN=dpXJ21-dpXJ22+dpXJ23+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN27=int(dpXN)
      dpXN=dpXJ12+dpXJ22-dpXJ32+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN29=int(dpXN)
      dpXN=dpXJ21+dpXJ22-dpXJ23+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN31=int(dpXN)
      dpXN=-dpXJ12-dpXJ21+dpXJ32+dpXJ23+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN33=int(dpXN)
      dpXN=dpXJ12+dpXJ22+dpXJ32+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN36=int(dpXN)+1
      dpXN=dpXJ21+dpXJ22+dpXJ23+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN38=int(dpXN)+1
      dpXN=-dpXJ13+dpXJ23+dpXJ33+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN41=int(dpXN)
      dpXN=-dpXJ13+dpXJ11+dpXJ12+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN44=int(dpXN)
      dpXN=dpXJ13-dpXJ23+dpXJ33+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN45=int(dpXN)
      dpXN=dpXJ13-dpXJ11+dpXJ12+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN48=int(dpXN)
      dpXN=dpXJ13+dpXJ23-dpXJ33+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN49=int(dpXN)
      dpXN=dpXJ13+dpXJ11-dpXJ12+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN52=int(dpXN)
      dpXN=-dpXJ23-dpXJ11+dpXJ33+dpXJ12+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN54=int(dpXN)
      dpXN=dpXJ13+dpXJ23+dpXJ33+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN56=int(dpXN)+1
      dpXN=dpXJ13+dpXJ11+dpXJ12+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN59=int(dpXN)+1
      if(iN9.lt.0.or.iN5.lt.0.or.iN1.lt.0.or.iN10.lt.0
     &   .or.iN6.lt.0.or.iN2.lt.0.or.iN29.lt.0.or.iN25.lt.0
     &   .or.iN21.lt.0.or.iN31.lt.0.or.iN27.lt.0.or.iN23.lt.0
     &   .or.iN49.lt.0.or.iN45.lt.0.or.iN41.lt.0.or.iN52.lt.0
     &   .or.iN48.lt.0.or.iN44.lt.0)then
       dpC9J=0.d0
       return
      end if
      iK=iN1+iN2+iN5+iN6+iN9+iN10+iN21+iN23+iN25+iN27+iN29+iN31
     &  +iN41+iN44+iN45+iN48+iN49+iN52
      if(iK.le.0)then
       dpC9J=1.d0
       return
      end if
      dpXN=2.d0*(dpXJ21-dpXJ32)+dpEPS1
      if(dpXN.lt.0.d0) dpXN=-(dpXN+dpEPS2)
      iJMIN=int(dpXN)
      dpXN=2.d0*(dpXJ11-dpXJ33)+dpEPS1
      if(dpXN.lt.0.d0)dpXN=-(dpXN+dpEPS2)
      iN=int(dpXN)
      if(iN.gt.iJMIN)iJMIN=iN
      dpXN=2.d0*(dpXJ12-dpXJ23)+dpEPS1
      if(dpXN.lt.0.d0)dpXN=-(dpXN+dpEPS2)
      iN=int(dpXN)
      if(iN.gt.iJMIN)iJMIN=iN
      dpXN=2.d0*(dpXJ21+dpXJ32)+dpEPS1
      iJMAX=int(dpXN)
      dpXN=2.d0*(dpXJ11+dpXJ33)+dpEPS1
      iN=int(dpXN)
      if(iN.lt.iJMAX)iJMAX=iN
      dpXN=2.d0*(dpXJ12+dpXJ23)+dpEPS1
      iN=int(dpXN)
      if(iN.lt.iJMAX)iJMAX=iN
      dpXJMIN=dfloat(iJMIN)/2.d0
      dpXJMAX=dfloat(iJMAX)/2.d0
      dpS=0.d0
      dpXJ=dpXJMIN
      dpXN=-dpXJ32+dpXJ21+dpXJ+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN3=int(dpXN)
      dpXN=-dpXJ11+dpXJ33+dpXJ+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN4=int(dpXN)
      dpXN=dpXJ32-dpXJ21+dpXJ+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN7=int(dpXN)
      dpXN=dpXJ11-dpXJ33+dpXJ+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN8=int(dpXN)
      dpXN=dpXJ32+dpXJ21-dpXJ+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN11=int(dpXN)
      dpXN=dpXJ11+dpXJ33-dpXJ+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN12=int(dpXN)
      dpXN=-dpXJ11-dpXJ32+dpXJ31+dpXJ+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN13=int(dpXN)
      dpXN=-dpXJ21-dpXJ33+dpXJ31+dpXJ+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN14=int(dpXN)
      dpXN=dpXJ32+dpXJ21+dpXJ+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN18=int(dpXN)+1
      dpXN=dpXJ11+dpXJ33+dpXJ+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN19=int(dpXN)+1
      dpXN=-dpXJ21+dpXJ+dpXJ32+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN22=int(dpXN)
      dpXN=-dpXJ12+dpXJ23+dpXJ+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN24=int(dpXN)
      dpXN=dpXJ21-dpXJ+dpXJ32+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN26=int(dpXN)
      dpXN=dpXJ12-dpXJ+dpXJ23+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN28=int(dpXN)
      dpXN=dpXJ21+dpXJ-dpXJ32+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN30=int(dpXN)
      dpXN=dpXJ12+dpXJ-dpXJ23+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN32=int(dpXN)
      dpXN=-dpXJ22-dpXJ+dpXJ32+dpXJ23+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN34=int(dpXN)
      dpXN=dpXJ12+dpXJ21+dpXJ22+dpXJ+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN35=int(dpXN)+1
      dpXN=dpXJ21+dpXJ+dpXJ32+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN37=int(dpXN)+1
      dpXN=dpXJ12+dpXJ+dpXJ23+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN39=int(dpXN)+1
      dpXN=-dpXJ+dpXJ11+dpXJ33+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN42=int(dpXN)
      dpXN=-dpXJ+dpXJ23+dpXJ12+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN43=int(dpXN)
      dpXN=dpXJ-dpXJ11+dpXJ33+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN46=int(dpXN)
      dpXN=dpXJ-dpXJ23+dpXJ12+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN47=int(dpXN)
      dpXN=dpXJ+dpXJ11-dpXJ33+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN50=int(dpXN)
      dpXN=dpXJ+dpXJ23-dpXJ12+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN51=int(dpXN)
      dpXN=-dpXJ13-dpXJ+dpXJ33+dpXJ12+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN53=int(dpXN)
      dpXN=dpXJ13+dpXJ+dpXJ23+dpXJ11+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN55=int(dpXN)+1
      dpXN=dpXJ+dpXJ11+dpXJ33+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN57=int(dpXN)+1
      dpXN=dpXJ+dpXJ23+dpXJ12+dpEPS1
      if(dpXN.lt.0.d0)dpXN=dpXN+dpEPS2
      iN58=int(dpXN)+1
      iStep=0
      do while(.true.)
       if(iStep.ne.0)then
        if(dpXJ.gt.dpXJMAX)exit
        iN3=iN3+1
        iN4=iN4+1
        iN7=iN7+1
        iN8=iN8+1
        iN11=iN11-1
        iN12=iN12-1
        iN13=iN13+1
        iN14=iN14+1
        iN18=iN18+1
        iN19=iN19+1
        iN22=iN22+1
        iN24=iN24+1
        iN26=iN26-1
        iN28=iN28-1
        iN30=iN30+1
        iN32=iN32+1
        iN34=iN34-1
        iN35=iN35+1
        iN37=iN37+1
        iN39=iN39+1
        iN42=iN42-1
        iN43=iN43-1
        iN46=iN46+1
        iN47=iN47+1
        iN50=iN50+1
        iN51=iN51+1
        iN53=iN53-1
        iN55=iN55+1
        iN57=iN57+1
        iN58=iN58+1
       end if
       iStep=1
       iK1=0.d0
       iL1=-iN13
       if(iL1.gt.iK1)iK1=iL1
       iL1=-iN14
       if(iL1.gt.iK1)iK1=iL1
       iL1=iN9
       if(iN10.lt.iL1)iL1=iN10
       if(iN11.lt.iL1)iL1=iN11
       if(iN12.lt.iL1)iL1=iN12
       if(iN15.lt.iL1)iL1=iN15
       dpF1=1.d0
       dpS1=1.d0
       i1=iK1+1
       do while(i1.le.iL1)
        i1M1=i1-1
        iNN1=(iN9-i1M1)*(iN10-i1M1)*(iN11-i1M1)*(iN12-i1M1)
        iND1=i1*(iN13+i1)*(iN14+i1)*(iN15-i1M1)
        dpF1=-dpF1*dfloat(iNN1)/dfloat(iND1)
        dpS1=dpS1+dpF1
        i1=i1+1
       end do
       iK2=0
       iL2=-iN33
       if(iL2.gt.iK2)iK2=iL2
       iL2=-iN34
       if(iL2.gt.iK2)iK2=iL2
       iL2=iN29
       if(iN30.lt.iL2)iL2=iN30
       if(iN31.lt.iL2)iL2=iN31
       if(iN32.lt.iL2)iL2=iN32
       if(iN35.lt.iL2)iL2=iN35
       dpF2=1.d0
       dpS2=1.d0
       i2=iK2+1
       do while(i2.le.iL2)
        i2M2=i2-1
        iNN2=(iN29-i2M2)*(iN30-i2M2)*(iN31-i2M2)*(iN32-i2M2)
        iND2=i2*(iN33+i2)*(iN34+i2)*(iN35-i2M2)
        dpF2=-dpF2*dfloat(iNN2)/dfloat(iND2)
        dpS2=dpS2+dpF2
        i2=i2+1
       end do
       iK3=0
       iL3=-iN53
       if(iL3.gt.iK3)iK3=iL3
       iL3=-iN54
       if(iL3.gt.iK3)iK3=iL3
       iL3=iN49
       if(iN50.lt.iL3)iL3=iN50
       if(iN51.lt.iL3)iL3=iN51
       if(iN52.lt.iL3)iL3=iN52
       if(iN55.lt.iL3)iL3=iN55
       dpF3=1.d0
       dpS3=1.d0
       i3=iK3+1
       do while(i3.le.iL3)
        i3M3=i3-1
        iNN3=(iN49-i3M3)*(iN50-i3M3)*(iN51-i3M3)*(iN52-i3M3)
        iND3=i3*(iN53+i3)*(iN54+i3)*(iN55-i3M3)
        dpF3=-dpF3*dfloat(iNN3)/dfloat(iND3)
        dpS3=dpS3+dpF3
        i3=i3+1
       end do
       dpS2N=dpaFLog(iN3+1)+dpaFLog(iN4+1)+dpaFLog(iN7+1)
     &      +dpaFLog(iN8+1)+dpaFLog(iN11+1)+dpaFLog(iN12+1)
     &      +dpaFLog(iN22+1)+dpaFLog(iN24+1)+dpaFLog(iN26+1)
     &      +dpaFLog(iN28+1)+dpaFLog(iN30+1)+dpaFLog(iN32+1)
     &      +dpaFLog(iN42+1)+dpaFLog(iN43+1)+dpaFLog(iN46+1)
     &      +dpaFLog(iN47+1)+dpaFLog(iN50+1)+dpaFLog(iN51+1)
       dpS2N=0.5d0*dpS2N
       dpS2D=dpaFLog(iN18+1)+dpaFLog(iN19+1)+dpaFLog(iN37+1)
     &      +dpaFLog(iN39+1)+dpaFLog(iN57+1)+dpaFLog(iN58+1)
       dpS2D=0.5d0*dpS2D
       iKM1=iK1-1
       iKP1=iK1+1
       iKM2=iK2-1
       iKP2=iK2+1
       iKM3=iK3-1
       iKP3=iK3+1
       dpS2N=dpS2N+dpaFLog(iN15-iKM1)+dpaFLog(iN35-iKM2)
     &      +dpaFLog(iN55-iKM3)
       dpS2D=dpS2D+dpaFLog(iKP1)+dpaFLog(iKP2)+dpaFLog(iKP3)
     &      +dpaFLog(iN9-iKM1)+dpaFLog(iN10-iKM1)+dpaFLog(iN11-iKM1)
     &      +dpaFLog(iN12-iKM1)+dpaFLog(iN13+iKP1)+dpaFLog(iN14+iKP1)
     &      +dpaFLog(iN29-iKM2)+dpaFLog(iN30-iKM2)+dpaFLog(iN31-iKM2)
     &      +dpaFLog(iN32-iKM2)+dpaFLog(iN33+iKP2)+dpaFLog(iN34+iKP2)
     &      +dpaFLog(iN49-iKM3)+dpaFLog(iN50-iKM3)+dpaFLog(iN51-iKM3)
     &      +dpaFLog(iN52-iKM3)+dpaFLog(iN53+iKP3)+dpaFLog(iN54+iKP3)
       dpF=dpS2D-dpS2N
       if(dpF.le.80.d0)then
        dpF=dpS2N/dpS2D
        if(.not.((dpF.lt.1.01d0).and.(dpF.gt.0.98d0)))then
         dpF=dpS1*dpS2*dpS3*dexp(dpS2N-dpS2D)
         dpF=dpF*(2.d0*dpXJ+1.d0)
         iL1=iK1+iK2+iK3
         iK1=iL1/2
         iK1=2*iK1
         if(iL1.ne.iK1)dpF=-dpF
         dpS=dpS+dpF
         dpXJ=dpXJ+1.d0
         cycle
        end if
       end if
       dpF=dpS1*dpS2*dpS3
       if(dpF.ne.0.d0)then
        if(dpF.lt.0.d0)then
         dpF=dlog(-dpF)
         dpF=-dexp(dpF+dpS2N-dpS2D)
        else if(dpF.gt.0.d0)then
         dpF=dlog(dpF)
         dpF=dexp(dpF+dpS2N-dpS2D)
        end if
        dpF=dpF*(2.d0*dpXJ+1.d0)
        iL1=iK1+iK2+iK3
        iK1=iL1/2
        iK1=2*iK1
        if(iL1.ne.iK1)dpF=-dpF
        dpS=dpS+dpF
       end if
       dpXJ=dpXJ+1.d0
      end do
      dpC2N=dpaFLog(iN1+1)+dpaFLog(iN2+1)+dpaFLog(iN5+1)
     &     +dpaFLog(iN6+1)+dpaFLog(iN9+1)+dpaFLog(iN10+1)
     &     +dpaFLog(iN21+1)+dpaFLog(iN23+1)+dpaFLog(iN25+1)
     &     +dpaFLog(iN27+1)+dpaFLog(iN29+1)+dpaFLog(iN31+1)
     &     +dpaFLog(iN41+1)+dpaFLog(iN44+1)+dpaFLog(iN45+1)
     &     +dpaFLog(iN48+1)+dpaFLog(iN49+1)+dpaFLog(iN52+1)
      dpC2N=0.5d0*dpC2N
      dpC2D=dpaFLog(iN16+1)+dpaFLog(iN17+1)+dpaFLog(iN36+1)
     &     +dpaFLog(iN38+1)+dpaFLog(iN56+1)+dpaFLog(iN59+1)
      dpC2D=0.5d0*dpC2D
      dpF=dpC2D-dpC2N
      if(dpF.le.80.d0)then
       dpF=dpC2N/dpC2D
       if(.not.((dpF.lt.1.01d0).and.(dpF.gt.0.98d0)))then
        dpC9J=dpS*dexp(dpC2N-dpC2D)
        iK=iN9+iN16+iN36+iN56-1
        iL=iK/2
        iL=2*iL
        if(iL.ne.iK)dpC9J=-dpC9J
        return
       end if
      end if
      if(dpS.lt.0.d0)then
       dpS=dlog(-dpS)
       dpC9J=-dexp(dpS+dpC2N-dpC2D)
      else if(dpS.eq.0.d0)then
       dpC9J=0.d0
       return
      else
       dpS=dlog(dpS)
       dpC9J=dexp(dpS+dpC2N-dpC2D)
      end if
      iK=iN9+iN16+iN36+iN56-1
      iL=iK/2
      iL=2*iL
      if(iL.ne.iK)dpC9J=-dpC9J

      return
      end


C     ================================================================
      subroutine SUBINTRAD_F1(i,iJ,i1,i2,i3,i4,iL1,iL4,dpSU)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_WF/dpaWF(IP_NNP,IP_NSP),dpaDWF(IP_NNP,IP_NSP)
      common/C_HF/iNMAX,iNOcc,iNUnocc,iNOrb
      common/C_MESHR/dpStep

      do iR=1,iNMAX
       dpR_i1=0.0d0
       dpR_i4=0.0d0
       call SUBOPRAD(iR,i,i1,iL1,dpR_i1)
       call SUBOPRAD(iR,iJ,i4,iL4,dpR_i4)
       dpSU=dpSU+dpR_i1*dpR_i4*dpaWF(iR,i2)*dpaWF(iR,i3)*dpStep
      end do

      return
      end


C     ================================================================
      subroutine SUBOPRAD(iR,iK,i,iL,dpR)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_WF/dpaWF(IP_NNP,IP_NSP),dpaDWF(IP_NNP,IP_NSP)
      common/C_MESHR/dpStep

      dpRadP=dpStep*iR
      dpR=(dpaDWF(iR,i)/dpRadP-1/(dpRadP**2)*dpaWF(iR,i)
     &   +if_POT((iK+1)/2)*(iL+(1-iK)/2)*dpaWF(iR,i)/(dpRadP**2))
     &   *dsqrt(iL+(iK+1.0d0)/2.0d0)

      return
      end


C     ================================================================
C     Set up A and B RPA matrices
C     ================================================================
      subroutine MATRIX(dpaA1,dpaB1)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_DIM/iNV,iNV1,iProt,iNeutr
      common/C_TODO/iRPA
      common/C_CTR/dpQBess,iISPIN,iORB,iPAR,iISFLIP,iRAD
      common/C_AB/iaNUM(IP_NCF)
      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)
      common/C_MATEL/dpaPHME(IP_NCF,IP_NCF),dpaPHMEEX(IP_NCF,IP_NCF)
      common/C_GFB/dpaECF(IP_NCF),dpaAKeep(IP_NCF)
      common/C_ENERG/dpaSPEP(IP_NSP),dpaSPEN(IP_NSP)

      dimension dpaA(IP_NCF,IP_NCF),dpaB(IP_NCF,IP_NCF)
      dimension dpaA1(IP_NCF,IP_NCF),dpaB1(IP_NCF,IP_NCF)

      do i1=1,iNV
       iPRTY1=(-1)**(iaLL(iaIPP(i1))+iaLL(iaINN(i1)))
       iQ1=iaIQ(iaIPP(i1))+iaIQ(iaINN(i1))
       do i2=i1,iNV
        iPRTY2=(-1)**(iaLL(iaIPP(i2))+iaLL(iaINN(i2)))
        dpaA(i1,i2)=0.d0
        dpaB(i1,i2)=0.d0
        if(iaJJ(i1).ne.iaJJ(i2))cycle
        if(iPRTY1.ne.iPRTY2)cycle
        dpEConf=0.d0
        if(i1.eq.i2)then
         if(iQ1.eq.2)dpEConf=dpaSPEP(iaIPP(i1))-dpaSPEP(iaINN(i1))
         if(iQ1.eq.0)dpEConf=dpaSPEN(iaIPP(i1))-dpaSPEN(iaINN(i1))
        end if
        dpSJP2=dfloat(iaLJ(iaIPP(i2)))*0.5d0
        dpSJN2=dfloat(iaLJ(iaINN(i2)))*0.5d0
        iJPMJN=int(dpSJP2-dpSJN2)
        dpaA(i1,i2)=dpEConf+dpaPHME(i1,i2)
        dpaB(i1,i2)=(-1)**(iaJJ(i1)+iJPMJN)*dpaPHMEEX(i1,i2)
        if(abs(iRPA).eq.1)dpaB(i1,i2)=0.d0
        dpaA(i2,i1)=dpaA(i1,i2)
        dpaB(i2,i1)=dpaB(i1,i2)
        if(i1.eq.i2)then
         dpaECF(i1)=dpEConf
         dpaAKeep(i1)=dpaA(i1,i1)-dpEConf
        end if
       end do
      end do
      i3=1
      i4=1
      do i1=1,iNV
       iPRTY1=(-1)**(iaLL(iaIPP(i1))+iaLL(iaINN(i1)))
       if(iaJJ(i1).ne.iISPIN)cycle
       iaNUM(i3)=i1
       do i2=1,iNV
        if(iaJJ(i2).ne.iISPIN)cycle
        dpaA1(i3,i4)=dpaA(i1,i2)
        dpaB1(i3,i4)=dpaB(i1,i2)
        i4=i4+1
       end do
       i4=1
       i3=i3+1
      end do
      iNV1=i3-1

      return
      end


C     ================================================================
C     Diagonalization of RPA matrices, see Ring & Schuck, Pag. 306
C     ================================================================
      subroutine SRPA(dpaA,dpaB,dpaD,iFlag)

      implicit integer(i)
      implicit double precision(d)

      include 'param.rpa'

      common/C_DIM/iNV,iNV1,iProt,iNeutr

      dimension dpaA(IP_NCF,IP_NCF),dpaB(IP_NCF,IP_NCF),dpaD(IP_NCF)
      dimension dpaTemp(IP_NCF,IP_NCF)
      dimension dpaE(IP_NCF),dpaC(IP_NCF,IP_NCF)
      dimension dpaXPY(IP_NCF,IP_NCF),dpaXMY(IP_NCF,IP_NCF)
      dimension dpaTemp2(IP_NCF,IP_NCF)

      iFlag=0
      do i=1,iNV1
       do iJ=1,iNV1
        dpaTemp(i,iJ)=dpaA(i,iJ)+dpaB(i,iJ)
       end do
      end do
      call TRED2(dpaTemp,iNV1,IP_NCF,dpaD,dpaE)
      call TQLI(dpaD,dpaE,iNV1,IP_NCF,dpaTemp)
      do i=1,iNV1
       if(dpaD(i).lt.0.d0)then
        write(2,*)' WARNING : NO SQUARE ROOT OF A+B !!!',dpaD(i)
        iFlag=1
        return
       end if
      end do
      iFlag=0
      do i=1,iNV1
       do iJ=1,iNV1
        dpaC(i,iJ)=0.d0
        dpaXPY(i,iJ)=0.d0
        dpaXMY(i,iJ)=0.d0
        dpaTemp(i,iJ)=dpaA(i,iJ)-dpaB(i,iJ)
       end do
      end do
      call TRED2(dpaTemp,iNV1,IP_NCF,dpaD,dpaE)
      call TQLI(dpaD,dpaE,iNV1,IP_NCF,dpaTemp)
      do i=1,iNV1
       if(dpaD(i).lt.0.d0)then
        write(2,*)' WARNING : NO SQUARE ROOT OF A-B !!!'
        iFlag=1
        return
       end if
       dpaD(i)=dsqrt(dpaD(i))
      end do
      do i=1,iNV1
       do iJ=1,iNV1
        do iK=1,iNV1
         dpaC(i,iJ)=dpaC(i,iJ)
     &             +dpaTemp(i,iK)*dpaD(iK)*dpaTemp(iJ,iK)
        end do
       end do
      end do
      do i=1,iNV1
       do iJ=1,iNV1
        dpaTemp2(i,iJ)=0.d0
        do iK=1,iNV1
         dpaTemp2(i,iJ)=dpaTemp2(i,iJ)
     &                 +dpaC(i,iK)*(dpaA(iK,iJ)+dpaB(iK,iJ))
        end do
       end do
      end do
      do i=1,iNV1
       do iJ=1,iNV1
        dpaTemp(i,iJ)=0.d0
        do iK=1,iNV1
         dpaTemp(i,iJ)=dpaTemp(i,iJ)
     &                +dpaTemp2(i,iK)*dpaC(iK,iJ)
        end do
       end do
      end do
      call TRED2(dpaTemp,iNV1,IP_NCF,dpaD,dpaE)
      call TQLI(dpaD,dpaE,iNV1,IP_NCF,dpaTemp)
      do i=1,iNV1
       dpSgn=dpaD(i)/dabs(dpaD(i))
       dpaD(i)=dpSgn*dsqrt(dabs(dpaD(i)))
      end do
      do i=1,iNV1
       do iJ=1,iNV1
        dpaXPY(i,iJ)=0.d0
        do iK=1,iNV1
         dpaXPY(i,iJ)=dpaXPY(i,iJ)
     &               +dpaC(i,iK)*dpaTemp(iK,iJ)/dsqrt(dabs(dpaD(iJ)))
        end do
       end do
      end do
      do i=1,iNV1
       do iJ=1,iNV1
        dpaXMY(i,iJ)=0.d0
        do iK=1,iNV1
         dpaXMY(i,iJ)=dpaXMY(i,iJ)+(dpaA(i,iK)
     &               +dpaB(i,iK))*dpaXPY(iK,iJ)/dpaD(iJ)
        end do
       end do
      end do
      do i=1,iNV1
       do iJ=1,iNV1
        dpaTemp(i,iJ)=0.5d0*(dpaXPY(i,iJ)+dpaXMY(i,iJ))
        dpaTemp2(i,iJ)=0.5d0*(dpaXPY(i,iJ)-dpaXMY(i,iJ))
       end do
      end do
      do i=1,iNV1
       do iJ=1,iNV1
        dpaA(i,iJ)=dpaTemp(i,iJ)
        dpaB(i,iJ)=dpaTemp2(i,iJ)
       end do
      end do

      return
      end


C     ================================================================
      subroutine TRED2(dpaA,iNV,iNP,dpaD,dpaE)

      implicit integer(i)
      implicit double precision(d)

      dimension dpaA(iNP,iNP),dpaD(iNP),dpaE(iNP)

      if(iNV.gt.1)then
       do i=iNV,2,-1
        iL=i-1
        dpH=0.d0
        dpScale=0.d0
        if(iL.gt.1)then
         do iK=1,iL
          dpScale=dpScale+dabs(dpaA(i,iK))
         end do
         if(dpScale.eq.0.d0)then
          dpaE(i)=dpaA(i,iL)
         else
          do iK=1,iL
           dpaA(i,iK)=dpaA(i,iK)/dpScale
           dpH=dpH+dpaA(i,iK)**2
          end do
          dpF=dpaA(i,iL)
          dpG=-dsign(dsqrt(dpH),dpF)
          dpaE(i)=dpScale*dpG
          dpH=dpH-dpF*dpG
          dpaA(i,iL)=dpF-dpG
          dpF=0.d0
          do iJ=1,iL
           dpaA(iJ,i)=dpaA(i,iJ)/dpH
           dpG=0.d0
           do iK=1,iJ
            dpG=dpG+dpaA(iJ,iK)*dpaA(i,iK)
           end do
           if(iL.gt.iJ)then
            do iK=iJ+1,iL
             dpG=dpG+dpaA(iK,iJ)*dpaA(i,iK)
            end do
           end if
           dpaE(iJ)=dpG/dpH
           dpF=dpF+dpaE(iJ)*dpaA(i,iJ)
          end do
          dpHH=dpF/(dpH+dpH)
          do iJ=1,iL
           dpF=dpaA(i,iJ)
           dpG=dpaE(iJ)-dpHH*dpF
           dpaE(iJ)=dpG
           do iK=1,iJ
            dpaA(iJ,iK)=dpaA(iJ,iK)-dpF*dpaE(iK)-dpG*dpaA(i,iK)
           end do
          end do
         end if
        else
         dpaE(i)=dpaA(i,iL)
        end if
        dpaD(i)=dpH
       end do
      end if
      dpaD(1)=0.d0
      dpaE(1)=0.d0
      do i=1,iNV
       iL=i-1
       if(dpaD(i).ne.0.d0)then
        do iJ=1,iL
         dpG=0.d0
         do iK=1,iL
          dpG=dpG+dpaA(i,iK)*dpaA(iK,iJ)
         end do
         do iK=1,iL
          dpaA(iK,iJ)=dpaA(iK,iJ)-dpG*dpaA(iK,i)
         end do
        end do
       end if
       dpaD(i)=dpaA(i,i)
       dpaA(i,i)=1.d0
       if(iL.ge.1)then
        do iJ=1,iL
         dpaA(i,iJ)=0.d0
         dpaA(iJ,i)=0.d0
        end do
       end if
      end do

      return
      end


C     ================================================================
      subroutine TQLI(dpaD,dpaE,iNV,iNP,dpaZ)

      implicit integer(i)
      implicit double precision(d)
      implicit logical(b)

      dimension dpaD(iNP),dpaE(iNP),dpaZ(iNP,iNP)

      if(iNV.gt.1)then
       do i=2,iNV
        dpaE(i-1)=dpaE(i)
       end do
       dpaE(iNV)=0.d0
       do iL=1,iNV
        iIter=0
        do while(.true.)
         bMax=.true.
         do iM=iL,iNV-1
          dpDD=dabs(dpaD(iM))+dabs(dpaD(iM+1))
          if(dabs(dpaE(iM))+dpDD.eq.dpDD)then
           bMax=.false.
           exit
          end if
         end do
         if(bMax.eqv..true.)iM=iNV
         if(iM.ne.iL)then
          if(iIter.eq.30)stop 'TOO MANY ITERATIONS'
          iIter=iIter+1
          dpG=(dpaD(iL+1)-dpaD(iL))/(2.d0*dpaE(iL))
          dpR=dsqrt(dpG**2+1.d0)
          dpG=dpaD(iM)-dpaD(iL)+dpaE(iL)/(dpG+sign(dpR,dpG))
          dpS=1.d0
          dpC=1.d0
          dpP=0.d0
          do i=iM-1,iL,-1
           dpF=dpS*dpaE(i)
           dpB=dpC*dpaE(i)
           if(dabs(dpF).ge.dabs(dpG))then
            dpC=dpG/dpF
            dpR=dsqrt(dpC**2+1.d0)
            dpaE(i+1)=dpF*dpR
            dpS=1.d0/dpR
            dpC=dpC*dpS
           else
            dpS=dpF/dpG
            dpR=dsqrt(dpS**2+1.d0)
            dpaE(i+1)=dpG*dpR
            dpC=1.d0/dpR
            dpS=dpS*dpC
           end if
           dpG=dpaD(i+1)-dpP
           dpR=(dpaD(i)-dpG)*dpS+2.d0*dpC*dpB
           dpP=dpS*dpR
           dpaD(i+1)=dpG+dpP
           dpG=dpC*dpR-dpB
           do iK=1,iNV
            dpF=dpaZ(iK,i+1)
            dpaZ(iK,i+1)=dpS*dpaZ(iK,i)+dpC*dpF
            dpaZ(iK,i)=dpC*dpaZ(iK,i)-dpS*dpF
           end do
          end do
          dpaD(iL)=dpaD(iL)-dpP
          dpaE(iL)=dpG
          dpaE(iM)=0.d0
         else
          exit
         end if
        end do
       end do
      end if

      return
      end


C     ================================================================
      subroutine EIGSRT(dpaD,dpaV,dpaW,iN,iNP)

      implicit integer(i)
      implicit double precision(d)

      dimension dpaD(iNP),dpaV(iNP,iNP),dpaW(iNP,iNP)

      do i=1,iN-1
       iK=i
       dpP=dpaD(i)
       do iJ=i+1,iN
        if(dpaD(iJ).le.dpP)then
         iK=iJ
         dpP=dpaD(iJ)
        end if
       end do
       if(iK.ne.i)then
        dpaD(iK)=dpaD(i)
        dpaD(i)=dpP
        do iJ=1,iN
         dpP=dpaV(iJ,i)
         dpaV(iJ,i)=dpaV(iJ,iK)
         dpaV(iJ,iK)=dpP
         dpP=dpaW(iJ,i)
         dpaW(iJ,i)=dpaW(iJ,iK)
         dpaW(iJ,iK)=dpP
        end do
       end if
      end do

      return
      end


C     ================================================================
C     Evaluate transition densities and reduced transition
C     probabilities B(EL)
C     ================================================================
      subroutine ELECTRO(dpaA,dpaB,dpaFreq,dpaBel_IS,dpaBel_IV,
     &                   dpaFractS0_IS,dpaFractS0_IV,dpaFractS1_IS,
     &                   dpaFractS1_IV)

      implicit integer(i)
      implicit double precision(d)
      implicit character*20(s)

      include 'param.rpa'

      common/C_DIM/iNV,iNV1,iProt,iNeutr
      common/C_UNITS/dpPI,dpHBARC,dpAMC2,dpHBDM,dpPMass,dpNMass
      common/C_CTR/dpQBess,iISPIN,iORB,iPAR,iISFLIP,iRAD
      common/C_AB/iaNUM(IP_NCF)
      common/C_ASIS/iaLev(IP_NSP),iaNN(IP_NSP),
     &       iaLL(IP_NSP),iaLJ(IP_NSP),iaIQ(IP_NSP),iaJJ(IP_NCF),
     &       iaIPP(IP_NCF),iaINN(IP_NCF)
      common/C_ASI/dpaFOcc(IP_NSP),dpaSPE(IP_NSP)
      common/C_HF/iNMAX,iNOcc,iNUnocc,iNOrb
      common/C_MESHR/dpStep
      common/C_WF/dpaWF(IP_NNP,IP_NSP),dpaDWF(IP_NNP,IP_NSP)
      common/C_IQI/iaQI(IP_NCF)
      common/C_NUC/dpANUCL,dpZNUCL
      common/C_GFB/dpaECF(IP_NCF),dpaAKeep(IP_NCF)
      common/C_RADII/dpaR2L(12)
      common/C_RC/dpEnCentrMin,dpEnCentrMax
      common/C_DIPO/dpaDipoKAP(7),dpaIPGL(7),dpaINBes(7)
      common/C_ENERG/dpaSPEP(IP_NSP),dpaSPEN(IP_NSP)
      common/C_BESSL/dpaSJ(0:10),dpaDJ(0:10)
      common/C_GAMMA/dpExIni,dpExFin,dpEStep,dpGamma

      dimension dpaFreq(IP_NCF),dpaA(IP_NCF,IP_NCF),dpaB(IP_NCF,IP_NCF)
      dimension dpaEConf(IP_NCF),dpaBEConf(IP_NCF)
      dimension dpaBel_IS(IP_NCF),dpaBel_IV(IP_NCF)
      dimension dpaHFBel(IP_NCF),dpaBBel(IP_NCF)
      dimension dpaFractS0_IS(IP_NCF),dpaFractS0_IV(IP_NCF)
      dimension dpaFractS1_IV(IP_NCF),dpaFractS1_IS(IP_NCF)
      dimension dpaRMAT(2,IP_NCF,IP_NCF)
      dimension dpaRhoN(IP_NCF,IP_NNP),dpaRhoP(IP_NCF,IP_NNP)
      dimension dpaRMAT1(2,IP_NCF,IP_NCF)
      dimension dpaBelElec(IP_NCF)

      open(unit=13,status='unknown',file='td.out')

      dpETA=5.d0*dpaR2L(2)/3.d0
      do i=1,iNV1
       iX=iaIPP(iaNUM(i))
       iY=iaINN(iaNUM(i))
       dpaRMAT1(1,iX,iY)=0.d0
       dpaRMAT1(2,iX,iY)=0.d0
       if(iISFLIP.eq.0)then
        dpaRMAT1(1,iX,iY)=dpf_YL(iaLL(iX),iaLJ(iX),iaLL(iY),
     &                           iaLJ(iY),iaJJ(iaNUM(i)))
       end if
       if(iISFLIP.eq.1)then
        dpaRMAT1(1,iX,iY)=dpf_TJL(iaLL(iX),iaLJ(iX),iaLL(iY),
     &                            iaLJ(iY),iISPIN,iORB)
       end if
       dpaRMAT1(2,iX,iY)=dpaRMAT1(1,iX,iY)*iaQI(iaNUM(i))
      end do
      do i1=1,iNV1
       do i=1,iNMAX
        dpaRhoN(i1,i)=0.d0
        dpaRhoP(i1,i)=0.d0
        dpRDist=dfloat(i)*dpStep
        do i2=1,iNV1
         iQ2=iaIQ(iaIPP(iaNUM(i2)))+iaIQ(iaINN(iaNUM(i2)))
         if(iaJJ(iaNUM(i1)).ne.iaJJ(iaNUM(i2)))cycle
         if(iaJJ(iaNUM(i1)).ne.iISPIN)then
          write(2,100)iaJJ(iaNUM(i1)),iISPIN
         end if
         i_RPA_FAC=1
         if(iISFLIP.eq.1)i_RPA_FAC=if_POT(iISPIN+iORB)
         dpXY=dpaA(i2,i1)+i_RPA_FAC*dpaB(i2,i1)
         iEXP=(-iaLL(iaIPP(iaNUM(i2)))+iaLL(iaINN(iaNUM(i2)))+iORB)/2
         iSGN=if_POT(iEXP)
         if(iQ2.eq.2)then
          dpaRhoP(i1,i)=dpaRhoP(i1,i)+iSGN
     &                 *dpaRMAT1(1,iaIPP(iaNUM(i2)),iaINN(iaNUM(i2)))
     &                 *dpXY*dpaWF(i,iaIPP(iaNUM(i2)))
     &                 *dpaWF(i,iaINN(iaNUM(i2)))/(dpRDist**2)
         else if(iQ2.eq.0)then
          dpaRhoN(i1,i)=dpaRhoN(i1,i)+iSGN
     &                 *dpaRMAT1(1,iaIPP(iaNUM(i2)),iaINN(iaNUM(i2)))
     &                 *dpXY*dpaWF(i,iaIPP(iaNUM(i2)))
     &                 *dpaWF(i,iaINN(iaNUM(i2)))/(dpRDist**2)
         end if
        end do
        dpaRhoP(i1,i)=-dpaRhoP(i1,i)/dsqrt(dfloat(2*iISPIN+1))
        dpaRhoN(i1,i)=-dpaRhoN(i1,i)/dsqrt(dfloat(2*iISPIN+1))
       end do
      end do
      dpS0_IS=0.d0
      dpS0_IV=0.d0
      dpHFS1=0.d0
      dpBS1=0.d0
      dpS1_IS=0.d0
      dpS1_IV=0.d0
      dpSM1_IS=0.d0
      dpSM1_IV=0.d0
      dpS3_IS=0.d0
      dpS3_IV=0.d0
      dpS0_IS_Part=0.d0
      dpS1_IS_Part=0.d0
      dpS0_IV_Part=0.d0
      dpS1_IV_Part=0.d0
      do i1=1,iNV1
       dpEX1=dpaSPE(iaIPP(iaNUM(i1)))
       dpEY1=dpaSPE(iaINN(iaNUM(i1)))
       dpRhoElec=0.d0
       dpRhoInt_IS=0.d0
       dpRhoInt_IV=0.d0
       dpRDInt=0.d0
       dpFAC1=dpZNucl/dpANucl
       dpFAC2=(dpANucl-dpZNucl)/dpANucl
       do i=1,iNMax
        dpRDist=dfloat(i)*dpStep
        if(iRAD.EQ.1)then
         dpXOP=dpRDist**(iORB+2)
        else if(iRAD.eq.2)then
         dpXOP=dpRDist**(iORB+4)
        else if(iRAD.eq.3)then
         dpXOP=dpRDist**(iORB+2) 
         if(iORB.EQ.1)then
          dpXOP_IS=dpRDist**(iORB+4)-dpEta*dpRDist**(iORB+2)
         end if
        else if(iRAD.EQ.4)then
         call SPHJ(dpQBess*dpRDist)
         dpXOP=dpaSJ(iORB)*dpRDist**2
        end if
        dpRhoElec=dpRhoElec+(dpaRhoP(i1,i))*dpXOP
        if(iORB.eq.1.and.iRAD.eq.3)then
         dpRhoInt_IS=dpRhoInt_IS+(dpaRhoP(I1,I)+dpaRhoN(I1,I))*dpXOP_IS
C     IV dipole
         dpRhoInt_IV=dpRhoInt_IV+(dpFAC1*dpaRhoN(i1,i)
     &              -dpFAC2*dpaRhoP(i1,i))*dpXOP
        else
         dpRhoInt_IS=dpRhoInt_IS+(dpaRhoP(i1,I)+dpaRhoN(i1,i))*dpXOP
         dpRhoInt_IV=dpRhoInt_IV+(dpaRhoP(i1,i)-dpaRhoN(i1,i))*dpXOP
        endif
        dpRDInt=dpRDInt+dpaWF(i,iaIPP(iaNUM(i1)))
     &         *dpaWF(i,iaINN(iaNUM(i1)))*dpXOP/(dpRDist**2)
       end do
       dpRhoElec=dpRhoElec*dpStep
       dpRhoInt_IS=dpRhoInt_IS*dpStep
       dpRhoInt_IV=dpRhoInt_IV*dpStep
       dpRDInt=dpRDInt*dpStep
       dpaRMAT(1,iaIPP(iaNUM(i1)),iaINN(iaNUM(i1)))
     &  =dpaRMAT1(1,iaIPP(iaNUM(i1)),iaINN(iaNUM(i1)))*dpRDInt
       dpaRMAT(2,iaIPP(iaNUM(i1)),iaINN(iaNUM(i1)))
     &  =dpaRMAT1(2,iaIPP(iaNUM(i1)),iaINN(iaNUM(i1)))*dpRDInt
       dpQRhoElec=dpRhoElec**2
       dpQRhoInt_IS=dpRhoInt_IS**2
       dpQRhoInt_IV=dpRhoInt_IV**2
       dpaBel_IS(i1)=dpQRhoInt_IS*(2*iISPIN+1)
       dpaBel_IV(i1)=dpQRhoInt_IV*(2*iISPIN+1)
       dpaBelElec(i1)=dpQRhoElec*(2*iISPIN+1)
       dpaHFBel(i1)=0.d0
       if(iaIPP(iaNUM(i1)).ne.iaINN(iaNUM(i1)))then
        dpaEConf(i1)=dabs(dpEX1-dpEY1)
        dpHFOCC=dsqrt((1.d0-dpaFOcc(iaIPP(iaNUM(i1))))
     &         *(dpaFOcc(iaINN(iaNUM(i1)))))
     &         +dsqrt((dpaFOcc(iaIPP(iaNUM(i1))))
     &         *(1.d0-dpaFOcc(iaINN(iaNUM(i1)))))
        dpaHFBel(i1)=(dpaRMAT(1,iaIPP(iaNUM(i1))
     &               ,iaINN(iaNUM(i1)))*dpHFOCC)**2
        iQ1=iaIQ(iaIPP(iaNUM(i1)))+iaIQ(iaINN(iaNUM(i1)))
        if(iQ1.eq.2)then
         dpaBEConf(i1)=dpaSPEP(iaIPP(iaNUM(i1)))
     &                -dpaSPEP(iaINN(iaNUM(i1)))
         dpaBBel(i1)=(dpaRMAT(1,iaIPP(iaNUM(i1)),iaINN(iaNUM(i1))))**2
        else if(iQ1.eq.0)then
         dpaBEConf(i1)=dpaSPEN(iaIPP(iaNUM(i1)))
     &                -dpaSPEN(iaINN(iaNUM(i1)))
         dpaBBel(i1)=(dpaRMAT(1,iaIPP(iaNUM(i1)),iaINN(iaNUM(i1))))**2
        end if
       end if
       if(iISPIN.ne.1.or.I1.ne.1.or.iRAD.ne.3)then
        dpHFS1=dpHFS1+dpaEConf(i1)*dpaHFBel(i1)
        dpBS1=dpBS1+dpaBEConf(i1)*dpaBBel(i1)
        dpS0_IS=dpS0_IS+dpaBEL_IS(i1)
        dpS1_IS=dpS1_IS+dpaFreq(i1)*dpaBEL_IS(i1)
        dpSM1_IS=dpSM1_IS+dpaBEL_IS(i1)/dpaFreq(i1)
        dpS3_IS=dpS3_IS+dpaBEL_IS(i1)*dpaFreq(i1)**3
        dpS0_IV=dpS0_IV+dpaBEL_IV(i1)
        dpS1_IV=dpS1_IV+dpaFreq(i1)*dpaBEL_IV(i1)
        dpSM1_IV=dpSM1_IV+dpaBEL_IV(i1)/dpaFreq(i1)
        dpS3_IV=dpS3_IV+dpaBEL_IV(i1)*dpaFreq(i1)**3
        if(dpaFreq(i1).gt.dpEnCentrMin.and.
     &     dpaFreq(i1).lt.dpEnCentrMax)then
         dpS0_IS_PART=dpS0_IS_PART+dpaBEL_IS(i1)
         dpS1_IS_PART=dpS1_IS_PART+dpaBEL_IS(i1)*dpaFreq(i1)
         dpS0_IV_Part=dpS0_IV_Part+dpaBel_IV(I1)
         dpS1_IV_Part=dpS1_IV_Part+dpaBel_IV(I1)*dpaFreq(I1)
        end if
       end if
      end do
      do i=1,iNV1
       if(dpS0_IS.ne.0.d0)then
        dpaFRACTS0_IS(i)=dpaBEL_IS(i)/dpS0_IS*100.d0
       else
        dpaFRACTS0_IS(i)=0.d0
       end if
       if(dpS1_IS.ne.0.d0)then
        dpaFRACTS1_IS(i)=dpaFreq(i)*dpaBEL_IS(i)/dpS1_IS*100.d0
       else
        dpaFRACTS1_IS(i)=0.d0
       end if
       if(dpS0_IV.ne.0.d0)then
        dpaFRACTS0_IV(i)=dpaBEL_IV(i)/dpS0_IV*100.d0
       else
        dpaFRACTS0_IV(i)=0.d0
       end if
       if(dpS1_IV.ne.0.d0)then
        dpaFRACTS1_IV(i)=dpaFreq(i)*dpaBEL_IV(i)/dpS1_IV*100.d0
       else
        dpaFRACTS1_IV(i)=0.d0
       end if
      end do
      if(iRAD.eq.1)write(2,101)
      if(iRAD.eq.2)write(2,102)
      if(iRAD.eq.3)write(2,103)
      if(iRAD.eq.4)write(2,104)
      write(2,110)
      do i=1,iNV1
       write(2,120)i,dpaFreq(i),dpaBEL_IS(i),
     &             dpaFRACTS0_IS(i),dpaBelElec(i),
     &             dpaBEL_IV(i),dpaFRACTS0_IV(i)
      end do
      iNV2=int(dpExFin-dpExIni)/dpEStep+1 
      sFileName='Plot_Bel_IS.dat'
      call AVERLOR(sFileName,dpaFreq,dpaBEL_IS,iNV2)
      sFileName='Plot_Bel_EM.dat'
      call AVERLOR(sFileName,dpaFreq,dpaBelElec,iNV2)
      sFileName='Plot_Bel_IV.dat'
      call AVERLOR(sFileName,dpaFreq,dpaBEL_IV,iNV2)
      write(2,*)''
      write(2,*)'SUM RULES FOR IS:'
      write(2,130)dpSM1_IS,dpS0_IS,dpS1_IS,dpS3_IS
      iCHO=iORB
      if(iRAD.EQ.1)then
       dpAM1_DC=dpHBDM*(dpANucl-1.d0)*iORB*(2*iORB+1)**2/4.d0/dpPI
     &         *dpaR2L(iCHO)
       write(2,210)dpAM1_DC,dpS1_IS/dpAM1_DC*100.d0
      else if(iRAD.eq.2)then
       dpAM1_DC=dpHBDM*(dpANucl-1.d0)*dpaR2L(2)/dpPI
       write(2,210)dpAM1_DC,dpS1_IS/dpAM1_DC*100.d0
      else if(iRAD.eq.3)then
       iCHO=iORB+2
       dpAM1_DC=dpHBDM*(dpANucl-1.d0)*(33.d0*dpaR2L(iCHO)
     &         -25.d0*dpaR2L(2)**2)/4.d0/dpPI
       write(2,210)dpAM1_DC,dpS1_IS/dpAM1_DC*100.d0
      else if(iRAD.eq.4)then
       dpAM1_DC=dpHBDM*(2*iORB+1)*dpaInBes(iORB+1)/4.d0/dpPI
     &         *(dpANucl-1.d0)/dpANucl
       write(2,210)dpAM1_DC,dpS1_IS/dpAM1_DC*100.d0
      end if
      write(2,*)''
      write(2,*)'SUM RULES FOR IV:'
      write(2,130)dpSM1_IV,dpS0_IV,dpS1_IV,dpS3_IV
      if(iRAD.eq.1)then
       dpAKappa=dpaDipoKAP(iORB+1)
       dpaKappaf=4.d0*(dpANucl-dpZNucl)*dpZNucl
     &          /dpANucl/(dpANucl-1.d0)
       dpAM1_DC=dpHBDM*(dpANucl-1.d0)*iORB*(2*iORB+1)**2/4.d0/dpPI
     &         *dpaR2L(ICHO)*(1+dpaKappaf*(dpAKappa))
       write(2,210)dpAM1_DC,dpS1_IV/dpAM1_DC*100.d0
      else if(iRAD.eq.2)then
       dpAKappa=dpaDipoKAP(iORB+1)
       dpaKappaf=4.d0*(dpANucl-dpZNucl)*dpZNucl
     &          /dpANucl/(dpANucl-1.d0)
       dpAM1_DC=dpHBDM*(dpANucl-1.d0)*dpAR2L(2)/dpPI
     &         *(1+dpaKappaf*(dpAKappa))
       write(2,210)dpAM1_DC,dpS1_IV/dpAM1_DC*100.d0
      else if(iRAD.eq.3)then
       dpAKappa=dpaDipoKAP(iORB+1)
       dpaKappaf=4.d0*(dpANucl-dpZNucl)*dpZNucl
     &          /dpANucl/(dpANucl-1.d0)
       dpTRK=14.85d0*(dpANucl-dpZNucl)*dpZNucl/dpANucl
       dpTRK=dpTRK*(dpANucl-1.d0)/dpANucl
       write(2,140)dpTRK
       dpAKappa=dpAKappa*dpANucl/(dpANucl-1.d0)
       dpAM1_DC=dpTRK*(1+dpAKappa)
       write(2,210)dpAM1_DC,dpS1_IV/dpAM1_DC*100.d0
      else if(iRAD.eq.4)then
       dpAKappa=dpaDipoKAP(iORB+1)
       dpaKappaf=4.d0*(dpANucl-dpZNucl)*dpZNucl/dpANucl/dpANucl
       dpAM1_DC=dpHBDM*(2*iORB+1)*dpaInBes(iORB+1)
     &         *(1+dpaKappaf*(dpAKappa))/4.d0/dpPI
     &         *(dpANucl-1.d0)/dpANucl
       write(2,210)dpAM1_DC,dpS1_IV/dpAM1_DC*100.d0
      end if
      write(2,*)''
      write(2,*)'CONSTRAINED, CENTROID, SCALING ENERGIES:'
      write(2,150)dsqrt(dpS1_IS/dpSM1_IS),dpS1_IS/dpS0_IS,
     &            dsqrt(dpS3_IS/dpS1_IS)
      write(2,*)''
      write(2,160)dsqrt(dpS1_IV/dpSM1_IV),dpS1_IV/dpS0_IV,
     &            dsqrt(dpS3_IV/dpS1_IV)
      if(dpEnCentrMin.gt.0.1d0.and.dpEnCentrMax.gt.0.1d0)then
       write(2,*)
       write(2,170)dpEnCentrMin,dpEnCentrMax
       write(2,180)dpS0_IS_Part,dpS1_IS_Part,dpS1_IS_Part/dpS0_IS_Part
       write(2,190)dpS0_IV_Part,dpS1_IV_Part,dpS1_IV_Part/dpS0_IV_Part
      end if
      write(2,*)
      iCount=0
      do i1=1,iNV1
       iK=i1
       dpX=dpaECF(i1)
       dpAAA=dpaAKeep(i1)
       do i2=i1,iNV1
        if(dpX-dpaECF(i2).gt.0.d0)then
         iK=i2
         dpX=dpaECF(i2)
         dpAAA=dpaAKeep(i2)
        end if
       end do
       dpaECF(iK)=dpaECF(i1)
       dpaAKeep(iK)=dpaAKeep(i1)
       dpaECF(i1)=dpX
       dpaAKeep(i1)=dpAAA
      end do
      do i1=1,iNV1
       iCount=iCount+1
       write(13,200)iCount,iORB,iISPIN,iISFLIP,dpaFreq(i1)
       write(13,*)dpaFractS0_IS(i1),dpaFractS0_IV(i1)
       do i=1,iNMAX
        dpRDist=dfloat(i)*dpStep
        write(13,*)dpRDist,dpaRhoP(i1,i),dpaRhoN(i1,i)
       end do
      end do

      close(13)

  100 format(/,10X,' WARNING: SPIN ERROR ! ',2X,'JT1=',I1,' ISPIN=',
     &       I1)
  101 format(//,' OPERATORS PROPORTIONAL TO r**J ARE EMPLOYED'/)
  102 format(//,' OPERATORS PROPORTIONAL TO r**(J+2) ARE EMPLOYED'/)
  103 format(//,' DIPOLE OPERATORS ARE EMPLOYED'/)
  104 format(//,' OPERATORS PROPORTIONAL TO BESSEL FUNCTIONS 
     &      ARE EMPLOYED'/)
  110 format(/,25('*'),' RPA STATES ',25('*'),//,
     &       '          E          B(IS)       %M0(IS)       B(EM)',
     &       '        B(IV)       %M0(IV)',/,
     &       '      ---------    ---------    ---------    ---------',
     &       '    ---------    ---------',/)
  120 format(I4,6(E12.5,1X))
  130 format('  M(-1)                 =',E12.6,/
     &       '  M(0)                  =',E12.6,/
     &       '  M(1)                  =',E12.6,/
     &       '  M(3)                  =',E12.6)
  140 format(' TRK SUM RULE FOR IVGDR =',E12.6)
  150 format('  SQRT[M(1)/M(-1)] (IS) =',E12.6,/
     &       '  M(1)/M(0)        (IS) =',E12.6,/
     &       '  SQRT[M(3)/M(1)]  (IS) =',E12.6)
  160 format('  SQRT[M(1)/M(-1)] (IV) =',E12.6,/
     &       '  M(1)/M(0)        (IV) =',E12.6,/
     &       '  SQRT[M(3)/M(1)]  (IV) =',E12.6)
  170 format('BETWEEN ',E12.6,' AND ',E12.6,' :')
  180 format('M(0),M(1),M(1)/M(0) IS = ',E12.6,2X,E12.6,2X,E12.6)
  190 format('M(0),M(1),M(1)/M(0) IV = ',E12.6,2X,E12.6,2X,E12.6)
  200 format(1X,4(I4,1X),E12.5)
  210 format('  M(1) D.C.             =',E12.6,/  
     &       '  % EXHAU.              =',E12.6,' %')

      return
      end


C     ================================================================
C     Reduced matrix elements of Y_L coupled to sigma (total spin J)
C     ================================================================
      function dpf_TJL(iLA,iJA,iLB,iJB,iJ,iL)

      implicit integer(i)
      implicit double precision(d)

      dpTJL=0.d0
      dpZ=0.d0
      dpQP=4.d0*3.14159265d0
      dpFJ=dfloat(iJ)
      dpFL=dfloat(iL)
      iJD=2*iJ
      iLP=iLA+iLB
      iLM=iabs(iLA-iLB)
      iLL=(iL-iLP)*(iL-iLM)
      if(iLL.le.0)then
       iLL=iLA+iLB+iL-2*((iLA+iLB+iL)/2)
       if(iLL.eq.0)then
        iJP=(iJA+iJB)/2
        iJJ=(iJ-iJP)*(iJ-iabs(iJA-iJB)/2)
        if(iJJ.le.0)then
         iLJP=iL+iJ
         iLJM=iabs(iL-iJ)
         iLJ=(1-iLJP)*(1-iLJM)
         if(iLJ.le.0)then
          if(iL.eq.iJ)then
           iIG=(iJA+iJB)/2+iJ
           dpCC=(iJD+1.d0)/(dpFJ*(dpFJ+1.d0))
           dpZ=-0.5d0*dsqrt(dpCC)*(iJB+1.d0+(1.d0-2.d0*mod(iIG,2))
     &        *(iJA+1.d0))
          else
           iIG=iLB+(iJB+1)/2
           dpCC=(iLA-0.5d0*iJA)*(iJA+1.d0)+(iLB-0.5d0*iJB)*(iJB+1.d0)
           dpSG=1.d0-2.d0*mod(iIG,2)
           if(iL.eq.(iJ+1))dpZ=dpSG*(dpCC+dpFL)/dsqrt(dpFL)
           if(iL.eq.(iJ-1))dpZ=dpSG*(dpCC-dpFL-1.d0)/dsqrt(dpFL+1.d0)
          end if
          dpC=(iJA+1.d0)*(iJB+1.d0)/((iJD+1.d0)*dpQP)
          dpXJA=iJA/2.d0
          dpXJB=iJB/2.d0
          dpXJD=iJ
          iJPHA=(iJA+iJB)/2+1
          iJPHA=1-2*mod(iJPHA,2)
          dpTJL=(1.d0-2.d0*mod(iLA,2))*dsqrt(dpC)*dpZ*iJPHA
     &         *dpf_COFCG(dpXJA,dpXJB,dpXJD,0.5d0,-0.5d0,0.d0)
         end if
        end if
       end if
      end if
      dpf_TJL=dpTJL

      return
      end


C     ================================================================
C     Reduced matrix elements of Y_L
C     ================================================================
      function dpf_YL(iLA,iJA,iLB,iJB,iL)

      implicit integer(i)
      implicit double precision(d)

      dpQP=4.d0*3.14159265d0
      dpYL=0.d0
      iLMinA=iabs(iJA-iJB)/2
      iLMaxA=(iJA+iJB)/2
      if(iL.ge.iLMinA)then
       if(iL.le.iLMaxA)then
        iPAR=1-2*mod(iLA+iLB+iL,2)
        if(iPAR.ge.0)then
         dpSG=1.d0-2.d0*mod(iL+(iJB-1)/2,2)
         dpXJA=iJA/2.d0
         dpXJB=iJB/2.d0
         dpXL=iL
         iJPHA=(iJA+iJB)/2+1
         iJPHA=1-2*mod(iJPHA,2)
         dpC=(iJA+1.d0)*(iJB+1.d0)/dpQP
         dpSG=dpSG*iJPHA
         dpYL=dpSG*dsqrt(dpC)
     &       *dpf_COFCG(dpXJA,dpXJB,dpXL,0.5d0,-0.5d0,0.d0)
        end if
       end if
      end if
      dpf_YL=dpYL

      return
      end


C     ================================================================
C     Clebsch-Gordan coefficients
C     ================================================================
      function dpf_COFCG(dpA,dpB,dpC,dpX,dpY,dpZ)

      implicit integer(i)
      implicit double precision(d)

      dpXRCNP=dabs(dpA)+dabs(dpB)+dabs(dpC)
      if(dpXRCNP.ge.200.d0)then
       write(9,*)'XJ1, XJ2, XJ3; A, B, C  AT POINT5 =',
     &           dpA,dpB,dpC,dpA,dpB,dpC
       stop
      end if
      call TROISJ(dpA,dpB,dpC,dpX,dpY,-dpZ,dpC3J)
      dpCOF=dabs(dpA-dpB+dpZ)+1.d-06
      iCOF=int(dpCOF)
      dpCOF=1.d0-2.d0*mod(iCOF,2)
      dpf_COFCG=dpC3J*dpCOF*dsqrt(2.d0*dpC+1.d0)

      return
      end


C     ================================================================
      function if_POT(iVal)

      implicit integer(i)
      implicit double precision(d)

      dpA=iVal
      i=iVal/2
      dpP=dpA/2.0
      dpQ=i
      if(dpP-dpQ.eq.0.d0)then
       if_POT=1
      else
       if_POT=-1
      end if

      return
      end


C     ================================================================
C     Purpose: Compute spherical Bessel functions yn(x) and
C              their derivatives
C     Input :  dpX --- Argument of yn(x) ( x ?0 )
C              iN --- Order of yn(x) ( n = 0,1,?)
C     Output:  dpaSJ(n) --- yn(x)
C              dpaDJ(n) --- yn'(x)
C              iNM --- Highest order computed
C     ================================================================
      SUBROUTINE SPHJ(dpX)

      implicit integer(i)
      implicit double precision(d)

      common/C_BESSL/dpaSJ(0:10),dpaDJ(0:10)

      iN=10
      iNM=iN
      if(dpX.lt.1.0d-60)then
       do iK=0,iN
        dpaSJ(iK)=-1.0d+300
        dpaDJ(iK)=1.0d+300
       end do
      else
       dpaSJ(0)=dsin(dpX)/dpX
       dpaSJ(1)=(dpaSJ(0)-dcos(dpX))/dpX
       dpF0=dpaSJ(0)
       dpF1=dpaSJ(1)
       do iK=2,iN
        dpF=(2.0d0*dfloat(iK)-1.0d0)*dpF1/dpX-dpF0
        dpaSJ(iK)=dpF
        if(dabs(dpF).ge.1.0d+300)exit
        dpF0=dpF1
        dpF1=dpF
       end do
       dpaDJ(0)=(dcos(dpX)-dsin(dpX)/dpX)/dpX
       do iK=1,iNM-1
        dpaDJ(iK)=(dfloat(iK)*dpaSJ(iK-1)-(dfloat(iK)+1.0d0)
     &           *dpaSJ(iK+1))/(2.d0*dfloat(iK)+1.d0)
       end do
      end if

      return
      end


C     ================================================================
C     AVERLOR is a short routine which makes an averaging with
C     Lorentzians of a given strength function
C     ================================================================
      subroutine AVERLOR(sFile,dpaEx,dpaBE,iNum)
      implicit integer(i)
      implicit double precision(d)
      implicit character*20(s)

      include 'param.rpa'

      common/C_CTR/dpQBess,iISPIN,iORB,iPAR,iISFLIP,iRAD
      common/C_UNITS/dpPI,dpHBARC,dpAMC2,dpHBDM,dpPMass,dpNMass
      common/C_GAMMA/dpExIni,dpExFin,dpEStep,dpGamma

      dimension dpaEx(iNum),dpaBE(iNum)
      dimension dpaSaf(iNum),dpaEx1(iNum)

      dpTPI=2.d0*dpPI
      iP_NS=int((dpExFin-dpExIni)/dpEStep)
      iNus=1
      if(iRAD.eq.3.and.iORB.eq.1)iNus=2
      if(iRAD.eq.4.and.iORB.eq.1)iNus=2
      do i=1,iP_NS
       dpaEx1(i)=dpExIni+dpEStep*dfloat(i-1)
       dpaSaf(i)=0.d0
       do ij=iNus,iNum
        dpFact=(dpGamma/dpTPI)*1.d0
     &        /((dpaEx1(i)-dpaEx(ij))**2+dpGamma**2/4.d0)
        dpaSaf(i)=dpaSaf(i)+dpaBE(ij)*dpFact
       end do
      end do
      open(unit=5,status='unknown',file=sFile)
      daint1=0.d0
      daint2=0.d0
      do i=1,iP_NS
       write(5,100)dpaEx1(i),dpaSaf(i)
       if(i.ne.1)then
        dpDeltaEn=dpaEx1(i)-dpaEx1(i-1)
        daint1=daint1+dpDeltaEn*(dpaSaf(i)+dpaSaf(i-1))/2.d0
        daint2=daint2+dpDeltaEn
     &        *(dpaSaf(i)*dpaEx1(i)+dpaSaf(i-1)*dpaEx1(i-1))/2.d0
       end if
      end do
      close(5)

 100  format(1X,5(E12.6,1X))

      return
      end