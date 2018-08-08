      module sascnv_mod
       implicit none
        private
        public :: SASCNV 

       contains

      SUBROUTINE SASCNV(KM,DELT,DEL,PRSL,PS,PHIL,QL,
     &       Q1,T1,U1,V1,CLDWRK,RN,KBOT,KTOP,KUO,SLIMSK,
     &       DOT,XKT2,ncloud)

      USE constants_mod, only: grav, CP => CP_AIR, HVAP => HLV 
     &,             RV => RVGAS, T0C => KELVIN, CLIQ => CP_WATER
     &,             RDGAS, CVAP => CP_VAPOR

      implicit none
      
      real, parameter :: FV = RV/RDGAS-1.
      real, parameter :: eps = RDGAS/RV
      real, parameter :: epsm1 = RDGAS/RV-1.
!
!
      integer            KM, ncloud,
     &                   KBOT, KTOP, KUO
      real DELT
      real PS,     DEL(KM),  PRSL(KM),
!     real             DEL(KM),  PRSL(KM),
     &                     QL(KM,2),Q1(KM),   T1(KM),
     &                     U1(KM),  V1(KM),  
     &                     CLDWRK, RN,      SLIMSK,
     &                     DOT(KM), XKT2,    PHIL(KM)
!
      integer              I, INDX, jmn, k, knumb, latd, lond, km1
!
      real adw,     alpha,   alphal,  alphas,
     &                     aup,     beta,    betal,   betas,
     &                     c0,      cpoel,   dellat,  delta,
     &                     desdt,   deta,    detad,   dg,
     &                     dh,      dhh,     dlnsig,  dp,
     &                     dq,      dqsdp,   dqsdt,   dt,
     &                     dt2,     dtmax,   dtmin,   dv1,
     &                     dv1q,    dv2,     dv2q,    dv1u,
     &                     dv1v,    dv2u,    dv2v,    dv3u,
     &                     dv3v,    dv3,     dv3q,    dvq1,
     &                     dz,      dz1,     e1,      edtmax,
     &                     edtmaxl, edtmaxs, el2orc,  elocp,
     &                     es,      etah,
     &                     evef,    evfact,  evfactl, fact1,
     &                     fact2,   factor,  fkm,
     &                     fuv,     g,       gamma,   onemf,
     &                     onemfu,  pdetrn,  pdpdwn,  pprime,
     &                     qc,      qlk,     qrch,    qs,
     &                     rain,    rfact,   shear,   tem1,
     &                     tem2,    terr,    val,     val1,
     &                     val2,    w1,      w1l,     w1s,
     &                     w2,      w2l,     w2s,     w3,
     &                     w3l,     w3s,     w4,      w4l,
     &                     w4s,     xdby,    xpw,     xpwd,
     &                     xqc,     xqrch,   xlambu,  mbdt,
     &                     tem
!
!
      integer              JMIN, KB, KBCON, KBDTR,
     &                     KT2,  KTCON, LMIN,
     &                     kbm,  kbmax, kmax
!
      real AA1,     ACRT,   ACRTFCT,
     &                     DELHBAR, DELQ,   DELQ2,
     &                     DELQBAR, DELQEV, DELTBAR,
     &                     DELTV,   DTCONV, EDT,
     &                     EDTO,    EDTX,   FLD,
     &                     HCDO,    HKBO,   HMAX,
     &                     HMIN,    HSBAR,  UCDO,
     &                     UKBO,    VCDO,   VKBO,
     &                     PBCDIF,  PDOT,   PO(KM),
     &                                  PWAVO,  PWEVO,
!    &                     PSFC,    PWAVO,  PWEVO,
     &                     QCDO,    QCOND,  QEVAP,
     &                     QKBO,    RNTOT,  VSHEAR,
     &                     XAA0,    XHCD,   XHKB,
     &                     XK,      XLAMB,  XLAMD,
     &                     XMB,     XMBMAX, XPWAV,
     &                     XPWEV,   XQCD,   XQKB
cc
C  PHYSICAL PARAMETERS
      PARAMETER(G=grav)
      PARAMETER(CPOEL=CP/HVAP,ELOCP=HVAP/CP,
     &          EL2ORC=HVAP*HVAP/(RV*CP))
      PARAMETER(TERR=0.,C0=.002,DELTA=fv)
      PARAMETER(FACT1=(CVAP-CLIQ)/RV,FACT2=HVAP/RV-FACT1*T0C)
C  LOCAL VARIABLES AND ARRAYS
      real PFLD(KM),    TO(KM),     QO(KM),
     &                     UO(KM),      VO(KM),     QESO(KM)
c  cloud water
      real QLKO_KTCON, DELLAL,    TVO(KM),
     &                     DBYO(KM),    ZO(KM),     SUMZ(KM),
     &                     SUMH(KM),    HEO(KM),    HESO(KM),
     &                     QRCD(KM),    DELLAH(KM), DELLAQ(KM),
     &                     DELLAU(KM),  DELLAV(KM), HCKO(KM),
     &                     UCKO(KM),    VCKO(KM),   QCKO(KM),
     &                     ETA(KM),     ETAU(KM),   ETAD(KM),
     &                     QRCDO(KM),   PWO(KM),    PWDO(KM),
     &                     RHBAR,      TX1
!
      LOGICAL TOTFLG, CNVFLG, DWNFLG, DWNFLG2, FLG
!
      real PCRIT(15), ACRITT(15), ACRIT(15)
cmy      SAVE PCRIT, ACRITT
      DATA PCRIT/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,
     &           350.,300.,250.,200.,150./
      DATA ACRITT/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,
     &           .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
C  GDAS DERIVED ACRIT
C     DATA ACRITT/.203,.515,.521,.566,.625,.665,.659,.688,
C    &            .743,.813,.886,.947,1.138,1.377,1.896/
cc
      real tf, tcr, tcrf
      parameter (TF=233.16, TCR=263.16, TCRF=1.0/(TCR-TF)) ! From Lord(1978)
!
!     parameter (tf=258.16, tcr=273.16, tcrf=1.0/(tcr-tf))
!
      real, parameter :: cons_0=0.0
c
c--------------------------------------------------------------------
!
      km1 = km - 1
C  INITIALIZE ARRAYS
C
        RN=0.
        KBOT=KM+1
        KTOP=0
!       KUO=0
        CNVFLG = .TRUE.
        DTCONV = 3600.
        CLDWRK = 0.
        PDOT = 0.
        KT2 = 0
        QLKO_KTCON = 0.
        DELLAL = 0.
!!
      DO K = 1, 15
        ACRIT(K) = ACRITT(K) * (975. - PCRIT(K))
      ENDDO
      DT2 = DELT
      val   =         1200.
      dtmin = max(dt2, val )
      val   =         3600.
      dtmax = max(dt2, val )
C  MODEL TUNABLE PARAMETERS ARE ALL HERE
      MBDT    = 10.
      EDTMAXl = .3
      EDTMAXs = .3
      ALPHAl  = .5
      ALPHAs  = .5
      BETAl   = .15
      betas   = .15
      BETAl   = .05
      betas   = .05
c     EVEF    = 0.07
      evfact  = 0.3
      evfactl = 0.3
      PDPDWN  = 0.
      PDETRN  = 200.
      xlambu  = 1.e-4
      val     =           1.
      fkm     = (float(km) / 28.) ** 2
      fkm     = max(fkm,val)
      W1l     = -8.E-3
      W2l     = -4.E-2
      W3l     = -5.E-3
      W4l     = -5.E-4
      W1s     = -2.E-4
      W2s     = -2.E-3
      W3s     = -1.E-3
      W4s     = -2.E-5
CCCCC IF(IM.EQ.384) THEN
        LATD  = 92
        lond  = 189
CCCCC ELSEIF(IM.EQ.768) THEN
CCCCC   LATD = 80
CCCCC ELSE
CCCCC   LATD = 0
CCCCC ENDIF
C
C  DEFINE TOP LAYER FOR SEARCH OF THE DOWNDRAFT ORIGINATING LAYER
C  AND THE MAXIMUM THETAE FOR UPDRAFT
C
        KBMAX = KM
        KBM   = KM
        KMAX  = KM
        TX1   = 1.0 / PS
!
      DO K = 1, KM
          IF (prSL(K)*tx1 .GT. 0.45) KBMAX = K + 1
          IF (prSL(K)*tx1 .GT. 0.70) KBM   = K + 1
          IF (prSL(K)*tx1 .GT. 0.04) KMAX  = K + 1
      ENDDO
        KBMAX = MIN(KBMAX,KMAX)
        KBM   = MIN(KBM,KMAX)
C
C   CONVERT SURFACE PRESSURE TO MB FROM CB
C
!!
      DO K = 1, KM
          if (K .le. kmax) then
            PFLD(k) = PRSL(K) * 10.0
            PWO(k)  = 0.
            PWDO(k) = 0.
            TO(k)   = T1(k)
            QO(k)   = Q1(k)
            UO(k)   = U1(k)
            VO(k)   = V1(k)
            DBYO(k) = 0.
            SUMZ(k) = 0.
            SUMH(k) = 0.
          endif
      ENDDO
C
C  COLUMN VARIABLES
C  P IS PRESSURE OF THE LAYER (MB)
C  T IS TEMPERATURE AT T-DT (K)..TN
C  Q IS MIXING RATIO AT T-DT (KG/KG)..QN
C  TO IS TEMPERATURE AT T+DT (K)... THIS IS AFTER ADVECTION AND TURBULAN
C  QO IS MIXING RATIO AT T+DT (KG/KG)..Q1
C
      DO K = 1, KM
          if (k .le. kmax) then
!jfe        QESO(k) = 10. * FPVS(T1(k))
!
            QESO(k) = 0.01 * fpvs(T1(K))      ! fpvs is in Pa
!
            QESO(k) = EPS * QESO(k) / (PFLD(k) + EPSM1*QESO(k))
            val1      =             1.E-8
            QESO(k) = MAX(QESO(k), val1)
            val2      =           1.e-10
            QO(k)   = max(QO(k), val2 )
c           QO(k)   = MIN(QO(k),QESO(k))
            TVO(k)  = TO(k) + DELTA * TO(k) * QO(k)
          endif
      ENDDO
C
C  HYDROSTATIC HEIGHT ASSUME ZERO TERR
C
      DO K = 1, KM
          ZO(k) = PHIL(k) / G
      ENDDO
C  COMPUTE MOIST STATIC ENERGY
      DO K = 1, KM
          if (K .le. kmax) then
!           tem       = G * ZO(k) + CP * TO(k)
            tem       = PHIL(k) + CP * TO(k)
            HEO(k)  = tem  + HVAP * QO(k)
            HESO(k) = tem  + HVAP * QESO(k)
C           HEO(k)  = MIN(HEO(k),HESO(k))
          endif
      ENDDO
C
C  DETERMINE LEVEL WITH LARGEST MOIST STATIC ENERGY
C  THIS IS THE LEVEL WHERE UPDRAFT STARTS
C
        HMAX = HEO(1)
        KB = 1
!!
      DO K = 2, KM
          if (k .le. kbm) then
            IF(HEO(k).GT.HMAX.AND.CNVFLG) THEN
              KB   = K
              HMAX = HEO(k)
            ENDIF
          endif
      ENDDO
C     DO K = 1, KMAX - 1
C         TOL(k) = .5 * (TO(k) + TO(k+1))
C         QOL(k) = .5 * (QO(k) + QO(k+1))
C         QESOL(k) = .5 * (QESO(k) + QESO(k+1))
C         HEOL(k) = .5 * (HEO(k) + HEO(k+1))
C         HESOL(k) = .5 * (HESO(k) + HESO(k+1))
C     ENDDO
      DO K = 1, KM1
          if (k .le. kmax-1) then
            DZ      = .5 * (ZO(k+1) - ZO(k))
            DP      = .5 * (PFLD(k+1) - PFLD(k))
!jfe        ES      = 10. * FPVS(TO(k+1))
!
            ES      = 0.01 * fpvs(TO(K+1))      ! fpvs is in Pa
!
            PPRIME  = PFLD(k+1) + EPSM1 * ES
            QS      = EPS * ES / PPRIME
            DQSDP   = - QS / PPRIME
            DESDT   = ES * (FACT1 / TO(k+1) + FACT2 / (TO(k+1)**2))
            DQSDT   = QS * PFLD(k+1) * DESDT / (ES * PPRIME)
            GAMMA   = EL2ORC * QESO(k+1) / (TO(k+1)**2)
            DT      = (G * DZ + HVAP * DQSDP * DP) / (CP * (1. + GAMMA))
            DQ      = DQSDT * DT + DQSDP * DP
            TO(k) = TO(k+1) + DT
            QO(k) = QO(k+1) + DQ
            PO(k) = .5 * (PFLD(k) + PFLD(k+1))
          endif
      ENDDO
!
      DO K = 1, KM1
          if (k .le. kmax-1) then
!jfe        QESO(k) = 10. * FPVS(TO(k))
!
            QESO(k) = 0.01 * fpvs(TO(K))      ! fpvs is in Pa
!
            QESO(k) = EPS * QESO(k) / (PO(k) + EPSM1*QESO(k))
            val1      =             1.E-8
            QESO(k) = MAX(QESO(k), val1)
            val2      =           1.e-10
            QO(k)   = max(QO(k), val2 )
c           QO(k)   = MIN(QO(k),QESO(k))
            HEO(k)  = .5 * G * (ZO(k) + ZO(k+1)) +
     &                  CP * TO(k) + HVAP * QO(k)
            HESO(k) = .5 * G * (ZO(k) + ZO(k+1)) +
     &                  CP * TO(k) + HVAP * QESO(k)
            UO(k)   = .5 * (UO(k) + UO(k+1))
            VO(k)   = .5 * (VO(k) + VO(k+1))
          endif
      ENDDO
c     k = kmax
c       HEO(k) = HEO(k)
c       hesol(k) = HESO(k)
c      IF(LAT.EQ.LATD.AND.lon.eq.lond.and.CNVFLG) THEN
c        PRINT *, '   HEO ='
c        PRINT 6001, (HEO(K),K=1,KMAX)
c        PRINT *, '   HESO ='
c        PRINT 6001, (HESO(K),K=1,KMAX)
c        PRINT *, '   TO ='
c        PRINT 6002, (TO(K)-273.16,K=1,KMAX)
c        PRINT *, '   QO ='
c        PRINT 6003, (QO(K),K=1,KMAX)
c        PRINT *, '   QSO ='
c        PRINT 6003, (QESO(K),K=1,KMAX)
c      ENDIF
C
C  LOOK FOR CONVECTIVE CLOUD BASE AS THE LEVEL OF FREE CONVECTION
C
        IF(CNVFLG) THEN
          INDX    = KB
          HKBO = HEO(INDX)
          QKBO = QO(INDX)
          UKBO = UO(INDX)
          VKBO = VO(INDX)
        ENDIF
        FLG    = CNVFLG
        KBCON  = KMAX
!!
      DO K = 1, KM
          if (k .le. kbmax) then
            IF(FLG.AND.K.GT.KB) THEN
              HSBAR   = HESO(k)
              IF(HKBO.GT.HSBAR) THEN
                FLG   = .FALSE.
                KBCON = K
              ENDIF
            ENDIF
          endif
      ENDDO
      
        IF(CNVFLG) THEN
          PBCDIF = -PFLD(KBCON) + PFLD(KB)
          PDOT   = 10.* DOT(KBCON)
          IF(PBCDIF.GT.150.)    CNVFLG = .FALSE.
          IF(KBCON.EQ.KMAX)  CNVFLG = .FALSE.
        ENDIF
!!
      TOTFLG = .TRUE.
        TOTFLG = TOTFLG .AND. (.NOT. CNVFLG)
      IF(TOTFLG) RETURN

 6001 FORMAT(2X,-2P10F12.2)
 6002 FORMAT(2X,10F12.2)
 6003 FORMAT(2X,3P10F12.2)
C
C  DETERMINE ENTRAINMENT RATE BETWEEN KB AND KBCON
C
        alpha = alphas
        if(SLIMSK.eq.1.) alpha = alphal
        IF(CNVFLG) THEN
          IF(KB.EQ.1) THEN
            DZ = .5 * (ZO(KBCON) + ZO(KBCON-1)) - ZO(1)
          ELSE
            DZ = .5 * (ZO(KBCON) + ZO(KBCON-1))
     &         - .5 * (ZO(KB) + ZO(KB-1))
          ENDIF
          IF(KBCON.NE.KB) THEN
            XLAMB = - LOG(ALPHA) / DZ
          ELSE
            XLAMB = 0.
          ENDIF
        ENDIF

C  DETERMINE UPDRAFT MASS FLUX
      DO K = 1, KM
          if (k .le. kmax .and. CNVFLG) then
            ETA(k)  = 1.
            ETAU(k) = 1.
          ENDIF
      ENDDO
      DO K = KM1, 2, -1
          if (k .le. kbmax) then
            IF(CNVFLG.AND.K.LT.KBCON.AND.K.GE.KB) THEN
              DZ        = .5 * (ZO(k+1) - ZO(k-1))
              ETA(k)  = ETA(k+1) * EXP(-XLAMB * DZ)
              ETAU(k) = ETA(k)
            ENDIF
          endif
      ENDDO
        IF(CNVFLG.AND.KB.EQ.1.AND.KBCON.GT.1) THEN
          DZ = .5 * (ZO(2) - ZO(1))
          ETA(1) = ETA(2) * EXP(-XLAMB * DZ)
          ETAU(1) = ETA(1)
        ENDIF
C
C  WORK UP UPDRAFT CLOUD PROPERTIES
C
        IF(CNVFLG) THEN
          INDX         = KB
          HCKO(INDX) = HKBO
          QCKO(INDX) = QKBO
          UCKO(INDX) = UKBO
          VCKO(INDX) = VKBO
          PWAVO     = 0.
        ENDIF
C
C  CLOUD PROPERTY BELOW CLOUD BASE IS MODIFIED BY THE ENTRAINMENT PROCES
C
      DO K = 2, KM1
          if (k .le. kmax-1) then
            IF(CNVFLG.AND.K.GT.KB.AND.K.LE.KBCON) THEN
              FACTOR = ETA(k-1) / ETA(k)
              ONEMF = 1. - FACTOR
              HCKO(k) = FACTOR * HCKO(k-1) + ONEMF *
     &                    .5 * (HEO(k) + HEO(k+1))
              UCKO(k) = FACTOR * UCKO(k-1) + ONEMF *
     &                    .5 * (UO(k) + UO(k+1))
              VCKO(k) = FACTOR * VCKO(k-1) + ONEMF *
     &                    .5 * (VO(k) + VO(k+1))
              DBYO(k) = HCKO(k) - HESO(k)
            ENDIF
            IF(CNVFLG.AND.K.GT.KBCON) THEN
              HCKO(k) = HCKO(k-1)
              UCKO(k) = UCKO(k-1)
              VCKO(k) = VCKO(k-1)
              DBYO(k) = HCKO(k) - HESO(k)
            ENDIF
          endif
      ENDDO
C  DETERMINE CLOUD TOP
        FLG = CNVFLG
        KTCON = 1
      DO K = 2, KM
          if (k .le. kmax) then
            IF(DBYO(k).LT.0..AND.FLG.AND.K.GT.KBCON) THEN
              KTCON = K
              FLG = .FALSE.
            ENDIF
          endif
      ENDDO
        IF(CNVFLG.AND.(PFLD(KBCON) - PFLD(KTCON)).LT.150.)
     &  CNVFLG = .FALSE.
      TOTFLG = .TRUE.
        TOTFLG = TOTFLG .AND. (.NOT. CNVFLG)
      IF(TOTFLG) RETURN
C
C  SEARCH FOR DOWNDRAFT ORIGINATING LEVEL ABOVE THETA-E MINIMUM
C
        HMIN = HEO(KBCON)
        LMIN = KBMAX
        JMIN = KBMAX
        DO K = KBCON, KBMAX
          IF(HEO(k).LT.HMIN.AND.CNVFLG) THEN
            LMIN = K + 1
            HMIN = HEO(k)
          ENDIF
        ENDDO
C
C  Make sure that JMIN is within the cloud
C
        IF(CNVFLG) THEN
          JMIN = MIN(LMIN,KTCON-1)
          XMBMAX = .1
          JMIN = MAX(JMIN,KBCON+1)
        ENDIF
C
C  ENTRAINING CLOUD
C
      do k = 2, km1
          if (k .le. kmax-1) then
            if(CNVFLG.and.k.gt.JMIN.and.k.le.KTCON) THEN
              SUMZ(k) = SUMZ(k-1) + .5 * (ZO(k+1) - ZO(k-1))
              SUMH(k) = SUMH(k-1) + .5 * (ZO(k+1) - ZO(k-1))
     &                  * HEO(k)
            ENDIF
          endif
      enddo
!!
        IF(CNVFLG) THEN
c         call random_number(XKT2)
c         call srand(fhour)
c         XKT2 = rand()
          KT2 = nint(XKT2*float(KTCON-JMIN)-.5)+JMIN+1
!         KT2 = nint(sqrt(XKT2)*float(KTCON-JMIN)-.5) + JMIN + 1
c         KT2 = nint(ranf() *float(KTCON-JMIN)-.5) + JMIN + 1
          tem1 = (HCKO(JMIN) - HESO(KT2))
          tem2 = (SUMZ(KT2) * HESO(KT2) - SUMH(KT2))
          if (abs(tem2) .gt. 0.000001) THEN
            XLAMB = tem1 / tem2
          else
            CNVFLG = .false.
          ENDIF
!         XLAMB = (HCKO(JMIN) - HESO(KT2))
!    &          / (SUMZ(KT2) * HESO(KT2) - SUMH(KT2))
          XLAMB = max(XLAMB,cons_0)
          XLAMB = min(XLAMB,2.3/SUMZ(KT2))
        ENDIF
!!
       DWNFLG  = CNVFLG
       DWNFLG2 = CNVFLG
       IF(CNVFLG) THEN
        if(KT2.ge.KTCON) DWNFLG = .false.
      if(XLAMB.le.1.e-30.or.HCKO(JMIN)-HESO(KT2).le.1.e-30)
     &  DWNFLG = .false.
        do k = JMIN, KT2
          if(DWNFLG.and.HEO(k).gt.HESO(KT2)) DWNFLG=.false.
        enddo
c       IF(CNVFLG.AND.(PFLD(KBCON)-PFLD(KTCON)).GT.PDETRN)
c    &     DWNFLG=.FALSE.
        IF(CNVFLG.AND.(PFLD(KBCON)-PFLD(KTCON)).LT.PDPDWN)
     &     DWNFLG2=.FALSE.
       ENDIF
!!
      DO K = 2, KM1
          if (k .le. kmax-1) then
            IF(DWNFLG.AND.K.GT.JMIN.AND.K.LE.KT2) THEN
              DZ        = .5 * (ZO(k+1) - ZO(k-1))
c             ETA(k)  = ETA(k-1) * EXP( XLAMB * DZ)
c  to simplify matter, we will take the linear approach here
c
              ETA(k)  = ETA(k-1) * (1. + XLAMB * dz)
              ETAU(k) = ETAU(k-1) * (1. + (XLAMB+xlambu) * dz)
            ENDIF
          endif
      ENDDO
!!
      DO K = 2, KM1
          if (k .le. kmax-1) then
c           IF(.NOT.DWNFLG.AND.K.GT.JMIN.AND.K.LE.KT2) THEN
            IF(.NOT.DWNFLG.AND.K.GT.JMIN.AND.K.LE.KTCON) THEN
              DZ        = .5 * (ZO(k+1) - ZO(k-1))
              ETAU(k) = ETAU(k-1) * (1. + xlambu * dz)
            ENDIF
          endif
      ENDDO

        if(DWNFLG) THEN
          KTCON = KT2
        ENDIF
C
C  CLOUD PROPERTY ABOVE CLOUD Base IS MODIFIED BY THE DETRAINMENT PROCESS
C
      DO K = 2, KM1
          if (k .le. kmax-1) then
            IF(CNVFLG.AND.K.GT.KBCON.AND.K.LE.KTCON) THEN
              FACTOR    = ETA(k-1) / ETA(k)
              ONEMF     = 1. - FACTOR
              fuv       = ETAU(k-1) / ETAU(k)
              onemfu    = 1. - fuv
              HCKO(k) = FACTOR * HCKO(k-1) + ONEMF *
     &                    .5 * (HEO(k) + HEO(k+1))
              UCKO(k) = fuv * UCKO(k-1) + ONEMFu *
     &                    .5 * (UO(k) + UO(k+1))
              VCKO(k) = fuv * VCKO(k-1) + ONEMFu *
     &                    .5 * (VO(k) + VO(k+1))
              DBYO(k) = HCKO(k) - HESO(k)
            ENDIF
          endif
      ENDDO
        if(CNVFLG.and.DWNFLG2.and.JMIN.le.KBCON)
     &     THEN
          CNVFLG = .false.
          DWNFLG = .false.
          DWNFLG2 = .false.
        ENDIF
!!
      TOTFLG = .TRUE.
        TOTFLG = TOTFLG .AND. (.NOT. CNVFLG)
      IF(TOTFLG) RETURN
!!
C
C  COMPUTE CLOUD MOISTURE PROPERTY AND PRECIPITATION
C
          AA1 = 0.
          RHBAR = 0.
      DO K = 1, KM
          if (k .le. kmax) then
            IF(CNVFLG.AND.K.GT.KB.AND.K.LT.KTCON) THEN
              DZ = .5 * (ZO(k+1) - ZO(k-1))
              DZ1 = (ZO(k) - ZO(k-1))
              GAMMA = EL2ORC * QESO(k) / (TO(k)**2)
              QRCH = QESO(k)
     &             + GAMMA * DBYO(k) / (HVAP * (1. + GAMMA))
              FACTOR = ETA(k-1) / ETA(k)
              ONEMF = 1. - FACTOR
              QCKO(k) = FACTOR * QCKO(k-1) + ONEMF *
     &                    .5 * (QO(k) + QO(k+1))
              DQ = ETA(k) * QCKO(k) - ETA(k) * QRCH
              RHBAR = RHBAR + QO(k) / QESO(k)
C
C  BELOW LFC CHECK IF THERE IS EXCESS MOISTURE TO RELEASE LATENT HEAT
C
              IF(DQ.GT.0.) THEN
                ETAH = .5 * (ETA(k) + ETA(k-1))
                QLK = DQ / (ETA(k) + ETAH * C0 * DZ)
                AA1 = AA1 - DZ1 * G * QLK
                QC = QLK + QRCH
                PWO(k) = ETAH * C0 * DZ * QLK
                QCKO(k) = QC
                PWAVO = PWAVO + PWO(k)
              ENDIF
            ENDIF
          endif
      ENDDO
        RHBAR = RHBAR / float(KTCON - KB - 1)
c
c  this section is ready for cloud water
c
      if(ncloud.gt.0) THEN
c
c  compute liquid and vapor separation at cloud top
c
        k = KTCON
        IF(CNVFLG) THEN
          GAMMA = EL2ORC * QESO(K) / (TO(K)**2)
          QRCH = QESO(K)
     &         + GAMMA * DBYO(K) / (HVAP * (1. + GAMMA))
          DQ = QCKO(K-1) - QRCH
C
C  CHECK IF THERE IS EXCESS MOISTURE TO RELEASE LATENT HEAT
C
          IF(DQ.GT.0.) THEN
            QLKO_KTCON = dq
            QCKO(K-1) = QRCH
          ENDIF
        ENDIF
      ENDIF
C
C  CALCULATE CLOUD WORK FUNCTION AT T+DT
C
      DO K = 1, KM
          if (k .le. kmax) then
            IF(CNVFLG.AND.K.GT.KBCON.AND.K.LE.KTCON) THEN
              DZ1 = ZO(k) - ZO(k-1)
              GAMMA = EL2ORC * QESO(k-1) / (TO(k-1)**2)
              RFACT =  1. + DELTA * CP * GAMMA
     &                 * TO(k-1) / HVAP
              AA1 = AA1 +
     &                 DZ1 * (G / (CP * TO(k-1)))
     &                 * DBYO(k-1) / (1. + GAMMA)
     &                 * RFACT
              val = 0.
              AA1=AA1+
     &                 DZ1 * G * DELTA *
     &                 MAX(val,(QESO(k-1) - QO(k-1)))
            ENDIF
          endif
      ENDDO
        IF(CNVFLG.AND.AA1.LE.0.) DWNFLG  = .FALSE.
        IF(CNVFLG.AND.AA1.LE.0.) DWNFLG2 = .FALSE.
        IF(CNVFLG.AND.AA1.LE.0.) CNVFLG  = .FALSE.
!!
      TOTFLG = .TRUE.
        TOTFLG = TOTFLG .AND. (.NOT. CNVFLG)
      IF(TOTFLG) RETURN
!!
        IF(CNVFLG) THEN
          VSHEAR = 0.
        ENDIF
      DO K = 1, KM
          if (k .le. kmax) then
            IF(K.GE.KB.AND.K.LE.KTCON.AND.CNVFLG) THEN
              shear= sqrt((UO(k+1)-UO(k)) ** 2
     &                          + (VO(k+1)-VO(k)) ** 2)
              VSHEAR = VSHEAR + SHEAR
            ENDIF
          endif
      ENDDO
        EDT = 0.
        IF(CNVFLG) THEN
          KNUMB = KTCON - KB + 1
          KNUMB = MAX(KNUMB,1)
          VSHEAR = 1.E3 * VSHEAR / (ZO(KTCON)-ZO(KB))
          E1=1.591-.639*VSHEAR
     &       +.0953*(VSHEAR**2)-.00496*(VSHEAR**3)
          EDT=1.-E1
          val =         .9
          EDT = MIN(EDT,val)
          val =         .0
          EDT = MAX(EDT,val)
          EDTO=EDT
          EDTX=EDT
        ENDIF
C  DETERMINE DETRAINMENT RATE BETWEEN 1 AND KBDTR
        KBDTR = KBCON
        beta = betas
        if(SLIMSK.eq.1.) beta = betal
        IF(CNVFLG) THEN
          KBDTR = KBCON
          KBDTR = MAX(KBDTR,1)
          XLAMD = 0.
          IF(KBDTR.GT.1) THEN
            DZ = .5 * ZO(KBDTR) + .5 * ZO(KBDTR-1)
     &         - ZO(1)
            XLAMD =  LOG(BETA) / DZ
          ENDIF
        ENDIF
C  DETERMINE DOWNDRAFT MASS FLUX
      DO K = 1, KM
          IF(k .le. kmax) then
            IF(CNVFLG) THEN
              ETAD(k) = 1.
            ENDIF
            QRCDO(k) = 0.
          endif
      ENDDO
      DO K = KM1, 2, -1
          if (k .le. kbmax) then
            IF(CNVFLG.AND.K.LT.KBDTR) THEN
              DZ        = .5 * (ZO(k+1) - ZO(k-1))
              ETAD(k) = ETAD(k+1) * EXP(XLAMD * DZ)
            ENDIF
          endif
      ENDDO
      K = 1
        IF(CNVFLG.AND.KBDTR.GT.1) THEN
          DZ = .5 * (ZO(2) - ZO(1))
          ETAD(k) = ETAD(k+1) * EXP(XLAMD * DZ)
        ENDIF
C
C--- DOWNDRAFT MOISTURE PROPERTIES
C
        PWEVO = 0.
        FLG = CNVFLG
        IF(CNVFLG) THEN
          JMN = JMIN
          HCDO = HEO(JMN)
          QCDO = QO(JMN)
          QRCDO(JMN) = QESO(JMN)
          UCDO = UO(JMN)
          VCDO = VO(JMN)
        ENDIF
      DO K = KM1, 1, -1
          if (k .le. kmax-1) then
            IF(CNVFLG.AND.K.LT.JMIN) THEN
              DQ = QESO(k)
              DT = TO(k)
              GAMMA      = EL2ORC * DQ / DT**2
              DH         = HCDO - HESO(k)
              QRCDO(k) = DQ+(1./HVAP)*(GAMMA/(1.+GAMMA))*DH
              DETAD      = ETAD(k+1) - ETAD(k)
              PWDO(k)  = ETAD(k+1) * QCDO -
     &                     ETAD(k) * QRCDO(k)
              PWDO(k)  = PWDO(k) - DETAD *
     &                    .5 * (QRCDO(k) + QRCDO(k+1))
              QCDO    = QRCDO(k)
              PWEVO   = PWEVO + PWDO(k)
            ENDIF
          endif
      ENDDO
C
C--- FINAL DOWNDRAFT STRENGTH DEPENDENT ON PRECIP
C--- EFFICIENCY (EDT), NORMALIZED CONDENSATE (PWAV), AND
C--- EVAPORATE (PWEV)
C
        edtmax = edtmaxl
        if(SLIMSK.eq.0.) edtmax = edtmaxs
        IF(DWNFLG2) THEN
          IF(PWEVO.LT.0.) THEN
            EDTO = -EDTO * PWAVO / PWEVO
            EDTO = MIN(EDTO,EDTMAX)
          ELSE
            EDTO = 0.
          ENDIF
        ELSE
          EDTO = 0.
        ENDIF
C
C
C--- DOWNDRAFT CLOUDWORK FUNCTIONS
C
C
      DO K = KM1, 1, -1
          if (k .le. kmax-1) then
            IF(DWNFLG2.AND.K.LT.JMIN) THEN
              GAMMA = EL2ORC * QESO(k+1) / TO(k+1)**2
              DHH=HCDO
              DT=TO(k+1)
              DG=GAMMA
              DH=HESO(k+1)
              DZ=-1.*(ZO(k+1)-ZO(k))
              AA1=AA1+EDTO*DZ*(G/(CP*DT))*((DHH-DH)/(1.+DG))
     &               *(1.+DELTA*CP*DG*DT/HVAP)
              val=0.
              AA1=AA1+EDTO*
     &        DZ*G*DELTA*MAX(val,(QESO(k+1)-QO(k+1)))
            ENDIF
          endif
      ENDDO

        IF(AA1.LE.0.) CNVFLG  = .FALSE.
        IF(AA1.LE.0.) DWNFLG  = .FALSE.
        IF(AA1.LE.0.) DWNFLG2 = .FALSE.
!!
      TOTFLG = .TRUE.
        TOTFLG = TOTFLG .AND. (.NOT. CNVFLG)
      IF(TOTFLG) RETURN
!!
C
C
C--- WHAT WOULD THE CHANGE BE, THAT A CLOUD WITH UNIT MASS
C--- WILL DO TO THE ENVIRONMENT?
C
      DO K = 1, KM
          IF(k .le. kmax .and. CNVFLG) THEN
            DELLAH(k) = 0.
            DELLAQ(k) = 0.
            DELLAU(k) = 0.
            DELLAV(k) = 0.
          ENDIF
      ENDDO

        IF(CNVFLG) THEN
          DP = 1000. * DEL(1)
          DELLAH(1) = EDTO * ETAD(1) * (HCDO
     &                - HEO(1)) * G / DP
          DELLAQ(1) = EDTO * ETAD(1) * (QCDO
     &                - QO(1)) * G / DP
          DELLAU(1) = EDTO * ETAD(1) * (UCDO
     &                - UO(1)) * G / DP
          DELLAV(1) = EDTO * ETAD(1) * (VCDO
     &                - VO(1)) * G / DP
        ENDIF
C
C--- CHANGED DUE TO SUBSIDENCE AND ENTRAINMENT
C
      DO K = 2, KM1
          if (k .le. kmax-1) then
            IF(CNVFLG.AND.K.LT.KTCON) THEN
              AUP = 1.
              IF(K.LE.KB) AUP = 0.
              ADW = 1.
              IF(K.GT.JMIN) ADW = 0.
              DV1= HEO(k)
              DV2 = .5 * (HEO(k) + HEO(k+1))
              DV3= HEO(k-1)
              DV1Q= QO(k)
              DV2Q = .5 * (QO(k) + QO(k+1))
              DV3Q= QO(k-1)
              DV1U= UO(k)
              DV2U = .5 * (UO(k) + UO(k+1))
              DV3U= UO(k-1)
              DV1V= VO(k)
              DV2V = .5 * (VO(k) + VO(k+1))
              DV3V= VO(k-1)
              DP = 1000. * DEL(K)
              DZ = .5 * (ZO(k+1) - ZO(k-1))
              DETA = ETA(k) - ETA(k-1)
              DETAD = ETAD(k) - ETAD(k-1)
              DELLAH(k) = DELLAH(k) +
     &            ((AUP * ETA(k) - ADW * EDTO * ETAD(k)) * DV1
     &        - (AUP * ETA(k-1) - ADW * EDTO * ETAD(k-1))* DV3
     &                    - AUP * DETA * DV2
     &                    + ADW * EDTO * DETAD * HCDO) * G / DP
              DELLAQ(k) = DELLAQ(k) +
     &            ((AUP * ETA(k) - ADW * EDTO * ETAD(k)) * DV1Q
     &        - (AUP * ETA(k-1) - ADW * EDTO * ETAD(k-1))* DV3Q
     &                    - AUP * DETA * DV2Q
     &       +ADW*EDTO*DETAD*.5*(QRCDO(k)+QRCDO(k-1))) * G / DP
              DELLAU(k) = DELLAU(k) +
     &            ((AUP * ETA(k) - ADW * EDTO * ETAD(k)) * DV1U
     &        - (AUP * ETA(k-1) - ADW * EDTO * ETAD(k-1))* DV3U
     &                     - AUP * DETA * DV2U
     &                    + ADW * EDTO * DETAD * UCDO
     &                    ) * G / DP
              DELLAV(k) = DELLAV(k) +
     &            ((AUP * ETA(k) - ADW * EDTO * ETAD(k)) * DV1V
     &        - (AUP * ETA(k-1) - ADW * EDTO * ETAD(k-1))* DV3V
     &                     - AUP * DETA * DV2V
     &                    + ADW * EDTO * DETAD * VCDO
     &                    ) * G / DP
            ENDIF
          endif
      ENDDO
C
C------- CLOUD TOP
C
        IF(CNVFLG) THEN
          INDX = KTCON
          DP = 1000. * DEL(INDX)
          DV1 = HEO(INDX-1)
          DELLAH(INDX) = ETA(INDX-1) *
     &                     (HCKO(INDX-1) - DV1) * G / DP
          DVQ1 = QO(INDX-1)
          DELLAQ(INDX) = ETA(INDX-1) *
     &                     (QCKO(INDX-1) - DVQ1) * G / DP
          DV1U = UO(INDX-1)
          DELLAU(INDX) = ETA(INDX-1) *
     &                     (UCKO(INDX-1) - DV1U) * G / DP
          DV1V = VO(INDX-1)
          DELLAV(INDX) = ETA(INDX-1) *
     &                     (VCKO(INDX-1) - DV1V) * G / DP
c
c  cloud water
c
          DELLAL = ETA(INDX-1) * QLKO_KTCON * g / dp
        ENDIF
C
C------- FINAL CHANGED VARIABLE PER UNIT MASS FLUX
C
      DO K = 1, KM
          if (k .le. kmax) then
            IF(CNVFLG.and.k.gt.KTCON) THEN
              QO(k) = Q1(k)
              TO(k) = T1(k)
              UO(k) = U1(k)
              VO(k) = V1(k)
            ENDIF
            IF(CNVFLG.AND.K.LE.KTCON) THEN
              QO(k) = DELLAQ(k) * MBDT + Q1(k)
              DELLAT = (DELLAH(k) - HVAP * DELLAQ(k)) / CP
              TO(k) = DELLAT * MBDT + T1(k)
              val   =           1.e-10
              QO(k) = max(QO(k), val  )
            ENDIF
          endif
      ENDDO
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C--- THE ABOVE CHANGED ENVIRONMENT IS NOW USED TO CALULATE THE
C--- EFFECT THE ARBITRARY CLOUD (WITH UNIT MASS FLUX)
C--- WOULD HAVE ON THE STABILITY,
C--- WHICH THEN IS USED TO CALCULATE THE REAL MASS FLUX,
C--- NECESSARY TO KEEP THIS CHANGE IN BALANCE WITH THE LARGE-SCALE
C--- DESTABILIZATION.
C
C--- ENVIRONMENTAL CONDITIONS AGAIN, FIRST HEIGHTS
C
      DO K = 1, KM
          IF(k .le. kmax .and. CNVFLG) THEN
!jfe        QESO(k) = 10. * FPVS(TO(k))
!
            QESO(k) = 0.01 * fpvs(TO(K))      ! fpvs is in Pa
!
            QESO(k) = EPS * QESO(k) / (PFLD(k)+EPSM1*QESO(k))
            val       =             1.E-8
            QESO(k) = MAX(QESO(k), val )
            TVO(k)  = TO(k) + DELTA * TO(k) * QO(k)
          ENDIF
      ENDDO
        IF(CNVFLG) THEN
          XAA0 = 0.
          XPWAV = 0.
        ENDIF
C
C--- MOIST STATIC ENERGY
C
      DO K = 1, KM1
          IF(k .le. kmax-1 .and. CNVFLG) THEN
            DZ = .5 * (ZO(k+1) - ZO(k))
            DP = .5 * (PFLD(k+1) - PFLD(k))
cjfe        ES = 10. * FPVS(TO(k+1))
!
            ES = 0.01 * fpvs(TO(K+1))      ! fpvs is in Pa
!
            PPRIME = PFLD(k+1) + EPSM1 * ES
            QS = EPS * ES / PPRIME
            DQSDP = - QS / PPRIME
            DESDT = ES * (FACT1 / TO(k+1) + FACT2 / (TO(k+1)**2))
            DQSDT = QS * PFLD(k+1) * DESDT / (ES * PPRIME)
            GAMMA = EL2ORC * QESO(k+1) / (TO(k+1)**2)
            DT = (G * DZ + HVAP * DQSDP * DP) / (CP * (1. + GAMMA))
            DQ = DQSDT * DT + DQSDP * DP
            TO(k) = TO(k+1) + DT
            QO(k) = QO(k+1) + DQ
            PO(k) = .5 * (PFLD(k) + PFLD(k+1))
          ENDIF
      ENDDO
      DO K = 1, KM1
          IF(k .le. kmax-1 .and. CNVFLG) THEN
cjfe        QESO(k) = 10. * FPVS(TO(k))
!
            QESO(k) = 0.01 * fpvs(TO(K))      ! fpvs is in Pa
!
            QESO(k) = EPS * QESO(k) / (PO(k) + EPSM1 * QESO(k))
            val1      =             1.E-8
            QESO(k) = MAX(QESO(k), val1)
            val2      =           1.e-10
            QO(k)   = max(QO(k), val2 )
c           QO(k)   = MIN(QO(k),QESO(k))
            HEO(k)   = .5 * G * (ZO(k) + ZO(k+1)) +
     &                    CP * TO(k) + HVAP * QO(k)
            HESO(k) = .5 * G * (ZO(k) + ZO(k+1)) +
     &                  CP * TO(k) + HVAP * QESO(k)
          ENDIF
      ENDDO
        k = kmax
        IF(CNVFLG) THEN
          HEO(k) = G * ZO(k) + CP * TO(k) + HVAP * QO(k)
          HESO(k) = G * ZO(k) + CP * TO(k) + HVAP * QESO(k)
c         HEO(k) = MIN(HEO(k),HESO(k))
        ENDIF
        IF(CNVFLG) THEN
          INDX = KB
          XHKB = HEO(INDX)
          XQKB = QO(INDX)
          HCKO(INDX) = XHKB
          QCKO(INDX) = XQKB
        ENDIF
C
C
C**************************** STATIC CONTROL
C
C
C------- MOISTURE AND CLOUD WORK FUNCTIONS
C
      DO K = 2, KM1
          if (k .le. kmax-1) then
C           IF(CNVFLG.AND.K.GT.KB.AND.K.LE.KBCON) THEN
            IF(CNVFLG.AND.K.GT.KB.AND.K.LE.KTCON) THEN
              FACTOR = ETA(k-1) / ETA(k)
              ONEMF = 1. - FACTOR
              HCKO(k) = FACTOR * HCKO(k-1) + ONEMF *
     &                    .5 * (HEO(k) + HEO(k+1))
            ENDIF
C           IF(CNVFLG.AND.K.GT.KBCON) THEN
C             HEO(k) = HEO(k-1)
C           ENDIF
          endif
      ENDDO
      DO K = 2, KM1
          if (k .le. kmax-1) then
            IF(CNVFLG.AND.K.GT.KB.AND.K.LT.KTCON) THEN
              DZ = .5 * (ZO(k+1) - ZO(k-1))
              GAMMA = EL2ORC * QESO(k) / (TO(k)**2)
              XDBY = HCKO(k) - HESO(k)
              val  =          0.
              XDBY = MAX(XDBY,val)
              XQRCH = QESO(k)
     &              + GAMMA * XDBY / (HVAP * (1. + GAMMA))
              FACTOR = ETA(k-1) / ETA(k)
              ONEMF = 1. - FACTOR
              QCKO(k) = FACTOR * QCKO(k-1) + ONEMF *
     &                    .5 * (QO(k) + QO(k+1))
              DQ = ETA(k) * QCKO(k) - ETA(k) * XQRCH
              IF(DQ.GT.0.) THEN
                ETAH = .5 * (ETA(k) + ETA(k-1))
                QLK = DQ / (ETA(k) + ETAH * C0 * DZ)
                XAA0 = XAA0 - (ZO(k) - ZO(k-1)) * G * QLK
                XQC = QLK + XQRCH
                XPW = ETAH * C0 * DZ * QLK
                QCKO(k) = XQC
                XPWAV = XPWAV + XPW
              ENDIF
            ENDIF
c           IF(CNVFLG.AND.K.GT.KBCON.AND.K.LT.KTCON) THEN
            IF(CNVFLG.AND.K.GT.KBCON.AND.K.LE.KTCON) THEN
              DZ1 = ZO(k) - ZO(k-1)
              GAMMA = EL2ORC * QESO(k-1) / (TO(k-1)**2)
              RFACT =  1. + DELTA * CP * GAMMA
     &                 * TO(k-1) / HVAP
              XDBY = HCKO(k-1) - HESO(k-1)
              XAA0 = XAA0
     &                + DZ1 * (G / (CP * TO(k-1)))
     &                * XDBY / (1. + GAMMA)
     &                * RFACT
              val=0.
              XAA0=XAA0+
     &                 DZ1 * G * DELTA *
     &                 MAX(val,(QESO(k-1) - QO(k-1)))
            ENDIF
          endif
      ENDDO
C
C------- DOWNDRAFT CALCULATIONS
C
C
C--- DOWNDRAFT MOISTURE PROPERTIES
C
        XPWEV = 0.
        IF(DWNFLG2) THEN
          JMN = JMIN
          XHCD = HEO(JMN)
          XQCD = QO(JMN)
          QRCD(JMN) = QESO(JMN)
        ENDIF
      DO K = KM1, 1, -1
          if (k .le. kmax-1) then
            IF(DWNFLG2.AND.K.LT.JMIN) THEN
              DQ = QESO(k)
              DT = TO(k)
              GAMMA    = EL2ORC * DQ / DT**2
              DH       = XHCD - HESO(k)
              QRCD(k)=DQ+(1./HVAP)*(GAMMA/(1.+GAMMA))*DH
              DETAD    = ETAD(k+1) - ETAD(k)
              XPWD     = ETAD(k+1) * QRCD(k+1) -
     &                   ETAD(k) * QRCD(k)
              XPWD     = XPWD - DETAD *
     &                 .5 * (QRCD(k) + QRCD(k+1))
              XPWEV = XPWEV + XPWD
            ENDIF
          endif
      ENDDO
C
        edtmax = edtmaxl
        if(SLIMSK.eq.0.) edtmax = edtmaxs
        IF(DWNFLG2) THEN
          IF(XPWEV.GE.0.) THEN
            EDTX = 0.
          ELSE
            EDTX = -EDTX * XPWAV / XPWEV
            EDTX = MIN(EDTX,EDTMAX)
          ENDIF
        ELSE
          EDTX = 0.
        ENDIF
C
C
C
C--- DOWNDRAFT CLOUDWORK FUNCTIONS
C
C
      DO K = KM1, 1, -1
          if (k .le. kmax-1) then
            IF(DWNFLG2.AND.K.LT.JMIN) THEN
              GAMMA = EL2ORC * QESO(k+1) / TO(k+1)**2
              DHH=XHCD
              DT= TO(k+1)
              DG= GAMMA
              DH= HESO(k+1)
              DZ=-1.*(ZO(k+1)-ZO(k))
              XAA0=XAA0+EDTX*DZ*(G/(CP*DT))*((DHH-DH)/(1.+DG))
     &                *(1.+DELTA*CP*DG*DT/HVAP)
              val=0.
              XAA0=XAA0+EDTX*
     &        DZ*G*DELTA*MAX(val,(QESO(k+1)-QO(k+1)))
            ENDIF
          endif
      ENDDO
C
C  CALCULATE CRITICAL CLOUD WORK FUNCTION
C
        ACRT = 0.
        IF(CNVFLG) THEN
C       IF(CNVFLG.AND.SLIMSK.NE.1.) THEN
          IF(PFLD(KTCON).LT.PCRIT(15))THEN
            ACRT=ACRIT(15)*(975.-PFLD(KTCON))
     &              /(975.-PCRIT(15))
          ELSE IF(PFLD(KTCON).GT.PCRIT(1))THEN
            ACRT=ACRIT(1)
          ELSE
            K =  int((850. - PFLD(KTCON))/50.) + 2
            K = MIN(K,15)
            K = MAX(K,2)
            ACRT=ACRIT(K)+(ACRIT(K-1)-ACRIT(K))*
     *           (PFLD(KTCON)-PCRIT(K))/(PCRIT(K-1)-PCRIT(K))
           ENDIF
C        ELSE
C          ACRT = .5 * (PFLD(KBCON) - PFLD(KTCON))
         ENDIF
      
        ACRTFCT = 1.
        IF(CNVFLG) THEN
          if(SLIMSK.eq.1.) THEN
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          ENDIF
c
c  modify critical cloud workfunction by cloud base vertical velocity
c
          IF(PDOT.LE.W4) THEN
            ACRTFCT = (PDOT - W4) / (W3 - W4)
          ELSEIF(PDOT.GE.-W4) THEN
            ACRTFCT = - (PDOT + W4) / (W4 - W3)
          ELSE
            ACRTFCT = 0.
          ENDIF
          val1    =             -1.
          ACRTFCT = MAX(ACRTFCT,val1)
          val2    =             1.
          ACRTFCT = MIN(ACRTFCT,val2)
          ACRTFCT = 1. - ACRTFCT
          DTCONV = DT2 + max((1800. - DT2),cons_0) *
     &                (PDOT - W2) / (W1 - W2)
          DTCONV = max(DTCONV,dtmin)
          DTCONV = min(DTCONV,dtmax)
 
        ENDIF
C
C--- LARGE SCALE FORCING
C
        FLG = CNVFLG
        IF(CNVFLG) THEN
C         F = AA1 / DTCONV
          FLD = (AA1 - ACRT * ACRTFCT) / DTCONV
          IF(FLD.LE.0.) FLG = .FALSE.
        ENDIF
        CNVFLG = FLG
        IF(CNVFLG) THEN
C         XAA0 = MAX(XAA0,0.)
          XK = (XAA0 - AA1) / MBDT
          IF(XK.GE.0.) FLG = .FALSE.
        ENDIF
C
C--- KERNEL, CLOUD BASE MASS FLUX
C
        CNVFLG = FLG
        IF(CNVFLG) THEN
          XMB = -FLD / XK
          XMB = MIN(XMB,XMBMAX)
        ENDIF
      TOTFLG = .TRUE.
        TOTFLG = TOTFLG .AND. (.NOT. CNVFLG)
      IF(TOTFLG) RETURN
c
c  restore t0 and QO to t1 and q1 in case convection stops
c
      do k = 1, km
          if (k .le. kmax) then
            TO(k) = T1(k)
            QO(k) = Q1(k)
!
            QESO(k) = 0.01 * fpvs(T1(K))      ! fpvs is in Pa
!
            QESO(k) = EPS * QESO(k) / (PFLD(k) + EPSM1*QESO(k))
            val     =             1.E-8
            QESO(k) = MAX(QESO(k), val )
          endif
      enddo
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C--- FEEDBACK: SIMPLY THE CHANGES FROM THE CLOUD WITH UNIT MASS FLUX
C---           MULTIPLIED BY  THE MASS FLUX NECESSARY TO KEEP THE
C---           EQUILIBRIUM WITH THE LARGER-SCALE.
C
        DELHBAR = 0.
        DELQBAR = 0.
        DELTBAR = 0.
        QCOND = 0.
      DO K = 1, KM
          if (k .le. kmax) then
            IF(CNVFLG.AND.K.LE.KTCON) THEN
              AUP = 1.
              IF(K.Le.KB) AUP = 0.
              ADW = 1.
              IF(K.GT.JMIN) ADW = 0.
              DELLAT = (DELLAH(k) - HVAP * DELLAQ(k)) / CP
              T1(k) = T1(k) + DELLAT * XMB * DT2
              Q1(k) = Q1(k) + DELLAQ(k) * XMB * DT2
              U1(k) = U1(k) + DELLAU(k) * XMB * DT2
              V1(k) = V1(k) + DELLAV(k) * XMB * DT2
              DP = 1000. * DEL(K)
              DELHBAR = DELHBAR + DELLAH(k)*XMB*DP/G
              DELQBAR = DELQBAR + DELLAQ(k)*XMB*DP/G
              DELTBAR = DELTBAR + DELLAT*XMB*DP/G
            ENDIF
          endif
      ENDDO
      DO K = 1, KM
          if (k .le. kmax) then
            IF(CNVFLG.AND.K.LE.KTCON) THEN
!jfe          QESO(k) = 10. * FPVS(T1(k))
!
              QESO(k) = 0.01 * fpvs(T1(K))      ! fpvs is in Pa
!
              QESO(k) = EPS * QESO(k)/(PFLD(k) + EPSM1*QESO(k))
              val     =             1.E-8
              QESO(k) = MAX(QESO(k), val )
c
c  cloud water
c
              if(ncloud.gt.0.and.cnvflg.and.k.eq.ktcon) then
                tem  = dellal * xmb * dt2
                tem1 = max(0.0, min(1.0, (tcr-t1(k))*tcrf))
                if (ql(k,2) .gt. -999.0) then
                  ql(k,2) = ql(k,2) + tem * tem1            ! ice
                  ql(k,1) = ql(k,1) + tem *(1.0-tem1)       ! water
                else
                  ql(k,1) = ql(k,1) + tem
                endif
                dp = 1000. * del(k)
                dellal = dellal * xmb * dp / g
              endif
!
!             if(ncloud.gt.0.and.CNVFLG.and.k.eq.KTCON) THEN
!               QL(k) = QL(k) + DELLAL * XMB * dt2
!               dp = 1000. * del(k)
!               DELLAL = DELLAL * XMB * dp / g
!             ENDIF
!
            ENDIF
          endif
      ENDDO
        RNTOT = 0.
        DELQEV = 0.
        DELQ2 = 0.
        FLG = CNVFLG
      DO K = KM, 1, -1
          if (k .le. kmax) then
            IF(CNVFLG.AND.K.LE.KTCON) THEN
              AUP = 1.
              IF(K.Le.KB) AUP = 0.
              ADW = 1.
              IF(K.GT.JMIN) ADW = 0.
              rain =  AUP * PWO(k) + ADW * EDTO * PWDO(k)
              RNTOT = RNTOT + rain * XMB * .001 * dt2
            ENDIF
          endif
      ENDDO
      DO K = KM, 1, -1
          if (k .le. kmax) then
            DELTV = 0.
            DELQ = 0.
            QEVAP = 0.
            IF(CNVFLG.AND.K.LE.KTCON) THEN
              AUP = 1.
              IF(K.Le.KB) AUP = 0.
              ADW = 1.
              IF(K.GT.JMIN) ADW = 0.
              rain =  AUP * PWO(k) + ADW * EDTO * PWDO(k)
              RN = RN + rain * XMB * .001 * dt2
            ENDIF
            IF(FLG.AND.K.LE.KTCON) THEN
              evef = EDT * evfact
              if(SLIMSK.eq.1.) evef=EDT * evfactl
c             if(SLIMSK.ne.1.) evef = 0.
              QCOND = EVEF * (Q1(k) - QESO(k))
     &                 / (1. + EL2ORC * QESO(k) / T1(k)**2)
              DP = 1000. * DEL(K)
              IF(RN.GT.0..AND.QCOND.LT.0.) THEN
                QEVAP = -QCOND * (1.-EXP(-.32*SQRT(DT2*RN)))
                QEVAP = MIN(QEVAP, RN*1000.*G/DP)
                DELQ2 = DELQEV + .001 * QEVAP * dp / g
              ENDIF
              if(RN.gt.0..and.QCOND.LT.0..and.
     &           DELQ2.gt.RNTOT) THEN
                QEVAP = 1000.* g * (RNTOT - DELQEV) / dp
                FLG = .false.
              ENDIF
              IF(RN.GT.0..AND.QEVAP.gt.0.) THEN
                Q1(k) = Q1(k) + QEVAP
                T1(k) = T1(k) - ELOCP * QEVAP
                RN = RN - .001 * QEVAP * DP / G
                DELTV = - ELOCP*QEVAP/DT2
                DELQ =  + QEVAP/DT2
                DELQEV = DELQEV + .001*dp*QEVAP/g
              ENDIF
              DELLAQ(k) = DELLAQ(k) + DELQ / XMB
              DELQBAR = DELQBAR + DELQ*DP/G
              DELTBAR = DELTBAR + DELTV*DP/G
            ENDIF
          endif
      ENDDO
C
C  PRECIPITATION RATE CONVERTED TO ACTUAL PRECIP
C  IN UNIT OF M INSTEAD OF KG
C
        IF(CNVFLG) THEN
C
C  IN THE EVENT OF UPPER LEVEL RAIN EVAPORATION AND LOWER LEVEL DOWNDRAF
C    MOISTENING, RN CAN BECOME NEGATIVE, IN THIS CASE, WE BACK OUT OF TH
C    HEATING AND THE MOISTENING
C
          if(RN.lt.0..and..not.FLG) RN = 0.
          IF(RN.LE.0.) THEN
            RN = 0.
          ELSE
            KTOP = KTCON
            KBOT = KBCON
            KUO = 1
            CLDWRK = AA1
          ENDIF
        ENDIF

      DO K = 1, KM
          if (k .le. kmax) then
            IF(CNVFLG.AND.RN.LE.0.) THEN
              T1(k) = TO(k)
              Q1(k) = QO(k)
            ENDIF
          endif
      ENDDO

      RETURN
      END SUBROUTINE 
      
      function fpvs(T)
        use sat_vapor_pres_mod, only : lookup_es
        real, intent(in) :: T
        real :: fpvs
       
        call lookup_es(T,fpvs)
    
        return
      end function fpvs 
      
      end module sascnv_mod 
