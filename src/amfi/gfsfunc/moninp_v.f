        module moninp_mod
        implicit none
        contains
      SUBROUTINE MONINP(KM,ntrac,DV,DU,TAU,RTG,
     &     U1,V1,T1,Q1,
     &     PSK,RBSOIL,FM,FH,QSS,HEAT,EVAP,STRESS,SPD1,KPBL,
     &     PRSI,DEL,PRSL,PRSLK,PHII,PHIL,DELTIM,
     &     DUSFC,DVSFC,DTSFC,DQSFC,HPBL,HGAMT,HGAMQ)

      USE constants_mod, only : grav, RD => RDGAS, CP => CP_AIR
     &,             HVAP => HLV, RVGAS

      implicit none
      real, parameter :: FV = RVGAS/RD - 1.
      integer KM, ntrac, KPBL
      real DELTIM
      real DV(KM),     DU(KM),
     &                     TAU(KM),    RTG(KM,ntrac),
     &                     U1(KM),     V1(KM),
     &                     T1(KM),     Q1(KM,ntrac),
     &                     PSK,       RBSOIL,
     &                     FM,        FH,
     &                     QSS,
     &                                    SPD1,
     &                     PRSI(KM+1), DEL(KM),
     &                     PRSL(KM),   PRSLK(KM),
     &                     PHII(KM+1), PHIL(KM),
     &                     DUSFC,
     &                     dvsfc,     dtsfc,
     &                     DQSFC,     HPBL,
     &                     HGAMT,     hgamq
      integer iprt,is,iun,k,kk,kmpbl,lond
      real :: dspheat
      real evap,  heat,    phih,
     &                     phim,  rbdn,    rbup,
     &                     the1,  stress,  beta,
     &                     the1v, thekv,   thermal,
     &                     thesv, ustar,   wscale
      real RDZT(KM-1),
     &                     ZI(KM+1),     ZL(KM),
     &                     DKU(KM-1),    DKT(KM-1),
     &                     AL(KM-1),     AD(KM),
     &                     AU(KM-1),     A1(KM),
     &                     A2(KM*ntrac), THETA(KM)
      logical              pblflg,   sfcflg, stable
      real aphi16,  aphi5,  bet1,   bvf2,
     &                     cfac,    conq,   cont,   conw,
     &                     conwrc,  dk,     dkmax,  dkmin,
     &                     dq1,     dsdz2,  dsdzq,  dsdzt,
     &                     dsig,    dt,     dthe1,  dtodsd,
     &                     dtodsu,  dw2,    dw2min, g,
     &                     gamcrq,  gamcrt, gocp,   gor, gravi,
     &                     hol,     pfac,   prmax,  prmin, prinv,
     &                     prnum,   qmin,   qtend,  rbcr,
     &                     rbint,   rdt,    rdz,
     &                     ri,      rimin,  rl2,    rlam,  rlamun,
     &                     rone,   rzero,   sfcfrac,
     &                     sflux,   shr2,   spdk2,  sri,
     &                     tem,     ti,     ttend,  tvd,
     &                     tvu,     utend,  vk,     vk2,
     &                     vpert,   vtend,  xkzo(km),   zfac,
     &                     zfmin,   zk,     tem1
      real :: xkzm=1.
      parameter (gravi=1.0/grav)
      PARAMETER(g=grav)
      PARAMETER(GOR=G/RD,GOCP=G/CP)
      PARAMETER(CONT=1000.*CP/G,CONQ=1000.*HVAP/G,CONW=1000./G)
      PARAMETER(RLAM=30.0,VK=0.4,VK2=VK*VK,PRMIN=1.0,PRMAX=4.)
      PARAMETER(DW2MIN=0.0001,DKMIN=0.0,DKMAX=1000.,RIMIN=-100.)
      PARAMETER(RBCR=0.25,CFAC=7.8,PFAC=2.0,SFCFRAC=0.1)
      PARAMETER(QMIN=1.E-8,         ZFMIN=1.E-8,APHI5=5.,APHI16=16.)
      PARAMETER(GAMCRT=3.,GAMCRQ=0., RLAMUN=150.0)
      PARAMETER(IUN=84)

      dspheat = 0.
      DT    = 2. * DELTIM
      RDT   = 1. / DT
      KMPBL = KM / 2
      do k=1,km
          zi(k) = phii(k) * gravi
          zl(k) = phil(k) * gravi
      enddo
      do k=1,kmpbl
          theta(k) = t1(k) * psk / prslk(k)
      enddo
      DO K = 1,KM-1
          RDZT(K) = 1.0 / (ZL(K+1) - ZL(K))
      ENDDO
         DUSFC = 0.
         DVSFC = 0.
         DTSFC = 0.
         DQSFC = 0.
         HGAMT = 0.
         HGAMQ = 0.
         WSCALE = 0.
         KPBL = 1
         HPBL = ZI(2)
         PBLFLG = .TRUE.
         SFCFLG = .TRUE.
         IF(RBSOIL.GT.0.0) SFCFLG = .FALSE.

         BETA  = DT / (zi(2)-zi(1))
         USTAR = SQRT(STRESS)
         THE1    = THETA(1)
         THE1V   = THE1*(1.+FV*MAX(Q1(1,1),QMIN))
         THERMAL = THE1V


         STABLE = .FALSE.
         RBUP = RBSOIL

      DO K = 2, KMPBL
          IF(.NOT.STABLE) THEN
             RBDN   = RBUP
             THEKV  = THETA(k)*(1.+FV*MAX(Q1(k,1),QMIN))
             SPDK2     = MAX((U1(k)**2+V1(k)**2),1.)
             RBUP   = (THEKV-THE1V)*(G*ZL(k)/THE1V)/SPDK2
             KPBL   = K
             STABLE = RBUP.GT.RBCR
          ENDIF
      ENDDO
         K = KPBL
         IF(RBDN.GE.RBCR) THEN
            RBINT = 0.
         ELSEIF(RBUP.LE.RBCR) THEN
            RBINT = 1.
         ELSE
            RBINT = (RBCR-RBDN)/(RBUP-RBDN)
         ENDIF
         HPBL = ZL(K-1) + RBINT*(ZL(K)-ZL(K-1))
         IF(HPBL.LT.ZI(KPBL)) KPBL = KPBL - 1
           HOL = MAX(RBSOIL*FM*FM/FH,RIMIN)
           IF(SFCFLG) THEN
              HOL = MIN(HOL,-ZFMIN)
           ELSE
              HOL = MAX(HOL,ZFMIN)
           ENDIF
           HOL = HOL*HPBL/ZL(1)*SFCFRAC
           IF(SFCFLG) THEN
              TEM  = 1.0 / (1. - APHI16*HOL)
              PHIH = SQRT(TEM)
              PHIM = SQRT(PHIH)
           ELSE
              PHIM = (1.+APHI5*HOL)
              PHIH = PHIM
           ENDIF
           WSCALE = USTAR/PHIM
           WSCALE = MIN(WSCALE,USTAR*APHI16)
           WSCALE = MAX(WSCALE,USTAR/APHI5)
         SFLUX  = HEAT + EVAP*FV*THE1
         IF(SFCFLG.AND.SFLUX.GT.0.0) THEN
           HGAMT   = MIN(CFAC*HEAT/WSCALE,GAMCRT)
           HGAMQ   = MIN(CFAC*EVAP/WSCALE,GAMCRQ)
           VPERT      = HGAMT + FV*THE1*HGAMQ
           VPERT      = MIN(VPERT,GAMCRT)
           THERMAL = THERMAL + MAX(VPERT,0.)
           HGAMT   = MAX(HGAMT,0.0)
           HGAMQ   = MAX(HGAMQ,0.0)
         ELSE
           PBLFLG = .FALSE.
         ENDIF
         IF(PBLFLG) THEN
            KPBL = 1
            HPBL = ZI(2)
         ENDIF
         IF(PBLFLG) THEN
            STABLE = .FALSE.
            RBUP = RBSOIL
         ENDIF
      DO K = 2, KMPBL
          IF(.NOT.STABLE.AND.PBLFLG) THEN
            RBDN   = RBUP
            THEKV  = THETA(k)*(1.+FV*MAX(Q1(k,1),QMIN))
            SPDK2     = MAX((U1(k)**2+V1(k)**2),1.)
            RBUP   = (THEKV-THERMAL)*(G*ZL(k)/THE1V)/SPDK2
            KPBL   = K
            STABLE = RBUP.GT.RBCR
          ENDIF
      ENDDO
         IF(PBLFLG) THEN
            K = KPBL
            IF(RBDN.GE.RBCR) THEN
               RBINT = 0.
            ELSEIF(RBUP.LE.RBCR) THEN
               RBINT = 1.
            ELSE
               RBINT = (RBCR-RBDN)/(RBUP-RBDN)
            ENDIF
            HPBL = ZL(K-1) + RBINT*(ZL(k)-ZL(K-1))
            IF(HPBL.LT.ZI(KPBL)) KPBL = KPBL - 1
            IF(KPBL.LE.1) PBLFLG = .FALSE.
         ENDIF

      DO K = 1,KM-1
          tem1      = 1.0 - prsi(k+1) / prsi(1)
          tem1      = tem1 * tem1 * 10.0
          xkzo(k) = xkzm * min(1.0, exp(-tem1))
      ENDDO
      DO K = 1, KMPBL
            IF(KPBL.GT.K) THEN
               PRINV = 1.0 / (PHIH/PHIM+CFAC*VK*.1)
               PRINV = MIN(PRINV,PRMAX)
               PRINV = MAX(PRINV,PRMIN)
               ZFAC = MAX((1.-(ZI(K+1)-ZL(1))/
     1                (HPBL-ZL(1))), ZFMIN)
               DKU(k) = XKZO(k) + WSCALE*VK*ZI(K+1)
     1                         * ZFAC**PFAC
               DKT(k) = DKU(k)*PRINV
               DKU(k) = MIN(DKU(k),DKMAX)
               DKU(k) = MAX(DKU(k),DKMIN)
               DKT(k) = MIN(DKT(k),DKMAX)
               DKT(k) = MAX(DKT(k),DKMIN)
            ENDIF
      ENDDO
      DO K = 1, KM-1
            IF(K.GE.KPBL) THEN
               TI   = 2.0 / (T1(k)+T1(K+1))
               RDZ  = RDZT(K)

               DW2  = (U1(k)-U1(K+1))**2
     &                      + (V1(k)-V1(K+1))**2
               SHR2 = MAX(DW2,DW2MIN)*RDZ*RDZ
               TVD  = T1(k)*(1.+FV*MAX(Q1(k,1),QMIN))
               TVU  = T1(K+1)*(1.+FV*MAX(Q1(K+1,1),QMIN))
               BVF2 = G*(GOCP+RDZ*(TVU-TVD)) * TI
               RI   = MAX(BVF2/SHR2,RIMIN)
               ZK   = VK*ZI(K+1)
               IF(RI.LT.0.) THEN ! UNSTABLE REGIME
                  RL2      = ZK*RLAMUN/(RLAMUN+ZK)
                  DK       = RL2*RL2*SQRT(SHR2)
                  SRI      = SQRT(-RI)
                  DKU(k) = XKZO(k) + DK*(1+8.*(-RI)/(1+1.746*SRI))
                  DKT(k) = XKZO(k) + DK*(1+8.*(-RI)/(1+1.286*SRI))
               ELSE             ! STABLE REGIME
                  RL2       = ZK*RLAM/(RLAM+ZK)
                  DK        = RL2*RL2*SQRT(SHR2)
                  DKT(k)  = XKZO(k) + DK/(1+5.*RI)**2
                  PRNUM     = 1.0 + 2.1*RI
                  PRNUM     = MIN(PRNUM,PRMAX)
                  DKU(k)  = (DKT(k)-XKZO(k))*PRNUM + XKZO(k)
               ENDIF
               DKU(k) = MIN(DKU(k),DKMAX)
               DKU(k) = MAX(DKU(k),DKMIN)
               DKT(k) = MIN(DKT(k),DKMAX)
               DKT(k) = MAX(DKT(k),DKMIN)
            ENDIF
      ENDDO
         AD(1) = 1.
         A1(1) = T1(1)   + BETA * HEAT
         A2(1) = Q1(1,1) + BETA * EVAP

      if(ntrac.ge.2) then
        do k = 2, ntrac
          is = (k-1) * km
            A2(1+is) = Q1(1,k)
        enddo
      endif
      DO K = 1,KM-1
          DTODSD = DT/DEL(K)
          DTODSU = DT/DEL(K+1)
          DSIG   = PRSL(K)-PRSL(K+1)
          RDZ    = RDZT(K)
          tem1   = DSIG * DKT(k) * RDZ
          IF(PBLFLG.AND.K.LT.KPBL) THEN
             tem   = 1.0 / HPBL
             DSDZT = tem1 * (GOCP-HGAMT*tem)
             DSDZQ = tem1 * (-HGAMQ*tem)
             A2(k)   = A2(k)+DTODSD*DSDZQ
             A2(k+1) = Q1(k+1,1)-DTODSU*DSDZQ
          ELSE
             DSDZT = tem1 * GOCP
             A2(k+1) = Q1(k+1,1)
          ENDIF
          DSDZ2     = tem1 * RDZ
          AU(k)   = -DTODSD*DSDZ2
          AL(k)   = -DTODSU*DSDZ2
          AD(k)   = AD(k)-AU(k)
          AD(k+1) = 1.-AL(k)
          A1(k)   = A1(k)+DTODSD*DSDZT
          A1(k+1) = T1(k+1)-DTODSU*DSDZT
      ENDDO

      if(ntrac.ge.2) then
        do kk = 2, ntrac
          is = (kk-1) * km
          do k = 1, km - 1
              A2(k+1+is) = Q1(k+1,kk)
          enddo
        enddo
      endif
      CALL TRIDIN(KM,ntrac,AL,AD,AU,A1,A2,AU,A1,A2)
      DO  K = 1,KM
            TTEND      = (A1(k)-T1(k))*RDT
            QTEND      = (A2(k)-Q1(k,1))*RDT
            TAU(k)   = TAU(k)+TTEND
            RTG(k,1) = RTG(k,1)+QTEND
            DTSFC   = DTSFC+CONT*DEL(K)*TTEND
            DQSFC   = DQSFC+CONQ*DEL(K)*QTEND
      ENDDO
      if(ntrac.ge.2) then
        do kk = 2, ntrac
          is = (kk-1) * km
          do k = 1, km 
              QTEND = (A2(K+is)-Q1(K,kk))*RDT
              RTG(K,kk) = RTG(K,kk)+QTEND
          enddo
        enddo
      endif
         AD(1) = 1.0 + BETA * STRESS / SPD1
         A1(1) = U1(1)
         A2(1) = V1(1)
      DO K = 1,KM-1
          DTODSD    = DT/DEL(K)
          DTODSU    = DT/DEL(K+1)
          DSIG      = PRSL(K)-PRSL(K+1)
          RDZ       = RDZT(K)
          DSDZ2     = DSIG*DKU(k)*RDZ*RDZ
          AU(k)   = -DTODSD*DSDZ2
          AL(k)   = -DTODSU*DSDZ2
          AD(k)   = AD(k)-AU(k)
          AD(k+1) = 1.-AL(k)
          A1(k+1) = U1(k+1)
          A2(k+1) = V1(k+1)
      ENDDO
      CALL TRIDI2(KM,AL,AD,AU,A1,A2,AU,A1,A2)
      DO K = 1,KM
            CONWRC = CONW
            UTEND = (A1(k)-U1(k))*RDT
            VTEND = (A2(k)-V1(k))*RDT
            DU(k)  = DU(k)+UTEND
            DV(k)  = DV(k)+VTEND
            DUSFC = DUSFC+CONWRC*DEL(K)*UTEND
            DVSFC = DVSFC+CONWRC*DEL(K)*VTEND
            TTEND=0.5*RDT*((U1(K)**2-A1(K)**2) +
     &           (V1(K)**2-A2(K)**2))/CP
            TAU(K)=TAU(K)+TTEND
            dspheat=dspheat+CONT*DEL(k)*TTEND
      ENDDO
      RETURN
      END SUBROUTINE



      SUBROUTINE TRIDI2(N,CL,CM,CU,R1,R2,AU,A1,A2)
      implicit none
      integer             k,n
      real fk

      real CL(2:N),CM(N),CU(N-1),R1(N),R2(N),
     &          AU(N-1),A1(N),A2(N)


        FK      = 1./CM(1)
        AU(1) = FK*CU(1)
        A1(1) = FK*R1(1)
        A2(1) = FK*R2(1)


      DO K=2,N-1
          FK      = 1./(CM(K)-CL(K)*AU(K-1))
          AU(K) = FK*CU(K)
          A1(K) = FK*(R1(K)-CL(K)*A1(K-1))
          A2(K) = FK*(R2(K)-CL(K)*A2(K-1))
      ENDDO

        FK      = 1./(CM(N)-CL(N)*AU(N-1))
        A1(N) = FK*(R1(N)-CL(N)*A1(N-1))
        A2(N) = FK*(R2(N)-CL(N)*A2(N-1))

      DO K=N-1,1,-1
          A1(K) = A1(K)-AU(K)*A1(K+1)
          A2(K) = A2(K)-AU(K)*A2(K+1)
      ENDDO

      RETURN
      END SUBROUTINE TRIDI2


      SUBROUTINE TRIDIN(N,nt,CL,CM,CU,R1,R2,AU,A1,A2)
      implicit none
      integer             is,k,kk,n,nt
      real fk
      real CL(2:N), CM(N), CU(N-1),
     &                     R1(N),   R2(N*nt),
     &                     AU(N-1), A1(N), A2(N*nt),
     &                     FKK(2:N-1)
        FK   = 1./CM(1)
        AU(1) = FK*CU(1)
        A1(1) = FK*R1(1)

      do k = 1, nt
        is = (k-1) * n
        a2(1+is) = fk * r2(1+is)
      enddo

      DO K=2,N-1
          FKK(K) = 1./(CM(K)-CL(K)*AU(K-1))
          AU(K)  = FKK(K)*CU(K)
          A1(K)  = FKK(K)*(R1(K)-CL(K)*A1(K-1))
      ENDDO

      do kk = 1, nt
        is = (kk-1) * n
        DO K=2,N-1
            A2(K+is) = FKK(K)*(R2(K+is)-CL(K)*A2(K+is-1))
        ENDDO
      ENDDO

        FK   = 1./(CM(N)-CL(N)*AU(N-1))
        A1(N) = FK*(R1(N)-CL(N)*A1(N-1))

      do k = 1, nt
        is = (k-1) * n
          A2(N+is) = FK*(R2(N+is)-CL(N)*A2(N+is-1))
      enddo

      DO K=N-1,1,-1
          A1(K) = A1(K) - AU(K)*A1(K+1)
      ENDDO

      do kk = 1, nt
        is = (kk-1) * n
        DO K=n-1,1,-1
            A2(K+is) = A2(K+is) - AU(K)*A2(K+is+1)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE TRIDIN

      end module moninp_mod
