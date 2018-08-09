      module gscond_mod
      implicit none
      contains
      subroutine gscond (km,dt,prsl,ps,q,cwm,t
     &,                  tp, qp, psp, tp1, qp1, psp1, u)
!
!     ******************************************************************
!     *                                                                *
!     *  subroutine for grid-scale condensation & evaporation          *
!     *  for the mrf model at ncep.                                    *
!     *                                                                *
!     ******************************************************************
!     *                                                                *
!     *  created by:   q.  zhao         jan. 1995                      *
!     *  modified by:  h.-l. pan        sep. 1998                      *
!     *  modified by:  s. moorthi       aug. 1999, 2000                *
!     *                                                                *
!     *  references:                                                   *
!     *                                                                *
!     ******************************************************************
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      use funcphys_mod, only : fpvs
      use physcons_mod, psat => con_psat, hvap => con_hvap, 
     &              grav => con_g
     &,             hfus => con_hfus, ttp => con_ttp, rd => con_rd
     &,             cp => con_cp, eps => con_eps, epsm1 => con_epsm1
     &,             rv => con_rv
      implicit none

      real  g,    h1,   h2, h1000
     &,                     d00,  d125,  d5,   elwv, eliv
     &,                     epsq, tm10,  eliw, arcp
     &,                     a1,   r,     cpr,  rcpr, rcp
      parameter (h1=1.e0,       h2=2.e0,               h1000=1000.0
     &,          d00=0.e0,      d125=.125e0,           d5=0.5e0
     &,          a1=psat  
     &,          elwv=hvap,     eliv=hvap+hfus, g=grav
     &,          epsq=2.e-12,   tm10=ttp-10.,   r=rd
     &,          cpr=cp*r,      rcpr=h1/cpr,    rcp=h1/cp)
!
      real, parameter :: cons_0=0.0, cons_m15=-15.0
!
      integer km
      real  q(km),    t(km),    cwm(km)
     &,                     prsl(km), ps, dt
     &,                     tp(km),   qp(km),   psp
     &,                     tp1(km),  qp1(km),  psp1
!
      real   qi, qint, u(km), ccrik, e0
     &,                      cond,   rdt, us, cclimit, climit
     &,                      u00b,   u00t, tmt0, tmt15, qik, cwmik
     &,                      ai, bi, qw, u00ik, tik, pres, pp0, fi 
     &,                      at, aq, ap, fiw, elv, qc, rqik
     &,                      rqikk, tx1, tx2, tx3, es, qs
     &,                      tsq, delq, condi, cone0, us00, ccrik1
     &,                      aa, ab, ac, ad, ae, af, ag
     &,                      el2orc, albycp, vprs           
      integer iw(km), k, iwik
      logical lprnt
!
!-----------------prepare constants for later uses-----------------
      lprnt = .false.
      el2orc = hvap*hvap / (rv*cp)
      albycp = hvap / cp

      rdt     = h1/dt
      us      = h1
      cclimit = 1.0e-3
      climit  = 1.0e-20

        iw(km) = d00

!  check for first time step

      if (tp(1) .lt. 1.) then
        do k = 1, km
            tp(k) = t(k)
            qp(k) = max(q(k),epsq)
            tp1(k) = t(k)
            qp1(k) = max(q(k),epsq)
        enddo
          psp  = ps
          psp1 = ps
      endif
c
c*************************************************************
c*******begining of  grid-scale condensation/evap. loop*******
c*************************************************************
c
!     do k = km-1,2,-1
      do k = km,1,-1
!       vprs(:) = 0.001 * fpvs(t(:,k))       ! fpvs in pa
c-----------------------------------------------------------------------
c------------------qw, qi and qint--------------------------------------
          tmt0  = t(k)-273.16                                                
          tmt15 = min(tmt0,cons_m15)                                            
          qik   = max(q(k),epsq)
          cwmik = max(cwm(k),climit)
!  the global qsat computation is done in cb
          pres    = prsl(k)
          qw      = min(pres, 0.001 * fpvs(t(k)))

          qw      = eps * qw / (pres + epsm1 * qw)
          qw      = max(qw,epsq)
          qi   = qw
          qint = qw

c-------------------ice-water id number iw------------------------------
          if(tmt0.lt.-15.0) then
            u00ik = u(k)
            fi    = qik - u00ik*qi    
            if(fi.gt.d00.or.cwmik.gt.climit) then                    
               iw(k) = 1                                                   
            else                                                           
              iw(k) = 0                                                   
            end if                                                         
          end if

          if(tmt0.ge.0.0) then
            iw(k) = 0
          end if

          if (tmt0 .lt. 0.0 .and. tmt0 .ge. -15.0) then
            iw(k) = 0
            if (k .lt. km) then
            if (iw(k+1) .eq. 1 .and. cwmik .gt. climit) iw(k) = 1
            endif
          end if
c--------------condensation and evaporation of cloud--------------------
c------------------------at, aq and dp/dt-------------------------------
          qik   = max(q(k),epsq)
          cwmik = max(cwm(k),climit)
          iwik  = iw(k)
          u00ik = u(k)
          tik   = t(k)
          pres  = prsl(k)   * h1000
          pp0   = (pres / ps) * psp
          at    = (tik-tp(k)) * rdt
          aq    = (qik-qp(k)) * rdt
          ap    = (pres-pp0)    * rdt
c----------------the satuation specific humidity------------------------
          fiw   = float(iwik)
          elv   = (h1-fiw)*elwv    + fiw*eliv
          qc    = (h1-fiw)*qint + fiw*qi
!     if (lprnt) print *,' qc=',qc,' qint=',qint,' qi=',qi
c----------------the relative humidity----------------------------------
          if(qc.le.1.0e-10) then
            rqik=d00 
          else
            rqik = qik/qc
          endif
c----------------cloud cover ratio ccrik--------------------------------
          if (rqik .lt. u00ik) then
             ccrik = d00
          elseif(rqik.ge.us) then
             ccrik = us
          else
             rqikk  = min(us,rqik)
             ccrik = h1-sqrt((us-rqikk)/(us-u00ik))
          endif
c-----------correct ccr if it is too small in large cwm regions--------
c         if(ccrik.ge.0.01.and.ccrik.le.0.2.and
c    &          .cwmik.ge.0.2e-3) then
c          ccrik=min(1.0,cwmik*1.0e3)
c         end if
c----------------------------------------------------------------------
!   if no cloud exists then evaporate any existing cloud condensate
c----------------evaporation of cloud water-----------------------------
          e0 = d00
          if (ccrik.le.cclimit.and.cwmik.gt.climit)  then 
!
!   first iteration - increment halved
!
            tx1 = tik
            tx3 = qik
!
            es   = min(pres, fpvs(tx1))
            qs   = u00ik * eps * es / (pres + epsm1*es)
            tsq  = tx1 * tx1
            delq = 0.5 * (qs - tx3) * tsq / (tsq + el2orc * qs)
!
            tx2   = delq
            tx1   = tx1 - delq * albycp
            tx3   = tx3 + delq
!
!   second iteration
!
            es   = min(pres, fpvs(tx1))
            qs   = u00ik * eps * es / (pres + epsm1*es)
            tsq  = tx1 * tx1
            delq = (qs - tx3) * tsq / (tsq + el2orc * qs)
!
            tx2  = tx2 + delq
            tx1  = tx1 - delq * albycp
            tx3  = tx3 + delq
!
!   third iteration
!
            es   = min(pres, fpvs(tx1))
            qs   = u00ik * eps * es / (pres + epsm1*es)
            tsq  = tx1 * tx1
            delq = (qs - tx3) * tsq / (tsq + el2orc * qs)
            tx2  = tx2 + delq
!
            e0   = max(tx2*rdt, cons_0)
            e0 = min(cwmik*rdt,   e0)
            e0 = max(cons_0,e0)
          end if
!   if cloud cover > 0.2 condense water vapor in to cloud condensate
c-----------the eqs. for cond. has been reorganized to reduce cpu------
          cond = d00
!         if (ccrik .gt. 0.20 .and. qc .gt. epsq) then
          if (ccrik .gt. cclimit .and. qc .gt. epsq) then
             us00   = us  - u00ik 
             ccrik1 = 1.0 - ccrik
             aa     = eps*elv*pres*qik
             ab     = ccrik*ccrik1*qc*us00
             ac     = ab + 0.5*cwmik
             ad     = ab * ccrik1
             ae     = cpr*tik*tik
             af     = ae * pres
             ag     = aa * elv
             ai     = cp * aa
             cond   = (ac-ad)*(af*aq-ai*at+ae*qik*ap)/(ac*(af+ag))
c-----------check & correct if over condensation occurs-----------------
             condi  = (qik   -u00ik   *qc*1.0)*rdt
             cond   = min(cond, condi)
             cond = max(cond, d00)
c-------------------update of t, q and cwm------------------------------
          end if
          cone0    = (cond-e0) * dt
          cwm(k) = cwm(k) + cone0
          t(k)   = t(k)   + elwv*rcp*cone0 !change elv to elwv for energy conservation: Prajeesh 
          q(k)   = q(k)   - cone0
      enddo                                    ! end of k-loop!

c*********************************************************************
c****************end of the condensation/evaporation loop*************
c*********************************************************************
c----------------store t, q, ps for next time step
      do k = 1, km
          tp(k)  = tp1(k)
          qp(k)  = qp1(k)
!
          tp1(k) = t(k)
          qp1(k) = max(q(k),epsq)
      enddo
        psp  = psp1
        psp1 = ps
c-----------------------------------------------------------------------
      return
      end subroutine
      end module
