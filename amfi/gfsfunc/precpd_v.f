       module precpd_mod
       implicit none
       
       contains 
       subroutine precpd (km,dt,del,prsl,ps,q,cwm,t,rn,sn, u00k)
!
!
!     ******************************************************************
!     *                                                                *
!     *           subroutine for precipitation processes               *
!     *           from suspended cloud water/ice                       *
!     *                                                                *
!     ******************************************************************
!     *                                                                *
!     *  originally created by  q. zhao                jan. 1995       *
!     *                         -------                                *    
!     *  modified and rewritten by shrinivas moorthi   oct. 1998       *
!     *                            -----------------                   *
!     *  and                       hua-lu pan                          *
!     *                            ----------                          *
!     *                                                                *
!     *  references:                                                   *
!     *                                                                *
!     *  zhao and carr (1997), monthly weather review (august)         *
!     *  sundqvist et al., (1989) monthly weather review. (august)     *
!     *                                                                *
!     ******************************************************************
!
!     in this code vertical indexing runs from surface to top of the
!     model
!
!     argument list:
!     --------------
!       im         : inner dimension over which calculation is made
!       ix         : maximum inner dimension
!       km         : number of vertical levels
!       dt         : time step in seconds
!       del(km)    : pressure layer thickness (bottom to top)
!       prsl(km)   : pressure values for model layers (bottom to top)
!       ps     : surface pressure (centibars)
!       q(km)   : specific humidity (updated in the code)
!       cwm(km) : condensate mixing ratio (updated in the code)
!       t(km)   : temperature       (updated in the code)
!       rn     : precipitation over one time-step dt (m/dt)
!       sr     : index (=-1 snow, =0 rain/snow, =1 rain)
!       tcw    : vertically integrated liquid water (kg/m**2)
!       cll(km) : cloud cover
!
      use funcphys_mod , only : fpvs
      use physcons_mod, grav => con_g, hvap => con_hvap, 
     &              hfus => con_hfus
     &,             ttp => con_ttp, cp => con_cp
     &,             eps => con_eps, epsm1 => con_epsm1
      implicit none
!     include 'constant.h'
!
      real  g,      h1,    h2,   h1000
     &,                     h1000g, d00,   d125, d5
     &,                     elwv,   eliv,  row
     &,                     epsq,   dldt,  tm10, eliw
     &,                     rcp,    rrow
       parameter (g=grav,         h1=1.e0,     h2=2.e0,     h1000=1000.0
     &,           h1000g=h1000/g, d00=0.e0,    d125=.125e0, d5=0.5e0
     &,           elwv=hvap,      eliv=hvap+hfus,   row=1.e3
     &,           epsq=2.e-12,    dldt=2274.e0,tm10=ttp-10.0
     &,           eliw=eliv-elwv, rcp=h1/cp,   rrow=h1/row)
!
      real, parameter :: cons_0=0.0,     cons_p01=0.01
     &,                                  cons_20=20.0
     &,                                  cons_m30=-30.0, cons_50=50.0
!
      integer km, lat
      real  q(km),   t(km),    cwm(km)
     &,                                 del(km),  prsl(km)
!    &,                     cll(km), del(km),  prsl(km)
     &,                     ps,     rn,      sr
     &,                     tcw,    dt, sn
!
!
      real  err,      ers,     precrl
     &,                     precsl,   precrl1, precsl1
     &,                     rq,       condt
     &,                     conde,    rconde,  tmt0
     &,                     wmin(km),  wmink,   pres
     &,                     wmini(km), ccr,     cclim(km)
     &,                     tt,       qq,      ww
     &,                     wfix(km),     u00k(km), es
     &,                     zaodt
!
      integer iw(km), iwl,     iwl1
!
       logical comput
       logical lprnt
!
      real  ke,   rdt,  us, cclimit, climit, cws, csm1
     &,                     crs1, crs2, cr, aa2,     dtcp,   c00, cmr
     &,                     tem,  c1,   c2, wwn
!    &,                     tem,  c1,   c2, u00b,    u00t,   wwn
     &,                     precrk, precsk, pres1,   qk,     qw,  qi
     &,                     ai,     bi, qint, fiw, wws, cwmk, expf
     &,                     psaut, psaci, amaxcm, tem1, tem2
     &,                     tmt0k, tmt15, psm1, psm2, ppr, pswi
     &,                     rprs,  erk,   pps, sid, rid, amaxps
     &,                     praut, pracw, fi, qc, amaxrq, rqkll
      integer i, k, ihpr, n
!
!-----------------------preliminaries ---------------------------------

      lprnt = .false.
!
      rdt     = h1 / dt
      ke      = 2.0e-5  ! commented on 09/10/99
      us      = h1
      cclimit = 1.0e-3
      climit  = 1.0e-20
      cws     = 0.025

      zaodt   = 800.0 * rdt

      csm1    = 5.0000e-8   * zaodt
      crs1    = 5.00000e-6  * zaodt
      crs2    = 6.66600e-10 * zaodt
      cr      = 5.0e-4      * zaodt
      aa2     = 1.25e-3     * zaodt

      ke      = ke * sqrt(rdt)
      dtcp    = dt * rcp
      c00 = 1.0e-4 * dt          !05/09/2000
      cmr = 1.0 / 3.0e-4
      c1  = 300.0
      c2  = 0.5
!
!
!--------calculate c0 and cmr using lc at previous step-----------------
!
      do k=1,km
          tem   = (prsl(k)*0.01)
          iw(k)    = 0.0
          wmin(k)  = 1.0e-5 * tem
          wmini(k) = 1.0e-5 * tem       ! testing for ras
      enddo
        iwl1    = 0
        precrl1 = d00
        precsl1 = d00
        comput  = .false.
        rn      = d00
        sn      = d00
        sr      = d00
        ccr     = d00
!------------select columns where rain can be produced--------------
      do k=1, km-1
          tem = min(wmin(k), wmini(k))
          if (cwm(k) .gt. tem) comput = .true.
      enddo

      ihpr = 0
      if (comput) ihpr =  1
!***********************************************************************
!-----------------begining of precipitation calculation-----------------
!***********************************************************************

      do k=km,1,-1
        if (ihpr>0) then
          precrl = precrl1
          precsl = precsl1
          err    = d00
          ers    = d00
          iwl    = 0

          tt     = t(k)
          qq     = q(k)
          ww     = cwm(k)
          wmink  = wmin(k)
          pres   = h1000 * prsl(k)

          precrk = max(cons_0,    precrl1)
          precsk = max(cons_0,    precsl1)
          wwn    = max(ww, climit)

          if (wwn .gt. climit .or. (precrk+precsk) .gt. d00) then
            comput = .true.
          else
            comput = .false.
          endif
        endif

        if (ihpr>0) then
          if (comput) then
            conde  = (h1000*dt/g) * del(k)
            condt  = conde * rdt
            rconde = h1 / conde
            qk        = max(epsq,  qq)
            tmt0   = tt - 273.16
            wwn       = max(ww, climit)


!  the global qsat computation is done in pa
            pres1   = pres 
            qw      = min(pres1, fpvs(tt))
            qw      = eps * qw / (pres1 + epsm1 * qw)
            qw      = max(qw,epsq)
            qi   = qw
            qint = qw


!-------------------ice-water id number iw------------------------------
            if(tmt0.lt.-15.) then
               fi = qk - u00k(k)*qi
               if(fi.gt.d00.or.wwn.gt.climit) then
                  iwl = 1
               else
                  iwl = 0
               endif
            elseif (tmt0.ge.0.) then
               iwl = 0
            else
              iwl = 0
              if(iwl1.eq.1.and.wwn.gt.climit) iwl=1
            endif

!----------------the satuation specific humidity------------------------
            fiw   = float(iwl)
            qc    = (h1-fiw)*qint + fiw*qi
!----------------the relative humidity----------------------------------
            if(qc .le. 1.0e-10) then
               rq = d00
            else
               rq = qk / qc
            endif
!----------------cloud cover ratio ccr----------------------------------
            if(rq.lt.u00k(k)) then
                   ccr=d00
            elseif(rq.ge.us) then
                   ccr=us
            else
                 rqkll=min(us,rq)
                 ccr= h1-sqrt((us-rqkll)/(us-u00k(k)))
            endif
          endif
        endif
!-------------------ice-water id number iwl------------------------------
!
!---   precipitation production --  auto conversion and accretion
!
        
        if (ihpr>0) then
          if (comput .and. ccr .gt. 0.0) then
            wws    = ww
            cwmk   = max(cons_0, wws)
            if (iwl .eq. 1) then                 !  ice phase
               amaxcm = max(cons_0, cwmk - wmini(k))
               expf      = dt * exp(0.025*tmt0)
               psaut     = min(cwmk, 4.0e-4*expf*amaxcm)
               ww     = ww - psaut
               cwmk      = max(cons_0, ww)
               psaci     = min(cwmk, aa2*expf*precsl1*cwmk)
               ww     = ww - psaci
               precsl = precsl + (wws - ww) * condt
               pswi = (wws - ww) * condt * rconde !Prajeesh: pswi is water -> ice
               tt = tt + dtcp * (eliw*pswi) 

            else                                    !  liquid water

!          for using sundqvist precip formulation of rain

               amaxcm    = cwmk
               tem1      = precsl1 + precrl1
               tem2      = min(max(cons_0, 268.0-tt), cons_20)
               tem       = (1.0+c1*sqrt(tem1*rdt)) * (1+c2*sqrt(tem2))
               tem2      = amaxcm * cmr * tem / max(ccr,cons_p01)
               tem2      = min(cons_50, tem2*tem2)
               praut     = c00  * tem * amaxcm * (1.0-exp(-tem2))
               praut     = min(praut, cwmk)
               ww     = ww - praut

!          below is for zhao's precip formulation (water)

               precrl = precrl + (wws - ww) * condt
            endif
          endif
        endif

!-----evaporation of precipitation-------------------------
!**** err & ers positive--->evaporation-- negtive--->condensation

        
        if (ihpr>0) then
          if (comput) then
            qk     = max(epsq,  qq)
            tmt0k  = max(cons_m30, tmt0)
            precrk = max(cons_0,    precrl)
            precsk = max(cons_0,    precsl)
            amaxrq = max(cons_0,    u00k(k)-rq) * conde
!----------------------------------------------------------------------
! increase the evaporation for strong/light prec
!----------------------------------------------------------------------
            ppr    = ke * amaxrq * sqrt(precrk)
            if (tmt0 .ge. 0.) then
              pps = 0.
            else
              pps = (crs1+crs2*tmt0k) * amaxrq * precsk / u00k(k)
            end if
!---------------correct if over-evapo./cond. occurs--------------------
            erk=precrk+precsk
            if(rq.ge.1.0e-10)  erk = amaxrq * qk * rdt / rq
            if (ppr+pps .gt. abs(erk)) then
               rprs   = erk / (precrk+precsk)
               ppr    = precrk * rprs
               pps    = precsk * rprs
            endif
            ppr       = min(ppr, precrk)
            pps       = min(pps, precsk)
            err    = ppr * rconde
            ers    = pps * rconde
            precrl = precrl - ppr
            precsl = precsl - pps
          endif
        endif
!--------------------melting of the snow--------------------------------

        if (ihpr>0) then
          if (comput) then
            if (tmt0 .gt. 0.) then
               amaxps = max(cons_0,    precsl)
               psm1   = csm1 * tmt0 * tmt0 * amaxps
               psm2   = cws * cr * max(cons_0, ww) * amaxps
               ppr    = (psm1 + psm2) * conde
               if (ppr .gt. amaxps) then
                 ppr  = amaxps
                 psm1 = amaxps * rconde
               endif
               precrl = precrl + ppr
               precsl = precsl - ppr
               psm1 = ppr * rconde !Prajeesh: psm1 is ice -> water
            else
               psm1 = d00
            endif

!---------------update t and q------------------------------------------
            tt = tt - dtcp * (elwv*err+eliv*ers+eliw*psm1) !Prajeesh
!            tt = tt - dtcp * (elwv*err+elwv*ers)
            qq = qq + dt * (err+ers)
          endif
        endif

        if (ihpr>0) then
          iwl1    = iwl
          precrl1 = max(cons_0, precrl)
          precsl1 = max(cons_0, precsl)
          t(k)     = tt
          q(k)     = qq
          cwm(k)   = ww
          iw(k)    = iwl
        endif

!  move water from vapor to liquid should the liquid amount be negative

          if (cwm(k) < 0.) then
            tem      = q(k) + cwm(k)
            if (tem >= 0.0) then
              q(k)   = tem
              t(k)   = t(k) - elwv * rcp * cwm(k)
              cwm(k) = 0.
            elseif (q(k) > 0.0) then
              cwm(k) = tem
              t(k)   = t(k) + elwv * rcp * q(k)
              q(k)   = 0.0
            endif
          endif
      enddo                               ! k loop ends here!
!**********************************************************************
!-----------------------end of precipitation processes-----------------
!**********************************************************************

      if (ihpr>0) then
        rn = precrl1 * rrow  ! precip at surface
        sn = precsl1 * rrow  ! precip at surface

!----sr=1 if sfc prec is rain ; ----sr=-1 if sfc prec is snow
!----sr=0 for both of them or no sfc prec

        rid = 0.
        sid = 0.
        if (precrl1 .ge. 1.e-13) rid = 1.
        if (precsl1 .ge. 1.e-13) sid = -1.
        sr = rid + sid  ! sr=1 --> rain, sr=-1 -->snow, sr=0 -->both
      endif

      return
      end subroutine
      end module 
