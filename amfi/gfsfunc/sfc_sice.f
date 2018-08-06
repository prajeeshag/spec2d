      module seaice_mod

        use mpp_mod, only: mpp_error, FATAL, WARNING

        use sat_vapor_pres_mod, only : compute_qs

        implicit none
        private

        integer, public, parameter :: kmi = 2        ! 2-layer of ice
        real, public, parameter :: tgice = 2.7120e+2     ! temp freezing sea     (K)
        real, public, parameter :: himin = 0.1      ! minimum ice thickness required
        real, parameter :: cimin = 0.      ! minimum ice thickness required

        character(len=512) :: msg

        public :: sfc_sice_drv

      contains

!--------------------------------------------------------------------------------   
        subroutine sfc_sice_drv(imax, pgr, wind, tgrs, qgrs, deltim, 
     &             sfcdlw, sfcnsw, sfcdsw, cd, cdq, prsl, prslki, mask,
     &             hice, tsice, hsnow_sice, fprecip, tice, 
     &             sice_snwdph, qss, sice_snowmt, gflux, cmm, chh, 
     &             zlvl, evap, hflx)
!--------------------------------------------------------------------------------   
            integer, intent(in) :: imax
            real, intent(in), dimension(imax) :: pgr, wind, tgrs, qgrs
            real, intent(in), dimension(imax) :: sfcdlw, sfcnsw, sfcdsw
            real, intent(in), dimension(imax) :: cd, cdq, prsl, prslki
            real, intent(in), dimension(imax) :: fprecip
            logical, intent(in), dimension(imax) :: mask
            real, intent(in) :: deltim
            real, intent(inout), dimension(imax) :: hice, tsice 
            real, intent(inout), dimension(imax) :: hsnow_sice
            real, intent(inout), dimension(kmi,imax) :: tice
            real, intent(out), dimension(imax) :: sice_snwdph, qss
            real, intent(out), dimension(imax) :: sice_snowmt, gflux
            real, intent(out), dimension(imax) :: cmm, chh, zlvl
            real, intent(out), dimension(imax) :: evap, hflx

            integer :: i
            
            do i = 1, imax
                call sfc_sice(pgr(i), wind(i), tgrs(i), qgrs(i), deltim, 
     &             sfcdlw(i), sfcnsw(i), sfcdsw(i), cd(i), cdq(i), 
     &             prsl(i), prslki(i), mask(i), hice(i), tsice(i), 
     &             hsnow_sice(i), fprecip(i), tice(:,i), sice_snwdph(i), 
     &             qss(i), sice_snowmt(i), gflux(i), cmm(i), chh(i), 
     &             zlvl(i), evap(i), hflx(i))

            enddo

        end subroutine sfc_sice_drv
     

!--------------------------------------------------------------------------------   
      subroutine sfc_sice                                               
!  ---  inputs:
     &     ( ps, wind, t1, q1, delt,                            
     &       dlwflx, sfcnsw, sfcdsw,                                    
     &       cm, ch, prsl1, prslki, flag,
!  ---  input/outputs:
     &       hice, tice, sheleg, tprcp, stc,                      
!  ---  outputs:
     &       snwdph, qsurf, snowmt, gflux, cmm, chh,                    
     &       zlvl, evap, hflx                                           
     &     )

      use constants_mod, only : hvap => HLV, sbc => STEFAN,
     &    cp => CP_AIR, RVGAS, RDGAS, GRAV, t0c => KELVIN

      implicit none

!  ---  constant parameters:
      real, parameter :: eps = RDGAS/RVGAS 
      real, parameter :: epsm1 = RDGAS/RVGAS-1.
      real, parameter :: rvrdm1 = RVGAS/RDGAS-1.
      real, parameter :: rd = RDGAS
      
      real, parameter :: cpinv = 1.0/cp
      real, parameter :: hvapi = 1.0/hvap
      real, parameter :: elocp = hvap/cp
      real, parameter :: himax = 8.0      ! maximum ice thickness allowed
      real, parameter :: hsmax = 2.0      ! maximum snow depth allowed
      real, parameter :: timin = 173.0    ! minimum temperature allowed for snow/ice
      real, parameter :: albfw = 0.06     ! albedo for lead
      real, parameter :: dsi   = 1.0/0.33

!  ---  inputs:
      real, intent(in) :: ps, wind,   
     &       t1, q1, dlwflx, sfcnsw, sfcdsw, cm, ch,   
     &       prsl1, prslki, tprcp

      logical, intent(in) :: flag


      real , intent(in) :: delt

!  ---  input/outputs:
      real , intent(inout) :: hice,      
     &       tice, sheleg

      real , dimension(kmi), intent(inout) :: stc

!  ---  outputs:
      real , intent(out) :: snwdph,      
     &       qsurf, snowmt, gflux, cmm, chh, zlvl, evap, hflx

!  ---  locals:
      real :: evapi,
     &       hflxi, sneti, qssi, hfd, hfi,           
     &       focn, snof, hi_save, hs_save, psurf, q0, qs1, rch, rho,    
     &       snowd, theta1, tv1, ps1

      real :: sfcemis=0.98 ! LW emissivity for sea-ice

      real :: t12, t14, tem, stsice(kmi)

      integer :: k
 

        snwdph = 0.; qsurf = 0.; snowmt = 0.; gflux = 0.; cmm = 0.
        chh = 0.; zlvl = 0.; evap = 0.; hflx = 0.

       if (.not.flag) then
           hice = himin; 
           tice = tgice 
           sheleg = 0.
           stc = tgice
           return
       endif

!  --- ...  snow-rain detection

       sheleg = sheleg + 1.e3 * tprcp * delt

!  --- ...  initialize variables. all units are supposedly m.k.s. unless specifie
!           psurf is in pascals, wind is wind speed, theta1 is adiabatic surface
!           temp from level 1, rho is density, qs1 is sat. hum. at level1 and qss
!           is sat. hum. at surface
!           convert slrad to the civilized unit from langley minute-1 k-4

          psurf = 1000.0 * ps
          ps1   = 1000.0 * prsl1

!         dlwflx has been given a negative sign for downward longwave
!         sfcnsw is the net shortwave flux (direction: dn-up)

          q0 = max(q1, 1.0e-8)
          theta1 = t1 * prslki
          tv1 = t1 * (1.0 + rvrdm1*q0)
          rho = ps1 / (rd*tv1)

          call compute_qs(t1,ps1,qs1) 
          qs1 = max(qs1, 1.e-8)
          q0  = min(qs1, q0)

          if (tice < timin) then
              tice = tgice
          endif

          call compute_qs(tice,psurf,qssi)

!  --- ...  snow depth in water equivalent is converted from mm to m unit

            snowd = sheleg * 0.001

!  --- ...  when snow depth is less than 1 mm, a patchy snow is assumed and
!           soil is allowed to interact with the atmosphere.
!           we should eventually move to a linear combination of soil and
!           snow under the condition of patchy snow.

!  --- ...  rcp = rho cp ch v

          rch = rho * cp * ch * wind
          cmm = cm * wind
          chh = rho * ch * wind
          zlvl = -rd * tv1 * log(ps1/psurf) / grav

!  --- ...  sensible and latent heat flux over open water & sea ice

          evapi = elocp * rch * (qssi - q0)

!  --- ...  update sea ice temperature

      do k = 1, kmi
            stsice(k) = stc(k)
      enddo

          sneti = sfcnsw

          t12 = tice * tice
          t14 = t12 * t12

!  --- ...  hfi = net non-solar and upir heat flux @ ice surface

          hfi = -dlwflx + sfcemis*sbc*t14 + evapi           
     &           + rch*(tice - theta1)
          if (hfi/=hfi) then
             print *, "hfi error=", hfi, dlwflx, evapi,
     &                 rch, tice, theta1, t1, prslki
             call mpp_error(fatal,'hfi error')
          endif
          hfd = 4.0*sfcemis*sbc*tice*t12                       
     &           + (1.0 + elocp*eps*hvap*qs1/(rd*t12)) * rch

          t12 = tgice * tgice
          t14 = t12 * t12

!  --- ...  hfw = net heat flux @ water surface (within ice)

          focn = 2.0     ! heat flux from ocean - should be from ocn model
          snof = 0.0     ! snowfall rate - snow accumulates in gbphys

          hice = max( min( hice, himax ), himin )
          snowd = min( snowd, hsmax )

          if (snowd > (2.0*hice)) then
            call mpp_error(WARNING,'Too much snow over '//
     &           'sea-ice, reseting!')
            snowd = 2.0 * hice
          endif

      call ice3lay
!  ---  inputs:                                                         !
!    &     ( im, kmi, flag, hfi, hfd, sneti, focn, delt,          !
!  ---  outputs:                                                        !
!    &       snowd, hice, stsice, tice, snof, snowmt, gflux )           !

          if (tice < timin) then
              write(msg,'(A,F10.5)') 
     &          'surface ice temp is too low; '// 
     &          'reseting it to: ', timin
              call mpp_error(WARNING,msg)
            tice = timin
          endif
        
          do k = 1, kmi
            if (stsice(k) < timin) then
              write(msg,'(A,I3.3,A,F10.5)') 
     &          'layer ', k, ' ice temp is too low; '// 
     &          'reseting it to: ', timin
              call mpp_error(WARNING,msg)
              stsice(k) = timin
            endif
          enddo

      do k = 1, kmi
        stc(k) = min(stsice(k), t0c)
      enddo

!  --- ...  calculate sensible heat flux (& evap over sea ice)

          hflxi = rch * (tice - theta1)
          hflx = hflxi
          evap = evapi

!  --- ...  the rest of the output

          qsurf = q1 + evap / (elocp*rch)

!  --- ...  convert snow depth back to mm of water equivalent

          sheleg = snowd * 1000.0
          snwdph = sheleg * dsi             ! snow depth in mm

          tem     = 1.0 / rho
          hflx = hflx * tem * cpinv
          evap = evap * tem * hvapi

      return

! =================
      contains
! =================


!-----------------------------------
      subroutine ice3lay
!...................................
!  ---  inputs:
!    &     ( im, kmi, fice, flag, hfi, hfd, sneti, focn, delt,          
!  ---  input/outputs:
!    &       snowd, hice, stsice, tice, snof,                           
!  ---  outputs:
!    &       snowmt, gflux                                              
!    &     )

!**************************************************************************
!                                                                         *
!            three-layer sea ice vertical thermodynamics                  *
!                                                                         *
! based on:  m. winton, "a reformulated three-layer sea ice model",       *
! journal of atmospheric and oceanic technology, 2000                     *
!                                                                         *
!                                                                         *
!        -> +---------+ <- tice - diagnostic surface temperature ( <= 0c )*
!       /   |         |                                                   *
!   snowd   |  snow   | <- 0-heat capacity snow layer                     *
!       \   |         |                                                   *
!        => +---------+                                                   *
!       /   |         |                                                   *
!      /    |         | <- t1 - upper 1/2 ice temperature; this layer has *
!     /     |         |         a variable (t/s dependent) heat capacity  *
!   hice    |...ice...|                                                   *
!     \     |         |                                                   *
!      \    |         | <- t2 - lower 1/2 ice temp. (fixed heat capacity) *
!       \   |         |                                                   *
!        -> +---------+ <- base of ice fixed at seawater freezing temp.   *
!                                                                         *
!  =====================  defination of variables  =====================  !
!                                                                         !
!  inputs:                                                         size   !
!     im, kmi  - integer, horiz dimension and num of ice layers      1    !
!     fice     - real, sea-ice concentration                         im   !
!     hfi      - real, net non-solar and heat flux @ surface(w/m^2)  im   !
!     hfd      - real, heat flux derivatice @ sfc (w/m^2/deg-c)      im   !
!     sneti    - real, net solar incoming at top  (w/m^2)            im   !
!     focn     - real, heat flux from ocean    (w/m^2)               im   !
!     delt     - real, timestep                (sec)                 1    !
!                                                                         !
!  input/outputs:                                                         !
!     snowd    - real, surface pressure                              im   !
!     hice     - real, sea-ice thickness                             im   !
!     stsice   - real, temp @ midpt of ice levels  (deg c)          im,kmi!     
!     tice     - real, surface temperature     (deg c)               im   !
!     snof     - real, snowfall rate           (m/sec)               im   !
!                                                                         !
!  outputs:                                                               !
!     snowmt   - real, snow melt during delt   (m)                   im   !
!     gflux    - real, conductive heat flux    (w/m^2)               im   !
!                                                                         !
!  locals:                                                                !
!     hdi      - real, ice-water interface     (m)                   im   !
!     hsni     - real, snow-ice                (m)                   im   !
!                                                                         !
! ======================================================================= !
!

!  ---  constant parameters: (properties of ice, snow, and seawater)
      real , parameter :: ds   = 330.0    ! snow (ov sea ice) density (kg/m^3)
      real , parameter :: dw   =1000.0    ! fresh water density  (kg/m^3)
      real , parameter :: t0c  =273.15    ! freezing temp of fresh ice (k)
      real , parameter :: ks   = 0.31     ! conductivity of snow   (w/mk)
      real , parameter :: i0   = 0.3      ! ice surface penetrating solar fraction
      real , parameter :: ki   = 2.03     ! conductivity of ice  (w/mk)
      real , parameter :: di   = 917.0    ! density of ice   (kg/m^3)
      real , parameter :: ci   = 2054.0   ! heat capacity of fresh ice (j/kg/k)
      real , parameter :: li   = 3.34e5   ! latent heat of fusion (j/kg-ice)
      real , parameter :: si   = 1.0      ! salinity of sea ice
      real , parameter :: mu   = 0.054    ! relates freezing temp to salinity
      real , parameter :: tfi  = -mu*si   ! sea ice freezing temp = -mu*salinity
      real , parameter :: tfw  = -1.8     ! tfw - seawater freezing temp (c)
      real , parameter :: tfi0 = tfi-0.0001
      real , parameter :: dici = di*ci 
      real , parameter :: dili = di*li 
      real , parameter :: dsli = ds*li 
      real , parameter :: ki4  = ki*4.0

!  ---  locals:
      real :: hdi, hsni, a, b, ip,      
     &      a1, b1, c1, a10, b10, k12, k32, h1, h2, dh, f1, tsf,        
     &      tmelt, bmelt

      real :: dt2, dt4, dt6, tmpt

      integer :: i
!
!===> ...  begin here
!
      dt2 = 2.0 * delt
      dt4 = 4.0 * delt
      dt6 = 6.0 * delt

          tmpt=tice
          snowd = snowd * dw / ds
          hdi = (ds*snowd + di*hice) / dw

          if (hice < hdi) then
            snowd = snowd + hice - hdi
            hsni  = (hdi - hice) * ds / di
            hice  = hice + hsni
          endif

          snof = snof * dw / ds
          tice = tice - t0c
          stsice(1) = min(stsice(1)-t0c, tfi0)     ! degc
          stsice(2) = min(stsice(2)-t0c, tfi0)     ! degc

          ip = i0 * sneti         ! ip +v (in winton ip=-i0*sneti as sol -v)
          if (snowd > 0.0) then
            tsf = 0.0
            ip  = 0.0
          else
            tsf = tfi
            ip  = i0 * sneti      ! ip +v here (in winton ip=-i0*sneti)
          endif
          tice = min(tice, tsf)

!  --- ...  compute ice temperature

          b = hfd
          a = hfi - sneti + ip - tice*b  ! +v sol input here
          k12 = ki4*ks / (ks*hice + ki4*snowd)
          k32 = 2.0 * ki / hice

          a10 = dici*hice/dt2 + k32*(dt4*k32 + dici*hice)
     &           / (dt6*k32 + dici*hice)
          b10 = -di*hice * (ci*stsice(1) + li*tfi/stsice(1))  
     &           / dt2 - ip                                          
     &           - k32*(dt4*k32*tfw + dici*hice*stsice(2))   
     &           / (dt6*k32 + dici*hice)

          a1 = a10 + k12*b / (k12 + b)
          b1 = b10 + a*k12 / (k12 + b)
          c1 = dili*tfi / dt2*hice

          stsice(1) = -(sqrt(b1*b1 - 4.0*a1*c1)           
     &                + b1)/(2.0*a1)
          tice = (k12*stsice(1) - a) / (k12 + b)

          if (tice > tsf) then
            a1 = a10 + k12
            b1 = b10 - k12*tsf
            stsice(1) = -(sqrt(b1*b1 - 4.0*a1*c1)         
     &                  + b1) / (2.0*a1)
            tice = tsf
            tmelt = (k12*(stsice(1) - tsf)                   
     &               - (a + b*tsf)) * delt
          else
            tmelt = 0.0
            snowd = snowd + snof*delt
          endif

          stsice(2) = (dt2*k32*(stsice(1) + 2.0*tfw)             
     &                + dici*hice*stsice(2))                       
     &                / (dt6*k32 + dici*hice)

          bmelt = (focn + ki4*(stsice(2) - tfw)/hice) * delt

!  --- ...  resize the ice ...

          h1 = 0.5 * hice
          h2 = 0.5 * hice

!  --- ...  top ...

          if (tmelt <= snowd*dsli) then
            snowmt = tmelt / dsli
            snowd  = snowd - tmelt/dsli
          else
            snowmt = snowd
            h1 = h1 - (tmelt - snowd*dsli)                  
     &            / (di * (ci - li/stsice(1)) * (tfi - stsice(1)))
            snowd = 0.0
          endif

!  --- ...  and bottom

          if (bmelt < 0.0) then
            dh = -bmelt / (dili + dici*(tfi - tfw))
            stsice(2) = (h2*stsice(2) + dh*tfw)               
     &                  / (h2 + dh)
            h2 = h2 + dh
          else
            h2 = h2 - bmelt / (dili + dici*(tfi - stsice(2)))
          endif

!  --- ...  if ice remains, even up 2 layers, else, pass negative energy back in snow

          hice = h1 + h2

          if (hice > 0.0) then
            if (h1 > 0.5*hice) then
              f1 = 1.0 - 2.0*h2 / hice
              stsice(2) = f1                                       
     &                    * (stsice(1) + li*tfi/(ci*stsice(1)))     
     &                    + (1.0 - f1)*stsice(2)

              if (stsice(2) > tfi) then
                hice = hice - h2*ci*(stsice(2) - tfi)        
     &                  / (li*delt)
                stsice(2) = tfi
              endif
            else
              f1 = 2.0 * h1 / hice
              stsice(1) = f1                                       
     &                    * (stsice(1) + li*tfi/(ci*stsice(1)))     
     &                    + (1.0 - f1)*stsice(2)
              stsice(1) = (stsice(1) - sqrt(stsice(1)*stsice(1) 
     &                    - 4.0*tfi*li/ci)) / 2.0
            endif

            k12 = ki4*ks / (ks*hice + ki4*snowd)
            gflux = k12 * (stsice(1) - tice)
          else
            snowd = snowd + (h1*(ci*(stsice(1) - tfi)        
     &               - li*(1.0 - tfi/stsice(1)))                      
     &               + h2*(ci*(stsice(2) - tfi) - li)) / li

            hice = max(0.0, snowd*ds/di)
            snowd = 0.0
            stsice(1) = tfw
            stsice(2) = tfw
            gflux = 0.0
          endif   ! end if_hice_block

          gflux = gflux
          snowmt = snowmt * ds / dw
          snowd = snowd * ds / dw
          tice = tice + t0c
          stsice(1) = stsice(1) + t0c
          stsice(2) = stsice(2) + t0c

      return
      end subroutine ice3lay

      end subroutine sfc_sice

      end module seaice_mod
