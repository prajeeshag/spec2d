      module seaice_mod

        use mpp_mod, only: mpp_error, fatal

        use sat_vapor_pres_mod, only : compute_qs

        implicit none
        private

        integer, public, parameter :: kmi = 2        ! 2-layer of ice
        real, public, parameter :: tgice = 2.7120e+2     ! temp freezing sea     (K)
        real, public, parameter :: himin = 0.1      ! minimum ice thickness required

        integer, parameter :: im = 1

        public :: sfc_sice_drv

      contains

!--------------------------------------------------------------------------------   
        subroutine sfc_sice_drv(imax, pgr, wind, tgrs, qgrs, deltim, 
     &             sfcdlw, sfcnsw, sfcdsw, cd, cdq, prsl, prslki, 
     &             hice, fice, tsice, hsnow_sice, fprecip, tice, 
     &             sice_snwdph, qss, sice_snowmt, gflux, cmm, chh, 
     &             zlvl, evap, hflx)
!--------------------------------------------------------------------------------   
            integer, intent(in) :: imax
            real, intent(in), dimension(im) :: pgr, wind, tgrs, qgrs
            real, intent(in), dimension(im) :: sfcdlw, sfcnsw, sfcdsw
            real, intent(in), dimension(im) :: cd, cdq, prsl, prslki
            real, intent(in), dimension(im) :: fprecip
            real, intent(in) :: deltim
            real, intent(inout), dimension(im) :: hice, fice, tsice 
            real, intent(inout), dimension(im) :: hsnow_sice
            real, intent(inout), dimension(kmi,im) :: tice
            real, intent(out), dimension(im) :: sice_snwdph, qss
            real, intent(out), dimension(im) :: sice_snowmt, gflux
            real, intent(out), dimension(im) :: cmm, chh, zlvl
            real, intent(out), dimension(im) :: evap, hflx

            integer :: i
            
            do i = 1, imax
                call sfc_sice(pgr(i), wind(i), tgrs(i), qgrs(i), deltim, 
     &             sfcdlw(i), sfcnsw(i), sfcdsw(i), cd(i), cdq(i), 
     &             prsl(i), prslki(i), hice(i), fice(i), tsice(i), 
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
     &       cm, ch, prsl1, prslki,                                     
!  ---  input/outputs:
     &       hice, fice, tice, sheleg, tprcp, stc,                      
!  ---  outputs:
     &       snwdph, qsurf, snowmt, gflux, cmm, chh,                    
     &       zlvl, evap, hflx                                           
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_sice                                                      !
!       inputs:                                                         !
!          ( ps, wind, t1, q1, delt,                          !
!            dlwflx, sfcnsw, sfcdsw, srflag,                   !
!            cm, ch, prsl1, prslki,                 !
!            flag_iter, mom4ice, lsm,                                   !
!       input/outputs:                                                  !
!            hice, fice, tice, sheleg, tprcp, stc, ep,           !
!       outputs:                                                        !
!            snwdph, qsurf, snowmt, gflux, cmm, chh,                    !
!            zlvl, evap, hflx )                                         !
!                                                                       !
!  subprogram called:  ice3lay.                                         !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im,  - integer, horiz dimension and num of soil layers   1    !
!     ps       - real, surface pressure                            im   !
!     u1, v1   - real, u/v component of surface layer wind         im   !
!     t1       - real, surface layer mean temperature ( k )        im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     delt     - real, time interval (second)                      1    !
!     sfcemis  - real, sfc lw emissivity ( fraction )              im   !
!     dlwflx   - real, total sky sfc downward lw flux ( w/m**2 )   im   !
!     sfcnsw   - real, total sky sfc netsw flx into ground(w/m**2) im   !
!     sfcdsw   - real, total sky sfc downward sw flux ( w/m**2 )   im   !
!     srflag   - real, snow/rain flag for precipitation            im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, surface layer mean pressure                 im   !
!     prslki   - real,                                             im   !
!     flag_iter- logical,                                          im   !
!     mom4ice  - logical,                                          im   !
!     lsm      - integer, flag for land surface model scheme       1    !
!                =0: use osu scheme; =1: use noah scheme                !
!                                                                       !
!  input/outputs:                                                       !
!     hice     - real, sea-ice thickness                           im   !
!     fice     - real, sea-ice concentration                       im   !
!     tice     - real, sea-ice surface temperature                 im   !
!     sheleg   - real, snow depth (water equiv)                    im   !
!     tprcp    - real, total precipitation                         im   !
!     stc      - real, soil temp (k)                              im, !
!     ep       - real, potential evaporation                       im   !
!                                                                       !
!  outputs:                                                             !
!     snwdph   - real, water equivalent snow depth (mm)            im   !
!     qsurf    - real, specific humidity at sfc                    im   !
!     snowmt   - real, snow melt (m)                               im   !
!     gflux    - real, soil heat flux (w/m**2)                     im   !
!     cmm      - real,                                             im   !
!     chh      - real,                                             im   !
!     zlvl     - real,                                             im   !
!     evap     - real, evaperation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!                                                                       !
! ===================================================================== !
!
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
      real , dimension(im), intent(in) :: ps, wind,   
     &       t1, q1, dlwflx, sfcnsw, sfcdsw, cm, ch,   
     &       prsl1, prslki, tprcp


      real , intent(in) :: delt

!  ---  input/outputs:
      real , dimension(im), intent(inout) :: hice,      
     &       fice, tice, sheleg

      real , dimension(im,kmi), intent(inout) :: stc

!  ---  outputs:
      real , dimension(im), intent(out) :: snwdph,      
     &       qsurf, snowmt, gflux, cmm, chh, zlvl, evap, hflx

!  ---  locals:
      real , dimension(im) :: ffw, evapi, evapw,        
     &       hflxi, hflxw, sneti, snetw, qssi, hfd, hfi, hfw,           
     &       focn, snof, hi_save, hs_save, psurf, q0, qs1, rch, rho,    
     &       snowd, theta1, tv1, ps1

      real :: sfcemis=0.98 ! LW emissivity for sea-ice

      real :: cimin, t12, t14, tem, stsice(im,kmi)

      integer :: i, k
 
      logical :: flag(im)
!
!  --- ...  set minimum ice concentration required

      cimin = 0.001          ! gfs only

!  --- ...  set flag for sea-ice
      flag(:) =.false.
      do i = 1, im
         flag(i) = fice(i)>=cimin
      enddo

!  --- ...  snow-rain detection

        do i = 1, im
          if (flag(i)) then
              sheleg(i) = sheleg(i) + 1.e3*tprcp(i)
          endif
        enddo

!  --- ...  initialize variables. all units are supposedly m.k.s. unless specifie
!           psurf is in pascals, wind is wind speed, theta1 is adiabatic surface
!           temp from level 1, rho is density, qs1 is sat. hum. at level1 and qss
!           is sat. hum. at surface
!           convert slrad to the civilized unit from langley minute-1 k-4

      do i = 1, im
        if (flag(i)) then
          psurf(i) = 1000.0 * ps(i)
          ps1(i)   = 1000.0 * prsl1(i)

!         dlwflx has been given a negative sign for downward longwave
!         sfcnsw is the net shortwave flux (direction: dn-up)

          q0(i) = max(q1(i), 1.0e-8)
          theta1(i) = t1(i) * prslki(i)
          tv1(i) = t1(i) * (1.0 + rvrdm1*q0(i))
          rho(i) = ps1(i) / (rd*tv1(i))

          call compute_qs(t1(i),ps1(i),qs1(i)) 
          qs1(i) = max(qs1(i), 1.e-8)
          q0 (i) = min(qs1(i), q0(i))

          ffw(i) = 1.0 - fice(i)
          if (fice(i) < cimin) then
            print *,'warning: ice fraction is low:', fice(i)
            fice(i) = cimin
            ffw (i) = 1.0 - fice(i)
            tice(i) = tgice
            print *,'fix ice fraction: reset it to:', fice(i)
          endif

          if (tice(i) < timin) then
              tice(i) = tgice
          endif

          call compute_qs(tice(i),psurf(i),qssi(i))

!  --- ...  snow depth in water equivalent is converted from mm to m unit

            snowd(i) = sheleg(i) * 0.001
!         flagsnw(i) = .false.

!  --- ...  when snow depth is less than 1 mm, a patchy snow is assumed and
!           soil is allowed to interact with the atmosphere.
!           we should eventually move to a linear combination of soil and
!           snow under the condition of patchy snow.
        endif
      enddo

!  --- ...  rcp = rho cp ch v

      do i = 1, im
        if (flag(i)) then
          rch(i) = rho(i) * cp * ch(i) * wind(i)
          cmm(i) = cm(i) * wind(i)
          chh(i) = rho(i) * ch(i) * wind(i)
          zlvl(i) = -rd * tv1(i) * log(ps1(i)/psurf(i)) / grav
        endif
      enddo

!  --- ...  sensible and latent heat flux over open water & sea ice

      do i = 1, im
        if (flag(i)) then
          evapi(i) = elocp * rch(i) * (qssi(i) - q0(i))
        endif
      enddo

!  --- ...  update sea ice temperature

      do k = 1, kmi
        do i = 1, im
          if (flag(i)) then
            stsice(i,k) = stc(i,k)
          endif
        enddo
      enddo

      do i = 1, im
        if (flag(i)) then
          snetw(i) = sfcdsw(i) * (1.0 - albfw)
          snetw(i) = min(3.0*sfcnsw(i)/(1.0+2.0*ffw(i)), snetw(i))
          sneti(i) = (sfcnsw(i) - ffw(i)*snetw(i)) / fice(i)

          t12 = tice(i) * tice(i)
          t14 = t12 * t12

!  --- ...  hfi = net non-solar and upir heat flux @ ice surface

          hfi(i) = -dlwflx(i) + sfcemis*sbc*t14 + evapi(i)           
     &           + rch(i)*(tice(i) - theta1(i))
          if (hfi(i)/=hfi(i)) then
             print *, "hfi error=", hfi(i), dlwflx(i), evapi(i),
     &                 rch(i), tice(i), theta1(i), t1(i), prslki(i)
             call mpp_error(fatal,'hfi error')
          endif
          hfd(i) = 4.0*sfcemis*sbc*tice(i)*t12                       
     &           + (1.0 + elocp*eps*hvap*qs1(i)/(rd*t12)) * rch(i)

          t12 = tgice * tgice
          t14 = t12 * t12

!  --- ...  hfw = net heat flux @ water surface (within ice)

!         hfw(i) = -dlwflx(i) + sfcemis(i)*sbc*t14 + evapw(i)           
!    &           + rch(i)*(tgice - theta1(i)) - snetw(i)

          focn(i) = 2.0     ! heat flux from ocean - should be from ocn model
          snof(i) = 0.0     ! snowfall rate - snow accumulates in gbphys

          hice(i) = max( min( hice(i), himax ), himin )
          snowd(i) = min( snowd(i), hsmax )

          if (snowd(i) > (2.0*hice(i))) then
            print *, 'warning: too much snow :',snowd(i)
            snowd(i) = 2.0 * hice(i)
            print *,'fix: decrease snow depth to:',snowd(i)
          endif
        endif
      enddo

      call ice3lay
!  ---  inputs:                                                         !
!    &     ( im, kmi, fice, flag, hfi, hfd, sneti, focn, delt,          !
!  ---  outputs:                                                        !
!    &       snowd, hice, stsice, tice, snof, snowmt, gflux )           !

      do i = 1, im
        if (flag(i)) then
          if (tice(i) < timin) then
            print *,'warning: snow/ice temperature is too low:',tice(i)
            tice(i) = timin
            print *,'fix snow/ice temperature: reset it to:',tice(i)
          endif

          if (stsice(i,1) < timin) then
            print *,'warning: layer 1 ice temp is too low:',stsice(i,1)
            stsice(i,1) = timin
            print *,'fix layer 1 ice temp: reset it to:',stsice(i,1)
          endif

          if (stsice(i,2) < timin) then
            print *,'warning: layer 2 ice temp is too low:',stsice(i,2)
            stsice(i,2) = timin
            print *,'fix layer 2 ice temp: reset it to:',stsice(i,2)
          endif
        endif
      enddo

      do k = 1, kmi
        do i = 1, im
          if (flag(i)) then
            stc(i,k) = min(stsice(i,k), t0c)
          else
            stc(i,k) = t0c
          endif
        enddo
      enddo

!  --- ...  calculate sensible heat flux (& evap over sea ice)

      do i = 1, im
        if (flag(i)) then
          hflxi(i) = rch(i) * (tice(i) - theta1(i))
          hflxw(i) = rch(i) * (tgice - theta1(i))
          hflx(i) = fice(i)*hflxi(i) + ffw(i)*hflxw(i)
          evap(i) = fice(i)*evapi(i) + ffw(i)*evapw(i)
        endif
      enddo

!  --- ...  the rest of the output

      do i = 1, im
        if (flag(i)) then
          qsurf(i) = q1(i) + evap(i) / (elocp*rch(i))

!  --- ...  convert snow depth back to mm of water equivalent

          sheleg(i) = snowd(i) * 1000.0
          snwdph(i) = sheleg(i) * dsi             ! snow depth in mm
        endif
      enddo

      do i = 1, im
        if (flag(i)) then
          tem     = 1.0 / rho(i)
          hflx(i) = hflx(i) * tem * cpinv
          evap(i) = evap(i) * tem * hvapi
        endif
      enddo

      do i = 1, im
        if (.not. flag(i)) then
          hice(i) = 0.0
          fice(i) = 0.0
        endif
      enddo

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
!     flag     - logical, ice mask flag                              1    !
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

!  ---  inputs:
!     integer, intent(in) :: im, kmi

!     real , dimension(im), intent(in) :: fice, hfi,    
!    &       hfd, sneti, focn

!     real , intent(in) :: delt

!     logical, dimension(im), intent(in) :: flag

!  ---  input/outputs:
!     real , dimension(im), intent(inout) :: snowd,     
!    &       hice, tice, snof

!     real , dimension(im,kmi), intent(inout) :: stsice

!  ---  outputs:
!     real , dimension(im), intent(out) :: snowmt,      
!    &       gflux

!  ---  locals:
      real , dimension(im) :: hdi, hsni, a, b, ip,      
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

      do i = 1, im
        if (flag(i)) then
          tmpt=tice(i)
          snowd(i) = snowd(i) * dw / ds
          hdi(i) = (ds*snowd(i) + di*hice(i)) / dw

          if (hice(i) < hdi(i)) then
            snowd(i) = snowd(i) + hice(i) - hdi(i)
            hsni (i) = (hdi(i) - hice(i)) * ds / di
            hice (i) = hice(i) + hsni(i)
          endif

          snof(i) = snof(i) * dw / ds
          tice(i) = tice(i) - t0c
          stsice(i,1) = min(stsice(i,1)-t0c, tfi0)     ! degc
          stsice(i,2) = min(stsice(i,2)-t0c, tfi0)     ! degc

          ip(i) = i0 * sneti(i)         ! ip +v (in winton ip=-i0*sneti as sol -v)
          if (snowd(i) > 0.0) then
            tsf(i) = 0.0
            ip (i) = 0.0
          else
            tsf(i) = tfi
            ip (i) = i0 * sneti(i)      ! ip +v here (in winton ip=-i0*sneti)
          endif
          tice(i) = min(tice(i), tsf(i))

!  --- ...  compute ice temperature

          b(i) = hfd(i)
          a(i) = hfi(i) - sneti(i) + ip(i) - tice(i)*b(i)  ! +v sol input here
          k12(i) = ki4*ks / (ks*hice(i) + ki4*snowd(i))
          k32(i) = 2.0 * ki / hice(i)

          a10(i) = dici*hice(i)/dt2 + k32(i)*(dt4*k32(i) + dici*hice(i))
     &           / (dt6*k32(i) + dici*hice(i))
          b10(i) = -di*hice(i) * (ci*stsice(i,1) + li*tfi/stsice(i,1))  
     &           / dt2 - ip(i)                                          
     &           - k32(i)*(dt4*k32(i)*tfw + dici*hice(i)*stsice(i,2))   
     &           / (dt6*k32(i) + dici*hice(i))

          a1(i) = a10(i) + k12(i)*b(i) / (k12(i) + b(i))
          b1(i) = b10(i) + a(i)*k12(i) / (k12(i) + b(i))
          c1(i) = dili*tfi / dt2*hice(i)

          stsice(i,1) = -(sqrt(b1(i)*b1(i) - 4.0*a1(i)*c1(i))           
     &                + b1(i))/(2.0*a1(i))
          tice(i) = (k12(i)*stsice(i,1) - a(i)) / (k12(i) + b(i))

          if (tice(i) > tsf(i)) then
            a1(i) = a10(i) + k12(i)
            b1(i) = b10(i) - k12(i)*tsf(i)
            stsice(i,1) = -(sqrt(b1(i)*b1(i) - 4.0*a1(i)*c1(i))         
     &                  + b1(i)) / (2.0*a1(i))
            tice(i) = tsf(i)
            tmelt(i) = (k12(i)*(stsice(i,1) - tsf(i))                   
     &               - (a(i) + b(i)*tsf(i))) * delt
          else
            tmelt(i) = 0.0
            snowd(i) = snowd(i) + snof(i)*delt
          endif

          stsice(i,2) = (dt2*k32(i)*(stsice(i,1) + 2.0*tfw)             
     &                + dici*hice(i)*stsice(i,2))                       
     &                / (dt6*k32(i) + dici*hice(i))

          bmelt(i) = (focn(i) + ki4*(stsice(i,2) - tfw)/hice(i)) * delt

!  --- ...  resize the ice ...

          h1(i) = 0.5 * hice(i)
          h2(i) = 0.5 * hice(i)

!  --- ...  top ...

          if (tmelt(i) <= snowd(i)*dsli) then
            snowmt(i) = tmelt(i) / dsli
            snowd (i) = snowd(i) - tmelt(i)/dsli
          else
            snowmt(i) = snowd(i)
            h1(i) = h1(i) - (tmelt(i) - snowd(i)*dsli)                  
     &            / (di * (ci - li/stsice(i,1)) * (tfi - stsice(i,1)))
            snowd(i) = 0.0
          endif

!  --- ...  and bottom

          if (bmelt(i) < 0.0) then
            dh(i) = -bmelt(i) / (dili + dici*(tfi - tfw))
            stsice(i,2) = (h2(i)*stsice(i,2) + dh(i)*tfw)               
     &                  / (h2(i) + dh(i))
            h2(i) = h2(i) + dh(i)
          else
            h2(i) = h2(i) - bmelt(i) / (dili + dici*(tfi - stsice(i,2)))
          endif

!  --- ...  if ice remains, even up 2 layers, else, pass negative energy back in snow

          hice(i) = h1(i) + h2(i)

          if (hice(i) > 0.0) then
            if (h1(i) > 0.5*hice(i)) then
              f1(i) = 1.0 - 2.0*h2(i) / hice(i)
              stsice(i,2) = f1(i)                                       
     &                    * (stsice(i,1) + li*tfi/(ci*stsice(i,1)))     
     &                    + (1.0 - f1(i))*stsice(i,2)

              if (stsice(i,2) > tfi) then
                hice(i) = hice(i) - h2(i)*ci*(stsice(i,2) - tfi)        
     &                  / (li*delt)
                stsice(i,2) = tfi
              endif
            else
              f1(i) = 2.0 * h1(i) / hice(i)
              stsice(i,1) = f1(i)                                       
     &                    * (stsice(i,1) + li*tfi/(ci*stsice(i,1)))     
     &                    + (1.0 - f1(i))*stsice(i,2)
              stsice(i,1) = (stsice(i,1) - sqrt(stsice(i,1)*stsice(i,1) 
     &                    - 4.0*tfi*li/ci)) / 2.0
            endif

            k12(i) = ki4*ks / (ks*hice(i) + ki4*snowd(i))
            gflux(i) = k12(i) * (stsice(i,1) - tice(i))
          else
            snowd(i) = snowd(i) + (h1(i)*(ci*(stsice(i,1) - tfi)        
     &               - li*(1.0 - tfi/stsice(i,1)))                      
     &               + h2(i)*(ci*(stsice(i,2) - tfi) - li)) / li

            hice(i) = max(0.0, snowd(i)*ds/di)
            snowd(i) = 0.0
            stsice(i,1) = tfw
            stsice(i,2) = tfw
            gflux(i) = 0.0
          endif   ! end if_hice_block

          gflux(i) = fice(i) * gflux(i)
          snowmt(i) = snowmt(i) * ds / dw
          snowd(i) = snowd(i) * ds / dw
          tice(i) = tice(i) + t0c
          stsice(i,1) = stsice(i,1) + t0c
          stsice(i,2) = stsice(i,2) + t0c
              
          
        endif   ! end if_flag_block
      enddo   ! end do_i_loop

      return
      end subroutine ice3lay

      end subroutine sfc_sice

      end module seaice_mod
