
      module sfc_ocean_mod
        use sat_vapor_pres_mod, only : compute_qs
        implicit none
        private

        public :: sfc_ocean

        contains

!-----------------------------------
      subroutine sfc_ocean                                              
!...................................
!  ---  inputs:
     &     ( im, ps, wind, t1, q1, tskin, cm, ch,                 
     &       prsl1, prslki, slimsk,                     
!  ---  outputs:
     &       qsurf, cmm, chh, evap, hflx                                
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_ocean                                                     !
!       inputs:                                                         !
!          ( im, ps, wind, t1, q1, tskin, cm, ch,                 !
!            prsl1, prslki, slimsk,                    !
!       outputs:                                                        !
!            qsurf, cmm, chh, evap, hflx )                              !
!                                                                       !
!                                                                       !
!  subprograms/functions called: fpvs                                   !
!                                                                       !
!                                                                       !
!  program history log:                                                 !
!         xxxx  --             created                                  !
!    oct  2006  -- h. wei      modified (need description)              !
!    apr  2009  -- y.-t. hou   modified to match the modified gbphys_v.f!
!                     rmoved unused variable from argument list.        !
!                     reformatted code and added program documentation. !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horizontal dimension                     1    !
!     ps       - real, surface pressure                            im   !
!     wind   - real, u/v component of surface layer wind         im   !
!     t1       - real, surface layer mean temperature ( k )        im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, surface layer mean pressure                 im   !
!     prslki   - real,                                             im   !
!     slimsk   - real, sea/land/ice mask (=0/1/2)                  im   !
!                                                                       !
!  outputs:                                                             !
!     qsurf    - real, specific humidity at sfc                    im   !
!     cmm      - real,                                             im   !
!     chh      - real,                                             im   !
!     evap     - real, evaperation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!                                                                       !
! ===================================================================== !

      use constants_mod, only : cp => CP_AIR, RDGAS, RVGAS, hvap => HLV

      implicit none
!
!  ---  constant parameters:

      real, parameter :: eps = RDGAS/RVGAS
      real, parameter :: epsm1 = RDGAS/RVGAS-1.
      real, parameter :: rvrdm1 = RVGAS/RDGAS-1.
      real, parameter :: rd = RDGAS
      
      real , parameter :: cpinv  = 1.0/cp
      real , parameter :: hvapi  = 1.0/hvap
      real , parameter :: elocp  = hvap/cp

!  ---  inputs:
      integer, intent(in) :: im

      real , dimension(im), intent(in) :: ps, wind,   
     &      t1, q1, tskin, cm, ch, prsl1, prslki

      logical, dimension(im), intent(in) :: slimsk

!  ---  outputs:
      real , dimension(im), intent(out) :: qsurf,       
     &       cmm, chh, evap, hflx

!  ---  locals:
      real , dimension(im) :: psurf, ps1, q0, qss,      
     &       rch, rho, theta1, tv1

      real  :: tem

      integer :: i

      logical :: flag(im)
!
!===> ...  begin here
!
!  --- ...  flag for open water
      do i = 1, im
         flag(i) = slimsk(i)
      enddo

!  --- ...  initialize variables. all units are supposedly m.k.s. unless specified
!           psurf is in pascals, wind is wind speed, theta1 is adiabatic surface
!           temp from level 1, rho is density, qss is sat. hum. at surface

      do i = 1, im
        if ( flag(i) ) then
          psurf(i) = 1000.0 * ps(i)
          ps1(i)   = 1000.0 * prsl1(i)
          theta1(i) = t1(i) * prslki(i)
          q0(i) = max( q1(i), 1.0e-8 )
          tv1(i) = t1(i) * (1.0 + rvrdm1*q0(i))
          rho(i) = ps1(i) / (rd*tv1(i))

          call compute_qs(tskin(i),psurf(i),qss(i))
        endif
      enddo

      do i = 1, im
        if ( flag(i) ) then
          evap(i) = 0.0
          hflx(i) = 0.0
        endif
      enddo

!  --- ...  rcp = rho cp ch v

      do i = 1, im
        if ( flag(i) ) then
          rch(i) = rho(i) * cp * ch(i) * wind(i)
          cmm(i) = cm(i) * wind(i)
          chh(i) = rho(i) * ch(i) * wind(i)
        endif
      enddo

!  --- ...  sensible and latent heat flux over open water

      do i = 1, im
        if ( flag(i) ) then
          hflx(i) = rch(i) * (tskin(i) - theta1(i))
          evap(i) = elocp*rch(i) * (qss(i) - q0(i))
          qsurf(i) = qss(i)
        endif
      enddo

      do i = 1, im
        if ( flag(i) ) then
          tem     = 1.0 / rho(i)
          hflx(i) = hflx(i) * tem * cpinv
          evap(i) = evap(i) * tem * hvapi
        endif
      enddo
!
      return
!...................................
      end subroutine sfc_ocean
!-----------------------------------

        end module sfc_ocean_mod
