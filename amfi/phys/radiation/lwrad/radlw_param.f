!!!!!  ==============================================================  !!!!!
!!!!!             lw-rrtm3 radiation package description               !!!!!
!!!!!  ==============================================================  !!!!!
!                                                                          !
!   this package includes ncep's modifications of the rrtm-lw radiation    !
!   code from aer inc.                                                     !
!                                                                          !
!    the rrtm3 package includes these parts:                               !
!                                                                          !
!       'radlw_rrtm3_param.f'                                              !
!       'radlw_rrtm3_datatb.f'                                             !
!       'radlw_rrtm3_main.f'                                               !
!                                                                          !
!    the 'radlw_rrtm3_param.f' contains:                                   !
!                                                                          !
!       'module_radlw_cntr_para'   -- control parameters set up            !
!       'module_radlw_parameters'  -- band parameters set up               !
!                                                                          !
!    the 'radlw_rrtm3_datatb.f' contains:                                  !
!                                                                          !
!       'module_radlw_avplank'     -- plank flux data                      !
!       'module_radlw_ref'         -- reference temperature and pressure   !
!       'module_radlw_cldprlw'     -- cloud property coefficients          !
!       'module_radlw_kgbnn'       -- absorption coeffients for 16         !
!                                     bands, where nn = 01-16              !
!                                                                          !
!    the 'radlw_rrtm3_main.f' contains:                                    !
!                                                                          !
!       'module_radlw_main'        -- main lw radiation transfer           !
!                                                                          !
!    in the main module 'module_radlw_main' there are only two             !
!    externally callable subroutines:                                      !
!                                                                          !
!       'lwrad'     -- main rrtm3 lw radiation routine                     !
!       'rlwinit'   -- to initialize rrtm3 lw radiation                    !
!                                                                          !
!    all the lw radiation subprograms become contained subprograms         !
!    in module 'module_radlw_rrtm' and many of them are not directly       !
!    accessable from places outside the module.                            !
!                                                                          !
!    compilation sequence is:                                              !
!                                                                          !
!       'radlw_rrtm3_param.f'                                              !
!       'radlw_rrtm3_datatb.f'                                             !
!       'radlw_rrtm3_main.f'                                               !
!                                                                          !
!    and all should be put in front of routines that use lw modules        !
!                                                                          !
!    ncep modifications history log:                                       !
!                                                                          !
!       see list in program "radlw_rrtm3_main.f"                           !
!                                                                          !
!!!!!  ==============================================================  !!!!!
!!!!!                         end descriptions                         !!!!!
!!!!!  ==============================================================  !!!!!


!========================================!
      module module_radlw_cntr_para      !
!........................................!
!
        implicit   none
!
        integer :: ilwrate, iaerlw, irgaslw, icfclw, iflagliq, iflagice

!
!  ---  set up control parameters for lw radiation
!
        parameter ( ilwrate=2 )     !===> ... lw heating rate unit selection
                                    ! =1: output in k/day
                        !(default)  ! =2: output in k/second

        parameter ( iaerlw=1 )      !===> ... control flag for aerosols
                                    ! =0: do not include aerosol effect
                        !(default)  ! =1: aeros opt prop are calc for each spectral band
                                    ! =2: broad band aeros opt prop are used for all bands

        parameter ( irgaslw=1 )     !===> ... control flag for rare gases (ch4,n2o,o2, etc.)
                                    ! =0: do not include rare gases
                        !(default)  ! =1: include all rare gases

        parameter ( icfclw=1  )     !===> ... control flag for halocarbon (cfc) gases
                                    ! =0: do not include cfc gases
                        !(default)  ! =1: include all cfc gases

        parameter ( iflagliq=1 )    !===> ... liq-cloud optical properties contrl flag
                                    ! =0: input cloud opt depth, ignor iflagice setting
                        !(default)  ! =1: input cwp rew, hu and stamnes(1993) method for liq cld

        parameter ( iflagice=3 )    !===> ... ice-cloud optical properties contrl flag
                                    !         used only when iflagliq>0, otherwise ignored
                                    ! =1: input cip rei, ebert and curry(1997) for ice clouds
                                    ! =2: input cip rei, streamer (1996) for ice clouds
                        !(default)  ! =3: input cip rei, fu (1998) for ice clouds


!
!........................................!
      end module module_radlw_cntr_para  !
!========================================!



!========================================!
      module module_radlw_parameters     !
!........................................!

      implicit none
!
      public
!
!  ---  define type construct for radiation fluxes at toa
!
      type :: topflw_type
        real :: upfxc         ! total sky upward flux at toa
        real :: upfx0         ! clear sky upward flux at toa
      end type
!
!  ---  define type construct for radiation fluxes at surface
!
      type :: sfcflw_type
        real :: upfxc         ! total sky upward flux at sfc
        real :: upfx0         ! clear sky upward flux at sfc
        real :: dnfxc         ! total sky downward flux at sfc
        real :: dnfx0         ! clear sky downward flux at sfc
      end type
!
!  ---  define type construct for optional radiation flux profiles
!
      type :: proflw_type
        real :: upfxc         ! level up flux for total sky
        real :: dnfxc         ! level dn flux for total sky
        real :: upfx0         ! level up flux for clear sky
        real :: dnfx0         ! level dn flux for clear sky
      end type
!
!  ---  parameter constants for lw band structures
!
      integer, parameter :: NBANDS = 16         ! num of total spectral bands
      integer, parameter :: NGPTLW = 140        ! num of total g-points
      integer, parameter :: NTBL   = 10000      ! lookup table dimension
      integer, parameter :: MAXGAS = 7          ! max num of absorbing gases
      integer, parameter :: MAXXSEC= 4          ! num of halocarbon gases
      integer, parameter :: NRATES = 6          ! num of ref rates of binary species
      integer, parameter :: NPLNK  = 181        ! dim for plank function table

      integer, parameter :: NBDLW  = NBANDS

!  ---  number of g-point in each band
      integer  :: NG01, NG02, NG03, NG04, NG05, NG06, NG07, NG08,       &
     &            NG09, NG10, NG11, NG12, NG13, NG14, NG15, NG16
      parameter (NG01=10, NG02=12, NG03=16, NG04=14, NG05=16, NG06=08,  &
     &           NG07=12, NG08=08, NG09=12, NG10=06, NG11=08, NG12=08,  &
     &           NG13=04, NG14=02, NG15=02, NG16=02)

!  ---  begining index of each band
      integer  :: NS01, NS02, NS03, NS04, NS05, NS06, NS07, NS08,       &
     &            NS09, NS10, NS11, NS12, NS13, NS14, NS15, NS16
      parameter (NS01=00, NS02=10, NS03=22, NS04=38, NS05=52, NS06=68,  &
     &           NS07=76, NS08=88, NS09=96, NS10=108, NS11=114,         &
     &           NS12=122, NS13=130, NS14=134, NS15=136, NS16=138)

!  ---  band indices for each g-point
      integer, dimension(NGPTLW) :: NGB
      data NGB(:) / 10*1, 12*2, 16*3, 14*4, 16*5,  8*6, 12*7,  8*8,     & ! band  1- 8
     &              12*9, 6*10, 8*11, 8*12, 4*13, 2*14, 2*15, 2*16 /      ! band  9-16

!  ---  band spectrum structures (wavenumber in cm**-1)
      real :: wvnlw1(NBANDS), wvnlw2(NBANDS)
      data wvnlw1  /                                                    &
     &         10.,  351.,  501.,  631.,  701.,  821.,  981., 1081.,    &
     &       1181., 1391., 1481., 1801., 2081., 2251., 2381., 2601. /
      data wvnlw2  /                                                    &
     &        350.,  500.,  630.,  700.,  820.,  980., 1080., 1180.,    &
     &       1390., 1480., 1800., 2080., 2250., 2380., 2600., 3250. /


!........................................!
      end module module_radlw_parameters !
!========================================!
