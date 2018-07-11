!!!!!  ==============================================================  !!!!!
!!!!!              sw-rrtm3 radiation package description              !!!!!
!!!!!  ==============================================================  !!!!!
!                                                                          !
!   this package includes ncep's modifications of the rrtm-sw radiation    !
!   code from aer inc.                                                     !
!                                                                          !
!   the sw-rrtm3 package includes these parts:                             !
!                                                                          !
!      'radsw_rrtm3_param.f'                                               !
!      'radsw_rrtm3_datatb.f'                                              !
!      'radsw_rrtm3_main.f'                                                !
!                                                                          !
!   the 'radsw_rrtm3_param.f' contains:                                    !
!                                                                          !
!      'module_radsw_cntr_para'   -- control parameters set up             !
!      'module_radsw_parameters'  -- band parameters set up                !
!                                                                          !
!   the 'radsw_rrtm3_datatb.f' contains:                                   !
!                                                                          !
!      'module_radsw_ref'         -- reference temperature and pressure    !
!      'module_radsw_cldprtb'     -- cloud property coefficients table     !
!      'module_radsw_sflux'       -- spectral distribution of solar flux   !
!      'module_radsw_kgbnn'       -- absorption coeffients for 14          !
!                                    bands, where nn = 16-29               !
!                                                                          !
!   the 'radsw_rrtm3_main.f' contains:                                     !
!                                                                          !
!      'module_radsw_main'        -- main sw radiation transfer            !
!                                                                          !
!   in the main module 'module_radsw_main' there are only two              !
!   externally callable subroutines:                                       !
!                                                                          !
!      'swrad'      -- main rrtm3 sw radiation routine                     !
!      'rswinit'    -- initialization routine                              !
!                                                                          !
!   all the sw radiation subprograms become contained subprograms          !
!   in module 'module_radsw_main' and many of them are not directly        !
!   accessable from places outside the module.                             !
!                                                                          !
!   compilation sequence is:                                               !
!                                                                          !
!      'radsw_rrtm3_param.f'                                               !
!      'radsw_rrtm3_datatb.f'                                              !
!      'radsw_rrtm3_main.f'                                                !
!                                                                          !
!   and all should be put in front of routines that use sw modules         !
!                                                                          !
!   ncep modifications history log:                                        !
!                                                                          !
!       see list in program "radsw_rrtm3_main.f"                           !
!                                                                          !
!!!!!  ==============================================================  !!!!!
!!!!!                         end descriptions                         !!!!!
!!!!!  ==============================================================  !!!!!



!========================================!
      module module_radsw_cntr_para      !
!........................................!
!
        implicit   none
!
        integer :: iswrate, iaersw, imodsw, irgassw
        integer :: iflagliq, iflagice

!
!  ---  set up control parameters for sw radiation
!
        parameter ( iswrate=2 )     !===> ... flag for heating rate unit
                                    ! =1: output in k/day
                        ! (default) ! =2: output in k/second
        parameter ( iaersw=1 )      !===> ... flag for aerosols
                        ! (default) ! =0: without aerosol effect
                                    ! =1: include aerosol effect

        parameter ( imodsw=2 )      !===> ... flag for 2-stream transfer scheme
                                    ! =1: delta-eddington    (joseph et al., 1976)
                        ! (default) ! =2: pifm               (zdunkowski et al., 1980)
                                    ! =3: discrete ordinates (liou, 1973)

        parameter ( irgassw=1 )     !===> ... control flag for rare gases (ch4,n2o,o2, etc.)
                                    ! =0: do not include rare gases
                        ! (default) ! =1: include all rare gases

        parameter ( iflagliq=1 )    !===> ... liq-cloud optical properties control flag
                                    ! =0: input cloud opt depth, ignor iflagice setting
                        !(default)  ! =1: input cwp,rew, hu and stamnes(1993) method for liq cld

        parameter ( iflagice=3 )    !===> ... ice-cloud optical properties control flag
                                    !         used only when iflagliq>0, otherwise ignored
                                    ! =0: not set up yet
                                    ! =1: input cip rei, ebert & curry (1992) for ice clouds
                                    ! =2: input cip rei, streamer v3 (2001) for ice clouds
                        !(default)  ! =3: input cip rei, fu (1996) for ice clouds


!
!........................................!
      end module module_radsw_cntr_para  !
!========================================!



!========================================!
      module module_radsw_parameters     !
!........................................!

      implicit   none
!
      public
!
!  ---  define type construct for radiation fluxes at toa
!
      type :: topfsw_type
        real :: upfxc         ! total sky upward flux at toa
        real :: dnfxc         ! total sky downward flux at toa
        real :: upfx0         ! clear sky upward flux at toa
      end type
!
!  ---  define type construct for radiation fluxes at surface
!
      type :: sfcfsw_type
        real :: upfxc         ! total sky upward flux at sfc
        real :: dnfxc         ! total sky downward flux at sfc
        real :: upfx0         ! clear sky upward flux at sfc
        real :: dnfx0         ! clear sky downward flux at sfc
      end type
!
!  ---  define type construct for optional radiation flux profiles
!
      type :: profsw_type
        real :: upfxc         ! total sky level upward flux
        real :: dnfxc         ! total sky level downward flux
        real :: upfx0         ! clear sky level upward flux
        real :: dnfx0         ! clear sky level downward flux
      end type
!
!  ---  define type construct for optional component downward fluxes at surface
!
      type :: cmpfsw_type
        real :: uvbfc         ! total sky downward uv-b flux at sfc
        real :: uvbf0         ! clear sky downward uv-b flux at sfc

        real :: nirbm         ! sfc downward nir direct beam flux
        real :: nirdf         ! sfc downward nir diffused flux
        real :: visbm         ! sfc downward uv+vis direct beam flx
        real :: visdf         ! sfc downward uv+vis diffused flux
      end type
!
!  ---  parameter constants for sw band structures
!
      integer, parameter :: NBLOW  = 16            ! band range lower limit
      integer, parameter :: NBHGH  = 29            ! band range upper limit
      integer, parameter :: NBANDS = NBHGH-NBLOW+1 ! num of spectral bands
      integer, parameter :: NGPTSW = 112           ! total num of g-point in all bands
      integer, parameter :: NGMAX  = 16            ! max num of g-point in one band
      integer, parameter :: MAXGAS = 7             ! max num of absorbing gases
      integer, parameter :: NTBMX  = 10000         ! indx upper lim of trans table

      integer, parameter :: NSWSTR = 1
      integer, parameter :: NSWEND = NBANDS
      integer, parameter :: NBDSW  = NBANDS

!  ---  number of g-point in each band
      integer  :: NG16, NG17, NG18, NG19, NG20, NG21, NG22,             &
     &            NG23, NG24, NG25, NG26, NG27, NG28, NG29
      parameter ( NG16=06, NG17=12, NG18=08, NG19=08, NG20=10,          &
     &            NG21=10, NG22=02, NG23=10, NG24=08, NG25=06,          &
     &            NG26=06, NG27=08, NG28=06, NG29=12)

      integer, dimension(NBLOW:NBHGH) :: NG
      data  NG / NG16, NG17, NG18, NG19, NG20, NG21, NG22,              &
     &           NG23, NG24, NG25, NG26, NG27, NG28, NG29  /

!  ---  starting index of each band
      integer  :: NS16, NS17, NS18, NS19, NS20, NS21, NS22,             &
     &            NS23, NS24, NS25, NS26, NS27, NS28, NS29
      parameter ( NS16=00,         NS17=NS16+NG16,  NS18=NS17+NG17,     &
     &            NS19=NS18+NG18,  NS20=NS19+NG19,  NS21=NS20+NG20,     &
     &            NS22=NS21+NG21,  NS23=NS22+NG22,  NS24=NS23+NG23,     &
     &            NS25=NS24+NG24,  NS26=NS25+NG25,  NS27=NS26+NG26,     &
     &            NS28=NS27+NG27,  NS29=NS28+NG28  )

      integer, dimension(NBLOW:NBHGH) :: NGS
      data  NGS / NS16, NS17, NS18, NS19, NS20, NS21, NS22,             &
     &            NS23, NS24, NS25, NS26, NS27, NS28, NS29  /

!  ---  band index for each g-point
      integer, dimension(NGPTSW)      :: NGB
      data NGB(:) / 16,16,16,16,16,16,                                  & ! band 16
     &              17,17,17,17,17,17,17,17,17,17,17,17,                & ! band 17
     &              18,18,18,18,18,18,18,18,                            & ! band 18
     &              19,19,19,19,19,19,19,19,                            & ! band 19
     &              20,20,20,20,20,20,20,20,20,20,                      & ! band 20
     &              21,21,21,21,21,21,21,21,21,21,                      & ! band 21
     &              22,22,                                              & ! band 22
     &              23,23,23,23,23,23,23,23,23,23,                      & ! band 23
     &              24,24,24,24,24,24,24,24,                            & ! band 24
     &              25,25,25,25,25,25,                                  & ! band 25
     &              26,26,26,26,26,26,                                  & ! band 26
     &              27,27,27,27,27,27,27,27,                            & ! band 27
     &              28,28,28,28,28,28,                                  & ! band 28
     &              29,29,29,29,29,29,29,29,29,29,29,29 /                 ! band 29

!  ---  band wavenumber intervals
      real , dimension(NBANDS):: wvnum1, wvnum2
      data wvnum1(:)    /                                               &
     &         2600.0, 3251.0, 4001.0, 4651.0, 5151.0, 6151.0, 7701.0,  &
     &         8051.0,12851.0,16001.0,22651.0,29001.0,38001.0,  820.0 /
      data wvnum2(:)    /                                               &
     &         3250.0, 4000.0, 4650.0, 5150.0, 6150.0, 7700.0, 8050.0,  &
     &        12850.0,16000.0,22650.0,29000.0,38000.0,50000.0, 2600.0 /

!
!........................................!
      end module module_radsw_parameters !
!========================================!
