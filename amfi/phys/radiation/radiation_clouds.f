!!!!!              module_radiation_clouds description             !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!    the 'radiation_clouds.f' contains:                                !
!                                                                      !
!       'module_radiation_clouds' ---  compute cloud related quantities!
!                for radiation computations                            !
!                                                                      !
!    the following are the externally accessable subroutines:          !
!                                                                      !
!       'cldinit'            --- initialization routine                !
!          inputs:                                                     !
!           (si, NLAY, np3d, icld, me)                          !
!          outputs:                                                    !
!           ( none )                                                   !
!                                                                      !
!       'progcld1'           --- zhao/moorthi prognostic cloud scheme  !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,qlyr,qstl,rhly,clw,                        !
!            xlat,slmsk,                                          !
!            IX, NLAY, NLP1, iovr,                              !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)                                    !
!                                                                      !
!    internal accessable only subroutines:                             !
!       'gethml'             --- get diagnostic hi, mid, low clouds    !
!                                                                      !
!                                                                      !
!    cloud array description:                                          !
!                ---  for prognostic cloud: icld=1  ---                !
!          clouds(:,:,1)  -  layer total cloud fraction                !
!          clouds(:,:,2)  -  layer cloud liq water path                !
!          clouds(:,:,3)  -  mean effective radius for liquid cloud    !
!          clouds(:,:,4)  -  layer cloud ice water path                !
!          clouds(:,:,5)  -  mean effective radius for ice cloud       !
!          clouds(:,:,6)  -  layer rain drop water path                !
!          clouds(:,:,7)  -  mean effective radius for rain drop       !
!   **     clouds(:,:,8)  -  layer snow flake water path               !
!          clouds(:,:,9)  -  mean effective radius for snow flake      !
!   ** fu's scheme need to be normalized by snow density (g/m**3/1.0e6)!
!                ---  for diagnostic cloud: icld=0  ---                !
!          clouds(:,:,1)  -  layer total cloud fraction                !
!          clouds(:,:,2)  -  layer cloud optical depth                 !
!          clouds(:,:,3)  -  layer cloud single scattering albedo      !
!          clouds(:,:,4)  -  layer cloud asymmetry factor              !
!                                                                      !
!    external modules referenced:                                      !
!                                                                      !
!       'module machine'             in 'machine.f'                    !
!       'module physcons'            in 'physcons.f'                   !
!       'module module_microphysics' in 'module_bfmicrophysics.f'      !
!                                                                      !
!                                                                      !
!    modification history log:                                         !
!                                                                      !
!       apr 2003,  yu-tai hou                                          !
!                  created 'module_rad_clouds' from combining the      !
!                  original subroutine 'cldjms', 'cldprp', and 'gcljms'!
!                                                                      !
!       may 2004,  yu-tai hou                                          !
!                  incorporate ferrier's cloud microphysics scheme.    !
!                                                                      !
!       apr 2005,  yu-tai hou                                          !
!                  modified cloud array and module structures.         !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!
      module module_radiation_clouds

      use constants_mod, only : con_pi => PI, con_g => GRAV
      use constants_mod, only : con_rd => RDGAS, RVGAS

      implicit none
      private

      integer, parameter, public :: NF_CLDS = 9   ! number of fields in cloud array

!  ---  pressure limits of cloud domain interfaces (low,mid,high) in mb (0.1kPa)
      real :: ptopc(4,2)

      data ptopc / 1050., 642., 350., 0.0,  1050., 750., 500., 0.0 /

      real, parameter :: climit = 0.001, climit2=0.05
      real, parameter :: ovcst  = 1.0 - 1.0e-8

!  ---  set default quantities as parameters (for prognostic cloud)

      real, parameter :: reliq_def = 10.0    ! default liq radius to 10 micron
      real, parameter :: reice_def = 50.0    ! default ice radius to 50 micron
      real, parameter :: rrain_def = 1000.0  ! default rain radius to 1000 micron
      real, parameter :: rsnow_def = 250.0   ! default snow radius to 250 micron

!  ---  set look-up table dimensions and other parameters (for diagnostic cloud)

      real, parameter :: cldssa_def = 0.99   ! default cld single scat albedo
      real, parameter :: cldasy_def = 0.84   ! default cld asymmetry factor

      real, parameter :: xlim=5.0

      real, parameter :: vvcld1= 0.0003e0
      real, parameter :: vvcld2=-0.0005e0

      integer         :: llyr

      real, parameter :: con_ttp = 2.7316e+2
      real, parameter :: con_fvirt = RVGAS/con_rd-1.

      public :: progcld1, cldinit

      contains

!--------------------------------------------------------------------------------   
      subroutine cldinit(si)
!--------------------------------------------------------------------------------   
      implicit none
      real, intent(in) :: si(:)

      integer :: k, kl, NLAY

        !  ---  compute llyr - the top of bl cld and is topmost non cld(low) layer
        !       for stratiform (at or above lowest 0.1 of the atmosphere)

        NLAY = size(si,1)

        kl = 2
        lab_do_k1 : do k = 1, NLAY
          kl = k
          if (si(k) < 0.9e0) exit lab_do_k1
        enddo  lab_do_k1

        llyr = kl - 1
      return
      end subroutine cldinit

!--------------------------------------------------------------------------------   
      subroutine progcld1 ( plyr,plvl,tlyr,qlyr,qstl,rhly,clw,
     &   xlat,slmsk,IX, NLAY, NLP1, iovr, sashal, crick_proof, ccnorm,
     &   clouds,clds,mtop,mbot)
!--------------------------------------------------------------------------------   

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld1    computes cloud related quantities using    !
!   zhao/moorthi's prognostic cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cldinit".                                       !
!                                                                       !
! program history log:                                                  !
!      11-xx-1992   y.h., k.a.c, a.k. - cloud parameterization          !
!         'cldjms' patterned after slingo and slingo's work (jgr,       !
!         1992), stratiform clouds are allowed in any layer except      !
!         the surface and upper stratosphere. the relative humidity     !
!         criterion may cery in different model layers.                 !
!      10-25-1995   kenneth campana   - tuned cloud rh curves           !
!         rh-cld relation from tables created using mitchell-hahn       !
!         tuning technique on airforce rtneph observations.             !
!      11-02-1995   kenneth campana   - the bl relationships used       !
!         below llyr, except in marine stratus regions.                 !
!      04-11-1996   kenneth campana   - save bl cld amt in cld(,5)      !
!      12-29-1998   s. moorthi        - prognostic cloud method         !
!      04-15-2003   yu-tai hou        - rewritten in frotran 90         !
!         modulized form, seperate prognostic and diagnostic methods    !
!         into two packages.                                            !
!                                                                       !
! usage:         call progcld1                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY) : layer cloud condensate amount                     !
!   xlat  (IX)      : grid latitude in radians                          !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!   iovr            : control flag for cloud overlap                    !
!                     =0 random overlapping clouds                      !
!                     =1 max/ran overlapping clouds                     !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1, iovr

      real, dimension(ix,nlay), intent(in) :: plyr,  
     &       tlyr, qlyr, qstl, rhly, clw

      real, dimension(ix,nlp1), intent(in) :: plvl

      real, dimension(ix),   intent(in) :: xlat
      real, dimension(ix), intent(in) :: slmsk
      logical, intent(in) :: sashal, crick_proof, ccnorm

!  ---  outputs
      real, dimension(ix,nlay,NF_CLDS), intent(out) :: clouds

      real, dimension(ix,5),   intent(out) :: clds

      integer, dimension(ix,3),   intent(out) :: mtop,mbot

!  ---  local variables:
      real, dimension(IX,NLAY) :: cldtot, cldcnv, 
     &       cwp, cip, crp, csp, rew, rei, res, rer, delp, tem2d, clwf

      real :: ptop1(IX,4)

      real :: clwmin, clwm, clwt, onemrh, value, 
     &       tem1, tem2, tem3

      integer, dimension(IX) :: kinver

      integer :: i, k, id, id1

      logical :: inversn(IX)

      clouds(:,:,:) = 0.0

      do k = 1, NLAY
        do i = 1, IX
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          clwf(i,k)   = 0.0
        enddo
      enddo
!
      if (crick_proof) then
        do i = 1, IX
          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay)
        enddo
        do k = 2, NLAY-1
          do i = 1, IX
            clwf(i,K) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
          enddo
        enddo
      else
        do k = 1, NLAY
          do i = 1, IX
            clwf(i,k) = clw(i,k)
          enddo
        enddo
      endif

!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = max( 0.0, 4.0*abs(xlat(i))/con_pi-1.0 )
          ptop1(i,id) = ptopc(id,1) + tem1*tem2
        enddo
      enddo

!  ---  compute liquid/ice condensate path in g/m**2

      tem1 = 1.0e+5 / con_g

        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
            clwt     = max(0.0, clwf(i,k)) * tem1 * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo

!  ---  effective liquid cloud droplet radius over land

      do i = 1, IX
        if (nint(slmsk(i)) == 1) then
          do k = 1, NLAY
            rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
          enddo
        endif
      enddo

!  ---  layer cloud fraction

        do i = 1, IX
          inversn(i) = .false.
          kinver (i) = NLAY
        enddo

        do k = 2, NLAY
          do i = 1, IX
            if (plyr(i,k) > 600.0 .and. (.not.inversn(i))) then
              tem1 = tlyr(i,k+1) - tlyr(i,k)

              if (tem1 > 0.1 .and. tlyr(i,k) > 278.0) then
                inversn(i) = .true.
                kinver(i)  = k
              endif
            endif
          enddo
        enddo

        clwmin = 0.0
        if (.not. sashal) then
        do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt .or.                                    &
     &         (inversn(i) .and. k <= kinver(i)) ) then

              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

              tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
              tem1  = 2000.0 / tem1

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
        enddo
        else
                do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt .or.                                    &
     &         (inversn(i) .and. k <= kinver(i)) ) then

              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
              tem1  = 100.0 / tem1

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
        enddo
        endif

      where (cldtot < climit)
        cldtot = 0.0
        cwp    = 0.0
        cip    = 0.0
        crp    = 0.0
        csp    = 0.0
      endwhere
!
      if (ccnorm) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!  ---  effective ice cloud droplet radius

      tem1 = con_g / con_rd
      do k = 1, NLAY
        do i = 1, IX
          tem2 = tlyr(i,k) - con_ttp

          if (cip(i,k) > 0.0) then
            tem3 = tem1 * cip(i,k) * ( plyr(i,k) / delp(i,k) )          &
     &           / (tlyr(i,k) * (1.0 + con_fvirt * qlyr(i,k)))

            if (tem2 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem2 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem2 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif
          endif
        enddo
      enddo
!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
!         clouds(i,k,6) = 0.0
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = 0.0
          clouds(i,k,9) = rei(i,k)
        enddo
      enddo


!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which is
!  ---  also used by the lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX,NLAY, iovr,   
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )
      return

      end subroutine progcld1


!--------------------------------------------------------------------------------   
      subroutine gethml                                                 &
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX,NLAY, iovr, 
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )

!  ===================================================================  !
!                                                                       !
! abstract: compute high, mid, low, total, and boundary cloud fractions !
!   and cloud top/bottom layer indices for model diagnostic output.     !
!   the three cloud domain boundaries are defined by ptopc.  the cloud  !
!   overlapping method is defined by control flag 'iovr', which is also !
!   used by lw and sw radiation programs.                               !
!                                                                       !
! program history log:                                                  !
!      04-29-2004   yu-tai hou        - separated to become individule  !
!         subprogram to calculate averaged h,m,l,bl cloud amounts.      !
!                                                                       !
! usage:         call gethml                                            !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   ptop1 (IX,4)    : pressure limits of cloud domain interfaces        !
!                     (sfc,low,mid,high) in mb (100Pa)                  !
!   cldtot(IX,NLAY) : total or straiform cloud profile in fraction      !
!   cldcnv(IX,NLAY) : convective cloud (for diagnostic scheme only)     !
!   IX              : horizontal dimention                              !
!   NLAY            : vertical layer dimensions                         !
!   iovr            : control flag for cloud overlap                    !
!                     =0 random overlapping clouds                      !
!                     =1 max/ran overlapping clouds                     !
!                                                                       !
! output variables:                                                     !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none!

!  ---  inputs:
      integer, intent(in) :: IX, NLAY, iovr

      real, dimension(:,:), intent(in) :: plyr, ptop1,
     &       cldtot, cldcnv

!  ---  outputs
      real, dimension(:,:), intent(out) :: clds

      integer,               dimension(:,:), intent(out) :: mtop, mbot

!  ---  local variables:
      real :: cl1(IX), cl2(IX)

      real :: pcur, pnxt, ccur, cnxt

      integer, dimension(IX):: idom, kbt1, kth1, kbt2, kth2

      integer :: i, k, id, id1, kstr, kend, kinc

!
!===> ... begin here
!
      do i = 1, IX
        clds(i,1) = 0.0
        clds(i,2) = 0.0
        clds(i,3) = 0.0
        clds(i,4) = 0.0
        clds(i,5) = 0.0
        mtop(i,1) = 1
        mtop(i,2) = 1
        mtop(i,3) = 1
        mbot(i,1) = 1
        mbot(i,2) = 1
        mbot(i,3) = 1
        cl1 (i) = 1.0
        cl2 (i) = 1.0
      enddo

!  ---  total and bl clouds, where cl1, cl2 are fractions of clear-sky view
!       layer processed from surface and up

        kstr = 1
        kend = NLAY
        kinc = 1

      if (iovr == 0) then                       ! random overlap

        do k = kstr, kend, kinc
          do i = 1, IX
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
            if (ccur >= climit) cl1(i) = cl1(i) * (1.0 - ccur)
          enddo

          if (k == llyr) then
            do i = 1, IX
              clds(i,5) = 1.0 - cl1(i)          ! save bl cloud
            enddo
          endif
        enddo

        do i = 1, IX
          clds(i,4) = 1.0 - cl1(i)              ! save total cloud
        enddo

      else                                      ! max/ran overlap

        do k = kstr, kend, kinc
          do i = 1, IX
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
            if (ccur >= climit) then             ! cloudy layer
              cl2(i) = min( cl2(i), (1.0 - ccur) )
            else                                ! clear layer
              cl1(i) = cl1(i) * cl2(i)
              cl2(i) = 1.0
            endif
          enddo

          if (k == llyr) then
            do i = 1, IX
              clds(i,5) = 1.0 - cl1(i) * cl2(i) ! save bl cloud
            enddo
          endif
        enddo

        do i = 1, IX
          clds(i,4) = 1.0 - cl1(i) * cl2(i)     ! save total cloud
        enddo

      endif                                     ! end_if_iovr

!  ---  high, mid, low clouds, where cl1, cl2 are cloud fractions
!       layer processed from one layer below llyr and up

        do i = 1, IX
          cl1 (i) = 0.0
          cl2 (i) = 0.0
          kbt1(i) = 1
          kbt2(i) = 1
          kth1(i) = 0
          kth2(i) = 0
          idom(i) = 1
        enddo

        do k = llyr+1, NLAY
          do i = 1, IX
            id = idom(i)
            id1= id + 1

            pcur = plyr(i,k)
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))

            if (k < NLAY) then
              pnxt = plyr(i,k+1)
              cnxt = min( ovcst, max( cldtot(i,k+1), cldcnv(i,k+1) ))
            else
              pnxt = -1.0
              cnxt = 0.0
            endif

            if (pcur < ptop1(i,id1)) then
              id = id + 1
              id1= id1 + 1
              idom(i) = id
            endif

            if (ccur >= climit) then
              if (kth2(i) == 0) kbt2(i) = k
              kth2(i) = kth2(i) + 1

              if (iovr == 0) then
                cl2(i) = cl2(i) + ccur - cl2(i)*ccur
              else
                cl2(i) = max( cl2(i), ccur )
              endif

              if (cnxt < climit .or. pnxt < ptop1(i,id1)) then
                kbt1(i) = nint( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i))       &
     &                  / (cl1(i) + cl2(i)) )
                kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i))       &
     &                  / (cl1(i) + cl2(i)) )
                cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)

                kbt2(i) = 1
                kth2(i) = 0
                cl2 (i) = 0.0
              endif     ! end_if_cnxt_or_pnxt
            endif       ! end_if_ccur

            if (pnxt < ptop1(i,id1)) then
              clds(i,id) = cl1(i)
              mtop(i,id) = max( kbt1(i), kbt1(i)+kth1(i)-1 )
              mbot(i,id) = kbt1(i)

              cl1 (i) = 0.0
              kbt1(i) = 1
              kth1(i) = 0
            endif     ! end_if_pnxt

          enddo       ! end_do_i_loop
        enddo         ! end_do_k_loop

      return
!...................................
      end subroutine gethml
!-----------------------------------

      end module module_radiation_clouds !

