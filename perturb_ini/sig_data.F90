module sig_data_mod

use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_init, FATAL, mpp_error, mpp_exit
use fms_mod, only : fms_init
use fms_io_mod, only : fms_io_exit, write_data

use mpp_domains_mod, only : mpp_define_domains, domain2d

use transforms_mod, only : init_transforms, get_wdecompa, spherical_to_grid, grid_to_spherical

use sigio_module
use sigio_r_module

implicit none

private

integer :: nft=1, iret

real, allocatable, target :: buff(:)
real, allocatable, target :: buff2(:,:)

type(domain2d) :: domain

integer :: nlat, nlon, trunc, nwaves, nwavesglob, nlev, nvar, ntrac, lnt2

integer, allocatable :: wdom(:,:), wdoma(:,:)

character (len=128) :: flnm='out.nc', fldnm='fld'

real, allocatable :: buffg(:,:,:)

complex, allocatable :: buffc(:,:,:)

integer :: i, j

logical :: init = .false.

public :: get_grid_data, write_grid_data_to_sig

contains

subroutine get_grid_data(sigfile,datp)
    character(len=*) :: sigfile
    real, pointer :: datp(:,:,:)
    type(sigio_head) :: head
    type(sigio_dbti) :: dati

    call mpp_init()
    call fms_init()

    nft = nft + 1

    call sigio_rropen(nft,trim(sigfile),iret)

    if (iret/=0) call mpp_error(FATAL, 'Error: could not open file :'//trim(sigfile))
    
    call sigio_alhead(head,iret)
    
    if (iret/=0) call mpp_error(FATAL, 'Error: could not allocate header '//trim(sigfile))
    
    call sigio_rrhead(nft,head,iret)

    if (iret/=0) call mpp_error(FATAL, 'Error: could not read header '//trim(sigfile))
    
    if (.not. init) then 
        nlat = head%latf
        nlon = head%lonf
        trunc = head%jcap
        nlev = head%levs
        ntrac = head%ntrac

        nvar = 1 + 3 * nlev + ntrac * nlev

        call mpp_define_domains([1,nlat,1,nlon],[1,1],domain,kxy=1)
        
        call init_transforms(domain,trunc,nwaves,Tshuffle=.false.)
    
        allocate(wdom(nwaves,2))
    
        call get_wdecompa(wdom,nwavesglob)
    
        lnt2 = (trunc+2)*(trunc+1)
        
        allocate(buff(lnt2))
        allocate(buff2(2,lnt2/2))
        allocate(buffc(nvar,nwaves,2))
        allocate(buffg(nvar,nlat,nlon))
        init = .true.
    endif
    
    do i = 1, nvar
        dati%i = 1 + i
        dati%f => buff 
        
        call sigio_rrdbti(nft,head,dati,iret)
        if (iret/=0) call mpp_error(FATAL, 'Error: could not read data '//trim(sigfile))
        
        buff2(:,:) = reshape(buff,[2,lnt2/2])
    
        do j = 1, nwaves 
            if(wdom(j,1)>0) buffc(i,j,1) = cmplx(buff2(1,wdom(j,1)),buff2(2,wdom(j,1)))
            if(wdom(j,2)>0) buffc(i,j,2) = cmplx(buff2(1,wdom(j,2)),buff2(2,wdom(j,2)))
        enddo
    
    enddo
    
    call sigio_rclose(nft,iret)

    if (iret/=0) call mpp_error(FATAL, 'Error: could not close file '//trim(sigfile))

    call spherical_to_grid(buffc,grid=buffg)

    if (associated(datp)) deallocate(datp) 
    allocate(datp(nvar,nlat,nlon))
    datp = buffg

end subroutine get_grid_data


subroutine write_grid_data_to_sig(sigfile,ofile,datg)
    character(len=*) :: sigfile
    character(len=*) :: ofile
    real, intent(in) :: datg(:,:,:)
    type(sigio_head) :: head
    type(sigio_dbti) :: dati

    integer, parameter :: ifid = 12, ofid = 13

    if (.not.init) call mpp_error(FATAL,'write_grid_data_to_sig: sig_data_mod: Not initialized')

    call sigio_rropen(ifid,trim(sigfile),iret)

    if (iret/=0) call mpp_error(FATAL, 'Error: could not open file :'//trim(sigfile))
    
    call sigio_alhead(head,iret)
    
    if (iret/=0) call mpp_error(FATAL, 'Error: could not allocate header '//trim(sigfile))
    
    call sigio_rrhead(ifid,head,iret)

    if (iret/=0) call mpp_error(FATAL, 'Error: could not read header '//trim(sigfile))


    call sigio_rwopen(ofid,trim(ofile),iret) 

    if (iret/=0) call mpp_error(FATAL, 'Error: could not open file from writing '//trim(ofile))

    call sigio_rwhead(ofid,head,iret)
    if (iret/=0) call mpp_error(FATAL, 'Error: could not write header '//trim(ofile))

    dati%i = 1
    dati%f => buff 
    call sigio_rrdbti(ifid,head,dati,iret)
    if (iret/=0) call mpp_error(FATAL, 'Error: could not read data '//trim(sigfile))
    
    call sigio_rwdbti(ofid,head,dati,iret)    
    if (iret/=0) call mpp_error(FATAL, 'Error: could not write data '//trim(ofile))

    call grid_to_spherical(datg,buffc) 
   
    do i = 1, size(datg,1) 
        do j = 1, nwaves
            if (wdom(j,1)>0) then
                buff2(1,wdom(j,1)) = real(buffc(i,j,1))
                buff2(2,wdom(j,1)) = aimag(buffc(i,j,1))
            endif
            if (wdom(j,2)>0) then
                buff2(1,wdom(j,2)) = real(buffc(i,j,2))
                buff2(2,wdom(j,2)) = aimag(buffc(i,j,2))
            endif
        enddo
        buff = reshape(buff2(:,:),[lnt2])
        
        dati%i = 1 + i
        dati%f => buff 
    
        call sigio_rwdbti(ofid,head,dati,iret)    
        if (iret/=0) call mpp_error(FATAL, 'Error: could not write data '//trim(ofile))
    enddo
    
end subroutine write_grid_data_to_sig
    
end module sig_data_mod
