program main

use, intrinsic :: iso_c_binding

use mpp_mod, only : mpp_init, FATAL, mpp_error, mpp_exit
use fms_mod, only : fms_init
use fms_io_mod, only : fms_io_exit, write_data

use mpp_domains_mod, only : mpp_define_domains, domain2d

use transforms_mod, only : init_transforms, get_wdecompa, spherical_to_grid

use sigio_module
use sigio_r_module

implicit none


integer :: nft=12, iret

character(len=128) :: cfile='sig_ini'

type(sigio_head) :: head
type(sigio_dbti) :: dati

real, allocatable, target :: buff(:)
real, allocatable, target :: buff2(:,:,:)

integer :: lnt2 = 4032

integer, parameter :: nlev = 64

integer, parameter :: nvar = nlev * 4 + 1

type(domain2d) :: domain

integer :: nlat = 94, nlon = 192, trunc = 62, nwaves, nwavesglob

integer, allocatable :: wdom(:,:), wdoma(:,:)

character (len=128) :: flnm='out.nc', fldnm='fld'

real, allocatable :: buffg(:,:,:)

complex, allocatable :: buffc(:,:,:)

integer :: i, j


call mpp_init()
call fms_init()

call mpp_define_domains([1,nlat,1,nlon],[1,1],domain,kxy=1)

call init_transforms(domain,trunc,nwaves,Tshuffle=.false.)

print *, nwaves

allocate(wdom(nwaves,2))

call get_wdecompa(wdom,nwavesglob)

print *, 'even = ', wdom(:,1)
print *, 'odd = ', wdom(:,2)

call getarg(1,cfile)

flnm = trim(cfile)//".nc"

lnt2 = (trunc+2)*(trunc+1)

print *, 'lnt2, nwavesglob =', lnt2, nwavesglob

call sigio_rropen(nft,trim(cfile),iret)

if (iret/=0) stop 'Error: could not open file'

call sigio_alhead(head,iret)

if (iret/=0) stop 'Error: could not allocate header'

call sigio_rrhead(nft,head,iret)

if (iret/=0) stop 'Error: could not read header'

allocate(buff(lnt2))
allocate(buff2(2,lnt2/2,nvar))
allocate(buffc(nvar,nwaves,2))
allocate(buffg(nvar,nlat,nlon))

do i = 1, nvar

    dati%i = 1 + i
    dati%f => buff 
    
    call sigio_rrdbti(nft,head,dati,iret)
    if (iret/=0) stop 'Error: could not read data'
    
    buff2(:,:,i) = reshape(buff,[2,lnt2/2])

    do j = 1, nwaves 
        if(wdom(j,1)>0) buffc(i,j,1) = cmplx(buff2(1,wdom(j,1),i),buff2(2,wdom(j,1),i))
        if(wdom(j,2)>0) buffc(i,j,2) = cmplx(buff2(1,wdom(j,2),i),buff2(2,wdom(j,2),i))
    enddo

enddo

call spherical_to_grid(buffc,grid=buffg)

call write_data(flnm,fldnm,buff2)

call write_data(flnm,trim(fldnm)//'g',buffg)

call fms_io_exit()
call mpp_exit()

end program main
