program main

use mpp_mod, only : mpp_init, FATAL, mpp_error
use fms_mod, only : fms_init

use mpp_domains_mod, only : mpp_define_domains, domain2d

use transforms_mod, only : init_transforms

use sigio_module
use sigio_r_module

implicit none


integer :: nft=12, iret

character(len=128) :: cfile='sig_ini'

type(sigio_head) :: head
type(sigio_dbti) :: dati

real, allocatable, target :: spec(:)

integer :: lnt2 = 4032

type(domain2d) :: domain

integer :: nlat = 94, nlon = 192, trunc = 62, nwaves 

call mpp_init()
call fms_init()

call mpp_define_domains([1,nlat,1,nlon],[1,1],domain,kxy=1,ishuff=0)

call init_transforms(domain,trunc,nwaves,Tshuffle=.false.)

call sigio_rropen(nft,trim(cfile),iret)

if (iret/=0) stop 'Error: could not open file'

call sigio_alhead(head,iret)

if (iret/=0) stop 'Error: could not allocate header'

call sigio_rrhead(nft,head,iret)

if (iret/=0) stop 'Error: could not read header'

allocate(spec(lnt2))

dati%i = 1
dati%f => spec

call sigio_rrdbti(nft,head,dati,iret)
if (iret/=0) stop 'Error: could not read data'

print *, shape(dati%f)

end program main
