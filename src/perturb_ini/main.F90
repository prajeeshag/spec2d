
program main

use sig_data_mod, only : get_grid_data, write_grid_data_to_sig

use mpp_mod, only : mpp_init, mpp_exit
use fms_io_mod, only : write_data, fms_io_exit

use mersenne_twister, only : random_gauss

implicit none

real, pointer :: dat1(:,:,:), dat2(:,:,:), dat3(:,:,:)

real, pointer :: epslon(:,:,:), dati(:,:,:), rand(:,:,:)

character(len=128) :: cfile, cfile1, ofile

integer :: i, k, j, nini, n, nsize

call mpp_init()

call getarg(1,cfile)
print *, trim(cfile)
call get_grid_data(cfile,dat1)
cfile=''

call getarg(2,cfile)
print *, trim(cfile)
cfile1 = cfile
call get_grid_data(cfile,dat2)
cfile=''

call getarg(3,cfile)
print *, trim(cfile)
call get_grid_data(cfile,dat3)
cfile=''

allocate(epslon(size(dat1,1),size(dat1,2),size(dat1,3)))
allocate(dati(size(dat1,1),size(dat1,2),size(dat1,3)))
allocate(rand(size(dat1,1),size(dat1,2),size(dat1,3)))

epslon = dat2 - (dat1+dat3)*0.5

nini = 30
nsize = size(dat1,1)*size(dat1,2)*size(dat1,3)

do n = 1, nini
    call random_gauss_prt(nsize,rand)
    dati = dat2 + rand * epslon
    write(ofile,'(A,I3.3)')trim(cfile1)//'_', n
    call write_grid_data_to_sig(cfile1,ofile,dati)
enddo    
   

deallocate(epslon, dati)
 
call fms_io_exit()
call mpp_exit()

contains 

subroutine random_gauss_prt(ns,randm)
    integer, intent(in) :: ns
    real :: randm(ns)

    call random_gauss(randm)
    return
end subroutine random_gauss_prt

end program main
