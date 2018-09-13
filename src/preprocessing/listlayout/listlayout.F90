program listlayout

use ocpack_mod

implicit none

type lstlay_type
    integer :: npes=0
    integer, allocatable :: xnpe(:), ynpe(:)
    integer :: npair=0
end type lstlay_type

integer :: nlat, i, npes, j, trunc
character(len=1024) :: msg
type(lstlay_type), allocatable :: lstlay(:)
integer :: nlstlay=0, n

nlstlay=0

write(*,*) "enter number of latitudes, truncation : "
read(*,*) nlat, trunc
write(*,*) nlat, trunc

if (nlat<=0) then
    stop "Error in number of latitudes"
endif

if (trunc<=0) then
    stop "Error in truncation"
endif

call init_ocpack(nlat,trunc,[1,1])

allocate(lstlay(nxpe*nype))

write(*,*)
write(*,*) 'valid npes_y valued are : ', npesy(1:nype)%npes
write(*,*)

write(*,*)
write(*,*) 'valid npes_x valued are : ', npesx(1:nxpe)%npes
write(*,*)



do j = 1, nype
    do i = 1, nxpe
        if (npesx(i)%npes<npesy(j)%npes) cycle
        call assign_lst(npesx(i)%npes, npesy(j)%npes)
    end do
end do

write(msg,*)nlat 
msg=trim(adjustl(msg))

call print_layout()

contains

subroutine assign_lst(npex, npey)
    integer, intent(in) :: npex, npey
    integer :: npexy, n, i

    npexy = npex*npey

    i = 0
    do n = 1, nlstlay
        if (npexy==lstlay(n)%npes) then
            i = n
            exit
        endif
    enddo

    if (i==0) then
        nlstlay = nlstlay+1
        i = nlstlay
        lstlay(i)%npes = npexy
        allocate(lstlay(i)%xnpe(npexy), lstlay(i)%ynpe(npexy))
        lstlay(i)%xnpe=-1; lstlay(i)%ynpe=-1
        lstlay(i)%npair = 0
    endif

    lstlay(i)%npair = lstlay(i)%npair+1

    lstlay(i)%xnpe(lstlay(i)%npair) = npex

    lstlay(i)%ynpe(lstlay(i)%npair) = npey
end subroutine

integer function ndigits(nn)
    integer, intent(in) :: nn
    integer :: n

    n = nn
    
    ndigits = 0

    do while(n >= 0)
        n = n/10
        ndigits = ndigits + 1
    end do

end function ndigits       

subroutine print_layout()
    integer :: pmin, loc
    character(len=1024) :: cly1, cly2
    integer, allocatable :: npes(:)

    allocate(npes(nlstlay))

    npes = lstlay(1:nlstlay)%npes

    pmin = 0

    write(*,*) "Its highly recommended to use a layout where npes_x >= npes_y"
    write(*,*) "All such valid processor layouts are listed below:"
    do n = 1, nlstlay
        loc = minloc(npes,1,mask=npes>pmin)
        pmin = lstlay(loc)%npes

        cly1 = ""; cly2 = "" 

        write(cly1,*) lstlay(loc)%npes, ' pes:-'
        cly2 = trim(cly2)//trim(adjustl(cly1))

        write(*,'(A)', advance='no') trim(cly2)

        do i = 1, lstlay(loc)%npair
             write(*,'(A)', advance='no') " ["
             write(cly1,*) lstlay(loc)%ynpe(i)
             write(*,'(A)', advance='no') trim(adjustl(cly1))//","
             write(cly1,*) lstlay(loc)%xnpe(i)
             write(*,'(A)', advance='no') trim(adjustl(cly1))//"]; "
        end do
        write(*,*)" " 
    end do
    
end subroutine print_layout 

end program listlayout
