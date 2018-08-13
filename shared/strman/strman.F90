module strman_mod
implicit none
private

public :: int2str

integer, parameter :: inttype = 1234
integer, parameter :: iclen = digits(inttype)+1

contains

function int2str(n)
    integer, intent(in) :: n
    character(len=iclen) :: int2str
    character(len=iclen) :: str

    write(str, *) n
    int2str = adjustl(str)

    return

end function int2str

end module strman_mod

#ifdef test_strman
program main

use strman_mod

print *, int2str(234)

end program main
#endif
