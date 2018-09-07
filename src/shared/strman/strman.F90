module strman_mod
implicit none
private

public :: int2str

integer, parameter :: inttype = 1234
integer, parameter :: iclen = digits(inttype)+1

interface int2str
	module procedure int2str_4, int2str_8
end interface int2str 

contains

#include <strman.inc>

end module strman_mod

#ifdef test_strman
program main

use strman_mod

print *, int2str(234)

end program main
#endif
