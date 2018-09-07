function int2str_(n)
    integer(kind=kind_), intent(in) :: n
    character(len=iclen) :: int2str_
    character(len=iclen) :: str

    write(str, *) n
    int2str_ = adjustl(str)

    return

end function int2str_

