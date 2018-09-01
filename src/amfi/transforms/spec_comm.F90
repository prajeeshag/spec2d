module spec_comm_mod

use mpp_mod, only : mpp_error, fatal, warning, note, mpp_init, mpp_npes, mpp_pe, &
        mpp_declare_pelist, mpp_get_current_pelist, mpp_set_current_pelist, mpp_root_pe, &
        mpp_gather, mpp_broadcast

use strman_mod, only : int2str

implicit none

include 'mpif.h'

private

public :: split_pelist, spec_comm_allgather, spec_comm_max, spec_comm_sumI, spec_comm_sumDC, &
          spec_comm_npes, spec_comm_pe, spec_comm_bcast

interface spec_comm_sum
    module procedure spec_comm_sumDC
    module procedure spec_comm_sumI
end interface spec_comm_sum

contains

subroutine spec_comm_allgather(buffin, buffout, commid)
    integer, intent(in) :: buffin(:)
    integer, intent(out) :: buffout(:)
    integer, intent(in), optional :: commid

    integer :: comm_id, n_in, n_out, ierr, npes

    if (present(commid)) then
        comm_id = commid
    else
        call mpp_get_current_pelist(commid=comm_id)
    endif
     
    npes = spec_comm_npes(comm_id)
    n_in = size(buffin)
    n_out = size(buffin)*npes
    
    if (size(buffout)/=n_out) call mpp_error(FATAL, 'spec_comm_allgather: buffout size error')
    
    call MPI_ALLGATHER(buffin, n_in, MPI_INTEGER, buffout, n_in, MPI_INTEGER, comm_id, ierr)

    if (ierr/=MPI_SUCCESS) call mpp_error(FATAL,'spec_comm_allgather: error '//int2str(ierr))

    return

end subroutine spec_comm_allgather


subroutine split_pelist(iamin, spes, snpes, commid)
    logical, intent(in) :: iamin
    integer, intent(inout) :: spes(:)
    integer, intent(inout) :: snpes
    integer, intent(inout), optional :: commid

    integer, allocatable :: allpes(:), parentpes(:), pes(:)
    integer :: pe, i, j, npes
    integer :: commidd

    pe = -1
    commidd = 0

    allocate(allpes(mpp_npes()),parentpes(mpp_npes()))

    call mpp_get_current_pelist(pelist=parentpes)

    if (iamin) pe = mpp_pe()

    call mpp_gather([pe], allpes)
    call mpp_broadcast(allpes,size(allpes),mpp_root_pe())
     
    npes = count(allpes>-1)
    allocate(pes(npes))

    pes = pack(allpes,allpes>-1)

    call mpp_declare_pelist(pes(1:npes))

    if (any(mpp_pe()==pes(1:npes))) then
        call mpp_set_current_pelist(pes,no_sync=.true.)
        call mpp_get_current_pelist(pelist=pes, commid=commidd)
        if (size(spes)<npes) call mpp_error(FATAL,'split_pelist: size of spes < the required size.')
        if (present(commid)) commid = commidd
        snpes = npes
        spes(1:npes) = pes(:) 
        call mpp_set_current_pelist(pelist=parentpes, no_sync=.true.)
    end if

    deallocate(allpes,parentpes)

    return
end subroutine split_pelist

integer function spec_comm_npes(commid)
    integer, intent(in) :: commid
    integer :: ierr, rank

    call MPI_COMM_SIZE(commid, rank, ierr)
    if (ierr/=MPI_SUCCESS) call mpp_error(FATAL,'spec_comm_npes: error '//int2str(ierr))
    spec_comm_npes = rank
    return
end function spec_comm_npes

integer function spec_comm_pe(commid)
    integer, intent(in) :: commid
    integer :: ierr, rank

    call MPI_COMM_RANK(commid, rank, ierr)
    if (ierr/=MPI_SUCCESS) call mpp_error(FATAL,'spec_comm_pe: error '//int2str(ierr))
    spec_comm_pe = rank
    return
end function spec_comm_pe

subroutine spec_comm_bcast(buff, n, pe, comm_id)
    integer, intent(in) :: n
    real, intent(inout) :: buff(n)
    integer, intent(in) :: pe
    integer, intent(in) :: comm_id

    integer :: ierr

    call MPI_BCAST(buff, n, MPI_REAL8, pe, comm_id, ierr)
    if (ierr/=MPI_SUCCESS) call mpp_error(FATAL,'spec_comm_bcast: error '//int2str(ierr))

    return
end subroutine spec_comm_bcast

!--------------------------------------------------------------------------------   
subroutine spec_comm_max(arr, n, comm_id)
!--------------------------------------------------------------------------------   
    integer, intent(in) :: n
    real(kind=8), intent(inout) :: arr(n) 
    integer, intent(in), optional :: comm_id
    real(kind=8) :: buff(n) 
    integer :: ierr, commID

    if (present(comm_id)) then
        commID = comm_id
    else
        call mpp_get_current_pelist(commid=commID)
    endif
    call MPI_ALLREDUCE(arr, buff, n, MPI_REAL8, MPI_MAX, commID, ierr)
    if (ierr/=MPI_SUCCESS) call mpp_error(FATAL,'spec_comm_max: error '//int2str(ierr))

    arr = buff

    return
end subroutine spec_comm_max

subroutine spec_comm_sumDC(arr, n, comm_id)
    integer, intent(in) :: n
    complex, intent(inout) :: arr(n) 
    integer, intent(in), optional :: comm_id
    complex :: buff(n) 
    integer :: ierr, commID

    if (present(comm_id)) then
        commID = comm_id
    else
        call mpp_get_current_pelist(commid=commID)
    endif

    call MPI_ALLREDUCE(arr, buff, n, MPI_DOUBLE_COMPLEX, MPI_SUM, commID, ierr)

    if (ierr/=MPI_SUCCESS) call mpp_error(FATAL,'spec_comm_sumDC: error '//int2str(ierr))

    arr = buff

    return
end subroutine spec_comm_sumDC

subroutine spec_comm_sumI(arr,n, comm_id)
    integer, intent(in) :: n
    integer, intent(inout) :: arr(n) 
    integer, intent(in), optional :: comm_id
    integer :: buff(n) 
    integer :: ierr, commID

    if (present(comm_id)) then
        commID = comm_id
    else
        call mpp_get_current_pelist(commid=commID)
    endif
    call MPI_ALLREDUCE(arr, buff, n, MPI_INTEGER, MPI_SUM, commID, ierr)
    if (ierr/=MPI_SUCCESS) call mpp_error(FATAL,'spec_comm_sum: error '//int2str(ierr))

    arr = buff

    return
end subroutine spec_comm_sumI

end  module spec_comm_mod

