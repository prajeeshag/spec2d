program main

use iso_c_binding

use diag_manager_mod, only : diag_manager_init, DIAG_OCEAN, DIAG_OTHER, DIAG_ALL
use diag_data_mod, only : files, num_files, mix_snapshot_average_fields, &
    VERY_LARGE_FILE_FREQ, output_fields, num_output_fields
use diag_util_mod, only : sync_file_times, diag_time_inc, get_time_string
use time_manager_mod
use fms_io_mod, only : open_file, close_file, file_exist
use mpp_mod, only : mpp_init, mpp_exit, mpp_error, FATAL, WARNING, NOTE

implicit none

interface 
    integer(C_INT) function nccp2r(narg,args) bind(C,name="nccp2r")
        import
        integer(C_INT), value :: narg
        type(C_PTR), dimension(*) :: args
    end function nccp2r
end interface

type filenm_type
    character(len=128), allocatable :: nm(:)
    integer :: n
end type 

type(C_PTR) :: argsptr(10)
character(len=512) :: args(10)=char(0)
integer :: nargs=2, i, ierr, unit=15, n, position
integer :: calendar_type, startdate(6)=0, enddate(6)=0
real :: totalwait=3600., swait=15., dt
logical :: fexist
character(len=1024) :: msg, suffix, filename, base_name
type(filenm_type), allocatable :: filenms(:)
type(time_type) :: starttime, endtime

dt = 0.
do while (dt <= totalwait)
    fexist=file_exist('run.timestamp')
    if (fexist) exit
    print '(A,F10.0,A)', 'waiting for file run.timestamp since ', dt, 'seconds'
    call wait_seconds(swait)
    dt = dt + swait
end do

if (.not.fexist) call mpp_error(FATAL,"run.timestamp does not exist after waiting")

open(unit,file='run.timestamp',status='old')
read(unit,*) calendar_type
read(unit,*) startdate
read(unit,*) enddate
close(unit,status='delete')

call set_calendar_type(calendar_type)
starttime=set_date(startdate(1),startdate(2),startdate(3), &
                  startdate(4),startdate(5),startdate(6))
endtime=set_date(enddate(1),enddate(2),enddate(3), &
                  enddate(4),enddate(5),enddate(6))

call print_date(starttime,'run_mppnccp2r: STARTTIME')
call print_date(endtime,'run_mppnccp2r: ENDTIME')

call diag_manager_init(DIAG_ALL)

if (num_files<=0) then
    call mpp_error(NOTE, "run_mppnccp2r: no files to process, exiting")
    call mpp_exit()
    stop 0
end if

allocate(filenms(num_files))

print*,'num_files=', num_files

do n = 1, num_files
    print *, trim(files(n)%name)
    CALL sync_file_times(n, starttime, err_msg=msg)
    if (trim(msg)/='') call mpp_error(FATAL,trim(msg))
    call set_filenames(n)
end do

do i = 1, 10
    argsptr(i) = C_LOC(args(i))
end do

args(1) = "nccp2r"//trim(char(0))

do n = 1, 2
    args(2) = trim(filenms(1)%nm(n))//trim(char(0))
    ierr = nccp2r(nargs,argsptr)
    if (ierr/=0) stop 1
end do


contains

subroutine set_filenames(n)
    integer, intent(in) :: n
    integer :: no, nfiles
    real :: nfilesr
    logical :: mid
    type(time_type) :: dt, time, dt_out, next_output, last_output, middle_time

    mid = .false.

    dt = set_time(0)
    dt = dt - diag_time_inc(dt, files(n)%new_file_freq, files(n)%new_file_freq_units)
    !call print_time(dt,"dt for "//trim(files(n)%name))

    dt_out = set_time(0)
    dt_out = dt_out - diag_time_inc(dt_out, files(n)%output_freq, files(n)%output_units)
    !call print_time(dt_out,"dt_out for "//trim(files(n)%name))

    nfilesr = (endtime-starttime)/dt
    nfiles = ceiling(nfilesr)
    !print *, "number of file for "//trim(files(n)%name)//" = ", nfiles
    nfiles = nfiles + 1

    allocate(filenms(n)%nm(nfiles))

    last_output = starttime
    next_output = starttime
    next_output = diag_time_inc(next_output, files(n)%output_freq, files(n)%output_units)

    do no = 1, num_output_fields
        if (output_fields(no)%output_file==n.and. &
                .not.output_fields(no)%static.and. &
                output_fields(no)%time_ops) then
            mid = .not.mix_snapshot_average_fields
            exit
        endif
    enddo

    IF (mid) THEN
       middle_time = (last_output+next_output)/2
    ELSE
       middle_time = next_output
    END IF

    base_name=files(n)%name
    IF ( files(n)%new_file_freq < VERY_LARGE_FILE_FREQ ) THEN
       position = INDEX(files(n)%name, '%')
       IF ( position > 0 )  THEN
          base_name = base_name(1:position-1)
       ELSE
          CALL mpp_error('run_mppnccp2r',&
               & 'file name '//TRIM(files(n)%name)// &
                ' does not contain % for time stamp string', FATAL) 
       END IF
       suffix = get_time_string(files(n)%name, middle_time)
    ELSE
       suffix = ''
    END IF

    nfiles = 1
    filenms(n)%nm(nfiles) = trim(base_name)//trim(suffix)//'.nc'

    time = starttime
    do while(time<endtime)
        time = time + dt_out
        last_output = next_output
        next_output = next_output + dt_out
    
        if (time >= files(n)%next_open) then
            IF (mid) THEN
               middle_time = (last_output+next_output)/2
            ELSE
               middle_time = next_output
            END IF

            base_name=files(n)%name
            IF ( files(n)%new_file_freq < VERY_LARGE_FILE_FREQ ) THEN
               position = INDEX(files(n)%name, '%')
               IF ( position > 0 )  THEN
                  base_name = base_name(1:position-1)
               ELSE
                  CALL mpp_error('run_mppnccp2r',&
                       & 'file name '//TRIM(files(n)%name)// &
                        ' does not contain % for time stamp string', FATAL) 
               END IF
               suffix = get_time_string(files(n)%name, middle_time)
            ELSE
               suffix = ''
            END IF
            nfiles = nfiles + 1
            filenms(n)%nm(nfiles) = trim(base_name)//trim(suffix)//'.nc'
            files(n)%next_open = files(n)%next_open + dt
        end if
    end do

    filenms(n)%n = nfiles

    return

end subroutine set_filenames


subroutine wait_seconds(seconds)
    real, intent(in) :: seconds
    real*8 :: rWait, rDT
    Integer :: iStart, iNew

    ! rWait: seconds that you want to wait for; 
    rWait = seconds; rDT = 0.d0
    call system_clock (iStart)
    do while (rDT <= rWait)
        call system_clock (iNew)
        rDT = floatj (iNew - iStart) / 10000.d0
    enddo
end subroutine wait_seconds


end program main

