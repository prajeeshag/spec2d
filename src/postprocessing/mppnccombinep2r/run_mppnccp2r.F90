program main

use iso_c_binding

use diag_manager_mod, only : diag_manager_init, DIAG_OCEAN, DIAG_OTHER, DIAG_ALL
use diag_data_mod, only : files, num_files, mix_snapshot_average_fields, &
    VERY_LARGE_FILE_FREQ, output_fields, num_output_fields
use diag_util_mod, only : sync_file_times, diag_time_inc, get_time_string
use time_manager_mod
use mpp_mod, only : mpp_init, mpp_exit, mpp_error, FATAL, WARNING, NOTE

implicit none

interface 
    integer(C_INT) function nccp2r(narg,args) bind(C,name="nccp2r")
        import
        integer(C_INT), value :: narg
        type(C_PTR), dimension(*) :: args
    end function nccp2r
end interface

integer, parameter :: maxarg=32

type filenm_type
    character(len=128), allocatable :: nm(:)
    integer :: total, done
    logical :: lowfreq=.false.
end type

integer :: maxrunlength=100 !years

type(C_PTR) :: args(maxarg)

integer :: nargs=1, i, ierr, unit=15, n, position, nf, stat
integer :: calendar_type=-11, startdate(6)=0, enddate(6)=0
logical :: fexist
character(len=1024) :: msg, fnm, fnm_next
type(filenm_type), allocatable :: filenms(:)
type(time_type) :: starttime, endtime

character(len=8) :: arg_r="-r"//char(0)
character(len=8) :: arg_v="-r"//char(0)
character(len=8) :: arg_vv="-vv"//char(0)
character(len=8) :: arg_n="-n"//char(0)
character(len=8) :: arg_n4="-n4"//char(0)
character(len=8) :: arg_u="-u"//char(0)
character(len=8) :: arg_ov="-ov"//char(0)

logical :: removein=.true.
integer :: startpe=0, nc4=4, run=2, atmpes=1, ocnpes=1
character(len=32) :: prgrm="nccp2r"//char(0)
character(len=1024) :: xgrid="INPUT/p_xgrd.nc", run_time_stamp='INPUT/atm.res'
character(len=64) :: cnc4, cstartpe
type(time_type) :: lowestfreq
real :: time1, time2, endwaittime=0., minendwaittime=10.
logical :: next_file_found=.false., end_check=.false.

namelist/opts_nml/removein, atmpes, ocnpes, nc4, xgrid, run, startpe, &
                  minendwaittime, startdate, calendar_type

call cpu_time(time1)

call mpp_init()

read(*,nml=opts_nml,iostat=stat)
if (stat/=0) call mpp_error(FATAL,"run_mppnccp2r: error while reading of options")

xgrid = trim(xgrid)//char(0)

nargs=1

args(nargs) = c_loc(prgrm)
nargs=nargs+1

if (removein) then
    args(nargs) = c_loc(arg_r)
    nargs=nargs+1
endif

if (startpe>0) then
    args(nargs) = c_loc(arg_n)
    nargs=nargs+1
    write(cstartpe,*) startpe 
    cstartpe=trim(adjustl(cstartpe))//char(0)
    args(nargs) = c_loc(cstartpe)
    nargs=nargs+1
endif 

if (nc4>0) then
    args(nargs) = c_loc(arg_n4)
    nargs=nargs+1
    write(cnc4,*) nc4 
    cnc4=trim(adjustl(cnc4))//char(0)
    args(nargs) = c_loc(cnc4)
    nargs=nargs+1
endif


args(nargs) = c_loc(arg_u)
nargs=nargs+1
args(nargs) = c_loc(xgrid)
nargs=nargs+1

if (all(startdate<=0).or.calendar_type==-11) then
    call mpp_error(NOTE,"startdate or calendar_type is not provided correctly "// &
                    "via STDIN, trying to read from run_time_stamp")
    fexist=file_exist(run_time_stamp)
    if (.not.fexist) call mpp_error(FATAL,"run_mppnccp2r: run_time_stamp file ("// &
                trim(run_time_stamp)//") do not exist")

    open(unit,file=trim(run_time_stamp),status='old')
    read(unit,*) calendar_type
    read(unit,*) startdate
    close(unit)
endif

call set_calendar_type(calendar_type)

starttime=set_date(startdate(1),startdate(2),startdate(3), &
                  startdate(4),startdate(5),startdate(6))

endtime=increment_date(starttime, years = maxrunlength)

lowestfreq = set_time(0) !lowest frequency of output

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

args(nargs) = c_loc(fnm)
end_check = .false.

do while (.true.) ! Infinite Loop

    next_file_found=.false.

    do nf = 1, num_files
        n = filenms(nf)%done+1
        fnm_next = trim(filenms(nf)%nm(n+1))
        if (.not.all_files_exist(trim(fnm_next),startpe,atmpes)) cycle

        next_file_found=.true. !Atleast one next file was found
        if (endwaittime<=0.) then
            if (filenms(nf)%lowfreq) then 
                !Find the lowest frequency output, set its frequency as
                !endwaittime in seconds
                call cpu_time(time2)
                endwaittime=max(time2-time1,minendwaittime) 
                write(msg,*) endwaittime
                call mpp_error(NOTE,"run_mppnccp2r: end waiting time is "//trim(adjustl(msg)))
            endif
        endif
            
        fnm = trim(filenms(nf)%nm(n))//char(0)
        ierr = nccp2r(nargs,args)
        if (ierr/=0) call mpp_error(FATAL,"nccpr failed for file "//trim(fnm))

        call mpp_error(NOTE,trim(fnm)//" done...")
        filenms(nf)%done=n
    end do

    if (end_check) then
        if (.not.next_file_found) then
            call mpp_error(NOTE,"No next file found after waiting time, exiting...")
            exit
        else
            end_check=.false.
        endif
    endif
             
    if (endwaittime>0.and..not.next_file_found) then
        call mpp_error(NOTE,"Found no next file for all files, waiting for sometime")
        call wait_seconds(endwaittime)
        end_check=.true.
    endif
    
end do

call mpp_error(NOTE,"run_mppnccp2r: Processing last files...")

do nf = 1, num_files
    n = filenms(nf)%done+1
    fnm = trim(filenms(nf)%nm(n))//char(0)
    ierr = nccp2r(nargs,args)
    if (ierr/=0) call mpp_error(FATAL,"nccpr failed for file "//trim(fnm))
    call mpp_error(NOTE,trim(fnm)//" done...")
    filenms(nf)%done=n
end do

call mpp_error(NOTE,"run_mppnccp2r: DONE")

contains

logical function all_files_exist(fnm,strt,cnt)
    character(len=*), intent(in) :: fnm
    integer, intent(in) :: strt, cnt
    character(len=len(fnm)+10) :: filnm
    character(len=10) :: suffix
    integer :: i
    
    all_files_exist = .true.

    do i = strt, strt+cnt-1 
        write(suffix,'(I4.4)') i
        suffix = "."//trim(adjustl(suffix))
        if(.not.file_exist(trim(fnm)//trim(suffix))) then
            all_files_exist=.false.
            return
        endif
    enddo
       
    return 
end function all_files_exist


subroutine set_filenames(n)
    integer, intent(in) :: n
    integer :: no, nfiles
    real :: nfilesr
    logical :: mid
    character(len=1024) :: base_name, suffix
    type(time_type) :: dt, time, dt_out, next_output, last_output, middle_time

    mid = .false.

    filenms(n)%lowfreq=.false.    

    dt = set_time(0)
    dt = dt - diag_time_inc(dt, files(n)%new_file_freq, files(n)%new_file_freq_units)
    if (dt>lowestfreq.and.files(n)%new_file_freq<VERY_LARGE_FILE_FREQ) then
        filenms(n)%lowfreq=.true.
        lowestfreq = dt
    endif
    !call print_time(dt,"dt for "//trim(files(n)%name))

    dt_out = set_time(0)
    dt_out = dt_out - diag_time_inc(dt_out, files(n)%output_freq, files(n)%output_units)
    !call print_time(dt_out,"dt_out for "//trim(files(n)%name))

    nfilesr = (endtime-starttime)/dt
    nfiles = ceiling(nfilesr)
    if (nfiles<=1) then
        filenms(n)%total = nfiles-1
        filenms(n)%done = 0
    endif
        
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

    filenms(n)%total = nfiles-1
    filenms(n)%done = 0

    return

end subroutine set_filenames


subroutine wait_seconds(seconds)
    real, intent(in) :: seconds
    real :: rWait, rDT
    integer :: iStart, iNew, count_rate
    integer :: sec1
    character(len=10) :: dtime
    sec1=0
    ! rWait: seconds that you want to wait for; 
    rWait = seconds; rDT = 0.d0
    call system_clock(iStart)
    call date_and_time(TIME=dtime)
    print *, dtime
    do while (rDT <= rWait)
        call system_clock(iNew,count_rate)
        rDT = float(iNew - iStart)/count_rate
    enddo
    call date_and_time(TIME=dtime)
    print *, dtime
end subroutine wait_seconds

logical function file_exist(flnm)
    character(len=*) :: flnm
    
    inquire(file=trim(flnm),exist=file_exist)

    return

end function file_exist


end program main

