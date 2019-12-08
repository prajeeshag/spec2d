program main

use iso_c_binding

use diag_manager_mod, only : diag_manager_init, DIAG_OCEAN, DIAG_OTHER, DIAG_ALL
use diag_data_mod, only : files, num_files, mix_snapshot_average_fields, &
    VERY_LARGE_FILE_FREQ, output_fields, num_output_fields
use diag_util_mod, only : sync_file_times, diag_time_inc, get_time_string
use time_manager_mod
use mpp_mod, only : mpp_init, mpp_exit, mpp_error, FATAL, WARNING, NOTE, &
        mpp_pe, mpp_root_pe, mpp_npes

implicit none

#ifdef use_libMPI
include 'mpif.h'
#endif

interface 
    integer(C_INT) function nccp2r(narg,args) bind(C,name="nccp2r")
        import
        integer(C_INT), value :: narg
        type(C_PTR), dimension(*) :: args
    end function nccp2r
    
    real(C_DOUBLE) function modtimediff(args) bind(C,name="modtimediff")
        import
        type(C_PTR), dimension(*), intent(in) :: args
    end function modtimediff

    real(C_DOUBLE) function modtime(args) bind(C,name="modtime")
        import
        type(C_PTR), dimension(*), intent(in) :: args
    end function modtime
    
    integer(C_INT) function show_status(args) bind(C,name="show_status")
        import
        real(C_DOUBLE), intent(in), value :: args
    end function show_status

    integer(C_INT) function rmfile(args) bind(C,name="rmfile")
        import
        type(C_PTR), dimension(*), intent(in) :: args
    end function rmfile
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

integer :: removein=0
integer :: startpe=0, nc4=0, atmpes=1, ocnpes=1, tfile=0
character(len=32) :: prgrm="nccp2r"//char(0)
character(len=1024) :: xgrid="INPUT/p_xgrd.nc", run_time_stamp='INPUT/atm.res'
character(len=64) :: cnc4, cstartpe
type(time_type) :: lowestfreq
integer :: lownf=0
real :: time1, time2, endwaittime=0., minendwaittime=30., maxwait=12*3600.
real :: mtime1, mtime2
logical :: next_file_found=.false., end_check=.false.
integer :: child_run=0, ov=1, verbose=0
logical :: no_files_found=.true.

call cpu_time(time1)
call mpp_init()
call read_options()

xgrid = trim(xgrid)//char(0)

nargs=1

args(nargs) = c_loc(prgrm)
nargs=nargs+1

if (removein/=0) then
    call mpp_error(NOTE,"Removein opiton ON")
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

lowestfreq = set_time(0,days=VERY_LARGE_FILE_FREQ) !lowest frequency of output

call print_date(starttime,'run_mppnccp2r: STARTTIME')
call print_date(endtime,'run_mppnccp2r: ENDTIME')

call diag_manager_init(DIAG_ALL)

if (num_files<=0) then
    call mpp_error(NOTE, "run_mppnccp2r: no files to process, exiting")
    call mpp_exit()
    stop 0
end if

allocate(filenms(num_files))

if (mpp_pe()==mpp_root_pe()) print*,'num_files=', num_files

do n = 1, num_files
    if (mpp_pe()==mpp_root_pe()) print *, trim(files(n)%name)
    CALL sync_file_times(n, starttime, err_msg=msg)
    if (trim(msg)/='') call mpp_error(FATAL,trim(msg))
    call set_filenames(n)
end do

if (lownf>0) filenms(lownf)%lowfreq=.true.

args(nargs) = c_loc(fnm)
end_check = .false.

call mpp_error(NOTE,"stage 1")
do nf = 1, num_files
    do n = 1, size(filenms(nf)%nm) 
        if (all_files_exist(trim(filenms(nf)%nm(n)),0,1)) then
            filenms(nf)%done=n-1
            exit
        endif
    end do
end do

if (mpp_pe()==mpp_root_pe()) then
    if (child_run/=0) then !Launched along with the of model run
        call mpp_error(NOTE,"stage 2")
        no_files_found=.true. 
        tfile=0
        mtime1 = 0.
        mtime2 = 0.
        do while (.true.) ! Infinite Loop
    
            next_file_found=.false.
    
            do nf = 1, num_files
    
                n = filenms(nf)%done+1
    
                fnm_next = trim(filenms(nf)%nm(n+1))
    
                if (.not.all_files_exist(trim(fnm_next),0,atmpes)) cycle
    
                next_file_found=.true. !Atleast one next file was found
    
                no_files_found=.false.
    
                if (filenms(nf)%lowfreq.and.endwaittime<=0.) then
                    !Find the lowest frequency output, set its frequency as
                    !endwaittime in seconds
                    if (mtime1<=0.) then
                        mtime1 = mod_time(trim(filenms(nf)%nm(n)),0)
                    elseif (mtime2<=0.) then
                        mtime2 = mod_time(trim(filenms(nf)%nm(n)),0)
                    else
                        endwaittime = (mtime2-mtime1)*3.
                        endwaittime=max(endwaittime,minendwaittime) 
                        write(msg,*) endwaittime
                        call mpp_error(NOTE,"run_mppnccp2r: end waiting time is " &
                                       //trim(adjustl(msg))//" seconds")
                    endif
                endif
                ierr = send_jobs(nf,n)
                if (ierr/=0) call quit_jobs()
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
    
            if (no_files_found) then 
                call cpu_time(time2)
                if (time2-time1 > maxwait) then
                    exit
                endif
            endif
    
        end do
        call mpp_error(NOTE,"run_mppnccp2r: Processing last files...")
    endif
    
    do nf = 1, num_files
        do n = filenms(nf)%done+1, filenms(nf)%total
            if (.not.all_files_exist(trim(filenms(nf)%nm(n)),0,atmpes)) exit
            ierr = send_jobs(nf,n)
            if (ierr/=0) call quit_jobs()
            filenms(nf)%done=n
        end do
    end do

    call all_done()
else
    call do_jobs()
endif

call mpp_error(NOTE,"run_mppnccp2r: DONE")

call mpp_exit()

contains

subroutine read_options()
    integer :: ierr, stat

    namelist/opts_nml/removein, atmpes, ocnpes, nc4, startpe, &
                      minendwaittime, startdate, calendar_type, child_run, &
                      ov, verbose
    if (mpp_pe()==mpp_root_pe()) then
        read(*,nml=opts_nml,iostat=stat)
        if (stat/=0) call mpp_error(FATAL,"run_mppnccp2r: error while reading of options")
    end if
    
    call mpi_bcast(removein, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(atmpes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(nc4, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ov, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(verbose, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(child_run, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(calendar_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(startdate, size(startdate), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

return
end subroutine read_options


integer function submit_processing(nf,n)
    integer, intent(in) :: nf, n
    integer :: ierr

    fnm = trim(filenms(nf)%nm(n))//char(0)
    if (ov/=0) then
        if (rm_file(fnm)==0) then
            print *, "Deleted old file "//trim(filenms(nf)%nm(n))
        endif
    endif
    ierr = nccp2r(nargs,args)
    if (ierr/=0) then
        print *, "ERROR: nccpr failed for file "//trim(fnm)
        submit_processing=-1
        return
    endif
     
    print *, trim(fnm)//" done..."
    submit_processing = 0
    return

end function submit_processing

subroutine do_jobs()
    integer :: jobids(2), done, tag
    INTEGER :: stat(MPI_STATUS_SIZE)
    integer :: nf, n, ierr

    tag = 0

    done = mpp_pe()
    do while (.true.)
        call MPI_SEND(done, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, ierr)
        call MPI_RECV(jobids, 2, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, stat, ierr)
        nf = jobids(1); n = jobids(2)
        if (nf==0) exit
        if(verbose>0)print *, "recieved file "//trim(filenms(nf)%nm(n))//" for processing" 
        ierr = submit_processing(nf,n)
        if (ierr/=0) then
            done = done * -1
        endif
    end do

end subroutine do_jobs

integer function send_jobs(nf,n)
    integer, intent(in) :: nf, n
    integer :: jobids(2), free_pe, tag, ierr
    INTEGER :: stat(MPI_STATUS_SIZE)
    character(len=512) :: msg, fnm

    tag = 0

    if (nf>0) then
        fnm = trim(filenms(nf)%nm(n))//char(0)
        if (ov/=0) then
            if (rm_file(fnm)==0) then
                print *, "Deleted old file "//trim(filenms(nf)%nm(n))
            endif
        endif
    endif

    if (mpp_npes()>1) then
        if (mpp_pe()==mpp_root_pe()) then
            jobids(1) = nf
            jobids(2) = n
            call MPI_RECV(free_pe, 1, MPI_INTEGER, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, stat, ierr)
            write(msg,*) abs(free_pe)
            if (free_pe<0) then
                jobids = 0
                free_pe=abs(free_pe)
                call MPI_SEND(jobids, 2, MPI_INTEGER, free_pe, tag, MPI_COMM_WORLD, ierr)
                send_jobs = -1
                return
            endif
            if(nf>0) then
                if(verbose>0) print *, "sending file "//trim(filenms(nf)%nm(n))//" for processing to pe "//trim(adjustl(msg))
            else
                print *, "sending end job msg to pe "//trim(adjustl(msg))
            endif
            call MPI_SEND(jobids, 2, MPI_INTEGER, free_pe, tag, MPI_COMM_WORLD, ierr)
        endif
    else
        ierr = submit_processing(nf,n)
        if (ierr/=0) then
            send_jobs = -1
            return
        endif
    endif

    send_jobs = 0
   
    return 
end function send_jobs

subroutine quit_jobs()
    integer :: n, ierr
    
    if (mpp_pe()/=mpp_root_pe()) return
    do n = 1, mpp_npes()-2
       ierr = send_jobs(0,0) 
    end do

    call mpp_error(FATAL,"ERROR")
    return
end subroutine quit_jobs

subroutine all_done()
    integer :: n
    
    if (mpp_pe()/=mpp_root_pe()) return
    do n = 1, mpp_npes()-1
       ierr = send_jobs(0,0) 
    end do
    return
end subroutine all_done

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

    base_name=files(n)%name

    dt = set_time(0)
    dt = dt - diag_time_inc(dt, files(n)%new_file_freq, files(n)%new_file_freq_units)
    !call print_time(dt,"dt for "//trim(files(n)%name))

    dt_out = set_time(0)
    dt_out = dt_out - diag_time_inc(dt_out, files(n)%output_freq, files(n)%output_units)
    !call print_time(dt_out,"dt_out for "//trim(files(n)%name))

    nfilesr = (endtime-starttime)/dt
    nfiles = ceiling(nfilesr)

    filenms(n)%lowfreq=.false.

    if (child_run/=0) then
        if (nfiles<3) then
            call mpp_error(NOTE,trim(base_name)//" does not satisfy minimum 3 file criteria, "// &
                                                   "won't process this file")
            filenms(n)%total = 0
            filenms(n)%done = 0
            return
        endif

        if (dt<lowestfreq) then
            lowestfreq = dt
            lownf = n
        endif
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
    integer :: percent1, percent2
    character(len=10) :: dtime

    ! rWait: seconds that you want to wait for; 
    rWait = seconds; rDT = 0.d0
    percent1=0; percent2=0
    call system_clock(iStart)
    print *, dtime
    do while (rDT <= rWait)
        call system_clock(iNew,count_rate)
        rDT = float(iNew - iStart)/count_rate
        percent2 = (rDT/rWait)*100
    enddo
end subroutine wait_seconds

logical function file_exist(flnm)
    character(len=*) :: flnm
    
    inquire(file=trim(flnm),exist=file_exist)

    return

end function file_exist


real function mod_time(file1,spe)
    character(len=*), intent(in) :: file1
    integer, intent(in), optional :: spe
    character(len=len(file1)+10) :: f1
    character(len=10) :: cpe
    type(C_PTR) :: ptr(1)

    cpe = ""

    if (present(spe)) then
        write(cpe,'(I4.4)')spe
        cpe = "."//trim(adjustl(cpe))
    endif

    f1 = trim(file1)//trim(cpe)//char(0)

    ptr(1) = c_loc(f1)

    mod_time = modtime(ptr)

    return 
end function mod_time

real function diff_modtime(file1,file2,spe)
    character(len=*), intent(in) :: file1, file2
    integer, intent(in), optional :: spe
    character(len=len(file1)+10) :: f1
    character(len=len(file2)+10) :: f2
    character(len=10) :: cpe
    type(C_PTR) :: ptr(2)

    cpe = ""

    if (present(spe)) then
        write(cpe,'(I4.4)')spe
        cpe = "."//trim(adjustl(cpe))
    endif

    f1 = trim(file1)//trim(cpe)//char(0)
    f2 = trim(file2)//trim(cpe)//char(0)

    ptr(1) = c_loc(f1)
    ptr(2) = c_loc(f2)

    diff_modtime = modtimediff(ptr)
    return 
end function diff_modtime

integer function rm_file(filename)
    character(len=*), intent(in) :: filename
    type(C_PTR) :: ptr(1)
    character(len=len(filename)+10) :: f1

    rm_file = 1

    if (mpp_pe()/=mpp_root_pe()) return

    f1 = trim(filename)//char(0)
    
    ptr(1) = c_loc(f1)

    rm_file = rmfile(ptr)

    return
end function rm_file

end program main

