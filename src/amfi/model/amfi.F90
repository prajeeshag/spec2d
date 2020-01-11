program amfi

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather
use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist

use mpp_io_mod, only : mpp_open, MPP_RDONLY, mpp_close

use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain

use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init
use fms_mod, only : file_exist
use fms_io_mod, only : fms_io_exit

use time_manager_mod, only : time_type, set_calendar_type, operator(-), &
    operator(+), set_date, set_time, days_in_month, operator(/), operator(>), & 
    print_date, increment_date, date_to_string, operator(>=), get_date

use diag_manager_mod, only: diag_manager_init, diag_manager_end
use diag_manager_mod, only: get_base_date, DIAG_OTHER

use atmos_mod, only : init_atmos, update_atmos, end_atmos, save_restart_atmos

implicit none

type(time_type) :: Time_init, Time, Time_start, Time_end, &
        Run_length, time_step, time_restart
integer :: months=0, days=0, hours=0, minutes=0, seconds=0
integer :: num_atmos_calls
character(len=32) :: resfile='INPUT/atm.res'
character(len=512) :: err_msg
integer :: date_init(6), date(6), calendar_type
integer :: restart_interval(6) = 0
integer :: dt_atmos=0, unit, m, n
character(len=32) :: timestamp

integer :: clck_atmos, clck_atmosinit1step

namelist/amfi_nml/months, days, hours, minutes, seconds, dt_atmos, restart_interval

call mpp_init()
call fms_init()

unit = open_namelist_file()
read(unit,nml=amfi_nml)
call close_file(unit)

if( file_exist(resfile) )then
    call mpp_open( unit, resfile, action=MPP_RDONLY )
    read(unit,*)calendar_type
    read(unit,*)date_init
    read(unit,*)date
else
    call mpp_error(FATAL,trim(resfile)//' does not exist!')
endif

call set_calendar_type(calendar_type, err_msg)
if(trim(err_msg) /= '') then
  call mpp_error(FATAL, 'ERROR in init_atmos: '//trim(err_msg))
endif

call diag_manager_init(DIAG_OTHER)

call get_base_date ( date_init(1), date_init(2), date_init(3), &
     date_init(4), date_init(5), date_init(6)  )

if ( date_init(1) == 0 ) date_init = date

Time_init = set_date (date_init(1), date_init(2), date_init(3), &
                      date_init(4), date_init(5), date_init(6))

Time = set_date (date(1), date(2), date(3),  &
                 date(4), date(5), date(6))

if (Time_init>Time) call mpp_error(FATAL,'Init Time > Current Time')

Time_start = Time

Time_end = Time

do m=1,months
   Time_end = Time_end + set_time(0,days_in_month(Time_end))
end do

Time_end   = Time_end + set_time(hours*3600+minutes*60+seconds, days)




if (all(restart_interval<=0)) then
    time_restart = increment_date(Time_end,1,0,0,0,0,0)
else
    time_restart = increment_date(Time, restart_interval(1), restart_interval(2), &
    restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
endif 

if (mpp_pe()==mpp_root_pe()) then
  call mpp_open(unit, 'time_stamp.out', nohdrs=.True.)
  call get_date(Time_start, date(1), date(2), date(3), &
    date(4), date(5), date(6))
  write(unit,'(6i4,2x)') date
  call get_date(Time_end, date(1), date(2), date(3), &
    date(4), date(5), date(6))
  write(unit,'(6i4,2x)') date
  call mpp_close(unit)
endif



Run_length = Time_end - Time

time_step = set_time(dt_atmos,0)

num_atmos_calls = Run_length / time_step

clck_atmos = mpp_clock_id('Atmos')
clck_atmosinit1step = mpp_clock_id('Atmos_init1step')

call mpp_clock_begin(clck_atmos)
call mpp_clock_begin(clck_atmosinit1step)

call init_atmos(Time,real(dt_atmos))

do n = 1, num_atmos_calls
    Time = Time + time_step
    call update_atmos(Time)
    if (Time >= time_restart) then
        timestamp = date_to_string(Time)
        time_restart = increment_date(Time, restart_interval(1), restart_interval(2), &
        restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
        call mpp_error(NOTE,'Intermediate restart file is written and '&
            //trim(timestamp)//' is appended as prefix to each restart file name')
        call save_restart_atmos(timestamp)
        call atm_restart(Time,timestamp) 
    endif
    if (n==1) call mpp_clock_end(clck_atmosinit1step)
enddo

call mpp_clock_end(clck_atmos)

call diag_manager_end(Time)

call end_atmos(Time)

call atm_restart(Time)

call fms_io_exit()

call mpp_exit()


contains

subroutine atm_restart(Time_run, time_stamp)
    type(time_type),   intent(in)           :: Time_run
    character(len=*), intent(in),  optional :: time_stamp
    character(len=128)                      :: file_run
    integer :: yr, mon, day, hr, min, sec, date(6), unit

    if (present(time_stamp)) then   
        file_run = 'RESTART/'//trim(time_stamp)//'.atm.res'
    else
        file_run = 'RESTART/atm.res'
    endif

    call get_date(Time_run, date(1), date(2), date(3),  &
                   date(4), date(5), date(6))
    if (mpp_pe()==mpp_root_pe()) then
        call mpp_open(unit, file_run, nohdrs=.True.)

        write( unit, '(i6,8x,a)' )calendar_type, &
             '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'
        write( unit, '(6i6,8x,a)' )date_init, &
             'Model start time:   year, month, day, hour, minute, second'
        write( unit, '(6i6,8x,a)' )date, &
             'Current model time: year, month, day, hour, minute, second'
        call mpp_close(unit)
    end if

    return
end subroutine atm_restart

end program amfi
