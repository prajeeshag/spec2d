module phys_mod

use mpp_mod, only : mpp_init, FATAL, WARNING, NOTE, mpp_error
use mpp_mod, only : mpp_npes, mpp_get_current_pelist, mpp_pe
use mpp_mod, only : mpp_exit, mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : mpp_sync, mpp_root_pe, mpp_broadcast, mpp_gather
use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist

use mpp_domains_mod, only : mpp_define_domains, domain2d, mpp_get_compute_domain

use fms_mod, only : read_data, write_data, open_namelist_file, close_file, fms_init
use fms_mod, only : file_exist

use time_manager_mod, only : time_type

use radiation_mod, only : init_radiation

use sfc_mod, only : init_sfc

implicit none
private

public :: init_phys

type(domain2D), pointer :: domain
real :: deltim
integer :: ntrac, nlev

logical :: initialized=.false.

contains

!--------------------------------------------------------------------------------   
subroutine init_phys(Time,deltim_in,domain_in,ntrac_in,nlev_in)
!--------------------------------------------------------------------------------   
    type(time_type), intent(in) :: Time
    type(domain2d), target :: domain_in
    real, intent(in) :: deltim_in 
    integer, intent(in) :: ntrac_in, nlev_in

    deltim = deltim_in

    domain => domain_in

    ntrac = ntrac_in
    
    nlev = nlev_in


    call init_radiation(Time,deltim,domain,ntrac,nlev) 

    call init_sfc(Time,deltim,domain)

    initialized = .true.

end subroutine init_phys


!!--------------------------------------------------------------------------------   
!subroutine phys(Time,tlyr,tr,p)
!!--------------------------------------------------------------------------------   
!
!end subroutine phys

end module phys_mod
