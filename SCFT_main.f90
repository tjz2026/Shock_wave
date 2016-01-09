!-------------------------------------------------------------------------
! project : AB-diblock-polymer-film-SCFT
! program : main
! source  : SCFT_main.f90
! type    : main program
! author  : Jiuzhou Tang (tangjiuzhou@iccas.ac.cn)
! history : 22/09/2013 by Jz Tang
! purpose : the main program of AB-diblock polymer SCFT
! input   : none
! output  : none
! status  : unstable
! comment :
!-------------------------------------------------------------------------
program main
use nrtype,only : DP
use global_para
use constants
use control
use mpi
use mmpi
use wall
implicit none
REAL*8 :: minFE,minx
REAL(DP),external :: cal_scft
REAL(DP),external :: brent


call initialize()


!minFE=brent(x1,x2,x3,cal_scft,tol,minx)

!minFE=cal_scft(minx)

minFE=cal_scft(x2)

call global_arrays_clean()

if(myid == master) then
call SCMFT_print_footer()
endif

call mp_barrier()
call mp_finalize()
stop
end

















