!-------------------------------------------------------------------------
! project : AB-diblock-polymer-film-SCFT
! program : global_para,mpi
! source  : global_parameters.f90
! type    : module
! author  : Jiuzhou Tang (tangjiuzhou@iccas.ac.cn)
! history : 19/09/2013 by Jz Tang
! purpose : define and declare all the global parameters and the global variables 
!           in the global_para module
! input   : none
! output  : none
! status  : unstable
! comment :
!-------------------------------------------------------------------------


module global_para
use nrtype,only:SP,DP
implicit none
!***********parameters for AB-diblock  chain*****
REAL(DP) :: NxAB   !flory-hugins para
REAL(DP) :: XNw_A,XNw_B   !flory-hugins para
REAL(DP) :: fa,fb 
INTEGER  :: Nmax,NA,NB !discretize grid number for polymer contour length L
REAL(DP) :: ds
REAL(DP) :: C !rho_0 /N
INTEGER :: N_kuhn
INTEGER :: N_lam
 

!*************parameters for spatial grid***************************
INTEGER,parameter ::  Ndim=3 !system dimention,1 for 1D,2 for 2D,3 for 3D
INTEGER,parameter :: Cell_dim=3
INTEGER ::  SIDEx,SIDEy,SIDEz !spatial grid number,make sure they are equal to 2^N for fftw 
INTEGER :: M_grid !M_grid=SIDEx*SIDEy*SIDEz for 3D case
INTEGER :: LOCAL_SIZE !locao_size=SIDEx*SIDEy*local_nz for 3D case
REAL(DP) :: Lx,Ly,Lz ! unit cell length for x,y,z
REAL(DP) :: dx,dy,dz
REAL(DP) :: x1,x2,x3 !box size range for brent
REAL(DP) :: M_v    !total volume
REAL(DP) :: density_ava !avaraged density 
REAL(DP) :: z_d  !wall depth

!!!para for free Energy
REAL(DP):: totden,pff_global,tempEE_global,FE_global,fe_wall

!*************parameters for SCFT *********************************
INTEGER,parameter :: MAXITS=50000 !maxmum iteration step number
INTEGER :: iter_method !iteration scheme,1 for simple mixing,2 for anderson_mixing
REAL(DP) :: lambda_simple 
REAL(DP),PARAMETER :: TOL=1.0e-2
REAL(DP),PARAMETER ::ksi=800.0
!!converge criterion
REAL(DP),PARAMETER :: CC=1.0e-5
REAL(DP),PARAMETER :: delta_t1=1.0
REAL(DP),PARAMETER :: delta_t2=0.12
INTEGER :: n_iter,n_iter_M
INTEGER :: wall_index
!**************variable for SCFT**************************
!*****for tensor like M,S
integer ,parameter :: N_dim_ddm=3
type node
DOUBLE PRECISION :: xx
DOUBLE PRECISION :: xy
DOUBLE PRECISION :: xz
!DOUBLE PRECISION :: yx
!DOUBLE PRECISION :: yy
DOUBLE PRECISION :: yz
!DOUBLE PRECISION :: zx
!DOUBLE PRECISION :: zy
DOUBLE PRECISION :: zz
end type

type node2
       INTEGER :: x
       INTEGER :: y
end type



REAL(DP),allocatable :: wA(:),wB(:),i_W_plus(:),W_minus(:)
REAL(DP),allocatable :: RA(:),RB(:),RWALL(:)
REAL(DP),allocatable :: R_half(:),R_end(:)
REAL(DP),allocatable :: RHOA(:,:),RHOB(:,:)




!*********for Anderson mixing 
integer,parameter :: n_r_WAB=10
INTEGER,parameter ::  Anderson_nim_AB=10
integer,parameter :: Num_step_simple_mixing_for_WAB=20
INTEGER,parameter ::  Anderson_nim_M=4
integer,parameter :: n_r_M=4
INTEGER,parameter ::  N_WABITER_PER_M=10
INTEGER,parameter ::  Num_simple_mixing_AB=10
INTEGER,parameter ::  Num_simple_mixing_M=10




REAL(DP),parameter :: lambda_WAB_anderson=0.10d0
REAL(DP),parameter :: lambda_M_anderson=0.10d0
!!!
REAL(DP),parameter :: lanbtwA=0.10d0
REAL(DP),parameter :: lanbtwB=0.10d0
REAL(DP),parameter :: lanbtM=0.080d0


REAL(DP) :: t1,t2,ta_diff,tb_diff
REAL(DP) :: lambda_WAB_anderson_const
REAL(DP),allocatable :: wA_out(:,:),wB_out(:,:)
REAL(DP),allocatable :: dw_A(:,:),dw_B(:,:)
REAL(DP),allocatable :: dA_anderson(:,:),dB_anderson(:,:)

!*****************************

!# if defined (Dim2)
!REAL(DP),allocatable :: phiA(:,:),phiB(:,:)
!REAL(DP),allocatable :: wA(:,:),wB(:,:)
!REAL(DP),allocatable :: yita(:,:)
!	type(node),allocatable :: Mr(:,:)
!	type(node),allocatable :: Sr(:,:)
!
!#endif  /* Dim2 */

# if defined (Dim3)
REAL(DP),allocatable :: phiA(:,:,:),phiB(:,:,:)
!REAL(DP),allocatable :: wA(:,:,:),wB(:,:,:)
!REAL(DP),allocatable :: yita(:,:,:)
#endif  /* Dim3 */

REAL(DP),allocatable :: qA(:,:) !balock 0~NA
REAL(DP),allocatable :: qB(:,:)
REAL(DP),allocatable :: qAstar(:,:)
REAL(DP),allocatable :: qBstar(:,:)
!!for pff calculation
REAL(DP),allocatable :: qBend(:)



type node3
      integer :: i
      integer :: j
      real(DP) :: value
end type
type node4
      integer :: i
      integer :: j
      integer :: k
      real(DP) :: value
end type

type node5
     integer :: i
     integer :: j
end type

type node6
     integer :: i
     integer :: j
     integer :: k
end type
type node7
     integer :: x
     integer :: y
     integer :: z
end type
type(node7),allocatable :: kD(:)

!type(node5),allocatable :: counter(:)
!different types to init the w field 
integer :: wAB_init_type  
 !6 read from input file named wfield.txt
 !0 for random
 !1 for FCC initial 3D
 !2 for  BCC initial 3D
 !3 for Gyroid initial 3D
 !4 for Cylinder initial 3D 
 !5 for Lamellar initial 3D 
 !6 for P4 initial 3D 



 REAL(DP),allocatable :: exp_K2(:)
end module global_para
   module mpi
  include 'mpif.h'
  end module mpi
