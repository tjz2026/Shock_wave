   
   module fftw3_mpi
   use, intrinsic :: iso_c_binding
   include 'fftw3-mpi.f03'
   end module fftw3_mpi

  module mpi_fftw3_operation
  use fftw3_mpi

  integer(C_INTPTR_T), public :: N1
  integer(C_INTPTR_T), public :: N2
  integer(C_INTPTR_T), public :: N3
  type(C_PTR),public :: fplan,bplan, cdata
  complex(C_DOUBLE_COMPLEX), pointer,public :: local_data(:)
  integer(C_INTPTR_T),public ::  alloc_local, local_N3, local_N3_offset
  integer :: local_nx,local_x_start
  contains

  subroutine fftw3_mpi_3d_dft_init()
  use mpi
  use fftw3_mpi
  use global_para,only :SIDEx,SIDEy,SIDEz
  USE control, ONLY: myid, nprocs
  implicit none
  N1=SIDEz   
  N2=SIDEy   
  N3=SIDEx   


    call fftw_mpi_init()
# if defined (Dim3)
!   get local data size and allocate (note dimension reversal)
  alloc_local = fftw_mpi_local_size_3d(N3, N2,N1, MPI_COMM_WORLD, &
                                       local_N3, local_N3_offset)
  cdata = fftw_alloc_complex(alloc_local)

  call c_f_pointer(cdata, local_data, [N1,N2,local_N3])

!   create MPI plan for in-place forward DFT (note dimension reversal)
  fplan = fftw_mpi_plan_dft_3d(N3,N2,N1, local_data, local_data, MPI_COMM_WORLD, &
                              FFTW_FORWARD, FFTW_EXHAUSTIVE)

  bplan = fftw_mpi_plan_dft_3d(N3,N2,N1, local_data, local_data, MPI_COMM_WORLD, &
                              FFTW_BACKWARD, FFTW_EXHAUSTIVE)

  write(*,*) "Local_N3=,local_N3_offset=",Local_N3,Local_N3_offset,"on",myid 

  local_nx=Local_N3
  Local_x_start=Local_N3_offset
# else /* Dim3 */
  alloc_local = fftw_mpi_local_size_2d(N3, N2, MPI_COMM_WORLD, &
                                       local_N3, local_N3_offset)
  cdata = fftw_alloc_complex(alloc_local)

  call c_f_pointer(cdata, local_data, [N2,local_N3])

!   create MPI plan for in-place forward DFT (note dimension reversal)
  fplan = fftw_mpi_plan_dft_2d(N3,N2, local_data, local_data, MPI_COMM_WORLD, &
                              FFTW_FORWARD, FFTW_EXHAUSTIVE)

  bplan = fftw_mpi_plan_dft_2d(N3,N2, local_data, local_data, MPI_COMM_WORLD, &
                              FFTW_BACKWARD, FFTW_EXHAUSTIVE)

  write(*,*) "Local_N3=,local_N3_offset=",Local_N3,Local_N3_offset,"on",myid 

  local_nx=Local_N3
  Local_x_start=Local_N3_offset
# endif /* Dim3 */


! compute transform (as many times as desired)
!  call fftw_mpi_execute_dft(fplan, data, data)
!  call fftw_mpi_execute_dft(bplan, data, data)
!  data=data/(32*32)
! write(*,*) "data(1,1)af fftw3",data(1,1),"on",myid

end subroutine fftw3_mpi_3d_dft_init

subroutine destroy_3d_mpi_fftw3()
  use fftw3_mpi
  implicit none
  call fftw_destroy_plan(fplan)
  call fftw_destroy_plan(bplan)
  call fftw_free(cdata)
end subroutine destroy_3d_mpi_fftw3

end module mpi_fftw3_operation




