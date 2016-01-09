! project : wormlike chain SCFT
! program : SCFT
! source  : wall_density.f90
! type    : module
! author  : Jiuzhou Tang (tangjiuzhou@iccas.ac.cn)
!purpose :: This module is for building the wall density profile according to 
!           Glenn's so called MASK skill.

!!!!  Warning, this module currently only works for 2D model
module wall
implicit none
     !integer,public :: periodic ! periodic sine like wall density for the bottom wall
     double precision,public :: periodic ! zhoupan change 
     double precision,public :: wall_thickness
     double precision,public :: zoverdelta
     double precision,public :: wave_intensity
     integer,public :: wall_at
     real*8,public,allocatable :: R_wall_top(:)
     real*8,public,allocatable :: R_wall_bottom(:)
     real*8,public,allocatable :: R_wall_sum(:)

  contains
         

       subroutine Wall_builder()
       USE nrtype,only :Pi
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
       USE mpi_fftw3_operation
       implicit none
       integer :: n,i,j,k,errorw
       real*8  :: arg,r,delta_pattern,temp_thickness
            
           
           !build the upper wall
           if (allocated(R_wall_top)) then
           else
           allocate(R_wall_top(0:local_size-1))
           endif

           if (allocated(R_wall_bottom)) then
           else
           allocate(R_wall_bottom(0:local_size-1))
           endif

           if (allocated(R_wall_sum)) then
           else
           allocate(R_wall_sum(0:local_size-1))
           endif

              write(*,*) "wall_thickness=",wall_thickness
              write(*,*) "zoverdelta",zoverdelta

               wall_at=1
              do i=0,LOCAL_SIZE-1
                r=(kD(i)%x-1)*dy
                arg=-1.0*wall_thickness
                arg=arg+r*(1.0-2.0*wall_at) 
                !arg=arg+wall_at*SIDEx*dx !yiqiande wall
                arg=arg+wall_at*SIDEy*dy!zhoupan changes the wall
                arg=arg*zoverdelta
                R_wall_top(i)=0.5*(1.0-tanh(arg))
               enddo
              
              wall_at=0
           !build the bottom sine like wall
              do i=0,LOCAL_SIZE-1
                r=(kD(i)%x-1)*dy 
                temp_thickness=wall_thickness
                temp_thickness=temp_thickness*(wave_intensity*cos(2*Pi*periodic* &
                ((kD(i)%y-1+local_x_start)/(SIDEx*1.0d0)))+1.0)
                !arg=-1.0*wall_thickness
                arg=-1.0*temp_thickness
                arg=arg+r*(1.0-2.0*wall_at) 
                arg=arg+wall_at*SIDEy*dy
                arg=arg*zoverdelta
                R_wall_bottom(i)=0.5*(1.0-tanh(arg))
                !R_wall_bottom(i)=R_wall_bottom(i)*(wave_intensity*sin(2*Pi*periodic* & 
                !                 ((kD(i)%y-1+local_x_start)/(SIDEx*1.0d0)))+1.0)
                !R_wall_top(i)=0   
                !R_wall_bottom(i)=0

                !R_wall_bottom(i)= R_wall_top(i) !zhoupan change to be plate
                R_wall_sum(i)=R_wall_top(i)+R_wall_bottom(i)
               enddo

           
         call wall_dump()

    end subroutine Wall_builder
    
    subroutine wall_dump()
       USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
       USE mpi_fftw3_operation
       implicit none
       integer :: i,j,k,aaa,bbb,K_i,k_j,K_k
        character(len=30)::aa
        character(len=30)::bb
        aaa=myid
       write(aa,*) aaa

         open(unit=23,file= 'Wall_top'//trim(adjustl(aa)) // '.dat',status='replace')
         k=0
        do K_i=0,local_nx-1
        do k_j=0,SIDEy-1
          write(23,'(2f10.6,f12.8)') (K_i+local_x_start)*dx,K_j*dy,R_wall_top(k)
          k=k+1
         enddo
         enddo
        close(23)

         open(unit=23,file= 'Wall_bottom'//trim(adjustl(aa)) // '.dat',status='replace')
         k=0
        do K_i=0,local_nx-1
        do k_j=0,SIDEy-1
          write(23,'(2f10.6,f12.8)') (K_i+local_x_start)*dx,K_j*dy,R_wall_bottom(k)
          k=k+1
         enddo
         enddo
        close(23)


   
         open(unit=23,file= 'Wall_sum'//trim(adjustl(aa)) // '.dat',status='replace')
         k=0
        do K_i=0,local_nx-1
        do k_j=0,SIDEy-1
          write(23,'(2f10.6,f12.8)') (K_i+local_x_start)*dx,K_j*dy,R_wall_sum(k)
          k=k+1
         enddo
          write(23,'(2f10.6,f12.8)') 
         enddo
        close(23)

     call mp_barrier()
end subroutine wall_dump
          
end module wall










     
