! project : wormlike chain SCFT
! program : SCFT
! source  : wall_confine_mod.f90
! type    : module
! author  : Jiuzhou Tang (tangjiuzhou@iccas.ac.cn)
!purpose :: This module is for building the wall density profile according to 
!           Glenn's so called MASK skill.Meanwhile,this module can be used generally.

module wall
implicit none
    type wall_confine
     integer,public :: confine_type ! what kind of confinement,1 for plane,2 for cylinder,3 for spher,4 for chemical pattern
     integer,public :: wall_dim  ! dimension of wall
     integer,public :: wall_number !how many walls
     real*8,public :: chemostripe ! the period of chemical pattern on X axis 
     real*8,public :: chemostriperatio ! the width ratio of filled and empty wall
     integer,public :: period_x,period_y ! hexagonal pattern on x and y axis
     real*8,public :: R_cylinder          ! the radius of cylinder,used both for confinement and pattern
     real*8,public :: R_spher         ! the radius of spehrical confinement wall
     double precision,public :: wall_thickness
     double precision,public :: zoverdelta
end type

     integer,public :: wall_axis
     integer,public :: wall_at
     real*8,public,allocatable :: R_wall(:,:)
     real*8,public,allocatable :: R_wall_ava(:)
     real*8,public,allocatable :: R_wall_sum(:)
     real*8,public :: Xw_A(6),Xw_B(6)

!    integer,public  :: confine_type  ! types of confinement,plane,cylinder,spher and other user defined types.
!    double precision,public ::  wall_thickness !  control parameter for wall thickness
!    double precision,public ::  zoverdeltd    ! control parameter for transition width
!    integer,public :: wall_axis   ! the axis that the wall build on, X,Y,or Z. especially for plane confinement
!    integer,public :: wall_at     ! two end of the wall, with value of either 0 or 1.
!    integer,public :: wall_number  ! number of different walls,currently only single wall is surpported. 
     type(wall_confine) :: wall_list(16)

  contains
         

       subroutine Wall_builder()
       USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
       USE mpi_fftw3_operation
       implicit none
       integer :: n,i,j,k,l,errorw
       real*8  :: arg,r,delta_pattern
       real*8,allocatable :: zero(:),R_temp(:) 
       ! for cylinder pattern only
       integer :: upper_x,upper_y,filled,box_id_x,box_id_y
       real*8 ::  center_x,center_y,center_z,hex_x,hex_y
       real*8 :: wall_thickness,zoverdelta,R_cylinder,R_spher 
       integer :: chemical_period
       integer :: period_x,period_y
       real*8 :: pattern_ratio
            
           
          call wall_info_init() 
          if (.not. allocated(R_wall_sum)) then
          allocate(R_wall_sum(0:local_size-1),stat=errorw)
          endif
          if(.not. allocated(zero)) then
          allocate(zero(0:local_size-1),stat=errorw)
          endif
          if(.not. allocated(R_temp)) then
          allocate(R_temp(0:local_size-1),stat=errorw)
          endif

          zero=0.0

          n=wall_list(wall_index)%wall_number
          wall_thickness=wall_list(wall_index)%wall_thickness
          zoverdelta=wall_list(wall_index)%zoverdelta
          R_cylinder=wall_list(wall_index)%R_cylinder
          R_spher=wall_list(wall_index)%R_spher
          period_x=wall_list(wall_index)%period_x
          period_y=wall_list(wall_index)%period_y
          
          if(.not. allocated (R_wall)) then
          allocate(R_wall(0:local_size-1,1:n),stat=errorw)
          allocate(R_wall_ava(0:local_size-1),stat=errorw)
          endif

      if(wall_list(wall_index)%confine_type==0) then
           R_wall(:,:)=0.0
           R_wall_sum(:)=0.0
      else if(wall_list(wall_index)%confine_type==1) then

           
          do j=1,n

           if(wall_index==1 .or. wall_index==4) then
             wall_axis=1
            else if(wall_index==2 .or. wall_index==5)then
             wall_axis=2
            else if(wall_index==3 .or. wall_index==6)then
             wall_axis=3
            else
             wall_axis=(j+mod(j,2))/2
            endif


             wall_at=1-mod(j,2)
             write(*,*) "wall axis is",wall_axis,"on",wall_at,"for",j,"wall"
             if(wall_axis==1) then
              do i=0,LOCAL_SIZE-1
                r=(kD(i)%z-1+local_x_start)*dx 
                arg=-1.0*wall_thickness
                arg=arg+r*(1.0-2.0*wall_at) 
                arg=arg+wall_at*SIDEx*dx
                arg=arg*zoverdelta
                R_wall(i,j)=0.5*(1.0-tanh(arg))
               enddo
              
             elseif(wall_axis==2) then
              do i=0,LOCAL_SIZE-1
                r=(kD(i)%y-1)*dy 
                arg=-1.0*wall_thickness
                arg=arg+r*(1.0-2.0*wall_at) 
                arg=arg+wall_at*SIDEy*dy
                arg=arg*zoverdelta
                R_wall(i,j)=0.5*(1.0-tanh(arg))
               enddo
             elseif(wall_axis==3) then
              do i=0,LOCAL_SIZE-1
                r=(kD(i)%x-1)*dz
                arg=-1.0*wall_thickness
                arg=arg+r*(1.0-2.0*wall_at) 
                arg=arg+wall_at*SIDEz*dz
                arg=arg*zoverdelta
                R_wall(i,j)=0.5*(1.0-tanh(arg))
               enddo
              endif
            enddo 

         !else if ....(! other confinement,cylinder,spher,chemical pattern)
         else if ( wall_list(wall_index)%confine_type==2) then
      ! 3D cylinder confinement,Z axis free.
      ! first ,locates the cylinder center,it is better that the XY box is square,ie.Lx=Ly,although not necessarily.
      center_x=SIDEx*dx/2.0
      center_y=SIDEy*dy/2.0
      do i=0,local_size-1
         r=(((kD(i)%z-1+local_x_start)*dx-center_x)**2 + ((kD(i)%y-1)*dy-center_y)**2)  
         r=sqrt(r)
         !note that here,R_cyinder is the radius of cylinder confinement wall
         if(r<R_cylinder) then 
         R_wall(i,1)=1.0
         else
         arg=(r-R_cylinder)-wall_thickness
         arg=arg*zoverdelta
         R_wall(i,1)=0.5*(1.0-tanh(arg))
         endif
         R_wall(i,1)=1.0-R_wall(i,1)
      enddo
      
         else if ( wall_list(wall_index)%confine_type==3) then
         !3D spehrical confinement,just like cylinderical confinement
         center_x=SIDEx*dx/2.0
         center_y=SIDEy*dy/2.0
         center_z=SIDEy*dy/2.0
      do i=0,local_size-1
         r=(((kD(i)%z-1+local_x_start)*dx-center_x)**2 + ((kD(i)%y-1)*dy-center_y)**2 + &
           ((kD(i)%x-1)*dz-center_z)**2)  
         r=sqrt(r)
         if(r<R_spher) then
         R_wall(i,1)=1.0
         else
         arg=(r-R_spher)-wall_thickness
         arg=arg*zoverdelta
         R_wall(i,1)=0.5*(1.0-tanh(arg))
         endif
         R_wall(i,1)=1.0-R_wall(i,1)
      enddo



         else if( wall_list(wall_index)%confine_type==4) then
              chemical_period=wall_list(wall_index)%chemostripe
              pattern_ratio=wall_list(wall_index)%chemostriperatio
             if(wall_index==13) then
             ! chmical pattern on X, confined on Y .
                 delta_pattern=(SIDEx*dx)/chemical_period
                R_wall(:,1)=0.0
                R_wall(:,2)=0.0
               do i=0,Local_size-1
                r=(kD(i)%z-1+local_x_start)*dx
             !warning,cruently only two patterns-wall is adopted,if wall_number>2,error may occur.  
                 
                 arg=((kD(i)%y-1)*dy)-wall_thickness
                 arg=arg*zoverdelta

                if((r/delta_pattern - int(r/delta_pattern))<pattern_ratio) then
                !this means we are inside the filled wall
                !Wall profile as we move away from Y=0
                 R_wall(i,1)=0.5*(1.0-tanh(arg))
                 else
                 R_wall(i,2)=0.5*(1.0-tanh(arg))
                 endif

                enddo

             else if(wall_index==14) then
             ! chmical pattern on X, wall on Z=0, 3D case .
                 delta_pattern=(SIDEx*dx)/(wall_list(wall_index)%chemostripe)
                  pattern_ratio=wall_list(wall_index)%chemostriperatio
                R_wall(:,1)=0.0
                R_wall(:,2)=0.0
               do i=0,Local_size-1
                r=(kD(i)%z-1+local_x_start)*dx
             !warning,cruently only two patterns-wall is adopted,if wall_number>2,error may occur.  
             
                 arg=((kD(i)%x-1)*dz)-wall_thickness
                 !arg=r-wall_thickness
                 arg=arg*zoverdelta
                if((r/delta_pattern - int(r/delta_pattern))<pattern_ratio) then
                 R_wall(i,1)=0.5*(1.0-tanh(arg))
                else
                 R_wall(i,2)=0.5*(1.0-tanh(arg))
                 endif
                enddo

             else if(wall_index==15) then
             ! chemical pattern on XY ,cylinder pattern,wall on z=0
             ! we will first generate a single unit of hexagonal cell,then copy this unit cell to fullfill the
             ! whole XY box.
             ! in order to determine the pattern,we need pattern period on both X and Y axis,namely,
             ! period_x,period_y,and the cylinder radius,R_cylinder.
             !in order to avoid the possible trouble, keep the period even ,like 2,4,8....
             upper_x=SIDEx/period_x -1
             upper_y=SIDEy/period_y -1
           ! OK, the first unit cell locates at (0,0),(0,upper_x),(upper_y,0),(upper_y,upper_x)
           ! there are five centers in a single unit cell of HEX
             R_wall(:,1)=0.0
             R_wall(:,2)=0.0
             ! first ,check out which unit cell are we in
                    
                 hex_x=SIDEx*dx/period_x
                 hex_y=SIDEy*dy/period_y

                 if(R_cylinder>0.5*Hex_y) then
                 write(*,*) "warning,R_cylinder too big,ready to exit!"
                 write(*,*) "R_cylinder",R_cylinder,"SIDEy*dy=",SIDEy*dy
                 stop
                 endif

                  
              do i=0,local_size-1
                 box_id_x=int((kD(i)%z-1+local_x_start)/(upper_x+1))
                 box_id_y=int((kD(i)%y-1)/(upper_y+1))

                 filled=0 
                 do l=1,5
                  if(l==1) then
                   center_x=0.0
                   center_y=0.0
                  else if(l==2) then
                   center_x=hex_x
                   center_y=0.0
                   else if(l==3) then
                   center_x=0.0
                   center_y=hex_y
                   else if(l==4) then
                   center_x=hex_x
                   center_y=hex_y
                    else 
                   center_x=hex_x/2
                   center_y=hex_y/2
                   endif
                   center_x=center_x+box_id_x*hex_x
                   center_y=center_y+box_id_y*hex_y
                   
                   
        r=sqrt((((kD(i)%z-1+local_x_start)*dx-center_x)**2+((kD(i)%y-1)*dy-center_y)**2))
                  if(r<R_cylinder) then
                  filled=1
                  exit
                  endif
               enddo  !enddo l=1,5

                 arg=((kD(i)%x-1)*dz)-wall_thickness
                 !arg=r-wall_thickness
                 arg=arg*zoverdelta
                  if(filled==1) then
                  R_wall(i,1)=0.5*(1.0-tanh(arg))
                  else 
                  R_wall(i,2)=0.5*(1.0-tanh(arg))
                  endif

              enddo !enddo i=0,local_size-1

              endif  !endif chemo pattern

         endif    ! endif all kinds of confinement   
                
                



        !correct the cross section of each walls,the earlier wall will eat out the later wall
        !if the wall density is larger than 1.0

        do i=1,wall_list(wall_index)%wall_number
            R_wall_sum(:)=0.0
            do j=1,i
             R_wall_sum(:)=R_wall_sum(:)+R_wall(:,j) 
            enddo
             R_wall_sum(:)=R_wall_sum(:)-1.0

             where (R_wall_sum<0.0) 
               R_wall_sum=zero
             end where

             !call zeroneg(R_wall_sum,Local_size)
             R_wall(:,i)=R_wall(:,i)-R_wall_sum(:)
             R_temp(:)=R_wall(:,i)

             where (R_temp<0.0) 
               R_temp=zero
             end where
             R_wall(:,i)=R_temp(:)

             !call zeroneg(R_wall(:,i),local_size)
         enddo
                      
            R_wall_sum(:)=0.0
        do i=1,wall_list(wall_index)%wall_number
            R_wall_sum(:)=R_wall_sum(:)+R_wall(:,i)
        enddo




         deallocate(R_temp)
         deallocate(zero)
           
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


       do i=1,wall_list(wall_index)%wall_number

        bbb=i
       write(bb,*) bbb
         open(unit=23,file= 'Wall' //trim(adjustl(bb)) // 'on' //trim(adjustl(aa)) // '.dat',status='replace')
         k=0
        do K_i=0,local_nx-1
        do k_j=0,SIDEy-1
         do K_k=0,SIDEz-1

          write(23,'(3f10.6,3f12.8)') (K_i+local_x_start)*dx,K_j*dy,K_k*dz,R_wall(k,i)
          k=k+1
         enddo
         enddo
        enddo
       close(23)
        enddo
     call mp_barrier()
   
         open(unit=23,file= 'Wall_sum'//trim(adjustl(aa)) // '.dat',status='replace')
         k=0
        do K_i=0,local_nx-1
        do k_j=0,SIDEy-1
         do K_k=0,SIDEz-1

          write(23,'(3f10.6,3f12.8)') (K_i+local_x_start)*dx,K_j*dy,K_k*dz,R_wall_sum(k)
          k=k+1
         enddo
         enddo
        enddo
        close(23)

     call mp_barrier()
end subroutine wall_dump
          





   subroutine wall_info_init()
       USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
       implicit none
  
       if(myid==0) then
         write(*,*) "wall confinement added,initialing wall profile ..."
       endif 
         
         SELECT CASE (wall_index) 
         case(1)
         write(*,*) "2D plane confinement on X axis "
         wall_list(1)%confine_type=1
         wall_list(1)%wall_dim=2
         wall_list(1)%wall_number=2
         case(2)
         write(*,*) "2D plane confinement on Y axis"
         wall_list(2)%confine_type=1
         wall_list(2)%wall_dim=2
         wall_list(2)%wall_number=2
         case(3)
         write(*,*) "2D plane confinement on XY axis"
         wall_list(3)%confine_type=1
         wall_list(3)%wall_dim=2
         wall_list(3)%wall_number=4
         case(4) 
         write(*,*) "3D plane confinement on X axis"
         wall_list(4)%confine_type=1
         wall_list(4)%wall_dim=3
         wall_list(4)%wall_number=2
         case(5) 
         write(*,*) "3D plane confinement on Y axis"
         wall_list(5)%confine_type=1
         wall_list(5)%wall_dim=3
         wall_list(5)%wall_number=2
         case(6) 
         write(*,*) "3D plane confinement on Z axis"
         wall_list(6)%confine_type=1
         wall_list(6)%wall_dim=3
         wall_list(6)%wall_number=2
         case(7) 
         write(*,*) "3D plane confinement on XY axis"
         wall_list(7)%confine_type=1
         wall_list(7)%wall_dim=3
         wall_list(7)%wall_number=4
         case(8) 
         write(*,*) "3D plane confinement on XZ axis"
         wall_list(8)%confine_type=1
         wall_list(8)%wall_dim=3
         wall_list(8)%wall_number=4
         case(9) 
         write(*,*) "3D plane confinement on YZ axis"
         wall_list(9)%confine_type=1
         wall_list(9)%wall_dim=3
         wall_list(9)%wall_number=4
         case(10) 
         write(*,*) "3D plane confinement on XYZ axis"
         wall_list(10)%confine_type=1
         wall_list(10)%wall_dim=3
         wall_list(10)%wall_number=6
         case(11) 
         write(*,*) "3D cylinder confinement "
         wall_list(11)%confine_type=2
         wall_list(11)%wall_dim=3
         wall_list(11)%wall_number=1
         case(12) 
         write(*,*) "3D spher confinement "
         wall_list(12)%confine_type=3
         wall_list(12)%wall_dim=3
         wall_list(12)%wall_number=1
         case(13) 
         write(*,*) "2D plane confinement,chemical pattern "
         !! to be done ,with two different walls on X,wall on Y=0
         wall_list(13)%confine_type=4
         wall_list(13)%wall_dim=2
         wall_list(13)%wall_number=2
         case(14) 
         write(*,*) "3D plane confinement,stripe chemical pattern on X "
         !! to be done ,with two different walls on X,wall on Z=0
         wall_list(14)%confine_type=4
         wall_list(14)%wall_dim=3
         wall_list(14)%wall_number=2
         case(15) 
         write(*,*) "3D plane confinement,cylinder chemical pattern "
         !! to be done ,with two different walls on XY,wall on Z=0
         wall_list(15)%confine_type=4
         wall_list(15)%wall_dim=3
         wall_list(15)%wall_number=2
         
         case(16) 
         write(*,*) "no wall "
         wall_list(16)%confine_type=0
         wall_list(16)%wall_dim=3
         wall_list(16)%wall_number=1

         case default 
        write(*,*) "oops,wall confinement type unknown,read to exit "
         stop
        end select
       
       
         if(wall_list(wall_index)%wall_dim/=Cell_dim) then
            if(myid==0) then
            write(*,*) "wall cell dim doesn't match Cell dim,exit"
            endif
            call mp_barrier()
            call mp_finalize()
            stop
         endif

      end subroutine wall_info_init
              
            
  end module wall   
          

         













     
