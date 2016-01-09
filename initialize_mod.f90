 subroutine initialize()
 USE nrtype,only :DP
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 !USE fftw !user defined module
 use fftw3_mpi
 use mpi_fftw3_operation
 use wall
 implicit none

 logical ::  exists
 integer :: istat1,i,j,k
  
 Integer*8 planff     		! forward plan for use of fft => xyz
 !logical ::  exists ! why this statement causes error ,find out later!!!

! initialize mpi envirnoment
# if defined (MPI)

! initialize the mpi execution environment
     call mp_init()

! determines the rank of the calling process in the communicator
     call mp_comm_rank(myid)

! determines the size of the group associated with a communicator
     call mp_comm_size(nprocs)

# endif  /* MPI */


# if defined (OPENMP)
call OMP_SET_NUM_THREADS(16)
# endif  /* OPENMP */


if(myid == master) then
call SCMFT_print_header()
endif

  if(myid==master) then
         exists = .false.

         inquire (file = 'input.txt', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then
      open(unit=mytmp,file='input.txt')
read(mytmp,*) NxAB             
read(mytmp,*) fA             ! 
read(mytmp,*) Nmax              ! number of bonds
read(mytmp,*) wAB_init_type              ! wAB initialized type
read(mytmp,*) N_lam              ! N_lam periodic of lamelar
read(mytmp,*) SIDEx    
read(mytmp,*) SIDEy    
read(mytmp,*) SIDEz    
read(mytmp,*) Lx    
read(mytmp,*) Ly    
read(mytmp,*) Lz   
read(mytmp,*) z_d   
read(mytmp,*) x1    
read(mytmp,*) x2    
read(mytmp,*) x3    
close(mytmp)
close(mytmp)
      else 
      write(mystd,*) "input file missing ,program stop!"
      stop
      endif
endif  !endif myid==master


  if(myid==master) then
         exists = .false.

         inquire (file = 'wall_para.txt', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then
      open(unit=mytmp,file='wall_para.txt')
       read(mytmp,*) wall_thickness
       read(mytmp,*) zoverdelta
       read(mytmp,*) wave_intensity
       read(mytmp,*) periodic
       read(mytmp,*) XNw_A
       read(mytmp,*) XNw_B
       close(mytmp)
       else 
       write(mystd,*) "wall input file missing ,program stop!"
       stop
       endif
        endif  



! since input parameters may be updated in master node, it is important
! to broadcast input parameters from root to all children processes

!-----------------------------------------------------------------------
!using the so called function overload scheme in the sofware engineering,
!we have the uniform subroutine mp_bcast for every type of data we want to
! bcast.
# if defined (MPI)
      call mp_barrier()
      call mp_bcast( NxAB, master )                                    
      call mp_bcast( fA, master )                                      
      call mp_bcast( Nmax, master )                                   
      call mp_bcast( wAB_init_type, master )                                
      call mp_bcast( N_lam, master )                                
      call mp_bcast( SIDEx, master )                                  
      call mp_bcast( SIDEy, master )         
      call mp_bcast( SIDEz, master )         
      call mp_bcast( lx , master )                                     
      call mp_bcast( ly , master )                                    
      call mp_bcast( lz , master )                                  
      call mp_bcast( z_d , master )                                  
      call mp_bcast( x1 , master )                                  
      call mp_bcast( x2 , master )                                  
      call mp_bcast( x3 , master )                                  
      call mp_barrier() 

      call mp_bcast( wall_thickness , master )                                  
      call mp_bcast( wave_intensity , master )                                  
      call mp_bcast( zoverdelta , master )                                  
      call mp_bcast( periodic , master )                                  
      call mp_bcast( XNw_A , master )                                  
      call mp_bcast( XNw_B , master )                                  
      call mp_barrier() 
# endif  /* MPI */
!write(*,*) "SIDEx,SIDEy,SIDEz=",SIDEx,SIDEy,SIDEz,"on",myid
!write(*,*) "labntWA,B,M",lanbtwA,lanbtwB,lanbtM,"on",myid




!set up some parameters: check out !
 fb=1.0-fa
 ds=1.0/Nmax
 NA=NINT(fA*Nmax)
 NB=Nmax-NA
 C=1.0
!brent parameters,may not needed
 x1=x1/SIDEx
 x2=x2/SIDEx
 x3=x3/SIDEx



!checking paras
write(*,*) "C=",C,"on",myid
write(*,*) "delta_t",delta_t1,"on",myid
write(*,*) "z_d",z_d,"on",myid
write(*,*) "fA",fA,"on",myid
if(myid==0) then
 if(mod(NA,2)==0 .and. mod(NB,2)==0) then
  write(*,*) "NA,NB checked!"
  else
  write(*,*) "NA,NB=",NA,NB
  stop
  endif
endif

# if defined (Dim3)
       M_grid=SIDEx*SIDEy*SIDEz
      call fftw3_mpi_3d_dft_init()
      local_size=local_nx*SiDEy*SIDEz
# else /* Dim3 */
       M_grid=SIDEx*SIDEy
      call fftw3_mpi_3d_dft_init()
      local_size=local_nx*SiDEy

# endif /* Dim3 */

      call mp_barrier()
!
!!!!***************************begin to intialize w and scft variables**

call init_arrays()
call init_k_index()

if(wAB_init_type<0) then
call readW()
else
call init_w()
endif


if(myid==0) then
write(*,*) "done init_w and wall density"
endif
call mp_barrier()


!!!!test region ,remove all these after the coding is done

call mp_barrier()
if(myid==0) then
write(*,*) "done test-simpson!"
endif



!!!test matrix inverse
!if(myid==0) then
!call test_matrix_inverse()
!endif




!!! init the dump file 
!***************************************************************
if(myid==0) then
open(unit=48,file='FreeEG.dat',status='new')
close(48)
open(unit=43,file='error.dat',status='new')
close(43)
endif

!***********************************************************
!checkout the input paras
if(myid==0) then
write(*,*) "wAB_init_type=",wAB_init_type,"on",myid
write(*,*) "SIDEy,SIZEz",SIDEy,SIDEz,"on",myid
write(*,*) "local_nx=",local_nx,"on",myid
write(*,*) "local_size=",local_size,"on",myid
write(*,*) "x1,x2,x3",x1,x2,x3,"on",myid
endif
 
call mp_barrier()
end subroutine initialize

!!!!test region ,remove after it is done

subroutine init_arrays()
 USE nrtype,only :DP
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 !USE fftw !user defined module
 implicit none
 integer :: istat2,i,j,k
 REAL(DP) :: temp,temp1,temp2,sumTHETA
 integer :: LS,MB,MA
 integer :: n_anderson
 integer :: nonzero_counter

 LS=LOCAL_SIZE-1

 !*****************note that N_dim_ddm==3!
  

 allocate(qA(0:NA,0:LS),stat=istat2)
 allocate(qB(0:NB,0:LS),stat=istat2)


 allocate(qAstar(0:NA,0:LS),stat=istat2)
 allocate(qBstar(0:NB,0:LS),stat=istat2)
 allocate(wA(0:LS),stat=istat2)
 allocate(wB(0:LS),stat=istat2)
 allocate(i_W_plus(0:LS),stat=istat2)
 allocate(W_minus(0:LS),stat=istat2)


 allocate(RA(0:LS),stat=istat2)
 allocate(RB(0:LS),stat=istat2)

 allocate(exp_K2(0:LS),stat=istat2)

if(istat2/=0) then
write(*,*) "allocate failed,exit!","on",myid
stop
else
write(*,*) "allocate arrays succeeded","on",myid

endif 


end subroutine init_arrays





          


!!!!!!!!!!!!!!!!!!!initialize the input field!!!!!!!!!!!!!!!!!!

       subroutine init_w()
          USE nrtype,only :PI
         USE global_para
         USE mpi
         USE control
         USE constants
         USE utility
         USE mmpi
         USE mpi_fftw3_operation
         !USE fftw !user defined module
         implicit none
         integer :: k_i,k_j,k_k
         integer :: KK,j,k,i
         integer :: initseed
         integer :: dt(8)
         
   if(wAB_init_type==0) then
          !random init
          call date_and_time(values=dt)
          initseed=(dt(8)-500)*12345*(myid+1)+43210
          
          KK=0  
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
# if defined (Dim3)
                 do k_i=0,SIDEz-1
# endif /* Dim3 */
                   WA(KK)=1.0*ran2(initseed)
                   WB(KK)=1.0*ran2(initseed)
                   KK=KK+1
                  enddo
             enddo
# if defined (Dim3)
          enddo
# endif /* Dim3 */

   else if(wAB_init_type==4) then
      !!!Cylinder initial 3D
     KK=0
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
# if defined (Dim3)
                 do k_i=0,SIDEz-1
# endif /* Dim3 */
                   WA(KK)=NXab*(1-fA*(1+0.8*cos(2.0*PI*(local_x_start+K_k+1)/SIDEx)* &
                         cos(2.0*PI*(K_j+1)/SIDEy)))

                   WB(KK)=NXab*fA*(1+0.8*cos(2.0*PI*(local_x_start+K_k+1)/SIDEx)* &
                         cos(2.0*PI*(K_j+1)/SIDEy))
                   KK=KK+1
                  enddo
             enddo
# if defined (Dim3)
          enddo
# endif /* Dim3 */

else if(wAB_init_type==5) then
      !!!Lamellar initial 3D
     KK=0
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
# if defined (Dim3)
                 do k_i=0,SIDEz-1
# endif /* Dim3 */
                   WA(KK)=NXab*(1-fA*(1+0.8*cos(2.0*N_lam*PI*(K_j)/SIDEy)))
                   WB(KK)=NXab*fA*(1+0.8*cos(2.0*N_lam*PI*(K_j/SIDEy)))
                   KK=KK+1
             enddo
          enddo
# if defined (Dim3)
          enddo
# endif /* Dim3 */


   else if(wAB_init_type==6) then
     KK=0
          do K_k=0,Local_nx-1
             do k_j=0,SIDEy-1
# if defined (Dim3)
                 do k_i=0,SIDEz-1
# endif /* Dim3 */
                   WA(KK)=NXab*(1-fA*(1+0.8*cos(2.0*N_lam*PI*(K_k+local_x_start)/SIDEx)))
                   WB(KK)=NXab*fA*(1+0.8*cos(2.0*PI*N_lam*(K_k+local_x_start)/SIDEx))
!                   WA(KK)=NXab*(1.0d0+cos(2.0*N_lam*PI*(K_k+local_x_start)/SIDEx))*0.5d0
!                   WB(KK)=NXab-WA(KK)
                   KK=KK+1
             enddo
          enddo
# if defined (Dim3)
          enddo
# endif /* Dim3 */




endif

    i_W_plus(:)=(WA(:)+WB(:))*0.5
    W_minus(:)=(WB(:)-WA(:))*0.5
   end subroutine init_w




  subroutine readW()
   USE nrtype,only :DP
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
  implicit none
  integer :: i,j,k,aaa
  integer :: ii,iitot
  REAL*8 :: x,y,z,sumaW,sumaWtot
  REAL*8 :: sumM00,sumM01,sumM02,sumM10,sumM11,sumM12,sumM20,sumM21,sumM22,sumMtot
  character(len=30)::aa

  aaa=myid
  ii=myid
  write(aa,*) aaa
  open(unit=23,file= 'W' //trim(adjustl(aa)) // '.dat',status='old')
  do K=0,LOCAL_SIZE-1
  read(23,*) i,WA(K),WB(K)
  
  enddo
  close(23)

    i_W_plus(:)=(WA(:)+WB(:))*0.5
    W_minus(:)=(WB(:)-WA(:))*0.5

 end subroutine readW





 subroutine global_arrays_clean()
 USE nrtype,only :DP
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 implicit none


 deallocate(qA)
 deallocate(qB)
 deallocate(qAstar)
 deallocate(qBstar)
 deallocate(wA)
 deallocate(wB)
 deallocate(i_W_plus)
 deallocate(W_minus)
 deallocate(RA)
 deallocate(RB)
 deallocate(exp_K2)



end subroutine global_arrays_clean

 



           







subroutine  exponents_Karray()
!!be careful,this subroutine is going to be called ever time 
!the box size is changed
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 implicit none
 integer :: i,j,k
 integer :: K_i,K_ii
 integer :: K_j,K_jj
 integer :: K_k,K_kk
 REAL(DP) ::k_x_factor,k_y_factor,k_z_factor


k_x_factor=2.0*PI/(SIDEx*dx)
k_y_factor=2.0*PI/(SIDEy*dy)
k_z_factor=2.0*PI/(SIDEz*dz)


!k counter
k=0
!!the first dimension must be z 

do K_i=0,local_nx-1 
!do K_i=0,0 
   
  if((k_i+local_x_start)<=(SiDEx/2)) then
   k_ii=k_i+local_x_start  
  else 
   k_ii=(k_i+local_x_start)-SIDEx 
  endif

    do k_j=0,SIDEy-1
    ! do k_j=1,1
        if(K_j<=(SIDEy/2)) then 
          k_jj=K_j
        else
          K_jj=K_j-SIDEy
        endif
# if defined (Dim3)
           do k_k=0,SIDEz-1
           !do k_k=0,0
             if(K_k<=(SIDEz/2)) then
             k_kk=K_k
             else
             k_kk=K_k-SIDEz
             endif
# endif /* Dim3 */

                    exp_K2(k)=exp(-1.0*ds*((k_ii*k_x_factor)**2 + &
                                     (k_jj*k_y_factor)**2 + &
                                      (k_kk*k_z_factor)**2))
                  k=k+1 

                 enddo
            enddo
# if defined (Dim3)
         enddo
# endif /* Dim3 */

         end subroutine  exponents_Karray


subroutine init_k_index()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 USE mpi_fftw3_operation
 implicit none
 integer :: istat,k,k_i,k_j,k_k
 allocate(kD(0:LOCAL_SIZE-1),stat=istat)
 if(istat/=0) then
  if(myid==0) then
  write(*,*) "error alocate k_index"
  call mp_barrier()
  stop
  endif
 endif
!!note that local_data is formed as (1:SIDEz,1:SIDEy,1:local_nx) not
!as (0:SIDEz-1,0:SIDEy-1,0:local_nx-1)
 k=0
 do k_i=0,local_nx-1
    do k_j=0,SIDEy-1
# if defined (Dim3)
      do k_k=0,SIDEz-1
# endif /* Dim3 */

# if defined (Dim3)
        kD(k)%x=k_k+1
        kD(k)%y=k_j+1
        kD(k)%z=k_i+1
# else /* Dim3 */
        kD(k)%x=k_j+1
        kD(k)%y=k_i+1
# endif /* Dim3 */
        k=k+1
      enddo
    enddo
# if defined (Dim3)
  enddo
# endif /* Dim3 */

end subroutine init_k_index











