       !note : assuming that only z direction is changed,more complex situation
       !will be considered later.
       subroutine Wfield_fill(SIDEz_old,Lz_old)
       USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
       USE mpi_fftw3_operation
       use wall
       implicit none
  integer :: i,j,k,aaa,K_i,k_j,K_k,k_prime,k_k_prime
  integer :: ii,iitot,SIDEz_old
  integer :: LOCAL_SIZE1
  real(DP),allocatable :: WA1(:),WB1(:)
  real(DP),allocatable :: zloc1(:),zloc(:)
  character(len=30)::aa
  integer :: initseed
  integer :: dt(8)
  real(DP) :: r,Lz_old
  type node10
     integer :: x
     integer :: y
     integer :: z
  end type
  type(node10),allocatable :: kD1(:)
         
          
 
  LOCAL_SIZE1=SIDEz_old*SIDEy*local_nx

  allocate(wA1(0:LOCAL_SIZE1-1)) 
  allocate(wB1(0:LOCAL_SIZE1-1)) 

  allocate(zloc1(0:LOCAL_SIZE1-1)) 
  allocate(kD1(0:LOCAL_SIZE1-1)) 

  allocate(zloc(0:LOCAL_SIZE-1)) 

  aaa=myid
  ii=myid
  write(aa,*) aaa
  ! better check out if the local_size1 matches the file,addtional code be added later.
  open(unit=23,file= 'W' //trim(adjustl(aa)) // '.dat',status='old')
  do K=0,LOCAL_SIZE1-1
  read(23,*) i,WA1(K),WB1(K)
  enddo
  close(23)
  
!!note that local_data is formed as (1:SIDEz,1:SIDEy,1:local_nx) not
!as (0:SIDEz-1,0:SIDEy-1,0:local_nx-1)
 k=0
 do k_i=0,local_nx-1
    do k_j=0,SIDEy-1
      do k_k=0,SIDEz_old-1
        kD1(k)%x=k_k+1
        kD1(k)%y=k_j+1
        kD1(k)%z=k_i+1
        zloc1(k)=((k_k)*1.0d0/(SIDEz_old))*Lz_old
        k=k+1
      enddo
    enddo
 enddo


!assuming that kD has alrready been calculated
 k=0
 do k_i=0,local_nx-1
    do k_j=0,SIDEy-1
      do k_k=0,SIDEz-1
        zloc(k)=k_k*dz
        k=k+1
      enddo
    enddo
 enddo

          !random init
          call date_and_time(values=dt)
          initseed=(dt(8)-500)*654321*(myid+1)+888888
          do K=1,1000
          r=ran2(initseed)
          enddo
          write(*,*) "check out initseed",initseed,ran2(initseed),"on",myid


do K=0,LOCAL_SIZE-1
 if(zloc(k)>Lz_old) then
 wA(K)=Nxab*(1-fA)*ran2(initseed)
 wB(K)=Nxab*(fA)*ran2(initseed)
 else if(zloc(k)<=(Lz_old/SIDEz_old)) then
     k_k_prime=0
     k_prime=0+(kD(k)%y-1)*SIDEz_old+(kD(k)%x-1)*SIDEz_old*SIDEy
     wA(K)=wA1(k_prime)
     wB(K)=wB1(k_prime)
 else    
      do ii=1,SIDEz_old-1
      if((ii+1)*(Lz_old/SIDEz_old)>=zloc(k)) then
      k_k_prime=ii
      exit
      endif
      enddo

     k_prime=k_k_prime+(kD(k)%y-1)*SIDEz_old+(kD(k)%x-1)*SIDEz_old*SIDEy
     wA(K)=wA1(k_prime)
     wB(K)=wB1(k_prime)
  endif 
 enddo 

  deallocate(wA1)
  deallocate(wB1)
  deallocate(zloc)
  deallocate(zloc1)
  deallocate(kD1)


end  subroutine Wfield_fill


