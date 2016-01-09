         subroutine SCFT_dump(converge)
         USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
       USE mpi_fftw3_operation
       use wall
       !USE fftw !user defined module
       implicit none
       logical,intent(in) :: converge
       integer :: k,aaa,K_i,K_j,K_k

         character(len=30)::aa
        aaa=myid
       write(aa,*) aaa

      if(converge) then
       
         if(myid==0) then

         open(unit=41,file='error.dat',status='old',position='append')
         write(41,*) n_iter,t1,t2
         close(41)

         open(unit=44,file='FreeEG.dat',status='old',position='append')
         write(44,*) n_iter,FE_global,pff_global,tempEE_global
         close(44)

        open(unit=45,file='result.txt',status='replace')  
        write(45,*) "FREE_ENERGY",FE_global,pff_global,tempEE_global,fe_wall
        write(45,*) "box size",SIDEx*dx,SIDEy*dy,SIDEz*dz
        write(45,*) "ava density",density_ava
         write(45,*) "converge",converge 
        close(45)

         endif

  open(unit=23,file= 'RHO_total' //trim(adjustl(aa)) // '.dat',status='replace')
  ! write(23,*) "local_x_start=",local_x_start
         k=0
        do K_i=0,local_nx-1
        do k_j=0,SIDEy-1
# if defined (Dim3)
         do K_k=0,SIDEz-1
          write(23,'(3f10.6,3f12.8)') (K_i+local_x_start)*dx,K_j*dy,K_k*dz,RA(k),RB(k),R_wall_sum(k)
# else /* Dim3 */
          write(23,'(2f10.6,3f12.8)') (K_i+local_x_start)*dx,K_j*dy,RA(k),RB(k),R_wall_sum(k)
# endif /* Dim3 */
            
          k=k+1
         enddo
         enddo
# if defined (Dim3)
        enddo
# endif /* Dim3 */
close(23)


       else


        if(myid==0) then
        
         open(unit=41,file='error.dat',status='old',position='append')
         write(41,*) n_iter,t1,t2
         close(41)
        endif






         if(myid==0) then
         open(unit=44,file='FreeEG.dat',status='old',position='append')
         write(44,*) n_iter,FE_global,pff_global,tempEE_global,fe_wall
         write(44,*) "dx,dy,dz",dx,dy,dz
        write(44,*) "box size",SIDEx*dx,SIDEy*dy,SIDEz*dz
         close(44)
         endif

  open(unit=23,file= 'RHO_total' //trim(adjustl(aa)) // '.dat',status='replace')
  ! write(23,*) "local_x_start=",local_x_start
         k=0
        do K_i=0,local_nx-1
        do k_j=0,SIDEy-1
# if defined (Dim3)
         do K_k=0,SIDEz-1
          write(23,'(3f10.6,3f12.8)') (K_i+local_x_start)*dx,K_j*dy,K_k*dz,RA(k),RB(k),R_wall_sum(k)
# else /* Dim3 */
          write(23,'(2f10.6,3f12.8)') (K_i+local_x_start)*dx,K_j*dy,RA(k),RB(k),R_wall_sum(k)
# endif /* Dim3 */
            
          k=k+1
         enddo
         enddo
# if defined (Dim3)
        enddo
# endif /* Dim3 */
close(23)

endif

 
        

  end subroutine SCFT_dump


  subroutine W_dump()
   USE nrtype,only :DP
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
  implicit none
  integer :: i,j,k,aaa
  REAL*8 :: x,y,z
  character(len=30)::aa
  aaa=myid
  write(aa,*) aaa
  open(unit=23,file= 'W' //trim(adjustl(aa)) // '.dat',status='replace')
  do K=0,LOCAL_SIZE-1
  write(23,*) K,WA(K),WB(K)
  enddo
  close(23)

  

 end subroutine W_dump



















  





















