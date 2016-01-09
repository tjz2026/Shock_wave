

subroutine MDE_q()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 !USE fftw !user defined module
 USE mpi_fftw3_operation
 implicit none
 integer :: i,j,s
 integer :: index_nonzero
 integer :: k
 integer :: istat2
 !for A
 real(DP),allocatable :: exp_WA(:)
 !for B
 real(DP),allocatable :: exp_WB(:)
 REAL(DP) :: M_grid_inv




!for A
 allocate(exp_WA(0:LOCAL_SIZE-1),stat=istat2)


!for B
 allocate(exp_WB(0:LOCAL_SIZE-1),stat=istat2)

  

 if(istat2/=0) then
 write(*,*) "allocate matrix in mde_q failed"
 stop
 endif


 M_grid_inv=1.0d0/M_grid




do k=0,LOCAL_SIZE-1
    qA(0,K)=1.0
enddo

do k=0,LOCAL_SIZE-1
   exp_WA(k)=exp(-ds*wA(k)*0.5)
   exp_WB(k)=exp(-ds*wB(k)*0.5)
enddo
       
 !!!!!
 do s=1,NA

!!!doing the forward FFT

        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(qA(s-1,k)*exp_WA(k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       !call fftwnd_f77_mpi(plancf,1,local_data,work,0,FFTW_NORMAL_ORDER)
       
        do k=0,LOCAL_SIZE-1
            local_data(k+1)=local_data(k+1)*exp_K2(k)
         enddo


       call fftw_mpi_execute_dft(bplan, local_data, local_data)
            
                do k=0,LOCAL_SIZE-1
                 qA(s,k)=real(local_data(k+1))*M_grid_inv*exp_WA(k)
                enddo
         
 enddo   !enddo s=1,NA


do k=0,LOCAL_SIZE-1
    qB(0,K)=qA(NA,k)
enddo

 do s=1+NA,Nmax


!!!doing the forward FFT

        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(qB(s-1-NA,k)*exp_WB(k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       !call fftwnd_f77_mpi(plancf,1,local_data,work,0,FFTW_NORMAL_ORDER)
       
        do k=0,LOCAL_SIZE-1
            local_data(k+1)=local_data(k+1)*exp_K2(k)
         enddo


       call fftw_mpi_execute_dft(bplan, local_data, local_data)
            
                do k=0,LOCAL_SIZE-1
                 qB(s-NA,k)=real(local_data(k+1))*M_grid_inv*exp_WB(k)
                enddo


  enddo  !enddo   s=1+NA,Nmax           


 deallocate(exp_WA)
 deallocate(exp_WB)


end subroutine MDE_q

      
subroutine MDE_qstar()
 USE nrtype,only :DP,PI
 USE global_para
 USE mpi
 USE control
 USE constants
 USE utility
 USE mmpi
 !USE fftw !user defined module
 USE mpi_fftw3_operation
 implicit none
 integer :: i,j,s
 integer :: index_nonzero
 integer :: k
 integer :: istat2
 !for A
 real(DP),allocatable :: exp_WA(:)
 !for B
 real(DP),allocatable :: exp_WB(:)
 REAL(DP) :: M_grid_inv




!for A
 allocate(exp_WA(0:LOCAL_SIZE-1),stat=istat2)


!for B
 allocate(exp_WB(0:LOCAL_SIZE-1),stat=istat2)

  

 if(istat2/=0) then
 write(*,*) "allocate matrix in mde_q failed"
 stop
 endif


 M_grid_inv=1.0d0/M_grid




do k=0,LOCAL_SIZE-1
    qBstar(0,K)=1.0
enddo

do k=0,LOCAL_SIZE-1
   exp_WA(k)=exp(-ds*wA(k)*0.5)
   exp_WB(k)=exp(-ds*wB(k)*0.5)
enddo
       
 !!!!!
 do s=1,NB

!!!doing the forward FFT

        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(qBstar(s-1,k)*exp_WB(k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       !call fftwnd_f77_mpi(plancf,1,local_data,work,0,FFTW_NORMAL_ORDER)
       
        do k=0,LOCAL_SIZE-1
            local_data(k+1)=local_data(k+1)*exp_K2(k)
         enddo


       call fftw_mpi_execute_dft(bplan, local_data, local_data)
            
                do k=0,LOCAL_SIZE-1
                 qBstar(s,k)=real(local_data(k+1))*M_grid_inv*exp_WB(k)
                enddo
         
 enddo   !enddo s=1,NB





do k=0,LOCAL_SIZE-1
    qAstar(0,K)=qBstar(NB,k)
enddo

 do s=1+NB,Nmax


!!!doing the forward FFT

        do k=0,LOCAL_SIZE-1
          local_data(k+1)=cmplx(qAstar(s-1-NB,k)*exp_WA(k),0.0)
         enddo
       
       
       call fftw_mpi_execute_dft(fplan, local_data, local_data)
       !call fftwnd_f77_mpi(plancf,1,local_data,work,0,FFTW_NORMAL_ORDER)
       
        do k=0,LOCAL_SIZE-1
            local_data(k+1)=local_data(k+1)*exp_K2(k)
         enddo


       call fftw_mpi_execute_dft(bplan, local_data, local_data)
            
                do k=0,LOCAL_SIZE-1
                 qAstar(s-NB,k)=real(local_data(k+1))*M_grid_inv*exp_WA(k)
                enddo


  enddo  !enddo   s=1+NB,Nmax           


 deallocate(exp_WA)
 deallocate(exp_WB)


 end subroutine MDE_qstar

