
       FUNCTION cal_scft(tz_scal)
       USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
       !USE fftw !user defined module
       USE mpi_fftw3_operation
       use wall
       implicit none
       REAL(DP),intent(inout) :: tz_scal
       logical :: converge
       REAL(DP) :: cal_scft


       if(abs(WAB_init_type)==4) then
       !6 and -6 for HEX ,-6 means input from file
       dx=tz_scal
       dy=dx*SIDEx/(sqrt(3.0)*SIDEy)
# if defined (Dim3)
       dz=dx*SIDEx/SIDEz
# endif /* Dim3 */
       else
       dx=Lx/SIDEx
       dy=Ly/SIDEy
     !   dx=tz_scal
     !   dy=tz_scal
# if defined (Dim3)
       dz=Lz/SIDEz
# endif /* Dim3 */
       endif

# if defined (Dim3)
       M_v=dx*dy*dz*SIDEx*SIDEy*SIDEz
# else /* Dim3 */
       M_v=dx*dy*SIDEx*SIDEy
# endif /* Dim3 */
       if(myid==0) then 
       write(*,*) "M_v=",M_v
       write(*,*) "dx,dy,dz",dx,dy,dz
       endif
       converge=.false.
       n_iter=1

           
      call  exponents_Karray()
      call init_W()


      call wall_builder()
      call wall_dump()
      call mp_barrier()


      if(myid==0) then
       write(*,*) "done init exp_K2"
      endif
      do while (converge/=.true. .and. n_iter<=MAXITS)
         if(myid==0) then
         write(*,*) "the",n_iter,"th iteration","on",myid
         endif
         call MDE_q()
         call MDE_qstar()
         call density()
         call free_energy()
         call w_update(converge)
        ! call simple_mixing(converge)
         if(mod(n_iter,50)==0) then
         call SCFT_dump(converge)
         endif

         !call mp_barrier()
         !if(n_iter>100) then
         !stop
         !endif
        ! call W_dump()
         
         n_iter=n_iter+1

      enddo

      

        cal_scft=FE_global 
        if(myid==0) then
        write(*,*) "exit the calscft","dx=",dx,"converge=",converge
        write(*,*) "FreeEG:",FE_global,"the",n_iter,"iteration"
        write(*,*) "t1:",t1,"t2",t2
        endif

       call SCFT_dump(converge)
       call w_dump()


      end FUNCTION cal_scft


