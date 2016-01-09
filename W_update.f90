       subroutine W_update(converge)
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
       integer :: K,j
       real(DP) :: a1,a2,sum_t1,sum_t2,sum_t1_tot,sum_t2_tot
       logical :: converge
       real(DP),allocatable :: i_W_plus_new(:),W_minus_new(:)
       real(DP) :: wall_force

       a1=2.0/(NXab)

       wall_force=0.0
       call mp_barrier()
       if(myid==0) then
       write(*,*) "before update w"
       endif
       call mp_barrier()
       sum_t1=0.0
       sum_t2=0.0  
       do K=0,LOCAL_SIZE-1
       sum_t1=sum_t1+abs(delta_t1*C*(RA(k)+RB(k)+(R_wall_sum(k)-1.0)* &
         (ksi/(ksi+Nxab))-(2/(Nxab+ksi)*i_W_plus(K))))
      i_W_plus(K)=i_W_plus(K)+delta_t1*C*(RA(k)+RB(k)+(R_wall_sum(k)-1.0)* &
         (ksi/(ksi+Nxab))-(2/(Nxab+ksi)*i_W_plus(K)))
        wall_force=(2.0*(XNw_A-XNw_B)/NXAB)*R_wall_bottom(K)
       W_minus(K)=W_minus(K)-delta_t2*C*(a1*W_minus(K)- wall_force -RA(K)+RB(K))
     
       sum_t2=sum_t2+abs(delta_t2*C*(a1*W_minus(K)- wall_force -RA(K)+RB(K)))
       enddo

       call mp_barrier()
       if(myid==0) then
       write(*,*) "af update w"
       endif
       call mp_barrier()
       

              sum_t1_tot=0.0
              sum_t2_tot=0.0
              call mp_barrier()
              call mp_allreduce(sum_t1,sum_t1_tot)
              call mp_allreduce(sum_t2,sum_t2_tot)
              call mp_barrier()
              sum_t1=sum_t1_tot/nprocs
              sum_t2=sum_t2_tot/nprocs
               t1=sum_t1
               t2=sum_t2
             if(myid==0) then
             write(*,*) "error t1,t2",sum_t1,sum_t2
             endif 
             if(max(abs(sum_t1),abs(sum_t2))<CC) then
                     converge=.true.
             else
                     converge=.false.
             endif


            WA(:)=i_W_plus(:)-W_minus(:)
            WB(:)=i_W_plus(:)+W_minus(:)

       end subroutine w_update


       subroutine simple_mixing(converge)
       USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
       USE mpi_fftw3_operation
       USE wall 
       implicit none
       integer :: K
       real(DP) :: a1,a2,sum_t1,sum_t2,sum_t1_tot,aaa
       logical :: converge
       real(DP) :: ta_diff_tot,tb_diff_tot,pressure_coeff,temp1,temp2


       ta_diff=0.0
       tb_diff=0.0


         !!!doing simple mixing for A,B
           do k=0,LOCAL_SIZE-1
               pressure_coeff=0.5*(WA(k)+WB(k)-R_wall_top(k)*(-NXab)- & 
                              R_wall_bottom(k)*(XNw_A+XNw_B-NXab))
               temp1=NXab*(RB(k))+XNw_A*R_wall_bottom(k)+pressure_coeff-WA(k)
               temp2=NXab*(RA(k))+XNw_B*R_wall_bottom(k)+pressure_coeff-WB(k)
               WA(k)=WA(k)+delta_t2*temp1
               WB(k)=WB(k)+delta_t2*temp2
                if(abs(temp1)>ta_diff) ta_diff=abs(temp1)
                if(abs(temp2)>tb_diff) tb_diff=abs(temp2)
           enddo

              ta_diff_tot=0.0
              tb_diff_tot=0.0
              call mp_barrier()
              call mp_allreduce(ta_diff,ta_diff_tot)
              call mp_allreduce(tb_diff,tb_diff_tot)
              call mp_barrier()

              ta_diff=ta_diff_tot/nprocs
              tb_diff=tb_diff_tot/nprocs


             if(myid==0) then
             write(*,*) "error ta_diff,tb_diff",ta_diff,tb_diff
             endif 
             if(ta_diff<CC .and. tb_diff<CC) then
                     converge=.true.
             else
                     converge=.false.
             endif

           i_W_plus(:)=(WA(:)+WB(:))*0.5
           W_minus(:)=(WB(:)-WA(:))*0.5


       end subroutine simple_mixing
