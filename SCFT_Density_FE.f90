       subroutine density
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
       integer :: i,j,s,K
       real(DP) :: pff_temp_global,pff_temp,totden_global,bigQ
       real(DP) :: RAtot,RAtot_global,density_ava_global
       REAL(DP),allocatable :: ar(:)
       REAL(DP),allocatable :: br(:)
       REAL(DP),allocatable :: denz(:)
       REAL(DP),allocatable :: c_A(:)
       REAL(DP),allocatable :: c_B(:)
       REAL(DP) :: sum_iA,sum_iB,sum_half,sum_end_B
       REAL(DP) :: nor_coeff
       REAL(DP) :: sumbr,sumbr_global


       allocate(ar(0:LOCAL_SIZE-1))
       
       do K=0,LOCAL_SIZE-1
         ar(K)=qB(NB,K)
       enddo

       RAtot=sum(ar)
       RAtot_global=0.0 
       density_ava_global=0.0
# if defined (Dim3)       
       pff_temp=simposon_3D_1D_mpi(SIDEz,SIDEy,local_nx,dz,dy,dx,ar)
       density_ava=simposon_3D_1D_mpi(SIDEz,SIDEy,local_nx,dz,dy,dx,R_wall_sum)
# else /* Dim3 */
       pff_temp=simposon_2D_1D_mpi(SIDEy,local_nx,dy,dx,ar)
       density_ava=simposon_2D_1D_mpi(SIDEy,local_nx,dy,dx,R_wall_sum)
# endif  /* Dim3 */




              density_ava_global=0.0
              pff_temp_global=0.0
              call mp_barrier()
              call mp_allreduce(pff_temp,pff_temp_global)
              call mp_allreduce(density_ava,density_ava_global)
              call mp_barrier()

        density_ava=density_ava_global/M_v
        density_ava=1.0-density_ava
        bigQ=pff_temp_global/M_v
        pff_global=bigQ

        if(myid==0) then
         write(*,*) "arsum",RAtot_global,"pff before log",pff_global
         endif    
    
        pff_global=log(pff_global)


        deallocate(ar)

        allocate(denz(0:LOCAL_SIZE-1))
        allocate(c_A(0:NA))
        allocate(c_B(0:NB))

       do K=0,LOCAL_SIZE-1
            do s=0,NA
              c_A(s)=qA(s,K)*qAstar(NA-s,K) 
            enddo
          RA(K)=simposon_1D_NR(0,NA,ds,c_A)*(density_ava/bigQ)

              
            do s=0,NB
              c_B(s)=qB(s,K)*qBstar(NB-s,K) 
            enddo
          RB(K)=simposon_1D_NR(0,NB,ds,c_B)*(density_ava/bigQ)
      
       enddo

  deallocate(denz)
  deallocate(c_A)
  deallocate(c_B)


     
   end subroutine density


   subroutine free_energy()
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
    integer :: j,i,K
    real*8 :: tempEE,NXAB_inv,Xw_XAB,wall_force,aa,bb,fe_wall_tot
    real*8,allocatable :: ar(:)
    real*8,allocatable :: temp(:)

   NXAB_inv=1.0/NXAB
   allocate(ar(0:LOCAL_SIZE-1))
   allocate(temp(0:LOCAL_SIZE-1))
    
  do K=0,LOCAL_SIZE-1
  wall_force=(2.0*(XNw_A-XNw_B)/NXab)*R_wall_bottom(K)       !zhoupan change
  aa=(-2.0*ksi/(Nxab+2*ksi))*(1.0-R_wall_sum(k))*i_W_plus(k)
  temp(k)=aa 
  bb=(1.0/(Nxab+ksi))*i_W_plus(k)**2
  ar(K)=NXab_inv*W_minus(k)**2 + wall_force*W_minus(K)+aa-bb
  enddo
# if defined (Dim3)   
   tempEE=simposon_3D_1D_mpi(SIDEz,SIDEy,local_nx,dz,dy,dx,ar)
   fe_wall=simposon_3D_1D_mpi(SIDEz,SIDEy,local_nx,dz,dy,dx,temp)
# else /* Dim3 */
   tempEE=simposon_2D_1D_mpi(SIDEy,local_nx,dy,dx,ar)
   fe_wall=simposon_2D_1D_mpi(SIDEy,local_nx,dy,dx,temp)
# endif /* Dim3 */

   tempEE_global=0.0
   fe_wall_tot=0.0d0
              call mp_barrier()
              call mp_allreduce(tempEE,tempEE_global)
              call mp_allreduce(fe_wall,fe_wall_tot)
              call mp_barrier()
        tempEE_global=tempEE_global/M_v
        fe_wall=fe_wall_tot/M_v
        tempEE=tempEE_global

 FE_global=-pff_global*density_ava+tempEE_global
 if(myid==0) then
 write(*,*) "FREE energy=",FE_global,pff_global*density_ava,tempEE_global,fe_wall
 endif
     deallocate(ar)


end subroutine free_energy    
     
    
    
     




 




                
                 
        








   
