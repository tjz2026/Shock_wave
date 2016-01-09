       subroutine Wall_density()
       USE nrtype,only :DP
       USE global_para
       USE mpi
       USE control
       USE constants
       USE utility
       USE mmpi
       USE mpi_fftw3_operation
       implicit none
       real(DP) :: dzw
       integer :: k,k_i,k_j,k_k


         k=0
        do K_i=0,local_nx-1
        do k_j=0,SIDEy-1
         do K_k=0,SIDEz-1
         
        dzw=min((K_k*dz-0.0),(dz*SIDEz-dz*K_k))

        RWALL(K)=0.5*(1.0+tanh(z_d*dzw))
        k=k+1
         enddo
         enddo
        enddo
        RWALL=RWALL-0.49999d0
        RWALL=RWALL*2.0d0
        RWALL=-1.0*RWALL+1.0

       end subroutine Wall_density
