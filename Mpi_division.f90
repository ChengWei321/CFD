subroutine Mpi_division()
use variables 
use mpi
implicit none



    ! 1-D for MPI division
    Zdv = (nz) / nproc
    Zr  = (nz) - Zdv * nproc 
    

    do i = 0,nproc-1
        if(i < Zr) then
            gstart(i) = 1 + i * (Zdv+1)
            gend0(i) = gstart(i) + Zdv
        else
            gstart(i) = 1 + i * Zdv + Zr
            gend0(i) = gstart(i) + Zdv - 1
        end if
        
        gcount(i) = gend0(i) - gstart(i) + 1
        gend(i) = gcount(i) + 2

        !write(*,*) i, gstart(i), gcount(i), gend(i),gend0(i) 
    end do


    

end subroutine Mpi_division

