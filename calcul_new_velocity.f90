subroutine calcul_new_velocity()
use variables
implicit none


!-------------------------------------------------!
!    Calculation of velocity field at t = dt*n+1  !
!-------------------------------------------------!

!!! In x direction !!!


do k=istart,iend
do j=1,ny; do i=1,nx

u1(i,j,k) = u_star(i,j,k) - dt*drho*(p(i+1,j,k)-p(i,j,k)) / Dxs(i)

enddo; enddo
enddo


!!! In y direction !!!


do k=istart,iend
do j=1,ny; do i=1,nx

v1(i,j,k) = v_star(i,j,k) - dt*drho*(p(i,j+1,k)-p(i,j,k)) / Dys(j)

enddo; enddo
enddo


!!! In w direction !!!


do k=istart,iend
do j=1,ny; do i=1,nx

w1(i,j,k) = w_star(i,j,k) - dt*drho*(p(i,j,k+1)-p(i,j,k)) / Dzs(k)

enddo; enddo
enddo




!-------------------------------------------------!
!    Calculation of velocity field for DFIB       !
!-------------------------------------------------!

!!! In x direction !!!


do k=istart,iend
do j=1,ny; do i=1,nx

u2(i,j,k) = ETA(i,j,k) *u_solid + (1-ETA(i,j,k)) * u1(i,j,k)
FX(i,j,k) = (u2(i,j,k) - u1(i,j,k)) / dt

enddo; enddo
enddo


!!! In y direction !!!


do k=istart,iend
do j=1,ny; do i=1,nx

v2(i,j,k) = ETA(i,j,k) *v_solid_y + (1-ETA(i,j,k)) * v1(i,j,k)
FY(i,j,k) = (v2(i,j,k) - v1(i,j,k)) / dt

enddo; enddo
enddo


!!! In w direction !!!


do k=istart,iend
do j=1,ny; do i=1,nx

w2(i,j,k) = ETA(i,j,k) *w_solid + (1-ETA(i,j,k)) * w1(i,j,k)
FZ(i,j,k) = (w2(i,j,k) - w1(i,j,k)) / dt

enddo; enddo
enddo




!----------data collect among nodes----------!
icount = gcount(myid)*(nx+4)*(ny+4)

!Send my results back to the master
if(myid>master)then

    
    itag = 350
    call MPI_SEND( FX(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
    itag = 360
    call MPI_SEND( FY(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
    itag = 370
    call MPI_SEND( FZ(-1,-1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )

end if


!Wait to receive results from each task
if(myid==master)then

    do i = 1, (nproc-1)

            itag = 350
            call MPI_RECV( FX(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 360
            call MPI_RECV( FY(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 370
            call MPI_RECV( FZ(-1,-1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            

    end do


end if



!icount = (nx+4)*(ny+4)*(nz+4)
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!call MPI_BCAST( FX(-1,-1,-1), icount, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
!call MPI_BARRIER(MPI_COMM_WORLD, ierr
!call MPI_BCAST( FY(-1,-1,-1), icount, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!call MPI_BCAST( FZ(-1,-1,-1), icount, MPI_REAL8, master, MPI_COMM_WORLD, ierr )
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!----------data collect among nodes----------!
      


!---------------------------------------------------------------------- -----!
!      COEFFICIENTS: CD and Center_CL,Distance,Velocity for solid            !
!----------------------------------------------------------------------------!

! only cpu 0 answer correct
call virtualForceIntegrator()


if (solidPos == 2) then

    call Runge_Kutta_VIV()
    
    if(myid==master)then
        open(88,file='distance_y.dat',position='append')
        write(88,*) time ,( (yc + dis_Y) - yc )
        close(88)
    end if

end if

!----------------------------------------------------------------------------------------!
!           Centering variables to get matrix (1:nx,1:ny,1:nz) for Tecplot               !
!----------------------------------------------------------------------------------------!

!----------data transformation among nodes----------!
icount = (nx+4)*(ny+4)

itag = 380
call MPI_SENDRECV( w2(-1,-1,iend), icount, MPI_REAL8, r_nbr, itag, &
                   w2(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
            
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!----------data transformation among nodes----------!


do k=istart,iend
do j=1,ny; do i=1,nx

uc(i,j,k) = 0.5*(u2(i,j,k)+u2(i-1,j,k))
vc(i,j,k) = 0.5*(v2(i,j,k)+v2(i,j-1,k))
wc(i,j,k) = 0.5*(w2(i,j,k)+w2(i,j,k-1))

enddo; enddo
enddo


end subroutine calcul_new_velocity










subroutine Updating_velocity()
use variables
implicit none
!---------------------------------------------------------!
!    loops to update the main u and v velocity fields     !
!---------------------------------------------------------!

!!! In x direction !!!


do k=istart,iend
do j=1,ny; do i=1,nx

    u(i,j,k) = u2(i,j,k)    

enddo; enddo
enddo


!!! In y direction !!!

do k=istart,iend
do j=1,ny; do i=1,nx

    v(i,j,k) = v2(i,j,k)  

enddo; enddo
enddo


!!! In w direction !!!

	
do k=istart,iend
do j=1,ny; do i=1,nx

    w(i,j,k) = w2(i,j,k) 

enddo; enddo
enddo 


end subroutine Updating_velocity


