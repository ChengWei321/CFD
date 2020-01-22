subroutine initial_conditions()
use variables
implicit none

!---------------------------------------------------!
!         Initial conditions calculation            !
!---------------------------------------------------!

do k=-1,nz+2; do j=-1,ny+2; do i=-1,nx+2
!    if (k>70 .and. ETA(i,j,k)<0.0001) then
!        u(i,j,k) = 0.5
!    else
        u(i,j,k) = 0
!    endif
    
   v(i,j,k) = 0

!   if (k<70 .and. ETA(i,j,k)<0.0001) then
!        w(i,j,k) = -0.5
!   else
        w(i,j,k) = 0
!   endif

   u1(i,j,k) = 0
   v1(i,j,k) = 0
   w1(i,j,k) = 0
   u2(i,j,k) = 0
   v2(i,j,k) = 0
   w2(i,j,k) = 0
   p(i,j,k) = 0
   u_star(i,j,k) = 0
   v_star(i,j,k) = 0
   w_star(i,j,k) = 0
   ETA(i,j,k) = 1


end do; end do; end do

end subroutine initial_conditions



