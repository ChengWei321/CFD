subroutine boundary_conditions()
  use variables
  implicit none
  
  
  !---------------------------------------------------!
  !         Boundary conditions calculation           !
  !---------------------------------------------------!
  
  
  !!!!!!!!Velocity boundary conditions!!!!!!!!
  
  !West vertical wall
    do k=-1,nz+2
  do j=-1,ny+2
    u(1,j,k) = 0
    u(0,j,k) = u(1,j,k)
    u(-1,j,k) = u(0,j,k)
    v(0,j,k) = -v(1,j,k)
    v(-1,j,k) = v(0,j,k)
    w(0,j,k) = -w(1,j,k)
    w(-1,j,k) = w(0,j,k)
  end do
    end do
  
  !East vertical wall
    do k=-1,nz+2
  do j=-1,ny+2
    u(nx,j,k) = 0
    u(nx+1,j,k) = u(nx,j,k)
    u(nx+2,j,k) = u(nx+1,j,k)
    v(nx+1,j,k) = -v(nx,j,k)
    v(nx+2,j,k) = v(nx+1,j,k)
    w(nx+1,j,k) = 2*w_east-w(nx,j,k)
    w(nx+2,j,k) = w(nx+1,j,k)
  end do
    end do
  
  !South horizontal wall
    do k=-1,nz+2
  do i=-1,nx+2
    u(i,0,k) =  -u(i,1,k)
    u(i,-1,k) = u(i,0,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    if (ETA(i,1,k)<0.5) then 
      v(i,1,k)=v_south*(1-ETA(i,1,k))
    else
      v(i,1,k)=0
    endif
    v(i,0,k) = v(i,1,k) 
    v(i,-1,k) = v(i,0,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    w(i,0,k) = 2*w_south-w(i,1,k)
    w(i,-1,k) = w(i,0,k)
  end do
    end do
  
  !North horizontal wall
    do k=-1,nz+2
  do i=-1,nx+2
    u(i,ny+1,k) = -u(i,ny,k)
    u(i,ny+2,k) = u(i,ny+1,k) 
    v(i,ny,k) = 0
    v(i,ny+1,k) = v(i,ny,k)
    v(i,ny+2,k) = v(i,ny+1,k)
    w(i,ny+1,k) = 2*w_north-w(i,ny,k)
    w(i,ny+2,k) = w(i,ny+1,k)
  end do
    end do
  
  !Back horizontal wall
    do j=-1,ny+2
  do i=-1,nx+2
    u(i,j,0) =  u(i,j,1)
    u(i,j,-1) = u(i,j,0)
    v(i,j,0) =  v(i,j,1)
    v(i,j,-1) = v(i,j,0)
    w(i,j,0) =  w(i,j,1)
    w(i,j,-1) = w(i,j,0)
  end do
    end do
  
  !Front horizontal wall
    do j=-1,ny+2
  do i=-1,nx+2
    u(i,j,nz+1) = -u(i,j,nz)
    u(i,j,nz+2) = u(i,j,nz+1)
    v(i,j,nz+1) = -v(i,j,nz)
    v(i,j,nz+2) = v(i,j,nz+1)
    w(i,j,nz) =  w_front
    w(i,j,nz+1) = w(i,j,nz)
    w(i,j,nz+2) = w(i,j,nz+1)
  end do
    end do
  
  
  
     do k=-1,nz+2 
    do j=-1,ny+2
  do i=-1,nx+2
  
    u1(i,j,k)=u(i,j,k)
    v1(i,j,k)=v(i,j,k)
    w1(i,j,k)=w(i,j,k)
  
    u2(i,j,k)=u(i,j,k)
    v2(i,j,k)=v(i,j,k)
    w2(i,j,k)=w(i,j,k)
    
    u_star(i,j,k)=u(i,j,k)
    v_star(i,j,k)=v(i,j,k)
    w_star(i,j,k)=w(i,j,k)
  
  
  end do
    end do
      end do
      
  !!!!!!!Pressure boundary conditions!!!!!!!
  
  !West vertical wall
    do k=-1,nz+2 
  do j=-1,ny+2
    p(0,j,k)=p(1,j,k)
    p(-1,j,k)=p(0,j,k)
  end do
    end do
  
  !East vertical wall
    do k=-1,nz+2
  do j=-1,ny+2
     p(nx+1,j,k)=p(nx,j,k)
     p(nx+2,j,k)=p(nx+1,j,k)
  end do
    end do
  
  !South horizontal wall
    do k=-1,nz+2
  do i=-1,nx+2
    p(i,0,k)=p(i,1,k)
    p(i,-1,k)=p(i,0,k)
  end do
    end do
    
  !North horizontal wall
    do k=-1,nz+2
  do i=-1,nx+2
    p(i,ny+1,k)=p(i,ny,k)
    p(i,ny+2,k)=p(i,ny+1,k)
  end do
    end do
  
  !Back horizontal wall
    do j=-1,ny+2
  do i=-1,nx+2
    if (ETA(i,j,1)<0.5)then
     p(i,j,0)=p_back
    else
     p(i,j,0)=p(i,j,1)
    endif
     p(i,j,-1)=p(i,j,0)
  end do
    end do
  
  !Front horizontal wall
    do j=-1,ny+2
  do i=-1,nx+2
     p(i,j,nz+1)=p(i,j,nz)
     p(i,j,nz+2)=p(i,j,nz+1)
  end do
    end do
  
  end subroutine boundary_conditions
