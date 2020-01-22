subroutine read_ETA()
  use variables
  implicit none
  
  real*8 ,dimension(nx*ny*nz):: A
  integer :: inumber,lk
  inumber = nx*ny*nz
  
  open(1,file='ETA.DAT')
  
  do i=1,inumber
  
    read(1,*) A(i)
    
  end do
  
  close(1)
  lk=1
  do k=1,nz;do j=1,ny;do i=1,nx
  
  ETA(i,j,k)=A(lk)
  lk=lk+1
  
  end do;end do;end do

  
end subroutine read_ETA

subroutine DFIB_Sphere()
use variables
implicit none
real*4 :: DISTANCE,SDIST
real*4 :: DIAGONAL
integer :: l ,m ,n
real*4 :: xi
real*4 :: dxg, dyg, dzg
real*4 ,dimension(1:nSubGrids+1) :: SX
real*4 ,dimension(1:nSubGrids+1) :: SY
real*4 ,dimension(1:nSubGrids+1) :: SZ

!open(88,file='ETA.dat')

do k=kBgnVOS,kEndVOS 
  do j=jBgnVOS,jEndVOS
    do i=iBgnVOS,iEndVOS

         DIAGONAL = sqrt( iDx(i)*iDx(i) + iDy(j)*iDy(j) + iDz(k)*iDz(k) ) / 2.D0
         DISTANCE = sqrt((((X(i)+X(i+1))/2.D0)-xc)**2.D0+ &
                         (((Y(j)+Y(j+1))/2.D0)- ( yc + dis_Y ) )**2.D0+ &
                         (((Z(k)+Z(k+1))/2.D0)-zc)**2.D0  )

           if( abs(DISTANCE - r) < DIAGONAL ) then
              
              do l=1,nSubGrids+1,1
              do m=1,nSubGrids+1,1
              do n=1,nSubGrids+1,1
                dxg = iDx(i) / nSubGrids
                dyg = iDy(j) / nSubGrids
                dzg = iDz(k) / nSubGrids
                SX(n)=X(i)+(n-1)*dxg
                SY(m)=Y(j)+(m-1)*dyg
                SZ(l)=Z(k)+(l-1)*dzg
              end do
              end do
              end do

              xi = 0.0
              do l=1,nSubGrids
                do m=1,nSubGrids
                  do n=1,nSubGrids
                    
                    SDIST=SQRT((((SX(n)+SX(n+1))/2.D0)-xc)**2.D0+ &
                               (((SY(m)+SY(m+1))/2.D0)- ( yc + dis_Y ) )**2.D0+ &
                               (((SZ(l)+SZ(l+1))/2.D0)-zc)**2.D0  )
                    if(SDIST <= r) then
                      xi = xi + 1
                    end if

                  end do
                end do
              end do   
              ETA(i,j,k) = xi / (nSubGrids*nSubGrids*nSubGrids)
              
              !print*,'ETA(i,j,k) = ',ETA(i,j,k),'xi = ',xi,i,j,k
           
           else if (DISTANCE <= r) then
             ETA(i,j,k)=1.D0
           else
             ETA(i,j,k)=0.D0
           end if
                
    end do
  end do
end do



end subroutine DFIB_Sphere

subroutine DFIB_Pipe()
  use variables
  implicit none
  real*8 :: DISTANCE,SDIST
  real*8 :: DIAGONAL
  integer :: l ,m ,n
  real*8 :: xi
  real*8 :: dxg, dyg, dzg
  real*8 ,dimension(1:nSubGrids+1) :: SX
  real*8 ,dimension(1:nSubGrids+1) :: SY
  real*8 ,dimension(1:nSubGrids+1) :: SZ
  open(88,file='ETA.dat')
  
  
  
  k=kBgnVOS
  do j=jBgnVOS,jEndVOS
  do i=iBgnVOS,iEndVOS
  
     DIAGONAL = sqrt(  iDx(i)*iDx(i) + iDy(j)*iDy(j) ) /2.D0
     DISTANCE = sqrt((((X(i)+X(i+1))/2.D0)-xc)**2.D0+(((Y(j)+Y(j+1))/2.D0)-yc)**2.D0)  
  
       if( abs(DISTANCE - r) < DIAGONAL ) then
        
        xi = 0.0
        
        do l=1,nSubGrids+1,1
        do m=1,nSubGrids+1,1
        do n=1,nSubGrids+1,1
          dxg = iDx(i) / nSubGrids
          dyg = iDy(j) / nSubGrids
          dzg = iDz(k) / nSubGrids
          SX(n)=X(i)+(n-1)*dxg
          SY(m)=Y(j)+(m-1)*dyg
          SZ(l)=Z(k)+(l-1)*dzg
        end do
        end do
        end do
        
  
        do l=1,nSubGrids,1
        do m=1,nSubGrids,1
        do n=1,nSubGrids,1
            SDIST=SQRT((((SX(n)+SX(n+1))/2.D0)-xc)**2.D0+(((SY(m)+SY(m+1))/2.D0)-yc)**2.D0)
            if(SDIST <= r) then
            xi = xi + 1
            end if
        end do
        end do
        end do   
        
        ETA(i,j,k) = 1 - xi / (nSubGrids*nSubGrids*nSubGrids)
        
        !print*,'ETA(i,j,k) = ',ETA(i,j,k),'xi = ',xi,i,j,k
       
       else if (DISTANCE > r) then
       ETA(i,j,k)=1.D0
       else
       ETA(i,j,k)=0.D0
       end if
        
  end do
  end do
  
  do k=kBgnVOS+1,kEndVOS
  do j=jBgnVOS,jEndVOS
  do i=iBgnVOS,iEndVOS
    ETA(i,j,k) = ETA(i,j,kBgnVOS)
  end do
  end do
  end do    

end subroutine DFIB_Pipe