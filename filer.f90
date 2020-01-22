subroutine filer_final()
use variables
implicit none

   open (1,file='information.dat',position='append')
   write(1,*) ' ' 
   write(1,*) '======================================================='
   write(1,*) 'total cost time = ',totalcosttime,'sec'
   write(1,*) 'final time = ',time,'sec'
   write(1,*) '======================================================='
!   write(1,*) '======================================================='
!   write(1,*) 'Pressure Gradient', GradientP
!   write(1,*) 'Pressure Gradient Difference', DifferenceGradientP
!   write(1,*) 'Max Velocity Difference', Diff_Max_U
!   write(1,*) 'Velocity Difference(RMS)', U_RMS
!   write(1,*) 'Velocity Difference (L2norm)', U_L2norm
!   write(1,*) '======================================================='
!   write(1,*) ' ' 
!   close(1)
!   open (25,file='Ans_U.dat',position='append')   
!   write(25,*) 'TITLE = ','"TIME = ',time,'"'
!   write(25,*) 'VARIABLES = "y","w_ans","w"'
!   write(25,*) ' ZONE ','j=', ny
!   write(25,*) 'DT = (DOUBLE, DOUBLE,DOUBLE)'
!   write(25,*) 'DATAPACKING = POINT'
!      do j=1,ny
!         write (25,'(E12.5,a,E12.5,a,E12.5)')  &
!         Ys(j),'      ', profileU(j),'      ', wc(nx/2,j,1)
!      enddo
!   close(25)
!   open (26,file='FullyDevelopedCheck.dat',position='append')
!   write(26,*) 'TITLE = ','"TIME = ',time,'"'
!   write(26,*) 'VARIABLES = "y","w0m","w1m","w2m","w3m","w4m","w5m","w6m","w7m","w8m","w9m","w10m"'
!   write(26,*) ' ZONE ','j=', ny
!   write(26,*) 'DT = (DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE)'
!   write(26,*) 'DATAPACKING = POINT'
!      do j=1,ny
!         write (26,'(E12.5,a,E12.5,a,E12.5,a,E12.5,a,E12.5,a,E12.5,a,E12.5,a,E12.5,a,E12.5,a,E12.5,a,E12.5,a,E12.5)')  &
!            Ya(j),'      ', wc(nx/2,j,1),'      ',wc(nx/2,j,20),'      ',&
!            wc(nx/2,j,40), '      ',wc(nx/2,j,60),'      ', wc(nx/2,j,80),&
!            '      ', wc(nx/2,j,100),'      ',wc(nx/2,j,120),'      ',&
!            wc(nx/2,j,140),'      ',wc(nx/2,j,160),'      ', wc(nx/2,j,180),&
!            '      ', wc(nx/2,j,200)  
!      enddo
!   close(26)

   do k=1,nz; do j=1,ny; do i=1,nx
      Qout(i,j,k,1)=pre(i,j,k)
      Qout(i,j,k,2)=uc(i,j,k)
      Qout(i,j,k,3)=vc(i,j,k)
      Qout(i,j,k,4)=wc(i,j,k)
      Qout(i,j,k,5)=ETA(i,j,k)
   enddo; enddo; enddo

	write(filename,'(I4.4)')num
 	fileformat = '.q'

   open (19,file='finaltimestep.q',position='append',form='unformatted')
   write(19) nblocks
   write(19) nx, ny, nz
   write(19) temp, temp, temp, temp
   write(19) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
   close(19)
end subroutine filer_final

subroutine filerInfo()
use variables
implicit none
   open (1,file='information.dat',position='append')

   write(1,*) '!-----------------------------------------------------------!'
   write(1,*) '!  Finite Volume Method by projection method with mpi       !'
   write(1,*) '!-----------------------------------------------------------!'
   write(1,*) ' '
   write(1,*) '======================================================='
   write(1,*) 'number of processor = ',nproc
   write(1,*) 'number of grid = ',nx*ny*nz
   write(1,*) 'SOR parameter = ',omega
   write(1,*) 'nx = ',nx
   write(1,*) 'ny = ',ny
   write(1,*) 'nz = ',nz
   write(1,*) 'nxSml = ',nxSml
   write(1,*) 'nySml = ',nySml
   write(1,*) 'nzSml = ',nzSml
   write(1,*) 'x-direction length = ',lx
   write(1,*) 'y-direction length = ',ly
   write(1,*) 'z-direction length = ',lz
   write(1,*) 'x-direction small length = ',lxSml
   write(1,*) 'y-direction small length = ',lySml
   write(1,*) 'z-direction small length = ',lzSml
   write(1,*) 'Large dx = ',dx
   write(1,*) 'Large dy = ',dy
   write(1,*) 'Large dz = ',dz
   write(1,*) 'Small dx = ',dxSml
   write(1,*) 'Small dy = ',dySml
   write(1,*) 'Small dz = ',dzSml
   write(1,*) 'Small grid xc = ',GridderXc
   write(1,*) 'Small grid yc = ',GridderYc
   write(1,*) 'Small grid zc = ',GridderZc
   write(1,*) 'lxSml  = ',lxSml
   write(1,*) 'lySml  = ',lySml
   write(1,*) 'lzSml  = ',lzSml
   write(1,*) 'DFIB xc = ',xc
   write(1,*) 'DFIB yc = ',yc
   write(1,*) 'DFIB zc = ',zc
   write(1,*) 'Kinematic Viscosity = ',nu
   write(1,*) 'density = ',rho
   write(1,*) 'Renolds Number = ',u_west*2*sqrt(2.218e-4*pi)/nu
   if (Q == 1) then; write(1,*) 'Velocity scheme = Upwind'
   else if (Q == 2) then; write(1,*) 'Velocity scheme = Centre'
   else if (Q == 3 .or. Q == 4) then; write(1,*) 'Velocity scheme = QUICK'
   end if
   if (solidPos == 1) then; write(1,*) 'solidPos = fixed'
   else if(solidPos == 2) then; write(1,*) 'solidPos = moving'
   else if(solidPos == 0) then; write(1,*) 'solidPos = no solid'
   end if
   if (steadiness == 1) then; write(1,*) 'steadiness = steady'
   else if(steadiness == 2) then; write(1,*) 'steadiness = unsteady'
   end if
   write(1,*) '======================================================='
   write(1,*) ' '
   write(1,*) '======================================================='
   if (Q1 == 1) then; write(1,*) 'LES mode on'
   else if(Q1 == 0) then; write(1,*) 'LES mode off'
   end if
   write(1,*) 'LES Cs = ',Cs
   write(1,*) '======================================================='
   write(1,*) ' '
   write(1,*) '======================================================='
   write(1,*) 'dt = ',dt
   write(1,*) 'max time step = ',nstep*dt,'sec'
   write(1,*) 'each time step = ',isto*dt,'sec'
   write(1,*) '======================================================='
   write(1,*) ' '
   write(1,*) '======================================================='
   write(1,*) 'p Residual =',zeta
   write(1,*) 'velocity Residual =',zeta_vel
   write(1,*) '======================================================='

   close(1)
end subroutine filerInfo

subroutine filerProcess()
use variables
implicit none
   real*8 :: px ,py ,pz
   real*8 :: px1,py1,pz1
   real*8 :: px2,py2,pz2
   real*8 :: px3,py3,pz3
   real*8 :: px4,py4,pz4
   
     
!      px = (Xc - lxSml/2.0) / dx + (lxSml/2) / dxSml
!      py = (Yc - lySml/2.0) / dy + (lySml/2) / dySml 
!      pz = (Zc - lzSml/2.0) / dz +  (1.8) / dzSml

   
   open (2,file='process.dat',position='append')
   write(2,*) time ,wc(86,76,13), p(86,76,13)
   close(2)

!   open (3,file='CD.dat',position='append')
!   write(3,*) time ,cDrag,cLift
!   close(3)



   if(steadiness==2 .AND. istep>=StartAvg)then
      open (4,file='CD_unsteady.dat',position='append')
      write(4,*) time ,cDragAvg,cLiftAvg,cDrag_rms,cLift_rms, &
                     totalFXAvg,totalFYAvg
      close(4)
   end if
   
end subroutine filerProcess

subroutine filereachtime()
use variables
implicit none

   do k=1,nz; do j=1,ny; do i=1,nx
      Qout(i,j,k,1)=pre(i,j,k)
      Qout(i,j,k,2)=uc(i,j,k)
      Qout(i,j,k,3)=vc(i,j,k)
      Qout(i,j,k,4)=wc(i,j,k)
      Qout(i,j,k,5)=ETA(i,j,k)
   enddo; enddo; enddo

	write(filename,'(I4.4)')num
 	fileformat = '.q'

   open (18,file=TRIM(filename)//fileformat,position='append',form='unformatted')
   write(18) nblocks
   write(18) nx, ny, nz
   write(18) temp, temp, temp, temp
   write(18) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
   close(18)
   
   num = num + 1
   
end subroutine filereachtime

subroutine filerweird()
   use variables
   implicit none
   
      do k=1,nz; do j=1,ny; do i=1,nx
         Qout(i,j,k,1)=pre(i,j,k)
         Qout(i,j,k,2)=uc(i,j,k)
         Qout(i,j,k,3)=vc(i,j,k)
         Qout(i,j,k,4)=wc(i,j,k)
         Qout(i,j,k,5)=ETA(i,j,k)
      enddo; enddo; enddo
      
      open (17,file='weird_steps.q',position='append',form='unformatted')
      write(17) nblocks
      write(17) nx, ny, nz
      write(17) temp, temp, temp, temp
      write(17) ( ( ( ( Qout(i,j,k,h), i = 1, nx), j = 1, ny), k = 1, nz), h = 1, 5 )
     
      close(17)
      
      num = num + 1
      
   end subroutine filerweird
