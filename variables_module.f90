MODULE variables
use mpi
implicit none 

!This module gathers all structures and variable declarations. 
!The different variables no need to be redefined and can this file can just be called in the subroutines.

! Properties
real*8 ,parameter :: mu = 1.87e-5 !Re = 1/nu 
real*8 ,parameter :: rho = 1.217
real*8 ,parameter :: nu = mu / rho !Re = 1/nu 
real*8 ,parameter :: dt = 1.0e-5

! respiratory volume flow rate
real*8 ,parameter :: FR = 261 * (1e-6) ! m^3/s   (ml/s * (m^3/cm^3) )
real*8 ,parameter :: inletA = 3.080168750000000e-04

! Grid definitions
integer ,parameter :: nx=151
integer ,parameter :: ny=116
integer ,parameter :: nz=230
real*8 ,parameter :: lx=151*0.25e-3
real*8 ,parameter :: ly=116*0.25e-3
real*8 ,parameter :: lz=230*0.25e-3
!------------------- OPENMP ------------------------!
integer :: nthreads
!------------------- OPENMP ------------------------!

!----------------------------MPI----------------------------!
integer :: nproc, myid, ierr, dest
integer ::  status(MPI_STATUS_SIZE)
integer ,parameter :: master=0
integer :: Zdv, Zr
integer ,dimension(1:40)::gstart ,gend, gend0, gcount
integer :: l_nbr, r_nbr, icount, iend, istart, itag, igcount
!----------------------------MPI----------------------------!


!--------------------- Unequal grid ---------------------!
real*8 ,parameter :: GridderD = 1
real*8 ,parameter :: GridderXc = 5
real*8 ,parameter :: GridderYc = 5
real*8 ,parameter :: GridderZc = 6.5 
integer ,parameter :: nxSml = 16
integer ,parameter :: nySml = 200
integer ,parameter :: nzSml = 200
real*8 ,parameter :: lxSml  = 4.0*GridderD
real*8 ,parameter :: lySml  = 4.0*GridderD
real*8 ,parameter :: lzSml  = 4.0*GridderD
!--------------------- Unequal grid ---------------------!

!--------------------- output ---------------------!
integer,parameter :: nblocks = 1
real,dimension(1:nx,1:ny,1:nz):: Xout, Yout, Zout
real,dimension(1:nx,1:ny,1:nz,5):: Qout
real :: temp = 1.0    ! mach, alpha, reyn, time 
integer :: h
!--------------------- output ---------------------!

!---------------------- input ---------------------!
integer :: inblocks
integer :: inx
integer :: iny
integer :: inz 
character(len=20) :: inputfile
!---------------------- input ---------------------!


integer :: ik, k, i, j, itmax, isto, istep, nstep, breaker
integer ,parameter :: ndim = nx*ny*nz,mdim=4
integer :: num
real*8 :: totalstarttime, totalfinaltime, totalcosttime
character(len=20) :: filename ,fileformat
integer :: solidPos,steadiness,Q,Q1
real*8 :: zeta, time, ddx, ddt, ddz, ddy, nu1, drho, sum ,VelocityDifference,zeta_vel
real*8 :: dxSml, dySml, dzSml
real*8 :: dx, dy, dz
real*8 ,dimension(-1:nx+2,-1:ny+2,-1:nz+2) :: p,u,v,w,u1,v1,w1,u_star,v_star,w_star, last_velocity
real*8 ,dimension(1:nx,1:ny,1:nz):: div,uc,vc,wc,pre

real*8 ,dimension(-1:nx+2,-1:ny+2,-1:nz+2) :: p_sum,u_sum,v_sum,w_sum
real*8 ,dimension(-1:nx+2,-1:ny+2,-1:nz+2) :: p_avg,u_avg,v_avg,w_avg
real*8 ,dimension(1:nx,1:ny,1:nz) :: uc_avg,vc_avg,wc_avg
integer :: StartAvg

!---------------------------Gauss Seidel-----------------------!
real*8 :: pNew, pChange, mChange, omega, pChangeMax, mChangeMax 
real*8 :: pChangeMax_
!---------------------------Gauss Seidel-----------------------!

!---------------------------B.C value-----------------------!
real*8 :: p_east=0 , p_west=0 , p_north=0 , p_south=0 , p_front=0 , p_back=0
real*8 :: u_east=0 , u_west=0 , u_north=0 , u_south=0 , u_front=0 , u_back=0 
real*8 :: v_east=0 , v_west=0 , v_north=0 , v_south=FR/inletA , v_front=0 , v_back=0 
real*8 :: w_east=0 , w_west=0 , w_north=0 , w_south=0 , w_front=0 , w_back=0.
!---------------------------B.C value-----------------------!

!---------------------------centre--------------------------!
real*8 :: u_tilde_x1, u_tilde_x2, u_tilde_y1, u_tilde_y2, u_tilde_z1, u_tilde_z2 
real*8 :: v_tilde_x1, v_tilde_x2, v_tilde_y1, v_tilde_y2, v_tilde_z1, v_tilde_z2
real*8 :: w_tilde_x1, w_tilde_x2, w_tilde_y1, w_tilde_y2, w_tilde_z1, w_tilde_z2
!---------------------------centre--------------------------!

!---------------------------QUICK---------------------------!
real*8 :: ue, uw, un, us, uf, ub, vnu, vsu, wfu, wbu
real*8 :: ve, vw, vn, vs, vf, vb, uev, uwv, wfv, wbv
real*8 :: we, ww, wn, ws, wf, wb, uew, uww, vnw, vsw
!---------------------------QUICK---------------------------!

!---------------------------AX=B----------------------------!
real*8 ,dimension(1:ndim,1:7) :: coef
integer, dimension(1:mdim):: jcoef  
real*8 ,dimension(1:ndim) :: div1,p_s,r_s,r2_s ,v_s, ss_s,t_s
real*8 ,dimension(1:ndim) :: x1
!---------------------------AX=B----------------------------!

!---------------------------LES-----------------------------!
real*8 ,parameter :: Cs = 0.18
real*8 :: nut,delta
!---------------------------LES-----------------------------!

!---------------------------Grid----------------------------!
!Initial grid coordinates for evaluating grid lengths
real*8 ,dimension (-1:nx+3) :: X
real*8 ,dimension (-1:ny+3) :: Y
real*8 ,dimension (-1:nz+3) :: Z 

!Actual grid cooridinates (with adjusted index)
real*8 ,dimension (1:nx+1) :: Xa
real*8 ,dimension (1:ny+1) :: Ya
real*8 ,dimension (1:nz+1) :: Za

!Grid lengths
real*8 ,dimension (-1:nx+2) :: iDx
real*8 ,dimension (-1:nx+2) :: Dxs
real*8 ,dimension (-1:ny+2) :: iDy
real*8 ,dimension (-1:ny+2) :: Dys
real*8 ,dimension (-1:nz+2) :: iDz
real*8 ,dimension (-1:nz+2) :: Dzs

!Midpoints of grid coordinates
real*8 ,dimension (1:nx) :: Xs
real*8 ,dimension (1:ny) :: Ys
real*8 ,dimension (1:nz) :: Zs
!---------------------------Grid----------------------------!

!---------------------------DFIB----------------------------!
real*8 ,dimension(-1:nx+2,-1:ny+2,-1:nz+2):: ETA

real*8 ,parameter :: d = 1
real*8 ,parameter :: r = d*0.5
real*8 ,parameter :: xc = lx/2
real*8 ,parameter :: yc = ly/2
real*8 ,parameter :: zc = lz/2
integer ,parameter :: nSubGrids = 200

integer ,parameter :: iBgnVOS = 1
integer ,parameter :: iEndVOS = nx
integer ,parameter :: jBgnVOS = 1
integer ,parameter :: jEndVOS = ny
integer ,parameter :: kBgnVOS = 1
integer ,parameter :: kEndVOS = nz

!---------------------------DFIB----------------------------!

!---------------------------VIV----------------------------!
real*8 ,parameter :: pi = 3.14159265358979323846
real*8 :: u_solid = 0
real*8 :: v_solid = 0
real*8 :: w_solid = 0
real*8 ,dimension(-1:nx+2,-1:ny+2,-1:nz+2) :: u2,v2,w2 
real*8 ,dimension(-1:nx+2,-1:ny+2,-1:nz+2) :: FX,FY,FZ
real*8 :: totalFX,totalFY,totalFZ
real*8 :: cDrag ,cLift

!---for unsteady---!
real*8 :: cDragSum ,cLiftSum
real*8 :: cDragSum_rms ,cLiftSum_rms
real*8 :: cDragAvg ,cLiftAvg
real*8 :: cDrag_rms ,cLift_rms
real*8 :: totalFXSum,totalFYSum
real*8 :: totalFXAvg,totalFYAvg
!---for unsteady---!

real*8 ,dimension(1:nx,1:ny) :: FXz ,FYz ,FZz
real*8 ,dimension(1:nx) :: FXy ,FYy ,FZy
real*8 :: v_solid_y = 0
real*8 :: dis_Y = 0
real*8 ,parameter :: U_RS = 5.568
real*8 ,parameter :: M_s = 0.148
real*8 ,parameter :: solid_zeta = 0.0012
!---------------------------VIV----------------------------!

!---------------------Analytical Cal---------------------!
!real*8 :: GradientP = 0
!real*8 ,dimension (1:ny) :: ProfileU
!real*8 :: DifferenceU = 0, U_RMS=0, U_L2norm=0
!real*8 :: DifferenceGradientP=0, Diff_Max_U=0
!real*8 :: tempp=0, axis_shift=0
!real*8 :: out_count=0, right_count=0, left_count=0
!---------------------Analytical Cal---------------------!

!--------------------endpoints---------------------!
real*8 :: p_in=0., p_in_count=0. , p_out=0. , p_out_count=0.
real*8 :: minADh=0., minAavgU=0.,  minARe=0., minA=242.3487e-6 , minAP=61.25e-3 
!---------------------Analytical Cal---------------------!



end module
