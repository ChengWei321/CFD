! Unit cube to be gridded.
!    y=1 ______________  
!       /             /|         
!      /             / |        
!     /____________ /  |        
!     |  |         |   |        
!     |  |         |   |        
!     |  | x=y=z=0 |   |        
!     |  |_________|___|x=1     
!     |  /         |  /
!     | /          | /
!     |/___________|/
!    z=1         
!

program main
   use variables
   use mpi
   use omp_lib
   implicit none

   !------------------- OPENMP ------------------------!
      !nthreads = 10
      !call omp_set_num_threads(nthreads)
   !------------------- OPENMP ------------------------!

   !------------------------ MPI ------------------------!
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call Mpi_division()
   !------------------------ MPI ------------------------!

   !--------for nz--------!
   istart = gstart(myid)  !
   iend = gend0(myid)     !
   igcount = gcount(myid) !
   !--------for nz--------!


   !-----------------Parameters for the simulation------------------!
      omega=1.8          ! Set value for SOR method
      zeta=1.e-5         ! zeta and itmax for solving pressure matrix
      zeta_vel=1.e-4     ! convergence condition for velocity field
      time = 0.          ! initialize time of simulation
      nstep = 200000       ! number of timesteps for the simulation
      isto = 1000          ! data stored every 'isto' steps
      Q1 =  1            ! the LES mode. 1 : on ; 0 : off
      solidPos = 1       ! fixed : 1 ; moving : 2 ; no solid : 0
      steadiness = 1     ! steady : 1 ; unsteady : 2 
      StartAvg = 500000  ! determine the time of beginning calculation of AvgFlow. if steadiness = 2 
      totalcosttime = 0  ! initialize wall time
!      inputfile = 'finaltimestep.q'
   !-----------------Parameters for the simulation------------------!


!
   call reading_data()           
   call gridder_equal()
   call initial_conditions()  
!   call reading_variables()  
   call read_ETA()
   call boundary_conditions()  
!   if(solidPos==1 .OR. solidPos==2)then; call DFIB_Cylinder(); end if

   if(myid==master)then

      call filereachtime()     
      call filerInfo()

   end if

   

   !----------------------------for sendrecv----------------------------!
   l_nbr = myid - 1                                                     !
   r_nbr = myid + 1                                                     !         
   if(myid == 0) then; l_nbr=MPI_PROC_NULL; end if                      !
   if(myid == (nproc-1)) then; r_nbr=MPI_PROC_NULL; end if              !
   !----------------------------for sendrecv----------------------------!


   
   write(*,*) myid, 'istart = ', istart, 'iend = ', iend, 'gcount = ', igcount



!-------------------------main loop on the timesteps----------------------!
   do istep=1,nstep  
      totalstarttime = MPI_WTIME()

      
      


!      if(solidPos == 2) then; call DFIB_Cylinder(); end if


      
      !----------data transformation among nodes----------!
      icount = 2*(nx+4)*(ny+4)
      itag = 110
      call MPI_SENDRECV( u(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         u(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 120
      call MPI_SENDRECV( v(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         v(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 130
      call MPI_SENDRECV( w(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         w(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )


      itag = 140
      call MPI_SENDRECV( u(-1,-1,iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         u(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 150
      call MPI_SENDRECV( v(-1,-1,iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         v(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 160
      call MPI_SENDRECV( w(-1,-1,iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         w(-1,-1,istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )



      icount = 2
      itag = 170
      call MPI_SENDRECV( iDz(istart), icount, MPI_REAL8, l_nbr, itag, &
                         iDz(iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 180
      call MPI_SENDRECV( Dzs(istart), icount, MPI_REAL8, l_nbr, itag, &
                         Dzs(iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )


      itag = 190
      call MPI_SENDRECV( iDz(iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         iDz(istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      itag = 200
      call MPI_SENDRECV( Dzs(iend-1), icount, MPI_REAL8, r_nbr, itag, &
                         Dzs(istart-2), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
                         
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!



      


      !------------------------------------------------------------------------------------------------------!
      call discretisation_QUICK_centre()
      !------------------------------------------------------------------------------------------------------!



      
      
      !----------data transformation among nodes----------!
      icount = (nx+4)*(ny+4)
      itag = 240
      call MPI_SENDRECV( w_star(-1,-1,iend), icount, MPI_REAL8, r_nbr, itag, &
                         w_star(-1,-1,istart-1), icount, MPI_REAL8, l_nbr, itag, MPI_COMM_WORLD, status, ierr )
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!


      !------------------------------------------------------------------------------------------------------!
      call gauss_seidel()       ! calculate pressure field
      !------------------------------------------------------------------------------------------------------!


      !----------data transformation among nodes----------!
      icount = (nx+4)*(ny+4)

      itag = 270
      call MPI_SENDRECV( p(-1,-1,istart), icount, MPI_REAL8, l_nbr, itag, &
                         p(-1,-1,iend+1), icount, MPI_REAL8, r_nbr, itag, MPI_COMM_WORLD, status, ierr )
                    
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data transformation among nodes----------!


      !------------------------------------------------------------------------------------------------------!
      call calcul_new_velocity() ! update velocity field 
      !------------------------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------------------------!
      call Updating_velocity()   ! update velocity field 
      !------------------------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------------------------!
      call boundary_conditions() ! recall boundary conditions to update them
      !------------------------------------------------------------------------------------------------------!
      









      time = time + dt
      !----------data collect among nodes for filer----------!
      icount = igcount*(nx)*(ny)
      !Send my results back to the master
      if(myid>master)then
         itag = 10
         call MPI_SEND( pre(1,1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 20
         call MPI_SEND( uc(1,1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 30
         call MPI_SEND( vc(1,1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
         itag = 40
         call MPI_SEND( wc(1,1,istart), icount, MPI_REAL8, master, itag, MPI_COMM_WORLD, ierr )
      end if
      !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      !Wait to receive results from each task
      if(myid==master)then
         do i = 1, (nproc-1)
            itag = 10
            call MPI_RECV( pre(1,1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 20
            call MPI_RECV( uc(1,1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 30
            call MPI_RECV( vc(1,1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
            itag = 40
            call MPI_RECV( wc(1,1,gstart(i)), icount, MPI_REAL8, i, itag, MPI_COMM_WORLD, status, ierr )
         end do
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !----------data collect among nodes for filer----------!



      if(myid==master)then

         write(*,*) 'time = ',time  
         call filerProcess()
               
      end if
      

       



      VelocityDifference = 1024
!      if (mod(istep,isto)==0) then ! write results if isto = istep
               
         if(myid==master)then

            call check_steady()
            
         end if 

!      end if


      if (mod(istep,isto*5)==0) then ! write results if isto = istep
               
         if(myid==master)then

            call filereachtime()

         end if 

      end if



      !----------calculate wall time----------!
      totalfinaltime = MPI_WTIME()
      totalcosttime = totalcosttime + (totalfinaltime-totalstarttime)
      if(myid==master)then; write(*,*) 'each cost = ', totalfinaltime-totalstarttime; end if
      if(myid==master)then; write(*,*) 'total cost = ', totalcosttime; end if
      !----------calculate wall time----------!


      !--------------exit time loop--------------!
      if ( VelocityDifference .lt. zeta_vel ) then
         exit 
      end if
      !--------------exit time loop--------------!
      


   end do
   istep = istep - 1
!-------------------------main loop on the timesteps----------------------!

   


   if(myid==master)then

      call endpoints()
      call filer_final()

   end if

   call MPI_FINALIZE(ierr) 

end program main
