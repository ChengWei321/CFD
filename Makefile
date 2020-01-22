# Makefile 
# Nom du compilateur
FC = mpiifort 

# Options de compilation: optimisation, debug etc...
OPT = -O3 -ipo  -xCORE-AVX2 -align array32byte -mcmodel=large

# -ipo -xCORE-AVX2 -align array32byte
# -mcmodel=large


# nom de l'executable
EXE = Flow_MPI
# Options de l'edition de lien..
LINKOPT = 

DIR_OBJ = ./obj
DIR_BIN = ./bin

# Defining the objects (OBJS) variables
OBJS =  \
    variables_module.o \
    main.o \
    reading_data.o \
    gridder.o \
    initial_conditions.o \
    boundary_conditions.o \
    DFIB.o \
    discretisation_QUICK_centre.o \
    calcul_new_velocity.o \
    check_steady.o \
    filer.o \
    VIV.o \
    Mpi_division.o \
    GS.o\
    endpoints.o
    
# Linking object files
exe :   $(OBJS)
	$(FC) $(OPT)   -o $(EXE) \
    $(OBJS) \
    $(OPT)

    
	@echo "       !-------------------------------------------------------------------------------------------------------!" 
	#      !    y=1 ______________                                                                                 ! 
	#      !       /             /|       |  Author:  Zi-Hsuan Wei                                                 ! 
	#      !      /             / |       |  Version: 1.5                                                          ! 
	#      !     /____________ /  |       |  Web:     http://smetana.me.ntust.edu.tw/                              ! 
	#      !     |  |         |   |                                                      ______________________    ! 
	#      !     |  |         |   |                       /*****\                       |                      |   ! 
	#      !     |  | x=y=z=0 |   |                      |  o o  |                      |    ____      ____    |   ! 
	#      !     |  |_________|___|x=1                   |  :_/  |                      |   |    |    |    |   |   ! 
	#      !     |  /         |  /                      //      | \                     |   |____|____|____|   |   ! 
	#      !     | /          | /                       (|  * * | )                     |      __|    |__      |   ! 
	#      !     |/___________|/                        /`\_ * _/`\                     |     |          |     |   ! 
	#      !    z=1                                     \___)=(___/                     |     |   ____   |     |   ! 
	#      !                                                                            |     |__|    |__|     |   ! 
	#      !                                                                            |                      |   ! 
	#      !                      ___   ___   ___   .       _     __                    |                      |   ! 
	#      !                     |     |___  |   \  |      /_\   |__\                   |______________________|   ! 
	#      !                     |___  |     |___/  |___  /   \  |__/                                              ! 
	#      !-------------------------------------------------------------------------------------------------------! 


%.o:%.f90
	$(FC) $(OPT) -c $<

main.o : variables_module.o

reading_data.o : variables_module.o

gridder.o : variables_module.o

initial_conditions.o : variables_module.o

boundary_conditions.o : variables_module.o

DFIB.o : variables_module.o

discretisation_upwind.o : variables_module.o

discretisation_QUICK_centre.o : variables_module.o

discretisation_QUICK.o : variables_module.o

calcul_new_velocity.o : variables_module.o

check_steady.o : variables_module.o

filer.o : variables_module.o

VIV.o : variables_module.o

Mpi_division.o : variables_module.o

GS.o : variables_module.o

endpoints.o : variables_module.o

# Removing object files
clean :
	/bin/rm -f $(OBJS) $(EXE)  *.mod

cleanall : 
	/bin/rm -f $(OBJS) $(EXE)  *.mod
	/bin/rm -f *.dat
	/bin/rm -f *.x
	/bin/rm -f *.q
    
config :
	if [ ! -d obj ] ; then mkdir obj ; fi
	if [ ! -d run ] ; then mkdir bin ; fi
