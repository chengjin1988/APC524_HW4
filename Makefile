FC=ifort
FFLAGS=-O3 -warn all
MPIFC=mpif90

all: heat_serial heat_omp heat_mpi
heat_serial: heat_serial.f90
	$(FC) $(FFLAGS) heat_serial.f90 -o heat_serial

heat_omp: heat_omp.f90
	$(FC) -openmp $(FFLAGS) heat_omp.f90 -o heat_omp

heat_mpi: heat_mpi.f90
	$(MPIFC) $(FFLAGS) heat_mpi.f90 -o heat_mpi

