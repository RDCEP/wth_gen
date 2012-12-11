
# export LD_LIBRARY_PATH := /autonfs/home/dmcinern/lib:$(LD_LIBRARY_PATH)

nc_wth_gen: nc_wth_gen.f90
#	gfortran -o nc_wth_gen nc_wth_gen.f90 -L/autonfs/home/dmcinern/lib -lnetcdf -lnetcdff -I/autonfs/home/dmcinern/include/ -mcmodel=medium
	gfortran -O0 -g -fbounds-check -fbacktrace -o nc_wth_gen nc_wth_gen.f90 -I$(NETCDF_DIR)/include -L$(NETCDF_DIR)/lib -lnetcdf -lnetcdff -mcmodel=medium

