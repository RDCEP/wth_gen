

# n_time_intervals = 5
# year_start_arr = 1950 1980 2010 2040 2070
# year_stop_arr =  1980 2010 2040 2070 2100
fname_arr = GENERIC1.WTH GENERIC2.WTH GENERIC3.WTH GENERIC4.WTH GENERIC5.WTH

in_dir = /scratch/local/isi-mip-input/wth_gen_input/HadGEM2-ES
out_dir = /scratch/local/isi-mip-input/grid/HadGEM2-ES

n_procs = 1

export LD_LIBRARY_PATH := /autonfs/home/dmcinern/lib:$(LD_LIBRARY_PATH)

nc_wth_gen: nc_wth_gen.f90
	gfortran -o nc_wth_gen nc_wth_gen.f90 -L/autonfs/home/dmcinern/lib -lnetcdf -lnetcdff -I/autonfs/home/dmcinern/include/ -mcmodel=medium

GENERIC1.WTH: nc_wth_gen
	./nc_wth_gen 1950 1980 $(in_dir) $(out_dir) $@ $(n_procs) 1 > log/GENERIC1.LOG

GENERIC2.WTH: nc_wth_gen
	./nc_wth_gen 1980 2010 $(in_dir) $(out_dir) $@ $(n_procs) 1 > log/GENERIC2.LOG

GENERIC3.WTH: nc_wth_gen
	./nc_wth_gen 2010 2040 $(in_dir) $(out_dir) $@ $(n_procs) 1 > log/GENERIC3.LOG

GENERIC4.WTH: nc_wth_gen
	./nc_wth_gen 2040 2070 $(in_dir) $(out_dir) $@ $(n_procs) 1 > log/GENERIC4.LOG

GENERIC5.WTH: nc_wth_gen
	./nc_wth_gen 2070 2099 $(in_dir) $(out_dir) $@ $(n_procs) 1 > log/GENERIC5.LOG

all: $(fname_arr)

.PHONY: $(fname_arr)
