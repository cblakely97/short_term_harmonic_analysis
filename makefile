make_prog:
	ifort -traceback -O0 short_term_lsq.F -I/afs/crc.nd.edu/x86_64_linux/n/netcdf/4.7.0/intel/18.0/include -L/afs/crc.nd.edu/x86_64_linux/n/netcdf/4.7.0/intel/18.0/lib -lnetcdf -lnetcdff -o short_term_lsq.o
