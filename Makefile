all: iso_server.f
	f90 $$F90FLAGS -o iso_server -s -w iso_server.f $$LINK_FNL
