 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
balance: $(objects) $(fft_lib) $(sourcedir)/post/balance.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/balance.f90 -o balance $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)

