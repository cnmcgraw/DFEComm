# prompt> make
#
CC     = g++ -std=c++11 -g     # the c compiler to use
MPICC  = mpig++-4.8.4 # the MPI cc compiler
CFLAGS = -O3 # optimize code
DFLAGS = -O0 # non-optimized version
DEPFILE	= .depends
DEPTOKEN	= '\# MAKEDEPENDS'
DEPFLAGS	= -Y -f $(DEPFILE) -s $(DEPTOKEN) -p 

# SRCS	 = Problem_Input.cc Quadrature.cc Cell.cc CellSet.cc Subdomain.cc Problem.cc ParallelComm.cc Task.cc Sweep_Mini_App.cc
 SRCS	= *.cc
 OBJS	 = $(SRCS:.cc=.o)
 OBJS_O	= $(foreach obj, $(OBJS), $(obj) )

 all: depend $(SRCS)
	$(MPICC) -c  $(CFLAGS) $(SRCS)
	$(MPICC) -o SimpleLD -L /usr/gapps/tamu/bgqos_0/gperftools-2.4-gcc-4.8.4/lib -ltcmalloc_minimal $(OBJS)
 depend:
	make $(DEPFILE)

 $(DEPFILE):
	@echo $(DEPTOKEN) > $(DEPFILE)
	makedepend $(DEPFLAGS) -- $(CFLAGS) -- $(SRC) >&/dev/null
	
 $(SRCS): 
 
 clean:
	rm -f $(DEPFILE)
	rm SimpleLD
	
 debug: depend $(SRCS)
	$(MPICC) -o SimpleLD -g $(DFLAGS) $(SRCS)
 

 # put this file in the last line of your Makefile
 sinclude $(DEPFILE)
