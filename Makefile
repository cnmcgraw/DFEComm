# prompt> make
#
CC     = g++ -std=c++11     # the c compiler to use
MPICC  = mpic++  # the MPI cc compiler
CFLAGS = -O3     # optimize code
DEPFILE	= .depends
DEPTOKEN	= '\# MAKEDEPENDS'
DEPFLAGS	= -Y -f $(DEPFILE) -s $(DEPTOKEN) -p 

# SRCS	 = Problem_Input.cc Quadrature.cc Cell.cc CellSet.cc Subdomain.cc Problem.cc Sweep_Mini_App.cc
 SRCS	= *.cc
 OBJS	 = $(SRCS:.c=.o)
 OBJS_O	= $(foreach obj, $(OBJS), $(obj) )

 all: depend $(SRCS)
	$(MPICC) $(CFLAGS) $(DFLAGS) $(SRCS)
	
 depend:
	rm -f $(DEPFILE)
	make $(DEPFILE)

 $(DEPFILE):
	@echo $(DEPTOKEN) > $(DEPFILE)
	makedepend $(DEPFLAGS) -- $(CFLAGS) -- $(SRC) >&/dev/null
	
 $(SRCS): 
 
 clean:
	rm -f $(DEPFILE)
	rm a.out

 # put this file in the last line of your Makefile
 sinclude $(DEPFILE)
