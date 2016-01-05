MPICXX = mpig++-4.8.4
#MPICXX  = mpixlcxx-12.1.0.13 
CXXFLAGS = -Ofast -ffast-math -ftree-vectorize # optimize code
#CXXFLAGS =  -O5 -qhot -qsimd=auto -qarch=qp -qlist
LIBS = -L /usr/gapps/tamu/bgqos_0/gperftools-2.4-gcc-4.8.4/lib -ltcmalloc_minimal
DFLAGS = -O0 # non-optimized version
DEPFILE	= .depends
DEPTOKEN	= '\# MAKEDEPENDS'
DEPFLAGS	= -Y -f $(DEPFILE) -s $(DEPTOKEN) -p 
SRCS	= *.cc
OBJS	 = $(SRCS:.cc=.o)
OBJS_O	= $(foreach obj, $(OBJS), $(obj) )

all: depend $(SRCS)
	$(MPICXX) $(CXXFLAGS) -c $(SRCS)
	$(MPICXX) $(CXXFLAGS) -o SimpleLD $(OBJS) $(LIBS)

depend:
	make $(DEPFILE)

$(DEPFILE):
	@echo $(DEPTOKEN) > $(DEPFILE)
	makedepend $(DEPFLAGS) -- $(CXXFLAGS) -- $(SRC) >&/dev/null
	
$(SRCS): 
 
clean:
	rm -f $(DEPFILE)
	rm -f SimpleLD
	
debug: depend $(SRCS)
	$(MPICXX) -o SimpleLD -g $(DFLAGS) $(SRCS)
 

# put this file in the last line of your Makefile
sinclude $(DEPFILE)