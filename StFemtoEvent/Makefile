# Define compiler
CXX = g++

# Define flags. -D_VANILLA_ROOT_ is needed to avoid StMessMgr confusion
CFLAGS = $(shell root-config --cflags) -O2 -fPIC -Wall -pipe -std=c++11 -D_VANILLA_ROOT_ -I.
LIBS = $(shell root-config --libs)
INCS = $(shell root-config --incdir)

# Define output library
STFEMTODST = libStFemtoDst.so

# Compile all *.cxx classes in the directory
SRC = $(shell find . -name "*.cxx")

all: $(STFEMTODST)

# $(SRC:.cc=.o)
$(STFEMTODST): $(SRC:.cxx=.o) StFemtoDst_Dict.C
	$(CXX) $(CFLAGS) -shared $^ -o $(STFEMTODST) $(LIBS)

%.o: %.cxx
	$(CXX) -fPIC $(CFLAGS) -c -o $@ $<

StFemtoDst_Dict.C: $(shell find . -name "*.h" ! -name "*LinkDef*")
	rootcint -f $@ -c -D_VANILLA_ROOT_ -DROOT_CINT -D__ROOT__ -I. -I$(INCS) $^ StFemtoDstLinkDef.h

clean:
	rm -vf *.o StFemtoDst_Dict*

distclean:
	rm -vf *.o StFemtoDst_Dict* $(STFEMTODST)
