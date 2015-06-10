ROOT_LIBS=`root-config --libs --glibs` -lRooFit -lRooFitCore
ROOT_INCS=`root-config --cflags`

INCDIR=./include/
SRCDIR=./src/
LIBDIR=./libs/
BINDIR=./bin/
VPATH=$(SRCDIR):$(INCDIR):$(LIBDIR)
CC = g++ $(CFLAGS) $(ROOT_INCS) -I$(INCDIR) -I$(LIBDIR)

all: bin/CutEfficiencies.exe bin/SimpleCuts.exe

clean: 
	rm bin/*.exe
##	rm libs/*

bin/%.exe: src/%.cpp 
	$(CC) -o $@ $(ROOT_LIBS) $(ROOT_INCS) $^

#### GENERAL RULE TO BUILD LIBRARIES ######
#added $(INCDIR)* so that make recompiles if there were changes in INCDIR
libs/%.o: %.cpp %.h $(INCDIR)*
	$(CC) -c $< -o $@
