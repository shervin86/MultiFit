##/usr/bin/g++
# compiler
CC=g++
# compilation options
CXXFLAGS=-Wall -std=c++11  -pedantic 

# addigionals for root
ROOT_LIB:=`root-config --libs --glibs`
ROOT_FLAGS:=`root-config --cflags --ldflags`
ROOT_INCLUDE:=`root-config --incdir`

# local directories
INCLUDEDIR=./include
SRCDIR=./src
BUILDDIR=./bin
OBJ_DIR=./lib

# Boost
#BOOST = /afs/cern.ch/cms/slc5_amd64_gcc434/external/boost/1.47.0
#BOOST=/afs/cern.ch/cms/slc5_amd64_gcc472/external/boost/1.50.0
SCRAMTOOL=$(shell which scram &>/dev/null || echo 1)

ifeq ($(SCRAMTOOL),1)
	ROOFIT_LIB="-lRooFit"
	ROOSTAT_LIB="-lRooStats"
	ROOFIT_INCLUDE="./"
else
	ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
	ROOFIT_LIB := -l$(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIB= | sed 's|LIB=||')
	ROOFIT_LIB += -l$(shell cd $(CMSSW_BASE); scram tool info roofit | grep LIB= | sed 's|LIB=||')
	ROOFIT_LIBDIR = -L$(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')
	ROOFIT_LIB+=$(ROOFIT_LIBDIR)
endif
ROOSTAT_LIB="-lRooStats"


#################
INCLUDE=-I$(INCLUDEDIR)  -isystem$(ROOT_INCLUDE)  -I$(ROOFIT_INCLUDE) #-I$(BOOST)/include
LIB=-L$(BOOST)/lib -L/usr/lib64 # -L/usr/lib 



#CFLAGS=-Iinterface/ $(ROOTCFLAGS) -I$(ROOTINCLUDES)  -c 
#LFLAGS=-I. -lm $(ROOTLIBS)

CSOURCES:=${wildcard src/*.cc}
ESOURCES:=${wildcard bin/*.cpp}
OBJECTS:=$(patsubst src/%.cc, lib/%.o, $(CSOURCES))
EXECS:=$(ESOURCES:.cpp=.exe)



#### Make the list of dependencies for a particular module

#------------------------------ MODULES (static libraries)

MAKEDEPEND = -MMD  -MT '$@ lib/$*.d'

# $<: first prerequisite -> put always the .cc as first 
#### General rule for classes (libraries) compiled statically
### Generate also a .d file with prerequisites
$(OBJECTS): $(OBJ_DIR)/%.o: $(SRCDIR)/%.cc
	@echo "--> Making $@"
	@$(COMPILE.cc) $(INCLUDE) $(MAKEDEPEND) $< -o $@

#this includes the statements contained in the .d files with all the dependencies
-include $(OBJECTS:.o=.d)


modules: $(OBJECTS)
all: modules $(EXECS)

###### Main program
default: test
#$(BUILDDIR)/study_Ra $(OBJECTS)

test: $(BUILDDIR)/test.exe
$(BUILDDIR)/test.exe: $(BUILDDIR)/test.cpp
	@$(CXX) $(CXXFLAGS) $(INCLUDE) $(MAKEDEPEND) -o $@ $< $(OBJECTS)  $(LIB) #$(ROOT_LIB) $(ROOFIT_LIB) $(ROOSTAT_LIB) $(ROOT_FLAGS) 


$(BUILDDIR)/study.exe:  $(BUILDDIR)/study.cpp $(MODULES)
	@echo "---> Making study $(COMPILE)"
	@g++ $(CXXFLAGS) $(INCLUDE) $(MAKEDEPEND) -o $@ $< $(MODULES)  $(LIB) $(ROOT_LIB) $(ROOFIT_LIB) $(ROOSTAT_LIB) $(ROOT_FLAGS) \
	#-lboost_program_options -lTreePlayer 

$(BUILDDIR)/study_Ra.exe:  $(BUILDDIR)/study_Ra.cpp $(OBJECTS)
	@echo "---> Making study $(COMPILE)"
	@g++ $(CXXFLAGS) $(INCLUDE) $(MAKEDEPEND) -o $@ $< $(OBJECTS)  $(LIB) $(ROOT_LIB) $(ROOFIT_LIB) $(ROOSTAT_LIB) $(ROOT_FLAGS) \
	#-lboost_program_options -lTreePlayer 

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(OBJ_DIR)/*.d
	rm -f $(BUILDDIR)/*.exe






##CFLAGS opzioni di compilazione
##LFLAGS opzioni di linking

##Primo passo: da .cc/.cpp a file oggetto
##Secondo passo: linking e creazione dell'eseguibile

.SUFFIXES:.exe
.MRPROPER:clean ##Mastro Lindo

