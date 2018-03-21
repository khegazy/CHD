BASEDIR=/reg/neh/home/khegazy/baseScripts
CXX=g++
CXXFLAGS=-std=c++11 -Wall -Wno-unused-local-typedefs $(shell root-config --cflags) -I/reg/g/psdm/sw/external/fftw/3.3.4/x86_64-rhel7-gcc48-opt/include -I/reg/g/psdm/sw/external/opencv2/2.4.11/x86_64-rhel7-gcc48-opt/include -I/reg/g/psdm/sw/external/boost/1.55.0-python2.7/x86_64-rhel7-gcc48-opt/include -I/reg/neh/home/khegazy/packages/Eigen/include
LINKLIBS=$(shell root-config --libs) -L/reg/g/psdm/sw/external/fftw/3.3.4/x86_64-rhel7-gcc48-opt/lib -lfftw3 -L/reg/g/psdm/sw/external/opencv2/2.4.11/x86_64-rhel7-gcc48-opt/lib -lopencv_core -lopencv_highgui -lopencv_contrib -lopencv_imgproc -lopencv_objdetect -L/reg/g/psdm/sw/external/boost/1.55.0-python2.7/x86_64-rhel7-gcc48-opt/lib -lboost_system -lboost_math_c99 -lboost_math_tr1 #-L${ROOTSYS}/lib -lHist -lMatrix -lTree -lNet -lProofPlayer -lProof -lRIO -lThread -lMathCore -lPhysics
#LINKLIBS=$(shell root-config --libs) -L/reg/g/psdm/sw/external/fftw/3.3.4/x86_64-rhel7-gcc48-opt/lib -lfftw3 -L/reg/g/psdm/sw/external/opencv2/2.4.11/x86_64-rhel7-gcc48-opt/lib -lopencv_core -lopencv_highgui -lopencv_contrib -lopencv_imgproc -lopencv_objdetect -L/reg/g/psdm/sw/external/boost/1.55.0-python2.7/x86_64-rhel7-gcc48-opt/lib #-L${ROOTSYS}/lib -lHist -lMatrix -lTree -lNet -lProofPlayer -lProof -lRIO -lThread -lMathCore -lPhysics
#CFLAGS=-I/reg/g/psdm/sw/external/root/5.30.06-python2.7/x86_64-rhel5-gcc41-opt/include/root  -I/reg/g/psdm/sw/external/fftw/3.3.4/x86_64-rhel5-gcc44-opt/include
#LINKLIBS=-L/reg/g/psdm/sw/external/root/5.30.06-python2.7/x86_64-rhel5-gcc41-opt/lib -L/reg/g/psdm/sw/external/fftw/3.3.4/x86_64-rhel5-gcc44-opt/lib -lfftw3 #-L${ROOTSYS}/lib -lHist -lMatrix -lTree -lNet -lProofPlayer -lProof -lRIO -lThread -lMathCore -lPhysics

DEP_HPP_FILES := $(wildcard $(BASEDIR)/*.h) preProcFuncts.h
DEP_CPP_FILES := $(wildcard $(BASEDIR)/*.cpp) preProcFuncts.cpp
DEPS := $(DEP_HPP_FILES) $(DEP_CPP_FILES)
#CPP_FILES := $(filter-out plotFunctions.cpp,$(ACPP_FILES))
CPP_FILES := preProcessing.cpp#$(wildcard *.cpp)
OBJ_FILES := $(CPP_FILES:%.cpp=%.o)
DEP_OBJ_FILES := $(DEP_CPP_FILES:%.cpp=%.o)
EXE_FILES := $(CPP_FILES:%.cpp=%.exe)
STD := -std=c++0x

all: $(EXE_FILES)

.PRECIOUS: %.o 
%.o: %.cpp $(DEPS) 
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.exe: %.o $(DEP_OBJ_FILES)
	$(CXX) -o $@ $^ $(LINKLIBS)

clean:
	rm -f *.o *.exe
