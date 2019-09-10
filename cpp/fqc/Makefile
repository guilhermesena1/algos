## name of C++ compiler
CXX = g++

## which C++ standard and warning level
CXXFLAGS = -std=c++11 -Wall

## optimization
CXXFLAGS += -O3

## targets that may be installed as programs
PROGRAMS = fqc

## compilation pre-processor flags
#CPPFLAGS =

all: $(PROGRAMS)

%: src/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(CPPFLAGS)
