PROGS = fqc

LIBS =  -lz -lhts
CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -std=c++14

all: $(PROGS)

install: $(PROGS)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin
	rm $(PROGS)

$(PROGS): OptionParser.o smithlab_utils.o

.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -r $(PROGS) *.o *.so *.a *~

.PHONY: clean
