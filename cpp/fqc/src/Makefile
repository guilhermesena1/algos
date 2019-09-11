ifndef SMITHLAB_CPP
	$(error Must define SMITHLAB_CPP varialbe)
endif

PROGS = fqc

INCLUDEDIRS = $(SMITHLAB_CPP)
#LIBS = 
INCLUDEARGS = $(addprefix -I, $(INCLUDEDIRS))
CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -std=c++11

all: $(PROGS)

install: $(PROGS)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin
	rm $(PROGS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, smithlab_os.o \
	smithlab_utils.o OptionParser.o)

.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	rm -r $(PROGS) *.o *.so *.a *~

.PHONY: clean
