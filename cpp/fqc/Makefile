PROGRAMS = fqc
CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -std=c++14

# optional zlib and hts
LIBS =
ifndef NO_ZLIB
	LIBS += -lz
	CXXFLAGS += -DUSE_ZLIB
endif

ifndef NO_HTS
	LIBS += -lhts
	CXXFLAGS += -DUSE_HTS
endif

# get path to binary at compile time to look for files
CXXFLAGS += "-D__MYFILE__=\"`pwd`\""

all: $(PROGRAMS)

install: $(PROGRAMS)
	@mkdir -p bin
	@install -m 755 $(PROGRAMS) bin
	rm $(PROGRAMS)

$(PROGRAMS): $(addprefix src/, HtmlMaker.o \
	                             FastqStats.o \
	                             StreamReader.o \
	                             Config.o \
	                             OptionParser.o \
															 smithlab_utils.o)

%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c -o src/$@ $<

%: src/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -r $(PROGRAMS) src/*.o

.PHONY: clean
