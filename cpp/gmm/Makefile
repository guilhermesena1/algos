CC=g++
CXXFLAGS=-lgsl -lgslcblas -O3 -fopenmp
INCLUDEARGS=/home/cmb-panasas2/desenabr/software/anaconda/include


%.o: %.cpp
	$(CC) $(CXXFLAGS) $(OPT) -c -o $@ $< -I$(INCLUDEARGS)

all: gmm.o
	$(CC) -o gmm gmm.o $(CXXFLAGS) $(O3FLAGS)
