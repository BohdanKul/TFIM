CC=g++
CFLAGS=-c -o3 -std=c++0x -fPIC
LDFLAGS=-L. -lboost_program_options -lboost_filesystem -lboost_system
SOURCES=main.cpp tfim.cpp communicator.cpp bond.cpp spin.cpp directedloop.cpp vertex.cpp helper.cpp randombase.cpp estimator.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main.e

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLE)

