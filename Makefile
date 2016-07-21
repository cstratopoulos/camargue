TARGET = PSEP_b
INCLUDES = -I /Users/christos/Applications/IBM/ILOG/CPLEX_Studio1261/cplex/include/ilcplex/
LIBS = /Users/christos/Applications/IBM/ILOG/CPLEX_Studio1261/cplex/lib/x86-64_osx/static_pic/libcplex.a
CONCORDE_LIBS = /Users/christos/Dropbox/school/research/programs/concorde/concorde.a
CC = g++
CFLAGS =  -g -Wall -pedantic -Wno-long-long -ggdb -std=c++11
#-O3
#-O3
.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) $(CONCORDE_LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)
