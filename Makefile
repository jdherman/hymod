# Hymod Makefile
# Jon Herman

TARGET = hymod
CC = g++
C_FLAGS = -O3

SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

all: $(SOURCES) $(TARGET)

.cpp.o:
	$(CC) -c $(CFLAGS) -DHYMOD $^ -o $@
	
$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) $(CFLAGS) -o $@ 

clean:
	rm -rf *.o $(TARGET)