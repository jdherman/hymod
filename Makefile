# Hymod Makefile
# Jon Herman

TARGET = hymod
CC = g++
C_FLAGS = -O0 -g

SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

all: $(SOURCES) $(TARGET)

.cpp.o:
	$(CC) -c $(C_FLAGS) $^ -o $@
	
$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) $(C_FLAGS) -o $@ 

clean:
	rm -rf *.o $(TARGET)