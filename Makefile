# Copyright (C) 2010-2013 Jon Herman, Josh Kollat, and others.

# Hymod is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Hymod is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with Hymod.  If not, see <http://www.gnu.org/licenses/>.

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
