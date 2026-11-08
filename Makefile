# C sources are compiled by gcc; C++ gravity sources by g++.
CC      = gcc
CXX     = g++

CFLAGS   = -O2 -Wall -I. -Iinclude
CXXFLAGS = -O2 -Wall -std=c++17 -I. -Iinclude

# Existing CSPICE installation in the project.
CSPICE = lib/cspice.a

# SHTOOLS C++ wrapper used by grav.cpp. Override on the make command line
# if the library has another name or is located elsewhere, e.g.
# make SHTOOLS_LIBS="path/to/libshtools_wrapper.a"
SHTOOLS_LIBS ?= -lSHTOOLS -lfftw3 -llapack -lblas -lgfortran
SHTOOLS_LDFLAGS ?= -L/usr/local/lib -Llib

TARGET = messenger_orbit
C_OBJECTS = dynamics.o integrator.o main.o
CXX_OBJECTS = mercury_gravity.o grav.o
OBJECTS = $(C_OBJECTS) $(CXX_OBJECTS)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(CSPICE) $(SHTOOLS_LDFLAGS) $(SHTOOLS_LIBS) -o $@ -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGET) orbit_diff.txt

.PHONY: all clean
