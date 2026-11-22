.DEFAULT_GOAL := all

CXX ?= g++

TARGET := messenger_pipeline
BUILD_DIR := build

CSPICE_LIB ?= lib/cspice.a
SOFA_LIB   ?= lib/libsofa_c.a

SHTOOLS_LIBS ?= -lSHTOOLS -lfftw3 -llapack -lblas -lgfortran
SHTOOLS_LDFLAGS ?= -Llib

QD_LIBS ?= -lqd

CPPFLAGS += -I. -Iinclude -Iinclude/orbit -Iinclude/light-time
CXXFLAGS += -std=c++17 -O2 -Wall -Wextra

ORBIT_OBJECTS := \
    build/orbit/dynamics.o \
    build/orbit/integrator.o \
    build/orbit/mercury_gravity.o \
    build/orbit/grav.o \
    build/orbit/orbit_runner.o

LIGHTTIME_OBJECTS := \
    build/lighttime/cubic_spline.o \
    build/lighttime/eop.o \
    build/lighttime/light_time.o \
    build/lighttime/lt_runner.o \
    build/lighttime/spacecraft_ephemeris.o \
    build/lighttime/spice.o \
    build/lighttime/stations.o \
    build/lighttime/time_scales.o

MAIN_OBJECT := $(BUILD_DIR)/main.o

.PHONY: all clean orbit-build lighttime-build run

all: $(TARGET)

orbit-build:
	$(MAKE) -f orbit/Makefile.orbit

lighttime-build:
	$(MAKE) -f light-time/Makefile.lighttime

$(MAIN_OBJECT): main.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(TARGET): orbit-build lighttime-build $(MAIN_OBJECT)
	$(CXX) -o $@ \
		$(MAIN_OBJECT) \
		$(ORBIT_OBJECTS) \
		$(LIGHTTIME_OBJECTS) \
		$(CSPICE_LIB) \
		$(SOFA_LIB) \
		$(SHTOOLS_LDFLAGS) $(SHTOOLS_LIBS) \
		$(QD_LIBS) \
		-lm

run: all
	mkdir -p results
	./$(TARGET) $(ARGS)

clean:
	$(MAKE) -f orbit/Makefile.orbit clean
	$(MAKE) -f light-time/Makefile.lighttime clean
	rm -rf $(BUILD_DIR) $(TARGET)