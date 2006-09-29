# Makefile for cxtpy module

#	$Id$


# SWIG is still experiencing rapid development--1.3.28 or better is required.
# A reasonably recent g++/libstdc++ may also be required.  Python 2.4 or
# better is assumed.

# Developed (on devel01) with swig 1.3.28, g++ 4.1.2, libstdc++.so.6 (ld 2.16.91)


.PHONY: all clean
.DELETE_ON_ERROR:

# for debugging (extra checking, slow)
#CXXFLAGS = -Wall -g3 -O0 -D_GLIBCXX_DEBUG

# for speed (fastest?, fewest checks)
#CXXFLAGS = -Wall -g -O3 -DNDEBUG

# reasonably fast
CXXFLAGS = -Wall -g -O2

SWIGCXXFLAGS = $(CXXFLAGS) -fno-strict-aliasing -Wno-unused-function -fPIC \
			-I$(PYTHON_I)


# Generally, this is where the 'Python.h' corresponding to your 'python' lives.
#PYTHON_I = /n/site/inst/Linux-i686/sys/include/python2.4
PYTHON_I = /usr/include/python2.4

MODULE = cxtpy

all :: _$(MODULE).so

$(MODULE)_wrap.cpp : $(MODULE).i $(MODULE).hpp
	swig -c++ -python -o $@ $<

$(MODULE)_wrap.o : $(MODULE)_wrap.cpp $(MODULE).hpp
	g++ $(SWIGCXXFLAGS) -c $<

$(MODULE).o : $(MODULE).cpp $(MODULE).hpp

_$(MODULE).so : $(MODULE).o $(MODULE)_wrap.o
	g++ $(CXXFLAGS) -shared $^ -o $@


clean::
	@rm -f $(MODULE).py $(MODULE)_wrap.cpp $(MODULE).o $(MODULE)_wrap.o \
		_$(MODULE).so *.py[co]
