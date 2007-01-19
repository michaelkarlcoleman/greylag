# Makefile for cgreylag module

#	$Id$


# SWIG is still experiencing rapid development--1.3.28 or better is required.
# A reasonably recent g++/libstdc++ may also be required.  Python 2.4 or
# better is assumed.

# Developed with swig 1.3.28, g++ 4.1.2, libstdc++.so.6 (ld 2.16.91)


.PHONY: all pycheck modsyms install clean tags
.DELETE_ON_ERROR:


MARCH = pentium3
# MARCH = pentium4
# MARCH = prescott
# MARCH = opteron
# MARCH = nocona


DEST = /n/site/inst/Linux-i686/bioinfo/greylag/

PYTHONVER=2.4

# Generally, this is where the 'Python.h' corresponding to your 'python' lives.
#PYTHON_I = /n/site/inst/Linux-i686/sys/include/python$(PYTHONVER)
PYTHON_I = /usr/include/python$(PYTHONVER)


# for debugging (extra checking, slow)
#CXXFLAGS = -Wall -g3 -O0 -D_GLIBCXX_DEBUG -march=$(MARCH)

# for speed (fastest?, fewest checks)
#CXXFLAGS = -Wall -g3 -O3 -DNDEBUG -ffast-math -mfpmath=sse -march=$(MARCH)
#CXXFASTFLAGS = -finline-limit=20000 --param inline-unit-growth=1000 --param large-function-growth=1000

# reasonably fast
CXXFLAGS = -Wall -g3 -O2 -ffast-math -mfpmath=sse -march=$(MARCH)

SWIGCXXFLAGS = $(CXXFLAGS) -fPIC -I$(PYTHON_I) -fno-strict-aliasing \
		-Wno-unused-function -Wno-uninitialized

PROGRAM = greylag
MODULE = cgreylag

all :: _$(MODULE).so

$(MODULE)_wrap.cpp : $(MODULE).i $(MODULE).hpp
	swig -c++ -python -o $@ $<

$(MODULE)_wrap.o : $(MODULE)_wrap.cpp $(MODULE).hpp
	g++ $(SWIGCXXFLAGS) -c $<

$(MODULE).o : $(MODULE).cpp $(MODULE).hpp
	g++ $(CXXFLAGS) $(CXXFASTFLAGS) -c $<

_$(MODULE).so : $(MODULE).o $(MODULE)_wrap.o
	g++ $(CXXFLAGS) $(CXXFASTFLAGS) -shared $^ -o $@


pycheck::
	PYTHONVER=$(PYTHONVER) pychecker --limit 1000 $(PROGRAM).py

# summary C++ modules symbols used by main script
modsyms::
	@sed -n -e 's/^.*\(cgreylag\.[a-zA-Z0-9_.]*\).*$$/\1/p' $(PROGRAM).py \
		| sort | uniq -c

tags :: TAGS
TAGS : $(MODULE).cpp $(MODULE).hpp
	etags --members $^

# FIX: we could compile the .py files here
install::
	[ -d $(DEST) ] || install -d $(DEST)
	install -p _$(MODULE).so $(DEST)
	install -p --mode=444 $(MODULE).py $(DEST)
	install -p $(PROGRAM).py $(DEST)

clean::
	-rm -f $(MODULE).py $(MODULE)_wrap.cpp $(MODULE).o $(MODULE)_wrap.o \
		_$(MODULE).so *.py[co] TAGS
