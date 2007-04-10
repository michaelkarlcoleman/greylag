# Makefile for cgreylag module

# SWIG is still experiencing rapid development--1.3.31 or later is required.
# Python 2.5 or later is required.
# A reasonably recent g++/libstdc++ may also be required.

# Developed with swig 1.3.31, g++ 3.4.6/4.1.2, libstdc++.so.6


.PHONY: all pycheck modsyms install clean tags check
.DELETE_ON_ERROR:


# MARCH = pentium3
# MARCH = pentium4
MARCH = prescott
# MARCH = opteron
# MARCH = nocona


DEST = /n/site/inst/Linux-i686/bioinfo/greylag/

PYTHONVER=2.5
PYTHONFLAGS = $(shell python$(PYTHONVER)-config --include)

CXXBASEFLAGS = -Wall -g3 -march=$(MARCH) -fPIC

# This makes it easy to compile different versions without editing this file.
# 0=debug, 2=fast, 3=faster and less safe
SPEED := 2

ifeq ($(SPEED),0)
# for debugging (extra checking, slow)
CXXFLAGS = $(CXXBASEFLAGS) -O0 -D_GLIBCXX_DEBUG
else
  ifeq ($(SPEED),2)
# reasonably fast
CXXFLAGS = $(CXXBASEFLAGS) -O2 -ffast-math -mfpmath=sse
  else
    ifeq ($(SPEED),3)
# for speed (fastest?, fewest checks)
CXXFLAGS = $(CXXBASEFLAGS) -O3 -DNDEBUG -ffast-math -mfpmath=sse
#CXXFASTFLAGS = -finline-limit=20000 --param inline-unit-growth=1000 --param large-function-growth=1000
    else
CXXFLAGS = ---INVALID-SPEED
    endif
  endif
endif

SWIGCXXFLAGS = $(CXXFLAGS) $(PYTHONFLAGS) -fno-strict-aliasing \
		-Wno-unused-function -Wno-uninitialized

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
	PYTHONVER=$(PYTHONVER) pychecker --limit 1000 greylag-grind.py

# summary C++ modules symbols used by main script
modsyms::
	@sed -n -e 's/^.*\(cgreylag\.[a-zA-Z0-9_.]*\).*$$/\1/p' greylag-grind.py \
		| sort | uniq -c

tags :: TAGS
TAGS : $(MODULE).cpp $(MODULE).hpp
	etags --members $^

# FIX: we could compile the .py files here
install::
	[ -d $(DEST) ] || install -d $(DEST)
	install -p _$(MODULE).so $(DEST)
	install -p --mode=444 $(MODULE).py $(DEST)
	install -p greylag-grind.py $(DEST)
	install -p greylag-mp.py $(DEST)

clean::
	-rm -f $(MODULE).py $(MODULE)_wrap.cpp $(MODULE).o $(MODULE)_wrap.o \
		_$(MODULE).so *.py[co] test/*.py[co] TAGS *~ .??*~ test/*~ \
		test/tmp*

check::
	nosetests --exe --with-doctest $(NOSEFLAGS)
