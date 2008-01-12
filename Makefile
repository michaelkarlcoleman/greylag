# Makefile for cgreylag module

# The SWIG Python/STL interface is pretty new, so SWIG 1.3.31 or later is
# required.
# Python 2.5 or later is required.
# A reasonably recent g++/libstdc++ may also be required.

# Developed with swig 1.3.31, g++ 3.4.6/4.1.2, libstdc++.so.6


.PHONY: all modsyms install clean tags check release
.DELETE_ON_ERROR:


DEST = /usr/local/lib/greylag/

PYTHONFLAGS = $(shell python-config --include)

CXXBASEFLAGS = -Wall -g3 -fPIC

# use to add flags from the make command-line (e.g., -march, -mtune)
CXXFASTFLAGS =

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
CXXFLAGS = $(CXXBASEFLAGS) -O3 -ffast-math -mfpmath=sse -DNDEBUG
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


# summary of C++ modules symbols used by main script
modsyms ::
	@sed -n -e 's/^.*\(cgreylag\.[a-zA-Z0-9_.]*\).*$$/\1/p' \
			greylag_grind.py \
		| sort | uniq -c

tags :: TAGS
TAGS : $(MODULE).cpp $(MODULE).hpp
	etags $^

# FIX: we could compile the .py files here
install :: all
	[ -d $(DEST) ] || install -d $(DEST)
	install -p _$(MODULE).so $(DEST)
	install -p --mode=444 $(MODULE).py $(DEST)
	install -p greylag_grind.py $(DEST)/greylag-grind
	install -p greylag_index_spectra.py $(DEST)/greylag-index-spectra
	install -p greylag_merge.py $(DEST)/greylag-merge
	install -p greylag_sqt.py $(DEST)/greylag-sqt
	install -p greylag_solo.py $(DEST)/greylag-solo
	install -p greylag_validate.py $(DEST)/greylag-validate

clean ::
	-rm -f $(MODULE).py $(MODULE)_wrap.cpp $(MODULE).o $(MODULE)_wrap.o \
		_$(MODULE).so *.py[co] TAGS *~ .??*~ \
		test/*.py[co] test/*.idx test/*.glw test/*-bad.sqt test/tmp* \
		test/*~

check :: all
	nosetests --exe --with-doctest $(NOSEFLAGS)


release :: all
	@echo "# did you update __version__ in greylag_grind.py?"
	git-tag -l "v$(VERSION)" || false # no such tag
	git-archive --format=tar --prefix=greylag-$(VERSION)/ v$(VERSION) \
	    | gzip -9 > ../greylag-$(VERSION).tgz

