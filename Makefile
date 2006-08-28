

# tool issues?  issues mixing g++, libstdc++, swig, python and python header versions?
# works on devel01 with swig 1.3.28, g++ 4.1.2, libstdc++.so.6 (ld 2.16.91)




.PHONY: all clean
.DELETE_ON_ERROR:


#CXXFLAGS = -Wall -g3 -O0 -D_GLIBCXX_DEBUG
CXXFLAGS = -Wall -g -O2

SWIGCXXFLAGS = -fno-strict-aliasing -fPIC -I$(PYTHON_I) $(CXXFLAGS)

# Generally, this is where the 'Python.h' corresponding to your 'python' lives
#PYTHON_I = /n/site/inst/Linux-i686/sys/include/python2.4
PYTHON_I = /usr/include/python2.4

MODULE = spectrum

all :: _$(MODULE).so test_spectrum

$(MODULE)_wrap.cpp : $(MODULE).i $(MODULE).h
	swig -c++ -python -shadow -o $@ $<

$(MODULE)_wrap.o : $(MODULE)_wrap.cpp $(MODULE).h
	g++ $(SWIGCXXFLAGS) -c $<

$(MODULE).o : $(MODULE).cpp $(MODULE).h

_$(MODULE).so : $(MODULE).o $(MODULE)_wrap.o
	g++ $(CXXFLAGS) -shared $^ -o $@


test_spectrum : test_spectrum.cpp spectrum.o
	g++ $(CXXFLAGS) $^ -o $@


clean::
	@rm -f $(MODULE)_wrap.cpp *.o *.so *.py[co] test_spectrum

