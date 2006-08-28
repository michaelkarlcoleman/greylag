



%module spectrum

%feature("autodoc");

%{
#include "spectrum.h"
%}

// catch STL exceptions, etc
%include "exception.i"
%exception {
    try {
        $action
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (...) {
	PyErr_SetString(PyExc_RuntimeError,"unknown C++ (swig) exception");
	return NULL;
    }
}


%include std_string.i
%apply const std::string & { std::string *name };

%include std_list.i
%include std_vector.i

%template(vector_double) std::vector<double>;
%template(vector_peak) std::vector<peak>;
%template(vector_vector_double) std::vector< std::vector<double> >;

%include file.i


%include "spectrum.h"
