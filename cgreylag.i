



%module cxtpy

%feature("autodoc");

%{
#include "cxtpy.h"
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

%template(vector_int) std::vector<int>;
%template(vector_double) std::vector<double>;
%template(vector_peak) std::vector<peak>;
%template(vector_spectrum) std::vector<spectrum>;
%template(vector_match) std::vector<match>;
%template(vector_vector_int) std::vector< std::vector<int> >;
%template(vector_vector_double) std::vector< std::vector<double> >;
%template(vector_vector_match) std::vector< std::vector<match> >;


%include std_pair.i
%include std_map.i
%include std_multimap.i

%template(pair_double_vector_size_type) 
    std::pair<double, std::vector<spectrum>::size_type>;
%template(pair_char_int) std::pair<char, int>;
%template(pair_char_double) std::pair<char, double>;

// These two lines seem to be required to make the multimap template work.
// The multimap is only exposed for debugging purposes, anyway, though.
%template() std::pair<swig::PyObject_ptr, swig::PyObject_ptr>;
%template(pymap) std::map<swig::PyObject_ptr, swig::PyObject_ptr>;

%template(multimap_double_vector_size_type)
    std::multimap<double, std::vector<spectrum>::size_type>;

%template(map_char_int) std::map<char, int>;
%template(map_char_double) std::map<char, double>;


%include "typemaps.i"
%apply int *OUTPUT { int *peak_count };


%include file.i


%include "cxtpy.h"
