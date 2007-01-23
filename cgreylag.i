// SWIG interface definition for module cgreylag

//	$Id$

//     Copyright (C) 2006-2007, Stowers Institute for Medical Research
//
//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License along
//     with this program; if not, write to the Free Software Foundation, Inc.,
//     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


%module cgreylag


// Declare this read-only, to suppress a warning about a possible memory leak.
%immutable mass_trace_item::description;


%feature("autodoc");

%{
#include "cgreylag.hpp"
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

// not sure these are useful--just avoiding SWIG warning
%rename(ion_increment) operator++(ion &);
%rename(ion_post_increment) operator++(ion &, int);


%include std_string.i
%apply const std::string & { std::string *name };

// %include std_list.i
%include std_vector.i

// Note: SWIG currently only exposes the outermost vector (in a vector of
// vector of X) as a modifiable object.  Inner vectors appear as tuples, and
// are thus unmodifiable.  They must therefore be assigned all at once.  This
// shortcoming will probably be fixed in a future version of SWIG.

%template(vector_int) std::vector<int>;
%template(vector_double) std::vector<double>;
%template(vector_sequence_run) std::vector<sequence_run>;
%template(vector_peak) std::vector<peak>;
%template(vector_spectrum) std::vector<spectrum>;
%template(vector_match) std::vector<match>;
%template(vector_mass_regime_parameters) std::vector<mass_regime_parameters>;
%template(vector_mass_trace_item) std::vector<mass_trace_item>;
%template(vector_vector_int) std::vector< std::vector<int> >;
%template(vector_vector_double) std::vector< std::vector<double> >;
%template(vector_vector_match) std::vector< std::vector<match> >;
%template(vector_vector_vector_double) std::vector< std::vector< std::vector<double> > >;


%include std_pair.i
%include std_map.i
%include std_multimap.i

%template(pair_double_vector_size_type) 
    std::pair<double, std::vector<spectrum>::size_type>;

// These two lines seem to be required to make the multimap template work.
// The multimap is only exposed for debugging purposes, anyway, though.
%template() std::pair<swig::PyObject_ptr, swig::PyObject_ptr>;
%template(pymap) std::map<swig::PyObject_ptr, swig::PyObject_ptr>;

%template(multimap_double_vector_size_type)
    std::multimap<double, std::vector<spectrum>::size_type>;

// currently unused
//%template(pair_char_int) std::pair<char, int>;
//%template(pair_char_double) std::pair<char, double>;

//%template(map_char_int) std::map<char, int>;
//%template(map_char_double) std::map<char, double>;


%include "typemaps.i"
%apply int *OUTPUT { int *peak_count };


%include file.i


%naturalvar search_context::pca_residues;

%include "cgreylag.hpp"
