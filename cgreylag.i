// SWIG interface definition for module cgreylag

//     greylag, a collection of programs for MS/MS protein analysis
//     Copyright (C) 2006-2008  Stowers Institute for Medical Research
//
//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//     Contact: Mike Coleman
//              Stowers Institute for Medical Research
//              1000 East 50th Street
//              Kansas City, Missouri  64110
//              USA


%module cgreylag


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
%template(vector_string) std::vector<std::string>;
%template(vector_vector_int) std::vector< std::vector<int> >;
%template(vector_vector_double) std::vector< std::vector<double> >;
%template(vector_vector_match) std::vector< std::vector<match> >;


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


%include file.i


%naturalvar search_context::pca_residues;

%include "cgreylag.hpp"
