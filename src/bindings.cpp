#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include "walker.h"

namespace py = pybind11;
using namespace walkers;

PYBIND11_MODULE(walkers, m)
{
	// Walker class
	py::class_<walker>(m, "Walker")
        .def(py::init<const uint&, const uint&>())
		.def(py::init<const std::string&>())
		.def(py::init<>())
		.def("get_coord", &walker::get_coord)
		.def("get_length", &walker::get_lengths)
        .def("get_rg", &walker::get_rg)
        //.def("get_com", &walker::get_com)
        .def("get_gene_length", py::overload_cast<> (&walker::get_gene_length))
        .def("get_gene_length", py::overload_cast<const vector_t> (&walker::get_gene_length),
                py::arg("lengths"))
        .def("get_looping_histogram", &walker::get_looping_histogram,
                py::arg("tol") = 0.0)
        .def("set_monomers", &walker::set_monomers,
                py::arg("n") = 0)
        .def("set_bond_length", &walker::set_bond_length,
                py::arg("l") = 0)
        .def("chain_growth", py::overload_cast<const f_type, const uint> (&walker::chain_growth),
                py::arg("tol") = 0.0, py::arg("max_trials") = 0)
		.def("pivot", &walker::pivot,
                py::arg("tol") = 0.0)
        .def("chain_growth", py::overload_cast<const vector_t> (&walker::chain_growth),
                py::arg("dists"))
        .def("dist", &walker::dist)
		.def("save", &walker::save);
	
}
