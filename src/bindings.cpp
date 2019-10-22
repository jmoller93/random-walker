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
		//.def(py::init<>(),py::arg())
		.def("get_coord", &walker::get_coord)
		.def("get_length", &walker::get_lengths)
        .def("get_rg", &walker::get_rg)
        .def("get_com", &walker::get_com)
        .def("get_gene_length", &walker::get_gene_length)
        .def("set_monomers", &walker::set_monomers,
                py::arg("n") = 0)
        .def("set_bond_length", &walker::set_bond_length,
                py::arg("l") = 0)
        .def("chain_growth", &walker::chain_growth,
                py::arg("tol") = 0.0, py::arg("max_trials") = 0)
        .def("dist", &walker::dist);
		//.def("save", &distance_matrix::save)
	
}
