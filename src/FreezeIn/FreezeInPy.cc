//This C++ program exposes the functions in the FreezeIn.h library to python
//using a light-weight header-only pybind11 library (included in this
//repository)

/********************/
/* FreezeIn Library */
/********************/

#include "FreezeIn.h"

/******************/
/* Pybind Library */
/******************/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

/************************************/
/* Exposing C++ functions to python */
/************************************/

PYBIND11_MODULE(FreezeIn, mod)
{
    mod.def("Read_gstar", &Read_gstar, R"pbdoc(
    Read gstar(S) data from various .tab files in the gstar folder.
    Set the variable gstarpath to the path to the gstar folder.
    Various choices:
    "standard": Gondolo Gelmini TQCD = 150 MeV (default)
    "HP_A": Hindmarsh-Philipsen equation of state A
    "HP_B": Hindmarsh-Philipsen equation of state B
    "HP_B2": Hindmarsh-Philipsen equation of state B2
    "HP_B3": Hindmarsh-Philipsen equation of state B3
    "HP_C": Hindmarsh-Philipsen equation of state C
    )pbdoc", py::arg("choice")="standard", py::arg("gstarpath")=GSTARPATH);

    mod.def("gstar", &gstar, R"pbdoc(
    g* as a function of temperature T in GeV.
    Evaluate Read_gstar(choice) with various input choices to use various g*(S)
    )pbdoc", py::arg("T"));
    
    mod.def("gstarS", &gstarS, R"pbdoc(
    g*S as a function of temperature T in GeV.
    Evaluate Read_gstar(choice) with various input choices to use various g*(S) 
    )pbdoc", py::arg("T"));

    mod.def("dlngstarSdlnT", &dlngstarSdlnT, R"pbdoc(
    dlng*S/dlnT as a function of temperature T in GeV.
    Evaluate Read_gstar(choice) with various input choices to use various g*(S) 
    )pbdoc", py::arg("T"));

    mod.def("RhoVisible", &RhoVisible, R"pbdoc(
    Energy density in the visible sector as a function of T
    )pbdoc", py::arg("T"));

    mod.def("EntropyVisible", &EntropyVisible, R"pbdoc(
    Comoving entropy in the visible sector as a function of T
    )pbdoc", py::arg("T"));

    mod.def("HubbleVisible", &HubbleVisible, R"pbdoc(
    Hubble rate in the visible sector as a function of T
    )pbdoc", py::arg("T"));

    mod.def("YieldEqVisible", &YieldEqVisible, R"pbdoc(
    Equilibrium yield as a function of T in the visible sector
    )pbdoc", py::arg("T"), py::arg("mchi"));

    mod.def("SigmaV_chi", &SigmaV_chi, R"pbdoc(
    Thermally-averaged cross section for SM SMbar -> chi chi process
    )pbdoc", py::arg("T"), py::arg("mchi"), py::arg("kappa"),
             py::arg("TQCD")=0.15);

    mod.def("kappa_FreezeIn", &kappa_FreezeIn, R"pbdoc(
    Computes kappa to freeze-in the observed relic abundance
    )pbdoc", py::arg("mchi"), py::arg("TQCD")=0.15);
    
    mod.def("SigmaDDe", &SigmaDDe, R"pbdoc(
    Direct detection cross section in cm^2
    )pbdoc", py::arg("mchi"), py::arg("kappa"));
    
};
