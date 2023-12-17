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
    //Read_gstar(choice, gstarpath)
    mod.def("Read_gstar", &Read_gstar, R"pbdoc(
    Read tabulated data for effective number of degrees of freedom from various
    .tab files in the gstar folder with three columns:
    {Temperature in GeV, g*S, g*}.

    Inputs
    ------

    choices: "standard": Gondolo-Gelmini (LambdaQCD = 150 MeV) (default)
             "HP_A": Hindmarsh-Philipsen equation of state A
             "HP_B": Hindmarsh-Philipsen equation of state B
             "HP_B2": Hindmarsh-Philipsen equation of state B2
             "HP_B3": Hindmarsh-Philipsen equation of state B3
             "HP_C": Hindmarsh-Philipsen equation of state C
             (These are taken directly from MicrOMEGAs package)
    
    gstarpath: By default, set to the path to the gstar folder provided with
               this package
    )pbdoc", py::arg("choice")="standard", py::arg("gstarpath")=GSTARPATH);
    
    //gstar(T)
    mod.def("gstar", &gstar, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    The effective number of degrees of freedom for energy density g*

    (By default uses the standard Gondolo-Gelmini g*(T). To use other choices
    for g*: evaluate Read_gstar(choice); see documentation for the function
    Read_gstar for more details.)
    )pbdoc", py::arg("T"));
    
    //gstarS(T)
    mod.def("gstarS", &gstarS, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    The effective number of degrees of freedom for entropy density g*S

    (By default uses the standard Gondolo-Gelmini g*S(T). To use other choices
    for g*S: evaluate Read_gstar(choice); see documentation for the function
    Read_gstar for more details.)
    )pbdoc", py::arg("T"));
    
    //dlngstarSdlnT(T)
    mod.def("dlngstarSdlnT", &dlngstarSdlnT, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    The derivative of log(g*S) with respect to log(T), where g*S is the
    effective number of degrees of freedom for entropy density

    (By default uses the standard Gondolo-Gelmini g*S(T). To use other choices
    for g*S: evaluate Read_gstar(choice); see documentation for the function
    Read_gstar for more details.)
    )pbdoc", py::arg("T"));
    
    //RhoVisible(T)
    mod.def("RhoVisible", &RhoVisible, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    Energy density in the visible sector
    )pbdoc", py::arg("T"));
    
    //EntropyVisible(T)
    mod.def("EntropyVisible", &EntropyVisible, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    Comoving entropy in the visible sector
    )pbdoc", py::arg("T"));

    //HubbleVisible(T)
    mod.def("HubbleVisible", &HubbleVisible, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    Hubble rate in the visible sector
    )pbdoc", py::arg("T"));

    //SigmaV_chi(T, mchi, kappa, LambdaQCD)
    mod.def("SigmaV_chi", &SigmaV_chi, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV
    mchi: mass of the dark matter in GeV
    kappa: portal coupling
    LambdaQCD: QCD confinement scale in GeV. Set to 0.15 GeV by default

    Returns
    -------

    Thermally-averaged cross section for SM SMbar -> chi chibar process
    )pbdoc", py::arg("T"), py::arg("mchi"), py::arg("kappa"),
             py::arg("LambdaQCD")=0.15);

    //kappa_FreezeIn(mchi, LambdaQCD)
    mod.def("kappa_FreezeIn", &kappa_FreezeIn, R"pbdoc(
    Inputs
    ------

    mchi: mass of the dark matter in GeV
    LambdaQCD: QCD confinement scale in GeV. Set to 0.15 GeV by default

    Returns
    -------

    Portal coupling kappa that reproduces the observed dark matter relic
    abundance for dark matter frozen-in via a light dark photon mediator
    )pbdoc", py::arg("mchi"), py::arg("LambdaQCD")=0.15);
    
};
