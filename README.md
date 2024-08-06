# FreezeIn v2.0

### Prudhvi N. Bhattiprolu, Robert McGehee, Evan Petrosky, and Aaron Pierce

This package computes the portal coupling $\kappa$ that reproduces the observed relic abundance, for dark matter $\chi$ frozen-in via a light dark photon mediator, as a function of its mass $m_\chi$, while accounting for plasmon decays (see Refs. [1, 2]). It also gives the corresponding direct detection cross section off electrons $\overline \sigma_e$.

The freeze-in model is a benchmark for ongoing direct detection experiments. In some cases, the literature for this benchmark has contained errors. Ref. [1] provides the corrected predictions and the corresponding code is available in this repository. Ref. [2] also accounted for plasmon decays in the predictions of the freeze-in model (based on Ref. [3]), relevant for sub-MeV dark matter masses.

## Data

The data for the freeze-in $\kappa$ and the corresponding $\overline \sigma_e$ as a function of $m_\chi$ in GeV, is included in this repository in the file [`Data/FreezeIn.txt`](https://github.com/prudhvibhattiprolu/FreezeIn/blob/main/Data/FreezeIn.txt). This data was generated using the standard Gondolo-Gelmini $g_{\ast(, s)} (T)$ located at [`gstar/std.tab`](https://github.com/prudhvibhattiprolu/FreezeIn/blob/main/gstar/std.tab). Furthermore, we have set the QCD scale $\Lambda_\text{QCD}$ to 0.15 GeV for these calculations.

[<img src="https://github.com/prudhvibhattiprolu/FreezeIn/blob/main/Data/KappaFI.png">](Data/KappaFI.png)

[<img src="https://github.com/prudhvibhattiprolu/FreezeIn/blob/main/Data/SigmaDDeFI.png">](Data/SigmaDDeFI.png)

The shaded band corresponds to variation of $\Lambda_\text{QCD} = 0.15 \pm 0.05$ GeV.

For sub-MeV dark matter masses, plasmons quickly become the dominant source of freezing-in dark matter.

[<img src="https://github.com/prudhvibhattiprolu/FreezeIn/blob/main/Data/SubMeVFI.png">](Data/SubMeVFI.png)

The dot-dashed curve corresponds to the freeze-in line neglecting the effect of plasmons.

## Code

The source code for this package, written in C++, is located at [`src/FreezeIn/FreezeIn.h`](https://github.com/prudhvibhattiprolu/FreezeIn/blob/main/src/FreezeIn/FreezeIn.h). The program located at [`src/FreezeIn/FreezeInPy.cc`](https://github.com/prudhvibhattiprolu/FreezeIn/blob/main/src/FreezeIn/FreezeInPy.cc) exposes the functions in the FreezeIn.h library to python using a light-weight header-only [pybind11](https://pybind11.readthedocs.io/en/stable/) library.

The installation requires Python3 (> 3.7) and a C++ compiler (e.g. gcc). Additional C++ headers from [boost](https://www.boost.org/) and [pybind11](https://pybind11.readthedocs.io/en/stable/) required for the installation are already included in this package (located at `extern/boost` and `extern/pybind11`), and need not be installed seperately.

### Automatic Installation

To automatically install FreezeIn as a python package using pip, run

```bash
python -m pip install "git+https://github.com/prudhvibhattiprolu/FreezeIn.git#egg=FreezeIn"
```

and the package can be readily imported from python.

### Manual Installation

Instead, if you want to edit the source code and/or build it into your own code/package (or if the automatic installation did not successfully go through) the package can also be manually cloned and installed:

```bash
# Git clone
git clone https://github.com/prudhvibhattiprolu/FreezeIn.git

# cd into the package
cd FreezeIn
```
Finally, to build the package:

```bash
# build the package manually by doing
python setup.py --quiet build_ext --inplace clean --all
```

If the installation is successful, a shared object file should appear at `src/FreezeIn/FreezeIn.*.so`. To start using the package in python, `cd src/` and launch python to import `FreezeIn` as a package.

### Usage

To use the package

```python
import FreezeIn
```

All the functions (along with the documentation) can be listed by using the Python help function

```python
help(FreezeIn)
```

The documentation for each function can be accessed by doing

```python
# To access the documentation of the function kappa_FreezeIn, e.g., we can either do
help(FreezeIn.kappa_FreezeIn)
# or
print(FreezeIn.kappa_FreezeIn.__doc__)
```

We can now compute the freeze-in $\kappa$ and the corresponding $\overline \sigma_e$ using the following functions:

`kappa_FreezeIn(mchi, LambdaQCD=0.15, IncludePlasmons=True)`:
computes the portal coupling, kappa, that reproduces the observed dark matter relic abundance for dark matter frozen-in via a light dark photon mediator as a function of the dark matter mass `mchi` in GeV. The QCD scale `LambdaQCD`, set to 0.15 GeV by default, can also be changed. Plasmon decay effects are included by default but can be excluded by setting `IncludePlasmons` to `False`.

`SigmaDDe(mchi, kappa)`:
computes the direct detection cross section (in cm^2) through the light dark photon mediator as a function of the dark matter mass `mchi` in GeV and the portal coupling `kappa`


In addition to the above functions, there are several other functions, e.g,

`SigmaV_chi(T, mchi, kappa, LambdaQCD=0.15, IncludePlasmons=True)`:
computes the effective thermally-averaged cross section for dark matter freeze-in $\langle \sigma v \rangle_\text{eff}$ (see Ref. [2]) as a function of the temperature in the visible sector `T` in GeV, the dark matter mass `mchi` in GeV, and the portal coupling `kappa`. The QCD scale `LambdaQCD` is set to 0.15 GeV by default. Plasmon decay effects are included by default but can be excluded by setting `IncludePlasmons` to `False`.

`RhoVisible(T)`:
computes the energy density in the visible sector as a function of the temperature `T` in GeV

`EntropyVisible(T)`:
computes the entropy density in the visible sector as a function of the temperature `T` in GeV

`Hubble(T)`:
computes the Hubble rate as a function of the temperature `T` in GeV

`gstar(T)`:
computes the effective number of degrees of freedom for energy density as a function of the temperature `T` in GeV

`gstarS(T)`:
computes the effective number of degrees of freedom for entropy density as a function of the temperature `T` in GeV

`dlngstarSdlnT(T)`:
computes the derivative of log(gstarS) with respect to log(T) as a function of the temperature `T` in GeV

All of the above functions use the standard Gondolo-Gelmini $g_{\ast(,s)}(T)$ by default. To use other choices, evaluate the following function:

`Read_gstar(choice="standard", gstarpath="<path to FreezeIn repository>/gstar")`:
Reads tabulated data for effective number of degrees of freedom from various .tab files in the `gstar/` folder with three columns: {Temperature in GeV, gstarS, gstar}. The `gstarpath` parameter is by default set to the path to the `gstar/` folder provided with this repository. The `choice` parameter can be set to one of the following:

* "standard": Gondolo-Gelmini (default)
* "HP_A": Hindmarsh-Philipsen equation of state A
* "HP_B": Hindmarsh-Philipsen equation of state B
* "HP_B2": Hindmarsh-Philipsen equation of state B2
* "HP_B3": Hindmarsh-Philipsen equation of state B3
* "HP_C": Hindmarsh-Philipsen equation of state C

(Various `gstar/*.tab` files are taken directly from the [MicrOMEGAs](https://lapth.cnrs.fr/micromegas/) package. See the documentation for MicrOMEGAS for more details on each of the above choices for $g_{\ast(,s)}$.

If you want to use your own .tab file for gstar(S), you can include it in the `gstar` folder once the repository is cloned. Then edit the `Read_gstar` function in `src/FreezeIn/FreezeIn.h` to include your own .tab file, and build the code manually as described above. Finally evaluate the `Read_gstar` function in python.)

#### Plasmon Properties

An implementation of the plasmon properties discussed in Ref. [2] is made available with v2.0 of this package. For example,

`Omegap(T)`:
Plasma frequency (typical plasmon mass) in GeV as a function of the temperature `T` in GeV

`vstar(T)`:
Typical electron velocity in the plasma as a function of the temperature `T` in GeV

[<img src="https://github.com/prudhvibhattiprolu/FreezeIn/blob/main/Data/PlasmonProperties.png">](Data/PlasmonProperties.png)

`PlasmonMasst(k, T)`:
Transverse plasmon mass in GeV as a function of the wavevector `k` in GeV and the temperature `T` in GeV

`PlasmonZt(k, T)`:
Transverse plasmon renormalization factor as a function of the wavevector `k` in GeV and the temperature `T` in GeV

[<img src="https://github.com/prudhvibhattiprolu/FreezeIn/blob/main/Data/TransversePlasmon.png">](Data/TransversePlasmon.png)

`PlasmonMassl(k, T)`:
Longitudinal plasmon mass in GeV as a function of the wavevector `k` in GeV and the temperature `T` in GeV

`PlasmonZl(k, T)`:
Longitudinal plasmon renormalization factor as a function of the wavevector `k` in GeV and the temperature `T` in GeV

[<img src="https://github.com/prudhvibhattiprolu/FreezeIn/blob/main/Data/LongitudinalPlasmon.png">](Data/LongitudinalPlasmon.png)

For faster evaluations of the plasmon properties mentioned above, the functions `Omegap(T)`, `vstar(T)`(, and `Omega1(T)`) are retrieved from the look-up table `gstar/OmegaInterpolation.tab` using 1D interpolation, with the default evaluation of the following function:

`Read_gstar(gstarpath="<path to FreezeIn repository>/gstar")`:
Read tabulated data for {Omega1, Omegap, vstar} as a function of temperature from `gstar/OmegaInterpolation.tab`. The `gstarpath` parameter is by default set to the path to the `gstar/` folder provided with this repository

## References

[1] P. N. Bhattiprolu, R. McGehee, A. Pierce, “A Dark Sink Enhances the Direct Detection of Freeze-in Dark Matter,” arXiv: 2312.14152 [hep-ph].

[2] P. N. Bhattiprolu, R. McGehee, E. Petrosky, A. Pierce, “Sub-MeV Dark Sink Dark Matter,” arXiv: 2408.***** [hep-ph].

[3] C. Dvorkin, T. Lin, K. Schutz, "Making dark matter out of light: freeze-in from plasma effects," arXiv: 1902.08623 [hep-ph]

