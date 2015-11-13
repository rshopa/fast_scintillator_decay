# Decay kinetics in fast scintillators

This is an attempt to evaluate the decay kinetic constant for the integral luminescence of BaF<sub>2</sub> nanoparticles upon high-energy X-ray excitation.

In case of fast scintillators like BaF<sub>2</sub>, luminescence decay times are short, of the same scale as the duration of excitation impulses. Therefore, the output signal, which is being measured, is composed of both decay and impulse. In order to exclude excitation one has to evaluate convolution equation *y*(*t*) = *f*(*t*) ∗ *g*(*t*), where '∗'-sign means convolution, *y*(*t*) - is the output (measured) signal, *g*(*t*) - excitation impulse, *f*(*t*) - unknown function of luminescence decay kinetics.

In order to calculate the 'true' decay kinetics *f*(*t*) curve, decay kinetic constant *t*<sub>0</sub> and its error, two approaches have been applied: nonlinear fitting (regression) and Fourier/inverse Fourier transforms with Tikhonov regularization. Python 3.3 along with numpy, scipy, matplotlib and lmfit extension packages were used as computational tools.

The environment for the code is crucial (especially about lmfit package version), so I specify it here:

> Windows 7 SP1 x64, Python 3.3.0

> Packages:
> * numpy 1.10.1
* matplotlib 1.5.0
* scipy 0.16.1
* lmfit 0.9.2 ('minimize' function works differently than in v.0.8.x!)

