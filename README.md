# Decay kinetics in fast scintillators

This is an attempt to evaluate the decay kinetic constant for the integral luminescence of __BaF<sub>2</sub>__ nanoparticles upon high-energy X-ray excitation.

In case of fast scintillators like __BaF<sub>2</sub>__, luminescence decay times are short, of the same scale as the duration of excitation impulses. Therefore, the output signal, which is being measured, is composed of both decay and impulse. In order to exclude excitation one has to evaluate convolution equation __*y*(*t*) = *f*(*t*) ∗ *g*(*t*)__, where '∗'-sign means convolution, __*y*(*t*)__ - is the output (measured) signal, __*g*(*t*)__ - excitation impulse, *f*(*t*) - unknown function of luminescence decay kinetics.

In order to calculate the 'true' decay kinetics *f*(*t*) curve, decay kinetic constant *t*<sub>0</sub> and its error, two approaches have been used: nonlinear fitting (regression) and Fourier/inverse Fourier transforms with Tikhonov regulatization. Python 3.3 with numpy, scipy, matplotlib and lmfit extension packages were used.
