# Solution 1: FFT/inverse FFT

The first solution method to derive the "pure" signal of the luminescence emitted from the fast scintillator is using [Fourier/inverse Fourier transform](https://en.wikipedia.org/wiki/Fourier_transform "Fourier transform").

Fourier transform decomposes an arbitrary function of time $$y(t)$$ into the corresponding function $$Y(\omega)$$ of frequencies $$\omega$$ which form $$y(t)$$. The useful hint in our case is that in the space of frequencies $$\omega$$ a convolution operation is replaced by simple multiplication. Hence the following transformation has to be made:

$$y(t)=f(t)*g(t) \rightarrow Y(\omega)=F(\omega)\cdot G(\omega)$$, where $$y,f,g$$ chars correspond to output signal, luminescence and excitation impulse, respectively.

From this equation, $$F(\omega)=Y(\omega)/G(\omega)$$. Therefore by performing the inverse operation called "inverse Fourier transform", one can obtain luminescence $$f(t)$$.

However, real measurements imply errors, so direct calculation will eventually provide inappropriate [noise-like results](https://github.com/rshopa/fast_scintillator_decay/blob/master/FFT%20method/FFT_no_regularization.png?raw=true "Figure 1"). This is also caused by applying simplified __fast Fourier transform (FFT)__ and its inverse analogue.

One of the possible solution to the problem is regularization. Here a modified [Tikhonov regularization model](http://www.ees.nmt.edu/outside/courses/GEOP505/Docs/deconv.pdf "Time Series/Data Processing and Analysis") has been used. The function $$F(\omega)$$ is multiplied and divided by conjugate $$G^*(\omega)$$: $$F(\omega)=\frac{Y(\omega)\cdot G^*(\omega)}{G(\omega)\cdot G^*(\omega)}$$. Later a small parameter $$\lambda$$ (regularization parameter) is added in order to remove noise:

$$F(\omega)=\frac{Y(\omega)\cdot G^*(\omega)}{G(\omega)\cdot G^*(\omega)+\lambda}$$.

By exploring different values of $$\lambda$$, one can obtain a value, good enough to obtain $$f(t)$$. From the curves in the next figure it is clear, that the best value is $$\lambda=10^8$$.

![figure 2](https://github.com/rshopa/fast_scintillator_decay/blob/master/FFT%20method/FFT_w_regularization.png?raw=true "Figure 2")

The simplest model for the luminescence decay is $$f(t)=I_0\exp(-t/t_0)$$. Log-scaled intensity behaves linearly in the $$33\div37$$ nm span. Therefore, just performing linear fitting (regression), the decay kinetic constant $$t_0$$ could be obtained.

The input data file "BaF2_78nm.dat" contains of three columns: time $$t$$ (ns), output $$y$$ and excitation impulse $$g$$. The code file "FFT_regularization.py" performs FFT/inverse FFT transform using Tikhonov regularization on the arbitrary file __FFT_regularization.py <file name>__ to obtain unknown integral luminescence decay function for polystyrene:BaF $$_2$$ composition (BaF $$_2$$ nanoparticles have mean size of 78 nm). It prints out four different functions $$f(t)$$ to file "Output_BaF2_78nm.dat" for different regularization steps $$\lambda$$. The program also performs linear regression in order to calculate the decay kinetic constant $$t_0$$ (printed to stdout), using log-scaled data from $$f(t)$$.