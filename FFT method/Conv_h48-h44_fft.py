import numpy as np

                # genfromtxt,fft,conj,log,isnan,

data = np.genfromtxt('BaF2_48_44.dat')

t = data[:,0]         # time - t
y_n = data[:,1]       # output signal - y(t)
g = data[:,2]         # impulse signal - g(t)

Y = np.fft.fft(y_n)      # fast Fourier transforms
G = np.fft.fft(g)        # for output and impulse (Y(w) and G(w))

lambdas = [0, 1e6, 1e7, 1e8, 1e9] # 4 6 8 11
                        #----------------------------------------------
                        # regularization parameter

f = []

for i in range(len(lambdas)):
    F = np.conj(G)*Y / (np.conj(G)*G + lambdas[i])
    decay = np.fft.ifft(F).real
    f.append(decay)

                # Fourier transform of decay function
                # and inverse Fourier transform of F
                # as an arrays for different regularizations
                 
#-----------------------------------------------------------------------
# write results in a .dat file
                
file1=open('Output_h48-h112_fft.dat','w')
file1.write('t'+'\t'+'f1'+'\t'+'f2'+'\t'+'f3'+'\t'+'f4'+'\n')

                # file with results, header and data for 4 lambdas [1-4]

for k in range(len(t)):
    file1.write(str(t[k])+'\t'+str(f[1][k])+'\t'+str(f[2][k])+'\t'+str(f[3][k])+'\t'+str(f[4][k])+'\n')
file1.close()


from pylab import polyfit, plot, subplot, hold, show, xlim, yscale, legend, fill_between

                     # a plot with 4 different
                     # regularization steps


subplot(2,1,1)
#yscale('log')
plot(t,y_n/max(y_n),'ko',label='Output')
plot(t,g/max(g),'mo',label='Impulse')
legend()

subplot(2,1,2)
yscale('log')

plot(t,f[1],'ro')
plot(t,f[2],'bo')
plot(t,f[3],'go')
plot(t,f[4],'mo')

show()

#-----------------------------------------------------------------------
# range used for regression - index i_fit, time t_fit and decay f_fit
# it's about 33 to 37 ns from the graph, and f[3] as the best option

i_fit = [i for i in range(len(t)) if t[i]>33 and t[i]<37]
t_fit = t[i_fit]
f_fit=np.log(f[3][i_fit])
    #----------------------
    # log scale will be used for linear regression f = (theta[0] + theta[1] * t)

theta,err = polyfit(t_fit[~np.isnan(f_fit)],f_fit[~np.isnan(f_fit)],1,cov=True)
                            #-----------------------------------
                            # NaN elements are being excluded

SD = np.sqrt(sum(((np.polyval(theta,t_fit)-f_fit)**2)[~np.isnan(f_fit)]/len(f_fit[~np.isnan(f_fit)])))

t_0 = 1/theta[0]
#-----------------------
# decay_kinetic_constant

se_theta0 = np.sqrt(np.diag(err)[0])  # standard error of theta[0] parameter
se_t_0 = se_theta0/theta[0]**2        # standard error for t_0

print('Decay_kinetic_constant t_0 =', "%.3g" % t_0,'+/-',"%.2g" % se_t_0)

#------plot fitting----------------------------

plot(t_fit,f_fit,'go')
hold('on')
plot(t_fit,np.polyval(theta,t_fit),'--r')
fill_between(t_fit,f_fit+SD,f_fit-SD,facecolor='cyan')

show()

