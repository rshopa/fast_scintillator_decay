import sys
from numpy import *         # yes, I know I had better specified 'em
from pylab import *
from scipy import integrate, signal, stats
from lmfit import minimize, Parameters, report_errors   #lmfit method

from numpy import random        # must override numpy (for random.choice)

#------------------------functions-------------------------------------

def Conv_residuals(params,x,y,g):     # computes residuals
                                    # y - convolution(f_model(t),g(t))
                                    
    a0 = params['I_0'].value          # initial values
    t0 = params['time0'].value      # of 2 parameters
    
    f_model = []                    # nonlinear regression model
    integrand = []                  # array of calculated y(t)
                                    # using f_model

    for j in range(len(x)):
        f_model.append(a0*exp(-(x[j])/t0))             
                                                        
    integrand = (convolve(f_model,g))[:len(x)]
    return y - integrand
                                
                                # an array of residuals is returned

#---------------------------------------------------------------------

def bootstrap(y,r):         # returns array with randomly sampled
                            # residuals added to y(t)
	rnd = []
	for z in random.choice(list(r),len(r),replace = True):
		rnd.append(z)

	return y + rnd

#----------------------------------------------------------------------

def main():

    #----------------------------------------------------------------------
    # import data
    #----------------------------------------------------------------------


    data = genfromtxt(sys.argv[1])

    t_c = data[2:,0]                    #time
    y_c = data[2:,1]                    #measured decay curve y(t)
    g_c = data[2:,2]                    #impulse g(t)

    #----------------------------------------------------------------------
    # lmfit section
    #----------------------------------------------------------------------

    params = Parameters()   
    params.add('I_0',value = 1.0)
    params.add('time0',value = 1.0)
                                    # minimize residuals using lmfit
                                    # with Levenberg-Marquardt method,
                                    # the first, "signal" fitting procedure

    rezult = minimize(Conv_residuals,
                      params,
                      args = (t_c,y_c,g_c),
                      method='leastsq')
        #--------------------------------------------------------------------
        # leastsq, nelder, lbfgsb, anneal, powell, cg, newton, cobyla, slsqp

    provisional_parameters = rezult      #will be printed at the end of the code

    residuals = Conv_residuals(rezult.params,t_c,y_c,g_c) #residuals
    y_lmfit = y_c + residuals                             #fitted line 


    #----------------------------------------------------------------------
    # bootstrapping section
    #----------------------------------------------------------------------


    I_0 = []                  # arrays of calculated parameters
    t_0 = []                # during bootstrapping

    iteration = int(input('How many bootstrapping iterations? '))

                                # number of boostrap procedures

    for i in range(iteration):
        
        params = Parameters()                # repeat lmfit minimize using 
        params.add('I_0',value = 10.0)         # Levenberg-Marquardt method
        params.add('time0',value = 1.0)      # with 
                                             # y(t) + randomly sampled residuals

        rezult = minimize(Conv_residuals,
                          params,
                          args = (t_c, bootstrap(y_c, residuals), g_c),
                          method ='leastsq')

        I_0.append(rezult.params['I_0'].value)
        t_0.append(rezult.params['time0'].value)


    #----------------------------------------------------------------------
    # quantiles/errors and means
    #----------------------------------------------------------------------

    err_I_0 = stats.mstats.mquantiles(I_0,[0.0, 1.0])
    err_t_0 = stats.mstats.mquantiles(t_0,[0.0, 1.0])

            # 0.25, 0.75 -quantiles may be used instead of the full span

    #----------------------------------------------------------------------
    # runs test
    #----------------------------------------------------------------------

    np = nm = 0     # number of positive and negative residuals, respectively
    nR = 1          # observed number of runs (changes of sign)


    if residuals[0] < 0:
        nm += 1

    for i in range(1,len(residuals)):       # loop for calculating
                                            # nm and nR
        if residuals[i] < 0:
            nm += 1

            if residuals[i-1] > 0:
                nR += 1

        elif residuals[i-1] < 0:
            nR += 1

    np = len(residuals) - nm                # np - number of positive residuals

    R = 1 + (2*np*nm)/(np + nm)     #expected number of runs

    sigma_R = sqrt(2*nm*np*(2*nm*np - np - nm)/((np + nm - 1)*(np + nm)**2))
                                    #variance of the expected number of runs

    if nR <= R:
        Z = (nR - R + 0.5)/sigma_R
    else:                               # estimated standard normal
        Z = (nR - R - 0.5)/sigma_R      # distribution (Z-score)


    #----------------------------------------------------------------------
    #report results of calculations
    #----------------------------------------------------------------------
        

    print('\nLMFIT report:\n')          # results from the 'signal' fit
    report_errors(provisional_parameters)

    print('\nBootstrapping report:\n\nI_0 =',"%.4f" % median(I_0),
          '\t(+',"%.4f" % (100*((median(I_0)-err_I_0[0])/median(I_0))),
          '% / -',"%.4f" % (100*((err_I_0[1]-median(I_0))/median(I_0))),'%)')

                        # NOTE! since the statistical approach
                        # has been used, medians are more relevant
                        # instead of means

    print('t_0 =',"%.4f" % median(t_0),
          '\t(+',"%.4f" % (100*((median(t_0)-err_t_0[0])/median(t_0))),
          '% / -',"%.4f" % (100*((err_t_0[1]-median(t_0))/median(t_0))),'%)\n')

    print('Runs test:\n\n Numbers of points:\n n_m =',nm,'\n n_p =',np,'\n\n'
          'Observed number of runs n_R =',nR,'\n'
          'Expected number of runs R =',"%.4f" % R,'+/-',"%.4f" % sigma_R,'\n'
          'The standard normal distribution score Z =',"%.4f" % Z)


    #----------------------------------------------------------------------
    # plot section
    #----------------------------------------------------------------------

    f_c=[]              # luminescence decay as exp-model

    for i in range(len(t_c)):
            f_c.append(median(I_0)*exp(-(t_c[i])/median(t_0)))


    suptitle(r'Decay kinetics of BaF$_2$ 78 nm nanoparticles',fontsize=18)

    subplot(211)
    plot(t_c,y_c/max(y_c),'bo',label = r'$y(t)$')       # all graphs are normalized
    plot(t_c,g_c/max(g_c),'ro',label = r'$g(t)$')
    plot(t_c,y_lmfit/max(y_lmfit),'m-',label = 'fitting curve')
    plot(t_c,f_c/max(f_c),'g.--',label = r'$f(t)$')

    xlabel('Time (ns)',fontsize = 15)
    ylabel('Intensity (a.u.)',fontsize = 16)
    legend(loc = 1)

    subplot(212)
    stem(t_c,residuals,linefmt='g--',markerfmt='bs',basefmt='r-')
    xlabel('Time (ns)',fontsize = 15)
    ylabel(r'Residuals $y - y_{model}$',fontsize = 16)

    subplots_adjust(hspace = 0.3, wspace = 0.3, right = 0.95, top = 0.92)
    show()

    #------------Histograms----------------------------

    suptitle(r'Decay kinetics of BaF$_2$ 78 nm nanoparticles',fontsize=18)

    subplot(121)
    hist(I_0, color = 'green')
    xlabel(r'$I_0$ (a.u.)',fontsize = 16)
    ylabel(r'Frequency',fontsize = 16)

    subplot(122)
    hist(t_0, color = 'green')
    xlabel(r'$t_0$ (ns)',fontsize = 16)
    ylabel('Frequency',fontsize = 16)

    subplots_adjust(hspace = 0.4, left = 0.1, right = 0.95, top = 0.92)
    show()

#-------------------------------------------------------------------

if __name__ == '__main__':
     main()

