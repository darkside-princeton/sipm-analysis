import numpy as np
from scipy.special import gamma,erf,erfc
from scipy.stats import chi2

def gauss(x,a,b,c):
    return a*np.exp(-(x-b)**2/(2*c**2))

def gauss_normalized(x,N,mu,sigma):
    return N*np.exp(-(x-mu)**2/(2*sigma**2))/np.sqrt(2*np.pi)/sigma

def line_simple(x,a,b):
    return a*x+b

def line(x,a,b):
    return a*(x-b)

def sipm_temp(x,V0,mu,sigma,tau,tau2):
    return V0/2.0 * erfc(1.0/np.sqrt(2.0) * (sigma/tau - (x-mu)/sigma)) * ( np.exp(0.5 * (sigma/tau)**2 - (x-mu)/tau) + np.exp(0.5 * (sigma/tau2)**2 - (x-mu)/tau2))

def pulse_jitter(t, a, tau, sigma, t0):
    return a*np.exp(sigma**2/2/tau**2)*np.exp(-(t-t0)/tau)*(1+erf((t-t0-sigma**2/tau)/sigma/np.sqrt(2)))/2

def compound_poisson(x,mu,p):
    k = [int(x_+0.5) for x_ in x]
    ans = []
    for k_ in k:
        if k_==0:
            ans.append(np.exp(-mu))
        else:
            ans_ = 0
            for i in range(1,k_+1):
                ans_ += gamma(k_+1)*gamma(k_)/gamma(i+1)/gamma(i)/gamma(k_-i+1)*(mu*(1-p))**i*p**(k_-i)
            ans.append(ans_*np.exp(-mu)/gamma(k_+1))
    return ans

def error_distance(df,sigma):
    return chi2.ppf(chi2.cdf(sigma**2,1),df)**0.5

