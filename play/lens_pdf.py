import numpy as np
from matplotlib import pyplot as plt

def sigma_gauss_func(mu_logn,sigma_logn):
    return np.log(1+sigma_logn**2/mu_logn**2)

def mu_gauss_func(mu_logn,sigma_gauss):
    return 0.5*np.log((mu_logn**2)/(1+sigma_gauss**2/mu_logn**2))

def p_kappa(kappa,mu_logn, sigma_logn, kappa_mean):
    sigma_gauss = sigma_gauss_func(mu_logn,sigma_logn)
    mu_gauss = mu_gauss_func(mu_logn,sigma_gauss)
    numerator = np.exp( -( np.log(kappa+mu_logn) - kappa_mean - mu_gauss)/(2*sigma_gauss) )
    return numerator/(np.sqrt(2*np.pi)*sigma_gauss*(kappa+mu_logn-kappa_mean))

kappa=np.arange(101)/50.-1.2

p_kappa_result = p_kappa(kappa, 1.0,1.0,1.0)

print(kappa)

plt.plot(kappa, p_kappa_result, 'r-o')

plt.show()



