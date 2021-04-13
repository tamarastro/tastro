import numpy as np
from matplotlib import pyplot as plt

def lognormal(x,mu,sigma):
    return np.exp( -(np.log(x)-mu)**2/(2*sigma**2) ) / (x*sigma*np.sqrt(2*np.pi))

def normal(x,mu,sigma):
    return np.exp( -(x-mu)**2/(2*sigma**2) ) / (sigma*np.sqrt(2*np.pi))

x = np.arange(0.01,500.01,0.01)
mu = 3.0
sigma_arr = [1.0,0.5,0.2]
color_arr = ['r','g','b']
#norm = normal(x,mu,sigma)

m = 10.

i=0
for sigma in sigma_arr:
    #v = sigma**2
    #mu = np.log(m/np.sqrt(1+v/(m**2)))
    mu = np.log(m)-sigma**2/2
    lognorm = lognormal(x,mu,sigma)
    mean = np.exp(mu+sigma**2/2)
    total = np.sum(lognorm)
    meanarr = np.sum(x*lognorm)/total
    print(mean,meanarr,total)
    if i==0: mean0=mean
    meandiff = mean-mean0
    var = (np.exp(sigma**2)-1)*np.exp(2*mu+sigma**2)

    #plt.plot(np.log(x),np.exp(norm))
    plt.loglog(x,lognorm,color_arr[i])
    plt.plot([mean0,mean0],[0.0001,np.max(lognorm)],color_arr[i],linestyle=':')
    i=i+1


plt.ylim(1.e-4,np.max(lognorm))
plt.ylabel('Probability distribution')
plt.xlabel('Mass')
plt.savefig('lognorm.png')
plt.show()
