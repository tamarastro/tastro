import numpy as np
from matplotlib import pyplot as plt


def geodesic_polar(ll,k,theta0):
    r=np.sqrt(ll**2+k**2)
    theta = np.arctan(ll/k)+theta0
    return r,theta

def rtheta(k,theta,theta0):
    return k/np.cos(theta-theta0)


larr = np.arange(101)/100.-1.0
k = 1.0# np.arange(3)*10+1
theta0 = np.pi/4
tarr = np.arange(101)/100.*theta0*2



for k in np.arange(3)+1:
#    r,theta = geodesic_polar(larr,k,theta0)
#    x = r*np.cos(theta)
#    y = r*np.sin(theta)
#    plt.plot(x,y)
#    length = np.sqrt((x[-1]-x[0])**2 + (y[-1]-y[0])**2)
#    print(length)
    r_theta = rtheta(k,tarr,theta0)
    plt.plot(tarr,r_theta)
plt.show()
#plt.show()
