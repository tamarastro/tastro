# Needs the following to have been imported
from scipy import integrate
import numpy as np


def ezinv(z,om=0.3,ox=0.7,w0=-1.0,wa=0.0,orr=0.0):
    ok = 1.-om-ox-orr
    a  = 1./(1.+z)
    f = a**(3*(1.+w0+wa)) * np.exp(3*wa*(1-a))
    ez = np.sqrt( orr/a**4 + om/a**3 + ok/a**2 + ox/f )
    return 1./ez

def xx(z,om=0.3,ox=0.7,w0=-1.0,wa=0.0,orr=0.0):
    xx,err = integrate.quad(ezinv,0.,z,args=(om,ox,w0,wa,orr))
    # Multiply by R_0 (i.e. c/H0 if flat or c/H0*sqrt(abs(1/ok)) to get in units of Mpc.
    # Remember to curvature correct if want proper distance.
    return xx

def vv(xx,c=299792.458): # c in km/s
    return xx*c

def vvapprox(z,q=-0.55,j=1.0,c=299792.):
    return c*z/(1+z)*( 1+0.5*(1-q)*z - (1./6)*(1-q-3*q**2+j)*z**2 )

def get_zp(vp):
    return np.sqrt((1+vp/c)/(1-vp/c))-1


def dist_curve_correct(xx, ok):
    # Corrects the comoving distance, x, for curvature, ok.
    # Result is curve-corrected distance / (c/H0)
    if ok < 0.0:
        dk = np.sin(np.sqrt(-ok)*xx)/np.sqrt(-ok)
    elif ok > 0.0:
        dk = np.sinh(np.sqrt(ok)*xx)/np.sqrt(ok)
    else:
        dk = xx
    return dk

def dist_lum(xx,ok,z):
    return dist_curve_correct(xx,ok)*(1+z)

def H0kmsmpc2Gyr(H0kmsmpc):
    # Convert to inverse seconds
    H0s = H0kmsmpc * 3.24e-20 #s-1
    # Convert to inverse Giga-years
    H0Gyr = H0s* 3.156e16 #Gyr-1
    return H0Gyr

def adotinv(a,om=0.3,ox=0.7,w0=-1.0,wa=0.0,orr=0.0,H0kmsmpc=70.):
    # Calculates 1/(a*dot(a)).
    ok = 1.-om-ox-orr
    H0Gyr=H0kmsmpc2Gyr(H0kmsmpc)
    adot = H0Gyr * a * np.sqrt ( orr/a**4 + om/a**3 + ok/a**2 + ox/( a**(3*(1+w0+wa))*np.exp(3*wa*(1-a)) ) )
    return 1./adot

def time(a0, a1, om=0.3,ox=0.7,w0=-1.0,wa=0.0,orr=0.0,H0kmsmpc=70.):
    # Calculates the time between two scalefactors.
    time = integrate.quad(adotinv,a0,a1,args=(om,ox,w0,wa,orr,H0kmsmpc))
    return time 



