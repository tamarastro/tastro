import numpy as np
from astropy.cosmology import FlatLambdaCDM as cosmo

cosmo1 = cosmo(H0=70, Om0=0.3)
print cosmo1

cosmo2 = cosmo(H0=72, Om0=0.28)
print cosmo2

print('But this might be better.... editing master.')