import numpy as np
from astropy.cosmology import FlatLambdaCDM as cosmo
cosmo1 = cosmo(H0=70, Om0=0.3)
print cosmo1