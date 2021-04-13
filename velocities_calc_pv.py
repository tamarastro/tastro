#################
# Program to calculate peculiar velocities given an observed redshift and observed distance.
# It assumes a cosmological model. 

import numpy as np
from matplotlib import pyplot as plt
from scipy import integrate
from numpy.polynomial import Polynomial as P

from functions import funcs_wa as wa

#################
# Set constants
c=299792.458  # km/s
H0fid = 70. # km/s/Mpc
omfid=0.3
oxfid=0.7
okfid = 1.0-omfid-oxfid

# Set array of redshifts,
nz = 1001
zcs = np.linspace(0.001,0.201,num=nz) # array of cosmological redshifts

# Choose a peculiar velocity 
vp = 300.0 #km/s
zp = np.sqrt((1+vp/c)/(1-vp/c))-1 #All peculiar velocities are 300km/s

# Calculate the array of observed redshifts
zos = (1+zcs)*(1+zp) -1 # array of observed redshifts

# Calculate distances (first comoving, then proper)
xzcs = np.zeros(nz)
for i, zc in enumerate(zcs): 
	xzcs[i]=wa.xx(zc,om=omfid,ox=oxfid,w0=-1.0,wa=0.0,orr=0.0)
Dzcs = (c/H0fid)*wa.dist_curve_correct(xzcs, okfid) 

# Make some Distances with uncertainties
nz_obs = 50
zcs_obs = np.linspace(0.002,0.20,num=nz_obs) # array of cosmological redshifts
zos_obs = (1+zcs_obs)*(1+zp) -1 # array of observed redshifts after peculiar velocities
xzcs_obs = np.zeros(nz_obs)
for i, zc in enumerate(zcs_obs): 
	xzcs_obs[i]=wa.xx(zc,om=omfid,ox=oxfid,w0=-1.0,wa=0.0,orr=0.0)
Dzcs_obs = (c/H0fid)*wa.dist_curve_correct(xzcs_obs, okfid)

add_error=False
if add_error:
	Dzcs_obs = Dzcs_obs *(1.0+np.random.randn(nz_obs)*0.01) # flat universe


###########################################################
# Now imagine that we observe Dzcs (Tully Fisher), and zos.  
# We go from Dzcs and calculate the redshift that should correspond to. 
# Then compare that to the observed redshift to get the zp.

# Test a single recovered z 
i=10 # Just choose one of the redshifts for now as a test
Dtest = Dzcs_obs[i]
zpred = np.interp(Dtest,Dzcs,zcs)  # Interpolate the distance-redshift relation to the chosen distance to get the predicted redshift
zp_meas = (1+zos[i])/(1+zpred)-1   # Calculate the peculiar redshift
vp_meas = c*((1+zp_meas)**2-1)/((1+zp_meas)**2+1)  # Calculate the peculiar velocity

print('Rel     z=',zp)
print('Rel zpmeas=',zp_meas)
print('Rel vpmeas=',vp_meas)

print('## Test non-rel approx ##')
print('Non-rel z=',300.0/c)  # Just have a look at how different it would be if we took the non-relativistic approx.
print('Non-rel v=',c*zp)
print('Non-rel vmeas=',c*zp_meas)

#plt.plot(zcs_fine,Dzcs_fine)
#plt.plot(zcs,Dzcs,':')
#plt.plot([zpred,zpred],[Dtest,Dtest],'o')
#plt.show()

##############################################
# Do for whole array instead of a single value
zpred = np.interp(Dzcs_obs,Dzcs,zcs)  # Interpolate the distance-redshift relation to the chosen distance to get the predicted redshift
zp_meas = (1+zos_obs)/(1+zpred)-1     # Measure the peculiar redshift using the observed redshift and redshift predicted from the observed distance
vp_meas = c*((1+zp_meas)**2-1)/((1+zp_meas)**2+1)
print('Array of measured peculiar velocities',vp_meas)

########################################
# Now test the velocity addition formula

# Use the observed redshift to predict the recession velocity, assuming zo = zc
xzos = np.zeros(nz)
for i, zc in enumerate(zcs): 
	xzos[i]=wa.xx(zos[i],om=omfid,ox=oxfid,w0=-1.0,wa=0.0,orr=0.0)
vt_from_zobs        = c*xzos  # This is instead of the common c*z
vt_from_zobs_approx = c*zos
h=0.7
vt_from_zobs_Springob = 99.939*(xzos*c/100)+0.00818*(xzos*c/100)**2

# Use the observed distance to predict the recession velocity, v=H0*D
vr_from_Dobs = H0fid * Dzcs

# Subtract one from the other to get peculiar velocity
vp_add        = vt_from_zobs         - vr_from_Dobs # This seems to be the wrong way around to me, but need to check.
vp_add_approx = vt_from_zobs_approx  - vr_from_Dobs
vp_Springob   = vt_from_zobs_Springob- vr_from_Dobs

# Plot the derived peculiar velocity vs redshift, for the various methods
plt.axhline(y=300.0, color="grey", linestyle="-",label='True')
plt.plot(zcs_obs,vp_meas,'.',label='(1+z) calculation')
plt.xlabel('Redshift')
plt.ylabel('Peculiar Velocity derived (km/s)')
plt.plot(zcs,vp_add       ,':' ,label='vt-HD first iteration')
plt.plot(zcs,vp_add_approx,'--',label='cz-HD approx')
plt.plot(zcs,vp_Springob,'-.',label='Springob-HD approx')
plt.legend(loc='lower right')
plt.ylim(250,400)
plt.xlim(0,0.2)
plt.savefig('calc_pv.png',bbox_inches='tight')
plt.close()


v_from_Springob_approx = 99.939*(h*Dzcs)+0.00818*(h*Dzcs)**2
plt.plot(zcs,vr_from_Dobs,label='true',color='black')
plt.plot(zcs,c*zcs,'-',label='cz approx',color='red')
plt.plot(zcs,v_from_Springob_approx,'--',label='Springob 2014 Eq16',color='orange',linewidth=2)
plt.plot(zcs,wa.vvapprox(zcs),'--',label='q0 j0 formula',color='green',linewidth=3)
plt.legend(loc='lower right')
plt.xlim(0.0,0.2)
plt.ylim(0.0,60000)
plt.xlabel('Redshift')
plt.ylabel('Velocity (km/s)')
plt.savefig('approx_pv.png')
plt.close()
#Dzo = (c/H0fid)*wa.dist_curve_correct(xzo, 1.0)
#vthry = c*wa.dist_curve_correct(xzobs,1.0)

#print('zp=',zp*c)

#vp = H0fid*Dobs - vthry 

#print('vp=',vp)

# Plot distance vs cosmological redshift, and vs observed redshift. 
#plt.plot(zcs,Dzcs-Dzcs_poly(zcs),'.')
#plt.plot(zcs,Dzcs_poly(zcs),'.')
#plt.plot(zos,Dobs,'.')
#plt.show()

# Scraps
#Dzcs_polynomial = P.fit(zcs,Dzcs,12)
#Dzcs_poly = np.poly1d(np.polyfit(zcs,Dzcs,12))
#zi = (Dzcs_polynomial-Dzcs[i]).roots()
#print(zcs[i],zi)
#exit()
