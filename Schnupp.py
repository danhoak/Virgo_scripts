#! /usr/bin/python

from numpy import *

import sys, os, time
sys.path.append('/users/swinkels/deploy/PythonVirgoTools/trunk/src')
from virgotools import getChannel

import matplotlib
matplotlib.use("Agg", warn=False)
import pylab
from scipy import signal
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.gridspec as gridspec
matplotlib.rcParams.update({'savefig.dpi':250})
from matplotlib.font_manager import FontProperties
from scipy import constants
from uncertainties import ufloat
from uncertainties import unumpy



## Parametric function: 'v' is the vector of parameters, 'x' the independent variable
def decay(t,A,m):
    return A + t*m




c = constants.c
lamb = 1064.0e-9





#### Code to demodulate data


def decimate(x, q, n=None, ftype='iir', axis=-1):
    if not isinstance(q, int):
        raise TypeError("q must be an integer")
    if n is None:
        if ftype == 'fir':
            n = 30
        elif q>32:
            n = 2
        else:
            n = 3
    if ftype == 'fir':
        b = signal.firwin(n + 1, 1. / q, window='hamming')
        a = 1.
    else:
        b, a = signal.cheby1(n, 0.05, 0.8 / q)
        #b, a = signal.butter(4, 0.8 / q, 'low', analog=False)

    #print b, a
    y = signal.filtfilt(b, a, x)
    sl = [slice(None)] * y.ndim
    sl[axis] = slice(None, None, q)
    return y[sl]


def bandpass(x,f1,f2,fs):
    nyq = fs/2
    b, a = signal.butter(4, [f1/nyq, f2/nyq], 'bandpass')
    y = signal.filtfilt(b, a, x)
    return y


def lowpass(x,fc,fs):
    nyq = fs/2
    b, a = signal.butter(4, fc/nyq, 'lowpass')
    y = signal.filtfilt(b, a, x)
    return y


def data_demod(t,xx,f0,Fs,duration):

    demod_I = zeros(len(t))
    demod_Q = zeros(len(t))

    # lowpass
    z = lowpass(xx,200.0,Fs)

    # decimation frequency
    Fs0 = 1000.0

    # decimate data
    w = decimate(z, int(Fs/Fs0))
    t_dec = arange(0,duration,1.0/Fs0)
    
    # first bandpass the data around f0
    zz = bandpass(w,f0-1,f0+1,Fs0)

    # demodulate and low-pass
    sinewave = sin(2*pi*f0*t_dec)
    cosinewave = cos(2*pi*f0*t_dec)

    Iphase = zz*sinewave
    Qphase = zz*cosinewave
    Ipass = lowpass(Iphase,0.5,Fs0)
    Qpass = lowpass(Qphase,0.5,Fs0)
    xdemod = sqrt(Ipass**2 + Qpass**2)
    
    # only use part of the time window to avoid filter transients
    i1 = Fs0*2
    i2 = Fs0*(duration-2)

    return xdemod




#gps_start = 1172947197
#gps_start = 1172948800
#gps_start = 1172949717
gps_start = 1173376317
duration = 200





fignum=0

fignum+=1
pylab.figure(fignum)


###############################################
### Grab the B1p data
chan = 'V1:LSC_B1p_56MHz_Q'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
ts = arange(0,duration,1./1000)
line = data_demod(ts,x.data,84.1,Fs,duration)

# Find the indices for each half of the data
idx1 = argmin(abs(ts-95.0))
idx2 = argmin(abs(ts-110.0))

# Fit each half
param_start = [1.0,1.0]
pfit1, fitCov = curve_fit(decay,ts[Fs:idx1],line[Fs:idx1],param_start)
perr = sqrt(diag(fitCov))
pfit2, fitCov = curve_fit(decay,ts[idx2:-1*Fs],line[idx2:-1*Fs],param_start)
perr = sqrt(diag(fitCov))

# Find the times of the zero crossings from each half-data fit
y1 = decay(ts,pfit1[0],pfit1[1])
id1 = argmin(abs(y1))
b1pA = ts[id1]
y2 = decay(ts,pfit2[0],pfit2[1])
id2 = argmin(abs(y2))
b1pB = ts[id2]
tB1p = ts
#################################################

fignum+=1
pylab.figure(fignum)

ax1 = pylab.subplot(2,1,1)
pylab.plot(ts,line,'r.',label='B1p 56MHz Q Demod',markersize=5)
pylab.plot(ts,decay(ts,pfit1[0],pfit1[1]),'k-')
pylab.plot(ts,decay(ts,pfit2[0],pfit2[1]),'k-')
pylab.grid(True)
pylab.xticks(visible=False)
#pylab.ylabel('DC Power [mW]',fontsize=10)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
pylab.ylim([0.0,1e-4])
pylab.xlim([2,198])
pylab.legend(fancybox=True,loc=2,prop={'size':8},numpoints=1)
ax1.get_yaxis().set_label_coords(-0.06,0.5)

pylab.title('NArm Schnupp Measurement')

###############################################
### Grab the B4 data
chan = 'V1:LSC_B4_56MHz_Q'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
ts = arange(0,duration,1./1000)
line = data_demod(ts,x.data,84.1,Fs,duration)

# Find the indices for each half of the data
idx1 = argmin(abs(ts-95.0))
idx2 = argmin(abs(ts-110.0))

# Fit each half
param_start = [1.0,1.0]
pfit1, fitCov = curve_fit(decay,ts[Fs:idx1],line[Fs:idx1],param_start)
perr = sqrt(diag(fitCov))
pfit2, fitCov = curve_fit(decay,ts[idx2:-1*Fs],line[idx2:-1*Fs],param_start)
perr = sqrt(diag(fitCov))

# Find the indices of the minima
y1 = decay(ts,pfit1[0],pfit1[1])
id1 = argmin(abs(y1))
b4A = ts[id1]
y2 = decay(ts,pfit2[0],pfit2[1])
id2 = argmin(abs(y2))
b4B = ts[id2]
tB4 = ts
#################################################


ax2 = pylab.subplot(2,1,2)

pylab.plot(ts,line,'r.',label='B4 56MHz Q Demod',markersize=5)
pylab.plot(ts,decay(ts,pfit1[0],pfit1[1]),'k-')
pylab.plot(ts,decay(ts,pfit2[0],pfit2[1]),'k-')
pylab.grid(True)
#pylab.xticks(visible=False)
#pylab.ylabel('DC Power [mW]',fontsize=10)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
pylab.ylim([0.0,1.2e-5])
pylab.legend(fancybox=True,loc=2,prop={'size':8},numpoints=1)
pylab.xlim([2,198])



###################################################

# Now get the phase data



print

chan = 'V1:SDB2_B1p_56MHz_phi0'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
t1 = arange(0,duration,1./Fs)
b1p_phi = x.data

chan = 'V1:SPRB_B4_56MHz_phi0'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
t1 = arange(0,duration,1./Fs)
b4_phi = x.data

yy = interp1d(t1,b1p_phi)
y1 = yy(b1pA)
y2 = yy(b1pB)
b1p_narm = ufloat((y1+y2)/2.,abs(y1-y2))
print 'NArm B1p:', b1p_narm

"""
yy = interp1d(t1,b4_phi)
y1 = yy(b4A)
y2 = yy(b4B)
b4_narm = ufloat((y1+y2)/2.,abs(y1-y2))
print 'NArm B4:', b4_narm


ax3 = pylab.subplot(3,1,3)

pylab.plot(t1,b1p_phi*180/pi + 37.0,'g--',label='B1p 56MHz Demod Phase')
pylab.plot(t1,b4_phi*180/pi - 3.0,'k--',label='B4 56MHz Demod Phase')
pylab.grid(True)
#pylab.xticks(visible=False)
#pylab.ylabel('Optical Gain [arb]',fontsize=10)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
pylab.ylim([-20.,20.])
#pylab.xlim([20,190])
pylab.legend(fancybox=True,loc=3,prop={'size':8},numpoints=1)
ax3.get_yaxis().set_label_coords(-0.06,0.5)
"""

pylab.xlabel('Time [sec]',fontsize=10)

pylab.savefig('Schnupp_narm.png')





# Now the west arm



#gps_start = 1172953457
#gps_start = 1172955727
gps_start = 1173378742
duration = 200



###############################################
### Grab the B1p data
chan = 'V1:LSC_B1p_56MHz_Q'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
ts = arange(0,duration,1./1000)
line = data_demod(ts,x.data,84.1,Fs,duration)

# Find the indices for each half of the data
idx1 = argmin(abs(ts-100.0))
idx2 = argmin(abs(ts-115.0))

# Fit each half
param_start = [1.0,1.0]
pfit1, fitCov = curve_fit(decay,ts[Fs:idx1],line[Fs:idx1],param_start)
perr = sqrt(diag(fitCov))
pfit2, fitCov = curve_fit(decay,ts[idx2:-0.3*Fs],line[idx2:-0.3*Fs],param_start)
perr = sqrt(diag(fitCov))

# Find the times of the zero crossings from each half-data fit
y1 = decay(ts,pfit1[0],pfit1[1])
id1 = argmin(abs(y1))
b1pA = ts[id1]
y2 = decay(ts,pfit2[0],pfit2[1])
id2 = argmin(abs(y2))
b1pB = ts[id2]
tB1p = ts
#################################################

fignum+=1
pylab.figure(fignum)

ax1 = pylab.subplot(2,1,1)
pylab.plot(ts,line,'r.',label='B1p 56MHz Q Demod',markersize=5)
pylab.plot(ts,decay(ts,pfit1[0],pfit1[1]),'k-')
pylab.plot(ts,decay(ts,pfit2[0],pfit2[1]),'k-')
pylab.grid(True)
pylab.xticks(visible=False)
#pylab.ylabel('DC Power [mW]',fontsize=10)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
pylab.ylim([0.0,2e-4])
pylab.xlim([2,198])
pylab.legend(fancybox=True,loc=2,prop={'size':8},numpoints=1)
ax1.get_yaxis().set_label_coords(-0.06,0.5)

pylab.title('WArm Schnupp Measurement')


###############################################
### Grab the B4 data
chan = 'V1:LSC_B4_56MHz_Q'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
ts = arange(0,duration,1./1000)
line = data_demod(ts,x.data,84.1,Fs,duration)

# Find the indices for each half of the data
idx1 = argmin(abs(ts-110.0))
idx2 = argmin(abs(ts-140.0))

# Fit each half
param_start = [1.0,1.0]
pfit1, fitCov = curve_fit(decay,ts[Fs:idx1],line[Fs:idx1],param_start)
perr = sqrt(diag(fitCov))
pfit2, fitCov = curve_fit(decay,ts[idx2:-0.3*Fs],line[idx2:-0.3*Fs],param_start)
perr = sqrt(diag(fitCov))

# Find the indices of the minima
y1 = decay(ts,pfit1[0],pfit1[1])
id1 = argmin(abs(y1))
b4A = ts[id1]
y2 = decay(ts,pfit2[0],pfit2[1])
id2 = argmin(abs(y2))
b4B = ts[id2]
tB4 = ts
#################################################


ax2 = pylab.subplot(2,1,2)

pylab.plot(ts,line,'r.',label='B4 56MHz Q Demod',markersize=5)
pylab.plot(ts,decay(ts,pfit1[0],pfit1[1]),'k-')
pylab.plot(ts,decay(ts,pfit2[0],pfit2[1]),'k-')
pylab.grid(True)
#pylab.xticks(visible=False)
#pylab.ylabel('DC Power [mW]',fontsize=10)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
pylab.ylim([0.0,2.3e-5])
pylab.legend(fancybox=True,loc=2,prop={'size':8},numpoints=1)
pylab.xlim([2,198])



###################################################

# Now get the phase data



print

chan = 'V1:SDB2_B1p_56MHz_phi0'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
t1 = arange(0,duration,1./Fs)
b1p_phi = x.data

chan = 'V1:SPRB_B4_56MHz_phi0'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
t1 = arange(0,duration,1./Fs)
b4_phi = x.data

yy = interp1d(t1,b1p_phi)
y1 = yy(b1pA)
y2 = yy(b1pB)
b1p_warm = ufloat((y1+y2)/2.,abs(y1-y2))
print 'WArm B1p:', b1p_warm

"""
yy = interp1d(t1,b4_phi)
y1 = yy(b4A)
y2 = yy(b4B)
b4_warm = ufloat((y1+y2)/2.,abs(y1-y2))
print 'WArm B4:', b4_warm



ax3 = pylab.subplot(3,1,3)

pylab.plot(t1,b1p_phi*180/pi + 71.0,'g--',label='B1p 56MHz Demod Phase')
pylab.plot(t1,b4_phi*180/pi + 28.0,'k--',label='B4 56MHz Demod Phase')
pylab.grid(True)
#pylab.xticks(visible=False)
#pylab.ylabel('Optical Gain [arb]',fontsize=10)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
pylab.ylim([-20.,20.])
#pylab.xlim([20,190])
pylab.legend(fancybox=True,loc=3,prop={'size':8},numpoints=1)
ax3.get_yaxis().set_label_coords(-0.06,0.5)
"""

pylab.xlabel('Time [sec]',fontsize=10)

pylab.savefig('Schnupp_warm.png')



print

dphi = abs(b1p_warm - b1p_narm)
print 'B1p DPHI:', dphi
fmod = 56436993.0

lam = constants.c / fmod
print 'B1p Schnupp:', dphi * lam / (4*pi)
print

"""
dphi = abs(b4_warm - b4_narm)
print 'B4 DPHI:', dphi
fmod = 56436993.0

lam = constants.c / fmod
print 'B4 Schnupp:', dphi * lam / (4*pi)
"""
