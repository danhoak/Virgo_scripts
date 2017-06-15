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

#### Fitting equation

def decay(t,A0,Q):
    tau = Q / (pi*f0)
    return A0*exp(-t/tau)



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


def data_demod(t,xx,f0,Fs,duration,lp=6000.0,Fs0=10000.0):

    # lowpass
    z = lowpass(xx,lp,Fs)

    # decimate data
    w = decimate(z, int(round(Fs/Fs0)))
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

    return xdemod[i1:i2], t_dec[i1:i2]



start_time = 1178210597
duration = 5000


# line frequency
f0 = 1256.585
#f0 = 1875.37

err_chan = 'V1:LSC_DARM_ERR'

# grab data in steps, 50% overlap
step_duration = 50

t_DARM = arange(0,duration,step_duration/2)
datapoints = zeros(len(t_DARM))

demod_I = zeros(len(t_DARM))
demod_Q = zeros(len(t_DARM))

for i in range(len(t_DARM)):

    gps_start = start_time + i*step_duration/2
    x = getChannel('raw',err_chan,gps_start,step_duration)
    Fs = x.fsample
    ts = arange(0,duration,1./Fs)
    pre_line, t_dec = data_demod(ts,x.data,f0,Fs,step_duration,lp=2500.0,Fs0=10000.0)

    datapoints[i] = median(pre_line)
    print i, gps_start, datapoints[i]


fignum=0

fignum+=1
pylab.figure(fignum)

gs = gridspec.GridSpec(2,1,height_ratios=[3,1])

# fit data

t = t_DARM
data = datapoints
param_start = [3e-9,1e6]
pfit, fitCov = curve_fit(decay,t,data,param_start)

perr = sqrt(diag(fitCov))

print pfit
print perr

Q = round(pfit[1],-3)
Qerr = round(perr[1],-2)

tau = pfit[1] / (pi*f0)
tau_err = perr[1] / (pi*f0)
print Q, '+/-', Qerr
print tau, '+/-', tau_err

ax1 = pylab.subplot(gs[0])

#plot with a 1-second pad to avoid filter transients
pylab.plot(t,data,'bo',markersize=6)
pylab.plot(t,decay(t,pfit[0],pfit[1]),'k--',linewidth=1.5)
pylab.ylabel('Demodulated DARM_ERR [arb units]')
pylab.grid(True)
pylab.legend(('Demodulated data, 1256.585 Hz','Exponential Fit, Q = %.0f +/- %.0f'%(Q, Qerr)),fancybox=True,prop={'size':10})
#pylab.legend(('Demodulated data, 1875.37 Hz','Exponential Fit, Q = %.0f +/- %.0f'%(Q, Qerr)),fancybox=True,prop={'size':10})
#pylab.title('BS Drum Mode Ringdown, GPS = %.0f'%start_time, fontsize=12)
pylab.title('BS Butterfly Ringdown, GPS = %.0f'%start_time, fontsize=12)
ax1.get_yaxis().set_label_coords(-0.06,0.5)

ax2 = pylab.subplot(gs[1])

pylab.plot(t,data-decay(t,pfit[0],pfit[1]),'bo',markersize=6)

pylab.ylabel('fit residuals', fontsize=12)
pylab.xlabel('time [sec]', fontsize=12)
pylab.grid(True)
ax2.get_yaxis().set_label_coords(-0.06,0.5)

#pylab.savefig('BS_drum_fit.png')
pylab.savefig('BS_butter_fit.png')

