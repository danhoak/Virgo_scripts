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
    zz = bandpass(w,f0-5,f0+5,Fs0)

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

    return xdemod, t_dec



gps_start = 1178209917
duration = 300





fignum=0

fignum+=1
pylab.figure(fignum)


###############################################
### Grab the PD data and demodulate
chan = 'V1:SSFS_B4_Err_pre_50kHz'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
ts = arange(0,duration,1./Fs)
pre_line, t_dec = data_demod(ts,x.data,3345.0,Fs,duration)

chan = 'V1:SSFS_B4_Err_post_50kHz'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
ts = arange(0,duration,1./Fs)
post_line, t_dec = data_demod(ts,x.data,3345.0,Fs,duration)

oltf = pre_line / post_line * 3354.0

ssfs_ugf = decimate(oltf, int(round(10000/50.0)))

#################################################

fignum+=1
pylab.figure(fignum)

ax1 = pylab.subplot(2,1,1)
pylab.plot(t_dec,oltf,'r.',label='SSFS OLTF (pre/post)',markersize=5)
pylab.grid(True)
#pylab.xticks(visible=False)
pylab.ylabel('SSFS UGF [Hz]',fontsize=10)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
#pylab.ylim([0.0,1e-4])
#pylab.xlim([2,198])
pylab.legend(fancybox=True,loc=2,prop={'size':8},numpoints=1)
#ax1.get_yaxis().set_label_coords(-0.06,0.5)

#pylab.title('NArm Schnupp Measurement')


### Grab the B4 data
chan = 'V1:SPRB_B4_112MHz_mag_50Hz'
x = getChannel('rds',chan,gps_start,duration)
Fs = x.fsample
ts = arange(0,duration,1./Fs)

ax2 = pylab.subplot(2,1,2)

pylab.plot(ts,x.data,'b.',label='B4 56MHz Q Demod',markersize=5)
pylab.grid(True)
#pylab.xticks(visible=False)
pylab.ylabel('Sideband power [mW]',fontsize=10)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
#pylab.ylim([0.0,1.2e-5])
pylab.legend(fancybox=True,loc=2,prop={'size':8},numpoints=1)
#pylab.xlim([2,198])

pylab.xlabel('Time [sec]',fontsize=10)

pylab.savefig('SSFS_gain_demod.png',bbox_inches='tight')





# Now try MICH


###############################################
### Grab the PD data and demodulate
chan = 'V1:LSC_MICH_ERR'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
ts = arange(0,duration,1./Fs)
pre_line, t_dec = data_demod(ts,x.data,54.7,Fs,duration,lp=100.0,Fs0=1000)

chan = 'V1:LSC_MICH'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
ts = arange(0,duration,1./Fs)
post_line, t_dec = data_demod(ts,x.data,54.7,Fs,duration,lp=100.0,Fs0=1000)

oltf = pre_line / post_line * 54.7

mich_ugf = decimate(oltf, int(round(1000/50.0)))
t = arange(0,duration,1./50.0)

#################################################

fignum+=1
pylab.figure(fignum)

ax1 = pylab.subplot(2,1,1)
pylab.plot(t_dec,oltf,'r.',label='MICH OLTF (pre/post)',markersize=5)
#pylab.plot(t,mich_ugf,'k.',label='MICH OLTF (decimate)',markersize=5)
pylab.grid(True)
#pylab.xticks(visible=False)
pylab.ylabel('MICH UGF [Hz]',fontsize=10)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
#pylab.ylim([0.0,1e-4])
#pylab.xlim([2,198])
pylab.legend(fancybox=True,loc=2,prop={'size':8},numpoints=1)
#ax1.get_yaxis().set_label_coords(-0.06,0.5)

#pylab.title('NArm Schnupp Measurement')


### Grab the B4 data
chan = 'V1:SPRB_B4_112MHz_mag_50Hz'
x = getChannel('rds',chan,gps_start,duration)
Fs = x.fsample
ts = arange(0,duration,1./Fs)

B4_112 = x.data

ax2 = pylab.subplot(2,1,2)

pylab.plot(ts,x.data,'b.',label='B4 56MHz Q Demod',markersize=5)
pylab.grid(True)
#pylab.xticks(visible=False)
pylab.ylabel('Sideband power [mW]',fontsize=10)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
#pylab.ylim([0.0,1.2e-5])
pylab.legend(fancybox=True,loc=2,prop={'size':8},numpoints=1)
#pylab.xlim([2,198])

pylab.xlabel('Time [sec]',fontsize=10)

pylab.savefig('MICH_gain_demod.png',bbox_inches='tight')




fignum+=1
pylab.figure(fignum)

x = arange(0.003,0.007,0.0002)
y = 12500 * x - 46.

pylab.plot(B4_112,mich_ugf,'r.',markersize=5, alpha=0.3)
pylab.plot(x,y,'k--')
pylab.grid(True)
pylab.ylabel('MICH UGF [Hz]',fontsize=10)
pylab.xlabel('B4 112MHz mag [mw]',fontsize=10)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
pylab.ylim([0.0,40])
pylab.xlim([0.0035,0.0065])
#pylab.legend(fancybox=True,loc=2,prop={'size':8},numpoints=1)
#ax1.get_yaxis().set_label_coords(-0.06,0.5)

#pylab.title('NArm Schnupp Measurement')

pylab.savefig('MICH_gain_b4112.png',bbox_inches='tight')



fignum+=1
pylab.figure(fignum)

y = 1500000 * x - 100.

pylab.plot(B4_112,ssfs_ugf,'r.',markersize=5, alpha=0.3)
pylab.plot(x,y,'k--')
pylab.grid(True)
pylab.ylabel('SSFS UGF [Hz]',fontsize=10)
pylab.xlabel('B4 112MHz mag [mw]',fontsize=10)
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
pylab.ylim([4500,10500])
pylab.xlim([0.0035,0.0065])
#pylab.legend(fancybox=True,loc=2,prop={'size':8},numpoints=1)
#ax1.get_yaxis().set_label_coords(-0.06,0.5)

#pylab.title('NArm Schnupp Measurement')

pylab.savefig('SSFS_gain_b4112.png',bbox_inches='tight')
