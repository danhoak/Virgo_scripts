#! /usr/bin/python

from numpy import *

import sys, os, time
sys.path.append('/users/swinkels/deploy/PythonVirgoTools/trunk/src')
from virgotools import getChannel

import matplotlib
matplotlib.use("Agg", warn=False)
import pylab
#from scipy import signal
#from scipy.optimize import leastsq
#from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
matplotlib.rcParams.update({'savefig.dpi':250})
from matplotlib.font_manager import FontProperties
from scipy import constants


# import segments

segments = genfromtxt('C8_segments.txt')
seg_start = segments[:,0]
seg_stop = segments[:,1]
duration = segments[:,2]


print len(seg_start)

# check if the ITF was unlocked in the second following the end of the segment

seg_stop_error = []
for i in range(len(seg_stop)):

    gps_start = seg_stop[i]-10
    duration = 20.0
    chan = 'LSC_ENABLE_50Hz'
    x = getChannel('rds',chan,gps_start,duration)
    Fs = x.fsample
    t = arange(0,duration,1./Fs)+gps_start

    if any(x.data==0):
        #print 'Segment was a lockloss!', seg_stop[i], sum(x.data)

        # get time of lockloss
        for j in range(len(x.data)):
            if x.data[j]==1 and x.data[j+1]<1:
                seg_stop_error.append(t[j+1] - seg_stop[i])
                print seg_stop[i], t[j+1], t[j+1] - seg_stop[i]
                continue
    else:
        print 'Segment was not a lockloss!', seg_stop[i], sum(x.data)


# histogram of unlock times

print
print seg_stop_error

fontP = FontProperties()
fontP.set_size('small')
matplotlib.rcParams.update({'savefig.dpi':250})

fignum=0
fignum=fignum+1
pylab.figure(fignum)

bins = arange(-10,10,1)
n, bins, patches = pylab.hist(seg_stop_error, bins, normed=0, histtype='bar', alpha=0.75)

pylab.xlabel('Time of lockloss relative to segment stop [sec]')
pylab.xlim(-10,5)
pylab.ylim(0,10)
pylab.grid(True, which='both', linestyle=':')
pylab.title('DQ_META_ITF_Mode is missing the end time of locks',fontsize=12)
pylab.savefig('seg_stop_error.png',bbox_inches='tight')
