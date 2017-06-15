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


################################################
#
# Find segments and save seglist
#
################################################

# C8 start and end times
gps_start = 1178028018 

# prior to 17:00 May 5 there were some times when DQ_META_ITF_Mode was not defined, start at 17:00
#gps_start = 1178038817

gps_stop = 1178262018
duration = gps_stop - gps_start

chan = 'DQ_META_ITF_Mode'
x = getChannel('rds',chan,gps_start,duration)
Fs = x.fsample

t = arange(0,duration,1./Fs)+gps_start

seg_start = []
seg_stop = []
seg_flag = False

for i in range(len(x.data)):

    # If a segment has just started
    if x.data[i]==1 and seg_flag==False:
        seg_start.append(t[i])
        seg_flag = True

    # If a segment has just finished
    elif x.data[i]<1 and seg_flag==True:
        seg_stop.append(t[i])
        seg_flag = False

    # If we are inside a segment
    elif x.data[i]==1 and seg_flag==True:
        continue

    # If we are outside a segment
    elif x.data[i]<1 and seg_flag==False:
        continue

    else:
        print 'Problem!', x.data[i], seg_flag, t[i]


# If we are still in a segment at the end of the data, close the final segment
if seg_flag==True:
        seg_stop.append(gps_stop)    



cumulative_time = 0
seg_durations = array(seg_stop)-array(seg_start)

#for i in range(len(seg_start)):
#    cumulative_time += seg_stop[i]-seg_start[i]
#    print i, seg_start[i], seg_stop[i], seg_stop[i]-seg_start[i]
    

print
print sum(seg_durations) / duration


output_data = zeros((len(seg_start),3))
for i in range(len(output_data)):
    output_data[i,0] = seg_start[i]
    output_data[i,1] = seg_stop[i]
    output_data[i,2] = seg_stop[i] - seg_start[i]

savetxt('C8_segments.txt',output_data,fmt='%.0f',delimiter='\t')



fontP = FontProperties()
fontP.set_size('small')
matplotlib.rcParams.update({'savefig.dpi':250})

fignum=0
fignum=fignum+1
pylab.figure(fignum)

dmin = 1
dmax = 55000

nbins = 16
logdmin=log10(dmin)
logdmax=log10(dmax)
dbinwidth = (logdmax-logdmin)/nbins
dbins=[0.0]*(nbins+1)
for i in range(len(dbins)):
   dbins[i] = pow(10,logdmin+i*dbinwidth)

n, b = histogram(seg_durations, dbins)

print n, b

plotbins = zeros((len(b)*2)) + 0.001
plotcounts = zeros((len(b)*2)) + 0.001

# Since we're plotting the logarithms, we need some non-zero values
# so that things will be defined.

for i in range(len(b)):
    plotbins[2*i] = b[i] - b[i]*0.0003
    plotbins[2*i+1] = b[i] + b[i]*0.0003

for i in range(len(n)-1):
    plotcounts[2*i+1] = n[i] + 0.001
    plotcounts[2*i+2] = n[i] + 0.001

print dbins
print plotbins, plotcounts


pylab.plot(plotbins,plotcounts,'k-',linewidth=1.0)

pylab.xscale('log')
pylab.yscale('linear')
pylab.xlabel('Duration [sec]',fontsize=14)
pylab.grid(True, which='major', linestyle=':')
pylab.grid(True, which='minor', linestyle=':')
pylab.xlim(dmin,dmax)
#pylab.ylim(0.7,1.2*max(plotcounts))
pylab.ylim(0.0,1.2*max(plotcounts))
pylab.title('C8 Segment Duration')
pylab.savefig('C8_segmentDurationHist.png',bbox_inches='tight')


