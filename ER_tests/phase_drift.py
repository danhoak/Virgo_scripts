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


# get the SSFS/MICH and PRCL phases for the science times
# only use segments that are >100sec in duration, and remove the last 10 seconds

PRCL_phi = array([])
SSFS_phi = array([])
PRCL_UGF = array([])
SSFS_UGF = array([])
DARM_UGF = array([])
MICH_UGF = array([])
B4_112mag = array([])
for i in range(len(seg_stop)):

    if duration[i] < 100:
        continue

    print i, seg_start[i], seg_stop[i], duration[i]

    chan = 'LSC_B2_8MHz_DPHI'
    x = getChannel('rds',chan,seg_start[i],duration[i]-10)
    PRCL_phi = append(PRCL_phi,x.data)

    chan = 'LSC_B4_56MHz_DPHI'
    x = getChannel('rds',chan,seg_start[i],duration[i]-10)
    SSFS_phi = append(SSFS_phi,x.data)

    chan = 'LSC_DARM_UGF'
    x = getChannel('rds',chan,seg_start[i],duration[i]-10)
    DARM_UGF = append(DARM_UGF,x.data)

    chan = 'LSC_SSFS_UGF'
    x = getChannel('rds',chan,seg_start[i],duration[i]-10)
    SSFS_UGF = append(SSFS_UGF,x.data)

    chan = 'LSC_MICH_UGF'
    x = getChannel('rds',chan,seg_start[i],duration[i]-10)
    MICH_UGF = append(MICH_UGF,x.data)

    chan = 'LSC_PRCL_UGF'
    x = getChannel('rds',chan,seg_start[i],duration[i]-10)
    PRCL_UGF = append(PRCL_UGF,x.data)

    chan = 'SPRB_B4_112MHz_mag_50Hz'
    x = getChannel('rds',chan,seg_start[i],duration[i]-10)
    B4_112mag = append(B4_112mag,x.data)


print median(SSFS_UGF)


fontP = FontProperties()
fontP.set_size('small')
matplotlib.rcParams.update({'savefig.dpi':250})

fignum=0
fignum=fignum+1
pylab.figure(fignum)

#bins = arange(-10,10,1)
n, bins, patches = pylab.hist(PRCL_phi*180/pi, 100, normed=1, histtype='bar', alpha=0.75)

pylab.xlabel('B2 8MHz DPHI [deg]', fontsize=10)
#pylab.xlim(-10,5)
#pylab.ylim(0,10)
pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)
pylab.grid(True, which='both', linestyle=':')
pylab.title('PRCL demod phase error in C8',fontsize=12)
pylab.savefig('PRCL_dphi_hist.png',bbox_inches='tight')


fignum=fignum+1
pylab.figure(fignum)

#bins = arange(-10,10,1)
n, bins, patches = pylab.hist(SSFS_phi*180/pi, 100, normed=1, histtype='bar', alpha=0.75)

pylab.xlabel('B4 56MHz DPHI [deg]', fontsize=10)
#pylab.xlim(-10,5)
#pylab.ylim(0,10)
pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)
pylab.grid(True, which='both', linestyle=':')
pylab.title('SSFS/MICH demod phase error in C8',fontsize=12)
pylab.savefig('SSFS_dphi_hist.png',bbox_inches='tight')



fignum=fignum+1
pylab.figure(fignum)

pylab.subplot(2,2,1)
n, bins, patches = pylab.hist(DARM_UGF, 100, normed=1, histtype='bar', alpha=0.75, facecolor='b',label='DARM')
pylab.xticks(fontsize=8)
pylab.yticks(fontsize=8)
pylab.grid(True, which='both', linestyle=':')
pylab.legend(loc='upper left',fancybox=True,prop={'size':8})

pylab.subplot(2,2,2)
n, bins, patches = pylab.hist(MICH_UGF, 100, normed=1, histtype='bar', alpha=0.75, facecolor='r',label='MICH')
pylab.xticks(fontsize=8)
pylab.yticks(fontsize=8)
pylab.grid(True, which='both', linestyle=':')
pylab.legend(loc='upper right',fancybox=True,prop={'size':8})

pylab.subplot(2,2,3)
n, bins, patches = pylab.hist(PRCL_UGF, 100, normed=1, histtype='bar', alpha=0.75, facecolor='g',label='PRCL')
pylab.xlabel('UGF [Hz]',fontsize=10)
pylab.xticks(fontsize=8)
pylab.yticks(fontsize=8)
pylab.grid(True, which='both', linestyle=':')
pylab.legend(loc='upper right',fancybox=True,prop={'size':8})

pylab.subplot(2,2,4)
n, bins, patches = pylab.hist(SSFS_UGF, 100, normed=1, histtype='bar', alpha=0.75, facecolor='k',label='SSFS')
pylab.xlabel('UGF [Hz]',fontsize=10)
pylab.xticks(fontsize=8)
pylab.yticks(fontsize=8)
pylab.grid(True, which='both', linestyle=':')
pylab.legend(loc='upper left',fancybox=True,prop={'size':8})

pylab.savefig('UGF_hist.png',bbox_inches='tight')



fignum=fignum+1
pylab.figure(fignum)

pylab.subplot(2,1,1)
n, bins, patches = pylab.hist(SSFS_UGF, 100, normed=1, histtype='bar', alpha=0.75, facecolor='k',label='SSFS UGF')
pylab.xlabel('UGF [Hz]',fontsize=10)
pylab.xticks(fontsize=8)
pylab.yticks(fontsize=8)
pylab.grid(True, which='both', linestyle=':')
pylab.legend(loc='upper left',fancybox=True,prop={'size':8})


pylab.subplot(2,1,2)
n, bins, patches = pylab.hist(B4_112mag, 100, normed=1, histtype='bar', alpha=0.75, facecolor='b',label='B4 112MHz Mag')

pylab.xlabel('Sideband power [mW]', fontsize=10)
pylab.xticks(fontsize=8)
pylab.yticks(fontsize=8)
pylab.grid(True, which='both', linestyle=':')
pylab.legend(loc='upper left',fancybox=True,prop={'size':8})

pylab.savefig('B4_112_hist.png',bbox_inches='tight')
