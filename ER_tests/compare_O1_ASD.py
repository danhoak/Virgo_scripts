#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A variation of the Rayleighgram script to handle long spans of data which may cover multiple science segments
# Data flagged as DMT-ANALYSIS_READY is used and CAT1 veto times are removed.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
from MCIV_tools import *
import matplotlib
matplotlib.use("Agg")
import pylab
import os, sys
from glue import segmentsUtils
import h5py

fig_width_pt = 600  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
#golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio
golden_mean = 0.6
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size = [fig_width,fig_height]

matplotlib.rcParams.update({'savefig.dpi':250,
                            'text.usetex':True,
                            'figure.figsize':fig_size,
                            'font.family':"serif",
                            'font.serif':["Times"],
                            'xtick.major.pad':'8'})



from matplotlib.mlab import *
from scipy import stats
from glue.segments import *


fmin = 8
fmax = 2050

res = '6sec'

# grab H1 data, L1 data

datafile = 'H1_segs/H1_O1_percentiles_' + res + '.hdf5'
H = h5py.File(datafile,'r')

dataset = H['p50']
H_p50 = dataset[...]

dataset = H['p95']
H_p95 = dataset[...]

dataset = H['p05']
H_p05 = dataset[...]

dataset = H['p997']
H_p997 = dataset[...]

dataset = H['p003']
H_p003 = dataset[...]

dataset = H['freq']
freq = dataset[...]

H.close()


datafile = 'L1_segs/L1_O1_percentiles_' + res + '.hdf5'
L = h5py.File(datafile,'r')

dataset = L['p50']
L_p50 = dataset[...]

dataset = L['p95']
L_p95 = dataset[...]

dataset = L['p05']
L_p05 = dataset[...]

dataset = L['p997']
L_p997 = dataset[...]

dataset = L['p003']
L_p003 = dataset[...]

L.close()


# Get percentiles for standard Rayleigh distribution with 1e6 samples
x = random.rayleigh(1.0,1000000)
xp50  = stats.scoreatpercentile(x,50)
xp95  = stats.scoreatpercentile(x,95.5)
xp997 = stats.scoreatpercentile(x,99.73)
xp05  = stats.scoreatpercentile(x,4.5)
xp003 = stats.scoreatpercentile(x,0.27)

xmp95  = xp95/xp50
xmp997 = xp997/xp50
xmp05  = xp05/xp50
xmp003 = xp003/xp50

fignum=0
fignum=fignum+1
pylab.figure(fignum)

ax = pylab.subplot(1,1,1)

pylab.loglog(freq,H_p50/sqrt(log(2)),'Red',linestyle="-",linewidth=1.0,label=r' H1 median noise')
pylab.loglog(freq,L_p50/sqrt(log(2)),'Blue',linestyle="-",linewidth=1.0,label=r' L1 median noise')

pylab.loglog(freq,H_p003/sqrt(log(2))/xmp003,'Fuchsia',linestyle="--",linewidth=1.0,label=r' H1 projected noise')
pylab.loglog(freq,L_p003/sqrt(log(2))/xmp003,'DarkTurquoise',linestyle="--",linewidth=1.0,label=r' L1 projected noise')

pylab.axis([fmin, fmax, 3e-24, 1e-19])
pylab.grid(True, which='both', linestyle=':')

pylab.xlabel('Frequency [Hz]',fontsize=12)
pylab.ylabel('Strain ASD [1/rt(Hz)]',fontsize=12)
pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)
#ax.get_yaxis().set_label_coords(-0.105,0.5)
pylab.legend(loc=1,prop={'size':10},fancybox=True)
pylab.title('Noise in O1 - Oct 15 to Dec 8 2015', fontsize=12)

pylab.savefig('O1_noise_' + res + '.png',bbox_inches='tight')
pylab.close()




fignum=fignum+1
pylab.figure(fignum)

ax = pylab.subplot(1,1,1)

pylab.semilogx(freq,H_p95/H_p50,'Fuchsia',linestyle="-",linewidth=1.4,label=r'H1 95$^{\mathrm{th}}$ Percentile')
pylab.semilogx(freq,H_p997/H_p50,'Firebrick',linestyle="-",linewidth=1.4,label=r'H1 99.7$^{\mathrm{th}}$ Percentile')

pylab.semilogx(freq,L_p95/L_p50,'DarkTurquoise',linestyle="-",linewidth=1.4,label=r'L1 95$^{\mathrm{th}}$ Percentile')
pylab.semilogx(freq,L_p997/L_p50,'steelblue',linestyle="-",linewidth=1.4,label=r'L1 99.7$^{\mathrm{th}}$ Percentile')

pylab.semilogx([fmin, fmax],[xmp95,xmp95],"y--",linewidth=1.8,label=r'Gaussian Noise: 95$^{\mathrm{th}}$ Percentile')
pylab.semilogx([fmin, fmax],[xmp997,xmp997],"g--",linewidth=1.8,label=r'Gaussian Noise: 99.7$^{\mathrm{th}}$ Percentile')

pylab.axis([fmin, fmax, 1, 10.0])
pylab.grid(True, which='both', linestyle=':')

pylab.xlabel('Frequency [Hz]',fontsize=12)
pylab.ylabel('Strain ASD, Normalized to Median',fontsize=12)
pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)
#ax.get_yaxis().set_label_coords(-0.08,0.5)
pylab.legend(loc=1,prop={'size':10},fancybox=True)

pylab.savefig('ASD_deviation_O1_' + res + '.png',bbox_inches='tight')
pylab.close()
