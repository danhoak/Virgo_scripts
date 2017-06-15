#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *

import sys, os, time
sys.path.append('/users/swinkels/deploy/PythonVirgoTools/trunk/src')
from virgotools import getChannel

import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab
from scipy import constants
from matplotlib.mlab import *
#from makewaveform import *
#from optimalSNR import *


fig_width_pt = 600  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
#golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio
golden_mean = 0.6
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size = [fig_width,fig_height]

matplotlib.rcParams.update({'savefig.dpi':350,
                            'text.usetex':True,
                            'figure.figsize':fig_size,
                            'font.family':"serif",
                            'font.serif':["Times"]})#,
#                            'xtick.major.pad':'8'})



def get_spectrum(gps_start, duration):
    chan = 'V1:Hrec_hoft_16384Hz'
    x = getChannel('raw',chan,gps_start,duration)
    Fs = x.fsample

    # Define parameters for FFT
    stride = 5.0   # FFT stride in seconds
    overlap = stride/2  # overlap in seconds (50%)
    
    Pxx, freq, t = specgram(x.data, NFFT=int(stride*Fs), Fs=int(Fs), noverlap=int(overlap*Fs))

    Axx = sqrt(Pxx)
    ASD = median(Axx,1)/sqrt(log(2))

    return ASD, freq



# Observed S6 H1 strain
#data1 = genfromtxt('lho4k_15May2010.txt')
data1 = genfromtxt('H1_S6_ASD.txt')
H1_freq_S6 = data1[:,0]
H1_S6 = data1[:,1]

# Observed S6 L1 strain
#data2 = genfromtxt('llo4k_31May2010.txt')
data2 = genfromtxt('L1_S6_ASD.txt')
L1_freq_S6 = data2[:,0]
L1_S6 = data2[:,1]

# Observed VSR4 strain
dataV = genfromtxt('VSR4_20110805.txt')
V1_freq = dataV[:,0]
V1_VSR4 = dataV[:,1]

"""
# aLIGO late
data5 = genfromtxt('ZERO_DET_high_P.txt')
late_freq = data5[:,0]
late_strain = data5[:,1]
"""

data = genfromtxt('3mm_e-3_advsens_steel.txt')
freq2 = data[:,0]
AdV_early_1e3 = data[:,1]


# Virgo early, 13W
data = genfromtxt('conf1.txt')
freq13W = data[0,:]
AdV_13W = data[1,:]

"""
data = genfromtxt('advsens_officialearly.txt')
freq3 = data[:,0]
AdV_early_fibers = data[:,1]


data = genfromtxt('adv_early_srtuned.txt')
freq4 = data[:,0]
AdV_early_fibers_SR = data[:,1]
"""

"""
# LIGO SRD
data3 = genfromtxt('SRD_strain_4k.txt')
SRD_freq = data3[:,0]
SRD_strain = data3[:,1]

# aLIGO early
data4 = genfromtxt('NO_SRM.txt')
early_freq = data4[:,0]
early_strain = data4[:,1]


# aLIGO H1
data9 = genfromtxt('CAL_DeltaL_7Jun.txt')
freq9 = data9[:,0]
H1_CAL = data9[:,1]/3994.5
H1_DARM = data9[:,2]

# aLIGO L1
data10 = genfromtxt('2015_06_06_L1_ER7.txt')
freq10 = data10[:,0]
L1_CAL = data10[:,1]
"""

"""
# H1 template waveform (from paper SVN)
# https://trac.ligo.caltech.edu/dac/browser/Papers/GW150914/matlab/timeseries/H1_template_1126259462_pyCBC.dat
tem_data = genfromtxt('H1_template_pyCBC.dat')
t = tem_data[:,0]
hp = tem_data[:,1]
hc = tem_data[:,2]

Fs = 16384.0
print len(t)/Fs
print 1/(t[1]-t[0])
print t[-1]-t[0]
CBC_PSD, freqCBC = matplotlib.mlab.psd(hp, NFFT=int(2*Fs), Fs=int(Fs))

# H1 best-fit SEOBNRv2ROM (from ben's page)
# https://ldas-jobs.ligo.caltech.edu/~bfarr/ER8/G184098/residuals/
tem_data = genfromtxt('H1_HOFT_SEOBNRv2ROM-8.dat')
tROM = tem_data[:,0]-1126259462.39
hROM = tem_data[:,1]

print len(tROM)/Fs
print 1/(tROM[2]-tROM[1])
SEOB_PSD, freqSEOB = matplotlib.mlab.psd(hROM, NFFT=int(16*Fs), Fs=int(Fs))

t_insp, hp_insp, hc_insp = makewaveform('inspiral','41~32~0~550',4.0,3.0,16384.0)

Fs = 1/(tROM[2]-tROM[1])
Fmax = 4000
print 'Sample rate for template', Fs


# aLIGO O1 medians
data9 = genfromtxt('H1_GW150914_ASD.txt')
freq9 = data9[:,0]
H1_O1 = data9[:,1]


# aLIGO O1 medians
data9 = genfromtxt('L1_GW150914_ASD.txt')
freq10 = data9[:,0]
L1_O1 = data9[:,1]

# aLIGO O1 medians
data = genfromtxt('O1_sensitivity.txt')
freq11 = data[:,0]
H1_O1 = data[:,1]
L1_O1 = data[:,2]

aligo_curves = genfromtxt('aligo_sensitivity.txt')
aLIGO_freq   = aligo_curves[:,0]
aLIGO_early  = aligo_curves[:,1]
aLIGO_mid    = aligo_curves[:,2]
aLIGO_late   = aligo_curves[:,3]
aLIGO_design = aligo_curves[:,4]
aLIGO_BNS    = aligo_curves[:,5]

avirgo_curves = genfromtxt('avirgo_sensitivity.txt')
aVIRGO_freq   = avirgo_curves[:,0]
aVIRGO_early1 = avirgo_curves[:,1]
aVIRGO_early2 = avirgo_curves[:,2]
aVIRGO_mid    = avirgo_curves[:,3]
aVIRGO_late1  = avirgo_curves[:,4]
aVIRGO_late2  = avirgo_curves[:,5]
aVIRGO_design = avirgo_curves[:,6]
aVIRGO_BNS    = avirgo_curves[:,7]
"""


data = genfromtxt('L1_O2_Sensitivity_strain_asd.txt')
freq_L1O2 = data[:,0]
L1_O2 = data[:,1]

data = genfromtxt('H1_O2_Sensitivity_strain_asd.txt')
freq_H1O2 = data[:,0]
H1_O2 = data[:,1]



# AdVirgo
data = genfromtxt('sensitivity.txt')
#data = genfromtxt('sensitivity_May3_1945.txt')
freq_NB = data[:,0]
AdV_NB = data[:,1]

fontP = FontProperties()
fontP.set_size('small')



fignum=0


fignum+=1
pylab.figure(fignum)

pylab.loglog(freq2,AdV_early_1e3,'k',linestyle='--',linewidth=1.6,label='AdVirgo, 25W, Steel SUS, no SR',alpha=1.0)
#pylab.loglog(freq13W,AdV_13W,'k',linestyle='--',linewidth=1.6,label='AdVirgo, 13W, Steel SUS, no SR',alpha=1.0)
#pylab.loglog(H1_freq_S6,H1_S6,'Fuchsia',linestyle='-',linewidth=0.8,label='H1 in S6 (May 15 2010)',alpha=0.7)
#pylab.loglog(L1_freq_S6,L1_S6,'DarkTurquoise',linestyle='-',linewidth=0.8,label='L1 in S6 (May 31 2010)',alpha=0.7)
pylab.loglog(V1_freq,V1_VSR4,'LimeGreen',linestyle='-',linewidth=0.7,label='VSR4 (5 Aug 2011)',alpha=0.7)
pylab.loglog(freq_H1O2,H1_O2,'r-',linewidth=0.7,label='H1 in O2',alpha=0.7)
pylab.loglog(freq_L1O2,L1_O2,'b-',linewidth=0.7,label='L1 in O2',alpha=0.7)


# Old noise budget data - calibration was wrong
"""
data = genfromtxt('sensitivity_May3_1945.txt')
freq_NB = data[:,0]
AdV_NB = data[:,1]
pylab.loglog(freq_NB,AdV_NB,linestyle='-',linewidth=1.0,label='3 May 2017',alpha=1.0)

#data = genfromtxt('sensitivity_May5_0108.txt')
data = genfromtxt('sensitivity_May5_1835.txt')
freq_NB = data[:,0]
AdV_NB = data[:,1]
pylab.loglog(freq_NB,AdV_NB,linestyle='-',linewidth=1.0,label='5 May 2017 - MICH subtraction',alpha=1.0)

data = genfromtxt('sensitivity.txt')
freq_NB = data[:,0]
AdV_NB = data[:,1]
pylab.loglog(freq_NB,AdV_NB,linestyle='-',linewidth=1.0,label='12 May 2017 - SDB1 LP, OMC2, more MICH',alpha=1.0)
"""

#### Grab Hrec data

#gps_start = 1179709240
#gps_start = 1179865427

gps_start = 1178655078
duration = 120
ASD, freq = get_spectrum(gps_start, duration)
pylab.loglog(freq,ASD,'DarkSlateBlue',linestyle='-',linewidth=0.8,label='V1 12 May 2017 (5.5 Mpc)',alpha=1.0)


gps_start = 1180914137
duration = 120
ASD, freq = get_spectrum(gps_start, duration)
pylab.loglog(freq,ASD,'DarkMagenta',linestyle='-',linewidth=0.8,label='V1 7 Jun 2017 (8 Mpc)',alpha=1.0)

gps_start = 1181231717
duration = 120
ASD, freq = get_spectrum(gps_start, duration)
pylab.loglog(freq,ASD,'Fuchsia',linestyle='-',linewidth=0.8,label='V1 11 Jun 2017 (9 Mpc)',alpha=1.0)

pylab.grid(True, which='both', linestyle=':',alpha=0.8)
pylab.ylim(3e-24,1e-18)
pylab.xlim(10,3000)
pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)
pylab.xlabel('Frequency [Hz]',fontsize=12)
pylab.ylabel('Strain ASD 1/rt[Hz]',fontsize=12)
pylab.legend(loc='upper right',fancybox=True,prop={'size':7})
#pylab.legend(bbox_to_anchor=(0.82,1.1),fancybox=True,prop={'size':9},ncol=2)
pylab.savefig('AdVirgo_progress.png',bbox_inches='tight')




"""
fignum+=1
pylab.figure(fignum)

#pylab.loglog(H1_freq_S6,H1_S6,'Fuchsia',linestyle='-',linewidth=0.8,label='H1 in S6 (May 15 2010)',alpha=0.7)
#pylab.loglog(L1_freq_S6,L1_S6,'DarkTurquoise',linestyle='-',linewidth=0.8,label='L1 in S6 (May 31 2010)',alpha=0.7)
pylab.loglog(V1_freq,V1_VSR4,'LimeGreen',linestyle='-',linewidth=0.8,label='V1 in VSR4 (Aug 5 2011)',alpha=0.7)

pylab.loglog(freq_H1O2,H1_O2,'r-',linewidth=0.7,label='H1 in O2',alpha=0.7)
pylab.loglog(freq_L1O2,L1_O2,'b-',linewidth=0.7,label='L1 in O2',alpha=0.7)

pylab.loglog(freq2,AdV_early_1e3,'k',linestyle='--',linewidth=1.6,label='AdVirgo, 25W, Steel SUS, no SR',alpha=1.0)
pylab.loglog(freq_NB,AdV_NB,'g',linestyle='-',linewidth=1.4,label='AdVirgo Preliminary',alpha=1.0)

pylab.grid(True, which='both', linestyle=':',alpha=0.8)
pylab.ylim(3e-24,3e-19)
pylab.xlim(10,3000)
pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)
pylab.xlabel('Frequency [Hz]',fontsize=12)
pylab.ylabel('Strain ASD 1/rt[Hz]',fontsize=12)
pylab.legend(bbox_to_anchor=(0.82,1.1),fancybox=True,prop={'size':9},ncol=2)

pylab.savefig('AdVirgo_strain_prelim.png',bbox_inches='tight')
"""




"""
fignum+=1
pylab.figure(fignum)

pylab.loglog(freq2,AdV_early_1e3,'k',linestyle='--',linewidth=1.6,label='AdVirgo, 25W, Steel SUS, no SR',alpha=1.0)
#pylab.loglog(freq13W,AdV_13W,'r',linestyle='--',linewidth=1.6,label='AdVirgo, 13W, Steel SUS, no SR',alpha=1.0)
#pylab.loglog(H1_freq_S6,H1_S6,'Fuchsia',linestyle='-',linewidth=0.8,label='H1 in S6 (May 15 2010)',alpha=0.7)
#pylab.loglog(L1_freq_S6,L1_S6,'DarkTurquoise',linestyle='-',linewidth=0.8,label='L1 in S6 (May 31 2010)',alpha=0.7)
#pylab.loglog(V1_freq,V1_VSR4,'LimeGreen',linestyle='-',linewidth=0.8,label='VSR4 (5 Aug 2011)',alpha=0.7)
#pylab.loglog(freq_H1O2,H1_O2,'r-',linewidth=0.7,label='H1 in O2',alpha=0.7)
#pylab.loglog(freq_L1O2,L1_O2,'b-',linewidth=0.7,label='L1 in O2',alpha=0.7)


data = genfromtxt('sensitivity.txt')
freq_NB = data[:,0]
AdV_NB = data[:,1]
pylab.loglog(freq_NB,AdV_NB,linestyle='-',linewidth=1.0,label='Noise budget - 25 May 2017',alpha=1.0)

pylab.loglog(V1freq,V1_ASD_25May,linestyle='-',linewidth=1.0,label='Hrec h(t) - 25 May 2017',alpha=1.0)


pylab.grid(True, which='both', linestyle=':',alpha=0.8)

pylab.ylim(5e-24,1e-18)
pylab.xlim(10,3000)

pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)

#pylab.title('Previous and Expected LIGO Sensitivies')
pylab.xlabel('Frequency [Hz]',fontsize=12)
pylab.ylabel('Strain ASD 1/rt[Hz]',fontsize=12)

pylab.legend(loc='upper right',fancybox=True,prop={'size':8})
#pylab.legend(bbox_to_anchor=(0.82,1.1),fancybox=True,prop={'size':9},ncol=2)

pylab.savefig('AdVirgo_cali_compare.png',bbox_inches='tight')
"""
