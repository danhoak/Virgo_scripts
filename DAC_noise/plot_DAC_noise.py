#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys, os, time
sys.path.append('/users/swinkels/deploy/PythonVirgoTools/trunk/src')
from virgotools import getChannel


from numpy import *
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab
from scipy import constants
from scipy.special import jn
import nbutils as nbu
from matplotlib.mlab import *
#from makewaveform import *
#from optimalSNR import *

fig_width_pt = 600  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
#golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio
golden_mean = 0.8
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size = [fig_width,fig_height]

matplotlib.rcParams.update({'savefig.dpi':350,
                            'text.usetex':True,
                            'figure.figsize':fig_size,
                            'font.family':"serif",
                            'font.serif':["Times"]})#,
#                            'xtick.major.pad':'8'})


# Define parameters for FFT
stride = 6.0   # FFT stride in seconds
overlap = 3.0  # overlap in seconds (50%)



#### Grab Hrec data
gps_start = 1180111697
duration = 100

chan = 'V1:Hrec_hoft_16384Hz'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
    
# for science segment, get PSD.
# 6-second FFTs, 50% overlap, hann window (default)
Pxx, V1freq, t = specgram(x.data, NFFT=int(stride*Fs), Fs=int(Fs), noverlap=int(overlap*Fs))

Axx = sqrt(Pxx)
V1_ASD_NIWI = median(Axx,1)/sqrt(log(2))




#### Grab Hrec data
gps_start = 1179865427
duration = 100

chan = 'V1:Hrec_hoft_16384Hz'
x = getChannel('raw',chan,gps_start,duration)
Fs = x.fsample
    
# for science segment, get PSD.
# 6-second FFTs, 50% overlap, hann window (default)
Pxx, V1freq, t = specgram(x.data, NFFT=int(stride*Fs), Fs=int(Fs), noverlap=int(overlap*Fs))

Axx = sqrt(Pxx)
V1_ASD_26May = median(Axx,1)/sqrt(log(2))




"""
# aLIGO O1 medians
data = genfromtxt('O1_sensitivity.txt')
freq11 = data[:,0]
H1_O1 = data[:,1]
L1_O1 = data[:,2]
"""

"""
# aLIGO O2 medians
data = genfromtxt('L1_O2_Sensitivity_strain_asd.txt')
freq11 = data[:,0]
L1_O1 = data[:,1]

# aLIGO O2 medians
data = genfromtxt('H1_O2_Sensitivity_strain_asd.txt')
freq11 = data[:,0]
H1_O1 = data[:,1]
"""

#data = genfromtxt('3mm_e-4_advsens_steel.txt')
#freq1 = data[:,0]
#AdV_early_1e4 = data[:,1]

data = genfromtxt('3mm_e-3_advsens_steel.txt')
freq2 = data[:,0]
AdV_early_1e3 = data[:,1]

data = genfromtxt('advsens_officialearly.txt')
freq3 = data[:,0]
AdV_early_fibers = data[:,1]

#data = genfromtxt('adv_early_srtuned.txt')
#freq4 = data[:,0]
#AdV_early_fibers_SR = data[:,1]


#data = genfromtxt('DARM_Mar28.txt')
#freq0328 = data[:,0]
#DARM_Mar28 = data[:,1]


avirgo_curves = genfromtxt('avirgo_sensitivity.txt')
aVIRGO_freq   = avirgo_curves[:,0]
aVIRGO_early1 = avirgo_curves[:,1]
aVIRGO_early2 = avirgo_curves[:,2]
aVIRGO_mid    = avirgo_curves[:,3]
aVIRGO_late1  = avirgo_curves[:,4]
aVIRGO_late2  = avirgo_curves[:,5]
aVIRGO_design = avirgo_curves[:,6]
aVIRGO_BNS    = avirgo_curves[:,7]


# Calculate aLIGO SUS actuator noise

esd_bias = 380 # V
esd_coeff = 2e-10 # N/V**2
#ff = logspace(-1, 5, 300)
ff = freq2

ligo_dac_noise = lambda ff: 300e-9*sqrt(1+(50./ff)**2)

def get_quad_tst_act(ff):
    zeros = [50, 50, 3250]
    poles = [2.2, 2.2, 152]
    gain = 2 * (2*esd_coeff*esd_bias)
    _, driver_tf = nbu.freqresp((zeros, poles, gain), ff)
    sus_resp = nbu.get_complex_interp('quadTstTstLL.txt', ff)
    return ff, abs(driver_tf * sus_resp) * ligo_dac_noise(ff)

def get_quad_bias_act(ff):
    zeros = []
    poles = [1.8, 1.8, 1590]
    gain = 40 * (2*esd_coeff*esd_bias)
    _, driver_tf = nbu.freqresp((zeros, poles, gain), ff)
    sus_resp = nbu.get_complex_interp('quadTstTstLL.txt', ff)
    return ff, abs(driver_tf * sus_resp) * ligo_dac_noise(ff)

def get_quad_pum_act(ff):
    # T1100378
    zeros = [6, 20, 1.35]
    poles = [0.5, 250, 110]
    gain = 0.268e-3*0.0309
    _, driver_tf = nbu.freqresp((zeros, poles, gain), ff)
    sus_resp = nbu.get_complex_interp('quadPumTstLL.txt', ff)
    return ff, 2 * abs(driver_tf * sus_resp) * ligo_dac_noise(ff)


def get_bs_pum_act(ff):
    # T1100479
    zeros = [10, 20, 10]
    poles = [1, 200, 100]
    gain = 0.32e-3*0.963
    _, driver_tf = nbu.freqresp((zeros, poles, gain), ff)
    sus_resp = nbu.get_complex_interp('bsfmM2M3LL.txt', ff)
    return ff, abs(driver_tf * sus_resp) * ligo_dac_noise(ff) / 280


ff, quadTSTnoise = get_quad_tst_act(ff)
ff, quadBiasnoise = get_quad_bias_act(ff)
ff, quadPUMnoise = get_quad_pum_act(ff)
ff, bsPUMnoise = get_bs_pum_act(ff)

LIGO_SUSnoise = sqrt(quadTSTnoise**2 + quadBiasnoise**2 + quadPUMnoise**2 + bsPUMnoise**2)


# Calculate AdV SUS actuator noise

virgo_dac_noise_lf = lambda ff: 800e-9*sqrt(1./ff)
virgo_dac_noise_hf = lambda ff: 100e-9

virgo_dac_noise = sqrt(virgo_dac_noise_lf(ff)**2 + virgo_dac_noise_hf(ff)**2)

# coil driver transfer function is flat out to 100kHz

# 10 ohm for resistors, 20 ohm for cable + coil (see Vincenzo's comment)

def get_TM_Mir_act(ff, res=10):
    zeros = []
    poles = [300]
    # 3.4mN/A per pair, 20 ohms (resistor+coil), noise adds incoherently
    gain = sqrt(2)*3.4e-3 / (res+20)
    _, driver_tf = nbu.freqresp((zeros, poles, gain), ff)
    sus_resp = nbu.get_complex_interp('TM_Mir_resp.txt', ff)
    return ff, abs(driver_tf * sus_resp) * virgo_dac_noise

def get_TM_Mar_act(ff):
    zeros = []
    poles = [300]
    # 40mN/A per pair, 30 ohms (resistor+coil), noise adds incoherently
    gain = sqrt(2)*40.0e-3 / 30
    _, driver_tf = nbu.freqresp((zeros, poles, gain), ff)
    sus_resp = nbu.get_complex_interp('TM_Mar_resp.txt', ff)
    return ff, abs(driver_tf * sus_resp) * virgo_dac_noise

def get_BS_Mir_act(ff):
    zeros = []
    poles = [300]
    # 1.7mN/A per pair, 30 ohms (resistor+coil), noise adds incoherently
    #gain = sqrt(2)*1.7e-3 / 30
    # factor of 6 increase for mirror from combination of magnets, coils, DAC actuation changes
    gain = sqrt(2)*1.7e-3*6 / 30
    _, driver_tf = nbu.freqresp((zeros, poles, gain), ff)
    sus_resp = nbu.get_complex_interp('BS_Mir_resp.txt', ff)
    return ff, abs(driver_tf * sus_resp) * virgo_dac_noise / 200

def get_BS_Mar_act(ff):
    zeros = []
    poles = [300]
    # 28mN/A per pair, 30 ohms (resistor+coil), noise adds incoherently
    gain = sqrt(2)*28e-3 / 30
    _, driver_tf = nbu.freqresp((zeros, poles, gain), ff)
    sus_resp = nbu.get_complex_interp('BS_Mar_resp.txt', ff)
    return ff, abs(driver_tf * sus_resp) * virgo_dac_noise / 200


ff, TM_Mir_noise = get_TM_Mir_act(ff,10)
ff, TM_Mar_noise = get_TM_Mar_act(ff)
ff, BS_Mir_noise = get_BS_Mir_act(ff)
ff, BS_Mar_noise = get_BS_Mar_act(ff)

ff, TM_Mir_noise_1k = get_TM_Mir_act(ff,1000)
ff, TM_Mir_noise_10k = get_TM_Mir_act(ff,10000)
ff, TM_Mir_noise_80k = get_TM_Mir_act(ff,80000)

Virgo_MARnoise = sqrt(4*TM_Mar_noise**2)
#Virgo_MIRnoise = sqrt(4*TM_Mir_noise**2)
Virgo_MIRnoise = sqrt(2*TM_Mir_noise**2)
Virgo_BSnoise = sqrt(BS_Mir_noise**2 + BS_Mar_noise**2)
Virgo_BS_MIRnoise = BS_Mir_noise

Virgo_MIRnoise_1k = sqrt(2*TM_Mir_noise_1k**2)
Virgo_MIRnoise_10k = sqrt(2*TM_Mir_noise_10k**2)
Virgo_MIRnoise_80k = sqrt(2*TM_Mir_noise_80k**2)



# estimate shot noise limited sensitivity for PRFPMI configuration

lam = 1064e-9
Pin = 25
g_cr2 = 35
finesse = 450.0
L = 2999.8
fc = 55.0
T_IFO = 0.878  # nominal=0.878

# calculate circulating power in the arms - input power * PRC gain * BS * arm gain * transmission into the IFO
armP = Pin * g_cr2 * 0.5 * 280 * T_IFO
#print armP

"""
detection bench losses:

faraday = 4%
1.5% to B1p
3% OMC losses from internal losses and alignment mismatch

current DCPDs are same as others - 88% QE - but can be replaced to 98%


ep = T_IFO * 0.86 * 0.96 * 0.98 # T_IFO * SDB1 losses * OMC mode matching * DCPD response

print 1/sqrt(ep)
h_shot = sqrt(1/ep)/(4*sqrt(g_cr2)*finesse*L) * sqrt(pi*constants.hbar*constants.c*lam/Pin) * abs(1+1j*freq2/fc)


# shot noise strain sensitivity for FPMI with DC readout

T_PR = 0.048
h_shot_FPMI = sqrt(1/ep)/(4*finesse*L) * sqrt(pi*constants.hbar*constants.c*lam/(T_PR * Pin)) * abs(1+1j*aVIRGO_freq/fc)


# calculation of displacement shot noise for RF readout of FPMI

# modulation depths, etc
J0 = jn(0,0.22)
J1 = jn(1,0.22)
gcr = 1

N = 8 * J0 * J1 * gcr
P_in = Pin * T_PR
P_DC = 0.18e-3

#d_shot_FPMI_RF = sqrt(1/0.015)/(N * P_in * finesse) * sqrt(constants.h*constants.c*lam*P_DC) * abs(1+1j*aVIRGO_freq/fc)

d_shot_FPMI_RF = h_shot_FPMI * L * sqrt(2)



AdV_early_losses = zeros(len(AdV_early_1e3))
for i in range(len(AdV_early_1e3)):
    AdV_early_losses[i] = max(AdV_early_1e3[i],h_shot[i])

"""


fignum=0

fontP = FontProperties()
fontP.set_size('small')
matplotlib.rcParams.update({'savefig.dpi':250})


fignum+=1
pylab.figure(fignum)

pylab.loglog(freq2,AdV_early_1e3,'k-',linewidth=1.2,label='AdV Early, Steel, 25W (VIR-0316D-16)',alpha=0.6)
#pylab.loglog(aVIRGO_freq,aVIRGO_design,'k--',linewidth=1.0,label='AdV Design (130 Mpc)')
pylab.loglog(freq3,AdV_early_fibers,'k--',linewidth=1.0,label='AdV Early, Fibers (VIR-0316D-16)')

#pylab.loglog(ff,LIGO_SUSnoise/3e3,'k-',linewidth=1.2,label='SUS Actuator Noise',alpha=0.6)
pylab.loglog(ff,Virgo_MIRnoise/2999.8,'b-',linewidth=2.0,label='DAC Noise Estimate, NI+WI',alpha=0.6)
#pylab.loglog(ff,Virgo_MIRnoise_1k/2999.8,'Teal',linestyle='-',linewidth=2.0,label='Test Mass DAC Noise, 1kOhm',alpha=0.6)
#pylab.loglog(ff,Virgo_MIRnoise_10k/2999.8,'Fuchsia',linestyle='-',linewidth=2.0,label='Test Mass DAC Noise, 10kOhm',alpha=0.6)
#pylab.loglog(ff,Virgo_MIRnoise_80k/2999.8,'r-',linewidth=2.0,label='Test Mass DAC Noise, 80kOhm',alpha=0.6)

#pylab.loglog(ff,Virgo_MARnoise/2999.8,'r-',linewidth=2.0,label='Test Mass DAC Noise (marionette)',alpha=0.6)
#pylab.loglog(ff,Virgo_BSnoise/2999.8,'g-',linewidth=2.0,label='BS DAC Noise (mirror+marionette)',alpha=0.6)
pylab.loglog(ff,Virgo_BS_MIRnoise/2999.8,'g-',linewidth=2.0,label='BS MIR DAC Noise Estimate (new magnets)',alpha=0.8)

#pylab.loglog(ff,(sqrt(AdV_early_1e3**2 + (Virgo_BS_MIRnoise/2999.8)**2)),'r--',linewidth=2.0,label='Steel + BS MIR DAC Noise',alpha=0.8)

pylab.loglog(V1freq,V1_ASD_26May,'DarkTurquoise',linewidth=1.0,label='V1 26 May 2017',alpha=0.8)
pylab.loglog(V1freq,V1_ASD_NIWI,'r',linewidth=1.0,label='V1 29 May, NI+WI MIR DAC enabled',alpha=0.8)

pylab.grid(True, which='both', linestyle=':', alpha=0.8)
pylab.ylim(1e-24,1e-17)
pylab.xlim(8,3000)
pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)
pylab.xlabel('Frequency [Hz]',fontsize=12)
pylab.ylabel('Strain ASD 1/rt[Hz]',fontsize=12)
pylab.legend(loc=1,fancybox=True,prop={'size':8})
pylab.savefig('AdV_DAC_noise.png',bbox_inches='tight')


