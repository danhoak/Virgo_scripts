#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Estimate optical cavity parameters of the Advanced Virgo interferometer (no signal recycling)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
from scipy.special import jn
#import matplotlib
#from matplotlib.font_manager import FontProperties
#matplotlib.use("Agg")
#import pylab

def cavity_params(r_i,r_e,RTL,FSR):

    finesse = pi*sqrt(r_i*r_e)/(1-r_i*r_e)

    cavity_pole = FSR/(2*finesse)

    cavity_gain = t_i/(1-r_i*r_e*sqrt(1-RTL))
    cavity_reflectivity = (r_e*sqrt(1-RTL)-r_i)/(1-r_i*r_e*sqrt(1-RTL))
    cavity_transmissivity = t_i*t_e/(1-r_i*r_e*sqrt(1-RTL))
    storage_time = -1/(FSR*log(r_i*r_e*sqrt(1-RTL)))

    return finesse, cavity_pole, cavity_gain, cavity_reflectivity, cavity_transmissivity, storage_time



# Note: all as-built numbers are for H1

c = 299792458.0
lamd = 1064.0e-9

#### Input mode cleaner (design specs)

r1 = sqrt(1-2500e-6)
r2 = sqrt(1-5.1e-6)
r3 = sqrt(1-2500e-6)

t1 = sqrt(1-r1**2)
t2 = sqrt(1-r2**2)
t3 = sqrt(1-r3**2)

L = 143.4259
FSR = c / (2*L)

IMC_FSR = FSR

finesse = pi*sqrt(r1*r2*r3)/(1-r1*r2*r3)

cavity_pole = FSR/(2*finesse)

cavity_gain = t1/(1-r1*r2*r3)
cavity_reflectivity = (r3-r1)/(1-r1*r2*r3)
cavity_transmittivity = t1*t3/(1-r1*r2*r3)
storage_time = -1/(FSR*log(r1*r2))

print
print 'Input Mode Cleaner'
print
print 'Intra-cavity gain (power):', cavity_gain**2
print 'Cavity reflection (power):', cavity_reflectivity**2
print 'Cavity transmission (power):', cavity_transmittivity**2
print
print 'FSR [Hz]:', FSR
print 'Finesse:', finesse
print 'Cavity pole [Hz]:', cavity_pole
print 'Storage time [sec]:', storage_time
print 2*finesse*L/(pi*c), 1/(4*finesse*L/(c))
print

# X-arm

L = 2999.8
FSR = c / (2*L)

# Test mass HR surface loss
# Table 1 page 8 of TDR
L_TM = 37.5e-6
#RTL = 2*L_TM
RTL = 44e-6

t_i = sqrt(0.01377)
t_e = sqrt(4.3e-6)

#r_i = sqrt(1 - t_i**2 - L_TM)
r_i = sqrt(1 - t_i**2)
r_e = sqrt(1 - t_e**2)



xfinesse, xcavity_pole, xcavity_gain, xcavity_reflectivity, xcavity_transmissivity, xstorage_time = cavity_params(r_i,r_e,RTL,FSR)

print
print 'X-arm'
print
print 'Intra-cavity gain (power):', xcavity_gain**2
print 'Cavity reflection (power):', xcavity_reflectivity**2
print 'Cavity transmission (power):', xcavity_transmissivity**2
print
print 'FSR [Hz]:', FSR
print 'Finesse:', xfinesse
print 'Cavity pole [Hz]:', xcavity_pole
print 'Storage time [sec]:', xstorage_time
print

print 'ITM reflectivity:', r_i**2


# Y-arm

L = 2999.8
FSR = c / (2*L)

L_TM = 37.5e-6
RTL = 2*L_TM

t_i = sqrt(0.01376)
t_e = sqrt(4.3e-6)

r_i = sqrt(1 - t_i**2)
r_e = sqrt(1 - t_e**2)

yfinesse, ycavity_pole, ycavity_gain, ycavity_reflectivity, ycavity_transmissivity, ystorage_time = cavity_params(r_i,r_e,RTL,FSR)

print
print 'Y-arm'
print
print 'Intra-cavity gain (power):', ycavity_gain**2
print 'Cavity reflection (power):', ycavity_reflectivity**2
print 'Cavity transmission (power):', ycavity_transmissivity**2
print
print 'FSR [Hz]:', FSR
print 'Finesse:', yfinesse
print 'Cavity pole [Hz]:', ycavity_pole
print 'Storage time [sec]:', ystorage_time
print


# Beamsplitter AR surface losses - for B5 beam
L_BS = 329e-6

t_BS = sqrt(0.5012)
r_BS = sqrt(1-0.5012)

# MICH reflectivity expression - carrier gets a sign flip on reflection from the arms
r_c = r_BS * ycavity_reflectivity * r_BS * (1-L_BS) + t_BS * xcavity_reflectivity * t_BS * (1-L_BS)
print 'MICH reflectivity (power, with arm+BS losses):', r_c**2

# MICH transmission expression
t_c = r_BS * ycavity_reflectivity * t_BS * (1-L_BS) - t_BS * xcavity_reflectivity * r_BS * (1-L_BS)
print 'MICH transmission (power, with arm+BS losses):', t_c**2

# PRC

# POP loss (one-way)
L_prc = 300e-6

# PRC length (one-way)
prcl = 11.952

t_r = sqrt(0.04835)
r_r = sqrt(1-t_r**2)

# Note the minus signs! For an anti-resonant recycling cavity
fcc = ycavity_pole*(1-r_r*r_c)/(1+r_r)

# Calculate recycling gain with PRC losses
pr_cavity_gain = t_r / (1 - r_r * r_c * sqrt(1-L_prc))

# IFO reflectivity
pr_cavity_reflectivity = (r_c*sqrt(1-L_prc) - r_r) / (1-r_r*r_c*sqrt(1-L_prc))

g_cr = pr_cavity_gain
r_PRC = pr_cavity_reflectivity

print
print 'PRC'
print
print 'Intra-cavity gain (power):', pr_cavity_gain**2
print 'Cavity reflection (power):', pr_cavity_reflectivity**2
print 'Coupled cavity pole [Hz]:', fcc
print


"""
t_r = sqrt(0.028)
r_r = sqrt(1-t_r**2)

# Note the minus signs! For an anti-resonant recycling cavity
fcc = ycavity_pole*(1-r_r*r_c)/(1+r_r)

# Calculate recycling gain with PRC losses
pr_cavity_gain = t_r / (1 - r_r * r_c * (1-L_prc)**0.5)

# IFO reflectivity
pr_cavity_reflectivity = (r_c*(1-L_prc)**0.5 - r_r)/(1-r_r*r_c*(1-L_prc)**0.5)

print
print 'PRC - with 2.8% PRM reflection'
print
print 'Intra-cavity gain (power):', pr_cavity_gain**2
print 'Cavity reflection (power):', pr_cavity_reflectivity**2
print 'Coupled cavity pole [Hz]:', fcc
print
"""

print
print '#### Sidebands ####'
print


# Sideband buildup

# Schnupp asymmetry is 23cm
schnupp = 0.23

f_RF = 6*IMC_FSR

# michelson reflectivity - assumes symmetric losses, contrast defect = 0
# uses ITM reflectivity (since sidebands don't resonate in the arms)
r_M = r_i * cos( 2*pi * f_RF * schnupp / c )

r_6M = cos( 2*pi * f_RF * schnupp / c )

# In the following formulae we assume the modulation frequency is resonant in the PRC
# If this is not the case, the buildup and reflectivity will change
# To check this, look at the detuning phase - this should be very small

print 'fMod1 is', f_RF
print 'Cavity detuning [deg] is', (pi-remainder(4*pi*f_RF*prcl/c,pi))*180/pi
print 'Michelson reflectivity for 6MHz is', r_M**2
print 'Michelson transmission for 6MHz is', 1-cos( 2*pi * f_RF * schnupp / c )**2

prc_gain_6MHz = t_r / (1 - r_r * r_M * (1-L_prc))
print 'PRC Cavity buildup for 6MHz is', abs(prc_gain_6MHz)**2

r_prc_6MHz = (r_M*(1-L_prc) - r_r)/(1-r_r*r_M*(1-L_prc))
print 'PRC Cavity reflectivity for 6MHz is', r_prc_6MHz**2
print



f_RF = 2*6*IMC_FSR

# michelson reflectivity - assumes symmetric losses, contrast defect = 0
# uses ITM reflectivity (since sidebands don't resonate in the arms)
r_M = r_i * cos( 2*pi * f_RF * schnupp / c )

# In the following formulae we assume the modulation frequency is resonant in the PRC
# If this is not the case, the buildup and reflectivity will change
# To check this, look at the detuning phase - this should be very small

print 'fMod1 2F is', f_RF
print 'Cavity detuning [deg] is', (pi-remainder(4*pi*f_RF*prcl/c,pi))*180/pi
print 'Michelson reflectivity for 12MHz is', r_M**2
print 'Michelson transmission for 12MHz is', 1-cos( 2*pi * f_RF * schnupp / c )**2

prc_gain_12MHz = t_r / (1 - r_r * r_M * (1-L_prc))
print 'PRC Cavity buildup for 12MHz is', abs(prc_gain_12MHz)**2

r_prc_12MHz = (r_M*(1-L_prc) - r_r)/(1-r_r*r_M*(1-L_prc))
print 'PRC Cavity reflectivity for 12MHz is', r_prc_12MHz**2
print


f_RF = 54*IMC_FSR
r_M = r_i * cos( 2*pi * f_RF * schnupp / c )
print 'fMod2 is', f_RF
print 'Cavity detuning [deg] is', (pi-remainder(4*pi*f_RF*prcl/c,pi))*180/pi
print 'Michelson reflectivity for 56MHz is', r_M**2
print 'Michelson transmission for 56MHz is', 1-cos( 2*pi * f_RF * schnupp / c )**2

prc_gain_56MHz = t_r / (1 - r_r * r_M * (1-L_prc))
print 'PRC Cavity buildup for 56MHz is', abs(prc_gain_56MHz)**2

r_prc_56MHz = (r_M*(1-L_prc) - r_r)/(1-r_r*r_M*(1-L_prc))
print 'PRC Cavity reflectivity for 56MHz is', r_prc_56MHz**2
print


f_RF = 8*IMC_FSR
#r_M = r_i * cos( 2*pi * f_RF * schnupp / c )
print 'fMod3 is', f_RF
print 'Cavity detuning [deg] is', (pi-remainder(4*pi*f_RF*prcl/c,pi))*180/pi
#print 'Michelson reflectivity for 8MHz is', r_M**2

#prc_gain_8MHz = t_r * exp(2*pi*1.0j*f_RF*prcl/c) / (1 - r_r * r_M * exp(4*pi*1.0j*f_RF*prcl/c) * (1-L_prc)**0.5)
#print 'PRC Cavity buildup for 8MHz is', abs(prc_gain_8MHz)**2
print


#f_RF = 126*IMC_FSR
f_RF = 114*IMC_FSR
r_M = r_i * cos( 2*pi * f_RF * schnupp / c )
print 'fMod4 is', f_RF
print 'Cavity detuning [deg] is', (pi-remainder(4*pi*f_RF*prcl/c,pi))*180/pi
print 'Michelson reflectivity for 119MHz is', r_M**2
print 'Michelson transmission for 119MHz is', 1-cos( 2*pi * f_RF * schnupp / c )**2

prc_gain_131MHz = t_r / (1 - r_r * r_M * (1-L_prc))
print 'PRC Cavity buildup for 119MHz is', abs(prc_gain_131MHz)**2

r_prc_131MHz = (r_M*(1-L_prc) - r_r)/(1-r_r*r_M*(1-L_prc)**0.5)
print 'PRC Cavity reflectivity for 119MHz is', r_prc_131MHz**2

fin131 = pi*sqrt(r_i*r_M)/(1-r_i*r_M)
print 'Finesse of 119MHz is', fin131

print





# Power in beams

# Assume that 8MHz and 131MHz are completely reflected by the PRC?

# Only the 6MHz and 56MHz beams have significant power at B2, B4

P_in = 10.0

Gamma6 = 0.22
Gamma56 = 0.1
Gamma119 = 0.12

J0 = jn(0,Gamma6)
J1 = jn(1,Gamma6)

# B1 gets PRC buildup and MICH transmission?

P6_in = P_in * J0 * J1

P6_PRC = P6_in * prc_gain_6MHz**2

P6_MICH_trans = P6_PRC * (1-r_6M**2)

print P6_in
print 'PRC 6MHz buildup', P6_PRC
#print (1-r_6M**2)
print 'MICH trans', P6_MICH_trans 
#print P6_MICH_trans / P6_in


#r_M = r_i * cos( 2*pi * 9e6 * 0.08 / c )
#print 1-cos( 2*pi * 9e6 * 0.08 / c )**2
#print 'LIGO 9MHz MICH transmission:', 1-r_M**2, r_i**2, r_M**2

"""
# B2 gets reflected carrier from PRC and reflected sidebands

b2_carrier = P_in   * (1-4*J1**2)*r_PRC**2
b2_6MHz_1F = 2*P_in * J0*J1*r_prc_6MHz*r_PRC
b2_6MHz_2F = 2*P_in * J1*J1*r_prc_6MHz**2

print J0**8, (1-4*J1**2)
print b2_carrier, b2_6MHz_1F, b2_6MHz_2F

# B4 is 300ppm

pop = 300e-6

b4_carrier = P_in*pop   * (1-4*J1**2)*g_cr**2
b4_6MHz_1F = 2*P_in*pop * J0*J1*prc_gain_6MHz*g_cr
b4_6MHz_2F = 2*P_in*pop * J1*J1*prc_gain_6MHz**2

print b4_carrier, b4_6MHz_1F, b4_6MHz_2F
"""
