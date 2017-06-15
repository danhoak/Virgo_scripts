#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
from scipy import constants
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab
#from scipy.special import jn



c = constants.c
lamb = 1064.0e-9
k = 2*pi/lamb

L_BSNI = 6.016
L_BSWI = 5.786

# PRC length to bring carrier on resonance
# assume this is locked as we change MICH fringe
L_PRC = 11.952
#L_PRC = 57.65575 # LIGO
#L_PRC = 57.656 # LIGO

L_PRBS = L_PRC - (L_BSNI + L_BSWI) / 2

schnupp = (L_BSNI - L_BSWI)


f0 = c / lamb

# MICH half-fringe offset - each arm gets this differentially
dL = arange(-lamb/3,lamb/3,1e-10)

Lx = L_BSNI + dL/2
Ly = L_BSWI - dL/2

# arm reflectivity with losses - change this assumption?
rM = sqrt(0.977)

tr = sqrt(0.04835)
rr = sqrt(1-0.04835)



fignum=0

matplotlib.rcParams.update({'savefig.dpi':250})



# Calculate PRCL error signal as a function of MICH reflectivity, for the 8MHz sideband

f_RF = 6270777.0
#f_RF = 8361036.0
#f_RF = 56436993.0
#f_RF = 119144763.0

#f_RF = 45497355.0 # LIGO

w0 = 2*pi*f0
w_RF = 2*pi*f_RF



def refl_coeff(r_i, r_e, w, L):
    # equation 4.1 from anamaria
    return (r_i - r_e * exp(-2j*w*L/c)) / (1.0 - r_i*r_e*exp(-2j*w*L/c))



fignum=fignum+1
pylab.figure(fignum)

ax1 = pylab.subplot(2,1,1)

L = L_PRC + dL
PD_I = imag( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
PD_Q = real( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
idx = argmax(gradient(PD_I))
x = L[idx]
pylab.plot((L-x)*1e9,PD_I,'r-',linewidth=1.2,label='B2_I')
pylab.plot((L-x)*1e9,PD_Q,'b-',linewidth=1.2,label='B2_Q')


delta = -1e-3
L = L_PRC + delta + dL
PD_I = imag( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
PD_Q = real( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
idx = argmax(gradient(PD_I))
x = L[idx]
pylab.plot((L-x)*1e9,PD_I,'r--',linewidth=1.2,label='B2_I, 10mm')
pylab.plot((L-x)*1e9,PD_Q,'b--',linewidth=1.2,label='B2_Q, 10mm')

delta = -5e-3
L = L_PRC + delta + dL
PD_I = imag( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
PD_Q = real( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
idx = argmax(gradient(PD_I))
x = L[idx]
pylab.plot((L-x)*1e9,PD_I,'m--',linewidth=1.2,label='B2_I, 30mm')
pylab.plot((L-x)*1e9,PD_Q,'c--',linewidth=1.2,label='B2_Q, 30mm')




pylab.legend(loc=1,prop={'size':8},fancybox=True)
pylab.grid(True,which='both')
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
pylab.xlim(-10,10)
pylab.ylabel(r'S$_{refl}$, Carrier Lock',fontsize=10)
ax1.get_yaxis().set_label_coords(-0.06,0.5)

#pylab.title('PRITF Beam Powers, Effect of Mode Matching + Lossy PRC',fontsize=12)


ax2 = pylab.subplot(2,1,2)


L = L_PRC + dL
PD_I = imag( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
PD_Q = real( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
idx = argmin(gradient(PD_I))
x = L[idx]
pylab.plot((L-x)*1e9,PD_I,'r-',linewidth=1.2,label='B2_I')
pylab.plot((L-x)*1e9,PD_Q,'b-',linewidth=1.2,label='B2_Q')


delta = -1e-3
L = L_PRC + delta + dL
PD_I = imag( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
PD_Q = real( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
idx = argmin(gradient(PD_I))
x = L[idx]
pylab.plot((L-x)*1e9,PD_I,'r--',linewidth=1.2,label='B2_I, 10mm')
pylab.plot((L-x)*1e9,PD_Q,'b--',linewidth=1.2,label='B2_Q, 10mm')

I_grad = mean(gradient(PD_I[idx-5:idx+5]))
Q_grad = mean(gradient(PD_Q[idx-5:idx+5]))

print Q_grad/I_grad
print 'Deg per 1mm:', arctan(Q_grad/I_grad)*180/pi



delta = -5e-3
L = L_PRC + delta + dL
PD_I = imag( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
PD_Q = real( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
idx = argmin(gradient(PD_I))
x = L[idx]
pylab.plot((L-x)*1e9,PD_I,'m--',linewidth=1.2,label='B2_I, 30mm')
pylab.plot((L-x)*1e9,PD_Q,'c--',linewidth=1.2,label='B2_Q, 30mm')

I_grad = mean(gradient(PD_I[idx-5:idx+5]))
Q_grad = mean(gradient(PD_Q[idx-5:idx+5]))

print Q_grad/I_grad
print 'Deg per 5mm:', arctan(Q_grad/I_grad)*180/pi


pylab.legend(loc=1,prop={'size':8},fancybox=True)
pylab.grid(True,which='both')
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
pylab.xlim(-10,10)
pylab.ylabel(r'S$_{refl}$, 6MHz Sideband Lock',fontsize=10)
ax2.get_yaxis().set_label_coords(-0.06,0.5)


pylab.xlabel(r'Cavity detuning [nm]',fontsize=10)

pylab.savefig('PRCL_length_signal.png',bbox_inches='tight')






fignum=fignum+1
pylab.figure(fignum)

ax1 = pylab.subplot(2,1,1)

delta = arange(-50e-3,50e-3,1e-3)
demod = zeros(len(delta))

for i in range(len(delta)):

    L = L_PRC + delta[i] + dL
    PD_I = imag( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
    PD_Q = real( refl_coeff(rr, rM, w0, L) * conj(refl_coeff(rr, rM, w0+w_RF, L)) - conj(refl_coeff(rr, rM, w0, L)) * refl_coeff(rr, rM, w0-w_RF, L) )
    idx = argmin(gradient(PD_I))
    x = L[idx]

    I_grad = mean(gradient(PD_I[idx-5:idx+5]))
    Q_grad = mean(gradient(PD_Q[idx-5:idx+5]))

    demod[i] = 180 - abs(arctan(Q_grad/I_grad)*180/pi)


pylab.plot(delta*1e3,demod,'b--',linewidth=0.8,markersize=4,label='B2_I, 30mm')

#pylab.legend(loc=1,prop={'size':8},fancybox=True)
pylab.grid(True,which='both')
pylab.yticks(fontsize=8)
pylab.xticks(fontsize=8)
pylab.xlim(-20,20)
pylab.ylabel(r'B2 6MHz Phase [deg]',fontsize=10)
ax1.get_yaxis().set_label_coords(-0.06,0.5)

pylab.xlabel(r'PRCL Length Detuning [mm]',fontsize=10)

pylab.savefig('PRCL_length_detune.png',bbox_inches='tight')
