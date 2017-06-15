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
golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size = [fig_width,fig_height]

matplotlib.rcParams.update({'savefig.dpi':250,
                            'text.usetex':True,
                            'figure.figsize':fig_size,
                            'font.family':"serif",
                            'font.serif':["Times"],})

from matplotlib.mlab import *
from scipy import stats
from glue.segments import *

seg_start = 1128470417  # Oct 10 2015 00:00
seg_stop  = 1128556817  # Oct 11 2015 00:00
#seg_stop  = 1130198417  # Oct 30 2015 00:00

epoch = segmentlist(segmentsUtils.segmentlist_range(seg_start, seg_stop, seg_stop-seg_start))

fmin = 10
fmax = 6000



##########################
#
# Begin H1 section
#
##########################

ifo = 'H1'

print
print 'Grabbing data for', ifo

# Get segments from segdb

seg_flag = ifo + ':DMT-ANALYSIS_READY:1'
xml_file = ifo + '_seg.xml'
txt_file = ifo + '_seg.txt'

seg_cmd = 'ligolw_segment_query_dqsegdb --query-segments --segment-url https://segments.ligo.org --gps-start-time ' + \
    str(seg_start) + ' --gps-end-time ' + str(seg_stop) + ' --include-segments ' + seg_flag + ' --output-file ' + xml_file
os.system(seg_cmd)

xml_cmd = "ligolw_print --table segment --column start_time --column end_time --delimiter ' ' " + xml_file +\
    " | awk '{print NR  \" \" $1 \" \" $2 \" \" $2-$1}' > " + txt_file
os.system(xml_cmd)

sciseg = segmentsUtils.fromsegwizard(open(txt_file)); sciseg.coalesce()
print 'Number of segments is', len(sciseg)

# get CAT1 veto segs and remove them from the scisegs
# figure out how much time is in cat1, typical duration of flag, whether it breaks up scisegs into tiny pieces...

cat1_file = ifo + '-VETOTIME_CAT1-1127271617-2111400.xml'
txt1_file  = ifo + '_cat1.txt'
cat2_file = ifo + '-VETOTIME_CAT2-1127271617-2111400.xml'
txt2_file  = ifo + '_cat2.txt'
xml_cmd = "ligolw_print --table segment --column start_time --column end_time --delimiter ' ' " + cat1_file +\
    " | awk '{print NR  \" \" $1 \" \" $2 \" \" $2-$1}' > " + txt1_file
os.system(xml_cmd)
#xml_cmd = "ligolw_print --table segment --column start_time --column end_time --delimiter ' ' " + cat2_file +\
#    " | awk '{print NR  \" \" $1 \" \" $2 \" \" $2-$1}' > " + txt2_file
#os.system(xml_cmd)

cat1seg = segmentsUtils.fromsegwizard(open(txt1_file)); cat1seg.coalesce()
#cat2seg = segmentsUtils.fromsegwizard(open(txt2_file)); cat1seg.coalesce()
print 'Number of CAT1 veto segments is', len(cat1seg & epoch)
#print 'Number of CAT2 veto segments is', len(cat2seg & epoch)

print 'Time lost to CAT1:', getSegmentSum(cat1seg & epoch) / float(getSegmentSum(sciseg))
#print 'Time lost to CAT2:', getSegmentSum(cat2seg & epoch) / float(getSegmentSum(sciseg))

#sci_good = sciseg - (cat1seg & epoch) - (cat2seg & epoch)
sci_good = sciseg - (cat1seg & epoch)

#print sciseg
#print cat1seg
#print sci_good
print 'Number of good segments is', len(sci_good)

# loop over segments and grab the data
# start arrays to hold the data
i=0
for segment in sci_good:
    i = i+1
    segment_start = int(segment[0])
    segment_stop = int(segment[1])
    duration = segment_stop - segment_start
    print i, segment_start, segment_stop, duration

    if duration < 128:
        print 'Segment is too short!'
        continue

    # Get the strain data
    frame_type = ifo + '_HOFT_C00'   # C00 frames used after Oct 9
    chans = [ifo + ':GDS-CALIB_STRAIN']
    Fs = 16384

    y = getChannelVector(ifo[0],str(segment_start),str(segment_stop),frame_type,chans,Fs,filter_rate=None)

    # Define parameters for FFT
    stride = 6.0   # FFT stride in seconds
    overlap = 3.0  # overlap in seconds (50%)
    
    # for science segment, get PSD.
    # 6-second FFTs, 50% overlap, hann window (default)
    Pxx, freq, t = specgram(y[:,0], NFFT=int(stride*Fs), Fs=int(Fs), noverlap=int(overlap*Fs))

    Axx = sqrt(Pxx)
    [frequencies] = shape(freq)

    if i==1:
       ASD = Axx
    else:
        ASD = hstack((ASD,Axx))

    print shape(Axx), shape(t), shape(freq)
    print shape(ASD)



datafile = ifo + '_Oct_ASDs.hdf5'
f = h5py.File(datafile,'w')

ASD_dset = f.create_dataset('ASD', data=ASD)

f.close()

"""

p50 = empty(frequencies)
p75 = empty(frequencies)
p90 = empty(frequencies)
p95 = empty(frequencies)
p975 = empty(frequencies)
p99 = empty(frequencies)
p995 = empty(frequencies)
p999 = empty(frequencies)

for i in range(frequencies):

    p50[i] = stats.scoreatpercentile(ASD[i,:],50)
    p75[i] = stats.scoreatpercentile(ASD[i,:],75)
    p90[i] = stats.scoreatpercentile(ASD[i,:],90)
    p95[i] = stats.scoreatpercentile(ASD[i,:],95)
    p975[i] = stats.scoreatpercentile(ASD[i,:],97.5)
    p99[i] = stats.scoreatpercentile(ASD[i,:],99)
    p995[i] = stats.scoreatpercentile(ASD[i,:],99.5)
    p999[i] = stats.scoreatpercentile(ASD[i,:],99.9)


# Get percentiles for standard Rayleigh distribution with 1e6 samples
x = random.rayleigh(1.0,1000000)
xp50 = stats.scoreatpercentile(x,50)
xp75 = stats.scoreatpercentile(x,75)
xp90 = stats.scoreatpercentile(x,90)
xp95 = stats.scoreatpercentile(x,95)
xp975 = stats.scoreatpercentile(x,97.5)
xp99 = stats.scoreatpercentile(x,99)
xp995 = stats.scoreatpercentile(x,99.5)
xp999 = stats.scoreatpercentile(x,99.9)

xmp75 = xp75/xp50
xmp90 = xp90/xp50
xmp95 = xp95/xp50
xmp975 = xp975/xp50
xmp99 = xp99/xp50
xmp995 = xp995/xp50
xmp999 = xp999/xp50


fignum=0
fignum=fignum+1
pylab.figure(fignum)

pylab.semilogy(freq,p50,"b-",linewidth=0.5)
pylab.semilogy(freq,p75,"g-",linewidth=0.5)
pylab.semilogy(freq,p90,"r-",linewidth=0.5)
pylab.semilogy(freq,p95,"c-",linewidth=0.5)
pylab.semilogy(freq,p975,"m-",linewidth=0.5)
pylab.semilogy(freq,p99,"y-",linewidth=0.5)
#pylab.semilogy(freq,p995,"k-",linewidth=0.5)
pylab.semilogy(freq,p999,"k-",linewidth=0.5)

pylab.axis([fmin, fmax, 3e-24, 3e-19])

pylab.grid(True, which='both', linestyle=':')

pylab.xlabel('Frequency [Hz]',fontsize=12)
pylab.ylabel('Strain ASD [1/rt(Hz)]',fontsize=12)
pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)

pylab.title('Distribution of detector noise for ' + ifo,fontsize=12)
#pylab.legend(('Median','75th percentile','90th percentile','95th percentile','97.5th percentile','99th percentile','99.5th percentile'),loc=1,prop={'size':10},fancybox=True)
pylab.legend(('Median','75th percentile','90th percentile','95th percentile','97.5th percentile','99th percentile','99.9th percentile'),loc=1,prop={'size':9},fancybox=True)

pylab.savefig('G184098_' + ifo + '_asd_rayleigh_180.png',bbox_inches='tight')
#pylab.savefig('G184098_' + ifo + '_asd_rayleigh.pdf',bbox_inches='tight')
pylab.close()




fignum=fignum+1
pylab.figure(fignum)

pylab.plot(freq,p50/p50,"b-",linewidth=0.5)
pylab.plot(freq,p75/p50,"g-",linewidth=0.5)
pylab.plot(freq,p90/p50,"r-",linewidth=0.5)
pylab.plot(freq,p95/p50,"c-",linewidth=0.5)
pylab.plot(freq,p975/p50,"m-",linewidth=0.5)
pylab.plot(freq,p99/p50,"y-",linewidth=0.5)
#pylab.plot(freq,p995/p50,"k-",linewidth=0.5)
pylab.plot(freq,p999/p50,"k-",linewidth=0.5)

pylab.plot([fmin, fmax],[1.0,1.0],"b--",linewidth=0.5)
pylab.plot([fmin, fmax],[xmp75,xmp75],"g--",linewidth=0.5)
pylab.plot([fmin, fmax],[xmp90,xmp90],"r--",linewidth=0.5)
pylab.plot([fmin, fmax],[xmp95,xmp95],"c--",linewidth=0.5)
pylab.plot([fmin, fmax],[xmp975,xmp975],"m--",linewidth=0.5)
pylab.plot([fmin, fmax],[xmp99,xmp99],"y--",linewidth=0.5)
#pylab.plot([fmin, fmax],[xmp995,xmp995],"k--",linewidth=0.5)
pylab.plot([fmin, fmax],[xmp999,xmp999],"k--",linewidth=0.5)

pylab.axis([fmin, fmax, 0.9, 10.0])

pylab.grid(True, which='both', linestyle=':')

pylab.xlabel('Frequency [Hz]',fontsize=12)
pylab.ylabel('Normalized Strain ASD',fontsize=12)
pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)

pylab.title('Median-normalized detector noise for ' + ifo,fontsize=12)
#pylab.legend(('Median','75th percentile','90th percentile','95th percentile','97.5th percentile','99th percentile','99.5th percentile'),loc=1,prop={'size':10},fancybox=True)
pylab.legend(('Median','75th percentile','90th percentile','95th percentile','97.5th percentile','99th percentile','99.9th percentile'),loc=1,prop={'size':9},fancybox=True)

pylab.savefig('G184098_' + ifo + '_asd_rayleigh_deviation_180.png',bbox_inches='tight')
#pylab.savefig('G184098_' + ifo + '_asd_rayleigh_deviation.pdf',bbox_inches='tight')
pylab.close()

"""
