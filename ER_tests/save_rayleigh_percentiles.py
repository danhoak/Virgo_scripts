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
import os, sys, time
from glue import segmentsUtils
import h5py

from matplotlib.mlab import *
from scipy import stats
from glue.segments import *

seg_start = 1128902417  # Oct 15 2015 00:00
#seg_stop  = 1129334417  # Oct 20 2015 00:00
seg_stop  = 1133568017  # Dec 08 2015 00:00
#seg_stop = 1134172817  # Dec 15 2015 00:00

epoch = segmentlist(segmentsUtils.segmentlist_range(seg_start, seg_stop, seg_stop-seg_start))


ifo = 'L1'
res = '6sec'

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


# get CAT1 veto segs from O1

ep = 'O1'
veto_gps_string = '1127271617-9982800'
cat1_file = '/home/dhoak/public_html/segs/' + ifo + '-VETOTIME_CAT1-' + veto_gps_string + '.xml'
O1_txt1_file  = ep + '_' + ifo + '_cat1.txt'
cat2_file = '/home/dhoak/public_html/segs/' + ifo + '-VETOTIME_CAT2-' + veto_gps_string + '.xml'
O1_txt2_file  = ep + '_' + ifo + '_cat2.txt'
xml_cmd = "ligolw_print --table segment --column start_time --column end_time --delimiter ' ' " + cat1_file +\
    " | awk '{print NR  \" \" $1 \" \" $2 \" \" $2-$1}' > " + O1_txt1_file
os.system(xml_cmd)

O1_cat1seg = segmentsUtils.fromsegwizard(open(O1_txt1_file)); O1_cat1seg.coalesce()

sci_good = sciseg - (O1_cat1seg & epoch)

#print sciseg
#print cat1seg
print 'Total science time is', getSegmentSum(sci_good), 'seconds'
print 'Number of good segments is', len(sci_good)

# loop over segments and grab the data
# start arrays to hold the data
i=0
sample_list = []
segment_list = []
for segment in sci_good:
    segment_start = int(segment[0])
    segment_stop = int(segment[1])
    duration = segment_stop - segment_start

    # counter to keep track of files
    i = i+1
    print i, segment_start, segment_stop, duration

    if duration < 512:
        print 'Segment is too short!'
        continue

    # Get the strain data
    #frame_type = ifo + '_HOFT_C00'   # C00 frames used after Oct 9
    #chans = [ifo + ':GDS-CALIB_STRAIN']

    segment_list.append(i)

    frame_type = ifo + '_HOFT_C01'
    chans = [ifo + ':DCS-CALIB_STRAIN_C01']

    #Fs = 16384
    Fs = 4096

    s1 = time.time()
    y = getChannelVector(ifo[0],str(segment_start),str(segment_stop),frame_type,chans,Fs,filter_rate=None)
    s2 = time.time()
    print 'Time to acquire data:', (s2-s1)/duration, 'seconds per sec of data'

    # Define parameters for FFT
    stride = 6.0   # FFT stride in seconds
    overlap = stride/2  # overlap in seconds (50%)
    
    # for science segment, get PSD.
    s1 = time.time()
    Pxx, freq, t = specgram(y[:,0], NFFT=int(stride*Fs), Fs=int(Fs), noverlap=int(overlap*Fs))
    s2 = time.time()
    print 'Time to transform data:', (s2-s1)/duration, 'seconds per sec of data'

    Axx = sqrt(Pxx)
    [frequencies] = shape(freq)

    ASD = Axx

    freqs, samples = shape(ASD)
    sample_list.append(samples)

    """
    s1 = time.time()
    if i==1:
       ASD = Axx
    else:
        ASD = hstack((ASD,Axx))
    s2 = time.time()
    print 'Time to concatenate data:', (s2-s1)/duration, 'seconds per sec of data'
    """

    datafile = ifo + '_segs/' + ifo + '_O1_ASDs_' + res + '_' + str(i) + '.hdf5'
    f = h5py.File(datafile,'w')

    ASD_dset = f.create_dataset('ASD', data=ASD)

    f.close()


# save a datafile with important numbers
datafile = ifo + '_segs/' + ifo + '_O1_ASDs_' + res + '_header.hdf5'
f = h5py.File(datafile,'w')

#segs_dset = f.create_dataset('segments', data=segments)
samples_dset = f.create_dataset('sample_list', data=sample_list)
segment_dset = f.create_dataset('segment_list', data=segment_list)
f_dset = f.create_dataset('freq', data=freq)

f.close()
