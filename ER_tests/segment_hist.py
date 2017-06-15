#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script plots a histogram of the segment length in a given science run epoch. 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
matplotlib.use("Agg")
import pylab
import commands

ifo = 'H1'

epoch = 'S6D'

start_time = '961545615'
end_time = '971654415'

cmd1 = 'ligolw_segment_query -dqa "' + ifo + ':DMT-SCIENCE" -s ' + start_time + ' -e ' + end_time + ' -o ' + 'segments.xml'
cmd2 = 'ligolw_print -t segment -c start_time -c end_time ' + 'segments.xml -d " "'


commands.getoutput(cmd1)
segment_list = commands.getoutput(cmd2)
cmd3 = 'rm segments.xml'
commands.getoutput(cmd3)

segments = segment_list.rsplit('\n')

numsegs = len(segments)

seglengths = zeros((numsegs))

i=0
for segment in segments:
    i+=1
    segment_start,segment_stop = segment.split()
    start = int(segment_start)
    stop = int(segment_stop)
    duration = stop-start

    seglengths[i-1] = duration

#    print start, stop, duration


print max(seglengths)

fignum=0
fignum=fignum+1
pylab.figure(fignum)

dmin = 1
dmax = 150000

nbins = 100
logdmin=log10(dmin)
logdmax=log10(dmax)
dbinwidth = (logdmax-logdmin)/nbins
dbins=[0.0]*(nbins+1)
for i in range(len(dbins)):
   dbins[i] = pow(10,logdmin+i*dbinwidth)

n, b = histogram(seglengths, dbins)

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

pylab.plot(plotbins,plotcounts,'k-',linewidth=1.0)

pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('Duration [sec]',fontsize=14)
pylab.grid(True, which='major', linestyle=':')
pylab.grid(True, which='minor', linestyle=':')
pylab.xlim(dmin,dmax)
pylab.ylim(0.7,1.2*max(plotcounts))
pylab.title(ifo + ' ' + epoch + ' ' + 'Segment Duration')
pylab.savefig('/home/dhoak/public_html/' + ifo + '_' + epoch + '_segmentDurationHist.png')


