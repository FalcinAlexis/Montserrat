#!/Users/thompsong/miniconda3/bin/python

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 13:02:50 2015

@author: thompsong
"""

from obspy import read
import os.path
import numpy as np
import datetime
import sys
import matplotlib.pyplot as plt
import struct
import os
import glob
import shutil
import obspy.core.utcdatetime as utcdt
#import scipy.stats as scistats


def sfilename(wavfile):
    spath = os.path.dirname(wavfile)
    spath = spath.replace('WAV','REA',1)
    wavbase = os.path.basename(wavfile)
    yyyy = wavbase[0:4]
    mm = wavbase[5:7]
    dd = wavbase[8:10]
    HH = wavbase[11:13]
    MI = wavbase[13:15]
    SS = wavbase[16:18]
    sbase = dd + "-" + HH + MI + "-" + SS + "L.S" + yyyy + mm
    sfile = os.path.join(spath, sbase)
    return sfile

def signaltonoise(tr):
    # Here we just make an estimate of the signal-to-noise ratio
    #
    # Normally the trace should be pre-processed before passing to this routine, e.g.
    # * remove ridiculously large values
    # * remove any sequences of 0 from start or end
    # * detrend
    # * bandpass filter
    #
    # Processing:
    #    1. ensure we still have at least 10 seconds
    #    2. take absolute values
    #    3. compute the maximum of each 1-s of data, call this time series M
    #    4. compute 95th and 5th percentile of M, call these M95 and M5
    #    5. estimate signal-to-noise ratio as M95/M5
    
    highval = -1
    lowval = -1
    snr = -1 
    a = tr.data

    fsamp = int(tr.stats.sampling_rate)
    npts = tr.stats.npts
    numseconds = int(npts/fsamp)
    if numseconds > 10:
        a = a[0:int(fsamp * numseconds - 1)]             # remove any fractional second from end of trace
        abs = np.absolute(a)                             # take absolute value of a        
        abs_resize = np.resize(abs, (fsamp, numseconds)) # resize so that each second is one row
        M = np.max(abs_resize,axis=0)                    # find max of each second / row
        highval = np.nanpercentile(M,95)                    # set highval to 95th percentile, to represent signal amplitude
        if highval < 1:
            highval = -1
            return (snr, highval, lowval)
        lowval = np.nanpercentile(M,5)                      # set lowval to 5th percentile, to represent noise amplitude
        snr = highval / lowval
        print(abs_resize.shape)
        print(M.shape)
    return (snr, highval, lowval,)

def compute_metrics(tr, antelopetop):
    # Here we compute simple metrics on each trace and write them to NET.STA.CHAN files. 
    # These metrics are:
    #     1. duration of signal
    #     2. signal amplitude
    #     3. noise amplitude
    #     4. signal-to-noise ratio

    duration = tr.stats.npts /  tr.stats.sampling_rate
    snr = signaltonoise(tr)
    print(snr)
    #dummy = input('any to continue')
    if snr[0] <= 1:
        return False

    starttime = tr.stats['starttime']
    fdir = os.path.join(antelopetop, 'stats', starttime.strftime('%Y/%m') ) 
    if not os.path.exists(fdir):
        os.makedirs(fdir)
    fname = os.path.join(fdir, tr.id + '.txt')
    print('Appending to %s' % fname)
    fileptr = open(fname,'a')
    fileptr.write('%s\t%7.2f\t%9d\t%9d\t%7.2f\n' % (starttime.strftime('%Y-%m-%d %H:%M:%S'),duration,snr[1], snr[2], snr[0],)) 
    fileptr.close()
    return True


def multiTraceSeisan2singleTraceMiniseed(PROJECT_TOP, file, antelopetop, db, flag_compute_metrics, last_endtime ):
    try:
        st = read(file)
    except:
        return

    this_starttime = st[0].stats.starttime
    this_endtime = st[0].stats.endtime
    max_overlap_seconds = 15
    if (this_starttime - last_endtime < -max_overlap_seconds):
        print('this file overlaps the last one')
        return this_endtime

    # Create PNG image file of this STREAM object, if it does not already exist
    pngpath = file + '.png'
    if not os.path.exists(pngpath):
        #10 pixels per second, and assume that 0.8 of plot width is data, therefore
        duration = st[0].stats.delta * st[0].stats.npts
        pixels_wide = (duration * 10) * 1.25
        pixels_high = len(st) * 100 + 50
        st.plot(size=(pixels_wide,pixels_high),outfile=pngpath,equal_scale=False)

    num_nonzero_traces = 0
    for tr in st:

        # ignore blank trace
        anyData = np.count_nonzero(tr.data)
        if anyData==0:
            continue 

############### PRE-PROCESSING #################       
        # remove absurdly large values
        AMP_LIMIT = 10000000
        a = tr.data
        np.clip(a, -AMP_LIMIT, AMP_LIMIT, out=a)
        np.where(a == AMP_LIMIT, 0, a)
        tr.data = a

        # For some SEISAN files from Montserrat - possibly from SEISLOG conversion,
        # the last value in the time series was always some absurdly large value
        # This sections was to change last value to be the mean of the rest
        #a = tr.data
        #a = np.delete(a,[np.size(a) - 1])
        #m = np.mean(a)
        #tr.data[-1] = m

        # remove linear trend
        tr.detrend(type='linear')


############### FIX NET, STA, LOC, CHAN CODES ###
        # fix the network, channel and location
        network = 'MV'
        tr.stats['network']=network
        sta = tr.stats['station'].strip()
        chan = tr.stats['channel'].strip()
        if chan=='PRS' or chan=='APS':
            chan='BDF'
        else:
            if chan[0]=='A':
                if tr.stats['location'] == 'J':
                    bandcode = 'S'
                #else:
                    #bandcode = 'B'
            else:
                if chan[1]=='B':
                    bandcode = 'B'
                else:
                    bandcode = chan[0]
                instrumentcode = 'H'
                if len(chan)==2:
                    orientationcode = tr.stats['location']
                else:
                    orientationcode = chan[2]

                chan = bandcode + instrumentcode + orientationcode

        if chan[0]=='A':
            print(tr.stats)
            print(chan)
            sys.exit('bugger!')
        tr.stats['channel'] = chan
        tr.stats['location']='--'

############### COMPUTE METRICS #################

        # Now high-pass filter, plot and compute_metrics on each trace
        tr.filter('highpass', freq=0.5, corners=2, zerophase=True)
        if flag_compute_metrics:
            flag_trace_okay = compute_metrics(tr, antelopetop)
            if not flag_trace_okay:
                next

############## WRITE SINGLE TRACES TO FILE ######

        # get sta, chan
        sta = tr.stats['station']
        chan = tr.stats['channel']

        # get time strings for this trace
        starttime = tr.stats['starttime']
        jjj = str(starttime.julday).zfill(3)
        yyyy = str(starttime.year).zfill(4)
        mm = str(starttime.month).zfill(2)
        dd = str(starttime.day).zfill(2)
        hh = str(starttime.hour).zfill(2)
        ii = str(starttime.minute).zfill(2)
        ss = str(starttime.second).zfill(2)

        # write out to SDS structure - but this only works for full day files
        if tr.stats.npts / tr.stats.sampling_rate == 86400:
            sdsdir = os.path.join("SDS", yyyy, network, sta, chan + "D")
            if not os.path.exists(sdsdir):
                os.makedirs(sdsdir)
            sdsfilename = tr.id + ".D." + yyyy + "." + jjj
            sdsfullpath = os.path.join(sdsdir, sdsfilename)
            if not os.path.exists(sdsfullpath):
                tr.write(sdsfullpath, format="MSEED")

        # Write Miniseed files to Antelope structure YYYY/JJJ/sta.chan.yyyy:jjj:hh:mm:ss
        antelopedir = os.path.join(antelopetop, db, yyyy, mm)
        if not os.path.exists(antelopedir):
            os.makedirs(antelopedir)
        antelopefile = sta + "." + chan + "." + yyyy + mm + dd + "_" + hh + ii + ss
        antelopefullpath = os.path.join(antelopedir, antelopefile)
        if not os.path.exists(antelopefullpath):
            tr.write(antelopefullpath, format="MSEED")

        # Plot seismogram to an PNG file
        pngpath = antelopefullpath + '.png'
        if not os.path.exists(pngpath):
            #10 pixels per second, and assume that 0.8 of plot width is data, therefore
            duration = tr.stats.delta * tr.stats.npts
            pixels_wide = (duration * 10) * 1.25
            tr.plot(size=(pixels_wide,300),outfile=pngpath)

        # Another good trace if we got this far 
        num_nonzero_traces += 1

    # Log number of good traces for this event
    f = open(os.path.join(PROJECT_TOP,'num_good_traces.txt'),"a")
    startt = st[0].stats.starttime
    starttfmt = startt.strftime('%Y-%m-%d %H:%M:%S')
    f.write("%s\t%d\n" % (starttfmt,num_nonzero_traces,))
    f.close()
    return this_endtime

def swap32(i):
    return struct.unpack("<i", struct.pack(">i", i))[0]

def main():
    PROJECT_TOP = '/Users/thompsong/Desktop/IPGP_Thompson_collaboration'
    antelopetop = os.path.join(PROJECT_TOP, 'miniseed')
    seisantop = os.path.join(PROJECT_TOP, 'seismo')
    db = 'MVOE_' #seisan database code
    last_endtime = utcdt.UTCDateTime(1970, 1, 1, 0, 0)
    flag_compute_metrics = False
    flag_reinitialize = True
    if flag_reinitialize:
        if os.path.exists('num_good_traces.txt'):
            os.remove('num_good_traces.txt')
        if os.path.exists(antelopetop):
            shutil.rmtree(antelopetop)

    # Start processing
    os.chdir(os.path.join(seisantop,'WAV',db))
    print(os.getcwd())
    years = sorted(glob.glob('[0-9]' * 4))
    for year in years:
        os.chdir(year)
        print(os.getcwd())
        months = sorted(glob.glob('[0-9]' * 2))
        for month in months:
            os.chdir(month)
            print(os.getcwd())

            # recompute stats - so blow away file before
            fdir = os.path.join(antelopetop, 'stats', year, month)
            fname = os.path.join(fdir, 'M*.txt')
            files = glob.glob(fname)
            if not files:
                flag_compute_metrics = True

            seisanfiles = sorted(glob.glob('*S.MVO*__0??'))
            print("Found %d Seislog files" % (len(seisanfiles)))
            for seisanfile in seisanfiles:
                last_endtime = multiTraceSeisan2singleTraceMiniseed(PROJECT_TOP, seisanfile, antelopetop, db, flag_compute_metrics, last_endtime)
                #pass
            os.chdir('..')
        os.chdir('..')

if __name__ == "__main__":
	main()
