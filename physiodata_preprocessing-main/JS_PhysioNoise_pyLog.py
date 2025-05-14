#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import sys
from numpy import *
from scipy.signal import cspline1d, cspline1d_eval, iirdesign, lfilter, lp2lp, freqz
from matplotlib import *
import pylab 
from pylab import plot, title, show , legend, figure, subplot, savefig, detrend_linear,  window_hanning,  slopes,  semilogy,  ylabel
import string
from scipy import interpolate
import scipy.linalg
from scipy.ndimage import spline_filter1d
# from scikits.timeseries.lib.moving_funcs import cmov_average
from matplotlib.mlab import psd

#Modified from peakdet.m by Eli Billauer, 3.4.05 (Explicitly not copyrighted).
#This function is released to the public domain; Any use is allowed.
#Accessed Sept 2007

def findextrema(TS, ChangeMag,TimeConstraint):
        MinMag = Inf
        MaxMag = -Inf
        FindMaxMag = 1
        MaxCounts=zeros(len(TS))
        MinCounts=zeros(len(TS))
        MaxCountsPeak=zeros(len(TS),float)
        MinCountsPeak=zeros(len(TS),float)
        maxcounter=0
        mincounter=0
        maxmincounter=0
        AllExtremaPeak=zeros(len(TS),float)
        AllExtremaPeakTime=zeros(len(TS))
        MaxPeak=0
        MaxPeakTime= -2*TimeConstraint
        MinPeak=0
        MinPeakTime= -2*TimeConstraint
        MinTime =  -2*TimeConstraint
        MaxTime = -2*TimeConstraint
        constraintindex=0
        AllExtremaPeakTime[0]= -2*TimeConstraint


        print "Finding Extrema"
        for i in xrange(len(TS)):
                TestMag = TS[i]
                if TestMag > MaxMag:
                        MaxMag = TestMag
                        MaxTime = i

                if TestMag < MinMag:
                        MinMag = TestMag
                        MinTime = i

                if FindMaxMag == 1:
                        if TestMag < ( MaxMag - ChangeMag) and abs(MaxTime-AllExtremaPeakTime[constraintindex])>=TimeConstraint:
                                MaxPeak = MaxMag
                                MaxPeakTime = MaxTime
                                MaxCounts[maxcounter] = MaxTime
                                MaxCountsPeak[maxcounter]= MaxMag
                                maxcounter= maxcounter + 1
                                AllExtremaPeak[maxmincounter]=MaxMag
                                AllExtremaPeakTime[maxmincounter]=MaxTime
                                constraintindex=maxmincounter
                                maxmincounter=maxmincounter+1
                                MinMag = TestMag
                                MinTime = i
                                FindMaxMag = 0

                elif TestMag > (MinMag + ChangeMag) and  abs(MinTime-AllExtremaPeakTime[constraintindex] )>=TimeConstraint:
                        MinPeak = MinMag
                        MinPeakTime = MinTime
                        MinCounts[mincounter]= MinTime
                        MinCountsPeak[mincounter]= MinMag
                        mincounter = mincounter + 1
                        AllExtremaPeak[maxmincounter]=MinMag
                        AllExtremaPeakTime[maxmincounter]=MinTime
                        constraintindex=maxmincounter
                        maxmincounter=maxmincounter+1
                        MaxMag = TestMag
                        MaxTime = i
                        FindMaxMag = 1
                        

        MaxCountsPeak = MaxCountsPeak[0:maxcounter]
        MaxCounts = MaxCounts[0:maxcounter]
        MinCountsPeak = MinCountsPeak[0:mincounter]
        MinCounts = MinCounts[0:mincounter]
        AllExtremaPeak= AllExtremaPeak[0:maxmincounter]
        AllExtremaPeakTime=AllExtremaPeakTime[0:maxmincounter]


        TopEnvTemp=zeros(size(TS))
        BotEnvTemp=zeros(size(TS))
        for counter in xrange(len(MaxCounts)):
                topms = MaxCounts[counter]
                TopEnvTemp[int(topms)] = MaxCountsPeak[counter]

        for counter in xrange(len(MinCounts)):
                botms = MinCounts[counter]
                BotEnvTemp[int(botms)] = MinCountsPeak[counter]

        return AllExtremaPeak, AllExtremaPeakTime,  MaxCountsPeak, MinCountsPeak, TopEnvTemp, BotEnvTemp, MaxCounts, MinCounts


#--------------------------------------------------------------------# FiltFilt
from numpy import vstack, hstack, eye, ones, zeros, linalg, \
newaxis, r_, flipud, convolve, matrix, array
from scipy.signal import lfilter

def lfilter_zi(b,a):
        #compute the zi state from the filter parameters. see [Gust96].

        #Based on:
        # [Gust96] Fredrik Gustafsson, Determining the initial states in forward-backward
        # filtering, IEEE Transactions on Signal Processing, pp. 988--992, April 1996,
        # Volume 44, Issue 4

        n=max(len(a),len(b))
        zin = (  eye(n-1) - hstack( (-a[1:n,newaxis],
        vstack((eye(n-2), zeros(n-2))))))
        zid=  b[1:n] - a[1:n]*b[0]
        zi_matrix=linalg.inv(zin)*(matrix(zid).transpose())
        zi_return=[]
        #convert the result into a regular array (not a matrix)

        for i in range(len(zi_matrix)):
                zi_return.append(float(zi_matrix[i][0]))
        return array(zi_return)

def filtfilt(b,a,x):
        #For now only accepting 1d arrays
        ntaps=max(len(a),len(b))
        edge=ntaps*3
        if x.ndim != 1:
                raise ValueError, "Filiflit is only accepting 1 dimension arrays."

        #x must be bigger than edge
        if x.size < edge:
                raise ValueError, "Input vector needs to be bigger than 3 * max(len(a),len(b)."
        if len(a) < ntaps:
                a=r_[a,zeros(len(b)-len(a))]

        if len(b) < ntaps:
                b=r_[b,zeros(len(a)-len(b))]
        zi=lfilter_zi(b,a)
        #Grow the signal to have edges for stabilizing
        #the filter with inverted replicas of the signal
        s=r_[2*x[0]-x[edge:1:-1],x,2*x[-1]-x[-1:-edge:-1]]
        #in the case of one go we only need one of the extrems
        # both are needed for filtfilt
        (y,zf)=lfilter(b,a,s,-1,zi*s[0])
        (y,zf)=lfilter(b,a,flipud(y),-1,zi*y[-1])
        return flipud(y[edge-1:-edge+1])

def ConnectTheDots(Counts,Hills,TotPoints,RepPt):
        if len(Counts)<=3:
                print "Not enough peaks detected"
                return zeros(len(TotPoints))
        tc=long(Counts[-1]-Counts[RepPt]+1)
        Tp=interpolate.interp1d(Counts[RepPt:],Hills)
        InterpTime=range(tc)
        for i in  range(tc):
                InterpTime[i]=Counts[RepPt] + i
        InterpEnv = Tp(InterpTime)
        FullInterpEnv=zeros(len(TotPoints))
        for i in  range(tc):
                FullInterpEnv[int(Counts[RepPt] + i) ]= InterpEnv[i]
        FullInterpEnv[0:int(Counts[RepPt])]=FullInterpEnv[int(Counts[RepPt])]
        FullInterpEnv[int(Counts[-1]):]=FullInterpEnv[int(Counts[-1])]
        return FullInterpEnv

def movingaverage(interval, window_size):
    window = ones(int(window_size))/float(window_size)
    return convolve(interval, window, 'same')






from optparse import OptionParser, OptionGroup
usage ="""
===========================================================================
Usage: %prog [options] arg1 arg2 arg3
Author: Dan Kelley
	Waisman Laboratory for Brain Imaging and Behavior
	University of Wisconsin-Madison
Rewritten by Joerg Stadler
    Leibniz Institute for Neurobiology
    Magdeburg
    Last Update 1.2.2010
    
This program:
1 Read Log files of physiological data and remove the trailing points
2 Removes noise with a low pass Butterworth filter with zero phase shift
3 Calculates the residual = Raw signal - filtered signal
4 Idenitifies dirty peaks with big residuals that need cleaning
5 Dirty peaks  +/- a window are removed to produce filtered, cleaned data
6 Calculates the best fitting spline through the cleaned data
7 Spline interpolates over the removed dirty peaks to produce a clean
respiratory waveform (RW) and cardiac waveform (CW)
This is fast on Mac OS. Other platforms use the --speed option.
8 Finds the respiratory peaks (RPpd) and cardiac peaks (CPpd).
Peak finding was slightly modified from the open source peakdet.m algorithm
to include magnitude and time thresholds and ported into python.
9 The top and bottom envelopes are generated using the RPpd. The
respiratory volumes over time (RVT) is their difference as in:
        Birn RM, Diamond JB, Smith MA, Bandettini PA.
        Separating respiratory-variation-related fluctuations
        from neuronal-activity-related fluctuations in fMRI.
        Neuroimage. 2006 Jul 15;31(4):1536-48. Epub 2006 Apr 24.
        PMID: 16632379
The same is done for the cardiac waveforms using the CPpd to produce CVT curves.
10 Cardiac rate time (CRT) courses based on the TTL (CRTttl)are generated as the
inverse change in time between CP centered points on the initial peak as in:
        Shmueli K, van Gelderen P, de Zwart JA, Horovitz SG, Fukunaga M,
        Jansma JM, Duyn JH.
        Low-frequency fluctuations in the cardiac rate as a source of variance in
        the resting-state fMRI BOLD signal.
        Neuroimage. 2007 Nov 1;38(2):306-20. Epub 2007 Aug 9.
        PMID: 17869543
  An RRT is also generated in the same manner using the RPpd.
11 Cardiac rate time (CRT) courses based on the third derivative (CRTd3)
 and the third derivative's R-wave estimate (CRTd3R) are generated as in:
        Chan GS, Middleton PM, Celler BG, Wang L, Lovell NH.
        Automatic detection of left ventricular ejection time from a finger
        photoplethysmographic pulse oximetry waveform: comparison with Doppler
        aortic measurement.
        Physiol Meas. 2007 Apr;28(4):439-52. Epub 2007 Mar 20.
        PMID: 17395998
12 PhysioNoise outputs the RW, CPd3, CPd3R, and CPttl
time courses as options for 3dretroicor and the RVT, RRT, CVT, CRTd3, and CRTd3R,
curves sampled on the TR and half TR as covariate options for neuroanalysis.

   usage: %prog [options] arg1 arg2

Example usage:
./JS_PhysioNoise_pylog.py -l logfile.txt 
./JS_PhysioNoise_pylog.py -l logfile.txt --numTR 360 --rmagthresh 80 --cmagthresh 10 -p PreIcor_Subj1 --speed

"""
epilogue ="""==========================================================================="""

parser = OptionParser(usage=usage, epilog=epilogue)
parser.add_option("-l", "--log",  action="store", type="string", dest="logfile",help="read  signals from FILE", metavar="FILE")
#parser.add_option("-e", "--ecg",  action="store", type="string", dest="ecgfile",help="read ecg signal from FILE", metavar="FILE")
#parser.add_option("-p", "--pulse",  action="store",type="string", dest="pulsefile",help="read pulse from FILE", metavar="FILE")
parser.add_option("--TR",  action="store",type="float", dest="TR",help="TR is TR seconds [2.0]", metavar="TR", default=2.0)
parser.add_option("--numTR",  action="store",type="int", dest="NumTR",help="Number of TRs [480]", metavar="TR", default=480)
parser.add_option("--inputHz",  action="store",type="int", dest="OrigSamp",help="input frequency in Hz[100]", metavar="Hz", default=100)
#parser.add_option("--outputHz",  action="store",type="int", dest="NewSamp",help="desired output frequency in Hz[40]", metavar="Hz", default=40)
parser.add_option("--fpass",  action="store",type="int", dest="fpass",help="Butterworth filter:Stop frequency of passband in Hz[1]", metavar="Hz", default=1)
parser.add_option("--fstop",  action="store",type="int", dest="fstop",help="Butterworth filter:Start frequency of the stopband  in Hz[10]", metavar="Hz", default=10)
parser.add_option("--trimwindow",  action="store",type="int", dest="TrimWindow",help="Number of points PTS to clean on either side of dirty residual peaks[20]", metavar="PTS", default=20)
parser.add_option("--plr",  action="store",type="float", dest="FwdRm",help="Percentage of points on the left of dirty peaks to clean on the right [1.2]", metavar="Percentage", default=1.2)
parser.add_option("--statthresh",  action="store",type="int", dest="FoldSTD",help="Mark points beyond FOLD standard deviation(s) as dirty residual peaks for cleaning  [6]", metavar="FOLD", default=6)
parser.add_option("--rmagthresh",  action="store",type="int", dest="RMagThresh",help="Respiratory Peak must be RMAG units away from nearest extrema  [100]", metavar="RMAG", default=100)
parser.add_option("--rtimethresh",  action="store",type="float", dest="RTimeThresh",help="Respratory Peak must be RTIME seconds away from the nearest extrema  [0]", metavar="RTIME", default=0)
parser.add_option("--cmagthresh",  action="store",type="int", dest="PulseMagThresh",help="Pulse Peak must be CMAG units away from nearest extrema  [20]", metavar="CMAG", default=20)
parser.add_option("--ctimethresh",  action="store",type="float", dest="PulseTimeThresh",help="Pulse Peak must be CTIME seconds away from the nearest extrema  [0]", metavar="CTIME", default=0)
parser.add_option("--plots",  action="store_true",dest="PlotAll",help="this turns on interactive plotting", default=False)
parser.add_option("--speed",  action="store_true",dest="Speed",help="This skips artifact detection and spline interpolation. Instead a cubic spline filter is used and does a decent job removing artifact. [artifact detect]", default=False)
parser.add_option("-o","--prefix",  action="store", type="string",dest="Prefix",help="Output prefix [PreIcor]", default='PreIcor')

options, args = parser.parse_args()

if '-help' in sys.argv:
        parser.print_help()
        raise SystemExit()

if not options.logfile  or '-help' in sys.argv:
        print "At least one physiological data file is required. Try --help "
        raise SystemExit()
        
# Daten Einlesen
if options.logfile:
    try: Logfile=fromfile(options.logfile,'f', sep=' ' )
    except IOError: sys.exit("Respiration File not found")
 
  # preprocess
    Logfile = reshape(Logfile, (4, -1))
    IResp = Logfile[0, ]
    IPulse = -Logfile[1, ]
    
    IResp=IResp  * 10000  # Da Kalkulation teilweise in INT skalieren
    IResp=detrend_linear(IResp) # linear trendremoval, entfernt auch mean
    
    IPulse=(IPulse - mean(IPulse)) * 10000
    
####
PlotAll=options.PlotAll
OrigSamp=options.OrigSamp
fpass=options.fpass #passband edge frequency
fstop=options.fstop  #stopband edge frequency
Nyquist=long(OrigSamp/2)
Speed=options.Speed
RespTime = range(len(IResp)) # x Achse in Punkten
PulseTime = range(len(IPulse))
RespMagThresh=options.RMagThresh
RespTimeThresh=options.RTimeThresh * OrigSamp
PulseMagThresh=options.PulseMagThresh 
PulseTimeThresh=options.PulseTimeThresh * OrigSamp
TotSamp = int(options.TR * options.NumTR * OrigSamp) # Totale Anzahl an Samples

####--------------------------------Butter Filter  Data
# Respiration

filter_b,filter_a=scipy.signal.iirdesign(float(fpass)/Nyquist,float(fstop)/Nyquist,gpass=5,gstop=60,analog=0, ftype="butter")

RespFilt=filtfilt(filter_b, filter_a, IResp)

w, h = freqz(filter_b,filter_a)
mag = abs(h)
phase = unwrap(angle(h))*180/pi
radsample=w/pi

x= int((log(OrigSamp)-log(0.01))/log(2))
NFFT=int(2**x)
nooverlap=int(0.5*(2**x))
pxx, freqs = psd( RespFilt, NFFT=NFFT , Fs=OrigSamp,noverlap=nooverlap , detrend=detrend_linear, window=window_hanning)
maxfreqpt = pxx.argmax()
maxfreq =freqs[maxfreqpt]
PeakPower=pxx[maxfreqpt]

if maxfreq==0:
        print "Respiratory power peak at zero. Threshold at",freqs[4],"Hz"
        Tempxx=pxx[5:]
        Tempfreqs=freqs[5:]
        maxfreqpt = Tempxx.argmax()
        maxfreq =Tempfreqs[maxfreqpt]
        PeakPower=Tempxx[maxfreqpt]

print "The spectral respiratory peak is:", maxfreq,"Hz "
print "The respiratory peak power density is",PeakPower
print "The spectral respiratory period is", 1./maxfreq,"seconds"

#Pulse

Pulsefilter_b,Pulsefilter_a=scipy.signal.iirdesign(float(fpass)/Nyquist,float(fstop)/Nyquist,gpass=5,gstop=60,analog=0, ftype="butter")
PulseFilt=filtfilt(Pulsefilter_b, Pulsefilter_a, IPulse)
Pulsepxx, Pulsefreqs = psd( PulseFilt, NFFT=NFFT , Fs=OrigSamp,noverlap=nooverlap , detrend=detrend_linear, window=window_hanning)
Pulsemaxfreqpt = Pulsepxx.argmax()
Pulsemaxfreq =Pulsefreqs[Pulsemaxfreqpt]
PulsePeakPower=Pulsepxx[Pulsemaxfreqpt]
if Pulsemaxfreq==0:
        print "Pulse power peak at zero. Threshold at",freqs[4],"Hz"
        PulseTempxx=Pulsepxx[5:]
        PulseTempfreqs=Pulsefreqs[5:]
        Pulsemaxfreqpt = PulseTempxx.argmax()
        Pulsemaxfreq =PulseTempfreqs [Pulsmaxfreqpt]
        PulsePeakPower=PulseTempxx[Pulsmaxfreqpt]

print "The spectral pulse peak is:", Pulsemaxfreq,"Hz "
print "The pulse peak power density is",PulsePeakPower
print "The spectral pulse period is", 1./Pulsemaxfreq,"seconds"

#---------------------------------------------------#Remove regions of high noise variance for interpolation

RespResid = IResp - RespFilt
PulseResid = IPulse - PulseFilt

if Speed:
        print "Spline Filtering"
        print "Respiratory"
        RespBestFit = spline_filter1d(RespFilt, order=3)
        print "Cardiac"
        PulseBestFit= spline_filter1d(PulseFilt, order=5)
else:
        print "Finding the regions of high variance"

        print "Start Respiration"
        FoldSTD=options.FoldSTD
        OutlierThresh =float(FoldSTD)*std(RespResid)
        RemoveRange= options.TrimWindow
        ExtremTime=0
        FwdRm=options.FwdRm
        RespAllOutliers=zeros(len(RespFilt))
        DirtyPoints= []
        SigThreshResid = abs(RespResid) - OutlierThresh
        idUglyPoints= nonzero(greater_equal(SigThreshResid,0))[0]
        for counter in idUglyPoints:
                RespAllOutliers[counter]=RespResid[counter]
                for i in (range(long( counter -RemoveRange),long(counter+FwdRm*RemoveRange))) :
                        DirtyPoints.extend([i])
        DirtyPoints= sort(unique(clip(unique(DirtyPoints),1,len(RespResid))))
        NoOutliers=zeros(len(RespFilt))
        NoOutlierTime=zeros(len(RespTime))
        NoOutliers=delete(RespFilt,DirtyPoints)
        NoOutlierTime=delete(RespTime,DirtyPoints)
        print("Respiration Outliers",  len(DirtyPoints))
        print "Spline interpolation over deleted acquisition artifact"
        print "Respiratory"
        SplineEst = scipy.interpolate.splrep(NoOutlierTime ,NoOutliers, k=3)
        RespBestFit = scipy.interpolate.splev(RespTime  , SplineEst)

        #Cardiac
        print "Start Pulse"
#        FoldSTD=options.FoldSTD
        OutlierThresh =float(FoldSTD)*std(PulseResid)
        RemoveRange= options.TrimWindow
        ExtremTime=0
        FwdRm=options.FwdRm
        PulseAllOutliers=zeros(len(PulseFilt))
        DirtyPoints= []
        SigThreshResid = abs(PulseResid) - OutlierThresh
        idUglyPoints= nonzero(greater_equal(SigThreshResid,0))[0]
        for counter in idUglyPoints:
                PulseAllOutliers[counter]=PulseResid[counter]
                for i in (range(long( counter -RemoveRange),long(counter+FwdRm*RemoveRange))) :
                        DirtyPoints.extend([i])
        DirtyPoints= sort(unique(clip(unique(DirtyPoints),1,len(PulseResid))))
        PulseNoOutliers=zeros(len(PulseFilt))
        PulseNoOutlierTime=zeros(len(PulseTime))
        PulseNoOutliers=delete(PulseFilt,DirtyPoints)
        PulseNoOutlierTime=delete(PulseTime,DirtyPoints)
        print("Pulse Outliers",  len(DirtyPoints))
        print "Pulse interpolation"
        PulseSplineEst = scipy.interpolate.splrep(PulseNoOutlierTime ,PulseNoOutliers, k=5)
        PulseBestFit = scipy.interpolate.splev(PulseTime  , PulseSplineEst)

# -------------------- Data is now clean --------------------------------
RespAllExtrema, RespAllExtremaTime,RespAllHills, RespAllValleys, RespTopEnv, RespBotEnv, RespTopCounts, RespBotCounts=findextrema(RespBestFit,RespMagThresh, RespTimeThresh)
PulseAllExtrema, PulseAllExtremaTime,PulseAllHills, PulseAllValleys, PulseTopEnv, PulseBotEnv, PulseTopCounts, PulseBotCounts=findextrema(PulseBestFit,PulseMagThresh, PulseTimeThresh)

#-----------------------------------------------------------------------------Find Volume over Time
print"Calculating Envelopes"
print "Respiratory"
RespTimeEnv=ConnectTheDots(RespTopCounts,diff(RespTopCounts) ,RespBestFit,1)
RespInterpTopEnv=ConnectTheDots(RespTopCounts,RespAllHills,RespBestFit,0)
RespInterpBotEnv=ConnectTheDots(RespBotCounts,RespAllValleys,RespBestFit,0)
RespRVTEnv=zeros(len(RespBestFit), 'f')
RespRVTEnv = OrigSamp*abs(RespInterpTopEnv-RespInterpBotEnv)/ (RespTimeEnv)

RespMeanRVT=zeros(len(RespBestFit), 'f')
RespMeanRVT = (    abs(RespInterpTopEnv)+   abs(RespInterpBotEnv))/ 2.

print "Pulse"
PulseTimeEnv=ConnectTheDots(PulseTopCounts,diff(PulseTopCounts) ,PulseBestFit,1)
PulseInterpTopEnv=ConnectTheDots(PulseTopCounts,PulseAllHills,PulseBestFit,0)
PulseInterpBotEnv=ConnectTheDots(PulseBotCounts,PulseAllValleys,PulseBestFit,0)
PulseRVTEnv=zeros(len(PulseBestFit), 'f')
PulseRVTEnv = OrigSamp*abs(PulseInterpTopEnv-PulseInterpBotEnv)/ (PulseTimeEnv)

PulseMeanRVT=zeros(len(PulseBestFit), 'f')
PulseMeanRVT = (abs(PulseInterpTopEnv)+abs(PulseInterpBotEnv))/ 2.

#-----------------------------Find third derivatives of Cardiac signal to create new trigger file
print"Finding Pulse Derivatives"
dPulse1=slopes(PulseTime,PulseBestFit)    * OrigSamp
dPulse2=slopes(PulseTime,dPulse1)   * OrigSamp
dPulse3=slopes(PulseTime,dPulse2)  * OrigSamp
#dPulseSmooth3=cmov_average(dPulse3,int(0.05*OrigSamp))
dPulseSmooth3=movingaverage(dPulse3,int(0.05*OrigSamp))
#Find large Troughs
dPulseAllExtrema, dPulseAllExtremaTime,dPulseAllHills, dPulseAllValleys, dPulseTopEnv, dPulseBotEnv, dPulseTopCounts, dPulseBotCounts=findextrema( dPulseSmooth3,PulseMagThresh, PulseTimeThresh)

idRise=nonzero(less_equal(dPulse1 ,0 ))
put(dPulseBotEnv,idRise,0)
CPd3=-1.*dPulseBotEnv
#CPd3 has location of the large deriv3 troughs
#use envelope of deriv3 peaks as a threshold
CPd3Threshd3 =ConnectTheDots(dPulseTopCounts,abs(dPulseAllHills),PulseBestFit,0)

print "Thresholding CPd3"

idRise=nonzero(less_equal(CPd3-CPd3Threshd3,0))[0]
#put(CPd3,idRise,0)
idRise=nonzero(less_equal(CPd3,0.5*mean(CPd3Threshd3)))[0]
#put(CPd3,idRise,0)

# now find R waves which are the deriv3 peaks that precede the trough which is in CPd3
CPd3Rwave=zeros(len(PulseBestFit), 'f')
CRwaveCounts=nonzero(greater(CPd3,0))[0]
for peaktime in CRwaveCounts:
    tempR=nonzero(less(dPulseTopCounts,peaktime))[0]
    CPd3Rwave[int(dPulseTopCounts[int(tempR[-1])])]= CPd3[int(peaktime)]

print "Length CPd3Rwave: " + str(len(CPd3Rwave))

#------------------Calculate respiratory (RRT) and cardiac rate per time (CRT)
print "Calculating Rates over Time"
#Calculate RRT
RespTopEnvTime=nonzero(greater(RespTopEnv,0))[0]
RespRate=OrigSamp * (diff(RespTopEnvTime)**(-1.))
RespTTopCounts=RespTopEnvTime[0:-1]
RRT=ConnectTheDots(RespTTopCounts,RespRate,RespBestFit,0)

#Calculate CRTd3 = Pulse rate von CPd3
PulseTopEnvTime=nonzero(greater(CPd3,0))[0]
PulseRate =OrigSamp * (diff(PulseTopEnvTime)**(-1.))
PulseTTopCounts=PulseTopEnvTime[0:-1]
CRTd3=ConnectTheDots(PulseTTopCounts,PulseRate,PulseBestFit,0)

#Calculate CRTd3R

RWaveTime=nonzero(greater(CPd3Rwave,0))[0]
RRate=OrigSamp * (diff(RWaveTime)**(-1.))
RWaveTopCounts=RWaveTime[ 0:-1]
CRTd3R=ConnectTheDots(RWaveTopCounts,RRate,PulseBestFit,0)

#--------------------------------------------#Plotten

if PlotAll:
        figure(1,figsize=(16, 12))
        subplot(311)
        semilogy(radsample*Nyquist,mag)
        ylabel('Magnitude')
        title('Bode Plot ' + os.path.basename(options.Prefix))
        pylab.ylim((10e-4, 10e0))
        pylab.xlim((0.0, fstop))

        figure(1)
        subplot(312)
        plot(freqs,pxx)
        pylab.ylabel('Respiratory Power')
        pylab.xlim((0.0, fstop))

        figure(1)
        subplot(313)
        plot(Pulsefreqs,Pulsepxx)
        pylab.xlabel('frequency(Hz)')
        pylab.ylabel('Pulse Power')
        pylab.xlim((0.0, fstop))
        
        savefig(options.Prefix + '_Fig_1.png',orientation='landscape')
        
        figure(2,figsize=(16, 12))
        subplot(211)
        plot(IResp, 'y', label='raw')
        plot(RespFilt, 'b', label='filtered')
        plot(RespResid, 'r',  label = 'residual')
        plot(RespBestFit, 'k', label='Spline')
        if not Speed:
            plot(RespTime, RespAllOutliers, 'go', label='outlier')
        legend()
        title('Respiratory ' + os.path.basename(options.Prefix))

    
        subplot(212)
        plot(IPulse, 'y',  label='raw')
        plot(PulseFilt, 'b', label='filtered')
        plot(PulseResid, 'r', label='residual')
        plot(PulseBestFit, 'k',  label='Spline')
        if not Speed:
            plot(PulseTime, PulseAllOutliers, 'go', label='outlier')
        legend()
        title('Puls')
        savefig(options.Prefix + '_Fig_2.png',orientation='landscape')

        figure(3,figsize=(16, 12))
        subplot(211)
        plot(RespBestFit, 'b', label='Spline')
        plot(RespTime, RespInterpTopEnv, 'r',  label='TopEnv')
        plot(RespTime, RespInterpBotEnv, 'g',  label='BotEnv')
        plot(RespTime, RespTopEnv, 'ro')
        plot(RespTime, RespBotEnv, 'go')
        plot(RespTime, RespRVTEnv, 'c', label='RVT')
        plot(RespTime, RespMeanRVT, 'y', label='MeanAbsEnv')
        legend()
        title('Respiratory data ' + os.path.basename(options.Prefix))
    
    
        subplot(212)
        plot(PulseBestFit, 'b', label='Spline')
        plot(PulseTime, PulseInterpTopEnv, 'r',  label='TopEnv')
        plot(PulseTime, PulseInterpBotEnv, 'g',  label='BotEnv')
        plot(PulseTime, PulseTopEnv, 'ro')
        plot(PulseTime, PulseBotEnv, 'go')
        plot(PulseTime, PulseRVTEnv, 'c', label='CVT')
        plot(PulseTime, PulseMeanRVT, 'y', label='MeanAbsEnv')
        legend()
        title('Pulse data')
        savefig(options.Prefix + '_Fig_3.png',orientation='landscape')
        
        figure(4,figsize=(16, 12))
        plot(PulseTime, PulseBestFit, 'k', label='Spline')
#            plot(PulseTime, dPulse1, 'b', label='deriv1')
#            plot(PulseTime, dPulse2, 'g', label='deriv2')
        plot(PulseTime, dPulse3, 'g', label='deriv3')
        plot(PulseTime, CPd3, 'b', label='CPd3')
        plot(PulseTime, CPd3Threshd3, 'y',  label='Threshold')
        plot(PulseTime,  CPd3Rwave, 'm', label='R Wave')
        legend()
        title(os.path.basename(options.Prefix))
        savefig(options.Prefix + '_Fig_4.png',orientation='landscape')        
        
        figure(333,figsize=(16, 12))
        subplot(321)
        plot(RRT)
        ylabel('RRT')
        title(os.path.basename(options.Prefix))
        subplot(322)
        plot(RespRVTEnv)
        ylabel('RVT')
        subplot(323)
        plot(CRTd3)
        ylabel('CRTd3')
        subplot(325)
        plot(CRTd3R)
        ylabel('CRTd3R')
        subplot(324)
        plot(PulseRVTEnv)
        ylabel('CVT')
        savefig(options.Prefix + '_Fig_333.png',orientation='landscape')
        

#-------------- calculate Fourier series-----------------
# difference between two R waves
ind=diff(nonzero(greater(CPd3Rwave, 0)))[0]
start=nonzero(greater(CPd3Rwave, 0)) # detect the first r wave
CardPhase=zeros(start[0][0])
print "CardPhase: " + str(CardPhase)
print "Length Pulse Phase: " + str(len(CardPhase))
for i in range(size(ind)):
    CardPhase=append(CardPhase, arange(ind[i])*2*pi/ind[i], axis=0)

print "Length Pulse Phase: " + str(len(CardPhase))
print "ind " + str(size(ind))
print "start " + str(size(start))

#-------- respiration
RespScale=RespBestFit + abs( min(RespBestFit))
HistSum=cumsum(histogram(RespScale,bins=100)[0])
ind=99*RespScale/max(RespScale) # scale Resp to 0-99
ind=ind.astype(int) # cast as int
# --- diff() hat 1 Datenpunkt weniger als ind !!! Daher ind[0:-1]
RespPhase = pi * HistSum[ind[0:-1]]/float(HistSum[-1]) * sign (diff(RespScale))

print "Length Resp Phase: " + str(len(RespPhase))

#------------------------- Trimmen fuer RETROICOR
if ( TotSamp > len(RespTime)):
    print "Time series is too short based on TR"
    print "Time series should be:", TotSamp
    print "Time series actually is:", len(RespTime)
    print "Was the TR actually",options.TR * len(RespTime)*(TotSamp)**-1,"not",options.TR,"seconds?"
#CardPhase=CardPhase[0:TotSamp]
#RespPhase=RespPhase[0:TotSamp]
RespBestFit = RespBestFit[0:TotSamp]
RespRVTEnv = RespRVTEnv[0:TotSamp]
PulseBestFit = PulseBestFit[0:TotSamp]
PulseRVTEnv = PulseRVTEnv[0:TotSamp]
CRTd3 = CRTd3[0:TotSamp]
CPd3 = CPd3[0:TotSamp]
CPd3Rwave = CPd3Rwave[0:TotSamp]

#---------------Calculating Stats
print "Calculating Statistics for", options.Prefix
OutputStat= zeros(54,'f')

OutputStat[0]=sum(greater(RespTopEnv ,0 ))
print  "Peakdet detected ",OutputStat[0]  ," positive respiratory peaks"
OutputStat[1]=( mean(diff(nonzero(greater( RespTopEnv ,0 ))))  )
print  "The peakdet average respiratory interpeak time was",  OutputStat[1]
OutputStat[2]=( std(diff(nonzero(greater( RespTopEnv ,0 ))))   )
print  "The peakdet std respiratory interpeak time was",  OutputStat[2]

OutputStat[3]=( maxfreq )
print  "The spectral respiratory peak is:",OutputStat[3]  ,"Hz "
OutputStat[4]=(PeakPower  )
print  "The respiratory peak power density is", OutputStat[4]
OutputStat[5]=(1./maxfreq  )
print  "The spectral respiratory period is",OutputStat[5]  ,"seconds"

OutputStat[6]=( Pulsemaxfreq )
print  "The spectral pulse peak is:",OutputStat[6]  ,"Hz "
OutputStat[7]=(PulsePeakPower  )
print  "The pules peak power density is", OutputStat[7]
OutputStat[8]=( 1./Pulsemaxfreq )
print  "The spectral pulse period is",OutputStat[8]  ,"seconds"

OutputStat[9]=sum(greater(CPd3Rwave ,0 ))
print  "CPd3R detected ",OutputStat[9]  ," positive cardiac peaks"
OutputStat[10]=( mean(diff(nonzero(greater( CPd3Rwave ,0 ))))  )
print  "The average cardiac CPd3R interpeak time was",  OutputStat[10]
OutputStat[11]=( std(diff(nonzero(greater( CPd3Rwave ,0 ))))  )
print  "The std  cardiac CPd3R interpeak time was", OutputStat[11]

OutputStat[12]=( std(CRTd3) )
print  "The CRTd3 std is:",  OutputStat[12]
OutputStat[13]=( std(CRTd3R) )
print  "The CRTd3R std is:",  OutputStat[13]

OutputStat[14]=( mean(CRTd3) )
print    "The CRTd3 mean is:",  OutputStat[14]
OutputStat[15]=( mean(CRTd3R)  )
print    "The CRTd3R mean is:",   OutputStat[15]

OutputStat[16]=( sum(greater(CPd3 ,0 )))
print  "CPd3 detected ", OutputStat[16] ," positive cardiac peaks"
OutputStat[17]=( mean(diff(nonzero(greater( CPd3 ,0 ))))  )
print  "The average cardiac CPd3 interpeak time was", OutputStat[17]
OutputStat[18]=( std(diff(nonzero(greater( CPd3 ,0 ))))  )
print  "The std  cardiac CPd3 interpeak time was", OutputStat[18]
OutputStat[19]=( mean(diff(nonzero(greater( CPd3Rwave ,0 ))))  )
print  "The avg  cardiac CPd3R peak time was", OutputStat[19]
OutputStat[20]=( std((nonzero(greater( CPd3 ,0 ))))  )
print  "The std  cardiac CPd3 peak time was", OutputStat[20]
OutputStat[21]=( std((nonzero(greater( CPd3Rwave ,0 ))))  )
print  "The std  cardiac CPd3R peak time was", OutputStat[21]

OutputStat[37]= sum(greater( CRTd3 ,mean(CRTd3)+3*std(CRTd3) ))
print  "The number of CRTd3 peaks > 3 std",OutputStat[37]

OutputStat[38]= len(nonzero(greater( CRTd3R ,mean(CRTd3R)+3*std(CRTd3R) )))
print  "The number of CRTd3R peaks > 3 std",OutputStat[38]

OutputStat[39]= len(RespTime)
print  "Length of Respiration Signal",OutputStat[39]

OutputStat[40]= len(PulseTime)
print  "Length of Pulse Signal",OutputStat[40]

OutputStat[41]= std(CRTd3[nonzero(greater( CRTd3 ,mean(CRTd3)+3*std(CRTd3) ))])
print  "The std of CRTd3 peaks > 3 std",OutputStat[41]

OutputStat[42]= std(CRTd3R[nonzero(greater( CRTd3R ,mean(CRTd3R)+3*std(CRTd3R) ))])
print  "The std of CRTd3R peaks > 3 std",OutputStat[42]

OutputStat[43]= sum(greater( CRTd3 ,mean(CRTd3)+1. ))
print  "The number of CRTd3 peaks > mean CRTd3 + 1 ",OutputStat[43]

OutputStat[44]= len(nonzero(greater( CRTd3R ,mean(CRTd3R)+1.)))
print  "The number of CRTd3R peaks > mean CRTd3R + 1 ",OutputStat[44]

OutputStat[45]= sum(greater( CRTd3 ,mean(CRTd3)-0.4))
print  "The number of CRTd3 peaks 0.4 Hz below mean ",OutputStat[45]

OutputStat[46]= len(nonzero(greater( CRTd3R ,mean(CRTd3R)-0.4)))
print  "The number of CRTd3R peaks 0.4 Hz below mean ",OutputStat[46]

OutputStat[47]=sum(greater(PulseTopEnv ,0 ))
print  "Peakdet detected ",OutputStat[47]  ," positive cardiac peaks"
OutputStat[48]=( mean(diff(nonzero(greater( PulseTopEnv ,0 ))))  )
print  "The  peakdet average cardiac interpeak time was",  OutputStat[48]
OutputStat[49]=( std(diff(nonzero(greater( PulseTopEnv ,0 ))))   )
print  "The peakdet std cardiac interpeak time was",  OutputStat[49]

OutputStat[50]= sum(greater( CRTd3 ,mean(CRTd3)+4*std(CRTd3) ))
print  "The number of CRTd3 peaks > 4 std",OutputStat[50]

OutputStat[51]= len(nonzero(greater( CRTd3R ,mean(CRTd3R)+4*std(CRTd3R) )))
print  "The number of CRTd3R peaks > 4 std",OutputStat[51]

OutputStat[52]= std(CRTd3[nonzero(greater( CRTd3 ,mean(CRTd3)+4*std(CRTd3) ))])
print  "The std of CRTd3 peaks > 4 std",OutputStat[52]

OutputStat[53]= std(CRTd3R[nonzero(greater( CRTd3R ,mean(CRTd3R)+4*std(CRTd3R) ))])
print  "The std of CRTd3R peaks > 4 std",OutputStat[53]
print OutputStat
"-----------------"

print "Saving files----------------------"

#-----------------------------------------------save files to disk

filename = options.Prefix +'_OutputStatistics'+".txt"
savetxt(filename, OutputStat, fmt='%8.6f')

#Respiration


filename = options.Prefix +'_Resp_Spline_'+ str(OrigSamp) + ".txt"
savetxt(filename, RespBestFit + abs( min(RespBestFit)), fmt='%8.6f')

filename = options.Prefix +'_RespPhase_'+ str(OrigSamp) + ".txt"
savetxt(filename, RespPhase, fmt='%8.6f')

#filename = options.Prefix +'_Resp_cos1_'+ str(OrigSamp) + ".txt"
#savetxt(filename, cos(RespPhase), fmt='%8.6f')
#
#filename = options.Prefix +'_Resp_cos2_'+ str(OrigSamp) + ".txt"
#savetxt(filename, cos(2*RespPhase), fmt='%8.6f')
#
#filename = options.Prefix +'_Resp_cos3_'+ str(OrigSamp) + ".txt"
#savetxt(filename, cos(3*RespPhase), fmt='%8.6f')
#
#filename = options.Prefix +'_Resp_sin1_'+ str(OrigSamp) + ".txt"
#savetxt(filename, sin(RespPhase), fmt='%8.6f')
#
#filename = options.Prefix +'_Resp_sin2_'+ str(OrigSamp) + ".txt"
#savetxt(filename, sin(2*RespPhase), fmt='%8.6f')
#
#filename = options.Prefix +'_Resp_sin3_'+ str(OrigSamp) + ".txt"
#savetxt(filename, sin(3*RespPhase), fmt='%8.6f')


#SigTRfull,SigTRPt, SigTRPtTime=TRSample(CleanSig,1.0,OrigSamp)
#filename = options.Prefix +'_TR_RW_Spline_'+ str(OrigSamp) + ".txt"
#savetxt(filename,  SigTRPt  + abs(min(SigTRPt)), fmt='%8.6f')
#
#HSigTRfull,HSigTRPt, HSigTRPtTime=TRSample(CleanSig,0.5,OrigSamp)
#filename = options.Prefix +'_HalfTR_RW_Spline_'+ str(OrigSamp) + ".txt"
#savetxt(filename,  HSigTRPt  + abs(min(HSigTRPt)), fmt='%8.6f'  )

#RVT
filename = options.Prefix +'_RVT_'+ str(OrigSamp) + ".txt"
savetxt(filename,  RespRVTEnv  + abs(min(RespRVTEnv)), fmt='%8.6f' )

#SigTRfull,SigTRPt, SigTRPtTime=TRSample( FullInterpRVTEnv,1.0,OrigSamp)
#filename = options.Prefix +'_TR_RVT_'+ str(OrigSamp) + ".txt"
#savetxt(filename,  SigTRPt  + abs(min(SigTRPt)), fmt='%8.6f')
#
#HSigTRfull,HSigTRPt, HSigTRPtTime=TRSample(FullInterpRVTEnv ,0.5,OrigSamp)
#filename = options.Prefix +'_HalfTR_RVT_'+ str(OrigSamp) + ".txt"
#savetxt(filename,  HSigTRPt  + abs(min(HSigTRPt)), fmt='%8.6f'  )

#RRT

filename = options.Prefix +'_RRT_'+ str(OrigSamp) + ".txt"
savetxt(filename, RRT  + abs(min(RRT)), fmt='%8.6f' )

#SigTRfull,SigTRPt, SigTRPtTime=TRSample( RRT,1.0,OrigSamp)
#filename = options.Prefix +'_TR_RRT_'+ str(OrigSamp) + ".txt"
#savetxt(filename,  SigTRPt  + abs(min(SigTRPt)), fmt='%8.6f')
#
#HSigTRfull,HSigTRPt, HSigTRPtTime=TRSample(RRT ,0.5,OrigSamp)
#filename = options.Prefix +'_HalfTR_RRT_'+ str(OrigSamp) + ".txt"
#savetxt(filename,  HSigTRPt  + abs(min(HSigTRPt)), fmt='%8.6f'  )

#Cardiac Peaks
filename = options.Prefix +'_CPd3_'+ str(OrigSamp) + ".txt"
savetxt(filename,  CPd3  , fmt='%8.6f' )

filename = options.Prefix +'_CPd3R_'+ str(OrigSamp) + ".txt"
savetxt(filename,  CPd3Rwave  , fmt='%8.6f' )


filename = options.Prefix +'_CardPhase_'+ str(OrigSamp) + ".txt"
savetxt(filename, CardPhase, fmt='%8.6f')

#filename = options.Prefix +'_Card_cos1_'+ str(OrigSamp) + ".txt"
#savetxt(filename, cos(CardPhase), fmt='%8.6f')
#
#filename = options.Prefix +'_Card_cos2_'+ str(OrigSamp) + ".txt"
#savetxt(filename, cos(2*CardPhase), fmt='%8.6f')
#
#filename = options.Prefix +'_Card_cos3_'+ str(OrigSamp) + ".txt"
#savetxt(filename, cos(3*CardPhase), fmt='%8.6f')
#
#filename = options.Prefix +'_Card_sin1_'+ str(OrigSamp) + ".txt"
#savetxt(filename, sin(CardPhase), fmt='%8.6f')
#
#filename = options.Prefix +'_Card_sin2_'+ str(OrigSamp) + ".txt"
#savetxt(filename, sin(2*CardPhase), fmt='%8.6f')
#
#filename = options.Prefix +'_Card_sin3_'+ str(OrigSamp) + ".txt"
#savetxt(filename, sin(3*CardPhase), fmt='%8.6f')

#CRTd3 trough

filename = options.Prefix +'_CRTd3_'+ str(OrigSamp) + ".txt"
savetxt(filename, CRTd3  + abs(min(CRTd3)), fmt='%8.6f' )

#SigTRfull,SigTRPt, SigTRPtTime=TRSample( CRTd3,1.0,OrigSamp)
#filename = options.Prefix +'_TR_CRTd3_'+ str(OrigSamp) + ".txt"
#savetxt(filename,  SigTRPt  + abs(min(SigTRPt)), fmt='%8.6f')
#
#HSigTRfull,HSigTRPt, HSigTRPtTime=TRSample(CRTd3 ,0.5,OrigSamp)
#filename = options.Prefix +'_HalfTR_CRTd3_'+ str(OrigSamp) + ".txt"
#savetxt(filename,  HSigTRPt  + abs(min(HSigTRPt)), fmt='%8.6f'  )

#CRTd3Rwave
filename = options.Prefix +'_CRTd3R_'+ str(OrigSamp) + ".txt"
savetxt(filename, CRTd3R  + abs(min(CRTd3R)), fmt='%8.6f' )

#SigTRfull,SigTRPt, SigTRPtTime=TRSample( CRTd3R,1.0,OrigSamp)
#filename = options.Prefix +'_TR_CRTd3R_'+ str(OrigSamp) + ".txt"
#savetxt(filename,  SigTRPt  + abs(min(SigTRPt)), fmt='%8.6f')

#HSigTRfull,HSigTRPt, HSigTRPtTime=TRSample(CRTd3R ,0.5,OrigSamp)
#filename = options.Prefix +'_HalfTR_CRTd3R_'+ str(OrigSamp) + ".txt"
#savetxt(filename,  HSigTRPt  + abs(min(HSigTRPt)), fmt='%8.6f'  )
print "Done"
if PlotAll:
    show()


