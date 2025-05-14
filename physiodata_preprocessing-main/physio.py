#!/usr/local/bin/python
#from __future__ import division

"""Author: Dan Kelley, University of Wisconsin - Madison
   Program: PhysioNoise.py
   Date: November 2007
"""


import string, sys
from numpy import *
import scipy as Sci
import scipy.linalg
import pytz
from pylab import *
from matplotlib import *
from scipy import linspace, polyval, polyfit, sqrt, stats, randn,signal, interpolate 
from pylab import plot, title, show , legend, figure
from scipy.signal import cspline1d, cspline1d_eval, iirdesign, lfilter, lp2lp, freqz
#import scipy.interpolate as interpolate
from scipy.ndimage import spline_filter1d 
from timeseries.lib.moving_funcs import cmov_average
from matplotlib.mlab import psd
from scipy.signal import cspline1d, cspline1d_eval
from scipy.signal.signaltools import hilbert
from matplotlib.mlab import slopes

import matplotlib.numerix.ma as M
######------------------------------------------------------------------------####### Find Extrema

#Modified from peakdet.m by Eli Billauer, 3.4.05 (Explicitly not copyrighted).
#This function is released to the public domain; Any use is allowed.
#Accessed Sept 2007 

def findextrema(TS, ChangeMag,TimeConstraint):
	MinMag = Inf
	MaxMag = -Inf
	FindMaxMag = 1
	MaxCounts=zeros(len(TS))
	MinCounts=zeros(len(TS))
	MaxCountsPeak=zeros(len(TS),Float)
	MinCountsPeak=zeros(len(TS),Float)
	maxcounter=0
	mincounter=0
	maxmincounter=0
	AllExtremaPeak=zeros(len(TS),Float)
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
		TopEnvTemp[topms] = MaxCountsPeak[counter]

	for counter in xrange(len(MinCounts)):
		botms = MinCounts[counter]
		BotEnvTemp[botms] = MinCountsPeak[counter]

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

#---------------------------------------------------------------------------Downsample			 


def DownSample(TS):
	SigDownsampleTime=range(long(len(TS)/SampRatio))
	SigDownsample=zeros(len(SigDownsampleTime ),'f')
	SigDownsampleFull=zeros(len(TS),'f')
	SigDownsampleFull[0]=TS[0]
	SigDownsample[0]=TS[0]
	SigDownsampleTime[0]=0
	for counter in range(1,len(SigDownsampleTime)) :
		ms =  (counter*SampRatio)-1
		SigDownsampleFull[ms] = TS[ms]
		SigDownsample[counter] = TS[ms]
		SigDownsampleTime[counter] = ms
	return SigDownsampleFull, SigDownsample, SigDownsampleTime
#---------------------------------------------------------------------------TRSample			 
def TRSample(TS,TRFreq, SampRate):
	SigDownsampleTime=range(long(NumTR/float(TRFreq)))
	SigDownsample=zeros(len(SigDownsampleTime ),'f')
	SigDownsampleFull=zeros(len(TS),'f')
	SigDownsampleFull[0]=TS[0]
	SigDownsample[0]=TS[0]
	SigDownsampleTime[0]=0
#	print "SigDownsampleTime", len(SigDownsampleTime)
#	print "NumTR", NumTR, "TRFreq", TRFreq, "No div T",long(NumTR/float(TRFreq)) 
#	print "length TS",len(TS)
	if (len(SigDownsampleTime)*TR *TRFreq *SampRate ) > len(TS):
		print "Time series is too short based on TR"
		print "Time series should be:", len(SigDownsampleTime)*TR *TRFreq *SampRate 
		print "Time series actually is:", len(TS)
		print "Was the TR actually",len(TS)*(len(SigDownsampleTime)*TRFreq *SampRate)**-1,"not",TR,"seconds?"
		return SigDownsampleFull, SigDownsample, SigDownsampleTime
	for counter in range(1,len(SigDownsampleTime)) :
		ms =  int(counter*TR *TRFreq *SampRate)-1
#		print "counter", counter, "ms", ms
		SigDownsampleFull[ms] = TS[ms]
		SigDownsample[counter] = TS[ms]
		SigDownsampleTime[counter] = ms
	return SigDownsampleFull, SigDownsample, SigDownsampleTime

#---------------------------------------------------------------------------End TRSample	
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
		FullInterpEnv[Counts[RepPt] + i ]= InterpEnv[i]
	FullInterpEnv[0:Counts[RepPt]]=FullInterpEnv[Counts[RepPt]]
	FullInterpEnv[Counts[-1]:]=FullInterpEnv[Counts[-1]]
	return FullInterpEnv
#------------------------------------------------------EndConnectTheDots

		 

from optparse import OptionParser, OptionGroup

usage ="""
===========================================================================
Usage: %prog [options] arg1 arg2 arg3
Author: Dan Kelley
	Waisman Laboratory for Brain Imaging and Behavior
	University of Wisconsin-Madison
	Last Update: November 1, 2007
This program:
1 Removes noise with a low pass Butterworth filter with zero phase shift
2 Calculates the residual = Raw signal - filtered signal 
3 Idenitifies dirty peaks with big residuals that need cleaning 
4 Dirty peaks  +/- a window are removed to produce filtered, cleaned data 
5 Calculates the best fitting spline through the cleaned data 
6 Spline interpolates over the removed dirty peaks to produce a clean
respiratory waveform (RW) and cardiac waveform (CW)
This is fast on Mac OS. Other platforms use the --speed option.
7 Finds the respiratory peaks (RPpd) and cardiac peaks (CPpd).
Peak finding was slightly modified from the open source peakdet.m algorithm 
to include magnitude and time thresholds and ported into python.
8 The top and bottom envelopes are generated using the RPpd. The 
respiratory volumes over time (RVT) is their difference as in:
 	Birn RM, Diamond JB, Smith MA, Bandettini PA.	 
 	Separating respiratory-variation-related fluctuations 
 	from neuronal-activity-related fluctuations in fMRI. 
 	Neuroimage. 2006 Jul 15;31(4):1536-48. Epub 2006 Apr 24. 
 	PMID: 16632379 
The same is done for the cardiac waveforms using the CPpd to produce CVT curves.
9 Cardiac rate time (CRT) courses based on the TTL (CRTttl)are generated as the 
inverse change in time between CP centered points on the initial peak as in:
   	Shmueli K, van Gelderen P, de Zwart JA, Horovitz SG, Fukunaga M, 
	Jansma JM, Duyn JH.	
  	Low-frequency fluctuations in the cardiac rate as a source of variance in 
	the resting-state fMRI BOLD signal.
	Neuroimage. 2007 Nov 1;38(2):306-20. Epub 2007 Aug 9.
	PMID: 17869543
  An RRT is also generated in the same manner using the RPpd. 
10 Cardiac rate time (CRT) courses based on the third derivative (CRTd3)
 and the third derivative's R-wave estimate (CRTd3R) are generated as in:
 	Chan GS, Middleton PM, Celler BG, Wang L, Lovell NH.	
	Automatic detection of left ventricular ejection time from a finger 
	photoplethysmographic pulse oximetry waveform: comparison with Doppler 
	aortic measurement.
	Physiol Meas. 2007 Apr;28(4):439-52. Epub 2007 Mar 20.
	PMID: 17395998
11 PhysioNoise outputs the RW (downsampled to 40Hz by default), CPd3, CPd3R, and CPttl 
time courses as options for 3dretroicor and the RVT, RRT, CVT, CRTd3, CRTd3R, and 
CRTttl curves sampled on the TR and half TR as covariate options for neuroanalysis. 
   
   usage: %prog [options] arg1 arg2

Example usage:
./PhysioNoise.py -c CardiacWaveform.txt -o CardiacTrigger.txt -r RespiratoryWaveform.txt
./PhysioNoise.py -c CardiacWaveform.txt -o CardiacTrigger.txt -r RespiratoryWaveform.txt --speed
./PhysioNoise.py -c CardiacWaveform.txt -o CardiacTrigger.txt -r RespiratoryWaveform.txt --numTR 360 --rmagthresh 80 --cmagthresh 10 -p PreIcor_Subj1 --speed

"""


epilogue ="""===========================================================================
"""

parser = OptionParser(usage=usage, epilog=epilogue)
parser.add_option("-r", "--resp",  action="store", type="string", dest="sigfile",help="read respiratory signal from FILE", metavar="FILE")
parser.add_option("-c", "--card",  action="store", type="string", dest="csigfile",help="read cardiac signal from FILE", metavar="FILE")
parser.add_option("-o", "--ox",  action="store",type="string", dest="ctrigfile",help="read cardiac pulse ox trigger from FILE", metavar="FILE")
parser.add_option("--TR",  action="store",type="float", dest="TR",help="TR is TR seconds [2.5]", metavar="TR", default=2.5)
parser.add_option("--numTR",  action="store",type="int", dest="NumTR",help="Number of TRs [480]", metavar="TR", default=480)
parser.add_option("--inputHz",  action="store",type="int", dest="OrigSamp",help="input frequency in Hz[1000]", metavar="Hz", default=1000)
parser.add_option("--outputHz",  action="store",type="int", dest="NewSamp",help="desired output frequency in Hz[40]", metavar="Hz", default=40)
parser.add_option("--fpass",  action="store",type="int", dest="fpass",help="Butterworth filter:Stop frequency of passband in Hz[2]", metavar="Hz", default=2)
parser.add_option("--fstop",  action="store",type="int", dest="fstop",help="Butterworth filter:Start frequency of the stopband  in Hz[10]", metavar="Hz", default=10)
parser.add_option("--trimwindow",  action="store",type="int", dest="TrimWindow",help="Number of points PTS to clean on either side of dirty residual peaks[650]", metavar="PTS", default=650)
parser.add_option("--plr",  action="store",type="float", dest="FwdRm",help="Percentage of points on the left of dirty peaks to clean on the right [1.2]", metavar="Percentage", default=1.2)
parser.add_option("--statthresh",  action="store",type="int", dest="FoldSTD",help="Mark points beyond FOLD standard deviation(s) as dirty residual peaks for cleaning  [6]", metavar="FOLD", default=6)
parser.add_option("--rmagthresh",  action="store",type="int", dest="RMagThresh",help="Respiratory Peak must be RMAG units away from nearest extrema  [100]", metavar="RMAG", default=100)
parser.add_option("--rtimethresh",  action="store",type="float", dest="RTimeThresh",help="Respratory Peak must be RTIME seconds away from the nearest extrema  [0]", metavar="RTIME", default=0)
parser.add_option("--cmagthresh",  action="store",type="int", dest="CMagThresh",help="Cardiac Peak must be CMAG units away from nearest extrema  [20]", metavar="CMAG", default=20)
parser.add_option("--ctimethresh",  action="store",type="float", dest="CTimeThresh",help="Cardiac Peak must be CTIME seconds away from the nearest extrema  [0]", metavar="CTIME", default=0)
parser.add_option("--plots",  action="store_true",dest="PlotAll",help="this turns on interactive plotting", default=False)
parser.add_option("--speed",  action="store_true",dest="Speed",help="This skips artifact detection and spline interpolation. Instead a cubic spline filter is used and does a decent job removing artifact. [artifact detect]", default=False)
parser.add_option("-p","--prefix",  action="store", type="string",dest="Prefix",help="Output prefix [PreIcor]", default='PreIcor')



options, args = parser.parse_args()



if '-help' in sys.argv:
	parser.print_help()
	raise SystemExit()
if not (options.sigfile and options.csigfile and options.ctrigfile) or '-help' in sys.argv:
	print "The respiratory signal file, cardiac signal file, and cardiac trigger file are required. Try --help "
	raise SystemExit()
        
if options.sigfile and options.csigfile and options.ctrigfile:
	Sig = fromfile( options.sigfile, 'f', sep='\n')      
	meanSig=float(mean(Sig[0:long(len(Sig)/2)]))
	Sig = Sig - meanSig
	SigT = fromfile( options.ctrigfile, 'f', sep='\n')      
	CSig = fromfile( options.csigfile, 'f', sep='\n')
	meanCSig=float(mean(CSig[0:long(len(CSig)/2)]))
	CSig = CSig - meanCSig
	CSigT = fromfile( options.ctrigfile, 'f', sep='\n')      

TR=options.TR 
OrigSamp=options.OrigSamp
NewSamp=options.NewSamp
SampRatio=long(OrigSamp/NewSamp)
SampRatioIndex=SampRatio-1
samplerate=OrigSamp
Nyquist=long(samplerate/2)
PlotAll=options.PlotAll	
Speed=options.Speed
MagThresh=options.RMagThresh
OrigTimeThresh=long(options.RTimeThresh*OrigSamp)
NewTimeThresh=long(options.RTimeThresh*NewSamp)
NumTR=options.NumTR
TotOutSampNew= long(TR*NumTR*NewSamp)
TotOutSampOrig= long(TR*NumTR*OrigSamp)

COrigSamp=OrigSamp
CNewSamp=NewSamp
CSampRatio=SampRatio
CSampRatioIndex=SampRatioIndex
Csamplerate=OrigSamp
CNyquist=long(Csamplerate/2)
CMagThresh=options.CMagThresh
COrigTimeThresh= long(options.CTimeThresh*OrigSamp)
CNewTimeThresh= long(options.CTimeThresh*NewSamp)
CTotOutSampNew= long(TR*NumTR*CNewSamp)
CTotOutSampOrig= long(TR*NumTR*COrigSamp)


#Trim
#Sig=Sig[0:TotOutSampOrig]
#CSig=CSig[0:TotOutSampOrig]

####--------------------------------Butter Filter  Data
	
fpass=options.fpass #passband edge frequency
fstop=options.fstop  #stopband edge frequency

filter_b,filter_a=scipy.signal.iirdesign(float(fpass)/Nyquist,float(fstop)/Nyquist,gpass=5,gstop=60,analog=0, ftype="butter")

FiltFiltSig=filtfilt(filter_b, filter_a, Sig)

w, h = freqz(filter_b,filter_a)
mag = abs(h)
phase = unwrap(angle(h))*180/pi
radsample=w/pi

x= int((log(OrigSamp)-log(0.01))/log(2))
NFFT=int(2**x)
nooverlap=int(0.5*(2**x))
pxx, freqs = psd( FiltFiltSig, NFFT=NFFT , Fs=OrigSamp,noverlap=nooverlap , detrend=detrend_linear, window=window_hanning)
maxfreqpt = pxx.argmax()
maxfreq =freqs[maxfreqpt]
TPeriod=1./maxfreq
PeakPower=pxx[maxfreqpt]

if maxfreq==0:
	print "Respiratory power peak at zero. Threshold at",freqs[4],"Hz"
	Tempxx=pxx[5:]
	Tempfreqs=freqs[5:]
	maxfreqpt = Tempxx.argmax()
	maxfreq =Tempfreqs[maxfreqpt]
	TPeriod=1./maxfreq
	PeakPower=Tempxx[maxfreqpt]

print "The spectral respiratory peak is:", maxfreq,"Hz "
print "The respiratory peak power density is",PeakPower
print "The spectral respiratory period is", 1./maxfreq,"seconds"


#Cardiac
Cfilter_b,Cfilter_a=scipy.signal.iirdesign(float(fpass)/Nyquist,float(fstop)/Nyquist,gpass=5,gstop=60,analog=0, ftype="butter")
CFiltFiltSig=filtfilt(Cfilter_b, Cfilter_a, CSig)
Cpxx, Cfreqs = psd( CFiltFiltSig, NFFT=NFFT , Fs=OrigSamp,noverlap=nooverlap , detrend=detrend_linear, window=window_hanning)
Cmaxfreqpt = Cpxx.argmax()
Cmaxfreq =Cfreqs[Cmaxfreqpt]
CTPeriod=1./Cmaxfreq
CPeakPower=pxx[maxfreqpt]
if Cmaxfreq==0:
	print "Cardiac power peak at zero. Threshold at",freqs[4],"Hz"
	CTempxx=Cpxx[5:]
	CTempfreqs=Cfreqs[5:]
	Cmaxfreqpt = CTempxx.argmax()
	Cmaxfreq =CTempfreqs[Cmaxfreqpt]
	CTPeriod=1./Cmaxfreq
	CPeakPower=CTempxx[Cmaxfreqpt]

print "The spectral cardiac peak is:", Cmaxfreq,"Hz "
print "The cardiac peak power density is",CPeakPower
print "The spectral cardiac period is", 1./Cmaxfreq,"seconds"

#FinalFreq=[]
#if Cmaxfreq==0:
#	CidPeaks= nonzero(greater_equal( Cpxx,100000))
#	Cmaxfreq=Cfreqs[CidPeaks]
#	print "The spectral cardiac power densities are:",Cpxx[CidPeaks]
#	print "The spectral cardiac frequency is:", Cmaxfreq,"Hz"
#	print "The spectral cardiac period is", 1./Cmaxfreq,"seconds"



if PlotAll:
	figure(1)
	subplot(311)
	semilogy(radsample*Nyquist,mag)
	ylabel('Magnitude')
	title('Bode Plot')
	pylab.ylim((10e-4, 10e0))
	pylab.xlim((0.0, fstop))
	
	figure(1)
	subplot(312) 
	plot(freqs,pxx)
	pylab.ylabel('Respiratory Power')
	pylab.xlim((0.0, fstop))

	figure(1)
	subplot(313) 
	plot(Cfreqs,Cpxx)
	pylab.xlabel('frequency(Hz)')
	pylab.ylabel('Cardiac Power')
	pylab.xlim((0.0, fstop))

#--------------------------------------------------------------------# Expand Trigger File      
FiltSig=FiltFiltSig 
SigTime_ms = range(len(Sig))
SigTPeak=zeros(len(Sig),'int')
SigTtips=zeros(len(SigT),'f')
SigTtipsTime=zeros(len(SigT),'int')

CFiltSig=CFiltFiltSig 
CSigTime_ms = range(len(CSig))
CSigTPeak=zeros(len(CSig),'int')
CSigTtips=zeros(len(CSigT),'f')
CSigTtipsTime=zeros(len(CSigT),'int')



for counter in xrange(len(CSigT)) :
	ms = CSigT[counter]
	CSigTPeak[ms] = CFiltSig[ms]
	CSigTtips[counter]=CFiltSig[ms]
	CSigTtipsTime[counter]= ms
CPttl=CSigTPeak
#---------------------------------------------------#Remove regions of high noise variance for interpolation

ResidSig=Sig - FiltSig
CResidSig=CSig - CFiltSig

if Speed:
	print "Spline Filtering"
	print "Respiratory"
	BestFit = spline_filter1d(FiltSig, order=3)
	print "Cardiac"
	CBestFit= spline_filter1d(CFiltSig, order=3)
else:
	print "Finding the regions of high variance"

	print "Start Resp"
	FoldSTD=options.FoldSTD
	OutlierThresh =float(FoldSTD)*std(ResidSig)
	RemoveRange= options.TrimWindow
	ExtremTime=0
	FwdRm=options.FwdRm
	AllOutliers=zeros(len(FiltSig))
	DirtyPoints= []
	SigThreshResid = abs(ResidSig) - OutlierThresh
	idUglyPoints= nonzero(greater_equal(SigThreshResid,0))
	for counter in idUglyPoints:
		AllOutliers[counter]=ResidSig[counter]
		for i in (range(long( counter -RemoveRange),long(counter+FwdRm*RemoveRange))) :
			DirtyPoints.extend([i])
	DirtyPoints= sort(unique(clip(unique(DirtyPoints),1,len(ResidSig))))
	NoOutliers=zeros(len(FiltSig))
	NoOutlierTime=zeros(len(SigTime_ms))
	NoOutliers=delete(FiltSig,DirtyPoints)
	NoOutlierTime=delete(SigTime_ms,DirtyPoints)

	#Cardiac
	print "Start Cardiac"
	CFoldSTD=options.FoldSTD
	COutlierThresh =float(FoldSTD)*std(CResidSig)
	CRemoveRange= options.TrimWindow
	CExtremTime=0
	CFwdRm=options.FwdRm
	CAllOutliers=zeros(len(CFiltSig))
	CDirtyPoints= []
	CSigThreshResid = abs(CResidSig) - COutlierThresh
	CidUglyPoints= nonzero(greater_equal(CSigThreshResid,0))
	for counter in CidUglyPoints:
		CAllOutliers[counter]=CResidSig[counter]
		for i in (range(long( counter -CRemoveRange),long(counter+CFwdRm*CRemoveRange))) :
			CDirtyPoints.extend([i])
	CDirtyPoints= sort(unique(clip(unique(CDirtyPoints),1,len(CResidSig))))
	CNoOutliers=zeros(len(CFiltSig))
	CNoOutlierTime=zeros(len(CSigTime_ms))
	CNoOutliers=delete(CFiltSig,CDirtyPoints)
	CNoOutlierTime=delete(CSigTime_ms,CDirtyPoints)





	print "Spline interpolation over deleted acquisition artifact"
	print "Respiratory"
	#Generate cubic spline....Too slow when not on Mac OS X but looks nice
	SplineEst = scipy.interpolate.splrep(NoOutlierTime ,NoOutliers, k=3)
	BestFit = scipy.interpolate.splev(SigTime_ms  , SplineEst)
	print "Cardiac"

	CSplineEst = scipy.interpolate.splrep(CNoOutlierTime ,CNoOutliers, k=3)
	CBestFit = scipy.interpolate.splev(CSigTime_ms  , CSplineEst)

#SplineEst = cspline1d(FiltSig)
#BestFit = scipy.signal.cspline1d_eval(SplineEst,SigTime_ms)
#print "Cardiac"
#CSplineEst = cspline1d(CFiltSig)
#CBestFit = scipy.signal.cspline1d_eval(CSplineEst,CSigTime_ms)




if PlotAll:
#	figure(2000)
#	plot(NoOutlierTime,NoOutliers,'bo',label='NoOutliers')	
#	pylab.legend()

	figure(2)
	subplot(211)
	plot(SigTime_ms,FiltSig,'b-', SigTime_ms,BestFit,'go', SigTime_ms, Sig, 'k', SigTime_ms, ResidSig,'r'  )	
	pylab.legend((r'Filtered', r'Spline', r'Raw', r'Residual'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'b')
	pylab.setp(ltext[1],  fontsize = 8, color = 'g')
	pylab.setp(ltext[2],  fontsize = 8, color = 'k')
	pylab.setp(ltext[3],  fontsize = 8, color = 'r')
	title('Respiratory')

	subplot(212)
	plot(CSigTime_ms,CFiltSig,'b-', CSigTime_ms,CBestFit,'go', CSigTime_ms, CSig, 'k', CSigTime_ms, CResidSig,'r'  )	
	pylab.legend((r'Filtered', r'Spline', r'Raw', r'Residual'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'b')
	pylab.setp(ltext[1],  fontsize = 8, color = 'g')
	pylab.setp(ltext[2],  fontsize = 8, color = 'k')
	pylab.setp(ltext[3],  fontsize = 8, color = 'r')
	title('Cardiac')


#---------------------------------------------Data is now clean
CleanSig=BestFit
CCleanSig=CBestFit	

AllExtrema, AllExtremaTime,AllHills, AllValleys, TopEnv, BotEnv, TopCounts, BotCounts=findextrema(CleanSig,MagThresh, OrigTimeThresh)
CAllExtrema, CAllExtremaTime,CAllHills, CAllValleys, CTopEnv, CBotEnv, CTopCounts, CBotCounts=findextrema(CCleanSig,CMagThresh, COrigTimeThresh)

"""
if PlotAll:
#Resp
	figure(3)
	subplot(321);plot(SigTime_ms,Sig,'b',SigTime_ms,FiltSig,'g',SigTime_ms, ResidSig,'r'  )
	ylabel('Amplitude')
	grid(True)
	pylab.legend((r'Raw', r'Filter', r'Residual'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'b')
	pylab.setp(ltext[1],  fontsize = 8, color = 'g')
	pylab.setp(ltext[2],  fontsize = 8, color = 'r')
	title('Respiratory')


	figure(3)
	subplot(323)
	plot(SigTime_ms, FiltSig, 'b', SigTime_ms, SigTPeak,'r',linewidth=1.0)
	pylab.legend((r'Filter', r'TTLTrigger'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'b')
	pylab.setp(ltext[1],  fontsize = 8, color = 'r')

	

	figure(3)
	subplot(325); plot(SigTime_ms,CleanSig,'k', SigTime_ms, FiltSig, 'b',SigTime_ms,TopEnv,'r', SigTime_ms, BotEnv, 'g');
	pylab.legend((r'Spline', r'Filtered', r'TopEnv', r'BotEnv'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'k')
	pylab.setp(ltext[1],  fontsize = 8, color = 'b')
	pylab.setp(ltext[2],  fontsize = 8, color = 'r')
	pylab.setp(ltext[3],  fontsize = 8, color = 'g')
	xlabel('Time (ms)')

#Cardiac	
	
	figure(3)
	subplot(322);plot(CSigTime_ms,CSig,'b',CSigTime_ms,CFiltSig,'g',CSigTime_ms, CResidSig,'r'  )
	ylabel('Amplitude')
	grid(True)
	pylab.legend((r'Raw', r'Filter', r'Residual'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'b')
	pylab.setp(ltext[1],  fontsize = 8, color = 'g')
	pylab.setp(ltext[2],  fontsize = 8, color = 'r')
	title('Cardiac')


	figure(3)
	subplot(324)
	plot(CSigTime_ms, CFiltSig, 'b', CSigTime_ms, CSigTPeak,'r',linewidth=1.0)
	pylab.legend((r'Filter', r'TTLTrigger'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'b')
	pylab.setp(ltext[1],  fontsize = 8, color = 'r')


	figure(3)
	subplot(326); plot(CSigTime_ms,CCleanSig,'k', CSigTime_ms, CFiltSig, 'b',CSigTime_ms,CTopEnv,'r', CSigTime_ms, CBotEnv, 'g');
	pylab.legend((r'Spline', r'Filtered', r'TopEnv', r'BotEnv'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'k')
	pylab.setp(ltext[1],  fontsize = 8, color = 'b')
	pylab.setp(ltext[2],  fontsize = 8, color = 'r')
	pylab.setp(ltext[3],  fontsize = 8, color = 'g')
	xlabel('Time (ms)')

#Resp	
	figure(4)
	subplot(211);
	plot(SigTime_ms,CleanSig,'k',SigTime_ms, FiltSig,'b',SigTime_ms, TopEnv,'c', SigTime_ms, BotEnv, 'm', SigTime_ms, AllOutliers,'go', SigTime_ms, Sig,'y', NoOutlierTime, NoOutliers,'ro',SigTime_ms, ResidSig,'r');
	pylab.legend((r'Spline', r'Filtered', r'TopEnv', r'BotEnv', r'Outliers', r'RawSignal', r'NoOutliers',r'Resid'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'k')
	pylab.setp(ltext[1],  fontsize = 8, color = 'b')
	pylab.setp(ltext[2],  fontsize = 8, color = 'c')
	pylab.setp(ltext[3],  fontsize = 8, color = 'm')
	pylab.setp(ltext[4],  fontsize = 8, color = 'g')
	pylab.setp(ltext[5],  fontsize = 8, color = 'y')
	pylab.setp(ltext[6],  fontsize = 8, color = 'r')
	pylab.setp(ltext[7],  fontsize = 8, color = 'r')
	title('Respiratory')

#Card
	figure(4)
	subplot(212);
	plot(CSigTime_ms,CCleanSig,'k',CSigTime_ms, CFiltSig,'b',CSigTime_ms, CTopEnv,'c', CSigTime_ms, CBotEnv, 'm', CSigTime_ms, CAllOutliers,'go', CSigTime_ms, CSig,'y', CNoOutlierTime, CNoOutliers,'ro',CSigTime_ms, CResidSig,'r');
	pylab.legend((r'Spline', r'Filtered', r'TopEnv', r'BotEnv', r'Outliers', r'RawSignal', r'NoOutliers',r'Resid'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'k')
	pylab.setp(ltext[1],  fontsize = 8, color = 'b')
	pylab.setp(ltext[2],  fontsize = 8, color = 'c')
	pylab.setp(ltext[3],  fontsize = 8, color = 'm')
	pylab.setp(ltext[4],  fontsize = 8, color = 'g')
	pylab.setp(ltext[5],  fontsize = 8, color = 'y')
	pylab.setp(ltext[6],  fontsize = 8, color = 'r')
	pylab.setp(ltext[7],  fontsize = 8, color = 'r')
	title('Cardiac')

"""

#-----------------------------------------------------------------------------Find Volume over Time
print"Calculating Envelopes"

print "Respiratory"

RWTopEnv=TopEnv
FullInterpTimeEnv=ConnectTheDots(TopCounts,diff(TopCounts) ,CleanSig,1)
FullInterpTopEnv=ConnectTheDots(TopCounts,AllHills,CleanSig,0)
FullInterpBotEnv=ConnectTheDots(BotCounts,AllValleys,CleanSig,0)
FullInterpRVTEnv=zeros(len(CleanSig), 'f')
FullInterpRVTEnv = OrigSamp*abs(FullInterpTopEnv-FullInterpBotEnv)/ (FullInterpTimeEnv)

MeanRVT=zeros(len(CleanSig), 'f')
MeanRVT = (    abs(FullInterpTopEnv)+   abs(FullInterpBotEnv))/ 2.

print "Cardiac"
CRWTopEnv=CTopEnv
CFullInterpTimeEnv=ConnectTheDots(CTopCounts,diff(CTopCounts) ,CCleanSig,1)
CFullInterpTopEnv=ConnectTheDots(CTopCounts,CAllHills,CCleanSig,0)
CFullInterpBotEnv=ConnectTheDots(CBotCounts,CAllValleys,CCleanSig,0)
CFullInterpRVTEnv=zeros(len(CCleanSig), 'f')
CFullInterpRVTEnv = COrigSamp*abs(CFullInterpTopEnv-CFullInterpBotEnv)/ (CFullInterpTimeEnv)

CMeanRVT=zeros(len(CCleanSig), 'f')
CMeanRVT = (abs(CFullInterpTopEnv)+abs(CFullInterpBotEnv))/ 2.



if PlotAll:
	figure(7)
	subplot(211)
	plot(SigTime_ms,CleanSig,'b'); hold(True);
	plot(SigTime_ms,FullInterpTopEnv, 'r', SigTime_ms , FullInterpBotEnv , 'g',SigTime_ms ,TopEnv, 'ro',SigTime_ms ,BotEnv, 'go',SigTime_ms , FullInterpRVTEnv  ,'c', SigTime_ms, MeanRVT, 'y')
	hold(False)
	pylab.legend((r'Spline', r'TopEnvInterpl', r'BotEnvInterpl', r'TopEnv', r'BotEnv', r'RVTinterpl',  r'MeanAbsEnv' ), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'b')
	pylab.setp(ltext[1],  fontsize = 8, color = 'r')
	pylab.setp(ltext[2],  fontsize = 8, color = 'g')
	pylab.setp(ltext[3],  fontsize = 8, color = 'r')
	pylab.setp(ltext[4],  fontsize = 8, color = 'g')
	pylab.setp(ltext[5],  fontsize = 8, color = 'c')
	pylab.setp(ltext[6],  fontsize = 8, color = 'y')
	title('Respiratory Data')

	figure(7)
	subplot(212)
	plot(CSigTime_ms,CCleanSig,'b'); hold(True);
	plot(CSigTime_ms,CFullInterpTopEnv, 'r', CSigTime_ms , CFullInterpBotEnv , 'g',CSigTime_ms ,CTopEnv, 'ro',CSigTime_ms ,CBotEnv, 'go',CSigTime_ms , CFullInterpRVTEnv  ,'c',  CSigTime_ms , CMeanRVT, 'y')
	hold(False)
	pylab.legend((r'Spline', r'TopEnvInterpl', r'BotEnvInterpl', r'TopEnv', r'BotEnv', r'CVTinterpl', r'MeanAbsEnv' ), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'b')
	pylab.setp(ltext[1],  fontsize = 8, color = 'r')
	pylab.setp(ltext[2],  fontsize = 8, color = 'g')
	pylab.setp(ltext[3],  fontsize = 8, color = 'r')
	pylab.setp(ltext[4],  fontsize = 8, color = 'g')
	pylab.setp(ltext[5],  fontsize = 8, color = 'c')
	pylab.setp(ltext[6],  fontsize = 8, color = 'y')
	title('Cardiac Data')
	

#-----------------------------Find third derivatives of Cardiac signal to create new trigger file
#CCSigBase=CCleanSig

CCSigBase= CCleanSig
print"Finding Cardiac Derivatives"
dCCleanSig1=slopes(CSigTime_ms,CCSigBase)    *(COrigSamp)
#figure(3333321)
#plot(CSigTime_ms,dCCleanSig1,'b',CSigTime_ms,CCSigBase,'r',CSigTime_ms,CCleanSig,'c')

dCCleanSig2=slopes(CSigTime_ms,dCCleanSig1)   *(COrigSamp)
dCCleanSig3=slopes(CSigTime_ms,dCCleanSig2)   *(COrigSamp)
CSmoothCCleanSig3=cmov_average(dCCleanSig3,int(0.05*COrigSamp))
#Find large Troughs
CAllExtrema, CAllExtremaTime,CAllHills, CAllValleys, CTopEnv, CBotEnv, CTopCounts, CBotCounts=findextrema( CSmoothCCleanSig3,MagThresh, NewTimeThresh)


idRise=nonzero(less_equal(dCCleanSig1 ,0 ))
put(CBotEnv,idRise,0)
CPd3=-1.*CBotEnv
#CPd3 has location of the large deriv3 troughs

#use envelope of deriv3 peaks as a threshold
CPd3Threshd3 =ConnectTheDots(CTopCounts,abs(CAllHills),CSig,0)

print "Thresholding CPd3"

idRise=nonzero(less_equal(CPd3-CPd3Threshd3,0))
put(CPd3,idRise,0)
idRise=nonzero(less_equal(CPd3,0.5*mean(CPd3Threshd3)))
put(CPd3,idRise,0)

# now find R waves which are the deriv3 peaks that precede the trough which is in CPd3
CPd3Rwave=zeros(len(CCSigBase), 'f')
CRwaveCounts=array(nonzero(greater(CPd3,0)))
for peaktime in CRwaveCounts:
	tempR=nonzero(less(CTopCounts,peaktime))
	CPd3Rwave[CTopCounts[tempR[-1]]]= CPd3[peaktime]
	
if PlotAll:

	figure(41)
	subplot(211)
	plot(CSigTime_ms,CCSigBase*(COrigSamp),'k', CSigTime_ms, CSigTPeak*(COrigSamp), 'r',CSigTime_ms,CPd3,'b', CSigTime_ms, dCCleanSig1*(COrigSamp), 'g', CSigTime_ms, dCCleanSig2, 'y',CSigTime_ms, CSmoothCCleanSig3, 'c');
	pylab.legend((r'Filter', r'CPttl', r'CPd3', r'Deriv1', r'Deriv2', r'Deriv3'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'k')
	pylab.setp(ltext[1],  fontsize = 8, color = 'r')
	pylab.setp(ltext[2],  fontsize = 8, color = 'b')
	pylab.setp(ltext[3],  fontsize = 8, color = 'g')
	pylab.setp(ltext[4],  fontsize = 8, color = 'y')
	pylab.setp(ltext[5],  fontsize = 8, color = 'c')
	xlabel('Time (ms)')
	subplot(212)
	plot(CSigTime_ms,CCSigBase*(COrigSamp),'k', CSigTime_ms, CSigTPeak*(COrigSamp), 'r',CSigTime_ms, CPd3 ,'b',CSigTime_ms, CSmoothCCleanSig3, 'c',CSigTime_ms, CPd3Threshd3, 'y',CSigTime_ms, CPd3Rwave, 'm');
	pylab.legend((r'Filter', r'CPttl', r'CPd3', r'Deriv3', r'Threshold', r'Rwave'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'k')
	pylab.setp(ltext[1],  fontsize = 8, color = 'r')
	pylab.setp(ltext[2],  fontsize = 8, color = 'b')
	pylab.setp(ltext[3],  fontsize = 8, color = 'c')
	pylab.setp(ltext[4],  fontsize = 8, color = 'y')
	pylab.setp(ltext[5],  fontsize = 8, color = 'm')
	xlabel('Time (ms)')

#------------------Calculate respiratory (RRT) and cardiac rate per time (CRT)
print "Calculating Rates over Time"
#Calculate RRT

RCPd3SigT=array(nonzero(greater(TopEnv,0)))
RCSigTHills=COrigSamp * (diff(RCPd3SigT)**(-1.))
RCSigTTopCounts=RCPd3SigT[0:-1]
RRT=ConnectTheDots(RCSigTTopCounts,RCSigTHills,Sig,0)


CSigTHills= (diff(CSigT)**(-1))*COrigSamp
CSigTTopCounts=CSigT[0:-1]
CRTttl=ConnectTheDots(CSigTTopCounts,CSigTHills,CSig,0)


#Calculate CRTd3

CPd3SigT=array(nonzero(greater(CPd3,0)))
CSigTHills=COrigSamp * (diff(CPd3SigT)**(-1.))
CSigTTopCounts=CPd3SigT[0:-1]
CRTd3=ConnectTheDots(CSigTTopCounts,CSigTHills,CSig,0)

#Calculate CRTd3R

CPd3RSigT=array(nonzero(greater(CPd3Rwave,0)))
CSigTHills=COrigSamp * (diff(CPd3RSigT)**(-1.))
CSigTTopCounts=CPd3RSigT[0:-1]
CRTd3R=ConnectTheDots(CSigTTopCounts,CSigTHills,CSig,0)

if PlotAll:
	figure(333)
	subplot(321)
	plot(CRTd3R)
	ylabel('CRTd3R')
	subplot(322)
	plot(CRTttl)
	ylabel('CRTttl')
	subplot(324)
	plot(CRTd3)
	ylabel('CRTd3')
	subplot(323)
	plot(RRT)
	ylabel('RRT')
	subplot(325)
	plot(FullInterpRVTEnv)
	ylabel('RVT')
	subplot(326)
	plot(CFullInterpRVTEnv)
	ylabel('CVT')

#-----------------------------------------------------------------------------DownSampling
print "Downsampling"
#Resp
SigDownSamp,SigDownSamplePt, SigDownSamplePtTime=DownSample(CleanSig)
#Cardiac
CSigDownSamp,CSigDownSamplePt, CSigDownSamplePtTime=DownSample(CCleanSig)
"""
if PlotAll:
	figure(5)
	subplot(211)
	plot(SigTime_ms, CleanSig,'b',SigTime_ms, Sig ,'k', SigTime_ms, SigDownSamp,'r')
	pylab.legend((r'Spline', r'Raw', r'Downsampled'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'b')
	pylab.setp(ltext[1],  fontsize = 8, color = 'k')
	pylab.setp(ltext[2],  fontsize = 8, color = 'r')
	title('Downsampled Respiratory Spline')

	figure(5)
	subplot(212)
	plot(CSigTime_ms, CCleanSig,'b',CSigTime_ms, CSig ,'k', CSigTime_ms, CSigDownSamp,'r')
	pylab.legend((r'Spline', r'Raw', r'Downsampled'), shadow = True, loc = 1)
	ltext = pylab.gca().get_legend().get_texts()
	pylab.setp(ltext[0],  fontsize = 8, color = 'b')
	pylab.setp(ltext[1],  fontsize = 8, color = 'k')
	pylab.setp(ltext[2],  fontsize = 8, color = 'r')
	title('Downsampled Cardiac Spline')
"""
#-----------------------------------------------Trim length for retroicor
print "Trim length"
SigDownSamplePt=SigDownSamplePt[0:TotOutSampNew]
RWTopEnv=TopEnv[0:TotOutSampOrig]
FullInterpRVTEnv= FullInterpRVTEnv[0:TotOutSampOrig]
RRT=RRT[0:TotOutSampOrig]

CleanSig=CleanSig[0:TotOutSampOrig]
CCleanSig=CCleanSig[0:CTotOutSampOrig]

CPd3=CPd3[0:TotOutSampOrig]
CPd3Rwave=CPd3Rwave[0:TotOutSampOrig]
CPttl=CPttl[0:TotOutSampOrig]

CRTd3=CRTd3[0:TotOutSampOrig]
CRTttl=CRTttl[0:TotOutSampOrig]
CRTd3R=CRTd3R[0:TotOutSampOrig]

CSigDownSamplePt=CSigDownSamplePt[0:CTotOutSampNew]
CRWTopEnv=CRWTopEnv[0:CTotOutSampOrig]
CFullInterpRVTEnv= CFullInterpRVTEnv[0:CTotOutSampOrig]

#---------------Calculating Stats
print "Calculating Statistics for", options.Prefix
OutputStat= zeros(65,'f')

OutputStat[0]=( len(nonzero(greater(RWTopEnv ,0 ))) )
print  "Peakdet detected ",OutputStat[0]  ," positive respiratory peaks"  		  
OutputStat[1]=( mean(diff(nonzero(greater( RWTopEnv ,0 ))))  )
print  "The RW peakdet average respiratory interpeak time was",  OutputStat[1]
OutputStat[2]=( std(diff(nonzero(greater( RWTopEnv ,0 ))))   )
print  "The RW peakdet std respiratory interpeak time was",  OutputStat[2]

OutputStat[3]=( maxfreq )
print  "The spectral respiratory peak is:",OutputStat[3]  ,"Hz "					  
OutputStat[4]=(PeakPower  )
print  "The respiratory peak power density is", OutputStat[4]						  	  
OutputStat[5]=(1./maxfreq  )
print  "The spectral respiratory period is",OutputStat[5]  ,"seconds"	  				  

OutputStat[6]=( Cmaxfreq )
print  "The spectral cardiac peak is:",OutputStat[6]  ,"Hz " 						  	  
OutputStat[7]=(CPeakPower  )
print  "The cardiac peak power density is", OutputStat[7] 					  
OutputStat[8]=( 1./Cmaxfreq )
print  "The spectral cardiac period is",OutputStat[8]  ,"seconds" 					  

OutputStat[9]=( len(nonzero(greater( CPttl ,0 ))) )
print  "The trigger detected ",OutputStat[9] ," positive cardiac peaks"    	  
OutputStat[10]=( len(nonzero(greater(CPd3Rwave ,0 ))) )
print  "CPd3R detected ",OutputStat[10]  ," positive cardiac peaks"  	  	  


OutputStat[11]=( mean(diff(nonzero(greater( CPttl ,0 )))) 	 )
print  "The average cardiac CPttl interpeak time was", 	 OutputStat[11] 
OutputStat[12]=( mean(diff(nonzero(greater( CPd3Rwave ,0 ))))  )
print  "The average cardiac CPd3R interpeak time was",  OutputStat[12]  
OutputStat[13]=( std(diff(nonzero(greater( CPttl ,0 )))) )
print  "The std  cardiac CPttl interpeak time was", OutputStat[13]	  	  
OutputStat[14]=( std(diff(nonzero(greater( CPd3Rwave ,0 ))))  )
print  "The std  cardiac CPd3R interpeak time was", OutputStat[14]	  

OutputStat[15]=( std(CRTttl) )
print  "The CRTttl std is:", 	OutputStat[15]  
OutputStat[16]=( std(CRTd3) )
print  "The CRTd3 std is:",  OutputStat[16]  
OutputStat[17]=( std(CRTd3R) )
print  "The CCRTd3R std is:",  OutputStat[17]	  

OutputStat[18]=( mean(CRTttl)  )
print    "The CRTttl mean is:",  OutputStat[18]  
OutputStat[19]=( mean(CRTd3) )
print    "The CRTd3 mean is:", 	OutputStat[19]  
OutputStat[20]=( mean(CRTd3R)  )
print    "The CCRTd3R mean is:",   OutputStat[20]

OutputStat[21]=( len(nonzero(greater(CPd3 ,0 ))) )
print  "CPd3 detected ", OutputStat[21] ," positive cardiac peaks"  	  	  
OutputStat[22]=( mean(diff(nonzero(greater( CPd3 ,0 ))))  )
print  "The average cardiac CPd3 interpeak time was", OutputStat[22]   
OutputStat[23]=( std(diff(nonzero(greater( CPd3 ,0 ))))  )
print  "The std  cardiac CPd3 interpeak time was", OutputStat[23]	  


OutputStat[24]=( mean((nonzero(greater( CPttl ,0 ))))  )
print  "The avg  cardiac CPttl peak time was", OutputStat[24]	  
OutputStat[25]=( mean((nonzero(greater( CPd3 ,0 ))))  )
print  "The avg  cardiac CPd3 peak time was", OutputStat[25]	  
OutputStat[26]=( mean((nonzero(greater( CPd3Rwave ,0 ))))  )
print  "The avg  cardiac CPd3R peak time was", OutputStat[26] 	  


OutputStat[27]=( std((nonzero(greater( CPttl ,0 ))))  )
print  "The std  cardiac CPttl peak time was", OutputStat[27] 	  
OutputStat[28]=( std((nonzero(greater( CPd3 ,0 ))))  )
print  "The std  cardiac CPd3 peak time was", OutputStat[28] 	  
OutputStat[29]=( std((nonzero(greater( CPd3Rwave ,0 ))))  )
print  "The std  cardiac CPd3R peak time was", OutputStat[29] 	  

# now find find time between CPttl and both CPd3 and CPd3R 
PostPeakd3=array(nonzero(greater(CPd3,0)))
PostPeakttl=array(nonzero(greater(CPttl,0)))
PeakBase=array(nonzero(greater(CPd3Rwave,0)))
InterPeakCountsd3=zeros(len(PeakBase), 'f')
InterPeakCountsttl=zeros(len(PeakBase), 'f')
IPd3Mttl=zeros(len(PeakBase), 'f')
counter=0
for peaktime in PeakBase[1:-1]:
	TempPostd3=nonzero(greater_equal(PostPeakd3,peaktime))
	TempPostttl=nonzero(greater_equal(PostPeakttl,peaktime))
	if len(TempPostd3)!=0  and len(TempPostttl)!=0  :
		InterPeakCountsd3[counter]=PostPeakd3[TempPostd3[0]]- peaktime
		InterPeakCountsttl[counter]=PostPeakttl[TempPostttl[0]]- peaktime
		IPd3Mttl[counter]=InterPeakCountsd3[counter]-InterPeakCountsttl[counter]
	counter=counter+1
OutputStat[30]= mean( InterPeakCountsd3 )
print  "The avg  time between CPd3R and CPd3 peaks was",OutputStat[30] 	  
OutputStat[31]= std( InterPeakCountsd3 )
print  "The std  time between CPd3R and CPd3 peaks was",OutputStat[31] 	  
OutputStat[32]= mean( InterPeakCountsttl )
print  "The avg  time between CPd3R and CPttl peaks was",OutputStat[32] 	  
OutputStat[33]= std( InterPeakCountsttl )
print  "The std  time between CPd3R and CPttl peaks was",OutputStat[33] 	  

OutputStat[34]= mean( IPd3Mttl )
print  "The avg  time between CPd3 and CPttl peaks was",OutputStat[34] 	  
OutputStat[35]= std( IPd3Mttl )
print  "The std  time between CPd3 and CPttl peaks was",OutputStat[35] 	  
#PostPeakBase=array(nonzero(greater(CPttl,0)))
#PeakBase=array(nonzero(greater(CPd3Rwave,0)))
#InterPeakCounts=zeros(len(PeakBase), 'f')
#counter=0
#for peaktime in PeakBase[1:-1]:
#	TempPost=nonzero(greater_equal(PrePeakBase,peaktime))
#	InterPeakCounts[counter]=PostPeakBase[TempPost[0]]- peaktime
#	counter=counter+1

OutputStat[36]= len(nonzero(greater( CRTttl ,mean(CRTttl)+3*std(CRTttl) )))
print  "The number of CRTttl peaks > 3 std",OutputStat[36] 	  

OutputStat[37]= len(nonzero(greater( CRTd3 ,mean(CRTd3)+3*std(CRTd3) )))
print  "The number of CRTd3 peaks > 3 std",OutputStat[37] 	  

OutputStat[38]= len(nonzero(greater( CRTd3R ,mean(CRTd3R)+3*std(CRTd3R) )))
print  "The number of CRTd3R peaks > 3 std",OutputStat[38] 	  

OutputStat[39]= len(Sig)
print  "Length of Sig",OutputStat[39] 	  

OutputStat[40]= len(CSig)
print  "Length of CSig",OutputStat[40] 	  

OutputStat[41]= std(CRTttl[nonzero(greater( CRTttl ,mean(CRTttl)+3*std(CRTttl) ))])
print  "The std of CRTttl peaks > 3 std",OutputStat[41] 	  

OutputStat[42]= std(CRTd3[nonzero(greater( CRTd3 ,mean(CRTd3)+3*std(CRTd3) ))])
print  "The std of CRTd3 peaks > 3 std",OutputStat[42] 	  

OutputStat[43]= std(CRTd3R[nonzero(greater( CRTd3R ,mean(CRTd3R)+3*std(CRTd3R) ))])
print  "The std of CRTd3R peaks > 3 std",OutputStat[43] 	  



OutputStat[44]= len(nonzero(greater( CRTttl ,mean(CRTttl)+1. )))
print  "The number of CRTttl peaks > mean CRTttl + 1 ",OutputStat[44] 	  

OutputStat[45]= len(nonzero(greater( CRTd3 ,mean(CRTd3)+1. )))
print  "The number of CRTd3 peaks > mean CRTd3 + 1 ",OutputStat[45] 	  

OutputStat[46]= len(nonzero(greater( CRTd3R ,mean(CRTd3R)+1.)))
print  "The number of CRTd3R peaks > mean CRTd3R + 1 ",OutputStat[46] 	  


OutputStat[47]= len(nonzero(greater( CRTttl ,mean(CRTttl)-0.4 )))
print  "The number of CRTttl peaks 0.4 Hz below mean ",OutputStat[47] 	  

OutputStat[48]= len(nonzero(greater( CRTd3 ,mean(CRTd3)-0.4)))
print  "The number of CRTd3 peaks 0.4 Hz below mean ",OutputStat[48] 	  

OutputStat[49]= len(nonzero(greater( CRTd3R ,mean(CRTd3R)-0.4)))
print  "The number of CRTd3R peaks 0.4 Hz below mean ",OutputStat[49] 	  

OutputStat[50]=( len(nonzero(greater(CRWTopEnv ,0 ))) )
print  "Peakdet detected ",OutputStat[50]  ," positive cardiac peaks"  		  
OutputStat[51]=( mean(diff(nonzero(greater( CRWTopEnv ,0 ))))  )
print  "The  peakdet average cardiac interpeak time was",  OutputStat[51]
OutputStat[52]=( std(diff(nonzero(greater( CRWTopEnv ,0 ))))   )
print  "The peakdet std cardiac interpeak time was",  OutputStat[52]

OutputStat[53]= len(nonzero(greater( CRTttl ,mean(CRTttl)+4*std(CRTttl) )))
print  "The number of CRTttl peaks > 4 std",OutputStat[53] 	  

OutputStat[54]= len(nonzero(greater( CRTd3 ,mean(CRTd3)+4*std(CRTd3) )))
print  "The number of CRTd3 peaks > 4 std",OutputStat[54] 	  

OutputStat[55]= len(nonzero(greater( CRTd3R ,mean(CRTd3R)+4*std(CRTd3R) )))
print  "The number of CRTd3R peaks > 4 std",OutputStat[55] 	  

OutputStat[56]= std(CRTttl[nonzero(greater( CRTttl ,mean(CRTttl)+4*std(CRTttl) ))])
print  "The std of CRTttl peaks > 4 std",OutputStat[56] 	  

OutputStat[57]= std(CRTd3[nonzero(greater( CRTd3 ,mean(CRTd3)+4*std(CRTd3) ))])
print  "The std of CRTd3 peaks > 4 std",OutputStat[57] 	  

OutputStat[58]= std(CRTd3R[nonzero(greater( CRTd3R ,mean(CRTd3R)+4*std(CRTd3R) ))])
print  "The std of CRTd3R peaks > 4 std",OutputStat[58] 	  
print OutputStat
"-----------------"	  

print "Saving files"

#-----------------------------------------------save files to disk

#from numpy import *
#from scipy.io import write_array
#write_array(filename , SigDownSamplePt , separator=' ', linesep='\n')
from pylab import save

filename = options.Prefix +'_OutputStatistics'+".txt"
save(filename, OutputStat, fmt='%8.6f')  

#RW
filename = options.Prefix +'_RW_Spline_'+ str(NewSamp) + ".txt"
save(filename, SigDownSamplePt + abs( min(SigDownSamplePt)), fmt='%8.6f')  


SigTRfull,SigTRPt, SigTRPtTime=TRSample(CleanSig,1.0,OrigSamp)
filename = options.Prefix +'_TR_RW_Spline_'+ str(OrigSamp) + ".txt"
save(filename,  SigTRPt  + abs(min(SigTRPt)), fmt='%8.6f') 

HSigTRfull,HSigTRPt, HSigTRPtTime=TRSample(CleanSig,0.5,OrigSamp)
filename = options.Prefix +'_HalfTR_RW_Spline_'+ str(OrigSamp) + ".txt"
save(filename,  HSigTRPt  + abs(min(HSigTRPt)), fmt='%8.6f'  )

#RVT
filename = options.Prefix +'_RVT_'+ str(OrigSamp) + ".txt"
save(filename,  FullInterpRVTEnv  + abs(min(FullInterpRVTEnv)), fmt='%8.6f' ) 

SigTRfull,SigTRPt, SigTRPtTime=TRSample( FullInterpRVTEnv,1.0,OrigSamp)
filename = options.Prefix +'_TR_RVT_'+ str(OrigSamp) + ".txt"
save(filename,  SigTRPt  + abs(min(SigTRPt)), fmt='%8.6f') 

HSigTRfull,HSigTRPt, HSigTRPtTime=TRSample(FullInterpRVTEnv ,0.5,OrigSamp)
filename = options.Prefix +'_HalfTR_RVT_'+ str(OrigSamp) + ".txt"
save(filename,  HSigTRPt  + abs(min(HSigTRPt)), fmt='%8.6f'  )


#RRT

filename = options.Prefix +'_RRT_'+ str(OrigSamp) + ".txt"
save(filename, RRT  + abs(min(RRT)), fmt='%8.6f' ) 

SigTRfull,SigTRPt, SigTRPtTime=TRSample( RRT,1.0,OrigSamp)
filename = options.Prefix +'_TR_RRT_'+ str(OrigSamp) + ".txt"
save(filename,  SigTRPt  + abs(min(SigTRPt)), fmt='%8.6f') 

HSigTRfull,HSigTRPt, HSigTRPtTime=TRSample(RRT ,0.5,OrigSamp)
filename = options.Prefix +'_HalfTR_RRT_'+ str(OrigSamp) + ".txt"
save(filename,  HSigTRPt  + abs(min(HSigTRPt)), fmt='%8.6f'  )

#Cardiac Peaks
filename = options.Prefix +'_CPd3_'+ str(OrigSamp) + ".txt"
save(filename,  CPd3  , fmt='%8.6f' ) 

filename = options.Prefix +'_CPd3R_'+ str(OrigSamp) + ".txt"
save(filename,  CPd3Rwave  , fmt='%8.6f' ) 


filename = options.Prefix +'_CPttl_'+ str(OrigSamp) + ".txt"
save(filename,  CPd3  , fmt='%8.6f' ) 



#CRTd3 trough

filename = options.Prefix +'_CRTd3_'+ str(OrigSamp) + ".txt"
save(filename, CRTd3  + abs(min(CRTd3)), fmt='%8.6f' ) 

SigTRfull,SigTRPt, SigTRPtTime=TRSample( CRTd3,1.0,OrigSamp)
filename = options.Prefix +'_TR_CRTd3_'+ str(OrigSamp) + ".txt"
save(filename,  SigTRPt  + abs(min(SigTRPt)), fmt='%8.6f') 

HSigTRfull,HSigTRPt, HSigTRPtTime=TRSample(CRTd3 ,0.5,OrigSamp)
filename = options.Prefix +'_HalfTR_CRTd3_'+ str(OrigSamp) + ".txt"
save(filename,  HSigTRPt  + abs(min(HSigTRPt)), fmt='%8.6f'  )

#CRTd3Rwave
filename = options.Prefix +'_CRTd3R_'+ str(OrigSamp) + ".txt"
save(filename, CRTd3R  + abs(min(CRTd3R)), fmt='%8.6f' ) 

SigTRfull,SigTRPt, SigTRPtTime=TRSample( CRTd3R,1.0,OrigSamp)
filename = options.Prefix +'_TR_CRTd3R_'+ str(OrigSamp) + ".txt"
save(filename,  SigTRPt  + abs(min(SigTRPt)), fmt='%8.6f') 

HSigTRfull,HSigTRPt, HSigTRPtTime=TRSample(CRTd3R ,0.5,OrigSamp)
filename = options.Prefix +'_HalfTR_CRTd3R_'+ str(OrigSamp) + ".txt"
save(filename,  HSigTRPt  + abs(min(HSigTRPt)), fmt='%8.6f'  )

#CRTttl

filename = options.Prefix +'_CRTttl_'+ str(OrigSamp) + ".txt"
save(filename, CRTttl  + abs(min(CRTttl)), fmt='%8.6f' ) 

SigTRfull,SigTRPt, SigTRPtTime=TRSample( CRTttl,1.0,OrigSamp)
filename = options.Prefix +'_TR_CRTttl_'+ str(OrigSamp) + ".txt"
save(filename,  SigTRPt  + abs(min(SigTRPt)), fmt='%8.6f') 

HSigTRfull,HSigTRPt, HSigTRPtTime=TRSample(CRTttl ,0.5,OrigSamp)
filename = options.Prefix +'_HalfTR_CRTttl_'+ str(OrigSamp) + ".txt"
save(filename,  HSigTRPt  + abs(min(HSigTRPt)), fmt='%8.6f'  )


#-----------------------------------------------------------------------------


if PlotAll:
	show()

