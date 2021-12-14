# Import libraries
import os
import sys
import time
import mne
import math
import numpy as np
import scipy.signal as ss
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import logging
import support_functions.spectralevents_functions as tse
import multiprocessing as mp
import warnings
from scipy import stats
from fooof import FOOOF
from fooof import FOOOFGroup
import random
from scipy import signal



def get_spectral_events(TFR, Fs, fVec, tVec,):

	# modified from Shin et al., eLife, 2017

    # Find transient spectral events based on TFR
    findMethod = 1
    thrFOM = 1 # This is irrelevant for this script b/c we dont take outlier events
    numTrials = TFR.shape[0] 
    classLabels = [1 for x in range(numTrials)]
    neighbourhood_size = (4,160)
    threshold = 0 # i.e., no thresholding setting 
    spectralEvents = tse.spectralevents_find(findMethod, thrFOM, tVec,
                fVec, TFR, classLabels, neighbourhood_size, threshold, Fs)
    df = pd.DataFrame(spectralEvents)

    #modify burst propterties
    allEvents = df.copy()
    allEvents['periods'] = allEvents['Event Duration']*allEvents['Peak Frequency']

    allEvents = allEvents.drop(['Trial', 'Hit/Miss', 'Normalized Peak Power', 'Outlier Event'], axis=1)

    return allEvents





def energyvec(f,s,Fs,width):
    
	# Modified from Shin et al., eLife, 2017

    # Return a vector containing the energy as a
    # function of time for frequency f. The energy
    # is calculated using Morlet's wavelets. 
    # s : signal
    # Fs: sampling frequency
    # width : width of Morlet wavelet (>= 5 suggested).


    dt = 1/Fs
    sf = f/width
    st = 1/(2 * np.pi * sf)

    t= np.arange(-3.5*st, 3.5*st, dt)
    m = morlet(f, t, width)

    y = np.convolve(s, m)
    y = (dt * np.abs(y))**2
    lowerLimit = int(np.ceil(len(m)/2))
    upperLimit = int(len(y)-np.floor(len(m)/2)+1)
    y = y[lowerLimit:upperLimit]

    return y

def morlet(f,t,width):

	# Modified from Shin et al., eLife, 2017

    # Morlet's wavelet for frequency f and time t. 
    # The wavelet will be normalized so the total energy is 1.
    # width defines the ``width'' of the wavelet. 
    # A value >= 5 is suggested.
    #
    # Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)

    sf = f/width
    st = 1/(2 * np.pi * sf)
    A = 1/np.sqrt((st/2 * np.sqrt(np.pi))) 
    y = A * np.exp(-t**2 / (2 * st**2)) * np.exp(1j * 2 * np.pi * f * t)

    return y



def interpolate_50Hz_notch(PSD, PSD_fVec):


	#linear interpolate over 50 Hz notch filter
	m = (PSD[0][0][205]-PSD[0][0][193])/(51.5-48.5)
	b = PSD[0][0][205] - m*51.5
	for i in range(193,206):
		PSD[0][0][i] = m*PSD_fVec[i]+b

	return PSD


def TFR_via_morlet_wavelet(epoch, fVec):

	# Adapted from spectralevents_ts2tfr function from Shin et al (2017)


	width = 10

	# get timecourse s(t) from epoch 
	# s should have dimensions: (#timepoints, 1) 
	s = np.squeeze(epoch.get_data(), axis=0)

	# obtain time vector (tVec) from timecourse (tVec starting with t=0s)
	numSamples = s.shape[1] 
	tVec = np.arange(numSamples)/epoch.info['sfreq']

	# find number of frequencies for convolution
	numFrequencies = len(fVec)

	# generate TFR row by row
	TFR = []
	B = np.zeros((numFrequencies, numSamples))
	# Frequency loop
	for j in np.arange(numFrequencies):
		B[j,:] = energyvec(fVec[j], signal.detrend(s[0,:]), epoch.info['sfreq'], width)
	TFR.append(B)

	return TFR, tVec








def find_PAPTO_bursts(epoch, outDir, region, subjectID):


	#################################################################
	# [1] compute TFR via Morlet wavelet convolution
	# technique based on: Tallon-Baudry et al., J. Neuroscience, 1997
	# code adapted from: Shin et al, eLife, 2017
	##################################################################
	
	fmin = 0.25
	fmax = 80
	fstep = 0.25
	fVec = np.arange(fmin, fmax, fstep)

	TFR, tVec = TFR_via_morlet_wavelet(epoch, fVec)


	###################################################################
	# [2] FOOOF modeling to obtain offset and exponent 
	####################################################################

	# PSD via welch method
	PSD, PSD_fVec = mne.time_frequency.psd_welch(epoch, fmin=fmin, fmax=fmax,  
								n_fft=1000, n_overlap=900, picks='all')

	PSD = interpolate_50Hz_notch(PSD, PSD_fVec)

	# FOOOF modeling
	fm = FOOOF(peak_width_limits=(2,10), max_n_peaks=4, aperiodic_mode='fixed', min_peak_height=0.05, peak_threshold=1.5)
	fm.fit(PSD_fVec, np.squeeze(PSD), [fmin,fmax])
	fm.save_report(outDir + subjectID + '/FOOOF_report' + region + '_' + subjectID + '.pdf')
	exponent = fm.get_params('aperiodic_params', 'exponent')
	offset = fm.get_params('aperiodic_params', 'offset')


	###################################################################
	# [3] calculate TFR normalization factor and apply it to TFR
	####################################################################

	# generate normalization factor n_c from aperiodic offset and exponent
	nc = (fVec**(exponent)) / (10**offset)

	# normalize the TFR
	TFR = TFR * nc[None,:,None] 


	###################################################################
	# [4] find events from normalized TFR
	####################################################################

	allEvents = get_spectral_events(TFR, epoch.info['sfreq'], fVec, tVec)









if __name__ == "__main__":

	
	mne.set_log_level(verbose='ERROR')
	
	# channel to analyze
	regions = [	'postcentral_lh', 'postcentral_rh', 
				'precentral_lh', 'precentral_rh']	

	# where to save files (one csv file per participant/region)
	outDir = 


	# Flist of subjectIDs to be analyzed
	subjectIDs = 

	############ HACK FOR TESTING ##############################################
	#subjectIDs = subjectIDs[0:3]
	#regions = regions[0:1]
	############################################################################



	for subjectID in subjectIDs:

		print(subjectID)

		for region in regions:

			# directory where source localized epoch.fif files are located
			epochFif = 
			
			if os.path.exists(epochFif):

				epochs = mne.read_epochs(epochFif)
				epoch = epochs.copy().pick(picks=region)

				# downsample the epoch
				epoch = epoch.resample(250, npad='auto') 

				find_PAPTO_bursts(epoch, outDir, region, subjectID)







			 


