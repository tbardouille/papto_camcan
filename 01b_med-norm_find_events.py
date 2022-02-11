'''
########################################################################################
################################### NOTES ##############################################

This script is designed find transient events using an adapted Shin, 2017 algorithm.
Reference: eLife 2017;6:e29086 doi: 10.7554/eLife.29086

Inputs (for each participant):
 - epoch.fif file (one per participant) with four timecourses (M1 & S1 - left & right)

Outputs (for each participant):
 - list (CSV) of transient events detect using the Shin 2017 algorithm.

########################################################################################
'''

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

def get_spectral_events(TFR, sfreq, tVec):

    # Find transient spectral events based on TFR
    findMethod = 1
    thrFOM = 1 # This is irrelevant for this script b/c we dont take outlier events
    numTrials = TFR.shape[0] 
    classLabels = [1 for x in range(numTrials)]
    neighbourhood_size = (4,160)
    threshold = 0.00
    spectralEvents = tse.spectralevents_find(findMethod, thrFOM, tVec,
                fVec, TFR, classLabels, neighbourhood_size, threshold, sfreq)
    df = pd.DataFrame(spectralEvents)


    allEvents = df.copy()
    allEvents['periods'] = allEvents['Event Duration']*allEvents['Peak Frequency']

    return allEvents


def get_TFR(S, sfreq, fmin, fmax, fstep):

    # Make the TFR for transient spectral event analysis 
    fVec = np.arange(fmin, fmax+1, fstep)
    width = 10
    TFR, tVec, fVec = tse.spectralevents_ts2tfr(np.expand_dims(S,1), fVec, sfreq, width) 

    return TFR, fVec, tVec


if __name__ == "__main__":
	
	mne.set_log_level(verbose='ERROR')

	regions = [	'postcentral_lh', 'postcentral_rh', 
				'precentral_lh', 'precentral_rh']	
	

	# re-sampling frequency in Hz
	resamp_freq = 250
	sfreq = 250

	# Analysis frequency range (for TFR generation and PSD analysis). measured in Hz.
	# It should be at least 1-60 Hz (wider ranges acceptable) for accurate PSD modeling
	fmin = 0.25
	fmax = 80
	fstep = 0.25


	# Rest times (s)
	tmin = -105
	tmax = 105

	# where to save epoch.fif files (one file per participant)
	outDir = 

	# Find subjects to be analysed 
	subjectIDs = 




	for subjectID in subjectIDs:
		
		epochFif = 

		if os.path.exists(epochFif):

			# load in cleaned epoch fif and extract timecourse data
			epochs = mne.read_epochs(epochFif)

			for region in regions:

				print(subjectID + '--' + region)

				epochs_copy = epochs.copy()
				region_epochs = epochs_copy.pick_channels([region])
				region_epochs = region_epochs.resample(resamp_freq, npad='auto')  
				epochData = np.squeeze(region_epochs.get_data())

				# make output file folder if it doesnt already exist
				subject_outDir =  outDir + subjectID + '/'
				if not os.path.exists(subject_outDir):
					os.makedirs(subject_outDir)

				# create TFR
				TFR, fVec, tVec = get_TFR(epochData, resamp_freq, fmin, fmax, fstep) #TFR shape = (1, 221, 210001)

				#create AverageTFR object and save
				TFR_obj = mne.time_frequency.AverageTFR(region_epochs.info, TFR, tVec, fVec, 1)

				
				# get TFR data
				TFR = TFR_obj.data
				fVec = TFR_obj.freqs
				tVec = TFR_obj.times
				
				#find spectral events
				allEvents = get_spectral_events(TFR, sfreq, tVec)
				subject_outDir =  outDir + subjectID + '/'
				allEvents.to_csv()


	

	







	




