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




def save_epoch_object(epochs, epochTimes, region, sfreq, subject_outDir, name):
	
	# adjust array to have apprpriate dimensions. shape should be (n_epochs, n_chans, n_times)
	epochs = np.expand_dims(epochs, axis=1)

	#Initialize the info structure
	info = mne.create_info(ch_names=[region], ch_types=['mag'], sfreq=sfreq)

	# what is the begining time of the epoch? (the end is inferred from sampling freq and n_samples)
	epoch_tmin = epochTimes[0]

	# Create an "events" array. shape: (n_events, 3). 
	# first column is index of event, second column is the length of the event, third column is event type	
	num_events = len(epochs) 
	c1 = np.arange(num_events)
	c2 = np.full(num_events, 1)
	c3 = c2
	events = np.stack((c1, c2, c3), axis=1)

	# Create epochs object
	subj_epochs = mne.EpochsArray(epochs, info, events, epoch_tmin)

	# save the epochs object
	subj_epochs.save()

	return subj_epochs 



def burst_triggered_average(S, sfreq, dfs, epoch_length_ms):

    # setting up waveform epochs    
    epochLength = epoch_length_ms/1000*sfreq
    epochTime = epochLength/sfreq
    epochTimes = np.arange(epochLength)/sfreq-(epochTime/2) #time in seconds
    peak_times = dfs['Peak Time'].tolist()
    peak_freqs = dfs['Peak Frequency'].tolist()

    # Trough finding using band pass filtered data
    troughIndexes = dict()
    i=0
    for pt in peak_times:


        # Narrow band-pass filter the time course (2Hz width around peak frequency)
        PF = peak_freqs[i]
        low_cut = PF - 1
        high_cut = PF + 1
        b,a = ss.butter(4, [low_cut, high_cut], 'bandpass', fs=sfreq) 
        S_BP = ss.filtfilt(b,a,np.expand_dims(S,0), padlen=100)

        peakIndex = int(pt*sfreq)
        firstIndex = int(peakIndex-epochLength/2)
        lastIndex = int(peakIndex+epochLength/2)
        thisData = S_BP.T[firstIndex:lastIndex] # shape = (500, 1)
        localMinIndex = ss.argrelextrema(-thisData, np.greater)[0]-epochLength/2
        localMinClosestToZero = min(localMinIndex, key=lambda x:abs(x))
        troughIndexes[pt] = peakIndex + localMinClosestToZero 
        i=i+1


    # Extract unfiltered time course for each event aligned to trough closest to peak time
    epochs = []
    for pt in troughIndexes:

        # get trough index (ti) for this event
        ti = troughIndexes[pt]

        # set up epoch bounds
        firstIndex = int(ti-epochLength/2)
        lastIndex = int(ti+epochLength/2)

        #
        thisData = S[firstIndex:lastIndex]
        epochs.append(thisData)

    J= np.asarray(epochs)
    

    return epochTimes, J




if __name__ == "__main__":

	
	mne.set_log_level(verbose='ERROR')
	
	sfreq = 1000

	# epoch length (in ms). Must be a whole number. Should be longer than 1000 ms.
	epoch_length_ms = 2000

	# channel to analyze
	regions = [	'postcentral_lh', 'postcentral_rh', 
				'precentral_lh', 'precentral_rh']	

	# where to save files (one csv file per participant/region)
	outDir = 

	# Find subjects to be analysed 
	subjectIDs = 




	for region in regions:

		# load it data and delete any event with peak power above 99 percentile (erronious events). also get rid of peak times within 1s of start

		df_all_PAPTO = pd.read_csv()
		df_all_PAPTO = df_all_PAPTO[['subjectID', 'Peak Time', 'Peak Frequency', 'Peak Power']]		
		df_all_PAPTO = df_all_PAPTO[df_all_PAPTO['Peak Time']>2]
		df_all_PAPTO = df_all_PAPTO[df_all_PAPTO['Peak Time']<208]
		NN_percentile = df_all_PAPTO['Peak Power'].mean() + 2.58*df_all_PAPTO['Peak Power'].std()
		df_all_PAPTO = df_all_PAPTO.drop(df_all_PAPTO[df_all_PAPTO['Peak Power'] > NN_percentile].index)


		df_all_Mnorm = pd.read_csv()
		df_all_Mnorm = df_all_Mnorm[['subjectID', 'Peak Time', 'Peak Frequency', 'Peak Power']]
		df_all_Mnorm = df_all_Mnorm[df_all_Mnorm['Peak Time']>2]
		df_all_Mnorm = df_all_Mnorm[df_all_Mnorm['Peak Time']<208]
		NN_percentile = df_all_Mnorm['Peak Power'].mean() + 2.58*df_all_Mnorm['Peak Power'].std()
		df_all_Mnorm = df_all_Mnorm.drop(df_all_Mnorm[df_all_Mnorm['Peak Power'] > NN_percentile].index)

		df_PAPTO_only = pd.read_csv()
		df_PAPTO_only = df_PAPTO_only[['subjectID', 'Peak Time', 'Peak Frequency', 'Peak Power']]
		df_PAPTO_only = df_PAPTO_only[df_PAPTO_only['Peak Time']>2]
		df_PAPTO_only = df_PAPTO_only[df_PAPTO_only['Peak Time']<208]
		NN_percentile = df_PAPTO_only['Peak Power'].mean() + 2.58*df_PAPTO_only['Peak Power'].std()
		df_PAPTO_only = df_PAPTO_only.drop(df_PAPTO_only[df_PAPTO_only['Peak Power'] > NN_percentile].index)

		df_Mnorm_only = pd.read_csv()
		df_Mnorm_only = df_Mnorm_only[['subjectID', 'Peak Time', 'Peak Frequency', 'Peak Power']]
		df_Mnorm_only = df_Mnorm_only[df_Mnorm_only['Peak Time']>2]
		df_Mnorm_only = df_Mnorm_only[df_Mnorm_only['Peak Time']<208]
		NN_percentile = df_Mnorm_only['Peak Power'].mean() + 2.58*df_Mnorm_only['Peak Power'].std()
		df_Mnorm_only = df_Mnorm_only.drop(df_Mnorm_only[df_Mnorm_only['Peak Power'] > NN_percentile].index)

		df_both = pd.read_csv()
		df_both = df_both[['subjectID', 'Peak Time', 'Peak Frequency', 'Peak Power']]
		df_both = df_both[df_both['Peak Time']>2]
		df_both = df_both[df_both['Peak Time']<208]		
		NN_percentile = df_both['Peak Power'].mean() + 2.58*df_both['Peak Power'].std()
		df_both = df_both.drop(df_both[df_both['Peak Power'] > NN_percentile].index)	


		df_all_PAPTO.name='all_PAPTO'
		df_all_Mnorm.name = 'all_Mnorm'
		df_PAPTO_only.name = 'PAPTO_only'
		df_Mnorm_only.name = 'Mnorm_only'
		df_both.name = 'both'


		dfs = [df_all_PAPTO, df_all_Mnorm, df_PAPTO_only, df_Mnorm_only, df_both]
		names = ['all_PAPTO','all_Mnorm', 'PAPTO_only', 'Mnorm_only', 'both']
	


		i=0
		for df in dfs:
			name = names[i]
			i=i+1
			for subjectID in subjectIDs:
				epochFif = 

				if os.path.exists(epochFif):

					print(subjectID + '       ' + region + '      ' + name)

					dfs = df[df['subjectID']==subjectID]

					if len(dfs)>0:
						
						# out directory for this subject
						subject_outDir = 

						epochs = mne.read_epochs(epochFif)
						epoch = epochs.copy().pick(picks=region)
						epochData = np.squeeze(epoch.get_data())
						S = epochData
						
						epochTimes, epochs = burst_triggered_average(S, sfreq, dfs, epoch_length_ms)


						# save all burst-level waveforms within one epochs object per participant. 
						subj_epochs_beta = save_epoch_object(epochs, epochTimes, region, sfreq, subject_outDir, name)
					else:
						print('no events :(')








			 


