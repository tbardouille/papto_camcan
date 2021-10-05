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
import random


def truncate(n, decimals):
	multiplier = 10 ** decimals
	return int(n * multiplier) / multiplier



def FOOOF_background(TFR, fVec, fmin, fmax, subjectID, beta_min, beta_max, subject_outDir, region):
	

	# Generate power spectral density (PSD) from TFR 
	TFR_PSD = np.mean(np.mean(TFR, axis=0), axis=1)/fVec 

	# Run FOOOF algorithm and fit to PSD
	fm = FOOOF(peak_width_limits=(1,25), max_n_peaks=4, aperiodic_mode='fixed', min_peak_height=0.0, peak_threshold=1.0)
	freq_range = [fmin,fmax]
	fm.fit(fVec, TFR_PSD, freq_range)

	# save FOOOF report
	fm.save_report(subject_outDir + '/FOOOF_report_'+region+'_loglog', plt_log=True)
	fm.save_report(subject_outDir + '/FOOOF_report_'+region+'_semilog', plt_log=False)

	# get aperiodic fit parameters (note: 'knee' parameter should be included for unfixed model)
	offset = fm.get_params('aperiodic_params', 'offset')
	exponent = fm.get_params('aperiodic_params', 'exponent')

	# get fitted parameter values for the peaks. Each row is a peak, as [center frequency, power, bandwidth]
	periodic_params = fm.get_params('peak_params') 

	# get gaussian fit parameters. These parameters define gaussian fit(s). Each row is a gaussian, as [mean, height, standard deviation]
	gaussian_params = fm.get_params('gaussian_params')

	# get R-squared of the fit between the input power spectrum and thee full model fit
	r_squared = fm.get_params('r_squared')
	
	# get error of the full model fit
	error = fm.get_params('error')

	# get number of peaks that were fit in the model
	n_peaks = fm.n_peaks_

	return offset, exponent, periodic_params, gaussian_params, r_squared, error, n_peaks


def bootstrappin(TFR_obj, subject_outDir, region):

	######################## 
	'''
	This function is used to estimate variance for F-test in script 01a

	#crop TFR to a 50s chunks base on random number


	'''
	########################


	# sub list of subjectIDs that are being analyzed
	bootstrap_parameters = pd.DataFrame(columns=['subjectID','Age','region', 'exponent','offset', 'r_squared', 'error'])

	#make epochs TFR from average TFR
	info = TFR_obj[0].info
	data = TFR_obj[0].data
	data = np.expand_dims(data, axis=0)
	fVec = TFR_obj[0].freqs
	tVec = TFR_obj[0].times
	TFR = mne.time_frequency.EpochsTFR(info, data, tVec, fVec)


	for i in range(0,25):

		TFR_copy = TFR[0].copy()
		r = random.random()
		tmin = r*150.0
		tmax = tmin+50.0
		TFR_chunk = TFR_copy.crop(tmin=tmin, tmax=tmax)

		TFR_data = TFR_chunk[0].data
		TFR_data = np.mean(TFR_data, axis=0)
		fVec = TFR_chunk[0].freqs
		tVec = TFR_chunk[0].times


		# find fit parameters via FOOOF 
		offset, exponent, periodic_params, gaussian_params, r_squared, error, n_peaks = FOOOF_background(TFR_data, fVec, fmin, fmax, subjectID, beta_min, beta_max, subject_outDir, region)	

		# append baseline parameters to output
		bootstrap_parameters = bootstrap_parameters.append({'subjectID':subjectID, 'Age':Ages_dict[subjectID], 'region':region,
															 'exponent':exponent, 
															 'offset':offset, 
															 'r_squared': r_squared,
															 'error':error},
															 ignore_index=True)
	
		bootstrap_parameters.to_csv(subject_outDir + '/' +region+ 'FOOOF_bootstrap_parameters.csv')




if __name__ == "__main__":
	
	mne.set_log_level(verbose='ERROR')
	
	# channel to analyze
	regions = [	'postcentral_lh', 'postcentral_rh', 
				'precentral_lh', 'precentral_rh']	
	
	# sampling frequency in Hz
	sfreq = 250

	# Analysis frequency range (for TFR generation and PSD analysis). measured in Hz.
	# It should be at least 1-60 Hz (wider ranges acceptable) for accurate PSD modeling
	fmin = 0.5
	fmax = 60
	fstep = 0.5

	beta_min=15
	beta_max=30

	# Rest times (s)
	tmin = -105
	tmax = 105

	# where to save files (one file per participant)
	outDir = '/media/NAS/bbrady/Burst_Triggered_Waveforms/anatomical_ROI_noise_corrected_results/'

	# Find subjects to be analysed and take only subjects with more than 55 task epochs
	homeDir = os.path.expanduser("~")
	dataDir = '/home/timb/camcan/'
	camcanCSV = dataDir + 'proc_data/oneCSVToRuleThemAll.csv'
	subjectData = pd.read_csv(camcanCSV)
	subjectData = subjectData[subjectData['numEpochs'] > 55]
	subjectIDs = subjectData['SubjectID'].tolist() 
	subjectAges = subjectData['Age_x'].tolist()
	Ages_dict = dict(zip(subjectIDs, subjectAges))

	############ HACK FOR TESTING ##############################################
	#subjectIDs = subjectIDs[0:10]
	############################################################################

	# sub list of subjectIDs that are being analyzed
	baseline_parameters = pd.DataFrame(columns=['subjectID','Age', 'region', 'exponent','offset'])
	periodic_parameters = pd.DataFrame(columns=['subjectID','Age', 'region', 'peak_num', 'CF','PW', 'BW', 'mean', 'height', 'standard_deviation'])
	fit_parameters = pd.DataFrame(columns=['subjectID','Age', 'region', 'r_squared', 'error'])

	for subjectID in subjectIDs:

		for region in regions:

			print(subjectID + '--__--' + region)
		

			# set TFR directory
			tfrFilename = 'rest-phase_'+ region + '_original-tfr.h5'
			tfr_file = os.path.join('/media/NAS/bbrady/Burst_Triggered_Waveforms/anatomical_ROI_noise_corrected_results/', subjectID, tfrFilename)
			subject_outDir = os.path.join('/media/NAS/bbrady/Burst_Triggered_Waveforms/anatomical_ROI_noise_corrected_results/', subjectID)

			if os.path.exists(tfr_file):

				# load in cleaned epoch fif and extract timecourse data
				TFR_obj = mne.time_frequency.read_tfrs(tfr_file)

				#bootstrap to estimate variance for F-test
				bootstrappin(TFR_obj, subject_outDir, region)

				
				# get TFR data
				TFR = TFR_obj[0].data
				fVec = TFR_obj[0].freqs
				tVec = TFR_obj[0].times
				
				# find fit parameters via FOOOF 
				offset, exponent, periodic_params, gaussian_params, r_squared, error, n_peaks = FOOOF_background(TFR, fVec, fmin, fmax, subjectID, beta_min, beta_max, subject_outDir, region)	
				

				# append baseline parameters to output
				baseline_parameters = baseline_parameters.append({'subjectID':subjectID, 'Age':Ages_dict[subjectID], 'region':region,
																	 'exponent':exponent, 
																	 'offset':offset}, ignore_index=True)


				
				# append periodic parameters to output
				for peak in range(n_peaks):
					periodic_parameters = periodic_parameters.append({'subjectID':subjectID, 'Age':Ages_dict[subjectID], 'region':region,
																		'peak_num':peak+1,
																		'CF':periodic_params[peak][0],
																		'PW':periodic_params[peak][1], 
																		'BW':periodic_params[peak][2], 
																		'mean':gaussian_params[peak][0], 
																		'height':gaussian_params[peak][1], 
																		'standard_deviation':gaussian_params[peak][2]}, ignore_index=True)
		


				# append fit parameters to output
				fit_parameters = fit_parameters.append({'subjectID':subjectID, 'Age':Ages_dict[subjectID], 'region':region,
															'r_squared':r_squared, 
															'error':error}, ignore_index=True)




	baseline_parameters.to_csv(outDir + 'aperiodic_parameters_by_subject.csv')
	periodic_parameters.to_csv(outDir + 'periodic_parameters_by_subject.csv')
	fit_parameters.to_csv(outDir + 'fit_parameters_by_subject.csv')








	




