import numpy as np
import matplotlib.pyplot as plt
import os
import mne
from mne.time_frequency import csd_morlet
from mne.minimum_norm import make_inverse_operator, apply_inverse
from mne.minimum_norm import apply_inverse_epochs, read_inverse_operator
from mne.minimum_norm import apply_inverse
import pandas as pd
import seaborn as sns
import multiprocessing as mp
import math
import statistics






def main_function(subjectID):

	# pick ROI
	regions = [	'postcentral_lh', 'postcentral_rh', 
				'precentral_lh', 'precentral_rh',]
	hemispheres = ['lh', 'rh', 'lh', 'rh']
	ROI_indices = [54,55,56,57]
	ROI_names = dict(zip(ROI_indices, regions))
	hemispheres = dict(zip(ROI_indices, hemispheres))

	# set some things
	sfreq=1000
	snr = 3.0
	lambda2 = 1.0 / snr**2
	method = 'dSPM'
	subjectsDir = '/home/timb/camcan/subjects/'

	# set rest data directories
	dataDir_1 = '/media/NAS/lpower/camcan/spectralEvents/rest/proc_data/'
	epochFifFilename = 'transdef_mf2pt2_rest_raw_rest_210s_cleaned-epo.fif'
	epochFif = os.path.join(dataDir_1, subjectID, epochFifFilename)

	# set some more directories
	transFif = subjectsDir + 'coreg/sub-' + subjectID + '-trans.fif'
	srcFif = subjectsDir + 'sub-' + subjectID + '/bem/sub-' + subjectID + '-5-src.fif'
	bemFif = subjectsDir + 'sub-' + subjectID + '/bem/sub-' + subjectID + '-5120-bem-sol.fif'
	emptyroomFif = '/media/NAS/lpower/BetaSourceLocalization/emptyroomData/' + subjectID + '/emptyroom_trans-epo.fif'
	
	if os.path.exists(epochFif) & os.path.exists(srcFif) & os.path.exists(emptyroomFif):	

		print(subjectID)

		# load in task data and average over trials 
		epochs = mne.read_epochs(epochFif)

		#load in source space
		src = mne.read_source_spaces(srcFif)

		# make forward solution 
		forward = mne.make_forward_solution(epochs.info,
		                                    trans=transFif, src=srcFif, bem=bemFif,
		                                    meg=True, eeg=False)

		# generate noise covariance matirx (306 x 306)
		epochs_emptyroom = mne.read_epochs(emptyroomFif)
		noise_cov= mne.compute_covariance(epochs_emptyroom) 

		# make inverse operator
		inverse_operator = make_inverse_operator(epochs.info, forward, noise_cov)

		#Pull the labels for all anatomically annotated regions & pick the parcellation of interest
		labels = mne.read_labels_from_annot('sub-'+subjectID, parc='aparc.a2009s', subjects_dir=subjectsDir)
		timecourses = []
		for ROI_index in ROI_indices:
			label = labels[ROI_index] 
			restricted_label = label.restrict(src)
			COMvertex = restricted_label.center_of_mass(subject='sub-'+subjectID, subjects_dir=subjectsDir, restrict_vertices=True)
			COMlabel = mne.Label([COMvertex], hemi=hemispheres[ROI_index], subject='sub-'+subjectID, name='COM_label')

			# compute the source estimate 
			ROI_VectorSourceEstimate = apply_inverse_epochs(epochs, inverse_operator, lambda2, method, 
															label, pick_ori='vector' , nave=1) 


			# Get source estimate object restricted to COMlabel
			ROI_COM_VectorSourceEstimate = ROI_VectorSourceEstimate[0].in_label(COMlabel)

			# Get source estimate normal to cortical surface
			NormalSourceEstimate = ROI_COM_VectorSourceEstimate.normal(src, use_cps=True) 

			# extract timecourse of source estimate normal to cortical surface
			COM_tc = NormalSourceEstimate.extract_label_time_course(COMlabel, src, mode='mean_flip', allow_empty=False)
			timecourses.append(np.expand_dims(COM_tc, axis=0))



		# create & save epoch object from COM timecourse
		events = np.array([[0,0,0]])
		info = mne.create_info(ch_names=regions, ch_types=['misc','misc','misc','misc'], sfreq=sfreq)
		tmin = -105
		#  data must be in shape of (n_epochs, n_channels, n_samples) = (1,4,210001)
		data = np.concatenate(timecourses, axis=1)
		epoch = mne.EpochsArray(data, info, events, tmin)


		# make output file folder if it doesnt already exist
		subject_outDir =  '/media/NAS/bbrady/Burst_Triggered_Waveforms/anatomical_ROI_noise_corrected_results/' + subjectID+'/'
		if not os.path.exists(subject_outDir):
			os.makedirs(subject_outDir)

		epoch.save(subject_outDir +subjectID+'_rest_anatomical_ROIs-epo.fif', overwrite=True)


		return






if __name__ == "__main__":

	mne.set_log_level(verbose='ERROR')

	# Find subjects to be analysed and take only subjects with more than 55 task epochs
	homeDir = os.path.expanduser("~")
	dataDir = homeDir + '/camcan/'
	camcanCSV = '/home/timb/camcan/proc_data/oneCSVToRuleThemAll.csv'
	subjectData = pd.read_csv(camcanCSV)
	subjectData = subjectData[subjectData['numEpochs'] > 55]
	subjectIDs = subjectData['SubjectID'].tolist() 

	############ HACK FOR TESTING ##############################################
	#subjectIDs = subjectIDs[0:5]
	############################################################################

	for subjectID in subjectIDs:
		main_function(subjectID)

