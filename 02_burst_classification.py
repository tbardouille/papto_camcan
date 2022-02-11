
'''
########################################################################################
################################### NOTES ##############################################

This script is designed to classify events as either papto-only, med-norm-only, or shared
(i.e. detected by both methods)

Inputs (for each participant):
- two events lists (csv): one papto events and one med-norm events

Outputs (for each participant):
 - same as the inputs but with an extra column indicating what category the event falls into

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
#import peakutils
from fooof import FOOOF
import random
from scipy import signal








if __name__ == "__main__":

	mne.set_log_level(verbose='ERROR')

	# list of subjects to be analysed 
	subjectIDs = 

	regions = [	'postcentral_lh', 'postcentral_rh', 
				'precentral_lh', 'precentral_rh']		

	i=0
	for subjectID in subjectIDs:
	
		epochFif = '/media/NAS/bbrady/PAPTO_manuscript/' + subjectID + '/' + subjectID + '_rest_anatomical_ROIs-epo.fif'
		
		if os.path.exists(epochFif):
			
			print(subjectID)
	

			for region in regions:

				print(region)

				# list of detected events for this participant/ROI
				csvFile_Mnorm = 
				csvFile_PAPTO = 
				
				df_PAPTO = pd.read_csv(csvFile_PAPTO)
				df_PAPTO = df_PAPTO[df_PAPTO['Peak Power']>8]
				df_PAPTO = df_PAPTO[df_PAPTO['Peak Frequency']>15]
				df_PAPTO = df_PAPTO[df_PAPTO['Peak Frequency']<30]
				df_PAPTO['burst category'] = 0.0

				df_Mnorm= pd.read_csv(csvFile_Mnorm)
				df_Mnorm = df_Mnorm[df_Mnorm['Normalized Peak Power']>5]
				df_Mnorm = df_Mnorm[df_Mnorm['Peak Frequency']>15]
				df_Mnorm = df_Mnorm[df_Mnorm['Peak Frequency']<30]
				df_Mnorm['burst category'] = 0.0

				for Pindex, Pburst in df_PAPTO.iterrows():				
					
					fmin = Pburst['Lower Frequency Bound']
					fmax = Pburst['Upper Frequency Bound']
					tmin = Pburst['Event Onset Time']
					tmax = Pburst['Event Offset Time']

					for Mindex, Mburst in df_Mnorm.iterrows():

						PT = Mburst['Peak Time']
						PF = Mburst['Peak Frequency']

						if ((PT<tmax) & (PT>tmin) & (PF<fmax) & (PF>fmin)):

							df_Mnorm.loc[Mindex, 'burst category'] = 1.0
							df_PAPTO.loc[Pindex, 'burst category'] = 1.0



				# save list of events with categorization
				df_PAPTO.to_csv()
				df_Mnorm.to_csv()

							



	
