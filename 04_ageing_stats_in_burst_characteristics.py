'''
########################################################################################
################################### NOTES ##############################################

This script is designed to generated ageing models for papto and med-norm burst characteristics

Inputs (for each participant):
- two events lists (csv): one papto events and one med-norm events

Outputs (for each participant):
 - ageing model stats in csv format
 - scatterplots (linear and quadratic) of age-related changes in each event characteristic

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
#import spectralevents_functions as tse
import multiprocessing as mp
import warnings
from scipy import stats
from fooof import FOOOF
import scipy.stats as stats
import statsmodels.formula.api as smf
import statistics


def calc_l_fit_coefficients(x, y, characteristic, weights):
	
	degree = 1

	coefs = np.polyfit(x, y, degree, w=weights)
	linear_parameters = ['characteristic', 'linear', 'intercept']
	df_linear_fit = pd.DataFrame(columns=linear_parameters)
	linear_values = [characteristic, coefs[0], coefs[1]]
	dictionary = dict(zip(linear_parameters, linear_values))
	df_linear_fit = df_linear_fit.append(dictionary, ignore_index=True)	

	return df_linear_fit


def calc_q_fit_coefficients(x, y, characteristic, weights):
	
	degree = 2
	coefs = np.polyfit(x, y, degree, w=weights)
	quadratic_parameters = ['characteristic', 'quadratic', 'linear', 'intercept']
	df_quadratic_fit = pd.DataFrame(columns=quadratic_parameters)
	quadratic_values = [characteristic, coefs[0], coefs[1], coefs[2]]
	dictionary = dict(zip(quadratic_parameters, quadratic_values))
	df_quadratic_fit = df_quadratic_fit.append(dictionary, ignore_index=True)		

	return df_quadratic_fit





def calc_chi_square_l(parameter, df, coef, region, variances):	
	
	chi_square = 0	
	for index, subject in df.iterrows():
		real = subject[c]
		pred = coef['linear'][0]*subject['age'] + coef['intercept'][0]
		subjectID = index
		var = variances[c][subjectID] 
		w = ((real-pred)**2)/var
		chi_square = w + chi_square

	return chi_square



def calc_chi_square_q(parameter, df, coef, region, variances):	
	
	chi_square = 0	
	for index, subject in df.iterrows():
		real = subject[c]
		pred = coef['quadratic'][0]*subject['age']**2 + coef['linear'][0]*subject['age'] + coef['intercept'][0]
		subjectID = index
		var = variances[c][subjectID] 
		w = ((real-pred)**2)/var
		chi_square = w + chi_square

	return chi_square








def calc_F_stat(N, chi_sq_l, chi_sq_q):

	'''
	###################
	NOTES:
	- see page 215 (eq. 13.5) of Bonamente, Statistics and Analysis of Scientific Data, 2nd ed. (2017)
	- Use f1 and f2 to get corresponding critical values of F in table A.8
	#####################
	'''

	m = 3 # m = num. of adjustable parameters in quad model (i.e. a,b,c)
	p = 1 # p = difference in num. of adjustable parameters between the two models (i.e. we fixed one parameter, a=0)
	f1 = p # degrees of freedom
	f2 = N-m # degrees of freedom 
	num = (chi_sq_l-chi_sq_q)/p 
	den = chi_sq_q/(N-m)
	F = num/den

	return F, f1, f2








if __name__ == "__main__":
	
	mne.set_log_level(verbose='ERROR')
	
	# channel to analyze
	regions = [	'postcentral_lh', 'postcentral_rh', 
				'precentral_lh', 'precentral_rh']	

	# Find subjects to be analysed a
	subjectIDs = subjectData['SubjectID'].tolist() 


	# where to save epoch.fif files (one file per participant)
	outDir = 

	columns = ['region', 'parameter', 'background_sub', 'F-stat', 'appropriate model', 'p_val_CORRECTED',
				'quad', 'lin', 'int']

	stat_tracking = pd.DataFrame(columns=columns)

	for region in regions:

		print(region)

		if region == 'postcentral_rh':
			color = 'lightsteelblue'
		if region == 'precentral_rh':
			color = 'lightsteelblue'
		if region == 'precentral_lh':
			color = 'lightsteelblue'
		if region == 'postcentral_lh':
			color = 'lightsteelblue'

		# load in participant-level characteristics (with=PAPTO, without=med-norm)
		df_with = 
		df_without = 

		df_with['age'] = df_with['subjectID'].map(Ages_dict)
		df_without['age'] = df_without['subjectID'].map(Ages_dict)
		
		# frequency resolution is 0.25Hz whereas event characteristic assumed it was 1Hz
		df_with['Frequency Span'] = df_with['Frequency Span']/4
		df_without['Frequency Span'] = df_without['Frequency Span']/4

		# get event duration in ms
		df_with['Event Duration'] = df_with['Event Duration']*1000
		df_without['Event Duration'] = df_without['Event Duration']*1000

		types = [df_with, df_without]
		characteristics = ['Peak Frequency', 'Frequency Span', 'Event Duration', 'periods']

		b=0
		for df in types:
			b=b+1
			if b==1:
				bs = 'PAPTO'
			if b==2:
				bs = 'Mnorm'

			for c in characteristics:

				print(c)
	
				df_c = df[['subjectID', 'age', c]]
				df_mean = df_c.groupby('subjectID').mean()
				df_var = df_c.groupby('subjectID').var()
				df_var[c] = df_var[c].replace({0:np.nan}) # replace 0 with nan
				df_var = df_var.fillna(df_var.mean()) # replace nans with mean
				df_std = df_var.apply(np.sqrt) 
				weights = 1/df_std[c]

				N = len(df_mean)

							
	
				# get fit coefficients
				df_linear_fit = calc_l_fit_coefficients(df_mean['age'], df_mean[c], c, weights)
				df_quadratic_fit = calc_q_fit_coefficients(df_mean['age'], df_mean[c], c, weights)								

				# F - test for model comparison
				chi_sq_l = calc_chi_square_l(c, df_mean, df_linear_fit, region, df_var)
				chi_sq_q = calc_chi_square_q(c, df_mean, df_quadratic_fit, region, df_var)



				F, f1, f2 = calc_F_stat(N, chi_sq_l, chi_sq_q)
				# test the null hypothesis
				if F>6.635000:
					appropriate_model = 'quadratic'
				if F<6.635000:
					appropriate_model = 'linear'
				

				'''
				
				H0: linear and quadratic are quivalent models to explain the data
				our question, at the 95% confidence level, is the additional model component requried

				NOTE: In our case: 
				f1 = 1    and    f2 = 594-3 = 591 (f2 is effectively infinite)
				Based on table A.8, this corresponds to:
				F_critical - condience level - prob that fit improvement is due by chance
				0.455 - 0.50 - 0.50
				0.708 - 0.60 - 0.40
				1.074 - 0.70 - 0.30
				1.642 - 0.80 - 0.20
				2.706 - 0.90 - 0.10
				3.842 - 0.95 - 0.05 --- This one
				6.635 - 0.99 - 0.01 
				>6.635 - >0.99 - <0.01

				high F-stat means the additional model component is important
				1-CL = probability that fit improvement is due by chance

				examples:
				ex. if F = 6.888, very very likely there is a real differnce in the models! There is <0.01 probability that difference is just by chance
				ex. if F = 2.888, likely there is a real differnce in the models. There is <0.10 but >0.05 probability that difference is just by chance
				ex. if F = 0.565, who knows if there is a difference. There is <0.50 but >0.40 probability that difference is just by chance
				'''
				

				###############################################################################################################
				# GET P VALS
				###############################################################################################################


				degree=1
				coefs = np.polyfit(df_mean['age'].to_numpy(), df_mean[c].to_numpy(), degree, w=weights.to_numpy())
				df_s=pd.DataFrame(columns=['y','x'])
				df_s['x'] = df_mean['age'].to_numpy()
				df_s['y'] = df_mean[c].to_numpy()
				model = np.poly1d(coefs) 
				results = smf.ols(formula ='y ~ model(x)', data=df_s).fit()
				p_val_lin = results.pvalues.loc['model(x)']


				degree=2
				coefs = np.polyfit(df_mean['age'].to_numpy(), df_mean[c].to_numpy(), degree, w=weights.to_numpy())
				df_s=pd.DataFrame(columns=['y','x'])
				df_s['x'] = df_mean['age'].to_numpy()
				df_s['y'] = df_mean[c].to_numpy()
				model = np.poly1d(coefs) 
				results = smf.ols(formula ='y ~ model(x)', data=df_s).fit()
				p_val_quad = results.pvalues.loc['model(x)']


				###############################################################################################################
				###############################################################################################################


				if appropriate_model == 'linear':
					quad = 0.0
					lin = df_linear_fit['linear'][0]
					intercept = df_linear_fit['intercept'][0]
					p_val = p_val_lin
				if appropriate_model == 'quadratic':
					quad = df_quadratic_fit['quadratic'][0]
					lin = df_quadratic_fit['linear'][0]
					intercept = df_quadratic_fit['intercept'][0]
					p_val = p_val_quad

				stat_tracking = stat_tracking.append({'region':region, 
														'parameter': c, 
														'background_sub':bs,
														'F-stat':F, 
														'appropriate model':appropriate_model, 
														'p_val_CORRECTED':p_val*26,
														'quad':quad,
														'lin':lin,
														'int':intercept}, ignore_index=True)



				
				sns.set(style='white', palette='deep' ,font_scale=1.5, rc={'figure.figsize':(3,5)})   # (2,5)
				sns.set_style('ticks', {'xtick.major.size':8})

				# keep only outlier events
				if c == 'Event Duration':
					ymin = 140
					ymax = 300
				if c == 'Frequency Span':
					ymin = 5
					ymax = 8
				if c == 'Peak Frequency':
					ymin = 19
					ymax = 25
				if c == 'periods':
					ymin = 4
					ymax = 5.5





				# generate quadratic fit plot - exponent
				sns.regplot(data=df_mean, x='age', y=c, ci=95, fit_reg=True, scatter=True, order=2, scatter_kws={'s':15, 'color':color}, line_kws={'linewidth':3.5, 'color':'dimgray'})#'limegreen'}) 
				axes = plt.gca()
				sns.despine()
				axes.set_ylim([ymin,ymax])
				axes.set_xlim([8,98])
				axes.set_xlabel('Age (years)')
				axes.set_ylabel(c)
				axes.set_xticks([20,50,80], minor=False)
				#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
				#axes.set_yticks([])
				axes.tick_params(direction='out')
				plt.tight_layout()

				plt.savefig()
				plt.close()


				# generate quadratic fit plot - exponent
				sns.regplot(data=df_mean, x='age', y=c, ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':15, 'color':color}, line_kws={'linewidth':3.5, 'color':'dimgray'})#'limegreen'}) 
				axes = plt.gca()
				sns.despine()
				axes.set_ylim([ymin,ymax])
				axes.set_xlim([8,98])
				axes.set_xlabel('Age (years)')
				axes.set_ylabel(c)
				axes.set_xticks([20,50,80], minor=False)
				#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
				#axes.set_yticks([])
				axes.tick_params(direction='out')
				plt.tight_layout()

				plt.savefig()
				plt.close()


	stat_tracking.to_csv(()








































