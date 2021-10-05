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
		real = subject[parameter]
		pred = coef['linear'][0]*subject['Age'] + coef['intercept'][0]
		var = subject['variance'] 
		w = ((real-pred)**2)/var
		chi_square = w + chi_square

	return chi_square



def calc_chi_square_q(parameter, df, coef, region, variances):	
	
	chi_square = 0	
	for index, subject in df.iterrows():
		real = subject[parameter]
		pred = coef['quadratic'][0]*subject['Age']**2 + coef['linear'][0]*subject['Age'] + coef['intercept'][0]
		var = subject['variance'] 
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
	regions = ['postcentral_rh']
	# Find subjects to be analysed and take only subjects with more than 55 task epochs
	homeDir = os.path.expanduser("~")
	dataDir = '/home/timb/camcan/'
	camcanCSV = dataDir + 'proc_data/oneCSVToRuleThemAll.csv'
	subjectData = pd.read_csv(camcanCSV)
	subjectData = subjectData[subjectData['numEpochs'] > 55]
	subjectIDs = subjectData['SubjectID'].tolist() 
	subjectAges = subjectData['Age_x'].tolist()
	Ages_dict = dict(zip(subjectIDs, subjectAges))


	# load in participant-level waveform characteristics
	df_aperiodic = pd.read_csv('/media/NAS/bbrady/Burst_Triggered_Waveforms/anatomical_ROI_noise_corrected_results/aperiodic_parameters_by_subject.csv')
	df_fit = pd.read_csv('/media/NAS/bbrady/Burst_Triggered_Waveforms/anatomical_ROI_noise_corrected_results/fit_parameters_by_subject.csv')	
	df_component_og = pd.concat([df_aperiodic, df_fit['r_squared'], df_fit['error']], axis=1)
	
	# total num. of participants to analyze
	N = len(df_component_og) 
	
	# where to save files (one file per participant)
	outDir = '/media/NAS/bbrady/Burst_Triggered_Waveforms/anatomical_ROI_noise_corrected_results/plots/'

	# set of FOOOF parameters to analyze
	parameters = ['exponent', 'offset']#, 'r_squared', 'error']


	columns = ['region', 'parameter', 'chi_sq_l', 'chi_sq_q', 'F-stat', 'appropriate model', 'p_val_lin', 'p_val_quad', 
				'quad-quad', 'quad-lin', 'quad-int',
				'lin-lin', 'lin-int']

	stat_tracking = pd.DataFrame(columns=columns)

	for region in regions:


		if region == 'postcentral_rh':
			color = 'navajowhite'
		if region == 'precentral_rh':
			color = 'navajowhite'
		if region == 'precentral_lh':
			color = 'lightsteelblue'
		if region == 'postcentral_lh':
			color = 'lightsteelblue'

		df_component_og_copy = df_component_og.copy()
		df_component = df_component_og_copy[df_component_og_copy['region'] == region]


		for parameter in parameters:


			print('region')
			print(region)

			print('parameter:')
			print(parameter)

			w=[]
			# load in variances
			for subjectID in subjectIDs:
				if os.path.exists('/media/NAS/bbrady/Burst_Triggered_Waveforms/anatomical_ROI_results/' + subjectID + '/' + region +'FOOOF_bootstrap_parameters.csv'):
					df_var = pd.read_csv('/media/NAS/bbrady/Burst_Triggered_Waveforms/anatomical_ROI_results/' + subjectID + '/' + region +'FOOOF_bootstrap_parameters.csv')	
					var = df_var[parameter].var()
					std = np.sqrt(var)
					weight = 1/std
					df_weight = pd.DataFrame({'subjectID':subjectID, 'weight':[weight], 'variance':[var]})
					w.append(df_weight)
				df_weights = pd.concat(w)
				weights = dict(zip(df_weights['subjectID'], df_weights['weight']))
				variances = dict(zip(df_weights['subjectID'], df_weights['variance']))
			df_component['weight'] = df_component['subjectID'].map(weights)
			df_component['variance'] = df_component['subjectID'].map(variances)

			# get fit coefficients
			df_linear_fit = calc_l_fit_coefficients(df_component['Age'].to_numpy(), df_component[parameter].to_numpy(), parameter, df_component['weight'].to_numpy())
			df_quadratic_fit = calc_q_fit_coefficients(df_component['Age'].to_numpy(), df_component[parameter].to_numpy(), parameter, df_component['weight'].to_numpy())
			
			# F - test for model comparison
			chi_sq_l = calc_chi_square_l(parameter, df_component, df_linear_fit, region, df_component)
			chi_sq_q = calc_chi_square_q(parameter, df_component, df_quadratic_fit, region, df_component)



			#########################################################################
			# Do Regression Statistics
			##########################################################################

			x = df_component['Age']
			y = df_component[parameter]
			df=pd.DataFrame(columns=['y','x'])
			df['x'] = x
			df['y'] = y


			#calc F-statistic and d.o.f
			F, f1, f2 = calc_F_stat(N, chi_sq_l, chi_sq_q)
			print('F-stat:')
			print(F)
			# test the null hypothesis
			if F>6.635000:
				appropriate_model = 'quadratic'
			if F<6.635000:
				appropriate_model = 'linear'
			print('appropriate model:')
			print(appropriate_model)


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
			3.842 - 0.95 - 0.05 ----- THIS ONE
			6.635 - 0.99 - 0.01
			>6.635 - >0.99 - <0.01

			high F-stat means the additional model component is important
			1-CL = probability that fit improvement is due by chance

			examples:
			ex. if F = 6.888, very very likely there is a real differnce in the models! There is <0.01 probability that difference is just by chance
			ex. if F = 2.888, likely there is a real differnce in the models. There is <0.10 but >0.05 probability that difference is just by chance
			ex. if F = 0.565, who knows if there is a difference. There is <0.50 but >0.40 probability that difference is just by chance

			'''






			#model evaluation to retireve p-value
			degree=1
			coefs = np.polyfit(df_component['Age'].to_numpy(), df_component[parameter].to_numpy(), degree, w=df_component['weight'].to_numpy())
			df=pd.DataFrame(columns=['y','x'])
			df['x'] = df_component['Age'].to_numpy()
			df['y'] = df_component[parameter].to_numpy()
			model = np.poly1d(coefs) 
			results = smf.ols(formula ='y ~ model(x)', data=df).fit()
			p_val_lin = results.pvalues.loc['model(x)']
			print('p-val lin:')
			print(p_val_lin)


			degree=2
			coefs = np.polyfit(df_component['Age'].to_numpy(), df_component[parameter].to_numpy(), degree, w=df_component['weight'].to_numpy())
			df=pd.DataFrame(columns=['y','x'])
			df['x'] = df_component['Age'].to_numpy()
			df['y'] = df_component[parameter].to_numpy()
			model = np.poly1d(coefs) 
			results = smf.ols(formula ='y ~ model(x)', data=df).fit()
			p_val_quad = results.pvalues.loc['model(x)']
			print('p-val quad:')
			print(p_val_quad)
			print('-------------------------------------')

			stat_tracking = stat_tracking.append({'region':region, 
													'parameter':parameter, 
													'chi_sq_l':chi_sq_l, 
													'chi_sq_q':chi_sq_q, 
													'F-stat':F, 
													'appropriate model':appropriate_model, 
													'p_val_lin':p_val_lin,
													'p_val_quad':p_val_quad,
													'quad-quad':df_quadratic_fit['quadratic'][0], 
													'quad-lin':df_quadratic_fit['linear'][0], 
													'quad-int':df_quadratic_fit['intercept'][0],
													'lin-lin':df_linear_fit['linear'][0],
													'lin-int':df_linear_fit['intercept'][0]}, ignore_index=True)


		'''
		subject_avgs = []

		#make average per participant for plotting purposes
		for subjectID in subjectIDs:
			df_copy = df_component.copy()
			df_subj = df_copy[df_copy['subjectID'] == subjectID]
			df_subj.loc['mean'] = df_subj.mean()
			df_subj = df_subj[-1:]
			df_subj['subjectID'] = subjectID
			subject_avgs.append(df_subj)

		df_avgs =pd.concat(subject_avgs) 
		'''

	#stat_tracking.to_csv(outDir + 'FOOOF_ageing_stats.csv')

		

	################################
	# Generate Plots
	################################

	
	
	sns.set(style='white', palette='deep' ,font_scale=1.5, rc={'figure.figsize':(2,5)})
	sns.set_style('ticks', {'xtick.major.size':8})


	# generate quadratic fit plot - exponent
	sns.regplot(data=df_aperiodic[df_aperiodic['region']==region], x='Age', y='exponent', ci=95, fit_reg=True, scatter=True, order=2, scatter_kws={'s':15, 'color':color}, line_kws={'linewidth':3.5, 'color':'dimgray'})#'limegreen'}) 
	axes = plt.gca()
	sns.despine()
	axes.set_ylim([0.45,1.55])
	axes.set_xlim([8,98])
	axes.set_xlabel('Age (years)')
	axes.set_ylabel('')
	axes.set_xticks([20,50,80], minor=False)
	#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
	axes.set_yticks([])
	axes.tick_params(direction='out')
	plt.tight_layout()
	plt.show()
	#plt.savefig(outDir + region + '_exponent_vs_age_quadratic.png', bbox_inches='tight')
	#plt.close()




	
	# generate quadratic fit plot - offset
	sns.regplot(data=df_aperiodic[df_aperiodic['region']==region], x='Age', y='offset', ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':15, 'color':color}, line_kws={'linewidth':3.5, 'color':'dimgray'})#'limegreen'})  
	axes = plt.gca()
	axes.set_ylim([-3.70,0.35])
	axes.set_xlim([8,98])
	axes.set_xlabel('Age (years)')
	axes.set_ylabel('')
	axes.set_xticks([20,50,80], minor=False)
	#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
	axes.set_yticks([])
	axes.tick_params(direction='out')
	sns.despine()
	plt.tight_layout()
	plt.show()
	#plt.savefig(outDir + region + '_offset_vs_age_quadratic.png', bbox_inches='tight')
	#plt.close()
	

	'''
		# generatelinear fit plot - offset
		sns.regplot(data=df_aperiodic[df_aperiodic['region']==region], x='Age', y='offset', ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':15, 'color':color}, line_kws={'linewidth':3.5, 'color':'dimgray'})#'limegreen'}) 
		axes = plt.gca()
		axes.set_ylim([-3.70,0.35])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('offset')
		axes.set_xticks([20,50,80], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='out')
		plt.savefig(outDir + region + '_offset_vs_age_linear.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_aperiodic[df_aperiodic['region']==region], x='Age', y='offset', ci=95, fit_reg=False, scatter=True, order=1, scatter_kws={'s':15, 'color':color}, line_kws={'linewidth':3.5, 'color':'dimgray'})#'limegreen'}) 
		axes = plt.gca()
		axes.set_ylim([-3.70,0.35])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('offset')
		axes.set_xticks([20,50,80], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='out')
		plt.savefig(outDir + region + '_offset_vs_age_no_fit.png', bbox_inches='tight')
		plt.close()





		# generate quadratic fit plot - offset
		sns.regplot(data=df_fit[df_fit['region']==region], x='Age', y='r_squared', ci=95, fit_reg=True, scatter=True, order=2, scatter_kws={'s':15, 'color':color}, line_kws={'linewidth':3.5, 'color':'dimgray'})#'limegreen'}) 
		axes = plt.gca()
		axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('r_squared')
		axes.set_xticks([20,50,80], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='out')
		plt.savefig(outDir + region + '_r_squared_vs_age_quadratic.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_fit[df_fit['region']==region], x='Age', y='r_squared', ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':15, 'color':color}, line_kws={'linewidth':3.5, 'color':'dimgray'})#'limegreen'}) 
		axes = plt.gca()
		axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('r_squared')
		axes.set_xticks([20,50,80], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='out')
		plt.savefig(outDir + region + '_r_squared_vs_age_linear.png', bbox_inches='tight')
		plt.close()


		# generatelinear fit plot - offset
		sns.regplot(data=df_fit[df_fit['region']==region], x='Age', y='r_squared', ci=95, fit_reg=False, scatter=True, order=1, scatter_kws={'s':15, 'color':color}, line_kws={'linewidth':3.5, 'color':'dimgray'})#'limegreen'}) 
		axes = plt.gca()
		axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('r_squared')
		axes.set_xticks([20,50,80], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='out')
		plt.savefig(outDir + region + '_r_squared_vs_age_no_fit.png', bbox_inches='tight')
		plt.close()



		# generate quadratic fit plot - offset
		sns.regplot(data=df_fit[df_fit['region']==region], x='Age', y='error', ci=95, fit_reg=True, scatter=True, order=2, scatter_kws={'s':15, 'color':color}, line_kws={'linewidth':3.5, 'color':'dimgray'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('error')
		axes.set_xticks([20,50,80], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='out')
		plt.savefig(outDir + region + '_error_vs_age_quadratic.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_fit[df_fit['region']==region], x='Age', y='error', ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':15, 'color':color}, line_kws={'linewidth':3.5, 'color':'dimgray'})
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('error')
		axes.set_xticks([20,50,80], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='out')
		plt.savefig(outDir + region + '_error_vs_age_linear.png', bbox_inches='tight')
		plt.close()

	'''

	











	'''
		# generate quadratic fit plot - offset
		sns.regplot(data=df_avgs, x='Age', y='CF', ci=95, fit_reg=True, scatter=True, order=2, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'orangered'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('center frequency (Hz)')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'beta_CF_vs_age_quadratic.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_avgs, x='Age', y='CF', ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('center frequency (Hz)')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'beta_CF_vs_age_linear.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_avgs, x='Age', y='CF', ci=95, fit_reg=False, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('center frequency (Hz)')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'beta_CF_vs_age_no_fit.png', bbox_inches='tight')
		plt.close()

		# generate quadratic fit plot - offset
		sns.regplot(data=df_avgs, x='Age', y='PW_scaled', ci=95, fit_reg=True, scatter=True, order=2, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'orangered'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('peak power')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'beta_PW_scaled_vs_age_quadratic.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_avgs, x='Age', y='PW_scaled', ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('peak power')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'beta_PW_scaled_vs_age_linear.png', bbox_inches='tight')
		plt.close()

		# generate quadratic fit plot - offset
		sns.regplot(data=df_avgs, x='Age', y='BW', ci=95, fit_reg=True, scatter=True, order=2, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'orangered'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('peak power')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'beta_BW_vs_age_quadratic.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_avgs, x='Age', y='BW', ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('peak power')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'beta_BW_vs_age_linear.png', bbox_inches='tight')
		plt.close()





		











		# generate quadratic fit plot - offset
		sns.regplot(data=df_high_beta_peaks, x='Age', y='CF', ci=95, fit_reg=True, scatter=True, order=2, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'orangered'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('center frequency (Hz)')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'high_beta_CF_vs_age_quadratic.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_high_beta_peaks, x='Age', y='CF', ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('center frequency (Hz)')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'high_beta_CF_vs_age_linear.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_high_beta_peaks, x='Age', y='CF', ci=95, fit_reg=False, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('center frequency (Hz)')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'high_beta_CF_vs_age_no_fit.png', bbox_inches='tight')
		plt.close()

		# generate quadratic fit plot - offset
		sns.regplot(data=df_high_beta_peaks, x='Age', y='PW', ci=95, fit_reg=True, scatter=True, order=2, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'orangered'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('peak power')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'high_beta_PW_vs_age_quadratic.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_high_beta_peaks, x='Age', y='PW', ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('peak power')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'high_beta_PW_vs_age_linear.png', bbox_inches='tight')
		plt.close()































		


		# generate quadratic fit plot - offset
		sns.regplot(data=df_mu_peaks, x='Age', y='CF', ci=95, fit_reg=True, scatter=True, order=2, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'orangered'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('center frequency (Hz)')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'mu_CF_vs_age_quadratic.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_mu_peaks, x='Age', y='CF', ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('center frequency (Hz)')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'mu_CF_vs_age_linear.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_mu_peaks, x='Age', y='CF', ci=95, fit_reg=False, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('center frequency (Hz)')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'mu_CF_vs_age_no_fit.png', bbox_inches='tight')
		plt.close()

		# generate quadratic fit plot - offset
		sns.regplot(data=df_mu_peaks, x='Age', y='PW', ci=95, fit_reg=True, scatter=True, order=2, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'orangered'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('peak power')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'mu_PW_vs_age_quadratic.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_mu_peaks, x='Age', y='PW', ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('peak power')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'mu_PW_vs_age_linear.png', bbox_inches='tight')
		plt.close()



		# generatelinear fit plot - offset
		sns.regplot(data=df_mu_peaks, x='Age', y='PW', ci=95, fit_reg=False, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('peak power')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'mu_PW_vs_age_no_fit.png', bbox_inches='tight')
		plt.close()

		











		
		# generate quadratic fit plot - offset
		sns.regplot(data=df_theta_peaks, x='Age', y='CF', ci=95, fit_reg=True, scatter=True, order=2, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'orangered'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('center frequency (Hz)')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'theta_CF_vs_age_quadratic.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_theta_peaks, x='Age', y='CF', ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('center frequency (Hz)')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'theta_CF_vs_age_linear.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_theta_peaks, x='Age', y='CF', ci=95, fit_reg=False, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('center frequency (Hz)')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'theta_CF_vs_age_no_fit.png', bbox_inches='tight')
		plt.close()

		# generate quadratic fit plot - offset
		sns.regplot(data=df_theta_peaks, x='Age', y='PW', ci=95, fit_reg=True, scatter=True, order=2, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'orangered'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('peak power')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'theta_PW_vs_age_quadratic.png', bbox_inches='tight')
		plt.close()

		# generatelinear fit plot - offset
		sns.regplot(data=df_theta_peaks, x='Age', y='PW', ci=95, fit_reg=True, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('peak power')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'theta_PW_vs_age_linear.png', bbox_inches='tight')
		plt.close()



		# generatelinear fit plot - offset
		sns.regplot(data=df_theta_peaks, x='Age', y='PW', ci=95, fit_reg=False, scatter=True, order=1, scatter_kws={'s':10}, line_kws={'linewidth':3.5, 'color':'limegreen'})#'limegreen'}) 
		axes = plt.gca()
		#axes.set_ylim([0.96,1.01])
		axes.set_xlim([8,98])
		axes.set_xlabel('subject age (y)')
		axes.set_ylabel('peak power')
		axes.set_xticks([15,30,45,60,75,90], minor=False)
		#axes.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
		axes.tick_params(direction='in')
		plt.savefig(outDir + 'theta_PW_vs_age_no_fit.png', bbox_inches='tight')
		plt.close()
		




	'''
























