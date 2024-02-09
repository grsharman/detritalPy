###############################################################
# Import required modules
###############################################################

from detritalpy import detritalFuncs as dFunc # Working with detrital geochron data
import numpy as np # For efficient manipulation of arrays of numbers
from matplotlib import pyplot as plt # For making plots
import pandas as pd # For working with pandas dataframes/ loading data
import time # For keeping track of model progress
from operator import itemgetter # TODO: Consider removing this
from sklearn.decomposition import NMF
import xlsxwriter

###############################################################
# Functions for nmf 
###############################################################

def nmf(ages, errors=None, dist_type='KDE', bw=20, bw_x=None, x1=0, x2=4500, nEMs=10, max_iter=10000, tol=1e-8, verbose=True):
	"""
	Code for conducting non-negative matrix factorization (NMF) using the scikit-learn library.

    Parameters
    ----------
    ages : list of arrays
    	list of arrays of ages for each sample or sample group. Output from sampleToData()
    errors : list of arrays, default: None
    	list of arrays of errors (1-sigma) for each sample or sample group. Output from sampleToData()
    dist_type : {'PDP','KDE'}, default: 'KDE'
    	type of distribution to model. Options: 'KDE' or 'PDP'
    bw : {'ISJ','optimizedFixed','optimizedVariable'}, int, or list of int, default: 20
    	type of bandwidth to use if dist_type is 'KDE'. Options are 'ISJ','optimizedFixed', 'optimizedVariable', or a number (bandwidth in Myr). Placing numbers in a list will result in use of a split KDE
    bw_x : None or list of int, default: None
    	list of x-axis split locations if multiple bw values are specified. Use None if not using a split bandwdith
    x1 : int or float, default: 0
    	beginning of range to compute the age distribution used for modeling (default = 0 Ma)
    x2 : int or float, default: 4500
    	beginning of range to compute the age distribution used for modeling (default = 4500 Ma)
	nEMs : int, default: 10
		number of end-members to model
	max_iter : int, default: 10000
		maximum number of iterations before the NMF algorithm times out. See sklearn.decomposition.NMF for more information.
	tol : float, default: 1e-8
		tolerance for the stopping condition of the NMF algorithm. See sklearn.decomposition.NMF for more information.
	verbose : bool, default: True
		Set to True to print a readout of progress

	Returns
	-------
	xAge : array
		array with x-axis values that correspond to end-member age distributions (X)
	X : array
		array with input sample age distributions (KDE or PDP) with shape (# of samples, # of x-axis values)
	Xmodeled : list of arrays
		list of arrays with modeled sample age distributions (i.e., W*H) with length nEMs. Each array has shape (# of samples, # of x-axis values)
	W : list of arrays
		list of arrays with raw sample weightings with length nEMs. Each array has shape (# of samples, # of end-members)
	Wnorm : list of arrays
		list of arrays with normalized sample weightings with length nEMs. Each array has shape (# of samples, # of end-members)
	Hnorm : list of arrays
		list of arrays with end-member age distributions with length nEMs. Each array has shape (# of end-members, # of x-axis values)
	modelRecstErr : list of float
		list of model reconstruction error (float). See sklearn.decomposition.NMF for more information.
	nint : list of int
		list of # of iterations for each model. See sklearn.decomposition.NMF for more information.
	in_r2_avg : float
		average r2 of input samples (excluding self-comparisons)
	in_r2_range : float
		range of r2 of input samples (excluding self-comparisons)
	in_Vmax_avg: float
		average Kuiper Vmax of input samples (excluding self-comparisons)
	in_Vmax_range : float
		range of Kuiper Vmax of input samples (excluding self-comparisons)
	out_r2 : list of arrays
		list of arrays with r2 values of comparisons of input samples with modeled samples with length nEMs. Each array has shape (# of samples,)
	out_Vmax : list of ararys
		list of arrays with Kuiper Vmax values of comparisons of input samples with modeled samples with length nEMs. Each array has shape (# of samples,)
	"""
	start = time.time() # Record the starting time of the process

	if not(dist_type in ['KDE','PDP']):
		raise Warning('dist_type must be either KDE or PDP')

	if dist_type == 'KDE':
		if bw == None:
			raise Warning('bw must not be None')
		xAge, X = dFunc.KDEcalcAges(ages, x1=x1, x2=x2, xdif=1, bw=bw, bw_x=bw_x, cumulative=False) # Use this for the KDE
	if dist_type == 'PDP':
		if errors == None:
			raise Warning('errors must not be None')
		xAge, X = dFunc.PDPcalcAges(ages, errors, x1=x1, x2=x2, xdif=1, cumulative=False) # Use this for the PDP

	# Create empty variables for storing results in
	W = []
	Wnorm = []
	H = []
	modelRecstErr = []
	nint = []

	# Execute NMF
	for i in range(nEMs): # One loop for each set of end-member calculations
		if verbose:
			print('------------------------')
			print('Starting end-member', i+1)
		start_i = time.time()
		model = NMF(n_components=i+1, solver='cd', max_iter=100000, init='nndsvda', tol=1e-12)#, init='random')#, random_state=0) # Note: can adjust max_iter and tol search parameters
		W.append(model.fit_transform(X))
		H.append(model.components_)
		modelRecstErr.append(model.reconstruction_err_)
		nint.append(model.n_iter_)
		if verbose:
			print('Finished! Time elapsed (s):',round(time.time()-start_i, 3))
	if verbose:
		print('Finished NMF....')

	# Normalize W and H so that they sums to one
	Wnorm = []
	for i in range(len(W)):
		WnormTemp = np.copy(W[i])
		for j in range(len(W[i])):
			WnormTemp[j] = W[i][j]/sum(W[i][j])
		Wnorm.append(WnormTemp)

	Hnorm = H.copy()
	for i in range(len(H)): # One loop for each end-member scenario
		for j in range(len(H[i])): # One loop for each end-member
			Hnorm[i][j] = H[i][j]/np.sum(H[i][j]) # Normalize so area under curve = 1

	if verbose:
		print('Finished normalization....')

	# Preallocate arrays to store comparisons in
	r2 = np.zeros(shape=(len(X),len(X)))
	Dmax = np.zeros(shape=(len(X),len(X)))
	Vmax = np.zeros(shape=(len(X),len(X)))

	for i in range(len(X)):
		for j in range(len(X)):
			r2[i,j] = dFunc.calcR2(X[i],X[j])
			Dmax[i,j] = dFunc.calcDmax(np.cumsum(X[i]),np.cumsum(X[j]))
			Vmax[i,j] = dFunc.calcVmax(np.cumsum(X[i]),np.cumsum(X[j]))

	d = np.diag_indices(len(X))

	r2_mask = np.ma.array(r2, mask=False)
	r2_mask.mask[d] = True
	in_r2_avg = np.average(r2_mask)
	in_r2_range = np.max(r2_mask)-np.min(r2_mask)

	Vmax_mask = np.ma.array(Vmax, mask=False)
	Vmax_mask.mask[d] = True
	in_Vmax_avg = np.average(Vmax_mask)
	in_Vmax_range = np.max(Vmax_mask)-np.min(Vmax_mask)

	Xmodeled, out_r2, out_Vmax = compute_child_fit(H=Hnorm, X=X, W=Wnorm)

	if verbose:
		print('------------------------')
		print('Finished! Total time elapsed (s):',round(time.time()-start, 3))
		print('Input sample r2 mean and range: ',round(in_r2_avg, 3),',', round(in_r2_range,3))
		print('Input sample Vmax mean and range: ',round(in_Vmax_avg, 3),',', round(in_Vmax_range,3))

	return xAge, X, Xmodeled, W, Wnorm, Hnorm, modelRecstErr, nint, in_r2_avg, in_r2_range, in_Vmax_avg, in_Vmax_range, out_r2, out_Vmax

def compute_child_fit(H, X, W):
	"""
	Function that computes goodness-of-fit between input samples and modeled samples (i.e., W*H)
	Outputs Kuiper Vmax and Cross-correlation (r2) coefficients

	Note: array components of H, X, and W should all be normalized such that they sum-to-one

    Parameters
    ----------

	H : list of arrays
		list of arrays with normalized end-member age distributions
	X : list of arrays
		list of arrays with input sample age distributions
	W : list of arrays
		list of arrays with normalized sample weightings

	Returns
	-------
	Xmodeled : list of arrays
		list of arrays with modeled sample age distributions (i.e., W*H) with length nEMs. Each array has shape (# of samples, # of x-axis values)
	out_r2 : list of arrays
		list of arrays with r2 values of comparisons of input samples with modeled samples with length nEMs. Each array has shape (# of samples,)
	out_Vmax : list of ararys
		list of arrays with Kuiper Vmax values of comparisons of input samples with modeled samples with length nEMs. Each array has shape (# of samples,)
	"""

	# First compute the relative distributions for each modeled daughter
	Xmodeled = []
	for i in range(len(H)): # One loop for each set of end-member scenarios
		XmodeledTemp= []
		for j in range(len(X)): # One loop for each observation (daughter sample)
			daughterModeled = np.zeros(shape=(H[0].shape[1],))
			for k in range(len(W[i][0])): # One loop for each end-member abundance
				toSum = [x * W[i][j][k] for x in H[i][k]]
				daughterModeled = [x + y for x, y in zip(daughterModeled, toSum)]
			XmodeledTemp.append(daughterModeled)
		Xmodeled.append(np.asarray(XmodeledTemp))

	# Compute similarity metrics between daughter and modeled daughter
	r2 = []
	Vmax = []
	for i in range(len(H)): # One loop for each set of end-member scenarios
		r2Temp = []
		VmaxTemp = []
		for j in range(len(X)): # One loop for each observation (daughter sample)
			r2Temp.append(dFunc.calcR2(Xmodeled[i][j],X[j]))
			Xmodeled_cumsum = np.cumsum(Xmodeled[i][j]/np.sum(Xmodeled[i][j]))
			X_cumsum = np.cumsum(X[j]/np.sum(X[j]))
			VmaxTemp.append(dFunc.calcVmax(Xmodeled_cumsum,X_cumsum))
		r2.append(np.asarray(r2Temp))
		Vmax.append(np.asarray(VmaxTemp))

	return Xmodeled, r2, Vmax


def nmf_to_excel(labels, numGrains, xAge, X, Wnorm, 
	H, modelRecstErr, nEMs, EMs_min, EMs_max, max_iter, nint, tol,
	in_r2_avg, in_r2_range, in_Vmax_avg, in_Vmax_range, out_r2, out_Vmax,
	dist_type, bw, bw_x=None, file_name='dPy-nmf_results.xlsx', version=None, verbose=True):
	"""
	Code for outputing NMF results to an Excel file

	Parameters
	----------

    labels : list
    	list of labels for each sample or sample group. Output from dFunc.sampleToData()
    numGrains : list
    	list of the number of analyses for each sample or sample group. Output from dFunc.sampleToData()
	xAge : array
		array with x-axis values that correspond to end-member age distributions (X). Output from nmf()
	X : array
		array with input sample age distributions (KDE or PDP) with shape (# of samples, # of x-axis values)
	Wnorm : list of arrays
		list of arrays with normalized sample weightings with length nEMs. Each array has shape (# of samples, # of end-members)
	H : list of arrays
		list of arrays with normalized end-member age distributions with length nEMs. Each array has shape (# of end-members, # of x-axis values)
	modelRecstErr : list of float
		list of model reconstruction error (float). See sklearn.decomposition.NMF for more information.
	nEMs : int, default: 10
		number of end-members to model
	EMs_min : int
		minimum number of end-members to report results for
	EMs_max : int
		maximum number of end-members to report results for
	max_iter : int
		maximum number of iterations before the NMF algorithm times out. See sklearn.decomposition.NMF for more information.
	nint : list of int
		list of # of iterations for each model. See sklearn.decomposition.NMF for more information.
	tol : float
		tolerance for the stopping condition of the NMF algorithm. See sklearn.decomposition.NMF for more information.
	in_r2_avg : float
		average r2 of input samples (excluding self-comparisons)
	in_r2_range : float
		range of r2 of input samples (excluding self-comparisons)
	in_Vmax_avg: float
		average Kuiper Vmax of input samples (excluding self-comparisons)
	in_Vmax_range : float
		range of Kuiper Vmax of input samples (excluding self-comparisons)
	out_r2 : list of arrays
		list of arrays with r2 values of comparisons of input samples with modeled samples with length nEMs. Each array has shape (# of samples,)
	out_Vmax : list of ararys
		list of arrays with Kuiper Vmax values of comparisons of input samples with modeled samples with length nEMs. Each array has shape (# of samples,)
    dist_type : {'PDP','KDE'}, default: 'KDE'
    	type of distribution to model. Options: 'KDE' or 'PDP'
    bw : {'ISJ','optimizedFixed','optimizedVariable'}, int, or list of int
    	type of bandwidth to use if dist_type is 'KDE'. Options are 'ISJ','optimizedFixed', 'optimizedVariable', or a number (bandwidth in Myr). Placing numbers in a list will result in use of a split KDE
    bw_x : None or list of int
    	list of x-axis split locations if multiple bw values are specified. Use None if not using a split bandwdith
    file_name : str, default: 'dPy-nmf_results.xlsx'
    	output file name to use. Should end with '.xlsx'
    version : str, default: None
    	version of detritalPy used in modeling. set to detritalpy.__version__
	verbose : bool, default: True
		set to True to print when export is completed

	Returns
	-------
	Excel workbook with NMF results
	"""
	workbook = xlsxwriter.Workbook(file_name)

	merge_format = workbook.add_format({'align': 'center'})
	cell_format = workbook.add_format({'align' : 'center'})
	percent_format = workbook.add_format({'num_format': '0.0%'})
	bold_format = workbook.add_format({'bold' : True})

	# Record model parameters
	worksheet = workbook.add_worksheet('Model_parameters')
	worksheet.write(0, 0, 'detritalPy version '+str(version))

	worksheet.write(2, 0, 'Input sample paramters', bold_format)
	worksheet.write(3, 0, 'Conducting NMF on the following samples: '+str(labels))
	worksheet.write(4, 0, 'Cross-correlation (r2) mean of samples: '+str(round(in_r2_avg,3)))
	worksheet.write(5, 0, 'Cross-correlation (r2) range of samples: '+str(round(in_r2_range,3)))
	worksheet.write(6, 0, 'Kuiper Vmax mean of samples: '+str(round(in_Vmax_avg,3)))
	worksheet.write(7, 0, 'Kuiper Vmax range of samples: '+str(round(in_Vmax_range,3)))

	worksheet.write(9, 0, 'Model parameters', bold_format)
	worksheet.write(10, 0, 'Number of samples: '+str(len(labels)))
	worksheet.write(11, 0, 'Total number of analyses: '+str(np.sum(numGrains)))
	worksheet.write(12, 0, 'Avg number of analyses per sample: '+str(round(np.average(numGrains),3)))
	worksheet.write(13, 0, 'Type of distribution used in modeling: '+dist_type)
	worksheet.write(14, 0, str('KDE bandwidth: '+str(bw)+' (only relevant if using a KDE for modeling)'))
	worksheet.write(15, 0, 'KDE bandwidth x-axis split: '+str(bw_x)+' (only relevant if using a KDE for modeling)')
	worksheet.write(16, 0, 'Maximum number of end-members modeled: '+str(nEMs))
	worksheet.write(17, 0, 'Selected number of end-members: '+str(EMs_max))
	worksheet.write(18, 0, 'Maximum number of iterations: '+str(max_iter))
	worksheet.write(19, 0, 'Tolerance: '+str(tol))

	worksheet = workbook.add_worksheet('Model_fit')
	worksheet.write(0, 0, 'Number of end-members')
	worksheet.write(0, 1, 'Final residual error')
	worksheet.write(0, 2, 'Number of iterations')
	worksheet.write(0, 3, 'r2 (avg)')
	worksheet.write(0, 4, 'r2 (st. dev.)')
	worksheet.write(0, 5, 'Vmax (avg)')
	worksheet.write(0, 6, 'Vmax (st. dev.)')

	for i in range(nEMs):
		worksheet.write(i+1, 0, i+1)
		worksheet.write(i+1, 1, modelRecstErr[i])
		worksheet.write(i+1, 2, nint[i])
		worksheet.write(i+1, 3, round(np.average(out_r2[i]),3))
		worksheet.write(i+1, 4, round(np.std(out_r2[i]),3))
		worksheet.write(i+1, 5, round(np.average(out_Vmax[i]),3))
		worksheet.write(i+1, 6, round(np.std(out_Vmax[i]),3))

	worksheet.write(nEMs+2, 0, 'r2 and Vmax values are comparisons of modeled vs actual sample age distributions')

	# Write worksheet with end-member age distributions
	worksheet = workbook.add_worksheet('EMs')

	worksheet.write(0, 0, 'Age (Ma)')
	for i in range(len(xAge)):
		worksheet.write(i+1, 0, xAge[i])
	c = 0 # Counter variable
	for i in range(EMs_min-1, EMs_max): # One loop for each end-member scenario
		for j in range(len(H[i])): # One loop for each end-member
			EMscenarioLabel = str(len(H[i])) + 'EM_EM' + str(j+1)
			worksheet.write(0, c+1, EMscenarioLabel)
			for k in range(len(H[i][j])):
				worksheet.write(k+1, c+1, H[i][j][k])
			c += 1

	# Write worksheet with end-member abundances
	worksheet = workbook.add_worksheet('Abundances')

	worksheet.write(0, 0, 'Sample')

	#EMlabels = []
	c = 0 # Counter variable
	for i in range(EMs_min-1, EMs_max): # One loop for each end-member scenario
		for j in range(len(H[i])): # One loop for each end-member
			EMscenarioLabel = str(len(H[i])) + 'EM_EM' + str(j+1)
			#EMlabels.append(EMscenarioLabel)
			worksheet.write(0, c+1, EMscenarioLabel)
			c += 1

	for i in range(len(labels)):
		worksheet.write(i+1, 0, labels[i])


	for i in range(len(labels)): # One loop per sample
		c = 0
		for j in range(EMs_min-1, EMs_max): # One loop for each end-member scenario
			for k in range(len(H[j])): # One loop for each EM
				worksheet.write(i+1, c+1, Wnorm[j][i][k])
				c += 1

	# Write worksheet with sample fit parameters
	worksheet = workbook.add_worksheet('Sample_fit')

	worksheet.write(0, 0, 'Sample')

	for i in range(EMs_min-1, EMs_max): # One loop for each end-member scenario
		R2_label = str(len(H[i])) + 'EM' + str('_R2')
		Vmax_label = str(len(H[i])) + 'EM' + str('_Vmax')
		worksheet.write(0, (2*i)+1, R2_label)
		worksheet.write(0, (2*i)+2, Vmax_label)

	for i in range(len(labels)): # One loop per sample
		worksheet.write(i+1, 0, labels[i])
		for j in range(EMs_min-1, EMs_max): # One loop for each end-member scenario
			worksheet.write(i+1, 2*j+1, out_r2[j][i])
			worksheet.write(i+1, 2*j+2, out_Vmax[j][i])

	workbook.close()

	if verbose:
		print('Workbook saved')

####################
# Plotting Functions
####################

def plot_reconstr_error(modelRecstErr, plot_width=8.0, plot_height=3.0):
	'''
	Function for plotting the model reconstruction error as a function of # of end-members modeled

	Parameters
	----------
	modelRecstErr : list of floats
		list of model reconstruction error (float). Output from nmf()
	plot_width : float
		plot width
	plot_height : float
		plot height

	Returns
	-------
	fig : matplotlib figure
	'''

	# Plot the reconstruction error as a function of the number of end-members
	fig, axs = plt.subplots(1,2, figsize=(plot_width,plot_height))
	x = np.arange(1,len(modelRecstErr),1)+1
	axs[0].plot(x, modelRecstErr[1:], 'ko-')
	axs[0].set_xlabel('Number of end-members')
	axs[0].set_ylabel('Reconstruction error')
	axs[0].set_xlim(2,len(modelRecstErr))

	sumSSR = SSR(x,modelRecstErr[1:])
	axs[1].plot(x, sumSSR, 'ko-')
	axs[1].set_xlabel('Number of end-members')
	axs[1].set_ylabel('SSR1 + SSR2')
	axs[1].set_xlim(2,len(modelRecstErr))

	# Plot minimum value
	minIndex = min(enumerate(sumSSR), key=itemgetter(1))[0] 
	axs[0].plot(x[minIndex], modelRecstErr[minIndex+1], 'o', color='red', markersize=10, label='Optimal')
	axs[1].plot(x[minIndex], sumSSR[minIndex], 'o', color='red', markersize=10)
	axs[0].legend(loc='best')

	y = modelRecstErr[1:]

	gxr_slope = np.polyfit(x[minIndex:],y[minIndex:],1)[0]
	gxr_yint = np.polyfit(x[minIndex:],y[minIndex:],1)[1]
	fxr_slope = np.polyfit(x[:(minIndex+1)],y[:(minIndex+1)],1)[0]
	fxr_yint = np.polyfit(x[:(minIndex+1)],y[:(minIndex+1)],1)[1]

	axs[0].plot((2, x[minIndex]), (fxr_slope*2+fxr_yint,fxr_slope*x[minIndex]+fxr_yint), '-', color='red', lw=0.5)
	axs[0].plot((x[minIndex], len(modelRecstErr)), (gxr_slope*x[minIndex]+gxr_yint,gxr_slope*len(modelRecstErr)+gxr_yint), '-', color='red', lw=0.5)

	plt.tight_layout()

	return fig


def plot_end_members(H, xAge, EMs_min, EMs_max, x_min=0, x_max=4000, plotLog=False, 
	plot_width_multiplier=1.0, plot_height_multiplier=1.0, c=2.0, w=3.0, agebins=None, agebinsc=None,
	EMcolors='Default', colorByAge=False, fillBetween=True):

	"""
	Function that plots end-member age distributions

	Parameters
	----------
	H : list of arrays
		list of arrays with normalized end-member age distributions with length nEMs. Each array has shape (# of end-members, # of x-axis values)
	xAge : array
		array with x-axis values that correspond to end-member age distributions (X)
	EMs_min : int
		the minimum number of end-members to include on the plot
	EMs_max : int
		the maximum number of end-members to include on the plot
	x_min : float or int, default: 0
		the minimum x-axis value in Ma
	x_max : float or int, default: 0
		the maximum x-axis value in Ma
	plotLog : bool, default: False
		Set to True to plot the x-axis in a log scale
	plot_width_multiplier : float, default: 1.0
		Parameter that scales the width of the plot
	plot_height_multiplier : float, default: 1.0
		Parameter that scales the height of the plot
	c : int, default: 2
		Height of the CDF plot in number of rows
	w : float, default: 3.0
		Width of individual subplots
	agebins : list, default: None
		array of bin edges in Myr. Format option 1: [age1, age2, age3, etc.]. Format option 2: [[bin1_min, bin1_max],[bin2_min, bin2_max],etc.]
    agebinsc : list, default: None
    	array of colors that correspond to age bins
	EMcolors : str or list, default: 'Default'
		'Default' or specify colors as a list or a single color as a string
	colorByAge : bool, default: False
		Set to True to color age distributions according to agebins and agebinsc
	fillBetween : bool, default: True
		Set to True to color age distributions according to EMcolors

	Returns
	-------
	fig : matplotlib figure
	"""


	if EMs_max > len(H):
		raise Warning('You may not plot more end-members than you have modeled')

	if type(EMs_max) != int:
		raise Warning('EMs_max must be an integer')

	if EMcolors == 'Default':
		EMcolors = [dFunc.colorMe(x) for x in np.arange(EMs_max)]

	if type(EMcolors) == str:
		EMcolors = [EMcolors] * int(EMs_max)

	EMlabels = ['EM'+str(x) for x in np.arange(EMs_min,EMs_max+1, 1)]

	fig = plt.figure(figsize=(plot_width_multiplier*(EMs_max-EMs_min+1)*w,plot_height_multiplier*c+(EMs_max)))

	mosaic = []
	for i in range(c): # Set up cumulative plots
	    mosaic_line = []
	    for j in range(EMs_max-EMs_min+1):
	        mosaic_line += ['0,{}'.format(j)]
	    mosaic.append(mosaic_line)

	for i in range(EMs_max):
	    mosaic_line = []
	    for j in range(EMs_max-EMs_min+1):
	        mosaic_line += ['{},{}'.format(i+1,j)]
	    mosaic.append(mosaic_line)

	axs = fig.subplot_mosaic(mosaic)
	fig.subplots_adjust(wspace=0)
	fig.subplots_adjust(hspace=0)

	for i in range(EMs_max-EMs_min+1):
		if i>0:
			axs['0,{}'.format(i)].set_yticklabels([])

	# First loop is for making CDF plots
	for ax_i, i in enumerate(range(EMs_min-1, EMs_max)): # One loop for each end-member scenario
		ax = axs['0,{}'.format(ax_i)]
		ax.set_xlim(x_min,x_max)
		ax.set_ylim(0,1.)
		for ax_j, j in enumerate(range(len(H[i]))): # One loop for each end-member in each scenario
			ax.plot(xAge, np.cumsum(H[i][j])/np.sum(H[i][j]), lw=2, color=EMcolors[ax_j], label = EMlabels[0:j+1][-1])
			if plotLog:
				ax.set_xscale('log')
			ax.get_xaxis().set_ticks([])
			if ax_i>0:
				ax.get_yaxis().set_ticks([])
		ax.set_title(str(len(H[i])) + 'EM')
		ax.legend(loc="lower right", prop={'size':12})
		if ax_i == 0:
			ax.set_ylabel('Cumulative')

	# Second loop is to set up relative plots
	for ax_i, i in enumerate(range(EMs_min-1, EMs_max)): # One loop for each end-member scenario
		for ax_j, j in enumerate(range(EMs_max)): # One loop for the max number of relative plots needed
			ax = axs['{},{}'.format(ax_j+1,ax_i)]
			if ax_j<=ax_i+EMs_min-1: # ax_j<=ax_i
				ax.plot(xAge, H[i][j], color='black', lw=1, label=('EM')+str(j+1))
				ax.legend(loc="upper right", prop={'size':12})
				if colorByAge:
					xdif = xAge[1]-xAge[0]
					if len(np.shape(agebins)) == 1:
						nage = len(agebins)-1
					if len(np.shape(agebins)) == 2:
						nage = len(agebins)
					for k in range(nage):
						if len(np.shape(agebins)) == 1:
							xage1 = agebins[k]
							xage2 = agebins[k+1]
						if len(np.shape(agebins)) == 2:
							xage1 = agebins[k][0]
							xage2 = agebins[k][1]
						xPart = np.arange(xage1, xage2+xdif, xdif)
						yPart = H[i][j][int(xage1/xdif):int((xage2+xdif)/xdif)]
						ax.fill_between(xPart, 0, yPart, alpha = 1, color=agebinsc[k])				
				if fillBetween:
					ax.fill_between(xAge, 0, H[i][j], color=EMcolors[ax_j])
				ax.set_xlim(x_min,x_max)
				ax.set_ylim(0)
				ax.get_yaxis().set_ticks([])
				if ax_j==ax_i+EMs_min-1:
					ax.set_xlabel('Age (Ma)')
				if ax_j<(ax_i+EMs_min-1):
					ax.get_xaxis().set_ticks([]) # Turn off x-ticks for all plots but the last one
				if plotLog:
					ax.set_xscale('log')

			else:
				ax.axis('off')

	return fig

def plot_child_goodness_of_fit(out_r2, out_Vmax, plot_width=15.0, plot_height=5.0):
	"""
	Function that makes a plot of r2 and Vmax that compares input samples with modeled
	samples (i.e., reconstructed samples, W*H)

	Parameters
	----------
	out_r2 : list of arrays
		list of arrays with r2 values of comparisons of input samples with modeled samples with length nEMs. Each array has shape (# of samples,)
	out_Vmax : list of ararys
		list of arrays with Kuiper Vmax values of comparisons of input samples with modeled samples with length nEMs. Each array has shape (# of samples,)
	plot_width : float, default: 8.0
		parameter that scales the total width of the plot
	plot_height : float, default: 5.0
		parameter that scales the total height of the plot

	Returns
	-------
	fig : matplotlib figure
	"""

	fig, axs = plt.subplots(1,2,figsize=(plot_width, plot_height))

	# Plot Cross-correlation (r2) results
	axs[0].plot(np.arange(1,len(out_r2)+1,1), out_r2, color='black')
	axs[0].plot(np.arange(1,len(out_r2)+1), np.average(out_r2, axis=1), color='red', lw=2, label='Avg')
	axs[0].set_title('r2')
	axs[0].set_xlabel('Number of end-members')
	axs[0].set_ylabel('Cross-correlation (r2)')
	axs[0].legend(fontsize='small')
	axs[1].plot(np.arange(1,len(out_r2)+1,1), out_Vmax, color='black')
	axs[1].plot(np.arange(1,len(out_r2)+1), np.average(out_Vmax, axis=1), color='red', lw=2, label='Avg')
	axs[1].set_title('Vmax')
	axs[1].set_xlabel('Number of end-members')
	axs[1].set_ylabel('Vmax')
	axs[1].legend(fontsize='small')

	plt.tight_layout()

	return fig

def plot_child_vs_modeled_CDF(labels, EMs_min, EMs_max, xAge, X, Xmodeled, out_r2, out_Vmax, 
	x_min=0, x_max=4000, fillBetween=False,
	plot_width_multiplier = 4.0, plot_height_multiplier=2.0, metric_to_plot='none'):

	"""
	Parameters
	----------
    labels : list
    	list of labels for each sample or sample group. Output from dFunc.sampleToData()
	EMs_min : int
		the minimum number of end-members to include on the plot
	EMs_max : int
		the maximum number of end-members to include on the plot
	xAge : array
		array with x-axis values that correspond to end-member age distributions (X)
	X : array
		array with input sample age distributions (KDE or PDP) with shape (# of samples, # of x-axis values)
	Xmodeled : list of arrays
		list of arrays with modeled sample age distributions (i.e., W*H) with length nEMs. Each array has shape (# of samples, # of x-axis values)
	out_r2 : list of arrays
		list of arrays with r2 values of comparisons of input samples with modeled samples with length nEMs. Each array has shape (# of samples,)
	out_Vmax : list of ararys
		list of arrays with Kuiper Vmax values of comparisons of input samples with modeled samples with length nEMs. Each array has shape (# of samples,)
	x_min : float or int, default: 0
		the minimum x-axis value in Ma
	x_max : float or int, default: 0
		the maximum x-axis value in Ma	
	fillBetween : bool, default: False
		Set to True to color between the CDFs
	plot_width_multiplier : float, default: 4.0
		Parameter that scales the width of the plot
	plot_height_multiplier : float, default: 2.0
		Parameter that scales the height of the plot
	metric_to_plot : {'r2', 'Vmax', 'none'}, default: 'none'
		Select 'r2' or 'Vmax' or 'none' to display the comparison metric between the two CDFs

	Returns
	-------
	fig : matplotlib figure
	"""

	# Make a plot that compares modeled with actual daughter
	fig, axs = plt.subplots(len(labels), EMs_max-EMs_min+1, figsize=(plot_width_multiplier*(EMs_max-EMs_min+1),plot_height_multiplier*len(labels)))

	c = 0 # Counter variable
	for i in range(len(labels)): # One loop for each sample
		for j_ax, j in enumerate(range(EMs_min-1, EMs_max)):
			ax = axs.flatten()[c]

			ax.plot(np.cumsum(Xmodeled[j][i])/np.sum(Xmodeled[j][i]),'--',label='modeled',color='black')
			ax.plot(np.cumsum(X[i]),label='observed',color='black')
			if fillBetween:
				ax.fill_between(x=xAge, y1=np.cumsum(Xmodeled[j][i])/np.sum(Xmodeled[j][i]), y2=np.cumsum(X[i]), color='gray')
			#ax[i,j].fill_between(x=x, y2=np.cumsum(Xmodeled[j][i])/np.sum(Xmodeled[j][i]), y1=np.cumsum(X[i]), color='darkblue')
			ax.set_xlim(x_min,x_max)
			ax.set_ylim(0,1.0)
			if metric_to_plot == 'Vmax':
				ax.text(x=(x_max-x_min)*0.97, y=0.07, ha='right', s='Vmax:'+str(np.round(out_Vmax[j][i],3)), size='small')
			if metric_to_plot == 'r2':
				ax.text(x=(x_max-x_min)*0.97, y=0.07, ha='right', s='R2:'+str(np.round(out_r2[j][i],3)), size='small')
			ax.get_yaxis().set_ticks([])
			if (i == 0) and (j_ax == 0):
				ax.legend(loc='upper left', fontsize='small')
			if i+1 < len(labels): # Only plot x-axis for last row
				ax.get_xaxis().set_ticks([])
			if i+1 == len(labels): # Plot x-axis label for last row
				ax.set_xlabel('Age (Ma)')
			if i == 0: # Plot titles for top row
				ax.set_title(str(j+1) + 'EM')
			if j_ax == 0: # Plot samples for left column
				ax.set_ylabel(labels[i], rotation=0, ha='right')
				#ax.text(x=(x_max-x_min)*-0.5, y=0.5, s=sampleList[i], size='x-small')
			c += 1 # Add 1 to counter variable

	fig.subplots_adjust(wspace=0)
	fig.subplots_adjust(hspace=0)

	return fig 

def plot_end_member_abundance(labels, Wnorm, EMcolors, EMs_min, EMs_max, plot_width_multiplier=2.0, plot_height_multiplier=0.25, bar_height=0.9):
	"""
	Function that plots end-member abundance as a horizontal bar plot
    
    Parameters
    ----------
    labels : list
    	list of labels for each sample or sample group. Output from dFunc.sampleToData()
	Wnorm : list of arrays
		list of arrays with normalized sample weightings with length nEMs. Each array has shape (# of samples, # of end-members)
	EMcolors : str or list, default: 'Default'
		'Default' or specify colors as a list or a single color as a string
	EMs_min : int
		the minimum number of end-members to include on the plot
	EMs_max : int
		the maximum number of end-members to include on the plot
	plot_width_multiplier : float, default: 2.0
		Parameter that scales the width of the plot
	plot_height_multiplier : float, default: 0.25
		Parameter that scales the height of the plot
	bar_height : float (0-1), default: 0.9
		Parameter that scales the relative height of the bars, such that a value of 1.0 indicates that they touch each othher

	Returns
	-------
	fig : matplotlib figure
	"""


	fig, axs = plt.subplots(1,EMs_max-EMs_min+1, figsize=(plot_width_multiplier*(EMs_max-EMs_min+1),plot_height_multiplier*len(labels)))

	y = np.arange(len(labels))
	y = np.flip(y)

	if EMcolors == 'Default':
		EMcolors = [dFunc.colorMe(x) for x in np.arange(EMs_max)]

	if type(EMcolors) == str:
		EMcolors = [EMcolors] * int(EMs_max)

	EMlabels = ['EM'+str(x) for x in np.arange(EMs_min,EMs_max+1, 1)]

	for h_ax, h in enumerate(range(EMs_min-1, EMs_max)): # One loop for each set of end-members
		if EMs_max == EMs_min:
			ax = axs
		else:
			ax = axs[h_ax]
		for i in range(len(Wnorm[h])): # One loop for each sample
			for j in range(len(Wnorm[h][0])): # One loop for each end-member
				#plt.barh(i, np.cumsum(WnormTest[i][j]), color=EMcolors[j])
				if j == 0:
					ax.barh(y[i], Wnorm[h][i][j], color=EMcolors[j], height=bar_height)
				else:
					ax.barh(y[i], Wnorm[h][i][j], left=np.cumsum(Wnorm[h][i])[j-1], color=EMcolors[j], height=bar_height)
		if h_ax == 0:
			ax.set_yticks(y)
			ax.set_yticklabels(labels, fontsize='x-small')
		else:
			ax.get_yaxis().set_ticks([])
		ax.set_xlim(0,1.)
		ax.xaxis.set_tick_params(labelsize='x-small')
		ax.set_title(EMlabels[h_ax], fontsize='medium')
		ax.set_ylim(-bar_height/2, y[0]+bar_height/2)

	return fig

def plot_end_member_abundance_as_pie(labels, Wnorm, EMcolors, EMs_min, EMs_max, plot_width_multiplier=2.0, plot_height_multiplier=0.25):
	"""
	Function that plots end-member abundance as pie plots
    
    Parameters
    ----------
    labels : list
    	list of labels for each sample or sample group. Output from dFunc.sampleToData()
	Wnorm : list of arrays
		list of arrays with normalized sample weightings with length nEMs. Each array has shape (# of samples, # of end-members)
	EMcolors : str or list, default: 'Default'
		'Default' or specify colors as a list or a single color as a string
	EMs_min : int
		the minimum number of end-members to include on the plot
	EMs_max : int
		the maximum number of end-members to include on the plot
	plot_width_multiplier : float, default: 2.0
		Parameter that scales the width of the plot
	plot_height_multiplier : float, default: 0.25
		Parameter that scales the height of the plot

	Returns
	-------
	fig : matplotlib figure
	"""

	fig, axs = plt.subplots(len(labels), EMs_max-EMs_min+1, figsize=(plot_width_multiplier*(EMs_max-EMs_min+1),plot_height_multiplier*len(labels)))

	y = np.arange(len(labels))
	y = np.flip(y)

	if EMcolors == 'Default':
		EMcolors = [dFunc.colorMe(x) for x in np.arange(EMs_max)]

	if type(EMcolors) == str:
		EMcolors = [EMcolors] * int(EMs_max)

	EMlabels = [str(x)+'EM' for x in np.arange(EMs_min,EMs_max+1, 1)]

	for h_ax, h in enumerate(range(EMs_min-1, EMs_max)): # One loop for each set of end-members
		for i in range(len(Wnorm[h])): # One loop for each sample
			axs[i,h_ax].pie(Wnorm[h][i], colors=EMcolors)
			if i == 0:
				axs[i, h_ax].set_title(EMlabels[h_ax])
			if h_ax == 0:
				axs[i, h_ax].set_ylabel(labels[i], rotation=0, fontsize='x-small', ha='right')

	return fig

##################
# Helper Functions
##################

def SSR(x, y):
	"""
	Function that calculates the sum of squared residuals of two linear regressions of reconstruction
	error vs number of end-members. For more details, please refer to Saylor et al. (2019): EPSL

	Parameters
	----------
	x : list or array
		end-members (exluding EM1 scenario), e.g., [2,3,4,5]
	y : list or array
		final model residuals

	Returns
	-------
	sumSSR : list of floats
		a list of the the sum of the squred residuals for both linear regressions, one value per end-member scenario

	"""

	# To supress rank warning because regression is poorly conditioned
	import warnings
	warnings.simplefilter('ignore', np.RankWarning)

	gxr_slope = []
	gxr_yint = []
	fxr_slope = []
	fxr_yint = []

	for i in range(len(x)):
		gxr_slope.append(np.polyfit(x[i:],y[i:],1)[0])
		gxr_yint.append(np.polyfit(x[i:],y[i:],1)[1])
		fxr_slope.append(np.polyfit(x[:(i+1)],y[:(i+1)],1)[0])
		fxr_yint.append(np.polyfit(x[:(i+1)],y[:(i+1)],1)[1])

	SSR2 = []
	for i in range(len(x)):
		if i == 0:
			SSR2.append(0 + (y[i]-(x[i]*fxr_slope[i]+fxr_yint[i]))**2)
		else:
			SSR2.append(SSR2[-1] + (y[i]-(x[i]*fxr_slope[i]+fxr_yint[i]))**2)

	SSR1 = []
	# Calculate SSR1 by interating through list in reverse
	for i, e in reversed(list(enumerate(x))):
		if i == len(x)-1:
			SSR1.append(0 + (y[i]-(x[i]*gxr_slope[i]+gxr_yint[i]))**2)
		else:
			SSR1.append(SSR1[-1] + (y[i]-(x[i]*gxr_slope[i]+gxr_yint[i]))**2)
	SSR1.reverse()

	sumSSR = []
	for i in range(len(x)):
		sumSSR.append(SSR1[i] + SSR2[i])
		
	return sumSSR