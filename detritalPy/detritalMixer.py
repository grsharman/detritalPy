# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 21:37:32 2020

@author: sjohnstone, grsharman

"""

from detritalpy import detritalFuncs as dFunc # Working with detrital geochron data
import numpy as np # For efficient manipulation of arrays of numbers
from matplotlib import pyplot as plt # For making plots
from scipy.optimize import minimize # For optimizations (e.g., finding best fits)
import pandas as pd # For working with pandas dataframes/ loading data
import time # For keeping track of model progress
from scipy.interpolate import interp1d

AVAILABLE_METRICS = ['dmax','vmax','similarity-pdp','likeness-pdp','r2-pdp', 'similarity-kde','likeness-kde','r2-kde', 'ss']

########################################################################################################################
## Helper functions
########################################################################################################################

def resample_datesErrors(dates, errors, n = None):
    '''
    This function resamples (with replacement) the dates and corresponding
    errors in the provided arrays. 

    Parameters
    ----------
    dates : np.ndarray
        An array of dates.
    errors : np.ndarray
        An array the same length as dates of the corresponding uncertainties.
    n : int, optional
        The number of items to draw in the resampling, by default will
        use the total number of dates available. The default is None.


    Returns
    -------
    dates_rs : np.ndarray
        An array of dates contained within the supplied array, but resampled
        with replacement.
    errors_rs : np.ndarray
        An array of dates contained within the supplied array, but resampled
        with replacement.

    '''
    
    #If the default value of n was given, assign the total number of entries as n
    if n is None:
      n = len(dates)
    
    #Get a random index of which entries to resample
    randIdx = np.random.randint(0,len(dates)-1,size = n)
    
    #resample the data
    dates_rs,errors_rs = dates[randIdx],errors[randIdx]
    
    return dates_rs,errors_rs

def resample_and_perturb_datesErrors(dates, errors, n=None):
    '''
    This function resamples (with replacement) the dates and corresponding
    errors in the provided arrays. In addition, it adds variability by perturbing
    the original dates by an amount randomly drawn from a normal distribution 
    with a standard deviation equal to the supplied error, which must be at the
    1-sigma level. If using 2-sigma uncertainties in the input data, you must first
    convert to 1-sigma prior to using this function.

    Parameters
    ----------
    dates : np.ndarray
        an array of dates.
    errors : np.ndarray
        an array the same length as dates of the corresponding 1-sigma uncertainties.
    n : int, optional
        the number of items to draw in the resampling, by default will
    use the total number of dates available. The default is None.

    Returns
    -------
    dates_rs : np.ndarray
        An array of dates contained within the supplied array, but resampled
        with replacement and perturbed based on assigned uncertainties.
    errors_rs : np.ndarray
        An array of dates contained within the supplied array, but resampled
        with replacement.

    '''
    
    dates_perturbed = perturb_by_normal_distribution(dates, errors)
    
    return resample_datesErrors(dates_perturbed,errors,n=n)

def perturb_by_normal_distribution(means,std_dev):
    '''
    This function perturbs a series of mean values by a random amount drawn
    from a normal distribution with the specified standard deviations

    Parameters
    ----------
    means : np.ndarray
        Array of mean values .
    std_dev : np.ndarray
        standard deviations associated with each mean.

    Returns
    -------
    means_pert : np.ndarray
        original mean values summed with random values drawn from a normal distribution with the standard deviation
        specified in std_dev.

    '''
    
    means_pert = means + std_dev*np.random.standard_normal(len(means))
    
    return means_pert


def lookup_functions_for_metric(objective_metric = 'Dmax',
                               x1=0, x2=4500, xdif = 1, bw = 2, bw_x=None):
    '''
    This function creates functions to compare and calculate distributions.

    Parameters
    ----------
    objective_metric : string, optional
        The name of the objective_metric to use. The default is 'Dmax'. For a full list of values call print(dMix.AVAILABLE_METRICS)
    x1 : float, optional
        Minimum value of age axis distributions are evaluated on (Ma). The default is 0.
    x2 : float, optional
        Maximum value of age axis distributions are evaluated on (Ma). The default is 4500.
    xdif : float, optional
        Spacing of age axis distributions are evaluated along, in Myr. The default is 1.
    bw : float, optional
        bandwidth used when calculating KDE. The default is 2.
    bw_x : float, optional
         X-axis location of bandwidth split (only used if multiple KDE values are specified). The default is None.

    Returns
    -------
    distribution_function : callable
        The function needed to create the distribution type (e.g., PDF vs CDF) expected by the objective_metric.
        This function accepts two inputs (ages: np.ndarray, errors: np.ndarray).
    comparison_function : callable
        The function that performs the comparison requested by objective_metric.
        This function accepts two inputs (distribution1: np.ndarray, distribution1: np.ndarray).
    areGoodFitsSmallValues : bool
        Some objective_metrics report good fits as small values, others do the opposite. This boolean keeps track of this for
        the function requested.

    '''
    
    if not(objective_metric.lower() in AVAILABLE_METRICS):
        raise Warning('Requested metric not found. The '+
                  'available metrics are: {}'.format(AVAILABLE_METRICS))
        
        distribution_function = None
        comparison_function = None
        
    elif objective_metric.lower() == 'dmax':
        
        distribution_function = lambda ages,errors: dFunc.CDFcalcAges(ages,
                                                                      x1,x2,
                                                                      xdif)
        comparison_function = dFunc.calcDmax
        areGoodFitsSmallValues = True
        
    elif objective_metric.lower() == 'vmax':
        
        distribution_function = lambda ages,errors: dFunc.CDFcalcAges(ages,
                                                                      x1,x2,
                                                                      xdif)
        comparison_function = dFunc.calcVmax
        areGoodFitsSmallValues = True

    elif objective_metric.lower() == 'ss':
        distribution_function = lambda ages,errors: dFunc.CDFcalcAges(ages,
                                                                      x1,x2,
                                                                      xdif)
        comparison_function = calc_ss
        areGoodFitsSmallValues = True        
        
    elif objective_metric.lower() == 'likeness-pdp':
        distribution_function = lambda ages,errors: dFunc.PDPcalcAges(ages,errors,
                                                                        x1,x2,
                                                                        xdif)
        comparison_function = dFunc.calcLikeness
        areGoodFitsSmallValues = False
        
    elif objective_metric.lower() == 'similarity-pdp':
        distribution_function = lambda ages,errors: dFunc.PDPcalcAges(ages,errors,
                                                                        x1,x2,
                                                                        xdif)
        comparison_function = dFunc.calcSimilarity
        areGoodFitsSmallValues = False
        
    elif objective_metric.lower() == 'r2-pdp':
        distribution_function = lambda ages,errors: dFunc.PDPcalcAges(ages,errors,
                                                                        x1,x2,
                                                                        xdif)
        comparison_function = dFunc.calcR2
        areGoodFitsSmallValues = False

    elif objective_metric.lower() == 'r2-kde':
        distribution_function = lambda ages,errors: dFunc.KDEcalcAges(ages,
                                                                        x1,x2,
                                                                        xdif,
                                                                        bw,bw_x)
        comparison_function = dFunc.calcR2
        areGoodFitsSmallValues = False
    
    elif objective_metric.lower() == 'likeness-kde':
        distribution_function = lambda ages,errors: dFunc.KDEcalcAges(ages,
                                                                        x1,x2,
                                                                        xdif,
                                                                        bw,bw_x)
        comparison_function = dFunc.calcLikeness
        areGoodFitsSmallValues = False

    elif objective_metric.lower() == 'similarity-kde':
        distribution_function = lambda ages,errors: dFunc.KDEcalcAges(ages,
                                                                        x1,x2,
                                                                        xdif,
                                                                        bw,bw_x)
        comparison_function = dFunc.calcSimilarity
        areGoodFitsSmallValues = False
    
    return distribution_function, comparison_function, areGoodFitsSmallValues

########################################################################################################################
## General statistical functions
########################################################################################################################

def permutation_test(sampleA,sampleB,comparison_function,n_perms = 10000, 
                     ageGoodFitsSmallValues = True):
    '''
    This is a function that performs permutation tests. This is a technique for assessing the null-hypothesis that the differences
    observed in two samples, as assessed by some comparison function, resulted from sampling a single population.
    
    Given a function 'comparison_function' that accepts two inputs (arrays of ages) and returns a value comparing those arrays,
    the permutation test provides a means of assessing if the value of that comparison function computed on the two supplied
    samples is distinct from a  random selection of two samples from the same distribution.

    The permutation test works by recomputing the comparison_function on many randomly sampled pseudo-samples drawn
    from a combined list of sampleA and sampleB.  If the comparison_function value calculated between the original
    sampleA and sampleB is similar to the comparison_function value calculated from pseudo-samples randomly selected from
    a pool of all the ages in sampleA and sampleB, then sampling effects alone might explain the differences observed 
    between sampleA and sampleB.

    Parameters
    ----------
    sampleA : np.ndarray
        An array of values (for example ages).
    sampleB : np.ndarray
        A different array of ages representing a second sample.
    comparison_function : function
        A function that accepts two inputs (each arrays of ages) and
        returns an number representing how those data compare.
        
        for example:
            def comparison_function(a,b):
                return np.mean(a) - np.mean(b)
        
    n_perms : TYPE, optional
        The number of permutations of the original data to use when evaluating
        the permutation test. The default is 10000.
        
    ageGoodFitsSmallValues : Boolean, optional
        A boolean describing whether the supplied objective
        function uses small values to describe good comparisons (True) 
        or uses large values to describe good comparisons (False).
        The default is True.

    Returns
    -------
    p_perm : TYPE
        The proportion of values from comparisons of randomly selected samples that are more different than the original comparison.
        In other words, the proportion of values under the null-hypothesis that are as different as sampleA and sampleB.
    orig_fun_val : TYPE
        The value of the comparison_function computed on the original sampleA and sampleB.
    perm_fun_vals : TYPE
        The values of the comparison_function computed for each of the random pseudo-samples drawn under the null-hypothesis.
    '''

    #Compute the original value of the comparison function
    orig_fun_val = comparison_function(sampleA,sampleB)
      
    #Preallocate space to store all the permuted values of the comparison function
    perm_fun_vals = np.zeros(n_perms)
      
    #Create a combined list of all observations
    all_data = np.hstack((sampleA,sampleB))
    
    #Create a list of labels (zeros and ones)
    all_labels = np.hstack((np.zeros(len(sampleA)), np.ones(len(sampleB))))
      
    #For each requested permutations
    for i in range(n_perms):
      
      #Shuffle the original labels - this modifies the data in place
      np.random.shuffle(all_labels) 
      
      #Get the random labelled datapoints
      rnd_A = all_data[all_labels == 0]
      rnd_B = all_data[all_labels == 1]
      
      #Calculate the comparison function for this permutation of labels
      perm_fun_vals[i] = comparison_function(rnd_A,rnd_B)
      
    #How often did the comparison function indicate a worse fit than random
    if ageGoodFitsSmallValues:
        #If good fits have small values, we are interested in how often random
        #shuffles are larger than the original function (e.g., a worse fit)
        p_perm = np.float(np.sum(perm_fun_vals > orig_fun_val)/n_perms)
    else:
        #If good fits have large values, we are interested in how often random
        #shuffles are smaller than the original function (e.g., a worse fit)
        p_perm = np.float(np.sum(perm_fun_vals < orig_fun_val)/n_perms)
        
    return p_perm, orig_fun_val, perm_fun_vals

def bootstrapped_self_comparisons_many_samples(main_byid_df,sampleList,
                                               doPerturbResampledAges=True,
                                               objective_metric='dmax',
                                               sigma='1sigma',
                                               x1=0, x2=4500,xdif=1,
                                               bw=2, bw_x=None, nBootstrapIterations=1000):
    """

    This function calls bootstrapped_self_comparison on many individual samples, listed in sampleList,
    and concatenates the results.

    bootstrapped_self_comparison is a function to assess the variability in an objective function that might result
    from random sampling.

    This function randomly resamples (with replacement) the age and uncertainty data from a single sample,
    and uses an objective function (e.g., Dmax) to compare the random resample to the original data. 

    This provides a way of assessing how objective function values calculated during mixture modelling compare
    to the variability that could be expected do to random sampling.

    Args:
        main_byid_df : dataframe:
             dataframe loaded using detritalPy's loadDataExcel function
        sampleList : list
            list of the samples to call bootstrapped_self_comparison on
        doPerturbResampledAges : Boolean, optional
            Should resampled ages be additionally randomized based on the assigned uncertainty. The default is True,
            which perturbs each age based on a draw from a normal distribution scaled to match the assigned uncertainty.
        objective_metric : String, optional
            The objective comparison function used to assess similarity between datasets. The default is 'dmax'. For a 
            full list of values call print(dMix.AVAILABLE_METRICS)
        nBootstrapIterations : int, optional
            The number of bootstrapping iterations performed. The default is 1000, for which the original sample will be compared
            to 1000 resampled versions of itself.
        x1 : float, optional
            The starting age for which the distribution of ages will be calculated and evaluated. The default is 0.
        x2 : float, optional
            The ending age for which the distribution of ages will be calculated and evaluated. The default is 4500.
        xdif : float, optional
            The spacing between values on the age axis along which distributions will be evaluated. The default is 1.
        bw : float, optional
            The bandwidth for kernels used in KDEs. The default is 2.
        bw_x : None or list, optional
            Set to None if not using a split KDE bw. Otherwise, set x-axis locations for bw split (in Ma) in a list, (e.g., bw_x=[300])

    Returns:
        selfCompMetrics_bs_set : list
            A list of arrays of length child_list. Each array of size nBootstrapIterations containing the value
            of objective_metric obtained on each randomly resampled comparison.
        childDists_bs_set: list
            A list of arrays of length child_list. Each array is of size nBootstrapIterations (rows) x the size of
            the distributions, containing each randomly resampled distribution.
    """
    
    
    ages,errors, nGrains, labels = dFunc.sampleToData(sampleList,
                                                     main_byid_df)
    
    selfCompMetrics_bs_set, childDists_bs_set = [],[]
    
    for i in range(len(sampleList)):
        res_i = bootstrapped_self_comparison(ages[i],errors[i],
                                             doPerturbResampledAges,
                                             objective_metric,
                                             x1,x2,xdif,bw,bw_x,nBootstrapIterations)
        
        selfCompMetrics_bs_set.append(res_i[0])
        childDists_bs_set.append(res_i[1])
        
    return selfCompMetrics_bs_set, childDists_bs_set

def bootstrapped_self_comparison(ages,errors,doPerturbResampledAges = True,
                                 objective_metric = 'dmax',x1=0, x2=4500,
                                 xdif = 1,bw = 2, bw_x = None, nBootstrapIterations = 1000):
    '''
    A function to assess the variability in an objective function that might result from random sampling.

    This function randomly resamples (with replacement) the age and uncertainty data from a single sample,
    and uses an objective function (e.g., Dmax) to compare the random resample to the original data. 

    This provides a way of assessing how objective function values calculated during mixture modelling compare
    to the variability that could be expected do to random sampling.

    Parameters
    ----------
    ages : np.ndarray
        A 1-D numpy array describing the ages observed for a single sample.
    errors : np.ndarray
        A 1-D numpy array describing the uncertainties associated with the ages observed for a single sample..
    doPerturbResampledAges : Boolean, optional
        Should resampled ages be additionally randomized based on the assigned uncertainty. The default is True,
        which perturbs each age based on a draw from a normal distribution scaled to match the assigned uncertainty.
    objective_metric : String, optional
        The objective comparison function used to assess similarity between datasets. The default is 'dmax'. For a 
        full list of values call print(dMix.AVAILABLE_METRICS)
    nBootstrapIterations : TYPE, optional
        The number of bootstrapping iterations performed. The default is 1000, for which the original sample will be compared
        to 1000 resampled versions of itself.
    x1 : float, optional
        The starting age for which the distribution of ages will be calculated and evaluated. The default is 0.
    x2 : float, optional
        The ending age for which the distribution of ages will be calculated and evaluated. The default is 4500.
    xdif : float, optional
        The spacing between values on the age axis along which distributions will be evaluated. The default is 1.
    bw : float, optional
        The bandwidth for kernels used in KDEs. The default is 2.
    bw_x : None or list, optional
        Set to None if not using a split KDE bw. Otherwise, set x-axis locations for bw split (in Ma) in a list, (e.g., bw_x=[300])

    Returns
    -------
    selfCompMetrics_bs : np.ndarray
        Array of size nBootstrapIterations containing the value of objective_metric obtained on each 
        randomly resampled comparison.
    childDists_bs : np.ndarray
        Array of size nBootstrapIterations (rows) x the size of
        the distributions, contains each randomly resampled distribution.
    '''

    #For the supplied distribution comparison function identified, 
    #get the functions needed to evaluate that metric
    dist_function,comp_function, areGoodFitsSmall = lookup_functions_for_metric(objective_metric,
                                                              x1,x2,xdif,bw,bw_x)

    if doPerturbResampledAges:
      #If we are perturbing the resampled ages use the appropriate function
      resampling_function = resample_and_perturb_datesErrors
    else:
      #If we are only resampling
      resampling_function = resample_datesErrors
    
    #First preallocate space to store the self-compared comparison metrics
    selfCompMetrics_bs = np.zeros(nBootstrapIterations)
    
    #Determine how much space needs to be preallocated by making distribution
    test_dist = dist_function([ages],[errors])[0]
    
    #Lets also preallocate space to store the resampled distributions
    childDists_bs = np.zeros((nBootstrapIterations,len(test_dist)))
    
    #Create the original child distribution
    dist_age_axis,child_dist = dist_function([ages],[errors])
    
    #For each bootstrap iterations
    for i in range(nBootstrapIterations):
    
        #Resample with replacement the child
        child_ages_rs,child_errors_rs = resampling_function(ages,
                                                            errors)
        
        #Calculate the distribution for this data
        child_dist_rs = dist_function([child_ages_rs],
                                              [child_errors_rs])[1]
        
        #Compare the two values and store the comparison
        selfCompMetrics_bs[i] = comp_function(child_dist[0],
                                                    child_dist_rs[0])
        
        #Store also this resampled distribution
        childDists_bs[i] = child_dist_rs[0]
      
    return selfCompMetrics_bs, childDists_bs

############################################################################
## Mixing functions
############################################################################

def mix_parents(parent_distributions,mixing_coefficients):
    '''
    Preforms a simple linear mixture of the distributions stored in parent distributions based on supplied mixing coefficients.

    Parameters
    ----------
    parent_distributions : np.ndarray
        Array of size len(mixing_coefficients) x len(distribution).
    mixing_coefficients : np.ndarray
        Mixing coefficients for each.

    Returns
    -------
    mixed_dist : np.ndarray
        Array of size len(distribution). The sum of the product of each mixing coefficient and its corresponding parent distribution.
    '''

    #Start off with an empty mixed CDF
    mixed_dist = np.zeros_like(parent_distributions[0])
  
    #Add together each of the mixed parents
    for P_i,phi_i in zip(parent_distributions,mixing_coefficients):
      mixed_dist += phi_i*P_i
  
    #Return the mixed CDF
    return mixed_dist

def find_best_fit_mix(parent_list, child_list, main_byid_df, sigma='1sigma',
                     objective_metric='dmax', x1=0, x2=4500, xdif=1, bw=2, bw_x=None, sampleLabel = 'Sample_ID',
                     verbose=False):

    '''
    Determine the best fitting mixing coefficients through minimization for a list of children. The goodness of fit for the
    minimization is assessed by the supplied objective_metric.  Distributions are calculated based on 
    the parameters x1,x2, xdif, bw, and bw_x

    Parameters
    ----------
    parent_list : list
        List of the samples to use as parents for mixture modelling.
    child_list : TYPE
        List of the samples to use as children for mixture modelling..
    main_byid_df : dataframe:
        dataframe loaded using detritalPy, e.g., with the loadDataExcel function
    sigma : string, optional
        Uncertainly level of input data (options: '1sigma' or '2sigma'). The default is '1sigma'.
    objective_metric : string, optional
        Name of objective_metric that will be used to compare distributions. The default is 'dmax'. For a 
        full list of values call print(dMix.AVAILABLE_METRICS)
    x1 : float, optional
        Minimum value of age axis distributions are evaluated on (Ma). The default is 0.
    x2 : float, optional
        Maximum value of age axis distributions are evaluated on (Ma). The default is 4500.
    xdif : float, optional
        Spacing of age axis distributions are evaluated along (Myr). The default is 1.
    bw : float, optional
        bandwidth used when calculating KDE. The default is 2.
    bw_x : None or list, optional
        Set to None if not using a split KDE bw. Otherwise, set x-axis locations for bw split (in Ma) in a list, (e.g., bw_x=[300])
    sampleLabel : string, optional
        The field which contains labels. Default = 'Sample_ID'

    Returns
    -------
    mix_coeffs_bf : list
        A list of arrays, each array contains the best fit mixing coefficients.
    obj_func_val : list
        A list of the values of the objective_metric calculated with the best fit mixing coefficients.
    best_mixed_dist : list
        A list of arrays, each array contains the best fitting mixed distribution.

    '''

    start = time.time()

    dist_function, comp_function, areGoodFitsSmallValues = lookup_functions_for_metric(objective_metric,
                                                                                          x1,x2,xdif,bw,bw_x)

    #And load the data corresponding to that parent
    ages_p, errors_p, numGrains_p, labels_p = dFunc.sampleToData(parent_list,
                                                         main_byid_df,
                                                         sigma = sigma,
                                                         sampleLabel=sampleLabel)

    #Load all the child data
    ages_c, errors_c, numGrains_c, labels_c = dFunc.sampleToData(child_list,
                                                     main_byid_df,
                                                     sigma = sigma,
                                                     sampleLabel=sampleLabel)

    mix_coeffs_bf = []
    obj_func_val = []
    best_mixed_dist = []

    for i in range(len(child_list)):
        if type(child_list[0])==tuple:
            print('Starting:',child_list[i][1])
        else:
            print('Starting:',child_list[i])
        mix_coeffs_bf_i, obj_func_val_i, best_mixed_dist_i = find_best_fit_parent_mixture(ages_p,errors_p,ages_c[i],errors_c[i],
                           dist_function, comp_function,
                           areGoodFitsSmallValues,
                           verbose=verbose)
        mix_coeffs_bf.append(mix_coeffs_bf_i)
        obj_func_val.append(obj_func_val_i)
        best_mixed_dist.append(best_mixed_dist_i)

    print('Completed mixing children.')
    print('Time elapsed:',str(round(time.time() - start, 3)),'seconds')

    return mix_coeffs_bf, obj_func_val, best_mixed_dist

def find_best_fit_parent_mixture(parentAges,parentErrors,childAges,childErrors,
                   distributionCreationFunction,objectiveFunction,
                   areGoodFitsSmallValues = True, verbose=False):
    '''
    Determine the best fitting mixing coefficients through minimization for a single child sample. The goodness of fit for the
    minimization is assessed by the supplied objective_metric.  Distributions are calculated based on 
    the parameters x1,x2, xdif, bw, and bw_x

    Parameters
    ----------
    parentAges : list
        a list of n arrays storing the dates sampled for each parent.
    parentErrors : list
        a list of n arrays storing the errors from each parent sample.
    childAges : np.ndarray
        an array of the dates measured for the child sample.
    childErrors : np.ndarray
        an array of the error of each measured date for the child sample.
    distributionCreationFunction : callable
        a function that takes as input a list of
        arrays of ages and errors and returns a a list of distributions necessary for
        the supplied objective function.
        Can be created with dMix.lookup_functions_for_metric(objective_metric,x1,x2,xdif,bw,bw_x)
    objectiveFunction : callable
        a function that compares a modeled child distribution to 
        an observed one. Needs to accept as input two distributions of the type
        created by the distribution creation function. 
        Can be created with dMix.lookup_functions_for_metric(objective_metric, x1,x2,xdif,bw,bw_x)
    areGoodFitsSmallValues : bool, optional
        a boolean describing whether the supplied objective
        function uses small values to describe good comparisons (True) or uses large
        values to describe good comparisons (False). The default is True.

    Returns
    -------
    mix_coeffs_bf : np.ndarray
        Array containing the best fit mixing coefficients.
    obj_func_val : float
        The value of the objective function for the best fit mixture.
    best_mixed_dist : TYPE
        The best fitting mixed distribution.

    '''

    
    ########################
    #First, prep by creating all the necessary distributions to mix
    
    #We will combine all the data into one list, to minimize the overhead
    #created by multiple calls to detritalPys distribution creation functions
    ########################
    
    #Create a list of the data with the child data at the end
    allAges = [p for p in parentAges]
    allErrors = [e for e in parentErrors]
    allAges.append(childAges)
    allErrors.append(childErrors)
    
    #Create the distributions to mix
    ageAxis,distributions = distributionCreationFunction(allAges,allErrors)
    
    #Split out the parent and child distributions
    parentDistributions = distributions[:-1]
    childDistribution = distributions[-1]
    
    ########################
    #Next, set up the optimization approach 
    
    #We will define our constraints, initial guess, and wrap our
    #objective function to accept a mixed distribution
    ########################
    
    # All values must sum to zero
    sumConst = ({'type': 'eq', 'fun': lambda x: 1 - sum(x)})
    
    # All mixing coefficients must be positive, b/w 0 & 1
    nParents = len(parentAges)
    boundingVals = [(0, 1) for i in range(nParents)]
    
    #Create the function to minimize, taking care to ensure that the minimization
    #approach identifies the best mixture (as opposed to the worse)
    if areGoodFitsSmallValues:
        objFun = lambda x: objectiveFunction(childDistribution,
                                            mix_parents(parentDistributions,x))
    else:
        objFun = lambda x: 1.0/objectiveFunction(childDistribution,
                                            mix_parents(parentDistributions,x))
    
    #Create an initial guess - that everything was mixed evenly
    initialGuess = np.ones(nParents)/nParents
    
    
    ########################
    #Finally, preform the actual minimization and store its predictions for return
    ########################
    
    #Minimize the function
    result = minimize(objFun, initialGuess, method='SLSQP',
                                bounds=boundingVals,
                              constraints=sumConst,
                              options={'disp' : verbose, 'maxiter' : 1e6})

    # result = minimize(objFun, initialGuess, method='trust-constr',
    #     bounds=boundingVals, constraints=sumConst, options={'maxiter' : 1e4})

    mix_coeffs_bf = result.x

    #print(result)

#    mix_coeffs_bf = minimize(objFun, initialGuess, method='Powell', bounds=boundingVals, constraints=sumConst).x#, tol=1e-20, options={'maxiter' : 1e6, 'disp' : False}).x
    
    #Calculate and store the best fitting mixture as well
    best_mixed_dist = mix_parents(parentDistributions,mix_coeffs_bf)
    
    #Find the best fitting value of the mixing coefficient
    obj_func_val = objectiveFunction(childDistribution,
                                     best_mixed_dist)
    
    #Return the mix. coefficients, obj. function value, and mixed distribution
    return mix_coeffs_bf, obj_func_val, best_mixed_dist

def bootstrap_find_best_mixtures(parentAges,parentErrors,
                                childAges,childErrors,
                                 distributionCreationFunction,
                                objectiveFunction,
                                areGoodFitsSmallValues,
                                nBootstrapIterations = 1000,
                                doPerturbResampledAges = True,
                                nGrainsToResample = None,
                                verbose = False,
                                start = None,
                                update_freq = None):
    '''
    This function determines a series of bootstrapped estimates for the mixing
    coefficients that best describe how the supplied parentAges and Errors combine
    to describe the supplied childAges. 

    Parameters
    ----------
    parentAges : list
        A list of n arrays storing the dates sampled for each parent.
    parentErrors : list
        A list of n arrays storing the errors from each parent sample.
    childAges : np.ndarray
        An array of the dates measured for a child sample.
    childErrors : np.ndarray
        An array of the error of each mearued date from the child sample.
    distributionCreationFunction : function
        A function that takes as input a list of
        arrays of ages and errors and returns a a list of distributions
        necessary for the supplied objective function.
    objectiveFunction : function
        A function that compares a modeled child distribution to 
        an observed one. Needs to accept as input two distributions of the
        type created by the distribution creation function.
    areGoodFitsSmallValues : Boolean
        A boolean describing whether the supplied objective
        function uses small values to describe good comparisons (True) 
        or uses large values to describe good comparisons (False).
    nBootstrapIterations : int, optional
        An integer specifying how many times to repeat the resampling with 
        replacement determination of mixing coefficients. The default is 1000.
    doPerturbResampledAges : Boolean, optional
        if True original dates are perturbed by  arandom amount, the magnitude
        of which is drawn from a normal distribution described by the
        corresponding error on that date. If False, the original dates are 
        resampled as-is. The default is True.
    nGrainsToResample : int, optional
        The number of grains to resample with each bootstrap iteration. The default is None,
        in which case the original number of grains is used.
    verbse : Boolean, optional
        Set to True to print code progress. Default = False.
    start : time.time() input or None, optional
        Input from the time.time() function (starting time of process). Default is None.
    update_freq : integer or None, optional
        Frequency of how often code progress is printed. Default is None.

    Returns
    -------
    mixing_coeffs_bs : np.ndarray
        np.ndarray of size nBootstrapIterations x len(parentAges). Contains the best fit mixing coefficients for each
        bootstrap iteration.
    obj_fun_vals_bs : np.ndarray
        np.ndarray of size nBootstrapIterations. Contains the values of the objective_metric calculated for each bootstrap iteration.
    best_fitting_dists_bs : np.ndarray
        np.ndarray of size nBootstrapIterations x len(distributions). Contains the best fitting mixed distribution for each
        bootstrap iteration.
    '''

    #Preallocate space for the outputs
    mixing_coeffs_bs = np.zeros((nBootstrapIterations,len(parentAges)))
    obj_fun_vals_bs = np.zeros(nBootstrapIterations)
    
    #Determine how big the output distribution should be
    dist_axis,dist_vals = distributionCreationFunction([childAges],
                                                       [childErrors])
    best_fitting_dists_bs = np.zeros((nBootstrapIterations,len(dist_axis)))
    
    #Which function should we use for resampling? We can execute this statement
    #outside the loop to simplify each iteration
    if doPerturbResampledAges:
        #If we are perturbing the resampled ages use the appropriate function
        resampling_function = resample_and_perturb_datesErrors
    else:
        #If we are only resampling
        resampling_function = resample_datesErrors
    
    #Now preform all the bootstrap iterations requested
    for i in range(nBootstrapIterations):
      
        #Resample each set of parent ages - first create a list to store the values
        parentAges_rs = []
        parentErrors_rs = []
        #For each parent's dates and errors
        for p,e in zip(parentAges,parentErrors):
            #Resample the values
            p_rs,e_rs = resampling_function(p,e,nGrainsToResample)
            #Add these to the list
            parentAges_rs.append(p_rs)
            parentErrors_rs.append(e_rs)
          
        #Resample the child dates and errors
        childDates_rs, childErrors_rs = resampling_function(childAges,
                                                            childErrors)
      
        #Find the best fitting values and store these
        res = find_best_fit_parent_mixture(parentAges_rs,parentErrors_rs,
                                           childDates_rs,childErrors_rs,
                                           distributionCreationFunction,
                                           objectiveFunction,
                                           areGoodFitsSmallValues)
      
        #Unpack the results into the respective arrays
        mixing_coeffs_bs[i,:],obj_fun_vals_bs[i],best_fitting_dists_bs[i,:] =  res

        if verbose:
            total_time = time.time() - start
            avg_freq = total_time/(i+1)
            if i%update_freq==0:
                print(round(i/nBootstrapIterations*100,4),'% completed,',
                    round(total_time/60,3),'minutes elapsed,',
                    round((nBootstrapIterations-i)*avg_freq/60, 3),'minutes remaining')
      
    #Return the results
    return mixing_coeffs_bs,obj_fun_vals_bs,best_fitting_dists_bs

def calc_model_fit(child_list, obj_func_val, obj_vals_all, selfCompMetrics_bs_set, objective_metric = 'dmax',
                     percentilesToSummarize=[50, 2.5, 97.5], formatString='{:.2f} ({:.2f} - {:.2f})', verbose=False,
                     alpha = 0.05):
    '''
    This function determines the goodness-of-fit of the mixture modeling by comparing  the variability in objective_metric
    values observed in bootstrapped mixture models to the variability in objective_metric values observed when comparing a sample
    to itself (after resampling).

    Parameters
    ----------
    child_list : list
        A list of each of the modeled children.
    obj_func_val : list
        A list of the values of the objective_metric calculated with the best fit mixing coefficients.
    obj_vals_all : list of arrays
        A list of len(child_list) containing arrays of the values of the objective_metric calculated during bootstrapped mixture
        modelling.
    selfCompMetrics_bs_set : list of arrays
        A list of len(child_list) containing arrays of the values of the objective_metric calculated by comparing each child sample
        to a resampled version of itself.
    objective_metric : TYPE, optional
        Name of objective_metric that will be used to compare distributions. The default is 'dmax'. For a 
        full list of values call print(dMix.AVAILABLE_METRICS)
    percentilesToSummarize : list, optional
        Percentile values to report of the obj. function distributions. Default = [50, 2.5, 97.5].
    formatString : string, optional
        How do we want to print the strings of the obj. function distributions. Default = '{:.2f} ({:.2f} - {:.2f})'.
    verbose : bool, optional
        Should summaries be printed out or just saved. Default = False.
    alpha : float (must be < 1)
        The decision criteria used for summarizing the maximum deviation between the model and observed child allowed. 
        For example, if alpha = 0.05, we will determine the critical objective function value based on the worst 5% of
        objective_metric values observed when comparing a sample to a resampled version of itself.

    Returns
    -------
    obj_fun_crit : np.ndarray
        'critical' value of the objective_metric for each value. This is calculated as the alpha-th percentile
        of self-compared objective metrics.
    worse_than_crit : float
        The proportion of bootstrapped mixture model results that have objective_metric values suggesting a greater deviation
        than obj_fun_crit.
    '''

    nBootstrapIterations = len(obj_vals_all[0])

    areGoodFitsSmallValues = lookup_functions_for_metric(objective_metric=objective_metric)[2]

    #Store summary values
    self_samples_comp_summs = np.zeros((len(child_list),len(percentilesToSummarize)))
    model_comp_summs = np.zeros_like(self_samples_comp_summs)
    obj_func_crit = np.zeros(shape=(len(child_list)))
    worse_than_crit = np.zeros(shape=(len(child_list)))

    ## Loop through all the children and report these results
    for i, child in enumerate(child_list):
        these_obj_vals_mod = obj_vals_all[i]
        these_obj_vals_selfComp = selfCompMetrics_bs_set[i]
        
        #Summarize distributions
        perc_obj_vals_mod = np.percentile(these_obj_vals_mod,percentilesToSummarize)
        perc_obj_vals_selfComp = np.percentile(these_obj_vals_selfComp,percentilesToSummarize)
        
        #Store these vals
        self_samples_comp_summs[i,:] = perc_obj_vals_selfComp
        model_comp_summs[i,:] = perc_obj_vals_mod
        
        #How do the model and observed distributions overlap?
        if areGoodFitsSmallValues:
            mod_obs_percentile = 100.0*np.sum(1.0*(these_obj_vals_mod > these_obj_vals_selfComp))/nBootstrapIterations
            obj_func_crit[i] = np.percentile(these_obj_vals_selfComp, 100 - (alpha*100))
            worse_than_crit[i] = np.sum(these_obj_vals_mod>obj_func_crit[i])/nBootstrapIterations
            if verbose:
                print(child)
                print('best-fit ',objective_metric+': ','{:.3f}'.format(obj_func_val[i]))
                print(objective_metric,'_crit: ','{:.3f}'.format(obj_func_crit[i]))
                if obj_func_crit[i] > obj_func_val[i]:
                    print('Best-fit mixture is better than bootstrapped child')
                else:
                    print('Best-fit mixture is worse than bootstrapped child')
        else:
            mod_obs_percentile = 100.0*np.sum(1.0*(these_obj_vals_mod < these_obj_vals_selfComp))/nBootstrapIterations
            obj_func_crit[i] = np.percentile(these_obj_vals_selfComp, alpha*100)
            worse_than_crit[i] = np.sum(these_obj_vals_mod<obj_func_crit[i])/nBootstrapIterations
            if verbose:
                print(child)
                print('best-fit ',objective_metric+': ','{:.3f}'.format(obj_func_val[i]))
                print(objective_metric,'_crit: ','{:.3f}'.format(obj_func_crit[i]))
                if obj_func_crit[i] < obj_func_val[i]:
                    print('Best-fit mixture is better than bootstrapped child')
                else:
                    print('Best-fit mixture is worse than bootstrapped child')

    if verbose:
        print('Finished!')
    return obj_func_crit, worse_than_crit

def bootstrap_solve_mixture_models(parent_list,child_list,main_byid_df,
                                    sigma = '1sigma',
                                    objective_metric = 'dmax',
                                    nBootstrapIterations = 1000, 
                                    doPerturbResampledAges = True,
                                    nGrainsToResample = None,
                                    x1=0, x2=4500,xdif = 1,bw = 2, bw_x=None,
                                    sampleLabel='Sample_ID',
                                    verbose=True,
                                    update_freq = 10):
    '''
    Determine the best fitting mixing coefficients for a suite of child samples, given the parents supplied in parent_list. 

    Rather than solve for the single 'best' solution, this function characterizes the uncertainty in mixture coefficients by bootstrapping.
    For a specified number of iterations (nBootstrapIterations) the parent and child samples and redrawn with replacement (and ages are 
    optionally perturbed by an amount dictated by their analytical uncertainty), and mixing coefficients are fit to these redrawn samples.
    This procedure provides an estimate of how sampling effects introduce uncertainty in our ability to infer the correct mixture.
    

    Parameters
    ----------
    parent_list : list or tuple
        List or tuple of the samples to use as parents for mixture modelling.
    child_list : list or tuple
        List or tuple of the samples to use as children for mixture modelling.
    main_byid_df : dataframe:
        dataframe loaded using detritalPy, e.g., with the loadDataExcel function
    sigma : string, optional
        Uncertainly level of input data (options: '1sigma' or '2sigma'). The default is '1sigma'.
    objective_metric : string, optional
        Name of objective_metric that will be used to compare distributions. The default is 'dmax'. For a 
        full list of values call print(dMix.AVAILABLE_METRICS)
    nBootstrapIterations : int, optional
        How many bootstrapping iterations (resampling and recomputing mixing coefficients) should be preformed. The default is 1000.
    doPerturbResampledAges : bool, optional
        Should resampled ages be additionally randomized by drawing a perturbation from a normal distribution scaled to the analytical
        uncertainty? The default is True.
    nGrainsToResample : int, optional
        How many grains should be resampled from each dataset. The default is None, in which case the original number of grains is used.
    x1 : float, optional
        Minimum value of age axis distributions are evaluated on. The default is 0.
    x2 : float, optional
        Maximum value of age axis distributions are evaluated on. The default is 4500.
    xdif : float, optional
        Spacing of age axis distributions are evaluated along. The default is 1.
    bw : float, optional
        bandwidth used when calculating KDE. The default is 2.
    bw_x : None or list, optional
        Set to None if not using a split KDE bw. Otherwise, set x-axis locations for bw split (in Ma) in a list, (e.g., bw_x=[300])
    verbose : bool
        Set to True to print out a readout of progress. Default = True.
    update_freq : integer, optional
        How frequent to update progress (only applies if verbose = True)

    Returns
    -------
    mix_coeffs_all : list
        The suite of bootstrapped mixing coefficients for each children. This is a list of len(child_list), consisting of arrays of size
        len(parent_list) x nBootstrapIterations.
    obj_vals_all : list
        A list of len(child_list) containing the suite of objective function values for each bootstrapped iteration.
    Cmod_all : TYPE
        A list of len(child_list) containing arrays of len(distributions) x nBootstrapIterations representing the best fit mixture modelled child
        distribution for each bootstrap iteration.
    '''

    #For the supplied distribution comparison function identified, 
    #get the functions needed to evaluate that metric
    dist_function, comp_function, areGoodFitsSmallValues = lookup_functions_for_metric(objective_metric,
                                                                                      x1,x2,xdif,bw,bw_x)
    
    #And load the data corresponding to that parent
    ages_p, errors_p, numGrains_p, labels_p = dFunc.sampleToData(parent_list,
                                                         main_byid_df,
                                                         sigma = sigma,
                                                         sampleLabel=sampleLabel)
    
    #Load all the child data
    ages_c, errors_c, numGrains_c, labels_c = dFunc.sampleToData(child_list,
                                                     main_byid_df,
                                                     sigma = sigma,
                                                     sampleLabel=sampleLabel)
    
    #Make lists for each of the outputs
    mix_coeffs_all = [None for c in child_list]
    obj_vals_all = [None for c in child_list]
    Cmod_all = [None for c in child_list]
    
    #For each of the requested children
    for i,(age_c,err_c) in enumerate(zip(ages_c,errors_c)):

        if verbose:
            print('-----------------------------------------------------------------')
            print('Starting :',labels_c[i], '(', i+1,'of',len(ages_c),'child samples)')
            start = time.time()
        else:
            start = None
            update_freq = None

        #Run all the bootstrap iterations
        res = bootstrap_find_best_mixtures(ages_p,errors_p,
                                           age_c,err_c,
                                           dist_function,
                                           comp_function,
                                           areGoodFitsSmallValues,
                                           nBootstrapIterations,
                                           doPerturbResampledAges,
                                           nGrainsToResample,
                                           verbose=verbose,
                                           start=start,
                                           update_freq = update_freq)
    
        #Get all the bootstrapped values from the result
        mix_coeffs_all[i],obj_vals_all[i],Cmod_all[i] = res

    if verbose:
        print('Finished bootstrapping!')

    return mix_coeffs_all,obj_vals_all,Cmod_all

############################################################################
## Plotting functions
############################################################################

def print_best_fit_mixture(parent_list, child_list, objective_metric, obj_func_val, mix_coeffs_bf):
    """
    A function that prints the best-fit results of a mixture model.

    Parameters
    ----------
    parent_list : list or tuple
        List or tuple of the samples to use as parents for mixture modelling.
    child_list : list or tuple
        List or tuple of the samples to use as children for mixture modelling.   
    objective_metric : string, optional
        Name of objective_metric that will be used to compare distributions. The default is 'dmax'. For a 
        full list of values call print(dMix.AVAILABLE_METRICS)
    obj_vals_all : list
        A list of len(child_list) containing the suite of objective function values for each bootstrapped iteration.
    mix_coeffs_bf : list
        A list of arrays, each array contains the best fit mixing coefficients.

    Returns
    -------
    A printout of best-fit mixture model results
    """
    if type(child_list[0]) == tuple:
        child_names = [x[1] for x in child_list]
    if type(child_list[0]) == str:
        child_names = [x for x in child_list]

    for i, child in enumerate(child_names):
        print('Sample:',child)
        print('Best-fit',objective_metric,':',np.round(obj_func_val[i],3))
        print('Mixing coefficients')
        for j in range(len(parent_list)):
            if type(parent_list[0]) == str:
                print('---',parent_list[j],':', np.round(mix_coeffs_bf[i][j],3))            
            if type(parent_list[0]) == tuple:
                print('---',parent_list[j][1],':', np.round(mix_coeffs_bf[i][j],3))

def plot_bootstrapped_mixturecoefficients_stratigraphy(parent_list,
                                                       child_list,
                                                       mix_coeffs_bf,
                                                       mix_coeffs_bs_set,
                                                       confidence_interval = 95,
                                                       ax = None,
                                                       plotWidth = 4.0,
                                                       plotHeight = 8.0,
                                                       parent_colors = 'Default',
                                                       do_plot_errorbars = True,
                                                       yAxisValues = None,
                                                       doFlipXY = False,
                                                       plot_alpha = 0.5,
                                                       separate_parents = False,
                                                       best_plot_type = 'best-fit'):
    """
    A function to create a plot that summarizes the bootstrapped distribution of mixing coefficients as shaded regions
    and how they change for a suite of children. For example, how the propotion of one parent changes with height in a stratigraphic section.

    For examples of this plot see, for example, Fig 8 from:
    Malkowski, M. A., et al. "Continental shelves as detrital mixers: U–Pb and Lu–Hf detrital zircon provenance of the
    Pleistocene–Holocene Bering Sea and its margins." The Depositional Record (2022).

    Parameters
    ----------
    parent_list : list or tuple
        List or tuple of the samples to use as parents for mixture modelling.
    child_list : list or tuple
        List or tuple of the samples to use as children for mixture modelling.
    mix_coeffs_bf : list
        A list of arrays, each array contains the best fit mixing coefficients.
    mix_coeffs_bs_set : list
        The list of bootstrapped mixing coefficients returned by bootstrap_solve_mixture_models
    confidence_interval : float (in the interval 0 - 100).
        The percentile of bootstrapped mixing coefficients to summarize. Defaults to 95, infills between the 2.5th and 97.5th percentile.
    ax : matplotlib axis:
        The matplotlib axis to plot on. Defaults to None, in which case an axis is created.
    plotWidth : float, optional
        The width of the plot. Default = 4.0
    plotHeight : float, optional
        The height of the plot. Default = 8.0
    parent_colors : str or list, optional
        The list of len(parent_names) of matplotlib colors to color mixing coefficients of each parent by.
        May also use a single color (e.g., 'black'). Default value is 'Default'
        in which case the detritalpy default colors are used.
    do_plot_errorbars : bool, optional
        Should error bars be added at each child location, or should the plot just consist of filled ranges. Defaults to True.
    yAxisValues : array of len(child_names)
        Array of the y axis position to plot children at. Defaults to None, in which case children are plotted with regular spacing.
    doFlipXY : bool, optional
        Should the x and y axis be swapped, so that the yAxisValues are plotted on the horizontal and mixing coefficients the vertical. Defaults to False.
    plot_alpha : float (0-1), optional
        Specifies the transparency of the parent shading
    separate_parents : Boolean, optional
        Set to True to plot each parent on its own subplot. Default = False.
    best_plot_type : str, optional
        Determins whether the best-fit mixture model or the average or median of all models are displayed.
        Options: 'best-fit', 'average', 'median'. Default = 'best-fit'

    Returns
    -------
        ax : matplotlib axis
            The axis on which the plot was created
    """
    #If no axis specified, create one:
    if ax is None:
        if separate_parents:
            if doFlipXY:
                f,axs = plt.subplots(len(parent_list), 1, figsize=(plotWidth, plotHeight))
            else:
                f,axs = plt.subplots(1, len(parent_list), figsize=(plotWidth, plotHeight))
        else:
            f,ax = plt.subplots(1,1, figsize=(plotWidth, plotHeight))

    if type(parent_list[0]) == tuple:
        parent_names = [parent[1] for parent in parent_list]
    if type(parent_list[0]) == str:
        parent_names = [parent for parent in parent_list]
    if type(child_list[0]) == tuple:
        child_names = [child[1] for child in child_list]
    if type(child_list[0]) == str:
        child_names = [child for child in child_list]
        
    #If no colors were specified, get them from detritalPy
    if parent_colors == 'Default':
        parent_colors = [dFunc.colorMe(x) for x in np.arange(len(parent_list))]

    if type(parent_colors) == str:
        parent_colors = [parent_colors] * len(parent_list)

    #Get the upper and lower percentile to report
    alpha = (100.0 - confidence_interval)/2.0
    
    #Preallocate some space to store the upper and lower bounds of each mix_coeff
    mix_coeff_UB = np.zeros((len(child_names),len(parent_names)))
    mix_coeff_LB = np.zeros_like(mix_coeff_UB)
    
    #If no y axis values were specified
    if yAxisValues is None:
        yAxisValues = 1+ np.arange(len(child_names))
    
    for i,(child_name, mix_coeffs_bs) in enumerate(zip(child_names,mix_coeffs_bs_set)):    
        for j, parent_name in enumerate(parent_names):
            if best_plot_type == 'median':
                percs = np.percentile(mix_coeffs_bs[:,j],[alpha, 50, 100 - alpha])
            if best_plot_type == 'average':
                percs = np.percentile(mix_coeffs_bs[:,j],[alpha, 100 - alpha])
                avg = np.average(mix_coeffs_bs[:,j])
                percs = np.insert(percs, 1, avg)
            if best_plot_type == 'best-fit':
                percs = np.percentile(mix_coeffs_bs[:,j],[alpha, 100 - alpha])
                bf = mix_coeffs_bf[i][j]
                percs = np.insert(percs, 1, bf)

            mix_coeff_UB[i,j] = percs[-1]
            mix_coeff_LB[i,j] = percs[0]
            
            if do_plot_errorbars:
                errbnds = np.zeros((2,1))
                errbnds[:,0] = [np.abs(percs[1] - percs[0]),np.abs(percs[2] - percs[1])]
                
                if separate_parents:
                    if not doFlipXY:
                        axs[j].errorbar(percs[1],yAxisValues[i],yerr = None, xerr = errbnds,
                                    fmt = 'ok',ecolor = 'k',elinewidth = 1,
                                    markerfacecolor = parent_colors[j],markeredgecolor = None)
                    else:
                        axs[j].errorbar(yAxisValues[i],percs[1],yerr = errbnds, xerr = None,
                                    fmt = 'ok',ecolor = 'k',elinewidth = 1,
                                    markerfacecolor = parent_colors[j],markeredgecolor = None)
                else:
                    if not doFlipXY:
                        ax.errorbar(percs[1],yAxisValues[i],yerr = None, xerr = errbnds,
                                    fmt = 'ok',ecolor = 'k',elinewidth = 1,
                                    markerfacecolor = parent_colors[j],markeredgecolor = None)
                    else:
                        ax.errorbar(yAxisValues[i],percs[1],yerr = errbnds, xerr = None,
                                    fmt = 'ok',ecolor = 'k',elinewidth = 1,
                                    markerfacecolor = parent_colors[j],markeredgecolor = None)
            else:
                if separate_parents:
                    if not doFlipXY:
                            axs[j].plot(percs[1], yAxisValues[i], 'o',
                                        markerfacecolor = parent_colors[j], markeredgecolor = 'k')
                else:
                    if not doFlipXY:
                        ax.plot(percs[1], yAxisValues[i], 'o',
                                    markerfacecolor = parent_colors[j], markeredgecolor = 'k')

    #Once all the bounds for each parent and child have been calc'd, plot them
    #as a shaded region
    
    #For each parent, plot the variations in mixing coefficients
    if separate_parents:
        if not doFlipXY:
            for j,parent_name in enumerate(parent_names):
                axs[j].fill_betweenx(yAxisValues,
                                 mix_coeff_LB[:,j], mix_coeff_UB[:,j],
                                 color = parent_colors[j],
                                 alpha = plot_alpha)
                axs[j].set_title(parent_name, fontsize='small')
                axs[j].set_xlabel(r'$\phi$')
                if j == 0: # Only do this for the leftmost plot
                    axs[j].set_yticks(yAxisValues)
                    axs[j].set_yticklabels(child_names,rotation = 0)
                else:
                    axs[j].get_yaxis().set_ticks([])
                axs[j].set_xlim(np.min(mix_coeff_LB),np.max(mix_coeff_UB))

        else:
            for j,parent_name in enumerate(parent_names):
                axs[j].fill_between(yAxisValues,
                                 mix_coeff_LB[:,j], mix_coeff_UB[:,j],
                                 color = parent_colors[j],label = parent_name,
                                 alpha = plot_alpha)
                axs[j].set_title(parent_name, fontsize='small')
                axs[j].set_ylabel(r'$\phi$')
                if j == len(parent_names)-1: # Only do this for the bottom plot
                    axs[j].set_xticks(yAxisValues)
                    axs[j].set_xticklabels(child_names,rotation = 45)
                else:
                    axs[j].get_xaxis().set_ticks([])
                axs[j].set_ylim(np.min(mix_coeff_LB),np.max(mix_coeff_UB))

        #axs[j].legend(loc=(1.02, 0))

    else:

        if not doFlipXY:
            for j,parent_name in enumerate(parent_names):
                ax.fill_betweenx(yAxisValues,
                                 mix_coeff_LB[:,j], mix_coeff_UB[:,j],
                                 color = parent_colors[j],label = parent_name,
                                 alpha = plot_alpha)
            ax.set_xlabel(r'$\phi$')
            ax.set_yticks(yAxisValues)
            ax.set_yticklabels(child_names,rotation = 0)
            ax.set_xlim(0,)
        else:
            for j,parent_name in enumerate(parent_names):
                ax.fill_between(yAxisValues,
                                 mix_coeff_LB[:,j], mix_coeff_UB[:,j],
                                 color = parent_colors[j],label = parent_name,
                                 alpha = plot_alpha)
            ax.set_ylabel(r'$\phi$')
            ax.set_xticks(yAxisValues)
            ax.set_xticklabels(child_names,rotation = 45)
            ax.set_xlim(0,)
        ax.legend(loc=(1.02, 0))

    return f, ax

def plot_bootstrapped_mixturecoefficients_yerror(parent_names,
                                                 mix_coeffs_bs,
                                                 confidence_interval = 95.0,
                                                 ax = None,
                                                 colors = None):
    """ Plot for a single child samples mixture model results that visualizes the median
    and percentile bounds of mixing coefficients as a point with error bars. 

    For an example of this plot see Figure 3 from:
    Gilbert, J. Clark, Zane R. Jobe, Samuel A. Johnstone, and Glenn R. Sharman. 2021. “Identifying Elusive
    Piercing Points along the North American Transform Margin Using Mixture Modeling of Detrital Zircon Data
    from Sedimentary Units and Their Crystalline Sources.” The Sedimentary Record 19 (2): 12–21.
    https://doi.org/10.2110/sedred.2021.2.3.

    Parameters
    ----------
    parent_names: list
        List of the names to be used to label parents in the legend.
    mix_coeffs_bs : np.array
        The array of bootstrapped mixing coefficients
    confidence_interval : float (in the interval 0 - 100).
        The percentile of bootstrapped mixing coefficients to summarize. Defaults to 95, infills between the 2.5th and 97.5th percentile.
    ax : matplotlib axis:
        The matplotlib axis to plot on. Defaults to None, in which case an axis is created.
    colors : list
        The list of len(parent_names) of matplotlib colors to color mixing coefficients of each parent by.  Defaults to None, 
        in which case the detritalpy default colors are used.

    Returns
    -------
    ax : matplotlib axis
        The axis on which the plot was created
    """
    
    #If no axis specified, create one:
    if ax is None:
        f,ax = plt.subplots(1,1)
    
    #If no colors were specified, get them from detritalPy
    if (colors is None) or not(len(colors) == len(parent_names)):
        colors = [dFunc.colorMe(i) for i in range(len(parent_names))]


    #Get the upper and lower percentile to report
    alpha = (100.0 - confidence_interval)/2.0

    #Loop through each parent and plot it
    for i in range(len(parent_names)):
        percs = np.percentile(mix_coeffs_bs[:,i],[alpha, 50, 100 - alpha])
        errbnds = np.zeros((2,1))
        errbnds[:,0] = [percs[1] - percs[0],percs[2] - percs[1]]
        ax.errorbar(i+1,percs[1],errbnds,fmt = 'ok',ecolor = 'k',
                          elinewidth = 1,markerfacecolor = colors[i],
                          markeredgecolor = None)
    
    ax.set_xticks(1+np.arange(len(parent_names)))
    ax.set_xticklabels(parent_names,rotation = 45)
    
    return ax

def plot_many_bootstrapped_mixturecoefficients_yerror(parent_names,
                                                      child_names,
                                                 mix_coeffs_set_bs,
                                                 confidence_interval = 95.0,
                                                 axs = None,
                                                 plotWidth = 4.0,
                                                 plotHeight = 8.0,
                                                 colors = None):
    """ Plot for a set of child samples mixture model results that visualizes the median
    and percentile bounds of mixing coefficients as a point with error bars. Each child's results
    are plotted on a new axis.

    For an example of this plot see Figure 3 from:
    Gilbert, J. Clark, Zane R. Jobe, Samuel A. Johnstone, and Glenn R. Sharman. 2021. “Identifying Elusive
    Piercing Points along the North American Transform Margin Using Mixture Modeling of Detrital Zircon Data
    from Sedimentary Units and Their Crystalline Sources.” The Sedimentary Record 19 (2): 12–21.
    https://doi.org/10.2110/sedred.2021.2.3.

    Args:
        parent_names: list
            List of the names to be used to label parents in the legend.
        child_names : list
            List of the names to be used to label the children along the y axis
        mix_coeffs_bs_set : list
            The list of bootstrapped mixing coefficients returned by bootstrap_solve_mixture_models
        confidence_interval : float (in the interval 0 - 100).
            The percentile of bootstrapped mixing coefficients to summarize. Defaults to 95, infills between the 2.5th and 97.5th percentile.
        axs : matplotlib axes:
            The matplotlib axes to plot on. Must be the same length as len(child_names) and len(mix_coefs_set_bs).
            Defaults to None, in which case an axis is created.
        plotWidth : float, optional
            The width of the plot. Default = 4.0
        plotHeight : float, optional
            The height of the plot. Default = 8.0
        colors : list
            The list of len(parent_names) of matplotlib colors to color mixing coefficients of each parent by.  Defaults to None, 
            in which case the detritalpy default colors are used.

    Returns:
        axs : matplotlib axes
            The matplotlib axes on which each child's plot was created
    """

    
    #If a suite of axis were provided, make them
    #if (axs is None) or not(len(axs) == len(child_names)):
    if axs is None:
        f,axs = plt.subplots(len(child_names),1, figsize=(plotWidth, plotHeight))
        
    
    for i, (mix_coeffs_bs,child_name) in enumerate(zip(mix_coeffs_set_bs,
                                                       child_names)):
        
        plot_bootstrapped_mixturecoefficients_yerror(parent_names,
                                                  mix_coeffs_bs,
                                                  confidence_interval,
                                                  axs[i],
                                                  colors)
        
        axs[i].set_title(child_name)
    
    return axs

def plot_bootstrapped_mixturecoefficients_violin(parent_names,
                                                 mix_coeffs_bs,
                                                 ax = None,
                                                 colors = None,
                                                 **violPlots_kwargs):
    
    """ Plot for a single child samples mixture model results that visualizes the distribution
    of mixing coefficients as a violin plots. 

    For an example of this plot see Figure 3 from:
    
    Gilbert, J. Clark, Zane R. Jobe, Samuel A. Johnstone, and Glenn R. Sharman. 2021. “Identifying Elusive
    Piercing Points along the North American Transform Margin Using Mixture Modeling of Detrital Zircon Data
    from Sedimentary Units and Their Crystalline Sources.” The Sedimentary Record 19 (2): 12–21.
    https://doi.org/10.2110/sedred.2021.2.3.

    Args:
        parent_names: list
            List of the names to be used to label parents in the legend.
        mix_coeffs_bs : np.array
            The array of bootstrapped mixing coefficients
        ax : matplotlib axis:
            The matplotlib axis to plot on. Defaults to None, in which case an axis is created.
        colors : list
            The list of len(parent_names) of matplotlib colors to color mixing coefficients of each parent by.  Defaults to None, 
            in which case the detritalpy default colors are used.

    Returns:
        ax : matplotlib axis
            The axis on which the plot was created
    """

    #If no axis specified, create one:
    if ax is None:
        f,ax = plt.subplots(1,1)
    
    #If no colors were specified, get them from detritalPy
    if (colors is None) or not(len(colors) == len(parent_names)):
        colors = [dFunc.colorMe(i) for i in range(len(parent_names))]
        

    #Plot the histograms of the mixing coefficients and the value of the mixture function
    violinParts = ax.violinplot(mix_coeffs_bs,**violPlots_kwargs)
    
    #Color the violin plots
    for i,pc in enumerate(violinParts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_edgecolor(colors[i])
        pc.set_alpha(1)

    ax.set_xticks(1+np.arange(len(parent_names)))
    ax.set_xticklabels(parent_names,rotation = 45)

    return ax

def plot_many_bootstrapped_mixturecoefficients_violin(parent_names,
                                                      child_names,
                                                 mix_coeffs_set_bs,
                                                 axs = None,
                                                 colors = None,
                                                 **violPlots_kwargs):
    
    """ Plot for a suite of child samples mixture model results that visualizes the distribution
    of mixing coefficients as a violin plots. 

    For an example of this plot see Figure 3 from:
    
    Gilbert, J. Clark, Zane R. Jobe, Samuel A. Johnstone, and Glenn R. Sharman. 2021. “Identifying Elusive
    Piercing Points along the North American Transform Margin Using Mixture Modeling of Detrital Zircon Data
    from Sedimentary Units and Their Crystalline Sources.” The Sedimentary Record 19 (2): 12–21.
    https://doi.org/10.2110/sedred.2021.2.3.

    Args:
        parent_names: list
            List of the names to be used to label parents in the legend.
        mix_coeffs_bs : np.array
            The array of bootstrapped mixing coefficients
        ax : matplotlib axis:
            The matplotlib axis to plot on. Defaults to None, in which case an axis is created.
        colors : list
            The list of len(parent_names) of matplotlib colors to color mixing coefficients of each parent by.  Defaults to None, 
            in which case the detritalpy default colors are used.

    Returns:
        ax : matplotlib axis
            The axis on which the plot was created
    """

    #If a suite of axis were provided, make them
    #if (axs is None) or not(len(axs) == len(child_names)):
    if axs is None: 
        f,axs = plt.subplots(len(child_names),1)
        
    for i, (mix_coeffs_bs,child_name) in enumerate(zip(mix_coeffs_set_bs,
                                                       child_names)):
        
        plot_bootstrapped_mixturecoefficients_violin(parent_names,
                                                 mix_coeffs_bs,
                                                 ax = axs[i],
                                                 colors = colors,
                                                 **violPlots_kwargs)
        
        axs[i].set_title(child_name)
    
    return axs

def plot_bootstrapped_mixturecoefficients_histograms(parent_names,
                                                     mix_coeffs_bs,
                                                     confidence_interval = 95.,
                                                     ax = None,
                                                     colors = None,
                                                     **hist_kwargs):

    """ Plot for a single child samples mixture model results that visualizes the distribution
    of mixing coefficients as a suit of histograms. 

    For an example of this plot see Figure 7 C from:
    
    Malkowski, M. A., et al. "Continental shelves as detrital mixers: U–Pb and Lu–Hf detrital zircon provenance of the
    Pleistocene–Holocene Bering Sea and its margins." The Depositional Record (2022).

    Args:
        parent_names: list
            List of the names to be used to label parents in the legend.
        mix_coeffs_bs : np.array
            The array of bootstrapped mixing coefficients
        confidence_interval: float (0 - 100)
            The percentile bound of bootstrapped results to summarize in the legend for each mixing coefficient.
        ax : matplotlib axis:
            The matplotlib axis to plot on. Defaults to None, in which case an axis is created.
        colors : list
            The list of len(parent_names) of matplotlib colors to color mixing coefficients of each parent by.  Defaults to None, 
            in which case the detritalpy default colors are used.
        **hist_kwargs :
            Additional keyword arguments to pass to matplotlib.pyplot.hist function . e.g. bins = 100 (to use 100 bins) or histtype = 'stepfilled'.

    Returns:
        ax : matplotlib axis
            The axis on which the plot was created
    """
    
    #If no axis specified, create one:
    if ax is None:
        f,ax = plt.subplots(1,1)
    
    #If no colors were specified, get them from detritalPy
    if (colors is None) or not(len(colors) == len(parent_names)):
        colors = [dFunc.colorMe(i) for i in range(len(parent_names))]
        
    #Get the upper and lower percentile to report
    alpha = (100.0 - confidence_interval)/2.0
    
    #For each of the parents, plot that histogram
    for i,parent_name in enumerate(parent_names):
        
        prcnt = np.percentile(mix_coeffs_bs[:,i],[alpha, 100-alpha])
        
        thisLabel = parent_name + '({:.1f}% - {:.1f}%)'.format(prcnt[0]*100.0,
                                                               prcnt[1]*100.0)
        
        ax.hist(mix_coeffs_bs[:,i],color = colors[i],label = thisLabel,
                 **hist_kwargs)
        
        ax.set_xlabel(r'$\Phi$')
        
    return ax
  
def plot_many_bootstrapped_mixturecoefficients_histograms(parent_names,
                                                          child_names,
                                                     mix_coeffs_set_bs,
                                                     confidence_interval = 95.,
                                                     axs = None,
                                                     colors = None,
                                                     **hist_kwargs):
    """ Plot for a suite of child samples mixture model results that visualizes the distribution
    of mixing coefficients as a suit of histograms. 

    For an example of this plot see Figure 7 C from:
    
    Malkowski, M. A., et al. "Continental shelves as detrital mixers: U–Pb and Lu–Hf detrital zircon provenance of the
    Pleistocene–Holocene Bering Sea and its margins." The Depositional Record (2022).

    Args:
        parent_names: list
            List of the names to be used to label parents in the legend.
        child_names : list
            List of the names to be used to label the children along the y axis
        mix_coeffs_bs : np.array
            The array of bootstrapped mixing coefficients
        confidence_interval: float (0 - 100)
            The percentile bound of bootstrapped results to summarize in the legend for each mixing coefficient.
        ax : matplotlib axis:
            The matplotlib axis to plot on. Defaults to None, in which case an axis is created.
        colors : list
            The list of len(parent_names) of matplotlib colors to color mixing coefficients of each parent by.  Defaults to None, 
            in which case the detritalpy default colors are used.
        **hist_kwargs :
            Additional keyword arguments to pass to matplotlib.pyplot.hist function . e.g. bins = 100 (to use 100 bins) or histtype = 'stepfilled'.

    Returns:
        ax : matplotlib axis
            The axis on which the plot was created
    """
    
    #If a suite of axis were provided, make them
#    if (axs is None) or not(len(axs) == len(child_names)):
    if axs is None:
        f,axs = plt.subplots(len(child_names),1)
        
    for i, (mix_coeffs_bs,child_name) in enumerate(zip(mix_coeffs_set_bs,
                                                       child_names)):
        
        plot_bootstrapped_mixturecoefficients_histograms(parent_names,
                                                         mix_coeffs_bs,
                                                         confidence_interval,
                                                         axs[i],
                                                         colors,
                                                         **hist_kwargs)
        
        axs[i].set_title(child_name)
    
    return axs

def plot_bootstrapped_distribution_bounds(distribution_x_axis,
                                          bootstrapped_distributions,
                                          ax = None,
                                          fill_color = 'k',fill_alpha = 0.5,
                                          confidence_interval = 95.0,
                                          label = ''):
    '''
    Plot summarizing the variability in bootstrapped distributions that examines the range in probabilities
    observed as a function of age. Plots a shaded region that bounds some percentile of bootstrapped distributions for a single chid sample.

    For an example of this plot see Figure 7 A from:
    
    Malkowski, M. A., et al. "Continental shelves as detrital mixers: U–Pb and Lu–Hf detrital zircon provenance of the
    Pleistocene–Holocene Bering Sea and its margins." The Depositional Record (2022).

    Parameters
    ----------
    distribution_x_axis : np.array
        The x axis (e.g., age) that the distributions were evaluated on.
    bootstrapped_distributions : np.ndarray 
        Array containing the distribution calculated for each bootstrap iteration.
    ax : matplotlib axis:
        The matplotlib axis to plot on. Defaults to None, in which case an axis is created.
    fill_color : string
        Matplotlib color string used to shade within the bounds of distributions. The default is 'k'.
    fill_alpha : float (0 - 1), optional
        The transparency level to assign for the shaded fill. The default is 0.5.
    confidence_interval: float (0 - 100)
        The percentile bound of bootstrapped results to summarize in the legend for each mixing coefficient.
    label : string, optional
        The label for this data in the plot legend. The default is ''.

    Returns
    -------
        ax : matplotlib axis
            The axis on which the plot was created

    '''
    
    #What percentile range do we want?
    alpha = (100.0 - confidence_interval)/2.0
    uprLwrBounds = np.percentile(bootstrapped_distributions,
                                 [alpha, 100.0-alpha],
                                 axis = 0)
    
    label+= ' ({:.1f}% - {:.1f}%)'.format(alpha,100-alpha)

    ax.fill_between(distribution_x_axis,uprLwrBounds[0],
                    uprLwrBounds[1],color = fill_color, alpha = fill_alpha,
                    label = label)
    
    ax.set_xlabel('Age (Ma)')
    
    
    return ax

def plot_many_bootstrapped_distribution_bounds(distribution_x_axis,
                                               sample_names,
                                               bootstrapped_distributions_set,
                                               axs = None,
                                               fill_colors = None,
                                               fill_alpha = 0.5,
                                               confidence_interval = 95.0):
    '''
    Plot summarizing the variability in bootstrapped distributions that examines the range in probabilities
    observed as a function of age. Plots a shaded region that bounds some percentile of bootstrapped distributions for a suite
    of samples.

    For an example of this plot see Figure 7 A from:
    
    Malkowski, M. A., et al. "Continental shelves as detrital mixers: U–Pb and Lu–Hf detrital zircon provenance of the
    Pleistocene–Holocene Bering Sea and its margins." The Depositional Record (2022).

    Parameters
    ----------
    distribution_x_axis : np.array
        The x axis (e.g., age) that the distributions were evaluated on.
    sample_names : list
        The list of names that will be used as labels in each legend.
    bootstrapped_distributions_set : list
        A list of arrays (one for each sample name) containing the distribution calculated for each bootstrap iteration.
    axs : matplotlib axes:
        An array of matplotlib axis to plot on. Defaults to None, in which case an axis is created.
    fill_color : string
        Matplotlib color string used to shade within the bounds of distributions. The default is 'k'.
    fill_alpha : float (0 - 1), optional
        The transparency level to assign for the shaded fill. The default is 0.5.
    confidence_interval: float (0 - 100)
        The percentile bound of bootstrapped results to summarize in the legend for each mixing coefficient.
    label : string, optional
        The label for this data in the plot legend. The default is ''.

    Returns
    -------
        ax : matplotlib axis
            The axis on which the plot was created

    '''

    #If a suite of axes were not provided, make them
    #if (axs is None) or not(len(axs) == len(sample_names)):
    if axs is None:
        f,axs = plt.subplots(len(sample_names),1)
        
    #If no fill colors were provided make some up
    if (fill_colors is None) or not(len(fill_colors) == len(sample_names)):
        fill_colors = [dFunc.colorMe(i) for i in range(len(sample_names))]
    
    if len(sample_names) == 1:
        plot_bootstrapped_distribution_bounds(distribution_x_axis,
                                              bootstrapped_distributions_set[0],
                                              ax = axs,
                                              fill_color = fill_colors[0],
                                              fill_alpha = 0.5,
                                              confidence_interval = 95.0,
                                              label = sample_names[0])

    else:
        for i, (sample_name,bootstrapped_distributions) in enumerate(zip(sample_names,bootstrapped_distributions_set)):
            plot_bootstrapped_distribution_bounds(distribution_x_axis,
                                                  bootstrapped_distributions,
                                                  ax = axs[i],
                                                  fill_color = fill_colors[i],
                                                  fill_alpha = 0.5,
                                                  confidence_interval = 95.0,
                                                  label = sample_name)
        
    return axs

def plot_child_bootstrappedmodel_distribution_comparison(main_byid_df,child_modelled_distributions,
                                                            child_list, xaxis_1=0, xaxis_2=4500, x1=0, x2=4500, xdif=1, bw=2, bw_x=None, objective_metric='dmax',
                                                            confidence_interval=95.0, fill_alpha = 0.5,
                                                            plot_self_comparisons = False,
                                                            childDists_bs_set = None,
                                                            axs = None, plotWidth = 4.0, subplotHeight = 3.0, child_colors = None, model_color = 'k',
                                                            sigma = '1sigma'):
    """
    
    Creates a plot to help visualize comparisons between the predictions of mixture modelling and the observed child distribution for a set of 
    modelled children. The 'best' fit mixture model is not necessarily a good fit; the modelled distributions may not reproduce important age
    modes or proportions observed in the child because the parents couldn't be combined to do that.

    Plot shows the percentile bounds of bootstrapped, modelled child distributions and the observed child distribution. The type of distribution thats
    plotted is controlled by the type of distribution used by the specified objective_metric.

    For an example of this plot see Figure 7 A from:
    
    Malkowski, M. A., et al. "Continental shelves as detrital mixers: U–Pb and Lu–Hf detrital zircon provenance of the
    Pleistocene–Holocene Bering Sea and its margins." The Depositional Record (2022).

    Args:

        main_byid_df : dataframe:
            dataframe loaded using detritalPy, e.g., with the loadDataExcel function
        child_modelled_distributions : list
            A list of arrays (one for each child name) containing the distribution calculated for each bootstrap iteration performed during mixture modelling.
        child_list : list
            List of the samples to use as children for mixture modelling.
        xaxis_1 : float, optional
            Minimum value of age axis to plot. The default is 0.
        xaxis_2 : float, optional
            Maximum value of age axis to plot. The default is 4500.
        xdif : float, optional
            Spacing of age axis distributions are evaluated along. The default is 1.
        bw : float, optional
            bandwidth used when calculating KDE. The default is 2.
        bw_x : None or list, optional
            Set to None if not using a split KDE bw. Otherwise, set x-axis locations for bw split (in Ma) in a list, (e.g., bw_x=[300])
        objective_metric : string, optional
            Name of objective_metric that will be used to compare distributions. The default is 'dmax'. For a 
            full list of values call print(dMix.AVAILABLE_METRICS)
        confidence_interval: float (0 - 100), optional
            The percentile bound of bootstrapped results to summarize in the legend for each mixing coefficient. Default = 95.0
        fill_alpha : float (0 - 1), optional
            The transparency level to assign for the shaded fill. The default is 0.5.
        self_comparison_bootstrapped_children : list, optional
             A list of the distributions of the child created by resampling it with replacement, created by bootstrapped_self_comparisons_many_samples or boostrapped_self_comparison.
             Defaults to None, in which only the observed child distribution is plotted.
        axs : matplotlib axes:
            An array of matplotlib axis to plot on. Defaults to None, in which case an axis is created
        plotWidth : float, optional
            Width of the plot. Default = 4.0
        subplotHeight : float, optional
            Width of each subplot. Default = 3.0
        child_colors  : str or list, optional
            The list of len(child_list) of matplotlib colors to color observed child distributions by.  Defaults to 'Default', 
            in which case the detritalpy default colors are used. May also specify a single color as a string (e.g., 'red')
        model_color : (str, optional)
            A matplotlib color to color the modelled distribution range by. Defaults to 'black'.
        sigma : string, optional
            Uncertainly level of input data (options: '1sigma' or '2sigma'). The default is '1sigma'.

    Returns:
        axs : matplotlib axes:
            An array of matplotlib axis that the plots were created on
    """

    #if (axs is None) or (len(axs) != len(child_list)):
    #    f,axs = plt.subplots(len(child_list),1,sharex = True,sharey = 'col')

    if axs is None:
        f,axs = plt.subplots(len(child_list),1,figsize = (plotWidth,subplotHeight*len(child_list)),sharex = True,sharey = 'col')

    if child_colors == 'Default':
        child_colors = [dFunc.colorMe(i) for i in range(len(child_list))]

    if type(child_colors) == str:
        child_colors = [child_colors] * len(child_list)

    #First, get the parent distributions to plot
    #This function gets the function that was used to generate distributions to mix (based on whats needed for the requested objective metric)
    distribution_function = lookup_functions_for_metric(objective_metric,x1=x1, x2=x2, xdif=xdif, bw=bw, bw_x=bw_x)[0]
 
    #Get the data for the children
    ages_c,errors_c,nGrains_c,labels_c = dFunc.sampleToData(child_list,
                                                        main_byid_df,
                                                        sigma = sigma)
    #Get the distributions of the children
    dist_axis,distributions_children = distribution_function(ages_c,errors_c)


    #Plot the bootstrapped model distributions
    plot_many_bootstrapped_distribution_bounds(dist_axis,[r'Model' for c in child_list],child_modelled_distributions,axs,
                                                fill_colors = ['k' for c in child_list],
                                                fill_alpha = fill_alpha, confidence_interval = confidence_interval)

    #Plot the bootstrapped child distributions
    if (plot_self_comparisons) and not(childDists_bs_set is None) and (len(childDists_bs_set) == len(child_list)):
        plot_many_bootstrapped_distribution_bounds(dist_axis,labels_c,childDists_bs_set,axs,
                                                    fill_colors = child_colors,fill_alpha = fill_alpha,
                                                    confidence_interval = confidence_interval)

    #Add the legends to the plots
    if len(child_list) == 1:
        axs.plot(dist_axis,distributions_children[0],'-',color = 'k',label = labels_c[0]) # color = child_colors[0]
        axs.legend(loc = 'best',fontsize = 'small')
        #axs.spines['right'].set_visible(False)
        #axs.spines['top'].set_visible(False)
        axs.set_xlim(xaxis_1, xaxis_2)
        #if i < (len(child_list)-1):
        #    axs.set_xlabel('')
        if (objective_metric == 'dmax') or (objective_metric == 'vmax'): # CDFs are plotted for Dmax and Vmax, relative age distributions for the others
            axs.set_ylim(0,1)
        else:
            axs.set_ylim(0,)

    else:
        for i in range(len(child_list)):
            axs[i].plot(dist_axis,distributions_children[i],'-', color = 'k', label = labels_c[i]) # color = child_colors[i]
            axs[i].legend(loc = 'best',fontsize = 'small')
            axs[i].set_xlim(xaxis_1, xaxis_2)
            if i < (len(child_list)-1):
                axs[i].set_xlabel('')
        if (objective_metric == 'dmax') or (objective_metric == 'vmax'): # CDFs are plotted for Dmax and Vmax, relative age distributions for the others
            axs[i].set_ylim(0,1)
        else:
            axs[i].set_ylim(0,)

    plt.subplots_adjust(hspace=0.0)
    plt.tight_layout()


    return axs


def plot_bootstrapped_metric_comparisons_model_observations(bs_modelled_metrics,
                                                            bs_self_compared_metrics,
                                                            objective_metric,
                                                            child_name,
                                                            obj_func_crit_i,
                                                            worse_than_crit_i,
                                                            obj_func_val_i,
                                                            modelled_color = 'k',
                                                            self_compared_color = 'r',
                                                            ax = None,
                                                            doAddSummaryTitle = True,
                                                            x_min=0,
                                                            x_max=1,
                                                            **histkwargs):
    """
    Creates a plot to help assess the model goodness of fit. This function creates plots that show histograms of the values of the 
    objective_metric computed during bootstrapped mixture modelling and values of the objective_metric observed when comparing the child
    to a resampled version of itself.  If the mixture modelling does a good job reproducing the observed child sample, then these distributions
    should have a large amount of overlap.

    For an example of this plot see Figure 7 B from:
    
    Malkowski, M. A., et al. "Continental shelves as detrital mixers: U–Pb and Lu–Hf detrital zircon provenance of the
    Pleistocene–Holocene Bering Sea and its margins." The Depositional Record (2022).

    Args:
        bs_modelled_metrics : np.array
            Array of the objective_metric values created by bootstrapped mixture modelling
        bs_self_compared_metrics : np.array
            Array of the objective_metric values created by bootstrapped comparisons of the child to a resampled version
            of itself.
        objective_metric : string, optional
            Name of objective_metric that will be used to compare distributions. The default is 'Dmax'. For a 
            full list of values call print(dMix.AVAILABLE_METRICS)
        child_name : string
            The name of the child
        modelled_color : (str, optional)
            Matplotlib color for coloring the histogram of model results. Defaults to 'k'.
        self_compared_color (str, optional)
            Matplotlib color for coloring the histogram of self comparison results. Defaults to 'r'.
        ax : matplotlib axis
            Axis to create the plot on. Defaults to None, in which case a new plot and axis are created.
        doAddSummaryTitle : bool
            Should we add a title summarizing the overlap in the distributions? Defaults to True.
        **hist_kwargs :
            Additional keyword arguments to pass to matplotlib.pyplot.hist function . e.g. bins = 100 (to use 100 bins) or histtype = 'stepfilled'.

    Returns:
        ax : matplotlib axis
            The axis the plot was created on.
    """
    if ax is None:
        f,ax = plt.subplots(1,1)
    
    #How does this objective metric scale?    
    areGoodFitsSmallValues = lookup_functions_for_metric(objective_metric)[-1]
    
    #How many iterations were preformed?
    nBootstrapIterations = np.min([len(bs_self_compared_metrics), len(bs_self_compared_metrics)])


    #Plot a histogram of the modelled comparisons
    ax.hist(bs_modelled_metrics,color = modelled_color, range=(x_min, x_max),
            label = 'Bootstrapped model comparisons',
            **histkwargs)


    #Plot a histogram of the self-compared values
    ax.hist(bs_self_compared_metrics, color = self_compared_color, range=(x_min, x_max),
            label = 'Bootstrapped self-comparisons',
            **histkwargs)

    ax.axvline(x=obj_func_crit_i, ymin=0, ymax=1, ls='--', color='red', label=objective_metric+'_crit')

    y1, y2 = ax.get_ylim()
    x1, x2 = ax.set_xlim()

    ax.arrow(x=obj_func_val_i, y=(y2-y1)*0.1, dx=0, dy=-(y2-y1)*0.1, head_width = 0.03*(x2-x1), head_length = 0.04*(y2-y1), length_includes_head=True, head_starts_at_zero=False, label='Best-fit', edgecolor='black', facecolor='white')

    if areGoodFitsSmallValues:
        ax.arrow(x=obj_func_crit_i, y=(y2-y1)/2, dx=-obj_func_crit_i*0.05, dy=0.0, head_width = 0.04*(y2-y1), head_length = 0.03*(x2-x1), length_includes_head=False, head_starts_at_zero=False, label='Good model fits are in this direction', edgecolor='black', facecolor='red')
    else:
        ax.arrow(x=obj_func_crit_i, y=(y2-y1)/2, dx=obj_func_crit_i*0.05, dy=0.0, head_width = 0.04*(y2-y1), head_length = 0.03*(x2-x1), length_includes_head=False, head_starts_at_zero=False, label='Good model fits are in this direction', edgecolor='black', facecolor='red')

    ax.legend(loc=('best'), fontsize='x-small')
    
    ax.set_xlabel(objective_metric)

    if doAddSummaryTitle:
        '''
        Lets also try and summarize this result, the key question is: how often
        are the resampled, self-comparisons worse than whats fitted by the models?
        '''
        
        if areGoodFitsSmallValues:
            #fraction_worse = np.sum(bs_self_compared_metrics[:nBootstrapIterations] > bs_modelled_metrics[:nBootstrapIterations])/nBootstrapIterations
            comparison_name = 'greater'
        else:
            #fraction_worse = np.sum(bs_self_compared_metrics[:nBootstrapIterations] < bs_modelled_metrics[:nBootstrapIterations])/nBootstrapIterations
            comparison_name = 'smaller'

        ax.set_title(child_name + '\n'+ '{:.1f}% of resampled observations have '.format(100*worse_than_crit_i)
                   + comparison_name + '\n' +
                  ' values of ' + objective_metric + ' (are worse than) than ' + objective_metric +'_crit',
                  fontsize = 'small')

    else:
        ax.set_title(child_name)

    return ax

def plot_many_bootstrapped_metric_comparisons_model_observations(bs_modelled_metrics_set,
                                                            bs_self_compared_metrics_set,
                                                            objective_metric,
                                                            child_list,
                                                            main_byid_df,
                                                            obj_func_crit,
                                                            worse_than_crit,
                                                            obj_func_val,
                                                            modelled_colors = None,
                                                            self_compared_colors = None,
                                                            axs = None,
                                                            plotWidth=4.0,
                                                            subplotHeight=3.0,
                                                            doAddSummaryTitle = True,
                                                            **histkwargs):
    
    """
    Creates a set of plots to help assess the model goodness of fit. This function creates plots for a suite of children that show histograms of the values of the 
    objective_metric computed during bootstrapped mixture modelling and values of the objective_metric observed when comparing the child
    to a resampled version of itself.  If the mixture modelling does a good job reproducing the observed child sample, then these distributions
    should have a large amount of overlap.

    For an example of one of these plot see Figure 7 B from:
    
    Malkowski, M. A., et al. "Continental shelves as detrital mixers: U–Pb and Lu–Hf detrital zircon provenance of the
    Pleistocene–Holocene Bering Sea and its margins." The Depositional Record (2022).

    Args:
        bs_modelled_metrics_set : list
            A list of arrays of the objective_metric values created by bootstrapped mixture modelling of children.
        bs_self_compared_metrics_set : np.array
            A list of arrays of the objective_metric values created by bootstrapped comparisons of the children to a resampled version
            of themselves.
        objective_metric : string, optional
            Name of objective_metric that will be used to compare distributions. The default is 'Dmax'. For a 
            full list of values call print(dMix.AVAILABLE_METRICS)
        child_list : list
            The list of child samples, either as a list or as a list of tuples
        modelled_colors : (list, optional)
           A list of matplotlib colors (of len(bs_modelled_metrics)) for coloring the histogram of model results. Defaults to None, in which case
           'k' is used for all examples.
        self_compared_colors : (str, optional)
            A list of matplotlib colors (of len(bs_modelled_metrics)) for coloring the histogram of model results. Defaults to None, in which case
           the default detritalpy colors are used.
        axs : matplotlib axes
            The matplotlib axes to create the plots on. Defaults to None, in which case a new plot and set of axes are created.
        doAddSummaryTitle : bool
            Should we add a title summarizing the overlap in the distributions to each plot? Defaults to True.

    Returns:
        axs : matplotlib axss
            The axss the plots were created on.
    """
    
    child_names = dFunc.sampleToData(child_list, main_byid_df)[3]

    x_min = np.min((np.min(obj_func_val),np.min(bs_modelled_metrics_set),np.min(bs_self_compared_metrics_set)))
    x_max = np.max((np.max(obj_func_val),np.max(bs_modelled_metrics_set),np.max(bs_self_compared_metrics_set)))

    if axs is None:
        f,axs = plt.subplots(len(bs_modelled_metrics_set),1,figsize = (plotWidth,subplotHeight*1),sharex = False);

    if modelled_colors is None:
        modelled_colors = ['k' for i in bs_modelled_metrics_set]

    if type(modelled_colors) == str:
        modelled_colors = [modelled_colors for i in bs_modelled_metrics_set]
        
    if self_compared_colors is None:
        self_compared_colors = [dFunc.colorMe(i) for i in range(len(bs_modelled_metrics_set))]
    
    if type(self_compared_colors) == str:
        self_compared_colors = [self_compared_colors for i in bs_modelled_metrics_set]

    if len(child_names) == 1: # If only plotting one child sample
        for i, (bs_modelled_metrics,bs_self_compared_metrics) in enumerate(zip(bs_modelled_metrics_set,bs_self_compared_metrics_set)):
            plot_bootstrapped_metric_comparisons_model_observations(bs_modelled_metrics,
                                                                bs_self_compared_metrics,
                                                                objective_metric, child_names[i],
                                                                obj_func_crit[i], worse_than_crit[i],
                                                                obj_func_val[i],
                                                                modelled_color = modelled_colors[i],
                                                                self_compared_color = self_compared_colors[i],
                                                                ax = axs,
                                                                doAddSummaryTitle = doAddSummaryTitle, x_min=x_min,
                                                                x_max=x_max,
                                                                **histkwargs)
        axs.set_xlim(x_min, x_max)
    else:
        for i, (bs_modelled_metrics,bs_self_compared_metrics) in enumerate(zip(bs_modelled_metrics_set,bs_self_compared_metrics_set)):
            plot_bootstrapped_metric_comparisons_model_observations(bs_modelled_metrics,
                                                                bs_self_compared_metrics,
                                                                objective_metric, child_names[i],
                                                                obj_func_crit[i], worse_than_crit[i],
                                                                obj_func_val[i],
                                                                modelled_color = modelled_colors[i],
                                                                self_compared_color = self_compared_colors[i],
                                                                ax = axs[i],
                                                                doAddSummaryTitle = doAddSummaryTitle, x_min=x_min,
                                                                x_max=x_max,
                                                                **histkwargs)
        
        if not(i+1 == len(modelled_colors)):
            axs[i].set_xlabel('')
        axs[i].tick_params(labelbottom=True)
        axs[i].legend(fontsize = 'x-small')
        axs[i].set_xlim(x_min, x_max)
    plt.tight_layout()
  
    
    return axs
    
def plotMix(main_byid_df, parent_list, child_list,
            plotType='KDE', bw=2, bw_x=None, x1=0, x2=4500, xdif=1,
            fillParent=True, parent_colors='Default', child_colors='Default', 
            color_by_age=False, agebins=None, agebinsc=None, agebinsc_alpha=1, 
            xaxis_1=0, xaxis_2=4000, w1=6, w2=4, c=4, plotPie = True, plotMixResults=False,
            best_plot_type = 'best-fit', plotResultType='line', violin_width=0.15, line_style=None,
            CDFlw_p=2, CDFlw_c=1, KDElw=1, PDPlw=1, 
           sigma='1sigma', mix_coeffs_all=None, mix_coeffs_bf=None, obj_func_val=None, best_mixed_dist=None,
           obj_vals_all=None, objective_metric=None):
    
    '''
    Creates a plot that summarizes mixture modeling results. Both children and parent samples are shown as a cumulative plot (CDF)
    above and as relative plots (KDE or PDP) below. Parent mixtures for each child sample are shown to the right of the corresponding sample.

    Parameters
    ----------
    main_byid_df : dataframe
        dataframe loaded using detritalPy, e.g., with the loadDataExcel function
    parent_list : list or tuple
        List or tuple of the samples to use as parents for mixture modelling.
    child_list : list or tuple
        List or tuple of the samples to use as children for mixture modelling..
    plotType : str, optional
        Specify the type of plot to make: options are 'KDE' or 'PDP'. Default = 'KDE'
    bw : float or string, optional
        bandwidth in Ma used when calculating KDE. Default = 2
    bw_x : float, optional
         X-axis location of bandwidth split (only used if multiple KDE values are specified). The default is None.
    xdif : float, optional
        Spacing of age axis distributions are evaluated along. The default is 1.
    fillParent : Boolean, optional
        Color the area under the parent KDE or PDE. Default = True
    parent_colors : str or list, optional
        The list of len(parent_list) of matplotlib colors to color mixing coefficients of each parent by.
        May also use a single color (e.g., 'black'). Default value is 'Default'
        in which case the detritalpy default colors are used.
    child_colors : str or list, optional
        The list of len(child_list) of matplotlib colors to color mixing coefficients of each child by.
        May also use a single color (e.g., 'black'). Default value is 'Default'
        in which case the detritalpy default colors are used.
    color_by_age : Boolean, optional
        Will color KDE or PDP according to age categories if set to True. Default = False
    plotLog : Boolean, optional
        Set to True to plot the x-axis on a log scale. Default = False.
    agebins : list, optional
        Array of bin edges in Myr. Format option 1: [age1, age2, age3, etc.]. Format option 2: [[bin1_min, bin1_max],[bin2_min, bin2_max],etc.]. Default = None.
    agebinsc : list, optional
        Array of colors that correspond to age bins, default = None.
    agebinsc_alpha : float (0 - 1), optional
        Alpha value for agebinsc color fill. Default = 1.
    xaxis_1 : float, optional
        Minimum x-axis value in Ma for plotting. Default = 0
    xaxis_2 : float, optional
        Maximum x-axis value in Ma for plotting. Default = 4500
    w1 : float, optional
        Width of the KDE or PDP plots. Default = 6
    w2 : float, optional
        Width of the mixture modeling plots. Default = 4
    c : float, optional
        Height of the CDF panel. Default = 4
    plotPie : Boolean, optional
        Set to True to plot mixture modeling results (best-fit or average) as pie diagrams. Default = True.
    plotMixResults : Boolean, optional
        Set to True to plot mixture modeling results to the right of the KDE or PDE plots. Default = True.
    best_plot_type : string, optional
        Specify what value should be plotted for the mixture. Options: 'best-fit', 'average', 'median'. Default = 'best-fit'
    plotResultType : string, optional
        Specify how the parent mixtures should be displayed. Options: 'line' or 'violin'. Default = 'line'
    violin_width : float, optional
        Width of violins if plotResultType = 'violin'. Default = 0.15
    line_style : list or None, optional
        List of matplotlib line formats. Default value used if None. Default = None.
    CDFlw_p : float, optional
        Line weight of the parent CDFs. Default = 2
    CDFlw_c : float, optional
        Line weight of the child CDFs. Default = 1
    KDElw : float, optional
        Line weight of the KDEs. Default = 1
    PDPlw : float, optional
        Line weight of the PDPs. Default = 1
    sigma : string, optional
        Specify whether errors are '1sigma' or '2sigma' errors. Default is '1sigma'.
    mix_coeffs_all : list, optional
        The suite of bootstrapped mixing coefficients for each children. The default is None.
    mix_coeffs_bf : list, optional
        A list of arrays, each array contains the best fit mixing coefficients. The default is None.
    obj_func_val : list, optional
        A list of the values of the objective_metric calculated with the best fit mixing coefficients. The default is None.
    best_mixed_dist : list, optional
        A list of arrays, each array contains the best fitting mixed distribution. The default is None.
    obj_vals_all : list of arrays, optional
        A list of len(child_list) containing arrays of the values of the objective_metric calculated during bootstrapped mixture
        modelling. Default is None.
    objective_metric : string, optional
        Name of objective_metric that was used to compare distributions. The default is None. For a 
        full list of values call print(dMix.AVAILABLE_METRICS)

    Returns
    -------
    Figure
    '''

    if parent_colors == 'Default':
        parent_colors = [dFunc.colorMe(x) for x in np.arange(len(parent_list))]

    if child_colors == 'Default':
        child_colors = ['black'] * len(child_list)

    if type(child_colors) == str:
        child_colors = [child_colors] * len(child_list)
    
    # Get the data for the parents
    ages_p,errors_p,nGrains_p,labels_p = dFunc.sampleToData(parent_list, main_byid_df, sigma = sigma)

    # Get the data for the children
    ages_c,errors_c,nGrains_c,labels_c = dFunc.sampleToData(child_list, main_byid_df, sigma = sigma)

    # Calculate the number of grains per child sample or sample group plotted
    numGrainsPlotted = np.zeros_like(nGrains_c)
    for i in range(len(child_list   )):
        if isinstance(xaxis_1, list) == False:
            numGrainsPlotted[i] = len([elem for elem in ages_c[i] if (elem < xaxis_2 and elem > xaxis_1)]) # Number of grains in plot
        else:
            numGrainsPlotted[i] = len([elem for elem in ages_c[i] if (elem < xaxis_2[-1] and elem > xaxis_1[0])]) # Number of grains in plot (note that this assumes no gaps in what you are plotting!!!)

    # Calculate the number of samples per child distribution
    N = np.zeros_like(nGrains_c)
    if type(child_list[0])==tuple:
        for i in range(len(child_list)):
            N[i] = len(child_list[i][0])
    else:
        N = N + 1
    
    if agebins is not None:
        nage = len(agebins)-1   
    n = len(child_list)

    # Set up the figure and axes
    num_rows = c+len(parent_list)+len(child_list)+2 # Number of rows to have on the plot
    fig = plt.figure(figsize=(w1+w2+1, num_rows))

    mosaic = []
    for i in range(c): # One loop for each CDF row
        mosaic_line = [0]
        mosaic_line += ['0,1']*w1
        mosaic_line += [0]*w2
        mosaic.append(mosaic_line)
        
    mosaic.append([0]*(w1+w2+1)) # Blank row

    for i in range(len(parent_list)):
        mosaic_line = ['{},0'.format(i+1)]
        mosaic_line += ['{},1'.format(i+1)]*w1
        mosaic_line += [0]*w2
        mosaic.append(mosaic_line)
        
    mosaic.append([0]*(w1+w2+1)) # Blank row

    for i in range(len(child_list)):
        mosaic_line = ['{},0'.format(len(parent_list)+i+1)]
        mosaic_line += ['{},1'.format(len(parent_list)+i+1)]*w1
        mosaic_line += ['{},2'.format(len(parent_list)+i+1)]*w2
        mosaic.append(mosaic_line)
        
    axs = fig.subplot_mosaic(mosaic, empty_sentinel=0)
    fig.subplots_adjust(wspace=0)
    fig.subplots_adjust(hspace=0)

    # CDF plot (for parent plotting)
    axs['0,1'].set_ylim(0,1)
    axs['0,1'].set_xlim(xaxis_1,xaxis_2)

    # CDF plot (for child plotting)
    axsC = axs['0,1'].twinx()
    axsC.set_ylim(0,1)
    axsC.get_yaxis().set_ticks([])

    for i in range(len(parent_list)):
        axs['{},0'.format(i+1)].axis('off') # Turn off the pie
        # Parent relative distribution plots
        axs['{},1'.format(i+1)].get_yaxis().set_ticks([])
        axs['{},1'.format(i+1)].set_xlim(xaxis_1,xaxis_2)
        if i <= len(parent_list)-2: # Don't remove the last parent x-axis
            axs['{},1'.format(i+1)].get_xaxis().set_ticks([])
    
    for i in range(len(child_list)):
        axs['{},0'.format(len(parent_list)+i+1)].axis('off') # Turn off the pie
        # Make the relative age distribution plot
        axs['{},1'.format(len(parent_list)+i+1)].get_yaxis().set_ticks([])
        axs['{},1'.format(len(parent_list)+i+1)].set_xlim(xaxis_1,xaxis_2)
        if i <= len(child_list)-2: # Don't remove the last parent x-axis
            axs['{},1'.format(len(parent_list)+i+1)].get_xaxis().set_ticks([])
        if plotPie and (mix_coeffs_all is not None):         
            if best_plot_type == 'average':
                axs['{},0'.format(len(parent_list)+i+1)].pie(np.average(mix_coeffs_all[i],axis=0), colors=parent_colors, normalize=True)
            if best_plot_type == 'median':
                axs['{},0'.format(len(parent_list)+i+1)].pie(np.median(mix_coeffs_all[i],axis=0), colors=parent_colors, normalize=True)
            if best_plot_type == 'best-fit':
                axs['{},0'.format(len(parent_list)+i+1)].pie(mix_coeffs_bf[i], colors=parent_colors, normalize=True)

    axs['1,1'].set_title('Parents', fontsize='x-large')
    axs['{},1'.format(len(parent_list)+1)].set_title('Children', fontsize='x-large')

    plt.subplots_adjust(hspace=0, wspace=0.75)

    # Make the CDF plot
    CDF_age, CDF_p = dFunc.CDFcalcAges(ages=ages_p, x1=x1, x2=x2, xdif=xdif)
    CDF_p_min = np.min(CDF_p, axis=0)
    CDF_p_max = np.max(CDF_p, axis=0)

    axs['0,1'].fill_between(CDF_age, CDF_p_min, CDF_p_max, color='lightgray', alpha=0.5)
    for i in range(len(parent_list)):
        axs['0,1'].plot(CDF_age, CDF_p[i], color=parent_colors[i], lw=CDFlw_p, label=labels_p[i])
        if plotPie:
            axs['{},0'.format(i+1)].pie([1], colors=[parent_colors[i]], normalize=True)
    axs['0,1'].legend(loc=(1.02, 0))

    CDF_age, CDF_c = dFunc.CDFcalcAges(ages=ages_c, x1=x1, x2=x2, xdif=xdif)
    for i in range(len(child_list)):
        axsC.plot(CDF_age, CDF_c[i], color=child_colors[i], lw=CDFlw_c, ls=lineStyleMe(i, line_style), label=labels_c[i])
    axsC.legend(loc='lower right')

    # Make the relative age distribution plot for the parents
    if plotType != 'KDE' and plotType != 'PDP':
        print('Warning: plotType should equal KDE or PDP')
    if plotType == 'KDE':
        KDE_age, KDE_p = dFunc.KDEcalcAges(ages=ages_p, x1=x1, x2=x2, xdif=xdif, bw=bw, bw_x=bw_x, cumulative=False)           
        for i in range(len(parent_list)):
            axs['{},1'.format(i+1)].plot(KDE_age, KDE_p[i], color='black', lw=KDElw, label=labels_p[i])
            axs['{},1'.format(i+1)].set_ylim(0,)
            axs['{},1'.format(i+1)].text(0.98,0.92, s=labels_p[i], transform=axs['{},1'.format(i+1)].transAxes,
                horizontalalignment='right', verticalalignment='top')
            if color_by_age:
                axs['{},1'.format(i+1)] = colorByAge(KDE=KDE_p[i], agebins=agebins, agebinsc=agebinsc,
                    agebinsc_alpha=agebinsc_alpha, xdif=xdif, ax=axs['{},1'.format(i+1)])
            if fillParent:
                axs['{},1'.format(i+1)].fill_between(KDE_age, KDE_p[i], color=parent_colors[i])

    if plotType =='PDP':
        PDP_age, PDP_p = dFunc.PDPcalcAges(ages=ages_p, errors=errors_p, x1=x1, x2=x2, xdif=xdif, cumulative=False)
        for i in range(len(parent_list)):
            axs['{},1'.format(i+1)].plot(PDP_age, PDP_p[i], color='black', lw=PDPlw, label=parent_colors[i])
            axs['{},1'.format(i+1)].set_ylim(0,)
            axs['{},1'.format(i+1)].text(0.98,0.92, s=labels_p[i], transform=axs['{},1'.format(i+1)].transAxes,
                horizontalalignment='right', verticalalignment='top')
            if color_by_age:
                axs['{},1'.format(i+1)] = colorByAge(KDE=PDP_p[i], agebins=agebins, agebinsc=agebinsc,
                    agebinsc_alpha=agebinsc_alpha, xdif=xdif, ax=axs['{},1'.format(i+1)])
            if fillParent:
                axs['{},1'.format(i+1)].fill_between(PDP_age, PDP_p[i], color=parent_colors[i])

    # Make the relative age distribution plot for the children
    if plotType == 'KDE':
        KDE_age, KDE_c = dFunc.KDEcalcAges(ages=ages_c, x1=x1, x2=x2, xdif=xdif, bw=bw, bw_x=bw_x, cumulative=False)               
        for i in range(len(child_list)):
            axs['{},1'.format(len(parent_list)+i+1)].plot(KDE_age, KDE_c[i], color=child_colors[i], lw=KDElw, label=labels_c[i])
            axs['{},1'.format(len(parent_list)+i+1)].text(0.98,0.92, s=labels_c[i], transform=axs['{},1'.format(len(parent_list)+i+1)].transAxes,
                horizontalalignment='right', verticalalignment='top')

            obj_med = np.round(np.percentile(obj_vals_all[i],50), 2) # Median of the objective function

            if best_plot_type == 'average':
                axs['{},1'.format(len(parent_list)+i+1)].text(0.01,0.92, s='Avg '+objective_metric+': '+str(np.round(np.average(obj_vals_all[i],axis=0),2)), transform=axs['{},1'.format(len(parent_list)+i+1)].transAxes,
                    horizontalalignment='left', verticalalignment='top')
            if best_plot_type == 'median':
                axs['{},1'.format(len(parent_list)+i+1)].text(0.01,0.92, s='Median '+objective_metric+': '+str(np.round(np.median(obj_vals_all[i],axis=0),2)), transform=axs['{},1'.format(len(parent_list)+i+1)].transAxes,
                    horizontalalignment='left', verticalalignment='top')            
            if best_plot_type == 'best-fit':
                axs['{},1'.format(len(parent_list)+i+1)].text(0.01,0.92, s='Best-fit '+objective_metric+': '+str(np.round(obj_func_val[i],2)), transform=axs['{},1'.format(len(parent_list)+i+1)].transAxes,
                    horizontalalignment='left', verticalalignment='top')

            if color_by_age:
                axs['{},1'.format(len(parent_list)+i+1)] = colorByAge(KDE=KDE_c[i], agebins=agebins, agebinsc=agebinsc,
                    agebinsc_alpha=agebinsc_alpha, xdif=xdif, ax=axs['{},1'.format(len(parent_list)+i+1)])
        # Plot the x-axis label for the last plot
        axs['{},1'.format(len(parent_list)+i+1)].set_xlabel('Age (Ma)')

    if plotType =='PDP':
        PDP_age, PDP_c = dFunc.PDPcalcAges(ages=ages_c, errors=errors_c, x1=x1, x2=x2, xdif=xdif, cumulative=False)
        for i in range(len(child_list)):
            axs['{},1'.format(len(parent_list)+i+1)].plot(PDP_age, PDP_c[i], color=child_colors[i], lw=PDPlw, label=labels_c[i])
            axs['{},1'.format(len(parent_list)+i+1)].set_ylim(0,)
            axs['{},1'.format(len(parent_list)+i+1)].text(0.98,0.92, s=labels_c[i], transform=axs['{},1'.format(len(parent_list)+i+1)].transAxes,
                horizontalalignment='right', verticalalignment='top')
            
            if best_plot_type == 'average':
                axs['{},1'.format(len(parent_list)+i+1)].text(0.01,0.92, s='Avg '+objective_metric+': '+str(np.round(np.average(obj_vals_all[i],axis=0),2)), transform=axs['{},1'.format(len(parent_list)+i+1)].transAxes,
                    horizontalalignment='left', verticalalignment='top')
            if best_plot_type == 'median':
                axs['{},1'.format(len(parent_list)+i+1)].text(0.01,0.92, s='Median '+objective_metric+': '+str(np.round(np.median(obj_vals_all[i],axis=0),2)), transform=axs['{},1'.format(len(parent_list)+i+1)].transAxes,
                    horizontalalignment='left', verticalalignment='top')
            if best_plot_type == 'best-fit':
                axs['{},1'.format(len(parent_list)+i+1)].text(0.01,0.92, s='Best-fit '+objective_metric+': '+str(np.round(obj_func_val[i],2)), transform=axs['{},1'.format(len(parent_list)+i+1)].transAxes,
                    horizontalalignment='left', verticalalignment='top')
            if color_by_age:
                axs['{},1'.format(len(parent_list)+i+1)] = colorByAge(KDE=PDP_c[i], agebins=agebins, agebinsc=agebinsc,
                    agebinsc_alpha=agebinsc_alpha, xdif=xdif, ax=axs['{},1'.format(len(parent_list)+i+1)])
        # Plot the x-axis label for the last plot
        axs['{},1'.format(len(parent_list)+i+1)].set_xlabel('Age (Ma)')    

    # Plot the mixing results for the children
    if plotMixResults and mix_coeffs_all is not None:
        y = np.arange(0,1,1/len(parent_list))+(1/len(parent_list)*1/2)
        for i in range(len(child_list)):
            # Calculate the average and median mixing coefficients
            averages = np.average(mix_coeffs_all[i],axis=0)
            medians = np.median(mix_coeffs_all[i],axis=0)

            # Plot the mixture result on top of the child sample
            if plotType == 'KDE':
                dist = np.zeros(shape=KDE_p[0].shape)
                for j in range(len(parent_list)):
                    if best_plot_type == 'average':
                        dist += averages[j]*KDE_p[j]
                    if best_plot_type == 'median':
                        dist += medians[j]*KDE_p[j]
                    if best_plot_type == 'best-fit':
                        dist += mix_coeffs_bf[i][j]*KDE_p[j]
                axs['{},1'.format(len(parent_list)+i+1)].plot(KDE_age, dist, '--', color=child_colors[i])
                axs['{},1'.format(len(parent_list)+i+1)].set_ylim(0,)

            if plotType == 'PDP':
                dist = np.zeros(shape=PDP_p[0].shape)
                for j in range(len(parent_list)):
                    if best_plot_type == 'average':
                        dist += averages[j]*PDP_p[j]
                    if best_plot_type == 'median':
                        dist += medians[j]*PDP_p[j]
                    if best_plot_type == 'best-fit':
                        dist += mix_coeffs_bf[i][j]*PDP_p[j]
                axs['{},1'.format(len(parent_list)+i+1)].plot(PDP_age, dist, '--', color=child_colors[i])
                axs['{},1'.format(len(parent_list)+i+1)].set_ylim(0,)

            # Plot parent mixtures to the right of each child sample
            for j in range(len(parent_list)):
                p025 = np.percentile(mix_coeffs_all[i][:,j], 2.5)
                p975 = np.percentile(mix_coeffs_all[i][:,j], 97.5)
                if best_plot_type == 'average':
                    err_plus = np.abs(p975-averages[j])
                    err_minus = np.abs(averages[j]-p025)
                if best_plot_type == 'median':
                    err_plus = np.abs(p975-medians[j])
                    err_minus = np.abs(medians[j]-p025)             
                if best_plot_type == 'best-fit':
                    err_plus = np.abs(p975-mix_coeffs_bf[i][j])
                    err_minus = np.abs(mix_coeffs_bf[i][j]-p025)
                if plotResultType == 'line':
                    if best_plot_type == 'average':
                        axs['{},2'.format(len(parent_list)+i+1)].errorbar(averages[j],y[j],xerr=[[err_minus],[err_plus]], ecolor='black', fmt='o', markerfacecolor=parent_colors[j], markeredgecolor='black', markersize=10)
                    if best_plot_type == 'median':
                        axs['{},2'.format(len(parent_list)+i+1)].errorbar(medians[j],y[j],xerr=[[err_minus],[err_plus]], ecolor='black', fmt='o', markerfacecolor=parent_colors[j], markeredgecolor='black', markersize=10)
                    if best_plot_type == 'best-fit':
                        axs['{},2'.format(len(parent_list)+i+1)].errorbar(mix_coeffs_bf[i][j], y[j], xerr=[[err_minus],[err_plus]], ecolor='black', fmt='o', markerfacecolor=parent_colors[j], markeredgecolor='black', markersize=10)
                if plotResultType == 'violin':
                    parts = axs['{},2'.format(len(parent_list)+i+1)].violinplot([x[j] for x in mix_coeffs_all[i]],positions=[y[j]], vert=False, widths=violin_width, showmeans = False, showextrema = False, showmedians = True)
                    for pc in parts['bodies']:
                        pc.set_color(parent_colors[j])
                        pc.set_edgecolor('black')
                        pc.set_alpha(1)
                    vp = parts['cmedians']
                    vp.set_edgecolor('black')
                axs['{},2'.format(len(parent_list)+i+1)].set_xlim(0,1)
                axs['{},2'.format(len(parent_list)+i+1)].set_ylim(0,1)
                axs['{},2'.format(len(parent_list)+i+1)].get_yaxis().set_ticks([])
            if i <= len(child_list)-2: # Don't remove the last parent x-axis
                axs['{},2'.format(len(parent_list)+i+1)].get_xaxis().set_ticks([])
            axs['{},2'.format(len(parent_list)+i+1)].set_ylim(0,1)    
        axs['{},2'.format(len(parent_list)+i+1)].set_xlabel('Mixing proportion')
        if plotResultType == 'violin':
            axs['{},2'.format(len(parent_list)+1)].set_title('Parent mixtures (violin plot)')
        else:
            if best_plot_type == 'average':
                axs['{},2'.format(len(parent_list)+1)].set_title('Parent mixtures (average)')
            if best_plot_type == 'median':
                axs['{},2'.format(len(parent_list)+1)].set_title('Parent mixtures (median)')
            if best_plot_type == 'best-fit':
                axs['{},2'.format(len(parent_list)+1)].set_title('Parent mixtures (best-fit)')
    else:
        for i in range(len(child_list)):
            axs['{},2'.format(len(parent_list)+i+1)].set_axis_off()
            axs['{},1'.format(len(parent_list)+i+1)].set_ylim(0,)

    return fig

def lineStyleMe(i, line_style=None):
    """
    Returns a line style for a given integer input. 
    
    Parameters
    ----------
    i : integer
    line_style : list or None, optional
        List of matplotlib linestyles. Default values will be used if None. Default = None
    
    Returns
    -------
    line_style : the name of a line style ('solid','dashed','dashdot','dotted'). Uses defualt values if None.
    
    Notes
    -----
    """  
    available_stlyes = ['default','solid','dashed','dashdot','dotted']
    
    if line_style == None:
        dashes = ['solid','dashed','dashdot','dotted',(0, (5, 5)),(0, (5, 1)),(0, (3, 5, 1, 5)),(0, (3, 5, 1, 5, 1, 5)),(0, (3, 1, 1, 1, 1, 1))]  
    else:
        if not(line_style.lower() in available_stlyes):
            dashes = ['solid','dashed','dashdot','dotted',(0, (5, 5)),(0, (5, 1)),(0, (3, 5, 1, 5)),(0, (3, 5, 1, 5, 1, 5)),(0, (3, 1, 1, 1, 1, 1))]  
        else:
            dashes = [line_style]

    return dashes[i%len(dashes)]

def colorByAge(KDE, agebins, agebinsc, agebinsc_alpha, xdif, ax=None):
    '''
    Function that colors the area under a KDE or PDP plot

    Parameters
    ----------
    KDE : array
        Array of KDE or PDP to plot
    agebins : list
        Array of bin edges in Myr. Format option 1: [age1, age2, age3, etc.]. Format option 2: [[bin1_min, bin1_max],[bin2_min, bin2_max],etc.]
    agebinsc : list
        Array of colors that correspond to age bins
    agebinsc_alpha : float (0 - 1)
        Alpha value for agebinsc color fill
    xdif : float
        Spacing of age axis distributions are evaluated along
    ax : axis
        Axis to use in plotting. One is created if none is supplied. Default = None.


    Returns
    -------
    ax : axis

    '''
    
    #If no axis specified, create one
    if ax is None:
        f,ax  = plt.subplots(1,1)

    if agebinsc != None: # Only define agebins_alpha if agebinsc is specified
        if agebinsc_alpha == None: # Set default alpha of 1 if not specified
            agebinsc_alpha = np.ones(len(agebinsc))
        else:
            if not isinstance(agebinsc_alpha, list):
                agebinsc_alpha = np.ones(len(agebinsc))*agebinsc_alpha
            else: # If agebins_alpha is specified as a list, make sure there are enough values
                if len(agebinsc_alpha) < len(agebinsc):
                    print('Warning: Not enough alpha values in agebinsc_alpha')
                    return None

    if len(np.shape(agebins)) == 1:
        nage = len(agebins)-1                    
        for k in range(nage):
            xage1 = agebins[k]
            xage2 = agebins[k+1]
            KDE_agePart = np.arange(xage1, xage2+xdif, xdif)        
            KDEpart = KDE[int(xage1/xdif):int((xage2+xdif)/xdif)]
            ax.fill_between(KDE_agePart, 0, KDEpart, color=agebinsc[k], lw=0, alpha = agebinsc_alpha[k])
    if len(np.shape(agebins)) ==  2:
        for k in range(len(agebins)):
            xage1 = agebins[k][0]
            xage2 = agebins[k][1]
            KDE_agePart = np.arange(xage1, xage2+xdif, xdif)
            KDEpart = KDE[int(xage1/xdif):int((xage2+xdif)/xdif)]
            ax.fill_between(KDE_agePart, 0, KDEpart[0], color=agebinsc[k], lw=0, alpha = agebinsc_alpha[k])

    return ax

##########################################################################################################################
### Functions for exporting
##########################################################################################################################

def export_results(parent_list, child_list, main_byid_df, objective_metric, xdif, bw, mix_coeffs_bf, obj_func_val, file_name='dMix_results.xlsx',
    verbose=True, version=None, bootstrap=False, nBootstrapIterations=None, doPerturbResampledAges=None, nGrainsToResample=None, mix_coeffs_all=None,
    obj_vals_all=None, obj_func_crit=None, worse_than_crit=None, selfCompMetrics_bs_set=None):
    """
    Export mixture modeling results as an Excel spreadsheet

    Parameters
    ----------
    parent_list : list or tuple
        List or tuple of the samples to use as parents for mixture modelling.
    child_list : list or tuple
        List or tuple of the samples to use as children for mixture modelling..
    main_byid_df : dataframe
        dataframe loaded using detritalPy, e.g., with the loadDataExcel function
    objective_metric : string, optional
        Name of objective_metric that will be used to compare distributions. The default is 'Dmax'. For a 
        full list of values call print(dMix.AVAILABLE_METRICS)
    xdif : float, optional
        Spacing of age axis distributions are evaluated along. The default is 1.
    bw : float or string
        bandwidth used when calculating KDE
    mix_coeffs_bf : list
        A list of arrays, each array contains the best fit mixing coefficients.
    obj_func_val : list
        A list of the values of the objective_metric calculated with the best fit mixing coefficients.
    best_mixed_dist : list
        A list of arrays, each array contains the best fitting mixed distribution.
    file_name : string, optional
        Name of output Excel file. Defaults to 'dMix_results.xlsx'.
    verbose : bool, optional
        Lets you know when the export if finished. Defaults to True
    version : string, optional
        Version of detritalPy used. Defaults to None.
    bootstrap : bool, optional
        Exports bootstrapping results if True. Defaults to False.
    nBootstrapIterations : int or None, optional
        The number of bootstrapping iterations performed. Default is None.
    doPerturbResampledAges : bool, optional
        Whether resampled ages were additionally randomized based on the assigned uncertainty. Default is None.
    nGrainsToResample : int, optional
        The number of grains to resample with each bootstrap iteration. The default is None.
    mix_coeffs_all : list, optional
        The suite of bootstrapped mixing coefficients for each children. The default is None.
    obj_vals_all : list of arrays, optional
        A list of len(child_list) containing arrays of the values of the objective_metric calculated during bootstrapped mixture
        modelling. Default is None.
    obj_fun_crit : np.ndarray, optional
        'critical' value of the objective_metric for each value. This is calculated as the alpha-th percentile
        of self-compared objective metrics. Default is None.
    worse_than_crit : float, optional
        The proportion of bootstrapped mixture model results that have objective_metric values suggesting a greater deviation
        than obj_fun_crit. Default is None.

    Returns
    -------
    Excel file with mixture modeling results
    """
    import xlsxwriter

    workbook = xlsxwriter.Workbook(file_name)

    areGoodFitsSmallValues = lookup_functions_for_metric(objective_metric=objective_metric)[2]

    if areGoodFitsSmallValues:
        comparison_name = 'smaller'
    else:
        comparison_name = 'greater'

    #And load the data corresponding to that parent
    ages_p, errors_p, numGrains_p, labels_p = dFunc.sampleToData(parent_list,
                                                         main_byid_df,
                                                         sigma = '1sigma',
                                                         sampleLabel='Sample_ID')
    #Load all the child data
    ages_c, errors_c, numGrains_c, labels_c = dFunc.sampleToData(child_list,
                                                     main_byid_df,
                                                     sigma = '1sigma',
                                                     sampleLabel='Sample_ID')

    merge_format = workbook.add_format({'align': 'center'})
    cell_format = workbook.add_format({'align' : 'center'})
    percent_format = workbook.add_format({'num_format': '0.0%'})

    # Record model parameters
    worksheet = workbook.add_worksheet('Model_parameters')
    worksheet.write(0, 0, 'detritalPy version '+str(version))
    worksheet.write(2, 0, 'Modeling the following samples: '+str(labels_c))
    worksheet.write(3, 0, 'as a mixture of the following sources: '+str(labels_p))
    worksheet.write(5, 0, 'Objective criteria: '+str(objective_metric))
    if bootstrap:
        worksheet.write(6, 0, 'Number of bootstrapping interations: '+str(nBootstrapIterations))
        worksheet.write(7, 0, 'Redrawn ages perturbed by uncertainty? '+str(doPerturbResampledAges))
        if nGrainsToResample == None:
            worksheet.write(8, 0, 'Bootstrapped estimates are drawn by sampling the observed number of analyses (default)')
        else:
            worksheet.write(8, 0, 'Number of analyses drawn per bootstrapped estimate :'+str(nGrainsToResample))
        c=4
    else:
        c=0
    worksheet.write(6+c, 0, 'Plotting parameters:')
    worksheet.write(7+c, 0, 'xdif: '+str(xdif))
    worksheet.write(8+c, 0, 'bw: '+str(bw)+' (only relevant if using a KDE-based objective criterion)')

    if bootstrap:
        # Make summary worksheet
        worksheet = workbook.add_worksheet('Summary')
        worksheet.write(1, 0, 'Sample')
        for i, parent in enumerate(labels_p):
            worksheet.merge_range(0, i*5+1, 0, i*5+5, parent, merge_format)
            worksheet.write(1, i*5+1, 'Best-fit', cell_format)
            worksheet.write(1, i*5+2, 'Average', cell_format)
            worksheet.write(1, i*5+3, 'Median', cell_format)
            worksheet.write(1, i*5+4, 'P2.5', cell_format)
            worksheet.write(1, i*5+5, 'P97.5', cell_format)

        for i, child in enumerate(labels_c):
            worksheet.write(i+2, 0, child)
            averages = np.average(mix_coeffs_all[i],axis=0)
            medians = np.percentile(mix_coeffs_all[i], 50, axis=0)
            p025 = np.percentile(mix_coeffs_all[i], 2.5, axis=0)
            p975 = np.percentile(mix_coeffs_all[i], 97.5, axis=0)
            for j in range(len(parent_list)):
                worksheet.write(i+2, j*5+1, mix_coeffs_bf[i][j])
                worksheet.write(i+2, j*5+2, averages[j])
                worksheet.write(i+2, j*5+3, medians[j])
                worksheet.write(i+2, j*5+4, p025[j])
                worksheet.write(i+2, j*5+5, p975[j])  

    # Best-fit worksheet
    worksheet = workbook.add_worksheet('Best-fit')
    worksheet.write(0, 0, 'Sample')
    for i, parent in enumerate(labels_p):
        worksheet.write(0, i+1, parent, cell_format)    
    for i, child in enumerate(labels_c):
        worksheet.write(i+1, 0, child)
        for j in range(len(parent_list)):
            worksheet.write(i+1, j+1, mix_coeffs_bf[i][j])   

    if bootstrap:
        # Average worksheet
        worksheet = workbook.add_worksheet('Average')
        worksheet.write(0, 0, 'Sample')
        for i, parent in enumerate(labels_p):
            worksheet.write(0, i+1, parent, cell_format)    
        for i, child in enumerate(labels_c):
            worksheet.write(i+1, 0, child)
            averages = np.average(mix_coeffs_all[i],axis=0)
            for j in range(len(parent_list)):
                worksheet.write(i+1, j+1, averages[j])

    if bootstrap:
        # Median worksheet
        worksheet = workbook.add_worksheet('Median')
        worksheet.write(0, 0, 'Sample')
        for i, parent in enumerate(labels_p):
            worksheet.write(0, i+1, parent, cell_format)    
        for i, child in enumerate(labels_c):
            worksheet.write(i+1, 0, child)
            medians = np.percentile(mix_coeffs_all[i], 50, axis=0)
            for j in range(len(parent_list)):
                worksheet.write(i+1, j+1, medians[j])

    # Write model goodness-of-fit
    worksheet = workbook.add_worksheet('Model_fit')
    worksheet.write(0, 0, 'Sample')
    worksheet.write(0, 1, 'Best-fit ('+str(objective_metric)+')')
    if bootstrap:
        worksheet.write(0, 2, str(objective_metric)+'_crit')
        worksheet.write(0, 3, '% bootstrapped best-fit mixtures worse than best-fit '+str(objective_metric))
        #worksheet.write(0, 4, '% resampled observations that yielded '+comparison_name+' values of '+objective_metric+' than mixture models')
        worksheet.write(0, 4, '% resampled observations that yielded better '+'('+comparison_name+') values of '+objective_metric+' than mixture models')
    for i, child in enumerate(labels_c):
        worksheet.write(i+1, 0, child)
        worksheet.write(i+1, 1, obj_func_val[i])
        if bootstrap:
            worksheet.write(i+1, 2, obj_func_crit[i])
            worksheet.write(i+1, 3, worse_than_crit[i], percent_format)
            if areGoodFitsSmallValues:
                fraction_worse = np.sum(selfCompMetrics_bs_set[i][:nBootstrapIterations] > obj_vals_all[i][:nBootstrapIterations])/nBootstrapIterations
                worksheet.write(i+1, 4, 1-fraction_worse, percent_format)
            else:
                fraction_worse = np.sum(selfCompMetrics_bs_set[i][:nBootstrapIterations] < obj_vals_all[i][:nBootstrapIterations])/nBootstrapIterations
                worksheet.write(i+1, 4, 1-fraction_worse, percent_format)
    if bootstrap:
        worksheet.write(len(child_list)+2, 0, 'Notes')
        if areGoodFitsSmallValues:
            worksheet.write(len(child_list)+3, 0, str(objective_metric)+'_crit is the 95th percentile of comparisons between the child and bootstrapped (sampled with replacement) versions of itself')
            worksheet.write(len(child_list)+4, 0, 'The '+str(objective_metric)+' of the best-fit mixture should be less than '+str(objective_metric)+'_crit')
        else:  
            worksheet.write(len(child_list)+4, 0, str(objective_metric)+'_crit is the 5th percentile of comparisons between the child and bootstrapped (sampled with replacement) versions of itself')
            worksheet.write(len(child_list)+5, 0, 'The '+str(objective_metric)+' of the best-fit mixture should be more than '+str(objective_metric)+'_crit')
        worksheet.write(len(child_list)+5, 0, 'Good mixture should yield a low % of bootstrapped best-fit mixtures that are worse than the best-fit '+str(objective_metric))
        worksheet.write(len(child_list)+6, 0, 'Samples that are indistingshable from mixture models should have ~50% of resampled observations that are better than the best-fit '+str(objective_metric))

    if bootstrap:
        # Write mixing coefficients and objective function values for each child sample
        for k, child in enumerate(labels_c):
            worksheet = workbook.add_worksheet(child)
            # Make header
            for i, parent in enumerate(labels_p):
                worksheet.write(0, i, parent)
            worksheet.write(0, len(parent_list), objective_metric)
            for i in range(nBootstrapIterations):
                for j in range(len(parent_list)):
                    worksheet.write(i+1, j, mix_coeffs_all[k][i][j])
                worksheet.write(i+1, len(parent_list), obj_vals_all[k][i])

    # Write ages & errors for parents
    worksheet = workbook.add_worksheet('Parents')
    for i, parent in enumerate(labels_p):
        worksheet.write(0, i*2, parent)
        worksheet.write(1, i*2, 'Age')
        worksheet.write(1, 2*i+1, 'Error')
        for j in range(len(ages_p[i])):
            worksheet.write(2+j, i*2, ages_p[i][j])
            worksheet.write(2+j, 2*i+1, errors_p[i][j])
            
    # Write ages & errors for children
    worksheet = workbook.add_worksheet('Children')
    for i, child in enumerate(labels_c):
        worksheet.write(0, i*2, child)
        worksheet.write(1, i*2, 'Age')
        worksheet.write(1, 2*i+1, 'Error')
        for j in range(len(ages_c[i])):
            worksheet.write(2+j, i*2, ages_c[i][j])
            worksheet.write(2+j, 2*i+1, errors_c[i][j])

    workbook.close()

    if verbose:
        print('Workbook saved')

# Experimental function. Proceed with caution.
def calc_ss(CDF1, CDF2):
    # Square root of the sum of squared residuals between the two CDFs
    return np.sqrt(np.sum((CDF1-CDF2)**2))

# Experimental function. Proceed with caution.
# def calc_w1(CDF1, CDF2):
#     # Simplified way of calculating first Wasserstein distance
#     return np.sum(np.abs(CDF1-CDF2))