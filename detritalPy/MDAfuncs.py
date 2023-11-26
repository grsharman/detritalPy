#########################################################################################################################################   
# Maximum depositional age calculations from Sharman and Malkowski: 2020: Earth-Science Reviews (doi.org/10.1016/j.earscirev.2020.103109)
#########################################################################################################################################   

import detritalpy.detritalFuncs as dFunc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pathlib

def YSG(ages, errors, sigma=1):
    """
    Calculate the youngest single grain (YSG), where the YSG is defined as the youngest grain age plus error of uncertainty sigma (default is 1 sigma).

    Parameters
    ----------
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    sigma : (optional) Options are 1 or 2 (default is 1)

    Returns
    -------
    YSG : a list of [youngest single grain (Ma), 2-sigma error of the youngest single grain (Ma)] such that len(YSG) = len(ages)

    """

    # Check to see if ages is a list of arrays or just a single list of ages
    if not hasattr(ages[0], '__len__'):
        ages = [ages]
        errors = [errors]

    YSG = []

    for i in range(len(ages)):
        if sigma == 1: 
            data_err1s = list(zip(ages[i], errors[i]))
            data_err1s.sort(key=lambda d: d[0] + d[1]) # Sort based on age + 1s error
            YSG.append([data_err1s[0][0],data_err1s[0][1]*2]) # Reporting 2-sigma error
        if sigma == 2:
            data_err2s = list(zip(ages[i], errors[i]*2))
            data_err2s.sort(key=lambda d: d[0] + d[1]) # Sort based on age + 2s error
            YSG.append([data_err2s[0][0],data_err2s[0][1]]) # Reporting 2-sigma error

    return YSG

def YC1s(ages, errors, min_cluster_size=2, contiguous=True):
    """
    Calculate the youngest grain cluster that overlaps at 1-sigma error (see Dickinson and Gehrels (2009): Earth and Planetary Science Letters and Sharman et al. (2018): The Depositional Record for an explanation).
    Note: The result obtained is dependent on how analyses are sorted and whether a cluster consists of 'contiguous' analyses

    Paramters
    ---------
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    min_cluster_size : (optional) minimum number of grains in the cluster (default = 2)
    contiguous : bool, set to True to require that analyses in a cluster be adjacent (default = True) 

    Returns
    -------
    YC1s : a list of [the weighted mean age of the youngest cluster, the 2-sigma uncertainty of the weighted mean age of the youngest cluster, the MSWD of the youngest cluster, the number of analyses in the youngest cluster] such that len(YC1s) = len(ages)

    """

    # Check to see if ages is a list of arrays or just a single list of ages
    if not hasattr(ages[0], '__len__'):
        ages = [ages]
        errors = [errors]

    YC1s = []

    for i in range(len(ages)):

        data_err1s = list(zip(ages[i], errors[i]))
        data_err1s.sort(key=lambda d: d[0] + d[1]) # Sort based on age + 1s error

        YC1s_cluster = find_youngest_cluster(data_err1s, min_cluster_size, sort_by='none', contiguous=contiguous)
        YC1s_WM = dFunc.weightedMean(np.array([d[0] for d in YC1s_cluster]), np.array([d[1] for d in YC1s_cluster]))

        # Return NaN if YC1s did not find a cluster
        if YC1s_WM[0] == 0.0:
            YC1s.append([np.nan,np.nan,np.nan,np.nan])
        else:
            YC1s.append([YC1s_WM[0], YC1s_WM[1], YC1s_WM[2], len(YC1s_cluster)])

    return YC1s

def YC2s(ages, errors, min_cluster_size=3, contiguous=True):
    """
    Calculate the youngest grain cluster that overlaps at 2-sigma error (see Dickinson and Gehrels (2009): Earth and Planetary Science Letters and Sharman et al. (2018): The Depositional Record for an explanation).
    Note: The result obtained is dependent on how analyses are sorted and whether a cluster consists of 'contiguous' analyses

    Paramters
    ---------
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    min_cluster_size : (optional) minimum number of grains in the cluster (default = 3)
    contiguous : bool, set to True to require that analyses in a cluster be adjacent (default = True)

    Returns
    -------
    YC2s : [the weighted mean age of the youngest cluster, the 2-sigma uncertainty of the weighted mean age of the youngest cluster, the MSWD of the youngest cluster, the number of analyses in the youngest cluster] such that len(YC2s) = len(ages)

    """

    # Check to see if ages is a list of arrays or just a single list of ages
    if not hasattr(ages[0], '__len__'):
        ages = [ages]
        errors = [errors]

    YC2s = []

    for i in range(len(ages)):

        data_err2s = list(zip(ages[i], errors[i]*2))
        data_err2s.sort(key=lambda d: d[0] + d[1]) # Sort based on age + 2s error

        YC2s_cluster = find_youngest_cluster(data_err2s, min_cluster_size, sort_by='none', contiguous=contiguous)
        YC2s_WM = dFunc.weightedMean(np.array([d[0] for d in YC2s_cluster]), np.array([d[1] for d in YC2s_cluster])/2.)

        # Return NaN if YC2s did not find a cluster
        if YC2s_WM[0] == 0.0:
            YC2s.append([np.nan,np.nan,np.nan,np.nan])
        else:
            YC2s.append([YC2s_WM[0], YC2s_WM[1], YC2s_WM[2], len(YC2s_cluster)])

    return YC2s

def YDZ(ages, errors, iterations=10000, chartOutput = False, bins=25):
    from scipy import stats
    """"

    NOTE: I Have not been able to reproduce Isoplot YDZ results using this adaptation - see Appendix B of Sharman and Malkowski (Earth-Science Reviews)
    Warning: Use this algorithm at your caution

    Calculate the youngest detrital zircon age based on the Monte Carlo approach of IsoPlot (Ludwig, 2012). The youngest analyses (i.e., within 5 sigma of the youngest analysis) are repeatedly resampled by a probability distribution defined by their age and uncertainty.
    
    The YDZ is defined as the mode of the resulting distribution of ages and the uncertainty is reported as the P2.5 and P97.5. The resulting range accounts for 95% of the total range in values.
    Note that the age mode is defined as the midpoint of the histogram bin with the highest value (following IsoPlot). The mode is thus not independent of the choice of how many bins to define.
    
    Paramters
    ---------
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    iterations : (optional) the number of Monte Carlo iterations performed (default = 10000)
    chartOutput : (optional) returns a figure showing the montecarlo distribution
    bins : (optional) the number of bins to use in the histogram

    Returns
    -------
    YDZ : [the mode of youngest ages, the positive error (2-sigma), and the negative error (2-sigma)] Note that the error distribution is not symmetrical because it is based on the P2.5 and P97.5 of the distribution
    """

    # Check to see if ages is a list of arrays or just a single list of ages
    if not hasattr(ages[0], '__len__'):
        ages = [ages]
        errors = [errors]

    YDZ = []

    for i in range(len(ages)):

        data_err1s = list(zip(ages[i], errors[i]))
        
        # Identify the youngest analysis
        YSG_age, YSG_err1s = YSG(ages[i], errors[i])[0]

        ageCutoff = YSG_age + YSG_err1s*5 # 5 for 5-sigma

        # Identify all analyses within 5 sigma of the youngest analysis
        data_err1s.sort(key=lambda d: d[0]) # Sort based on age
        filtered = list(filter(lambda x: x[0] < ageCutoff, data_err1s)) # Filter out ages too old

        minAges = []
        for i in range(iterations):
            newAge_Ma = []
            for analysis in filtered:
                newAge_Ma.append(np.random.normal(loc = analysis[0], scale=analysis[1]))
            minAges.append(min(newAge_Ma))
    
        # Find the mode of the minimum ages
        binIndex, binAge = np.histogram(minAges, bins=bins)
        binMaxIndex = np.argmax(binIndex)
        binMaxAge = binAge[binMaxIndex]
        mode = binMaxAge + (binAge[binMaxIndex+1] - binMaxAge)/2

        YDZ.append([mode, np.percentile(minAges, 97.5)-mode, mode-np.percentile(minAges, 2.5)])

        if chartOutput:
            #KDE_age, KDE = dFunc.KDEcalcAges_2([minAges], bw=1, xdif=0.1)
            fig, ax = plt.subplots(1)
            ax.set_xlim(int(min(minAges))-1,int(max(minAges))+1,0.5)
            #ax.plot(KDE_age, KDE[0])
            ax.hist(minAges, bins=bins)
            ax.axvline(mode,color='black')
            ax.axvline(np.percentile(minAges,2.5),linestyle='--',color='black')
            ax.axvline(np.percentile(minAges,97.5),linestyle='--',color='black')
            ax.set_xlabel('Age (Ma)')

    return YDZ

def Y3Za(ages, errors):
    """
    Calculates the weighted mean average of the youngest three zircons, regardless of whether they overlap within error (see discussion in Coutts et al. (2019): Geoscience Frontiers)

    Parameters
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups

    Returns
    -------
    Y3Za : [the weighted mean age, the 2-sigma uncertainty, and the MSWD]

    """

    # Check to see if ages is a list of arrays or just a single list of ages
    if not hasattr(ages[0], '__len__'):
        ages = [ages]
        errors = [errors]

    Y3Za = []

    for i in range(len(ages)):
        data_err1s_ageSort = (list(zip(ages[i], errors[i])))
        data_err1s_ageSort.sort(key=lambda d: d[0]) # Sort based on age
        Y3Za_WM, Y3Za_WMerr2s, Y3Za_WM_MSWD = dFunc.weightedMean([x[0] for x in data_err1s_ageSort[:3]], [x[1] for x in data_err1s_ageSort[:3]])
        if len(ages[i]) < 3: # Return nulls if the samples has less than 3 analyses
            Y3Za.append([np.nan,np.nan,np.nan])
        else:
            Y3Za.append([Y3Za_WM, Y3Za_WMerr2s, Y3Za_WM_MSWD])

    return Y3Za


def Y3Zo(ages, errors, sigma=2, contiguous=True):
    """
    Calculates the weighted mean average of the youngest three zircons that overlap within uncertainty of sigma (default is 2-sigma) (see discussion in Coutts et al. (2019): Geoscience Frontiers)

    Parameters
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    sigma : (optional) level of uncertainty to evaluate overlap (default is 2-sigma)
    contiguous : bool, set to True to require that analyses in a cluster be adjacent (default = True)

    Returns
    -------
    Y3Zo : [the weighted mean age, the 2-sigma uncertainty, and the MSWD]

    """

    # Check to see if ages is a list of arrays or just a single list of ages
    if not hasattr(ages[0], '__len__'):
        ages = [ages]
        errors = [errors]

    Y3Zo = []

    for i in range(len(ages)):

        if sigma == 1:
            data_err1s = list(zip(ages[i], errors[i]))
            data_err1s.sort(key=lambda d: d[0] + d[1]) # Sort based on age + 1s error     
        if sigma == 2:
            data_err2s = list(zip(ages[i], errors[i]*2))
            data_err2s.sort(key=lambda d: d[0] + d[1]) # Sort based on age + 2s error

        if sigma == 1:
            Y3Zo_cluster = find_youngest_cluster(data_err1s, min_cluster_size=3, sort_by='none', contiguous=contiguous)
            Y3Zo_WM, Y3Zo_WM_err2s, Y3Zo_WM_MSWD = dFunc.weightedMean(np.array([d[0] for d in Y3Zo_cluster[:3]]), np.array([d[1] for d in Y3Zo_cluster[:3]]))
        if sigma == 2:
            Y3Zo_cluster = find_youngest_cluster(data_err2s, min_cluster_size=3, sort_by='none', contiguous=contiguous)
            Y3Zo_WM, Y3Zo_WM_err2s, Y3Zo_WM_MSWD = dFunc.weightedMean(np.array([d[0] for d in Y3Zo_cluster[:3]]), np.array([d[1]/2 for d in Y3Zo_cluster[:3]]))
        
        # Return NaN if Y3Zo did not find a cluster
        if Y3Zo_WM == 0.0:
            Y3Zo.append([np.nan,np.nan,np.nan])
        else:
            Y3Zo.append([Y3Zo_WM, Y3Zo_WM_err2s, Y3Zo_WM_MSWD])

    return Y3Zo

def YPP(ages, errors, min_cluster_size=2, thres=0.01, minDist=1, xdif=0.1):
    """
    Calculates the youngest graphical peak (after Dickinson and Gehrels (2009): EPSL) that has at least min_cluster_size analyses.

    Parameters
    ----------
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    min_cluster_size : (optional) the minimum number of analyses that fall within error of the youngest graphical peak (default = 2)
    thres : (optional) threshold of what constitues a peak (from 0 to 1). Default = 0.01
    minDist : (optional) minimum distance (Myr) between adjacent peaks. Default = 1
    xdif : (optional) bin size to compute PDP (default = 0.1 Ma)

    Returns
    -------
    YPP : [the age of the youngest graphical peak, rounded to the nearest 0.1 Ma]

    """



    import peakutils

    # Check to see if ages is a list of arrays or just a single list of ages
    if not hasattr(ages[0], '__len__'):
        ages = [ages]
        errors = [errors]

    # Calculate the PDP - note that a small xdif may be desired for increased precision
    PDP_age, PDP = dFunc.PDPcalcAges(ages, errors, xdif=xdif)

    YPP = []
    for i in range(len(ages)):

        # Calculate peak indexes
        indexes = list(peakutils.indexes(PDP[i], thres=thres, min_dist=minDist))
        # Peak ages
        peakAges = PDP_age[indexes]
        # Number of grains per peak
        peakAgeGrain = dFunc.peakAgesGrains([peakAges], [ages[i]], [errors[i]])[0]
        # Zip peak ages and grains per peak
        peakAgesGrains = list(zip(peakAges, peakAgeGrain))
        # Filter out peaks with less than min_cluster_size grains
        peakAgesGrainsFiltered = list(filter(lambda x: x[1] >= min_cluster_size, peakAgesGrains))

        # Check if a YPP was found, and if not return NaN
        if len(peakAgesGrainsFiltered) > 0:
            YPP.append(np.round(np.min([x[0] for x in peakAgesGrainsFiltered]),1))
        else:
            YPP.append(np.nan)

    return YPP

def YSP(ages, errors, min_cluster_size=2, MSWD_threshold=1, sort_by='age'):
    """
    Calculates the youngest statistical population after Coutts et al. (2019): Geoscience Frontiers. The YSP is the weighted average of the youngest group of 2 or more analyses that have a MSWD close to 1.
    If the youngest two analyses have an MSWD < MSWD_threshold, then a third grain is added and so forth.
    The final analyses to be included in the weighted average is the one with the closest MSWD value to 1.

    Warning: There is subjectivity in how the YSP algorithm is interpreted and implemented, including what threshold to use for what MSWD is
    sufficiently 'close to 1' and how analyses are sorted.

    Parameters
    ----------
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    min_cluster_size : (optional) the minimum number of analyses to calculate a MSWD from (default = 2)
    MSWD_threshold : (optional) the MSWD threshold of the youngest two analyses to begin searching for the cluster of analyses a MSWD close to 1
    sort_by : (optional) how to sort the analyses prior to looping through (options: 'age', 'age+err1s','age+err2s')

    Returns
    -------
    YSP : [the weighted mean age in Ma, the 2-sigma uncertainty of the weighted mean age, the MSWD of the weighted mean age, and the number of analyses included in the weighted mean age]

    """ 

    # Check to see if ages is a list of arrays or just a single list of ages
    if not hasattr(ages[0], '__len__'):
        ages = [ages]
        errors = [errors]   

    YSP = []

    for i in range(len(ages)): # One loop for each sample or sample group
        
        # Zip ages and errors and sort
        if sort_by == 'age':
            data_err1s_ageSort = list(zip(ages[i], errors[i]))
            data_err1s_ageSort.sort(key=lambda d: d[0]) # Sort based on age
        if sort_by == 'age+err1s':
            data_err1s_ageSort = list(zip(ages[i], errors[i]))
            data_err1s_ageSort.sort(key=lambda d: d[0] + d[1]) # Sort based on age + 1s errors
        if sort_by == 'age+err2s':
            data_err1s_ageSort = list(zip(ages[i], errors[i]))
            data_err1s_ageSort.sort(key=lambda d: d[0] + d[1]*2) # Sort based on age + 2s errors         

        for j in range(len(data_err1s_ageSort)): # One loop for each analysis. Loop repeated if MSWD of the first pair is not <1.
            # Creat list of MSWD
            MSWD = []
            for k in range(len(data_err1s_ageSort)):
                MSWD.append(dFunc.weightedMean(np.array([d[0] for d in data_err1s_ageSort[:(k+2)]]), np.array([d[1] for d in data_err1s_ageSort[:(k+2)]]))[2])

            # Add MSWD to the ages & errors tuple   
            data_err1s_MSWD = []
            for k in range(len(data_err1s_ageSort)):
                if k == 0: # Assign the first age an MSWD of 0 (so it is always included in the MSWD)
                    data_err1s_MSWD.append((data_err1s_ageSort[k][0], data_err1s_ageSort[k][1], 0))
                else: # Assign analyses the MSWD of the previous analysis, such that the filtering returns the correct analyses
                    data_err1s_MSWD.append((data_err1s_ageSort[k][0], data_err1s_ageSort[k][1], MSWD[k-1]))

            # Need to exit the algorithm if no YSP is found
            if j == len(ages[i])-1:
                YSP.append([float('nan'), float('nan'), float('nan'), float('nan')])
                break

            # Find the index of the analysis with an MSWD closest to 1
            idx = (np.abs(np.array([d[2] for d in data_err1s_MSWD][1:])-1)).argmin()+1 # Need to add 1 because we excluded the first one that had an assigned MSWD of 0

            # Filter analyses beyond the one which has a MSWD closest to 1
            agesFiltered = data_err1s_MSWD[0:idx+1]

            YSP_WM, YSP_WM_err2s, YSP_WM_MSWD = dFunc.weightedMean(np.array([d[0] for d in agesFiltered]), np.array([d[1] for d in agesFiltered]))

            if (agesFiltered[1][2] < MSWD_threshold and len(agesFiltered) >= min_cluster_size): # The first one is excluded because the MSWD is made to be 0. The second youngest analysis must have a MSWD < MSWD_threshold to proceed. The minimum cluster size must also be met or exceeded.
                YSP.append([YSP_WM, YSP_WM_err2s, YSP_WM_MSWD, len(agesFiltered)])
                break
            else:
                del data_err1s_ageSort[0] # Delete the first analysis, which was no use at all, and try again

    return YSP

def tauMethod(ages, errors, min_cluster_size=3, thres=0.01, minDist=1, xdif=0.1, x1=0, x2=4000):
    """
    Calculates the tau parameter, which is the mean weighted average of analyses that fall between probability minima (troughs) of a PDP plot (after Barbeau et al. (2009): EPSL)

    Note: results are sensitive to even minor troughs (probability minima)

    Parameters
    ----------
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    min_cluster_size : (optional) the minimum number of analyses to calculate mean weighted average (default = 3)
    thres : (optional) threshold of what constitues a peak (from 0 to 1). Default = 0.01
    minDist : (optional) minimum distance (Myr) between adjacent peaks. Default = 1
    xdif : (optional) bin size to compute PDP (default = 1 Ma)
    x1 : (optional) minimum x-axis value (default = 0 Ma)
    x2 : (optional) maximum x-axis value (default = 4000 Ma)

    Returns
    -------
    tauMethod : [the weighted mean age in Ma, the 2-sigma uncertainty of the weighted mean age, the MSWD of the weighted mean age, and the number of analyses included in the weighted mean age]

    """

    import peakutils

    # Check to see if ages is a list of arrays or just a single list of ages
    if not hasattr(ages[0], '__len__'):
        ages = [ages]
        errors = [errors]

    # Calculate the PDP - note that a small xdif may be desired for increased precision
    PDP_age, PDP = dFunc.PDPcalcAges(ages, errors, xdif=xdif, x1=x1, x2=x2)

    tauMethod = []
    for i in range(len(ages)):  

        # Calculate peak indexes
        peakIndexes = list(peakutils.indexes(PDP[i], thres=thres, min_dist=minDist))
        # Peak ages
        peakAges = PDP_age[peakIndexes]

        # Number of grains per peak
        #peakAgeGrain = dFunc.peakAgesGrains([peakAges], [ages[i]], [errors[i]])[0]

        # Calculate trough indexes
        troughIndexes = list(peakutils.indexes(PDP[i]*-1, thres=thres, min_dist=minDist))
        # Trough ages
        troughAges = [0] + list(PDP_age[troughIndexes]) + [4500] # Append a 0 because there is no trough on the young size of the youngest peak and no trough on the old side of the oldest peak

        # Identify analyses within each peak (trough-to-trough)
        ages_in_troughs = [] # Will be a list of lists
        for j in range(len(troughAges)-1):
            ages_in_troughs.append([x for x in ages[i] if troughAges[j] <= x <= troughAges[j+1]]) # Return ages within each trough

        # Count the number of analyses in eac peak (trough-to-trough)
        n_per_trough = np.asarray([len(x) for x in ages_in_troughs])

        # Determine the indexes of True values (i.e., peaks with sufficient numbers of analyses)
        true_indx = [i for i, x in enumerate(n_per_trough >= min_cluster_size) if x] # Returns indexes of True values

        # Stop the loop if no peaks are present with the min_cluster_size
        if true_indx == []:
            tauMethod.append([np.nan, np.nan, np.nan, np.nan])
            continue

        # Select the nearest trough that is younger than the youngest peak with at least min_cluster_size analyses
        troughYoung = troughAges[true_indx[0]] # Minimum bound
        troughOld = troughAges[true_indx[0]+1] # Maximum bound

        # Select ages and errors that fall between troughYoung and troughOld
        ages_errors1s = list(zip(ages[i], errors[i]))
        ages_errors1s_filtered = list(filter(lambda x: x[0] <= troughOld and x[0] >= troughYoung, ages_errors1s))

        tauMethod_WM, tauMethod_WM_err2s, tauMethod_WM_MSWD = dFunc.weightedMean(np.array([d[0] for d in ages_errors1s_filtered]), np.array([d[1] for d in ages_errors1s_filtered]))

        tauMethod.append([tauMethod_WM, tauMethod_WM_err2s, tauMethod_WM_MSWD, len(ages_errors1s_filtered)])

    return tauMethod

##### Helper Functions #####

def find_youngest_cluster(data_err, min_cluster_size, sort_by = 'age', contiguous=True):
    """
    Finds the youngest cluster of analyses that overlap within uncertainty
    (updated 6 November 2023)

    Parameters
    ----------
    data_err : array of tuples [(age1, error1), (age2, error2), etc.]
    min_cluster_size : integer, minimum number of analyses in the cluster
    sort_by : type, options: 'age+error', 'age', or 'none' (default = 'age'). Note: use 'none' if inputting a pre-sorted tuple.
    contiguous : bool, set to True to require that analyses in a cluster be adjacent (default = True) 

    Returns
    -------
    result : array of tuples, youngest ages & errors of analyses that overlap

    """
    
    if not sort_by in ['age+error','age','none']:
        print('Warning: options for sort_by kwarg are "age+error" or "age". Results may be incorrect.')

    if sort_by == 'age+error':
        data_err.sort(key=lambda d: d[0] + d[1]) # Sort based on age + error
    if sort_by == 'age':
        data_err.sort(key=lambda d: d[0]) # Sort based on age
    
    result = []
    tops = np.asarray(data_err)[:,0]+np.asarray(data_err)[:,1] # Age plus uncertainty
    bottoms = np.asarray(data_err)[:,0]-np.asarray(data_err)[:,1] # Age minus uncertainty

    for i in range(len(data_err)): # One loop for each analysis
        overlaps = bottoms[i:]<tops[i] # Boolean
        if not contiguous:
            if sum(overlaps) >= min_cluster_size: # Any analysis, regardless of order
                result = np.asarray(data_err[i:])[overlaps]
                break
        else: 
            false_indx = [i for i, x in enumerate(overlaps) if not x] # Returns indexes of False values
            if false_indx == []: # The case of all analyses overlapping or there only being 1 analysis
                if len(data_err[i:]) >= min_cluster_size: # If all analyses overlap, need to make sure there are enough of them
                    result = np.asarray(data_err[i:])
                    break
                else:
                    continue # Keep checking until no more analyses to check
            elif false_indx[0] >= min_cluster_size: # Analyses must be contiguous and overlapping
                result = np.asarray(data_err[i:(i+false_indx[0])])
                break

    return result