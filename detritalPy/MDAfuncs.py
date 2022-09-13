#########################################################################################################################################   
# Maximum depositional age calculations from Sharman and Malkowski: 2020: Earth-Science Reviews (doi.org/10.1016/j.earscirev.2020.103109)
#########################################################################################################################################   

import detritalpy.detritalFuncs as dFunc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pathlib

def YSG(ages, errors, sigma=2):
    """
    Calculate the youngest single grain (YSG), where the YSG is defined as the youngest grain age plus error of uncertainty sigma (default is 2 sigma).

    Parameters
    ----------
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    sigma : (optional) Options are 1 or 2 (default is 2)
=
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

def YC1s(ages, errors, min_cluster_size=2):
    """
    Calculate the youngest grain cluster that overlaps at 1-sigma error (see Dickinson and Gehrels (2009): Earth and Planetary Science Letters and Sharman et al. (2018): The Depositional Record for an explanation).

    Paramters
    ---------
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    min_cluster_size : (optional) minimum number of grains in the cluster (default = 2)

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
        data_err1s_ageSort = list(zip(ages[i], errors[i]))
        data_err1s_ageSort.sort(key=lambda d: d[0]) # Sort based on age
        data_err1s.sort(key=lambda d: d[0] + d[1]) # Sort based on age + 1s error


        YC1s_cluster, YC1s_imax = find_youngest_cluster(data_err1s, min_cluster_size)
        YC1s_WM = dFunc.weightedMean(np.array([d[0] for d in YC1s_cluster]), np.array([d[1] for d in YC1s_cluster]))

        # Return NaN if YC1s did not find a cluster
        if YC1s_WM[0] == 0.0:
            YC1s.append([np.nan,np.nan,np.nan,np.nan])
        else:
            YC1s.append([YC1s_WM[0], YC1s_WM[1], YC1s_WM[2], len(YC1s_cluster)])

    return YC1s

def YC2s(ages, errors, min_cluster_size=3):
    """
    Calculate the youngest grain cluster that overlaps at 2-sigma error (see Dickinson and Gehrels (2009): Earth and Planetary Science Letters and Sharman et al. (2018): The Depositional Record for an explanation).

    Paramters
    ---------
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    min_cluster_size : (optional) minimum number of grains in the cluster (default = 3)

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
        data_err2s_ageSort = list(zip(ages[i], errors[i]*2))
        data_err2s_ageSort.sort(key=lambda d: d[0]) # Sort based on age
        data_err2s.sort(key=lambda d: d[0] + d[1]) # Sort based on age + 2s error

        YC2s_cluster, YC2s_imax = find_youngest_cluster(data_err2s, min_cluster_size)
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


def Y3Zo(ages, errors, sigma=2):
    """
    Calculates the weighted mean average of the youngest three zircons that overlap within uncertainty of sigma (default is 2-sigma) (see discussion in Coutts et al. (2019): Geoscience Frontiers)

    Parameters
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    sigma : (optional) level of uncertainty to evaluate overlap (default is 2-sigma)

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
            Y3Zo_cluster, Y3Zo_imax = find_youngest_cluster(data_err1s, 3)
            Y3Zo_WM, Y3Zo_WM_err2s, Y3Zo_WM_MSWD = dFunc.weightedMean(np.array([d[0] for d in Y3Zo_cluster[:3]]), np.array([d[1] for d in Y3Zo_cluster[:3]]))
        if sigma == 2:
            Y3Zo_cluster, Y3Zo_imax = find_youngest_cluster(data_err2s, 3)
            Y3Zo_WM, Y3Zo_WM_err2s, Y3Zo_WM_MSWD = dFunc.weightedMean(np.array([d[0] for d in Y3Zo_cluster[:3]]), np.array([d[1]/2 for d in Y3Zo_cluster[:3]])/2.)
        
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

def YSP(ages, errors, min_cluster_size=2, MSWD_threshold=1):
    """
    Calculates the youngest statistical population after Coutts et al. (2019): Geoscience Frontiers. The YSP is the weighted average of the youngest group of 2 or more analyses that have a MSWD close to the MSWD_threshold (default=1),
    where the the MSWD of the youngest two analyses is less than the MSWD_threshold. The algorithm first considers the youngest two analyses. If they have an MSWD < 1, then a third grain is added and so forth.
    The final analyses to be included in the weighted average is the one with the closest value to MSWD_threshold (default of 1).

    Parameters
    ----------
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    min_cluster_size : (optional) the minimum number of analyses to calculate a MSWD from (default = 2)
    MSWD_threshold : (optional) the MSWD threshold from which to select analyses from

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
        
        # Zip ages and errors and sort by age
        data_err1s_ageSort = list(zip(ages[i], errors[i]))
        data_err1s_ageSort.sort(key=lambda d: d[0]) # Sort based on age
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
                else: # Assign analyses the MSWD of the previos analysis, such that the filtering returns the correct analyses
                    data_err1s_MSWD.append((data_err1s_ageSort[k][0], data_err1s_ageSort[k][1], MSWD[k-1]))

            # Need to exit the algorithm if no YSP is found
            if j == len(ages[i])-1:
                YSP.append([float('nan'), float('nan'), float('nan'), float('nan')])
                break

            # Find the index of the analysis with an MSWD closest to 1
            idx = (np.abs(np.array([d[2] for d in data_err1s_MSWD][1:])-MSWD_threshold)).argmin()+1 # Need to add 1 because we excluded the first one that had an assigned MSWD of 0

            # Filter analyses beyond the one which has a MSWD closest to MSWD_threshold
            agesFiltered = data_err1s_MSWD[0:idx+1]

            YSP_WM, YSP_WM_err2s, YSP_WM_MSWD = dFunc.weightedMean(np.array([d[0] for d in agesFiltered]), np.array([d[1] for d in agesFiltered]))

            if (agesFiltered[1][2] < 1 and len(agesFiltered) >= min_cluster_size): # The first one is excluded because the MSWD is made to be 0. The second youngest analysis must have a MSWD < 1 to proceed. The minimum cluster size must also be met or exceeded.
                YSP.append([YSP_WM, YSP_WM_err2s, YSP_WM_MSWD, len(agesFiltered)])
                break
            else:
                del data_err1s_ageSort[0] # Delete the first analysis, which was no use at all, and try again

    return YSP


def tauMethod(ages, errors, min_cluster_size=3, thres=0.01, minDist=1, xdif=1, chartOutput = False, x1=0, x2=4000):
    """
    Calculates the tau parameter, which is the mean weighted average of analyses that fall between probability minima (troughs) of a PDP plot (after Barbeau et al. (2009): EPSL)

    Parameters
    ----------
    ages : a 2-D array of ages, len(ages)=number of samples or sample groups
    errors : a 2-D array of 1-sigma errors for each sample or sample group, len(errors)=number of samples or sample groups
    min_cluster_size : (optional) the minimum number of analyses to calculate mean weighted average (default = 3)
    thres : (optional) threshold of what constitues a peak (from 0 to 1). Default = 0.01
    minDist : (optional) minimum distance (Myr) between adjacent peaks. Default = 1
    xdif : (optional) bin size to compute PDP (default = 1 Ma)
    chartOutput : (optional) set to True to create plots
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
    PDP_age, PDP = dFunc.PDPcalcAges(ages, errors, xdif)

    tauMethod = []
    for i in range(len(ages)):  

        # Calculate peak indexes
        peakIndexes = list(peakutils.indexes(PDP[i], thres=thres, min_dist=minDist))
        # Peak ages
        peakAges = PDP_age[peakIndexes]
        # Number of grains per peak
        peakAgeGrain = dFunc.peakAgesGrains([peakAges], [ages[i]], [errors[i]])[0]

        # Calculate trough indexes
        troughIndexes = list(peakutils.indexes(PDP[i]*-1, thres=thres, min_dist=minDist))
        # Trough ages
        troughAges = [0] + list(PDP_age[troughIndexes]) + [4500] # Append a 0 because there is no trough on the young size of the youngest peak and no trough on the old side of the oldest peak

        # Zip peak ages and grains per peak
        peakAgesGrains = list(zip(peakAges, peakAgeGrain))
        # Filter out peaks with less than min_cluster_size grains (default is 3, following Barbeau et al., 2009: EPSL)
        peakAgesGrainsFiltered = list(filter(lambda x: x[1] >= min_cluster_size, peakAgesGrains))

        # Stop the loop if no peaks are present with the min_cluster_size
        if peakAgesGrainsFiltered == []:
            tauMethod.append([np.nan, np.nan, np.nan, np.nan])
            continue

        # Select the nearest trough that is younger than the youngest peak with at least min_cluster_size analyses
        troughYoung = np.max(list(filter(lambda x: x < peakAgesGrainsFiltered[0][0], troughAges)))

        # Select the nearest trough that is older than the youngest peak with at least min_cluster_size analyses
        troughOld = np.min(list(filter(lambda x: x > peakAgesGrainsFiltered[0][0], troughAges)))

        # Select ages and errors that fall between troughYoung and troughOld
        ages_errors1s = list(zip(ages[i], errors[i]))
        ages_errors1s_filtered = list(filter(lambda x: x[0] < troughOld and x[0] > troughYoung, ages_errors1s))

        tauMethod_WM, tauMethod_WM_err2s, tauMethod_WM_MSWD = dFunc.weightedMean(np.array([d[0] for d in ages_errors1s_filtered]), np.array([d[1] for d in ages_errors1s_filtered]))

        tauMethod.append([tauMethod_WM, tauMethod_WM_err2s, tauMethod_WM_MSWD, len(ages_errors1s_filtered)])

        if chartOutput:
            fig, ax = plt.subplots(1)
            # Creates a plot output to check results
            ax.plot(PDP_age, PDP[i])
            ax.plot(PDP_age[peakIndexes], PDP[i][peakIndexes],'o')
            ax.plot(PDP_age[troughIndexes], PDP[i][troughIndexes],'o')
            ax.plot(tauMethod_WM,0,'s')
            ax.plot(ages_errors1s_filtered,np.zeros_like(ages_errors1s_filtered),'s')
            #ax.plot(tauMethod_WM-tauMethod_WM_err2s,0,'s')     
            ax.set_xlim(0,300)

    return tauMethod

##### Helper Functions #####

def find_youngest_cluster(data_err, min_cluster_size):
    """
    Finds the youngest cluster of analyses that overlap 

    Parameters
    ----------
    data_err : array of tuples [(age1, error1), (age2, error2), etc.]

    Returns
    -------
    Array of tuples of youngest cluster of analyses that overlap
    Number of analyses in the youngest cluster

    """
    i_min = 0
    i_max = 0
    for i in range(1, len(data_err)):
        top = data_err[i_min][0] + data_err[i_min][1]
        bottom = data_err[i][0] - data_err[i][1]
        if (top >= bottom):
            i_max = i
        elif i_max - i_min + 1 >= min_cluster_size:
            break
        else:
            i_min = i
            i_max = i
    return data_err[i_min: i_max + 1] if i_min < i_max else [], i_max