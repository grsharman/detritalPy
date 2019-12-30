# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 07:43:18 2017

@author: glennrsharman, jonathanpsharman, zoltansylvester
"""
#%%

###############################################################
# Import required modules
###############################################################

import numpy as np
import matplotlib.pyplot as plt
import pathlib
import pandas as pd

###############################################################
# Functions for loading a dataset and selecting samples 
###############################################################

def loadData(samples, analyses, ID_col = 'Sample_ID'):
    """
    Creates a database from a table of samples and analyses
    
    Parameters
    ----------
    samples : A table that contains sample information.
    analyses : A table that contains analysis information
    ID_col : (optional) The name of the column that contains unique sample identifiers. Default is 'Sample_ID'
    
    Returns
    -------
    main_byid_df : the database indexed by sample name
    
    Notes
    -----
    Portions of this script were provided by Kevin Befus (University of Wyoming)

    """
    main_byid_df = None
    main_byid_df = samples.copy()
    for sample_ind in range(main_byid_df.shape[0]):
        active_sample_id = main_byid_df.iloc[sample_ind][ID_col]
        if active_sample_id in analyses.index: # Allows samples to exist without analyses data
            active_UPb_data = analyses.loc[active_sample_id]
            for colname in active_UPb_data:
                if colname not in [ID_col]:
                # Make colname if not already in dataframe
                    if colname not in main_byid_df.columns:
                        main_byid_df[colname] = (np.nan*np.empty(shape=(len(main_byid_df),1))).tolist()
                        #main_byid_df[colname] = main_byid_df[colname].astype(np.ndarray)
                        main_byid_df[colname] = np.asarray(main_byid_df[colname])                       
                    main_byid_df.at[active_sample_id,colname] = active_UPb_data[colname].values
    return main_byid_df

def loadDataExcel(dataToPlot, mainSheet = 'Samples', dataSheet = 'ZrUPb', ID_col = 'Sample_ID'):
    """
    Loads an Excel file containing detrital geochronologic and/or thermochronologic data
    
    Parameters
    ----------
    dataToPlot : An array with the full filePath of each data file to be loaded, including the directory, file name, and extension
    mainSheet : (optional) The name of the Excel worksheet that contains sample information. The default name is 'Samples'
    dataSheet : (optional) The name of the Excel worksheet that contains grain analysis information. The default name is 'ZrUPb'
    ID_col : (optional) The name of the column that contains unique sample identifiers. Default is 'Sample_ID'

    Returns
    -------
    main_df : the database
    main_byid_df : the database indexed by sample name
    
    Notes
    -----
    Portions of this script were provided by Kevin Befus (University of Wyoming)
    """    
    
    obj1 = []
    obj2 = []
    obj3 = []
    obj4 = []
    for i in range(len(dataToPlot)):
        dfs = pd.read_excel(dataToPlot[i],sheet_name=None)
        main_df = None
        main_df = dfs[mainSheet]
        samples_df = main_df.copy()
        analyses_df = dfs[dataSheet]
    
        for sample_ind in range(main_df.shape[0]): # loop through entries in main_df
            active_sample_id = main_df.ix[sample_ind,ID_col]
            active_UPb_data = dfs[dataSheet].ix[dfs[dataSheet][ID_col].isin([active_sample_id]),:]
            for colname in active_UPb_data:
                if colname not in [ID_col]:
                    # Make colname if not already in dataframe
                    if colname not in main_df.columns:
                        main_df[colname] = (np.nan*np.empty(shape=(len(main_df),1))).tolist()
                        #main_df[colname] = (np.nan*np.ones(main_df.shape)).tolist()
                        main_df[colname] = np.asarray(main_df[colname])
                        #main_df[colname] = main_df[colname].astype(np.ndarray)
                    main_df.at[sample_ind,colname] = active_UPb_data[colname].values
    
        # Make a copy of the dataset and set the sample ID as index
        main_byid_df = main_df.copy()
        main_byid_df.set_index(ID_col,inplace=True,drop=False)
        obj1.append(main_df)
        obj2.append(main_byid_df)
        obj3.append(samples_df)
        obj4.append(analyses_df)
    main_df = pd.concat(obj1, sort=False)
    main_byid_df = pd.concat(obj2, sort=False)
    samples_df = pd.concat(obj3, sort=False)
    analyses_df = pd.concat(obj4, sort=False)
    
    return main_df, main_byid_df, samples_df, analyses_df

def plotSampleDist(main_byid_df, ID_col = 'Sample_ID', bestAge = 'BestAge', numBins = 25):
    """
    Plots a histogram displaying the distribution of sample size (number of analyses per sample) within the database.
    
    Parameters
    ----------
    main_df : the database from loadData()
    ID_col : (optional) The name of the column that contains unique sample identifiers. Default is 'Sample_ID'
    bestAge : (optional) The name of the column that contains the 'Best U-Pb Age' of the analysis. Default is 'BestAge'
    numBins : (optional) The number of histogram bins to plot. Default is 25
    
    Returns
    -------
    A histogram plot
    
    Notes
    -----
    """        
    numGrains = np.zeros_like(main_byid_df[ID_col])
    numSamples = len(main_byid_df[ID_col])
    c = 0 # Counter variable
    for sample in main_byid_df[ID_col]:
        numGrains[c] = len(main_byid_df.loc[sample,bestAge])
        c += 1
    fig, ax = plt.subplots(1,1)
    ax.hist(list(numGrains), numBins)
    ax.set_xlabel('Number of analyses')
    ax.set_ylabel('Frequency')
    ax.text(0.7, 0.8, s=("Samples:"+str(numSamples)), horizontalalignment='left',verticalalignment='center',transform=ax.transAxes)
    ax.text(0.7, 0.9, s=("Analyses: "+str(sum(numGrains))), horizontalalignment='left',verticalalignment='center',transform=ax.transAxes)
    
    return
                
def sampleToData(sampleList, main_byid_df, sampleLabel='Sample_ID', bestAge='BestAge', bestAgeErr='BestAge_err', sigma='1sigma'):
    """
    Returns arrays of single grain ages, 1 sigma errors, the number of analyses, and labels for
    individual samples or groups of samples
    
    Parameters
    ----------
    sampleList : array of sample IDs
        Must be in form for individual samples: ['Sample1', 'Sample2', . . . , etc.]
        Must be in the form for groups of samples: [(['Sample1','Sample2', . . . , etc.], 'Group 1 name'),
                                                    (['Sample1','Sample2', . . . , etc.], 'Group 2 name')]
    main_byid_df : DZ database, output from loadData()
    sampleLabel : (optional) The name of the column that contains the desired sample label. Default is 'Sample_ID'
    bestAge : (optional) The name of the column that contains the 'Best U-Pb Age' of the analysis. Default is 'BestAge'
    bestAgeErr : (optional) The name of the column that contains the 'Best U-Pb Age' uncertainty of the analysis. Default is 'BestAge_err'
    sigma : (optional) Specify whether bestAgeErr are 1-sigma or 2-sigma errors. Default is '1sigma', but '2sigma' can also be specified.

    Returns
    -------
    ages : array of ages for each sample or sample group
    errors : array of 1 sigma errors for each sample or sample group
    numGrains : array of numbers of analyses within each sample or sample group
    labels : array of labels for samples or sample groups
    Notes
    -----
    """
    N = len(sampleList)
    ages = []    
    errors = []
    numGrains = []
    labels = []
    if type(sampleList[0])==tuple:
        for i in range(N):
            samples = sampleList[i][0]
            sampleAges = []
            sampleErrors = []
            for sample in samples:                             
                sampleAges = np.append(sampleAges, main_byid_df.ix[sample, bestAge])
                if sigma == '2sigma':
                    sampleErrors = np.append(sampleErrors, main_byid_df.ix[sample, bestAgeErr]/2)
                else:
                    sampleErrors = np.append(sampleErrors, main_byid_df.ix[sample, bestAgeErr])                
            ages.append(sampleAges)
            errors.append(sampleErrors)
            numGrains.append(len(sampleAges))
            labels.append(sampleList[i][1])
    else:
        for sample in sampleList:
            ages.append(main_byid_df.ix[sample, bestAge])
            if sigma == '2sigma':
                errors.append(main_byid_df.ix[sample, bestAgeErr]/2.)
            else:
                errors.append(main_byid_df.ix[sample, bestAgeErr])
            numGrains.append(len(main_byid_df.ix[sample, bestAge]))
            labels.append(main_byid_df.ix[sample,sampleLabel])

    return ages, errors, numGrains, labels

def sampleToVariable(sampleList, main_byid_df, variableName):
    """
    Returns arrays of a grain variable (e.g., U_ppm) for
    an array of sample IDs or a tuple of sample IDs and group name
    
    Parameters
    ----------
    sampleList : array of sample IDs
        Must be in form for individual samples: ['Sample1', 'Sample2', . . . , etc.]
        Must be in the form for groups of samples: [(['Sample1','Sample2', . . . , etc.], 'Group 1 name'),
                                                    (['Sample1','Sample2', . . . , etc.], 'Group 2 name')]
    main_byid_df : Detrital database. Output from loadData()
    variableName : Column name of the desired variable

    Returns
    -------
    variable : array of specified variable for each sample or sample group
    
    Notes
    -----
    """
    N = len(sampleList)
    variable = []
    if type(sampleList[0])==tuple:
        for i in range(N):
            samples = sampleList[i][0]
            sampleVariables = []
            for sample in samples:
                sampleVariables = np.append(sampleVariables, main_byid_df.ix[sample, variableName])
            variable.append(sampleVariables)
    else:
        for sample in sampleList:
            variable.append(main_byid_df.ix[sample, variableName])

    return variable


###############################################################
# Functions for plotting and analyzing detrital data
###############################################################

def plotAll(sampleList, ages, errors, numGrains, labels, whatToPlot='both', separateSubplots=True, plotCDF=True, plotCPDP=False, plotCKDE=False, plotDKW=False,
 normPlots=False, plotKDE=False, colorKDE=False, colorKDEbyAge=False, plotPDP=True, colorPDP=True, colorPDPbyAge=False, plotColorBar=False, plotHist=False,
  plotLog=False, plotPIE=False, x1=0, x2=4000, b=25, bw=10, xdif=1, agebins=None, agebinsc=None, w=10, c=4, h=5, plotAgePeaks=False, agePeakOptions=None,
   CDFlw=3, KDElw=1, PDPlw=1, plotDepoAge = False, depoAge = [0], plotAgesOnCDF = False, plotHeatMap = False, heatMapType = None, heatMap = 'inferno_r'):
    """
    Creates a plot of detrital age distributions using a variety of the most common data visualization approaches. The plotting function is divided into a cumulative distribution plot and a relative distribution plot. When both are plotted together, the cumulative distribution is shown on top and the relative distribution for each sample or group of samples is shown below.

    Parameters
    ----------
    sampleList : array of sample IDs.
        Must be in form for individual samples: ['Sample1', 'Sample2', . . . , etc.].
        Must be in the form for groups of samples: [(['Sample1','Sample2', . . . , etc.], 'Group 1 name'),
                                                    (['Sample1','Sample2', . . . , etc.], 'Group 2 name')]
    ages : array of ages for each sample or sample group. Output from sampleToData()
    errors : array of errors for each sample or sample group. Output from sampleToData()
    numGrains : array of number of analyses per sample or sample group. Output from sampleToData()
    labels : array of labels for each sample or sample group. Output from sampleToData()
    whatToPlot : specify plot parameters. Options: 'cumulative', 'relative', or 'both'
    plotCDF : set to True to plot a discretized CDF
    plotCPDP : set to True to plot the cumulative PDP
    plotCKDE : set to True to plot the cumulative KDE
    normPlots : set to True to normalize the PDP and/or KDE
    plotKDE : set to True to plot a KDE
    colorKDE : set to True to color the KDE according to same coloration as used in CDF plotting
    colorKDEbyAge : set to True to color the KDE according to user-specified age populations
    plotPDP : set to True to plot a PDP
    colorPDP : set to True to color the PDP according to same coloration as used in CDF plotting
    colorPDPbyAge : set to True to color the PDP according to user-specified age populations
    plotColorBar : set to True to color age categories as vertical bars
    plotHist : set to True to plot a histogram
    plotLog : set to True to plot x-axis and y-axis on a log scale    
    plotPie : set to True to plot a pie diagram, using user-specified age populations
    x1 : lower limit (Myr) of x-axis
    x2 : upper limit (Myr) of x-axis
    b : histogram bin size (Myr)
    bw : KDE bandwidth. Options are 'optimizedFixed', 'optimizedVariable', or a number (bandwidth in Myr)
    xdif : interval (Myr) over which distributions are calculated
    agebins : array of bin edges in Myr. Format option 1: [age1, age2, age3, etc.]. Format option 2: [[bin1_min, bin1_max],[bin2_min, bin2_max],etc.]
    agebinsc : array of colors that correspond to age bins
    w : width of the plot
    c : height of the plot
    CDFlw : (optional) weight of CDF line
    KDElw : (optional) weight of KDE line
    PDPlw : (optional) weight of PDP line
    plotDepoAge: (optional) set to True to plot the depositional age of samples (only available if separateSubplots == True)
    depoAge: (optional) List of depositional age (Ma) of samples. If len(depoAge) == 1, then the one single age will be used throughout.
    plotAgesOnCDF: (optional) set to True to plot individual analyses with error on the CDF plot

    Returns
    -------
    fig : a figure with the plotted age distribution(s)
    
    Notes
    -----
    """

    if separateSubplots:
        fig = plotAll_1(sampleList, ages, errors, numGrains, labels, whatToPlot, plotCDF, plotCPDP, plotCKDE, plotDKW, normPlots, plotKDE, 
            colorKDE, colorKDEbyAge, plotPDP, colorPDP, colorPDPbyAge, plotColorBar, plotHist, plotLog, plotPIE, x1, x2, b, bw, xdif, agebins, 
            agebinsc, w, c, plotAgePeaks, agePeakOptions, CDFlw, KDElw, PDPlw, plotDepoAge, depoAge, plotAgesOnCDF, plotHeatMap, heatMapType, heatMap)
    else:
        fig = plotAll_2(sampleList, ages, errors, numGrains, labels, whatToPlot, plotCDF, plotCPDP, plotCKDE, plotDKW, normPlots, plotKDE, 
            colorKDE, colorKDEbyAge, plotPDP, colorPDP, colorPDPbyAge, plotColorBar, plotHist, plotLog, plotPIE, x1, x2, b, bw, xdif, agebins, 
            agebinsc, w, c, h, CDFlw, KDElw, PDPlw, plotAgesOnCDF)
    return fig

def plotAll_1(sampleList, ages, errors, numGrains, labels, whatToPlot, plotCDF, plotCPDP, plotCKDE, plotDKW, normPlots, plotKDE, colorKDE, 
    colorKDEbyAge, plotPDP, colorPDP, colorPDPbyAge, plotColorBar, plotHist, plotLog, plotPIE, x1, x2, b, bw, xdif, agebins, agebinsc, w, c, 
    plotAgePeaks, agePeakOptions, CDFlw, KDElw, PDPlw, plotDepoAge, depoAge, plotAgesOnCDF, plotHeatMap, heatMapType, heatMap):


    if (plotLog and x1 == 0):
        x1 = 0.1 # Ensures that 0 will not be plotted on a log scale

    # Calculate the number of grains per sample or sample group plotted
    numGrainsPlotted = np.zeros_like(numGrains)
    for i in range(len(sampleList)):
        numGrainsPlotted[i] = len([elem for elem in ages[i] if (elem < x2 and elem > x1)]) # Number of grains in plot
        
    # Number of samples per plotted distribution
    N = np.zeros_like(numGrains)
    if type(sampleList[0])==tuple:
        for i in range(len(sampleList)):
            N[i] = len(sampleList[i][0])
    else:
        N = N + 1
    
    if agebins is not None:
        nage = len(agebins)-1   
    n = len(sampleList)
    
    # Sets figure font options
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 14}
    plt.rc('font', **font)
    
    # Sets the matplotlib figure and axes structure
    if whatToPlot == 'cumulative':
        fig, axs = plt.subplots(c,w, figsize=(w,c))
        axs[0,0] = plt.subplot2grid((c,w),(0,0),rowspan=c) # empty subplot
        axs[0,0].axis('off') # delete axis from empty subplot
        axs[0,1] = plt.subplot2grid((c,w),(0,1),rowspan=c,colspan=w-1) # panel for CDF plot
        if plotLog:
            axs[0,0].set_xscale('log')
            axs[0,1].set_xscale('log')
    if (whatToPlot == 'both' or whatToPlot == 'relative'):
        if whatToPlot == 'relative':
            c = 0
        if (n == 1 and whatToPlot == 'relative'):
            c = c+1 # To avoid an index error when only one sample or sample group is plotted
        fig, axs = plt.subplots(n+c,w, figsize=(w,n+c))
        if c > 0:
            axs[0,0] = plt.subplot2grid((n+c,w),(0,0),rowspan=c) # empty subplot
            axs[0,1] = plt.subplot2grid((n+c,w),(0,1),rowspan=c,colspan=w-1) # panel for CDF plot
        if c == 0:
            axs[0,0] = plt.subplot2grid((n,w),(0,0)) # empty subplot
            axs[0,1] = plt.subplot2grid((n,w),(0,1),colspan=w-1)
        axs[0,0].axis('off') # delete axis from empty subplot
        axs[0,1].get_xaxis().set_ticks([])
        axs[0,1].tick_params(direction='out')
        axs[0,1].yaxis.set_ticks_position('left')
        if plotLog:
            axs[0,0].set_xscale('log')
            axs[0,1].set_xscale('log')
        if (n == 1 and whatToPlot == 'relative'):
            axs[0,1].axis('off') # Hide this axis if n = 1
        for i in range(n):
            axs[c+i,0] = plt.subplot2grid((n+c,w),(c+i,0)) # subplots for pie plots
            axs[c+i,0].get_xaxis().set_ticks([])
            axs[c+i,0].get_yaxis().set_ticks([])
            if not plotPIE:
                axs[c+i,0].axis('off') # delete axis from pie subplot if pies not plotted
            axs[c+i,1] = plt.subplot2grid((n+c,w),(c+i,1),colspan=w-1) # subplots for KDE plots
            axs[c+i,1].get_yaxis().set_ticks([])
            if i<n-1: # This insures that only the last plot will have an x-axis
                axs[c+i,1].get_xaxis().set_ticks([])
                if plotLog:
                    axs[c+i, 1].set_xscale('log')
            else:
                axs[c+i,1].tick_params(direction='out')
                axs[c+i,1].xaxis.set_ticks_position('bottom')
                axs[c+i,1].set_xlabel('Age (Ma)')
                if plotLog:
                    axs[c+i, 1].set_xscale('log')
        fig.subplots_adjust(wspace=0)
        fig.subplots_adjust(hspace=0)
        
    # Plot the cumulative distribution
    if (whatToPlot == 'cumulative' or whatToPlot == 'both'):
        if plotDKW: # Calculate the Kvoretsky-Kiefer-Wolfowitz inequality following Anderson et al. (2018): Basin Research
            alpha = 0.05 # Default is 95% confidence interval
            epsilon = np.empty(shape=(len(numGrains),1))
            for i in range(len(sampleList)):
                epsilon[i] = np.sqrt(np.log(2./alpha)/(2.*numGrains[i]))       
        if plotCDF:
            CDF_age, CDF = CDFcalcAges(ages=ages, x1=0, x2=4500, xdif=xdif)
            for i in range(len(sampleList)):
                axs[0,1].plot(CDF_age, CDF[i], color=colorMe(i), alpha=1, lw=CDFlw, label=labels[i]+(', N=(%d' % N[i])+(', %d' % numGrainsPlotted[i])+('/%d' % numGrains[i])+(')'))
                if plotDKW:
                    DFWmin, DFWmax = calcDFW(CDF[i], epsilon[i])
                    axs[0,1].plot(CDF_age, DFWmax, '--', color=colorMe(i))
                    axs[0,1].plot(CDF_age, DFWmin, '--', color=colorMe(i))
                    axs[0,1].fill_between(CDF_age, DFWmax, DFWmin, color=colorMe(i), alpha = 0.5)
                if plotAgesOnCDF:
                    agesErrorsSort = pd.DataFrame({'Ages' : ages[i],'Errors' : errors[i]})
                    agesErrorsSort = agesErrorsSort.sort_values(by=['Ages'])
                    agesErrorsSortY = np.arange(0,1.,1/len(agesErrorsSort['Ages']))
                    axs[0,1].errorbar(x=agesErrorsSort['Ages'], y=agesErrorsSortY, xerr=agesErrorsSort['Errors'],linestyle='',marker='o',ecolor='black',capsize=2,markerfacecolor=colorMe(i),markeredgecolor='black')        # Need to add functionality to additional options below
        if plotCPDP:
            CPDP_age, CPDP = PDPcalcAges(ages=ages, errors=errors, x1=0, x2=4500, xdif=xdif, cumulative=True)        
            for i in range(len(sampleList)):
                axs[0,1].plot(CPDP_age, CPDP[i], color=colorMe(i), alpha=1, lw=CDFlw, label=labels[i]+(', N=(%d' % N[i])+(', %d' % numGrainsPlotted[i])+('/%d' % numGrains[i])+(')'))
                if plotDKW:
                    DFWmin, DFWmax = calcDFW(CPDP[i], epsilon[i])
                    axs[0,1].plot(CPDP_age, DFWmax, '--', color=colorMe(i))
                    axs[0,1].plot(CPDP_age, DFWmin, '--', color=colorMe(i))
                    axs[0,1].fill_between(CPDP_age, DFWmax, DFWmin, color=colorMe(i), alpha = 0.5)
                if plotAgesOnCDF:
                    agesErrorsSort = pd.DataFrame({'Ages' : ages[i],'Errors' : errors[i]})
                    agesErrorsSort = agesErrorsSort.sort_values(by=['Ages'])
                    agesErrorsSortY = np.arange(0,1.,1/len(agesErrorsSort['Ages']))
                    axs[0,1].errorbar(x=agesErrorsSort['Ages'], y=agesErrorsSortY, xerr=agesErrorsSort['Errors'],linestyle='',marker='o',ecolor='black',capsize=2,markerfacecolor=colorMe(i),markeredgecolor='black')
        if plotCKDE:
            if bw == 'optimizedFixed':
                CKDE_age, CKDE = KDEcalcAgesLocalAdapt(ages=ages, x1=0, x2=4500, xdif=xdif, cumulative=True)
            if bw == 'optimizedVariable':
                CKDE_age, CKDE = KDEcalcAgesGlobalAdapt(ages=ages, x1=0, x2=4500, xdif=xdif, cumulative=True)
            if type(bw) != str:
                CKDE_age, CKDE = KDEcalcAges_2(ages=ages, x1=0, x2=4500, xdif=xdif, bw=bw, cumulative=True)
            for i in range(len(sampleList)):
                axs[0,1].plot(CKDE_age, CKDE[i], color=colorMe(i), alpha=1, lw=CDFlw, label=labels[i]+(', N=(%d' % N[i])+(', %d' % numGrainsPlotted[i])+('/%d' % numGrains[i])+(')'))
                if plotDKW:
                    DFWmin, DFWmax = calcDFW(CKDE[i], epsilon[i])
                    axs[0,1].plot(CKDE_age, DFWmax, '--', color=colorMe(i))
                    axs[0,1].plot(CKDE_age, DFWmin, '--', color=colorMe(i))
                    axs[0,1].fill_between(CKDE_age, DFWmax, DFWmin, color=colorMe(i), alpha = 0.5)
                if plotAgesOnCDF:
                    agesErrorsSort = pd.DataFrame({'Ages' : ages[i],'Errors' : errors[i]})
                    agesErrorsSort = agesErrorsSort.sort_values(by=['Ages'])
                    agesErrorsSortY = np.arange(0,1.,1/len(agesErrorsSort['Ages']))
                    axs[0,1].errorbar(x=agesErrorsSort['Ages'], y=agesErrorsSortY, xerr=agesErrorsSort['Errors'],linestyle='',marker='o',ecolor='black',capsize=2,markerfacecolor=colorMe(i),markeredgecolor='black')
        axs[0,1].set_ylabel("Cumulative Distribution")
        axs[0,1].legend(loc="lower right", prop={'size':8})
        axs[0,1].set_xlim(x1, x2)
        axs[0,1].set_ylim(0, 1.)
        if whatToPlot == 'both':
            axs[0,1].get_xaxis().set_ticks([])
            axs[0,1].get_yaxis().set_visible(True)
        else:
            axs[0,1].set_xlabel('Age (Ma)')
        if plotColorBar:
            if len(np.shape(agebins)) == 1:
                for j in range(nage):
                    axs[0,1].axvspan(xmin=agebins[j],xmax=agebins[j+1], color = agebinsc[j])
            if len(np.shape(agebins)) == 2:
                for j in range(len(agebins)):
                    axs[0,1].axvspan(xmin=agebins[j][0],xmax=agebins[j][1], color = agebinsc[j])

        # Plot depositional age as a vertical line, if selected
        if plotDepoAge:
            if len(depoAge) == 1:
                axs[0,1].axvline(x=depoAge, color='darkred')
            else:
                for i in range(len(depoAge)):
                    axs[0,1].axvline(x=depoAge[i], color=colorMe(i))       
    
    # Plot the relative distribution (PDP and/or KDE)
    if (whatToPlot == 'both' or whatToPlot == 'relative'):
        # Cycle through each sample for normalized plots            
        if plotKDE or (plotHeatMap and heatMapType == 'KDE'):
            if bw == 'optimizedFixed':
                KDE_age, KDE = KDEcalcAgesLocalAdapt(ages=ages, x1=0, x2=4500, xdif=xdif, cumulative=False)
            if bw == 'optimizedVariable':
                KDE_age, KDE = KDEcalcAgesGlobalAdapt(ages=ages, x1=0, x2=4500, xdif=xdif, cumulative=False)
            if type(bw) != str:
                KDE_age, KDE = KDEcalcAges_2(ages=ages, x1=0, x2=4500, xdif=xdif, bw=bw, cumulative=False)
        if plotPDP or (plotHeatMap and heatMapType == 'PDP'):
            PDP_age, PDP = PDPcalcAges(ages=ages, errors=errors, x1=0, x2=4500, xdif=xdif, cumulative=False)
        if plotAgePeaks:
            if (agePeakOptions[0] == 'KDE' and plotKDE):
                peakAges, indexes, peakAgesGrains = peakAge(KDE_age, KDE, ages, errors, thres=agePeakOptions[1], minDist=agePeakOptions[2], minPeakSize=agePeakOptions[3])
            if (agePeakOptions[0] == 'PDP' and plotPDP):
                peakAges, indexes, peakAgesGrains = peakAge(PDP_age, PDP, ages, errors, thres=agePeakOptions[1], minDist=agePeakOptions[2], minPeakSize=agePeakOptions[3])

        # Calculated if plotting a histogram heat map
        if plotHeatMap and heatMapType == 'hist':
            dist = []
            for i in range(len(sampleList)):
                dist.append(np.histogram(ages[i], bins = np.arange(x1, x2+xdif, xdif), range=(x1, x2))[0])

        for i in range(len(sampleList)):
            # Plot the KDE as a heat map
            if plotHeatMap:
                axHeat = axs[c+i,1].twinx()
                axHeat.set_xlim([x1, x2])
                axHeat.get_yaxis().set_visible(False)
                if heatMapType == 'KDE':
                    axHeat.imshow([KDE[i]], aspect='auto', cmap=heatMap, interpolation='none', extent=[0-xdif*0.5, 4500-xdif*0.5, 0, 1])
                if heatMapType == 'PDP':
                    axHeat.imshow([PDP[i]], aspect='auto', cmap=heatMap, interpolation='none', extent=[0-xdif*0.5, 4500-xdif*.5, 0, 1])
                if heatMapType == 'hist':
                    axHeat.imshow([dist[i]], aspect='auto', cmap=heatMap, interpolation='none', extent=[x1-xdif*0.5, x2-xdif*0.5, 0, 1])
            
            # KDE plot
            if plotKDE:
                axKDE = axs[c+i,1].twinx()
                # Plot depositional age as a vertical line, if selected                
                if plotDepoAge:
                    if len(depoAge) == 1:
                        axKDE.axvline(x=depoAge, color='darkred')
                    else:
                        axKDE.axvline(x=depoAge[i], color=colorMe(i))
                # Plot KDE as a line
                axKDE.plot(KDE_age, KDE[i], color='black', lw=KDElw, label=labels[i])
                # Plot age peaks
                if (plotAgePeaks and agePeakOptions[0] == 'KDE'):
                    axKDE.plot(KDE_age[indexes[i]],KDE[i][indexes[i]], '|', color='black')
                    if agePeakOptions[4]:
                        for j in range(len(peakAges[i])):
                            if (peakAges[i][j]>x1 and peakAges[i][j]<x2): # Only plot the peak age if within plotting range
                                axKDE.text(x=KDE_age[indexes[i][j]],y=KDE[i][indexes[i][j]], s=peakAges[i][j], size='x-small')
                    exportPeakAge(labels, peakAges, peakAgesGrains)
                # Fill the KDE      
                if colorKDE:
                    axKDE.fill_between(KDE_age, 0, KDE[i], alpha = 1, color=colorMe(i), lw=0)
                if colorKDEbyAge:
                    if len(np.shape(agebins)) == 1:
                        nage = len(agebins)-1                    
                        for k in range(nage):
                            xage1 = agebins[k]
                            xage2 = agebins[k+1]
                            KDE_agePart = np.arange(xage1, xage2+xdif, xdif)        
                            KDEpart = KDE[i][int(xage1/xdif):int((xage2+xdif)/xdif)]
                            axKDE.fill_between(KDE_agePart, 0, KDEpart, alpha = 1, color=agebinsc[k], lw=0)
                    if len(np.shape(agebins)) ==  2:
                        for k in range(len(agebins)):
                            xage1 = agebins[k][0]
                            xage2 = agebins[k][1]
                            KDE_agePart = np.arange(xage1, xage2+xdif, xdif)
                            KDEpart = KDE[i][int(xage1/xdif):int((xage2+xdif)/xdif)]
                            axKDE.fill_between(KDE_agePart, 0, KDEpart, alpha = 1, color=agebinsc[k], lw=0)
                axKDE.set_xlim(x1, x2)
                axKDE.legend(loc="upper right", prop={'size':8})
                # Adjust the y-axis scale, depending on normalization
                if normPlots:
                    kdeMax = 0
                    for k in range(len(sampleList)):
                        if max(KDE[k]) > kdeMax:
                            kdeMax = max(KDE[k])
                    axKDE.set_ylim(0, kdeMax)
                else:
                    axKDE.set_ylim([0, max(KDE[i])+max(KDE[i])*0.05])                
                axKDE.get_yaxis().set_visible(False)
    
            # PDP plot
            if plotPDP:
                axPDP = axs[c+i,1].twinx() # to allow the PDP to plot on a different scale
                # Plot depositional age as a vertical line, if selected                
                if plotDepoAge:
                    if len(depoAge) == 1:
                        axs[c+i,1].axvline(x=depoAge, color='darkred')
                    else:
                        axs[c+i,1].axvline(x=depoAge[i], color=colorMe(i))
                axPDP.plot(PDP_age, PDP[i], color='black', ls='-', alpha=1, lw=PDPlw, label=labels[i])
                if not plotKDE: # Only print the label if the KDE is not already plotted
                    axPDP.legend(loc="upper right", prop={'size':8})
                # Plot age peaks
                if (plotAgePeaks and agePeakOptions[0] == 'PDP'):
                    axPDP.plot(PDP_age[indexes[i]],PDP[i][indexes[i]], '|', color='black')
                    if agePeakOptions[4]:
                        for j in range(len(peakAges[i])):
                            if (peakAges[i][j]>x1 and peakAges[i][j]<x2): # Only plot the peak age if within plotting range
                                axPDP.text(x=PDP_age[indexes[i][j]],y=PDP[i][indexes[i][j]], s=peakAges[i][j], size='x-small')
                    exportPeakAge(labels, peakAges, peakAgesGrains)
                if colorPDP:
                    axPDP.fill_between(PDP_age, PDP[i], alpha = 1, color=colorMe(i))
                if colorPDPbyAge:
                    if len(np.shape(agebins)) == 1:
                        nage = len(agebins)-1
                        for j in range(nage):                
                            xage1 = agebins[j]
                            xage2 = agebins[j+1]
                            if (xage2 > x2 and xage1 <= x2): # Avoids a problem that would otherwise occur if any age bins are greater than x2
                                xage2 = x2
                            if (xage2 > x2 and xage1 >= x2):
                                break
                            PDP_agePart = np.arange(xage1, xage2+xdif, xdif)
                            PDPpart = PDP[i][int(xage1/xdif):int((xage2+xdif)/xdif)]
                            axPDP.fill_between(PDP_agePart, 0, PDPpart, alpha = 1, color=agebinsc[j])
                    if len(np.shape(agebins)) == 2:
                        for j in range(len(agebins)):
                            xage1 = agebins[j][0]
                            xage2 = agebins[j][1]
                            if (xage2 > x2 and xage1 <= x2): # Avoids a problem that would otherwise occur if any age bins are greater than x2
                                xage2 = x2
                            if (xage2 > x2 and xage1 >= x2):
                                break
                            PDP_agePart = np.arange(xage1, xage2+xdif, xdif)
                            PDPpart = PDP[i][int(xage1/xdif):int((xage2+xdif)/xdif)]
                            axPDP.fill_between(PDP_agePart, 0, PDPpart, alpha = 1, color=agebinsc[j])
                axPDP.set_xlim([x1, x2])
                if normPlots:
                    pdfMax = 0
                    for k in range(len(sampleList)):
                        if max(PDP[k]) > pdfMax:
                            pdfMax = max(PDP[k]) 
                    axPDP.set_ylim([0, pdfMax])
                else:
                    axPDP.set_ylim([0, max(PDP[i])+max(PDP[i])*0.05])
                axPDP.get_yaxis().set_visible(False)
                if plotKDE:
                    axPDP.get_xaxis().set_visible(False) # Do not plot the x-axis if it has already been plotted
                
            # Histogram plot
            if plotHist:
                axHist = axs[c+i,1].twinx() # to allow the histogram to plot on a different scale
                # Plot depositional age as a vertical line, if selected                
                if plotDepoAge:
                    if len(depoAge) == 1:
                        axs[c+i,1].axvline(x=depoAge, color='darkred')
                    else:
                        axs[c+i,1].axvline(x=depoAge[i], color=colorMe(i))
                bin_array = np.arange(x1, x2+xdif, b)
                axHist.hist(ages[i], bins=bin_array, color='black', fill=None, alpha=1, histtype='bar', density=False)
                axHist.set_xlim([x1, x2]) # Use this code to set the x-axis scale
                if normPlots:
                    histMax = 0
                    for k in range(len(sampleList)):
                        if max(np.histogram(ages[k], bins=bin_array)[0]) > histMax:
                            histMax = max(np.histogram(ages[k], bins=bin_array)[0]) 
                    axHist.set_ylim([0, histMax])
                axHist.get_yaxis().set_visible(True) # This makes the y-axis numbers invisible
                if (plotPDP or plotKDE):
                    axHist.get_xaxis().set_visible(False) # Do not plot the x-axis if it has already been plotted
                
            # Pie plot
            if plotPIE:
                if len(np.shape(agebins)) == 1:
                    hist = np.histogram(ages[i], agebins)[0]
                if len(np.shape(agebins)) == 2:
                    hist = []
                    for j in range(len(agebins)):
                        hist.append(np.histogram(ages[i],agebins[j])[0][0])
                axs[c+i,0].pie(hist, colors=agebinsc, startangle=90, counterclock=False, radius=0.75)
            
            # Plot colored vertical bars, if selected
            if plotColorBar:
                if len(np.shape(agebins)) == 1:
                    for j in range(nage):
                        axs[c+i,1].axvspan(xmin=agebins[j],xmax=agebins[j+1], color = agebinsc[j])
                if len(np.shape(agebins)) == 2:
                    for j in range(len(agebins)):
                        axs[c+i,1].axvspan(xmin=agebins[j][0],xmax=agebins[j][1], color = agebinsc[j])

    return fig

def plotAll_2(sampleList, ages, errors, numGrains, labels, whatToPlot, plotCDF, plotCPDP, plotCKDE, plotDKW, normPlots, plotKDE, 
    colorKDE, colorKDEbyAge, plotPDP, colorPDP, colorPDPbyAge, plotColorBar, plotHist, plotLog, plotPIE, x1, x2, b, bw, xdif, agebins, 
    agebinsc, w, c, h, CDFlw, KDElw, PDPlw, plotAgesOnCDF):
    # Reverse sample order, to make plotting order consistent with plotAll_1()
    sampleList = sampleList[::-1]
    ages = ages[::-1]
    errors = errors[::-1]
    numGrains = numGrains[::-1]
    labels = labels[::-1]
    
    if (plotLog and x1 == 0):
        x1 = 0.1 # Ensures that 0 will not be plotted on a log scale

    # Calculate the number of grains per sample or sample group plotted
    numGrainsPlotted = np.zeros_like(numGrains)
    for i in range(len(sampleList)):
        numGrainsPlotted[i] = len([elem for elem in ages[i] if (elem < x2 and elem > x1)]) # Number of grains in plot
    
    # Number of samples per plotted distribution
    N = np.zeros_like(numGrains)
    if type(sampleList[0])==tuple:
        N = N + 1
    else:
        for i in range(len(sampleList)):
            N[i] = len(sampleList[i][0])
    nage = len(agebins)-1
    n = len(sampleList)
    if h == 'auto':
        h = n
    
    # Sets figure font options
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 14}
    plt.rc('font', **font)

    # Sets the matplotlib figure and axes structure
    if whatToPlot == 'cumulative':
        fig, axs = plt.subplots(c,w, figsize=(w,c))
        axs[0,0] = plt.subplot2grid((c,w),(0,0),rowspan=c,colspan=w) # panel for CDF plot
        if plotLog:
            axs[0,0].set_xscale('log')
    if (whatToPlot == 'both' or whatToPlot == 'relative'):
        if whatToPlot == 'relative':
            c = 0
        if (n == 1 and whatToPlot == 'relative'):
            c = c+1 # To avoid an index error when only one sample or sample group is plotted            
        fig, axs = plt.subplots(h+c,w, figsize=(w+2,h+c))
        if c > 0:
            axs[0,0] = plt.subplot2grid((h+c,w),(0,0),rowspan=c,colspan=w) # panel for CDF plot
        if c == 0:
            axs[0,0] = plt.subplot2grid((h,w),(0,0),colspan=w)
        axs[0,0].get_xaxis().set_ticks([])
        axs[0,0].tick_params(direction='out')
        axs[0,0].yaxis.set_ticks_position('left')
        axs[0,0].set_xlim(x1,x2)
        if plotLog:
            axs[0,0].set_xscale('log')
        if (n == 1 and whatToPlot == 'relative'):
            axs[0,0].axis('off')
        axs[c,0] = plt.subplot2grid((h+c,w),(c,0),rowspan=h,colspan=w) # panel for the relative plot
        axs[c,0].get_yaxis().set_ticks([])
        axs[c,0].tick_params(direction='out')
        axs[c,0].xaxis.set_ticks_position('bottom')
        axs[c,0].set_xlabel('Age (Ma)')
        if plotLog:
            axs[c,0].set_xscale('log')
        fig.subplots_adjust(wspace=0)
        fig.subplots_adjust(hspace=0)        
        
    # Plot the cumulative distribution
    if (whatToPlot == 'cumulative' or whatToPlot == 'both'):
        if plotDKW: # Calculate the Kvoretsky-Kiefer-Wolfowitz inequality following Anderson et al. (2018): Basin Research
            alpha = 0.05 # Default is 95% confidence interval
            epsilon = np.empty(shape=(len(numGrains),1))
            for i in range(len(sampleList)):
                epsilon[i] = np.sqrt(np.log(2./alpha)/(2.*numGrains[i]))            
        if plotCDF:
            CDF_age, CDF = CDFcalcAges(ages=ages, x1=0, x2=4500, xdif=xdif)
            for i in range(len(sampleList)):
                axs[0,0].plot(CDF_age, CDF[i], color=colorMe(i), alpha=1, lw=CDFlw, label=str(labels[i])+(', N=(%d' % N[i])+(', %d' % numGrainsPlotted[i])+('/%d' % numGrains[i])+(')'))
                if plotDKW:
                    DFWmin, DFWmax = calcDFW(CDF[i], epsilon[i])
                    axs[0,0].plot(CDF_age, DFWmax, '--', color=colorMe(i))
                    axs[0,0].plot(CDF_age, DFWmin, '--', color=colorMe(i))
                    axs[0,0].fill_between(CDF_age, DFWmax, DFWmin, color=colorMe(i), alpha = 0.5)
                if plotAgesOnCDF:
                    agesErrorsSort = pd.DataFrame({'Ages' : ages[i],'Errors' : errors[i]})
                    agesErrorsSort = agesErrorsSort.sort_values(by=['Ages'])
                    agesErrorsSortY = np.arange(0,1.,1/len(agesErrorsSort['Ages']))
                    axs[0,0].errorbar(x=agesErrorsSort['Ages'], y=agesErrorsSortY, xerr=agesErrorsSort['Errors'],linestyle='',marker='o',ecolor='black',capsize=2,markerfacecolor=colorMe(i),markeredgecolor='black')
        if plotCPDP:
            CPDP_age, CPDP = PDPcalcAges(ages=ages, errors=errors, x1=0, x2=4500, xdif=xdif, cumulative=True)        
            for i in range(len(sampleList)):
                axs[0,0].plot(CPDP_age, CPDP[i], color=colorMe(i), alpha=1, lw=CDFlw, label=str(labels[i])+(', N=(%d' % N[i])+(', %d' % numGrainsPlotted[i])+('/%d' % numGrains[i])+(')'))
                if plotDKW:
                    DFWmin, DFWmax = calcDFW(CPDP[i], epsilon[i])
                    axs[0,0].plot(CPDP_age, DFWmax, '--', color=colorMe(i))
                    axs[0,0].plot(CPDP_age, DFWmin, '--', color=colorMe(i))
                    axs[0,0].fill_between(CPDP_age, DFWmax, DFWmin, color=colorMe(i), alpha = 0.5)
                if plotAgesOnCDF:
                    agesErrorsSort = pd.DataFrame({'Ages' : ages[i],'Errors' : errors[i]})
                    agesErrorsSort = agesErrorsSort.sort_values(by=['Ages'])
                    agesErrorsSortY = np.arange(0,1.,1/len(agesErrorsSort['Ages']))
                    axs[0,0].errorbar(x=agesErrorsSort['Ages'], y=agesErrorsSortY, xerr=agesErrorsSort['Errors'],linestyle='',marker='o',ecolor='black',capsize=2,markerfacecolor=colorMe(i),markeredgecolor='black')
        if plotCKDE:
            if bw == 'optimizedFixed':
                CKDE_age, CKDE = KDEcalcAgesLocalAdapt(ages=ages, x1=0, x2=4500, xdif=xdif, cumulative=True)
            if bw == 'optimizedVariable':
                CKDE_age, CKDE = KDEcalcAgesGlobalAdapt(ages=ages, x1=0, x2=4500, xdif=xdif, cumulative=True)
            if type(bw) != str:
                CKDE_age, CKDE = KDEcalcAges_2(ages=ages, x1=0, x2=4500, xdif=xdif, bw=bw, cumulative=True)
            for i in range(len(sampleList)):
                axs[0,0].plot(CKDE_age, CKDE[i], color=colorMe(i), alpha=1, lw=CDFlw, label=str(labels[i])+(', N=(%d' % N[i])+(', %d' % numGrainsPlotted[i])+('/%d' % numGrains[i])+(')'))
                if plotDKW:
                    DFWmin, DFWmax = calcDFW(CKDE[i], epsilon[i])
                    axs[0,0].plot(CKDE_age, DFWmax, '--', color=colorMe(i))
                    axs[0,0].plot(CKDE_age, DFWmin, '--', color=colorMe(i))
                    axs[0,0].fill_between(CKDE_age, DFWmax, DFWmin, color=colorMe(i), alpha = 0.5)
                if plotAgesOnCDF:
                    agesErrorsSort = pd.DataFrame({'Ages' : ages[i],'Errors' : errors[i]})
                    agesErrorsSort = agesErrorsSort.sort_values(by=['Ages'])
                    agesErrorsSortY = np.arange(0,1.,1/len(agesErrorsSort['Ages']))
                    axs[0,0].errorbar(x=agesErrorsSort['Ages'], y=agesErrorsSortY, xerr=agesErrorsSort['Errors'],linestyle='',marker='o',ecolor='black',capsize=2,markerfacecolor=colorMe(i),markeredgecolor='black')
        axs[0,0].set_ylabel("Cumulative Distribution")
        axs[0,0].legend(loc="lower right", prop={'size':8})
        axs[0,0].set_xlim(x1, x2)
        axs[0,0].set_ylim(0, 1.)
        if whatToPlot == 'both':
            axs[0,0].get_xaxis().set_ticks([])
            axs[0,0].get_yaxis().set_visible(True)
        else:
            axs[0,0].set_xlabel('Age (Ma)')
        if plotColorBar:
            if len(np.shape(agebins)) == 1:
                for j in range(nage):
                    axs[0,0].axvspan(xmin=agebins[j],xmax=agebins[j+1], color = agebinsc[j])
            if len(np.shape(agebins)) == 2:
                for j in range(len(agebins)):
                    axs[0,0].axvspan(xmin=agebins[j][0],xmax=agebins[j][1], color = agebinsc[j])

    # Plot the relative distribution (PDP and/or KDE)
    if (whatToPlot == 'both' or whatToPlot == 'relative'):
        # Cycle through each sample for normalized plots            
        if plotKDE:
            if bw == 'optimizedFixed':
                KDE_age, KDE = KDEcalcAgesLocalAdapt(ages=ages, x1=0, x2=4500, xdif=xdif, cumulative=False)
            if bw == 'optimizedVariable':
                KDE_age, KDE = KDEcalcAgesGlobalAdapt(ages=ages, x1=0, x2=4500, xdif=xdif, cumulative=False)
            if type(bw) != str:
                KDE_age, KDE = KDEcalcAges_2(ages=ages, x1=0, x2=4500, xdif=xdif, bw=bw, cumulative=False)
            # Determine the maximum KDE value for each sample or sample group
            KDEmax = np.empty(shape=(n,1))
            np.zeros_like(numGrains)
            for i in range(len(sampleList)):
                KDEmax[i] = np.max(KDE[i])

        if plotPDP:
            PDP_age, PDP = PDPcalcAges(ages=ages, errors=errors, x1=0, x2=4500, xdif=xdif, cumulative=False)  
            # Determine the maximum PDP value for each sample or sample group
            PDPmax = np.empty(shape=(n,1))
            for i in range(len(sampleList)):
                PDPmax[i] = np.max(PDP[i])
        
        # Calculate a list of y-axis displacements
        if (plotPDP and plotKDE):
            distMax = list(map(max, zip(KDEmax, PDPmax))) 
        else:
            if plotKDE:
                distMax = KDEmax.copy()
            if plotPDP:
                distMax = PDPmax.copy()
        distMaxCumSum = [0] + list(np.cumsum(distMax))       
        
        # Create KDE and/or PDP plots
        for i in range(len(sampleList)):
            if plotKDE:
                # Shift KDE y-values                
                KDEshift = KDE.copy()
                KDEshift[i] = KDEshift[i]+distMaxCumSum[i]
                axs[c,0].plot(KDE_age, KDEshift[i], color='black', lw=KDElw)
                if colorKDE:
                    axs[c,0].fill_between(KDE_age, np.min(KDEshift[i]), KDEshift[i], alpha = 1, color=colorMe(i), lw=0)

                if colorKDEbyAge:
                    if len(np.shape(agebins)) == 1:
                        for k in range(nage):
                            xage1 = agebins[k]
                            xage2 = agebins[k+1]
                            KDE_agePart = np.arange(xage1, xage2+xdif, xdif)        
                            KDEshiftPart = KDEshift[i][int(xage1/xdif):int((xage2+xdif)/xdif)]
                            axs[c,0].fill_between(KDE_agePart, np.min(KDEshift[i]), KDEshiftPart, alpha = 1, color=agebinsc[k], lw=0)
                    if len(np.shape(agebins)) ==  2:
                        for k in range(len(agebins)):
                            xage1 = agebins[k][0]
                            xage2 = agebins[k][1]
                            KDE_agePart = np.arange(xage1, xage2+xdif, xdif)
                            KDEshiftPart = KDEshift[i][int(xage1/xdif):int((xage2+xdif)/xdif)]
                            axs[c,0].fill_between(KDE_agePart, np.min(KDEshift[i]), KDEshiftPart, alpha = 1, color=agebinsc[k], lw=0)
            if plotPDP:
                # Shift PDP y-values
                PDPshift = PDP.copy()
                PDPshift[i] = PDPshift[i]+distMaxCumSum[i]
                axs[c,0].plot(PDP_age, PDPshift[i], color='black', lw=PDPlw)
                if colorPDP:
                    axs[c,0].fill_between(PDP_age, np.min(PDPshift[i]), PDPshift[i], alpha = 1, color=colorMe(i), lw=0)
                if colorPDPbyAge:
                    if len(np.shape(agebins)) == 1:        
                        for k in range(nage):
                            xage1 = agebins[k]
                            xage2 = agebins[k+1]
                            PDP_agePart = np.arange(xage1, xage2+xdif, xdif)        
                            PDPshiftPart = PDPshift[i][int(xage1/xdif):int((xage2+xdif)/xdif)]
                            axs[c,0].fill_between(PDP_agePart, np.min(PDPshift[i]), PDPshiftPart, alpha = 1, color=agebinsc[k], lw=0)
                    if len(np.shape(agebins)) == 2:
                        for k in range(len(agebins)):
                            xage1 = agebins[k][0]
                            xage2 = agebins[k][1]
                            PDP_agePart = np.arange(xage1, xage2+xdif, xdif)        
                            PDPshiftPart = PDPshift[i][int(xage1/xdif):int((xage2+xdif)/xdif)]
                            axs[c,0].fill_between(PDP_agePart, np.min(PDPshift[i]), PDPshiftPart, alpha = 1, color=agebinsc[k], lw=0)                          
            axs[c,0].text(x2+(x2-x1)*0.01, distMaxCumSum[i], s=labels[i], size='x-small')
        
        # Plot colored vertical bars, if selected
        if plotColorBar:
            if len(np.shape(agebins)) == 1:
                for j in range(nage):
                    axs[c,0].axvspan(xmin=agebins[j],xmax=agebins[j+1], color = agebinsc[j])
            if len(np.shape(agebins)) == 2:
                for j in range(len(agebins)):
                    axs[c,0].axvspan(xmin=agebins[j][0],xmax=agebins[j][1], color = agebinsc[j])

        axs[c,0].set_ylim(0)
        axs[c,0].set_xlim(x1,x2)
   
    return fig             

def plotRimsVsCores(main_byid_df, sampleList, ages, errors, labels, x1=0, x2=4000, y1=0, y2=4000, plotLog=False, plotError=False, w=8, c=8, grainIDcol='Grain_ID', rimCoreCol='RimCore', rimID='Rim', coreID='Core', bestAge='BestAge',bestAgeErr='BestAge_err'):
    """
    Creates a plot of rims versus cores.

    Parameters
    ----------
    main_byid_df : Detrital database. Output from loadData()
    sampleList : array of sample IDs.
        Must be in form for individual samples: ['Sample1', 'Sample2', . . . , etc.].
        Must be in the form for groups of samples: [(['Sample1','Sample2', . . . , etc.], 'Group 1 name'),
                                                    (['Sample1','Sample2', . . . , etc.], 'Group 2 name')]
    ages : array of ages for each sample or sample group. Output from sampleToData()
    errors : array of errors for each sample or sample group. Output from sampleToData()
    labels : array of labels for each sample or sample group. Output from sampleToData()
    x1 : lower limit (Myr) of x-axis
    x2 : upper limit (Myr) of x-axis
    y1 : lower limit (Myr) of y-axis
    y2 : upper limit (Myr) of y-axis
    plotLog : will plot x-axis and y-axis on a log scale if set to True
    plotError : will plot error bars if set to True
    w : width of the plot
    c : height of the plot
    grainIDcol : (optional) The name of the column that contains a unique identifer for each grain. Default is 'Grain_ID'
    rimCoreCol : (optional) The name of the column that identifies the analysis as a rim or core. Default is 'RimCore'
    rimID : (optional) The identifier of rim analyses. Default is "Rim"
    coreID : (optional) The identifier of core analyses. Default is "Core"
    bestAge : (optional) The name of the column that contains the 'Best U-Pb Age' of the analysis. Default is 'BestAge'
    bestAgeErr : (optional) The name of the column that contains the 'Best U-Pb Age' uncertainty of the analysis. Default is 'BestAge_err'

    Returns
    -------
    figRimCore : a figure with the plotted rim(s) versus core(s)
    
    Notes
    -----
    """    
    from matplotlib.path import Path
    import matplotlib.patches as patches

    # Sets figure font options
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 14}
    plt.rc('font', **font)

    # Assign parameters for plotting shading
    verts = [
            (x1, x1), #left, bottom
            (x1, y2), #left, top
            (x2, y2),
            (x2, x2), #right, top
            (x1, x1)  # ignored
            ]
    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY
             ]
    path = Path(verts, codes)

    figRimCore, ax = plt.subplots(1,1, figsize=(w,c))

    if plotLog:
        ax.set_xscale('log')
        ax.set_yscale('log')
    
    patch = patches.PathPatch(path, facecolor='lightgray', lw=0)
    ax.add_patch(patch)

    for i in range(len(sampleList)): # One loop for each sample or sample group
        rimVsCore = {} # Create an empty dictionary
        if type(sampleList[0])==tuple:
            c = 0 # Counter variable used to ensure that the label is only plotted once per sample group
            for j in range(len(sampleList[i][0])): # One loop for each sample in each group
                for k in range(len(main_byid_df.ix[sampleList[i][0][j],bestAge])): # One loop for each analysis in each sample
                    grainID = main_byid_df.ix[sampleList[i][0][j],grainIDcol][k]
                    if not grainID in rimVsCore:
                        rimVsCore[grainID] = [None, None, None, None]
                    if main_byid_df.ix[sampleList[i][0][j],rimCoreCol][k] == rimID:
                        rimVsCore[grainID][0] = main_byid_df.ix[sampleList[i][0][j],bestAge][k]
                        rimVsCore[grainID][1] = main_byid_df.ix[sampleList[i][0][j],bestAgeErr][k]
                    elif main_byid_df.ix[sampleList[i][0][j],rimCoreCol][k] == coreID:
                        rimVsCore[grainID][2] = main_byid_df.ix[sampleList[i][0][j],bestAge][k]
                        rimVsCore[grainID][3] = main_byid_df.ix[sampleList[i][0][j],bestAgeErr][k]            
                for grain in list(rimVsCore):
                    if rimVsCore[grain][0] is None or rimVsCore[grain][2] is None:
                        del rimVsCore[grain]
                for grain in rimVsCore:
                    if plotError:
                        ax.errorbar(rimVsCore[grain][2], rimVsCore[grain][0], rimVsCore[grain][3], rimVsCore[grain][1], fmt='s', color=colorMe(i), ecolor='gray', capthick=2, label=labels[i] if c == 0 else "")
                    else:
                        ax.plot(rimVsCore[grain][2],rimVsCore[grain][0],'s',color=colorMe(i),label=labels[i] if c == 0 else "")
                    c = c+1        
        else:
            for j in range(len(ages[i])): # One loop for analysis
                grainID = main_byid_df.ix[sampleList[i],grainIDcol][j]
                if not grainID in rimVsCore:
                    rimVsCore[grainID] = [None, None, None, None]
                if main_byid_df.ix[sampleList[i],rimCoreCol][j] == rimID:
                    rimVsCore[grainID][0] = main_byid_df.ix[sampleList[i],bestAge][j]
                    rimVsCore[grainID][1] = main_byid_df.ix[sampleList[i],bestAgeErr][j]
                elif main_byid_df.ix[sampleList[i],rimCoreCol][j] == coreID:
                    rimVsCore[grainID][2] = main_byid_df.ix[sampleList[i],bestAge][j]
                    rimVsCore[grainID][3] = main_byid_df.ix[sampleList[i],bestAgeErr][j]            
            for grain in list(rimVsCore):
                if rimVsCore[grain][0] is None or rimVsCore[grain][2] is None:
                    del rimVsCore[grain]
            c = 0 # Counter variable used to ensure that the label is only plotted once per sample
            for grain in rimVsCore:
                if plotError:
                    ax.errorbar(rimVsCore[grain][2], rimVsCore[grain][0], rimVsCore[grain][3], rimVsCore[grain][1], fmt='s', 
                            color=colorMe(i), ecolor='gray', capthick=2, label=labels[i] if c == 0 else "")
                else:
                    ax.plot(rimVsCore[grain][2],rimVsCore[grain][0],'s',color=colorMe(i),label=labels[i] if c == 0 else "")
                c = c+1        
    ax.set_xlim(x1,x2)
    ax.set_xlabel('Core age (Ma)')
    ax.set_ylim(y1,y2)
    ax.set_ylabel('Rim age (Ma)')
    ax.yaxis.set_label_position('right')    
    ax.yaxis.set_ticks_position('right')
    ax.legend(loc="upper left", prop={'size':8})

    return figRimCore

def plotDouble(sampleList, main_byid_df, ages, errors, numGrains, labels, variableName, plotError, variableError, normPlots, 
    plotKDE, colorKDE, colorKDEbyAge, plotPDP, colorPDP, colorPDPbyAge, plotHist, x1, x2, 
    autoScaleY, y1, y2, b, bw, xdif, agebins, agebinsc, w, t, l, plotLog, plotColorBar, 
    plotMovingAverage, windowSize, KDElw=1, PDPlw=1, averageType = 'Mean'):
    """
    Creates a figure where a numeric variable is plotted above detrital age distributions for each sample or sample group. Examples could include the uranium concentration (U_ppm), the thorium to uranium ratio (Th_U), the epsilon hafnium value (eHf), or the concentration of a trace element.

    Parameters
    ----------
    sampleList : array of sample IDs.
        Must be in form for individual samples: ['Sample1', 'Sample2', . . . , etc.].
        Must be in the form for groups of samples: [(['Sample1','Sample2', . . . , etc.], 'Group 1 name'),
                                                    (['Sample1','Sample2', . . . , etc.], 'Group 2 name')]
    main_byid_df : Detrital database. Output from loadData()
    ages : array of ages for each sample or sample group. Output from sampleToData()
    errors : array of errors for each sample or sample group. Output from sampleToData()
    numGrains : array of number of analyses per sample or sample group. Output from sampleToData()
    labels : array of labels for each sample or sample group. Output from sampleToData()
    variableName : label of the data column to plot (e.g., 'Th_U')
    plotError : set to True to plot error bars
    variableError : select the label of the data column containing errors (e.g., 'Th_U_err'), or specify the error as a percentage (e.g., 0.05)
    normPlots : set to True to normalize the PDP and/or KDE    
    plotKDE : set to True to plot a KDE
    colorKDE : set to True to color the KDE according to same coloration as used in CDF plotting
    colorKDEbyAge : set to True to color the KDE according to user-specified age populations
    plotPDP : set to True to plot a PDP
    colorPDP : set to True to color the PDP according to same coloration as used in CDF plotting
    colorPDPbyAge : set to True to color the PDP according to user-specified age populations
    plotHist : set to True to plot a histogram    
    x1 : lower limit (Myr) of x-axis
    x2 : upper limit (Myr) of x-axis
    autoScaleY : set to True to automatically select the y-axis scale
    y1 : lower limit of the y-axis (if autoScaleY is set to False)
    y2 : upper limit of the y-axis (if autoScaleY is set to False)
    b : histogram bin size (Myr)
    bw : KDE bandwidth. Options are 'optimizedFixed', 'optimizedVariable', or a number (bandwidth in Myr)
    xdif : interval (Myr) over which distributions are calculated
    agebins : array of bin edges (Myr)
    agebinsc : array of colors that correspond to age bins
    w : width of the plot
    t : height of the top panel
    l : height of the bottom panel
    plotLog : set to True to plot x-axis and y-axis on a log scale    
    plotColorBar : set to True to color age categories as vertical bars
    plotMovingAverage : set to True to plot a moving average
    windowSize : specify the window of the moving average
    KDElw : (optional) weight of KDE line
    PDPlw : (optional) weight of PDP line

    Returns
    -------
    fig : a figure with the plotted variable and age distribution(s)
    
    Notes
    -----
    """ 
    variables = sampleToVariable(sampleList, main_byid_df, variableName)

    if type(variableError) is str:
        variableErrorToPlot = sampleToVariable(sampleList, main_byid_df, variableError)

    if (plotLog and x1 == 0):
        x1 = 0.1 # Ensures that 0 will not be plotted on a log scale

    # Calculate the number of grains per sample or sample group plotted
    numGrainsPlotted = np.zeros_like(numGrains)
    for i in range(len(sampleList)):
        numGrainsPlotted[i] = len([elem for elem in ages[i] if (elem < x2 and elem > x1)]) # Number of grains in plot
    
    # Number of samples per plotted distribution
    N = np.zeros_like(numGrains)
    if type(sampleList[0])==tuple:
        for i in range(len(sampleList)):
            N[i] = len(sampleList[i][0])
    else:
        N = N + 1
    
    nage = len(agebins)-1   
    n = len(sampleList)
    
    # Sets figure font options
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 14}
    plt.rc('font', **font)
    
    # Set up the figure grid
    fig, axs = plt.subplots(n*(t+l),w, figsize=(w,n*(t+l)))
    c = 0 # counter variable
    for i in range(n):
        axs[c,0] = plt.subplot2grid((n*(t+l),w),(c,0),colspan=w, rowspan=t) # panel for top plot
        axs[c+t,0] = plt.subplot2grid((n*(t+l),w),(c+t,0),colspan=w, rowspan=l) # panel for bottom plot (age distribution)
        axs[c,0].set_ylabel(variableName)
        axs[c,0].get_xaxis().set_ticks([])
        axs[c+t,0].get_yaxis().set_ticks([])
        if i<n-1: # This insures that only the last plot will have an x-axis
            axs[c+t,0].get_xaxis().set_ticks([])            
        else:
            axs[c+t,0].tick_params(direction='out')
            axs[c+t,0].xaxis.set_ticks_position('bottom')
            axs[c+t,0].set_xlabel('Age (Ma)')
        if plotLog:
            axs[c,0].set_xscale('log')
            axs[c+t,0].set_xscale('log')            
        fig.subplots_adjust(wspace=0)
        fig.subplots_adjust(hspace=0)
        c = c+t+l
    
    # Plot
    c = 0 # counter variable
    for i in range(n):
        if autoScaleY:
            # Determine the minimum and maximum variable in the dataset
            y1 = min(variables[i])
            y2 = max(variables[i])
        # Create the upper scatter plot
        # Note that only data points within y1-y2 range will be plotted
        axs[c,0].scatter(ages[i][np.where((variables[i] > y1) & (variables[i] < y2))],variables[i][np.where((variables[i] > y1) & (variables[i] < y2))], color='white', edgecolor='black', marker='s', s=10, zorder=2)
        axs[c,0].set_xlim(x1, x2)
        axs[c,0].set_ylim(y1, y2)

        if plotError:
            if type(variableError) is not str:
                axs[c,0].errorbar(ages[i][np.where((variables[i] > y1) & (variables[i] < y2))], variables[i][np.where((variables[i] > y1) & (variables[i] < y2))],
                 xerr = errors[i][np.where((variables[i] > y1) & (variables[i] < y2))], yerr = variables[i][np.where((variables[i] > y1) & (variables[i] < y2))]*variableError,
                  fmt='none', color='black', ecolor='black', capthick=2, zorder=1)
            else:
                axs[c,0].errorbar(ages[i][np.where((variables[i] > y1) & (variables[i] < y2))], variables[i][np.where((variables[i] > y1) & (variables[i] < y2))],
                 xerr = errors[i][np.where((variables[i] > y1) & (variables[i] < y2))], yerr = variableErrorToPlot[i][np.where((variables[i] > y1) & (variables[i] < y2))],
                  fmt='none', color='black', ecolor='black', capthick=2, zorder=1)
        if plotMovingAverage:
            def isNaN(num):
                return num != num
            zipAgeVariable = list(zip(ages[i], variables[i]))
            zipAgeVariableF = [(x,y) for x, y in zipAgeVariable if not isNaN(y)] # Exclude grains with no data
            zipAgeVariableF.sort(key=lambda d: d[0]) # Sort based on age
            ageSort = np.asarray(zipAgeVariableF)[:,0]
            variableSort = np.asarray(zipAgeVariableF)[:,1]

            if averageType == 'Mean':
                # Moving average option #1
                def movingAverage(interval, windowSize):
                    window = np.ones(int(windowSize))/float(windowSize)
                    return np.convolve(interval, window, 'same')
                # Moving average option #2
                def running_mean(x, N):
                    cumsum = np.cumsum(np.insert(x, 0, 0)) 
                    return (cumsum[N:] - cumsum[:-N]) / float(N)
                xAvg = running_mean(ageSort,windowSize)
                yAvg = running_mean(variableSort,windowSize)
            if averageType == 'Median':
                def RunningMedian(x,N):
                    idx = np.arange(N) + np.arange(len(x)-N+1)[:,None]
                    b = [row[row>0] for row in x[idx]]
                    return list(map(np.median,b))
                    #return np.array(map(np.median,b))
                xAvg = RunningMedian(ageSort, windowSize)
                yAvg = RunningMedian(variableSort,windowSize)
            axs[c,0].plot(xAvg,yAvg,color='darkred',lw=2)
                
        # Plot the relative distribution (PDP and/or KDE)
        if plotKDE:
            if bw == 'optimizedFixed':
                KDE_age, KDE = KDEcalcAgesLocalAdapt(ages=ages, x1=0, x2=4500, xdif=xdif, cumulative=False)
            if bw == 'optimizedVariable':
                KDE_age, KDE = KDEcalcAgesGlobalAdapt(ages=ages, x1=0, x2=4500, xdif=xdif, cumulative=False)
            if type(bw) != str:
                KDE_age, KDE = KDEcalcAges_2(ages=ages, x1=0, x2=4500, xdif=xdif, bw=bw, cumulative=False)
        if plotPDP:
            PDP_age, PDP = PDPcalcAges(ages=ages, errors=errors, x1=x1, x2=x2, xdif=xdif, cumulative=False)
            
        # KDE plot
        if plotKDE:
            # Plot KDE as a line
            axs[c+t,0].plot(KDE_age, KDE[i], color='black', lw=KDElw, label=labels[i])
            # Fill the KDE      
            if colorKDE:
                axs[c+t,0].fill_between(KDE_age, 0, KDE[i], alpha = 1, color=colorMe(i), lw=0)
            if colorKDEbyAge:
                if len(np.shape(agebins)) == 1:
                    for k in range(nage):
                        xage1 = agebins[k]
                        xage2 = agebins[k+1]
                        KDE_agePart = np.arange(xage1, xage2+xdif, xdif)        
                        KDEpart = KDE[i][int(xage1/xdif):int((xage2+xdif)/xdif)]
                        axs[c+t,0].fill_between(KDE_agePart, 0, KDEpart, alpha = 1, color=agebinsc[k], lw=0)
                if len(np.shape(agebins)) == 2:
                    for k in range(len(agebins)):
                        xage1 = agebins[k][0]
                        xage2 = agebins[k][1]
                        KDE_agePart = np.arange(xage1, xage2+xdif, xdif)        
                        KDEpart = KDE[i][int(xage1/xdif):int((xage2+xdif)/xdif)]
                        axs[c+t,0].fill_between(KDE_agePart, 0, KDEpart, alpha = 1, color=agebinsc[k], lw=0)
            axs[c+t,0].set_xlim(x1, x2)
            axs[c+t,0].legend(loc="upper right", prop={'size':8})
            # Adjust the y-axis scale, depending on normalization
            if normPlots:
                kdeMax = 0
                for k in range(len(sampleList)):
                    if max(KDE[k]) > kdeMax:
                        kdeMax = max(KDE[k])
                axs[c+t,0].set_ylim(0, kdeMax)
            else:
                axs[c+t,0].set_ylim([0, max(KDE[i])+max(KDE[i])*0.05])                
            axs[c+1,0].get_yaxis().set_visible(False)
    
        # PDP plot
        if plotPDP:
            axPDP = axs[c+t,0].twinx() # to allow the PDP to plot on a different scale
            axPDP.plot(PDP_age, PDP[i], color='black', ls='-', alpha=1, lw=PDPlw, label=labels[i])
            if not plotKDE: # Only print the label if the KDE is not already plotted
                axPDP.legend(loc="upper right", prop={'size':8})
            if colorPDP:
                axPDP.fill_between(PDP_age, PDP[i], alpha = 1, color=colorMe(i))
            if colorPDPbyAge:
                if len(np.shape(agebins)) == 1:
                    nage = len(agebins)-1
                    for j in range(nage):                
                        xage1 = agebins[j]
                        xage2 = agebins[j+1]
                        if (xage2 > x2 and xage1 <= x2): # Avoids a problem that would otherwise occur if any age bins are greater than x2
                            xage2 = x2
                        if (xage2 > x2 and xage1 >= x2):
                            break
                        PDP_agePart = np.arange(xage1, xage2+xdif, xdif)
                        PDPpart = PDP[i][int(xage1/xdif):int((xage2+xdif)/xdif)]
                        axPDP.fill_between(PDP_agePart, 0, PDPpart, alpha = 1, color=agebinsc[j])
                if len(np.shape(agebins)) == 2:
                    for j in range(len(agebins)):
                        xage1 = agebins[j][0]
                        xage2 = agebins[j][1]
                        if (xage2 > x2 and xage1 <= x2): # Avoids a problem that would otherwise occur if any age bins are greater than x2
                            xage2 = x2
                        if (xage2 > x2 and xage1 >= x2):
                            break
                        PDP_agePart = np.arange(xage1, xage2+xdif, xdif)
                        PDPpart = PDP[i][int(xage1/xdif):int((xage2+xdif)/xdif)]
                        axPDP.fill_between(PDP_agePart, 0, PDPpart, alpha = 1, color=agebinsc[j])
            axPDP.set_xlim([x1, x2])
            if normPlots:
                pdfMax = 0
                for k in range(len(sampleList)):
                    if max(PDP[k]) > pdfMax:
                        pdfMax = max(PDP[k]) 
                axPDP.set_ylim([0, pdfMax])
            else:
                axPDP.set_ylim([0, max(PDP[i])+max(PDP[i])*0.05])
            axPDP.get_yaxis().set_visible(False)
            if plotKDE:
                axPDP.get_xaxis().set_visible(False) # Do not plot the x-axis if it has already been plotted
                
        # Histogram plot
        if plotHist:
            axHist = axs[c+t,0].twinx() # to allow the histogram to plot on a different scale
            bin_array = np.arange(x1, x2+xdif, b)
            axHist.hist(ages[i], bins=bin_array, color='black', fill=None, alpha=1, histtype='bar', density=False)
            axHist.set_xlim([x1, x2]) # Use this code to set the x-axis scale
            if normPlots:
                histMax = 0
                for k in range(len(sampleList)):
                    if max(np.histogram(ages[k], bins=bin_array)[0]) > histMax:
                        histMax = max(np.histogram(ages[k], bins=bin_array)[0]) 
                axHist.set_ylim([0, histMax])
            axHist.get_yaxis().set_visible(True) # This makes the y-axis numbers invisible
            if (plotPDP or plotKDE):
                axHist.get_xaxis().set_visible(False) # Do not plot the x-axis if it has already been plotted
            
        # Plot colored vertical bars, if selected
        if plotColorBar:
            if len(np.shape(agebins)) == 1:
                for j in range(nage):
                    axs[c+1,0].axvspan(xmin=agebins[j],xmax=agebins[j+1], color = agebinsc[j])
            if len(np.shape(agebins)) == 2:
                for j in range(len(agebins)):
                    axs[c,0].axvspan(xmin=agebins[j][0],xmax=agebins[j][1], color = agebinsc[j])
        c = c+t+l
        
    return fig
    
def ageProportionsCSV(ages, errors, numGrains, labels, agebins, fileName):
    """
    Calculates the number and percentage of analyses in user-specified age bins for each sample or sample group. Exports a CSV file.

    Parameters
    ----------
    ages : array of ages for each sample or sample group. Output from sampleToData()
    errors : array of errors for each sample or sample group. Output from sampleToData()
    numGrains : array of number of analyses per sample or sample group. Output from sampleToData()
    labels : array of labels for each sample or sample group. Output from sampleToData()
    agebins : array of bin edges (Myr)
    fileName : name of the output file that will be created
    
    Notes
    -----
    """ 
    import csv

    if len(np.shape(agebins)) == 1:
        numCategories = len(agebins)-1
    if len(np.shape(agebins)) == 2:
        numCategories = len(agebins)
    agecategory = []
    agecategoryp = []
    rowLabels = ['Sample','numGrains']

    pathlib.Path('Output').mkdir(parents=True, exist_ok=True) # Recursively creates the directory and does not raise an exception if the directory already exists 
    with open(pathlib.Path('Output/') / fileName, 'w', newline='') as f: #Select the CSV file to save data to , 'wb'
        writer = csv.writer(f)

        # Create category names (age ranges)
        for i in range(numCategories):
            if len(np.shape(agebins)) == 1:
                agecategory = np.append(agecategory, ('%d' % agebins[i])+('-%d' % agebins[i+1])+(' Ma'))
            if len(np.shape(agebins)) == 2:
                agecategory = np.append(agecategory, ('%d' % agebins[i][0])+('-%d' % agebins[i][1])+(' Ma'))
            rowLabels = np.append(rowLabels, agecategory[i])
        for i in range(numCategories):
            if len(np.shape(agebins)) == 1:
                agecategoryp = np.append(agecategoryp, ('%d' % agebins[i])+('-%d' % agebins[i+1])+(' Ma (%)'))
            if len(np.shape(agebins)) == 2:
                agecategoryp = np.append(agecategoryp, ('%d' % agebins[i][0])+('-%d' % agebins[i][1])+(' Ma (%)'))
            rowLabels = np.append(rowLabels, agecategoryp[i])

        writer.writerow((rowLabels))

        for i in range(len(ages)):
            data_row = [labels[i], float(numGrains[i])]
            if len(np.shape(agebins)) == 1:
                hist, histp = np.histogram(ages[i], agebins)
                histp = hist/float(numGrains[i])
            if len(np.shape(agebins)) == 2:
                hist = []
                for j in range(len(agebins)):
                    hist.append(np.histogram(ages[i],agebins[j])[0][0])
                histp = hist/np.sum(hist)          
            for j in range(numCategories):
                data_row = np.append(data_row, hist[j])
            for j in range(numCategories):
                data_row = np.append(data_row, histp[j])
            writer.writerow(data_row)

def plotBar(width, height, overlap, main_byid_df, sampleList, ages, numGrains, labels, agebins, agebinsc, separateGroups, savePlot):
    """
    Creates a figure with the age bin proportions are displayed as bar graphs for samples or groups of samples.

    Parameters
    ----------
    width : width of the plot
    height : height of the plot
    overlap : overlap between adjacent bars (a value of 1.0 results in no gap)
    main_byid_df : Detrital database. Output from loadData()    
    sampleList : array of sample IDs.
        Must be in form for individual samples: ['Sample1', 'Sample2', . . . , etc.].
        Must be in the form for groups of samples: [(['Sample1','Sample2', . . . , etc.], 'Group 1 name'),
                                                    (['Sample1','Sample2', . . . , etc.], 'Group 2 name')]
    ages : array of ages for each sample or sample group. Output from sampleToData()
    numGrains : array of number of analyses per sample or sample group. Output from sampleToData()
    labels : array of labels for each sample or sample group. Output from sampleToData()
    agebins : array of bin edges (Myr)
    agebinsc : array of colors that correspond to age bins
    separateGroups : set to True to separate sample groups into individual samples

    Returns
    -------
    fig : a figure age bin proportions displayed as bar graphs
    
    Notes
    -----
    """
    if type(sampleList[0])==tuple:
        if separateGroups:
            for i in range(len(sampleList)): # Loop for each group
                N = len(sampleList[i][0]) # Number of samples in the group
                figBar, ax = plt.subplots(1, figsize=(2*width,N*(height+1-overlap)))           
                y = np.arange(N,0.,-1.)
                for j in range(N): # Loop for each sample
                    ages = sampleToData(sampleList[i][0], main_byid_df)[0]
                    if len(np.shape(agebins)) == 1:
                        hist = np.histogram(ages[j], agebins)[0]
                    if len(np.shape(agebins)) == 2:
                        hist = []
                        for k in range(len(agebins)):
                            hist.append(np.histogram(ages[j],agebins[k])[0][0])
                        hist = np.asarray(hist)
                    hist = hist/len(ages[j])
                    left = 0
                    for k in range(len(hist)): # Loop for each age population category
                        ax.barh(y=y[j],width=hist[k],color=agebinsc[k], left=left, height=overlap)
                        left += hist[k]
                ax.set_xlabel('Proportion')
                ax.set_title(sampleList[i][1]+' age proportions')
                ax.set_ylabel('Sample ID')
                ax.set_yticks(y)
                ax.set_yticklabels((sampleList[i][0])) 
                ax.set_xlim(0.,1.)
                if savePlot:
                    pathlib.Path('Output').mkdir(parents=True, exist_ok=True) # Recursively creates the directory and does not raise an exception if the directory already exists 
                    figBar.savefig(pathlib.Path('Output/') / (('Bar_%s' %sampleList[i][1])+('.pdf')))                      
        else:
            N = len(sampleList) # Number of groups
            figBar, ax = plt.subplots(1, figsize=(2*width,N*(height+1-overlap)))            
            y = np.arange(N,0.,-1.)
            for i in range(len(sampleList)):
                if len(np.shape(agebins)) == 1:
                    hist = np.histogram(ages[i], agebins)[0]
                if len(np.shape(agebins)) == 2:
                    hist = []
                    for j in range(len(agebins)):
                        hist.append(np.histogram(ages[i],agebins[j])[0][0])
                    hist = np.asarray(hist)
                hist = hist/len(ages[i])
                left = 0
                for j in range(len(hist)):
                    ax.barh(y=y[i],width=hist[j],color=agebinsc[j], left=left, height=overlap)
                    left += hist[j]
            ax.set_xlabel('Proportion')
            ax.set_title('Age proportions')
            ax.set_ylabel('Group')
            ax.set_yticks(y)               
            ax.set_yticklabels((labels))             
            ax.set_xlim(0.,1.) 
            if savePlot:
                pathlib.Path('Output').mkdir(parents=True, exist_ok=True) # Recursively creates the directory and does not raise an exception if the directory already exists 
                figBar.savefig(pathlib.Path('Output/') / (('BarGroups')+('.pdf')))    
    else:
        N = len(sampleList) # Number of samples in the list
        figBar, ax = plt.subplots(1, figsize=(2*width,N*height))        
        y = np.arange(N,0.,-1.)
        for i in range(len(sampleList)):
            if len(np.shape(agebins)) == 1:
                hist = np.histogram(ages[i], agebins)[0]
            if len(np.shape(agebins)) == 2:
                hist = []
                for j in range(len(agebins)):
                    hist.append(np.histogram(ages[i],agebins[j])[0][0])
                hist = np.asarray(hist)
            hist = hist/len(ages[i])
            left = 0
            for j in range(len(hist)):
                ax.barh(y=y[i],width=hist[j],color=agebinsc[j], left=left, height=overlap)
                left += hist[j]
        ax.set_xlabel('Proportion')
        ax.set_title('Age proportions')
        ax.set_ylabel('Sample ID')
        ax.set_yticks(y)
        ax.set_yticklabels((labels)) 
        ax.set_xlim(0.,1.)
        if savePlot:
            pathlib.Path('Output').mkdir(parents=True, exist_ok=True) # Recursively creates the directory and does not raise an exception if the directory already exists 
            figBar.savefig(pathlib.Path('Output/') / (('BarSamples')+('.pdf')))    
    return figBar

def plotFoliumMap(sampleList, main_byid_df, ages, errors, numGrains, plotMapKDE, plotMapPDP, plotCumulative, 
    x2, bw, mapType, exportKML, descrpt, stickyPopups = False, width=400, height=100):
    """
    Displays sample locations on an interactive map.

    Parameters
    ----------
    sampleList : array of sample IDs.
        Must be in form for individual samples: ['Sample1', 'Sample2', . . . , etc.].
        Must be in the form for groups of samples: [(['Sample1','Sample2', . . . , etc.], 'Group 1 name'),
                                                    (['Sample1','Sample2', . . . , etc.], 'Group 2 name')]
    main_byid_df : Detrital database. Output from loadData()
    ages : array of ages for each sample or sample group. Output from sampleToData()
    errors : array of errors for each sample or sample group. Output from sampleToData()
    numGrains : array of number of analyses per sample or sample group. Output from sampleToData()
    plotMapKDE : set to True to enable KDE plots when samples are clicked
    plotMapPDP : set to True to enable PDP plots when samples are clicked
    plotCumulative : set to True to plot the KDE and/or PDP as a cumulative distribution when samples are clicked
    x2 : upper limit (Myr) of x-axis    
    bw : KDE bandwidth. Options are 'optimizedFixed', 'optimizedVariable', or a number (bandwidth in Myr)
    mapType : select map background type. Options are 'NatGeo_World_Map', 'World_Street_Map', 'World_Topo_Map', 'World_Light_Gray', 'World_Shaded_Relief', 'World_Terrain_Base', 'World_Hillshade', 'World_Physical_Map'
    
    Returns
    -------
    fig : a map with sample locations
    
    Notes
    -----
    Requires an internet connection
    """ 
    import bisect
    import folium # Must be installed: pip install folium
    import vincent # Must be installed: pip install vincent
    import simplekml # Must be installed: pip install simplekml

    def isNaN(num):
        return num != num

    if exportKML:
        kml = simplekml.Kml()
        style = simplekml.Style()
        style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle_highlight.png'
    
    # Find information about the spatial extent (lat/long) of the dataset
    if type(sampleList[0])==tuple:
        # Number of samples per plotted distribution
        N = np.zeros_like(numGrains)
        numSamples = 0
        for i in range(len(sampleList)):
            N[i] = len(sampleList[i][0])
            numSamples = numSamples + N[i]
        latitudes = np.empty(shape=(numSamples,1))
        longitudes = np.empty(shape=(numSamples,1))
        count = 0
        for i in range(len(sampleList)):
            for j in range(len(sampleList[i][0])):
                if not isNaN(main_byid_df.ix[sampleList[i][0][j],'Latitude']):
                    latitudes[count] = main_byid_df.ix[sampleList[i][0][j],'Latitude']
                    latitudes.astype(float)
                    longitudes[count] = float(main_byid_df.ix[sampleList[i][0][j],'Longitude'])
                    longitudes.astype(float)
                    count = count+1
    else:
        latitudes = np.empty(shape=(len(sampleList),1))
        longitudes = np.empty(shape=(len(sampleList),1))
        for i in range(len(sampleList)):
            if not isNaN(main_byid_df.ix[sampleList[i],'Latitude']):
                latitudes[i] = main_byid_df.ix[sampleList[i],'Latitude']
                latitudes.astype(float)
                longitudes[i] = float(main_byid_df.ix[sampleList[i],'Longitude'])
                longitudes.astype(float)
    latRange = latitudes.max()-latitudes.min()
    longRange = longitudes.max()-longitudes.min()
    print(latRange)
    print(longRange)
    
    latRef = [
        (0, 12),
        (0.1, 11),
        (0.5, 9),
        (1, 8),
        (2.5, 7),
        (5, 6),
        (10, 5),
        (21, 4),
        (60, 3),
        (90, 2),
        (125, 2),
        (360, 2)
        ]
    latRef.sort() # list must be sorted

    longRef = [
        (0, 12),
        (0.1, 10),
        (0.5, 9),
        (1, 8),
        (2.5, 7),
        (5, 6),
        (10, 5),
        (21, 4),
        (45, 3),
        (90, 2),
        (125, 2),
        (360, 2)
        ]
    latRef.sort() # list must be sorted

    latPos = bisect.bisect_right(latRef, (latRange,))
    lat_zoom_start = latRef[latPos][1]
    longPos = bisect.bisect_right(longRef, (longRange,))
    long_zoom_start = longRef[longPos][1]

    print('latitude %s -> %s' % (latRange, latRef[latPos]))
    print('longitude %s -> %s' % (longRange, longRef[longPos]))

    if lat_zoom_start < long_zoom_start:
        latlong_zoom_start = lat_zoom_start
    else:
        latlong_zoom_start = long_zoom_start

    # Basemap options
    if mapType == 'NatGeo_World_Map':
        url_base = 'http://server.arcgisonline.com/ArcGIS/rest/services/'
        service = 'NatGeo_World_Map/MapServer/tile/{z}/{y}/{x}'
    if mapType == 'World_Street_Map':
        url_base = 'http://server.arcgisonline.com/ArcGIS/rest/services/'
        service = 'World_Street_Map/MapServer/tile/{z}/{y}/{x}'
    if mapType == 'World_Topo_Map':
        url_base = 'http://server.arcgisonline.com/ArcGIS/rest/services/'
        service = 'World_Topo_Map/MapServer/tile/{z}/{y}/{x}'        
    if mapType == 'Ocean_Basemap':
        url_base = 'http://server.arcgisonline.com/ArcGIS/rest/services//'
        service = 'Ocean_Basemap/MapServer/tile/{z}/{y}/{x}'
    if mapType == 'World_Physical_Map':
        url_base = 'http://server.arcgisonline.com/ArcGIS/rest/services//'
        service = 'World_Physical_Map/MapServer/tile/{z}/{y}/{x}'        
    if mapType == 'World_Shaded_Relief':
        url_base = 'http://server.arcgisonline.com/ArcGIS/rest/services//'
        service = 'World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}'        
    if mapType == 'World_Terrain_Base':
        url_base = 'http://server.arcgisonline.com/ArcGIS/rest/services//'
        service = 'World_Terrain_Base/MapServer/tile/{z}/{y}/{x}'
    if mapType == 'World_Hillshade':
        url_base = 'http://server.arcgisonline.com/ArcGIS/rest/services/Elevation//'
        service = 'World_Hillshade/MapServer/tile/{z}/{y}/{x}'        
          
    tileset = url_base + service

    m = folium.Map(location=[latitudes.mean(),longitudes.mean()], zoom_start=latlong_zoom_start, control_scale=True, 
        tiles=tileset,attr='NG style') # Another option for a basemap tiles='Stamen Terrain'
#    m.add_child(folium.LatLngPopup()) # Allows lat long to generated by clicking on map

    if type(sampleList[0])==tuple:
        for i in range(len(sampleList)): # Loop for each group
            feature_group = folium.FeatureGroup(name=sampleList[i][1])
            ages = sampleToData(sampleList[i][0], main_byid_df)[0]
            errors = sampleToData(sampleList[i][0], main_byid_df)[1]
            if (plotMapKDE and not plotCumulative):
                if bw == 'optimizedFixed':
                    dist = KDEcalcAgesLocalAdapt(ages=ages, x1=0, x2=x2, xdif=1, cumulative=False)
                if bw == 'optimizedVariable':
                    dist = KDEcalcAgesGlobalAdapt(ages=ages, x1=0, x2=x2, xdif=1, cumulative=False)            
                if type(bw) != str:
                    dist = KDEcalcAges_2(ages, x1=0, x2=x2, xdif=1, bw=bw)
            if (plotMapKDE and plotCumulative):
                dist = KDEcalcAges_2(ages, x1=0, x2=x2, bw=bw, xdif=1, cumulative=True)
            if (plotMapPDP and not plotCumulative):
                dist = PDPcalcAges(ages, errors, x1=0, x2=x2, xdif=1)
            if (plotMapPDP and plotCumulative):
                dist = PDPcalcAges(ages, errors, x1=0, x2=x2, xdif=1, cumulative=True)
            if (not plotMapKDE and not plotMapPDP and plotCumulative):
                dist = CDFcalcAges(ages, x1=0, x2=x2, xdif=1)        
            for j in range(len(sampleList[i][0])): # Loop for each sample within each group
                if not isNaN(main_byid_df.ix[sampleList[i][0][j],'Latitude']):
                    if ((plotMapKDE or plotMapPDP) and not plotCumulative):
                        distArea = vincent.Area(dist[1][j].tolist(), width=width, height=height)
                        distArea.axis_titles(x='Age (Ma)', y='')
                        distArea.legend(title=sampleList[i][0][j])
                        if stickyPopups:
                            popup = folium.Popup(max_width=600, sticky=True)
                        else:
                            popup = folium.Popup(max_width=600)
                        folium.Vega(distArea, height=int(height*1.5), width=int(width*1.25)).add_to(popup)
                        folium.RegularPolygonMarker([main_byid_df.ix[sampleList[i][0][j],'Latitude'],main_byid_df.ix[sampleList[i][0][j],
                                                'Longitude']], fill_color=colorMe(i), radius=6, popup=popup).add_to(feature_group)
                    if (plotCumulative):
                        distLine = vincent.Line(dist[1][j].tolist(), width=width, height=height)
                        distLine.axis_titles(x='Age (Ma)', y='')
                        distLine.legend(title=sampleList[i][0][j])
                        if stickyPopups:
                            popup = folium.Popup(max_width=600, sticky=True)
                        else:
                            popup = folium.Popup(max_width=600)
                        folium.Vega(distLine, height=int(height*1.5), width=int(width*1.25)).add_to(popup)
                        folium.RegularPolygonMarker([main_byid_df.ix[sampleList[i][0][j],'Latitude'],main_byid_df.ix[sampleList[i][0][j],
                                                'Longitude']], fill_color=colorMe(i), radius=6, popup=popup).add_to(feature_group)
                    if (not (plotMapKDE or plotMapPDP or plotCumulative)):
                        folium.RegularPolygonMarker([main_byid_df.ix[sampleList[i][0][j],'Latitude'],main_byid_df.ix[sampleList[i][0][j],
                                                'Longitude']], fill_color=colorMe(i), radius=6, popup=sampleList[i][0][j]).add_to(feature_group)
                    feature_group.add_to(m)
                    if exportKML:
                        pnt = kml.newpoint(name=sampleList[i][0][j], description=main_byid_df.ix[sampleList[i][0][j],descrpt],coords=[(main_byid_df.ix[sampleList[i][0][j],'Longitude'],main_byid_df.ix[sampleList[i][0][j],'Latitude'])])
                        pnt.style = style
        folium.LayerControl().add_to(m)
    else:
        if (plotMapKDE and not plotCumulative):
            if bw == 'optimizedFixed':
                dist = KDEcalcAgesLocalAdapt(ages=ages, x1=0, x2=x2, xdif=1, cumulative=False)
            if bw == 'optimizedVariable':
                dist = KDEcalcAgesGlobalAdapt(ages=ages, x1=0, x2=x2, xdif=1, cumulative=False)            
            if type(bw) != str:
                dist = KDEcalcAges_2(ages, x1=0, x2=x2, xdif=1, bw=bw)
        if (plotMapKDE and plotCumulative):
            if bw == 'optimizedFixed':
                dist = KDEcalcAgesLocalAdapt(ages=ages, x1=0, x2=x2, xdif=1, cumulative=True)
            if bw == 'optimizedVariable':
                dist = KDEcalcAgesGlobalAdapt(ages=ages, x1=0, x2=x2, xdif=1, cumulative=True)               
            if type(bw) != str:
                dist = KDEcalcAges_2(ages, x1=0, x2=x2, bw=bw, xdif=1, cumulative=True)
        if (plotMapPDP and not plotCumulative):
            dist = PDPcalcAges(ages, errors, x1=0, x2=x2, xdif=1, cumulative=False)
        if (plotMapPDP and plotCumulative):
            dist = PDPcalcAges(ages, errors, x1=0, x2=x2, xdif=1, cumulative=True)
        if ((not plotMapKDE and not plotMapPDP) and plotCumulative):
            dist = CDFcalcAges(ages, x1=0, x2=x2, xdif=1)
        feature_group = folium.FeatureGroup(name='Samples')
        for i in range(len(sampleList)):
            if not isNaN(main_byid_df.ix[sampleList[i],'Latitude']):
                if ((plotMapKDE or plotMapPDP) and not plotCumulative):
                    distArea = vincent.Area(dist[1][i].tolist(), width=width, height=height)
                    distArea.axis_titles(x='Age (Ma)', y='')
                    distArea.legend(title=sampleList[i])
                    if stickyPopups:
                        popup = folium.Popup(max_width=600, sticky=True)
                    else:
                        popup = folium.Popup(max_width=600)
                    folium.Vega(distArea, height=int(height*1.5), width=int(width*1.25)).add_to(popup)
                    folium.RegularPolygonMarker([main_byid_df.ix[sampleList[i],'Latitude'],main_byid_df.ix[sampleList[i],
                                            'Longitude']], fill_color=colorMe(0), radius=6, popup=popup).add_to(feature_group)
                if (plotCumulative):
                    distLine = vincent.Line(dist[1][i].tolist(), width=width, height=height)
                    distLine.axis_titles(x='Age (Ma)', y='')
                    distLine.legend(title=sampleList[i])
                    if stickyPopups:
                        popup = folium.Popup(max_width=600, sticky=True)
                    else:
                        popup = folium.Popup(max_width=600)
                    folium.Vega(distLine, height=int(height*1.5), width=int(width*1.25)).add_to(popup)
                    folium.RegularPolygonMarker([main_byid_df.ix[sampleList[i],'Latitude'],main_byid_df.ix[sampleList[i],
                                            'Longitude']], fill_color=colorMe(0), radius=6, popup=popup).add_to(feature_group)
                if (not (plotMapKDE or plotMapPDP or plotCumulative)):
                    folium.RegularPolygonMarker([main_byid_df.ix[sampleList[i],'Latitude'],main_byid_df.ix[sampleList[i],
                                            'Longitude']], fill_color=colorMe(0), radius=6, popup=sampleList[i]).add_to(feature_group)
                if exportKML:
                    pnt = kml.newpoint(name=sampleList[i], description=main_byid_df.ix[sampleList[i],descrpt],coords=[(main_byid_df.ix[sampleList[i],'Longitude'],main_byid_df.ix[sampleList[i],'Latitude'])])
                    pnt.style = style
        feature_group.add_to(m)
        folium.LayerControl().add_to(m)
    if exportKML:
        pathlib.Path('Output').mkdir(parents=True, exist_ok=True) # Recursively creates the directory and does not raise an exception if the directory already exists 
        ####
        kml.save(str(pathlib.Path('Output/') / (('samples')+('.kml'))))
    return m

def MDAtoCSV(sampleList, ages, errors, numGrains, labels, fileName, sortBy, barWidth, plotWidth, plotHeight, ageColors, alpha, makePlot, fillMDACalcs):
    """
    Calculate the maxmium depositional age (MDA) for a sample or group of samples. Results will be exported to a CSV file. Individual plots can be made for each sample or group of samples showing the youngest grains and different calculations of the maxmimum depositional age.

    Parameters
    ----------
    sampleList : array of sample IDs.
        Must be in form for individual samples: ['Sample1', 'Sample2', . . . , etc.].
        Must be in the form for groups of samples: [(['Sample1','Sample2', . . . , etc.], 'Group 1 name'),
                                                    (['Sample1','Sample2', . . . , etc.], 'Group 2 name')]
    ages : array of ages for each sample or sample group. Output from sampleToData()
    errors : array of errors for each sample or sample group. Output from sampleToData()
    numGrains : array of number of analyses per sample or sample group. Output from sampleToData()
    labels : array of labels for each sample or sample group. Output from sampleToData()
    fileName : specify the output CSV file name
    sortBy : specify how grains are sorted: by best age, best age + 1 sigma error, or best age + 2 sigma error. Options: 'mean', '1sigma', '2sigma'
    barWidth : specify the width of bars
    plotWidth : specify the width of the plot
    plotHeight : specify the height of the plot
    ageColors : an array of colors of the horizontal bars for YSG, YC1S(2+), and YC2s(3+), respectively
    alpha : specify the transparency of filled 1 sigma confidence limits  
    makePlot : set to True to create a plot of youngest grain analyses and MDA calculations
    fillMDACalc : set to True to shade each MDA calculation within 1 sigma confidence limits
    sigma : (optional) Specify whether bestAgeErr are 1-sigma or 2-sigma errors. Default is '1sigma', but '2sigma' can also be specified.

    Returns
    -------
    fig : (optional) a figure with youngest grain analyses and MDA calculations
    
    Notes
    -----
    """ 
    import csv
    import matplotlib.pyplot as plt
    
    pathlib.Path('Output').mkdir(parents=True, exist_ok=True) # Recursively creates the directory and does not raise an exception if the directory already exists 
    #####
    with open(str(pathlib.Path('Output/')) + fileName, 'w', newline='') as f: #Select the CSV file to save data to
        writer = csv.writer(f)
        writer.writerow(('Sample','N', 'YSG', 'YSG_err1s', 'YC1S WM', 'YC1S WM 1serr', 'YC1S WM MSWD', 'YC1S cluster size', 'YC2S WM', 'YC2S WM 1serr', 'YC2S WM MSWD', 'YC2S cluster size'))
        N = len(sampleList) # Number of samples or groups of samples
        if makePlot:
            figMDA, ax = plt.subplots(N,1, figsize=(plotWidth,N*plotHeight)) #figsize=(width,height)
            for i in range(N):
                if N > 1:
                    ax[i] = plt.subplot2grid((N,1),(i,0))
                else:
                    break
        for i in range(N):
            if makePlot:
                if N > 1:
                    ax_i = ax[i]
                else:
                    ax_i = ax
            error1s = errors[i]
            error2s = errors[i]*2
            data_err1s_ageSort = list(zip(ages[i], error1s))            
            data_err1s = list(zip(ages[i], error1s))
            data_err2s = list(zip(ages[i], error2s))
            data_err1s_ageSort.sort(key=lambda d: d[0]) # Sort based on age
            data_err1s.sort(key=lambda d: d[0] + d[1]) # Sort based on age + 1s error
            data_err2s.sort(key=lambda d: d[0] + d[1]) # Sort based on age + 2s error
            YSG = data_err1s[0][0]
            YSG_err1s = data_err1s[0][1]

            def find_youngest_cluster_1s(min_cluster_size):
                i_min = 0
                i_max = 0
                for i in range(1, len(data_err1s)):
                    top = data_err1s[i_min][0] + data_err1s[i_min][1]
                    bottom = data_err1s[i][0] - data_err1s[i][1]
                    if (top >= bottom):
                        i_max = i
                    elif i_max - i_min + 1 >= min_cluster_size:
                        break
                    else:
                        i_min = i
                        i_max = i
                return data_err1s[i_min: i_max + 1] if i_min < i_max else [], i_max

            def find_youngest_cluster_2s(min_cluster_size):
                i_min = 0
                i_max = 0
                for i in range(1, len(data_err1s)):
                    top = data_err2s[i_min][0] + data_err2s[i_min][1]
                    bottom = data_err2s[i][0] - data_err2s[i][1]
                    if (top >= bottom):
                        i_max = i
                    elif i_max - i_min + 1 >= min_cluster_size:
                        break
                    else:
                        i_min = i
                        i_max = i
                return data_err2s[i_min: i_max + 1] if i_min < i_max else [], i_max

            YC1S, YC1S_imax = find_youngest_cluster_1s(2)
            YC1S_WM = weightedMean(np.array([d[0] for d in YC1S]), np.array([d[1] for d in YC1S]))
            YC2S, YC2S_imax = find_youngest_cluster_2s(3)
            YC2S_WM = weightedMean(np.array([d[0] for d in YC2S]), np.array([d[1] for d in YC2S])/2)

            writer.writerow((labels[i],numGrains[i], YSG, YSG_err1s, YC1S_WM[0], YC1S_WM[1]/2, YC1S_WM[2], len(YC1S), YC2S_WM[0], YC2S_WM[1]/2, YC2S_WM[2], len(YC2S)))
            
            if makePlot:
                # Specify plot characteristics and axes    
                if N > 1:
                    ax_i.set_ylabel("Age (Ma)")
                    ax_i.get_xaxis().set_ticks([])
                else:
                    ax.set_ylabel("Age (Ma)")
                    ax.get_xaxis().set_ticks([])
                
                # Choose which way to sort the ages
                if sortBy == 'mean':
                    ageErrors = data_err1s_ageSort
                if sortBy == '1sigma':
                    ageErrors = data_err1s
                if sortBy == '2sigma':
                    ageErrors = data_err2s
                
                # Determine how many grains to plot per sample or sample group
                toPlot = np.zeros_like(np.empty(shape=(len(ageErrors),1)))
                YC1S_max = data_err1s[YC1S_imax][0]+data_err1s[YC1S_imax][1]
                YC2S_max = data_err2s[YC2S_imax][0]+data_err2s[YC2S_imax][1]
                if YC1S_max >= YC2S_max:
                    plotMax = YC1S_max
                else:
                    plotMax = YC2S_max
                for j in range(len(ageErrors)):
                    if (sortBy == 'mean' or sortBy == '1sigma'):
                        if ageErrors[j][0]+ageErrors[j][1]*2 > plotMax:
                            toPlot[j] = 0
                        else:
                            toPlot[j] = 1
                    if sortBy == '2sigma':
                        if ageErrors[j][0]+ageErrors[j][1] > plotMax:
                            toPlot[j] = 0
                        else:
                            toPlot[j] = 1                    
                count = np.max(np.nonzero(toPlot == 1))+1
                for j in range(count):
                    barx = np.arange(0,count,1.)
                    barymin = ageErrors[j][0]-ageErrors[j][1]
                    barymax = ageErrors[j][0]+ageErrors[j][1]
                    height = barymax - barymin
                    gridlines = ax_i.get_ygridlines()
                    for line in gridlines:
                        line.set_linestyle('-')
                    if fillMDACalcs:
                    # Plot the age span (age +- 1 sigma error) of the three MDA calculations
                        ax_i.axhspan(ymin=YSG-YSG_err1s, ymax=YSG+YSG_err1s, xmin=np.min(barx), xmax=np.max(barx), color = ageColors[0], alpha=alpha, zorder=1)                
                        ax_i.axhspan(ymin=YC1S_WM[0]-YC1S_WM[1]/2, ymax=YC1S_WM[0]+YC1S_WM[1]/2, xmin=np.min(barx), xmax=np.max(barx), color = ageColors[1], alpha=alpha, zorder=2)
                        ax_i.axhspan(ymin=YC2S_WM[0]-YC2S_WM[1]/2, ymax=YC2S_WM[0]+YC2S_WM[1]/2, xmin=np.min(barx), xmax=np.max(barx), color = ageColors[2], alpha=alpha, zorder=3)
                    if not fillMDACalcs:
                        # Plot the mean and 1 sigma error limits of the three MDA calculations
                        ax_i.hlines(y=YSG, xmin=np.min(barx), xmax=np.max(barx), color = ageColors[0], lw=1)
                        ax_i.hlines(y=YSG+YSG_err1s, xmin=np.min(barx), xmax=np.max(barx), color = ageColors[0], lw=1, linestyles='dashed')
                        ax_i.hlines(y=YSG-YSG_err1s, xmin= np.min(barx), xmax = np.max(barx), color = ageColors[0], lw=1, linestyles='dashed')
                        ax_i.hlines(y=YC1S_WM[0], xmin= np.min(barx), xmax = np.max(barx), color = ageColors[1], lw=1)
                        ax_i.hlines(y=YC1S_WM[0]+YC1S_WM[1]/2, xmin= np.min(barx), xmax = np.max(barx), color = ageColors[1], lw=1, linestyles='dashed')
                        ax_i.hlines(y=YC1S_WM[0]-YC1S_WM[1]/2, xmin= np.min(barx), xmax = np.max(barx), color = ageColors[1], lw=1, linestyles='dashed')      
                        ax_i.hlines(y=YC2S_WM[0], xmin= np.min(barx), xmax = np.max(barx), color = ageColors[2], lw=1)
                        ax_i.hlines(y=YC2S_WM[0]+YC2S_WM[1]/2, xmin= np.min(barx), xmax = np.max(barx), color = ageColors[2], lw=1, linestyles='dashed')
                        ax_i.hlines(y=YC2S_WM[0]-YC2S_WM[1]/2, xmin= np.min(barx), xmax = np.max(barx), color = ageColors[2], lw=1, linestyles='dashed')     
                    # Plot the mean age, 1 sigma error, and 2 sigma error for each grain
                    ax_i.hlines(y = ageErrors[j][0], xmin=barx[j]-barWidth/2, xmax=barx[j]+barWidth/2, color='white', zorder=11)          
                    if (sortBy == 'mean' or sortBy == '1sigma'):
                        ax_i.bar(x = barx[j], height = height, width=barWidth, bottom = ageErrors[j][0]-ageErrors[j][1], color='black', zorder=10)
                        ax_i.bar(x = barx[j], height = height*2, width=barWidth, bottom = ageErrors[j][0]-ageErrors[j][1]*2, color='gray', zorder=9)
                    if sortBy == '2sigma':
                        ax_i.bar(x = barx[j], height = height/2, width=barWidth, bottom = ageErrors[j][0]-ageErrors[j][1]/2, color='black', zorder=10)
                        ax_i.bar(x = barx[j], height = height, width=barWidth, bottom = ageErrors[j][0]-ageErrors[j][1], color='gray', zorder=9)                    
                ax_i.grid(True)
                xmin, xmax = ax_i.get_xlim()
                ymin, ymax = ax_i.get_ylim()
                xdif = xmax-xmin
                ydif = ymax-ymin
                ax_i.text(xmin+xdif*0.02, ymax+ydif*0.02, s=labels[i])
                ax_i.text(xmin+xdif*0.02, ymax+ydif*-0.07, s=("YSG: "+str(round(YSG,1))+"$\pm$"+str(round(YSG_err1s,2))+" Ma"), fontsize='small', backgroundcolor='white')
                ax_i.text(xmin+xdif*0.02, ymax+ydif*-0.14, s=("YC1"+r'$\sigma$'+"(2+): "+str(round(YC1S_WM[0],1))+"$\pm$"+str(round(YC1S_WM[1]/2,2))+" Ma"), fontsize='small', backgroundcolor='white')
                ax_i.text(xmin+xdif*0.02, ymax+ydif*-0.21, s=("YC2"+r'$\sigma$'+"(3+): "+str(round(YC2S_WM[0],1))+"$\pm$"+str(round(YC2S_WM[1]/2,2))+" Ma"), fontsize='small', backgroundcolor='white')
                ax_i.text(xmin+xdif*plotWidth*0.05, ymax+ydif*0.02, s=("Calculated uncertainties are reported at the 1"+r'$\sigma$'+" level"), fontsize='small')

        if makePlot:
            return figMDA
                
def MDS(ages, errors, labels, sampleList, metric=False, plotWidth='10', plotHeight='8', plotPie=False, pieSize=0.05, agebins=None, agebinsc=None, criteria='Dmax', bw='optimizedFixed', color='Default', main_byid_df=None, plotLabels=True):
    """
    Create a multi-dimensional scaling (MDS) plot for individual samples or groups of samples.

    Parameters
    ----------
    ages : array of ages for each sample or sample group. Output from sampleToData()
    labels : array of labels for each sample or sample group. Output from sampleToData()
    sampleList : array of sample IDs.
        Must be in form for individual samples: ['Sample1', 'Sample2', . . . , etc.].
        Must be in the form for groups of samples: [(['Sample1','Sample2', . . . , etc.], 'Group 1 name'),
                                                    (['Sample1','Sample2', . . . , etc.], 'Group 2 name')]
    metric : set to False for non-metric MDS
    plotWidth : specify the width of the plot
    plotHeight : specify the height of the plot
    plotPie : set to True to plot data points as pies
    pieSize : specify the size of pie plots
    agebins : array of bin edges in Myr. Format option 1: [age1, age2, age3, etc.]. Format option 2: [[bin1_min, bin1_max],[bin2_min, bin2_max],etc.]
    agebinsc : array of colors that correspond to age bins
    criteria : (optional) similiarty metric used in the MDS calculation. Options: 'Dmax' (default), 'Vmax', 'R2-PDP', 'R2-KDE'
    bw : (optional) KDE bandwidth. Options are 'optimizedFixed', 'optimizedVariable', or a number (bandwidth in Myr)
    color : (optional) if set to equal a column name in Samples, will color by this category
    main_byid_df : (optional) required if color <> 'Default'. Set equal to main_byid_df.

    Returns
    -------
    fig : a figure with a MDS plot
    
    Notes
    -----
    """     
    from scipy import stats
    from sklearn import manifold
    import math

    matrix = np.empty(shape=(len(ages),len(ages)))

    figMDS, ax = plt.subplots(1, figsize=(plotWidth,plotHeight))

    # Calculate the CDF matrix
    if criteria == 'Dmax' or criteria == 'Vmax':
        CDF = CDFcalcAges(ages)[1]
    if criteria == 'R2-PDP':
        CDF = PDPcalcAges(ages=ages, errors=errors, x1=0, x2=4500, xdif=1, cumulative=True)[1]
    if criteria == 'R2-KDE':
        if bw == 'optimizedFixed':
            CDF = KDEcalcAgesLocalAdapt(ages=ages, x1=0, x2=4500, xdif=1, cumulative=True)[1]
        if bw == 'optimizedVariable':
            CDF = KDEcalcAgesGlobalAdapt(ages=ages, x1=0, x2=4500, xdif=1, cumulative=True)[1]
        if type(bw) != str:
            CDF = KDEcalcAges_2(ages=ages, x1=0, x2=4500, xdif=1, bw=bw, cumulative=True)[1]
    for i in range(len(ages)):
        for j in range(len(ages)):
            if criteria == 'Dmax':
                matrix[i,j] = stats.ks_2samp(ages[i],ages[j])[0]
            if criteria == 'Vmax':
                matrix[i,j] = calcVmax(CDF[i], CDF[j])
            if criteria == 'R2-PDP' or criteria == 'R2-KDE':
                matrix[i,j] = calcComplR2(CDF[i], CDF[j])     
    mds = manifold.MDS(random_state=1, dissimilarity='precomputed', n_init=1)
    pos = mds.fit(matrix).embedding_
    posStress = mds.fit(matrix).stress_     
    nmds = manifold.MDS(metric=False, random_state=1, dissimilarity='precomputed', n_init=1)
    npos = nmds.fit_transform(matrix, init=pos)
    nposStress = mds.fit(matrix).stress_ 
    if metric:
        m = pos
        stress = posStress
    else:
        m = npos
        stress = nposStress

    # For coloring by category
    if color != 'Default':
        values = list(set(main_byid_df[color]))
        dicts = {}
        c = 0
        for value in values:
            dicts[value] = colorMe(c)
            c += 1

    for i in range(len(npos)):
        if plotPie:
            if len(np.shape(agebins)) == 1:
                hist = np.histogram(ages[i],agebins)[0]
                histP = np.cumsum([0]+list(hist/np.sum(hist)))
                for j in range(len(hist)): # One loop for each bin
                    x = [0] + np.cos(np.linspace(2*math.pi*histP[j], 2*math.pi*histP[j+1], 100)).tolist()
                    y = [0] + np.sin(np.linspace(2*math.pi*histP[j], 2*math.pi*histP[j+1], 100)).tolist()
                    ax.fill(np.array(x)*pieSize+m[i][0],np.array(y)*pieSize+m[i][1],facecolor=agebinsc[j])
            if len(np.shape(agebins)) == 2:
                hist = [0]
                for j in range(len(agebins)):
                    hist.append(np.histogram(ages[i],agebins[j])[0][0])
                histP = np.cumsum(list(hist/np.sum(hist)))
                for j in range(len(hist)-1): # One loop for each bin
                        x = [0] + np.cos(np.linspace(2*math.pi*histP[j], 2*math.pi*histP[j+1], 100)).tolist()
                        y = [0] + np.sin(np.linspace(2*math.pi*histP[j], 2*math.pi*histP[j+1], 100)).tolist()
                        ax.fill(np.array(x)*pieSize+m[i][0],np.array(y)*pieSize+m[i][1],facecolor=agebinsc[j])

            if plotLabels:
                ax.text(m[i][0]+pieSize/1.5,m[i][1]+pieSize/1.5,labels[i])
            ax.set_aspect('equal')
        else:
            if color == 'Default':
                ax.plot(m[i][0],m[i][1],'o',label=sampleList[i],color=colorMe(i))
                if plotLabels:
                    ax.text(m[i][0]+0.01,m[i][1]+0.01,labels[i])
            else:
                ax.plot(m[i][0],m[i][1],'o',label=sampleList[i],color=dicts[main_byid_df.loc[sampleList[i],color]])
                if plotLabels:
                    ax.text(m[i][0]+0.01,m[i][1]+0.01,labels[i])

    return figMDS, stress

def plotDoubleDating(main_byid_df, sampleList, x1, x2, y1, y2, plotKDE, colorKDE, colorKDEbyAge, plotPDP, colorPDP, colorPDPbyAge, plotHist, b, bw, xdif, width, height, savePlot, agebins, agebinsc, coolingAge='ZHe_Age', coolingAgeErr='ZHe_Age_err'):
    """
    Creates a figure where detrital cooling ages are plotted against detrital crystallization ages.

    Parameters
    ----------
    main_byid_df : the database indexed by sample name
    sampleList : array of sample IDs.
        Must be in form for individual samples: ['Sample1', 'Sample2', . . . , etc.].
        Must be in the form for groups of samples: [(['Sample1','Sample2', . . . , etc.], 'Group 1 name'),
                                                    (['Sample1','Sample2', . . . , etc.], 'Group 2 name')]
    x1 : lower limit (Myr) of x-axis
    x2 : upper limit (Myr) of x-axis
    y1 : lower limit (Myr) of y-axis
    y2 : upper limit (Myr) of y-axis
    plotKDE : set to True to plot a KDE
    colorKDE : set to True to color the KDE according to same coloration as used in CDF plotting
    colorKDEbyAge : set to True to color the KDE according to user-specified age populations
    plotPDP : set to True to plot a PDP
    colorPDP : set to True to color the PDP according to same coloration as used in CDF plotting
    colorPDPbyAge : set to True to color the PDP according to user-specified age populations
    plotHist : set to True to plot a histogram
    b : histogram bin size (Myr)
    bw : KDE bandwidth. Options are 'optimizedFixed', 'optimizedVariable', or a number (bandwidth in Myr)
    xdif : interval (Myr) over which distributions are calculated
    width : width of the plot
    height : height of the plot
    savePlot : set to True to save the plot as a PDF file
    agebins : array of bin edges in Myr. Format option 1: [age1, age2, age3, etc.]. Format option 2: [[bin1_min, bin1_max],[bin2_min, bin2_max],etc.]
    agebinsc : array of colors that correspond to age bins
    coolingAge : (optional) label for column with cooling ages. Default: 'ZHe_Age'
    coolingAgeErr : (optional) label for column with cooling age errors. Default: 'ZHe_Age_Err'

    Returns
    -------
    fig : the figure
    
    Notes
    -----
    """
    from matplotlib.path import Path
    import matplotlib.patches as patches

    # Assign x and y data (ages and errors)
    for i in range(len(sampleList)):
        x, xerr, xNumGrains, xLabels = sampleToData(sampleList, main_byid_df)
        y, yerr, yNumGrains, yLabels = sampleToData(sampleList, main_byid_df, bestAge=coolingAge, bestAgeErr=coolingAgeErr)
        
        # Filter out analyses with no (U-Th)/He data
        def isNaN(num):
            return num != num
        
        def ageFilter(list1, list2):
            zipped = zip(list1, list2)
            filtered = filter(lambda pair: not isNaN(pair[1]), zipped)
            unzipped = (list(x) for x in zip(*filtered))
            return unzipped
                              
        xF, yF = ageFilter(x[i],y[i])
        xerrF, yerrF = ageFilter(xerr[i], yerr[i])

        # Assign parameters for plotting shading
        verts = [
                (x1, x1), #left, bottom
                (x1, y2), #left, top
                (x2, y2),
                (x2, x2), #right, top
                (x1, x1)  # ignored
                ]
        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY
                 ]
        path = Path(verts, codes)

        # Specify plot characteristics and axes
        fig, axs = plt.subplots(5,5, figsize=(width,height))

        # Scatterplot characteristics
        axs[0,0] = plt.subplot2grid((5,5),(0,0),rowspan=4, colspan=4) # scatterplot subplot
        axs[0,0].xaxis.set_ticks_position('top')
        axs[0,0].plot([x1,x2],[x1,x2], color='black', lw=1)
        axs[0,0].errorbar(x = xF, y = yF, xerr = xerrF, yerr = yerrF, fmt='s', color='black', ecolor='gray', capthick=2)
        patch = patches.PathPatch(path, facecolor='lightgray', lw=0)
        axs[0,0].add_patch(patch)
        axs[0,0].xaxis.set_label_position('top') 
        axs[0,0].set_xlim(x1, x2)
        axs[0,0].set_ylim(y1, y2)

        # (U-Th)/He plot characteristics
        axs[0,1] = plt.subplot2grid((5,5),(0,4), rowspan=4) # panel for ZHe distribution
        axs[0,1].get_xaxis().set_ticks([])
        axs[0,1].set_ylabel('Zircon (U-Th)/He Age (Ma)')
        axs[0,1].yaxis.set_label_position('right')    
        axs[0,1].yaxis.set_ticks_position('right')

        # U-Pb plot characteristics
        axs[1,0] = plt.subplot2grid((5,5),(4,0), colspan=4) # panel for the ZUPb distribution
        axs[1,0].get_yaxis().set_ticks([])
        axs[1,0].set_xlabel('Zircon U-Pb Age (Ma)')

        # Emply plot in the lower right corner
        axs[1,1] = plt.subplot2grid((5,5),(4,4)) # empty corner
        axs[1,1].axis('off')

        # Zircon U-Pb age distribution
        if plotHist:
            bin_array = np.arange(x1, x2+xdif, b)
            axsHist = axs[1,0].twinx()
            axsHist.hist([xF], bins=bin_array, color='black', fill=None, alpha=1, histtype='bar', density=False)
            axsHist.set_xlim(x1, x2)
            axsHist.set_xlabel('Zircon U-Pb Age (Ma)')
        if plotPDP:
            xPDP = PDPcalcAges(ages=[xF],errors=[xerrF], xdif=xdif)
            axsPDP = axs[1,0].twinx()
            axsPDP.plot(xPDP[0],xPDP[1][0], color='black', lw=1)
            axsPDP.set_xlim(x1, x2)
            axsPDP.set_ylim(0,np.max(xPDP[1][0])*1.1)
            axsPDP.get_yaxis().set_ticks([])
            if colorPDP:
                axsPDP.fill_between(xPDP[0], xPDP[1][0], alpha = 1, color=colorMe(i))

            if colorPDPbyAge:
                if len(np.shape(agebins)) == 1:
                    nage = len(agebins)-1
                    for j in range(nage):                
                        xage1 = agebins[j]
                        xage2 = agebins[j+1]
                        PDPpart = xPDP[1][0][int(xage1/xdif):int((xage2+xdif)/xdif)]
                        PDP_agePart = np.arange(xage1, xage1+len(PDPpart), xdif)
                        axsPDP.fill_between(PDP_agePart, 0, PDPpart, alpha = 1, color=agebinsc[j])
                if len(np.shape(agebins)) == 2:
                    for j in range(len(agebins)):
                        xage1 = agebins[j][0]
                        xage2 = agebins[j][1]
                        PDPpart = xPDP[1][0][int(xage1/xdif):int((xage2+xdif)/xdif)]
                        PDP_agePart = np.arange(xage1, xage2+xdif, xdif)
                        PDP_agePart = np.arange(xage1, xage1+len(PDPpart), xdif)
                        axsPDP.fill_between(PDP_agePart, 0, PDPpart, alpha = 1, color=agebinsc[j])
        if plotKDE:
            if bw == 'optimizedFixed':
                xKDE = KDEcalcAgesLocalAdapt(ages=x, xdif=xdif)
            if bw == 'optimizedVariable':
                xKDE = KDEcalcAgesGlobalAdapt(ages=x, xdif=xdif)
            if type(bw) != str:
                xKDE = KDEcalcAges_2(ages=x, xdif=xdif, bw=bw)
            #xKDE = KDEcalcAges_2(ages=x, bw=bw, xdif=xdif)
            axsKDE = axs[1,0].twinx()
            axsKDE.plot(xKDE[0],xKDE[1][0], color='black', lw=1)
            axsKDE.set_xlim(x1, x2)
            axsKDE.set_ylim(0,np.max(xKDE[1][0])*1.1)
            axsKDE.get_yaxis().set_ticks([])
            if colorKDE:
                axsKDE.fill_between(xKDE[0], xKDE[1][0], alpha = 1, color=colorMe(i))
            if colorKDEbyAge:
                if len(np.shape(agebins)) == 1:
                    nage = len(agebins)-1
                    for j in range(nage):
                        xage1 = agebins[j]
                        xage2 = agebins[j+1] 
                        KDEpart = xKDE[1][0][int(xage1/xdif):int((xage2+xdif)/xdif)]
                        KDE_agePart = np.arange(xage1, xage1+len(KDEpart), xdif)
                        axsKDE.fill_between(KDE_agePart, 0, KDEpart, alpha = 1, color=agebinsc[j])
                if len(np.shape(agebins)) == 2:
                     for j in range(len(agebins)):
                        xage1 = agebins[j][0]
                        xage2 = agebins[j][1]
                        KDEpart = xKDE[1][0][int(xage1/xdif):int((xage2+xdif)/xdif)]
                        KDE_agePart = np.arange(xage1, xage2+xdif, xdif)
                        KDE_agePart = np.arange(xage1, xage1+len(KDEpart), xdif)
                        axsKDE.fill_between(KDE_agePart, 0, KDEpart, alpha = 1, color=agebinsc[j])                   

        # Zircon (U-Th)/He age distribution   
        if plotHist:
            axsHist = axs[0,1].twiny()
            axsHist.invert_xaxis()
            bin_array = np.arange(x1, x2+xdif, b)
            axsHist.hist(yF, orientation='horizontal', bins=bin_array, color='black', fill=None, alpha=1, histtype='bar', density=False)
            axsHist.set_ylim(y1, y2)
            axsHist.set_ylabel('Zircon (U-Th)/He Age (Ma)')
        if plotPDP:
            yPDP = PDPcalcAges(ages=[yF],errors=[yerrF], xdif=xdif)
            axsPDP = axs[0,1].twiny()
            axsPDP.plot(yPDP[1][0],yPDP[0], color='black', lw=1)
            axsPDP.set_ylim(y1, y2)
            axsPDP.set_xlim(np.max(yPDP[1][0])*1.1,0)
            axsPDP.get_xaxis().set_ticks([])
            if colorPDP:
                axsPDP.fill_betweenx(yPDP[0], 0, yPDP[1][0], alpha = 1, color=colorMe(i))

            if colorPDPbyAge:
                if len(np.shape(agebins)) == 1:
                    nage = len(agebins)-1
                    for j in range(nage):                
                        xage1 = agebins[j]
                        xage2 = agebins[j+1]
                        PDPpart = yPDP[1][0][int(xage1/xdif):int((xage2+xdif)/xdif)]
                        PDP_agePart = np.arange(xage1, xage1+len(PDPpart), xdif)
                        axsPDP.fill_betweenx(PDP_agePart, 0, PDPpart, alpha = 1, color=agebinsc[j])
                if len(np.shape(agebins)) == 2:
                    for j in range(len(agebins)):
                        xage1 = agebins[j][0]
                        xage2 = agebins[j][1]
                        PDPpart = yPDP[1][0][int(xage1/xdif):int((xage2+xdif)/xdif)]
                        PDP_agePart = np.arange(xage1, xage2+xdif, xdif)
                        PDP_agePart = np.arange(xage1, xage1+len(PDPpart), xdif)
                        axsPDP.fill_betweenx(PDP_agePart, 0, PDPpart, alpha = 1, color=agebinsc[j])
        if plotKDE:
            if bw == 'optimizedFixed':
                yKDE = KDEcalcAgesLocalAdapt(ages=[yF], xdif=xdif)
            if bw == 'optimizedVariable':
                yKDE = KDEcalcAgesGlobalAdapt(ages=[yF], xdif=xdif)
            if type(bw) != str:
                yKDE = KDEcalcAges_2(ages=[yF], xdif=xdif, bw=bw)            
            axsKDE = axs[0,1].twiny()
            axsKDE.plot(yKDE[1][0],yKDE[0], color='black', lw=1)
            axsKDE.set_ylim(y1, y2)
            axsKDE.set_xlim(np.max(yKDE[1][0])*1.1,0)
            axsKDE.get_xaxis().set_ticks([])
            if colorKDE:
                axsKDE.fill_betweenx(yKDE[0], 0, yKDE[1][0], alpha = 1, color=colorMe(i))
            if colorKDEbyAge:
                if len(np.shape(agebins)) == 1:
                    nage = len(agebins)-1
                    for j in range(nage):
                        xage1 = agebins[j]
                        xage2 = agebins[j+1] 
                        KDEpart = yKDE[1][0][int(xage1/xdif):int((xage2+xdif)/xdif)]
                        KDE_agePart = np.arange(xage1, xage1+len(KDEpart), xdif)
                        axsKDE.fill_betweenx(KDE_agePart, 0, KDEpart,  alpha = 1, color=agebinsc[j])
                if len(np.shape(agebins)) == 2:
                    for j in range(len(agebins)):
                        xage1 = agebins[j][0]
                        xage2 = agebins[j][1]
                        KDEpart = yKDE[1][0][int(xage1/xdif):int((xage2+xdif)/xdif)]
                        KDE_agePart = np.arange(xage1, xage2+xdif, xdif)
                        KDE_agePart = np.arange(xage1, xage1+len(KDEpart), xdif)
                        axsKDE.fill_betweenx(KDE_agePart, 0, KDEpart, alpha = 1, color=agebinsc[j])                    

        # Label the upper left hand corner of the plot with the sample ID
        if type(sampleList[0])==tuple:
            axs[0,0].text(x=x1+(x2-x1)*0.02, y=y2+(y2-y1)*-0.05, s=sampleList[i][1])
        else:
            axs[0,0].text(x=x1+(x2-x1)*0.02, y=y2+(y2-y1)*-0.05, s=sampleList[i])
        fig.subplots_adjust(wspace=0)
        fig.subplots_adjust(hspace=0)
        if savePlot:
            pathlib.Path('Output').mkdir(parents=True, exist_ok=True) # Recursively creates the directory and does not raise an exception if the directory already exists 
            if type(sampleList[0])==tuple:
                fig.savefig(('Output/DoubleDating_%s' %sampleList[i][1])+('.pdf'))
            else:
                fig.savefig(('Output/DoubleDating_%s' %sampleList[i])+('.pdf'))    
    return fig
            
def exportDist(ages, errors, labels, exportType, cumulative, x1, x2, xdif, bw, fileName, normalize):
    """
    Creates a CSV file with raw age distribution data. Distribution types supported are cumulative density functions (CDF), probability density plots (PDPs), and kernal density estimations (KDEs).

    Parameters
    ----------
    ages : array of ages for each sample or sample group. Output from sampleToData()
    errors : array of errors for each sample or sample group. Output from sampleToData()
    labels : array of labels for each sample or sample group. Output from sampleToData()
    exportType : specify the type of distribution to export. Options: 'CDF', 'PDP', or 'KDE'
    x1 : lower limit (Myr) of x-axis
    x2 : upper limit (Myr) of x-axis
    bw : KDE bandwidth. Options are 'optimizedFixed', 'optimizedVariable', or a number (bandwidth in Myr)
    xdif : interval (Myr) over which distributions are calculated
    normalize : set to True to require distributions to sum to 1
    fileName : specify name of output CSV file
    
    Notes
    -----
    """    
    import csv
    
    if exportType == 'CDF':
        distAge, dist = CDFcalcAges(ages, x1, x2, xdif)
    if exportType == 'KDE':
        if bw == 'optimizedFixed':
            distAge, dist = KDEcalcAgesLocalAdapt(ages, x1, x2, xdif, cumulative)
        if bw == 'optimizedVariable':
            distAge, dist = KDEcalcAgesGlobalAdapt(ages, x1, x2, xdif, cumulative)
        if type(bw) != str:
            distAge, dist = KDEcalcAges_2(ages, x1, x2, xdif, bw, cumulative)       
    if exportType == 'PDP':
        distAge, dist = PDPcalcAges(ages, errors, x1, x2, xdif, cumulative)
    
    for i in range(len(dist)):
        for j in range(len(dist[i])):
            dist[i][j] = round(dist[i][j], 5)
        if normalize:
            dist[i] = dist[i]/sum(dist[i])      
    
    pathlib.Path('Output').mkdir(parents=True, exist_ok=True) # Recursively creates the directory and does not raise an exception if the directory already exists 
    with open(str(pathlib.Path('Output/')) + fileName, 'w', newline='') as f: #Select the CSV file to save data to
        writer = csv.writer(f)

        dataRow = ['Age'] # create an empty array    
        for j in range(len(distAge)):
            dataRow = np.append(dataRow, distAge[j])
        writer.writerow(dataRow)

        for i in range(len(ages)):
            dataRow = [labels[i]] # clear the dataRow        
            for k in range (len(distAge)):
                dataRow = np.append(dataRow, dist[i][k])
            writer.writerow(dataRow)
    
def agesErrorsCSV(ages, errors, sampleList, fileName):
    """
    Creates a CSV file with sample or sample group ages and 1 sigma errors reported in the two-column format that is used by other plotting and analysis software (e.g., DZstats (Saylor and Sundel, 2016) and Arizona LaserChron Center in-house excel macros).

    Parameters
    ----------
    ages : array of ages for each sample or sample group. Output from sampleToData()
    errors : array of errors for each sample or sample group. Output from sampleToData()
    sampleList : array of sample IDs.
    Must be in form for individual samples: ['Sample1', 'Sample2', . . . , etc.].
    Must be in the form for groups of samples: [(['Sample1','Sample2', . . . , etc.], 'Group 1 name'),
                                                (['Sample1','Sample2', . . . , etc.], 'Group 2 name')]
    fileName : specify name of output CSV file
    
    Notes
    -----
    """        
    import csv
    import itertools

    agesErrors = [sorted(list(zip(ages[i], errors[i])), key=lambda d: d[0]) for i in range(len(ages))]
    maxNumGrains = max([len(l) for l in agesErrors])

    pathlib.Path('Output').mkdir(parents=True, exist_ok=True) # Recursively creates the directory and does not raise an exception if the directory already exists 
    with open(str(pathlib.Path('Output/')) + fileName, 'w', newline='') as f: #Select the CSV file to save data to
        writer = csv.writer(f)

        # Write sample name column headings. The extra blank is to take up two columns per sample heading.
        sampleHeaders = list(itertools.chain.from_iterable([[sample, ''] for sample in sampleList]))
        writer.writerow(sampleHeaders)
        
        # Write Age and Error column headings for each sample.
        ageErrorHeaders = ['Age', 'Error'] * len(sampleList)
        writer.writerow(ageErrorHeaders)

        # Write data rows.
        for i in range(maxNumGrains):
            row = []
            for j in range(len(sampleList)):
                if len(agesErrors[j]) > i:
                    row.extend(agesErrors[j][i])
                else:
                    # Skip sample if doesn't have at least i grains. Write blanks as filler.
                    row.extend(['', ''])
            writer.writerow(row)    
    
def calcComparisonCSV(ages, errors, numGrains, labels, sampleList, calculateSimilarity, calculateLikeness, calculateKS, calculateKuiper, 
                  calculateR2, fileName, distType, bw):
    """
    Creates matricies of sample comparisons using a number of different metrics (see Saylor and Sundell, 2016). Similiarity, likness, Kolgomorov-Smirnov statistic (Dmax and p-value), Kuiper statistic (Vmax and p-value), and cross-correlation of relative probability density functions. Similiarty, likeness, and cross-correlation values are computed based on either the probability density plot (PDP) or kernal density estimation (KDE).

    Parameters
    ----------
    ages : array of ages for each sample or sample group. Output from sampleToData()
    errors : array of errors for each sample or sample group. Output from sampleToData()
    numGrains : array of numbers of analyses within each sample or sample group. Output from sampleToData()
    labels : array of labels for each sample or sample group. Output from sampleToData()
    sampleList : array of sample IDs.
    Must be in form for individual samples: ['Sample1', 'Sample2', . . . , etc.].
    Must be in the form for groups of samples: [(['Sample1','Sample2', . . . , etc.], 'Group 1 name'),
                                                (['Sample1','Sample2', . . . , etc.], 'Group 2 name')]
    
    calculateSimilarity : set to True to calculate the similarity index
    calculateLikeness : set to True to caclulate the likeness index
    calculateKS : set to True to calculate Dmax and the Kolgomorov-Smirnov statistic
    calculateKuiper : set to True to calculate Vmax and the Kuiper statistic
    calcualteR2 : set to True to calculate the cross-correlation coefficient
    fileName : specify name of output CSV file
    distType : specify the type of relative distribution to compare for similarity, likeness, and cross-correlation calculations. Options: 'PDP' or 'KDE'
    bw : KDE bandwidth. Options are 'optimizedFixed', 'optimizedVariable', or a number (bandwidth in Myr)
    
    Notes
    -----
    """     
    import csv
    import numpy as np
    from scipy import stats

    # Determine what type of relative distribution to use
    if (calculateSimilarity or calculateLikeness or calculateR2):
        if distType == 'PDP':
            dist = PDPcalcAges(ages, errors)[1]
        if distType == 'KDE':
            if bw == 'optimizedFixed':
                dist = KDEcalcAgesLocalAdapt(ages)[1]
            if bw == 'optimizedVariable':
                dist = KDEcalcAgesGlobalAdapt(ages)[1]
            if type(bw) != str:
                dist = KDEcalcAges_2(ages, bw=bw)[1]
    
    pathlib.Path('Output').mkdir(parents=True, exist_ok=True) # Recursively creates the directory and does not raise an exception if the directory already exists 
    with open(str(pathlib.Path('Output/')) + fileName, 'w', newline='') as f: #Select the CSV file to save data to , 'wb'
        writer = csv.writer(f)
        labelRow = ['Label','n']
        for i in range(len(sampleList)):
            labelRow = np.append(labelRow, labels[i])

        # Calculate the similiary matrix
        if calculateSimilarity:
            writer.writerow("")
            simMatrix = np.empty(shape=(len(sampleList),len(sampleList)))
            writer.writerow((['Similarity']))
            writer.writerow(labelRow)
            for i in range(len(sampleList)):
                dataRow = [labels[i],numGrains[i]]
                for j in range(len(sampleList)):
                    simMatrix[i][j] = calcSimilarity(dist[i],dist[j])
                    dataRow = np.append(dataRow,simMatrix[i][j])
                writer.writerow(dataRow)

        # Calculate the likeness matrix
        if calculateLikeness:
            writer.writerow("")
            likeMatrix = np.empty(shape=(len(sampleList),len(sampleList)))
            writer.writerow((['Likeness']))
            writer.writerow(labelRow)        
            for i in range(len(sampleList)):
                dataRow = [labels[i],numGrains[i]]
                for j in range(len(sampleList)):
                    likeMatrix[i][j] = calcLikeness(dist[i],dist[j])
                    dataRow = np.append(dataRow,likeMatrix[i][j])
                writer.writerow(dataRow)

        # Calculate the K-S Dmax matrix
        if calculateKS:
            writer.writerow("")
            KSdMatrix = np.empty(shape=(len(sampleList),len(sampleList)))

            writer.writerow((['Kolgomorov-Smirnov Dmax']))
            writer.writerow(labelRow)
            # Calculate the K-S Dmax and p-value matrices
            for i in range(len(sampleList)):
                dataRow = [labels[i],numGrains[i]]
                for j in range(len(sampleList)):
                    KSdMatrix[i][j] = stats.ks_2samp(ages[i],ages[j])[0]
                    dataRow = np.append(dataRow,KSdMatrix[i][j])
                writer.writerow(dataRow)
            writer.writerow("")
            KSpMatrix = np.empty(shape=(len(sampleList),len(sampleList)))        
            writer.writerow((['Kolgomorov-Smirnov p-value']))  
            writer.writerow(labelRow)
            for i in range(len(sampleList)):
                dataRow = [labels[i],numGrains[i]]
                for j in range(len(sampleList)):
                    KSpMatrix[i][j] = stats.ks_2samp(ages[i],ages[j])[1]
                    dataRow = np.append(dataRow,KSpMatrix[i][j])
                writer.writerow(dataRow)

        # Calculate the Kuiper Vmax matrix
        if calculateKuiper:
            writer.writerow("")        
            kupierVMatrix = np.empty(shape=(len(sampleList),len(sampleList)))
            writer.writerow((['Kuiper Vmax']))
            writer.writerow(labelRow)        
            CDF = CDFcalcAges(ages)
            for i in range(len(sampleList)):
                dataRow = [labels[i],numGrains[i]]
                for j in range(len(sampleList)):
                    kupierVMatrix[i][j] = calcKuiper(CDF[1][i],CDF[1][j],len(ages[i]),len(ages[j]))[0]    
                    dataRow = np.append(dataRow,kupierVMatrix[i][j])
                writer.writerow(dataRow)

        # Calculate the cross-correlation of a relative distribution function
        if calculateR2:
            writer.writerow("")
            r2Matrix = np.empty(shape=(len(sampleList),len(sampleList)))
            if distType == 'PDP':
                writer.writerow((['Cross-correlation of probability density plot']))
            if distType == 'KDE':
                writer.writerow((['Cross-correlation of kernal density estimation']))
            writer.writerow(labelRow)        
            for i in range(len(sampleList)):
                dataRow = [labels[i],numGrains[i]]
                for j in range(len(sampleList)):
                    r2Matrix[i][j] = calcR2(dist[i],dist[j])
                    dataRow = np.append(dataRow, r2Matrix[i][j])
                writer.writerow(dataRow)            

def weightedMeanCSV(ages, errors, numGrains, labels, fileName='weightedMean.csv'):
    """
    Creates a CSV file with weighted mean calculations.

    Parameters
    ----------
    ages : array of ages for each sample or sample group. Output from sampleToData()
    errors : array of errors for each sample or sample group. Output from sampleToData()
    numGrains : array of numbers of analyses within each sample or sample group. Output from sampleToData()    
    labels : array of labels for each sample or sample group. Output from sampleToData()    
    
    Notes
    -----
    """    
    import csv
    
    pathlib.Path('Output').mkdir(parents=True, exist_ok=True) # Recursively creates the directory and does not raise an exception if the directory already exists 
    with open(str(pathlib.Path('Output/')) + fileName, 'w', newline='') as f: #Select the CSV file to save data to
        writer = csv.writer(f)
        labelRow = ['Sample','Ngrains','Mean','2s internal error','MSWD']
        writer.writerow(labelRow)
        for i in range(len(ages)):
            dataRow = [] # create an empty array
            dataRow = np.append(dataRow, labels[i])
            dataRow = np.append(dataRow, numGrains[i])
            dataRow = np.append(dataRow, weightedMean(ages[i], errors[i])[0])
            dataRow = np.append(dataRow, weightedMean(ages[i], errors[i])[1])
            dataRow = np.append(dataRow, weightedMean(ages[i], errors[i])[2])
            writer.writerow(dataRow)
                
###############################################################
# Helper functions 
###############################################################            

def calcDFW(CDF, epsilon):
    """
    Computes the Kvoretsky-Kiefer-Wolfowitz inequality for a cumulative distribution function.
    
    Parameters
    ----------
    CDF : array of a cumulative distribution function for a single sample or sample group
    episilon : the half-width of the simultaneous confidence interval
    
    Returns
    -------
    DFWmin : array that defines the minimum bound on the Kvoretsky-Kiefer-Wolfowitz inequality
    DFWmax : that defines the maximum bound on the Kvoretsky-Kiefer-Wolfowitz inequality
    
    Notes
    -----
    Based on Anderson et al. (2018): Basin Research, v. 30, p. 132-147, doi: 10.1111/bre.12245
    """    
    DFWmin = CDF-epsilon
    DFWmax = CDF+epsilon
    DFWmin[DFWmin < 0] = 0
    DFWmax[DFWmax > 1] = 1
    
    return DFWmin, DFWmax
                

def CDFcalcAges(ages, x1=0, x2=4501, xdif=1):
    """
    Computes the CDF for an array of samples.
    
    Parameters
    ----------
    ages : array of ages, len(ages)=number of samples or sample groups
    x1 : (optional) beginning of range to compute CDF (default = 4500 Ma)
    x2 : (optional) end of range to compute CDF (default = 4500 Ma)
    xdif : (optional) bin size to compute CDF (default = 1 Ma)
    
    Returns
    -------
    CDF_age : array of ages that CDF is computed over    
    CDF : array of CDF functions
    
    Notes
    -----
    """
    CDF_age = np.arange(x1, x2+xdif, xdif)
    CDF = np.empty(shape=(len(ages),len(CDF_age)))
    for i in range(len(ages)):
        for j in range(len(CDF_age)):
            CDF[i][j] = np.sum(ages[i] <= CDF_age[j])
        CDF[i] = CDF[i]*1.0/len(ages[i])
    return CDF_age, CDF

def PDPcalcAges(ages, errors, x1=0, x2=4500, xdif=1, cumulative=False):    
    """
    Computes the PDP for an array of ages.
    
    Parameters
    ----------
    ages : array of ages, len(ages)=number of samples or sample groups
    errors : array of 1s errors
    x1 : (optional) beginning of range to compute PDP (default = 0 Ma)
    x2 : (optional) end of range to compute PDP (default = 4500 Ma)
    xdif : (optional) bin size to compute PDP (default = 1 Ma)
    cumulative : (optional) If True, will compute a cumulative PDP (CPDP)
    
    Returns
    -------
    PDP_age : array of ages that PDP is computed over    
    PDP : array of PDP functions
    
    Notes
    -----
    """
    from scipy.stats import norm
    PDP_age = np.arange(0, 4500+xdif, xdif) # Ensures that the PDP is calculated over all of geologic time
    PDPe = np.empty(shape=(len(ages),len(PDP_age))) # Create an empty array of appropriate shape
    PDP = np.zeros_like(PDPe) # Copy the array, but with zeros
    for i in range(len(ages)):
        data = ages[i]
        data_err = errors[i]
        pdf_cum = 0 # Creates an empty variable        
        for j in range(len(data)):      
            age = data[j]
            error = data_err[j]
            pdf = norm.pdf(PDP_age, age, error)
            pdf_cum = pdf_cum + pdf 
        pdf_cum = np.around(pdf_cum/np.trapz(pdf_cum), decimals=10)
        if cumulative:
            pdf_cum = np.cumsum(pdf_cum)        
        PDP[i] = pdf_cum
    PDP_age = PDP_age[int(x1/xdif):int((x2+xdif)/xdif)] # Only select the values within the specified plotting age range
    PDPportionRange = np.arange(x1, x2+xdif, xdif)
    PDPportionEmpty = np.empty(shape=(len(ages),len(PDPportionRange)))
    PDPportion = np.zeros_like(PDPportionEmpty) # Copy the array, but with zeros
    for i in range(len(ages)):
        PDPportion[i] = PDP[i][int(x1/xdif):int((x2+xdif)/xdif)] # Only select the values within the specified plotting age range
    return PDP_age, PDPportion

def KDEcalcAges_2(ages, x1=0, x2=4500, xdif=1, bw=2.5, cumulative=False):
    """
    Computes the KDE for an array of ages.
    
    Parameters
    ----------
    ages : array of ages, len(ages)=number of samples or sample groups
    x1 : (optional) beginning of range to compute KDE (default = 0 Ma)
    x2 : (optional) end of range to compute KDE (default = 4500 Ma)
    xdif : (optional) bin size to compute KDE (default = 1 Ma)
    bw : (optional) bandwidth used in KDE calculation (default = 2.5 Ma)
    cumulative : (optional) If True, will compute a cumulative KDE (CKDE) (default = False)
    
    Returns
    -------
    KDE_age : array of ages that KDE is computed over    
    KDE : array of KDE functions
    
    Notes
    -----
    """    
    from statsmodels.nonparametric.kde import KDEUnivariate

    KDE_age = np.arange(x1, x2+xdif, xdif)
    KDE = np.zeros(shape=(len(ages),len(KDE_age)))
    for i in range(len(ages)):
        kde = KDEUnivariate(ages[i])
        kde.fit(bw=bw)
        kde = kde.evaluate(KDE_age)
        if cumulative:
            kde = np.cumsum(kde)*xdif # Seems like this must be multiplied by xdif for it to work
        KDE[i,:] = kde
    return KDE_age, KDE    

def KDEcalcAgesLocalAdapt(ages, x1=0, x2=4500, xdif=1, cumulative=False):
    """
    Computes the KDE for an array of ages using an optimized, fixed kernal bandwidth following Shimazaki and Shinomoto, 2010: Journal of Computational Neuroscience, v. 29, p. 171-182, doi:10.1007/s10827-009-0180-4.
    
    Parameters
    ----------
    ages : array of ages, len(ages)=number of samples or sample groups
    x1 : (optional) beginning of range to compute KDE (default = 0 Ma)
    x2 : (optional) end of range to compute KDE (default = 4500 Ma)
    xdif : (optional) bin size to compute KDE (default = 1 Ma)
    bw : (optional) bandwidth used in KDE calculation (default = 2.5 Ma)
    cumulative : (optional) If True, will compute a cumulative KDE (CKDE)
    
    Returns
    -------
    KDE_age : array of ages that KDE is computed over    
    KDE : array of KDE functions
    
    Notes
    -----
    """
    import detritalpy.adaptiveKDE as akde

    KDE_age = np.arange(0, 4500+xdif, xdif) # Ensures that the KDE is calculated over all of geologic time
    KDE = np.zeros(shape=(len(ages),len(KDE_age)))
    for i in range(len(ages)):
        kde = akde.sskernel(np.asarray(ages[i]), tin=KDE_age, nbs=int(1e3))[0]
        if cumulative:
            kde = np.cumsum(kde)*xdif # Seems like this must be multiplied by xdif for it to work
        KDE[i,:] = kde
    KDE_age = KDE_age[int(x1/xdif):int((x2+xdif)/xdif)] # Only select the values within the specified plotting age range
    KDEportionRange = np.arange(x1, x2+xdif, xdif)
    KDEportionEmpty = np.empty(shape=(len(ages),len(KDEportionRange)))
    KDEportion = np.zeros_like(KDEportionEmpty) # Copy the array, but with zeros
    for i in range(len(ages)):
        KDEportion[i] = KDE[i][int(x1/xdif):int((x2+xdif)/xdif)] # Only select the values within the specified plotting age range
    return KDE_age, KDEportion
       
def KDEcalcAgesGlobalAdapt(ages, x1=0, x2=4500, xdif=1, cumulative=False):
    """
    Computes the KDE for an array of ages using an optimized, variable kernal bandwidth following Shimazaki and Shinomoto, 2010: Journal of Computational Neuroscience, v. 29, p. 171-182, doi:10.1007/s10827-009-0180-4.
    
    Parameters
    ----------
    ages : array of ages, len(ages)=number of samples or sample groups
    x1 : (optional) beginning of range to compute KDE (default = 0 Ma)
    x2 : (optional) end of range to compute KDE (default = 4500 Ma)
    xdif : (optional) bin size to compute KDE (default = 1 Ma)
    bw : (optional) bandwidth used in KDE calculation (default = 2.5 Ma)
    cumulative : (optional) If True, will compute a cumulative KDE (CKDE)
    
    Returns
    -------
    KDE_age : array of ages that KDE is computed over    
    KDE : array of KDE functions
    
    Notes
    -----
    Requires installation of the adaptivekde library.
    """
    import detritalpy.adaptiveKDE as akde

    KDE_age = np.arange(x1, x2+xdif, xdif)
    KDE = np.zeros(shape=(len(ages),len(KDE_age)))
    for i in range(len(ages)):
        kde = akde.ssvkernel(np.asarray(ages[i]), tin=KDE_age, nbs=int(1e2))[0]
        if cumulative:
            kde = np.cumsum(kde)*xdif # Seems like this must be multiplied by xdif for it to work
        KDE[i,:] = kde
    return KDE_age, KDE
    
def colorMe(i):
    """
    Returns a color for a given integer input
    
    Parameters
    ----------
    i : integer
    
    Returns
    -------
    color : the name of a color
    
    Notes
    -----
    """    
    colors = ['darkgreen', 'firebrick', 'darkblue', 'goldenrod', 'darkviolet', 'olive', 'lightcoral', 'gray', 'saddlebrown', 'dodgerblue', 'darkkhaki', 'mediumpurple', 'black', 'darkred', 'darkolivegreen']
    return colors[i%len(colors)]

def weightedMean(ages,error1s,conf=0.95):
    """
    Calculates the weighted mean, its 2-sigma uncertainty, and MSWD

    Paramters
    ---------
    ages : a 1D array of ages
    errors : an array of 1-sigma errors
    conf : (optional) confidence level

    Returns
    -------
    Twm : weighted mean age
    sm : 2-sigma uncertainty
    MSWD : Mean Square of the Weighted Deviation

    """


    from scipy import stats
    
    w=np.array(error1s)**(-2)/np.sum(np.array(error1s)**(-2)) # weight
    Twm=np.sum(w*np.array(ages)) # weight mean of age
    S=np.sum((np.array(ages)-Twm)**2/np.array(error1s)**2) # S
    N=len(ages)
    MSWD=S/(N-1) # Mean Square of the Weighted Deviation
    # Standard deviation of the weighted mean (2 sigma)
    sm=stats.norm.ppf(conf+(1-conf)/2.)*np.sqrt(1./np.sum(np.array(error1s)**(-2)))
    
    return(Twm,sm,MSWD)

def calcDmax(CDF1, CDF2):
    """
    Calculates Dmax (maximum separation between two CDFs)
    
    Parameters
    ----------
    CDF1 : First CDF array
    CDF2 : Second CDF array
        
    Returns
    -------
    Dmax : Value (float) of Dmax
    
    Notes
    -----
    CDF1 and CDF2 must have the same length.
    """      
    Dmax = max(np.abs(CDF1-CDF2))
    return Dmax

def calcSimilarity(dist1, dist2):
    """
    Computes the similarity metric on 2 samples.
    
    Parameters
    ----------
    dist1, dist2 : Relative probability distributions (e.g., PDP, KDE)
        
    Returns
    -------
    sim : float
        similarity metric
    
    Notes
    -----
    dist1 and dist2 must have the same length.
    Based on Saylor and Sundell, 2016: Geosphere, v. 12, doi:10.1130/GES01237.1    
    """
    sim = np.sum((dist1*dist2)**0.5)
    return sim

def calcLikeness(dist1, dist2):
    """
    Computes the likeness metric on 2 samples.
    
    Parameters
    ----------
    dist1, dist2 : Relative probability distributions (e.g., PDP, KDE)
        
    Returns
    -------
    sim : float
        likeness metric
    
    Notes
    -----
    dist1 and dist2 must have the same length.
    Based on Saylor and Sundell, 2016: Geosphere, v. 12, doi:10.1130/GES01237.1
    """
    likeness = 1-((np.sum(abs(dist1-dist2))/2))
    return likeness

def calcVmax(CDF1, CDF2):
    """
    Computes Vmax
    
    Parameters
    ----------
    CDF1, CDF2 : Cumulative probability distributions
        
    Returns
    -------
    V : Vmax
    
    Notes
    -----
    Based on Saylor and Sundell, 2016: Geosphere, v. 12, doi:10.1130/GES01237.1    
    """
    deltaCDF1 = CDF2 - CDF1
    maxDeltaCDF1 = np.max(deltaCDF1)
    deltaCDF2 = CDF1 - CDF2
    maxDeltaCDF2 = np.max(deltaCDF2)
    V = maxDeltaCDF1 + maxDeltaCDF2
    
    return V

def calcKuiper(CDF1, CDF2, nCDF1, nCDF2):
    """
    Computes the Kuiper statistic and Vmax on 2 samples.
    
    Parameters
    ----------
    CDF1, CDF2 : Cumulative probability distributions
    nCDF1, nCDF2 : Number of analyses in each cumulative probability distribution
        
    Returns
    -------
    V : Vmax
    p : p-value
    
    Notes
    -----
    Based on Saylor and Sundell, 2016: Geosphere, v. 12, doi:10.1130/GES01237.1
    """    
    deltaCDF1 = CDF2 - CDF1
    maxDeltaCDF1 = np.max(deltaCDF1)
    deltaCDF2 = CDF1 - CDF2
    maxDeltaCDF2 = np.max(deltaCDF2)
    V = maxDeltaCDF1 + maxDeltaCDF2
    ne = (nCDF1*nCDF2)/(nCDF1+nCDF2)
    Lambda = (np.sqrt(ne) + 0.155 + 0.25/np.sqrt(ne)) * V  
    if Lambda < 0.4:
        p = 1
    j = np.arange(1.,101.)
    pare = 4*Lambda**2*(j**2)-1
    expo = np.exp(-2*Lambda**2*(j**2))
    argo = pare*expo
    p = 2*np.sum(argo)
                  
    return V, p
                  
def calcR2(dist1, dist2):
    """
    Computes the cross-correlation coefficient on 2 samples.
    
    Parameters
    ----------
    dist1, dist2 : Relative probability distributions (e.g., PDP, KDE)
        
    Returns
    -------
    r2 : float
        cross-correlation coefficient
    
    Notes
    -----
    dist1 and dist2 must have the same length.
    Based on Saylor and Sundell, 2016: Geosphere, v. 12, doi:10.1130/GES01237.1
    """    
    lenDist1 = len(dist1)
    lenDist2 = len(dist2)
    meanDist1 = np.mean(dist1)
    meanDist2 = np.mean(dist2)             
    xcov = np.empty(shape=(lenDist1,1))
    ycov = np.empty(shape=(lenDist2,1))
    for i in range(lenDist1):
        xcov[i] = dist1[i] - meanDist1
    for i in range(lenDist2):
        ycov[i] = dist2[i] - meanDist2
    xcov = np.transpose(xcov)
    ycov = np.transpose(ycov)            
    mult = xcov*ycov
    numerator = np.sum(mult)     
    xcov2 = xcov*xcov
    sumxcov2 = np.sum(xcov2)
    ycov2 = ycov*ycov
    sumycov2 = np.sum(ycov2)
    multi2 = sumxcov2*sumycov2
    denominator = np.sqrt(multi2)     
    r = numerator/denominator
    r2 = r*r
                  
    return r2

def calcComplR2(dist1, dist2):
    """
    Computes the complement of the cross-correlation coefficient.
    
    Parameters
    ----------
    dist1, dist2 : Relative probability distributions (e.g., PDP, KDE)
        
    Returns
    -------
    cR2 : complement of the cross-correlation coefficient
    
    Notes
    -----
    dist1 and dist2 must have the same length.
    Based on Saylor et al., 2017: Basin Research, doi:10.1111/bre.12264
    """
    r2 = calcR2(dist1, dist2)
    cR2 = 1-r2
    return cR2
    

def peakAge(DF_age, DF, ages, errors, thres=0.05, minDist=5, minPeakSize=3):
    """
    Identifies peak ages in a given relative probability distribution (e.g., PDP, KDE)
    
    Parameters
    ----------
    DF_age, DF : Arrays of ages and probability values of a relative probability distribution (e.g., PDP, KDE)
    thres : Threshold of what constitues a peak (from 0 to 1). Default = 0.05
    minDist : Minimum distance (Myr) between adjacent peaks. Default = 5
        
    Returns
    -------
    peakAges : array of peak ages (Myr) for each sample or sample group
    indexes : array of the indexes where peak ages occur
    
    Notes
    -----
    Requires the package peakutils to be installed
    pip install peakutils
    """     
    import peakutils
    peakAges = []
    indexes = []
    for i in range(len(DF)):
        indexes = indexes + [list(peakutils.indexes(DF[i], thres=thres, min_dist=minDist))]
        peakAges = peakAges + [list(DF_age[indexes[i]])]
    peakAgeGrain = peakAgesGrains(peakAges, ages, errors)
    
    # Remove peaks with fewer grains than the minimum peak size
    for i in reversed(range(len(peakAges))):
        for j in reversed(range(len(peakAgeGrain[i]))):
            if peakAgeGrain[i][j] < minPeakSize:
                del peakAges[i][j]
                del peakAgeGrain[i][j]
                del indexes[i][j]
    return peakAges, indexes, peakAgeGrain

def peakAgesGrains(peakAges, ages, errors, sigma=2):
    import copy
    peakAgeGrain = copy.deepcopy(peakAges)
    for i in range(len(peakAges)): # One loop per sample or sample group
        for j in range(len(peakAges[i])): # One loop per peak
            c = 0 # counter variable
            for k in range(len(ages[i])): # Loop through each analysis
                if (ages[i][k]-errors[i][k]*sigma <= peakAges[i][j] <= ages[i][k]+errors[i][k]*sigma):
                    c = c+1
            peakAgeGrain[i][j] = c
    return peakAgeGrain   

def exportPeakAge(labels, peakAges, peakAgesGrains, fileName='peakAges.csv'):
    import csv
    with open(fileName, 'w', newline='') as f: #Select the CSV file to save data to
        writer = csv.writer(f)
        writer.writerow((['Peak Age']))
        writer.writerow((''))
        for i in range(len(peakAges)):
            writer.writerow(([labels[i]]))
            writer.writerow((peakAges[i][j]) for j in range(len(peakAges[i])))
            writer.writerow((peakAgesGrains[i][j]) for j in range(len(peakAgesGrains[i])))

def calcMDA(ages, errors):

    for i in range(len(ages)):
        error1s = errors[i]
        error2s = errors[i]*2
        data_err1s_ageSort = list(zip(ages[i], error1s))            
        data_err1s = list(zip(ages[i], error1s))
        data_err2s = list(zip(ages[i], error2s))
        data_err1s_ageSort.sort(key=lambda d: d[0]) # Sort based on age
        data_err1s.sort(key=lambda d: d[0] + d[1]) # Sort based on age + 1s error
        data_err2s.sort(key=lambda d: d[0] + d[1]) # Sort based on age + 2s error
        YSG = data_err1s[0][0]
        YSG_err1s = data_err1s[0][1]

        def find_youngest_cluster_1s(min_cluster_size):
            i_min = 0
            i_max = 0
            for i in range(1, len(data_err1s)):
                top = data_err1s[i_min][0] + data_err1s[i_min][1]
                bottom = data_err1s[i][0] - data_err1s[i][1]
                if (top >= bottom):
                    i_max = i
                elif i_max - i_min + 1 >= min_cluster_size:
                    break
                else:
                    i_min = i
                    i_max = i
            return data_err1s[i_min: i_max + 1] if i_min < i_max else [], i_max

        def find_youngest_cluster_2s(min_cluster_size):
            i_min = 0
            i_max = 0
            for i in range(1, len(data_err1s)):
                top = data_err2s[i_min][0] + data_err2s[i_min][1]
                bottom = data_err2s[i][0] - data_err2s[i][1]
                if (top >= bottom):
                    i_max = i
                elif i_max - i_min + 1 >= min_cluster_size:
                    break
                else:
                    i_min = i
                    i_max = i
            return data_err2s[i_min: i_max + 1] if i_min < i_max else [], i_max

        YC1S, YC1S_imax = find_youngest_cluster_1s(2)
        YC1S_WM = weightedMean(np.array([d[0] for d in YC1S]), np.array([d[1] for d in YC1S]))
        YC2S, YC2S_imax = find_youngest_cluster_2s(3)
        YC2S_WM = weightedMean(np.array([d[0] for d in YC2S]), np.array([d[1] for d in YC2S])/2)

    return YSG, YSG_err1s, YC1S_WM[0], YC1S_WM[1]/2, YC1S_WM[2], len(YC1S), YC2S_WM[0], YC2S_WM[1]/2, YC2S_WM[2], len(YC2S)
