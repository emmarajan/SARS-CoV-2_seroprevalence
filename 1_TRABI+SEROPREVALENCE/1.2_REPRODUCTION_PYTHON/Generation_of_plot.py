#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 12:15:04 2020

@author: rjacquat

Script defining function for :
    generating the df fiving new files.
    loading existing df
    generate QDA
    comparing QDA model
    ploting prevalence
    ploting histogram (age pyramid)

__Example_plot__()
    generate with the default parameter
    load andcreat all plot possible
"""
import pandas as pd
import os
from glob import glob
import numpy as np
from math import lgamma
from load_data_extract_df import load_all_data, filter_USZ
from verify_age import sort_blood
from gaussian_estimation import filter_cohort, getPosteriorFromMVNFull
from scipy.stats import multivariate_normal as mvn
from library_of_plot_df import plot_violin
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.utils import resample
from gaussian_estimation import Probability_field_FMW_bootstraps
from gaussian_estimation import Probability_field_LDA_bootstraps


def dependancy():
    main_folder = path_to_main_folder()
    name_file_adhoc = "*ptn_adhoc_20200914_annot*.xlsx"
    name_file_adhoc_pos = "*ptn_adhoc_20200926_seropositi*.xlsx"
    name_QC = "*Overview_data_transfer*.csv"
    name_files_BDS = ("NEED to be update", "", "", "")
    list_paths = (main_folder,
                  name_file_adhoc,
                  name_file_adhoc_pos,
                  name_QC,
                  name_files_BDS)
    return list_paths


def setup_existing_df(df_path,
                      df_BDS,
                      df_USZ):
    """
    Load the different dataframes
    retuns
    ------
    df, df_BDS, df_USZ
    """
    df = pd.read_csv(df_path, sep="\t")
    df_BDS = pd.read_csv(df_BDS, sep="\t")
    df_USZ = pd.read_csv(df_USZ, sep="\t")
    df = df.set_index("Patient_or_Control_ID")
    df_BDS = df_BDS.set_index("Patient_or_Control_ID")
    df_USZ = df_USZ.set_index("Patient_or_Control_ID")
    return df, df_BDS, df_USZ


def generate_dataframe(mainfolder_path=False,
                       adhoc_path="*ptn_adhoc_20200914_annot*.xlsx",
                       QC_path="*Overview_data_transfer*.csv",
                       blood_path='../USZ_SARS/Meta_files/BDS/RunScript/',
                       saving_path=None):
    """
    Load all the EC50 values from all plate.
    Filter from QC and annotate with adhoc_path
    save 3 differents files in .tsv files and return their links

    Parameters
    ----------
    mainfolder_path : string
        folder pointing on the arborescence of all subfolder containing data
    adhoc_path : string
        link to adhoc file containing all additional info data
    QC_path : string
        link to QC file alowing to QC
    saving_path : string
        link where to save the different files
    Returns
    -------
    3 link to df, df_BDS, df_USZ
    df_path, df_BDS_path, df_USZ_path
    """

    dic_cohort_sub = {
            "hospital patient < 2019": ("hospital patient", -1),
            "hospital patient": ("hospital patient", 0),
            "blood donation negative": ("blood donation", -1),
            "serum donation positive": ("blood donation", 1),
            "blood donation": ("blood donation", 0),
            "Infectious_Disease": ("Infectious_Disease", 0),
            "Oxford known negative": ("oxford", -1),
            "Oxford acute positive": ("oxford", 1),
            "Oxford convalescent positive": ("oxford", 1),
            "Oxford low positive serial dilution": ("oxford", 1),
            "Oxford other resp acute": ("oxford", 0),
            "Oxford other coronavirus conv": ("oxford", 0),
            "Oxford high positive serial dilution": ("oxford", 1),
            "Oxford other coronavirus acute": ("oxford", 0),
            "Oxford medium positive serial dilution": ("oxford", 1),
            "Oxford_neg": ("oxford", -1),
            "KIM_pos": ("hospital patient", 1),
            "KIM_neg": ("hospital patient", -1),
            "Only number": ("dermatology", 0)
    }

    """
    Load data
    """
    df, df_pos = load_all_data(mainfolder_path, QC_path)
    df.loc[:, "hasCovid"] = df.loc[:, "cohort"]
    df.loc[:, "subcohort"] = df.loc[:, "cohort"]

    """
    creat the subcohort and define if known covid positive 1, no info 0 and
    known covid negative -1
    """
    for name in df.loc[:, "subcohort"].drop_duplicates().values:
        mask = df.loc[:, "subcohort"] == name
        df.loc[mask, "hasCovid"] = dic_cohort_sub[name][1]
        df.loc[mask, "cohort"] = dic_cohort_sub[name][0]

    """
    Remove any plate which does not contain three EC50 result.
    """
    df = df[~np.isnan(df.RBD_pEC50)]
    df = df[~np.isnan(df.NC_pEC50)]
    df = df[~np.isnan(df.Spike_pEC50)]

    """
    Generate the df_USZ using the files path
    """
    ad_hoc_file = adhoc_path
    df_USZ = filter_USZ(df, ad_hoc_file)

    df_blood = sort_blood(df, blood_path)
    df_blood = df_blood.loc[df_blood[
                                     ['NC_pEC50', 'Spike_pEC50', 'RBD_pEC50']
                                    ].dropna().index]
    """
    Save data into savingpath+name of df
    if savingpath is None, then saving in current folder
    """
    if saving_path is None:
        saving_path = os.path.dirname(os.path.abspath(''))

    df_USZ_path = os.path.join(saving_path, "df_USZ.tsv")
    df_USZ.to_csv(df_USZ_path, sep="\t")
    df_BDS_path = os.path.join(saving_path, "df_BDS.tsv")
    df_blood.to_csv(df_BDS_path, sep="\t")
    df_path = os.path.join(saving_path, "df.tsv")
    df.to_csv(df_path, sep="\t")

    return df_path, df_BDS_path, df_USZ_path


def dmvnorm(x, mu, Sigma, df, log):
    '''
    Multivariate t-student density. Returns the density
    of the function at points specified by x.

    input:
        x = parameter (n x d numpy array)
        mu = mean (d dimensional numpy array)
        Sigma = scale matrix (d x d numpy array)
        df = degrees of freedom
        log = log scale or not

    '''
    p_3 = 3  # Dimensionality
    dec = np.linalg.cholesky(Sigma)
    R_x_m = np.linalg.solve(dec,
                            np.transpose(np.asarray(x)-np.asarray(mu)))
    rss = np.power(R_x_m, 2).sum(axis=0)
    logretval = (lgamma(1.0 * (p_3 + df) / 2) -
                 (lgamma(1.0 * df / 2) + np.sum(np.log(dec.diagonal()))
                 + p_3/2 * np.log(np.pi * df)) -
                 0.5 * (df + p_3) * np.log1p((rss/df)))
    if log is False:
        return(np.exp(logretval))
    else:
        return logretval


def Calculate_LDA(df,
                  model,
                  list_pEC50=["NC_pEC50", "Spike_pEC50", "RBD_pEC50"],
                  seed=111):
    """
    Parameters
    ----------
    df : TYPE
        DESCRIPTION.
    model : TYPE
        DESCRIPTION.
    list_pEC50 : TYPE, optional
        DESCRIPTION. The default is ["NC_pEC50", "Spike_pEC50", "RBD_pEC50"].
    seed : TYPE, optional
        DESCRIPTION. The default is 111.

    Returns
    -------
    None.
    """
    np.random.seed(seed)
    if model == "USZ":
        USEBLOOD = False
    elif model == "BDS":
        USEBLOOD = True
    else:
        RuntimeError(f"Please fill model either BDS or USZ not : {model}")

    castTblUsed, castTblUnused, knownPos, knownNeg = filter_cohort(df,
                                                                   USEBLOOD)
    posIndexVec = np.random.choice(np.repeat(np.arange(0, 10),
                                   np.ceil(len(knownPos)/10))[0:len(knownPos)],
                                   size=len(knownPos),
                                   replace=True,
                                   p=None)
    negIndexVec = np.random.choice(np.repeat(np.arange(0, 10),
                                   np.ceil(len(knownNeg)/10))[0:len(knownNeg)],
                                   size=len(knownNeg),
                                   replace=True,
                                   p=None)
    # mymat = castTblUsed.loc[:, list_pEC50]
    allVals = castTblUsed.loc[:, "Spike_pEC50"]
    for myi in np.arange(0, 10):
        removedPos = knownPos.iloc[posIndexVec == myi]
        removedNeg = knownNeg.iloc[negIndexVec == myi]
        torm = pd.concat([removedPos, removedNeg])
        mymat2 = castTblUsed.loc[~castTblUsed.index.isin(torm.index),
                                 list_pEC50]
        mat_pos_remove = knownPos.loc[~knownPos.index.isin(removedPos.index),
                                      list_pEC50]
        mat_neg_remove = knownNeg.loc[~knownNeg.index.isin(removedNeg.index),
                                      list_pEC50]

        fit = estimateLDAFrac(mymat2,
                              mat_pos_remove,
                              mat_neg_remove,
                              nround=50)
        myposPosterior = getPosteriorFromMVNFull(removedPos.loc[:, list_pEC50],
                                                 fit)
        mynegPosterior = getPosteriorFromMVNFull(removedNeg.loc[:, list_pEC50],
                                                 fit)
        allVals[allVals.index.isin(knownPos.iloc[posIndexVec ==
                                                 myi].index)] = myposPosterior
        allVals[allVals.index.isin(knownNeg.iloc[negIndexVec ==
                                                 myi].index)] = mynegPosterior

    fullfit = estimateLDAFrac(castTblUsed.loc[:, list_pEC50],
                              knownPos.loc[:, list_pEC50],
                              knownNeg.loc[:, list_pEC50],
                              nround=50)
    unknown = ~(allVals.index.isin(knownPos.index) |
                allVals.index.isin(knownNeg.index))
    allVals[unknown] = fullfit[1][unknown]

    allVals.to_frame().rename({"Spike_pEC50": "PosteriorProb"})

    df["LDA_prob"] = allVals
    return df, fullfit[2]


def estimateLDAFrac(valueMat, knownPos, knownNeg, nround=50):
    """
    Generate posterior value using S,RBD,NC pEC50 values
    using LDA model

    Parameters
    ----------
    valueMat : 3D Array
        S, RBD, NC pEC50 values to evaluate
    knownPos : 3D Array
        known positive values to train data
    knownNeg : 3D Array
        known negative values to train data
    nround : int
        number of round to train the data. The default is 50.

    Returns
    -------
    outl : myp, myt, paraml
        myp is the final result where 1
        myt is the final result where
        paraml parameter for QDA
    """
    use1dim = True
    use1dimCor = True
    mean1 = np.nanmean(knownPos, axis=0)
    mean2 = np.nanmean(knownNeg, axis=0)
    cov2 = np.cov(np.transpose(knownNeg))
    mean_cov1 = mvn(mean1, cov2)
    mean_cov2 = mvn(mean2, cov2)
    mask_pos = valueMat.index.isin(knownPos.index)
    mask_neg = valueMat.index.isin(knownNeg.index)
    myp = 0.02
    linearComb = np.linalg.solve(cov2,mean1-mean2) #  solve(cov2)%*%(mean1-mean2)
    tstat = valueMat.dot(linearComb) #  valueMat%*%linearComb
    linVar = np.transpose(mean1-mean2).dot(linearComb) #  (t(mean1-mean2)%*%solve(cov2)%*%(mean1-mean2))[1]
    linMean1 = np.transpose(mean1-mean2).dot(np.linalg.solve(cov2,mean1)) #  (t(mean1-mean2)%*%solve(cov2)%*%(mean1))[1]
    linMean2 = np.transpose(mean1-mean2).dot(np.linalg.solve(cov2,mean2)) #  (t(mean1-mean2)%*%solve(cov2)%*%(mean2))[1]
    sd1 = np.sqrt(np.var(tstat[knownPos.index]))
    sd2 = np.sqrt(np.var(tstat[knownNeg.index]))
    for i in np.arange(0, nround):
        if use1dim:
            density1_lin = norm.pdf(tstat, linMean1, np.sqrt(linVar)) #  dnorm(tstat,mean=linMean1,sqrt(linVar))
            density2_lin = norm.pdf(tstat, linMean2, np.sqrt(linVar))

            if use1dimCor:
                density1 = norm.pdf(tstat,
                                    np.mean(tstat[knownPos.index]),
                                    sd1)
                density2 = norm.pdf(tstat,
                                    np.mean(tstat[knownNeg.index]),
                                    sd2)

            density1 = density1_lin
            density2 = density2_lin
        else:

            density1 = mean_cov1.pdf(valueMat)
            density2 = mean_cov2.pdf(valueMat)
            density1 = dmvnorm(valueMat, mu=mean1, Sigma=cov2, df=1, log=False)
            density2 = dmvnorm(valueMat, mu=mean2, Sigma=cov2, df=1, log=False)

        myt = (myp * density1) / ((myp * density1) + ((1 - myp) * density2))
        myt[mask_pos] = 1
        myt[mask_neg] = 0
        myp = np.mean(myt[(mask_pos == False) & (mask_neg == False)])

    cov1 = np.cov(np.transpose(knownPos))
    paraml = (mean1, cov1, mean2, cov2, linearComb,
              linMean1, linMean2, sd1, sd2)
    outl = (myp, myt, paraml, tstat)
    return(outl)


def Calculate_QDA(df,
                  model,
                  list_pEC50=["NC_pEC50", "Spike_pEC50", "RBD_pEC50"],
                  seed=111):
    """
    Parameters
    ----------
    df : Pandas data frame
        data frame containing logEC50 of NC, S, RBD
    model : String
        Can be "USZ", "BDS"
        Knowing "USZ" incorporate the positive, "BDS" incorporate the negative
    list_pEC50: array of string
        array containing columns name of the model
    seed: int
        seed of the random number generated (default 111)
    Returns
    -------
    df : Pandas data frame
        return data frame with value QDA in QDA_prob
    param_fit : list of vector
        Paarameter for the fit (mean_pos, cov1, mean_neg, cov2)
        mean_pos = 3D vectors (NC, S, RBD)
        cov1 = matrix of the covariance parameter used
        mean neg = 3D vectors (NC, S, RBD)
        cov2 = matrix of the covariance parameter used
    """
    np.random.seed(seed)
    if model == "USZ":
        USEBLOOD = False
    elif model == "BDS":
        USEBLOOD = True
    else:
        RuntimeError(f"Please fill model either BDS or USZ not : {model}")

    castTblUsed, castTblUnused, knownPos, knownNeg = filter_cohort(df,
                                                                   USEBLOOD)

    posIndexVec = np.random.choice(np.repeat(np.arange(0, 10),
                                   np.ceil(len(knownPos) /
                                           10))[0:len(knownPos)],
                                   size=len(knownPos),
                                   replace=True,
                                   p=None)
    negIndexVec = np.random.choice(np.repeat(np.arange(0, 10),
                                   np.ceil(len(knownNeg) /
                                           10))[0:len(knownNeg)],
                                   size=len(knownNeg),
                                   replace=True,
                                   p=None)
    allVals = castTblUsed.loc[:, "Spike_pEC50"]
    for myi in np.arange(0, 10):
        removedPos = knownPos.iloc[posIndexVec == myi]
        removedNeg = knownNeg.iloc[negIndexVec == myi]
        torm = pd.concat([removedPos, removedNeg])
        mymat2 = castTblUsed.loc[~castTblUsed.index.isin(torm.index),
                                 list_pEC50]
        mat_pos_remove = knownPos.loc[~knownPos.index.isin(removedPos.index),
                                      list_pEC50]
        mat_neg_remove = knownNeg.loc[~knownNeg.index.isin(removedNeg.index),
                                      list_pEC50]
        fit = estimateQDAFrac(mymat2, mat_pos_remove, mat_neg_remove)
        myposPosterior = getPosteriorFromMVNFull(removedPos.loc[:, list_pEC50],
                                                 fit)
        mynegPosterior = getPosteriorFromMVNFull(removedNeg.loc[:, list_pEC50],
                                                 fit)
        allVals[allVals.index.isin(
            knownPos.iloc[posIndexVec == myi].index)] = myposPosterior
        allVals[allVals.index.isin(
            knownNeg.iloc[negIndexVec == myi].index)] = mynegPosterior

    fullfit = estimateQDAFrac(castTblUsed.loc[:, list_pEC50],
                              knownPos.loc[:, list_pEC50],
                              knownNeg.loc[:, list_pEC50],
                              nround=50)
    unknown = ~(allVals.index.isin(knownPos.index) |
                allVals.index.isin(knownNeg.index))
    allVals[unknown] = fullfit[1][unknown]
    allVals.to_frame().rename({"Spike_pEC50": "PosteriorProb"})

    df["QDA_prob"] = allVals

    return df, fullfit[2]


def estimateQDAFrac(valueMat, knownPos, knownNeg, nround=50):
    """
    Generate posterior value using S,RBD,NC pEC50 values

    Parameters
    ----------
    valueMat : 3D Array
        S, RBD, NC pEC50 values to evaluate
    knownPos : 3D Array
        known positive values to train data
    knownNeg : 3D Array
        known negative values to train data
    nround : int
        number of round to train the data. The default is 50.

    Returns
    -------
    outl : myp, myt, paraml
        myp is the final result where 1
        myt is the final result where
        paraml parameter for QDA

    """
    mean_pos = np.mean(knownPos, axis=0)
    mean_neg = np.mean(knownNeg, axis=0)
    cov1 = np.cov(np.transpose(knownPos))
    cov2 = np.cov(np.transpose(knownNeg))
    mask_pos = valueMat.index.isin(knownPos.index)
    mask_neg = valueMat.index.isin(knownNeg.index)
    myp = 0.005
    for i in np.arange(1, nround):
        mean_cov1 = mvn(mean_pos, cov1)
        mean_cov2 = mvn(mean_neg, cov2)
        density1 = mean_cov1.pdf(valueMat)
        density2 = mean_cov2.pdf(valueMat)
        myt = (myp * density1) / ((myp * density1) + ((1 - myp) * density2))
        myt[mask_pos] = 1
        myt[mask_neg] = 0
        myp = np.nanmean(myt[(~mask_pos) & (~mask_neg)])
        # print(myp)

    # print("done")
    paraml = (mean_pos, cov1, mean_neg, cov2)
    outl = (myp, myt, paraml)
    return outl


def QDAvalue_system(valueMat, prob_mean, parameter):
    """
    Calculate the QDA value from a 3D vector of pEC50 in NC, S, RBD, given a
    specific already trained model parameter (mean probability, cov matrix
                                              positif, mean_positive value,
                                              cov matrix negative)

    Parameters
    ----------
    valueMat : 3D array (Pandas df of 3 columns (NC, S, RBD))
        Value which want to be tested
    prob_mean : float value (between 0-1)
        mean probability of the population (month or entire set)
        which need to be evaluated
    parameter : list of (mean_pos, cov1, mean_neg, cov2)
        mean_pos list of 3 mean of the known positif NC, S, RBD
        cov1 matrix of the known positif 3D vector from pEC50 NC, S, RBD
        mean_pos list of 3 mean of the known negatif NC, S, RBD
        cov2 matrix of the known negatif 3D vector from pEC50 NC, S, RBD


    Returns
    -------
    myt : list of float (value between 0,1)
        list of float representing the same order as the index appear
        in ValueMat

    """
    myp = prob_mean
    (mean_pos, cov1, mean_neg, cov2) = parameter
    mean_cov1 = mvn(mean_pos, cov1)
    mean_cov2 = mvn(mean_neg, cov2)
    density1 = mean_cov1.pdf(valueMat)
    density2 = mean_cov2.pdf(valueMat)
    myt = (myp * density1) / ((myp * density1) + ((1 - myp) * density2))
    return myt


def plot_prevalence(df, columns):
    """
    Function which plot prevalence in function of the year and the month
    Becarefull does not save the figure in a file, and also
    is problematic for separate cohort/subcohort with same year and month.
    For example the positive serum donation donnor.
    Parameters
    ----------
    df : pandas dataframe
        df which you want to plot.
    columns : string
        the name of the column that will be plot

    Returns
    -------
    bool
        return True if stuff work.

    """
    for year in df.year.drop_duplicates():
        for month in df.month.drop_duplicates():
            mask_year = df.loc[:, "year"] == year
            mask_month = df.loc[:, "month"] == month
            plt.bar(month, np.nanmean(df.loc[mask_month & mask_year,
                                             columns]))
    return True


def plot_prevalence_withknownPos_USZ(df, columns):
    """
    Plot prevalence in function of the parameter "x axis violin plot columns"
    This plot can only be made on USZ after the ad_hoc_file is given,
    The columns 'delta_day_PCRpos' should be resent as well as hasCovid

    Parameters
    ----------
    df : Pandas dataframe
        data frame set you want to plot.
    columns : string
        Name of the columns you want to plot
    color : 3D or 4D array of float (between 0-1)
        color of the bar plot which will be displaied
    centered : string, optional
        Can be either "c", "r", "l". if the plot is offcentred on the right or
        on the leftThe default is "c".

    Returns
    -------
    bool
        if function ex3ecute correctly return True

    """

    list_column = ["year", "month"]

    temp_df = df.loc[df["x axis violin name"].drop_duplicates().index]
    temp_df = temp_df.sort_values(by=list_column)
    list_name = temp_df.sort_values(by=list_column).loc[:,
                                                        "x axis violin name"]
    list_mean = []
    list_mean_pcrpos = []
    list_mean_pos = []
    mask_0days = df.loc[:, 'delta_day_PCRpos'] < 0
    # mask_14days = df.loc[:, 'delta_day_PCRpos'] < -14
    mask_hasCovid = df.loc[:, "hasCovid"] == 1
    for name in list_name:
        mask_name = df.loc[:, "x axis violin name"] == name
        posterior_probability = df.loc[:,
                                       columns]
        posterior_probability[df.loc[:, columns].isna()] = 0
        list_mean.append(np.nanmean(df.loc[mask_name,
                                           columns]))
        total_len = len(df.loc[mask_name,
                               columns])
        value_temp = np.nansum(df.loc[mask_name & mask_0days,
                                      columns])
        list_mean_pcrpos.append(value_temp / total_len)

        value_temp = np.nansum(df.loc[((mask_name & mask_hasCovid) &
                                      mask_0days),
                                      columns])
        list_mean_pos.append(value_temp / total_len)

    dx = 0
    barWidth = 0.8

    x1 = np.arange(0, len(list_name))
    plt.bar(x1+dx,
            list_mean,
            width=barWidth,
            label=columns)
    plt.bar(x1+dx,
            list_mean_pcrpos,
            width=barWidth,
            label="PCR_detected")
    plt.bar(x1+dx,
            list_mean_pos,
            width=barWidth,
            label="hasCovid")
    plt.xticks(x1, list_name)
    plt.legend()

    return x1+dx, list_mean, list_mean_pcrpos, list_mean_pos


def plot_prevalence_name(df, model, columns, color, centered="c"):
    """
    Plot prevalence in function of the parameter "x axis violin plot columns"

    Parameters
    ----------
    df : Pandas dataframe
        data frame set you want to plot.
    model : string
        Can be ether "USZ" or "BDS", depending where to go
    columns : string
        Name of the columns you want to plot
    color : 3D or 4D array of float (between 0-1)
        color of the bar plot which will be displaied
    centered : string, optional
        Can be either "c", "r", "l". if the plot is offcentred on the right or
        on the leftThe default is "c".

    Returns
    -------
    bool
        if function ex3ecute correctly return True

    """
    if model == "USZ":
        list_column = ["year", "month"]
    elif model == "BDS":
        list_column = ["hasCovid", "year", "month"]
    else:
        RuntimeError(f"Please enter model easer USZ or BDS not {model}")
    temp_df = df.loc[df["x axis violin name"].drop_duplicates().index]
    temp_df = temp_df.sort_values(by=list_column)
    list_name = temp_df.sort_values(by=list_column).loc[:,
                                                        "x axis violin name"]
    list_mean = []
    for name in list_name:
        mask_name = df.loc[:, "x axis violin name"] == name
        posterior_probability = df.loc[:,
                                       columns]
        posterior_probability[df.loc[:, columns].isna()] = 0
        list_mean.append(np.nanmean(df.loc[mask_name,
                                           columns]))
    if centered == "c":
        dx = 0
        barWidth = 0.8
    elif centered == "r":
        dx = -0.2
        barWidth = 0.4
    elif centered == "l":
        dx = 0.2
        barWidth = 0.4

    x1 = np.arange(0, len(list_name))
    plt.bar(x1+dx,
            list_mean,
            color=color,
            width=barWidth,
            label=columns)
    plt.xticks(x1, list_name)
    return x1+dx, list_mean


def plot_prevalence_USZ_vs_BDS(df_1, df_2, model1, model2):
    """
    Plot the prevalence of one model given the second one, and compare it with
    its own model.

    Parameters
    ----------
    df_1 : pandas dataframes
        DESCRIPTION.
    df_2 : pandas dataframes
        DESCRIPTION.
    model1 : String
        give the model of the dataframe 1 ("USZ" or "BDS").
    model2 : String
        give the model of the dataframe 2 ("USZ" or "BDS").

    Returns
    -------
    bool
        if function ex3ecute correctly return True

    """
    list_pEC50 = ["NC_pEC50", "Spike_pEC50", "RBD_pEC50"]
    df_2, parameter_2 = Calculate_QDA(df_2, model2, seed=111)
    df_1, parameter_1 = Calculate_QDA(df_1, model1, seed=111)
    castTblUsed, castTblUnused, knownPos, knownNeg = filter_cohort(df_2,
                                                                   True)
    df_2.loc[castTblUsed.index,
             f"QDA_from_{model1}_model"] = QDAvalue_system(
                                                    castTblUsed.loc[:,
                                                                    list_pEC50
                                                                    ],
                                                    np.nanmean(df_1.QDA_prob),
                                                    parameter_1
                                                           )
    plt.figure()
    plt.title(f"QDA model of {model2}")
    plot_prevalence_name(df_2,
                         model2,
                         "QDA_prob",
                         [0, 0.2, 0.8, 0.5],
                         "r")
    plot_prevalence_name(df_2,
                         model2,
                         f"QDA_from_{model1}_model",
                         [1, 0, 0, 0.5],
                         "l")
    plt.legend()
    plt.savefig("QDA_{model1}vs{model2}")
    return True


def age_pyramide(df_USZ,
                 df_BDS,
                 name_csv="*KANTON_ZUERICH_bevoelkerung_1jahresklassen.csv",
                 path_folder=False):
    """
    Plot Age pyramid from both cohort and compare with the age pyramid of
    Kanton of Zurich. The default Zurich age pyramid is provided by the federal
    office of the statistic and give 2018 data.

    Parameters
    ----------
    df_USZ : pandas data frame
        USZ data frame which has been generated by the function
        generate_dataframe( xxxx ).
    df_BDS : pandas data frame
        BDS data frame which has been generated by the function
        generate_dataframe( xxxx ).
    name_csv : String, optional
        Name of the file containing stat of Zurich canton population.
        The default is "*KANTON_ZUERICH_bevoelkerung_1jahresklassen.csv".
    path_folder : string, optional
        path where the .csv file stand. if fault return the current folder
        The default is False.

    Returns
    -------
    bool
        if no error return True.

    """

    if path_folder is False:
        path_folder = path_to_main_folder()
    path_csv = glob(os.path.join(path_folder,
                                 "Code",
                                 "*",
                                 name_csv),
                    recursive=True)
    Zurich = pd.read_csv(path_csv[0],
                         sep=";")
    mask_men_ZH = Zurich["GESCHLECHT"] == "Mann"
    mask_women_ZH = Zurich["GESCHLECHT"] == "Frau"
    mask_year_ZH = Zurich["JAHR"] == 2019
    men = np.zeros(len(Zurich["ALTERSKLASSE_CODE"].drop_duplicates()))
    women = np.zeros(len(Zurich["ALTERSKLASSE_CODE"].drop_duplicates()))
    list_ages = Zurich["ALTERSKLASSE_CODE"].drop_duplicates().sort_values()
    for age in list_ages:
        mask_age_ZH = Zurich["ALTERSKLASSE_CODE"] == age
        men[age] = np.sum(Zurich.loc[mask_men_ZH & mask_year_ZH & mask_age_ZH,
                                     "ANZAHL_PERSONEN"])
        women[age] = np.sum(Zurich.loc[((mask_women_ZH &
                                         mask_year_ZH) &
                                        mask_age_ZH),
                                       "ANZAHL_PERSONEN"])

    df_USZ, parameter = Calculate_QDA(df_USZ,
                                      "USZ",
                                      ["NC_pEC50", "Spike_pEC50", "RBD_pEC50"],
                                      seed=111)
    df_BDS, parameter = Calculate_QDA(df_BDS,
                                      "BDS",
                                      ["NC_pEC50", "Spike_pEC50", "RBD_pEC50"],
                                      seed=111)

    plt.figure()
    df_2 = df_USZ
    maks_copendemic = df_2["year"] > 2018
    df_2 = df_2.loc[maks_copendemic]
    maskW = df_USZ.sex == "female"
    span = np.arange(np.min(df_2.loc[:, "age"]),
                     np.max(df_2.loc[:, "age"]), 2)
    x_M, y_M = np.histogram(df_2.loc[~maskW, "age"], span)
    x_W, y_W = np.histogram(df_2.loc[maskW, "age"], span)
    df_pyramid = pd.DataFrame({'Age': span[:-1],
                               'Male': -x_M,
                               'Female': x_W})
    plt.bar(df_pyramid["Age"]+1,
            df_pyramid["Male"] / 40,
            label="Male density (counts/40)",
            width=1.8)
    plt.bar(df_pyramid["Age"]+1,
            df_pyramid["Female"] / 40,
            label="Female density (counts/40)",
            width=1.8)

    maskP = df_USZ.loc[:, "QDA_prob"] >= 0.5
    x_M, y_M = np.histogram(df_2.loc[~maskW & maskP, "age"], span)
    x_W, y_W = np.histogram(df_2.loc[maskW & maskP, "age"], span)
    df_pyramid = pd.DataFrame({'Age': span[:-1],
                               'Male': -x_M,
                               'Female': x_W})

    plt.bar(df_pyramid["Age"]+1,
            df_pyramid["Male"],
            label="Male detected positive",
            width=1.8)
    plt.bar(df_pyramid["Age"]+1,
            df_pyramid["Female"],
            label="Female detected positive",
            width=1.8)
    plt.xlabel("Age")
    plt.ylabel("nbr")
    plt.plot(list_ages,
             -men / (np.sum(men) + np.sum(women)) * 1000,
             c="#165782ff",
             label="zurich overall population for thousand")
    plt.plot(list_ages,
             women / (np.sum(men) + np.sum(women)) * 1000,
             c="#a04f08ff")
    plt.legend()

    plt.figure()
    # df_blood = extract_date_Hnumber(df_blood)
    df_blood_2 = df_BDS[df_BDS.loc[:, "Geschlecht"].notna()]
    df_blood_2.loc[:, "age"] = (2020 -
                                pd.DatetimeIndex(df_blood_2.loc[:,
                                                                'Geburtsdatum']
                                                 ).year)
    mask_pos = df_blood_2["x axis violin name"] == "blood donation\n positif"
    df_2 = df_blood_2[~mask_pos]
    maskW = df_2.Geschlecht == "W"
    span = np.arange(np.min(df_2.loc[:, "age"]),
                     np.max(df_2.loc[:, "age"]),
                     2)
    x_M, y_M = np.histogram(df_2.loc[~maskW, "age"], span)
    x_W, y_W = np.histogram(df_2.loc[maskW, "age"], span)
    df_pyramid = pd.DataFrame({'Age': span[:-1],
                               'Male': -x_M,
                               'Female': x_W})
    plt.bar(df_pyramid["Age"] + 1,
            df_pyramid["Male"] / 5,
            label="Male density (counts/5)",
            width=1.8)
    plt.bar(df_pyramid["Age"] + 1,
            df_pyramid["Female"] / 5,
            label="Female density (counts/5)",
            width=1.8)

    maskP = df_2.loc[:, "QDA_prob"] >= 0.5
    x_M, y_M = np.histogram(df_2.loc[~maskW & maskP, "age"], span)
    x_W, y_W = np.histogram(df_2.loc[maskW & maskP, "age"], span)
    df_pyramid = pd.DataFrame({'Age': span[:-1],
                               'Male': -x_M,
                               'Female': x_W})
    plt.bar(df_pyramid["Age"] + 1,
            df_pyramid["Male"],
            label="Male detected positive",
            width=1.8)
    plt.bar(df_pyramid["Age"] + 1,
            df_pyramid["Female"],
            label="Female detected positive",
            width=1.8)
    plt.xlabel("Age")
    plt.ylabel("nbr")
    plt.plot(list_ages,
             -men / (np.sum(men) + np.sum(women)) * 1000,
             c="#165782ff",
             label="zurich overall population for thousand")
    plt.plot(list_ages,
             women / (np.sum(men)+np.sum(women)) * 1000,
             c="#a04f08ff")
    plt.legend()
    return True


def plot_clinic_repartission(df_USZ):
    """
    Plot clinic repartission of the blood sample, as well as the blood which
    has been detected as positive (QDA_prob >= 0.5 ).
    Work only for specific USZ data frame

    Parameters
    ----------
    df_USZ : pandas data frame
        USZ data frame which has been generated by the function
        generate_dataframe( xxxx ).

    Returns
    -------
    bool
        Return true if no error.

    """
    repartission = df_USZ.pivot_table(index=['clinic'], aggfunc='size')
    repartission = repartission.sort_values(ascending=False)
    list_clinic = repartission.index
    hist_clinic = repartission.values
    plt.figure()
    plt.semilogy()
    plt.xticks(rotation=45, ha='right')
    plt.bar(list_clinic, hist_clinic)
    try:
        mask_positive = df_USZ.QDA_prob >= 0.5
    except: # TODO find a way to ask if df_USZ.QDA_prob exist (like except User.DoesNotExist:)
        RuntimeError("QDA_prob should be calculate first")
    repartission = df_USZ.loc[mask_positive]
    repartission = repartission.pivot_table(index=['clinic'], aggfunc='size')
    repartission = repartission.sort_values(ascending=False)
    list_clinic = repartission.index
    hist_clinic = repartission.values
    plt.bar(list_clinic, hist_clinic)
    plt.ylabel("Occurance")
    return True


def plot_Origin_canton(df_USZ,
                       func_generate_file=dependancy):
    """
    Plot the origin of the canton,
    Becarefull the adhoc file is big and take time to open just to read
    the kanton of origin

    Parameters
    ----------
    df_USZ : pandas data frame
        USZ data frame which has been generated by the function
        generate_dataframe( xxxx ).
    func_generate_file : function, optional
        function which return the different files dependancy.
        The default is the function dependancy.
        but should be replace by a function which take value from the graphic
        interface

    Returns
    -------
    bool
        True

    """
    (main_folder,
     name_file_adhoc,
     name_file_adhoc_pos,
     name_QC,
     name_files_BDS) = func_generate_file()
    path_adhoc = glob(os.path.join(main_folder,
                                   "**",
                                   name_file_adhoc),
                      recursive=True)
    xls = pd.ExcelFile(path_adhoc[0])
    serie_kanton = pd.read_excel(xls, 'kantone')
    path_adhoc_pos = glob(os.path.join(main_folder,
                                       "**",
                                       name_file_adhoc_pos),
                          recursive=True)
    xls = pd.ExcelFile(path_adhoc_pos[0])
    serie_kanton_pos = pd.read_excel(xls, 'kantone')

    plt.figure()
    plt.semilogy()
    plt.xticks(rotation=45, ha='right')
    plt.bar(serie_kanton.canton, serie_kanton.sample_count)
    plt.bar(serie_kanton_pos.canton, serie_kanton_pos.sample_count)
    plt.ylabel("Occurence")
    return True


def plot_scatterQDA(df,
                    model,
                    IgG1="Spike_pEC50",
                    IgG2="RBD_pEC50",
                    color_parameter="QDA_prob"):
    """
    Plot the scatter plot of a specific dataframe following
    two columns

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.
    model : TYPE
        DESCRIPTION.
    IgG1 : String name of columns df, optional
        Can be a columns usualy NC_pEC50,S_pEC50 and RBD_pEC50.
        The default is "Spike_pEC50".
    IgG2 : String name of columns df, optional
        Can be a columns usualy NC_pEC50,S_pEC50 and RBD_pEC50.
        The default is "Spike_pEC50".
    color_parameter : String name of columns df, optional
        Can be a columns usualy QDA_prob but can be LDA_prob or other pEC50.
        The default is "QDA_prob".

    Returns
    -------
    bool
        True.

    """

    mask_pos = df.loc[:, "hasCovid"] == 1
    mask_neg = df.loc[:, "hasCovid"] == -1

    plt.figure()
    plt.scatter(df.loc[:, IgG1],
                df.loc[:, IgG2],
                c=df.loc[:, color_parameter],
                label="all value")
    plt.colorbar()
    plt.scatter(df.loc[mask_neg, IgG1],
                df.loc[mask_neg, IgG2],
                marker=",",
                color=[0.8, 0.8, 0.8, 0.3])
    plt.scatter(df.loc[mask_pos, IgG1],
                df.loc[mask_pos, IgG2],
                marker="s",
                color=[1, 0.3, 0.3, 0.8])
    plt.xlabel(IgG1)
    plt.ylabel(IgG2)
    plt.title(f"scatter plot of the {model} data," +
              f" colorbar are {color_parameter}")
    return True


def bootstrap(df, model, model_fonction, n_iterations=1000,
              path_save="QDA_error.csv", seed=111):
    """
    Calculate and save all the different prevalence per x axis name.
    helping to have an idea of the error of our graph. The bootstrap method
    is used.

    Parameters
    ----------
    df : dataframe of QDA (BDS, or USZ)
        DESCRIPTION.
    model : string
        Should be either "USZ" or "BDS".
    model_fonction : string
        Should be either "QDA" or "LDA".
    n_iterations : int
        number of iteration.
    path_save : string
        Name of the file
    seed : int
        seed used for setting the random generator.

    Returns
    -------
    tuple of 2 array
    First is a 1D array string name x-axis array
    second is a 2D array (size(x-axis),n_iteration) of the prevalence found

    """
    if model_fonction == "QDA":
        probability_field = Probability_field_FMW_bootstraps
    if model_fonction == "LDA":
        probability_field = Probability_field_LDA_bootstraps
    np.random.seed(seed)
    statistics = []
    subcohort = []
    for i in np.arange(0, n_iterations):  # bootstraps:
        castTblUsed, castTblUnused, knownPos, knownNeg = filter_cohort(df,
                                                                       True)
        mask_pos = castTblUsed["hasCovid"] == 1
        mask_neg = castTblUsed["hasCovid"] == -1
        mask_neutral = castTblUsed["hasCovid"] == 0
        sample_pos = resample(castTblUsed[mask_pos],
                              replace=True, random_state=i)
        sample_neg = resample(castTblUsed[mask_neg],
                              replace=True, random_state=i)
        sample_neutral = resample(castTblUsed[mask_neutral],
                                  replace=True, random_state=i)
        sample = pd.concat([sample_pos, sample_neg, sample_neutral])
        stat = bootstrap_prevalence(sample, castTblUnused,
                                    knownPos, knownNeg, model, df,
                                    probability_field)
        statistics.append(stat[1])
        subcohort.append(stat[0])
        print(f"bootstrap is at {i}")
        print(f"mean = {np.mean(np.asarray(statistics)[:,1])}, percentile 25-75. {np.percentile(np.asarray(statistics)[:,1],2.5)}-{np.percentile(np.asarray(statistics)[:,1],97.5)}")
        np.savetxt(path_save, statistics)
    return stat[0], np.asarray(statistics)


def bootstrap_prevalence(castTblUsed, castTblUnused,
                         knownPos, knownNeg, model, df, probability_field):
    """
    Calculate QDA of a resample df using .

    Parameters
    ----------
    castTblUsed : array
        DESCRIPTION.
    castTblUnused : array
        Should be either "USZ" or "BDS".
    knownPos : array
        number of iteration.
    knownNeg : array
        Name of the file
    model : string
        Should be either "USZ" or "BDS".
    df : dataframe of QDA (BDS, or USZ)
        used to know xaxis order
    probability_field : fonction for bootstrap QDA or LDA

    Returns
    -------
    tuple of 2 array
    First is a 1D array string name x-axis array
    second is a 2D array (size(x-axis),n_iteration) of the prevalence found

    """
    if model == "USZ":
        df = df.sort_values(["hasCovid", "year", "month"])
        mask_usz = df["cohort"] == "hospital patient"
        mask_KIM = ((df["subcohort"] == "KIM_neg") |
                    (df["subcohort"] == "KIM_pos"))
        list_xaxis = df.loc[mask_usz & ~mask_KIM,
                            "x axis violin name"].drop_duplicates()
    elif model == "BDS":
        df = df.sort_values(["hasCovid", "year", "month"])
        list_xaxis = df.loc[:, "x axis violin name"].drop_duplicates()
    else:
        RuntimeError(f"Please fill model either BDS or USZ not : {model}")
    intermediate_values = probability_field(castTblUsed,
                                            castTblUnused,
                                            knownPos,
                                            knownNeg)

    allVals, castTblUsed, castTblUnused, knownPos = intermediate_values
    prevalence = np.zeros(len(list_xaxis))
    prevalence_true = np.zeros(len(list_xaxis))
    for i, name in enumerate((list_xaxis)):
        mask_name = castTblUsed.loc[:, "x axis violin name"] == name

        mask_Unused = castTblUnused.loc[:, "x axis violin name"] == name
        y = allVals.loc[allVals.index.isin(castTblUsed.loc[mask_name].index)]
        pos_count = knownPos.loc[knownPos.index.isin(
                                castTblUsed.loc[mask_name].index)]
        prevalence[i] = np.sum(y.values) / (len(y)+np.sum(mask_Unused))
        prevalence_true[i] = len(pos_count) / (len(y)+np.sum(mask_Unused))
    return (list_xaxis, prevalence)


def plot_prevalence_blood(df_BDS):
    mask_bloodA = df_BDS.bloodGroup == "A"
    mask_bloodB = df_BDS.bloodGroup == "B"
    mask_blood0 = df_BDS.bloodGroup == "0"
    mask_bloodAB = df_BDS.bloodGroup == "AB"

    df_blood = df_BDS[mask_bloodA | mask_bloodB | mask_blood0 | mask_bloodAB]
    mask_hasCovid = df_blood.hasCovid == 0
    QDA_prob_mean = np.nanmean(df_blood.loc[mask_hasCovid,
                                            "QDA_prob"])

    # prevalenceA = np.nanmean(df_blood.loc[mask_hasCovid & mask_bloodA,
    #                                       "QDA_prob"])
    prevsampleA = sampling_riskfactor(df_blood,
                                      mask_hasCovid & mask_bloodA, 1000)
    prevsampleB = sampling_riskfactor(df_blood,
                                      mask_hasCovid & mask_bloodB, 1000)
    prevsample0 = sampling_riskfactor(df_blood,
                                      mask_hasCovid & mask_blood0, 1000)
    prevsampleAB = sampling_riskfactor(df_blood,
                                       mask_hasCovid & mask_bloodAB, 1000)
    dict_blood = {"A": (mask_bloodA, prevsampleA),
                  "B": (mask_bloodB, prevsampleB),
                  "0": (mask_blood0, prevsample0),
                  "AB": (mask_bloodAB, prevsampleAB)}
    mean = np.zeros(4)
    percentil_2p5 = np.zeros(4)
    percentil_97p5 = np.zeros(4)
    for i, name in enumerate(dict_blood):
        mean[i] = np.nanmean(df_blood.loc[mask_hasCovid & dict_blood[name][0],
                                          "QDA_prob"])
        percentil_2p5[i] = np.percentile(dict_blood[name][1], 2.5)
        percentil_97p5[i] = np.percentile(dict_blood[name][1], 97.5)
        plt.bar(name,
                mean[i]/QDA_prob_mean,
                label=name)
    plt.errorbar(np.arange(0, 4),
                 mean/QDA_prob_mean,
                 yerr=[np.abs(percentil_2p5-mean),
                       np.abs(percentil_97p5-mean)] / QDA_prob_mean,
                 marker=" ",
                 linestyle='None',
                 color=[0, 0, 0])
    plt.plot([-0.5, 3.5],
             [QDA_prob_mean/QDA_prob_mean,
              QDA_prob_mean/QDA_prob_mean], "--",
             c=[0.6, 0.6, 0.6], label="overall risk")
    plt.legend()
    print(f"Risk ratio for {name} = {QDA_prob_mean}")

    plt.legend()
    plt.figure()
    plt.title("blood Rh risk ratio")
    mask_bloodRhp = df_BDS.bloodRh == "+"
    mask_bloodRhm = df_BDS.bloodRh == "-"
    prevsampleRhp = sampling_riskfactor(df_blood,
                                        mask_hasCovid & mask_bloodRhp,
                                        1000)
    prevsampleRhm = sampling_riskfactor(df_blood,
                                        mask_hasCovid & mask_bloodRhm,
                                        1000)
    dict_blood = {"Rh+": (mask_bloodRhp, prevsampleRhp),
                  "Rh-": (mask_bloodRhm, prevsampleRhm)}
    mean = np.zeros(2)
    percentil_2p5 = np.zeros(2)
    percentil_97p5 = np.zeros(2)
    for i, name in enumerate(dict_blood):
        mean[i] = np.nanmean(df_blood.loc[mask_hasCovid & dict_blood[name][0],
                                          "QDA_prob"])
        percentil_2p5[i] = np.percentile(dict_blood[name][1], 2.5)
        percentil_97p5[i] = np.percentile(dict_blood[name][1], 97.5)
        plt.bar(name,
                mean[i]/QDA_prob_mean,
                label=name)

    plt.errorbar(np.arange(0, 2),
                 mean/QDA_prob_mean,
                 yerr=[np.abs(percentil_2p5-mean),
                       np.abs(percentil_97p5-mean)] / QDA_prob_mean,
                 marker=" ",
                 linestyle='None',
                 color=[0, 0, 0])
    plt.plot([-0.5, 1.5],
             [QDA_prob_mean/QDA_prob_mean,
              QDA_prob_mean/QDA_prob_mean], "--",
             c=[0.6, 0.6, 0.6], label="overall risk")
    plt.legend()
    return True


def sampling_riskfactor(df, mask_I, n_iteration):
    """
    Return the prevalence of a given sample
    sampling n_iteration times in order to calculate the risk factor given a
    control group mask_C
    """
    prevalenceSampleA = np.zeros(n_iteration)
    for i in np.arange(0, n_iteration):
        sampling_A = resample(df.loc[mask_I, "QDA_prob"],
                              replace=True, random_state=i)
        prevalenceSampleA[i] = np.nanmean(sampling_A)
    return prevalenceSampleA


def table_creation(df, model, name_save):
    df.loc[df.loc[:, "QDA_prob"].isna(), "QDA_prob"] = 0
    if model == "BDS":
        name_save = "BDS_" + name_save
        df.loc[:, 'Geburtsdatum'] = pd.to_datetime(df.loc[:,
                                                          'Geburtsdatum'])
        df.loc[:, 'Entnahmedatum'] = pd.to_datetime(df.loc[:,
                                                           'Entnahmedatum'])
        df.loc[:, 'age'] = (df.loc[:, 'Entnahmedatum'] -
                            df.loc[:, 'Geburtsdatum']).astype('<m8[Y]')
        mask_nan = df.Geschlecht.isna()
        df = df.loc[~mask_nan]
        mask_woman = df.Geschlecht == "W"
        mask_none = np.invert(mask_nan)
        mask_positive = df.QDA_prob > 0.5
        for i, mask in enumerate([mask_none, mask_positive]):
            median = np.median(df.loc[mask, "age"])
            median_w = np.median(df.loc[mask_woman & mask, "age"])
            median_m = np.median(df.loc[(~mask_woman) & mask, "age"])
            p2p5 = np.percentile(df.loc[mask, "age"], 2.5)
            p2p5_w = np.percentile(df.loc[mask & mask_woman, "age"], 2.5)
            p2p5_m = np.percentile(df.loc[mask & ~mask_woman, "age"], 2.5)
            p97p5 = np.percentile(df.loc[mask, "age"], 97.5)
            p97p5_w = np.percentile(df.loc[mask & mask_woman, "age"], 97.5)
            p97p5_m = np.percentile(df.loc[mask & ~mask_woman, "age"], 97.5)
            name = [" ", "median", "median_w", "median_m",
                    "p2p5", "p2p5_w", "p2p5_m",
                    "p97p5", "p97p5_w", "p97p5_m"]
            value = ["Age [year]", median, median_w, median_m,
                     p2p5, p2p5_w, p2p5_m,
                     p97p5, p97p5_w, p97p5_m]
            numbers = ["Numbers sample",
                       np.sum(mask),
                       np.sum(mask_woman & mask),
                       np.sum(np.invert(mask_woman) & mask),
                       np.sum(mask), np.sum(mask_woman & mask),
                       np.sum(np.invert(mask_woman) & mask),
                       np.sum(mask), np.sum(mask_woman & mask),
                       np.sum(np.invert(mask_woman) & mask)]
            if i == 0:
                print("None")
                np.savetxt("SampleAll_"+name_save,
                           np.asarray([name, value, numbers], str),
                           delimiter="\t",
                           fmt='%s')
            else:
                print("Positive")
                np.savetxt("SamplePositive_"+name_save,
                           np.asarray([name, value, numbers], str),
                           delimiter="\t", fmt='%s')

        df = df.loc[df.Spendernummer.drop_duplicates().index]
        name_save = "BDS_" + name_save
        df.loc[:, 'Geburtsdatum'] = pd.to_datetime(df.loc[:,
                                                          'Geburtsdatum'])
        df.loc[:, 'Entnahmedatum'] = pd.to_datetime(df.loc[:,
                                                           'Entnahmedatum'])
        df.loc[:, 'age'] = (df.loc[:, 'Entnahmedatum'] -
                            df.loc[:, 'Geburtsdatum']).astype('<m8[Y]')
        mask_nan = df.Geschlecht.isna()
        df = df.loc[~mask_nan]
        mask_woman = df.Geschlecht == "W"
        mask_none = np.invert(mask_nan)
        mask_positive = df.QDA_prob > 0.5
        for i, mask in enumerate([mask_none, mask_positive]):
            median = np.median(df.loc[mask, "age"])
            median_w = np.median(df.loc[mask_woman & mask, "age"])
            median_m = np.median(df.loc[(~mask_woman) & mask, "age"])
            p2p5 = np.percentile(df.loc[mask, "age"], 2.5)
            p2p5_w = np.percentile(df.loc[mask & mask_woman, "age"], 2.5)
            p2p5_m = np.percentile(df.loc[mask & ~mask_woman, "age"], 2.5)
            p97p5 = np.percentile(df.loc[mask, "age"], 97.5)
            p97p5_w = np.percentile(df.loc[mask & mask_woman, "age"], 97.5)
            p97p5_m = np.percentile(df.loc[mask & ~mask_woman, "age"], 97.5)
            name = [" ", "median", "median_w", "median_m",
                    "p2p5", "p2p5_w", "p2p5_m",
                    "p97p5", "p97p5_w", "p97p5_m"]
            value = ["Age [year]", median, median_w, median_m,
                     p2p5, p2p5_w, p2p5_m,
                     p97p5, p97p5_w, p97p5_m]
            numbers = ["Numbers sample",
                       np.sum(mask),
                       np.sum(mask_woman & mask),
                       np.sum(np.invert(mask_woman) & mask),
                       np.sum(mask), np.sum(mask_woman & mask),
                       np.sum(np.invert(mask_woman) & mask),
                       np.sum(mask), np.sum(mask_woman & mask),
                       np.sum(np.invert(mask_woman) & mask)]
            if i == 0:
                print("None")
                np.savetxt("Individual_All_"+name_save,
                           np.asarray([name, value, numbers], str),
                           delimiter="\t",
                           fmt='%s')
            else:
                print("Positive")
                np.savetxt("Individual_Positive_"+name_save,
                           np.asarray([name, value, numbers], str),
                           delimiter="\t",
                           fmt='%s')

    if model == "USZ":
        name_save = "USZ_" + name_save
        mask_nan = df.age.isna()
        df = df.loc[~mask_nan]
        mask_woman = df.sex == "female"
        mask_none = np.invert(mask_nan)
        mask_positive = df.QDA_prob > 0.5
        for i, mask in enumerate([mask_none, mask_positive]):
            median = np.median(df.loc[mask, "age"])
            median_w = np.median(df.loc[mask_woman & mask, "age"])
            median_m = np.median(df.loc[(~mask_woman) & mask, "age"])
            p2p5 = np.percentile(df.loc[mask, "age"], 2.5)
            p2p5_w = np.percentile(df.loc[mask & mask_woman, "age"], 2.5)
            p2p5_m = np.percentile(df.loc[mask & ~mask_woman, "age"], 2.5)
            p97p5 = np.percentile(df.loc[mask, "age"], 97.5)
            p97p5_w = np.percentile(df.loc[mask & mask_woman, "age"], 97.5)
            p97p5_m = np.percentile(df.loc[mask & ~mask_woman, "age"], 97.5)
            name = [" ", "median", "median_w", "median_m",
                    "p2p5", "p2p5_w", "p2p5_m",
                    "p97p5", "p97p5_w", "p97p5_m"]
            value = ["Age [year]", median, median_w, median_m,
                     p2p5, p2p5_w, p2p5_m,
                     p97p5, p97p5_w, p97p5_m]
            numbers = ["Numbers sample",
                       np.sum(mask),
                       np.sum(mask_woman & mask),
                       np.sum(np.invert(mask_woman) & mask),
                       np.sum(mask), np.sum(mask_woman & mask),
                       np.sum(np.invert(mask_woman) & mask),
                       np.sum(mask), np.sum(mask_woman & mask),
                       np.sum(np.invert(mask_woman) & mask)]
            if i == 0:
                print("None")
                np.savetxt("Sample_All_"+name_save,
                           np.asarray([name, value, numbers], str),
                           delimiter="\t",
                           fmt='%s')
            else:
                print("Positive")
                np.savetxt("Sample_Positive_"+name_save,
                           np.asarray([name, value, numbers], str),
                           delimiter="\t",
                           fmt='%s')
        df = df.loc[df.unique_PatientID.drop_duplicates().index]
        mask_nan = df.age.isna()
        df = df.loc[~mask_nan]
        mask_woman = df.sex == "female"
        mask_none = np.invert(mask_nan)
        mask_positive = df.QDA_prob > 0.5
        for i, mask in enumerate([mask_none, mask_positive]):
            median = np.median(df.loc[mask, "age"])
            median_w = np.median(df.loc[mask_woman & mask, "age"])
            median_m = np.median(df.loc[(~mask_woman) & mask, "age"])
            p2p5 = np.percentile(df.loc[mask, "age"], 2.5)
            p2p5_w = np.percentile(df.loc[mask & mask_woman, "age"], 2.5)
            p2p5_m = np.percentile(df.loc[mask & ~mask_woman, "age"], 2.5)
            p97p5 = np.percentile(df.loc[mask, "age"], 97.5)
            p97p5_w = np.percentile(df.loc[mask & mask_woman, "age"], 97.5)
            p97p5_m = np.percentile(df.loc[mask & ~mask_woman, "age"], 97.5)
            name = [" ", "median", "median_w", "median_m",
                    "p2p5", "p2p5_w", "p2p5_m",
                    "p97p5", "p97p5_w", "p97p5_m"]
            value = ["Age [year]", median, median_w, median_m,
                     p2p5, p2p5_w, p2p5_m,
                     p97p5, p97p5_w, p97p5_m]
            numbers = ["Numbers sample",
                       np.sum(mask),
                       np.sum(mask_woman & mask),
                       np.sum(np.invert(mask_woman) & mask),
                       np.sum(mask),
                       np.sum(mask_woman & mask),
                       np.sum(np.invert(mask_woman) & mask),
                       np.sum(mask),
                       np.sum(mask_woman & mask),
                       np.sum(np.invert(mask_woman) & mask)]
            if i == 0:
                print("None")
                np.savetxt("Individual_All_"+name_save,
                           np.asarray([name, value, numbers], str),
                           delimiter="\t",
                           fmt='%s')
            else:
                print("Positive")
                np.savetxt("Individual_Positive_"+name_save,
                           np.asarray([name, value, numbers], str),
                           delimiter="\t",
                           fmt='%s')

    return True


def creat_save_table(df_BDS, df_USZ, name_save):
    table_USZ = table_Marc_USZ(df_USZ)
    np.savetxt("table_USZ_Marc.tsv",
               np.asarray(table_USZ).T,
               delimiter="\t",
               fmt='%s')
    table_BDS = table_Marc_BDS(df_BDS)
    np.savetxt("table_BDS_Marc.tsv",
               np.asarray(table_BDS).T,
               delimiter="\t",
               fmt='%s')
    table_USZ = table_Marc_age_USZ(df_USZ)
    np.savetxt("table_USZ_Marc_age.tsv",
               np.asarray(table_USZ).T,
               delimiter="\t",
               fmt='%s')
    table_BDS = table_Marc_age_BDS(df_BDS)
    np.savetxt("table_BDS_Marc_age.tsv",
               np.asarray(table_BDS).T,
               delimiter="\t",
               fmt='%s')
    return True


def table_Marc_BDS(df):
    df.loc[df.loc[:, "QDA_prob"].isna(), "QDA_prob"] = 0
    df = df.loc[df.hasCovid == 0]
    # name_save = "BDS_"
    df.loc[:, 'Geburtsdatum'] = pd.to_datetime(df.loc[:,
                                                      'Geburtsdatum'])
    df.loc[:, 'Entnahmedatum'] = pd.to_datetime(df.loc[:,
                                                       'Entnahmedatum'])
    df.loc[:, 'age'] = (df.loc[:, 'Entnahmedatum'] -
                        df.loc[:, 'Geburtsdatum']).astype('<m8[Y]')
    mask_nan = df.Geschlecht.isna()
    df = df.loc[~mask_nan]
    mask_woman = df.Geschlecht == "W"
    mask_none = np.invert(mask_nan)
    mask_positive = df.QDA_prob > 0.5
    df_2 = df.loc[df.Spendernummer.drop_duplicates().index]
    mask_nan_2 = df_2.Geschlecht.isna()
    mask_woman_2 = df_2.Geschlecht == "W"
    mask_none_2 = np.invert(mask_nan_2)
    mask_positive_2 = df_2.QDA_prob > 0.5
    oldvalue = []
    name_i = {0: "total", 1: "positif"}
    mask2 = {0: mask_none_2, 1: mask_positive_2}
    table = [" ", "number samples", "number individual patientsORspenders",
             "median age samples", "age IQR samples",
             "median age individual", "age IQR individual",
             "Seroprevalence_QDA", "Stdev_QDA", "Seroprevalence_LDA",
             "Stdev_LDA",
             "number female samples",
             "number female individual patientsORspenders",
             "median age female samples", "age IQR  female samples",
             "median age  female individual", "age IQR  female individual",
             "Seroprevalence_QDA female", "Stdev_QDA female",
             "Seroprevalence_LDA female", "Stdev_LDA female",
             "number male samples",
             "number male individual patientsORspenders",
             "median age male samples",
             "age IQR male samples", "median age male individual",
             "age IQR male individual",
             "Seroprevalence_QDA male",
             "Stdev_QDA male", "Seroprevalence_LDA male", "Stdev_LDA male",
             "Positive ", "Positive number samples",
             "Positive number individual patients",
             "Positive median age samples",
             "Positive age IQR samples", "Positive median age individual",
             "Positive age IQR individual",
             "number Positive female samples",
             "number Positive female individual patientsORspenders",
             "Positive median age female samples",
             "Positive age IQR  female samples",
             "Positive median age  female individual",
             "Positive age IQR  female individual",
             "number Positive female samples",
             "number Positive female individual patientsORspenders",
             "Positive median age male samples",
             "Positive age IQR male samples",
             "Positive median age male individual",
             "Positive age IQR male individual"]
    for year in [3020]:
        print(year)
        mask_year = df.year == year
        # mask_year_2 = df_2.year == year
        for month in df.loc[mask_year].month.drop_duplicates():
            mask_month = df.month == month
            # mask_month_2 = df_2.month == month
            for i, mask in enumerate([mask_none, mask_positive]):
                mask = mask_year & mask_month & mask
                mask_2 = mask2[i] & mask_year & mask_month
                median = np.median(df.loc[mask, "age"])
                median_2 = np.median(df_2.loc[mask_2, "age"])
                median_w = np.median(df.loc[mask_woman & mask, "age"])
                median_w_2 = np.median(df_2.loc[mask_woman_2 & mask_2, "age"])
                median_m = np.median(df.loc[(~mask_woman) & mask, "age"])
                median_m_2 = np.median(df_2.loc[(~mask_woman_2) & mask_2,
                                                "age"])
                p2p5 = np.percentile(df.loc[mask, "age"], 2.5)
                if np.sum(mask_2) > 0:
                    p2p5_2 = np.percentile(df_2.loc[mask_2, "age"], 2.5)
                else:
                    p2p5_2 = np.nan
                p2p5_w = np.percentile(df.loc[mask & mask_woman, "age"], 2.5)
                if np.sum(np.sum(mask_2 & mask_woman_2)) > 0:
                    p2p5_w_2 = np.percentile(df_2.loc[mask_2 & mask_woman_2,
                                                      "age"], 2.5)
                else:
                    p2p5_w_2 = np.nan
                p2p5_m = np.percentile(df.loc[mask & ~mask_woman, "age"], 2.5)
                if np.sum(np.sum(mask_2 & ~mask_woman_2)) > 0:
                    p2p5_m_2 = np.percentile(df_2.loc[mask_2 & ~mask_woman_2,
                                                      "age"], 2.5)
                else:
                    p2p5_m_2 = np.nan
                p97p5 = np.percentile(df.loc[mask, "age"], 97.5)
                if np.sum(mask_2) > 0:
                    p97p5_2 = np.percentile(df_2.loc[mask_2, "age"], 97.5)
                else:
                    p97p5_2 = np.nan
                p97p5_w = np.percentile(df.loc[mask & mask_woman, "age"], 97.5)
                if np.sum(np.sum(mask_2 & mask_woman_2)) > 0:
                    p97p5_w_2 = np.percentile(df_2.loc[mask_2 & mask_woman_2,
                                                       "age"], 97.5)
                else:
                    p97p5_w_2 = np.nan
                p97p5_m = np.percentile(df.loc[mask & ~mask_woman, "age"],
                                        97.5)
                if np.sum(np.sum(mask_2 & ~mask_woman_2)) > 0:
                    p97p5_m_2 = np.percentile(df_2.loc[mask_2 & ~mask_woman_2,
                                                       "age"], 97.5)
                else:
                    p97p5_m_2 = np.nan
                if i == 0:
                    Seroprevalence_QDA = np.mean(df.loc[mask, "QDA_prob"])
                    Stdev_QDA = np.std(df.loc[mask, "QDA_prob"])
                    Seroprevalence_LDA = np.mean(df.loc[mask, "LDA_prob"])
                    Stdev_LDA = np.std(df.loc[mask, "LDA_prob"])

                    Seroprevalence_QDA_w = np.mean(df.loc[mask & mask_woman,
                                                          "QDA_prob"])
                    Stdev_QDA_w = np.std(df.loc[mask & mask_woman, "QDA_prob"])
                    Seroprevalence_LDA_w = np.mean(df.loc[mask & mask_woman,
                                                          "LDA_prob"])
                    Stdev_LDA_w = np.std(df.loc[mask & mask_woman, "LDA_prob"])

                    Seroprevalence_QDA_m = np.mean(df.loc[mask & ~mask_woman,
                                                          "QDA_prob"])
                    Stdev_QDA_m = np.std(df.loc[mask & ~mask_woman,
                                                "QDA_prob"])
                    Seroprevalence_LDA_m = np.mean(df.loc[mask & ~mask_woman,
                                                          "LDA_prob"])
                    Stdev_LDA_m = np.std(df.loc[mask & ~mask_woman,
                                                "LDA_prob"])

                    value = [f"{month}.{year-1000} {name_i[i]}",
                             np.sum(mask),
                             np.sum(mask_2),
                             median, f"{p2p5}-{p97p5}",
                             median_2, f"{p2p5_2}-{p97p5_2}",
                             Seroprevalence_QDA,
                             Stdev_QDA,
                             Seroprevalence_LDA,
                             Stdev_LDA,
                             np.sum(mask & mask_woman),
                             np.sum(mask_2 & mask_woman_2),
                             median_w,
                             f"{p2p5_w}-{p97p5_w}",
                             median_w_2,
                             f"{p2p5_w_2}-{p97p5_w_2}",
                             Seroprevalence_QDA_w,
                             Stdev_QDA_w,
                             Seroprevalence_LDA_w,
                             Stdev_LDA_w,
                             np.sum(mask & ~mask_woman),
                             np.sum(mask_2 & ~mask_woman_2),
                             median_m, f"{p2p5_m}-{p97p5_m}",
                             median_m_2, f"{p2p5_m_2}-{p97p5_m_2}",
                             Seroprevalence_QDA_m, Stdev_QDA_m,
                             Seroprevalence_LDA_m, Stdev_LDA_m]
                else:
                    value = [f"{month}.{year-1000} {name_i[i]}",
                             np.sum(mask), np.sum(mask_2),
                             median, f"{p2p5}-{p97p5}",
                             median_2, f"{p2p5_2}-{p97p5_2}",
                             np.sum(mask & mask_woman),
                             np.sum(mask_2 & mask_woman_2),
                             median_w, f"{p2p5_w}-{p97p5_w}",
                             median_w_2, f"{p2p5_w_2}-{p97p5_w_2}",
                             np.sum(mask & ~mask_woman),
                             np.sum(mask_2 & ~mask_woman_2),
                             median_m, f"{p2p5_m}-{p97p5_m}",
                             median_m_2, f"{p2p5_m_2}-{p97p5_m_2}"]

                value_2 = np.concatenate((np.asarray(oldvalue),
                                          np.asarray(value)))
                oldvalue = value
            table = np.vstack((np.asarray(table), np.asarray(value_2)))
        return table


def table_Marc_USZ(df):
    df.loc[df.loc[:, "QDA_prob"].isna(), "QDA_prob"] = 0
    # df = df.loc[df.hasCovid==0]
    # name_save = "USZ_"
    mask_nan = df.sex.isna()
    df = df.loc[~mask_nan]
    mask_woman = df.sex == "female"
    mask_none = np.invert(mask_nan)
    mask_positive = df.QDA_prob > 0.5
    df_2 = df.loc[df.unique_PatientID.drop_duplicates().index]
    mask_nan_2 = df_2.sex.isna()
    mask_woman_2 = df_2.sex == "female"
    mask_none_2 = np.invert(mask_nan_2)
    mask_positive_2 = df_2.QDA_prob > 0.5
    oldvalue = []
    name_i = {0: "total", 1: "positif"}
    mask2 = {0: mask_none_2, 1: mask_positive_2}
    table = [" ", "number samples", "number individual patientsORspenders",
             "median age samples", "age IQR samples", "median age individual",
             "age IQR individual",
             "Seroprevalence_QDA", "Stdev_QDA", "Seroprevalence_LDA",
             "Stdev_LDA",
             "number female samples",
             "number female individual patientsORspenders",
             "median age female samples", "age IQR  female samples",
             "median age  female individual", "age IQR  female individual",
             "Seroprevalence_QDA female", "Stdev_QDA female",
             "Seroprevalence_LDA female", "Stdev_LDA female",
             "number male samples",
             "number male individual patientsORspenders",
             "median age male samples", "age IQR male samples",
             "median age male individual", "age IQR male individual",
             "Seroprevalence_QDA male", "Stdev_QDA male",
             "Seroprevalence_LDA male", "Stdev_LDA male",
             "Positive ", "Positive number samples",
             "Positive number individual patients",
             "Positive median age samples", "Positive age IQR samples",
             "Positive median age individual", "Positive age IQR individual",
             "number Positive female samples",
             "number Positive female individual patientsORspenders",
             "Positive median age female samples",
             "Positive age IQR  female samples",
             "Positive median age  female individual",
             "Positive age IQR  female individual",
             "number Positive female samples",
             "number Positive female individual patientsORspenders",
             "Positive median age male samples",
             "Positive age IQR male samples",
             "Positive median age male individual",
             "Positive age IQR male individual"]
    for year in [2019, 2020]:
        print(year)
        mask_year = df.year == year
        # mask_year_2 = df_2.year == year
        for month in df.loc[mask_year].month.drop_duplicates():
            mask_month = df.month == month
            # mask_month_2 = df_2.month == month
            for i, mask in enumerate([mask_none, mask_positive]):
                mask = mask_year & mask_month & mask
                mask_2 = mask2[i] & mask_year & mask_month
                median = np.median(df.loc[mask, "age"])
                median_2 = np.median(df_2.loc[mask_2, "age"])
                median_w = np.median(df.loc[mask_woman & mask, "age"])
                median_w_2 = np.median(df_2.loc[mask_woman_2 & mask_2, "age"])
                median_m = np.median(df.loc[(~mask_woman) & mask, "age"])
                median_m_2 = np.median(df_2.loc[(~mask_woman_2) & mask_2,
                                                "age"])
                if np.sum(mask) > 0:
                    p2p5 = np.percentile(df.loc[mask, "age"], 2.5)
                    p97p5 = np.percentile(df.loc[mask, "age"], 97.5)
                else:
                    p2p5 = np.nan
                    p97p5 = np.nan
                if np.sum(mask_2) > 0:
                    p2p5_2 = np.percentile(df_2.loc[mask_2, "age"], 2.5)
                else:
                    p2p5_2 = np.nan
                if np.sum(mask & mask_woman):
                    p2p5_w = np.percentile(df.loc[mask & mask_woman,
                                                  "age"], 2.5)
                    p97p5_w = np.percentile(df.loc[mask & mask_woman,
                                                   "age"], 97.5)
                else:
                    p2p5_w = np.nan
                    p97p5_w = np.nan
                if np.sum(np.sum(mask_2 & mask_woman_2)) > 0:
                    p2p5_w_2 = np.percentile(df_2.loc[mask_2 & mask_woman_2,
                                                      "age"], 2.5)
                else:
                    p2p5_w_2 = np.nan
                if np.sum(mask & ~mask_woman):
                    p2p5_m = np.percentile(df.loc[mask & ~mask_woman,
                                                  "age"], 2.5)
                    p97p5_m = np.percentile(df.loc[mask & ~mask_woman,
                                                   "age"], 97.5)
                else:
                    p2p5_m = np.nan
                    p97p5_m = np.nan
                if np.sum(np.sum(mask_2 & ~mask_woman_2)) > 0:
                    p2p5_m_2 = np.percentile(df_2.loc[mask_2 & ~mask_woman_2,
                                                      "age"], 2.5)
                else:
                    p2p5_m_2 = np.nan
                if np.sum(mask_2) > 0:
                    p97p5_2 = np.percentile(df_2.loc[mask_2, "age"], 97.5)
                else:
                    p97p5_2 = np.nan
                if np.sum(np.sum(mask_2 & mask_woman_2)) > 0:
                    p97p5_w_2 = np.percentile(df_2.loc[mask_2 & mask_woman_2,
                                                       "age"], 97.5)
                else:
                    p97p5_w_2 = np.nan
                if np.sum(np.sum(mask_2 & ~mask_woman_2)) > 0:
                    p97p5_m_2 = np.percentile(df_2.loc[mask_2 & ~mask_woman_2,
                                                       "age"], 97.5)
                else:
                    p97p5_m_2 = np.nan
                if i == 0:
                    Seroprevalence_QDA = np.mean(df.loc[mask, "QDA_prob"])
                    Stdev_QDA = np.std(df.loc[mask, "QDA_prob"])
                    Seroprevalence_LDA = np.mean(df.loc[mask, "LDA_prob"])
                    Stdev_LDA = np.std(df.loc[mask, "LDA_prob"])

                    Seroprevalence_QDA_w = np.mean(df.loc[mask & mask_woman,
                                                          "QDA_prob"])
                    Stdev_QDA_w = np.std(df.loc[mask & mask_woman, "QDA_prob"])
                    Seroprevalence_LDA_w = np.mean(df.loc[mask & mask_woman,
                                                          "LDA_prob"])
                    Stdev_LDA_w = np.std(df.loc[mask & mask_woman, "LDA_prob"])

                    Seroprevalence_QDA_m = np.mean(df.loc[mask & ~mask_woman,
                                                          "QDA_prob"])
                    Stdev_QDA_m = np.std(df.loc[mask & ~mask_woman,
                                                "QDA_prob"])
                    Seroprevalence_LDA_m = np.mean(df.loc[mask & ~mask_woman,
                                                          "LDA_prob"])
                    Stdev_LDA_m = np.std(df.loc[mask & ~mask_woman,
                                                "LDA_prob"])

                    value = [f"{month}.{year} {name_i[i]}",
                             np.sum(mask), np.sum(mask_2),
                             median, f"{p2p5}-{p97p5}", median_2,
                             f"{p2p5_2}-{p97p5_2}",
                             Seroprevalence_QDA, Stdev_QDA,
                             Seroprevalence_LDA, Stdev_LDA,
                             np.sum(mask & mask_woman),
                             np.sum(mask_2 & mask_woman_2),
                             median_w, f"{p2p5_w}-{p97p5_w}",
                             median_w_2, f"{p2p5_w_2}-{p97p5_w_2}",
                             Seroprevalence_QDA_w,
                             Stdev_QDA_w, Seroprevalence_LDA_w, Stdev_LDA_w,
                             np.sum(mask & ~mask_woman),
                             np.sum(mask_2 & ~mask_woman_2),
                             median_m, f"{p2p5_m}-{p97p5_m}", median_m_2,
                             f"{p2p5_m_2}-{p97p5_m_2}",
                             Seroprevalence_QDA_m, Stdev_QDA_m,
                             Seroprevalence_LDA_m, Stdev_LDA_m]
                else:
                    value = [f"{month}.{year-1000} {name_i[i]}",
                             np.sum(mask), np.sum(mask_2),
                             median, f"{p2p5}-{p97p5}", median_2,
                             f"{p2p5_2}-{p97p5_2}",
                             np.sum(mask & mask_woman),
                             np.sum(mask_2 & mask_woman_2),
                             median_w, f"{p2p5_w}-{p97p5_w}",
                             median_w_2, f"{p2p5_w_2}-{p97p5_w_2}",
                             np.sum(mask & ~mask_woman),
                             np.sum(mask_2 & ~mask_woman_2),
                             median_m, f"{p2p5_m}-{p97p5_m}",
                             median_m_2, f"{p2p5_m_2}-{p97p5_m_2}"]
                value_2 = np.concatenate((np.asarray(oldvalue),
                                          np.asarray(value)))
                oldvalue = value
            table = np.vstack((np.asarray(table),np.asarray(value_2)))
    return table


def table_Marc_age_USZ(df, range_ages):
    df.loc[df.QDA_prob.isna(),"QDA_prob"]=0
    name_save = "USZ_"
    mask_nan = df.sex.isna()
    df = df.loc[~mask_nan]
    mask_woman = df.sex == "female"
    mask_none = np.invert(mask_nan)
    mask_positive = df.QDA_prob > 0.5
    # df = df.loc[df.unique_PatientID.drop_duplicates().index]
    mask_nan_2 = df.sex.isna()
    mask_woman_2 = df.sex == "female"
    mask_none_2 = np.invert(mask_nan_2)
    mask_positive_2 = df.QDA_prob > 0.5
    oldvalue = []
    name_i = {0: "total", 1: "positif"}
    mask2 = {0: mask_none_2, 1: mask_positive_2}
    if range_ages==0:
        list_ages = ((0, 18), (18, 31), (31, 41), (41, 51), (51, 61), (61, 71),
                     (71, 81), (81, 91), (91, 101), (101, 200))
        name = (
            "-18 seroprevalence",
        "-18 cumulative pos",
        "-18 # total",
        "-18 # ratio pos/tot",
        "18-30 seroprevalence",
        "18-30 cumulative pos",
        "18-30 # total",
        "18-30 # ratio pos/tot",
        "31-40 seroprevalence",
        "31-40 cumulative pos",
        "31-40 # total",
        "31-40 # ratio pos/tot",
        "41-50 seroprevalence",
        "41-50 cumulative pos",
        "41-50 # total",
        "41-50 # ratio pos/tot",
        "51-60 seroprevalence",
        "51-60 cumulative pos",
        "51-60 # total",
        "51-60 # ratio pos/tot",
        "61-70 seroprevalence",
        "61-70 cumulative pos",
        "61-70 # total",
        "61-70 # ratio pos/tot",
        "71-80 seroprevalence",
        "71-80 cumulative pos",
        "71-80 # total",
        "71-80 # ratio pos/tot",
        "81-90 seroprevalence",
        "81-90 cumulative pos",
        "81-90 # total",
        "81-90 # ratio pos/tot",
        "91-100 seroprevalence",
        "91-100 cumulative pos",
        "91-100 # total",
        "91-100 # ratio pos/tot",
        "100+ seroprevalence",
        "100+ cumulative pos",
        "100+ # total",
        "100+ # ratio pos/tot")
    if range_ages==1:
        list_ages = ((0, 18), (18, 30), (30, 40), (40, 50), (50, 60),
                     (60, 70), (70, 80), (80,200))
        name = (
            "-18 seroprevalence",
        "-18 cumulative pos",
        "-18 # total",
        "-18 # ratio pos/tot",
        "18-30 seroprevalence",
        "18-30 cumulative pos",
        "18-30 # total",
        "18-30 # ratio pos/tot",
        "30-40 seroprevalence",
        "30-40 cumulative pos",
        "30-40 # total",
        "30-40 # ratio pos/tot",
        "40-50 seroprevalence",
        "40-50 cumulative pos",
        "40-50 # total",
        "40-50 # ratio pos/tot",
        "50-60 seroprevalence",
        "50-60 cumulative pos",
        "50-60 # total",
        "50-60 # ratio pos/tot",
        "60-70 seroprevalence",
        "60-70 cumulative pos",
        "60-70 # total",
        "60-70 # ratio pos/tot",
        "70-80 seroprevalence",
        "70-80 cumulative pos",
        "70-80 # total",
        "70-80 # ratio pos/tot",
        "80+ seroprevalence",
        "80+ cumulative pos",
        "80+ # total",
        "80+ # ratio pos/tot")
    else:
        list_ages = ((0, 18), (18, 35), (35, 50), (50, 65), (65,200))
        name = (
            "-18 seroprevalence",
        "-18 cumulative pos",
        "-18 # total",
        "-18 # ratio pos/tot",
        "18-35 seroprevalence",
        "18-35 cumulative pos",
        "18-35 # total",
        "18-35 # ratio pos/tot",
        "35-50 seroprevalence",
        "35-50 cumulative pos",
        "35-50 # total",
        "35-50 # ratio pos/tot",
        "50-65 seroprevalence",
        "50-65 cumulative pos",
        "50-65 # total",
        "50-65 # ratio pos/tot",
        "65+ seroprevalence",
        "65+ cumulative pos",
        "65+ # total",
        "65+ # ratio pos/tot")

    oldvalue = []
    name = np.concatenate(([ nam + "_total" for nam in name],
                       [ nam + "_woman" for nam in name],
                       [ nam + "_man" for nam in name]))
    pandemic_time = np.asarray(df.loc[df.loc[:,"hasCovid"]!=-1,
                          ["year","month"]].drop_duplicates())
    for i, date in enumerate(pandemic_time):
        oldvalue = []
        for mask in [mask_none_2, mask_none_2 & mask_woman_2, mask_none_2 & ~mask_woman_2]:
            for age_range in list_ages:
                mask_age = df.age >= age_range[0]
                mask_age = mask_age & (df.age < age_range[1])
                mask_age = mask_age & (df.year==date[0])
                mask_age = mask_age & (df.month==date[1])

                seroprevalence = np.mean(df[mask & mask_age].QDA_prob)
                cumulative_pos = np.sum(df[mask & mask_age].QDA_prob>0.5)
                size = len(df[mask & mask_age].QDA_prob)
                ratio = cumulative_pos / size

                value = [seroprevalence, cumulative_pos, size, ratio]
                oldvalue = np.concatenate((np.asarray(oldvalue),
                                           np.asarray(value)))

        date_name = str(date[0]) + "-" + str(date[1])
        oldvalue = np.append(np.asarray(date_name), oldvalue)
        if i == 0:
            name = np.append(np.asarray(" "), name)
        name = np.vstack((np.asarray(name),np.asarray(oldvalue)))

    oldvalue = []
    for mask in [mask_none_2, mask_none_2 & mask_woman_2, mask_none_2 & ~mask_woman_2]:
        for age_range in list_ages:
            mask_age = df.age >= age_range[0]
            mask_age = mask_age & (df.age < age_range[1])
            mask_age = mask_age & (df.hasCovid ==-1) # prepandemic mask


            seroprevalence = np.mean(df[mask & mask_age].QDA_prob)
            cumulative_pos = np.sum(df[mask & mask_age].QDA_prob>0.5)
            size = len(df[mask & mask_age].QDA_prob)
            ratio = cumulative_pos / size

            value = [seroprevalence, cumulative_pos, size, ratio]
            oldvalue = np.concatenate((np.asarray(oldvalue),
                                       np.asarray(value)))

    date_name = "prepandemic"
    oldvalue = np.append(np.asarray(date_name), oldvalue)
    name = np.vstack((np.asarray(name),np.asarray(oldvalue)))

    oldvalue = []
    for mask in [mask_none_2, mask_none_2 & mask_woman_2, mask_none_2 & ~mask_woman_2]:
        for age_range in list_ages:
            mask_age = df.age >= age_range[0]
            mask_age = mask_age & (df.age < age_range[1])


            seroprevalence = np.mean(df[mask & mask_age].QDA_prob)
            cumulative_pos = np.sum(df[mask & mask_age].QDA_prob>0.5)
            size = len(df[mask & mask_age].QDA_prob)
            ratio = cumulative_pos / size

            value = [seroprevalence, cumulative_pos, size, ratio]
            oldvalue = np.concatenate((np.asarray(oldvalue),
                                       np.asarray(value)))

    date_name = "Overall data"
    oldvalue = np.append(np.asarray(date_name), oldvalue)

    name = np.vstack((np.asarray(name),np.asarray(oldvalue)))

    return name


def table_Marc_age_BDS(df, range_ages):
    df = df.loc[df.hasCovid==0]
    df.loc[df.QDA_prob.isna(),"QDA_prob"]=0
    name_save = "BDS_"
    df.loc[:, 'Geburtsdatum'] = pd.to_datetime(df.loc[:,
                                                      'Geburtsdatum'])
    df.loc[:, 'Entnahmedatum'] = pd.to_datetime(df.loc[:,
                                                       'Entnahmedatum'])
    df.loc[:, 'age'] = (df.loc[:, 'Entnahmedatum'] -
                            df.loc[:, 'Geburtsdatum']).astype('<m8[Y]')
    mask_nan = df.Geschlecht.isna()
    df = df.loc[~mask_nan]
    mask_woman = df.Geschlecht == "W"
    mask_none = np.invert(mask_nan)
    mask_positive = df.QDA_prob > 0.5
    # df = df.loc[df.Spendernummer.drop_duplicates().index]
    mask_nan_2 = df.Geschlecht.isna()
    mask_woman_2 = df.Geschlecht == "W"
    mask_none_2 = np.invert(mask_nan_2)
    mask_positive = df.QDA_prob > 0.5
    name_i = {0: "total", 1: "positif"}
    if range_ages==0:
        list_ages = ((0, 18), (18, 31), (31, 41), (41, 51), (51, 61), (61, 71),
                     (71, 81), (81, 91), (91, 101), (101, 200))
        name = (
            "-18 seroprevalence",
        "-18 cumulative pos",
        "-18 # total",
        "-18 # ratio pos/tot",
        "18-30 seroprevalence",
        "18-30 cumulative pos",
        "18-30 # total",
        "18-30 # ratio pos/tot",
        "31-40 seroprevalence",
        "31-40 cumulative pos",
        "31-40 # total",
        "31-40 # ratio pos/tot",
        "41-50 seroprevalence",
        "41-50 cumulative pos",
        "41-50 # total",
        "41-50 # ratio pos/tot",
        "51-60 seroprevalence",
        "51-60 cumulative pos",
        "51-60 # total",
        "51-60 # ratio pos/tot",
        "61-70 seroprevalence",
        "61-70 cumulative pos",
        "61-70 # total",
        "61-70 # ratio pos/tot",
        "71-80 seroprevalence",
        "71-80 cumulative pos",
        "71-80 # total",
        "71-80 # ratio pos/tot",
        "81-90 seroprevalence",
        "81-90 cumulative pos",
        "81-90 # total",
        "81-90 # ratio pos/tot",
        "91-100 seroprevalence",
        "91-100 cumulative pos",
        "91-100 # total",
        "91-100 # ratio pos/tot",
        "100+ seroprevalence",
        "100+ cumulative pos",
        "100+ # total",
        "100+ # ratio pos/tot")
    if range_ages==1:
        list_ages = ((0, 18), (18, 30), (30, 40), (40, 50), (50, 60),
                     (60, 70), (70, 80), (80,200))
        name = (
            "-18 seroprevalence",
        "-18 cumulative pos",
        "-18 # total",
        "-18 # ratio pos/tot",
        "18-30 seroprevalence",
        "18-30 cumulative pos",
        "18-30 # total",
        "18-30 # ratio pos/tot",
        "30-40 seroprevalence",
        "30-40 cumulative pos",
        "30-40 # total",
        "30-40 # ratio pos/tot",
        "40-50 seroprevalence",
        "40-50 cumulative pos",
        "40-50 # total",
        "40-50 # ratio pos/tot",
        "50-60 seroprevalence",
        "50-60 cumulative pos",
        "50-60 # total",
        "50-60 # ratio pos/tot",
        "60-70 seroprevalence",
        "60-70 cumulative pos",
        "60-70 # total",
        "60-70 # ratio pos/tot",
        "70-80 seroprevalence",
        "70-80 cumulative pos",
        "70-80 # total",
        "70-80 # ratio pos/tot",
        "80+ seroprevalence",
        "80+ cumulative pos",
        "80+ # total",
        "80+ # ratio pos/tot")
    else:
        list_ages = ((0, 18), (18, 35), (35, 50), (50, 65), (65,200))
        name = (
            "-18 seroprevalence",
        "-18 cumulative pos",
        "-18 # total",
        "-18 # ratio pos/tot",
        "18-35 seroprevalence",
        "18-35 cumulative pos",
        "18-35 # total",
        "18-35 # ratio pos/tot",
        "35-50 seroprevalence",
        "35-50 cumulative pos",
        "35-50 # total",
        "35-50 # ratio pos/tot",
        "50-65 seroprevalence",
        "50-65 cumulative pos",
        "50-65 # total",
        "50-65 # ratio pos/tot",
        "65+ seroprevalence",
        "65+ cumulative pos",
        "65+ # total",
        "65+ # ratio pos/tot")

    oldvalue = []
    name = np.concatenate(([ nam + "_total" for nam in name],
                       [ nam + "_woman" for nam in name],
                       [ nam + "_man" for nam in name]))
    pandemic_time = np.asarray(df.loc[df.loc[:,"year"]>=3019,
                          ["year","month"]].drop_duplicates())
    for i, date in enumerate(pandemic_time):
        oldvalue = []
        for mask in [mask_none_2, mask_none_2 & mask_woman_2, mask_none_2 & ~mask_woman_2]:
            for age_range in list_ages:
                mask_age = df.age >= age_range[0]
                mask_age = mask_age & (df.age < age_range[1])
                mask_age = mask_age & (df.year==date[0])
                mask_age = mask_age & (df.month==date[1])

                seroprevalence = np.mean(df[mask & mask_age].QDA_prob)
                cumulative_pos = np.sum(df[mask & mask_age].QDA_prob>0.5)
                size = len(df[mask & mask_age].QDA_prob)
                ratio = cumulative_pos / size

                value = [seroprevalence, cumulative_pos, size, ratio]
                oldvalue = np.concatenate((np.asarray(oldvalue),
                                           np.asarray(value)))

        date_name = str(date[0]) + "-" + str(date[1])
        oldvalue = np.append(np.asarray(date_name), oldvalue)
        if i == 0:
            name = np.append(np.asarray(" "), name)
        name = np.vstack((np.asarray(name),np.asarray(oldvalue)))

    oldvalue = []
    for mask in [mask_none_2, mask_none_2 & mask_woman_2, mask_none_2 & ~mask_woman_2]:
        for age_range in list_ages:
            mask_age = df.age >= age_range[0]
            mask_age = mask_age & (df.age < age_range[1])


            seroprevalence = np.mean(df[mask & mask_age].QDA_prob)
            cumulative_pos = np.sum(df[mask & mask_age].QDA_prob>0.5)
            size = len(df[mask & mask_age].QDA_prob)
            ratio = cumulative_pos / size

            value = [seroprevalence, cumulative_pos, size, ratio]
            oldvalue = np.concatenate((np.asarray(oldvalue),
                                       np.asarray(value)))

    date_name = "Overall data on the table"
    oldvalue = np.append(np.asarray(date_name), oldvalue)

    name = np.vstack((np.asarray(name),np.asarray(oldvalue)))

    return name


def path_to_main_folder():
    return os.path.dirname(os.path.abspath(''))


def __Example_plot__(mainfolder_path=False,
                     adhoc_path="*ptn_adhoc_202012*_annot*.xlsx",
                     QC_path="*Overview_data_transfer*.csv",
                     blood_path='../USZ_SARS/Meta_files/BDS/RunScript/',
                     saving_path=None):
    """
    Test of different plot of prevalence But does not save as pdf figure()

    Returns
    -------
    bool
        return True if no error

    """
    generate_dataframe(mainfolder_path,
                        adhoc_path,
                        QC_path,
                        blood_path,
                        saving_path)
    df, df_BDS, df_USZ = setup_existing_df(
        *('/Users/rjacquat/Documents/USZ/df.tsv',
          '/Users/rjacquat/Documents/USZ/df_BDS.tsv',
          '/Users/rjacquat/Documents/USZ/df_USZ.tsv'))
    df_USZ, parameter = Calculate_QDA(df_USZ, "USZ", seed=111)
    df_BDS, parameter = Calculate_QDA(df_BDS, "BDS", seed=111)
    df_USZ, parameter = Calculate_LDA(df_USZ, "USZ", seed=111)
    df_BDS, parameter = Calculate_LDA(df_BDS, "BDS", seed=111)

    print("Value of extremly negative case pEC50 smaller than -3")
    nbr_extneg_USZ_or = np.sum(np.logical_or(np.logical_or(df_USZ.NC_pEC50<-3,df_USZ.RBD_pEC50<-3),df_USZ.Spike_pEC50<-3))
    nbr_extneg_BDS_or = np.sum(np.logical_or(np.logical_or(df_BDS.NC_pEC50<-3,df_BDS.RBD_pEC50<-3),df_BDS.Spike_pEC50<-3))
    print(f"Nbr of sample concerned by pEC50 smaller than -3: {nbr_extneg_BDS_or+ nbr_extneg_USZ_or}")
    nbr_extneg_USZ_and = np.sum(np.logical_and(np.logical_and(df_USZ.NC_pEC50<-3,df_USZ.RBD_pEC50<-3),df_USZ.Spike_pEC50<-3))
    nbr_extneg_BDS_and = np.sum(np.logical_and(np.logical_and(df_BDS.NC_pEC50<-3,df_BDS.RBD_pEC50<-3),df_BDS.Spike_pEC50<-3))
    print(f"Nbr of sample with three value triger as low by pEC50 smaller than -3: {nbr_extneg_BDS_and+ nbr_extneg_USZ_and}")
    print(f"Nbr of Spike pEC50 smaller than -3: {np.sum(df_BDS.Spike_pEC50<-3)+ np.sum(df_USZ.Spike_pEC50<-3)}")
    print(f"Nbr of RBD pEC50 smaller than -3: {np.sum(df_BDS.RBD_pEC50<-3)+ np.sum(df_USZ.RBD_pEC50<-3)}")
    print(f"Nbr of NC pEC50 smaller than -3: {np.sum(df_BDS.NC_pEC50<-3)+ np.sum(df_USZ.NC_pEC50<-3)}")
    # a, stat = bootstrap(df_BDS, "BDS", "QDA", n_iterations = 1000,
    #                      path_save = "bds_QDA_error_order_2.csv", seed = 111)

    plt.figure()
    plt.title("QDA Prevalence of the USZ")
    plot_prevalence_name(df_USZ,
                         "USZ",
                         'QDA_prob',
                         [0.5, 0.6, 0.9],
                         centered="c")
    print("plot Prevalence USZ done")
    plt.figure()
    plt.title("QDA Prevalence of the BDS")
    x_axis, list_mean = plot_prevalence_name(df_BDS,
                                             "BDS",
                                             'QDA_prob',
                                             [0.9, 0.15, 0.15],
                                             centered="c")
    statistics = np.loadtxt("bds_QDA_error_order_2.csv")
    plt.errorbar(x_axis,
                 list_mean,
                 yerr=[-np.percentile(statistics,
                                      2.5, axis=0)+list_mean,
                       np.percentile(statistics,
                                     97.5, axis=0)-list_mean],
                 marker=" ",
                 linestyle='None',
                 color=[0, 0, 0])
    print("plot Prevalence BDS done")

    plt.figure()
    plt.title("LDA Prevalence of the USZ")
    x_axis, list_mean = plot_prevalence_name(df_USZ,
                                             "USZ",
                                             "LDA_prob",
                                             [0.5, 0.6, 0.9],
                                             centered="c")
    statistics = np.loadtxt("usz_LDA_error_order.csv")
    plt.errorbar(x_axis,
                  list_mean,
                  yerr=[-np.percentile(statistics,
                                      2.5, axis=0)+list_mean,
                        np.percentile(statistics,
                                      97.5, axis=0)-list_mean],
                  marker=" ",
                  linestyle='None',
                  color=[0, 0, 0])
    print("plot Prevalence USZ LDA done")

    plt.figure()
    plt.title("LDA Prevalence of the BDS")
    x_axis, list_mean = plot_prevalence_name(df_BDS,
                                             "BDS",
                                             "LDA_prob",
                                             [0.9, 0.15, 0.15],
                                             centered="c")
    statistics = np.loadtxt("bds_LDA_error_order.csv")
    plt.errorbar(x_axis,
                  list_mean,
                  yerr=[-np.percentile(statistics,
                                      2.5, axis=0)+list_mean,
                        np.percentile(statistics,
                                      97.5, axis=0)-list_mean],
                  marker=" ",
                  linestyle='None',
                  color=[0, 0, 0])
    print("plot Prevalence BDS LDA done")

    plot_prevalence_USZ_vs_BDS(df_USZ, df_BDS, "USZ", "BDS")
    plot_prevalence_USZ_vs_BDS(df_BDS, df_USZ, "BDS", "USZ")
    print("plot Prevalence using training other training dataset done")

    plt.figure()
    values = plot_prevalence_withknownPos_USZ(df_USZ, "QDA_prob")
    x_axis, list_mean, list_mean_pcrpos, list_mean_pos = values
    statistics = np.loadtxt("usz_QDA_error.csv")
    plt.errorbar(x_axis,
                 list_mean,
                 yerr=[-np.percentile(statistics, 2.5, axis=0)+list_mean,
                       np.percentile(statistics, 97.5, axis=0)-list_mean],
                 marker=" ",
                 linestyle='None',
                 color=[0, 0, 0])
    plt.title("prevalence add PCR pos")

    plt.figure()
    age_pyramide(df_USZ,
                 df_BDS,
                 name_csv="*KANTON_ZUERICH_bevoelkerung_1jahresklassen.csv",
                 path_folder=False)
    print("plot age pyramid done")

    plot_Origin_canton(df_USZ,
                       dependancy)
    print("plot origin of canton done")

    plot_clinic_repartission(df_USZ)
    print("plot clinic origin done")

    plot_scatterQDA(df_USZ, "USZ", "Spike_pEC50", "RBD_pEC50", "QDA_prob")
    plot_scatterQDA(df_BDS, "BDS", "Spike_pEC50", "RBD_pEC50", "QDA_prob")
    print("plot_scatter done")

    df.loc[df_BDS.index, "QDA_prob"] = df_BDS.loc[:, "QDA_prob"]
    df.loc[df_USZ.index, "QDA_prob"] = df_USZ.loc[:, "QDA_prob"]
    df_pos = df[df.loc[:, "QDA_prob"] >= 0.5]
    plt.figure()
    for method in ["NC", "Spike", "RBD"]:
        color = [0, 0.367, 0.656]
        plot_violin(df, df_pos, 'hospital patient', method, arg=(color))
        plt.xlim((0, 10.5))
        plt.ylim((-3, 6.5))

        color = [1, 0, 0]
        plot_violin(df, df_pos, 'blood donation', method, arg=(color))
        plt.xlim((0, 12.5))
        plt.ylim((-3, 4.5))
    # TODO ROC curve Plot
    # violin plot

    plt.figure()
    plt.title("Spendernumer Evolution covid level")
    df_NoDuplicate = df_BDS[~df_BDS.Spendernummer.isna()]
    mask_hascovid = df_NoDuplicate.hasCovid == 0
    for name in df_NoDuplicate.Spendernummer.drop_duplicates():
        QDA_prob = df_NoDuplicate.loc[((df_NoDuplicate.Spendernummer == name)
                                       & mask_hascovid)]
        QDA_prob.loc[QDA_prob.QDA_prob.isna(),"QDA_prob"]=0
        if len(QDA_prob) > 1:
            plt.plot(QDA_prob.month,QDA_prob.QDA_prob,"x-")
    plt.xlabel("month")
    plt.ylabel("QDA probability")

    plt.title("Spendernumer Evolution covid level")
    df_NoDuplicate = df_BDS[~df_BDS.Spendernummer.isna()]
    mask_hascovid = df_NoDuplicate.hasCovid == 0
    for name in df_NoDuplicate.Spendernummer.drop_duplicates():
        QDA_prob = df_NoDuplicate.loc[((df_NoDuplicate.Spendernummer == name)
                                       & mask_hascovid)]
        QDA_prob.loc[QDA_prob.QDA_prob.isna(),"QDA_prob"]=0
        if len(QDA_prob) > 1:
            plt.plot(QDA_prob.month,QDA_prob.QDA_prob,"x-")
    plt.xlabel("month")
    plt.ylabel("QDA probability")

    plt.figure()
    plt.title("blood group prevalence")
    plot_prevalence_blood(df_BDS)
    print("plot prevalence blood done")

    print("Start make table")
    table = table_Marc_age_USZ(df_USZ, 2)
    np.savetxt("table_mark_USZ_age_65plus.tsv",
               np.asarray(table).T, delimiter="\t", fmt='%s')

    table = table_Marc_age_USZ(df_USZ, 1)
    np.savetxt("table_mark_USZ_age_80plus.tsv",
               np.asarray(table).T,delimiter="\t", fmt='%s')

    table = table_Marc_age_BDS(df_BDS, 1)
    np.savetxt("table_mark_BDS_age_65plus.tsv",
               np.asarray(table).T,delimiter="\t", fmt='%s')

    table = table_Marc_age_BDS(df_BDS, 1)
    np.savetxt("table_mark_BDS_age_80plus.tsv",
               np.asarray(table).T,delimiter="\t", fmt='%s')

    print("END make table")

    return True
