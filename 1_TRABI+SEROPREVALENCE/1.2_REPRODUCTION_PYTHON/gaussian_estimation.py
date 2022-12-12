# -*- coding: utf-8 -*-
"""
Created on Tue May 19 16:41:02 2020

@author: Raphael Jacquat
"""
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import multivariate_normal as mvn
import os
from glob import glob
import scipy.stats
import math
from math import lgamma
from sklearn.utils import resample

list_pEC50 = ["NC_pEC50", "Spike_pEC50", "RBD_pEC50"]
CHOOSE_dop = False


def estimateFrac_gaussian(valueVec,mean1,sd1,mean2,sd2):
    pInit = 0.02
    density1 = norm.pdf(valueVec, mean1, sd1)
    density2 = norm.pdf(valueVec, mean2, sd2)
    myp = pInit
    for i in np.arange(1, 100):
        myt = ((myp * density1) /
               ((myp * density1) + ((1 - myp) * density2)))
        myp = np.nanmean(myt)
    return (myp, myt)


def plot_density_plot(df_USZ):
    # %%
    import scipy.stats as st
    from matplotlib import cm
    df_USZ_copy = df_USZ.copy()
    allVals, castTblUnused, knownPos = Probability_field_FMW(df_USZ_copy, True)
    df_USZ.loc[allVals.index,"Posterior"]= allVals
    x_neg = df_USZ.loc[df_USZ["hasCovid"] == -1, "RBD_pEC50"]
    y_neg = df_USZ.loc[df_USZ["hasCovid"] == -1, "Spike_pEC50"]
    x_0 = df_USZ.loc[:, "RBD_pEC50"]
    y_0 = df_USZ.loc[:, "Spike_pEC50"]
    colour_0 = df_USZ.loc[:, "Posterior"]
    colour_0.loc[colour_0.isna()] = 0
    x_pos = df_USZ.loc[df_USZ["hasCovid"] == 1, "RBD_pEC50"]
    y_pos = df_USZ.loc[df_USZ["hasCovid"] == 1, "Spike_pEC50"]
    knownPos = df_USZ.loc[df_USZ["hasCovid"] == 1,
                          ["RBD_pEC50", "Spike_pEC50"]]
    knownPos[knownPos < -3] = np.nan
    knownPos = knownPos.dropna(how='any')
    # mean1 = np.nanmean(knownPos, axis=0)
    # cov1 = np.cov(np.transpose(knownPos))
    # mean_cov1 = mvn(mean1, cov1)
    # knownNeg = df_USZ.loc[df_USZ["hasCovid"] == -1,
    #                       ["RBD_pEC50", "Spike_pEC50"]]
    x = df_USZ.loc[df_USZ["Posterior"] > 0.5, "RBD_pEC50"]
    y = df_USZ.loc[df_USZ["Posterior"] > 0.5, "Spike_pEC50"]
    values = np.vstack([x, y])
    kernel_pos = st.gaussian_kde(values)

    x = df_USZ.loc[df_USZ["Posterior"] < 0.01, "RBD_pEC50"]
    y = df_USZ.loc[df_USZ["Posterior"] < 0.01, "Spike_pEC50"]
    xmax = 10
    xmin = -3
    ymin = -3
    ymax = 10
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)


    # f_pos = np.reshape(kernel_pos(positions).T, xx.shape)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.gca()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.title('2D Gaussian Kernel density estimation')

    # viridis = cm.get_cmap('viridis', 1000)
    colormap = cm.viridis(colour_0)

    plt.scatter(x_0, y_0, marker="o", c=colormap, s=50)
    plt.scatter(x_pos, y_pos, marker="+", color=[0, 0, 0], linewidths=3, s=80)
    plt.scatter(x_neg, y_neg, marker='+', color=[1, 1, 1], linewidths=3, s=80)
    # cfset = ax.contourf(xx, yy, f, alpha=0.1, cmap='RdPu')
    # cfset = ax.contourf(xx, yy, f_pos,alpha=0.1, cmap='RdPu_r')

    ax.imshow(np.rot90(f), cmap='RdPu', extent=[xmin, xmax, ymin, ymax])
    # cset = ax.contour(xx, yy, f, colors='k')
    # cset = ax.contour(xx, yy, f-f_pos/5, colors='k')

    # ax.clabel(cset, inline=1, fontsize=10)
# %%
    return True


def plot_scatter(df_USZ):
    import scipy.stats as st
    from matplotlib import cm

    df_USZ_copy = df_USZ.copy()
    allVals, castTblUnused, knownPos = Probability_field_FMW(df_USZ_copy)
    df_USZ.loc[allVals.index, "Posterior"] = allVals
    x = df_USZ.loc[df_USZ["hasCovid"] == -1, "RBD_pEC50"]
    y = df_USZ.loc[df_USZ["hasCovid"] == -1, "Spike_pEC50"]
    x_0 = df_USZ.loc[df_USZ["hasCovid"] == 0, "RBD_pEC50"]
    y_0 = df_USZ.loc[df_USZ["hasCovid"] == 0, "Spike_pEC50"]
    colour_0 = df_USZ.loc[df_USZ["hasCovid"] == 0, "Posterior"]
    x_pos = df_USZ.loc[df_USZ["hasCovid"] == 1, "RBD_pEC50"]
    y_pos = df_USZ.loc[df_USZ["hasCovid"] == 1, "Spike_pEC50"]
    knownPos = df_USZ.loc[df_USZ["hasCovid"] == 1,
                          ["RBD_pEC50", "Spike_pEC50"]]
    knownPos[knownPos < -3] = np.nan
    knownPos = knownPos.dropna(how='any')
    mean1 = np.nanmean(knownPos, axis=0)
    cov1 = np.cov(np.transpose(knownPos))
    mean_cov1 = mvn(mean1, cov1)
    knownNeg = df_USZ.loc[df_USZ["hasCovid"] == -1,
                          ["RBD_pEC50", "Spike_pEC50"]]

    knownNeg[knownNeg < -3] = np.nan
    knownNeg = knownNeg.dropna(how='any')
    mean2 = np.nanmean(knownNeg, axis=0)
    cov2 = np.cov(np.transpose(knownNeg))
    mean_cov2 = mvn(mean2, cov2)

    myp = 0.9
    xmax = 5
    xmin = -1
    ymin = -1
    ymax = 5
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    # values = np.vstack([x, y])
    # kernel_neg = st.gaussian_kde(values)
    values_pos = np.vstack([x_pos, y_pos])
    # kernel_pos = st.gaussian_kde(values_pos)

    # kernel_pos = mean_cov1.pdf(knownPos)
    f = ((np.reshape((myp * mean_cov1.pdf(positions.T) /
                      (myp * mean_cov1.pdf(positions.T) +
                       (1 - myp) * mean_cov2.pdf(positions.T))).T,
                     xx.shape)))
    fig = plt.figure(figsize=(8, 8))
    ax = fig.gca()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    # cfset = ax.contourf(xx, yy, f, cmap='Purples')
    ax.imshow(np.rot90(f), cmap='Purples', extent=[xmin, xmax, ymin, ymax])
    cset = ax.contour(xx, yy, f, colors='k')
    ax.clabel(cset, inline=1, fontsize=10)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.title('2D Gaussian Kernel density estimation')
    plt.scatter(x, y, color=[0.2, 0.2, 0.2, 0.6])
    # viridis = cm.get_cmap('viridis', 1000)
    plt.scatter(x_0, y_0, marker="d", color=cm.viridis(colour_0 * 1000))
    plt.scatter(x_pos, y_pos, marker="x", color=[1, 1, 1])

    return True


def plot_scatter3D(df_USZ):
    # %%

    import scipy.stats as st
    from matplotlib import cm

    df_USZ_copy = df_USZ.copy()
    allVals, castTblUnused, knownPos = Probability_field_FMW(df_USZ_copy)
    df_USZ.loc[allVals.index, "Posterior"] = allVals
    x = df_USZ.loc[df_USZ["hasCovid"] == -1, "RBD_pEC50"]
    y = df_USZ.loc[df_USZ["hasCovid"] == -1, "Spike_pEC50"]
    x_0 = df_USZ.loc[df_USZ["hasCovid"] == 0, "RBD_pEC50"]
    y_0 = df_USZ.loc[df_USZ["hasCovid"] == 0, "Spike_pEC50"]
    colour_0 = df_USZ.loc[df_USZ["hasCovid"] == 0, "Posterior"]
    x_pos = df_USZ.loc[df_USZ["hasCovid"] == 1, "RBD_pEC50"]
    y_pos = df_USZ.loc[df_USZ["hasCovid"] == 1, "Spike_pEC50"]
    knownPos = df_USZ.loc[df_USZ["hasCovid"] == 1,
                          ["RBD_pEC50", "Spike_pEC50", "NC_pEC50"]]
    knownPos[knownPos < -3] = np.nan
    knownPos = knownPos.dropna(how='any')
    mean1 = np.nanmean(knownPos, axis=0)
    cov1 = np.cov(np.transpose(knownPos))
    mean_cov1 = mvn(mean1, cov1)
    knownNeg = df_USZ.loc[df_USZ["hasCovid"] == -1,
                          ["RBD_pEC50", "Spike_pEC50", "NC_pEC50"]]

    knownNeg[knownNeg < -3] = np.nan
    knownNeg = knownNeg.dropna(how='any')
    mean2 = np.nanmean(knownNeg, axis=0)
    cov2 = np.cov(np.transpose(knownNeg))
    mean_cov2 = mvn(mean2, cov2)
    myp = 0.1

    xmax = 5
    xmin = -1
    ymin = -1
    ymax = 5
    zmin = -1
    zmax = 5
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    x3D, y3D, zz = np.mgrid[xmin:xmax:100j, ymin:ymax:100j, zmin:zmax:100j]
    # positions = np.vstack([xx.ravel(), yy.ravel()])
    positions3D = np.vstack([x3D.ravel(), y3D.ravel(), zz.ravel()])
    # values = np.vstack([x, y])
    # kernel_neg = st.gaussian_kde(values)
    # values_pos = np.vstack([x_pos, y_pos])
    # kernel_pos = st.gaussian_kde(values_pos)

    # kernel_pos = mean_cov1.pdf(knownPos)

    f = (np.mean(np.reshape((myp * mean_cov1.pdf(positions3D.T) /
                             (myp*mean_cov1.pdf(positions3D.T) +
                              (1-myp)*mean_cov2.pdf(positions3D.T))).T,
                            zz.shape),
                 axis=2))
    fig = plt.figure(figsize=(10, 10))
    ax = fig.gca()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.imshow(np.rot90(f), cmap='Purples', extent=[xmin, xmax, ymin, ymax])
    cset = ax.contour(xx, yy, f, colors='k')
    ax.clabel(cset, inline=1, fontsize=10)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.title('2D Gaussian Kernel density estimation')
    plt.scatter(x, y, color=[0.2, 0.2, 0.2, 0.6])
    # viridis = cm.get_cmap('viridis', 1000)
    colormap = cm.viridis(colour_0 * 1000)
    colormap[:, -1] = np.log(colour_0 + 1) + 0.4
    colormap[colormap[:, -1] > 1, -1] = 1
    plt.scatter(x_0, y_0, marker="d", c=colormap)
    plt.scatter(x_pos, y_pos, marker="x", color=[1, 1, 1])
    # %%
    return True


def estimateFracTest():
    mean1 = 0.7
    mean2 = -0.1
    sd1 = 0.5
    sd2 = 0.3
    xx1 = np.random.normal(mean1, sd1, 200)
    xx2 = np.random.normal(mean2, sd2, 800) + np.random.normal(-12, sd2, 800)
    mean = np.mean(xx2)
    std = np.std(xx2)
    xx2 = np.append(np.random.normal(mean2, sd2, 400),
                    np.random.normal(-12, sd2, 400))
    valueVec = np.append(xx1,xx2)
    out = estimateFrac_gaussian(valueVec, mean1, sd1, mean, std)
    return out


def estimateFracFull(valueVec, knownPos, knownNeg):
    mean1 = np.nanmean(valueVec[knownPos])
    mean2 = np.nanmean(valueVec[knownNeg])
    sd1 = np.std(valueVec[knownPos])
    sd2 = np.std(valueVec[knownNeg])
    myp = 0.02
    for i in np.arange(1, 50):
        density1 = norm.pdf(valueVec, mean1, sd1)
        density2 = norm.pdf(valueVec, mean2, sd2)
        myt = ((myp * density1) /
               ((myp * density1) + ((1 - myp) * density2)))
        myt[knownPos] = 1
        myt[knownNeg] = 0
        myp = np.nanmean(myt[(knownPos == False) & (knownNeg == False)])
        mean1 = np.sum(myt * valueVec) / np.sum(myt)
        mean2 = np.sum((1 - myt) * valueVec) / np.sum(1 - myt)
        var1 = np.sum(myt*((valueVec - mean1) ** 2)) / np.sum(myt)
        var2 = np.sum((1 - myt) * ((valueVec - mean2) ** 2)) / np.sum(1-myt)
        sd1 = np.sqrt(var1)
        sd2 = np.sqrt(var2)

    paraml = (mean1, sd1, mean2, sd2)
    outl = (myp, myt, paraml)
# %%
    return(outl)


def estimateFracFullTest():
    mean1 = 0.7
    mean2 = -0.1
    sd1 = 0.2
    sd2 = 0.2
    xx1 = np.random.normal(mean1, sd1, 10000)
    xx2 = np.random.normal(mean2, sd2, 100000)
    valueVec = np.append(xx1, xx2)
    knownPos = np.repeat(False, len(valueVec))
    knownPos[:600] = True
    knownNeg = np.repeat(False, len(valueVec))
    knownNeg[-2000:] = True
    out = estimateFracFull(valueVec, knownPos, knownNeg)
    return(out)


def estimateFrac_gaussian(valueVec, mean1, sd1, mean2, sd2):
    pInit = 0.02
    density1 = norm.pdf(valueVec, mean1, sd1)
    density2 = norm.pdf(valueVec, mean2, sd2)
    myp = pInit
    for i in np.arange(1, 100):
        myt = ((myp * density1) /
               ((myp * density1) + ((1 - myp) * density2)))
        myp = np.nanmean(myt)
    return (myp, myt)


def estimateFracTest_NDim():
    mean1 = (1, 0.7, 0.3)
    mean2 = (0.2, -0.1, 0.0)
    cov1 = [[3, 0, 1], [0, 2, 0], [1, 0, 1]]
    cov2 = [[3, 0, 1], [0, 2, 0], [1, 0, 1]]
    xx1 = np.random.multivariate_normal(mean1, cov1, 800)
    xx2 = np.random.multivariate_normal(mean2, cov2, 200)
    valueVec = np.append(xx1, xx2, axis=0)
    out = estimateFrac_gaussianNDim(valueVec, mean1, cov1, mean2, cov2)
    return out


def estimateFrac_gaussianNDim(valueVec,mean1,cov1,mean2,cov2):
    pInit = 0.02
    mean_cov1 = mvn(mean1, cov1)
    mean_cov2 = mvn(mean2, cov2)
    density1 = mean_cov1.pdf(valueVec)
    density2 = mean_cov2.pdf(valueVec)
    myp = pInit
    for i in np.arange(1, 100):
        myt = ((myp * density1) /
               ((myp * density1) + ((1 - myp) * density2)))
        myp = np.nanmean(myt)
    return (myp, myt)


def estimateFrac_3D(df_tot):
    mask_pos = df_tot.loc[:, "hasCovid"] == +1
    mask_neg = df_tot.loc[:, "hasCovid"] == -1
    cov_pos = np.cov((df_tot.loc[mask_pos, "NC_pEC50"],
                      df_tot.loc[mask_pos, "Spike_pEC50"],
                      df_tot.loc[mask_pos, "RBD_pEC50"]))
    cov_neg = np.cov((df_tot.loc[mask_neg, "NC_pEC50"],
                      df_tot.loc[mask_neg, "Spike_pEC50"],
                      df_tot.loc[mask_neg, "RBD_pEC50"]))
    mean_pos = (np.nanmean(df_tot.loc[mask_pos, "NC_pEC50"]),
                np.nanmean(df_tot.loc[mask_pos, "Spike_pEC50"]),
                np.nanmean(df_tot.loc[mask_pos, "RBD_pEC50"]))
    mean_neg = (np.nanmean(df_tot.loc[mask_neg, "NC_pEC50"]),
                np.nanmean(df_tot.loc[mask_neg, "Spike_pEC50"]),
                np.nanmean(df_tot.loc[mask_neg, "RBD_pEC50"]))
    valueVec = (df_tot.loc[~(mask_pos | mask_neg), "NC_pEC50"],
                df_tot.loc[~(mask_pos | mask_neg), "Spike_pEC50"],
                df_tot.loc[~(mask_pos | mask_neg), "RBD_pEC50"])
    out = estimateFrac_gaussianNDim(valueVec,
                                    mean_pos,
                                    cov_pos,
                                    mean_neg,
                                    cov_neg)
    return out


def QDA(valueMat, knownPos, knownNeg):
    mean1 = np.nanmean(knownPos, axis=0)
    mean2 = np.nanmean(knownNeg, axis=0)
    cov1 = np.cov(np.transpose(knownPos))
    cov2 = np.cov(np.transpose(knownNeg))
    mean_cov1 = mvn(mean1, cov1)
    mean_cov2 = mvn(mean2, cov2)
    ldensity1 = np.log(mean_cov1.pdf(valueMat))
    ldensity2 = np.log(mean_cov2.pdf(valueMat))
    myt = ldensity1-ldensity2
    paraml = (mean1, cov1, mean2, cov2)
    outl = (myt, paraml)
    return outl


def getQDAFromFit(valueMat, fittedModel):
    mean1 = fittedModel[1][0]
    cov1 = fittedModel[1][1]
    mean2 = fittedModel[1][2]
    cov2 = fittedModel[1][3]
    mean_cov1 = mvn(mean1, cov1)
    mean_cov2 = mvn(mean2, cov2)
    ldensity1 = np.log(mean_cov1.pdf(valueMat))
    ldensity2 = np.log(mean_cov2.pdf(valueMat))
    myt = ldensity1 - ldensity2
    return myt


def estimateMVNFrac(valueMat, knownPos, knownNeg, nround=50):
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

    paraml = (mean_pos, cov1, mean_neg, cov2)
    outl = (myp, myt, paraml)
    return outl


def getPosteriorFromMVNFull(valueMat, fittedModel):
    myp = fittedModel[0]
    paraml = fittedModel[2]

    mean_pos = paraml[0]
    mean_neg = paraml[2]
    cov1 = paraml[1]
    cov2 = paraml[3]
    mean_cov1 = mvn(mean_pos, cov1)
    mean_cov2 = mvn(mean_neg, cov2)
    density1 = mean_cov1.pdf(valueMat)
    density2 = mean_cov2.pdf(valueMat)
    myt = (myp * density1) / ((myp * density1) + ((1 - myp) * density2))
    return myt


def Probability_field(castTblUsed, knownPos, knownNeg):
    # np.random.seed(seed=1011)
    posIndexVec = np.random.choice(np.repeat(np.arange(0, 10),
                                   np.ceil(len(knownPos) /
                                           10))[0:len(knownPos)],
                                   size=len(knownPos),
                                   replace=True,
                                   p=None)
    negIndexVec = np.random.choice(np.repeat(np.arange(0, 10),
                                   np.ceil(len(knownPos) /
                                           10))[0:len(knownNeg)],
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
        fit = QDA(mymat2, mat_pos_remove, mat_neg_remove)
        myposPosterior = getQDAFromFit(removedPos.loc[:, list_pEC50], fit)
        mynegPosterior = getQDAFromFit(removedNeg.loc[:, list_pEC50], fit)
        allVals[allVals.index.isin(knownPos.iloc[posIndexVec == myi].index)] = myposPosterior
        allVals[allVals.index.isin(knownNeg.iloc[negIndexVec == myi].index)] = mynegPosterior

    fullfit = QDA(castTblUsed.loc[:, list_pEC50],
                  knownPos.loc[:, list_pEC50],
                  knownNeg.loc[:, list_pEC50])
    unknown = ~(allVals.index.isin(knownPos.index) |
                allVals.index.isin(knownNeg.index))
    allVals[unknown] = fullfit[0][unknown]
    return allVals


def Probability_field_FMW(df_USZ, USEBLOOD=False):
    # np.random.seed(seed=111)
    castTblUsed, castTblUnused, knownPos, knownNeg  =  filter_cohort(df_USZ, USEBLOOD)

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
        fit = estimateMVNFrac(mymat2, mat_pos_remove, mat_neg_remove)
        myposPosterior = getPosteriorFromMVNFull(removedPos.loc[:, list_pEC50],
                                                 fit)
        mynegPosterior = getPosteriorFromMVNFull(removedNeg.loc[:, list_pEC50],
                                                 fit)
        allVals[allVals.index.isin(knownPos.iloc[posIndexVec == myi].index)] = myposPosterior
        allVals[allVals.index.isin(knownNeg.iloc[negIndexVec == myi].index)] = mynegPosterior

    fullfit = estimateMVNFrac(castTblUsed.loc[:, list_pEC50],
                              knownPos.loc[:, list_pEC50],
                              knownNeg.loc[:, list_pEC50],
                              nround=50)
    unknown = ~(allVals.index.isin(knownPos.index) |
                allVals.index.isin(knownNeg.index))
    allVals[unknown] = fullfit[1][unknown]
    allVals.to_frame().rename({"Spike_pEC50": "PosteriorProb"})
    return allVals, castTblUnused, knownPos


def Probability_field_FMW_given_prevalence(df, df_USZ,
                                           prevalence, USEBLOOD=False):

    #extract pos and neg for creation of covarience
    castTblUsed, castTblUnused, knownPos, knownNeg = filter_cohort(df_USZ,
                                                                   USEBLOOD)

    valueMat = df.loc[:, list_pEC50]
    valueMat[valueMat < -3] = -3
    mean_pos = np.mean(knownPos.loc[:, list_pEC50], axis=0)
    mean_neg = np.mean(knownNeg.loc[:, list_pEC50], axis=0)
    cov1 = np.cov(np.transpose(knownPos.loc[:, list_pEC50]))
    cov2 = np.cov(np.transpose(knownNeg.loc[:, list_pEC50]))
    # mask_pos = valueMat.index.isin(knownPos.loc[:, list_pEC50].index)
    # mask_neg = valueMat.index.isin(knownNeg.loc[:, list_pEC50].index)
    myp = prevalence

    mean_cov1 = mvn(mean_pos, cov1)
    mean_cov2 = mvn(mean_neg, cov2)
    density1 = mean_cov1.pdf(valueMat)
    density2 = mean_cov2.pdf(valueMat)
    myt = (myp * density1) / ((myp * density1) + ((1 - myp) * density2))
    df.loc[:, "posterior_fromUSZ"] = myt
    return df


def Probability_field_FMW_given_prevalence_2(df, sub_df,
                                             prevalence, method="USZ"):

    #e xtract pos and neg for creation of covarience
    castTblUsed, castTblUnused, knownPos, knownNeg = filter_cohort_2(sub_df,
                                                                     method)

    valueMat = df.loc[:, list_pEC50]
    valueMat[valueMat < -3] = -3
    mean_pos = np.mean(knownPos.loc[:, list_pEC50], axis=0)
    mean_neg = np.mean(knownNeg.loc[:, list_pEC50], axis=0)
    cov1 = np.cov(np.transpose(knownPos.loc[:, list_pEC50]))
    cov2 = np.cov(np.transpose(knownNeg.loc[:, list_pEC50]))
    mask_pos = valueMat.index.isin(knownPos.loc[:, list_pEC50].index)
    mask_neg = valueMat.index.isin(knownNeg.loc[:, list_pEC50].index)
    myp = prevalence

    mean_cov1 = mvn(mean_pos, cov1)
    mean_cov2 = mvn(mean_neg, cov2)
    density1 = mean_cov1.pdf(valueMat)
    density2 = mean_cov2.pdf(valueMat)
    for i in range(0, 50):
        myt = (myp * density1) / ((myp * density1) + ((1 - myp) * density2))
        myt[mask_pos] = 1
        myt[mask_neg] = 0
        myp = np.nanmean(myt[(~mask_pos) & (~mask_neg)])
    print(myp)
    df.loc[:, "posterior_fromUSZ"] = myt
    return df


def filter_cohort(df, USEBLOOD=False):
    df = df.loc[~df.RBD_pEC50.isnull()]
    for name in list_pEC50:
        mask = df.loc[:, name] < -3
        df.loc[mask, name] = -np.inf

    remove = df.loc[:, name] == "Take the entire as False first"
    for name in list_pEC50:
        remove = (df.loc[:, name] == -np.inf) | remove

    castTblUsed = df.loc[~remove]
    castTblUnused = df.loc[remove]

    if USEBLOOD:
        df.loc[df["x axis violin name"] == "BDS 2019\nDec", "hasCovid"] = -1
        df.loc[df["x axis violin name"] == "BDS 2020\nJan", "hasCovid"] = -1
        knownNeg = df[(~remove) & (df.loc[:, "hasCovid"] == -1)]
        knownPos = df[(~remove) & (df.loc[:, "hasCovid"] == 1)]
    else:
        if CHOOSE_dop:
            path_to_main_folder = os.path.dirname(os.path.abspath(''))
            paths_meta = glob(os.path.join(path_to_main_folder,
                                           "**",
                                           "*20200520_COVID_Spike_reduced*.xlsx"),
                              recursive=True)
            meta = pd.read_excel(paths_meta[0])
            myvals = meta.loc[meta.loc[:, "CoV2 time since onset"] == ">14"]

            mask_pos_dop14 = df.index.isin(myvals.PatientIDList)
            mask_pos_PCR = np.repeat(False, np.shape(df)[0])
        else:
            mask_pos_dop14 = np.repeat(False, np.shape(df)[0])
            mask_pos_PCR = df.loc[:, "hasCovid"] == 1
        knownPos = df.loc[(~remove) & ((mask_pos_dop14) | (mask_pos_PCR))]
        knownNeg = df.loc[(~remove) & (df.loc[:,"year"] < 2019)]
    return castTblUsed, castTblUnused, knownPos, knownNeg


def filter_cohort_2(df, method="USZ"):
    df = df.loc[~df.RBD_pEC50.isnull()]
    for name in list_pEC50:
        mask = df.loc[:, name] < -3
        df.loc[mask, name] = -np.inf
    remove = df.loc[:, name] == "Take the entire as False first"
    for name in list_pEC50:
        remove = (df.loc[:, name] == -np.inf) | remove

    castTblUsed = df.loc[~remove]
    castTblUnused = df.loc[remove]
    if method == "USZ":
        return filter_cohort(df, USEBLOOD=False)
    if method == "BDS":
        return filter_cohort(df, USEBLOOD=True)
    if method == "USZneg_BDSpos":
        mask_neg = df.loc[:, "subcohort"] == "hospital patient < 2019"
        df.loc[~mask_neg, "hasCovid"] = 0
        mask_pos = df.loc[:, "subcohort"] == "serum donation positive"
        df.loc[mask_pos, "hasCovid"] = 1
    knownNeg = df[(~remove) & (df.loc[:, "hasCovid"] == -1)]
    knownPos = df[(~remove) & (df.loc[:, "hasCovid"] == 1)]
    return castTblUsed, castTblUnused, knownPos, knownNeg


def estimateMVMixedFrac(valueMat, knownPos, knownNeg, cohort, nround=50):
    # gg=fit.mst(valueMat[knownNeg,,drop=FALSE], method = "BFGS")
    if cohort == "blood donation":
        mydf = 5.670253
    elif cohort == "hospital patient":
        mydf = 6.31327
    cov1 = np.cov(np.transpose(knownPos))
    cov2 = np.cov(np.transpose(knownNeg))
    sigma2 = (mydf - 2) / mydf * cov2
    mean1 = np.mean(knownPos, axis=0)
    mean2 = np.mean(knownNeg, axis=0)
    mean_cov1 = mvn(mean1, cov1)
    myp = 0.02
    mask_pos = valueMat.index.isin(knownPos.index)
    mask_neg = valueMat.index.isin(knownNeg.index)
    for i in np.arange(1, nround):
        density1 = mean_cov1.pdf(valueMat)
        density2 = dmvt(valueMat.values, mean2, sigma2, mydf, log=False)
        myt = (myp * density1) / ((myp * density1) + ((1 - myp) * density2))
        myt[mask_pos] = 1
        myt[mask_neg] = 0
        myp = np.nanmean(myt[(~mask_pos) & (~mask_neg)])
#        print(myp)

#    print("done")
    paraml = (mean1, cov1, mean2, sigma2, mydf)
    outl = (myp, myt, paraml)
    return outl


def getPosteriorFromVMixedFull(valueMat, fittedModel):
    myp = fittedModel[0]
    paraml = fittedModel[2]
    mean1 = paraml[0]
    cov1 = paraml[1]
    mean2 = paraml[2]
    sigma2 = paraml[3]
    mydf = paraml[4]
    mean_cov1 = mvn(mean1, cov1)
    density1 = mean_cov1.pdf(valueMat)
    density2 = dmvt(valueMat, mean2, sigma2, mydf, log=False)
    myt = (myp * density1) / ((myp * density1) + ((1 - myp) * density2))
    return myt


def dmvt(x, mu, Sigma, df, log):
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
    R_x_m = np.linalg.solve(dec, np.transpose(np.asarray(x) - np.asarray(mu)))
    rss = np.power(R_x_m, 2).sum(axis=0)
    logretval = (lgamma(1.0 * (p_3 + df) / 2) - (lgamma(1.0 * df / 2) +
                 np.sum(np.log(dec.diagonal())) + p_3/2 * np.log(math.pi * df))
                 - 0.5 * (df + p_3) * np.log1p((rss/df)))
    if log is False:
        return(np.exp(logretval))
    else:
        return logretval


def Probability_field_VMixedFull(df_USZ, USEBLOOD=False):
    if USEBLOOD is False:
        cohort = "hospital patient"
    else:
        cohort = "blood donation"
    castTblUsed, castTblUnused, knownPos, knownNeg  =  filter_cohort(df_USZ,
                                                                     USEBLOOD)
    posIndexVec = np.random.choice(np.repeat(np.arange(0, 10),
                                   np.ceil(len(knownPos)/10))[0:len(knownPos)],
                                   size=len(knownPos),
                                   replace=True,
                                   p=None)
    negIndexVec = np.random.choice(np.repeat(np.arange(0, 10),
                                   np.ceil(len(knownPos) /
                                           10))[0:len(knownNeg)],
                                   size=len(knownNeg),
                                   replace=True,
                                   p=None)
    # mymat = castTblUsed.loc[:, list_pEC50]
    allVals = castTblUsed.loc[:, "Spike_pEC50"]
    for myi in np.arange(0,10):
        removedPos = knownPos.iloc[posIndexVec == myi]
        removedNeg = knownNeg.iloc[negIndexVec == myi]
        torm = pd.concat([removedPos, removedNeg])
        mymat2 = castTblUsed.loc[~castTblUsed.index.isin(torm.index),
                                 list_pEC50]
        mat_pos_remove = knownPos.loc[~knownPos.index.isin(removedPos.index),
                                      list_pEC50]
        mat_neg_remove = knownNeg.loc[~knownNeg.index.isin(removedNeg.index),
                                      list_pEC50]
        fit=estimateMVMixedFrac(mymat2,
                                mat_pos_remove,
                                mat_neg_remove,
                                cohort,
                                nround=50)
        myposPosterior = getPosteriorFromVMixedFull(removedPos.loc[:,
                                                                   list_pEC50],
                                                    fit)
        mynegPosterior = getPosteriorFromVMixedFull(removedNeg.loc[:,
                                                                   list_pEC50],
                                                    fit)
        allVals[allVals.index.isin(knownPos.iloc[posIndexVec==myi].index)] = myposPosterior
        allVals[allVals.index.isin(knownNeg.iloc[negIndexVec==myi].index)] = mynegPosterior

    fullfit = estimateMVMixedFrac(castTblUsed.loc[:, list_pEC50],
                                  knownPos.loc[:, list_pEC50],
                                  knownNeg.loc[:, list_pEC50],
                                  cohort,
                                  nround=50)
    unknown = ~(allVals.index.isin(knownPos.index) |
                allVals.index.isin(knownNeg.index))
    allVals[unknown] = fullfit[1][unknown]

    return allVals, castTblUnused, knownPos


def Probability_field_FMW_bootstraps(castTblUsed,
                                     castTblUnused,
                                     knownPos,
                                     knownNeg):

    castTblUsed = castTblUsed.reset_index()
    castTblUnused = castTblUnused.reset_index()
    knownPos = castTblUsed[castTblUsed.loc[:,
                                           "Patient_or_Control_ID"].isin(knownPos.index)]
    knownNeg = castTblUsed[castTblUsed.loc[:,
                                           "Patient_or_Control_ID"].isin(knownNeg.index)]

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
        fit = estimateMVNFrac(mymat2, mat_pos_remove, mat_neg_remove)
        myposPosterior = getPosteriorFromMVNFull(removedPos.loc[:, list_pEC50],
                                                 fit)
        mynegPosterior = getPosteriorFromMVNFull(removedNeg.loc[:, list_pEC50],
                                                 fit)
        allVals.iloc[removedPos.index] = myposPosterior
        allVals.iloc[removedNeg.index] = mynegPosterior

    fullfit = estimateMVNFrac(castTblUsed.loc[:, list_pEC50],
                              knownPos.loc[:, list_pEC50],
                              knownNeg.loc[:, list_pEC50],
                              nround=50)
    unknown = ~(allVals.index.isin(knownPos.index) |
                allVals.index.isin(knownNeg.index))
    allVals[unknown] = fullfit[1][unknown]
#    print(fullfit[0])
#    print("Newseed")
    return allVals, castTblUsed, castTblUnused, knownPos


def estimateFrac_fix_Prior(df, prevalence, USEBLOOD=False):
    castTblUsed, castTblUnused, knownPos, knownNeg = filter_cohort(df,
                                                                   USEBLOOD)
    # mymat2 = castTblUsed.loc[:, list_pEC50]
    mat_pos_remove = knownPos.loc[:,
                                  list_pEC50]
    mat_neg_remove = knownNeg.loc[:,
                                  list_pEC50]
    valueMat = castTblUsed.loc[:,
                               list_pEC50]
    mean_pos = np.mean(mat_pos_remove, axis=0)
    mean_neg = np.mean(mat_neg_remove, axis=0)
    cov1 = np.cov(np.transpose(mat_pos_remove))
    cov2 = np.cov(np.transpose(mat_neg_remove))
    mask_pos = valueMat.index.isin(knownPos.index)
    mask_neg = valueMat.index.isin(knownNeg.index)
    myp = prevalence

    mean_cov1 = mvn(mean_pos, cov1)
    mean_cov2 = mvn(mean_neg, cov2)
    density1 = mean_cov1.pdf(valueMat)
    density2 = mean_cov2.pdf(valueMat)
    myt = (myp * density1) / ((myp * density1) + ((1 - myp) * density2))
    myp = np.nanmean(myt[(~mask_pos) & (~mask_neg)])

    paraml = (mean_pos, cov1, mean_neg, cov2)
    outl = (myp, myt, paraml)
    return outl


def estimateLDAFrac(valueMat,
                    knownPos,
                    knownNeg,
                    nround=50,
                    use1dim=True,
                    use1dimCor=True):
    mean1 = np.nanmean(knownPos, axis=0)
    mean2 = np.nanmean(knownNeg, axis=0)
    # cov1 = np.cov(np.transpose(knownPos))
    cov2 = np.cov(np.transpose(knownNeg))
    mean_cov1 = mvn(mean1, cov2)
    mean_cov2 = mvn(mean2, cov2)
    mask_pos = valueMat.index.isin(knownPos.index)
    mask_neg = valueMat.index.isin(knownNeg.index)
    myp = 0.02

    linearComb = np.linalg.solve(cov2, mean1-mean2)  # solve(cov2)%*%(mean1-mean2)
    tstat = valueMat.dot(linearComb)  # valueMat%*%linearComb
    linVar = np.transpose(mean1-mean2).dot(linearComb)  # (t(mean1-mean2)%*%solve(cov2)%*%(mean1-mean2))[1]
    linMean1 = np.transpose(mean1-mean2).dot(np.linalg.solve(cov2, mean1))  # (t(mean1-mean2)%*%solve(cov2)%*%(mean1))[1]
    linMean2 = np.transpose(mean1-mean2).dot(np.linalg.solve(cov2, mean2))  # (t(mean1-mean2)%*%solve(cov2)%*%(mean2))[1]
    sd1 = np.sqrt(np.var(tstat[knownPos.index]))
    sd2 = np.sqrt(np.var(tstat[knownNeg.index]))
    for i in np.arange(0, nround):
        if use1dim:
            density1_lin = norm.pdf(tstat, linMean1, np.sqrt(linVar)) #dnorm(tstat,mean=linMean1,sqrt(linVar))
            density2_lin = norm.pdf(tstat, linMean2, np.sqrt(linVar))

            if use1dimCor:
                    density1 = norm.pdf(tstat, np.mean(tstat[knownPos.index]),
                                        sd1)
                    density2 = norm.pdf(tstat, np.mean(tstat[knownNeg.index]),
                                        sd2)

            density1 = density1_lin
            density2 = density2_lin
        else:

           density1 = mean_cov1.pdf(valueMat)
           density2 = mean_cov2.pdf(valueMat)
           density1 = dmvnorm(valueMat,mean=mean1,sigma=cov2,log=FALSE)
           density2 = dmvnorm(valueMat,mean=mean2,sigma=cov2,log=FALSE)

        myt = (myp * density1) / ((myp * density1) + ((1 - myp) * density2))
        myt[mask_pos] = 1
        myt[mask_neg] = 0
        myp = np.mean(myt[(mask_pos == False) & (mask_neg == False)])

    cov1 = np.cov(np.transpose(knownPos)) #Cov1 is not used here but still send back for debug
    paraml = (mean1, cov1, mean2, cov2, linearComb, linMean1, linMean2, sd1, sd2)
    outl = (myp, myt, paraml, tstat)
    return(outl)


def Probability_field_LDA(df_USZ, USEBLOOD=False):
    np.random.seed(seed=111)
    castTblUsed, castTblUnused, knownPos, knownNeg = filter_cohort(df_USZ,
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
                              nround=50,
                              use1dim=True,
                              use1dimCor=True)
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
                              nround=50,
                              use1dim=True,
                              use1dimCor=True)
    unknown = ~(allVals.index.isin(knownPos.index) |
                allVals.index.isin(knownNeg.index))
    allVals[unknown] = fullfit[1][unknown]

    return allVals, castTblUnused, knownPos


def Probability_field_LDA_bootstraps(castTblUsed,
                                     castTblUnused,
                                     knownPos,
                                     knownNeg):
    castTblUsed = castTblUsed.reset_index()
    castTblUnused = castTblUnused.reset_index()
    knownPos = castTblUsed[
            castTblUsed.loc[:, "Patient_or_Control_ID"].isin(knownPos.index)]
    knownNeg = castTblUsed[
            castTblUsed.loc[:, "Patient_or_Control_ID"].isin(knownNeg.index)]

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
                              nround=50,
                              use1dim=True,
                              use1dimCor=True)
        myposPosterior = getPosteriorFromMVNFull(removedPos.loc[:, list_pEC50],
                                                 fit)
        mynegPosterior = getPosteriorFromMVNFull(removedNeg.loc[:, list_pEC50],
                                                 fit)
        allVals[allVals.index.isin(
                        knownPos.iloc[posIndexVec == myi].index)
                ] = myposPosterior
        allVals[allVals.index.isin(
                        knownNeg.iloc[negIndexVec == myi].index)
                ] = mynegPosterior

    fullfit = estimateLDAFrac(castTblUsed.loc[:, list_pEC50],
                              knownPos.loc[:, list_pEC50],
                              knownNeg.loc[:, list_pEC50],
                              nround=50,
                              use1dim=True,
                              use1dimCor=True)

    unknown = ~(allVals.index.isin(knownPos.index) |
                allVals.index.isin(knownNeg.index))
    unknown_index = allVals[~(allVals.index.isin(knownPos.index) |
                              allVals.index.isin(knownNeg.index))].index
    allVals[unknown_index] = fullfit[1][unknown]


    return allVals, castTblUsed, castTblUnused, knownPos


def mask_pcr_fct(df_2):
    # Not only having the covid before 14th day
    # but just PCR positive once before
    path_to_main_folder = os.path.dirname(
                                     os.path.abspath(''))
    path_annotation = glob(os.path.join(path_to_main_folder,
                           "**", "*ptn_adhoc_20200914_annot*.xlsx"),
                           recursive=True)
    df_annotation = pd.read_excel(path_annotation[0])
    df_annotation = df_annotation.set_index("j_number")
    df_annotation = df_annotation[~pd.isna(df_annotation.research_id)]
    df_annotation = df_annotation.loc[df_2.index]
    df_2.loc[:,
             'once_covid_positive'] = df_annotation.loc[:,
                                                        'once_covid_positive']
    df_2.loc[:,
             'datediff_days_lab_request_and_covid_result'] = df_annotation.loc[:, 'datediff_days_lab_request_and_covid_result']
    mask = df_2.loc[:, "once_covid_positive"] == "Ja"
    mask = mask & (df_2.loc[:, 'datediff_days_lab_request_and_covid_result']
                   <= 0)
    mask_pos_2 = (df_annotation["Clinically_manifest_COVID_"]
                  == "Yes") & (df_2.loc[:,
                                        'datediff_days_lab_request_and_covid_result']
                               < -14)
    mask = mask & np.invert(mask_pos_2)
    return df_2[mask].index

def calculate_prevalence_USZ_LDA_biweekly_bootstrap(df,
                                                    allVals,
                                                    castTblUnused,
                                                    knownPos,
                                                    knownNeg,
                                                    plot=False):
    allVals, castTblUsed, castTblUnused, knownPos = Probability_field_LDA_bootstraps(allVals, castTblUnused, knownPos, knownNeg)
    castTblUsed["posterior_Prob"] = allVals
    cast = castTblUsed.set_index("Patient_or_Control_ID")
    df.loc[cast.index, "LDA_prob"] = cast.loc[:, "posterior_Prob"]
    castTblUsed = pd.concat([castTblUsed, castTblUnused])
    castTblUsed = castTblUsed.reset_index()
    monthly = True
    if monthly:
        list_month_day = [(0, (0, np.inf), "2016-2018", "2016-2018"),
                          (12, (0, np.inf), "2019\nDec", "2019.12"),
                          (1, (0, np.inf), "2020\nJan", "2020.01"),
                          (2, (0, np.inf), "2020\nFeb", "2020.02"),
                          (3, (0, np.inf), "2020\nMar", "2020.03"),
                          (4, (0, np.inf), "2020\nApr", "2020.04"),
                          (5, (0, np.inf), "2020\nMay", "2020.05"),
                          (6, (0, np.inf), "2020\nJun", "2020.06"),
                          (7, (0, np.inf), "2020\nJul", "2020.07"),
                          (8, (0, np.inf), "2020\nAug", "2020.08"),
                          (9, (0, np.inf), "2020\nSep", "2020.09"),
                          (10, (0, np.inf), "2020\nOct", "2020.10"),
                          (11, (0, np.inf), "2020\nNov", "2020.11"),
                          (12, (0, np.inf), "2020\nDec", "2020.12"),]
    else:
        list_month_day = [(0, (0, np.inf), "2016-2018", "2016-2018"),
                          (12, (0, np.inf), "2019\nDec", "2019\nDec"),
                          (1, (0, np.inf), "2020\nJan", "2020.01-1st"),
                          (1, (16, np.inf), "2020\nJan", "2020.01-2nd"),
                          (2, (0, 16), "2020\nFeb", "2020.02-1st"),
                          (2, (16, np.inf), "2020\nFeb", "2020.02-2nd"),
                          (3, (0, 16), "2020\nMar", "2020.03-1st"),
                          (3, (16, np.inf), "2020\nMar", "2020.03-2nd"),
                          (4, (0, 16), "2020\nApr", "2020.04-1st"),
                          (4, (16, np.inf), "2020\nApr", "2020.04-2nd"),
                          (5, (0, 16), "2020\nMay", "2020.05-1st"),
                          (5, (16, np.inf), "2020\nMay", "2020.05-2nd"),
                          (6, (0, 16), "2020\nJun", "2020.06-1st"),
                          (6, (16, np.inf), "2020\nJun", "2020.06-2nd"),
                          (7, (0, 16), "2020\nJul", "2020.07-1st"),
                          (7, (16, np.inf), "2020\nJul", "2020.07-2nd"),
                          (8, (0, 16), "2020\nAug", "2020.07-1st"),
                          (8, (16, np.inf), "2020\nAug", "2020.07-2nd")]
    prevalence = np.zeros(len(list_month_day))
    prevalence_true = np.zeros(len(list_month_day))
    prevalence_pcr = np.zeros(len(list_month_day))
    list_xaxis = np.repeat("no_axis_name_found", len(list_month_day))
    for i, name in enumerate((list_month_day)):
        mask_name = castTblUsed["x axis violin name"] == name[2]
        mask_day = (castTblUsed["day"] >=
                    name[1][0]) & (castTblUsed["day"] < name[1][1])
        if name[0] == 0:
            mask_month = castTblUsed["month"] != "I want True everywhere"
        else:
            mask_month = castTblUsed["month"] == name[0]
        prevalence[i] = (np.nansum(castTblUsed[((mask_month &
                                                mask_day) &
                                               mask_name)].posterior_Prob) /
                         len(castTblUsed[((mask_month &
                                         mask_day) &
                                         mask_name)]))

        df_knownpos = castTblUsed.loc[((mask_month &
                                        mask_day) &
                                       mask_name)]
        df_pos_detected = df_knownpos[df_knownpos.posterior_Prob >= 0.5]
        pos_count = knownPos.loc[knownPos.index.isin(df_pos_detected.index)]
        # pos_count = knownPos.loc[
        #     knownPos.index.isin(castTblUsed.loc[((mask_month &
        #                                           mask_day) &
        #                                          mask_name)].index)
        #                         ]
        prevalence_true[i] = (len(pos_count) /
                              (len(castTblUsed[((mask_month &
                                                 mask_day) &
                                                mask_name)
                                               ]))
                              )
        list_xaxis[i] = name[3]
        if plot:
            df_mask = castTblUsed.loc[mask_month & mask_day & mask_name]
            df_mask = df_mask.set_index("Patient_or_Control_ID")
            index_pcr = mask_pcr_fct(df_mask)
            prevalence_pcr[i] = ((np.nansum(df_mask.loc[index_pcr,
                                                        "posterior_Prob"]) +
                                 len(pos_count)) /
                                 (len(castTblUsed[((mask_month &
                                                  mask_day) &
                                                  mask_name)])))
        print((f"{np.nansum(castTblUsed[mask_month & mask_day & mask_name].posterior_Prob)} / "+
              f"{len(castTblUsed[mask_month & mask_day & mask_name])}"))

    if plot:
        plt.figure()
        plt.bar(list_xaxis, prevalence * 100, label="USZ unkown")
        plt.bar(list_xaxis, prevalence_pcr * 100, label="USZ pcr detected")
        plt.bar(list_xaxis, prevalence_true * 100, label="USZ known")
        plt.ylabel("prevalence")
    return (list_xaxis, prevalence)


def bootstrap_USZ_LDA_biweekly(df, df_USZ, n_iterations=1000):
    np.random.seed(seed=111)
    statistics = []
    subcohort = []

    for i in np.arange(0, n_iterations):  # bootstraps:
        castTblUsed, castTblUnused, knownPos, knownNeg = filter_cohort(df_USZ,
                                                                       False)
        mask_pos = castTblUsed["hasCovid"] == 1
        mask_neg = castTblUsed["hasCovid"] == -1
        mask_neutral = castTblUsed["hasCovid"] == 0
        sample_pos = resample(castTblUsed[mask_pos],
                              replace=True,
                              random_state=i)
        sample_neg = resample(castTblUsed[mask_neg],
                              replace=True,
                              random_state=i)
        sample_neutral = resample(castTblUsed[mask_neutral],
                                  replace=True,
                                  random_state=i)
        sample = pd.concat([sample_pos,
                            sample_neg,
                            sample_neutral])
        stat = calculate_prevalence_USZ_LDA_biweekly_bootstrap(df,
                                                               sample,
                                                               castTblUnused,
                                                               knownPos,
                                                               knownNeg)
        statistics.append(stat[1])
        subcohort.append(stat[0])
        print(f"bootstrap is at {i}")
        print((f"mean = {np.mean(np.asarray(statistics)[:,-1])}, " +
               "percentile 2.5-97.5: " +
               f"{np.percentile(np.asarray(statistics)[:,-1],2.5)}" +
               f"-{np.percentile(np.asarray(statistics)[:,-1],97.5)}"))
    return np.asarray(statistics)


def calculate_prevalence_BDS_LDA_biweekly_bootstrap(df,
                                                    allVals,
                                                    castTblUnused,
                                                    knownPos,
                                                    knownNeg,
                                                    plot=False,
                                                    df_blood=False):
    allVals, castTblUsed, castTblUnused, knownPos = Probability_field_LDA_bootstraps(allVals, castTblUnused, knownPos, knownNeg)
    castTblUsed["posterior_Prob"] = allVals
    cast = castTblUsed.set_index("Patient_or_Control_ID")
    df.loc[cast.index, "LDA_prob"] = cast.loc[:, "posterior_Prob"]
    castTblUsed = castTblUsed.reset_index()
    list_xaxis = ['blood donation\n 2015',
                  'BDS 2019\nDec',
                  'BDS 2020\nJan',
                  'BDS 2020\nFeb',
                  'BDS 2020\nMar',
                  'BDS 2020\nApr',
                  'BDS 2020\nMay',
                  'BDS 2020\nJun',
                  'BDS 2020\nJul',
                  'BDS 2020\nAug',
                  'BDS 2020\nSep',
                  'blood donation\n positif']
    prevalence = np.zeros(len(list_xaxis))
    prevalence_true = np.zeros(len(list_xaxis))
    for i, name in enumerate((list_xaxis)):
        mask_name = castTblUsed["x axis violin name"] == name
        mask_name_unused = castTblUnused["x axis violin name"] == name
        if len(castTblUsed[mask_name]) > 0:
            prevalence[i] = (np.nansum(castTblUsed[mask_name].posterior_Prob) /
                             (len(castTblUsed[mask_name]) +
                              len(castTblUnused[mask_name_unused])))
            pos_count = knownPos.loc[knownPos.index.isin(
                                        castTblUsed.loc[mask_name].index)]
            prevalence_true[i] = len(pos_count) / (len(castTblUsed[mask_name]))
        else:
            prevalence[i] = np.nan
            prevalence_true[i] = np.nan
        print(f"BDS {name} :{np.nansum(castTblUsed[mask_name].posterior_Prob)} / {len(castTblUsed[mask_name])}")
    if plot:
        plt.bar(list_xaxis, prevalence * 100, label="USZ BDS")
        plt.ylabel("prevalence")
    return (list_xaxis, prevalence)


def bootstrap_BDS_LDA_biweekly(df, df_blood, n_iterations = 1000):
    np.random.seed(seed=111)
    statistics = []
    subcohort = []

    for i in np.arange(0, n_iterations):  # bootstraps:
        castTblUsed, castTblUnused, knownPos, knownNeg = filter_cohort(df_blood,
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
        stat = calculate_prevalence_BDS_LDA_biweekly_bootstrap(df,
                                                               sample,
                                                               castTblUnused,
                                                               knownPos,
                                                               knownNeg)
        statistics.append(stat[1])
        subcohort.append(stat[0])
    return np.asarray(statistics)


def calculate_prevalence_BDS_biweekly_bootstrap(df,
                                                allVals,
                                                castTblUnused,
                                                knownPos,
                                                knownNeg):
    allVals, castTblUsed, castTblUnused, knownPos= Probability_field_FMW_bootstraps(allVals, castTblUnused, knownPos, knownNeg)
    castTblUsed["posterior_Prob"] = allVals
    castTblUsed = pd.concat([castTblUsed, castTblUnused])
    castTblUsed = castTblUsed.reset_index()
    list_xaxis = ['blood donation\n 2015',
                  'BDS 2019\nDec',
                  'BDS 2020\nJan',
                  'BDS 2020\nFeb',
                  'BDS 2020\nMar',
                  'BDS 2020\nApr',
                  'BDS 2020\nMay',
                  'BDS 2020\nJun',
                  'BDS 2020\nJul',
                  'BDS 2020\nAug',
                  'BDS 2020\nSep',
                  'blood donation\n positif']
    prevalence = np.zeros(len(list_xaxis))
    prevalence_true = np.zeros(len(list_xaxis))
    for i, name in enumerate((list_xaxis)):
        mask_name = castTblUsed["x axis violin name"] == name
        if len(castTblUsed[mask_name]) > 0:
            prevalence[i] = (np.nansum(castTblUsed[mask_name].posterior_Prob) /
                             len(castTblUsed[mask_name]))
            pos_count = knownPos.loc[knownPos.index.isin(castTblUsed.loc[mask_name].index)]
            prevalence_true[i] = len(pos_count) / (len(castTblUsed[mask_name]))
        else:
            prevalence[i] = np.nan
            prevalence_true[i] = np.nan
    return (list_xaxis, prevalence)


def bootstrap_BDS_biweekly(df, df_blood, n_iterations = 1000):
    np.random.seed(seed=111)
    statistics = []
    subcohort = []

    for i in np.arange(0,n_iterations):  # bootstraps:
        castTblUsed, castTblUnused, knownPos, knownNeg = filter_cohort(df_blood,
                                                                       True)
        mask_pos = castTblUsed["hasCovid"] == 1
        mask_neg = castTblUsed["hasCovid"] == -1
        mask_neutral = castTblUsed["hasCovid"] == 0
        sample_pos = resample(castTblUsed[mask_pos],
                              replace=True,
                              random_state=i)
        sample_neg = resample(castTblUsed[mask_neg],
                              replace=True,
                              random_state=i)
        sample_neutral = resample(castTblUsed[mask_neutral],
                                  replace=True,
                                  random_state=i)
        sample = pd.concat([sample_pos, sample_neg, sample_neutral])
        stat = calculate_prevalence_BDS_biweekly_bootstrap(df,
                                                           sample,
                                                           castTblUnused,
                                                           knownPos,
                                                           knownNeg)
        statistics.append(stat[1])
        subcohort.append(stat[0])
    return np.asarray(statistics)