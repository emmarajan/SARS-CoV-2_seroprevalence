# -*- coding: utf-8 -*-
"""
Created on Sun May  3 20:44:40 2020

@author: rapha
"""
import pandas as pd
import numpy as np
import re
from glob import glob
import os
from gaussian_estimation import Probability_field_FMW
import matplotlib.pyplot as plt

path_to_main_folder = os.path.dirname(os.path.abspath(''))
paths_REKO = glob("*REKO.xlsx")
reko = pd.read_excel(paths_REKO[0])
reko1 = [rek[:13] for rek in reko.Patient_or_Control_ID]
reko["Patient_or_Control_ID"] = reko1
paths_oxford_filter = glob(os.path.join(path_to_main_folder,
                                        "**",
                                        "*Positive_Negative Plate Annotations*.xlsx"),
                           recursive=True)
oxford_df = pd.read_excel(paths_oxford_filter[0])


paths_truepositive = glob((os.path.join(path_to_main_folder,
                                        "**",
                                        "*COVID19_samples_Dropbox.csv")),
                          recursive=True)
df_list_TP = pd.read_csv(paths_truepositive[0])
df_list_TP = df_list_TP[np.isnan(df_list_TP.Age) == False]
df_list_TP.rename(columns={'Pseudo_ID': 'Patient_or_Control_ID'},
                  inplace=True)

IgGnames = ['NC_IgG', 'Spike_IgG', 'RBD_IgG']
sample_dic = {"Patient": ("[{J}}]{1}[\d]{1,}",
                          "[{H}}]{1}[\d]{1,}",
                          "[{C}}]{1}[\d]{1,}",
                          "Negative[\d]{1,}",
                          "COVID19_[\d]{1,}",
                          "BDS[\d]{1,}",
                          "^[\d]{1,}"),
              "Negative_serum": "neg",
              "Positive_serum": "Positive_serum",
              "colour": "colour[\d]{1,}"}


def specific_data_load(path_csvfile, sample_value):
    """
    Load a dataframe from fits_out.csv, from a specific type of sample put
    inside a dictionary
    Patient, Negative_serum, Positive_serum, colour
    example: specific_data_load("*fits_out.csv", ("Patient",))
             specific_data_load("*fits_out.csv", ("Patient","Negative_serum"))

    Parameters
    ----------
    path_csvfile : String finishing by fits_out.csv
        serve to load csv with the following column
        Columns: [Patient_or_Control_ID, Plate_ID, target, Ab_conc, KD,
        IC50, Max_fittable_IC50, Error]
    sample_value : list of str of needed value to load using the following dic
        sample_dic = {"Patient":"J[\d]{1,}",
              "Negative_serum":"neg",
              "Positive_serum":"Positive_serum",
              "colour":"colour[\d]{1,}"}

    Returns
    -------
    entire_data : dataframe of sample value wanted
        dataframe with the following info:
        Columns: [Patient_or_Control_ID, Plate_ID, target, Ab_conc, KD,
        IC50, Max_fittable_IC50, Error]

    """
    data = pd.read_csv(path_csvfile, sep="\t")
    #Problem is depending of type of data you need to  use sep="\t" or not
    #this is why the if condition is used, as well as 2 different name:
    #'Patient_or_Control_ID' or "patient_or_control_id" for same effect in df
    if 'Patient_or_Control_ID' in data:
        for name in sample_value:
            entire_data = pd.DataFrame()
            if type(sample_dic[name]) is tuple:
                for name_sample in sample_dic[name]:
                    mask = data["Patient_or_Control_ID"].str.count(name_sample) == 1
                    load = data[mask]
                    entire_data = entire_data.append(load)
            else:
                mask = data["Patient_or_Control_ID"].str.count(sample_dic[name]) == 1
                load = data[mask]
                entire_data = entire_data.append(load)
        return entire_data

    else:
        data = pd.read_csv(path_csvfile)
        entire_data = pd.DataFrame()
        for name in sample_value:
            mask = data["patient_or_control_id"].str.count(sample_dic[name]) == 1
            load = data[mask]
            entire_data = entire_data.append(load)
        return entire_data


def get_date_from_id(serum_id):
    """
    Return the date of an identifyer of serum

    Parameters
    ----------
    serum_id : string
        following the code JAMMJJ.

    Returns
    -------
    year : string
    month : string
    day : string

    """
    if serum_id[0] == " ":
        serum_id = serum_id[1:]
    regex = r"[(\w{J,H})]?(\d{1})?(\d{2})?(\d{2})?"
    matches = re.match(regex, serum_id)
    if matches[0][0] == 'J':
        cohort = "hospital patient"
        if matches[1] == '0':
            year = "202" + matches[1]
        else:
            year = "201" + matches[1]
        month = matches[2]
        day = matches[3]
        if int(year) < 2019:
            cohort = "hospital patient < 2019"
    elif matches[0][0] == 'H':
        if serum_id[:13] in reko1:
            cohort = "serum donation positive"
            year = "3020"
            month = "13"
            day = "0"
        else:
            # year, month, day = extract_date_Hnumber(serum_id)
            cohort = "blood donation"
            year = "3020"
            month = "04"
            day = "0"
    elif re.match("COVID19", serum_id):
        if serum_id in df_list_TP["Patient_or_Control_ID"].values:
            cohort = "KIM_pos"
            year = "6020"
            month = "00"
            day = "0"
        else:
            cohort = "KIM_neg"
            year = "6020"
            month = "00"
            day = "0"
    elif re.match("[{C}}]{1}[\d]{1,}", serum_id, flags=0):
        if serum_id in oxford_df["barcode"].values:
            cohort = ("Oxford " +
                      oxford_df[
                          oxford_df["barcode"] == serum_id
                              ]["sample_cat"].values[0]
                      )
            year = "4020"
            month = "00"
            day = "0"
        else:
            print(serum_id)
            cohort = "Oxford Nomatch"
            year = "4020"
            month = "00"
            day = "0"
    elif re.match("Negative[\d]{1,}", serum_id, flags=0):
        cohort = "Oxford_neg"
        year = "4020"
        month = "13"
        day = "0"
    elif re.match("BDS", serum_id):
        matches = re.match("BDS(\d{4})?(\d{2})?", serum_id, flags=0)
        if int(matches[1]) <= 2018:
            cohort = "blood donation negative"
            year = str(int(matches[1]) + 1000)
            month = matches[2]
            day = "0"
        else:
            cohort = "blood donation"
            year = str(int(matches[1]) + 1000)
            month = matches[2]
            day = "0"
    elif re.match("\d{1,}", serum_id):
        cohort = "Only number"
        year = str(0)
        month = str(0)
        day = "0"
    else:
        raise RuntimeError("Name of serum does not start " +
                           f"with J or H: \n {serum_id}")
    return (year, month, day, cohort)


def filter_for_renameing_plot(df):
    """
    insert a name as string for different cohorte easier to
    plot the data more efficiently.

    Parameters
    ----------
    df : Panda dataframe, containing year month and cohort.

    Returns
    -------
    dataframe Panda with new axis name.

    """
    dict_month = {0: "NA", 1: "Jan", 2: "Feb", 3: "Mar", 4: "Apr",
                  5: "May",6: "Jun", 7: "Jul", 8: "Aug", 9: "Sep",
                  10: "Oct", 11: "Nov", 12: "Dec", 13: "NA"}
    axisname = np.zeros(df.shape[0])
    axisname = np.asarray(axisname, dtype=str)
    axisname[axisname == "0.0"] = ("No cathegorie yet \n please" +
                                   " look at filter_for_renaming_plot")

    axisname[((df.year <= 2018)
             & (df.cohort == "hospital patient < 2019"))] = "2016-2018"
    axisname[df.cohort == "blood donation negative"] = "blood donation\n 2015"
    axisname[df.cohort == "blood donation"] = "blood donation\n 2020"
    for mm in range(1, 13):
        axisname[(df.year == 2019) & (df.cohort == "hospital patient") &
                 (df.month == mm)] = ("2019\n" + dict_month[mm])
        axisname[(df.year == 2020) & (df.cohort == "hospital patient") &
                 (df.month == mm)] = ("2020\n" + dict_month[mm])
        axisname[(df.year == 3019) & (df.cohort == "blood donation") &
                 (df.month == mm)] = ("BDS 2019\n" + dict_month[mm])
        axisname[(df.year == 3020) & (df.cohort == "blood donation") &
                 (df.month == mm)] = ("BDS 2020\n" + dict_month[mm])

    axisname[df.cohort == "KIM_pos"] = "KIM_pos"
    axisname[df.cohort == "KIM_neg"] = "KIM_neg"

    axisname[df.cohort == "Oxford_neg"] = "Oxford_neg"
    axisname[df.cohort == "Oxford known negative"] = "Oxford_neg"
    axisname[df.cohort == "Oxford acute positive"] = "Oxford_pos"
    axisname[df.cohort == "Oxford convalescent positive"] = "Oxford_pos"
    axisname[df.cohort == "Oxford low positive serial dilution"] = "Oxford_pos"
    axisname[df.cohort == "Oxford other coronavirus acute"] = "Oxford"
    axisname[df.cohort == "Oxford other coronavirus conv"] = "Oxford"
    axisname[df.cohort == "Oxford high positive serial dilution"] = "Oxford_pos"
    axisname[df.cohort == "Oxford other resp acute"] = "Oxford"
    axisname[df.cohort == "Oxford medium positive serial dilution"] = "Oxford_pos"

    axisname[df.cohort == "serum donation positive"] = "blood donation\n positif"
    df.insert(1, "x axis violin name", axisname)
    return df

def method_use(df, method, arg=None):
    if method == "average_SpikeRBD":
        df.loc[:, "methodpEC50"] = df.loc[:,
                                          ["Spike_pEC50", "RBD_pEC50"]
                                          ].mean(axis=1)
    elif method == "min_SpikeRBD":
        df.loc[:, "methodpEC50"] = df.loc[:,
                                          ["Spike_pEC50", "RBD_pEC50"]
                                          ].min(axis=1)
    elif method == "NC":
        df.loc[:, "methodpEC50"] = df.NC_pEC50
    elif method == "Spike":
        df.loc[:, "methodpEC50"] = df.Spike_pEC50
    elif method == "RBD":
        df.loc[:, "methodpEC50"] = df.RBD_pEC50
    elif method == "posterior":
        mask_blood = df["cohort"] == "blood donation"
    else:
        raise RuntimeError("Method does not exists")
    return df


def find_sourcewell(df, name, methods):
    path_to_main_folder = os.path.dirname(os.path.abspath(''))
    paths_csvfilesdata = glob(os.path.join(path_to_main_folder,
                                           "**",
                                           "*hts*.csv"),
                              recursive=True)
    df_data = pd.DataFrame()
    for path_csvfile in paths_csvfilesdata:
        data = pd.read_csv(path_csvfile)
        df_data = df_data.append(data)
    df_data = df_data.rename(columns={"patient_or_control_id": "Patient_or_Control_ID"})
    df_data = df_data.set_index("Patient_or_Control_ID")
    df3 = pd.merge(df,
                   df_data,how = "inner",
                   left_index=True,
                   right_index=True)
    df3 = df3.loc[~df3.index.duplicated(keep='first')]
    df3.loc[:, "methods"] = np.full(len(df3), methods)
    df3.to_csv(name, index=True)
    return df3


def filter_oxford_data(df):
    path_to_main_folder = os.path.dirname(os.path.abspath(''))
    paths_csvfilesdata = glob(os.path.join(path_to_main_folder,
                                           "**",
                                           "*Positive_Negative Plate Annotations (1).xlsx"),
                              recursive=True)
    csv_to_keep = pd.read_excel(paths_csvfilesdata[0])
    mask_oxford = df["cohort"] == "oxford"
    df_oxford = df.loc[mask_oxford]
    x = pd.merge(csv_to_keep, df_oxford, left_on="barcode", right_index=True)
    x.to_excel("oxford.xlsx")


def Age_gender_usz(df, allVals):
    paths_KIM = glob((os.path.join(path_to_main_folder,
                                   "**",
                                   "*COVID19_samples_Dropbox*.csv")),
                     recursive=True)
    df_KIM = pd.read_csv(paths_KIM[0])
    paths_description = glob((os.path.join(path_to_main_folder,
                                           "**",
                                           "*20200524_descrp_stats*.xlsx")),
                             recursive=True)
    paths_clinique = glob(os.path.join(path_to_main_folder,
                                       "**",
                                       "*Organisationseinheit*.xlsx"),
                          recursive=True)
    paths_selected_person = glob(os.path.join(path_to_main_folder,
                                              "**",
                                              "*20200530_descrp_stats_full_list_2*.xlsx"),
                                 recursive=True)
    df_selected = pd.read_excel(paths_selected_person[0])
    df_selected = df_selected.rename(columns={"PatientIDList": "Patient_or_Control_ID"})
    df_selected = df_selected[df_selected.Patient_or_Control_ID.isin(df.index)]
    list_patient = df_selected.Patient_or_Control_ID.drop_duplicates()
    index = np.zeros(len(list_patient))
    for i, pat_id in enumerate(list_patient):
        mask_patient = df_selected.Patient_or_Control_ID == pat_id
        list_date = df_selected.loc[mask_patient, "Eindat"]
        if np.int(pat_id[1]) == 0:
            minimum_store = +np.inf
            for date in list_date:
                minimum = (np.abs(2020 - np.int(date[0:4])) * 365 +
                           np.abs(np.int(pat_id[2:4]) - np.int(date[5:7])) *
                           31 +
                           np.abs(np.int(pat_id[4:6]) - np.int(date[8:10])))
                if minimum < minimum_store:
                    print(f"list patient : {i}: {minimum:3.2f} {minimum_store:3.2f}")
                    minimum_store = minimum
                    mask_index = date == df_selected.loc[:, "Eindat"]
                    index[i] = df_selected.loc[mask_patient & mask_index, "Eindat"].index[0]
        if np.int(pat_id[1]) != 0:
            minimum_store = +np.inf
            for date in list_date:
                minimum = (np.abs(2010 + np.int(pat_id[1]) -
                                  np.int(date[0:4])) * 365 +
                           np.abs(np.int(pat_id[2:4]) -
                                  np.int(date[5:7])) * 31 +
                           np.abs(np.int(pat_id[4:6]) - np.int(date[8:10])))
                if minimum < minimum_store:
                    print(f"list patient : {i}: {minimum:3.2f} {minimum_store:3.2f}")
                    minimum_store = minimum
                    mask_index = date == df_selected.loc[:, "Eindat"]
                    index[i] = df_selected.loc[mask_patient & mask_index,
                                               "Eindat"].index[0]

    usz_description = pd.read_excel(paths_description[0])
    df_clinique_number = pd.read_excel(paths_clinique[0])
    # Comparison
    df = df[~np.isnan(df["RBD_pEC50"])]
    np.sum(df.index.isin(usz_description["PatientIDList"]))
    usz_description = usz_description.rename(columns={"PatientIDList":
                                                      "Patient_or_Control_ID"})
    usz_description = usz_description.set_index("Patient_or_Control_ID")
    usz_description = usz_description[usz_description.index.isin(df.index)]
    usz_description["Age"] = 2020 - usz_description.loc[:, "YearBirth"]
    mask_Men = usz_description.loc[:, "Gender"] == "M"
    mask_Woman = usz_description.loc[:, "Gender"] == "W"
    span = np.arange(np.min(usz_description.loc[:, "Age"]),
            np.max(usz_description.loc[:, "Age"]), 2)
    x_M, y_M = np.histogram(usz_description.loc[mask_Men, "Age"], span)
    x_W, y_W = np.histogram(usz_description.loc[mask_Woman, "Age"], span)
    df_pyramid = pd.DataFrame({'Age': span[:-1],
                               'Male': -x_M,
                               'Female': x_W})
    plt.bar(df_pyramid["Age"],
            df_pyramid["Male"],
            label="Male",
            width=1.8)
    plt.bar(df_pyramid["Age"],
            df_pyramid["Female"],
            label="Female",
            width=1.8)
    plt.legend()
    plt.xlim([10, 100])
    plt.ylim([-300, 300])
    plt.xlabel("Age")
    plt.ylabel("population")

    # repartission par clinique
    plt.figure()
    index = np.loadtxt("index_big.csv")
    index_2 = index > 0
    index_2 = index[index_2]
    index_2 = np.asarray(index_2, dtype="int")
    clinique_occurance = df_selected.iloc[index_2]
    clinique_occurance = clinique_occurance["Orgfa"].value_counts()
    # clinique_occurance = usz_description["Orgfa"].value_counts()
    clinique_occurance = pd.DataFrame(clinique_occurance)
    clinique_occurance = clinique_occurance.rename(columns={"Orgfa": "Counts"})
    clinique_occurance["Orgfa"] = clinique_occurance.index
    clinique_occurance = pd.merge(clinique_occurance,
                                  df_clinique_number,
                                  on="Orgfa")
    plt.bar(clinique_occurance["clinique"],
            clinique_occurance["Counts"],
            bottom=0.1)
    plt.xticks(rotation=90)

    # Add infected person
    mask_usz = df["cohort"] == "hospital patient"
    mask_KIM = (df["subcohort"] == "KIM_neg") | (df["subcohort"] == "KIM_pos")
    allVals, castTblUnused, knownPos = Probability_field_FMW(df.loc[mask_usz
                                                                    & ~mask_KIM])
    castTblUnused.loc[:, "Spike_pEC50"] = 0.001
    allVals = pd.concat([allVals, castTblUnused.Spike_pEC50])
    df_experiment = pd.merge(df, allVals, on="Patient_or_Control_ID")
    df_experiment.Spike_pEC50_y
    x = pd.merge(df_selected.iloc[index_2], df_experiment,
                 on="Patient_or_Control_ID")

    mask = (x.Spike_pEC50_y > 0.5)
    x = x.loc[mask]
    clinique_occurance_infected = x["Orgfa"].value_counts()
    clinique_occurance_infected = pd.DataFrame(clinique_occurance_infected)
    clinique_occurance_infected = clinique_occurance_infected.rename(columns={"Orgfa":
                                                                              "Counts"})
    clinique_occurance_infected["Orgfa"] = clinique_occurance_infected.index
    clinique_occurance_infected = pd.merge(clinique_occurance_infected,
                                           df_clinique_number,
                                           on="Orgfa")
    plt.bar(clinique_occurance_infected["clinique"],
            clinique_occurance_infected["Counts"],
            bottom=0.1)
    plt.xticks(rotation=90)

# %%
    plt.figure()
    x = pd.merge(usz_description, df_experiment, on="Patient_or_Control_ID")
    mask_M = x.loc[:, "Gender"] == "M"
    mask_W = x.loc[:, "Gender"] == "W"
    mask_P = (x.Spike_pEC50_y > 0.5)
    mask_noP = (x.Spike_pEC50_y <= 0.5)
    plt.scatter(x.loc[mask_M, "Age"],
                x.loc[mask_M, "Spike_pEC50_y"],
                marker="X",
                label="Male",
                color=np.asarray([31, 119, 180]) / 255.)
    plt.scatter(x.loc[mask_W, "Age"],
                x.loc[mask_W, "Spike_pEC50_y"],
                marker="P",
                label="Female",
                color=np.asarray([14, 176, 255]) / 255.)
    plt.legend()
    plt.xlabel("Age [year]")
    plt.ylabel("Posterior probability associated")

    plt.figure()
    span = np.arange(np.min(x.loc[:, "Age"]),
            np.max(x.loc[:, "Age"]), 2)
    x_M, y_M = np.histogram(x.loc[mask_M & mask_noP, "Age"], span)
    x_W, y_W = np.histogram(x.loc[mask_W & mask_noP, "Age"], span)
    df_pyramid = pd.DataFrame({'Age': span[:-1],
                               'Male': -x_M,
                               'Female': x_W})
    plt.bar(df_pyramid["Age"],
            df_pyramid["Male"] / 40,
            label="Male count/40",
            width=2)
    plt.bar(df_pyramid["Age"],
            df_pyramid["Female"] / 40,
            label="Female count/40",
            width=2)

    x_M, y_M = np.histogram(x.loc[mask_M & mask_P, "Age"], span)
    x_W, y_W = np.histogram(x.loc[mask_W & mask_P, "Age"], span)
    df_pyramid = pd.DataFrame({'Age': span[:-1],
                               'Male': -x_M,
                               'Female': x_W})
    plt.bar(df_pyramid["Age"],
            df_pyramid["Male"], label="Male P>0.5", width=1.8)
    plt.bar(df_pyramid["Age"],
            df_pyramid["Female"], label="Female P>0.5", width=1.8)
    plt.legend()
    plt.xlabel("Age")
    plt.ylabel("population")
# %%
    mask_M = x.loc[:, "Gender"] == "M"
    mask_W = x.loc[:, "Gender"] == "W"
    mask_P = (x.Spike_pEC50_y > 0.5)
    mask_noP = np.isnan(x.Spike_pEC50_y)
    mask_noP = (x.Spike_pEC50_y <= 0.5) | mask_noP
    print(f"M P>0.5: number: {np.sum(mask_M & mask_P)}")
    print(f"W P>0.5: number: {np.sum(mask_W & mask_P)}")
    print(f"tot P>0.5: number: {np.sum(mask_P)}")
    print(f"M P<=0.5: number: {np.sum(mask_M & mask_noP)}")
    print(f"W P<=0.5: number: {np.sum(mask_W & mask_noP)}")
    print(f"tot P<=0.5: number: {np.sum(mask_noP)}")
    print(f"M : number: {np.sum(mask_M)}")
    print(f"W : number: {np.sum(mask_W)}")
    print(f"tot : number: {np.sum(mask_W | mask_M)}")
    print(f"M : Median age :{np.median(x.Age.loc[mask_M]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_M],25):2.0f}-{np.percentile(x.Age.loc[mask_M],75):2.0f}")
    print(f"M(P>0.5) : Median age :{np.median(x.Age.loc[mask_M & mask_P]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_M & mask_P],25):2.0f}-{np.percentile(x.Age.loc[mask_M & mask_P],75):2.0f}")
    print(f"M(P<=0.5) : Median age :{np.median(x.Age.loc[mask_M & mask_noP]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_M & mask_noP],25):2.0f}-{np.percentile(x.Age.loc[mask_M & mask_noP],75):2.0f}")
    print(f"W : Median age :{np.median(x.Age.loc[mask_W]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_W],25):2.0f}-{np.percentile(x.Age.loc[mask_W],75):2.0f}")
    print(f"W(P>0.5) : Median age :{np.median(x.Age.loc[mask_W & mask_P]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_W & mask_P],25):2.0f}-{np.percentile(x.Age.loc[mask_W & mask_P],75):2.0f}")
    print(f"W(P<=0.5) : Median age :{np.median(x.Age.loc[mask_W & mask_noP]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_W & mask_noP],25):2.0f}-{np.percentile(x.Age.loc[mask_W & mask_noP],75):2.0f}")
    print(f"tot : Median age :{np.median(x.Age):2.0f}, IQR:{np.percentile(x.Age,25):2.0f}-{np.percentile(x.Age,75):2.0f}")
    print(f"tot(P>0.5) : Median age :{np.median(x.Age.loc[mask_P]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_P],25):2.0f}-{np.percentile(x.Age.loc[mask_P],75):2.0f}")
    print(f"tot(P<=0.5) : Median age :{np.median(x.Age.loc[mask_noP]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_noP],25):2.0f}-{np.percentile(x.Age.loc[mask_noP],75):2.0f}")
#%%
    plt.figure()
    x = pd.merge(usz_description, df_experiment, on="Patient_or_Control_ID")
    x = pd.merge(x, df_clinique_number, on="Orgfa")
    mask_M = x.loc[:, "Gender"] == "M"
    mask_W = x.loc[:, "Gender"] == "W"
    plt.scatter(x.loc[mask_M, "clinique"],
                x.loc[mask_M, "Spike_pEC50_y"],
                marker="X",
                label="Male",
                color=np.asarray([31, 119, 180]) / 255.)
    plt.scatter(x.loc[mask_W, "clinique"],
                x.loc[mask_W, "Spike_pEC50_y"],
                marker="P",
                label="Female",
                color=np.asarray([14, 176, 255]) / 255.)
    plt.legend()
    plt.ylabel("Posterior probability associated")
    plt.xticks(rotation=90)

    nbr_patient = usz_description["HashID_KSIM"].value_counts()
    nbr_patient = pd.DataFrame(nbr_patient)
    nbr_patient = nbr_patient["HashID_KSIM"].value_counts()
    print("Number of time a patient give the analysed blood")
    for i, patient_nbre in enumerate(nbr_patient):
        print(f"{patient_nbre} patients give blood {nbr_patient.index[i]}x ")

# %%
    x = pd.merge(nbr_patient,
                 df_experiment,
                 n="Patient_or_Control_ID")
    mask_duplicate = x["HashID_KSIM"].drop_duplicates()
    x = x.loc[mask_duplicate.index]
    mask_M = x.loc[:, "Gender"] == "M"
    mask_W = x.loc[:, "Gender"] == "W"
    mask_P = (x.Spike_pEC50_y > 0.5)
    mask_noP = np.isnan(x.Spike_pEC50_y)
    mask_noP = (x.Spike_pEC50_y <= 0.5) | mask_noP
    print("NUMBER OF PATIENT!")
    print(f"M P>0.5: number: {np.sum(mask_M & mask_P)}")
    print(f"W P>0.5: number: {np.sum(mask_W & mask_P)}")
    print(f"tot P>0.5: number: {np.sum(mask_P)}")
    print(f"M P<=0.5: number: {np.sum(mask_M & mask_noP)}")
    print(f"W P<=0.5: number: {np.sum(mask_W & mask_noP)}")
    print(f"tot P<=0.5: number: {np.sum(mask_noP)}")
    print(f"M : number: {np.sum(mask_M)}")
    print(f"W : number: {np.sum(mask_W)}")
    print(f"tot : number: {np.sum(mask_W | mask_M)}")
    print(f"M : Median age :{np.median(x.Age.loc[mask_M]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_M],25):2.0f}-{np.percentile(x.Age.loc[mask_M],75):2.0f}")
    print(f"M(P>0.5) : Median age :{np.median(x.Age.loc[mask_M & mask_P]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_M & mask_P],25):2.0f}-{np.percentile(x.Age.loc[mask_M & mask_P],75):2.0f}")
    print(f"M(P<=0.5) : Median age :{np.median(x.Age.loc[mask_M & mask_noP]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_M & mask_noP],25):2.0f}-{np.percentile(x.Age.loc[mask_M & mask_noP],75):2.0f}")
    print(f"W : Median age :{np.median(x.Age.loc[mask_W]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_W],25):2.0f}-{np.percentile(x.Age.loc[mask_W],75):2.0f}")
    print(f"W(P>0.5) : Median age :{np.median(x.Age.loc[mask_W & mask_P]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_W & mask_P],25):2.0f}-{np.percentile(x.Age.loc[mask_W & mask_P],75):2.0f}")
    print(f"W(P<=0.5) : Median age :{np.median(x.Age.loc[mask_W & mask_noP]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_W & mask_noP],25):2.0f}-{np.percentile(x.Age.loc[mask_W & mask_noP],75):2.0f}")
    print(f"tot : Median age :{np.median(x.Age):2.0f}, IQR:{np.percentile(x.Age,25):2.0f}-{np.percentile(x.Age,75):2.0f}")
    print(f"tot(P>0.5) : Median age :{np.median(x.Age.loc[mask_P]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_P],25):2.0f}-{np.percentile(x.Age.loc[mask_P],75):2.0f}")
    print(f"tot(P<=0.5) : Median age :{np.median(x.Age.loc[mask_noP]):2.0f}, IQR:{np.percentile(x.Age.loc[mask_noP],25):2.0f}-{np.percentile(x.Age.loc[mask_noP],75):2.0f}")

