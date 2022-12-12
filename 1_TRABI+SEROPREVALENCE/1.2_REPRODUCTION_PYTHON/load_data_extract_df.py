# -*- coding: utf-8 -*-
"""
Created on Fri May 15 08:26:59 2020

@author: rapha
"""
from glob import glob
import pandas as pd
import numpy as np
import os
import re
from filtering_serumid_script import *


def load_all_data(path_folder=False,
                  QC_fname="*Overview_data_transfer*.csv"):
    normalisation = False
    IgGnames = ['NC_IgG', 'Spike_IgG', 'RBD_IgG']
    sample_dic = {"Patient": ("[{J}}]{1}[\d]{1,}",
                             "[{H}}]{1}[\d]{1,}",
                             "COVID19_[\d]{1,}"),
                  "Negative_serum": "neg",
                  "Positive_serum": "Positive_serum",
                  "colour": "colour[\d]{1,}"}
    #  Code source on a folder separate from one to the entire data
    #  access to the similar commun folder with the data:
    if path_folder is False:
        path_to_main_folder = os.path.dirname(
                                         os.path.abspath(''))
    else:
        path_to_main_folder = path_folder

    #  look at all the csv files present in the subfolder
    paths_csvfiles = glob((os.path.join(path_to_main_folder,
                                        "**",
                                        "*fits_out.csv")),
                          recursive=True)
    # paths_csvfiles = [i for i in paths_csvfiles if '20200519' not in i]
    #
    paths_csvfiles = [i for i in paths_csvfiles if "IgA" not in i]
    Not_selected_SARS = ['NC_SARS_IgG', 'Spike_SARS_IgG', 'RBD_SARS_IgG',
                         "HKU1_IgG", "OC43_IgG"]
    for name in Not_selected_SARS:
        paths_csvfiles = [i for i in paths_csvfiles if name not in i]
    if QC_fname is not False:
        paths_QC = glob((os.path.join(path_to_main_folder,
                                      "**",
                                      QC_fname)),
                        recursive=True)
        QC = pd.read_csv(paths_QC[0], sep=";")
        QC = QC[QC["QCUser"] == 0]

    struct_df = [pd.DataFrame(), pd.DataFrame(), pd.DataFrame()]
    struct_df_pos = [pd.DataFrame(), pd.DataFrame(), pd.DataFrame()]
    for index_igG, IgGname in enumerate(IgGnames):
        total_db_serum = pd.DataFrame()
        total_positive = pd.DataFrame()
        for paths_csvfile in paths_csvfiles:
            matchObj = re.search(IgGname, paths_csvfile, flags=0)
            if matchObj:
                data = pd.read_csv(paths_csvfile, sep="\t")
                serum = specific_data_load(paths_csvfile, ("Patient",))
                dates = [get_date_from_id(name)
                         for name in serum["Patient_or_Control_ID"]]
                dates = np.asarray(dates, dtype=str)

                serum.insert(4, "cohort", dates[:, 3])
                date = np.asarray(dates[:, 0:3], dtype=int)
                serum.insert(1, "year", date[:, 0])
                serum.insert(2, "month", date[:, 1])
                serum.insert(3, "day", date[:, 2])
                total_db_serum = total_db_serum.append(serum)

                positive = specific_data_load(paths_csvfile,
                                              ("Positive_serum",))

                positive.insert(4, "cohort", dates[0:positive.shape[0], 3])
                date = np.asarray(date[:, 0:3], dtype=int)
                date = np.mean(date, axis=0)
                date = np.ones((positive.shape[0], 3)) * date
                positive.insert(1, "year", date[:, 0])
                positive.insert(2, "month", date[:, 1])
                positive.insert(3, "day", date[:, 2])
                positive.insert(1, "IgG_use", IgGname)
                total_positive = total_positive.append(positive)

        if QC_fname is not False:
            total_db_serum = total_db_serum[~total_db_serum["Plate_ID"].isin(QC["Destination plate"])]
            total_positive = total_positive[~total_positive["Plate_ID"].isin(QC["Destination plate"])]
        total_db_serum = filter_for_renameing_plot(total_db_serum)
        total_db_serum = total_db_serum.sort_values(by=["year", "month"],
                                                    ascending=(True, True))
        total_db_serum = total_db_serum.sort_values("Plate_ID")
        total_db_serum = total_db_serum.reset_index(drop=True)
        noduplicate = total_db_serum["Patient_or_Control_ID"].drop_duplicates(keep='first')

        idx_noduplicate = noduplicate.index
        total_db_serum = total_db_serum.iloc[idx_noduplicate]
        tricky_violin_plototal_db_serum = total_db_serum.sort_values(by=["year", "month"],
                                                                     ascending=(True, True))
        total_db_serum = total_db_serum.reset_index(drop=True)
        struct_df[index_igG] = total_db_serum
        struct_df_pos[index_igG] = total_positive
        struct_df_pos[index_igG] = struct_df_pos[index_igG].set_index("Plate_ID")
        struct_df[index_igG] = struct_df[index_igG].set_index("Patient_or_Control_ID")

    df = struct_df[0]
    df.loc[:, "NC_EC50"] = struct_df[0].loc[:, "IC50"]
    df.loc[:, "Spike_EC50"] = struct_df[1].loc[:, "IC50"]
    df.loc[:, "RBD_EC50"] = struct_df[2].loc[:, "IC50"]
    df.loc[:, "NC_pEC50"] = -np.log10(struct_df[0].loc[:, "IC50"])
    df.loc[:, "Spike_pEC50"] = -np.log10(struct_df[1].loc[:, "IC50"])
    df.loc[:, "RBD_pEC50"] = -np.log10(struct_df[2].loc[:, "IC50"])
    df.loc[:, "NC_Plate_ID"] = struct_df[0].loc[:, "Plate_ID"]
    df.loc[:, "Spike_Plate_ID"] = struct_df[1].loc[:, "Plate_ID"]
    df.loc[:, "RBD_Plate_ID"] = struct_df[2].loc[:, "Plate_ID"]
    df = df.drop(['Plate_ID', 'target', 'Ab_conc', 'KD', 'IC50',
                  'Max_fittable_IC50', 'Error', ],
                 axis=1)
    df.loc[:, "Normalise"] = False
    df_pos = struct_df_pos[0]
    df_pos.loc[:, "NC_EC50"] = struct_df_pos[0].loc[:, "IC50"]
    df_pos.loc[:, "Spike_EC50"] = struct_df_pos[1].loc[:, "IC50"]
    df_pos.loc[:, "RBD_EC50"] = struct_df_pos[2].loc[:, "IC50"]
    df_pos.loc[:, "NC_pEC50"] = -np.log10(struct_df_pos[0].loc[:, "IC50"])
    df_pos.loc[:, "Spike_pEC50"] = -np.log10(struct_df_pos[1].loc[:, "IC50"])
    df_pos.loc[:, "RBD_pEC50"] = -np.log10(struct_df_pos[2].loc[:, "IC50"])
    df_pos = df_pos.drop(['target', 'Ab_conc', 'KD', 'IC50',
                          'Max_fittable_IC50', 'Error', ],
                         axis=1)
    return df, df_pos


def filter_USZ(df, adhoc_file):
    mask_usz = df["cohort"] == "hospital patient"
    mask_KIM = (df["subcohort"] == "KIM_neg") | (df["subcohort"] == "KIM_pos")
    df_USZ = df.loc[mask_usz & ~mask_KIM]
    df_USZ = df_USZ.loc[df_USZ[
                                ['NC_pEC50', 'Spike_pEC50', 'RBD_pEC50']
                                ].dropna().index]

    path_annotation = glob((os.path.join(path_to_main_folder,
                                         "**",
                                         adhoc_file)),
                           recursive=True)
    df_annotation = pd.read_excel(path_annotation[0])
    df_annotation = df_annotation.set_index("j_number")
    df_annotation = df_annotation[~pd.isna(df_annotation.research_id)]
    df_USZ = df_USZ[df_USZ.index.isin(df_annotation.index)]

    mask_pos = df_annotation.datediff_days_lab_request_and_covid_result < -14

    mask_clinicalconfirm = df_annotation.Clinically_manifest_COVID_ == "Yes"

    df_USZ.loc[df_annotation.loc[mask_pos & mask_clinicalconfirm].index,
               "hasCovid"] = 1
    allVals, castTblUnused, knownPos = Probability_field_FMW(df_USZ)
    df_USZ.loc[allVals.index, "posterior"] = allVals

    df_USZ.loc[:, "age"] = df_annotation.age
    df_USZ.loc[:, "sex"] = df_annotation.sex
    df_USZ.loc[:, "clinic"] = df_annotation.treatment_unit_clinic
    df_USZ.loc[:, "unique_PatientID"] = df_annotation.research_id
    df_USZ.loc[:, "Once_PCR_positive"] = df_annotation.once_covid_positive
    df_USZ.loc[:, "delta_day_PCRpos"] = df_annotation.datediff_days_lab_request_and_covid_result
    return df_USZ
