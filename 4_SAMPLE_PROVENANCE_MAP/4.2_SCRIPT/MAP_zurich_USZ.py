#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 10:51:43 2021
@author: rjacquat

This map can be execute directly with the correct path for the
plz_nach_monaten.xlsx, PLZO_PLZ.dbf, Postleitzahlen-Schweiz.xlsx, and
map_geo_admin_ch_KML_20210104060626.kml

plz_nach_monaten is confidential, containing info of antibody detected against
sarscov2, residence and time when blood has been selected.
PLZO_PLZ.dbf is a file of the swiss governenement free of use thanks to the
Open Governement Data that swiss federation used
https://shop.swisstopo.admin.ch/en/free-geodata
Postleitzahlen-Schweiz.xlsx is a file given by Swiss post which allow to linked
PLZ/ZIP code to which canton it refere to exclude patient that come from a
canton outside the canton of Zurich. (be aware that some PLZ are share between
                                      several cantons, the map is therefore a
                                      bit different that the canton of zurich)
map_geo_admin_ch_KML_20210104060626.kml is a home made files using
map.geo.admin.ch which show the 3 main lake boundary.of Zurich's canton

As I was not familiar with map creation, the construction of such map has been
possible using help from different source like:
https://rosenfelder.ai/create-maps-with-python/
https://github.com/interactivethings/swiss-maps
https://medium.com/@v.brusylovets/your-population-density-map-of-switzerland-d108301f3dac
"""
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

import tilemapbase
import warnings
import matplotlib.cbook
from matplotlib import cm
import shapely.speedups


border_ZIP = 'PLZO_SHP_LV03/PLZO_PLZ.dbf'
link_plz_namekanton = 'Postleitzahlen-Schweiz.xlsx'
fn_handdraw_lake = 'map_geo_admin_ch_KML_20210104060626.kml'
fn_USZ_patient_ZIP_monate = 'plz_nach_monaten.xlsx'
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)
shapely.speedups.enable()

cm_name = "YlGn"  # set "viridis"

PLZ = gpd.read_file(border_ZIP)
gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'

#  For creation of city border in orange need to merge ZIP/PLZ of the city
city = PLZ.loc[np.logical_and((PLZ.PLZ > 8000), (PLZ.PLZ < 8064))]
city = city.dissolve(by="STATUS")  # merge city of Zurich zipcode
city = city.to_crs("epsg:3857")
#  For boundary of the main 3 lakes in Zurich, the shape as been made by hand
#  in map.geo.admin.ch
lake_homemade = gpd.read_file(fn_handdraw_lake,
                              driver="KML")
lake_homemade = lake_homemade.to_crs({"init": "EPSG:3857"})
plz_list = pd.read_excel(link_plz_namekanton)
zh_canton_number = plz_list.loc[plz_list.Canton == "Zurich",
                                'Postleitzahl / Code Postal / Codice Postale']

seroprevalence_place = pd.read_excel(fn_USZ_patient_ZIP_monate)
MinimumNumberPatientperPlace = 50
Saturation_value_percent = 12
windows_size_month = 4
perarea = False

PLZ = PLZ.set_index("PLZ")
PLZ.loc[:, "PLZ"] = PLZ.index

PLZ = PLZ.to_crs("epsg:3857")
list_month = np.asarray((np.arange(0, 14-windows_size_month),
                         np.arange(windows_size_month, 14))).T
list_month = ((3, 7), (9, 13))
list_month = ((1, 7), (7, 13))
vmax_coef = 100 / Saturation_value_percent


tilemapbase.start_logging()
# # # Don't need if you have run before; DB file will already exist.
tilemapbase.init(create=True)
# # # Use open street map
t = tilemapbase.tiles.build_OSM()
KANTON = PLZ[PLZ.PLZ.isin(zh_canton_number)].dissolve(by="STATUS")
for months in list_month:
    PLZ = gpd.read_file(border_ZIP)

    fig, ax = plt.subplots(figsize=(10, 10),
                           dpi=120,
                           facecolor='w',
                           edgecolor='k')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    mask_year = seroprevalence_place.jahr == 2020
    mask_month = seroprevalence_place.monat >= months[0]
    mask_month = mask_month & (seroprevalence_place.monat < months[1])
    place = seroprevalence_place.loc[mask_year & mask_month]
    a = place.groupby("zip_code")
    sommes = a.sum()
    number_sample = sommes.loc[:, "sample_count_total"]
    sommes = sommes.loc[number_sample > MinimumNumberPatientperPlace]
    ratio = (sommes.loc[:, "sample_count_seropos"] /
             sommes.loc[:, "sample_count_total"])

    ZH_mask = PLZ.loc[PLZ.PLZ.isin(ratio.index)]

    colormap = cm.get_cmap(cm_name, 256)
    color = colormap(ratio[ZH_mask.PLZ] * vmax_coef)
    values = ratio[ZH_mask.PLZ] * vmax_coef

    KANTON.plot(ax=ax,
                facecolor=[0.5, 0.5, 0.5])

    my_office = (8.65, 47.43)
    degree_range = 0.4
    extent = tilemapbase.Extent.from_lonlat(my_office[0] - degree_range,
                                            my_office[0] + degree_range,
                                            my_office[1] - degree_range,
                                            my_office[1] + degree_range)
    extent = extent.to_aspect(1.0)
    extent = extent.to_project_3857()

    ZH_mask = ZH_mask.to_crs("epsg:3857")
    ZH_mask = ZH_mask.loc[ZH_mask.PLZ.isin(zh_canton_number)]
    cm_name = "YlGn"
    colormap = cm.get_cmap(cm_name, 256)
    color_number = (ratio[ZH_mask.PLZ])
    color_number = (ratio[ZH_mask.PLZ]) * vmax_coef
    color = colormap(color_number)

    ZH_mask.plot(ax=ax,
                 edgecolor="black",
                 linewidth=1,
                 facecolor=color)
    plotter = tilemapbase.Plotter(extent, t, width=600)
    plotter.plot(ax, alpha=0.4)
    ZH_mask.plot(ax=ax,
                 edgecolor="black",
                 linewidth=1,
                 facecolor=color)

    ZH_mask.plot(ax=ax,
                 edgecolor="black",
                 linewidth=1,
                 facecolor=[0, 0, 0, 0])

    city.plot(ax=ax,
              edgecolor=[1, 0.3, 0],
              linewidth=2.5,
              facecolor=[0, 0, 0, 0])
    KANTON.plot(ax=ax,
                edgecolor=[1, 0.6, 0.3],
                linewidth=2.5,
                facecolor=[0, 0, 0, 0])

    lake_homemade = lake_homemade.to_crs("epsg:3857")
    lake_homemade.plot(ax=ax,
                       edgecolor="black",
                       linewidth=1,
                       facecolor=[0.4, 0.5, 1])

    norm = cm.colors.Normalize(vmax=1 * 100 / vmax_coef,  #
                               vmin=np.min(ratio[ZH_mask.PLZ]))
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm,
                                          cmap=colormap),
                        ax=ax)
    cbar.set_label('Prevalence %', rotation=270)
    plt.title(f"Overall data between month {months[0]}-{months[1]}")
    print(f"Number PLZ bigger than 2% {np.sum(ratio[ZH_mask.PLZ]>0.02)}")

    cm_name = "RdPu"

    colormap = cm.get_cmap(cm_name, 256)
    number_sample = sommes.loc[:, "sample_count_total"]

    ZH_mask_plz = ZH_mask.set_index("PLZ")

    if perarea is True:
        color_number = number_sample[ZH_mask.PLZ] / (ZH_mask_plz.area / 10**6)
    else:
        PLZ = gpd.read_file(border_ZIP)
        PLZ = PLZ.set_index("PLZ")
        ZH_mask = PLZ.loc[place.loc[place.zip_code.isin(zh_canton_number),
                                    "zip_code"]]
        a = place.groupby("zip_code")
        number_sample = a.sum().loc[:, "sample_count_total"]
        color_number = number_sample[ZH_mask.index]
    color = colormap(color_number / np.nanmax(color_number))

    tilemapbase.start_logging()

    # Don't need if you have run before; DB file will already exist.
    tilemapbase.init(create=True)

    # Use open street map
    t = tilemapbase.tiles.build_OSM()

    my_office = (8.65, 47.43)
    degree_range = 0.4
    extent = tilemapbase.Extent.from_lonlat(my_office[0] - degree_range,
                                            my_office[0] + degree_range,
                                            my_office[1] - degree_range,
                                            my_office[1] + degree_range)
    extent = extent.to_aspect(1.0)
    extent = extent.to_project_3857()

    fig, ax = plt.subplots(figsize=(10, 10),
                           dpi=120,
                           facecolor='w',
                           edgecolor='k')

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    ZH_mask = ZH_mask.to_crs("epsg:3857")

    KANTON.plot(ax=ax,
                facecolor=[0.5, 0.5, 0.5])

    ZH_mask.plot(ax=ax,
                 edgecolor="black",
                 linewidth=1,
                 facecolor=color)

    ZH_mask.plot(ax=ax,
                 edgecolor="black",
                 linewidth=1,
                 facecolor=[0, 0, 0, 0])

    city.plot(ax=ax,
              edgecolor=[1, 0.3, 0],
              linewidth=2.5,
              facecolor=[0, 0, 0, 0])
    KANTON.plot(ax=ax,
                edgecolor=[1, 0.6, 0.3],
                linewidth=2.5,
                facecolor=[0, 0, 0, 0])
    plotter = tilemapbase.Plotter(extent, t, width=600)
    plotter.plot(ax, alpha=0.4)

    lake_homemade = lake_homemade.to_crs("epsg:3857")
    lake_homemade.plot(ax=ax,
                       edgecolor="black",
                       linewidth=1,
                       facecolor=[0.4, 0.5, 1])

    norm = cm.colors.Normalize(vmax=np.nanmax(color_number),
                               vmin=0)
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm,
                                          cmap=cm.get_cmap(cm_name, 256)),
                        ax=ax)
    if perarea is True:
        cbar.set_label('# of sample per km$^2$', rotation=270)
    else:
        cbar.set_label('# of sample', rotation=270)
    plt.title(f"Number of Sample collected, month {months[0]}-{months[1]}")
