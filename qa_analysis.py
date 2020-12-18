from IPython.display import display
import json
import pickle
import copy
import os
import pandas as pd
import recordlinkage
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
import pydicom as dicom
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import pandas as pd
from pathlib import Path
from collections import defaultdict 
import dicom_contour.contour as dcm
from scipy import ndimage, sparse
import medpy.metric as metric
from dicom_contour.contour import get_contour_file,get_roi_names, coord2pixels, cfile2pixels, plot2dcontour, slice_order, get_contour_dict, get_data,  create_image_mask_files, fill_contour, get_data_new, get_contour_files
from Metrics_Analysis import metrics_organ, loop_metrics_organs
from NameMatching import generate_mask, generate_physician_mask, Rename_fuzzy, geometry_relation, dir_similirity, angle_3axl, loop_patients_geoStat, process_geo_relation, check_stdName
from NameMatching import df_stat
from rename_file import rename_ct


def create_database(parent_path):
    geo_relation = loop_patients_geoStat(parent_path, standard_list)
    geo_stat = process_geo_relation(geo_relation, 3)
    return geo_relation, geo_stat


def loop_test_patients(parent_path, standard_list, geo_relation,geo_stat):
    subdirs = [os.path.join(parent_path, o) for o in os.listdir(parent_path) if os.path.isdir(os.path.join(parent_path,o))]
    df = pd.DataFrame(columns=['patient','orignal name','matched name'])
    i = 1
    for path in subdirs:
        patient_name = path.split('\\')[-1]
        print('work on patient', path)
        test_patient = generate_physician_mask(path)
        Rename_fuzzy(test_patient, 'Physician', standard_list, 20)
        
        ### random switch some of the names
        
        data_newName = check_stdName(test_patient['Physician'], standard_list, geo_relation, 99, geo_stat)
        
        for organ in data_newName.keys():
            if organ=='image' or organ =='zxy_dimention':
                continue
        
            if('PTV' in organ) or ('GTV' in organ) or ('CTV' in organ) or ('ITV' in organ):
                continue
            
            if(data_newName[organ]['stdName'][-1] =='*'):
                df.loc[i] = [patient_name, organ, data_newName[organ]['stdName']]
                i = i+1
           
    display(df)
    return df
                                 

