from IPython.display import display
import numpy as np
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
import pandas as pd
import sys
import glob
from pathlib import Path
from collections import defaultdict 
import dicom_contour.contour as dcm
from scipy import ndimage, sparse
import medpy.metric as metric
from dicom_contour.contour import get_contour_file,get_roi_names, coord2pixels, cfile2pixels, plot2dcontour, slice_order, get_contour_dict, get_data,  create_image_mask_files, fill_contour, get_data_new, get_contour_files
from Metrics_Analysis import metrics_organ, loop_metrics_organs
from NameMatching import generate_mask, generate_physician_mask, Rename_fuzzy, geometry_relation, dir_similirity, angle_3axl, loop_patients_geoStat, process_geo_relation, check_stdName
from NameMatching import df_stat
from Rename_analysis import loop_test_patients
from rename_file import rename_ct

standard_list =["Tongue_Base", "Tongue_Oral", "Trachea", "Bone", "SpinalCord", "Scar", "Retinas", "Retina_R", "Retina_L", "Lens_R", "Lens_L", "Parotid_R","Parotid_L","Parotids","OpticNrv_R","OpticNrv_L", "OpticNrv","OpticChiasm", "Musc_Constrict_S", "Musc_Constrict_M", "Musc_Constrict_I", "Musc_Constrict", "Lungs", "Lips","Larnx_SG","Larnyx","Glnd_Submands","Glnd_Submand_L","Glnd_Submand_R", "Eyes", "Eye_R", "Eye_L","Esophagus","Ear_Internal_R", "Ear_Internal_L","Cochlea_R","Cochlea_L","Cochlea","Brainstem","Brain","Brachialplexs","BrachialPlex_R", "BrachialPlex_L", "Bone_Mandible","Body","Mouth_Floor",'Cavity_Oral','Pitutary','Thyroid','Carotid','Carotid_L','Carotid_R',"Arytenoid", "Arytenoid_L", "Arytenoid_R", "Ethmoid", "Occipital","Parietal_L","Parietal_R","Parietal","Sphenoid","Temporal", "Temporal_L", "Temporal_R", "Cavity_Nasal", "Cerebrum", "Cerebrum_L","Cerebrum_R","Glnd_Lacrimal","Glnd_Lacrimal_L","Glnd_Lacrimal_R","Glottis", "Nasopharynx"]


def create_database(parent_path):
    geo_relation = loop_patients_geoStat(parent_path, standard_list)
    geo_stat = process_geo_relation(geo_relation, 3)
    return geo_relation, geo_stat

def test(parent_path):
    loop_test_patients(parent_path, standard_list, geo_relation,geo_stat)
    




