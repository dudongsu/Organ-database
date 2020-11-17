"""
Geometry analysis for two matrics
"""
#%load_ext autoreload
#%autoreload 2
#%matplotlib notebook
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
from pathlib import Path
from collections import defaultdict 
import dicom_contour.contour as dcm
from scipy import ndimage, sparse
import medpy.metric as metric
from dicom_contour.contour import get_contour_file,get_roi_names, coord2pixels, cfile2pixels, plot2dcontour, slice_order, get_contour_dict, get_data,  create_image_mask_files, fill_contour, get_data_new, get_contour_files

#### Analyze the metrics for single organ ##############################

def metrics_organ(mask_A, mask_B, zxy_dimention):
    # mask_A is the result and mask_B is the reference 
    if(mask_A.shape != mask_B.shape):
        raise Exception("matrix is not in same shape")
    # calculate the dice coefficient
    metrics = {}
    dice = metric.binary.dc(mask_A, mask_B)
    Jaccard = metric.binary.jc(mask_A, mask_B)
    ASD = metric.binary.asd(mask_A, mask_B, zxy_dimention)
    Hausdorff = metric.binary.hd(mask_A, mask_B, zxy_dimention)
    mask_diff = np.subtract(mask_A,mask_B)
    maskA_minus_maskB = np.where(mask_diff>0,mask_diff,0)
    maskB_minus_maskA = np.where(mask_diff<0,mask_diff,0)
    if 0 == np.count_nonzero(maskA_minus_maskB):
        HD95_1 = 0
    else:
        HD95_1 = metric.binary.hd95(maskA_minus_maskB, mask_B, zxy_dimention, connectivity=1)
    if 0 == np.count_nonzero(maskB_minus_maskA):
        HD95_2 = 0
    else:
        HD95_2 = metric.binary.hd95(maskB_minus_maskA, mask_A, zxy_dimention, connectivity=1)
    HD95 = max(HD95_1, HD95_2)
    metrics['dice']=dice
    metrics['Jaccard'] = Jaccard
    metrics['ASD'] = ASD
    metrics['Hausdorff'] = Hausdorff
    metrics['HD95'] = HD95
    return metrics


######################################## loop over one patient all organs ##################

# loop over all patient to get the baseline geometry relation between 
def loop_metrics_organs(cu1_data, cu2_data):
    
    zxy_dimention = [0,0,0]
    metrics_all_organ =  {}
    for organ in cu1_data.keys():
        if(organ !='zxy_dimention'):
            continue
        zxy_dimention = cu1_data['zxy_dimention']
    
    if zxy_dimention == [0,0,0]:
        raise Exception("no dimention found for this patient datasets ")
    
    for organ in cu1_data.keys():
        if(organ=='image'  or organ=='zxy_dimention' or ('stdName' not in curr_data[organ].keys()) ):
            continue
            
        if(curr_data[organ]['stdName']=='NA'):
            continue
            
        stdName_1 = curr_data[organ]['stdName']
            
        for organ_ref in curr_data.keys():    
            if(organ_ref=='image'  or organ_ref=='zxy_dimention' or ('stdName' not in curr_data[organ_ref].keys())):
                continue
            if(curr_data[organ_ref]['stdName']!=stdName_1):
                continue
            metrics = metrics_organ(curr_data[organ]['mask'], curr_data[organ_ref]['mask'], zxy_dimention)
            metrics_all_organ[stdName_1] = metrics
        
    return metrics_all_organ

