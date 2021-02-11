"""

~~~~~~~~~~

Important functions for structure re-naming
"""
#%load_ext autoreload
#%autoreload 2
#%matplotlib notebook
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

#### generate mask function

def generate_mask(patient_path):
    path = patient_path
    contour_files = get_contour_files(path)
    all_data = defaultdict()
    for contour_file in contour_files:
        if path[-1] != '/': path += '/'
        f = dicom.dcmread(path + contour_file)
        RS_name = f.StructureSetLabel
        if "MedMind" in RS_name:
            RS_name = "MedMind"
        elif "DLCExpert" in RS_name:
            RS_name = "Mirada"
        else:
            RS_name = "Physician"
        
      #  print("work on RS structure ",RS_name,'+++++++++++++++++++++++++++')
        all_data[RS_name] = {}
        all_data[RS_name]["zxy_dimention"] =[]
        for s in os.listdir(path):
            img = dicom.dcmread(path + '/' + s)
            if hasattr(img, 'pixel_array'):  # to ensure not to read contour file
                img_arr = img.pixel_array
                # physical distance between the center of each pixel
                x_spacing, y_spacing = float(img.PixelSpacing[0]), float(img.PixelSpacing[1])
                slice_thickness = float(img.SliceThickness)
                all_data[RS_name]["zxy_dimention"] = [slice_thickness, x_spacing, y_spacing]
                break
        
        roi_seq_names = [roi_seq.ROIName for roi_seq in list(f.StructureSetROISequence)]
        for i in range(0,len(roi_seq_names)):
            
            roi_name = roi_seq_names[i]
       #     print("work on organ", roi_name,'_____________')
            ### convert a organ name to known matching, PC --- Musc_Constrict
            if RS_name == "Physician":
                if len(roi_name)>2 and roi_name[0:2]=='PC'  :
                    roi_name = 'Musc_Constrict' + roi_name[3:len(roi_name)] 
       #             print("the roi name converted to ", roi_name,'============')
                if len(roi_name)>3 and roi_name[0:3]=='SMG':
                    roi_name = 'Glnd_Submand_' + roi_name[3:len(roi_name)] 
                   # print("the roi name converted to ", roi_name,'============')
            
            #
            
            ## check if this contour is empty contour, if empty, skip it
            RTV_temp = f.ROIContourSequence[i]
            if not hasattr(RTV_temp, 'ContourSequence'):
                continue
            all_data[RS_name][roi_name] = {}
            img_voxel, mask_voxel = get_data_new(path, contour_file,i)
            all_data[RS_name][roi_name]["mask"] = mask_voxel
            all_data[RS_name][roi_name]["index"] = i
            all_data[RS_name]["image"] = img_voxel
    
    return all_data


############################# generate mask for physician contour only
def generate_physician_mask(patient_path):
    path = patient_path
    contour_files = get_contour_files(path)
    all_data = defaultdict()
    for contour_file in contour_files:
        if path[-1] != '/': path += '/'
        f = dicom.dcmread(path + contour_file)
        RS_name = f.StructureSetLabel
        if "MedMind" in RS_name:
            RS_name = "MedMind"
        elif "DLCExpert" in RS_name:
            RS_name = "Mirada"
        else:
            RS_name = "Physician"
        
        if RS_name != "Physician":
            continue
        print("work on RS structure ",RS_name,'+++++++++++++++++++++++++++')
        all_data[RS_name] = {}
        all_data[RS_name]["zxy_dimention"] =[]
        for s in os.listdir(path):
            img = dicom.dcmread(path + '/' + s)
            if hasattr(img, 'pixel_array'):  # to ensure not to read contour file
                img_arr = img.pixel_array
                # physical distance between the center of each pixel
                x_spacing, y_spacing = float(img.PixelSpacing[0]), float(img.PixelSpacing[1])
                slice_thickness = float(img.SliceThickness)
                all_data[RS_name]["zxy_dimention"] = [slice_thickness, x_spacing, y_spacing]
                break
        
        roi_seq_names = [roi_seq.ROIName for roi_seq in list(f.StructureSetROISequence)]
        for i in range(0,len(roi_seq_names)):
            
            roi_name = roi_seq_names[i]
          #  print("work on organ", roi_name,'_____________')
            ### convert a organ name to known matching, PC --- Musc_Constrict
            if RS_name == "Physician":
                if len(roi_name)>2 and roi_name[0:2]=='PC'  :
                    roi_name = 'Musc_Constrict' + roi_name[3:len(roi_name)] 
           #         print("the roi name converted to ", roi_name,'============')
                if len(roi_name)>3 and roi_name[0:3]=='SMG':
                    roi_name = 'Glnd_Submand_' + roi_name[3:len(roi_name)] 
           #         print("the roi name converted to ", roi_name,'============')
            
            #
            
            ## check if this contour is empty contour, if empty, skip it
            RTV_temp = f.ROIContourSequence[i]
            if not hasattr(RTV_temp, 'ContourSequence'):
                continue
            all_data[RS_name][roi_name] = {}
            img_voxel, mask_voxel = get_data_new(path, contour_file,i)
            all_data[RS_name][roi_name]["mask"] = mask_voxel
            all_data[RS_name][roi_name]["index"] = i
            all_data[RS_name]["image"] = img_voxel
    
    return all_data


############################ Rename the structure's names #############################################

def Rename_fuzzy(all_data, RS_name, standard_list, thre):
    if RS_name not in all_data:
        print(RS, "is not in all_data!")
    curr_data = all_data[RS_name]
    name_list = list(curr_data.keys())
    scores = {}
    standard_list_lower = standard_list.copy()
    
    for i in range(0,len(standard_list)):
        standard_list_lower[i] = standard_list_lower[i].lower()
    
    standard_map = dict()
    
    for i in range(0,len(standard_list)):
        standard_map[standard_list_lower[i]] = standard_list[i]
    
    
    for organ in name_list:
        if(organ=='image' or organ=='zxy_dimention'):
            continue
        
        if('PTV' in organ) or ('GTV' in organ) or ('CTV' in organ) or ('ITV' in organ):
            print(organ, ' is  a target structure')
            continue
  #      print('start work on organ')
        
        cand = process.extractOne(organ.lower(), standard_list_lower, scorer=fuzz.ratio)
        if(cand[1]>=thre):
            if(cand[0] in scores.keys()):
                score1 = fuzz.ratio(organ.lower(), cand[0])
                score2 = fuzz.ratio(scores[cand[0]].lower(), cand[0])
                if(score1>score2):
                    curr_data[scores[cand[0]]]['stdName'] = "NA"
                    curr_data[organ]['stdName'] = standard_map[cand[0]]
                    scores[cand[0]] = organ
                    print(organ," --> ", curr_data[organ]["stdName"], ' score is ', cand[1])
                else:
                    curr_data[organ]['stdName'] = "NA"
            else:
                curr_data[organ]['stdName'] = standard_map[cand[0]]
                scores[cand[0]] = organ
                print(organ," --> ", curr_data[organ]["stdName"], ' score is ', cand[1])
        else:
          #  print(curr_data[organ])
            curr_data[organ]['stdName'] = "NA"
     #       print(organ," --> ", curr_data[organ]["stdName"])
    for organ in name_list:
        if(organ=='image'):
            continue
        
    return curr_data


###########################  geometry_relation ###################################################################

## calculate relative position and distance for two organs
def geometry_relation(mask_A, mask_B, vec_ref=[[1,0,0],[0,1,0],[0,0,1]]):
    """Get the geometry relationship between two strucures, we use this function to find the L to R Parotid direction.
    
    """
    
    if mask_A.shape != mask_B.shape:
        raise ValueError('Please resample the data to same position is!')
    result = {}
    center_A = ndimage.measurements.center_of_mass(mask_A) # [SUP-INF, ANT-POST, Right-Left]
    center_B = ndimage.measurements.center_of_mass(mask_B)
    vector = [center_A[0]-center_B[0], center_A[1]-center_B[1], center_A[2]-center_B[2]]
    distance = np.linalg.norm(vector)
    vector = vector /distance
    result['distance'] = distance
    result['direction'] = vector
    result['angle'] = angle_3axl(vector, vec_ref)
    return result


############################# dis_similiairty ########################################################################

def dir_similirity(vec1, vec2):
    """
    
    """
    
    
    if len(vec1)!=3 or len(vec2)!=3:
        raise ValueError('need to be 3d vector!')
    
    vector_test = test_metric['direction']
    vector_ref = ref_metric['direction']
    unit_vector_1 = vector_test / np. linalg. norm(vector_test)
    unit_vector_2 = vector_ref / np. linalg. norm(vector_ref)
    dot_product = np. dot(unit_vector_1, unit_vector_2)
    angle3d = np. arccos(dot_product)
    vector_test_2d = vector_test[1:3]
    vector_ref_2d = vector_ref[1:3]
    unit_vector_1_2d = vector_test_2d / np. linalg. norm(vector_test_2d)
    unit_vector_2_2d = vector_ref_2d / np. linalg. norm(vector_ref_2d)
    dot_product = np. dot(unit_vector_1_2d, unit_vector_2_2d)
    angle2d = np. arccos(dot_product)
    result = {}
    result['angle3d'] = angle3d
    result['angle2d'] = angle2d
    return result


#####################################################################################

def angle_3axl(vec1, ref):
    """ get the vector angle to reference z, x, y axis
    
    """
    
    if len(vec1)!=3:
        raise ValueError('need to be 3d vector!')
    
    # vect1 direction is [SUP-INF, ANT-POST, Right-Left]
    
    z_axis = ref[0]
    x_axis = ref[1]
    y_axis = ref[2]
    vec1_2D = [0,vec1[1],vec1[2]]
    
    # get z axis angle first
    unit_vector_1 = vec1 / np. linalg. norm(vec1)
    dot_product_z = np. dot(unit_vector_1, z_axis)
    angle_z = np. arccos(dot_product_z)
    
    # get x and y axis angles
    unit_vector_1 = vec1_2D / np. linalg. norm(vec1_2D)
    dot_product_x = np. dot(unit_vector_1, x_axis)
    angle_x = np. arccos(dot_product_x)
    dot_product_y = np. dot(unit_vector_1, y_axis)
    angle_y = np. arccos(dot_product_y)
    result = [angle_z, angle_y, angle_x]
    
    return result

################################################################################################

def find_axis(mask):
    """ find the axis of a cylindrical shape in sup-inf direction
    
    """
    
    index = np.where(mask==1)
    top = max(index[0])
    bottom = min(index[0])
    center_top = ndimage.measurements.center_of_mass(mask[top,:,:])
    center_bottom = ndimage.measurements.center_of_mass(mask[bottom,:,:])
    vec =  [top-bottom, center_top[0]-center_bottom[0], center_top[1]-center_bottom[1]]
    dis = np.linalg.norm(vec)
    vec = vec /dis
    return vec
        


#################################################################################################
def find_ref_vec(curr_data):
    """
        find reference axis, z--Brainstem, x
    
    """
    ref_LR = np.array([0,0,0])
    ref_SI = np.array([0,0,0])
    for organ in curr_data.keys():
        if(organ=='image'  or organ=='zxy_dimention' or ('stdName' not in curr_data[organ].keys()) ):
            continue
        if(curr_data[organ]['stdName']!='Brainstem'):
            continue
        mask = curr_data[organ]['mask']
        ref_SI = find_axis(mask)
        
    for organ_R in curr_data.keys():
        if(organ_R=='image'  or organ_R=='zxy_dimention' or ('stdName' not in curr_data[organ_R].keys()) ):
            continue
            
        if(curr_data[organ_R]['stdName']!='Parotid_R'):
            continue
            
        for organ_L in curr_data.keys():
            if(organ_L=='image'  or organ_L=='zxy_dimention' or ('stdName' not in curr_data[organ_L].keys()) ):
                continue
            if(curr_data[organ_L]['stdName']!='Parotid_L'):
                continue
            result = geometry_relation(curr_data[organ_R]['mask'], curr_data[organ_L]['mask'],vec_ref=[[1,0,0],[0,1,0],[0,0,1]])
            ref_LR =  result['direction']
        
    #print('ref_LR is', ref_LR)
    #print('ref_SI is', ref_SI)
    if  (ref_LR==np.array([0,0,0])).all() or (ref_SI==np.array([0,0,0])).all():
        return [np.array([0,0,0]),np.array([0,0,0]) ,np.array([0,0,0])]
    ref_AP = np.cross(ref_SI, ref_LR) 
    ref_vec = [ref_SI, ref_AP, ref_LR]
    return ref_vec


##################################################################################
# loop over all patient to get the baseline geometry relation between 
def loop_patients_geoStat(parent_path, standard_list):
    subdirs = [os.path.join(parent_path, o) for o in os.listdir(parent_path) if os.path.isdir(os.path.join(parent_path,o))]
    
    geo_relation =  {}
    for path in subdirs:
        print('work on patient', path)
        all_data = generate_physician_mask(path)
        if 'Physician' not in all_data.keys():
            print("patient ", path, "doesn't have physician contour")
            continue
        RS_name = 'Physician'
        curr_data = Rename_fuzzy(all_data, RS_name, standard_list, 70)
        
        ## Have to get the reference vector direction for each of the organ relative positions######################
        ref_LR = np.array([0,0,0])
        ref_SI = np.array([0,0,0])
        for organ in curr_data.keys():
            if(organ=='image'  or organ=='zxy_dimention' or ('stdName' not in curr_data[organ].keys()) ):
                continue
            if(curr_data[organ]['stdName']!='Brainstem'):
                continue
            mask = curr_data[organ]['mask']
            ref_SI = find_axis(mask)
        
        for organ_R in curr_data.keys():
            if(organ_R=='image'  or organ_R=='zxy_dimention' or ('stdName' not in curr_data[organ_R].keys()) ):
                continue
            
            if(curr_data[organ_R]['stdName']!='Parotid_R'):
                continue
            
            for organ_L in curr_data.keys():
                if(organ_L=='image'  or organ_L=='zxy_dimention' or ('stdName' not in curr_data[organ_L].keys()) ):
                    continue
                if(curr_data[organ_L]['stdName']!='Parotid_L'):
                    continue
                result = geometry_relation(curr_data[organ_R]['mask'], curr_data[organ_L]['mask'],vec_ref=[[1,0,0],[0,1,0],[0,0,1]])
                ref_LR =  result['direction']
        
        print('ref_LR is', ref_LR)
        print('ref_SI is', ref_SI)
        if  (ref_LR==np.array([0,0,0])).all() or (ref_SI==np.array([0,0,0])).all():
            continue
        ref_AP = np.cross(ref_SI, ref_LR) 
        
        ref_vec = [ref_SI, ref_AP, ref_LR] 
        #############################################################################
        for organ in curr_data.keys():
            
            if(organ=='image'  or organ=='zxy_dimention' or ('stdName' not in curr_data[organ].keys()) ):
                continue
            
            if(curr_data[organ]['stdName']=='NA'):
                continue
            
            stdName_1 = curr_data[organ]['stdName']
            if stdName_1 not in geo_relation.keys():
                geo_relation[stdName_1] = {}
            
            for organ_target in curr_data.keys():
                
                if(organ_target=='image'  or organ_target=='zxy_dimention' or ('stdName' not in curr_data[organ_target].keys())):
                    continue
                if(curr_data[organ_target]['stdName']=='NA'):
                    continue
                stdName_2 = curr_data[organ_target]['stdName']
                if stdName_2 not in geo_relation[stdName_1].keys():
                    geo_relation[stdName_1][stdName_2] = []
                
                print('geometry relation start to work on organs', organ, ' and ', organ_target)
                result = geometry_relation(curr_data[organ]['mask'], curr_data[organ_target]['mask'], ref_vec)    
                geo_relation[stdName_1][stdName_2].append(result)
        del all_data
        
    return geo_relation

##########################################################

def process_geo_relation(geo_relation, thre):
    
    geo_stat = {}
    print('initial run, get the relative direction for all organs, get the each organ distance')
    distance = {}
    for organ_1 in geo_relation.keys():
        geo_stat[organ_1]={}
        distance[organ_1]={}
        for organ_2 in geo_relation[organ_1].keys():
            if organ_2 == organ_1:
                continue
            
            
            all_ref = []
            
            all_vecs = geo_relation[organ_1][organ_2]
            
            dirs = [c['direction'] for c in all_vecs]
            dist_cur = [c['distance'] for c in all_vecs]
            angles = [c['angle'] for c in all_vecs]
            distance[organ_1][organ_2] =  np.mean([c for c in dist_cur])
            
            if(len(all_vecs)>=thre):
                print('working on relation ', organ_1,  ' to ', organ_2)
                geo_stat[organ_1][organ_2]={}
                ref_angle3 = []
                for d1 in angles:
                    ref_angle3.append(d1)
                geo_stat[organ_1][organ_2]['mean']=[0,0,0]
                geo_stat[organ_1][organ_2]['std'] = [0,0,0]
                geo_stat[organ_1][organ_2]['mean'][0] = np.mean([c[0] for c in ref_angle3])
                geo_stat[organ_1][organ_2]['mean'][1] = np.mean([c[1] for c in ref_angle3])
                geo_stat[organ_1][organ_2]['mean'][2] = np.mean([c[2] for c in ref_angle3])
                geo_stat[organ_1][organ_2]['std'][0] = np.std([c[0] for c in ref_angle3])
                geo_stat[organ_1][organ_2]['std'][1] = np.std([c[1] for c in ref_angle3])
                geo_stat[organ_1][organ_2]['std'][2] = np.std([c[2] for c in ref_angle3])
   #     sorted(distance[organ_1].items(), key = lambda kv:(kv[1], kv[0]))
        print('sorted distance is ', distance[organ_1])
    
    near_map = {}
    serial_organs = ['Esophagus', 'SpinalCord', 'Brainstem']
    ## for each organ already has std, assign it to multiple organs that close to it
    print("start organ mapping phase -------------------------")
    for organ_1 in geo_stat.keys():
        if(organ_1 in serial_organs): 
            continue 
        print('working on ', organ_1, '   ^^^  ')
        for organ_2 in geo_relation[organ_1].keys():
            if(organ_2==organ_1): continue
            if(organ_2 not in near_map.keys()):
                near_map[organ_2] = organ_1
            else:
                if distance[organ_1][organ_2]<distance[near_map[organ_2]][organ_2]:
                    near_map[organ_2] = organ_1 
    
    
    print('second run to detect the std for lack of data geo relations')
    for organ_1 in geo_relation.keys():
        for organ_2 in geo_relation[organ_1].keys():
            if organ_2 == organ_1:
                continue
            all_vecs = geo_relation[organ_1][organ_2]
            dirs = [c['direction'] for c in all_vecs]
            angles = [c['angle'] for c in all_vecs]
            if(len(all_vecs)<thre):
                print('working on relation ', organ_1,  ' to ', organ_2)
                geo_stat[organ_1][organ_2]={}
                ref_angle3 = []
                for d1 in angles:
                    ref_angle3.append(d1)
                geo_stat[organ_1][organ_2]['mean']=[0,0,0]
                geo_stat[organ_1][organ_2]['std'] = [0,0,0]
                geo_stat[organ_1][organ_2]['mean'][0] = np.mean([c[0] for c in ref_angle3])
                geo_stat[organ_1][organ_2]['mean'][1] = np.mean([c[1] for c in ref_angle3])
                geo_stat[organ_1][organ_2]['mean'][2] = np.mean([c[2] for c in ref_angle3])
                if (near_map[organ_2] in geo_stat[organ_1]):
                    geo_stat[organ_1][organ_2]['std'] = geo_stat[organ_1][near_map[organ_2]]['std']
                elif (near_map[organ_1] in geo_stat[organ_2]):
                    geo_stat[organ_1][organ_2]['std'] = geo_stat[organ_2][near_map[organ_1]]['std']
                elif (near_map[organ_1] in geo_stat[near_map[organ_2]]):
                    geo_stat[organ_1][organ_2]['std'] = geo_stat[near_map[organ_2]][near_map[organ_1]]['std']
    
    return geo_stat


####################################################################################

# check the name for specific RS dataset 
def check_stdName(cu_data, standard_list, geo_relation, threshold, geo_stat):
    
    std_cnt = 3
    # start to check if the reference organ detected in cur_data
    curr_data = copy.deepcopy(cu_data)
    available_ref = {}
    for organ in curr_data.keys():
        if organ=='image' or organ =='zxy_dimention':
            continue
        if('PTV' in organ) or ('GTV' in organ) or ('CTV' in organ) or ('ITV' in organ):
            continue
        if curr_data[organ]['stdName']=="Parotid_R" or curr_data[organ]['stdName']=="Parotid_L" or curr_data[organ]['stdName']=="Cochlea_L" or curr_data[organ]['stdName']=="Cochlea_L":
            if fuzz.ratio(organ.lower(), curr_data[organ]['stdName'].lower())>threshold:
                available_ref[curr_data[organ]['stdName']] = organ
    print("the available reference organs are ", available_ref.keys())
    
    ## for some organ, the reference geo-ralation is too few, that we have to get the standard deviation angels from other organs
    
    ref_vec = find_ref_vec(curr_data)
    if (ref_vec[0]==np.array([0,0,0])).all():
        print(" Couldn't find the reference axis!!")
        return None
    
    # start to check if there is any reference for each organ to check and check their relative position
    for organ in curr_data.keys():
        if organ=='image' or organ =='zxy_dimention':
            continue
        
        if('PTV' in organ) or ('GTV' in organ) or ('CTV' in organ) or ('ITV' in organ):
       #     print(organ, ' is  a target structure')
            continue
        
        
        print('Start to work on organ', organ, '+++++++++++++++++++++++++++++++++++++++' )
        if curr_data[organ]['stdName']=='NA':
            print('The organ', organ, " does not have initial stndard Name")
            continue
        cand = curr_data[organ]['stdName']
        fuzz_score = fuzz.ratio(organ.lower(), cand.lower())
        print('The name mathcing score is ', fuzz_score)
        if  fuzz.ratio(organ.lower(), cand.lower())>threshold:
            continue
    
        if cand not in geo_relation.keys():
            print("the organ", cand, "is not in geo_relation libary")
            curr_data[organ]['stdName'] = curr_data[organ]['stdName']+'*'
            continue
        
        
        # only consider Parotid and cocolea as reference organ
        if ("Parotid_R" in geo_stat[cand].keys()) and ("Parotid_L" in geo_stat[cand].keys()) and ("Parotid_R" in available_ref.keys()) and ("Parotid_L" in available_ref.keys()):
            ref1 = "Parotid_R"
            ref2 = "Parotid_L"
        elif ("Cochlea_R" in geo_stat[cand].keys()) and ("Cochlea_L" in geo_stat[cand].keys()) and ("Cochlea_L" in available_ref.keys()) and ("Cochlea_R" in available_ref.keys()):
            ref1 = "Cochlea_R"
            ref2 = "Cochlea_L"
        else: 
            print("the organ", cand, "is not in geo_relation libary with reference organs")
            curr_data[organ]['stdName'] = curr_data[organ]['stdName']+'*'
            continue
        
        ref1_stat = {}
        ref2_stat = {}
        ref1_stat['mean'] = geo_stat[cand][ref1]['mean']
        ref1_stat['std'] = geo_stat[cand][ref1]['std']
        ref2_stat['mean'] = geo_stat[cand][ref2]['mean']
        ref2_stat['std'] = geo_stat[cand][ref2]['std']
        
        target1_dir_dist = geometry_relation(curr_data[organ]['mask'], curr_data[available_ref[ref1]]['mask'], ref_vec)
        target1_dir = target1_dir_dist['direction']
        target2_dir_dist = geometry_relation(curr_data[organ]['mask'], curr_data[available_ref[ref2]]['mask'], ref_vec)
        target2_dir = target2_dir_dist['direction']
        
        target_ref1_angles = target1_dir_dist['angle']
        target_ref2_angles = target2_dir_dist['angle']
        print('The geo stat range is', ref1_stat['mean'], ' +- ', ref1_stat['std']) 
        print("The target and ref1 direction ", target_ref1_angles)
        print('The geo stat range is', ref2_stat['mean'], ' +- ', ref2_stat['std']) 
        print("The target and ref2 direction ", target_ref2_angles)
        
        
        serial_organs = ['Esophagus', 'SpinalCord', 'Brainstem']
        # depends on the target organ, we will ignore the angle in 3rd direction comparision comparison
    
        if target_ref1_angles[1] > ref1_stat['mean'][1]+std_cnt* ref1_stat['std'][1] or target_ref1_angles[1] < ref1_stat['mean'][1]-std_cnt* ref1_stat['std'][1]:
            curr_data[organ]['stdName'] = curr_data[organ]['stdName']+'*'
            print("angle 1 y-axis out of range for reference organ 1 ")
            continue
        
        if target_ref1_angles[2] > ref1_stat['mean'][2]+std_cnt* ref1_stat['std'][2] or target_ref1_angles[2] < ref1_stat['mean'][2]-std_cnt* ref1_stat['std'][2]:
            curr_data[organ]['stdName'] = curr_data[organ]['stdName']+'*'
            print("angle 2 y-axis out of range for reference organ 1")
            continue
        if organ not in serial_organs:
            if target_ref1_angles[0] > ref1_stat['mean'][0]+std_cnt* ref1_stat['std'][0] or target_ref1_angles[0] < ref1_stat['mean'][0]-std_cnt* ref1_stat['std'][0]:
                curr_data[organ]['stdName'] = curr_data[organ]['stdName']+'*'
                print("angle 3 z-axis out of range for reference organ 1")
                continue
       
        ### start to check the geometry for ref 2
        if target_ref2_angles[1] > ref2_stat['mean'][1]+std_cnt* ref2_stat['std'][1] or target_ref2_angles[1] < ref2_stat['mean'][1]-std_cnt* ref2_stat['std'][1]:
            curr_data[organ]['stdName'] = curr_data[organ]['stdName']+'*'
            print("angle 1 y-axis out of range for reference organ 2 ")
            continue
        
        if target_ref2_angles[2] > ref2_stat['mean'][2]+std_cnt* ref2_stat['std'][2] or target_ref2_angles[2] < ref2_stat['mean'][2]-std_cnt* ref2_stat['std'][2]:
            curr_data[organ]['stdName'] = curr_data[organ]['stdName']+'*'
            print("angle 2 y-axis out of range for reference organ 2")
            continue
        if organ not in serial_organs:
            if target_ref2_angles[0] > ref2_stat['mean'][0]+std_cnt* ref2_stat['std'][0] or target_ref2_angles[0] < ref2_stat['mean'][0]-std_cnt* ref2_stat['std'][0]:
                curr_data[organ]['stdName'] = curr_data[organ]['stdName']+'*'
                print("angle 3 z-axis out of range for reference organ 2")
                continue
    
    print("Final result -------------------------------------------------------------------------------")
    
    for organ in curr_data.keys():
        if organ=='image' or organ =='zxy_dimention':
            continue
        
        if('PTV' in organ) or ('GTV' in organ) or ('CTV' in organ) or ('ITV' in organ):
           # print(organ, ' is  a target structure')
            continue
        
        print(organ, '-------->', curr_data[organ]['stdName'])
    
    return curr_data



def df_stat(geo_relation):
    df = pd.DataFrame(np.nan, index=list(geo_relation.keys()), columns=list(geo_relation.keys()))
    for key1 in geo_relation.keys():
        for key2 in geo_relation[key1].keys():
            df.loc[key1][key2] = len(geo_relation[key1][key2])
            
    display(df)
    
    return df
    
    









