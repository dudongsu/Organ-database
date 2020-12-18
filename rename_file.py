import numpy as np
import json
import pickle
import copy
import os
import pydicom as dicom


def rename_ct(file_path):
    
    ct_files = [os.path.join(file_path, o) for o in os.listdir(file_path) if not os.path.isdir(os.path.join(file_path,o))]
    
    for ct_file in ct_files:
       
        f = dicom.dcmread(ct_file)
        if(f.Modality != 'CT'):
            continue
        ct_name = f.SOPInstanceUID
        old_file = ct_file
        new_file = file_path + '\CT.'+ ct_name+'.dcm'
        print(old_file)
        print(new_file)
        os.rename(old_file,new_file)
    