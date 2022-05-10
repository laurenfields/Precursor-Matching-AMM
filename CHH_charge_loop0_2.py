# -*- coding: utf-8 -*-
"""
Created on Tue May  3 10:01:01 2022

@author: lawashburn
"""

import os
import csv
import pandas as pd
import numpy as np
from datetime import datetime
now = datetime.now()

spectra_import = r"C:\Users\lawashburn\Documents\ForLauren3\ForLauren3\NhuData\A\A_ms2_output_list.txt" #path to directory containing spectra
fragment_list_import = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Sample_input_data\by_cz_ions_CHH_QEHF_FusionLumos_LF.csv" #path to directory with list of target fragment ions
target_list_import = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Sample_input_data\precursor_list_CHH_QEHF_Lumos.csv"#path to directory with inclusion lists
final_dir =r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Lauren_Results\zero_charge_assign" #path to directory for final, processed data

sample_name = 'A'

error_precursor = 20 #+/- ppm, for precursor
error_fragment = 0.02 #+/- Da, for fragment ion, charge state 1

precursor_charges = [7]
fragment_charges = [1,2,3,4,5,6,7,8]

h_mass = 1.00784

print('loading files', datetime.now())

spectra_read= pd.read_csv(spectra_import, sep=" ",skiprows=[0], names= ["m/z", "resolution", "charge", "intensity","MS2",'scan_number','precursor_charge','empty','empty2'])
spectra_value = spectra_read.rename(columns={'m/z': 'Fragment m/z', 'resolution': 'Fragment resolution','charge':'Fragment ion charge','intensity':
                                             'Fragment intensity','MS2':'Actual precursor','scan_number':'Scan #','precursor_charge':'Precursor actual charge'})

spectra_value = spectra_value.drop(spectra_value[spectra_value['Precursor actual charge']==0].index) #remove charges equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['Fragment ion charge']!=0].index) #remove charges not equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['Fragment intensity']<100].index) #remove intensity values less than 100

target_list = pd.read_csv(target_list_import)

zero_precursor_matching = pd.DataFrame()

for b in precursor_charges:
        precursor_target = target_list[str(b)].values.tolist()
        spectra_value_filtered = spectra_value[spectra_value['Precursor actual charge'] == b]
        for c in precursor_target:
            spectra_value_filtered['Precursor error (ppm)'] =  ((abs(spectra_value_filtered['Actual precursor'] - c))/c) * 1E6
            spectra_value_filtered_err = spectra_value_filtered[spectra_value_filtered['Precursor error (ppm)'] <= error_precursor]
            spectra_value_filtered_err['Theoretical precursor'] = c
            spectra_value_filtered_err['Precursor theoretical charge'] = b
            if len(spectra_value_filtered_err) > 0:
                zero_precursor_matching = zero_precursor_matching.append(spectra_value_filtered_err)

file_name = sample_name + '_zero_precursor_matches.csv'
file_out_path = final_dir + '\\' + file_name
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        zero_precursor_matching.to_csv(filec,index=False)