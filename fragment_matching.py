# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 13:34:48 2022

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
precursor_matches = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Lauren_Results\precursor_matching\A_precursor_matches.csv"#path to directory with inclusion lists
final_dir =r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Lauren_Results\fragment_matching" #path to directory for final, processed data

sample_name = 'A'

error_precursor = 20 #+/- ppm, for precursor
error_fragment = 0.02 #+/- Da, for fragment ion, charge state 1

precursor_charges = [7]
fragment_charges = [1,2,3,4,5,6,7,8]

h_mass = 1.00784

precursor_matches = pd.read_csv(precursor_matches)

fragment_database = pd.read_csv(fragment_list_import)

spectra_read= pd.read_csv(spectra_import, sep=" ",skiprows=[0], names= ["m/z", "resolution", "charge", "intensity","MS2",'scan_number','precursor_charge','empty','empty2'])

spectra_value = spectra_read.rename(columns={'m/z': 'Fragment m/z', 'resolution': 'Fragment resolution','charge':'Fragment ion charge','intensity':
                                             'Fragment intensity','MS2':'Actual precursor','scan_number':'Scan #','precursor_charge':'Precursor actual charge'})

spectra_value = spectra_value.drop(spectra_value[spectra_value['Precursor actual charge']==0].index) #remove charges equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['Fragment ion charge']==0].index) #remove charges equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['Fragment intensity']<100].index) #remove intensity values less than 100

peptides = []
peptides_dups = fragment_database['peptide'].values.tolist()
for jj in peptides_dups:
    if jj not in peptides:
        peptides.append(jj)

fragment_matches = pd.DataFrame()
for jk in peptides:
    fragment_database_filtered = fragment_database[fragment_database['peptide'] == jk]
    
    for a in fragment_charges: 
        theo_fragment = fragment_database_filtered[str(a)].values.tolist() #theoretical fragments for each charge
        seq_exp = precursor_matches[precursor_matches['Fragment ion charge'] == a] #actual fragments for each charge
        for u in theo_fragment:
            theo_frag_M = (u * a) - (h_mass * a)
            seq_exp['Actual fragment M'] = (seq_exp['Fragment m/z'] * a) - (h_mass * a)
            seq_exp['Theoretical fragment m/z'] = u
            seq_exp['Peptide'] = jk
            seq_exp['Theoretical fragment M'] = theo_frag_M
            seq_exp['Fragment error (Da)'] = abs(seq_exp['Actual fragment M'] - theo_frag_M)
            seq_exp_filter = seq_exp.sort_values(by='Fragment error (Da)')
            seq_exp_filter = seq_exp_filter[seq_exp_filter['Fragment error (Da)'] <= error_fragment]    
            
            if len(seq_exp_filter) > 0:
                fragment_matches = fragment_matches.append(seq_exp_filter) 
file_name = sample_name + '_fragment_matches.csv'
file_out_path = final_dir + '\\' + file_name
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        fragment_matches.to_csv(filec,index=False)