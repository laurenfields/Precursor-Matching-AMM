# -*- coding: utf-8 -*-
"""
Created on Tue May  3 10:36:01 2022

@author: lawashburn
"""

import os
import csv
import pandas as pd
import numpy as np
from datetime import datetime
now = datetime.now()

#spectra_import = r"C:\Users\lawashburn\Documents\ForLauren3\ForLauren3\NhuData\A\A_ms2_output_list.txt" #path to directory containing spectra
fragment_list_import = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Sample_input_data\by_cz_ions_CHH_QEHF_FusionLumos_LF.csv" #path to directory with list of target fragment ions
target_list_import = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Sample_input_data\precursor_list_CHH_QEHF_Lumos.csv"#path to directory with precursor list
working_directory = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Lauren_Results\zero_charge_assign"
#final_dir =r"C:\Users\lawashburn\Documents\ForLauren3\ForLauren3\NhuData\A_Lauren_Run\Zero_Fragment_Matches\FD" #path to directory for final, processed data
precursor_matches = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Lauren_Results\zero_charge_assign\A_zero_precursor_matches.csv"

sample_name = 'A'

error_precursor = 20 #+/- ppm, for precursor
error_fragment = 0.02 #+/- Da, for fragment ion, charge state 1

precursor_charges = [7]
fragment_charges = [1,2,3,4,5,6,7,8]

h_mass = 1.00784

precursor_matches = pd.read_csv(precursor_matches)
fragment_database = pd.read_csv(fragment_list_import)
first_charge = fragment_charges[0]
       
fragment_matches = pd.DataFrame()



prelim_matches = len(precursor_matches)

new_charge = []

for x in range(0,prelim_matches):
    new_charge.append(first_charge)

precursor_matches['New charge assignment'] = new_charge

fragment_types_dups = fragment_database['peptide'].values.tolist()

fragment_types = []
for a in fragment_types_dups:
    if a not in fragment_types:
        fragment_types.append(a)

for b in fragment_types:
    fragment_database_filtered = fragment_database[fragment_database['peptide'] == b]   
    theo_fragment = fragment_database_filtered[str(first_charge)].values.tolist() #theoretical fragments for each charge
    for u in theo_fragment:
            theo_frag_M = (u * first_charge) - (h_mass * first_charge)
            precursor_matches['Actual fragment M'] = (precursor_matches['Fragment m/z'] * precursor_matches['New charge assignment']) - (h_mass * precursor_matches['New charge assignment'])
            precursor_matches['Theoretical fragment m/z'] = u
            precursor_matches['Theoretical fragment M'] = theo_frag_M
            precursor_matches['Peptide'] = b
            precursor_matches['Fragment error (Da)'] = abs(precursor_matches['Actual fragment M'] - theo_frag_M)
            precursor_matches_filter = precursor_matches.sort_values(by='Fragment error (Da)')
            precursor_matches_filter = precursor_matches_filter[precursor_matches_filter['Fragment error (Da)'] <= error_fragment]           
            if len(precursor_matches_filter) > 0:
                fragment_matches = fragment_matches.append(precursor_matches_filter)
fragment_matches = fragment_matches.drop_duplicates()
fragment_matches = fragment_matches.drop(columns=['empty','empty2'])

file_name = sample_name + '_' + str(first_charge) + '_zero_reassign_fragment_matches.csv'
file_out_path = working_directory + '\\' + file_name
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        fragment_matches.to_csv(filec,index=False)

