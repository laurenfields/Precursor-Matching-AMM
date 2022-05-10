# -*- coding: utf-8 -*-
"""
Created on Tue May  3 11:03:17 2022

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
remaining_charges = fragment_charges[1:]

def get_file_names_with_strings(str_list):
    full_list = os.listdir(working_directory)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

#input info here for loop 1, charge1

all_fragment_matches = pd.DataFrame()
unavailable_spectra = pd.DataFrame()
all_used_fragments = pd.DataFrame()
for zz in remaining_charges:
    zs = str(zz)
    file_query = '_zero_reassign_fragment_matches.csv'
    fragment_list = (get_file_names_with_strings([file_query]))
    for a in fragment_list:
        fragment_list_path = working_directory + '\\' + a
        fragment_read = pd.read_csv(fragment_list_path)
        all_used_fragments = all_used_fragments.append(fragment_read)

    available_spectra = pd.merge(precursor_matches, all_used_fragments,how='left', left_on=['Fragment m/z', 'Scan #','Fragment resolution','Fragment ion charge',
                                                                                            'Fragment intensity','Actual precursor','Precursor actual charge','Precursor error (ppm)',
                                                                                            'Theoretical precursor','Precursor theoretical charge'],
                                 right_on=['Fragment m/z','Scan #','Fragment resolution','Fragment ion charge','Fragment intensity','Actual precursor',
                                           'Precursor actual charge','Precursor error (ppm)','Theoretical precursor','Precursor theoretical charge'])
    

    available_spectra_filtered = available_spectra[available_spectra['Fragment error (Da)'].isna()]                                                                                                                                              
    
    true_working_spectra = available_spectra_filtered
    true_working_spectra['New charge assignment'] = zz
    #true_working_spectra = true_working_spectra.rename(columns={"resolution_x": "resolution", "intensity_x": "intensity", "Precursor_x": "Precursor", "Precursor_Charge_x": "Precursor_Charge"})
    true_working_spectra = true_working_spectra.drop(['empty','empty2'], axis=1)
    true_working_spectra2 = true_working_spectra.drop_duplicates()
    
    fragment_types_dups = fragment_database['peptide'].values.tolist()

    fragment_types = []
    for a in fragment_types_dups:
        if a not in fragment_types:
            fragment_types.append(a)
    fragment_matches = pd.DataFrame()
    for b in fragment_types:
        fragment_database_filtered = fragment_database[fragment_database['peptide'] == b]   
        theo_fragment = fragment_database_filtered[zs].values.tolist() #theoretical fragments for each charge
        for u in theo_fragment:
                theo_frag_M = (u * zz) - (h_mass * zz)
                true_working_spectra2['Actual fragment M'] = (true_working_spectra2['Fragment m/z'] * true_working_spectra2['New charge assignment']) - (h_mass * true_working_spectra2['New charge assignment'])
                true_working_spectra2['Theoretical fragment m/z'] = u
                true_working_spectra2['Theoretical fragment M'] = theo_frag_M
                true_working_spectra2['Peptide'] = b
                true_working_spectra2['Fragment error (Da)'] = abs(true_working_spectra2['Actual fragment M'] - theo_frag_M)
                fragment_matches_filter = true_working_spectra2.sort_values(by='Fragment error (Da)')
                fragment_matches_filter = fragment_matches_filter[fragment_matches_filter['Fragment error (Da)'] <= error_fragment]           
                if len(fragment_matches_filter) > 0:
                    fragment_matches = fragment_matches.append(fragment_matches_filter)
    fragment_matches = fragment_matches.drop_duplicates()
    file_name = sample_name + '_' + zs + '_zero_reassign_fragment_matches.csv'
    file_out_path = working_directory + '\\' + file_name
    with open(file_out_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            fragment_matches.to_csv(filec,index=False)

