# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 11:54:58 2022

@author: lawashburn
"""

import pandas as pd
import csv
import os

zero_dir = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Lauren_Results\zero_charge_assign"
normal_df = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Lauren_Results\fragment_matching\A_fragment_matches.csv"
output_dir = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Lauren_Results\all_fragments"
sample_name = 'A'

def get_file_names_with_strings(str_list):
    full_list = os.listdir(zero_dir)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

file_query = '_zero_reassign_fragment_matches.csv'
fragment_list = (get_file_names_with_strings([file_query]))

all_fragments = pd.DataFrame()

for a in fragment_list:
    a_path = zero_dir + '\\' + a
    a_df = pd.read_csv(a_path)
    all_fragments = all_fragments.append(a_df)


normal_df = pd.read_csv(normal_df)
normal_df_restructure = pd.DataFrame()
normal_df_restructure['Fragment m/z'] = normal_df['Fragment m/z']
normal_df_restructure['Fragment resolution'] = normal_df['Fragment resolution']
normal_df_restructure['Fragment ion charge'] = normal_df['Fragment ion charge']
normal_df_restructure['Fragment intensity'] = normal_df['Fragment intensity']
normal_df_restructure['Actual precursor'] = normal_df['Actual precursor']
normal_df_restructure['Scan #'] = normal_df['Scan #']
normal_df_restructure['Precursor actual charge'] = normal_df['Precursor actual charge']
normal_df_restructure['Precursor error (ppm)'] = normal_df['Precursor error (ppm)']
normal_df_restructure['Theoretical precursor'] = normal_df['Theoretical precursor']
normal_df_restructure['Precursor theoretical charge'] = normal_df['Theoretical precursor charge']
normal_df_restructure['New charge assignment'] = normal_df['Fragment ion charge']
normal_df_restructure['Actual fragment M'] = normal_df['Actual fragment M']
normal_df_restructure['Theoretical fragment m/z'] = normal_df['Theoretical fragment m/z']
normal_df_restructure['Theoretical fragment M'] = normal_df['Theoretical fragment M']
normal_df_restructure['Peptide'] = normal_df['Peptide']
normal_df_restructure['Fragment error (Da)'] = normal_df['Fragment error (Da)']


all_fragments = all_fragments.append(normal_df_restructure)
all_fragments = all_fragments.drop_duplicates()
#all_fragments = all_fragments.iloc[:,:16]

file_name = sample_name + '_all_fragments.csv'
file_out_path = output_dir + '\\' + file_name
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        all_fragments.to_csv(filec,index=False)