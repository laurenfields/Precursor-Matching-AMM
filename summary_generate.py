# -*- coding: utf-8 -*-
"""
Created on Tue May  3 12:04:08 2022

@author: lawashburn
"""

import os
import csv
import pandas as pd
import numpy as np
from datetime import datetime
now = datetime.now()

fragment_matches = pd.read_csv(r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Lauren_Results\all_fragments\A_all_fragments.csv")
final_dir = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220502_distribution\20220504\Lauren_Results\report_summary"
sample_name = 'A'

all_peptides = fragment_matches['Peptide'].values.tolist()

peptides = []
for a in all_peptides:
    if a not in peptides:
        peptides.append(a)
        
for b in peptides:
    peptide_fragments = fragment_matches[fragment_matches['Peptide'] == b]
    
    file_name = sample_name + '_' + b + '_fragments.csv'
    file_out_path = final_dir + '\\' + file_name
    with open(file_out_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            peptide_fragments.to_csv(filec,index=False)    
    
    scans_present = peptide_fragments['Scan #'].values.tolist()
    
    scans_filtered = []
    scan_summary_no = []
    scan_summary_instances = []
    
    for u in scans_present:
        if u not in scans_filtered:
            scans_filtered.append(u)
    
    for v in scans_filtered:
        scan_count = scans_present.count(v)
        scan_summary_no.append(v)
        scan_summary_instances.append(scan_count)

    scan_report = pd.DataFrame()
    scan_report['Scan #'] = scan_summary_no
    scan_report['# instances'] = scan_summary_instances

    file_name = sample_name + '_' + b + '_scan_summary.csv'
    file_out_path = final_dir + '\\' + file_name
    with open(file_out_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            scan_report.to_csv(filec,index=False)    