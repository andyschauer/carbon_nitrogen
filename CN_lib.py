#!/usr/bin/env python3

"""
Functions, constants, and other useful stuff for carbon and nitrogen isotopic data processing.

    Version 1.2 mod date 2024-06-09 => spacing for readability
    Version 2.0 mod date 2024-06-20 => changed how I refer to the standards: calibration_standards, etc; combined isolab_lib.py with shrekCN_lib.py and saved as CN_lib.py
"""


__author__ = "Andy Schauer"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024-06-20"
__version__ = "2.0"
__copyright__ = "Copyright 2024, Andy Schauer"
__license__ = "Apache 2.0"



# ---------- IMPORTS ----------
import json
import numpy as np
import os
import re



# ---------- LISTS ----------
calibration_standards = ['SA', 'GA1', "GA2"]
reference_material_list = calibration_standards[:] + ['Cocoa', 'NIST1547', 'USGS40', 'USGS41', 'MAL', 'DSM']
non_samples_list = reference_material_list[:] + ['zero', 'qtycal', 'blank', 'emptytin']

meta_headers = ['Amount', 'Analysis', 'Comment', 'Date', 'Identifier1', 'Identifier2', 'Information', 'Line', 'Method', 'Row', 'Time']

N_headers = ['Ampl28', 'Ampl29', 'AreaAll', 'Area28', 'Area29', 'BGD28', 'BGD29', 'Gasconfiguration', 'R15N14N', 'd15N14N', 'PeakNr', 'Start', 'Width']

C_headers = ['Ampl44', 'Ampl45', 'Ampl46', 'AreaAll', 'Area44', 'Area45', 'Area46', 'BGD44', 'BGD45', 'BGD46', 'Gasconfiguration',
             'PeakNr', 'R13C12C', 'd13C12C', 'Start', 'Width']

supp_headers = ['file', 'trust', 'notes', 'peak_center', 'pyversions', 'empty']

numlist = ['Amount', 'Analysis', 'Line','Row', 'trust', 'peak_center',
           'N_wg_PeakNr', 'N_wg_Start', 'N_wg_Width', 'N_wg_Ampl28', 'N_wg_Ampl29', 'N_wg_Area28', 'N_wg_Area29', 'N_wg_AreaAll', 
           'N_wg_BGD28', 'N_wg_BGD29', 'N_wg_R15N14N', 'N_wg_d15N14N', 
           'N_sam_PeakNr', 'N_sam_Start', 'N_sam_Width', 'N_sam_Ampl28', 'N_sam_Ampl29', 'N_sam_Area28', 'N_sam_Area29', 'N_sam_AreaAll',
           'N_sam_BGD28', 'N_sam_BGD29', 'N_sam_R15N14N', 'N_sam_d15N14N', 
           'C_sam_PeakNr', 'C_sam_Start', 'C_sam_Width', 'C_sam_Ampl44', 'C_sam_Ampl45', 'C_sam_Ampl46', 'C_sam_Area44', 'C_sam_Area45', 'C_sam_Area46', 'C_sam_AreaAll', 
           'C_sam_BGD44', 'C_sam_BGD45', 'C_sam_BGD46', 'C_sam_R13C12C', 'C_sam_d13C12C',           
           'C_wg_PeakNr', 'C_wg_Start', 'C_wg_Width', 'C_wg_Ampl44', 'C_wg_Ampl45', 'C_wg_Ampl46', 'C_wg_Area44', 'C_wg_Area45', 'C_wg_Area46', 'C_wg_AreaAll', 
           'C_wg_BGD44', 'C_wg_BGD45', 'C_wg_BGD46', 'C_wg_R13C12C', 'C_wg_d13C12C']


# ---------- CN-specific 'standards' ----------
blank = {'names': ['blank', 'void'],
         'material': None,
         'notes': 'no material dropped into EA'}
emptytin = {'names': ['empty_tin', 'Empty Tin'],
            'material': 'tin cups',
            'percentN': None,
            'percentC': None}
qtycal = {'names': ['qtycal_GA1', 'qtycal_ga1', 'qtycal.GA1'],
          'material': 'GA1',
          'index': [],
          'fractionN': 0.0952,
          'fractionC': 0.4082,
          'notes': 'material weighed across a range to calibrate peak area to quantity'}
zero = {'names': ['zero'],
        'material': 'reference gas peaks treated as unknowns'}



# ---------- EXPORT DATA FILE COLUMN HEADERS ----------
meta_data = {}
N_wg_data = {}
N_sam_data = {}
C_sam_data = {}
C_wg_data = {}
supp_data = {}

shrekCN_analysis_log_headers = []
shrekCN_analysis_log_headers.extend(meta_headers)
shrekCN_analysis_log_headers.extend([f"N_wg_{i}" for i in N_headers])
shrekCN_analysis_log_headers.extend([f"N_sam_{i}" for i in N_headers])
shrekCN_analysis_log_headers.extend([f"C_sam_{i}" for i in C_headers])
shrekCN_analysis_log_headers.extend([f"C_wg_{i}" for i in C_headers])
shrekCN_analysis_log_headers.extend(supp_headers)

data_to_write = []
data_to_write.extend([f"meta_data['{i}'][ii]" for i in meta_headers])
data_to_write.extend([f"N_wg_data['{i}'][ii]" for i in N_headers])
data_to_write.extend([f"N_sam_data['{i}'][ii]" for i in N_headers])
data_to_write.extend([f"C_sam_data['{i}'][ii]" for i in C_headers])
data_to_write.extend([f"C_wg_data['{i}'][ii]" for i in C_headers])
data_to_write.extend([f"supp_data['{i}'][ii]" for i in supp_headers])
data_to_write = str(data_to_write).replace("\"", "")





# ---------- FUNCTIONS ----------
def get_outliers(data, sigma):
    m = np.nanmean(data)
    s = np.nanstd(data) * sigma
    o = [i for i, e in enumerate(data) if e > m + s or e < m - s]
    return m, s, o



def get_path(instrument, desired_path):
    """Make your life easier with this section. These are the paths that seem to change depending on the computer we are working on."""
    path_file = os.path.join(os.getcwd(), f'{instrument}_path.txt')
    if os.path.isfile(path_file):
        # print(' :-) Using existing shrekCN path file for a warm and fuzzy experience. (-:')
        with open(path_file, 'r') as ppf:
            home, python_path, project_path, standards_path = ppf.readline().split(',')
            python_path = home + python_path
            project_path = home + project_path
            standards_path = home + standards_path

    else:
        python_path_check = False
        project_path_check = False
        standards_path_check = False
        print(' )-: Happy path file does not exist yet. :-(')
        print(" Let's make one... :-| ")
        while python_path_check is False:
            python_path = input(f'Enter the current path to the python scripts. Perhaps it is {os.getcwd()}. ')
            if os.path.isdir(python_path):
                python_path_check = True
                if python_path[-1] != '/':
                    python_path += '/'
            else:
                print(f'oops, try typing that in again (you typed {python_path}): ')

        while project_path_check is False:
            project_path = input('Enter the current path to your projects: ')
            if os.path.isdir(project_path):
                project_path_check = True
                if project_path[-1] != '/':
                    project_path += '/'
            else:
                print(f'oops, try typing that in again (you typed {project_path}): ')

        while standards_path_check is False:
            standards_path = input('Enter the current path and filename to your reference materials file: ')
            if os.path.isfile(standards_path):
                standards_path_check = True
            else:
                print(f'oops, try typing that in again (you typed {standards_path}): ')

        with open(path_file, 'w') as ppf:
            ppf.write(f'{python_path},{project_path},{standards_path}')

    if desired_path == "project":
        return project_path
    elif desired_path == "python":
        return python_path
    elif desired_path == "standards":
        return standards_path
    elif desired_path == "refmat_meta":
        with open(standards_path) as file:
            refmat = json.load(file)
            return f"{refmat['file_meta_data']['file']} - {refmat['file_meta_data']['modification_date']}"
    else:
        unknown_path = input('Unknown path, please enter the path to your project: ')
        return unknown_path



def make_file_list(directory, filetype):
    """Create and return a list of files contained within a directory
    of file type."""
    filelist = []
    initial_list = os.listdir(directory)
    for file in initial_list:
        if re.search(filetype, file):
            filelist.append(file)
    return filelist



def read_file(file_to_import, delim=None, header_row=1):
    """Read in a delimited text file containing a single header row
    followed by data and return those headers as a list and the data
    as a dictionary."""
    with open(file_to_import, 'r') as f:
        if header_row > 1:
            f.readline()
        headers = f.readline().split(delim)

        # remove unwanted characters from headers using a regular expression
        p = re.compile(r'[./\s()%]')  # list of characters to match
        for ii in range(len(headers)):
            m = p.findall(headers[ii])
            for em in m:
                headers[ii] = headers[ii].replace(em, '')

        data = {}
        for h in headers:
            data[h] = []

        for line in f:
            row = line.split(delim)
            if delim is None:
                if len(row) < len(headers):
                    row.append(0)
            # populate dictionary with all data in all rows
            for h, v in zip(headers, row):
                v = v.replace('\n', '').replace('1.#IO', '').replace('1.#INF000', '')
                if v == '':
                    v = None
                data[h].append(v)

    return headers, data
