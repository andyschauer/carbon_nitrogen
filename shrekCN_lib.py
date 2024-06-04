#!/usr/bin/env python3
"""
Random items needed by both shrekCN.py and shrekCNcalibrate.py
"""

__author__ = "Andy Schauer"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024-06-03"
__version__ = "1.1"
__copyright__ = "Copyright 2024, Andy Schauer"
__license__ = "Apache 2.0"


refmat_list = ['NIST1547', 'GA1', 'GA2', 'SA', 'MAL', 'DSM']
knowns_list = refmat_list[:] + ['zero', 'qtycal', 'blank', 'emptytin', 'USGS40', 'USGS41']

meta_headers = ['Amount', 'Analysis', 'Comment', 'Date', 'Identifier1', 'Identifier2', 'Information', 'Line', 'Method', 'Row', 'Time']

N_headers = ['AT15N14N', 'Ampl28', 'Ampl29', 'AreaAll', 'Area28', 'Area29', 'BGD28', 'BGD29', 'Gasconfiguration', 'R15N14N', 'R29N228N2', 'd15N14N',
             'd29N228N2', 'PeakNr', 'rArea28', 'rArea29', 'rR29N228N2', 'rd29N228N2', 'Start', 'Width']

C_headers = ['AT13C12C', 'Ampl44', 'Ampl45', 'Ampl46', 'AreaAll', 'Area44', 'Area45', 'Area46', 'BGD44', 'BGD45', 'BGD46', 'Gasconfiguration',
             'PeakNr', 'R13C12C', 'R45CO244CO2', 'R46CO244CO2', 'd13C12C', 'd45CO244CO2', 'd46CO244CO2',
             'rArea45', 'rR45CO244CO2', 'rR46CO244CO2', 'rd45CO244CO2', 'rd46CO244CO2', 'Start', 'Width']

supp_headers = ['file', 'trust', 'notes', 'peak_center', 'pyversions', 'empty']

numlist = ['Amount', 'Analysis', 'Line','Row', 'N_wg_AT15N14N', 'N_wg_Ampl28', 'N_wg_Ampl29', 'N_wg_AreaAll', 'N_wg_Area28', 'N_wg_Area29', 'N_wg_BGD28', 'N_wg_BGD29',
           'N_wg_R15N14N', 'N_wg_R29N228N2', 'N_wg_d15N14N', 'N_wg_d29N228N2', 'N_wg_PeakNr', 'N_wg_rArea28', 'N_wg_rArea29', 'N_wg_rR29N228N2', 'N_wg_rd29N228N2', 'N_wg_Start',
           'N_wg_Width', 'N_sam_AT15N14N', 'N_sam_Ampl28', 'N_sam_Ampl29', 'N_sam_AreaAll', 'N_sam_Area28', 'N_sam_Area29', 'N_sam_BGD28', 'N_sam_BGD29', 'N_sam_R15N14N',
           'N_sam_R29N228N2', 'N_sam_d15N14N', 'N_sam_d29N228N2', 'N_sam_PeakNr', 'N_sam_rArea28', 'N_sam_rArea29', 'N_sam_rR29N228N2', 'N_sam_rd29N228N2', 'N_sam_Start',
           'N_sam_Width', 'C_sam_AT13C12C', 'C_sam_Ampl44', 'C_sam_Ampl45', 'C_sam_Ampl46', 'C_sam_AreaAll', 'C_sam_Area44', 'C_sam_Area45', 'C_sam_Area46', 'C_sam_BGD44',
           'C_sam_BGD45', 'C_sam_BGD46', 'C_sam_PeakNr', 'C_sam_R13C12C', 'C_sam_R45CO244CO2', 'C_sam_R46CO244CO2', 'C_sam_d13C12C', 'C_sam_d45CO244CO2',
           'C_sam_d46CO244CO2', 'C_sam_rArea45', 'C_sam_rR45CO244CO2', 'C_sam_rR46CO244CO2', 'C_sam_rd45CO244CO2', 'C_sam_rd46CO244CO2', 'C_sam_Start', 'C_sam_Width',
           'C_wg_AT13C12C', 'C_wg_Ampl44', 'C_wg_Ampl45', 'C_wg_Ampl46', 'C_wg_AreaAll', 'C_wg_Area44', 'C_wg_Area45', 'C_wg_Area46', 'C_wg_BGD44', 'C_wg_BGD45', 'C_wg_BGD46',
           'C_wg_PeakNr', 'C_wg_R13C12C', 'C_wg_R45CO244CO2', 'C_wg_R46CO244CO2', 'C_wg_d13C12C', 'C_wg_d45CO244CO2', 'C_wg_d46CO244CO2', 'C_wg_rArea45', 'C_wg_rR45CO244CO2',
           'C_wg_rR46CO244CO2', 'C_wg_rd45CO244CO2', 'C_wg_rd46CO244CO2', 'C_wg_Start', 'C_wg_Width', 'trust', 'peak_center']

qtycal = {'names': ['qtycal_GA1', 'qtycal_ga1', 'qtycal.GA1'],
          'material': 'GA1',
          'index': [],
          'fractionN': 0.0952,
          'fractionC': 0.4082,
          'notes': 'material weighed across a range to calibrate peak area to quantity'}

blank = {'names': ['blank'],
         'material': None,
         'notes': 'no material dropped into EA'}

emptytin = {'names': ['empty_tin', 'Empty Tin'],
            'material': 'tin cups',
            'percentN': None,
            'percentC': None}

zero = {'names': ['zero'],
        'material': 'reference gas peaks treated as unknowns'}



meta_data = {}
N_wg_data = {}
N_sam_data = {}
C_sam_data = {}
C_wg_data = {}
supp_data = {}

# prepare export file column headers
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
