#!/usr/bin/env python3
"""
This is the calibration script for carbon and nitrogen stable isotope elemental analysis data that came
from a ThermoFinnigan 253 with Eurovector Elemental Analyzer. It summarizes the run(s) for quality control
and then calibrates d13C and d15N to the VPDB and AirN2 scales, respectively. Ideally, CN.py was run
prior to the present script.

Version 0.1 mod date 2021-03-03 => created
Version 0.2 mod date 2021-10-08 => updating nomenclature and adding figures
Version 0.3 mod date 2022-01-18 => slacker didn't make any notes
Version 0.4 mod date 2022-02-10 => trying to make this useful.
Version 0.5 mod date 2022-03-23 => saved from shrekCNcalibrate.py
Version 0.6 mod date 2022-05-19 => adding more to figures and generally trying to move this along
Version 0.7 mod date 2022-10-04 => combined shrekCNlog with shrekCNcalibrate. use "--verbose" argument when calling to get more log type info. Without, it will output a more sample data focused report (although, this report is not yet made)
Version 1.0 mod date 2024-04-10 => trying to get this finished up to a version 1 level and upload to github
Version 2.0 mod date 2024-06-03 => removed all specific standard references in favor of the automatically chosen type; code now separates out individual runs to be calibrated individually; also creates a summary file but this file still requires a fair bit of manual work to get it presentable
Version 2.1 mod date 2024-06-22 => made instrument a variable, added unify argument, started updating for bokeh deprecations, renamed to be CN_calibrate.py, touched up figures a bit
Version 2.2 mod date 2024-06-23 => mistake in Nqty and Cqty calculation found, needed to use the blank corrected peak areas, fixed now
Version 2.3 mod date 2024-07-13 => removed std_1, std_2, std_3 picking and now ask the user to enter n-1 reference materials to correct to, all others are used as qaqc
"""

__author__ = "Andy Schauer"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024-07-13"
__version__ = "2.3"
__copyright__ = "Copyright 2025, Andy Schauer"
__license__ = "Apache 2.0"
__acknowledgements__ = "Shrek"



# ---------- imports ---------- 
import argparse
from bokeh.io import output_file, show
from bokeh.layouts import row, column, grid
from bokeh.plotting import figure
from bokeh.palettes import Category20
from bokeh.resources import CDN, INLINE
from bokeh.embed import file_html
from CN_lib import *
import csv
import datetime as dt
import dateutil.parser
import json
import matplotlib.pyplot as pplt
import numpy as np
import os
import re
import shutil
import sys
import time
import webbrowser



# ---------- FUNCTIONS ---------- 
def add_calculation_note(note):
    calculation_notes.append(note)



# ---------- ARGUMENTS ----------
argument_string = '' 
parser = argparse.ArgumentParser()
parser.add_argument("--verbose", help="Include exhaustive diagnostic information and figures in report.", action="store_true")
parser.add_argument("--unify", help="Calibrate the entire log file as a single unified run.", action="store_true")
args = parser.parse_args()
if args.verbose:
    verbose = True
    argument_string += 'verbose, '
else:
    verbose = False
if args.unify:
    unify = True
    argument_string += 'unified calibration, '
else:
    unify = False
print(f'\nArguments: {argument_string}')


# ---------- LOAD CONFIGURATION ----------
with open("CN_config.json", 'r') as f:
    config = json.load(f)


# ---------- SETUP ---------- 
version = os.path.basename(__file__) + ' - ' + time.ctime(os.path.getctime(__file__))

home_directory = config["local_directories"]["home"]
python_directory = f"{home_directory}{config['local_directories']['python']}"
method_directory = f"{home_directory}{config['local_directories']['method_data_directory']}"
reference_materials_file = f"{home_directory}{config['local_directories']['standards']}"
new_data_directory = 'rawdata_new'
archive_data_directory = 'rawdata_archive'
junk_data_directory = 'rawdata_junk'

python_scripts = {'CN_lib.py': '', 'CN.py': '', 'CN_calibrate.py': ''}
python_scripts = {key: (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(f'{python_directory}{key}')))) for key, value in python_scripts.items()}

CN_log_file_list = make_file_list(method_directory, '_analysis_log.csv')

print('\nWhat analysis log file to you wish to process?\n')

[print(f'    {i}') for i in CN_log_file_list]
identified_file = 0
while identified_file == 0:
    CN_log_file_search = input('\nEnter a project analysis log file from above: ')
    isfile = [CN_log_file_search[0: len(CN_log_file_search)] in x for x in CN_log_file_list]
    if len(np.where(isfile)[0]) == 1:
        identified_file = 1
        log_file_name = CN_log_file_list[np.where(isfile)[0][0]]
        print(f'    Processing CN log file {log_file_name}...')
    else:
        print('\n** More than one file found. **\n')

if os.path.isdir(method_directory) is False:
    print('Method directory does not exist...exiting....')
    sys.exit()

project_directory = os.path.join(method_directory, f'{log_file_name[:-17]}/')
if os.path.exists(project_directory)==False:
    os.mkdir(project_directory)

archive_directory = os.path.join(project_directory, 'archive/')
if os.path.exists(archive_directory)==False:
    os.mkdir(archive_directory)


summary_data_filename = f'{log_file_name[:-17]}_summary.csv'
summary_data_file = os.path.join(project_directory, summary_data_filename)


# ---------- load reference material information ----------
ref_mat = {}

with open(reference_materials_file, 'r') as f:
    reference_materials = json.load(f)

non_samples_list = []
ref_mat['keys'] = list(reference_materials['organics'].keys())
for i in ref_mat['keys']:
    globals()[i] = reference_materials['organics'][i]
    globals()[i]['index'] = np.empty(0, dtype="int16")
    non_samples_list.append(i)

for i in config['corrective_measurements'].keys():
    globals()[i] = config['corrective_measurements'][i]
    globals()[i]['index'] = np.empty(0, dtype="int16")
    non_samples_list.append(i)


# ---------- get data ---------- 
print('Reading in data...')
headers, data = read_file(os.path.join(method_directory, log_file_name), ',')
entire_data_set = data.copy()
strlist = set(headers) - set(numlist)



# ---------- Get indices of standards, samples, and anything else we may need to know the location in the array ----------

if unify:
    runs = {'files': log_file_name}
    runs['indices'] = [np.where(data['Analysis'])[0]]
    runs['names'] = [dateutil.parser.parse(data['Date'][i[0]]).strftime("%Y%m%d") for i in runs['indices']]
else:
    runs = {'files': list(set(data['file']))}
    runs['indices'] = [np.where(np.asarray(data['file'])==i)[0] for i in runs['files']]
    runs['names'] = [dateutil.parser.parse(data['Date'][i[0]]).strftime("%Y%m%d") for i in runs['indices']]



for current_run_index, current_run_name in enumerate(runs['names']):

    print(f"Run = {runs['names'][current_run_index]}")

    current_data_set = data.copy()

    # remove trust 0 analyses from data
    trust0_indices = [i for i, e in enumerate(data['trust']) if int(e) == 0]
    trust1_indices = [i for i, e in enumerate(data['trust']) if int(e) == 1]
    current_indices = runs['indices'][current_run_index]
    trust1_currrun_indices = np.intersect1d(current_indices, trust1_indices)

    for header in headers[:-1]:
        current_data_set[header] = [current_data_set[header][index] for index in trust1_currrun_indices]

    for i in numlist:
        globals()[i] = np.asarray(current_data_set[i], dtype=float)

    for i in strlist:
        globals()[i] = np.asarray(current_data_set[i])

    calculation_notes = []


    # ------------------- Indices --------------------------
    all_indices = [i for i, _ in enumerate(Analysis)]

    non_samples_indices = []
    for i in non_samples_list:
        eval(i)['index'] = [j for j, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in eval(i)['names'])]
        non_samples_indices.extend(eval(i)['index'])

    # included_isotope_standards = list(set([i for i in Identifier1 if i in calibration_standards]))
    sample_indices = list(set(all_indices) - set(non_samples_indices))




    ref_mat['indices'] = []
    for i in list(ref_mat['keys']):
        eval(i)['index'] = [j for j, e in enumerate(Identifier1) if str(e).lower() in (name.lower() for name in eval(i)['names'])]
        ref_mat['indices'].extend(eval(i)['index'])

    ref_mat['id1_set'] = list(set(Identifier1[ref_mat['indices']]))
    print(f'\nThese reference materials were included in your run / set:')
    [print(f'    {i}') for i in ref_mat['id1_set']]
    print('Choose reference materials from the list above you wish to normalize to.')
    ref_mat['chosen'] = input(f"Enter at least 2 and at most {len(ref_mat['id1_set']) - 1} (e.g., std1, std2, std3): ")
    if ',' in ref_mat['chosen']:
        ref_mat['chosen'] = ref_mat['chosen'].split(',')
    else:
        ref_mat['chosen'] = ref_mat['chosen'].split()
    ref_mat['chosen'] = [i.strip().upper() for i in ref_mat['chosen']]
    ref_mat['qaqc'] = list(np.setdiff1d(ref_mat['id1_set'], ref_mat['chosen']))



    # ---------- Blanks and blank corrections ----------
    #    The blank correction is a simple two-source mixing model as follows:
    #        (d15N_measured * Size_measured) = (d15N_sample * Size_sample) - (d15N_blank * Size_blank)
    if blank['index']:
        if np.any(~np.isnan(N_sam_AreaAll[blank['index']])):
            blank['size_Vs'] = np.nanmean(N_sam_AreaAll[blank['index']])
        else:
            blank['size_Vs'] = np.nan
        if np.any(~np.isnan(N_sam_d15N14N[blank['index']])):
            blank['d15N'] = np.nanmean(N_sam_d15N14N[blank['index']])
        else:
            blank['d15N'] = np.nan
        if np.all(~np.isnan([blank['size_Vs'], blank['d15N']])):
            N_sam_AreaAll_blank_corr = N_sam_AreaAll - blank['size_Vs']
            d15N_blank_corr = ((N_sam_d15N14N * N_sam_AreaAll) - (blank['d15N'] * blank['size_Vs'])) / (N_sam_AreaAll - blank['size_Vs'])
            add_calculation_note("nitrogen blank correction applied")
        else:
            N_sam_AreaAll_blank_corr = N_sam_AreaAll
            d15N_blank_corr = N_sam_d15N14N
            add_calculation_note("nitrogen blank correction NOT applied")
    else:
        N_sam_AreaAll_blank_corr = N_sam_AreaAll
        d15N_blank_corr = N_sam_d15N14N
        add_calculation_note("nitrogen blank correction NOT applied")


    # ToDo: add carbon blank correction which would come from tins and only noticable during undiluted carbon runs
    d13C_blank_corr = C_sam_d13C12C[:]
    C_sam_AreaAll_blank_corr = C_sam_AreaAll



    # ---------- N and C Quantity ----------
    qtycal['Nqty'] = Amount[qtycal['index']] * qtycal['fractionN']
    qtycal['Nfit'] = np.polyfit(qtycal['Nqty'], N_sam_AreaAll_blank_corr[qtycal['index']], 1)

    qtycal['Cqty'] = Amount[qtycal['index']] * qtycal['fractionC']
    qtycal['Cfit'] = np.polyfit(qtycal['Cqty'], C_sam_AreaAll_blank_corr[qtycal['index']], 1)

    Nqty = (N_sam_AreaAll_blank_corr - qtycal['Nfit'][1]) / qtycal['Nfit'][0]
    Cqty = (C_sam_AreaAll_blank_corr - qtycal['Cfit'][1]) / qtycal['Cfit'][0]

    with np.errstate(divide='ignore', invalid='ignore'):
        PercentN = np.nan_to_num(Nqty / Amount * 100)
        PercentC = np.nan_to_num(Cqty / Amount * 100)

    add_calculation_note("nitrogen and carbon quantities calculated from peak area using quantity calibration standards (qtycal) as knowns")



    # ---------- Isotope Calibration Setup ----------

    for i in ref_mat['chosen']:
        eval(i)['purpose'] = 'd15N calibration; d13C calibration'
    for i in ref_mat['qaqc']:
        eval(i)['purpose'] = 'd15N quality control; d13C quality control'


    # ---------- Quantity Residuals ----------
    qtycal['Nresidual'] = N_sam_AreaAll_blank_corr[qtycal['index']] - (qtycal['Nfit'][0] * qtycal['Nqty'] + qtycal['Nfit'][1])
    qtycal['Cresidual'] = C_sam_AreaAll_blank_corr[qtycal['index']] - (qtycal['Cfit'][0] * qtycal['Cqty'] + qtycal['Cfit'][1])

    Nqty_residual = np.asarray([item for subarray in [(Nqty[eval(i)['index']] - Amount[eval(i)['index']] * eval(i)['fractionN']) * 1000 for i in ref_mat['chosen']] for item in subarray])
    Cqty_residual = np.asarray([item for subarray in [(Cqty[eval(i)['index']] - Amount[eval(i)['index']] * eval(i)['fractionC']) * 1000 for i in ref_mat['chosen']] for item in subarray])


    # ----------------- good data index - gdi ---------------------------

    gdi = np.asarray(all_indices)
    # gdi = np.where(Nqty>0.050)[0] 

    for i in ref_mat['chosen']:
        eval(i)['gdi'] = np.intersect1d(eval(i)['index'], gdi)

    for i in ref_mat['chosen']:
        eval(i)['d15N_residual'] = N_sam_d15N14N[eval(i)['index']] - np.nanmean(N_sam_d15N14N[np.intersect1d(eval(i)['index'], gdi)])
        eval(i)['d13C_residual'] = C_sam_d13C12C[eval(i)['index']] - np.nanmean(C_sam_d13C12C[eval(i)['index']])

    ref_mat['d15N_residual_std'] = np.std([item for subarray in [eval(i)['d15N_residual'] for i in ref_mat['chosen']] for item in subarray])
    ref_mat['d13C_residual_std'] = np.std([item for subarray in [eval(i)['d13C_residual'] for i in ref_mat['chosen']] for item in subarray])


    # ---------- Isotope Drift Calculation ----------
    Ndrift_fit = np.polyfit(np.asarray([item for subarray in [Analysis[eval(i)['index']] for i in ref_mat['chosen']] for item in subarray]),
                            np.asarray([item for subarray in [eval(i)['d15N_residual'] for i in ref_mat['chosen']] for item in subarray]),
                            1)
    Ndrift_corrfac = Ndrift_fit[0] * Analysis + Ndrift_fit[1]
    d15N_drift_corr = d15N_blank_corr - Ndrift_corrfac

    Cdrift_fit = np.polyfit(np.asarray([item for subarray in [Analysis[eval(i)['index']] for i in ref_mat['chosen']] for item in subarray]),
                            np.asarray([item for subarray in [eval(i)['d13C_residual'] for i in ref_mat['chosen']] for item in subarray]),
                            1)
    Cdrift_corrfac = Cdrift_fit[0] * Analysis + Cdrift_fit[1]
    d13C_drift_corr = d13C_blank_corr - Cdrift_corrfac

    add_calculation_note("d15N and d13C were corrected for drift")



    # ---------- Isotope Calibration ----------

    ref_mat['d13Cacc'] = [eval(i)['d13C'] for i in ref_mat['chosen']]
    ref_mat['d13Cmeas'] = [np.nanmean(d13C_drift_corr[eval(i)['index']]) for i in ref_mat['chosen']]
    ref_mat['d13C_fit'] = np.polyfit(ref_mat['d13Cmeas'], ref_mat['d13Cacc'], 1)
    d13C_VPDB = np.asarray(ref_mat['d13C_fit'][0] * d13C_drift_corr + ref_mat['d13C_fit'][1])

    ref_mat['d15Nacc'] = [eval(i)['d15N'] for i in ref_mat['chosen']]
    ref_mat['d15Nmeas'] = [np.nanmean(d15N_drift_corr[eval(i)['index']]) for i in ref_mat['chosen']]
    ref_mat['d15N_fit'] = np.polyfit(ref_mat['d15Nmeas'], ref_mat['d15Nacc'], 1)
    d15N_AirN2 = np.asarray(ref_mat['d15N_fit'][0] * d15N_drift_corr + ref_mat['d15N_fit'][1])

    add_calculation_note("d15N and d13C were normalized to AirN2 and VPDB, respectively")



    # ---------- Post-normalization residual calculation ----------
    for i in ref_mat['chosen']:
        eval(i)['d15N_AirN2_residual'] = d15N_AirN2[eval(i)['index']] - np.nanmean(d15N_AirN2[np.intersect1d(eval(i)['index'], gdi)])
        eval(i)['d13C_VPDB_residual'] = d13C_VPDB[eval(i)['index']] - np.nanmean(d13C_VPDB[eval(i)['index']])

    ref_mat['d15N_AirN2_residual_std'] = np.std([item for subarray in [eval(i)['d15N_AirN2_residual'] for i in ref_mat['chosen']] for item in subarray])
    ref_mat['d13C_VPDB_residual_std'] = np.std([item for subarray in [eval(i)['d13C_VPDB_residual'] for i in ref_mat['chosen']] for item in subarray])



    # populate current run dictionary 
    runs[current_run_name] = {'non_samples_indices': non_samples_indices,
                              'sample_indices': sample_indices,
                              'Identifier1': Identifier1,
                              'Date': Date,
                              'Analysis': Analysis,
                              'Amount': Amount,
                              'Nqty': Nqty,
                              'd15N_AirN2': d15N_AirN2,
                              'Cqty': Cqty,
                              'd13C_VPDB': d13C_VPDB}


    # ---------- FIGURES ---------- 
    print('Making figures...')

    figures = {}
    fig_n = 1
    font_size = "24pt"
    symbols = ['circle', 'square', 'triangle', 'diamond', 'inverted_triangle', 'asterisk', 'cross', 'x', 'hex', 'y']
    colors = Category20[20]

    for i,j in enumerate(ref_mat['id1_set']):
        eval(j)['symbol_color'] = [colors[i]]


    if verbose:

        figures[fig_n] = {}
        figures[fig_n]['cap'] = f"""Figure {fig_n}. The sample peak start time in chronological order can show how the chromatography is changing
                                    with time. Ideally, the sample peak start time does not change. Fluctuations are expected depending on the size
                                    of the peak but directional drift in the start time generally indicates something is changing that we don't want
                                    to change. For example, the GC column is getting wet, the magnesium perchlorate is become saturated and clogged,
                                    or the temperature in the room or of the GC column is changing."""
        figures[fig_n]['fig'] = figure(title="Nitrogen- and carbon-sample-peak start time vs analysis number", width=1200, height=600, background_fill_color="#fafafa")
        figures[fig_n]['fig'].scatter(Analysis, N_sam_Start, legend_label="Nitrogen", marker='triangle', size=8, color="blue", alpha=0.8)
        figures[fig_n]['fig'].scatter(Analysis, C_sam_Start, legend_label="Carbon", marker='square', size=8, color="black", alpha=0.8)
        figures[fig_n]['fig'].yaxis.axis_label = "Sample peak start time (seconds)"
        figures[fig_n]['fig'].xaxis.axis_label = "Analysis number"
        figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].title.text_font_size = font_size
        figures[fig_n]['fig'].legend.label_text_font_size = font_size

        fig_n += 1

        figures[fig_n] = {}
        figures[fig_n]['cap'] = f"""Figure {fig_n}. The separation between the end of the sample nitrogen peak and the start of the carbon sample peak is also
                                    indicative of the chromatography conditions. Generally, if the peaks are moving closer to one another, the GC column
                                    is getting wet and needs to be baked out."""
        figures[fig_n]['fig'] = figure(title="Nitrogen and carbon sample peak separation vs analysis number", width=1200, height=600, background_fill_color="#fafafa")
        figures[fig_n]['fig'].scatter(Analysis, C_sam_Start-(N_sam_Start+N_sam_Width), marker='circle', size=8, color="red", alpha=0.8)
        figures[fig_n]['fig'].yaxis.axis_label = "Sample peak separation (seconds)"
        figures[fig_n]['fig'].xaxis.axis_label = "Analysis number"
        figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].title.text_font_size = font_size

        fig_n += 1

        figures[fig_n] = {}
        figures[fig_n]['cap'] = f"""Figure {fig_n}. The sample peak width, as with the above figure, if changing in a directional manner, may be
                                    indicative of a problem with the chromatography."""
        figures[fig_n]['fig'] = figure(title="Nitrogen- and carbon-sample-peak width vs analysis number", width=1200, height=600, background_fill_color="#fafafa")
        figures[fig_n]['fig'].scatter(Analysis, N_sam_Width, legend_label="Nitrogen", marker='triangle', size=8, color="blue", alpha=0.8)
        figures[fig_n]['fig'].scatter(Analysis, C_sam_Width, legend_label="Carbon", marker='square', size=8, color="black", alpha=0.8)
        figures[fig_n]['fig'].yaxis.axis_label = "Sample peak width (seconds)"
        figures[fig_n]['fig'].xaxis.axis_label = "Analysis number"
        figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].title.text_font_size = font_size
        figures[fig_n]['fig'].legend.label_text_font_size = font_size

        fig_n += 1

        figures[fig_n] = {}
        figures[fig_n]['cap'] = f"""Figure {fig_n}. The peak height can help us observe several bits. One, the working gas peak height shows 
                                    mass spectrometer source sensitivity consistency, assuming the working gas pressure is not changing. The
                                    sample peak height gives you a glimps at the size of your samples. If you see any at the maximum value, 
                                    you know you have weighed too much material."""
        figures[fig_n]['fig'] = figure(title="Nitrogen- and carbon-sample-peak height vs analysis number", width=1200, height=600, background_fill_color="#fafafa")
        figures[fig_n]['fig'].scatter(Analysis, N_wg_Ampl28, legend_label="Nitrogen Working Gas", marker='triangle', size=5, line_color='blue', alpha=0.8)
        figures[fig_n]['fig'].scatter(Analysis, N_sam_Ampl28, legend_label="Nitrogen Sample", marker='triangle', size=8, color='blue', alpha=0.8)
        figures[fig_n]['fig'].scatter(Analysis, C_sam_Ampl44, legend_label="Carbon Sample", marker='square', size=8, color='black', alpha=0.8)
        figures[fig_n]['fig'].scatter(Analysis, C_wg_Ampl44, legend_label="Carbon Working Gas", marker='square', size=8, line_color='black', alpha=0.8)
        figures[fig_n]['fig'].yaxis.axis_label = 'Peak amplitude (mV)'
        figures[fig_n]['fig'].xaxis.axis_label = 'Analysis number'
        figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].title.text_font_size = font_size
        figures[fig_n]['fig'].legend.label_text_font_size = font_size

        fig_n += 1


    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. Nitrogen amount for each sample is calculated from the N2 peak area and from a linear equation
                                generated as a least squares fit from peak area and measured amount of the quantity calibration standards (qtycal).
                                The nitrogen amount of the other standards is also considered known and plotted here as such. We are
                                assuming, however, that we do not know the nitrogen amount for the samples."""
    figures[fig_n]['fig'] = figure(title="Peak Area vs Nitrogen quantity", width=1200, height=600, background_fill_color="#fafafa")
    figures[fig_n]['fig'].scatter(qtycal['Nqty'], N_sam_AreaAll_blank_corr[qtycal['index']], legend_label="qtycal", marker='circle', size=12, fill_color='yellow', line_color='black', alpha=0.8)
    for i in ref_mat['id1_set']:
        figures[fig_n]['fig'].scatter(Amount[eval(i)['index']] * eval(i)['fractionN'], N_sam_AreaAll_blank_corr[eval(i)['index']], legend_label=eval(i)['names'][0], marker='circle', size=8, color=eval(i)['symbol_color'][0], line_color='black', alpha=0.8)
    figures[fig_n]['fig'].scatter(Nqty[sample_indices], N_sam_AreaAll_blank_corr[sample_indices], legend_label="samples", marker='triangle', size=6, color='black', alpha=0.8)
    figures[fig_n]['fig'].line(([np.nanmin(N_sam_AreaAll_blank_corr), np.nanmax(N_sam_AreaAll_blank_corr)] - qtycal['Nfit'][1]) / qtycal['Nfit'][0], [np.nanmin(N_sam_AreaAll_blank_corr), np.nanmax(N_sam_AreaAll_blank_corr)], line_width=3, color='black')
    figures[fig_n]['fig'].legend.location = 'top_left'
    figures[fig_n]['fig'].yaxis.axis_label = 'Peak Area (Vs)'
    figures[fig_n]['fig'].xaxis.axis_label = 'Nitrogen Quantity (mg)'
    figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].title.text_font_size = font_size
    figures[fig_n]['fig'].legend.label_text_font_size = font_size

    fig_n += 1

    if verbose:

        figures[fig_n] = {}
        figures[fig_n]['cap'] = f"""Figure {fig_n}. Here we are able to see how the standards vary around the least squares line of best fit from above. The quality
                                    of the data used in the fit may be assessed by the residual standard deviation (<strong>2-&sigma;={np.round(np.std(qtycal['Nresidual'])*2, 1)} &micro;g</strong>). The 
                                    isotope reference materials have a residual <strong>2-&sigma;={np.round(np.std(Nqty_residual)*2, 1)} &micro;g</strong>."""
        figures[fig_n]['fig'] = figure(title="Nitrogen quantity residual vs Analysis number", width=1200, height=600, background_fill_color="#fafafa")
        figures[fig_n]['fig'].scatter(Analysis[qtycal['index']], qtycal['Nresidual'], legend_label="qtycal", marker='circle', size=12, fill_color='yellow', line_color='black', alpha=0.8)
        for i in ref_mat['id1_set']:
            figures[fig_n]['fig'].scatter(Analysis[eval(i)['index']], (Nqty[eval(i)['index']] - Amount[eval(i)['index']] * eval(i)['fractionN']) * 1000, legend_label=eval(i)['names'][0], marker='circle', size=8, color=eval(i)['symbol_color'][0], line_color='black', alpha=0.8)    
        figures[fig_n]['fig'].line([np.min(Analysis), np.max(Analysis)], [0, 0], line_width=3, color='black')
        figures[fig_n]['fig'].legend.location = 'top_left'
        figures[fig_n]['fig'].yaxis.axis_label = 'Nitrogen quantity residual (&micro;g)'
        figures[fig_n]['fig'].xaxis.axis_label = 'Analysis number'
        figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].title.text_font_size = font_size
        figures[fig_n]['fig'].legend.label_text_font_size = font_size

        fig_n += 1


    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. Carbon amount for each sample is calculated from the CO2 peak area and from a linear equation
                                generated as a least squares fit from peak area and measured amount of the quantity calibration standards (qtycal).
                                The carbon amount of the other standards is also considered known and plotted here as such. We are
                                assuming, however, that we do not know the carbon amount for the samples."""
    figures[fig_n]['fig'] = figure(title="Peak Area vs Carbon quantity", width=1200, height=600, background_fill_color="#fafafa")
    figures[fig_n]['fig'].scatter(qtycal['Cqty'], C_sam_AreaAll[qtycal['index']], legend_label="qtycal", marker='circle', size=12, fill_color='yellow', line_color='black', alpha=0.8)
    for i in ref_mat['id1_set']:
        figures[fig_n]['fig'].scatter(Amount[eval(i)['index']] * eval(i)['fractionC'], C_sam_AreaAll[eval(i)['index']], legend_label=eval(i)['names'][0], marker='circle', size=8, color=eval(i)['symbol_color'][0], line_color='black', alpha=0.8)
    figures[fig_n]['fig'].scatter(Cqty[sample_indices], C_sam_AreaAll[sample_indices], legend_label="samples", marker='triangle', size=6, color='black', alpha=0.8)
    figures[fig_n]['fig'].line(([np.nanmin(C_sam_AreaAll), np.nanmax(C_sam_AreaAll)] - qtycal['Cfit'][1]) / qtycal['Cfit'][0], [np.nanmin(C_sam_AreaAll), np.nanmax(C_sam_AreaAll)], line_width=3, color='black')
    figures[fig_n]['fig'].legend.location = 'top_left'
    figures[fig_n]['fig'].yaxis.axis_label = 'Peak Area (Vs)'
    figures[fig_n]['fig'].xaxis.axis_label = 'Carbon Quantity (mg)'
    figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].title.text_font_size = font_size
    figures[fig_n]['fig'].legend.label_text_font_size = font_size

    fig_n += 1


    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}. Nitrogen isotope composition versus carbon isotope composition of all reference materials and samples. This
                                figure allows you to see where the reference materials are in relation to your samples and assess the appropriateness of those
                                reference materials to your samples."""
    figures[fig_n]['fig'] = figure(title="d15N vs d13C", width=1200, height=600, background_fill_color="#fafafa")
    for i in ref_mat['id1_set']:
        figures[fig_n]['fig'].scatter(d15N_AirN2[eval(i)['index']], d13C_VPDB[eval(i)['index']], legend_label=eval(i)['names'][0], marker='circle', size=8, color=eval(i)['symbol_color'][0], line_color='black', alpha=0.8)
    figures[fig_n]['fig'].scatter(d15N_AirN2[sample_indices], d13C_VPDB[sample_indices], legend_label="samples", marker='triangle', size=6, color='black', alpha=0.8)
    figures[fig_n]['fig'].legend.location = 'top_left'
    figures[fig_n]['fig'].yaxis.axis_label = 'd15N vs Air-N2 (permil)'
    figures[fig_n]['fig'].xaxis.axis_label = 'd13C vs VPDB (permil)'
    figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].title.text_font_size = font_size
    figures[fig_n]['fig'].legend.label_text_font_size = font_size

    fig_n += 1

    if verbose:

        figures[fig_n] = {}
        figures[fig_n]['cap'] = f"""Figure {fig_n}. Here again we are able to see how the standards vary around the least squares line of best fit from above. The quality
                                    of the data used in the fit may be assessed by the residual standard deviation (<strong>2-&sigma;={np.round(np.std(qtycal['Cresidual'])*2, 1)} &micro;g</strong>). The 
                                    isotope reference materials have a residual <strong>2-&sigma;={np.round(np.std(Cqty_residual)*2, 1)} &micro;g</strong>."""
        figures[fig_n]['fig'] = figure(title="Carbon quantity residual vs Analysis number", width=1200, height=600, background_fill_color="#fafafa")
        figures[fig_n]['fig'].scatter(Analysis[qtycal['index']], qtycal['Cresidual'], legend_label="qtycal", marker='circle', size=12, fill_color='yellow', line_color='black', alpha=0.8)
        for i in ref_mat['id1_set']:
            figures[fig_n]['fig'].scatter(Analysis[eval(i)['index']], (Cqty[eval(i)['index']] - Amount[eval(i)['index']] * eval(i)['fractionC']) * 1000, legend_label=eval(i)['names'][0], marker='circle', size=8, color=eval(i)['symbol_color'][0], line_color='black', alpha=0.8)
        figures[fig_n]['fig'].line([np.min(Analysis), np.max(Analysis)], [0, 0], line_width=3, color='black')
        figures[fig_n]['fig'].legend.location = 'top_left'
        figures[fig_n]['fig'].yaxis.axis_label = 'Carbon quantity residual (&micro;g)'
        figures[fig_n]['fig'].xaxis.axis_label = 'Analysis number'
        figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].title.text_font_size = font_size
        figures[fig_n]['fig'].legend.label_text_font_size = font_size

        fig_n += 1

        figures[fig_n] = {}
        figures[fig_n]['cap'] = f"""Figure {fig_n}."""
        figures[fig_n]['fig'] = figure(title="d15N residual", width=1200, height=600, background_fill_color="#fafafa")
        for i in ref_mat['chosen']:
            figures[fig_n]['fig'].scatter(Analysis[eval(i)['index']], eval(i)['d15N_residual'], legend_label=eval(i)['names'][0], marker='circle', size=8, color=eval(i)['symbol_color'][0], line_color='black', alpha=0.8)
        figures[fig_n]['fig'].yaxis.axis_label = 'd15N vs AirN2 residual (permil)'
        figures[fig_n]['fig'].xaxis.axis_label = 'Analysis Number'
        figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].title.text_font_size = font_size
        figures[fig_n]['fig'].legend.label_text_font_size = font_size

        fig_n += 1

        figures[fig_n] = {}
        figures[fig_n]['cap'] = f"""Figure {fig_n}."""
        figures[fig_n]['fig'] = figure(title="d15N residual", width=1200, height=600, background_fill_color="#fafafa")
        for i in ref_mat['chosen']:
            figures[fig_n]['fig'].scatter(N_sam_AreaAll[eval(i)['index']], eval(i)['d15N_AirN2_residual'], legend_label=eval(i)['names'][0], marker='circle', size=8, color=eval(i)['symbol_color'][0], line_color='black', alpha=0.8)        
        figures[fig_n]['fig'].yaxis.axis_label = 'd15N vs AirN2 residual (permil)'
        figures[fig_n]['fig'].xaxis.axis_label = 'Peak Area (Vs)'
        figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].title.text_font_size = font_size
        figures[fig_n]['fig'].legend.label_text_font_size = font_size

        fig_n += 1

        figures[fig_n] = {}
        figures[fig_n]['cap'] = f"""Figure {fig_n}."""
        figures[fig_n]['fig'] = figure(title="d13C residual", width=1200, height=600, background_fill_color="#fafafa")
        for i in ref_mat['chosen']:
            figures[fig_n]['fig'].scatter(Analysis[eval(i)['index']], eval(i)['d13C_VPDB_residual'], legend_label=eval(i)['names'][0], marker='circle', size=8, color=eval(i)['symbol_color'][0], line_color='black', alpha=0.8)
        figures[fig_n]['fig'].yaxis.axis_label = 'd13C VPDB residual (permil)'
        figures[fig_n]['fig'].xaxis.axis_label = 'Analysis Number'
        figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].title.text_font_size = font_size
        figures[fig_n]['fig'].legend.label_text_font_size = font_size

        fig_n += 1

        figures[fig_n] = {}
        figures[fig_n]['cap'] = f"""Figure {fig_n}."""
        figures[fig_n]['fig'] = figure(title="d13C residual vs Peak Area", width=1200, height=600, background_fill_color="#fafafa")
        for i in ref_mat['chosen']:
            figures[fig_n]['fig'].scatter(C_sam_AreaAll[eval(i)['index']], eval(i)['d13C_VPDB_residual'], legend_label=eval(i)['names'][0], marker='circle', size=8, color=eval(i)['symbol_color'][0], line_color='black', alpha=0.8)
        figures[fig_n]['fig'].yaxis.axis_label = 'd13C vs VPDB residual (permil)'
        figures[fig_n]['fig'].xaxis.axis_label = 'Peak Area (Vs)'
        figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
        figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
        figures[fig_n]['fig'].title.text_font_size = font_size
        figures[fig_n]['fig'].legend.label_text_font_size = font_size

        fig_n += 1




    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(title="d15N vs Nqty", width=1200, height=600, background_fill_color="#fafafa")
    for i in ref_mat['id1_set']:
        figures[fig_n]['fig'].scatter(Nqty[eval(i)['index']], d15N_AirN2[eval(i)['index']], legend_label=eval(i)['names'][0], marker='circle', size=8, color=eval(i)['symbol_color'][0], line_color='black', alpha=0.8)
    figures[fig_n]['fig'].scatter(Nqty[sample_indices], d15N_AirN2[sample_indices], legend_label="samples", marker='triangle', size=8, color='black', alpha=0.8)
    figures[fig_n]['fig'].yaxis.axis_label = 'd15N vs AirN2 (permil)'
    figures[fig_n]['fig'].xaxis.axis_label = 'Nitrogen amount (mg)'
    figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].title.text_font_size = font_size
    figures[fig_n]['fig'].legend.label_text_font_size = font_size

    fig_n += 1

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(title="d13C vs Cqty", width=1200, height=600, background_fill_color="#fafafa")
    for i in ref_mat['id1_set']:
        figures[fig_n]['fig'].scatter(Cqty[eval(i)['index']], d13C_VPDB[eval(i)['index']], legend_label=eval(i)['names'][0], marker='circle', size=8, color=eval(i)['symbol_color'][0], line_color='black', alpha=0.8)
    figures[fig_n]['fig'].scatter(Cqty[sample_indices], d13C_VPDB[sample_indices], legend_label="samples", marker='triangle', size=8, color='black', alpha=0.8)
    figures[fig_n]['fig'].yaxis.axis_label = 'd13C vs VPDB (permil)'
    figures[fig_n]['fig'].xaxis.axis_label = 'Carbon amount (mg)'
    figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].title.text_font_size = font_size
    figures[fig_n]['fig'].legend.label_text_font_size = font_size

    fig_n += 1

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(title="d15N vs PercentN", width=1200, height=600, background_fill_color="#fafafa")
    figures[fig_n]['fig'].scatter(Nqty[sample_indices]/Amount[sample_indices]*100, d15N_AirN2[sample_indices], legend_label="samples", marker='triangle', size=8, color='black', alpha=0.8)
    figures[fig_n]['fig'].yaxis.axis_label = 'd15N vs AirN2 (permil)'
    figures[fig_n]['fig'].xaxis.axis_label = 'Percent Nitrogen (%)'
    figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].title.text_font_size = font_size
    figures[fig_n]['fig'].legend.label_text_font_size = font_size

    fig_n += 1

    figures[fig_n] = {}
    figures[fig_n]['cap'] = f"""Figure {fig_n}."""
    figures[fig_n]['fig'] = figure(title="d13C vs PercentC", width=1200, height=600, background_fill_color="#fafafa")
    figures[fig_n]['fig'].scatter(Cqty[sample_indices]/Amount[sample_indices]*100, d13C_VPDB[sample_indices], legend_label="samples", marker='triangle', size=8, color='black', alpha=0.8)
    figures[fig_n]['fig'].yaxis.axis_label = 'd13C vs VPDB (permil)'
    figures[fig_n]['fig'].xaxis.axis_label = 'Percent Carbon (%)'
    figures[fig_n]['fig'].xaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.axis_label_text_font_size = font_size
    figures[fig_n]['fig'].xaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].yaxis.major_label_text_font_size = font_size
    figures[fig_n]['fig'].title.text_font_size = font_size
    figures[fig_n]['fig'].legend.label_text_font_size = font_size



    # -------------------- REPORT SETUP -------------------

    report_directory = os.path.join(project_directory, f"{current_run_name}_report")
    if os.path.exists(report_directory):
        shutil.move(report_directory, os.path.join(archive_directory, f"{report_directory}_archive_{int(dt.datetime.utcnow().timestamp())}"))
    os.mkdir(report_directory)
    os.mkdir(os.path.join(report_directory, "data/"))
    os.mkdir(os.path.join(report_directory, "python/"))
    shutil.copy2(os.path.join(python_directory, 'CN_report.css'), report_directory)
    [shutil.copy2(os.path.join(python_directory, script), os.path.join(report_directory, f"python/{script}_REPORT_COPY")) for script in python_scripts]
    shutil.copy2(os.path.join(method_directory, log_file_name), os.path.join(report_directory, 'data/'))
    report_page = os.path.join(report_directory, f'{current_run_name}_calibration_summary.html')



    # -------------------- CALIBRATED DATA FILE --------------------
    print(f'\n    Creating calibrated data file.')
    calibrated_data_filename = f'{current_run_name}_CN_calibrated_data.csv'
    calibrated_data_file = os.path.join(method_directory, calibrated_data_filename)
    calibrated_file_headers = ['Sample ID', 'Date', 'Analysis Number', 'Total Mass (mg)', 'Nitrogen mass (mg)', 'd15N vs AirN2 (permil)', 'Carbon mass (mg)', 'd13C vs VPDB (permil)']
    data_to_write = '[Identifier1[ii], Date[ii], int(Analysis[ii]), Amount[ii], round(Nqty[ii], 3), round(d15N_AirN2[ii], 2), round(Cqty[ii], 3), round(d13C_VPDB[ii], 2)]'
    data_to_write = str(data_to_write).replace("'", "")
    with open(calibrated_data_file, 'w', newline='') as csvfile:
        datawriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
        datawriter.writerow(calibrated_file_headers)
        for ii in non_samples_indices:
            datawriter.writerow(eval(data_to_write))
        datawriter.writerow([])
        for ii in sample_indices:
            datawriter.writerow(eval(data_to_write))


    # -------------------- SUMMARY OF RUNS --------------------
    summary_file_headers = ['Run', 'Sample ID', 'Date', 'Analysis Number', 'Total Mass (mg)', 'Nitrogen mass (mg)', 'd15N vs AirN2 (permil)', 'Carbon mass (mg)', 'd13C vs VPDB (permil)']
    summary_data_to_write = '[current_run_name, Identifier1[ii], Date[ii], int(Analysis[ii]), Amount[ii], round(Nqty[ii], 3), round(d15N_AirN2[ii], 2), round(Cqty[ii], 3), round(d13C_VPDB[ii], 2)]'
    summary_data_to_write = str(summary_data_to_write).replace("'", "")


    if os.path.exists(summary_data_file):
        summary_file_exists = True
    else:
        summary_file_exists = False

    with open(summary_data_file, 'a', newline='') as csvfile:
        datawriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
        if summary_file_exists is False:
            datawriter.writerow(summary_file_headers)
        for ii in non_samples_indices:
            datawriter.writerow(eval(summary_data_to_write))
        for ii in sample_indices:
            datawriter.writerow(eval(summary_data_to_write))



    # ---------- REPORT BITS ---------- 
    print('Making html page...')
    header = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <!-- py by Andrew Schauer -->
            <meta http-equiv="Content-Type" content="text/html charset=UTF-8" />
            <meta name="viewport" content="width=device-width,initial-scale=1">
            <link rel="stylesheet" type="text/css" href="CN_report.css">
            <title>CN Calibration Report</title>
        </head>"""

    calculation_notes_block = str([f"<li>{i}</li>" for i in calculation_notes]).replace("[", "").replace("'", "").replace("]", "").replace(", ", "")
    refmat_block = str([f"<tr><td>{eval(i)['names'][0]}</td><td>{eval(i)['material']}</td><td>{eval(i)['d15N']}</td><td>{eval(i)['fractionN']}</td><td>{eval(i)['d13C']}</td><td>{eval(i)['fractionC']}</td><td>{eval(i)['purpose']}</td></tr>" for i in ref_mat['id1_set']]).replace("[", "").replace("'", "").replace("]", "").replace(", ", "")


    data_quality_block_1 = ""
    for i in ref_mat['qaqc']:
        data_quality_block_1 += f"""<tr><td>&delta;<sup>15</sup>N</td><td>{np.round(np.std(d15N_AirN2[eval(i)['index']])*2,3)} &permil;</td>
                                       <td>{np.round(np.mean(d15N_AirN2[eval(i)['index']]) - eval(i)['d15N'], 3)} &permil;</td><td>{i}</td></tr>
                                   <tr><td>&delta;<sup>13</sup>C</td><td>{np.round(np.std(d13C_VPDB[eval(i)['index']])*2,3)} &permil;</td>
                                       <td>{np.round(np.mean(d13C_VPDB[eval(i)['index']]) - eval(i)['d13C'],3)} &permil;</td><td>{i}</td></tr>
                                 <tr><td>N quantity</td>
                                     <td>{np.round(np.nanstd([Nqty[eval(i)['index']] - eval(i)['fractionN']*Amount[eval(i)['index']]])*1000)*2} &micro;g</td>
                                     <td>{np.round(np.nanmean([Nqty[eval(i)['index']] - eval(i)['fractionN']*Amount[eval(i)['index']]]))*1000} &micro;g</td>
                                     <td>{i}</td></tr>
                                 <tr><td>C quantity</td>
                                     <td>{np.round(np.nanstd([Cqty[eval(i)['index']] - eval(i)['fractionC']*Amount[eval(i)['index']]])*1000)*2} &micro;g</td>
                                     <td>{np.round(np.nanmean([Cqty[eval(i)['index']] - eval(i)['fractionC']*Amount[eval(i)['index']]]))*1000} &micro;g</td>
                                     <td>{i}</td></tr>"""

    data_quality_block_2 = f"""<tr><td><br></td><td> </td><td> </td><td> </td></tr>
                               <tr><td>&delta;<sup>15</sup>N</td><td>{np.round(ref_mat['d15N_AirN2_residual_std']*2, 3)} &permil;</td>
                                   <td> </td><td>all isotope reference materials</td></tr>
                               <tr><td>&delta;<sup>13</sup>C</td><td>{np.round(ref_mat['d13C_VPDB_residual_std']*2, 3)} &permil;</td>
                                   <td> </td><td>all isotope reference materials</td></tr>"""

    body = f"""
        <body>
        <div class="entire_page">
        <h2>CN Stable Isotope Analysis and Calibration Report</h2>
        <div class="created-date">Created - {str(dt.datetime.now())}</div>
        <h2>Introduction</h2>
        <div class="text-indent">
            <p>This report is meant to be a stand-alone collection of methods,
            data, scripts, and notes related to carbon and nitrogen isotopic
            analysis. Your samples were analyzed on a {config['methods']['instrumentation']},
            for d13C and d15N of solid material. You can read more
            about our implementation of this method on our website
            <a href="{config['methods']['procedure_link']}">
            {config['methods']['procedure']}</a>.</p>

            <p>The data and python scripts used to generate this page are linked
            and described in the <a href="#refs">References</a> section below. If you wish
            to save this report, <a href="report.zip">save the zip file</a> or copy and paste
            or download the entire 'report' directory to a place of your choosing and all
            html, images, data files, and python scripts will be saved. <strong>
            &iexcl; <a href="report.zip">Save a copy if you are finished analyzing your samples</a> !</strong></p>
        </div>

        <h2>My data</h2>
        <div class="text-indent">
            <p>This technical stuff is fine and all but where are <strong><a href="data/{calibrated_data_filename}">my data</a></strong>?
            This calibrated data file contains sample IDs, dates of analyses, unique analysis numbers, total mass weighed for analysis,
            mass of nitrogen, d15N versus Air N2, mass of carbon, d13C vs VPDB, flags, and notes. Each section of data is separated
            by an empty row. The first section of data are the trusted reference materials; the second section of data are trusted samples;
            the third section of data are untrusted. Under the "trust" heading, "1" indicates good, trusted data while "0" indicates poor
            quality data that should probably be distrusted. Up to you if you want to use it. Untrusted data are given the reason for
            distrust. <strong>&iexcl; <a href="report.zip">If you are done analyzing samples, save a copy of the entire report directory elsewhere,
            not just a copy of your data file</a> !</strong></p>
        </div>

        <h2>Data operations</h2>
        <div class="text-indent">
            <p>A suite of mathmatical operations were completed on these data prior to claiming they are final. Here are notes associated with
            these calculations:
            <ul>
            {calculation_notes_block}
            </ul>
            </p>
        </div>

        <h2>Run inventory</h2>
        <div class="text-indent">
            <table>
                <tr><td>Total number of project runs</td><td>{len(set(file))}</td></tr>
                <tr><td>Total number of analyses in run</td><td>{len(set(Analysis))}</td></tr>
                <tr><td>Total number of standards in run</td><td>{len(non_samples_indices)}</td></tr>
                <tr><td>Total number of samples in run</td><td>{len(sample_indices)}</td></tr>
                <tr><td><br></td></tr>
                <tr><td>Number of blanks in run</a></td><td>{len(blank['index'])}</td></tr>
                <tr><td>Number of qtycal in run</a></td><td>{len(qtycal['index'])}</td></tr>
                <tr><td>Number of <a href="#excluded">excluded analyses in run</a></td><td>{len(trust0_indices)}</td></tr>
            </table>
        </div>

        <h2>Reference materials</h2>
        <div class="text-indent">
            <p>All internationally recognized reference material accepted values can be found at the CIAAW (http://www.ciaaw.org/).
            All IsoLab in-house reference material accepted values can be found at http://isolab.ess.washington.edu/isolab/reference-materials#solid.
            For this particular analysis, the accepted values are from {reference_materials['file_meta_data']['file']} - {reference_materials['file_meta_data']['modification_date']}.
            The reference materials and their accepted values, normalized to the Air-N2 and VPDB scales, included in this run / set are:</p>
            <table>
                <tr><th>Reference<br>name</th><th>Reference<br>material</th><th>d15N<br>accepted<br>(permil)</th><th>Percent<br>Nitrogen<br>(%)</th><th>d13C<br>accepted<br>(permil)</th><th>Percent<br>Carbon<br>(%)</th><th>Purpose</th></tr>
                {refmat_block}
            </table>
        </div>

        <h2>Data quality</h2>
        <div class="text-indent"><p>Precision and accuracy estimates are derived from quality control reference material(s) with known &delta;<sup>15</sup>N and
            &delta;<sup>13</sup>C values but treated as an unknown sample. Additional precision estimates come directly from repeated unknown sample analyses.
            Precision is <strong>two standard deviations</strong> over all replicates of the quality control reference material. Accuracy is
            the difference of the mean of all replicates of the quality control reference material from the accepted value. Another more conservative
            estimate of precision uses a pooled residual 2-&sigma;. Typical &delta;<sup>15</sup>N long-term 2-&sigma; precision is 0.4 &permil;. 
            Typical &delta;<sup>13</sup>C long-term 2-&sigma; precision is 0.2 &permil;. Accuracy should be with the precision estimate of zero permil; if
            accuracy is outside greater or less than the precision estimate from zero permil, something is likely wrong with the run, calibration, or data processing.</p>
            <table>
                <tr><th>Parameter</th><th>Precision</th><th>Accuracy</th><th>Reference<br>Material(s)<br>Used</th></tr>
                {data_quality_block_1}
                {data_quality_block_2}
            </table><br>
            <p><strong>N2 blank</strong>: 
            <ul>
                <li>Mean peak area: {np.round(np.nanmean(N_sam_AreaAll[blank['index']]),3)} Vs</li>
                <li>Mean nitrogen quantity: {np.round(np.nanmean(N_sam_AreaAll[blank['index']]/qtycal['Nfit'][0]*1000),2)} &micro;g</li>
                <li>Mean &delta;<sup>15</sup>N: {np.round(np.nanmean(N_sam_d15N14N[blank['index']]),2)} permil</li>
                <li>n: {len(blank['index'])}</li>
            </ul>
            </p>
        </div>

        <h2>Figures</h2>"""

    figure_block = [f"""<div class="clear-both">{file_html(figures[i]['fig'], CDN)}{figures[i]['cap']}<hr></div>""" for i in figures.keys()]

    python_scripts_block = str([f'<li><a href="python/{key}_REPORT_COPY">{key}</a> - {value}</li>' for key, value in python_scripts.items()]).replace("[", "").replace("'", "").replace("]", "").replace(", ", "")

    footer = f"""
        <h2 id="refs">References</h2>
        <div class="references">
        <ul>
            <li>Python scripts - modification date:
                <ul>
                    {python_scripts_block}
                </ul>
            <li><a href="https://github.com/andyschauer/carbon-nitrogen">github repository</a></li>
            <li>Data files - <a href="data/{calibrated_data_filename}">{calibrated_data_filename}</a></li>
            <li><a href="https://isolab.ess.washington.edu/laboratory/solid-CN.php">IsoLab's carbon and nitrogen analysis overiew.</a></li>
            <li><a href="https://isolab.ess.washington.edu/SOPs/shrek-cn.php">IsoLab's water carbon and nitrogen analysis method.</a></li>
            <li><a href="report.zip">Zip file of entire report directory.</a></strong>.</li>
        </ul>
        </div>
        </body></div></html>"""



    # -------------------- WRITE REPORT --------------------
    with open(report_page, 'w') as report:
        report.write(header)
        report.write(body)
        [report.write(i) for i in figure_block]
        report.write(footer)
        report.close()
    webbrowser.open(report_page)



    # -------------------- REPORT ZIP --------------------
    shutil.make_archive('report', 'zip', report_directory)
    shutil.move('report.zip', os.path.join(report_directory, 'report.zip'))




