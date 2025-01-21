#!/usr/bin/env python3
"""
This script opens raw carbon and nitrogen stable isotope elemental analysis data
files from shrek, a ThermoFinnigan 253 with Eurovector Elemental Analyzer, and
writes salient data to shrekCN_analysis_log.csv

Version 0.1 mod date 2021-02-03 => created
Version 0.2 mod date 2021-02-25 => working version
Version 0.3 mod date 2021-03-01 => added Nonly and Conly, moves raw file to either archive or unread directories,
        tested it with 2008 to 2018 data
Version 0.4 mod date 2021-03-03 => changed Height to Ampl28 and Ampl44
Version 0.5 mod date 2021-10-08 => updating nomenclature
Version 0.6 mod date 2022-02-17 => updating, linting, more nomenclature
Version 0.7 mod date 2022-03-23 => export data and headers comes from shrekCN_standards.py; script tries to detect
        computer you are working on and adjust directories accordingly;
Version 0.8 mod date 2022-12-31 => updated, linted, tested. I think this now works after a complete reprocess of all shrek CN data.
Version 1.0 mod date 2024-04-08 => changed lab to isolab_lib; changed flag to trust, some refactoring for readability,
        function to handle sample notestrying to get this finished up to a version 1 level and upload to github
Version 1.1 mod date 2024-04-23 => found possible bug in joining notes. Used a comma. I think it needed to be a semi-colon or something other than a comma since we are writing a csv file
Version 1.2 mod date 2024-06-03 => added project level analysis log file and ability to choose one or create one. The script still writes to an exhaustive project-independent analysis log file
Version 2.0 mod date 2024-06-21 => removed isolab_lib in favor of single CN_lib.py file
Version 2.1 mod date 2024-06-23 => changed shrekCN to CN throughout
Version 2.2 mod date 2024-11-22 => found bug in GasConfiguration names and fixed
"""

__author__ = "Andy Schauer"
__email__ = "aschauer@uw.edu"
__last_modified__ = "2024-11-22"
__version__ = "2.2"
__copyright__ = "Copyright 2025, Andy Schauer"
__license__ = "Apache 2.0"
__acknowledgements__ = "Shrek"




# -------------------- imports --------------------
import csv
from numpy import where
import os
import re
from CN_lib import *
import sys
import time



# -------------------- functions --------------------
def append_meta_data():
    for i in meta_data:
        meta_data[i].append(data[i][index]) if i in data else meta_data[i].append(None)


def append_supp_data():
    supp_data['file'].append(file)
    supp_data['trust'].append(trust)
    supp_data['notes'].append(note)
    supp_data['pyversions'].append(version)
    supp_data['empty'].append('')

    if 'Information' in data:
        if data['Information'][index] is None:
            supp_data['peak_center'].append(None)
        elif re.match('Peak Center found at', data['Information'][index]) is None:
            supp_data['peak_center'].append(None)
        else:
            m = re.findall(r'\d+', data['Information'][index])
            supp_data['peak_center'].append(m[0])
    else:
        supp_data['peak_center'].append(None)


def append_N_wg_data():  # put data for nitrogen refence peak into list
    for i in N_wg_data:
        N_wg_data[i].append(data[i][index + peak_number_offset]) if i in data else N_wg_data[i].append(None)


def append_N_sam_data():  # put data for nitrogen sample peak into list
    for i in N_sam_data:
        N_sam_data[i].append(data[i][index + peak_number_offset]) if i in data else N_sam_data[i].append(None)
    if int(data['Ampl28'][index + peak_number_offset]) > 49950 or int(data['Ampl29'][index + peak_number_offset]) > 49950:
        global trust
        trust = 0
        sample_note('N2 cup saturated')


def append_C_sam_data():  # put data for carbon sample peak into list
    for i in C_sam_data:
        C_sam_data[i].append(data[i][index + peak_number_offset]) if i in data else C_sam_data[i].append(None)
    if int(data['Ampl44'][index + peak_number_offset]) > 49950 or int(data['Ampl45'][index + peak_number_offset]) > 49950 or int(data['Ampl46'][index + peak_number_offset]) > 49950:
        global trust
        trust = 0
        sample_note('CO2 cup saturated')


def append_C_wg_data():  # put data for carbon refernce peak into list
    for i in C_wg_data:
        C_wg_data[i].append(data[i][index + peak_number_offset]) if i in data else C_wg_data[i].append(None)


def N_wg_none():  # sets all nitrogen reference peak data to None
    for i in N_wg_data:
        N_wg_data[i].append(None)


def N_sam_none():  # sets all nitrogen sample peak data to None
    for i in N_sam_data:
        N_sam_data[i].append(None)


def C_sam_none():  # sets all carbon sample peak data to None
    for i in C_sam_data:
        C_sam_data[i].append(None)


def C_wg_none():  # sets all carbon reference peak data to None
    for i in C_wg_data:
        C_wg_data[i].append(None)


def sample_note(currnote):
    if trust == 0:
        print(f"\033[91m Analysis {str(data['Analysis'][index])} in file {file} - {currnote}\033[0m")
    else:
        print(f" Analysis {str(data['Analysis'][index])} in file {file} - {currnote}")
    note.append(currnote)



# ---------- LOAD CONFIGURATION ----------
with open("CN_config.json", 'r') as f:
    config = json.load(f)


# -------------------- SETUP --------------------

version = os.path.basename(__file__) + ' - ' + time.ctime(os.path.getctime(__file__))

home_directory = config["local_directories"]["home"]
python_directory = f"{home_directory}{config['local_directories']['python']}"
method_directory = f"{home_directory}{config['local_directories']['method_data_directory']}"
reference_materials_file = f"{home_directory}{config['local_directories']['standards']}"

new_data_directory = 'rawdata_new'
archive_data_directory = 'rawdata_archive'
junk_data_directory = 'rawdata_junk'
exhaustive_log_file_name = 'CN_exhaustive_analysis_log.csv'

if os.path.isdir(method_directory) is False:
    print('directory does not exist...exiting....')
    sys.exit()

CN_log_file_list = make_file_list(method_directory, '_analysis_log.csv')

print('\nWhere do you want all this raw data to go?\n')

[print(f'    {i}') for i in CN_log_file_list]
identified_file = 0
while identified_file == 0:
    CN_log_file_search = input('\nEnter a project analysis log file from above that you wish to append raw data to or leave blank to create a new one: ')
    if CN_log_file_search:
        isfile = [CN_log_file_search[0: len(CN_log_file_search)] in x for x in CN_log_file_list]
        if len(where(isfile)[0]) == 1:
            identified_file = 1
            project_log_file_name = CN_log_file_list[where(isfile)[0][0]]
            print(f'    Appending to CN log file {project_log_file_name}...')
        else:
            print('\n** More than one file found. **\n')
    else:
        print("\n\nCreate a project analysis log file.")
        project_log_file_name = input("Enter a project key word: ").replace(' ', '_')
        project_log_file_name += "_analysis_log.csv"
        identified_file = 1


# -------------------- get list of files --------------------
filelist = make_file_list(os.path.join(method_directory, new_data_directory), 'csv')
if not filelist:
    print('No files in raw data directory.')
else:
    filelist = sorted(filelist)



# -------------------- Main loop --------------------
for file in filelist:
    for i in meta_headers:
        meta_data[i] = []

    for i in N_headers:
        N_wg_data[i] = []
        N_sam_data[i] = []

    for i in C_headers:
        C_wg_data[i] = []
        C_sam_data[i] = []

    for i in supp_headers:
        supp_data[i] = []

    print(f"\nReading in file {file}...")
    headers, data = read_file(os.path.join(method_directory, new_data_directory, file), ',')  # read file and return headers (headers) and data (data)

    # test for file problems
    if 'Analysis' not in headers:
        print('problem with file ' + file)
        os.rename(os.path.join(method_directory, new_data_directory, file), os.path.join(method_directory, junk_data_directory, file))  # done with datafile, put it in the unread directory
        print('file ' + file + ' was moved to the junk folder')
        continue

    for n in data['Analysis']:
        try:
            data['Analysis'] = [int(index) for index in data['Analysis']]
        except ValueError:
            print(f'File {file} contains strings in the Analysis column')
            os.rename(os.path.join(method_directory, new_data_directory, file), os.path.join(method_directory, junk_data_directory, file))  # done with datafile, put it in the unread directory
            print(f'File {file} was moved to the junk folder.')
            continue

    # sample_index_first_row = first row of each sample
    analysis_difference = [j - index for index, j in zip(data['Analysis'][: -1], data['Analysis'][1:])]
    sample_index_first_row = [index for index, e in enumerate(analysis_difference) if e != 0]
    sample_index_first_row = [index + 1 for index in sample_index_first_row]
    sample_index_first_row.append(0)
    sample_index_first_row = sorted(sample_index_first_row)

    # rows_per_sample
    rows_per_sample = [j - index for index, j in zip(sample_index_first_row[:-1], sample_index_first_row[1:])]
    rows_per_sample.append(len(data['Analysis']) - sample_index_first_row[len(sample_index_first_row) - 1])  # adds last sample number of rows

    # sample_index_last_row = last row of each sample
    sample_index_last_row = [j + index for index, j in zip(sample_index_first_row, rows_per_sample)]

    # run type - N, C, CN
    gas_config_names_list = ['GasConfiguration', 'Gasconfiguration']
    gas_config_name = [gas_config_name for gas_config_name in gas_config_names_list if gas_config_name in list(data.keys())]
    if gas_config_name:
        gas_configuration = set(data[gas_config_name[0]])
    else:
        gas_configuration = 'undefined'
    if any(index in ['N2', 'CO2'] for index in gas_configuration):
        if all(index in gas_configuration for index in ['N2', 'CO2']):
            run_type = 'CN'
        elif 'N2' not in gas_configuration:
            run_type = 'C'
        elif 'CO2' not in gas_configuration:
            run_type = 'N'
    else:
        run_type = 'unk'

    if any([run_type == 'CN', run_type == 'N', run_type == 'C']):  # only run this code if the run_type is CN, N, or C

        for index, rows in zip(sample_index_first_row, rows_per_sample):
            append_meta_data()
            note = []
            trust = 1

            # Identify peaks based on run type and number of peaks per analysis
            if run_type == 'CN':

                if rows == 1:  # 1 peak in CN mode means problems, so don't export any data and let user know
                    trust = 0
                    sample_note('Only one peak detected; ')
                    peak_number_offset = 0
                    N_wg_none()
                    N_sam_none()
                    C_sam_none()
                    C_wg_none()

                elif rows == 2:  # 2 peaks in a CN run could mean Nref and Cref, or a problem happened and there are no carbon peaks, index can't think of a scenario in a CN run where you would have no nitrogen peaks but would have carbon peaks
                    # PeakNr1 is Nref
                    peak_number_offset = 0
                    append_N_wg_data()

                    # PeakNr2 could be Nsam or Cref or some other problem
                    peak_number_offset = 1
                    if data[gas_config_name[0]][index + peak_number_offset] == 'N2':  # Nsam is present but Csam is missing
                        sample_note('No carbon peaks')
                        append_N_sam_data()
                        C_sam_none()
                        C_wg_none()
                    elif data[gas_config_name[0]][index + peak_number_offset] == 'CO2':  # Nsam and Csam are missing but Cref is present
                        sample_note('No nitrogen or carbon sample peaks')
                        N_sam_none()
                        C_sam_none()
                        append_C_wg_data()

                elif rows == 3:  # second most common scenario where peaks are either Nref, Nsam, Cref or Nref, Csam, Cref
                    # PeakNr1 is Nref
                    peak_number_offset = 0
                    append_N_wg_data()

                    # PeakNr2
                    peak_number_offset = 1
                    if data[gas_config_name[0]][index + peak_number_offset] == 'N2':  # Nsam is present but Csam is missing
                        append_N_sam_data()
                        C_sam_none()
                    elif data[gas_config_name[0]][index + peak_number_offset] == 'CO2':  # Nsam is missing but Csam is present
                        N_sam_none()
                        append_C_sam_data()

                    # PeakNr3 is Cref
                    peak_number_offset = 2
                    append_C_wg_data()

                elif rows == 4:  # most common scenario where peaks are Nref, Nsam, Csam, Cref
                    # PeakNr1 is Nref
                    peak_number_offset = 0
                    append_N_wg_data()

                    # PeakNr2 is Nsam
                    peak_number_offset = 1
                    append_N_sam_data()

                    # PeakNr3 is Csam
                    peak_number_offset = 2
                    append_C_sam_data()

                    # PeakNr4 is Cref
                    peak_number_offset = 3
                    append_C_wg_data()

                elif rows == 5:
                    # determine which gas config has extra peak then assign peaks and exclude extra but make note of it
                    if data[gas_config_name[0]][index:index + rows].count('N2') == 2:  # CO2 has extra peak
                        sample_note('extra CO2 peak detected')
                        # PeakNr1 is Nref
                        peak_number_offset = 0
                        append_N_wg_data()

                        # PeakNr2 is Nsam
                        peak_number_offset = 1
                        append_N_sam_data()

                        # PeakNr 3 or 4 is Csam
                        test = data['AreaAll'][index + peak_number_offset + 1] > data['AreaAll'][index + peak_number_offset + 2] if type(data['AreaAll'][index + peak_number_offset + 1]) == int and type(data['AreaAll'][index + peak_number_offset + 2]) == int else None
                        if test is True or test is False:  # if peak 3 is larger than 4, peak 3 is sample
                            if test is True:
                                peak_number_offset = 2
                                sample_note('using larger first of two sample peaks')
                            elif test is False:
                                peak_number_offset = 3
                                sample_note('using larger second of two sample peaks')
                            append_C_sam_data()

                        else:
                            trust = 0
                            sample_note('problem with AreaAll data')
                            C_sam_none()

                        # PeakNr5 is Cref
                        peak_number_offset = 4
                        append_C_wg_data()

                    elif data[gas_config_name[0]][index:index + rows].count('CO2') == 2:  # N2 has extra peak
                        N_wg_none()
                        N_sam_none()
                        # PeakNr3 is Csam
                        peak_number_offset = 3
                        append_C_sam_data()

                        # PeakNr4 is Cref
                        peak_number_offset = 4
                        append_C_wg_data()
                        sample_note('5 peaks in sample and N2 has extra')
                    else:
                        N_wg_none()
                        N_sam_none()
                        C_sam_none()
                        C_wg_none()
                        trust = 0
                        sample_note('5 peaks in sample but something unexpected happened')

                elif rows == 6:
                    N_wg_none()
                    N_sam_none()
                    C_sam_none()
                    C_wg_none()
                    trust = 0
                    sample_note('6 peaks in sample')

                else:
                    N_wg_none()
                    N_sam_none()
                    C_sam_none()
                    C_wg_none()
                    trust = 0
                    sample_note('too many peaks causes anxiety')

            elif run_type == 'N':
                C_sam_none()
                C_wg_none()
                if rows == 4:
                    # PeakNr2 and 5 are Nref but index need to build in averaging these two peaks, for now it is PeakNr2
                    peak_number_offset = 1
                    append_N_wg_data()
                    N_sam_none()

                elif rows == 5:
                    # PeakNr2 and 5 are Nref but index need to build in averaging these two peaks, for now it is PeakNr2
                    peak_number_offset = 1
                    append_N_wg_data()

                    # PeakNr4 is Nsam
                    peak_number_offset = 4
                    append_N_sam_data()
                    sample_note('No nitrogen sample peaks')

                else:
                    N_wg_none()
                    N_sam_none()
                    trust = 0
                    sample_note('Number of peaks suggests an error')

            elif run_type == 'C':
                N_wg_none()
                N_sam_none()
                if rows == 4:
                    # PeakNr2 and 5 are Nref but index need to build in averaging these two peaks, for now it is PeakNr2
                    peak_number_offset = 1
                    append_C_wg_data()
                    C_sam_none()

                elif rows == 5:
                    # PeakNr2 and 5 are Nref but index need to build in averaging these two peaks, for now it is PeakNr2
                    peak_number_offset = 1
                    append_C_wg_data()

                    # PeakNr4 is Csam
                    peak_number_offset = 4
                    append_C_sam_data()
                    sample_note('No carbon sample peaks')

                else:
                    C_wg_none()
                    C_sam_none()
                    trust = 0
                    sample_note('Number of peaks suggests and error')

            else:  # run type is not known or something about identifying it went wrong
                N_wg_none()
                N_sam_none()
                C_sam_none()
                C_wg_none()
                trust = 0
                sample_note('Run type unknown')

            note = "; ".join(note)
            append_supp_data()


        # write data to the exhaustive analysis log
        if os.path.isfile(os.path.join(method_directory, exhaustive_log_file_name)) is False:
            with open(os.path.join(method_directory, exhaustive_log_file_name), 'w', newline='') as csvfile:  # if the log file has not been created, create it with column headers and data
                datawriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
                datawriter.writerow(CN_analysis_log_headers)
                for ii in range(len(meta_data['Analysis'])):
                    datawriter.writerow(eval(data_to_write))

        else:
            with open(os.path.join(method_directory, exhaustive_log_file_name), 'a', newline='') as csvfile:
                datawriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
                for ii in range(len(meta_data['Analysis'])):
                    datawriter.writerow(eval(data_to_write))

        # write data to the project analysis log
        if os.path.isfile(os.path.join(method_directory, project_log_file_name)) is False:
            with open(os.path.join(method_directory, project_log_file_name), 'w', newline='') as csvfile:  # if the log file has not been created, create it with column headers and data
                datawriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
                datawriter.writerow(CN_analysis_log_headers)
                for ii in range(len(meta_data['Analysis'])):
                    datawriter.writerow(eval(data_to_write))

        else:
            with open(os.path.join(method_directory, project_log_file_name), 'a', newline='') as csvfile:
                datawriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
                for ii in range(len(meta_data['Analysis'])):
                    datawriter.writerow(eval(data_to_write))

        os.rename(os.path.join(method_directory, new_data_directory, file), os.path.join(method_directory, archive_data_directory, file))  # done with datafile, put it in the archive directory

    else:  # not a C or N or CN data file
        print(f"{file} is not recognized as a C or N or CN data file; putting it in the junk folder")
        os.rename(os.path.join(method_directory, new_data_directory, file), os.path.join(method_directory, junk_data_directory, file))  # done with datafile, put it in the junk directory


print("\nDone")
print("\n\nYou may now wish to run CN_calibrate.py by typing 'python3 CN_calibrate.py' without the quotes. You may also type 'python3 -i CN_calibrate.py --verbose' for more figures.\n\n\n")
