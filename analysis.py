#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 5 13:17:31 2024

@author: mattmcdonald
"""

import pandas as pd
import numpy as np
import os
import csv
import scipy.signal as sp
from pyteomics import mzxml
from matplotlib import pyplot as plt
import yaml
from pymcr.mcr import McrAR
from pymcr.regressors import OLS
from pymcr.constraints import ConstraintNonneg
from scipy.signal import savgol_filter
import pybaselines.whittaker as baselines
import subprocess
import threading
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

plotting = False # True will generate two plots for each sample, if running many samples may crash matplotlib

""" UPDATE model_dir AND plates_root TO MATCH THE LOCATION OF THE MODEL AND DATA ON LOCAL MACHINE """
model_dir = os.path.join("mcdonald_et_al_data", "Deep4Chem_chemprop")
plates_root = "mcdonald_et_al_data"
plates = ["reaction_set_1", "reaction_set_2", "simulated_reactions"]

results_file = "figure_3_and_4_data.csv"


def log_e_chemprop_predict(wellplate: dict, file_name="log_e", solvents=['O','CC#N']):
    """
    Run chemprop log_e predictions on a wellplate. Dump results into file_name. Can also give solvents, acceptable
    solvents are (based on model training data) water ('O'), ACN ('CC#N'), DCM ('ClCCl'), and EtOH ('CCO'). 
    Chemprop is run as a subprocess the same way it would be run from the command line in Chemprop v1.X.X
    Starting in Chemprop v2.0 prediction commands can be made directly from python, this function should be
    updated to reflect that change for better performance and stability
    Parameters
    ----------
    smiles : dict
        contents dictionary
    file_name : str
        file_name to be used to generate and retrieve predictions
    solvents : list, optional
        list of solvents to predict in. The default is ['O','CC#N'].

    Returns
    -------
    The name of the file that the chemprop predictions will be available in
    """
    
    def chemprop_thread_func(command):
        process = subprocess.run(command, shell = True, capture_output = True)    
        formatted_stderr = []
        for line in process.stderr.decode().split('\r'):
            if 'WARNING' in line:
                continue
            clean_line = line.replace('\r', '').strip('\n')
            if clean_line.strip() == '':
                continue
            if '#' not in clean_line:
                continue
            formatted_stderr.append(clean_line)
        print(process.returncode)
        print('Finishing!: %s' % command)
    
    # generate csv file to be fed to chemprop
    to_calc_file = f"{file_name}_to_predict.csv"
    preds_file = f"{file_name}_predictions.csv"
    if os.path.exists(preds_file):
        return preds_file
    mol_solv_list = []
    for well in wellplate['contents'].values():
        for solv in solvents:
            mol_solv = [well['target_product'][0][0], solv]
            if mol_solv not in mol_solv_list:
                mol_solv_list.append(mol_solv)
    header = ['SMILES', 'solvent']
    with open(to_calc_file, 'w', newline='', encoding='utf-8') as new_csv:
        write = csv.writer(new_csv)
        write.writerow(header)
        write.writerows(mol_solv_list)
        
    # generate chemprop command string
    command_l1 = f"chemprop_predict --test_path '{to_calc_file}' --preds_path '{preds_file}' --checkpoint_dir '{model_dir}'"
    command_l2 = " --ensemble_variance --number_of_molecules 2"
    command = command_l1 + command_l2
    chemprop_thread = threading.Thread(target=chemprop_thread_func, args=(command,), daemon=True)
    chemprop_thread.start()
    print("Working on chemprop predictions")
    chemprop_thread.join()
    
    #with open(preds_file, 'w', newline='', encoding='utf-8') as new_csv:
    #    write = csv.writer(new_csv)
    #    write.writerow(['Results pending...'])
        
    return preds_file


def lc_wellplate_analysis(wellplate: dict, chemprop_file: str, data_folder: str):
    """
    Start analysis of well plate of reactions

    Parameters
    ----------
    wellplate : dict
        well plate dictionary as formatted in platform operation (examples provided in mcdonald_et_al_data).
    chemprop_file : str
        CSV file where the chemprop predictions have been saved.
    data_folder : str
        Directory containing the mass spec (.mzXML) and PDA (.txt) data, as exported from Shimadzu LabSolutions.

    Returns
    -------
    rxn_info : dict
        A dictionary of information regarding the yield prediction.

    """
    pda_file_table = preprocess_pda_files(data_folder)
    
    logE = {}  # predicted ext coef
    UV = {}    # predicted UV max
    em = {}    # Exact mass, to be calculated but for now pull from sheet
    with open(chemprop_file, 'r') as csv_in:
        reader = csv.reader(csv_in)
        next(reader) #skip header
        for line in reader:
            if line[0] in logE.keys():
                logE[line[0]][line[1]] = float(line[4])
            else:
                logE[line[0]] = {line[1]: float(line[4])}
            if line[0] in UV.keys():
                UV[line[0]][line[1]] = float(line[2])
            else:
                UV[line[0]] = {line[1]: float(line[2])}
            em[line[0]] = rdMolDescriptors.CalcExactMolWt(Chem.MolFromSmiles(line[0]))
    
    rxn_info = {}
    for well, content in wellplate['contents'].items():
        try:
            smiles = content['target_product'][0][0]
        except TypeError:
            smiles = content['reaction_smiles']
        # long term should calculate the exact mass with rdkit, but harder to debug from rdkit env
        if plotting:
            fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
            ret_time, best_mz = ms_data_analysis(em[smiles], 
                                                 os.path.join(data_folder, f"{well}.mzXML"),
                                                 ax1)
        else:
            ret_time, best_mz = ms_data_analysis(em[smiles], 
                                                 os.path.join(data_folder, f"{well}.mzXML"))
        print(f"Analyzing well {well}")
        if ret_time > 0:
            if plotting:
                rxn_info[well] = pda_data_analysis(logE[smiles], UV[smiles], ret_time, pda_file_table[well], ax2)
                plt.title(well)
                plt.show()
            else:
                rxn_info[well] = pda_data_analysis(logE[smiles], UV[smiles], ret_time, pda_file_table[well])
            rxn_info[well]['smiles'] = smiles
        else:
            continue
    return rxn_info

# Associate PDA data (.txt) to the corresponding MS data (.mzXML) and well in well plate
def preprocess_pda_files(folder: str, ):
    pda_file_table = {}
    for file in os.listdir(folder):
        if file.split('.')[-1] == 'txt' and file[0] != '.':
            with open(os.path.join(folder, file), 'r') as fin:
                for line in fin:
                    if line[0:4] == "Data":
                        well = line[-20:-1].split('\\')[-1].split('.')[0]
                        pda_file_table[well] = os.path.join(folder, file)
                        break
    return pda_file_table

# For processing blank chromatograms for background subtraction
def process_blank(files: list, ):
    collated = []
    for file in files:
        with open(file, 'r') as data_in:
            lines = data_in.readlines()
        for i, line in enumerate(lines):
            if len(line) > 200:
                break
        data = pd.read_csv(file, header=i-2, index_col=0)
        collated.append(data)     
    average = np.zeros(collated[0].shape)
    for sample in collated:
        average += sample/len(collated)
    return average

# Extract the retention time corresponding to the desired product, if reaction has any yield at all
def ms_data_analysis(exact_mass: float, mzxml_file: str, ax=None):
    # peak identification parameters: number of consecutive points m std devs above the baseline noise
    n_con_pts = 10
    m_sigma = 5
    def baseline_helper(data, n_segs, n_toss):
        # helper function to find standard deviation of the baseline
        # randomly remove points to resize data into n_segs equal sized segments
        data = np.delete(data, np.random.permutation(len(data))[0:(len(data) % n_segs)])
        data = data.reshape((-1, n_segs), order='F')
        mus = np.mean(data, axis=0)
        idx = np.argsort(mus)
        data = np.delete(data, idx[-n_toss:], 1)
        return np.mean(data), np.std(data)

    # look for adducts [+proton, -proton, +ACN+H]
    mws = np.array([exact_mass + 1, 1 - exact_mass, exact_mass + 42])

    experiment = mzxml.read(mzxml_file)
    spec = experiment.next()
    # rows ought to be the total length of the experiment divided by 2... 
    rows = int(480)
    cols = int(2500)
    upsample_pts = np.linspace(spec['startMz'], spec['endMz'], cols)

    # even scans are positive polarity, odd scans are negative polarity
    flipper = 0
    rt_pos = np.ndarray(shape=(rows, 1), dtype=np.float32)
    rt_neg = np.ndarray(shape=(rows, 1), dtype=np.float32)
    data_pos = np.ndarray(shape=(cols, rows), dtype=np.float32)
    data_neg = np.ndarray(shape=(cols, rows), dtype=np.float32)
    i = 0
    rt_pos[i] = spec['retentionTime']
    for spec in experiment:
        mzs = spec['m/z array']
        ints = spec['intensity array']
        rt = spec['retentionTime']
        if flipper:
            data_pos[:, i] = np.interp(upsample_pts, mzs, ints)
            flipper -= 1
            rt_pos[i] = rt
        else:
            data_neg[:, i] = np.interp(upsample_pts, mzs, ints)
            flipper += 1
            rt_neg[i] = rt
            i += 1

    ret_time = -1.0
    max_intensity = 0
    intensities = [] # list of [summed_intensity, ret_time] pairings
    best_mz = 0
    # indentify peaks based on N consecutive points > M*sigma baseline
    try:
        for mz in mws:
            if mz > 0:
                mic_low = np.where(upsample_pts > (mz - 0.5))[0][0]
                mic_high = np.where(upsample_pts < (mz + 0.5))[0][-1]
                mic_pos = np.sum(data_pos[mic_low:mic_high, :], 0)
                mic_pos = mic_pos / max(mic_pos)
                mu, sigma = baseline_helper(mic_pos, 10, 2)
                idx = np.array(np.where(mic_pos > (mu + m_sigma * sigma)))
                for i in range(0, n_con_pts - 1):
                    idx = idx[np.diff(idx, append=-1) == 1]
                j = 0 
                for i in idx:
                    if i == j+1:
                        intensities[-1][0] += ((mic_pos[i] - mu)/sigma)
                    else:
                        intensities.append([0, float(rt_pos[i])])
                    j = i
                
                if intensities:
                    ret_time = max(intensities)[1]
                    if max(intensities)[0] > max_intensity:
                        max_intensity = max(mic_pos)
                        best_mz = mz
                        chromatogram = mic_pos
            else:
                mic_low = np.where(upsample_pts > (abs(mz) - 0.5))[0][0]
                mic_high = np.where(upsample_pts < (abs(mz) + 0.5))[0][-1]
                mic_neg = np.sum(data_neg[mic_low:mic_high, :], 0)
                mic_neg = mic_neg / max(mic_neg)
                mu, sigma = baseline_helper(mic_neg, 10, 2)
                idx = np.array(np.where(mic_neg > (mu + m_sigma * sigma)))
                for i in range(0, n_con_pts - 1):
                    idx = idx[np.diff(idx, append=-1) == 1]
                if idx.size:
                    if max(mic_neg) > max_intensity:
                        ret_time = rt_neg[idx[0]]
                        max_intensity = max(mic_neg)
                        best_mz = mz
                j = 0 
                for i in idx:
                    if i == j+1:
                        intensities[-1][0] += ((mic_neg[i] - mu)/sigma)
                    else:
                        intensities.append([0, float(rt_neg[i])])
                    j = i
                
                if intensities:
                    ret_time = max(intensities)[1]
                    if max(intensities)[0] > max_intensity:
                        max_intensity = max(mic_pos)
                        best_mz = mz
                        chromatogram = mic_neg
                        
    except IndexError:
        print(f"Issue with {mzxml_file}")
        ret_time = -1
        best_mz = 0
    
    if (0.6 > ret_time > 0) or (ret_time > 7.7):
        print(f"Peak detected too early/late to isolate ({ret_time})... possibly unrealistic ")
        ret_time = -1.0
    
    try:
        if ax:
            if ret_time > 0:
                ax.plot(rt_pos, chromatogram)
                ax.scatter(ret_time, max_intensity)
            else:
                plt.close()
    except:
        print("weird error, continuing")

    # time delay between PDA detector and MS is 0.05 minutes on our machine
    return (ret_time+0.05), best_mz

# Integrate the area of the PDA peak at the retention time corresponding to the target product
def pda_data_analysis(logE: dict, UV: float, ret_time: float, pda_txt: str, ax=None):
    """

    Parameters
    ----------
    logE : dict
        {'O': logE_h2o,
         'CC#N': logE_acn}
        (eventually change to [low_logE, hi_logE] and estimate C.I.)
    UV : float
        UVmax as predicted by chemprop
    ret_time : float
        retention time as calculated by MS
    pda_txt : str
        PDA text file of data

    Returns
    -------
    conc_info : dict
        {'conc': float,
         'chromatogram': np.array}
    """
    flow_rate = 0.5e-3/60
    inj_vol = 1e-6
    path_length = 1
    
    with open(pda_txt, 'r') as data_in:
        lines = data_in.readlines()
    for i, line in enumerate(lines):
        if len(line) > 200:
            break
    data = pd.read_csv(pda_txt, header=i-2, index_col=0)
    wl = [float(i)/100 for i in list(data.columns)]
    #remove baseline
    
    #baseline = np.ones(data.shape)
    for i, column in enumerate(data.columns):
        z = baselines.arpls(np.array(data[column]), lam=1e3)
        data[column] = data[column] - z[0]
    """
        baseline[:,i] = z[0]
    baseline = savgol_filter(baseline, 7, 2, axis=1, mode='nearest')
    for i, column in enumerate(data.columns):
        data[column] = data[column] - baseline[:,i]
    """
        
    # Assume that logE is the weighted avg of logE in pure water and ACN
    ### this is hard coded version of gradient! y = fraction acetonitrile
    y = min(max(95/6*(ret_time-0.5)+5,5),100)/100
    ext_coef = 10**(y*logE['CC#N'] + (1-y)*logE['O'])
    
    if ret_time < 0:
        print(f"No target peak detected in {pda_txt}")
        conc_info = {'conc': 0.0,
                     'chromatogram': np.array([wl, pd.DataFrame.sum(data[1:])]),
                     }
        return conc_info
    
    # just out of curiosity, compare the predicted UVvis to the measured one at the retention time
    # get max wavelength with prominence of XXXX at MS retention time
    ms_rt_idx = np.where(data.index<ret_time)[0][-1]   # get the index of the MS retention time
    wl_peaks, _ = sp.find_peaks(data.iloc[ms_rt_idx,:], prominence=2e4)   # get the max wavelengths at that ret time
    if not wl_peaks.any():
        wl_peaks = np.array([0])
    pda_peaks, _ = sp.find_peaks(data[data.columns[wl_peaks[-1]]], prominence=1e3)   
    
    # bound the hit peak at the wavelength of max absorption
    time = pd.Series(data.index)
    wl_sig = data[data.columns[wl_peaks[-1]]]
    match = np.argmin(abs(pda_peaks-ms_rt_idx))
    peak_width = sp.peak_widths(data[data.columns[wl_peaks[-1]]], pda_peaks, rel_height=0.98)
    left_bound = int(peak_width[2][match])
    right_bound = int(peak_width[3][match])+1
    
    # Calculate the baseline for the peak
    half_width = int((right_bound - left_bound)/2)
    try:
        left_bline = np.argmin(wl_sig.iloc[(left_bound-half_width):(left_bound)]) + left_bound - half_width
        right_bline = np.argmin(wl_sig.iloc[(right_bound):(right_bound+half_width)]) + right_bound
    except ValueError:
        print("Value error")
        left_bline = left_bound
        right_bline = right_bound
    slope = (wl_sig.iloc[right_bline] - wl_sig.iloc[left_bline])/(right_bline-left_bline)
    bline = [wl_sig.iloc[left_bline]+(i*slope)+(left_bound-left_bline)*slope 
             for i in range(len(time[left_bound:right_bound]))]
    signal_to_integrate = data[data.columns[wl_peaks[-1]]].iloc[left_bound:right_bound] - bline
    
    # calculate important values to print
    area = np.trapz(signal_to_integrate.values, signal_to_integrate.index)*60
    ist_area = internal_standard(data)    
    abs_units = area*1e-6
    conc = abs_units * flow_rate / (ext_coef * path_length * inj_vol)
    info_to_print = {'conc': conc,
                     'AU*s': area,
                     'logE': np.log10(ext_coef),
                     'istd': ist_area,
                     'Abs_max': data.columns[wl_peaks[-1]],
                     'ret_time': ret_time}
    
    # plot the peak that's being integrated
    if ax:
        plt.plot(time, wl_sig)
        plt.scatter(time.iloc[pda_peaks], wl_sig.iloc[pda_peaks])
        plt.fill_between(time[left_bound:right_bound], wl_sig.iloc[left_bound:right_bound], bline)
    
    # peak deconvolution
    peak = np.array(data.iloc[left_bound:right_bound, :])
    peak = peak#/peak.max()
    #init_spec = np.array((data.iloc[pda_peaks[match]], data.iloc[left_bound]))#/peak.max()
    sample_points = [pda_peaks[match], left_bound, right_bound, int(left_bound/2 + right_bound/2)]
    best_error = np.inf
    for i in range(len(sample_points)):
        init_spec = np.array((data.iloc[sample_points[0:i+1]]))
        mcrar = McrAR(max_iter=100, st_regr='NNLS', c_regr=OLS(), 
            c_constraints=[ConstraintNonneg(),])
        mcrar.fit(peak.reshape((-1, np.size(data.columns))), ST=init_spec)
        if mcrar.err[-1] < best_error:
            best_error = mcrar.err[-1]
            #print(f"error = {best_error}")
            comps = mcrar.C_opt_
            spectra = mcrar.ST_opt_
        else:
            break
    
    smoothed_spectra = savgol_filter(spectra, 19, 2, axis=1, mode='nearest')
    spec_lambda_max, _ = sp.find_peaks(smoothed_spectra[0,:], prominence=2e4)
    if not spec_lambda_max.any():
        spec_lambda_max = [0]
    new_area = np.trapz(comps[:,0], time.iloc[left_bound:right_bound])*spectra[0,spec_lambda_max[-1]]*60
    
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(wl, smoothed_spectra.T)
    ax1.scatter(wl[spec_lambda_max[-1]], max(spectra.T[:,0]))
    ax2.plot(time.iloc[left_bound:right_bound], comps)
    info_to_print['MCR_adjusted_area'] = new_area
    
    return info_to_print


def internal_standard(data):
    """
    This function just analyzes the last minute of the chromatogram to return the internal standard peak area
    """
    end = data.loc[7:8,:]
    
    index_of_257nm = np.argmin(abs(end.columns.astype(int)-25700))
    wl_sig = end.iloc[:,index_of_257nm]
    wl_sig = wl_sig - min(wl_sig)
    ist_peak, _ = sp.find_peaks(wl_sig, prominence=1e3)
    
    try:
        peak_width = sp.peak_widths(wl_sig, ist_peak, rel_height=0.98)
        left_bound = int(peak_width[2][0])
        right_bound = int(peak_width[3][0])
        signal_to_integrate = wl_sig.iloc[left_bound:right_bound]
        area = np.trapz(signal_to_integrate.values, signal_to_integrate.index)*60
    except IndexError:
        print("Unable to indetify internal standard peak, assuming no internal standard used!")
        area = 0
        
    return area


all_reaction_data = {}
for plate in plates:
    with open(os.path.join(plates_root, plate, f"{plate}.yaml"), 'r') as yaml_in:
        wellplate = yaml.load(yaml_in, Loader=yaml.Loader)
    chemprop_file = log_e_chemprop_predict(wellplate['reaction_plate'], plate)
    all_reaction_data[plate] = lc_wellplate_analysis(wellplate['reaction_plate'], 
                                                     chemprop_file, 
                                                     os.path.join(plates_root, plate))

#TODO! Turn this into an easy to read csv because YAML hates numpy
with open(results_file, 'w') as yaml_out:
    yaml.dump(all_reaction_data, yaml_out)
    