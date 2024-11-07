#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 08:30:19 2021

@author: mattmcdonald
"""

import glob
import xmltodict
from collections import OrderedDict
import collections
import sys
import math
import os
import csv
from pprint import pprint
from rdkit import Chem


# use to extract UVvis data as well
def get_eac_val(entry):
    if type(entry) == list:
        strengths_list = []
        for item in entry:
            try:
                if type(item) in [dict, collections.OrderedDict]:
                    try:
                        eac_val = float(item['UV.EAC']['#text'])
                    except:
                        eac_val = float(item['UV.EAC'])
                else:
                    eac_val = float(item['UV.EAC'])
                strengths_list.append(eac_val)
            except KeyError:
                try:
                    loge_val = float(item['UV.LOGE'])
                    strengths_list.append(10**loge_val)
                except KeyError:
                    strengths_list.append(-1)
        return_val = strengths_list
    elif type(entry) in [dict, collections.OrderedDict]:
        strengths_list = []
        try:
            if type(entry['UV.EAC']) in [dict, collections.OrderedDict]:
                eac_val = float(entry['UV.EAC']['#text'])
            else:
                eac_val = float(entry['UV.EAC'])
            strengths_list.append(eac_val)
        except KeyError:
            try:
                loge_val = float(entry['UV.LOGE'])
                strengths_list.append(10**loge_val)
            except KeyError:
                strengths_list.append(-1)
        return_val = strengths_list
    else:
        return_val = [-1]
        
    return return_val


def get_abs_max(entry):
    try:
        if type(entry) == list:
            ab_max_list = []
            for item in entry:
                if type(item) in [dict, collections.OrderedDict]:
                    try:
                        ab_max_list.append(float(item['UV.AM']['#text']))
                    except:
                        ab_max_list.append(float(item['UV.AM']))
                else:
                    ab_max_list.append(float(item['UV.AM']))
            return_val = ab_max_list
        elif type(entry) in [dict, collections.OrderedDict]:
            if type(entry['UV.AM']) in [dict, collections.OrderedDict]:
                return_val = [float(entry['UV.AM']['#text'])]
            else:
                return_val = [float(entry['UV.AM'])]
        else:
            return_val = [-1]
    except KeyError:
        raise
    return return_val


def recursive_ordered_dict_to_dict(ordered_dict):
    simple_dict = {}
    for key, value in ordered_dict.items():
        if isinstance(value, OrderedDict):
            simple_dict[key] = recursive_ordered_dict_to_dict(value)
        else:
            simple_dict[key] = value
    return simple_dict


num_args = len(sys.argv)
if num_args != 3:
    print(f"Incorrect number of arguments, must provide xml_folder and output.csv, but {num_args-1} args found")
    print(sys.argv)
    sys.exit()

folder = sys.argv[1]
output = sys.argv[2]
files = glob.glob(os.path.join(folder,'*.xml'))

# Dictionary to convert solvents (as referenced by reaxys) into smiles for chemprop
solvents = {'acetonitrile': 'CC#N',
            'MeCN': 'CC#N',
            'ethanol': 'CCO',
            'water': 'O',
            'H2O': 'O',
            'dichloromethane': 'ClCCl',
            'toluene': 'Cc1ccccc1',
            'methanol': 'CO',
            '2-methyl-propan-2-ol': 'CC(O)(C)C',
            'isopropyl alcohol': 'CC(O)C',
            'benzene': 'c1ccccc1',
            '1,2-dichloro-benzene': 'Clc1c(Cl)cccc1',
            '1,4-dioxane': 'C1OCCOC1',
            'dimethyl sulfoxide': 'CS(=O)C',
            'dimethylsulfoxide': 'CS(=O)C',
            'chloroform': 'ClC(Cl)Cl',
            'CHCl3': 'ClC(Cl)Cl',
            'dimethylformamide': 'CN(C=O)C',
            'N,N-dimethyl-formamide': 'CN(C=O)C',
            '1-methyl-pyrrolidin-2-one': 'CN1CCCC1=O',
            'hexane': 'CCCCCC',
            'cyclohexane': 'C1CCCCC1',
            'ethyl acetate': 'O=C(OCC)C',
            }

molecule_dict = {}
molecule_index = 0
file_index = 1
print("Importing Data")
for file in files:
    print('Working on %s (%s of %s)' % (file, file_index, len(files)))
    file_index += 1
    with open(file, 'r') as infile:
        xml_dict = xmltodict.parse(infile.read())
    xml_dict['xf'].pop('response', None)
    
    for item in xml_dict['xf']['substances']['substance']:
        if molecule_index not in molecule_dict.keys():
            molecule_dict[molecule_index] = {}
        sub_length = len(list(molecule_dict[molecule_index].keys()))
        molecule_dict[molecule_index][sub_length] = recursive_ordered_dict_to_dict(item)
        molecule_index += 1

rdkit_mdl_header = '\n     RDKit          2D\n'

molecule_data_export = {}
repeats = {}
index = 1
for molecule in molecule_dict.keys():
    for mol_entry in molecule_dict[molecule].keys():
        try:
            mdlfile_string = molecule_dict[molecule][mol_entry]['YY'][1]['YY.STR']['#text']
        except KeyError:
            try:
                mdlfile_string = molecule_dict[molecule][mol_entry]['YY']['YY.STR']['#text']
            except KeyError:
                pprint(molecule)
                continue
        temp_split = mdlfile_string.split('\n')
        cleaned_lines = []
        for element in temp_split:
            if 'M  ' in element:
                if 'M  END' in element:
                    cleaned_lines.append(element)
            else:
                cleaned_lines.append(element)
        
        mdlfile_string = mdlfile_string.replace('HDR', rdkit_mdl_header)
        mdlmol = Chem.MolFromMolBlock(mdlfile_string, sanitize = True)
        try:
            mdlsmiles = Chem.MolToSmiles(mdlmol)
        except:
            continue

        molecule_data_export[mdlsmiles] = {} #{'fingerprint': fingerprint_string}

        try:
            reaxys_id = molecule_dict[molecule][mol_entry]['IDE']['IDE.XRN']
        except KeyError:
            print('Molecule does not have reaxys id?')

        try:
            uv_data_entries = molecule_dict[molecule][mol_entry]['UV']
        except KeyError:
            pprint(f"No UV data found for {molecule_dict[molecule][mol_entry]}")
            uv_data_entries = 'Not Found'
        
        if uv_data_entries != 'Not Found':
            uv_dataset = []
            if type(uv_data_entries) in [dict, collections.OrderedDict]:
                try:
                    abs_max = get_abs_max(uv_data_entries['UV01'])
                except KeyError:
                    abs_max = ['None']
                    continue
                abs_strengths = get_eac_val(uv_data_entries['UV01'])
                try:
                    if type(uv_data_entries['UV.SOL']) in [dict, collections.OrderedDict]:
                        solvent = [uv_data_entries['UV.SOL']['hi']]
                    else:
                        solvent = [uv_data_entries['UV.SOL']]
                except KeyError:
                    solvent = ['None']
                uv_dataset.append([abs_max, abs_strengths, solvent, reaxys_id])
                molecule_data_export[mdlsmiles]['UV'] = uv_dataset
                
            elif type(uv_data_entries) == list:
                for uv_entry in uv_data_entries:
                    try:
                        abs_max = get_abs_max(uv_entry['UV01'])
                    except KeyError:
                        abs_max = ['None']
                        continue
                    abs_strengths = get_eac_val(uv_entry['UV01'])
                    try:
                        if type(uv_entry['UV.SOL']) in [dict, collections.OrderedDict]:
                            solvent = [uv_entry['UV.SOL']['#text']]
                        else:
                            solvent = [uv_entry['UV.SOL']]
                    except KeyError:
                        solvent = ['None']
                    uv_dataset.append([abs_max, abs_strengths, solvent, reaxys_id])
                molecule_data_export[mdlsmiles]['UV'] = uv_dataset 
            else:
                print(uv_data_entries)
                sys.exit()
    index += 1

# pare down to only data we want in csv
ext_coef_dict = {}
uv_max_dict = {}
for smiles in molecule_data_export.keys():
    try:
        smiles_split = smiles.split('.')
        if len(smiles_split) > 2:
            pass
        elif len(smiles_split) == 2:
            if len(smiles_split[0]) > len(smiles_split[1]):
                counter_ion = smiles_split[1]
                mol = Chem.MolFromSmiles(counter_ion)
                charge = Chem.GetFormalCharge(mol)
                if charge == 1:
                    counter_ion = ".[Na+]"
                elif charge == -1:
                    counter_ion = ".[Cl-]"
                else:
                    counter_ion = ""
                target = smiles_split[0]+counter_ion
            else:
                counter_ion = smiles_split[0]
                mol = Chem.MolFromSmiles(counter_ion)
                charge = Chem.GetFormalCharge(mol)
                if charge == 1:
                    counter_ion = ".[Na+]"
                elif charge == -1:
                    counter_ion = ".[Cl-]"
                else:
                    counter_ion = ""
                target = smiles_split[1]+counter_ion
        else:
            target = smiles
        
        try:
            solvent_name = molecule_data_export[smiles]['UV'][0][2][0]
            if type(solvent_name) is list:
                # these measuremtnes made in mixed solvents!
                pass
            else:
                solvent_smiles = solvents[solvent_name]
        except KeyError:
            # These measurements were taken in odd solvents!
            continue
        max_uv = max(molecule_data_export[smiles]['UV'][0][0])
        max_index = molecule_data_export[smiles]['UV'][0][0].index(max_uv)
        ext_coef = molecule_data_export[smiles]['UV'][0][1][max_index]
        """
        reaxys_id = molecule_data_export[smiles]['UV'][0][3]
        if (ext_coef > 1e1) and (ext_coef < 1e7) and (max_uv < 1000):
            if smiles not in ext_coef_dict.keys():
                ext_coef_dict[smiles] = [[target, solvent_smiles, math.log10(ext_coef), max_uv, reaxys_id]]
            else:
                ext_coef_dict[smiles].append([target, solvent_smiles, math.log10(ext_coef), max_uv, reaxys_id])
        """
        if (ext_coef > 1e1) and (ext_coef < 1e7) and (max_uv < 1000):
            if smiles not in ext_coef_dict.keys():
                ext_coef_dict[smiles] = [[target, solvent_smiles, math.log10(ext_coef), max_uv]]
            else:
                ext_coef_dict[smiles].append([target, solvent_smiles, math.log10(ext_coef), max_uv])
    except KeyError:
        print(f"{smiles} doesn't have UV data")


#data_columns = ['SMILES', 'SOLVENT', 'logE', 'UVmax', 'Reaxys Registry Number']
data_columns = ['SMILES', 'SOLVENT', 'logE', 'UVmax']
with open(output, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(data_columns)
    for entries in ext_coef_dict.values():
        for entry in entries: 
            writer.writerow(entry)

print(f"\n\n\n\tCompleted extracting and validating {len(ext_coef_dict)} molecules from {len(files)} files\n\n\n")
