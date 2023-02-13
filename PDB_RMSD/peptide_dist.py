#!/ usr / bin / env python
import math
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
import pandas as pd

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Similarity by Residue of Different SH2-Paptide Combos')
parser.add_argument('-l', required=True, help='List of PDB name, SH2 residues, peptide residues (line 1 is reference)')
parser.add_argument('-f', required=True, help= 'File extension for output')

#Import Arguments
args = parser.parse_args()
File_input = args.l
if File_input.split('.')[-1] != 'txt': #Add file extension if not in input
    File_input = File_input + '.txt'
file_ext = args.f

#Load reference
input_data = open(File_input, 'r').readlines()
[ref, res_sh2_ref, res_pep_minus_ref, res_pep_plus_ref] = input_data[0].split()
res_pep_ref = int(res_pep_plus_ref) + int(res_pep_minus_ref) + 1
ref_pdb = md.load_pdb('../PDB/' + ref + '.pdb')

#Identify location of PTR
PTR_ref = int(res_sh2_ref) + int(res_pep_minus_ref)
PTR_ref_check = ref_pdb.topology.select('resid ' + str(PTR_ref) + ' and resname PTR')
if len(PTR_ref_check) == 0:
    sys.exit('Error reference PTR not identified properly')

#Set up alignment indices
#align_ind_ref = ref_pdb.topology.select('(resid > ' + str(res_sh2_ref_start) + ' and resid < ' + str(res_sh2_ref_end) + ' and backbone) or resid ' + str(PTR_ref))
align_ind_ref = ref_pdb.topology.select('resid ' + str(PTR_ref) + ' and backbone')

#Declare empty dataframe
df = pd.DataFrame()

#Loop through input file
Dis_all = []
for line in range(1, len(input_data)):
    #Load data from line
    [pdb, res_sh2, res_pep_minus, res_pep_plus] = input_data[line].split()

    #Load pdb file
    target_pdb = md.load_pdb('1qg1_align/' + pdb + '_align.pdb')
    res_pep = int(res_pep_minus)+int(res_pep_plus)+1

    #Identify location of PTR
    PTR = int(res_sh2) + int(res_pep_minus)
    
    #Check PTR is identified properly
    PTR_check = target_pdb.topology.select('resid ' + str(PTR) + ' and resname PTR')
    if len(PTR_check) == 0:
        sys.exit('Error PTR not identified properly for ' + pdb)

    #Set up alignment indices    
#    align_ind_target = target_pdb.topology.select('(resid > ' + str(res_sh2_start) + ' and resid < ' + str(res_sh2_end) + ' and backbone) or resid ' + str(PTR))
    align_ind_target = target_pdb.topology.select('resid ' + str(PTR) + ' and backbone')

    #Align on SH2 domain and PTR
    target_align = target_pdb.superpose(ref_pdb, atom_indices = align_ind_target, ref_atom_indices = align_ind_ref)
#    ref_pdb.save('temp_ref.pdb')
#    target_align.save('temp.pdb')

    #Set residue array for +/- in target
    res = np.linspace(PTR - int(res_pep_minus), PTR + int(res_pep_plus), res_pep)
    res_ref = np.linspace(PTR_ref - int(res_pep_minus), PTR_ref + int(res_pep_plus), res_pep)

    #Loop through SH2 residues to determine RMSD for each
    RMSD = []
    for i in range(int(res_pep)):
        pep_ref = ref_pdb.topology.select('resid ' + str(res_ref[i]) + ' and backbone')
        pep_target = target_align.topology.select('resid ' + str(res[i]) + ' and backbone')
        
        #Get coordinates of backbone atoms for reference
        ref_corr = ref_pdb.xyz[0, pep_ref]

        #Determine frames in target PDB
        frames = target_align.n_frames
        res_rmsd = []
        for frame in range(frames):
            target_corr = target_align.xyz[0, pep_target]
        
            #Compute squared displacment for all backbone atoms
            displacment = np.zeros(4)
            for n in range(4):
                displacment[n] = (ref_corr[n][0] - target_corr[n][0])**2 + (ref_corr[n][1] - target_corr[n][1])**2 + (ref_corr[n][2] - target_corr[n][2])**2
            res_rmsd.append(math.sqrt(np.mean(displacment)))
        RMSD.append(res_rmsd)

    if len(res_rmsd) == 1:
        if int(res_pep_plus) < int(res_pep_plus_ref):
            for i in range(int(res_pep_plus_ref) - int(res_pep_plus)):
                RMSD.append(-1)
        if int(res_pep_minus) < int(res_pep_minus_ref):
            for i in range(int(res_pep_minus_ref) - int(res_pep_minus)):
                RMSD.insert(0, -1)
        df[pdb] = RMSD
    else:
        #Reformat RMSD to np array
        RMSD_reformat = np.full((len(RMSD[0]), res_pep_ref), -1, dtype=float)
        n = int(res_pep_minus_ref) - int(res_pep_minus)
        for i in range(n, len(RMSD)):
            RMSD_reformat[:,i] = RMSD[i]
        for i in range(len(RMSD[0])):
            df[pdb + ' ' + str(i)] = RMSD_reformat[i,:]

df.index = ['-4', '-3', '-2', '-1', '0', '+1', '+2', '+3', '+4', '+5', '+6', '+7', '+8']
df.to_csv('rmsd_' + file_ext + '.csv')
