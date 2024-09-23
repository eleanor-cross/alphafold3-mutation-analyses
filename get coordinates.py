import pymol
from pymol import cmd
import os 
import itertools
from itertools import product
import numpy as np 
from Bio import PDB
from ase import Atoms
from ase.io import read 
import pickle 

directory = "C:/Users/ezraa/OneDrive/Desktop/PURA/PURA folded mutations/"
os.chdir(directory)
pymol.finish_launching()

for dirpath, dirnames, filenames in os.walk(directory):
    # Iterate over all files in the current directory
    for filename in filenames:
        # Check if the file ends with .cif
        if filename.endswith(".cif"):
            name =  filename.replace(".cif","")
            subdir = os.path.join(dirpath,name)
            try:
                os.mkdir(subdir)
            except FileExistsError:
                pass
                # print(f"Directory '{directory}' already exists")

            ## load in everything
            pymol.finish_launching()

            cif_file_path = os.path.join(dirpath, filename)
            cmd.load(cif_file_path,name)

            subdir = os.path.normpath(subdir)
            # print(f'path = {subdir}')
            
            # mutant
            cmd.load(r"C:\Users\ezraa\OneDrive\Desktop\PURA\fold_wt_pura\fold_wt_pura_model_0.cif")
            # model 0 for wt 
            cmd.load(r"C:\Users\ezraa\OneDrive\Desktop\PURA\8cht.cif", "8cht") 
            #Crystal structure of human PURA (fragment Glu57-Glu212, PUR repeat I and II)
            # cmd.get_chains('8cht')
            cmd.select('chain_to_keep','8cht and chain C')
            cmd.remove('8cht and not chain_to_keep')
            # this structure has 4 identical chains for some reason 
            # also some have multiple alpha carbons for a given residue, which shouldn't be happening 
            # select chain C because it has the least errors
            cmd.load(r"C:\Users\ezraa\OneDrive\Desktop\PURA\8chw.cif", "8chw") 
            #Crystal structure of human PURA (fragment Pro216-Lys280, PUR repeat III)

            ## alignments 
            cmd.align("fold_wt_pura_model_0",name)
            cmd.align("8chw", name)
            cmd.align("8cht", name)

            # print(count_resi)

            all_arrays = {
                '_8chw_array' : [f"{atom.name},{atom.resi},{atom.coord[0]},{atom.coord[1]},{atom.coord[2]}" for atom in cmd.get_model('8chw').atom],
                '_8cht_array' : [f"{atom.name},{atom.resi},{atom.coord[0]},{atom.coord[1]},{atom.coord[2]}" for atom in cmd.get_model('8cht').atom],
                'fold_wt_pura_model_0_array' : [f"{atom.name},{atom.resi},{atom.coord[0]},{atom.coord[1]},{atom.coord[2]}" for atom in cmd.get_model('fold_wt_pura_model_0').atom],
                f'{name}_array' : [f"{atom.name},{atom.resi},{atom.coord[0]},{atom.coord[1]},{atom.coord[2]}" for atom in cmd.get_model(name).atom]
            }

            # print(subdir)

            with open(f'{subdir}\{name}.pkl','wb') as file:
                pickle.dump(all_arrays, file)
                # print(f'wrote {subdir}{name}.pkl')

            cmd.delete(name)
            cmd.delete('8chw')
            cmd.delete('8cht')
            cmd.delete('fold_wt_pura_model_0')
            


