
import os 
import itertools
from itertools import product
import numpy as np 
import scipy
from scipy.spatial import distance
import pickle
import pandas as pd

directory = "C:/Users/ezraa/OneDrive/Desktop/PURA/PURA folded mutations/"
os.chdir(directory)

def get_df(name):
    array = loaded_named_lists[name]
    df =  pd.DataFrame([s.split(',') for s in array], columns = [f"atom","resi", "x", "y","z"])
    df.iloc[:, 2:6] = df.iloc[:,2:6].apply(pd.to_numeric)
    # print(df.iloc[:,2:6])
    df = df[df[f'atom']=='CA']
    df['resi'] = df['resi'].astype(float)
    # df['original'] = name
    if df['resi'].duplicated().any():
        df = df.groupby(['resi','atom']).mean().reset_index()
        # print(f'{name} has duplicates')
    df.set_index('resi', inplace = True, drop = False)
    return(df)

def get_mutant_df(name):
    array = loaded_named_lists[name]
    df =  pd.DataFrame([s.split(',') for s in array], columns =[f"atom","resi", "x", "y","z"])
    df.iloc[:, 2:6] = df.iloc[:,2:6].apply(pd.to_numeric)
    # print(df.iloc[:,2:6])
    df = df[df[f'atom']=='CA']
    if df.shape[0] <  322:
        rows_to_add = max(0, 322 - df.shape[0])
        empty_rows = pd.DataFrame(np.nan, index=range(rows_to_add), columns=df.columns)
        df = pd.concat([df, empty_rows], ignore_index=True)
    # df['original'] = name
    df['resi'] = df['resi'].astype(float)
    df.set_index('resi', inplace = True, drop = False)
    return(df)

def row_function(row,to_compare):
    mutant_coords = row[['x','y','z']].values.flatten()
    residue = row['resi']
    distances = []

    if residue in to_compare['resi'] and not np.isnan(residue):
            to_compare_row = to_compare[to_compare['resi']==residue]
            # print(to_compare_row)
            to_compare_coords = to_compare_row[['x','y','z']].values.flatten()
            # print(to_compare_coords)
            dist = distance.euclidean(mutant_coords, to_compare_coords)
            distances.append(dist)
    else:
        distances.append('')
            # print(f'{residue} missing')
    return(distances)

for root, dirs, files in os.walk(directory):
    # Ensure root is correctly normalized
    root = os.path.normpath(root)
    for file in files:
        # Construct full file path
        # print(file)
        if file.endswith(".pkl"):
            file_path = os.path.join(root, file)
            # print(f'path = {file_path}, file = {file}')

            name = os.path.splitext(file)[0]

            with open(file_path, 'rb') as file:
                loaded_named_lists = pickle.load(file)
            # for item in loaded_named_lists:
                # print(item)

            _8chw_df = get_df('_8chw_array')
            _8cht_df = get_df('_8cht_array')
            fold_wt_pura_model_0_df = get_df('fold_wt_pura_model_0_array')
            mutant_df= get_mutant_df(f"{name}_array")   

            _8chw_comparison = mutant_df.apply(lambda row: row_function(row, _8chw_df), axis=1)
            _8cht_comparison = mutant_df.apply(lambda row: row_function(row, _8cht_df), axis=1)
            fold_wt_pura_model_0_comparison = mutant_df.apply(lambda row: row_function(row, to_compare = fold_wt_pura_model_0_df), axis=1)

            mutant_df['_8chw_distances'] = _8chw_comparison
            mutant_df['_8cht_distances'] = _8cht_comparison
            mutant_df['fold_wt_pura_model_0_distances'] = fold_wt_pura_model_0_comparison

            csv_name = os.path.join(os.path.dirname(file_path),f'{name}.csv')
            mutant_df.to_csv(csv_name, index=False)

            print(name)

            # print(mutant_df)




            
            


