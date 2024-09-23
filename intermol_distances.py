
import os 
import itertools
from itertools import product
import numpy as np 
import scipy
from scipy.spatial import distance
import pickle
import pandas as pd
import json

directory = "C:/Users/ezraa/OneDrive/Desktop/PURA/PURA folded dimers/"
os.chdir(directory)

def get_df(name):
    array = loaded_named_lists[name]
    df =  pd.DataFrame([s.split(',') for s in array], columns = [f"atom","resi", "x", "y","z","chain"])
    df.iloc[:, 2:5] = df.iloc[:,2:5].apply(pd.to_numeric)
    # print(df.iloc[:,2:6])
    df = df[df[f'atom']=='CA']
    df['resi'] = df['resi'].astype(float)
    # df['original'] = name
    if df['resi'].duplicated().any():
        df = df.groupby(['resi','atom']).mean().reset_index()
        # print(f'{name} has duplicates')
    df.set_index('resi', inplace = True, drop = False)

    # print(f'got df {name}')

    return(df)

def get_mutant_df(name):
    array = loaded_named_lists[name]
    df =  pd.DataFrame([s.split(',') for s in array], columns =[f"atom","resi", "x", "y","z","chain"])
    df.iloc[:, 2:5] = df.iloc[:,2:5].apply(pd.to_numeric)
    df = df[df[f'atom']=='CA']
    # if df.shape[0] <  322:
    #     rows_to_add = max(0, 322 - df.shape[0])
    #     empty_rows = pd.DataFrame(np.nan, index=range(rows_to_add), columns=df.columns)
    #     df = pd.concat([df, empty_rows], ignore_index=True)
    # df['original'] = name
    df['resi'] = df['resi'].astype(float)
    df['xyz'] =  df[['x', 'y', 'z']].apply(lambda row: row.values.flatten(), axis=1)
    df.set_index('resi', inplace = True, drop = False)

    # print(f'got df {name}')
    
    return(df)

def row_function(row):
    # distances = []
    if row['chain'] == 'A':
        mutant_coords = row[['x','y','z']].values.flatten()
        residue = row['resi']

        if  not np.isnan(residue):
                to_compare_row = mutant_df[(mutant_df['resi']==residue) &( mutant_df['chain']=='B')]
                # print(to_compare_row)
                to_compare_coords = to_compare_row[['x','y','z']].values.flatten()
                # print(to_compare_coords)

                if to_compare_coords.size > 0:
                    dist = distance.euclidean(mutant_coords, to_compare_coords)
                
                else:
                    dist = ''
            
                # distances.append(dist)
    else:
        dist = ''
                # print(f'{residue} missing')
    return(dist)

i = 0

for root, dirs, files in os.walk(directory):
    # Ensure root is correctly normalized
    root = os.path.normpath(root)
    for file in files:
        # Construct full file path
        # print(file)
        if file.endswith(".pkl"):
            # print(file)
            # print(f'starting {file}')
            file_path = os.path.join(root, file)
            # print(f'path = {file_path}, file = {file}')

            name = os.path.splitext(file)[0]

            with open(file_path, 'rb') as file:
                loaded_named_lists = pickle.load(file)
            
            mutant_df= get_mutant_df(f"{name}_array")   

            a = mutant_df.loc[mutant_df['chain']=='A',['x','y','z']].to_numpy(dtype=float)

            b = mutant_df.loc[mutant_df['chain']=='B',['x','y','z']].to_numpy(dtype=float)

            dist = distance.cdist(a,b)

            csv_name = os.path.join(os.path.dirname(file_path),f'{name}_intermol_distances.csv')
            df = pd.DataFrame(dist, index=[f'Point_A{i}' for i in range(len(a))],
                  columns=[f'Point_B{i}' for i in range(len(b))])
            df.to_csv(csv_name, index=False)

            i = i + 1
            print(f'{i}')

            print(name)





            
            


