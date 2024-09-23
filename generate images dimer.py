import pymol
from pymol import cmd
import os 
import itertools
from itertools import product
import numpy as np 
import glob

directory = r"C:\Users\ezraa\OneDrive\Desktop\PURA\small folder dimer"

os.chdir(directory)
# pymol.finish_launching()

cmd.bg_color("white")
cmd.clip('near',1000000)
cmd.set('fog', 0)

view = (-0.23348447680473328,
 0.9462260007858276,
 -0.22391505539417267,
 0.6936360597610474,
 0.3234599828720093,
 0.6436133980751038,
 0.6814343929290771,
 -0.005041631869971752,
 -0.7318618893623352,
 0.0,
 0.0,
 -407.23876953125,
 -0.5424041748046875,
 0.19318389892578125,
 -0.38320159912109375,
 321.07012939453125,
 493.40740966796875,
 -20.0)

A_base_color = "0xECCBAE"
A_disordered_color = "0xF98400"
A_I_color = "0xFF0000"
A_II_color = "0x046C9A"
A_III_color = "0x00A08A"

B_base_color = "0xF6F0EC"
B_disordered_color = "0xF0DCC5"
B_I_color = "0xF0C7C7"
B_II_color = "0xB2D9E9"
B_III_color = "0xB1EAE3"

c = {
    'A':[A_base_color,A_disordered_color,A_I_color,A_II_color,A_III_color],
    'B':[B_base_color,B_disordered_color,B_I_color,B_II_color,B_III_color]
}

def color_function(name,chain):
    cmd.color(c[chain][0], f'{name} and chain {chain}')
    count_resi = pymol.cmd.count_atoms(f"{name} and polymer and chain {chain} and name CA")
    if count_resi > 54:
                cmd.select("disordered1", f"{name} and chain {chain} and resi 1-55")
                cmd.color(c[chain][1], "disordered1")
    else:
                cmd.select("disordered1", f"{name} and chain {chain} and resi 1-{count_resi}")
                cmd.color(c[chain][1], "disordered1")
            
    if count_resi > 321:
                cmd.select("disordered2", f"{name} and chain {chain} and resi 295-322")
                cmd.color(c[chain][1], "disordered2")
    elif count_resi > 295:
                cmd.select("disordered2", f"{name} and chain {chain} and resi 295-{count_resi}")
                cmd.color(c[chain][1], "disordered2")

    if count_resi > 124:
                cmd.select("PUR_I", f"{name} and chain {chain} and resi 60-125")
                cmd.color(c[chain][2], "PUR_I")
    elif count_resi > 60:
                cmd.select("PUR_I", f"{name} and chain {chain} and resi 60-{count_resi}")
                cmd.color(c[chain][2], "PUR_I")
            
    if count_resi > 212:
                cmd.select("PURII",f"{name} and chain {chain} and resi 142-213")
                cmd.color(c[chain][3], "PURII")
    elif count_resi > 142:
                cmd.select("PURII", f"{name} and chain {chain} and resi 142-{count_resi}")
                cmd.color(c[chain][3], "PURII")
            
    if count_resi > 280:
                cmd.select("PUR_III", f"{name} and chain {chain} and resi 215-281")
                cmd.color(c[chain][4], "PUR_III")
    elif count_resi > 215:
                cmd.select("PUR_III", f"{name} and chain {chain} and resi 215-{count_resi}")
                cmd.color(c[chain][4], "PUR_III")
    cmd.deselect()


for dirpath, dirnames, filenames in os.walk(directory):
    # Iterate over all files in the current directory
    for filename in filenames:
        # Check if the file ends with .cif
        if filename.endswith(".cif"):
            cif_file_path = os.path.join(dirpath, filename)
            name = os.path.splitext(os.path.basename(cif_file_path))[0]
            subdir = os.path.join(dirpath,name)

            if os.path.exists(os.path.join(subdir,f"{name}.png")):
                    print(f"{os.path.join(subdir,f"{name}.png")} exists")
            else:
                try:
                    os.mkdir(os.path.join(dirpath,name))
                except:
                    pass

                print(name)

                # print(count_resi)

                cmd.load(cif_file_path)
                cmd.load(r"C:\Users\ezraa\OneDrive\Desktop\PURA\fold_pura_wt_wt\fold_pura_wt_wt_model_0.cif")

                cmd.align(name, "fold_pura_wt_wt_model_0")
                cmd.set_view(view)
                cmd.delete("fold_pura_wt_wt_model_0")

                color_function(name,'A')
                color_function(name,'B')


                cmd.png(os.path.join(subdir,f"{name}.png"), width=1080, height=1080, dpi=300, ray=1)
                
                # print(os.path.join(subdir,f"{name}.png"))

                cmd.delete(name)


            