import pymol
from pymol import cmd
import os 
import itertools
from itertools import product
import numpy as np 
import glob

directory = "C:/Users/ezraa/OneDrive/Desktop/PURA/PURA folded mutations/"

os.chdir(directory)
pymol.finish_launching()

view =  (-0.148281127,    0.963571966,    0.222553030,
     0.866201699,    0.017958345,    0.499366224,
     0.477183700,    0.266826838,   -0.837314010,
    -0.000097036,    0.000026226, -407.238769531,
    30.573596954,   -7.856406212,   12.066558838,
   321.070007324,  493.407287598,  -20.000000000 )

cmd.bg_color("white")
cmd.clip('near',1000000)
cmd.set('fog', 0)

base_color = "0xECCBAE"
disordered_color = "0xF98400"
I_color = "0xFF0000"
II_color = "0x046C9A"
III_color = "0x00A08A"


for dirpath, dirnames, filenames in os.walk(directory):
    # Iterate over all files in the current directory
    for filename in filenames:
        # Check if the file ends with .cif
        if filename.endswith(".cif"):
            cif_file_path = os.path.join(dirpath, filename)
            name = os.path.splitext(os.path.basename(cif_file_path))[0]
            print(name)

            # print(count_resi)

            cmd.load(cif_file_path)
            
            cmd.load(r"C:\Users\ezraa\OneDrive\Desktop\PURA\fold_wt_pura\fold_wt_pura_model_0.cif")
            cmd.set_view(view)

            cmd.align(name, "fold_wt_pura_model_0")
            cmd.delete("fold_wt_pura_model_0")

            count_resi = pymol.cmd.count_atoms(f"{name} and polymer and name CA")

            cmd.color(base_color, name)

            if count_resi > 54:
                cmd.select("disordered1", "resi 1-55")
                cmd.color(disordered_color, "disordered1")
            else:
                cmd.select("disordered1", f"resi 1-{count_resi}")
                cmd.color(disordered_color, "disordered1")
            
            if count_resi > 321:
                cmd.select("disordered2", "resi 295-322")
                cmd.color(disordered_color, "disordered2")
            elif count_resi > 295:
                cmd.select("disordered2", f"resi 295-{count_resi}")
                cmd.color(disordered_color, "disordered2")

            if count_resi > 124:
                cmd.select("PUR_I", "resi 60-125")
                cmd.color(I_color, "PUR_I")
            elif count_resi > 60:
                cmd.select("PUR_I", f"resi 60-{count_resi}")
                cmd.color(I_color, "PUR_I")
            
            if count_resi > 212:
                cmd.select("PURII","resi 142-213")
                cmd.color(II_color, "PURII")
            elif count_resi > 142:
                cmd.select("PURII", f"resi 142-{count_resi}")
                cmd.color(II_color, "PURII")
            
            if count_resi > 280:
                cmd.select("PUR_III", "resi 215-281")
                cmd.color(III_color, "PUR_III")
            elif count_resi > 215:
                cmd.select("PUR_III", f"resi 215-{count_resi}")
                cmd.color(III_color, "PUR_III")
            
            cmd.deselect()

            cmd.png(f"{cif_file_path}.png", width=1080, height=1080, dpi=300, ray=1)

            cmd.delete(name)





