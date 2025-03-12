import numpy as np
from ase import Atoms
from ase.visualize import view
import pickle
#this script visualize data

# Load data
data_root="./mcDataAll/Nx4_Ny4_Nz4/T100/P0.05/U_dist_dataFiles/"
coord_folder=data_root+"/coord/"

coord_init=coord_folder+"/init.coord.pkl"

raw_init=data_root+"/raw.pkl"
box_file=data_root+"/box.pkl"


with open(raw_init,"rb") as fptr:
    raw_arr=pickle.load(fptr)

with open(coord_init,"rb") as fptr:
    coord_arr=pickle.load(fptr)
with open(box_file,"rb") as fptr:
    box_arr=pickle.load(fptr)

n_atoms = len(raw_arr)

type_map=["O","H"]

symbols = [type_map[t] for t in raw_arr]
positions = coord_arr.reshape(n_atoms, 3)  # Critical reshape
atoms = Atoms(symbols=symbols,
              positions=positions,
              cell=box_arr.reshape(3,3),
              pbc=True)

view(atoms)