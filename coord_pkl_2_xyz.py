import numpy as np
import pickle
from decimal import Decimal, getcontext
import glob
import re
#this script loads pkl data for coord
# and converts to pkl file for visualization


def format_using_decimal(value, precision=8):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

N=4
Nx=N
Ny=N
Nz=N


T=100

P=0.1

TStr=format_using_decimal(T)

PStr=format_using_decimal(P)
total_atom_xyz_components_num=Nx*Ny*Nz*3*3
total_atom_num=Nx*Ny*Nz*3
total_molecule_num=Nx*Ny*Nz
#read all or some of the pkl files
n_dig=8
def sort_data_files_by_flushEnd(oneDir):
    dataFilesAll=[]
    flushEndAll=[]
    for oneDataFile in glob.glob(oneDir+"/flushEnd*.pkl"):
        dataFilesAll.append(oneDataFile)
        matchEnd=re.search(r"flushEnd(\d+)",oneDataFile)
        if matchEnd:
            indTmp=int(matchEnd.group(1))
            flushEndAll.append(indTmp)

    endInds=np.argsort(flushEndAll)
    sortedDataFiles=[dataFilesAll[i] for i in endInds]
    return sortedDataFiles

def read_one_pkl_file(one_pkl_file_name):
    with open(one_pkl_file_name,"rb") as fptr:
        arr=np.array(pickle.load(fptr))
    arr=arr.reshape((-1,total_atom_xyz_components_num))
    return arr

coord_pkl_path=f"./mcDataAll/Nx{N}_Ny{N}_Nz{N}/T{TStr}/P{PStr}/U_dist_dataFiles/coord/"
sorted_coord_pkl_files_all=sort_data_files_by_flushEnd(coord_pkl_path)

#each row of arr0 is a frame
arr0=read_one_pkl_file(sorted_coord_pkl_files_all[0])
name_vec=["O","H","H"]*total_molecule_num
name_vec=np.array(name_vec)
def row_2_block(row_num,one_row,name_vec):
    """
    :param row_num: row number
    :param one_row: a row of atoms, corresponding to one frame in xyz
    :param name_vec: name of atoms
    :return: a block of atom coordinates data, in xyz file
    """
    out_text=[
        f"{total_atom_num}\n",
        f"frame{row_num}\n",
    ]
    # one_row=list(one_row)
    str_one_row=np.array([f"{elem:.{n_dig}f}" for elem in one_row],dtype=str)
    str_one_row=str_one_row.reshape((-1,3))
    # print(str_one_row.shape)
    # print(name_vec.shape)
    combined=np.column_stack((name_vec,str_one_row))
    # print(combined)
    for row in combined:
        # Join the elements with 4 spaces and add a newline at the end
        formatted_line = "    ".join(row) + "\n"
        out_text.append(formatted_line)


    return out_text



# print(arr0[0][0])
out_text=row_2_block(0,arr0[0],name_vec)


out_xyz_file="./coord.xyz"
with open(out_xyz_file,"w+") as fptr:
    for file_ind in range(0,len(sorted_coord_pkl_files_all)):
        file_name_tmp=sorted_coord_pkl_files_all[file_ind]
        arr_tmp=read_one_pkl_file(file_name_tmp)
        row_num,_=arr_tmp.shape
        content=[]
        for j in range(0,row_num):
            out_text_tmp=row_2_block(j,arr_tmp[j],name_vec)

            fptr.writelines(out_text_tmp)
