import sys
import glob
import re
import json
from decimal import Decimal, getcontext
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
import pickle
import random

#this script loads previous data

numArgErr=4
valErr=5
if (len(sys.argv)!=3):
    print("wrong number of arguments.")
    exit(numArgErr)


jsonDataFromConf =json.loads(sys.argv[1])
jsonFromSummary=json.loads(sys.argv[2])

U_dist_dataDir=jsonFromSummary["U_dist_dataDir"]
startingFileInd=jsonFromSummary["startingFileInd"]

Nx=int(jsonDataFromConf["Nx"])
Ny=int(jsonDataFromConf["Ny"])
Nz=int(jsonDataFromConf["Nz"])

box_x_init=float(jsonDataFromConf["box_x_init"])
box_y_init=float(jsonDataFromConf["box_y_init"])
box_z_init=float(jsonDataFromConf["box_z_init"])
eps=1e-1
half_bond_angle=52/180*np.pi

#Angstrom
bond_length_approx=0.95
#search flushEnd
pklFileList=[]
flushEndAll=[]
#read U files
for file in glob.glob(U_dist_dataDir+"/U/flushEnd*.pkl"):
    pklFileList.append(file)
    matchEnd=re.search(r"flushEnd(\d+)",file)
    if matchEnd:
        flushEndAll.append(int(matchEnd.group(1)))


flushLastFile=-1

def gen_p2(nx,ny,nz,theta,p1,q1,r1,q2):
    p2Val=nz*np.cos(2*theta)/(p1*nz-r1*nx)-(q1*nz-r1*ny)/(p1*nz-r1*nx)*q2
    return p2Val
def gen_r2(nx,ny,nz,p2,q2):
    r2Val=-p2*nx/nz-q2*ny/nz

    return r2Val

def alpha_beta_to_directions(alpha,beta,theta):
    nx=np.sin(alpha)*np.cos(beta)

    ny=np.sin(alpha)*np.sin(beta)
    nz=np.cos(alpha)

    p1=np.random.uniform(-1+eps,1-eps)

    q1=np.sqrt(1-eps-p1**2)

    r1=-p1*nx/nz-q1*ny/nz

    while True:
        q2=np.random.uniform(-1+eps,1-eps)

        p2=gen_p2(nx,ny,nz,theta,p1,q1,r1,q2)

        if q2**2+p2**2<=1:
            break

    r2=gen_r2(nx,ny,nz,p2,q2)

    H1_direction=np.array([p1,q1,r1])

    H2_direction=np.array([p2,q2,r2])

    return H1_direction,H2_direction





# def create_init_atom_positions(U_dist_dataDir,Nx,Ny,Nz,box_x,box_y,box_z):
#     """
#     create initial atom positions, in the sequence O,H,H,O,H,H,...,O,H,H
#     :param U_dist_dataDir:
#     :param Nx:
#     :param Ny:
#     :param Nz:
#     :param box_x:
#     :param box_y:
#     :param box_z:
#     :return:
#     """
#     coord_init=[]
#     type_init=[]
#     sec_num_x=Nx+1
#
#     sec_num_y=Ny+1
#
#     sec_num_z=Nz+1
#
#     cell_length_x=box_x/sec_num_x
#
#     cell_length_y=box_y/sec_num_y
#
#     cell_length_z=box_z/sec_num_z
#     Ox_Oy_Oz_all=[]
#     Hx_Hy_Hz_all=[]
#     for a in range(1,sec_num_x):
#         for b in range(1,sec_num_y):
#             for c in range(1,sec_num_z):
#
#                 O_x=a*cell_length_x
#
#                 O_y=b*cell_length_y
#
#                 O_z=c*cell_length_z
#                 Ox_Oy_Oz_all.append(O_x)
#                 Ox_Oy_Oz_all.append(O_y)
#                 Ox_Oy_Oz_all.append(O_z)
#                 alpha=np.random.uniform(0,np.pi/2-eps)
#                 beta=np.random.uniform(0,2*np.pi)
#                 theta=np.random.uniform(half_bond_angle-eps/10,half_bond_angle+eps/10)
#                 d=np.random.uniform(bond_length_approx-eps/2,bond_length_approx+eps/2)
#                 H1_direction,H2_direction=alpha_beta_to_directions(alpha,beta,theta)
#
#                 r_O=np.array([O_x,O_y,O_z])
#
#                 H1_position=r_O+d*H1_direction
#                 H2_position=r_O+d*H2_direction
#
#                 H1_x=H1_position[0]
#
#                 H1_y=H1_position[1]
#
#                 H1_z=H1_position[2]
#
#                 H2_x=H2_position[0]
#
#                 H2_y=H2_position[1]
#
#                 H2_z=H2_position[2]
#                 Hx_Hy_Hz_all+=[H1_x,H1_y,H1_z,
#                                H2_x,H2_y,H2_z]
#                 coord_init+=[
#                             O_x,O_y,O_z,
#                              H1_x,H1_y,H1_z,
#                             H2_x,H2_y,H2_z
#                              ]
#                 type_init+=[0,1,1]
#
#     coord_init=np.array(coord_init)
#     type_init=np.array(type_init)
#     out_coord_dir=U_dist_dataDir+"/coord/"
#     Path(out_coord_dir).mkdir(exist_ok=True,parents=True)
#
#     out_coord_file=out_coord_dir+"/init.coord.pkl"
#     with open(out_coord_file,"wb") as fptr:
#         pickle.dump(coord_init,fptr)
#
#     out_type_file=U_dist_dataDir+"/raw.pkl"
#     with open(out_type_file,"wb") as fptr:
#         pickle.dump(type_init,fptr)
#
#     out_box_dir=U_dist_dataDir+"/box/"
#     Path(out_box_dir).mkdir(exist_ok=True,parents=True)
#     out_box_file=out_box_dir+"/init.box.pkl"
#     box=np.array([box_x,0,0,0,box_y,0,0,0,box_z])
#
#     with open(out_box_file,"wb") as fptr:
#         pickle.dump(box,fptr)
#     return Ox_Oy_Oz_all,Hx_Hy_Hz_all
#
#

def create_init_atom_positions(U_dist_dataDir, Nx, Ny, Nz, box_x, box_y, box_z):
    """
    Create initial atom positions, in the sequence O,H,H,O,H,H,...,O,H,H.

    The function ensures that both hydrogen coordinates (H)
    are within the simulation box defined by box_x, box_y, box_z.

    Parameters:
        U_dist_dataDir (str): Base directory for distribution data.
        Nx, Ny, Nz (int): Number of grid cells in x, y, and z directions.
        box_x, box_y, box_z (float): Dimensions of the simulation box.

    Returns:
        tuple: (Ox_Oy_Oz_all, Hx_Hy_Hz_all)
    """
    coord_init = []
    type_init = []
    sec_num_x = Nx + 1
    sec_num_y = Ny + 1
    sec_num_z = Nz + 1

    cell_length_x = box_x / sec_num_x
    cell_length_y = box_y / sec_num_y
    cell_length_z = box_z / sec_num_z

    Ox_Oy_Oz_all = []
    Hx_Hy_Hz_all = []

    for a in range(1, sec_num_x):
        for b in range(1, sec_num_y):
            for c in range(1, sec_num_z):
                # Compute oxygen coordinate
                O_x = a * cell_length_x
                O_y = b * cell_length_y
                O_z = c * cell_length_z
                Ox_Oy_Oz_all.extend([O_x, O_y, O_z])

                # Use a do-while style loop to generate valid hydrogen positions.
                while True:
                    # Generate random angles and d for hydrogen placement
                    alpha = np.random.uniform(0, np.pi / 2 - eps)
                    beta = np.random.uniform(0, 2 * np.pi)
                    theta = np.random.uniform(half_bond_angle - eps/10, half_bond_angle + eps/10)
                    d = np.random.uniform(bond_length_approx - eps/2, bond_length_approx + eps/2)

                    H1_direction, H2_direction = alpha_beta_to_directions(alpha, beta, theta)
                    r_O = np.array([O_x, O_y, O_z])

                    H1_position = r_O + d * H1_direction
                    H2_position = r_O + d * H2_direction

                    H1_x, H1_y, H1_z = H1_position
                    H2_x, H2_y, H2_z = H2_position

                    # Check if hydrogen positions are within the simulation box.
                    # Adjust the boundary conditions if your box origin is not (0,0,0)
                    if (0 <= H1_x <= box_x and 0 <= H1_y <= box_y and 0 <= H1_z <= box_z and
                            0 <= H2_x <= box_x and 0 <= H2_y <= box_y and 0 <= H2_z <= box_z):
                        break  # Accept these positions if the condition is satisfied

                # Append valid hydrogen positions to the lists.
                Hx_Hy_Hz_all.extend([H1_x, H1_y, H1_z, H2_x, H2_y, H2_z])
                coord_init.extend([O_x, O_y, O_z, H1_x, H1_y, H1_z, H2_x, H2_y, H2_z])
                type_init.extend([0, 1, 1])

    coord_init = np.array(coord_init)
    type_init = np.array(type_init)

    # Create output directories and save the data.
    out_coord_dir = U_dist_dataDir + "/coord/"
    Path(out_coord_dir).mkdir(exist_ok=True, parents=True)
    out_coord_file = out_coord_dir + "/init.coord.pkl"
    with open(out_coord_file, "wb") as fptr:
        pickle.dump(coord_init, fptr)

    out_type_file = U_dist_dataDir + "/raw.pkl"
    with open(out_type_file, "wb") as fptr:
        pickle.dump(type_init, fptr)

    out_box_dir = U_dist_dataDir + "/box/"
    Path(out_box_dir).mkdir(exist_ok=True, parents=True)
    out_box_file = out_box_dir + "/init.box.pkl"
    box = np.array([box_x, 0, 0, 0, box_y, 0, 0, 0, box_z])
    with open(out_box_file, "wb") as fptr:
        pickle.dump(box, fptr)

    return Ox_Oy_Oz_all, Hx_Hy_Hz_all
def create_loadedJsonData(flushLastFileVal):

    initDataDict={

        "flushLastFile":str(flushLastFileVal)
    }
    # print(initDataDict)
    return json.dumps(initDataDict)




#if no data found, return flush=-1
if len(pklFileList)==0:
    Ox_Oy_Oz_all,Hx_Hy_Hz_all=create_init_atom_positions(U_dist_dataDir,Nx,Ny,Nz,box_x_init, box_y_init, box_z_init)
    # with open("min_O.txt","w") as fptr:
    #     fptr.write(str(np.min(Ox_Oy_Oz_all)))
    # with open("min_H.txt","w") as fptr:
    #     fptr.write(str(min(Hx_Hy_Hz_all)))
    out_U_path=U_dist_dataDir+"/U/"
    Path(out_U_path).mkdir(exist_ok=True,parents=True)
    loadedJsonDataStr=create_loadedJsonData(flushLastFile)
    loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
    print(loadedJsonData_stdout)
    exit(0)


#if found pkl data with flushEndxxxx
sortedEndInds=np.argsort(flushEndAll)
sortedflushEnd=[flushEndAll[ind] for ind in sortedEndInds]
loadedJsonDataStr=create_loadedJsonData(sortedflushEnd[-1])
loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
print(loadedJsonData_stdout)
exit(0)

