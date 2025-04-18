from pathlib import Path
from decimal import Decimal, getcontext
#This script creates directories and conf files for mc


def format_using_decimal(value, precision=15):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

#number of Oxygen atoms in each direction
Nx=4
Ny=4
Nz=4
default_flush_num=1
#T vals, unit is K
#P vals, unit is GPa

TVals=[100,200]

PVals=[5e-2,7e-2,9e-2,1e-1,1.2e-1,1.3e-1,1.4e-1]
#P directories are within 1 T directory
TStrAll=[]
PStrAll=[]
for j in range(0,len(TVals)):
    T=TVals[j]

    TStr=format_using_decimal(T)

    TStrAll.append(TStr)



for k in range(0,len(PVals)):
    P=PVals[k]
    PStr=format_using_decimal(P)
    PStrAll.append(PStr)

model_file="./se_e2_a/compressed_model_water.pth"
effective_data_num_required=200
sweep_to_write=5
sweep_multiple=7
observable_name="U_dist"
dataRoot="./mcDataAll/"
# Path(dataRoot).mkdir(exist_ok=True,parents=True)
dataOutDir=dataRoot
box_x_init=12.5
box_y_init=12.5
box_z_init=12.5
box_upper_bound=14
h=box_x_init/(Nx+1)/5
def contents_to_conf(j,k):
    """

    :param j: index of T
    :param k: index of P
    :return:
    """
    contents=[
        "#This is the configuration file for mc computations\n",
        "#System has H2O\n",
        "\n" ,
        "#Temperature, K\n",
        f"T={TStrAll[j]}\n",
        "\n",
        "#pressure, GPa\n",
        f"P={PStrAll[k]}\n",
        "\n",
        f"model_file={model_file}\n",
        "\n",
        f"Nx={Nx}\n",
        "\n",
        f"Ny={Ny}\n",
        "\n",
        f"Nz={Nz}\n",
        "\n",
        "#box size, A\n"
        f"box_x_init={box_x_init}\n",
        "\n",
        f"box_y_init={box_y_init}\n",
        "\n",
        f"box_z_init={box_z_init}\n",
        "\n",
        "#this is the data number in each pkl file, i.e., in each flush\n"
        f"sweep_to_write={sweep_to_write}\n",
        "\n",
        "#within each flush,  sweep_to_write*sweep_multiple mc computations are executed\n",
        "\n",
        f"default_flush_num={default_flush_num}\n",
        "\n",
        f"observable_name={observable_name}\n",
        "\n",
        "#coordinate step length, A\n"
        f"h={h}\n",
        "#the configurations of the system are saved to file if the sweep number is a multiple of sweep_multiple\n",
        "\n",
        f"sweep_multiple={sweep_multiple}\n",
        "\n",
        f"effective_data_num_required={effective_data_num_required}\n",
        f"box_upper_bound={box_upper_bound}\n"
        ]
    outDir=dataOutDir+f"/Nx{Nx}_Ny{Ny}_Nz{Nz}/T{TStrAll[j]}/P{PStrAll[k]}/"
    Path(outDir).mkdir(exist_ok=True,parents=True)
    outConfName=outDir+f"/run_T{TStrAll[j]}_P{PStrAll[k]}"+".mc.conf"

    with open(outConfName,"w+") as fptr:
        fptr.writelines(contents)


for j in range(0,len(TStrAll)):
    for k in range(0,len(PStrAll)):
        contents_to_conf(j,k)