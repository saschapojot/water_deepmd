import re
import sys

import json
import os

# this script parse conf file and return the parameters as json data

fmtErrStr = "format error: "
fmtCode = 1
valueMissingCode = 2
paramErrCode = 3
fileNotExistErrCode = 4

if (len(sys.argv) != 2):
    print("wrong number of arguments.")
    exit(paramErrCode)
inConfFile = sys.argv[1]

def removeCommentsAndEmptyLines(file):
    """

    :param file: conf file
    :return: contents in file, with empty lines and comments removed
    """
    with open(file, "r") as fptr:
        lines = fptr.readlines()

    linesToReturn = []
    for oneLine in lines:
        oneLine = re.sub(r'#.*$', '', oneLine).strip()
        if not oneLine:
            continue
        else:
            linesToReturn.append(oneLine)
    return linesToReturn


def parseConfContents(file):
    """

    :param file: conf file
    :return:
    """
    file_exists = os.path.exists(file)
    if not file_exists:
        print(file + " does not exist,")
        exit(fileNotExistErrCode)

    linesWithCommentsRemoved = removeCommentsAndEmptyLines(file)

    TStr = ""
    PStr=""
    NxStr=""
    NyStr=""
    NzStr=""
    box_xStr=""
    box_yStr=""
    box_zStr=""
    obs_name = ""
    effective_data_num_required = ""
    sweep_to_write = ""
    default_flush_num = ""
    hStr = ""
    swp_multiplyStr = ""
    for oneLine in linesWithCommentsRemoved:
        matchLine = re.match(r'(\w+)\s*=\s*(.+)', oneLine)
        if matchLine:
            key = matchLine.group(1).strip()
            value = matchLine.group(2).strip()

            # match T
            if key == "T":
                match_TValPattern = re.match(r"T\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$", oneLine)
                if match_TValPattern:
                    TStr = match_TValPattern.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)
            # match Nx
            if key=="Nx":
                match_Nx_pattern=re.match(r"Nx\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$", oneLine)
                if match_Nx_pattern:
                    NxStr=match_Nx_pattern.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)

            # match Ny
            if key=="Ny":
                match_Ny_pattern=re.match(r"Ny\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$", oneLine)
                if match_Ny_pattern:
                    NyStr=match_Ny_pattern.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)

            # match Nz
            if key=="Nz":
                match_Nz_pattern=re.match(r"Nz\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$", oneLine)
                if match_Nz_pattern:
                    NzStr=match_Nz_pattern.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)
            # match box_x
            if key=="box_x":
                match_box_x_pattern=re.match(r"box_x\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$", oneLine)
                if match_box_x_pattern:
                    box_xStr=match_box_x_pattern.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)
            # match box_y
            if key=="box_y":
                match_box_y_pattern=re.match(r"box_y\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$", oneLine)
                if match_box_y_pattern:
                    box_yStr=match_box_y_pattern.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)

            # match box_z
            if key=="box_z":
                match_box_z_pattern=re.match(r"box_z\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$", oneLine)
                if match_box_z_pattern:
                    box_zStr=match_box_z_pattern.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)


            # match observable_name
            if key == "observable_name":
                # if matching a non word character
                if re.search(r"[^\w]", value):
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)

                obs_name = value

            # match sweep_to_write
            if key == "sweep_to_write":
                if re.search(r"[^\d]", value):
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)
                sweep_to_write = value

            # match default_flush_num
            if key == "default_flush_num":
                if re.search(r"[^\d]", value):
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)
                default_flush_num = value

            # match effective_data_num_required
            if key == "effective_data_num_required":
                if re.search(r"[^\d]", value):
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)
                effective_data_num_required = value

            # match h
            if key == "h":
                match_h = re.match(r'([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)', value)
                # print(value)
                if match_h:
                    hStr = match_h.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)
            # match sweep_multiply
            if key == "sweep_multiple":
                match_swpMultiply = re.match(r"(\d+)", value)
                if match_swpMultiply:
                    swp_multiplyStr = match_swpMultiply.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)
        else:
            print("line: " + oneLine + " is discarded.")
            continue
    if TStr == "":
        print("T not found in " + str(file))
        exit(valueMissingCode)

    if NxStr=="":
        print("Nx not found in " + str(file))
        exit(valueMissingCode)
    if NyStr=="":
        print("Ny not found in " + str(file))
        exit(valueMissingCode)

    if NzStr=="":
        print("Nz not found in " + str(file))
        exit(valueMissingCode)
    if box_xStr=="":
        print("box_x not found in "+str(file))
        exit(valueMissingCode)
    if box_yStr=="":
        print("box_y not found in "+str(file))
        exit(valueMissingCode)
    if box_zStr=="":
        print("box_z not found in "+str(file))
        exit(valueMissingCode)


    if effective_data_num_required == "":
        print("effective_data_num_required not found in " + str(file))
        exit(valueMissingCode)

    if sweep_to_write == "":
        print("sweep_to_write not found in " + str(file))
        exit(valueMissingCode)

    if default_flush_num == "":
        print("default_flush_num not found in " + str(file))
        exit(valueMissingCode)
    if hStr == "":
        print("h not found in " + str(file))
        exit(valueMissingCode)
    if obs_name == "":
        print("observable_name not found in " + str(file))
        exit(valueMissingCode)
    if swp_multiplyStr == "":
        swp_multiplyStr = "1"

    dictTmp = {
        "T": TStr,
        "Nx":NxStr,
        "Ny":NyStr,
        "Nz":NzStr,
        "box_x":box_xStr,
        "box_y":box_yStr,
        "box_z":box_zStr,
        "observable_name": obs_name,
        "effective_data_num_required": effective_data_num_required,
        "sweep_to_write": sweep_to_write,
        "default_flush_num": default_flush_num,
        "confFileName": file,
        "h": hStr,
        "sweep_multiple": swp_multiplyStr,

    }
    return dictTmp


jsonDataFromConf=parseConfContents(inConfFile)

confJsonStr2stdout = "jsonDataFromConf=" + json.dumps(jsonDataFromConf)

print(confJsonStr2stdout)