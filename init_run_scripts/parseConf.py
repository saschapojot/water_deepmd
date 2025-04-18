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
    model_fileStr=""
    NxStr=""
    NyStr=""
    NzStr=""
    box_x_initStr=""
    box_y_initStr=""
    box_z_initStr=""
    obs_nameStr = ""
    effective_data_num_requiredStr = ""
    sweep_to_writeStr = ""
    default_flush_numStr = ""
    hStr = ""
    swp_multiplyStr = ""
    box_upper_boundStr=""
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
            # match P
            if key=="P":
                match_PValPattern=re.match(r"P\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$", oneLine)
                if match_PValPattern:
                    PStr=match_PValPattern.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)
            #match model_file
            if key=="model_file":
                match_model_file_Pattern=re.match(r"model_file\s*=\s*(.+)",oneLine)
                if match_model_file_Pattern:
                    model_fileStr=match_model_file_Pattern.group(1)

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
            if key=="box_x_init":
                match_box_x_pattern=re.match(r"box_x_init\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$", oneLine)
                if match_box_x_pattern:
                    box_x_initStr=match_box_x_pattern.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)
            # match box_y
            if key=="box_y_init":
                match_box_y_pattern=re.match(r"box_y_init\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$", oneLine)
                if match_box_y_pattern:
                    box_y_initStr=match_box_y_pattern.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)

            # match box_z
            if key=="box_z_init":
                match_box_z_pattern=re.match(r"box_z_init\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$", oneLine)
                if match_box_z_pattern:
                    box_z_initStr=match_box_z_pattern.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)


            # match observable_name
            if key == "observable_name":
                # if matching a non word character
                if re.search(r"[^\w]", value):
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)

                obs_nameStr = value

            # match sweep_to_write
            if key == "sweep_to_write":
                if re.search(r"[^\d]", value):
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)
                sweep_to_writeStr = value

            # match default_flush_num
            if key == "default_flush_num":
                if re.search(r"[^\d]", value):
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)
                default_flush_numStr = value

            # match effective_data_num_required
            if key == "effective_data_num_required":
                if re.search(r"[^\d]", value):
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)
                effective_data_num_requiredStr = value

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
            # box_upper_bound
            if key=="box_upper_bound":
                match_box_upper_bound=re.match(r"([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",value)
                if match_box_upper_bound:
                    box_upper_boundStr=match_box_upper_bound.group(1)
                else:
                    print(fmtErrStr + oneLine)
                    exit(fmtCode)

        else:
            print("line: " + oneLine + " is discarded.")
            continue
    if TStr == "":
        print("T not found in " + str(file))
        exit(valueMissingCode)
    if PStr=="":
        print("P not found in " + str(file))
        exit(valueMissingCode)
    if model_fileStr=="":
        print("model_file not found in " + str(file))
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
    if box_x_initStr=="":
        print("box_x not found in "+str(file))
        exit(valueMissingCode)
    if box_y_initStr=="":
        print("box_y not found in "+str(file))
        exit(valueMissingCode)
    if box_z_initStr=="":
        print("box_z not found in "+str(file))
        exit(valueMissingCode)


    if effective_data_num_requiredStr == "":
        print("effective_data_num_required not found in " + str(file))
        exit(valueMissingCode)

    if sweep_to_writeStr == "":
        print("sweep_to_write not found in " + str(file))
        exit(valueMissingCode)

    if default_flush_numStr == "":
        print("default_flush_num not found in " + str(file))
        exit(valueMissingCode)
    if hStr == "":
        print("h not found in " + str(file))
        exit(valueMissingCode)
    if obs_nameStr == "":
        print("observable_name not found in " + str(file))
        exit(valueMissingCode)
    if swp_multiplyStr == "":
        swp_multiplyStr = "1"
    if box_upper_boundStr=="":
        print("box_upper_bound not found in " + str(file))
        exit(valueMissingCode)

    dictTmp = {
        "T": TStr,
        "P":PStr,
        "model_file":model_fileStr,
        "Nx":NxStr,
        "Ny":NyStr,
        "Nz":NzStr,
        "box_x_init":box_x_initStr,
        "box_y_init":box_y_initStr,
        "box_z_init":box_z_initStr,
        "observable_name": obs_nameStr,
        "effective_data_num_required": effective_data_num_requiredStr,
        "sweep_to_write": sweep_to_writeStr,
        "default_flush_num": default_flush_numStr,
        "confFileName": file,
        "h": hStr,
        "sweep_multiple": swp_multiplyStr,
        "box_upper_bound":box_upper_boundStr

    }
    return dictTmp


jsonDataFromConf=parseConfContents(inConfFile)

confJsonStr2stdout = "jsonDataFromConf=" + json.dumps(jsonDataFromConf)

print(confJsonStr2stdout)