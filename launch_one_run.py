import re
import subprocess
import sys

import json
argErrCode=2
if (len(sys.argv)!=2):
    print("wrong number of arguments")
    print("example: python launch_one_run.py /path/to/mc.conf")
    exit(argErrCode)


confFileName=str(sys.argv[1])
invalidValueErrCode=1
summaryErrCode=2
loadErrCode=3
confErrCode=4

#################################################
#parse conf, get jsonDataFromConf
confResult=subprocess.run(["python3", "./init_run_scripts/parseConf.py", confFileName], capture_output=True, text=True)
confJsonStr2stdout=confResult.stdout
# print(confJsonStr2stdout)
if confResult.returncode !=0:
    print("Error running parseConf.py with code "+str(confResult.returncode))
    # print(confResult.stderr)
    exit(confErrCode)
match_confJson=re.match(r"jsonDataFromConf=(.+)$",confJsonStr2stdout)
if match_confJson:
    jsonDataFromConf=json.loads(match_confJson.group(1))
else:
    print("jsonDataFromConf missing.")
    exit(confErrCode)
# print(jsonDataFromConf)

##################################################

##################################################

#read summary file, get jsonFromSummary
parseSummaryResult=subprocess.run(["python3","./init_run_scripts/search_and_read_summary.py", json.dumps(jsonDataFromConf)],capture_output=True, text=True)
# print(parseSummaryResult.stdout)
if parseSummaryResult.returncode!=0:
    print("Error in parsing summary with code "+str(parseSummaryResult.returncode))
    # print(parseSummaryResult.stdout)
    # print(parseSummaryResult.stderr)
    exit(summaryErrCode)

match_summaryJson=re.match(r"jsonFromSummary=(.+)$",parseSummaryResult.stdout)
if match_summaryJson:
    jsonFromSummary=json.loads(match_summaryJson.group(1))
# print(jsonFromSummary)

##################################################


###############################################
#load previous data, to get paths
#get loadedJsonData
# print("before entering load")
loadResult=subprocess.run(["python3","./init_run_scripts/load_previous_data.py", json.dumps(jsonDataFromConf), json.dumps(jsonFromSummary)],capture_output=True, text=True)

if loadResult.returncode!=0:
    print("Error in loading with code "+str(loadResult.returncode))
    print(loadResult.stdout)
    print(loadResult.stderr)
    exit(loadErrCode)

match_loadJson=re.match(r"loadedJsonData=(.+)$",loadResult.stdout)
if match_loadJson:
    loadedJsonData=json.loads(match_loadJson.group(1))
else:
    print("loadedJsonData missing.")
    exit(loadErrCode)

# print(f"loadedJsonData={loadedJsonData}")
###############################################

###############################################
#construct parameters that are passed to mc

TStr=jsonDataFromConf["T"]
PStr=jsonDataFromConf["P"]
model_fileStr=jsonDataFromConf["model_file"]
sweep_to_write=jsonDataFromConf["sweep_to_write"]
flushLastFile=loadedJsonData["flushLastFile"]

newFlushNum=jsonFromSummary["newFlushNum"]
TDirRoot=jsonFromSummary["TDirRoot"]
U_dist_dataDir=jsonFromSummary["U_dist_dataDir"]

hStr=jsonDataFromConf["h"]
sweep_multipleStr=jsonDataFromConf["sweep_multiple"]

NxStr=jsonDataFromConf["Nx"]

NyStr=jsonDataFromConf["Ny"]

NzStr=jsonDataFromConf["Nz"]
box_x_init=jsonDataFromConf["box_x_init"]

box_y_init=jsonDataFromConf["box_y_init"]
box_z_init=jsonDataFromConf["box_z_init"]
box_upper_bound=jsonDataFromConf["box_upper_bound"]
params2cppInFile=[
    TStr+"\n",
    PStr+"\n",
    model_fileStr+"\n",
    sweep_to_write+"\n",
    flushLastFile+"\n",
    newFlushNum+"\n",
    TDirRoot+"\n",
    U_dist_dataDir+"\n",
    hStr+"\n",
    sweep_multipleStr+"\n",
    NxStr+"\n",
    NyStr+"\n",
    NzStr+"\n",
    # box_x_init+"\n",
    # box_y_init+"\n",
    # box_z_init+"\n",
    box_upper_bound+"\n"
    ]

cppInParamsFileName=TDirRoot+"/cppIn.txt"
with open(cppInParamsFileName,"w+") as fptr:
    fptr.writelines(params2cppInFile)