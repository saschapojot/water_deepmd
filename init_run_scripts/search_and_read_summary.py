import re
from decimal import Decimal
import json
import sys
from pathlib import Path
import os
import shutil


import numpy as np

#this script initializes the parameters for mc computation by reading summary file
#and also creates/erases necessary folders

invalidValueErrCode=1
mcErrCode=2
pathErrCode=3
numArgErr=4

if (len(sys.argv)!=2):
    print("wrong number of arguments.")
    exit(numArgErr)


jsonDataFromConf=json.loads(sys.argv[1])
#read json
T=float(jsonDataFromConf["T"])
if T<=0:
    print("invalid temperature: "+str(T))
    exit(invalidValueErrCode)


effective_data_num_required=int(jsonDataFromConf["effective_data_num_required"])
sweep_to_write=int(jsonDataFromConf["sweep_to_write"])
default_flush_num=int(jsonDataFromConf["default_flush_num"])


sweep_multiple=int(jsonDataFromConf["sweep_multiple"])

confFileName=jsonDataFromConf["confFileName"]
TDirRoot=os.path.dirname(confFileName)
U_dist_dataDir=TDirRoot+"/U_dist_dataFiles/"
#create dataDir if not exists
Path(U_dist_dataDir).mkdir(exist_ok=True, parents=True)

swpNumInOneFlush=sweep_to_write*sweep_multiple


#parameters to guide mc computation
lag=-1
startingFileInd=-1
newDataPointNum=-1
newMcStepNum=swpNumInOneFlush*default_flush_num
newFlushNum=default_flush_num

def create_jsonFromSummary(startingFileIndVal,newMcStepNumVal,
                           newDataPointNumVal,newFlushNumVal,TDirRootStr,U_dist_dataDirStr):

    """

    :param startingFileIndVal:
    :param newMcStepNumVal:
    :param newDataPointNumVal:
    :param newFlushNumVal:
    :param TDirRootStr:
    :param U_dist_dataDirStr:
    :return: jsonFromSummary as string
    """
    outDict={
        "startingFileInd":str(startingFileIndVal),
        "newMcStepNum":str(newMcStepNumVal),
        "newDataPointNum":str(newDataPointNumVal),
        "newFlushNum": str(newFlushNumVal),
        "TDirRoot":str(TDirRootStr),
        "U_dist_dataDir": str(U_dist_dataDirStr),
    }

    return json.dumps(outDict)

obs_name=jsonDataFromConf["observable_name"]
summaryFileName=TDirRoot+"/summary_"+obs_name+".txt"
summaryFileExists= os.path.isfile(summaryFileName)

#if summary file does not exist, return -1, sweep_to_write*default_flush_num, then exit with code 0
if summaryFileExists==False:
    jsonFromSummaryStr=create_jsonFromSummary(startingFileInd,newMcStepNum,
                                              newDataPointNum,newFlushNum,TDirRoot,U_dist_dataDir)
    jsonFromSummary_stdout="jsonFromSummary="+jsonFromSummaryStr
    print(jsonFromSummary_stdout)
    exit(0)


#parse summary file
with open(summaryFileName,"r") as fptr:
    linesInSummaryFile= fptr.readlines()


for oneLine in linesInSummaryFile:
    matchErr=re.search(r"error",oneLine)
    #if "error" is matched
    if matchErr:
        print("error in previous computation, please re-run.")
        exit(mcErrCode)

    #if "continue" is matched
    matchContinue=re.search(r"continue",oneLine)
    if matchContinue:
        jsonFromSummaryStr=create_jsonFromSummary(startingFileInd,newMcStepNum,
                                                  newDataPointNum,newFlushNum,TDirRoot,U_dist_dataDir)
        jsonFromSummary_stdout="jsonFromSummary="+jsonFromSummaryStr
        print(jsonFromSummary_stdout)
        exit(0)
    #if "high" is matched
    matchHigh=re.search(r"high",oneLine)
    if matchHigh:
        jsonFromSummaryStr=create_jsonFromSummary(startingFileInd,newMcStepNum,
                                                  newDataPointNum,newFlushNum,TDirRoot,U_dist_dataDir)
        jsonFromSummary_stdout="jsonFromSummary="+jsonFromSummaryStr
        print(jsonFromSummary_stdout)
        exit(0)

    #the rest of the cases is "equilibrium"
    matchEq=re.search(r"equilibrium",oneLine)
    if matchEq:
        continue

    #match lag
    matchLag=re.match(r"lag\s*=\s*(\d+)",oneLine)
    if matchLag:
        lag=int(matchLag.group(1))

    #match newDataPointNum
    matchNew=re.match(r"newDataPointNum\s*=\s*(\d+)",oneLine)
    if matchNew:
        newDataPointNum=int(matchNew.group(1))
    #match startingFileInd
    matchStartingFileInd=re.match(r"startingFileInd\s*=\s*(\d+)",oneLine)
    if matchStartingFileInd:
        startingFileInd=int(matchStartingFileInd.group(1))


newMcStepNum=lag*newDataPointNum

newFlushNum=int(np.ceil(newMcStepNum/(sweep_to_write)))

jsonFromSummaryStr=create_jsonFromSummary(startingFileInd,newMcStepNum,
                                          newDataPointNum,newFlushNum,TDirRoot,U_dist_dataDir)
jsonFromSummary_stdout="jsonFromSummary="+jsonFromSummaryStr
print(jsonFromSummary_stdout)
exit(0)