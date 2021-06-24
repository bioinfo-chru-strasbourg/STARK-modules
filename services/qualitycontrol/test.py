import os
import fnmatch
import argparse
import re
import time
import glob

def parseArgs():
    args = myoptions()
    #args.nonregScript, args.nonregArgs
    val = dict()
    parse = args.nonregArgs
    for item in parse.split(" "):
       ar = item.rsplit("=",1)
       val[ar[0]] = ar[-1]
    return val


       #

       #id = item.split("=")[0]
       #i = re.sub(r'[^A-Za-z0-9 ]+', '', id)
       #value = item.split("=")[1]
       #print(value)
       #val[i] = item.split("=")[1]
       
def checkResults():
    rg = parseArgs()
    keys = ['-r', '--results']
    output = ""
    for val in keys:
        if val not in rg:
            continue
        else:
            output = rg[val]
    if output == "":
        print("#ERROR: it appears that the output folder option was not specified and results could be located in the same folder as the script path")
        print("#[INFO] End of nonregression")
        exit()
    else:
        print("#[INFO] Results folder: "+output)
 
        #output variable is results path
        for file in os.listdir(output):
            folder = output+"/"+file
            if re.search("VALIDATION", file):
                #print("FILE="+file)
                #print(folder)
                for files in glob.iglob(folder+"/validation*.txt"):
                    #print(files)
                    with open(files) as f:
                        if 'correctly generated' in f.read():
                            print("#ERROR issues in STARK analysis, VCF not correclty generated "+ files +"!")
                        else:
                            print("#[INFO] Validation folder: "+os.path.dirname(files))
        print("#[INFO] End of nonregression")
        exit()
    
def myoptions():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--script", type = str, default = " ", help = "Path of the non regression script: <PATH_NON-REG_SCRIPT>", dest = 'nonregScript')
    parser.add_argument("-a", "--args", type = str, default = " ", help = "Args of the non regression script", dest = 'nonregArgs')
    return parser.parse_args()

if __name__ == "__main__":
    #input = " --set=TEST --folder=/app/databases --cov=30 --results=/home1/data/STARK/data/TESTJB"
    #checkResults(input)
    args = myoptions()
    parseArgs()
    checkResults()
    #./nonregression/cli/dockerfile/app/bin/docker_VALIDATION_SETS_JB.sh --set=TEST --folder=/home1/data/STARK/databases/validation --cov=30 --results=/home1/data/STARK/data/TESTJB

#python3 test.py -s='/home1/bin/STARK-modules/dev/services/qualitycontrol/nonregression/cli/dockerfile/app/bin/docker_VALIDATION_SETS_JB.sh' -a=' --set=TEST --folder=/home1/data/STARK/databases/validation --cov=30 --results=/home1/data/STARK/data/TESTJB'

#python3 nonreg.py -s='/home1/bin/STARK-modules/dev/services/qualitycontrol/nonregression/cli/dockerfile/app/bin/docker_VALIDATION_SETS_JB.sh' -a=' --set=TEST --folder=/home1/data/STARK/databases/validation --cov=30 --results=/home1/data/STARK/data/TESTJB'