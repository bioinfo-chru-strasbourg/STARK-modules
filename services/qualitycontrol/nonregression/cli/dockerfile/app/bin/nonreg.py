
import sys 
import os
import re
import subprocess
import glob
import time
import argparse
import fnmatch
import logging
import time
import datetime as dt

#from alive_progress import alive_bar

#############
#          python3 nonreg.py -s "docker_VALIDATION_SETS_JB.sh" -a "--set=CPSGEN --folder=/home1/data/STARK/databases/validation --results=/home1/data/STARK/data/TESTJB"
###############

'''
This script is set to run non-regression for all version of STARK,
it will be needed to assure that STARK is still able to retrieve at least the same number of variants than the previous version
Author: Jean-Baptiste LAMOUCHE

Prerequisite:
- STARK will be launch in the STARK-CLI-modify, need to be up and image build, use this image and 
- assure bcftools, tabix, bgzip, samtools, vcftools are available in the CLI
- assure you have enough memory especially for JAVA based tools
- SETS folder (containing Raw data for each datasets) are available in the container, rights are OK
- output folder is the good one 


'''

#"/STARK/services/qualitycontrol/nonregression/cli"
log = 'test.ronregression.log'

def launch(scriptname="nonregression", release="0.1", date="20201219", author="Jean-Baptiste Lamouche", group="UF7363 Bioinfo"):
    #Header of the script
    
    print("#############################################")
    print("## Name               "+scriptname)
    print("## Release            "+release+"/"+date)
    print("## Author             "+author)
    print("## Copyright          "+group)
    print("#############################################")
    print("\n")
    print("[#INFO] LAUNCH NONREGRESSION")

def setupAnalysis(cli, tools):
    
    #Check Folder and input files and container existence 
    #run(["docker",  ])
    print("[#INFO] Check configuration")
    stock = []
    #tools = ['bcftools', 'samtools', 'htslib', 'vcftools']
    T = subprocess.Popen(["docker ps -f name=stark --no-trunc --format 'table {{.Image}}\t{{.Names}}\t{{.Status}}' | grep -w "+cli], stdout=subprocess.PIPE, shell=True, executable="/bin/bash")

    docker = T.communicate()
    if docker:
        #print(T.stdout.readlines())
        for items in list(docker):
            #print("ITEMS="+items.strip())
            if items.decode('utf-8') != 0:
                print("[#INFO] Containers infos: "+items.decode('utf-8')+"\n")
                condi = subprocess.Popen(["ls -1 /STARK/tools 2>>"+log], stdout=subprocess.PIPE, shell=True, executable="/bin/bash"  )
                for elems in condi.stdout.readlines():
                    stock.append(elems.decode('utf-8').strip())
                stock.sort()
                tools.sort()
                for i in tools:
                    if i in stock:
                        continue
                    else:
                        print("ERROR "+i+" is not available in the container, you shoud rebuild it")
                        exit()
                print("#[INFO] Tools available: ", tools)
                checkRam()

                #time.sleep(10)
                #global search
                search = subprocess.Popen(['find /app/ -type d -name SETS_STARK -exec bash -c "ls -1 {} | wc -l \; 2>>'+log], stdout=subprocess.PIPE, shell=True,  executable='/bin/bash')
                for u in search.stdout.readlines():
                    #Number of datasets in SETS folder
                    #print(u.decode('utf-8').strip())               
                    if int(u.decode('utf-8').strip()) > 0:
                        sets = subprocess.Popen(['find /app/ -type d -name SETS -exec bash -c "ls -C {}" \; 2>>'+log], stdout=subprocess.PIPE, shell=True,  executable='/bin/bash') 
                        for items in sets.stdout.readlines():
                            print("\n[#INFO] Datasets available "+items.decode('utf-8').strip())
                        print("[#INFO] Setup | tools | config OK")
                    

            else:
                print("ERROR SETS are not available in the container, you should rebuild it or move SETS folder in the CLI")
                exit()
                
                #if u is None:
                #    continue
                #else:
                #    print(u.decode('utf-8').strip())
    #else:
    #    print("#[INFO] STARK cli not UP")
    #    print("#[INFO] Looking for STARK services")
    #    runCommand(cmd="/home1/bin/STARK-modules/dev/services/services.sh", args=" --modules=qualitycontrol --submodules=nonregression --command=up")
            
    
def runCommand(cmd, args):
    #Allowed to run bash function in python
    #print(cmd+" "+args+" 2>>"+log)
    subprocess.call(cmd+" "+args+" 2>>"+log, shell=True)
    #V = subprocess.Popen([cmd, args], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    #V.stdout.readlines()

def checkRam():
    val= subprocess.Popen(["free -m | awk '/^Mem/ {print$7}' 2>>"+log], stdout=subprocess.PIPE, shell=True, executable='/bin/bash')
    for items in val.stdout.readlines():
        if int(items.decode('utf-8').strip()) < 10000:
            print("\nWARNING: Not recommanded to run STARK if you have less than 10Go of RAM available")
            runCommand("free", " -t -m -h")
        else:
            print("[#INFO] Ram available: OK")

def parseArgs():
    '''
    #####DEPRECATED
    It return a dict with bash nonregression script, 
    sets: CPSGEN, HUSHEMATO etc, list of sets, if this value is empty take all set in setfolder
    sets_path: path of the set folder
    results: path of STARK results 
    '''

    args = myoptions()
    #args.nonregScript, args.nonregArgs
    val = dict()
    parse = args.nonregArgs
    for item in parse.split(" "):
        ar = item.rsplit("=",1)
        if ar[0] == '-s' or ar[0] == '--set':
            val['sets'] = ar[1].split(",")
        
        elif ar[0] == '-f' or ar[0] == '--folder':     
            val['sets_path'] = ar[1]

        elif ar[0] == '-r' or ar[0]  == '--results':      
            val['results'] = ar[1]
        #val[ar[0]] = ar[-1]
    #print(val)
    return val


def checkResults(output):
    #TODO
    if output == "":
        print("ERROR: it appears that the output folder option was not specified and results could be located in the same folder as the script path")
        print("[#INFO] End of nonregression")
        exit()
    elif os.path.exists(output):
        print("[#INFO] Results folder: "+output)
        # check if files finish is located in results folder
        
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
                            print("ERROR Validation folder: "+os.path.dirname(files)+" issues in STARK analysis, VCF not correclty generated ")
                        else:
                            print("[#INFO] Validation folder: "+os.path.dirname(files))
        print("[#INFO] End of nonregression")
        exit()

    else:
        print("ERROR: results folder does not exist")
        print("[#INFO] End of nonregression")
        exit()

def compressAndIndex(input, output):
    '''
    Compress with bgzip and index vcf file with tabix input and output should contain full path
    '''
    print("[#INFO] Compress and Index "+input)
    logging.info("RUN bgzip tabix input: "+input+" output: "+output)
    runCommand("/home1/TOOLS/tools/bgzip/current/bin/bgzip", "-c "+input+" > "+output)
    runCommand("/home1/TOOLS/tools/tabix/current/bin/tabix", "-p vcf "+output)


def checkVcf(inputFolder, outputFolder):
    '''
    Check in folder for all vcf file if they are compress in an output folder, otherwise file will be process
    '''
    listVcf = []
    for vcf in os.listdir(inputFolder):
        if vcf.endswith(".vcf"):
            listVcf.append(vcf)

    for v in listVcf:
        if os.path.exists(outputFolder+"/"+v+'.gz') and os.path.getsize(outputFolder+"/"+v+'.gz') != 0:
            continue
        else:
            try:
                os.remove(outputFolder+"/"+v+'.gz')
            except:
                print("[#INFO] File "+outputFolder+"/"+v+'.gz does not exist, Creation of file')            
            compressAndIndex(inputFolder+"/"+v, outputFolder+"/"+v+".gz") 


def checkTagfile(output):
    ''' 
    Create or modify $sample.tag file if it already exists, this file specify from which reference it belong
    '''
    results_path = output
    #day = "VALIDATION_"+dt.date.today().strftime("%Y")+dt.date.today().strftime("%m")+dt.date.today().strftime("%d")
    #print(day)
    for res in results_path:
        if "CPSGEN" in res:
            print("RUN="+res)
            for sample in os.listdir(os.path.join(results_path, res)):
                #print(sample)
                if os.path.isdir(os.path.join(results_path, res,sample)):
                    print("SAMPLE="+sample)
                    #print("echo CQI#"+sample+" > "+os.path.join(results_path, res, sample, sample+".tag"))
                    #runCommand("echo","OK")
                    with open(os.path.join(results_path, res, sample, sample+'.tag')) as f:
                        if not 'CQI' in f.read():
                            ag = 'CQI#'+sample+' > '+os.path.join(results_path, res, sample, sample+'.tag')
                            #print(ag)
                            runCommand('echo',ag)
                            #print("echo "+ag)
                        else:
                            print("[#INFO] TAG OK "+os.path.join(results_path, res, sample, sample+'.tag'))


        

def createDictRef(reference, datasets):
    #print(paramDict)
    #sets = []
    #set_folder = []
    #for key, items in paramDict.items():
    #    if key == '-s' or key == '--set':
    #        sets.append(items)
    #    elif key == '-f' or key == '--folder':
    #        set_folder.append(items)
    sets_path = reference+'/SETS_STARK/'
    #
    #if datasets:
    #    true_sets = datasets.split(',')
    #else:
    #    true_sets = []
    #    [true_sets.append(j) for j in os.listdir(sets_path)]
    #
    #[print("[#INFO] Dataset: "+i) for i in true_sets]
    
    #Create CQI dict to calculate metrics:
    for datas in datasets:
        fdata = os.path.join(sets_path, datas)
        #print(type(fdata))
        file = fdata+"/"+datas+"_CQI_dict.json"
        print("[#INFO] JSON: "+file)       
        if not os.path.exists(fdata+"/tmp.CQI"): 
            os.mkdir(fdata+"/tmp.CQI")
        

        #Write to JSON 
        cqilist = []
        if os.path.exists(file):
            os.remove(file)
        with open(file, "a+") as f:
            checkVcf(fdata, fdata+"/tmp.CQI")
            for vcf in os.listdir(fdata+"/tmp.CQI/"):
                if vcf.endswith(".vcf.gz"):
                    dict_json = {}
                    print("[#INFO] convert "+datas+" "+vcf+" into json")
                    dict_json["name"] = os.path.basename(os.path.splitext(os.path.splitext(vcf)[0])[0])
                    dict_json["VCF"] = fdata+"/tmp.CQI/"+vcf
                    cqilist.append(dict_json)

            f.write("{ \
                        CQI"+str(cqilist)+"\
                            }")
            f.close()
    

   
def main_all(cli, tools, reference, datasets, output, script):
    now = dt.datetime.now()
    ##print("[#INFO] "+now.strftime("%Y-%m-%d %H:%M:%S"))
    #launch()
    #logging.basicConfig(filename=log,format='%(levelname)s:%(asctime)s %(message)s',\
    #     datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)
    #logging.info(" START Nonregression")
    #
#
    #setupAnalysis(cli, tools)
    ##LAUNCH STARK
    #
    #print("[#INFO] "+now.strftime("%Y-%m-%d %H:%M:%S")+" STARK analysis")
    #runCommand(script, "-s "+datasets+" -f "+reference+" -r "+output)
    #folder = glob.glob(output+'/results/VALIDATION*') 
    #latest = max(folder, key=os.path.getctime)
    #
    #if datasets == 'ALL':
    #    datasets = os.listdir(os.path.join(reference, "SETS_STARK"))
    print("[#INFO] SETS used: ", datasets)
    for group in datasets:
        print("[#INFO] SET: "+group+" check for results folder")
        last = []
        for res in os.listdir(os.path.join(output,'results')):
            #print(os.listdir(os.path.join(output,'results')))
            if group in res:
                last.append(os.path.join(output,'results',res))
                #folder = glob.glob(output'/results/VALIDATION*') 
        if len(last) == 1:
            latest = last[0]
        else:
            latest = max(last, key=os.path.getctime)
        
        res = os.path.join(output,'results',latest)
        print("[#INFO] Results: "+res)
        while not os.path.exists(latest+'/STARKComplete.txt'):
            time.sleep(1)
            if os.path.isfile(res):
                createDictRef(reference, group)
                checkTagfile(res)
                print("[#INFO] "+now.strftime("%Y-%m-%d %H:%M:%S")+" CQI analysis "+res)
                #Run CQI
                runCommand('docker run', "--rm --name=nonregression_"+latest+" stark/stark-module-qualitycontrol-submodule-cqi-service-cli:1.0 /bin/bash -c 'source activate variant && /app/bin/CQI_1.3.sh --genome=/STARK/databases/genomes/current/hg19.fa  --json="+reference+"/SETS_STARK/"+group+"/"+group+"_CQI_dict.json"+" --run="+res+"'")
            elif os.path.isfile(res) and os.path.isfile(res+'/CQIComplete.txt'):
                print("[#INFO] "+now.strftime("%Y-%m-%d %H:%M:%S")+" CQI already done")

    print("[#INFO] "+now.strftime("%Y-%m-%d %H:%M:%S")+" Exit nonregression")


        
        #Search if results are ok otherwise disp errors
        
        #stark = Path()
        #checkResults()

def main_nr(reference, datasets, output, script):
        print("[#INFO] Run: "+output)
        while not os.path.exists(output+'/STARKComplete.txt'):
            time.sleep(1)
            if os.path.isfile(output):
                createDictRef(reference, datasets)
                checkTagfile(output)
                print("[#INFO] "+now.strftime("%Y-%m-%d %H:%M:%S")+" CQI analysis "+output)
                #Run CQI
                runCommand('docker run', "--rm --name=nonregression_"+os.path.basename+" stark/stark-module-qualitycontrol-submodule-cqi-service-cli:1.0 /bin/bash -c 'source activate variant && /app/bin/CQI_1.3.sh --genome=/STARK/databases/genomes/current/hg19.fa  --json="+reference+"/SETS_STARK/"+datasets+"/"+datesets+"_CQI_dict.json"+" --run="+res+"'")
            elif os.path.isfile(res) and os.path.isfile(res+'/CQIComplete.txt'):
                print("[#INFO] "+now.strftime("%Y-%m-%d %H:%M:%S")+" CQI already done")

    print("[#INFO] "+now.strftime("%Y-%m-%d %H:%M:%S")+" Exit nonregression")
        



def outputCopy():
    '''
    Get application from file in row data, 
    VALIDATION_$date_SET_$setname_APP_$appname
    Defaut results are in run 
    '''
    args = myoptions()
    if args.output:
        outputFolder = []
        parse = parseArgs()
        for key, ars in parse.items():
            if key == '-r' or key =='--results':
                outputFolder.append(ars)



def myoptions():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    #ALL analysis
    parser_all = subparsers.add_parser('all')
    parser_all.add_argument("-s", "--script", type = str, default = " ", help = "PATH of the non regression script: <PATH_NON-REG_SCRIPT>", dest = 'script')
    parser_all.add_argument("-d", "--datasets", type = str, nargs='+', default = "ALL", help = "name of Datasets in SETS_STARK folder conf reference arg", dest = 'datasets')
    parser_all.add_argument("-c", "--cov", type = int, default = 30, help = "Threshold coverage to keep variants in stat")
    parser_all.add_argument("-o", "--output", type = str, default = "", help = "PATH of STARK results, it will contains folder results, demultiplexing, log and tmp", dest = 'output')
    parser_all.add_argument("-r", "--reference", type = str, default = " ", help = "PATH of reference datasets Example: /app/databases   should contains folder SETS_STARK", dest = 'reference')
    parser_all.add_argument("-x", "--xternal", type = str, default = "", help = "Path of output folder for metrics", dest = 'xternal')
    parser_all.set_defaults(mode=main_all)

    #Only nonregression on run
    parser_nr = subparsers.add_parser('nr')
    parser_nr.add_argument("-o", "--output", type = str, default = "", help = "PATH of STARK run", dest = 'output')
    parser_nr.set_defaults(mode=main_nr)
  
    
    return parser.parse_args()
    
if __name__ == "__main__":
    #runCommand("figlet", "STARK -f Speed -w100 -d /figlet-fonts/ | lolcat")
    #runCommand("figlet", "'    nonregression' -f Speed -w100 -d /figlet-fonts/ | lolcat") 
    print("\n\n") 
    runCommand('bash STARK_nonregression.sh', ' --help')
    args = myoptions()
    args.mode(args)
    
    cli="stark-module-qualitycontrol-submodule-nonregression-service-cli"
    tools=['bcftools', 'samtools', 'htslib', 'vcftools'] 
    main(cli, tools, args.reference, args.datasets, args.output, args.script)   
    
    #setupAnalysis(cli=cli, tools=tools)
    
    #checkRam()
    
    
    #python3 nonreg.py -s='/app/bin/docker_VALIDATION_SETS_JB.sh' -a=' --set=TEST --folder=/app/databases --cov=30 --results=/STARK/data/TESTJB'

    #python3 nonreg.py -s="docker exec -it stark-module-qualitycontrol-submodule-nonregression-service-cli /app/bin/docker_VALIDATION_SETS_JB.sh" -a=" --set=TEST --folder=/STARK/databases/validation --cov=30 --results=/STARK/data/TESTJB"

    #python3 nonreg.py -s="/app/bin/docker_VALIDATION_SETS_JB.sh" -a=" --set=HUSHEMATO --folder=/STARK/databases/validation --cov=30 --results=/STARK/data/TESTJB" &