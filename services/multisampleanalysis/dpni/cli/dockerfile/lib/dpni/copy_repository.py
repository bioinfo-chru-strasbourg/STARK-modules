import subprocess
import os
from os.path import join as osj
import json
from datetime import datetime
import sys

"""
Copy result of DPNI module into folder

example:
python copy_repository.py <resultpath> <configjson>
"""


date = datetime.now().strftime("%Y%m%d-%H%M%S")

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

print(sys.argv[1])
input = [osj(sys.argv[1], items) for items in os.listdir(sys.argv[1]) if items.endswith(".final.vcf")]
cg = sys.argv[2]
with open(cg, "r") as f:
	config = json.load(f)


ps = {
	"bgzip": config["tools"]["bgzip"],
	"tabix": config["tools"]["tabix"],
	"bcftools": config["tools"]["bcftools"]
}

params = Struct(**ps)

print(params.bgzip)
print(config)
print(input)

whole = ["FATHER", "MOTHER", "FETUS"]


def systemcall(command):
    '''
    *passing command to the shell*
    *return list containing stdout lines*
    command - shell command (string)
    '''
    p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
    return p.stdout.read().decode('utf8').strip().split('\n')

def copy_repository():
    if os.path.exists(config["env"]["repository"]):
        print("#[INFO] Copying into repository "+config["env"]["repository"])
        for sample in os.listdir(config["env"]["repository"]):
            if not os.path.isdir(osj(config["env"]["repository"], sample)):
                continue #not a sample
            print("#[INFO] Sample ", sample)
            if os.path.isdir(osj(osj(config["env"]["repository"], sample))):
                #In case of first analysis
                if not os.path.exists(osj(config["env"]["repository"], sample, 'DPNI')):
                    print("#[INFO] First analysis create DPNI folder")
                    os.mkdir(osj(config["env"]["repository"], sample, 'DPNI'))
                else:
                    if not os.path.exists(osj(config["env"]["repository"], sample, "DPNI", "DPNI."+date)):
                        os.mkdir(osj(config["env"]["repository"], sample, "DPNI", "DPNI."+date))
                    #move old analysis files if they exists in DPNI date folder
                    systemcall("mv "+osj(config["env"]["repository"], sample, "DPNI/*.vcf.gz*")+" "+osj(config["env"][  "repository"], sample, "DPNI", "DPNI."+date))
                    #Same with log folder
                    systemcall("mv "+osj(config["env"]["repository"], sample, "DPNI/.log")+" "+osj(config["env"]["repository"], sample, "DPNI", "DPNI."+date))
                #Copy new vcf analysis and index into repo sample
                print("#[INFO] Copying sample "+sample)
                systemcall("rsync "+osj(config["env"]["output"], sample+".final.vcf.gz*")+" "+osj(config["env"]["repository"], sample, 'DPNI/'))
                #Copying new logs
                systemcall("rsync -ar "+osj(config["env"]["output"], ".log")+" "+osj(config["env"]["repository"], sample, 'DPNI/'))
            if config["env"]["depository"]:
                print("#[INFO] Copying into depository, sample ", sample)
                systemcall("rsync -ar "+osj(config["env"]["repository"], sample, 'DPNI')+" "+osj(config["env"]["depository"], sample+"/"))
    else:
        print("ERROR "+config["env"]["repository"]+" path does not exist EXIT")

#Compress and Index (same cmd in bash above)
for vcf in input:
	print(vcf)
	alias = os.path.basename(vcf).split(".", 1)[0]
	if alias in ["MOTHER", "FATHER"]:
		print(params.bgzip+" --force "+vcf+" -c > "+vcf+".tmp.gz && "+params.tabix+" -p vcf "+vcf+".tmp.gz")
		systemcall(params.bgzip+" --force "+vcf+" -c > "+vcf+".tmp.gz && "+params.tabix+" -p vcf "+vcf+".tmp.gz")
		id = config["family"][alias]["affinity"][0]
		print(alias)
		print(id)
		systemcall("echo "+id+" > "+osj(config["env"]["output"], alias))
		systemcall(params.bcftools+" reheader -s "+osj(config["env"]["output"], alias)+" "+vcf+".tmp.gz -o "+osj(config["env"]["output"], id+".final.vcf.gz")+" && "+params.tabix+" -p vcf "+osj(config["env"]["output"], id+".final.vcf.gz"))
	else:
		systemcall(params.bgzip+" --force "+vcf+" && "+params.tabix+" -p vcf "+vcf+".gz")
	copy_repository()
