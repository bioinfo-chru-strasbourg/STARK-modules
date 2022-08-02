import json
import os
import subprocess
from os.path import join as osj


class config:
    def __init__(
        self, run, trio, genome, tools, output, repository, depository, config_path
    ):
        self.run = run
        self.trio = trio
        self.genome = genome
        self.tools = tools
        self.output = output
        self.repository = repository
        self.depository = depository
        self.config_path = config_path

    def configure(self):
        trio = self.trio.split(",")
        family = {}
        if type(self.tools) == dict:
            family["tools"] = self.tools
        else:
            with open(self.tools, "r") as js:
                self.tools = json.load(js)
        family["family"] = {}
        for samples in os.listdir(self.run):
            if samples in trio and os.path.isdir(osj(self.run, samples)):
                print(samples)
                samp = {}
                tagfile = self.systemcall(
                    "find "
                    + osj(self.run, samples)
                    + ' -maxdepth 2  -name "'
                    + samples
                    + '.tag"'
                )[0]
                if not tagfile:
                    # TODO search ped file
                    print("WARNING no tag file in sample " + samples)
                else:
                    with open(tagfile, "r") as t:
                        for lines in t:
                            whole = lines.strip().split("!")
                            print(whole)
                            for tag in whole:
                                if "SEX" in tag and "_" in tag:
                                    sex = tag.split("_")[-1]
                                    samp["sex"] = sex
                                elif "SEX" in tag and "#" in tag:
                                    sex = tag.split("#")[-1]
                                    samp["sex"] = sex
                                    # family['family'][samples]['sex'] = sex
                                # MOTHER, FOETUS, FATHER
                                elif "FETAL" in tag:
                                    fam = tag.split("#")[-1]
                                    if fam == "FOETUS":
                                        family["foetus"] = samp
                                        family["foetus"]["name"] = samples
                                    else:
                                        samp["affinity"] = []
                                        samp["affinity"].extend([samples, fam])
                                        # family['family'] = samp
                # Looking for bam files
                bamfile = self.systemcall(
                    "find "
                    + osj(self.run, samples)
                    + ' -maxdepth 3 -name "*.bwamem.bam" ! -name "*.validation.*"'
                )[0]
                samp["bam"] = bamfile
                if not bamfile:
                    print("ERROR no bam file in sample " + samples + " exit !")
                    exit()
                vcfile = self.systemcall(
                    "find "
                    + osj(self.run, samples)
                    + ' -maxdepth 3 -name "'
                    + samples
                    + '.final.vcf.gz"'
                )[0]
                samp["vcf"] = vcfile
                if not vcfile:
                    print("ERROR no vcf file in sample " + samples + " exit !")
                    exit()
                if "affinity" in samp:
                    s = {}
                    s[samples] = {}
                    s[samples]["bam"] = bamfile
                    family["family"][samp["affinity"][1]] = samp

        bed = self.getbed()

        family["env"] = {}
        family["env"]["output"] = self.output
        family["env"]["genome"] = self.genome
        family["env"]["bed"] = bed
        family["env"]["run"] = self.run
        family["env"]["repository"] = self.repository
        family["env"]["depository"] = self.depository

        config = self.writejson(family, self.config_path)
        # print(json.dumps(family, indent=2))
        return config

    @staticmethod
    def writejson(d, f):
        """
        Read dict, return a json file
        """
        with open(f, "w+") as outfile:
            json.dump(d, outfile, indent=4)
        print("#[INFO] Create " + f + " configfile")
        return f

    @staticmethod
    def systemcall(command):
        """
        *passing command to the shell*
        *return list containing stdout lines*
        command - shell command (string)
        """
        p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
        return p.stdout.read().decode("utf8").strip().split("\n")

    def getbed(self):
        for samples in os.listdir(self.run):
            if os.path.isdir(osj(self.run, samples)):
                bed = self.systemcall(
                    "find " + osj(self.run, samples) + ' -name "' + samples + '.bed"'
                )[0]
                if bed:
                    return bed
        print("ERROR no bed in run " + self.run + " exit !")
        exit()


# run = "/STARK/data/users/jb/run_bbs"
# trio = ["ADN2104278", "ADN2104279", "ASG184957"]
# genome = "/STARK/databases/genomes/current/hg19.fa"
# output = "/app/lib/snakemake/test_analysis.json"


# configure(run, trio, genome, tools, output)
