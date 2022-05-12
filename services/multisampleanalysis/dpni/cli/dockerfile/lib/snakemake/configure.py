import argparse
import json
import os
import subprocess
from os.path import join as osj

tools = {
    "java": "java",
    "gatk": "/STARK/tools/gatk/3.8-1-0/bin/GenomeAnalysisTK.jar",
    "VaRank": "/home1/TOOLS/tools/varank/VaRank_1.4.3",
    "Alamut-batch": "/home1/TOOLS/tools/alamut_batch/alamut-batch-standalone-1.11",
    "howard": "/STARK/tools/howard/current/bin/HOWARD",
    "bcftools": "/STARK/tools/bcftools/current/bin/bcftools",
    "bgzip": "/STARK/tools/htslib/current/bin/bgzip",
    "tabix": "/STARK/tools/htslib/current/bin/tabix",
    "threads": 6,
}


def configure(run, trio_tmp, genome, tools, output):
    trio = trio_tmp.split(",")
    family = {}
    family["tools"] = tools
    family["family"] = {}
    for samples in os.listdir(run):
        if samples in trio and os.path.isdir(osj(run, samples)):
            samp = {}
            tagfile = systemcall(
                "find "
                + osj(run, samples)
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
                            elif "FOETAL" in tag:
                                fam = tag.split("#")[-1]
                                if fam == "FOETUS":
                                    family["foetus"] = samp
                                else:
                                    samp["affinity"] = []
                                    samp["affinity"].extend([samples, fam])
                                    # family['family'] = samp
            # Looking for bam files
            bamfile = systemcall(
                "find "
                + osj(run, samples)
                + ' -maxdepth 3 -name "*.bwamem.bam" ! -name "*.validation.*"'
            )[0]
            samp["bam"] = bamfile
            if not bamfile:
                print("ERROR no bam file in sample " + samples + " exit !")
                exit()
            vcfile = systemcall(
                "find "
                + osj(run, samples)
                + ' -maxdepth 3 -name "'
                + samples
                + '.final.vcf"'
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

    bed = getbed(run)

    family["env"] = {}
    family["env"]["output"] = "/app/res"
    family["env"]["genome"] = genome
    family["env"]["bed"] = bed

    config = writejson(family, output)
    # print(json.dumps(family, indent=2))
    return config


def writejson(d, f):
    """
    Read dict, return a json file
    """
    with open(f, "w") as outfile:
        json.dump(d, outfile, indent=4)


def systemcall(command):
    """
    *passing command to the shell*
    *return list containing stdout lines*
    command - shell command (string)
    """
    p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
    return p.stdout.read().decode("utf8").strip().split("\n")


def getbed(run):
    for samples in os.listdir(run):
        if os.path.isdir(osj(run, samples)):
            bed = systemcall(
                "find " + osj(run, samples) + ' -name "' + samples + '.bed"'
            )[0]
            if bed:
                return bed
    print("ERROR no bed in run " + run + " exit !")
    exit()


def parseargs():
    """
    parse arguments
    """
    parser = argparse.ArgumentParser(
        # formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-r", "--run", type=str, help="Absolute path of run folder repository format"
    )
    parser.add_argument(
        "-o", "--output", type=str, help="Path of the output config.json"
    )
    parser.add_argument(
        "-g",
        "--genome",
        type=str,
        default="/STARK/databases/genomes/current/hg19.fa",
        help="Absolute path of genome assemby fasta file",
    )
    parser.add_argument(
        "--tools",
        type=dict,
        default=tools,
        help="Dict containing path of tools executable",
    )
    parser.add_argument(
        "-t", "--trio", type=str, help="List of sample representing example s1,s2,s3"
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parseargs()
    configure(args.run, args.trio, args.genome, args.tools, args.output)


# run = "/STARK/data/users/jb/run_bbs"
# trio = ["ADN2104278", "ADN2104279", "ASG184957"]
# genome = "/STARK/databases/genomes/current/hg19.fa"
# output = "/app/lib/snakemake/test_analysis.json"


# configure(run, trio, genome, tools, output)
