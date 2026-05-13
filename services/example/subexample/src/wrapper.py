"""
@Goal: Start the subexamples pipeline and include its result properly in a STARK repository

@Author: Samuel Nicaise (2020)
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import time

from os.path import join as osj

from subexamples import main as subexample_main


def find_any_samplesheet(run_dir: str, from_res_dir=False) -> str:
    """
    Input:
            run_dir: String ; path/to/run where you want to find a samplesheet
            from_res_dir: boolean ; if True the run dir is organized as a STARK result dir,
                                                    else as a STARK repository dir (i.e. patient dirs contains a STARK dir)
    Output:
            String with samplesheet path or "NO_SAMPLESHEET_FOUND"

    Method:
    1) look up recursively all files named SampleSheet.csv in the runDir
    2) check if file path follows an expected samplesheet name and location
            (the latter depends on if we're in a STARK result or repository dir,
            defined by the bool fromResDir)
    3) first correct file path is returned
    """
    p = subprocess.Popen(
        "find -L " + run_dir + " -maxdepth 3 -name *SampleSheet.csv",
        stdout=subprocess.PIPE,
        shell=True,
    )
    out = p.stdout.readlines()
    for ss in out:
        ss = ss.decode("utf-8").strip()
        if from_res_dir:
            r = re.match(run_dir.rstrip("/") + "/(.*)/(.*).SampleSheet.csv", ss)
        else:
            r = re.match(run_dir.rstrip("/") + "/(.*)/STARK/(.*).SampleSheet.csv", ss)
        if r is None:
            continue
        elif r.group(1) == r.group(2):  # checks if (.*) == (.*)
            return ss
    return "NO_SAMPLESHEET_FOUND"


def get_sample_list_from_samplesheet(samplesheet_path: str) -> list[str]:
    """
    Returns a python list containing all sample names in a samplesheet.
    """
    assert (
        samplesheet_path != "NO_SAMPLESHEET_FOUND"
    ), "[ERROR] find_any_samplesheet() couldn't find any samplesheet. Check if the --fromResultDir argument is set correctly."
    in_data_table = False
    sample_list = []
    with open(samplesheet_path, "r") as f:
        for l in f:
            if not in_data_table:
                if l.startswith("Sample_ID,"):
                    in_data_table = True
            else:
                if "," in l:
                    sample_list.append(l.strip().split(",")[0])
    # if there are spaces in samplesheet names, change them to "_" because that's what demultiplexing.sh will do
    # otherwise the fastq won't be found when looking in the DEM dir
    for i in range(len(sample_list)):
        if " " in sample_list[i]:
            sample_list[i] = sample_list[i].replace(" ", "_")
    return sample_list


def clean_list(str) -> list[str]:
    return [x.strip() for x in str.split(",")]


def is_valid_sample(run_dir: str, sample: str, exclude: list[str]) -> bool:
    with open(osj(run_dir, sample, "STARK", sample + ".tag"), "r") as f:
        for l in f:
            for item in exclude:
                if item in l:
                    return False
    return True


def tag_field_to_dict(tag_string: str) -> dict[str, list[str]]:
    """
    >>> tag_field_to_dict("APP#HUS")
    {'APP': ['HUS']}
    >>> tag_field_to_dict("APP#HUSTUMSOL.XTHS!#test!ceci est une description")
    {'': ['test'], 'APP': ['HUSTUMSOL.XTHS']}
    >>> tag_field_to_dict("SEX#F!APP#DIAG.DI#subexample!")
    {'APP': ['DIAG.DI', 'subexample'], 'SEX': ['F']}

    Trailing "!" are authorized and must be dealt with.
    Fields without a # are not returned.
    """
    tag_dict = {}
    tag_string = tag_string.rstrip("\r\n")  # as it will usually come from a samplesheet
    if "#" not in tag_string:
        return {}
    if "!" not in tag_string:
        v = tag_string.split("#")
        if len(v) == 2:
            tag_dict[v[0]] = [v[1]]
        elif len(v) > 2:
            tag_dict[v[0]] = v[1:]
    else:
        for ts in tag_string.split("!"):
            if ts == "":  # deals with trailing "!"
                continue
            v = ts.split("#")
            if len(v) == 2:
                tag_dict[v[0]] = [v[1]]
            elif len(v) > 2:
                tag_dict[v[0]] = v[1:]
    return tag_dict


def is_subexample(run_dir: str, sample: str, sex_tag: str) -> bool:
    """
    checks if sample is a subexample of samples of given sex according to sample tag file
    """
    is_subexample = False
    correct_sex = False
    with open(osj(run_dir, sample, "STARK", sample + ".tag"), "r") as f:
        for l in f:
            tags_dict = tag_field_to_dict(l)
            for key in tags_dict.keys():
                if key == "PLUGAPP":
                    if "subexample" in tags_dict[key]:
                        is_subexample = True
                elif key == "APP":
                    if "subexample" in tags_dict[key]:
                        is_subexample = True
                if key == "SEX":
                    if tags_dict[key] == tag_field_to_dict(sex_tag)[key]:
                        correct_sex = True
    if is_subexample and correct_sex:
        return True
    else:
        return False


def get_subexample_dict(sample_list: list[str], run_dir: str) -> dict[str, list[str]]:
    subexample_dict = {}
    for s in sample_list:
        with open(osj(run_dir, s, "STARK", s + ".tag"), "r") as f:
            for l in f:
                if "subexample" in tag_field_to_dict(l).keys():
                    subexample_index = tag_field_to_dict(l)["subexample"]
                    subexample_index.sort()
                    if not "#".join(subexample_index) in subexample_dict.keys():
                        subexample_dict["#".join(subexample_index)] = []
                    subexample_dict["#".join(subexample_index)].append(s)
    if not subexample_dict:
        subexample_list = []
        for s in sample_list:
            with open(osj(run_dir, s, "STARK", s + ".tag"), "r") as f:
                for l in f:
                    if is_subexample(run_dir, s, "SEX#F"):
                        subexample_list.append(s)
                    if is_subexample(run_dir, s, "SEX#M"):
                        subexample_list.append(s)
        subexample_list.sort()
        subexample_dict["#".join(subexample_list)] = sample_list
    print("getsubexampleDict: ", subexample_dict)
    return subexample_dict


def create_sample_repository(run_dir: str, sample: str) -> tuple[str, str]:
    res_dir = osj(run_dir, sample, "subexample")
    log_dir = osj(res_dir, "logs")
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    return (res_dir, log_dir)


def write_error_log(sample: str, run_dir: str, error_msg: str) -> None:
    res_dir, log_dir = create_sample_repository(run_dir, sample)
    with open(osj(res_dir, "ERROR.log"), "w") as f:
        f.write(error_msg)


def copy_results(sample_list: list[str], run_dir: str, docker_output_dir: str) -> None:
    print("hello copyResults")
    print("sampleList:", sample_list)
    print("runDir:", run_dir)
    print("dockerOutputDir:", docker_output_dir)
    for s in sample_list:
        if os.path.exists(osj(docker_output_dir, s + ".final.vcf.gz")):
            res_dir, log_dir = create_sample_repository(run_dir, s)
            shutil.copyfile(
                osj(docker_output_dir, s + ".final.vcf.gz"),
                osj(res_dir, s + ".final.vcf.gz"),
            )
        else:
            print(
                f"[ERROR] Missing final VCF file for sample {s} in docker output directory."
            )
            write_error_log(s, run_dir, "[ERROR] Missing final VCF file.")


def launch_analysis(
    sample_list: list[str], key: str, run_dir: str, work_dir: str, genome: str
) -> str | None:
    vcf_list = []
    subexample_F_str = "init"
    subexample_M_str = "init"
    bed = "init"
    for s in sample_list:
        vcf = osj(run_dir, s, "STARK", s + ".reports", s + ".final.vcf")
        vcf_list.append(vcf)
    for subexample in key.split("#"):
        if os.path.exists(osj(run_dir, subexample, "STARK", subexample + ".tag")):
            if is_subexample(run_dir, subexample, "SEX#M"):
                assert (
                    subexample_M_str == "init"
                ), "[ERROR] More than one sample is named subexample_([A-Z]*)_M_([0-9]*) in the samplesheet"
                vcf = osj(
                    run_dir, subexample, "STARK", subexample + ".reports", subexample + ".final.vcf"
                )
                bam = osj(run_dir, subexample, "STARK", subexample + ".bwamem.bam")
                vcf_list.append(vcf)
                subexample_M_str = ":".join([subexample, vcf, bam])
            elif is_subexample(run_dir, subexample, "SEX#F"):
                assert (
                    subexample_F_str == "init"
                ), "[ERROR] More than one sample is named subexample_([A-Z]*)_F_([0-9]*) in the samplesheet"
                vcf = osj(
                    run_dir, subexample, "STARK", subexample + ".reports", subexample + ".final.vcf"
                )
                bam = osj(run_dir, subexample, "STARK", subexample + ".bwamem.bam")
                vcf_list.append(vcf)
                subexample_F_str = ":".join([subexample, vcf, bam])
            else:
                for s in sample_list:
                    write_error_log(
                        s, run_dir, "[ERROR] Can't find " + subexample + " sample's sex."
                    )
                return "[ERROR] Can't find " + subexample + " sample's sex."
        else:
            for s in sample_list:
                write_error_log(
                    s, run_dir, "[ERROR] Specified " + subexample + " sample not found."
                )
            return "[ERROR] Specified " + subexample + " sample not found."
    bed = osj(run_dir, subexample, "STARK", subexample + ".bed")
    if bed == "init":
        for s in sample_list:
            write_error_log(s, run_dir, "[ERROR] Missing bed for subexample analysis.")
        return "[ERROR] Missing bed for subexample analysis."
    # TODO: be able to fetch a subexample from a different run ?
    # cmd = 'python /app/lib/subexample/subexample.py sample -o '+dockerOutputDir+' -s "'+','.join(vcfList)+'" -p "'+subexampleFStr+","+subexampleMStr+'" -b '+bed+' -g '+genome
    print("Launching")
    subexample_main(work_dir, ",".join(vcf_list), subexample_F_str, subexample_M_str, bed, genome)
    print()
    # subprocess.call(cmd, shell=True)
    copy_results(sample_list, run_dir, work_dir)


def main(args: argparse.Namespace) -> None:
    run_work_dir = osj(args.workDir, os.path.basename(args.runDir))
    if not os.path.exists(run_work_dir):
        os.mkdir(run_work_dir)

    # get only samples that will be analysed
    sample_list = get_sample_list_from_samplesheet(find_any_samplesheet(args.runDir))
    removed_samples = []
    for s in sample_list:
        if not is_valid_sample(args.runDir, s, clean_list(args.exclude)):
            removed_samples.append(s)
        elif osj(args.runDir, s) not in glob.glob(osj(args.runDir, "*")):
            print(
                "[WARNING] Sample "
                + s
                + " from samplesheet not in repository, ignoring it. Could be normal if it belongs to a different application"
            )
            removed_samples.append(s)
    for removed in removed_samples:
        sample_list.remove(removed)
    # create a dictionary with subexample_ID1#subexample_ID2 as key, and samples as values, depending on tag subexample#subexample_ID1#subexample_ID2
    subexample_dict = get_subexample_dict(sample_list, args.runDir)
    for key in subexample_dict:
        launch_analysis(subexample_dict[key], key, args.runDir, run_work_dir, args.genome)

    shutil.rmtree(run_work_dir)

    with open(osj(args.runDir, "subexampleComplete.txt"), "w") as f:
        f.write(time.ctime())
    if os.path.exists(osj(args.runDir, "subexampleRunning.txt")):
        os.remove(osj(args.runDir, "subexampleRunning.txt"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="wrapper to launch routine analysis")
    parser.set_defaults(mode=main)
    parser.add_argument(
        "-i",
        "--runDir",
        type=str,
        help="path to run in a STARK 0.9.18 repository",
        required=True,
    )
    parser.add_argument(
        "-g", "--genome", help="genome file", type=str, dest="genome", required=True
    )
    parser.add_argument(
        "-e",
        "--exclude",
        help="list of tags for which samples are removed if have them",
        type=str,
        dest="exclude",
        default="CQI#",
    )
    parser.add_argument(
        "-w",
        "--workDir",
        help="work directory",
        type=str,
        dest="workDir",
        default="/dev/shm/subexamples",
    )

    args = parser.parse_args()
    args.mode(args)
