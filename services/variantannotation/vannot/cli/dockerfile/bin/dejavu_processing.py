import glob
from os.path import join as osj
import os
import logging as log
import shutil
import subprocess
import commons
import time
from time import sleep

import howard_launcher
import howard_processing


def convert_vcf_parquet(run_informations, args):
    threads = commons.get_threads("threads_dejavu")
    memory = commons.get_memory("memory_dejavu")
    run_dict = {}

    log.info(
        f"Generating new dejavu database for {run_informations["run_platform_application"]}"
    )
    if "run" in args:
        if not os.path.isdir(run_informations["parquet_db_run_folder"]):
            os.makedirs(run_informations["parquet_db_run_folder"], 0o775)
        else:
            shutil.rmtree(run_informations["parquet_db_run_folder"])
            os.makedirs(run_informations["parquet_db_run_folder"], 0o775)

    if not os.path.isdir(run_informations["tmp_analysis_folder"]):
        os.makedirs(run_informations["tmp_analysis_folder"], 0o775)
    else:
        shutil.rmtree(run_informations["tmp_analysis_folder"])
        os.makedirs(run_informations["tmp_analysis_folder"], 0o775)

    if "run" in args:
        vcf_files = glob.glob(osj(run_informations["archives_run_folder"], "*.vcf.gz"))
    elif "dejavu" in args:
        vcf_files = glob.glob(
            osj(run_informations["archives_project_folder"], "VCF", "*", "*.vcf.gz")
        )

    for vcf_file in vcf_files:
        run = vcf_file.split("/")[-2]
        vcf_name = vcf_file.split("/")[-1].split(".")[0]
        if run in run_dict.keys():
            previous_value = run_dict[run]
            previous_value.append(vcf_name)
            run_dict[run] = previous_value
        else:
            run_dict[run] = [vcf_name]

    for vcf_file in vcf_files:
        subprocess.call(
            ["rsync", "-rvt", vcf_file, run_informations["tmp_analysis_folder"]],
            universal_newlines=True,
        )

    if "dejavu" in args:
        for vcf_file in glob.glob(
            osj(run_informations["tmp_analysis_folder"], "*.vcf.gz")
        ):
            howard_processing.unmerge_vcf(vcf_file, run_informations)

    for vcf_file in glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz")):
        tmp_output = osj(os.path.dirname(vcf_file), "tmp_" + os.path.basename(vcf_file))
        cmd = ["bcftools", "annotate", "-x", "INFO", vcf_file]
        with open(tmp_output, "w") as writefile:
            subprocess.call(cmd, universal_newlines=True, stdout=writefile)
        os.remove(vcf_file)
        os.rename(tmp_output, vcf_file)

        exact_time = time.time() + 7200
        local_time = time.localtime(exact_time)
        actual_time = time.strftime("%H%M%S", local_time)
        start = actual_time
        container_name = f"VANNOT_dejavu_{start}_{run_informations["run_name"]}_{os.path.basename(vcf_file).split(".")[0]}"
        vcf_minimalize = osj(
            run_informations["tmp_analysis_folder"],
            f"{os.path.basename(vcf_file).split(".")[0]}.mini.parquet",
        )
        output_parquet = osj(
            run_informations["tmp_analysis_folder"],
            f"{os.path.basename(vcf_file).split(".")[0]}.parquet",
        )

        launch_minimalize_arguments = [
            "minimalize",
            "--input",
            vcf_file,
            "--output",
            vcf_minimalize,
            "--minimalize_info",
            "--minimalize_samples",
            "--threads",
            threads,
            "--memory",
            memory,
        ]
        launch_parquet_arguments = [
            "process",
            "--input",
            vcf_minimalize,
            "--output",
            output_parquet,
            "--calculations=BARCODE",
            "--explode_infos",
            "--explode_infos_fields=barcode",
            '--query=SELECT "#CHROM", POS, ID, REF, ALT, QUAL, FILTER, INFO, barcode FROM variants',
            "--threads",
            threads,
            "--memory",
            memory,
        ]
        log.info("Minimalizing vcfs")
        howard_launcher.launch(container_name, launch_minimalize_arguments)
        log.info("Converting to parquet format")
        howard_launcher.launch(container_name, launch_parquet_arguments)
        os.remove(vcf_minimalize)
        os.remove(vcf_minimalize + ".hdr")

    lockfile = osj(
        run_informations["tmp_analysis_folder"],
        "dejavu_lock_" + run_informations["run_platform_application"],
    )

    if os.path.isfile(lockfile):
        while os.path.isfile(lockfile):
            sleep(60)

    if not os.path.isfile(lockfile):
        with open(lockfile, "w") as write_file:
            pass

    for key, values in run_dict.items():
        for sample_name in values:
            parquet_file = osj(
                run_informations["tmp_analysis_folder"], sample_name + ".parquet"
            )
            parquet_file_hdr = osj(
                run_informations["tmp_analysis_folder"], parquet_file + ".hdr"
            )
            sample_dejavu_db_folder = osj(
                run_informations["parquet_db_project_folder"],
                f"RUN={key}",
                f"SAMPLE={sample_name}",
                "",
            )
            if not os.path.isdir(sample_dejavu_db_folder):
                os.makedirs(sample_dejavu_db_folder)
            subprocess.call(
                ["rsync", "-rvt", parquet_file, sample_dejavu_db_folder],
                universal_newlines=True,
            )
            subprocess.call(
                ["rsync", "-rvt", parquet_file_hdr, sample_dejavu_db_folder],
                universal_newlines=True,
            )

    # shutil.rmtree(run_informations["tmp_analysis_folder"])


def calculate_dejavu(run_informations):
    threads = commons.get_threads("threads_dejavu")
    memory = commons.get_memory("memory_dejavu")
    lockfile = osj(
        run_informations["tmp_analysis_folder"],
        "dejavu_lock_" + run_informations["run_platform_application"],
    )

    container_name = f"VANNOT_dejavu_{run_informations["run_name"]}"
    parquet_db_project = run_informations["parquet_db_howard_folder"]
    project = run_informations["run_application"]

    inner_dejavu_root_folder = os.path.dirname(parquet_db_project)
    inner_dejavu_output_parquet = osj(
        inner_dejavu_root_folder,
        f"dejavu.{run_informations["run_application"]}.parquet",
    )

    output_root_folder = os.path.dirname(run_informations["parquet_db_folder"])
    dejavu_output_parquet_hdr = osj(
        output_root_folder, f"dejavu.{run_informations["run_application"]}.parquet.hdr"
    )
    dejavu_output_parquet = osj(
        output_root_folder, f"dejavu.{run_informations["run_application"]}.parquet"
    )
    dejavu_previous_output_parquet_hdr = osj(
        output_root_folder,
        f"previous.dejavu.{run_informations["run_application"]}.parquet.hdr",
    )
    dejavu_previous_output_parquet = osj(
        output_root_folder,
        f"previous.dejavu.{run_informations["run_application"]}.parquet",
    )

    if os.path.isfile(dejavu_previous_output_parquet):
        os.remove(dejavu_previous_output_parquet)
    if os.path.isfile(dejavu_previous_output_parquet_hdr):
        os.remove(dejavu_previous_output_parquet_hdr)
    if os.path.isfile(dejavu_output_parquet):
        os.rename(dejavu_output_parquet, dejavu_previous_output_parquet)
    if os.path.isfile(dejavu_output_parquet_hdr):
        os.rename(dejavu_output_parquet_hdr, dejavu_previous_output_parquet_hdr)

    sample_count = len(
        glob.glob(
            osj(run_informations["parquet_db_project_folder"], "*", "*", "*.parquet")
        )
    )

    log.info("Calculating new frequencies")
    allelecount = "ALLELECOUNT"
    hetcount = "HETCOUNT"
    homcount = "HOMCOUNT"
    allelefreq = "ALLELEFREQ"
    query = f'SELECT "#CHROM", POS, REF, ALT, sum(CAST(barcode AS INT)) AS {allelecount}, count(barcode) FILTER(barcode=1) AS {hetcount}, count(barcode) FILTER(barcode=2) AS {homcount}, sum(CAST(barcode AS INT))/({sample_count}*2) AS {allelefreq} FROM variants WHERE PROJECT=\'{project}\' GROUP BY "#CHROM", POS, REF, ALT'
    # query = f"SELECT \"#CHROM\", POS, ANY_VALUE(ID) AS ID, REF, ALT, ANY_VALUE(QUAL) AS QUAL, ANY_VALUE(FILTER) AS FILTER, ANY_VALUE(INFO) AS INFO, sum(CAST(barcode AS INT)) AS {allelecount}, count(barcode) FILTER(barcode=1) AS {hetcount}, count(barcode) FILTER(barcode=2) AS {homcount}, sum(CAST(barcode AS INT))/({sample_count}*2) AS {allelefreq} FROM variants WHERE PROJECT='{project}' GROUP BY \"#CHROM\", POS, REF, ALT"

    launch_query_arguments = [
        "query",
        "--input",
        parquet_db_project,
        "--query",
        query,
        "--output",
        inner_dejavu_output_parquet,
        "--threads",
        threads,
        "--memory",
        memory,
    ]

    howard_launcher.launch(container_name, launch_query_arguments)
    os.remove(dejavu_output_parquet_hdr)
    with open(dejavu_output_parquet_hdr, "w") as writefile:
        writefile.write("##fileformat=VCFv4.2\n")
        writefile.write(
            '##INFO=<ID=ALLELECOUNT,Number=.,Type=Float,Description="allele count annotation">\n'
        )
        writefile.write(
            '##INFO=<ID=HETCOUNT,Number=.,Type=Float,Description="heterozygot count annotation">\n'
        )
        writefile.write(
            '##INFO=<ID=HOMCOUNT,Number=.,Type=Float,Description="homozygot count annotation">\n'
        )
        writefile.write(
            '##INFO=<ID=ALLELEFREQ,Number=.,Type=Float,Description="allele frequency annotation">\n'
        )
        writefile.write(
            f"#CHROM\tPOS\tREF\tALT\t{allelecount}\t{hetcount}\t{homcount}\t{allelefreq}\n"
        )
    os.remove(lockfile)
