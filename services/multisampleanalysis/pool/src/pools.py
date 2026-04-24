"""
This script replaces the POOL module.
previously launched with
cmd = 'python /app/lib/pool/pool.py sample -o '+dockerOutputDir+' -s "'+','.join(vcfList)+'" -p "'+poolFStr+","+poolMStr+'" -b '+bed+' -g '+genome
"""

from dataclasses import dataclass
import logging as log
import os
from pathlib import Path

from cyvcf2 import cyvcf2
from typing import TypedDict
import subprocess
from multiprocessing import Pool as ProcessPool  # too many things are named Pool here


class AnnotationData(TypedDict, total=False):
    GT: str
    DP: int
    base_counts: str


DataByVariantKey = dict[str, AnnotationData]
DataByChromosome = dict[str, DataByVariantKey]


@dataclass
class Pool:
    name: str
    vcf_path: Path
    bam_path: Path

    @staticmethod
    def from_string(pool_str: str) -> "Pool | None":
        if pool_str == "init":
            return None
        parts = pool_str.split(":")
        if len(parts) != 3:
            raise ValueError(
                f"Invalid pool string format. Expected: either 'init' (no pool) or '<pool_id>:<path_to_vcf>:<path_to_bam>' ; Got instead: {pool_str}"
            )
        return Pool(name=parts[0], vcf_path=Path(parts[1]), bam_path=Path(parts[2]))


def get_variant_key(variant: cyvcf2.Variant) -> str:
    return f"{variant.CHROM}:{variant.POS}:{variant.REF}:{','.join(variant.ALT)}"


def get_pool_data(
    pool: Pool, work_dir: str, all_variants: DataByChromosome
) -> DataByChromosome:
    """
    Get the necessary data from a Pool.

    Returns:
    A dictionary where the keys are variant IDs and the values are lists of annotations for that variant in the pool VCF.
    Relevant data: GT, depth, base counts
    """
    # data = {}
    if pool is not None:
        # With gts012=True, gt_types will be 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN. If False, 3 and 2 are flipped.
        # Here, it is expected to be set to True.
        vcf = cyvcf2.VCF(pool.vcf_path, gts012=True)
        if len(vcf.samples) != 1:
            raise ValueError(
                f"Expected exactly one sample in the pool VCF, but found {len(vcf.samples)} in {pool.vcf_path}"
            )
        for variant in vcf:
            if variant.CHROM not in all_variants:
                raise ValueError(
                    f"Variant chromosome {variant.CHROM} not found in the merged variants from the sample VCFs. This should not happen since the merged variants are coming from the sample VCFs and the pool VCF. Variant: {get_variant_key(variant)}"
                )
            key = get_variant_key(variant)
            gt = str(variant.gt_types[0])
            # the previous pipeline fetched DP data from GATK too, so let's keep using it instead of the VCF
            # dp = int(variant.format("DP")[0][0])
            all_variants[variant.CHROM][key] = AnnotationData(GT=gt)

    data = add_base_counts_to_data(pool, work_dir, all_variants)

    return data


def add_base_counts_to_data(
    pool: Pool, work_dir: str, data: DataByChromosome
) -> DataByChromosome:
    # 1) get gakt DepthOfCoverage data
    # 2) parse results and add data to the data dict if coordinates match

    for chrom, variants in data.items():
        coverage_file = f"{work_dir}/{pool.name}/depth_of_coverage_{chrom}"

        coverage_data = {}
        with open(coverage_file, "r") as f:
            for line in f:
                if line.startswith(f"{chrom}"):
                    columns = line.strip().split(",")
                    pos = columns[0].split(":")[1]
                    depth = int(columns[3])
                    base_counts = columns[4]
                    # replace spaces with pipes to keep current formating
                    base_counts = base_counts.replace(" ", "|")
                    coverage_data[pos] = (depth, base_counts)

        # print(f"Coverage data for chromosome {chrom} in pool {pool.name if pool is not None else 'None'}: {dict(list(coverage_data.items())[:10])}")

        for key in variants.keys():
            _, pos, _, _ = key.split(":")
            if pos in coverage_data:
                depth, base_counts = coverage_data[pos]
                data[chrom][key]["DP"] = depth
                data[chrom][key]["base_counts"] = base_counts

    print("helloooo there @@@@@@@@@@@")
    try:
        print("checking", data["chr1"]["chr1:1454424:C:T"])
    except:
        pass
    # print(f"Data for pool {pool.name if pool is not None else 'None'}: {dict(list(data.items())[:3])}")
    return data


def run_gatk(args) -> None:
    interval_file, work_dir, pool_name, bam_path, genome = args
    chrom = Path(interval_file).stem
    output_prefix = f"{work_dir}/{pool_name}/depth_of_coverage_{chrom}"
    log_file = Path(work_dir) / pool_name / f"gatk_{chrom}.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)

    command = (
        f"docker compose -f /app/src/docker-compose.yml --env-file /app/.env run --rm "
        f"stark-module-multisampleanalysis-submodule-pool-service-gatk gatk DepthOfCoverage "
        f"-I {bam_path} "
        f"-L {interval_file} "
        f"-O {output_prefix} "
        f"-R {genome} "
        f"--print-base-counts "
        f"--include-deletions "
        f"--omit-interval-statistics "
        f"--interval-padding 50 "
        f"--read-filter MappingQualityReadFilter "
        f"--minimum-mapping-quality 10"
    )
    with open(log_file, "a") as log:
        subprocess.run(command, shell=True, check=True, stdout=log, stderr=log)


def generate_cov_data(pool: Pool, work_dir: str, bed: str, genome: str) -> None:
    """
    Run GATK DepthOfCoverage, parallelized by chromosome, to get the base counts on every position within the bed file.

    Generate "{pool.name}_coverage_done.txt" in the work_dir when done. This is used to check if the data was already generated for a given pool, since the wrapper can launch an analysis containing the same pool twice. This happens when some samples are linked to only one pool, and others are linked to both parental pools, leading to two separate analysis.
    """
    intervals_dir = f"{work_dir}/intervals"
    os.makedirs(intervals_dir, exist_ok=True)

    # Create interval files for each chromosome
    interval_files = []
    with open(bed, "r") as bed_file:
        for line in bed_file:
            chrom, start, end = line.strip().split()[:3]
            interval_file = f"{intervals_dir}/{chrom}.list"
            if interval_file not in interval_files:
                interval_files.append(interval_file)
            with open(interval_file, "a") as f:
                f.write(f"{chrom}:{start}-{end}\n")

    with ProcessPool() as pool_executor:
        pool_executor.map(
            run_gatk,
            [
                (interval_file, work_dir, pool.name, pool.bam_path, genome)
                for interval_file in interval_files
            ],
        )

    done_file_path = Path(work_dir) / f"{pool.name}_coverage_done.txt"
    with open(done_file_path, "w") as done_file:
        done_file.write("DONE\n")


def create_output_vcf(
    index_vcf_path: Path,
    output_vcf_path: Path,
    pool_f_data: DataByChromosome,
    pool_m_data: DataByChromosome,
) -> None:
    input_vcf = cyvcf2.VCF(str(index_vcf_path), gts012=True)
    if len(input_vcf.samples) != 1:
        raise ValueError(
            f"Expected exactly one sample, but found {len(input_vcf.samples)} in {index_vcf_path}"
        )

    input_vcf.add_info_to_header(
        {
            "ID": "BARCODE",
            "Description": "POOL Barcode ordered as: index, maternal pool, paternal pool",
            "Type": "String",
            "Number": "1",
        }
    )
    input_vcf.add_info_to_header(
        {
            "ID": "POOL_F_Depth",
            "Description": "Depth of the variant in the maternal pool",
            "Type": "Integer",
            "Number": "1",
        }
    )
    input_vcf.add_info_to_header(
        {
            "ID": "POOL_F_BASE_COUNTS",
            "Description": "Base counts of the variant in the maternal pool",
            "Type": "String",
            "Number": "1",
        }
    )
    input_vcf.add_info_to_header(
        {
            "ID": "POOL_M_Depth",
            "Description": "Depth of the variant in the paternal pool",
            "Type": "Integer",
            "Number": "1",
        }
    )
    input_vcf.add_info_to_header(
        {
            "ID": "POOL_M_BASE_COUNTS",
            "Description": "Base counts of the variant in the paternal pool",
            "Type": "String",
            "Number": "1",
        }
    )

    output_vcf = cyvcf2.Writer(str(output_vcf_path), input_vcf)
    for variant in input_vcf:
        key = get_variant_key(variant)
        gt = str(variant.gt_types[0])

        # delete all the initial INFO annotation to save space
        info_keys = [k[0] for k in variant.INFO]
        for info_key in info_keys:
            del variant.INFO[info_key]

        if pool_f_data != {}:
            variant_data = pool_f_data[variant.CHROM].get(key, {})
            gt_pool_f = variant_data.get("GT", "0")
            variant.INFO["POOL_F_Depth"] = variant_data.get("DP", "")
            variant.INFO["POOL_F_BASE_COUNTS"] = variant_data.get("base_counts", "")
        else:
            gt_pool_f = "0"

        if pool_m_data != {}:
            variant_data = pool_m_data[variant.CHROM].get(key, {})
            gt_pool_m = variant_data.get("GT", "0")
            variant.INFO["POOL_M_Depth"] = variant_data.get("DP", "")
            variant.INFO["POOL_M_BASE_COUNTS"] = variant_data.get("base_counts", "")
        else:
            gt_pool_m = "0"

        barcode = f"{gt}{gt_pool_f}{gt_pool_m}"
        variant.INFO["BARCODE"] = barcode

        output_vcf.write_record(variant)
    output_vcf.close()


def merge_variants(vcf_list: list[str]) -> DataByChromosome:
    """
    Merge the variants from all the VCF paths in input. Used to have a list of variants from which coverage data should be fetched.
    Fetching coverage data from variants present in pools only is not enough, as some variants can be present in index samples but not in pools.
    """
    merged_data = {}
    for vcf_path in vcf_list:
        vcf = cyvcf2.VCF(vcf_path, gts012=True)
        if len(vcf.samples) != 1:
            raise ValueError(
                f"Expected exactly one sample in the VCF, but found {len(vcf.samples)} in {vcf_path}"
            )
        for variant in vcf:
            chrom = variant.CHROM
            key = get_variant_key(variant)
            if chrom not in merged_data:
                merged_data[chrom] = {}
            if key not in merged_data[chrom]:
                merged_data[chrom][key] = {}
    return merged_data


def main(
    output_dir: str,
    vcf_list_as_str: str,
    pool_f_str: str,
    pool_m_str: str,
    bed: str,
    genome: str,
) -> None:
    """
    This function is called by the wrapper for each pool analysis to be performed.

    In a given run, there can be multiple pools, and each sample can be attached to one or two pools.
    With 2 pairs of pools in a run there can be up to 6 pool analysis: with one or the other parent, or with the two parents, for each pair of pool.

    Args:
    output_dir: in-container output dir. This is in fact more of a work dir, and the wrapper will copy the results to the actual output dir.
    vcf_list_as_str: comma-separated list of VCF paths of samples attached to the pool(s) being analyzed.
    pool_f_str: identify the maternal pool being analyzed. Each string has the following form:
        <pool_id>:<path_to_vcf>:<path_to_bam>
        For example:
        POOL_F:/STARK/data/260217_NB551000_0001_AHLCWGZZZZ/POOL_F_130/STARK/POOL_F_130.reports/POOL_F_130.final.vcf:/STARK/data/260217_NB551000_0001_AHLCWGZZZZ/POOL_F_130/STARK/POOL_F_130.bwamem.bam
        If there is no maternal pool, this string is set as "init" by the wrapper.
    pool_m_str: identify the paternal pool being analyzed. Same format as pool_f_str. Same "init" behavior if there no paternal pool.
    bed: run bed path
    genome: STARK genome path

    The following results are copied by the wrapper to the run repository, hence are expected:
        for s in sampleList:
                        osj(dockerOutputDir, s+".final.vcf.gz")
                        osj(dockerOutputDir, s+".final.vcf.gz.tbi")
        + log files that I will have to do
    """
    vcf_list = vcf_list_as_str.split(",")

    pool_f = Pool.from_string(pool_f_str)
    pool_m = Pool.from_string(pool_m_str)

    for p in [pool_f, pool_m]:
        if p is not None:
            done_file_path = Path(output_dir) / f"{p.name}_coverage_done.txt"
            if not done_file_path.exists():
                generate_cov_data(p, output_dir, bed, genome)

    if pool_f is not None:
        all_variants = merge_variants(vcf_list + [str(pool_f.vcf_path)])
        pool_f_data = get_pool_data(pool_f, output_dir, all_variants)
    else:
        pool_f_data = {}

    if pool_m is not None:
        all_variants = merge_variants(vcf_list + [str(pool_m.vcf_path)])
        pool_m_data = get_pool_data(pool_m, output_dir, all_variants)
    else:
        pool_m_data = {}

    for vcf_path in vcf_list:
        print(f"Processing sample VCF: {vcf_path}")
        sample_vcf = Path(vcf_path)
        output_vcf = Path(output_dir) / f"{sample_vcf.stem}.vcf.gz"
        create_output_vcf(sample_vcf, output_vcf, pool_f_data, pool_m_data)

    print("Pool analysis completed successfully.")


# TODO: split base counts with pipes
# TODO: explain empty data on variant SGT2000780 chr1	1454424
