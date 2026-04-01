#!/usr/bin/env python3
"""
CQI - Quality control analysis
Author : Jean-Baptiste Lamouche
Copyright: HUS
Licence: GNU-GPL
Release: 2.0
Date: 2026
"""

import argparse
import json
import logging
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

DATEFILE = datetime.now().strftime("%Y%m%d-%H%M%S")
RELEASE_NOTES = """
0.9.18 : Script creation
1.0    : Full refactor
2.0    : Python rewrite using vcftoolz
"""

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging(log_path: Path) -> logging.Logger:
    log = logging.getLogger("CQI")
    log.setLevel(logging.DEBUG)
    fmt = logging.Formatter("[%(levelname)s] %(message)s")
    # file handler
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    log.addHandler(fh)
    # stderr handler (errors only)
    sh = logging.StreamHandler(sys.stderr)
    sh.setLevel(logging.ERROR)
    sh.setFormatter(fmt)
    log.addHandler(sh)
    return log

# ---------------------------------------------------------------------------
# Shell helpers
# ---------------------------------------------------------------------------

def run(cmd: str, log: logging.Logger, check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command, stream stderr to log file."""
    log.debug(f"CMD: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.stderr:
        log.debug(result.stderr.strip())
    if check and result.returncode != 0:
        log.error(f"Command failed (rc={result.returncode}): {cmd}")
        raise RuntimeError(f"Command failed: {cmd}")
    return result

def index_vcf(vcf_gz: Path, log: logging.Logger):
    run(f"tabix -f -p vcf {vcf_gz}", log)

def compress_and_index(infile: Path, outfile: Path, log: logging.Logger):
    if infile.is_file() and infile.stat().st_size > 0:
        run(f"bgzip -c {infile} > {outfile}", log)
    else:
        # Empty VCF placeholder
        empty = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        run(f"echo -e '{empty}' | bgzip -c > {outfile}", log)
    index_vcf(outfile, log)

def norm_and_index(infile: Path, outfile: Path, genome: Path, log: logging.Logger) -> Path:
    """Normalize VCF; returns outfile path on success, infile path on failure."""
    run(f"bcftools norm -f {genome} {infile} -O z -o {outfile}", log, check=False)
    if outfile.is_file() and outfile.stat().st_size > 0:
        index_vcf(outfile, log)
        return outfile
    return infile

def count_variants(vcf: Path, log: logging.Logger) -> int:
    if not vcf.is_file() or vcf.stat().st_size == 0:
        return 0
    r = run(f"bcftools view -H {vcf} 2>/dev/null | wc -l", log)
    return int(r.stdout.strip())

def pct(num: int, den: int) -> str:
    if den == 0:
        return "0.00"
    return f"{(num / den) * 100:.2f}"

# ---------------------------------------------------------------------------
# Interval helpers
# ---------------------------------------------------------------------------

def prepare_intervals(bed: str, cqi_dir: Path, log: logging.Logger) -> Path:
    if "," in bed:
        out = cqi_dir / f"CQI.{DATEFILE}.intervals.genes.bed"
        bed_files = " ".join(bed.split(","))
        run(
            f"cat {bed_files} | sort -k1,1V -k2,2n | bedtools merge -i stdin > {out}",
            log
        )
        if not out.is_file() or out.stat().st_size == 0:
            log.error(f"Gene file empty: {out}")
            sys.exit(2)
        return out
    return Path(bed)

# ---------------------------------------------------------------------------
# BCF preparation pipeline
# ---------------------------------------------------------------------------

def prepare_bcf(invcf: Path, genes: Path, genome: Path, log: logging.Logger, out_bcf: Path):
    """reheader -> intersect -> sort -> BCF"""
    run(
        f"bcftools view {invcf} | "
        f"bcftools reheader --fai {genome}.fai | "
        f"bedtools intersect -a stdin -b {genes} -header -u | "
        f"bcftools sort -O b -o {out_bcf}",
        log
    )

def split_bcf(bcf: Path, prefix: str, snv_dir: Path, indel_dir: Path, log: logging.Logger):
    """Split one BCF into SNV + INDEL vcf.gz."""
    for folder, flag in [(snv_dir, "--types snps"), (indel_dir, "--exclude-types snps")]:
        out = folder / f"{prefix}.vcf.gz"
        run(f"bcftools view {flag} -O z -o {out} {bcf}", log)
        index_vcf(out, log)

# ---------------------------------------------------------------------------
# vcftoolz comparison (replaces bcftools isec)
# ---------------------------------------------------------------------------

def run_vcftoolz(cqi_vcf: Path, sample_vcf: Path, out_dir: Path, log: logging.Logger) -> dict:
    """
    Use vcftoolz compare to intersect two VCFs.
    Returns dict with keys: expected, found, positive, missing, noise
    vcftoolz compare <truth> <query> outputs:
        - only in truth  -> missing  (0000.vcf equivalent)
        - only in query  -> noise    (0001.vcf equivalent)
        - in both        -> positive (0002.vcf / 0003.vcf equivalent)
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    missing_vcf  = out_dir / "missing.vcf"
    noise_vcf    = out_dir / "noise.vcf"
    positive_vcf = out_dir / "positive.vcf"

    run(
        f"vcftoolz compare {cqi_vcf} {sample_vcf} "
        f"--unique-to-first {missing_vcf} "
        f"--unique-to-second {noise_vcf} "
        f"--intersection {positive_vcf}",
        log
    )

    def _count(f: Path) -> int:
        return count_variants(f, log) if f.is_file() else 0

    pos = _count(positive_vcf)
    exp = count_variants(cqi_vcf, log)
    fnd = count_variants(sample_vcf, log)

    return {
        "expected"  : exp,
        "found"     : fnd,
        "positive"  : pos,
        "missing"   : exp - pos,
        "noise"     : fnd - pos,
        "missing_vcf": missing_vcf,
    }

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def write_report_header(out: Path, run_dir: Path, sample: str):
    with open(out, "w") as f:
        f.write("##########################\n")
        f.write("### CQI vcf comparison \n")
        f.write(f"### RUN: {run_dir.name}\n")
        f.write(f"### CQI: {sample}\n")
        f.write("##########################\n\n")

def calculate_metrics(
    folder: Path, cqi: Path, sample: Path,
    out: Path, log: logging.Logger,
    bed: Path, raw_cqi_n: int, raw_sam_n: int,
    report: Path
):
    stats = run_vcftoolz(cqi, sample, folder / "vcftoolz", log)

    exp = stats["expected"]
    fnd = stats["found"]
    pos = stats["positive"]
    mis = stats["missing"]
    noi = stats["noise"]

    tot = 0
    if bed.is_file() and bed.stat().st_size > 0:
        r = run(f"awk -F'\\t' '{{S+=$3-$2}}END{{print S}}' {bed}", log)
        tot = int(r.stdout.strip() or 0)

    lines = [
        f"Number of variants VCF REF: {raw_cqi_n}",
        f"Number of variants VCF INPUT: {raw_sam_n}",
        "",
    ]
    if bed.name:
        lines += [
            f"Number of variants after filtering on {bed.name} VCF REF: {exp}",
            f"Number of variants after filtering on {bed.name} VCF INPUT: {fnd}",
            "",
        ]
    lines += [
        "#################",
        "###  Metrics",
        "#################",
        "",
        f"# TYPES:       {folder.name}",
        f"# EXPECTED:    {exp}",
        "#              ALL ",
        f"# FOUND:       {fnd}",
        f"# POSITIVE:    {pos}",
        f"# MISSING:     {mis}",
        f"# NOISE:       {noi}",
        f"# SENSITIVITY: {pct(pos, exp)}%",
        f"# PPV:         {pct(pos, fnd)}%     Positive Predictive Value",
        f"# SPECIFICITY: {pct(tot - exp, tot - exp + noi)}%",
        "#",
    ]

    with open(out, "a") as f:
        f.write("\n".join(lines) + "\n")

    # Append to master report
    with open(report, "a") as rep, open(out) as src:
        rep.write(src.read())

    return stats["missing_vcf"]

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="CQI - Quality control analysis",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("-r", "--run",      required=True,  help="RUN directory")
    p.add_argument("-g", "--genes",    required=True,  help="BED file(s), comma-separated")
    p.add_argument("-a", "--archives", default="",     help="ARCHIVES option")
    p.add_argument("-j", "--json",     default="/databases/CQI/latest/REF.json", help="JSON config")
    p.add_argument("-o", "--genome",   required=True,  help="Reference genome FASTA")
    p.add_argument("-n", "--release",  action="store_true", help="Print release notes")
    return p.parse_args()

def main():
    print("#######################################")
    print(f"# CQI [2.0-{DATEFILE}]")
    print("# Quality control analysis")
    print("# Jean-Baptiste Lamouche @ HUS © GNU-GPL")
    print("#######################################")

    args = parse_args()
    if args.release:
        print(RELEASE_NOTES)
        sys.exit(0)

    run_dir  = Path(args.run)
    genome   = Path(args.genome)
    json_cfg = Path(args.json)

    # Temporary LOG before we have a sample-level log
    logging.basicConfig(level=logging.ERROR)
    base_log = logging.getLogger("CQI")

    if not json_cfg.is_file():
        base_log.error("No VCF JSON file")
        sys.exit(2)

    if (run_dir / "CQIComplete.txt").is_file():
        base_log.info("CQI already completed.")
        sys.exit(0)

    with open(json_cfg) as f:
        cfg = json.load(f)
    cqi_index = {entry["name"]: entry["VCF"] for entry in cfg.get("CQI", [])}

    for sample_dir in sorted(run_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        cqi_sample = sample_dir.name
        tag_file   = sample_dir / f"{cqi_sample}.tag"

        if not tag_file.is_file() or "CQI" not in tag_file.read_text():
            continue

        tag_content = tag_file.read_text()
        cqi_block   = next(
            (line for line in tag_content.splitlines() if "CQI" in line), ""
        )
        tag_str  = cqi_block.split("CQI#", 1)[-1].split("!")[0]
        full_tags = [t for t in tag_str.split("#") if t]

        for tag in full_tags:
            jfile_raw = cqi_index.get(tag)
            if not jfile_raw:
                base_log.error(f"{tag} not in list EXIT")
                sys.exit(2)

            jfile = Path(jfile_raw.lstrip("/STARK"))
            res   = next(sample_dir.rglob(f"{cqi_sample}.final.vcf.gz"), None)

            cqi_dir = run_dir / cqi_sample / "CQI" / tag
            cqi_dir.mkdir(parents=True, exist_ok=True)

            log_path    = cqi_dir / f"{cqi_sample}.analysis.{DATEFILE}.report.log"
            report_path = cqi_dir / f"{cqi_sample}.analysis.{DATEFILE}.report.tsv"
            log_path.touch(); report_path.touch()
            log = setup_logging(log_path)

            if not res or not res.is_file():
                log.error("INPUT VCF not found"); continue
            if not jfile.is_file():
                log.error(f"REF VCF not found: {jfile}"); continue

            log.info(f"SAMPLE={cqi_sample} TAG={tag}")

            genes = prepare_intervals(args.genes, cqi_dir, log)
            log.info(f"Filter interval on: {genes}")

            # Pre-count raw variants once
            raw_cqi_n = count_variants(jfile, log)
            raw_sam_n = count_variants(res, log)

            # Parallel BCF preparation
            bcf_cqi = cqi_dir / ".tmp_CQI.bcf"
            bcf_sam = cqi_dir / ".tmp_SAM.bcf"

            from concurrent.futures import ThreadPoolExecutor
            with ThreadPoolExecutor(max_workers=2) as ex:
                f1 = ex.submit(prepare_bcf, jfile, genes, genome, log, bcf_cqi)
                f2 = ex.submit(prepare_bcf, res,   genes, genome, log, bcf_sam)
                f1.result(); f2.result()

            # Split both BCFs into SNV + INDEL in one shared loop
            snv_dir   = cqi_dir / "SNV";   snv_dir.mkdir(exist_ok=True)
            indel_dir = cqi_dir / "INDEL"; indel_dir.mkdir(exist_ok=True)

            for prefix, bcf in [("CQI_VCF", bcf_cqi), ("SAMPLE_VCF", bcf_sam)]:
                split_bcf(bcf, prefix, snv_dir, indel_dir, log)
            bcf_cqi.unlink(missing_ok=True)
            bcf_sam.unlink(missing_ok=True)

            # Per-type analysis
            for type_name, folder in [("SNV", snv_dir), ("INDEL", indel_dir)]:
                output  = folder / f"{cqi_sample}.metrics_{DATEFILE}.tsv"
                cqi_pro = folder / "CQI_VCF.vcf.gz"
                sam_pro = folder / "SAMPLE_VCF.vcf.gz"

                write_report_header(output, run_dir, cqi_sample)

                if type_name == "INDEL":
                    log.info("Normalization INDEL...")
                    sam_pro = norm_and_index(sam_pro, folder / "SAMPLE_norm.vcf.gz", genome, log)
                    cqi_pro = norm_and_index(cqi_pro, folder / "CQI_norm.vcf.gz",   genome, log)

                missing_vcf = calculate_metrics(
                    folder, cqi_pro, sam_pro, output, log,
                    genes, raw_cqi_n, raw_sam_n, report_path
                )

                if type_name == "SNV":
                    compress_and_index(
                        missing_vcf,
                        cqi_dir / f"{cqi_sample}.SNV.missing.{DATEFILE}.vcf.gz",
                        log
                    )

                # Clean up type subfolder
                import shutil
                shutil.rmtree(folder, ignore_errors=True)

if __name__ == "__main__":
    main()