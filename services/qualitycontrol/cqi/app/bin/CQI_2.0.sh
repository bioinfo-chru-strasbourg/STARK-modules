#!/bin/bash
# shellcheck disable=SC2086
SCRIPT_NAME="CQI"
SCRIPT_DESCRIPTION="Quality control analysis "
SCRIPT_RELEASE="1.4"
SCRIPT_DATE="2026"
SCRIPT_AUTHOR="Jean-Baptiste Lamouche"
SCRIPT_COPYRIGHT="HUS"
SCRIPT_LICENCE="GNU-GPL"

RELEASE_NOTES="# 0.9.18: Script creation\n# 1.0: Full refactor\n# 2.0\n"
DATEFILE="$(date +'%Y%m%d-%H%M%S')"

function header() {
    echo "#######################################"
    echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]"
    echo "# $SCRIPT_DESCRIPTION"
    echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © $SCRIPT_LICENCE"
    echo "#######################################"
}
function release() { echo -e "# RELEASE NOTES:\n$RELEASE_NOTES"; }
function usage() {
    echo "#USAGE: $(basename "$0") --run=<RUN> --genes=<BED> --json=<JSON> --genome=<GENOME> [options]"
    echo "# -r/--run -g/--genes -a/--archives -j/--json -o/--genome -n/--release -h/--help"
}

function log_info()  { echo "[#INFO] $*" | tee -a "$LOG"; }
function log_error() { echo "[ERROR] $*" | tee -a "$LOG" >&2; }

# --- Helpers ---

function index_vcf() { tabix -f -p vcf "$1" 2>>"$2"; }

function compress_and_index() {
    local in="$1" out="$2" log="$3"
    if [ -s "$in" ]; then
        bgzip -c "$in" > "$out" 2>>"$log"
    else
        echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | bgzip -c > "$out"
    fi
    index_vcf "$out" "$log"
}

function try_norm() {
    local in="$1" out="$2" genome="$3" log="$4"
    local -n _ref="$5"
    bcftools norm -f "$genome" "$in" -O z -o "$out" 2>>"$log"
    if [ -s "$out" ]; then
        index_vcf "$out" "$log"
        _ref="$out"
    fi
}

function count_variants() {
    [ -s "$1" ] && bcftools view -H "$1" 2>/dev/null | wc -l || echo 0
}

function pct() {
    local num="$1" den="$2"
    [ "$den" -ne 0 ] && echo "scale=4; ($num/$den)*100" | bc | awk '{printf "%.2f\n", $0}' || echo "0.00"
}

function prepare_intervals() {
    local dir="$1" log="$2" out="$BED"
    if [[ "$BED" == *","* ]]; then
        out="$dir/CQI.$DATEFILE.intervals.genes.bed"
        cat ${BED//,/ } | sort -k1,1V -k2,2n | bedtools merge -i stdin > "$out" 2>>"$log"
        [ ! -s "$out" ] && { log_error "Gene file empty: $out"; exit 2; }
    fi
    echo "$out"
}

function prepare_bcf() {
    local in="$1" genes="$2" genome="$3" log="$4" out="$5"
    bcftools view "$in" 2>>"$log" | \
    bcftools reheader --fai "${genome}.fai" 2>>"$log" | \
    bedtools intersect -a stdin -b "$genes" -header -u 2>>"$log" | \
    bcftools sort -O b -o "$out" 2>>"$log"
}

function write_report_header() {
    local out="$1" run="$2" sample="$3"
    {
        echo "##########################"
        echo "### CQI vcf comparison "
        echo "### RUN: $(basename "$run")"
        echo "### CQI: $sample"
        echo "##########################"
        echo ""
    } > "$out"
}

function calculate_metrics() {
    local folder="$1" cqi="$2" sample="$3" out="$4" log="$5" bed="$6"
    local raw_cqi_n="$7" raw_sam_n="$8" report="$9"
    local tot=0

    local cqi_n sam_n
    cqi_n=$(count_variants "$cqi")
    sam_n=$(count_variants "$sample")

    {
        echo "Number of variants VCF REF: $raw_cqi_n"
        echo "Number of variants VCF INPUT: $raw_sam_n"
        echo ""
    } >> "$out"

    if [ -n "$bed" ]; then
        {
            echo "Number of variants after filtering on $(basename "$bed") VCF REF: $cqi_n"
            echo "Number of variants after filtering on $(basename "$bed") VCF INPUT: $sam_n"
            echo ""
        } >> "$out"
    fi

    bcftools isec -c none "$cqi" "$sample" -p "$folder/isec" 2>>"$log"

    local exp="$cqi_n" fnd="$sam_n" pos mis noi
    pos=$(count_variants "$folder/isec/0002.vcf")
    mis=$((exp - pos)); noi=$((fnd - pos))
    [ -s "$bed" ] && tot=$(awk -F'\t' '{S+=$3-$2}END{print S}' "$bed")

    {
        echo "#################"
        echo "###  Metrics"
        echo "#################"
        echo ""
        echo -e "# TYPES:       $(basename "$folder")"
        echo -e "# EXPECTED:    $exp"
        echo -e "#              ALL "
        echo -e "# FOUND:       $fnd"
        echo -e "# POSITIVE:    $pos"
        echo -e "# MISSING:     $mis"
        echo -e "# NOISE:       $noi"
        echo -e "# SENSITIVITY: $(pct "$pos" "$exp")%"
        echo -e "# PPV:         $(pct "$pos" "$fnd")%     Positive Predictive Value"
        echo -e "# SPECIFICITY: $(pct $((tot - exp)) $((tot - exp + noi)))%"
        echo -e "#"
    } >> "$out"

    cat "$out" >> "$report"
}

# --- Main ---

function main() {
    header

    local ARGS RUN BED ARCHIVES JSON GENOME
    ARGS=$(getopt -o "r:g:a:j:o:nh" --long "run:,genes:,archives:,json:,genome:,release,help" -- "$@" 2>/dev/null)
    eval set -- "$ARGS"
    while true; do
        case "$1" in
            -r|--run)      RUN="$2";      shift 2 ;;
            -g|--genes)    BED="$2";      shift 2 ;;
            -a|--archives) ARCHIVES="$2"; shift 2 ;;
            -j|--json)     JSON="$2";     shift 2 ;;
            -o|--genome)   GENOME="$2";   shift 2 ;;
            -n|--release)  release; exit 0 ;;
            -h|--help)     usage;   exit 0 ;;
            --) shift; break ;;
            *)  usage; exit 1 ;;
        esac
    done

    [ -z "$JSON" ]                && JSON="/databases/CQI/latest/REF.json"
    [ ! -f "$JSON" ]              && { log_error "No VCF JSON file"; exit 2; }
    [ -f "$RUN/CQIComplete.txt" ] && { log_info  "CQI already completed."; exit 0; }

    for SAMPLE in "$RUN"/*/; do
        [ ! -d "$SAMPLE" ] && continue
        CQI_SAMPLE=$(basename "$SAMPLE")
        grep -q "CQI" "$SAMPLE/$CQI_SAMPLE.tag" 2>/dev/null || continue

        TAG_FILE=$(grep "CQI" "$SAMPLE/$CQI_SAMPLE.tag" | sed 's/.*CQI#\(.*\)*/\1/')
        IFS='#' read -ra FULL_TAG <<< "${TAG_FILE%%!*}"

        for TAG in "${FULL_TAG[@]}"; do

            JFILE=$(jq -cr --arg TAG "$TAG" '.CQI[] | select(.name==$TAG) | .VCF' "$JSON")
            [ -z "$JFILE" ] && { log_error "$TAG not in list EXIT"; exit 2; }
            JFILE="${JFILE#/STARK}"

            RES=$(find "$SAMPLE" -name "$CQI_SAMPLE.final.vcf.gz" -print -quit)
            [ ! -f "$RES" ]   && { log_error "INPUT VCF not found";         continue; }
            [ ! -f "$JFILE" ] && { log_error "REF VCF not found: $JFILE"; continue; }

            CQI="$RUN/$CQI_SAMPLE/CQI/$TAG"; mkdir -p "$CQI"
            LOG="$CQI/${CQI_SAMPLE}.analysis.$DATEFILE.report.log";    touch "$LOG"
            REPORT="$CQI/${CQI_SAMPLE}.analysis.$DATEFILE.report.tsv"; touch "$REPORT"

            log_info "SAMPLE=$CQI_SAMPLE TAG=$TAG"

            GENES=$(prepare_intervals "$CQI" "$LOG")
            log_info "Filter interval on: $GENES"

            RAW_CQI_N=$(count_variants "$JFILE")
            RAW_SAM_N=$(count_variants "$RES")

            local bcf_cqi="$CQI/.tmp_CQI.bcf" bcf_sam="$CQI/.tmp_SAM.bcf"
            prepare_bcf "$JFILE" "$GENES" "$GENOME" "$LOG" "$bcf_cqi" &
            prepare_bcf "$RES"   "$GENES" "$GENOME" "$LOG" "$bcf_sam" &
            wait

            mkdir -p "$CQI/SNV" "$CQI/INDEL"

            for type_cfg in "SNV --types snps" "INDEL --exclude-types snps"; do
                local type="${type_cfg%% *}" flag="${type_cfg#* }"
                for pair in "CQI_VCF $bcf_cqi" "SAMPLE_VCF $bcf_sam"; do
                    local prefix="${pair%% *}" bcf="${pair#* }"
                    local out="$CQI/$type/${prefix}.vcf.gz"
                    bcftools view $flag -O z -o "$out" "$bcf" 2>>"$LOG"
                    index_vcf "$out" "$LOG"
                done
            done
            rm -f "$bcf_cqi" "$bcf_sam"

            for TYPE in SNV INDEL; do
                FOLDER="$CQI/$TYPE"
                OUTPUT="$FOLDER/${CQI_SAMPLE}.metrics_$DATEFILE.tsv"
                CQI_PRO="$FOLDER/CQI_VCF.vcf.gz"
                SAM_PRO="$FOLDER/SAMPLE_VCF.vcf.gz"

                write_report_header "$OUTPUT" "$RUN" "$CQI_SAMPLE"

                if [[ "$TYPE" == "INDEL" ]]; then
                    try_norm "$SAM_PRO" "$FOLDER/SAMPLE_norm.vcf.gz" "$GENOME" "$LOG" SAM_PRO
                    try_norm "$CQI_PRO" "$FOLDER/CQI_norm.vcf.gz"   "$GENOME" "$LOG" CQI_PRO
                fi

                calculate_metrics "$FOLDER" "$CQI_PRO" "$SAM_PRO" "$OUTPUT" "$LOG" "$GENES" \
                                   "$RAW_CQI_N" "$RAW_SAM_N" "$REPORT"

                [ "$TYPE" == "SNV" ] && compress_and_index "$FOLDER/isec/0000.vcf" \
                                        "$CQI/${CQI_SAMPLE}.SNV.missing.$DATEFILE.vcf.gz" "$LOG"
                rm -rf "$FOLDER"
            done
        done
    done
}

main "$@"