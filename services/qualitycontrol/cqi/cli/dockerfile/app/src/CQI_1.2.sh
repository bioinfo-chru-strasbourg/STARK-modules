#!/bin/bash
#################################
##
## VCF comparison
##
#################################

SCRIPT_NAME="STARK comparisonVCF"
SCRIPT_DESCRIPTION="Quality control for STARK analysis "
SCRIPT_RELEASE="0.9.18"
SCRIPT_DATE="25/06/2020"
SCRIPT_AUTHOR="Jean-Baptiste Lamouche"
SCRIPT_COPYRIGHT="HUS"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.18 date: Script creation\n";

# Header
function header () {
        echo "#######################################";
        echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
        echo "# $SCRIPT_DESCRIPTION ";
        echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © $SCRIPT_LICENCE";
        echo "#######################################";
}

# Release
function release () {
        echo "# RELEASE NOTES:";
        echo -e $RELEASE_NOTES
}

# Usage
function usage {
        echo "#USAGE: $(basename $0) --run =<RUN>|--cqi=<CQI_SAMPLE>|--genes=<GENES>|--vaf=<VAF,FLOAT,:ALLELE>|--cqigz=<CQI_VCF>|--set=<SET>|--json=<JSON>|--genome=<GENOME> [options...]";
		echo "# -c/--cqi              CQI option";
        echo "# -r/--run              RUN option";
        echo "# -g/--genes            GENES option";
		echo "# -a/--vaf              VAF option";
		echo "# -q --cqigz            CQI_VCF option";
		echo "# -s/--set              SETS option";
		echo "# -j/--json             JSON option";
        echo "# -o/--genome           GENOME option";
		echo "# -v/--verbose          VERBOSE option";
        echo "# -d/--debug            DEBUG option";
        echo "# -n/--release          RELEASE option";
        echo "# -h/--help             HELP option";
        echo "#";
}

# header
header;

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "r:c:g:a:q:s:j:o:vdnh" --long "run:,cqi:,genes:,vaf:,cqigz:,set:,json:,genome:,verbose,debug,release,help" -- "$@" 2> /dev/null)
if [ $? -ne 0 ]; then
        :
fi;
PARAM=$@


eval set -- "$ARGS"
while true
do
        case "$1" in
                -r|--run)
                        RUN="$2"
                        shift 2
                        ;;
                -c|--cqi)
                        CQI_SAMPLE="$2"
                        shift 2
                        ;;
                -g|--genes)
                        GENES="$2"
                        shift 2
                        ;;
                -a|--vaf)
                        VAF="$2"
                        shift 2
                        ;;
                -q|--cqigz)
                        CQI_VCF="$2"
                        shift 2
                        ;;
				-s|--set)
                        SETS="$2"
                        shift 2
						;;
				-j|--json)
                        JSON="$2"
                        shift 2
						;;
				-o|--genome)
                        GENOME="$2"
                        shift 2
						;;
                -v|--verbose)
                        VERBOSE=1
                        shift 1
                        ;;
                -d|--debug)
                        VERBOSE=1
                        DEBUG=1
                        shift 1
                        ;;
                -n|--release)
                        release;
                        exit 0
                        ;;
                -h|--help)
                        usage
                        exit 0
                        ;;
                --) shift
                        break
                        ;;
                *)      echo "# Option $1 is not recognized. " "Use -h or --help to display the help." && \
                        exit 1
                        ;;
        esac
done

#source variables

#NONREG
#NONREG=false 

echo -e "\n"
echo "Parameters provided: $PARAM"
echo "CQI_SAMPLE=$CQI_SAMPLE"
echo "RUN=$RUN"
echo "GENES=$GENES"
echo "VAF=$VAF"
echo "CQI_VCF=$CQI_VCF"
echo "SETS=$SETS"
echo "JSON=$JSON"
echo "GENOME=$GENOME"
echo -e "\n"


#date=$(date '+%Y-%m-%d')
DATE="$(date +'%c')"


#FOLDER run $DATA only for BBS/RP cuz of bad vcf header 
DATA=/home1/L_PROD/NGS/BAS/lamouchj/input
JSON=/app/VCFDB/current/REF.json

#NON regression SETS abosulte path of SETS ref
#SETS=/home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB/SETS/CPSGEN


#"CORIELL,HORIZON,ACCROMETRIX"
#CQI=$DATA/$RUN/$CQI_SAMPLE/CQI

#CQI=$RUN/$CQI_SAMPLE/CQI  avant
#echo "CQI=$CQI"

#NONREGRESSION
REG=false

#SAMPLE_VCF=($($CQI_SAMPLE | tr -d '[],'))
#CQI_VCF=($($CQI_VCF | tr -d '[],'))
#echo " CQI_VCF_ALL=${CQI_VCF[@]}"
#echo " SAMPLE_VCF_ALL=${CQI_SAMPLE[@]}"

#For non regression $RUN will represent SET folder FULL PATH
if [ ! -f $RUN/list_vcf_ref.txt ] && [ "$REG" = true ]; then 
	find $SETS -maxdepth 1 -name '*.vcf' > $RUN/list_vcf_ref.txt
	#find $SETS -maxdepth 1 -mindepth 1 -type d > $RUN/list_vcf_ref.txt
fi;

######REF JSON
if [ ! -f $JSON ]; then 
	echo "[ERROR] No VCF JSON file"
	exit 2
fi;

#According to listener features a text file containing all VCF of reference is located in run's root to perform this analysis 

#touch $RUN/CQIRunning.txt && echo -e "$DATE" >> $RUN/CQIRunning.txt

#if CQIComplete here do not process files
if [ ! -f $RUN/CQIComplete.txt ]; then
	for SAMPLE in $(find $RUN -mindepth 1 -maxdepth 1 -type d); do
		CQI_SAMPLE=$(basename "$SAMPLE")
		if grep -q "CQI" $SAMPLE/$CQI_SAMPLE.tag 2>/dev/null; then
			echo "[#INFO] SAMPLE as CQI: $CQI_SAMPLE"
			TAG_FILE=$(grep "CQI" $SAMPLE/$CQI_SAMPLE.tag)
			IFS='!' read -ra FULL_TAG <<< "$(echo $TAG_FILE | cut -d '#' -f2)"
			for TAG in "${FULL_TAG[@]}"; do
				echo "[#INFO] TAG=$TAG"
				# NAME OF SAMPLE = $(basename -s .vcf ${LINE})
				#SAMPLE_VCF is the CQI in STARKrun
				#SAMPLE_VCF=$DATA/$RUN/$CQI_SAMPLE/$CQI_SAMPLE.final.vcf
				
				#Standard RUN CQI
				if [ "$REG" = false ]; then
					#TOOL to parse JSON in bash    > JQ
					if cat $JSON | jq -cr '.CQI | .[] | .name' ; then
						head $JSON
						echo "[#INFO] CQI available in JSON"
						cat $JSON | jq -cr '.CQI | .[] | .name'
						if [ $? -eq 0 ]; then
							echo "[#INFO] Parsing JQ Extract VCF"
							JFILE=$(cat $JSON | jq -cr --arg TAG "$TAG" '.CQI | . [] | select( .name==$TAG ) | .VCF')
							echo "JFILE=$JFILE"
						else 
							echo "[ERROR] $TAG not in list EXIT"
							exit 2
						fi;
					else
						echo "[ERROR] No $JSON file or JQ issues"
						exit 2
					fi;
					
					#Modif 1er Octobre
					#CQI=$RUN/$CQI_SAMPLE/CQI_$TAG
					CQI=$RUN/$CQI_SAMPLE/CQI/$TAG
					mkdir -p $CQI
					chmod 777 $CQI
					RES=$(find $RUN/$CQI_SAMPLE -name $CQI_SAMPLE.final.vcf.gz -print -quit)
					
					gunzip -c $RES > $CQI/SAMPLE_VCF.vcf 2>/dev/null
					SAMPLE_VCF=$CQI/SAMPLE_VCF.vcf
					
					gunzip -c $JFILE > $CQI/CQI_VCF.vcf 2>/dev/null
					CQI_VCF=$CQI/CQI_VCF.vcf
					
					
					#for RP/BBS ONLY tests (wrong vcf header comming from STARK17)
					#if [ $CQI_SAMPLE = "TeCoriell" ]; then
					#	CQI=$RUN/$CQI_SAMPLE/CQI_$TAG
					#	mkdir -p $CQI
					#	RES=$DATA/CQI_REF/correctif_vcf_header/$CQI_SAMPLE.final.vcf.gz
					#	gunzip -c $RES > $CQI/SAMPLE_VCF.vcf 2>/dev/null
					#	SAMPLE_VCF=$CQI/SAMPLE_VCF.vcf
					#	
					#	gunzip -c $JFILE > $CQI/CQI_VCF.vcf 2>/dev/null
					#	CQI_VCF=$CQI/CQI_VCF.vcf
					#	
					#else
					#	CQI=$RUN/$CQI_SAMPLE/CQI_$TAG
					#	mkdir -p $CQI
					#	RES=$(find $RUN/$CQI_SAMPLE -name $CQI_SAMPLE.final.vcf.gz -print -quit)
					#	gunzip -c $RES > $CQI/SAMPLE_VCF.vcf 2>/dev/null
					#	SAMPLE_VCF=$CQI/SAMPLE_VCF.vcf
					#	
					#	gunzip -c $JFILE > $CQI/CQI_VCF.vcf 2>/dev/null
					#	CQI_VCF=$CQI/CQI_VCF.vcf
					#fi;
					
					
					
				#only case in SETS CQI_VCF are simple vcf (NON REG)
				else
					#Try to put here a while tru instruct to perform NON REG
					echo "[#INFO] NON REGRESSION"
					RES=$(find $RUN/$(basename -s .vcf ${LINE}) -name $(basename -s .vcf ${LINE}).final.vcf)
					#If no vcf final in STARK continue
					if [ ! -z "${RES}" ]; then 
						CQI=$RUN/$(basename -s .vcf ${LINE})/CQI
						mkdir -p $CQI
						cp $RES $CQI/SAMPLE_VCF.vcf
						SAMPLE_VCF=$CQI/SAMPLE_VCF.vcf
						cp ${LINE} $CQI/CQI_VCF.vcf 
						CQI_VCF=$CQI/CQI_VCF.vcf
					else
						echo "[#INFO] No vcf correctly generated for $(basename -s .vcf ${LINE}).final.vcf"
						continue 
					fi;
				fi;
				
				#If no bed (python issues or whatever) find in run
				if [ "$REG" = true ]; then 
					#GENES=$(find $RUN/$(basename -s .vcf ${LINE}) -name $(basename -s .vcf ${LINE}).bed)
					echo "GENES=$GENES"
				else 
					if [[ -n $(find $RUN/$CQI_SAMPLE -name $CQI_SAMPLE.*.manifest.*bed.genes -print -quit) ]]; then 
						GENES=$(find $RUN/$CQI_SAMPLE -name $CQI_SAMPLE.*.manifest.*bed.genes -print -quit)
					else 
						GENES=$(find $RUN/$CQI_SAMPLE -name $CQI_SAMPLE.*.genes -print -quit)
					fi;
					echo "[#INFO Filter interval on: $(basename $GENES)]"
				fi;
				
				
				touch $CQI/CQI.log
				LOG=$CQI/CQI.log
				
				#NONREGRESSION
				#SAMPLE_VCF=($($CQI_SAMPLE | tr -d '[],'))
				#first item ${SAMPLE_VCF[0]} etc     all items ${OUTPUT[@]}
				
				#Variable set in Listener, depending on CQI sample 
				
				echo "SAMPLE_VCF=$SAMPLE_VCF"
				echo "CQI_VCF=$CQI_VCF"
				
				#With bedtools intersect
				#No bed leeds to exit the script, Options filter on VAF, taking all variants with a float value min/max VAF conf bcftools options
				AR=($CQI_VCF $SAMPLE_VCF)
				for VCF in "${AR[@]}"; do 
					#Filter on bed (.genes) / check pour le faire avec le .genes 
					if [ ! -z "${GENES}" ]; then 
						#miss some header informations
						bedtools intersect -a $VCF -b $GENES -header -u > $CQI/$(basename -s .vcf $VCF)_bed.vcf 2>>$LOG
						#vcftools --vcf $VCF --bed $GENES --out $CQI/$(basename -s .vcf $VCF)_bed.vcf --recode --keep-INFO-all
						#CQI_VCF_OR=$CQI/TeCoriell_ref.vcf
						if [ ! -z "${VAF}" ]; then
							FILTER=$(echo $VAF | cut -f1 -d,)
							VALUE=$(echo $VAF | cut -f2 -d,)
							ALLELE=$(echo $VAF | cut -f3 -d,)
							bcftools view --$FILTER $VALUE$ALLELE $CQI/$(basename -s .vcf $VCF)_bed.vcf -O v -o $CQI/$(basename -s .vcf $VCF)_bed_vaf.vcf 2>>$LOG
							if [[ "$VCF" = "$CQI_VCF" ]]; then
								CQI_VCF_OR=$CQI/$(basename -s .vcf $VCF)_bed_vaf.vcf
							else
								SAMPLE_VCF_OR=$CQI/$(basename -s .vcf $VCF)_bed_vaf.vcf
								#echo "Sample_or=$SAMPLE_VCF_OR"
							fi;
						else
							if [[ "$VCF" = "$CQI_VCF" ]]; then
								CQI_VCF_OR=$CQI/$(basename -s .vcf $VCF)_bed.vcf
							else 
								SAMPLE_VCF_OR=$CQI/$(basename -s .vcf $VCF)_bed.vcf
							fi;
						fi;
					else 
						echo -e "[Warning] No Bed provided!\n" 
						CQI_VCF_OR=$CQI/CQI_VCF.vcf
						SAMPLE_VCF_OR=$CQI/SAMPLE_VCF.vcf
						##Display help
						#usage;
						#exit 1
					fi;
				done;
				#echo "CQI_VCF_OR=$CQI_VCF_OR"
				#echo "SAMPLE_VCF_OR=$SAMPLE_VCF_OR"
					
				ARRAY=($CQI_VCF_OR $SAMPLE_VCF_OR)
				#SPLITTING SNV/INDELS
				for VCF in "${ARRAY[@]}"; do
					mkdir -p $CQI/SNV $CQI/INDEL
					vcf-sort -c $VCF > $VCF.sorted.vcf 2>>$LOG
					bgzip -c $VCF.sorted.vcf > $VCF.gz
					tabix -p vcf $VCF.gz
					#vcftools --vcf $VCF --remove-indels --recode --recode-INFO-all --out $CQI/SNV/$(basename -s .recode.vcf $VCF) 2>>$LOG
					#mv $CQI/SNV/*.recode.vcf $CQI/SNV/$(basename -s .recode.vcf $VCF)
					#vcftools --vcf $VCF --keep-only-indels --recode --recode-INFO-all --out $CQI/INDEL/$(basename -s .recode.vcf $VCF) 2>>$LOG
					#mv $CQI/INDEL/*.recode.vcf $CQI/INDEL/$(basename -s .recode.vcf $VCF)
					
					#Given that VCF ref and VCF produce by STARK have the same process, it doesn't matter if a few sample are consider INDEL or SNV ( it appears only in nonreg sample cuz of old vcf I guess)
					bcftools view --types snps -O v -o $CQI/SNV/$(basename -- $VCF) $VCF.gz
					bcftools view --exclude-types snps -O v -o $CQI/INDEL/$(basename -- $VCF) $VCF.gz
				done;
				
				#SNV / INDELS folder
				SNV_FOLDER=$CQI/SNV
				INDEL_FOLDER=$CQI/INDEL
				ARRAY=($SNV_FOLDER $INDEL_FOLDER)
				
				#Do analysis for SNV then for INDELS
				for FOLDER in "${ARRAY[@]}"; do 
					OUTPUT=$FOLDER/VCF_comparison.txt
					
					CQI_VCF=$FOLDER/$(basename -- $CQI_VCF_OR)
					SAMPLE_VCF=$FOLDER/$(basename -- $SAMPLE_VCF_OR)
				
					echo "SAMPLE_VCF=$SAMPLE_VCF"
					echo "CQI_VCF=$CQI_VCF"
					echo "OUTPUT=$OUTPUT"
					
					echo "##########################" >> $OUTPUT
					echo "### CQI vcf comparison " >> $OUTPUT
					echo "### RUN: $(basename -- $RUN)" >> $OUTPUT
					echo "### CQI: $CQI_SAMPLE" >> $OUTPUT
					echo "##########################" >> $OUTPUT
					printf "\n" >> $OUTPUT
					
					echo "Number of variants VCF REF: $(cat $CQI/CQI_VCF.vcf | grep -v "^#" | wc -l)" >> $OUTPUT
					echo "Number of variants VCF STARK: $(cat $CQI/SAMPLE_VCF.vcf | grep -v "^#" | wc -l)" >> $OUTPUT
					printf "\n" >> $OUTPUT
					
					if [ ! -z "${GENES}" ]; then
						echo "Number of variants after filtering on $(basename $GENES) VCF REF: $(cat $CQI_VCF_OR | grep -v "^#" | wc -l)" >> $OUTPUT
						echo "Number of variants after filtering on $(basename -- $GENES) VCF STARK: $(cat $SAMPLE_VCF_OR | grep -v "^#" | wc -l)" >> $OUTPUT
						printf "\n" >> $OUTPUT
					fi;
				
					#Sort vcf
					vcf-sort -c $SAMPLE_VCF > $FOLDER/SAMPLE_VCF_sorted.vcf 2>>$LOG
					vcf-sort -c $CQI_VCF > $FOLDER/CQI_VCF_sorted.vcf 2>>$LOG
					bgzip -c $FOLDER/SAMPLE_VCF_sorted.vcf > $FOLDER/SAMPLE_VCF_sorted.vcf.gz 2>>$LOG
					tabix -f $FOLDER/SAMPLE_VCF_sorted.vcf.gz 2>>$LOG
					bgzip -c $FOLDER/CQI_VCF_sorted.vcf > $FOLDER/CQI_VCF_sorted.vcf.gz 2>>$LOG
					tabix -f $FOLDER/CQI_VCF_sorted.vcf.gz 2>>$LOG
					
					#variables
					if [[ $(basename -- $FOLDER) == "INDEL" ]]; then
						echo "[#INFO] NORM INDEL"
						bcftools norm -f $GENOME $FOLDER/SAMPLE_VCF_sorted.vcf.gz -O z -o $FOLDER/SAMPLE_VCF_sorted_norm.vcf.gz 2>>$LOG
						tabix -f $FOLDER/SAMPLE_VCF_sorted_norm.vcf.gz 2>>$LOG
						
						bcftools norm -f $GENOME $FOLDER/CQI_VCF_sorted.vcf.gz -O z -o $FOLDER/CQI_VCF_sorted_norm.vcf.gz 2>>$LOG
						tabix -f $FOLDER/CQI_VCF_sorted_norm.vcf.gz 2>>$LOG
						
						SAMPLE_VCF=$FOLDER/SAMPLE_VCF_sorted_norm.vcf.gz
						CQI_VCF=$FOLDER/CQI_VCF_sorted_norm.vcf.gz
					else 
						SAMPLE_VCF=$FOLDER/SAMPLE_VCF_sorted.vcf.gz
						CQI_VCF=$FOLDER/CQI_VCF_sorted.vcf.gz
					fi;
					
					#total number of variants
					#echo "Number of variants CQI | $(basename -- $FOLDER): $(zcat $CQI_VCF | grep -v "^#" | wc -l)" >> $OUTPUT
					#echo "Number of variants SAMPLE | $(basename -- $FOLDER): $(zcat $SAMPLE_VCF | grep -v "^#" | wc -l)" >> $OUTPUT
					#printf "\n" >> $OUTPUT
					
					# number of variants after filtering
					#echo "Total number of variants after filtering (unnormalised variants): $(bcftools view -v snps $CQI/SAMPLE_VCF_sorted.filter.vcf.gz | grep -v "^#" | wc -l)" >> $RES_VCF_SAMPLE 
					#echo "Total number of variants after filtering (unnormalised variants): $(bcftools view -v snps $CQI/CQI_VCF_sorted.filter.vcf.gz | grep -v "^#" | wc -l)" >> $RES_VCF_CQI
					
					#vcf-compare $CQI_VCF $SAMPLE_VCF >> $OUTPUT 2>>$LOG
					
					bcftools isec -c none $CQI_VCF $SAMPLE_VCF -p $FOLDER/bcftools_isec 2>>$LOG
					
					gunzip -c $SAMPLE_VCF > $FOLDER/final_SAMPLE_VCF.vcf 2>>$LOG
					gunzip -c $CQI_VCF > $FOLDER/final_CQI_VCF.vcf  2>>$LOG
					
					# calculate Jaccard index ( Area of overlap / area of union )
					#echo "###################" >> $OUTPUT 
					#echo "### Jaccard index" >> $OUTPUT
					#echo "###################" >> $OUTPUT
					#printf "\n" >> $OUTPUT
					#
					##bedtools sort -chrThenSizeA -header -i $FOLDER/final_CQI_VCF.vcf > $FOLDER/final_sorted_CQI_VCF.vcf 2>>$LOG
					##bedtools sort -chrThenSizeA -header -i $FOLDER/final_SAMPLE_VCF.vcf > $FOLDER/final_sorted_SAMPLE_VCF.vcf 2>>$LOG
					#echo "$(bedtools jaccard -a $FOLDER/final_CQI_VCF.vcf -b $FOLDER/final_SAMPLE_VCF.vcf)" >> $OUTPUT 2>>$LOG
					#printf "\n" >> $OUTPUT
					
					# I use the -u parameter, which reports the variant in -a if an overlap was found
					#echo "Variants in common: $(bedtools intersect -u -a $FOLDER/final_CQI_VCF.vcf -b $FOLDER/final_SAMPLE_VCF.vcf | wc -l)" >> $OUTPUT 2>>$LOG
					#printf "\n" >> $OUTPUT
					
					#vcf file 
					#SnpSift concordance -v $FOLDER/final_CQI_VCF.vcf $FOLDER/final_SAMPLE_VCF.vcf > $FOLDER/snp_concordance.txt 2>>$LOG
					#mv $RUN/concordance* $FOLDER 2>>$LOG
					#mv /concordance* $FOLDER 2>>$LOG
					
					#### output for sample detail  expected without bed filtering
					#NB_EXPECTED=$(zcat ${LINE} | grep -v "^#" | wc -l)
					
					NB_EXPECTED=$(cat $FOLDER/final_CQI_VCF.vcf | grep -v "^#" | wc -l)
					
					NB_FOUND=$(cat $FOLDER/final_SAMPLE_VCF.vcf | grep -v "^#" | wc -l)
					#NB_FOUND=$(find $RUN/$(basename -s .vcf.gz ${LINE}) -name "$(basename -s .vcf.gz ${LINE}).final.vcf.gz" -exec zcat {} \; | grep -v "^#" | wc -l)
					NB_POSITIVE=$(cat $FOLDER/bcftools_isec/0002.vcf | grep -v "^#" | wc -l)
					NB_MISSING=$(($NB_EXPECTED-$NB_POSITIVE))
					NB_NOISE=$(($NB_FOUND-$NB_POSITIVE))
					
					## calcul of SENSITIVITY for the sample
					if [ $NB_EXPECTED -eq 0 ]; then 
						SEN=0
					else
						SEN=$(echo "scale=4; ($NB_POSITIVE/$NB_EXPECTED)*100" | bc | sed s/00$// | awk '{printf "%.2f\n", $0}')
					fi;
					
				
					# calcul of SPECIFICITY for the sample
					# ALL
					if [ $NB_FOUND -eq 0 ]; then 
						SPE=0;
					else
						SPE=$(echo "scale=4; ($NB_POSITIVE/$NB_FOUND)*100" | bc | sed s/00$// | awk '{printf "%.2f\n", $0}')
					fi;
					
					echo "#################" >> $OUTPUT 
					echo "### VCF-compare" >> $OUTPUT
					echo "#################" >> $OUTPUT
					printf "\n" >> $OUTPUT
					echo -e "# TYPES:       $(basename -- $FOLDER)" >> $OUTPUT
					echo -e "# EXPECTED:    $NB_EXPECTED" >> $OUTPUT
					echo -e "#              ALL " >> $OUTPUT
					echo -e "# FOUND:       $NB_FOUND" >> $OUTPUT
					echo -e "# POSITIVE:    $NB_POSITIVE" >> $OUTPUT
					echo -e "# MISSING:     $NB_MISSING" >> $OUTPUT
					echo -e "# NOISE:       $NB_NOISE" >> $OUTPUT
					echo -e "# SENSITIVITY: $SEN%" >> $OUTPUT
					echo -e "# SPECIFICITY: $SPE%" >> $OUTPUT
					echo -e "#" >> $OUTPUT
					
					#Files used by RTG-tools
					bgzip -c $FOLDER/final_SAMPLE_VCF.vcf > $FOLDER/final_SAMPLE_VCF.vcf.gz
					tabix -p vcf $FOLDER/final_SAMPLE_VCF.vcf.gz
								
					bgzip -c $FOLDER/final_CQI_VCF.vcf > $FOLDER/final_CQI_VCF.vcf.gz
					tabix -p vcf $FOLDER/final_CQI_VCF.vcf.gz
					
					#RES folder for biologist
					mkdir -p $CQI/RES && touch $CQI/RES/all_miss.vcf
					zcat final_CQI_VCF.vcf.gz | grep "#CHROM" > $CQI/RES/all_miss_sorted.vcf && grep -v "^#" $FOLDER/bcftools_isec/0000.vcf >> $CQI/RES/all_miss.vcf
					cp $FOLDER/VCF_comparison.txt $CQI/RES/$(basename -- $FOLDER)_comparison.txt
					
					#ALL
					SUM_NB_EXPECTED=$(($SUM_NB_EXPECTED+$NB_EXPECTED))
					SUM_NB_FOUND=$(($SUM_NB_FOUND+$NB_FOUND))
					SUM_NB_POSITIVE=$(($SUM_NB_POSITIVE+$NB_POSITIVE))
					SUM_NB_MISSING=$(($SUM_NB_MISSING+$NB_MISSING))
					SUM_NB_NOISE=$(($SUM_NB_NOISE+$NB_NOISE))
				
				done; #SNV/ INDELS
				
				#SORT Missing variants
				sort -k1,1V -k2,2n $CQI/RES/all_miss.vcf >> $CQI/RES/all_miss_sorted.vcf
			done; #IF multiple tag CQI for one sample 
			#remove intermediate vcf
			rm -f $CQI/CQI_VCF* $CQI/SAMPLE* $CQI/*.vcf.gz $CQI/*/CQI_VCF* $CQI/*/SAMPLE_VCF* $CQI/*/final*.vcf  #/home1/BAS/lamouchj/CQI/bin/concordance*.txt
		else #SAMPLE which are not CQI 
			continue
		fi;
	done; #Loop over sample in run 
else #CQIComplete.txt
	echo "[#INFO] CQI already done"
fi;


#NONREGRESSION bilan globale TODO trouble with name of GLOB
if [ "$REG" = true ]; then
	GLOB=$RUN/validation_$(echo $(basename -- $RUN) | cut -d '_' -f 3-6).txt 
	for SAMPLE in $(find $RUN -name "CQI" -type d); do 
		printf "\n" >> $GLOB
		echo "#####################################" >> $GLOB
		echo "####    SAMPLE: $(basename -- $(dirname "$SAMPLE"))" >> $GLOB
		echo "#####################################" >> $GLOB
		printf "\n" >> $GLOB
		cat $SAMPLE/SNV/VCF_comparison.txt >> $GLOB
		printf "\n" >> $GLOB
		echo "### MISSING SNV VARIANTS" >> $GLOB
		cat $SAMPLE/SNV/bcftools_isec/0000.vcf | grep -v "^#" >> $GLOB 
		printf "\n" >> $GLOB
		cat $SAMPLE/INDEL/VCF_comparison.txt >> $GLOB
		printf "\n" >> $GLOB
		echo "### MISSING INDEL VARIANTS" >> $GLOB
		cat $SAMPLE/INDEL/bcftools_isec/0000.vcf | grep -v "^#" >> $GLOB
		
	done;
		### CALCUL for APP
	if [ $SUM_NB_EXPECTED -eq 0 ]; then 
		SEN=0;
	else
		SEN=$(echo "scale=4; ($SUM_NB_POSITIVE/$SUM_NB_EXPECTED)*100" | bc | sed s/00$// | awk '{printf "%.2f\n", $0}')
	fi;
	
	# ALL
	if [ $SUM_NB_FOUND -eq 0 ]; then 
		SPE=0;
	else
		SPE=$(echo "scale=4; ($SUM_NB_POSITIVE/$SUM_NB_FOUND)*100" | bc | sed s/00$// | awk '{printf "%.2f\n", $0}');
	fi;
				
	### output for APP detail
	echo -e "\n" >> $GLOB
	echo -e "######################################" >> $GLOB
	echo -e "### VALIDATION STARK $SCRIPT_RELEASE" >> $GLOB
	echo -e "### SET $(basename -- $SETS) [$(cat $RUN/list_vcf_ref.txt | wc -l) SAMPLES]" >> $GLOB
	echo -e "### APP $(echo $(basename -- $RUN) | cut -d '_' -f6)" >> $GLOB
	echo -e "######################################" >> $GLOB
	echo -e "# EXPECTED:    $SUM_NB_EXPECTED" >> $GLOB
	echo -e "#              ALL " >> $GLOB
	echo -e "# FOUND:       $SUM_NB_FOUND" >> $GLOB
	echo -e "# POSITIVE:    $SUM_NB_POSITIVE" >> $GLOB
	echo -e "# MISSING:     $SUM_NB_MISSING" >> $GLOB
	echo -e "# NOISE:       $SUM_NB_NOISE" >> $GLOB
	echo -e "# SENSITIVITY: $SEN%" >> $GLOB
	echo -e "# SPECIFICITY: $SPE% " >> $GLOB
	echo -e "######################################" >> $GLOB
fi;

#STARK_FOLDER
touch $RUN/CQIComplete.txt && chmod 777 $RUN/CQIComplete.txt && echo -e "$DATE" >> $RUN/CQIComplete.txt
rm -f $RUN/CQIRunning.txt

### Copy to depository check with Antony
#for RES in $(find $RUN -name "CQI" -type d); do 
#	rsync -rvt $RES /home1/data/current/output/depository/$(dirname $(echo $RES | cut -d '/' -f7-))
#done;