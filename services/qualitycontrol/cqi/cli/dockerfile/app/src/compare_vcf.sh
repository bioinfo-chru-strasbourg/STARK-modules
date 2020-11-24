#!/bin/bash
#################################
##
## VCF comparison
##
#################################

SCRIPT_NAME="STARK_comparisonVCF"
SCRIPT_DESCRIPTION="STARK Comparison 2 VCF"
SCRIPT_RELEASE="0.9.18"
SCRIPT_DATE="25/05/2020"
SCRIPT_AUTHOR="Jean-Baptiste Lamouche"
SCRIPT_COPYRIGHT="IRC"
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
        echo "#USAGE: $(basename $0) --run =<RUN>|--cqi=<CQI_SAMPLE> [options...]";
        echo "# -c/--cqi              CQI option";
		#echo "# USAGE: $(basename $0) --vcf=<VCF>|--ref=<VCF_REF> [options...]";
        #echo "# -f/--vcf            VCF to compare with the referecen VCF";
        #echo "# -r/--ref            Reference VCF";
        echo "# -r/--run              RUN option";
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
ARGS=$(getopt -o "r:c:vdnh" --long "run:,cqi:,verbose,debug,release,help" -- "$@" 2> /dev/null)
if [ $? -ne 0 ]; then
        :
fi;
PARAM=$@


eval set -- "$ARGS"
while true
do
        case "$1" in
                -r|--run)
                        RUN=$2
                        shift 2
                        ;;
                -c|--cqi)
                        CQI_SAMPLE="$2"
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



#DATE=$(date +"%Y%m%d")
#OUTPUT=$DATE"_CQI".txt
#RESULTS="run_folder"
DATA=/home1/BAS/lamouchj/CQI
#"CORIELL,HORIZON,ACCROMETRIX"
#CQI_SAMPLE="CORIELL,HORIZON,ACROMETRIX"
CQI=$DATA/$RUN/CQI
#MERGE_VCF=$CQI/sample_vcf_merge

#Loop pour merge tous les sample vcf si ils représente le CQI 
#
#if [ -f $MERGE_VCF ]; then
#	echo "File is already created"
#else
#	#In results folder not the same PATH
#	for SAMPLE_VCF in $DATA/$RUN/*/*.final.vcf.gz; do
#	CQI_VCF=$DATA/$RUN/$CQI_SAMPLE/$CQI_SAMPLE.reports/*.final.vcf.gz
#		if [ $(basename $SAMPLE_VCF .final.vcf.gz) == $CQI_SAMPLE ] || [ $SAMPLE_VCF == "*-*" ]; then 
#			echo "CQI" 
#		continue
#		fi;
#		echo $SAMPLE_VCF
#		#echo "basenamevar=$(basename $SAMPLE_VCF)"
#		echo $SAMPLE_VCF >> $MERGE_VCF
#	done;
#fi;
#bcftools merge --file-list $MERGE_VCF -Oz -o $CQI/SAMPLE_VCF.vcf.gz

#if VCF_COMPARE.complete.txt is present redo only rtg-tools
if [ ! -f $CQI/VCF_COMPARE.complete.txt ]; then
	mkdir -p $CQI
	#echo "$CQI/SNV.complete.txt"
	
	#Juste la à modif en fait TODO, SAMPLE_VCF is the CQI in STARKrun
	#SAMPLE_VCF=$DATA/$RUN/$CQI_SAMPLE/$CQI_SAMPLE.final.vcf
	SAMPLE_VCF_OR=/home1/BAS/lamouchj/CQI/$RUN/correctif_vcf_header/$CQI_SAMPLE.final.vcf
	
	#PATH results again, CQI ref, variants known
	CQI_VCF_OR=/home1/BAS/lamouchj/Coriell/TeCoriell_all_bedfilter.vcf
	
	
	ARRAY=($CQI_VCF_OR $SAMPLE_VCF_OR)


	#SPLITTING SNV/INDELS
	for VCF in "${ARRAY[@]}"; do
		mkdir -p $CQI/SNV $CQI/INDEL
		vcftools --vcf $VCF --remove-indels --recode --recode-INFO-all --out $CQI/SNV/$(basename -s .recode.vcf $VCF)
		mv $CQI/SNV/*.recode.vcf $CQI/SNV/$(basename -s .recode.vcf $VCF)
		vcftools --vcf $VCF --keep-only-indels --recode --recode-INFO-all --out $CQI/INDEL/$(basename -s .recode.vcf $VCF)
		mv $CQI/INDEL/*.recode.vcf $CQI/INDEL/$(basename -s .recode.vcf $VCF)
	done;

	#SNV / INDELS folder
	SNV_FOLDER=$CQI/SNV
	INDEL_FOLDER=$CQI/INDEL
	ARRAY=($SNV_FOLDER $INDEL_FOLDER)
	
	#Do analysis for SNV then for INDELS
	for FOLDER in "${ARRAY[@]}"; do 
		RES_VCF_SAMPLE=$FOLDER/$RUN"_VCF_SAMPLE".txt
		RES_VCF_CQI=$FOLDER/$RUN"_VCF_CQI".txt
		OUTPUT=$FOLDER/VCF_comparison.txt

		
		CQI_VCF=$FOLDER/$(basename -- $CQI_VCF_OR)
		SAMPLE_VCF=$FOLDER/$(basename -- $SAMPLE_VCF_OR)
	
		echo "SAMPLE_VCF=$SAMPLE_VCF"
		echo "CQI_VCF=$CQI_VCF"
		echo "OUTPUT=$OUTPUT"
		
		echo "##########################" >> $OUTPUT
		echo "### CQI vcf comparison " >> $OUTPUT
		echo "### RUN: $RUN" >> $OUTPUT
		echo "### CQI: $CQI_SAMPLE" >> $OUTPUT
		echo "##########################" >> $OUTPUT
		printf "\n" >> $OUTPUT
		
		#Sort vcf
		vcf-sort -c $SAMPLE_VCF > $FOLDER/SAMPLE_VCF_sorted.vcf
		vcf-sort -c $CQI_VCF > $FOLDER/CQI_VCF_sorted.vcf
		bgzip -c $FOLDER/SAMPLE_VCF_sorted.vcf > $FOLDER/SAMPLE_VCF_sorted.vcf.gz
		tabix -f $FOLDER/SAMPLE_VCF_sorted.vcf.gz
		bgzip -c $FOLDER/CQI_VCF_sorted.vcf > $FOLDER/CQI_VCF_sorted.vcf.gz
		tabix -f $FOLDER/CQI_VCF_sorted.vcf.gz
		
		#variables
		SAMPLE_VCF=$FOLDER/SAMPLE_VCF_sorted.vcf.gz
		CQI_VCF=$FOLDER/CQI_VCF_sorted.vcf.gz
		
		#total number of variants
		echo "Total number of variants SAMPLE STARK: $(zcat $SAMPLE_VCF | grep -v "^#" | wc -l)" >> $OUTPUT
		echo "Total number of variants $CQI_SAMPLE: $(zcat $CQI_VCF | grep -v "^#" | wc -l)" >> $OUTPUT
		printf "\n"
		
		# number of variants after filtering
		#echo "Total number of variants after filtering (unnormalised variants): $(bcftools view -v snps $CQI/SAMPLE_VCF_sorted.filter.vcf.gz | grep -v "^#" | wc -l)" >> $RES_VCF_SAMPLE
		#echo "Total number of variants after filtering (unnormalised variants): $(bcftools view -v snps $CQI/CQI_VCF_sorted.filter.vcf.gz | grep -v "^#" | wc -l)" >> $RES_VCF_CQI
		
		printf "\n" >> $OUTPUT
		echo "#################" >> $OUTPUT 
		echo "### VCF-compare" >> $OUTPUT
		echo "#################" >> $OUTPUT
		printf "\n" >> $OUTPUT
		
		vcf-compare $CQI_VCF $SAMPLE_VCF >> $OUTPUT
		
		bcftools isec  $CQI_VCF $SAMPLE_VCF -p $FOLDER/bcftools_isec
		
		gunzip -c $SAMPLE_VCF > $FOLDER/final_SAMPLE_VCF.vcf
		gunzip -c $CQI_VCF > $FOLDER/final_CQI_VCF.vcf 
		
		# calculate Jaccard index ( Area of overlap / area of union )
		printf "\n" >> $OUTPUT
		echo "###################" >> $OUTPUT 
		echo "### Jaccard index" >> $OUTPUT
		echo "###################" >> $OUTPUT
		printf "\n" >> $OUTPUT
		echo "$(bedtools jaccard -a $CQI_VCF -b $SAMPLE_VCF)" >> $OUTPUT
		
		# I use the -u parameter, which reports the variant in -a if an overlap was found
		echo "Variants in common: $(bedtools intersect -u -a $FOLDER/final_CQI_VCF.vcf -b $FOLDER/final_SAMPLE_VCF.vcf | wc -l)" >> $OUTPUT
		
		#vcf file 
		SnpSift concordance -v $FOLDER/final_CQI_VCF.vcf $FOLDER/final_SAMPLE_VCF.vcf > $FOLDER/snp_concordance.txt
		mv /home1/BAS/lamouchj/CQI/concordance_final_CQI_VCF_final_SAMPLE_VCF.by_sample.txt /home1/BAS/lamouchj/CQI/concordance_final_CQI_VCF_final_SAMPLE_VCF.summary.txt $FOLDER
		
		#rm -f $CQI*.sorted.vcf*
		touch $CQI/VCF_COMPARE.complete.txt
		
	done; #SNV/ INDELS
else 
	echo "#VCF_COMPARE already OK"
fi;


#############################################
################# RTG-TOOLS
################# 

function RTG () {

GENOME=/home1/TOOLS/genomes/hg19/hg19.fa
#SAMPLE_VCF=$DATA/$RUN/$CQI_SAMPLE/$CQI_SAMPLE.final.vcf
SAMPLE_VCF=SAMPLE_VCF_sorted.vcf.gz
CQI_VCF=CQI_VCF_sorted.vcf.gz
ARRAY_VCF=($CQI_VCF $SAMPLE_VCF)

#Normalize left trimming
#bcftools norm -f $HG -O z > $INDEL/$CQI_SAMPLE.vcf.gz

#SNV / INDELS folder
SNV_FOLDER=$CQI/SNV
INDEL_FOLDER=$CQI/INDEL
ARRAY_FOLDER=($SNV_FOLDER $INDEL_FOLDER)



#Create .sdf needed for analysis
if [ ! -d $CQI/hg19.sdf ]; then
	rtg format -o $CQI/hg19.sdf $GENOME
else 
	echo "#GENOME SDF already OK"
fi;


#check if analysis was performed before
if [ ! -f $CQI/RTG.complete.txt ]; then
	#Do analysis for SNV then for INDELS
	for FOLDER in "${ARRAY_FOLDER[@]}"; do
		mkdir -p $FOLDER/RTG-TOOLS
		#for CQIref and CQISTARK18	
		for VCF in "${ARRAY_VCF[@]}"; do 
			gunzip -c $FOLDER/$VCF > $FOLDER/RTG-TOOLS/${VCF%.*} 
			rtg vcfstats $FOLDER/RTG-TOOLS/${VCF%.*} > $FOLDER/RTG-TOOLS/${VCF%.*}.stats
		done;
		rtg vcfeval -b $FOLDER/$CQI_VCF -c $FOLDER/$SAMPLE_VCF -t $CQI/hg19.sdf -o $FOLDER/RTG-TOOLS/vcfeval
	done;
	touch $CQI/RTG.complete.txt
else
	echo "Test_CQI already OK"
fi;
}

#remove intermediate vcf
rm -f $CQI/*/CQI_VCF*.vcf $CQI/*/SAMPLE_VCF*.vcf $CQI/*/TeCoriell*.vcf $CQI/*/final*

#TODO,dossier CQI dans l'échantillon CQI du run, rm files to gain some place, check le docker-compose, check all values 