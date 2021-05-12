#!/bin/bash
#############################################
################# RTG-TOOLS
################# 

# Usage
function usage {
        echo "#USAGE: $(basename $0) --run=<RUN> | --cqi =<CQI>[options...]";
        echo "# -c/--cqi              CQI option";
		echo "# -r/--run              RUN option";

}

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "c:r:vdnh" --long "cqi:,run:,verbose,debug,release,help" -- "$@" 2> /dev/null)
if [ $? -ne 0 ]; then
        :
fi;
PARAM=$@


eval set -- "$ARGS"
while true
do
        case "$1" in
                -c|--cqi)
                        CQI="$2"
                        shift 2
                        ;;
                -r|--run)
                        RUN="$2"
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



GENOME=$GENOME
#SAMPLE_VCF=$DATA/$RUN/$CQI_SAMPLE/$CQI_SAMPLE.final.vcf
SAMPLE_VCF=final_SAMPLE_VCF.vcf.gz
CQI_VCF=final_CQI_VCF.vcf.gz
ARRAY_VCF=($CQI_VCF $SAMPLE_VCF)
#CQI=$RUN/$CQI/CQI old
DATE="$(date +'%c')"

#Normalize left trimming
#bcftools norm -f $HG -O z > $INDEL/$CQI_SAMPLE.vcf.gz

#NONREGRESSION loop over all sample which have COMPAREComplete.txt in run root
for SAMPLE in $(find $RUN -maxdepth 1 -mindepth 1 -type d); do
	if [ -d "$SAMPLE/CQI" ]; then
		CQI=$SAMPLE/CQI
		#SNV / INDELS folder
		SNV_FOLDER=$CQI/SNV
		INDEL_FOLDER=$CQI/INDEL
		ARRAY_FOLDER=($SNV_FOLDER $INDEL_FOLDER)
		
		touch $CQI/rtg_compare.log
		LOG=$CQI/rtg_compare.log
		
		#Create .sdf needed for analysis / at terms keep the sdf in a dedicated folder to skip redoing the sdf genome
		
		if [ ! -d "$CQI/hg19.sdf" ]; then
			rtg format -o $CQI/hg19.sdf $GENOME 2>>$LOG
			SDF=$CQI/hg19.sdf
		else 
			echo "#GENOME SDF already OK"
		fi;
		
		#check if analysis was performed before
		if [ ! -f $RUN/CQIComplete.txt ]; then
			#Do analysis for SNV then for INDELS
			for FOLDER in "${ARRAY_FOLDER[@]}"; do
				mkdir -p $FOLDER/RTG-TOOLS
				#for CQIref and CQISTARK18 / stats about vcf
				for VCF in "${ARRAY_VCF[@]}"; do 
					gunzip -c $FOLDER/$VCF > $FOLDER/RTG-TOOLS/${VCF%.*} 2>>$LOG
					rtg vcfstats $FOLDER/RTG-TOOLS/${VCF%.*} > $FOLDER/RTG-TOOLS/${VCF%.*}.stats 2>>$LOG
				done;
				rtg vcfeval --all-records -b $FOLDER/$CQI_VCF -c $FOLDER/$SAMPLE_VCF -t $SDF -o $FOLDER/RTG-TOOLS/vcfeval 2>>$LOG
				
				#GLOBAL file results for biologist
				#for item in $CQI/SNV/VCF_comparison.txt $CQI/SNV/RTG-TOOLS/final_CQI_VCF.vcf.stats $CQI/SNV/RTG-TOOLS/final_SAMPLE_VCF.vcf.stats $CQI/SNV/RTG-TOOLS/vcfeval/summary.txt; do (cat "${item}"; echo) >> $CQI/SNV_report.txt; done
				mkdir -p $CQI/RES/$(basename -- $FOLDER)
				for item in $FOLDER/VCF_comparison.txt $FOLDER/RTG-TOOLS/final_CQI_VCF.vcf.stats $FOLDER/RTG-TOOLS/final_SAMPLE_VCF.vcf.stats $FOLDER/RTG-TOOLS/vcfeval/summary.txt; do 
					(cat "${item}"; echo) >> $CQI/RES/$(basename -- $FOLDER)/$(basename -- $FOLDER)_report.txt; 
				done
				sed -i '0,/.*CQI_VCF.*/s//####   VCF REF   ####\n/' $CQI/RES/$(basename -- $FOLDER)/$(basename -- $FOLDER)_report.txt 
				sed -i '0,/.*SAMPLE_VCF.*/s//####   VCF STARK   ####\n/' $CQI/RES/$(basename -- $FOLDER)/$(basename -- $FOLDER)_report.txt 
				cp $FOLDER/bcftools_isec/0000.vcf $CQI/RES/$(basename -- $FOLDER)/$(basename -- $FOLDER)_miss_no_filter.vcf
				#gunzip -c $FOLDER/RTG-TOOLS/vcfeval/fn.vcf.gz > $CQI/RES/$(basename -- $FOLDER)/$(basename -- $FOLDER)_miss.vcf
				done;
		else #analysis of this script are alread performed
			echo "RTG already OK"
		fi;		
	else #first part of CQI analysis
		echo "[#INFO] Need to do execute first part to generate intermediate" 
	fi;
done #SAMPLE in RUN
touch $RUN/CQIComplete.txt && echo -e "$DATE" >> $RUN/CQIComplete.txt
rm -f $RUN/COMPAREComplete.txt
