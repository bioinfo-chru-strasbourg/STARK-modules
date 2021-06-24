#!/bin/bash
if [ -z $STARK ]; then 
	STARK=/STARK/tools/stark/current/
	
fi;

##### Tools needed during test, need to be mounted check versions in bin/env_tools.sh from STARK17
TOOL="/STARK/tools"
BCFTOOLS="$TOOL/bcftools/current/bin/bcftools"
TABIX="$TOOL/htslib/current/bin/tabix"
VCFTOOLS="$TOOLS/vcftools/0.1.14/bin"
BGZIP="$TOOL/htslib/current/bin/bgzip"
SAMTOOLS="$TOOL/samtools/current/bin/samtools"

##### VERSIONS tools STARK17 nonregression
#VCFTOOLS_VERSION=0.1.13
#BCFTOOLS_VERSION=1.8 
#BGZIP and TABIX htslib 1.8
#BGZIP_VERSION=1.8
#TABIX_VERSION=1.8

##### STARK18
STARK_docker=("/STARK/tools/stark/current/bin/STARK")

##### Folders
#DATA="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" #/home/TOOLS/data
#TODO test in personnal folder

#DATA, location of SETS_STARK
#DATA=$DATA

#VALIDATION=$VALIDATION

#ANALYSIS_PATTERN="VALIDATION_"$VALIDATION"_SET_" #$(date +"%Y%m%d")
#ANALYSIS_PATTERN=$ANALYSIS_PATTERN

#RESULTS="/home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB/validation_test_res"
#RESULTS=$RESULTS

#Location of apps in container
#APPS=$APPS

#APPS="/STARK/config/myapps"
#RESULTS=$DATA"/VALIDATION"


# Usage
function usage {
	echo "# USAGE: $(basename $0) --set=<SET_NAME>|--folder=<DATA>|--cov=<COVERAGE>|--results=<RESULTS>--output=<OUTPUT> ";
	echo "# -s/--set 		SET name ";
	echo "# -f/--folder 		SET folder path ";
	echo "# -c/--cov 		filter coverage";
	echo "# -r/--results 		output path";
	echo "# -o/--output 	output path for copy to externel server";
	echo "# -h/--help 		HELP option";
	echo "#";
}

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "s:f:c:r:o:h" --long "set:,folder:,cov:,results:,output:,help" -- "$@" 2> /dev/null)
if [ $? -ne 0 ]; then
	: 
fi;
PARAM=$@

eval set -- "$ARGS"
while true
do
	case "$1" in
		-s|--set)
			SET_INPUT="$2"
			shift 2 
			;;
		-f|--folder)
			DATA="$2"
			shift 2
			;;
		-c|--cov)
			COVERAGE="$2"
			shift 2
			;;
		-r|--results)
			RESULTS="$2"
			shift 2
			;;
		-o|--output)
			OUTPUT="$2"
			shift 2 
			;;
		-h|--help)
			usage
			exit 0
			;;
		--) shift
			break 
			;;
		*) 	echo "# Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done

echo "SET_INPUT=$SET_INPUT"
echo "SET_FOLDER=$DATA"
echo "COVERAGE=$COVERAGE"
echo "RESULTS=$RESULTS"
echo "OUTPUT=$OUTPUT"


#TEST
#DATA="/home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB"
#DATA="/STARK/data/TESTJB"
VALIDATION=$(date +"%Y%m%d") # 20171202
#VALIDATION="20210104"
ANALYSIS_PATTERN="VALIDATION_"$VALIDATION"_SET_" #$(date +"%Y%m%d")
#RESULTS="/home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB/stark_0.9.18d"
#RESULTS="/home1/L_PROD/NGS/VALIDATION/VALIDATION_STARK/stark_0.9.18.2b/VALIDATION_20200809"
#RESULTS="/STARK/data/TESTJB/output"
APPS="/STARK/config/myapps"

#source les outils VCFTOOLS et compagnie
#source /home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB/scripts/VALIDATION/.env

# # PARAMS
# SET_INPUT=$1 #"CPSGEN" #"HUSHEMATO" #"HUSTUMSOL" #"BORDEAUX" #"ROUEN"
# COVERAGE=$2

# STARK Analysis
EXEC_STARK=1 # by default set to 0 mean not to do stark analysis

## copy to externel server emclabo VALIDATION folder
# OUTPUT=/home1/L/NGS/VALIDATION/VALIDATION_STARK
COPY_SERVER=0 # by default set to 0 mean not to copy to externel server 


mkdir -p $RESULTS

# FOREACH SET
for S in $DATA/SETS_STARK/*; do

	SET_NAME=$(basename $S)
	
	if [ "$SET_NAME" == "$SET_INPUT" ] || [ "$SET_INPUT" == "" ]; then
	
		## foreach application selected in application file
		for APP in $(cat $S/application); do
			

			#echo $APP;
			### if application exist we continu if not break
			if [ ! -e $APPS/*/$APP.app ] && [ ! -e $STARK/config/apps/$APP.app ] && [ ! -e $STARK/config/apps/myapps/$APP.app ]; then
				#|| [ ! -e $APPS/DIAG/$APP.app ] || [ ! -e $APPS/ONCO/$APP.app ]; then
				echo "$APP.app not found !"
				continue
			fi;


			#source configs de chaque app
			#source $STARK/config/apps/$APP.app 
			#source $(find /home1/data/STARK/config/myapps -name '$APP".app".*')
			
			RUN=$ANALYSIS_PATTERN$SET_NAME"_APP_"$APP

			VALIDATION_OUTPUT=$RESULTS/$RUN/validation_$SET_NAME"_APP_"$APP.txt
			mkdir -p $RESULTS/$RUN

			echo "######################################" >> $VALIDATION_OUTPUT
			echo "### VALIDATION '$VALIDATION' STARK '$STARK_VERSION'" >> $VALIDATION_OUTPUT
			echo "### SET '$SET_NAME'" >> $VALIDATION_OUTPUT
			echo "### APP '$APP'" >> $VALIDATION_OUTPUT
			echo "### COVERAGE '$COVERAGE X'" >> $VALIDATION_OUTPUT
			echo "######################################" >> $VALIDATION_OUTPUT
			echo "#" >> $VALIDATION_OUTPUT
		
		
			E=$APP
	
			R1_LIST="";
			R2_LIST="";
			N_LIST="";
			D_LIST="";
			V_LIST="";
	

			nb_s=0
			NB_SAMPLE=0
			
			## if bed file exist
			[ -f $DATA/SETS_STARK/$SET_NAME/*.bed ] && BED=$(ls $DATA/SETS_STARK/$SET_NAME/*.bed | head -n1)	

			for SAMPLE_VCF in $DATA/SETS_STARK/$SET_NAME/*vcf; do
		
				nb_s=$(($nb_s+1))
		
				SAMPLE_NAME=$(basename $SAMPLE_VCF | sed s/.vcf$//gi)
				R1=$DATA/SETS_STARK/$SET_NAME/$SAMPLE_NAME"_R1.fastq.gz"
				R2=$DATA/SETS_STARK/$SET_NAME/$SAMPLE_NAME"_R2.fastq.gz"
				
				## if manifest file exist
				[ -f $DATA/SETS_STARK/$SET_NAME/*.manifest ] && MANIFEST=$(ls $DATA/SETS_STARK/$SET_NAME/*.manifest | head -n1)
		
				### created list of R1,R2, sample name to apply STARK by samples 
				if [ -e $R1 ] && [ -e $R2 ]; then
		
					NB_SAMPLE=$(($NB_SAMPLE+1))
					
					N=$SAMPLE_NAME
					B=$BED
					M=$MANIFEST
					V=$SAMPLE_VCF
					F="$RESULTS/RES/$RUN/$N/$N.reports/$N.final.vcf"

					if [ -e "$M" ]; then ### select manifest 
						D_LIST="$M,$D_LIST"
						R1_LIST="$R1,$R1_LIST"
						R2_LIST="$R2,$R2_LIST"
						N_LIST="$N,$N_LIST"
					
					elif [ -e "$B" ]; then ### select bed 
						D_LIST="$B,$D_LIST"
						R1_LIST="$R1,$R1_LIST"
						R2_LIST="$R2,$R2_LIST"
						N_LIST="$N,$N_LIST"
					else
						R1_LIST="$R1,$R1_LIST"
						R2_LIST="$R2,$R2_LIST"
						N_LIST="$N_LIST$N,"
						echo "[ERROR] No MANIFEST/BED file";
						exit 0;
					fi;
		
				fi;

					### created a list of R1 , R2, sample_name
					R1_LIST=$(echo $R1_LIST | sed s/,$//gi)
					R2_LIST=$(echo $R2_LIST | sed s/,$//gi)
					N_LIST=$(echo $N_LIST | sed s/,$//gi)
					D_LIST=$(echo $D_LIST | sed s/,$//gi)
		
			done;
		
		
	
			R1_LIST=$(echo $R1_LIST | sed s/,$//gi)
			R2_LIST=$(echo $R2_LIST | sed s/,$//gi)
			N_LIST=$(echo $N_LIST | sed s/,$//gi)
			D_LIST=$(echo $D_LIST | sed s/,$//gi)
			V_LIST=$(echo $V_LIST | sed s/,$//gi)

	
			#Care paired end reads1, reads2
			# echo "$STARK/STARK -e $E -f $R1_LIST -q $R2_LIST -s $N_LIST -b $D_LIST -r $RUN -o $RESULTS --threads=$THREADS;"
			if (($EXEC_STARK)); then
				echo "$STARK_docker --application=$E --reads=$R1_LIST --reads2=$R2_LIST --sample=$N_LIST --design=$D_LIST --analysis=$RUN --output=$RESULTS if $EXEC_STARK;"

				#older version 0.9.18b $STARK_docker -e $E -f $R1_LIST -q $R2_LIST -s $N_LIST -b $D_LIST -r $RUN -o $RESULTS --threads=2 if $EXEC_STARK;
				"${STARK_docker[@]}" --application=$E --reads=$R1_LIST --reads2=$R2_LIST --sample=$N_LIST --design=$D_LIST --analysis_name=$RUN --output=$RESULTS;


			fi;
		
			### set all variable to 0
			SUM_NB_EXPECTED=0
			SUM_NB_FOUND=0
			SUM_NB_POSITIVE=0
			SUM_NB_MISSING=0
			SUM_NB_NOISE=0
			

			NB_SAMPLE=$(ls $DATA/SETS_STARK/$SET_NAME/*.vcf | wc -l)
	
			for SAMPLE_VCF in $DATA/SETS_STARK/$SET_NAME/*.vcf; do
				echo "SAMPLE_VCF=$SAMPLE_VCF"
				#SAMPLE_NAME=$(basename $SAMPLE_VCF | sed s/.vcf$//gi | sed s/.20200421-090423//gi)
				SAMPLE_NAME=$(basename $SAMPLE_VCF | sed s/.vcf$//gi)
				echo "SAMPLE_NAME=$SAMPLE_NAME"
				
				#date mise manuellement comme j'avais sorti les résultat le jouor d'avant sinon se fait auto
				#RUN="VALIDATION_20200424_SET_HUSHEMATO_APP_HUSHEMATO.TSOMYELOID"
				
				#manque dossier results en plus du fait qu'aucun VCF n'est généré, erreur de manifest on dirait vu les logs 
				#/home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB/validation_test_res/results
				VCF_FINAL=$RESULTS/results/$RUN/$SAMPLE_NAME/$SAMPLE_NAME.reports/$SAMPLE_NAME.final.vcf
				VCF_missing_variant=$RESULTS/results/$RUN/$SAMPLE_NAME/$SAMPLE_NAME.reports/$SAMPLE_NAME.VCF_missing_variant.vcf
				VCF_SAMPLE_LOW_COV_FILTER=$RESULTS/results/$RUN/$SAMPLE_NAME/$SAMPLE_NAME.reports/$SAMPLE_NAME.VCF_SAMPLE_LOW_COV_FILTER.vcf
				cat $SAMPLE_VCF | grep '#' >> $VCF_SAMPLE_LOW_COV_FILTER								
				VCF_SAMPLE_filter_variant=$RESULTS/results/$RUN/$SAMPLE_NAME/$SAMPLE_NAME.reports/$SAMPLE_NAME.VCF_SAMPLE_filter_variant.vcf	
				VCF_SAMPLE_filter=$RESULTS/results/$RUN/$SAMPLE_NAME/$SAMPLE_NAME.reports/$SAMPLE_NAME.VCF_SAMPLE_filter.vcf		
				cat $SAMPLE_VCF | grep '#' >> $VCF_SAMPLE_filter
				VCF_SAMPLE_HIGH_COV_FILTER=$RESULTS/results/$RUN/$SAMPLE_NAME/$SAMPLE_NAME.reports/$SAMPLE_NAME.VCF_SAMPLE_HIGH_COV_FILTER.vcf		
				cat $SAMPLE_VCF | grep '#' >> $VCF_SAMPLE_HIGH_COV_FILTER

				BAM=$RESULTS/results/$RUN/$SAMPLE_NAME/$SAMPLE_NAME.bwamem.bam
				echo "BAM=$BAM"

				if [ -s $VCF_FINAL ]; then 

					# CHECK generated VCF
					#TODO mmodifs tools for the other tools
					#/STARK/tools/stark/current/toolbox/check_vcf.sh -r $SAMPLE_VCF -f $VCF_FINAL
					#ADD script comparaison de vcf
					#$STARK/toolbox/check_vcf.sh -r $SAMPLE_VCF -f $VCF_FINAL
					
					#docker run -dti --name=non-regression-test --env-file=/home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB/scripts/VALIDATION/bin/var -v /home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB:/home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB cqi:latest bash -c '/home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB/scripts/VALIDATION/bin/compare_vcf_non-regression_1.1.sh --set=$SET_NAME'
					
					#docker run -dti --name=VALIDATION_STARK_COMP --entrypoint "/bin/bash '/home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB/scripts/VALIDATION/bin/docker_VALIDATION_SETS_STARK_JB.sh --cov=${STARK_COV} --set=${SET_NAME}'" env-file=/home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB/scripts/VALIDATION/bin/var -v /home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB:/home1/L_PROD/NGS/BIO_INFO/BIO_INFO_JB cqi:latest
					
					### fichier sortie STARK18 deja indexé par tabix donc a voir pour la comparaison bcftools
					### sort and gzip SAMPLE_VCF and  VCF_FINAL
					echo "$VCFTOOLS/vcf-sort"
					$VCFTOOLS/vcf-sort -c $SAMPLE_VCF > $SAMPLE_VCF.sorted
					$BGZIP -c $SAMPLE_VCF.sorted > $SAMPLE_VCF.gz
					$TABIX -f $SAMPLE_VCF.gz
					
					$VCFTOOLS/vcf-sort -c $VCF_FINAL > $VCF_FINAL.sorted 
					$BGZIP -c $VCF_FINAL.sorted > $VCF_FINAL.gz 
					$TABIX -f $VCF_FINAL.gz
					#meme avec erreur résultat cohérent jusqu'ici
					### comparison between VCF_FINAL and SAMPLE_VCF to find all missing variations from SAMPLE_VCF
					$BCFTOOLS isec -c all -C $SAMPLE_VCF.gz $VCF_FINAL.gz >> $VCF_missing_variant


					#### to find missing variation with depth low/high than $COVERAGE
					#### to find missing variation with depth low/high than $COVERAGE
					count=0
					count_HIGH=0
					while read -r line ; do
						POS=$(echo $line | awk '{print $2}')
						CHROM=$(echo $line | awk '{ print $1}')
						DEPTH=$($SAMTOOLS depth -r $CHROM:$POS-$POS $BAM | awk '{print $3}')
						#DEPTH=$(/STARK/samtools/current/bin/samtools depth -r $CHROM:$POS-$POS $BAM | awk '{print $3}')
						### if depth of missing variation is lower than $COVERAGE we don't include for $NB_EXPECTED calculation 
						if [[ $DEPTH -lt $COVERAGE ]]; then  
							let count++
  							cat $SAMPLE_VCF | grep $CHROM | grep $POS >> $VCF_SAMPLE_LOW_COV_FILTER

  						else
  							cat $SAMPLE_VCF | grep $CHROM | grep $POS >> $VCF_SAMPLE_HIGH_COV_FILTER 
  							let count_HIGH++
  						fi;
					done < $VCF_missing_variant

					### sort and gzip VCF_FINAL_FILTER_FAILED
					$VCFTOOLS/vcf-sort -c $VCF_SAMPLE_LOW_COV_FILTER > $VCF_SAMPLE_LOW_COV_FILTER.sorted 
					$BGZIP -c $VCF_SAMPLE_LOW_COV_FILTER.sorted > $VCF_SAMPLE_LOW_COV_FILTER.gz 
					$TABIX -f $VCF_SAMPLE_LOW_COV_FILTER.gz

					### created VCF_SAMPLE_filter_variant with only variation who pass the coverage filter 
					$BCFTOOLS isec -c all -C $SAMPLE_VCF.gz $VCF_SAMPLE_LOW_COV_FILTER.gz >> $VCF_SAMPLE_filter_variant	
					while read -r line ; do
						POS=$(echo $line | awk '{print $2}')
						CHROM=$(echo $line | awk '{ print $1}')
						# echo -e "CHROM=$CHROM \t POS=$POS" 
						cat $SAMPLE_VCF | grep "$CHROM" | grep "$POS" >> $VCF_SAMPLE_filter
					done < $VCF_SAMPLE_filter_variant
					$VCFTOOLS/vcf-sort -c $VCF_SAMPLE_filter > $VCF_SAMPLE_filter.sorted 
					$BGZIP -c $VCF_SAMPLE_filter.sorted > $VCF_SAMPLE_filter.gz 
					$TABIX -f $VCF_SAMPLE_filter.gz

				
					### quality calculation
					NB_MISSING_real=$(cat $VCF_missing_variant | wc -l)
					NB_EXPECTED_real=$(grep ^# -cv $SAMPLE_VCF)
					NB_EXPECTED=$($BGZIP -dc $VCF_SAMPLE_filter.gz | grep ^# -cv)
					NB_FOUND=$($BGZIP -dc $VCF_FINAL.gz | grep ^# -cv)
					NB_MISSING=$($BCFTOOLS isec -c all -C $VCF_SAMPLE_filter.gz $VCF_FINAL.gz 2>/dev/null | wc -l) 
					NB_POSITIVE=$(($NB_EXPECTED-$NB_MISSING))
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
				
		
					#### output for sample detail
					echo -e "\n" >> $VALIDATION_OUTPUT
					echo -e "### Sample '$SAMPLE_NAME'" >> $VALIDATION_OUTPUT
					echo -e "# EXPECTED: $NB_EXPECTED" >> $VALIDATION_OUTPUT
					echo -e "#              ALL " >> $VALIDATION_OUTPUT
					echo -e "# FOUND:       $NB_FOUND" >> $VALIDATION_OUTPUT
					echo -e "# POSITIVE:    $NB_POSITIVE" >> $VALIDATION_OUTPUT
					echo -e "# MISSING:     $NB_MISSING" >> $VALIDATION_OUTPUT
					echo -e "# NOISE:       $NB_NOISE" >> $VALIDATION_OUTPUT
					echo -e "# SENSITIVITY: $SEN%" >> $VALIDATION_OUTPUT
					echo -e "# SPECIFICITY: $SPE%" >> $VALIDATION_OUTPUT
					echo -e "#" >> $VALIDATION_OUTPUT
					
					#### print missing variants lower than coverage
					if (($count)); then
						echo "# MISSING VARIANTS WITH COVERAGE < $COVERAGE X :" >> $VALIDATION_OUTPUT
						cat $VCF_SAMPLE_LOW_COV_FILTER | grep -v '#' >> $VALIDATION_OUTPUT
					fi;					
					
					#### print missing variants higher than coverage
					if (($NB_MISSING)); then
						echo "# MISSING VARIANTS WITH COVERAGE > $COVERAGE X :" >> $VALIDATION_OUTPUT
						cat $VCF_SAMPLE_HIGH_COV_FILTER | grep -v '#' >> $VALIDATION_OUTPUT
						# $TOOL/$BCFTOOLS isec -c all -C $VCF_SAMPLE_filter.gz $VCF_FINAL.gz 2>/dev/null >> $VALIDATION_OUTPUT
					fi;
				

					# count low/high coverage variant 
					SUM_COUNT_LOW_COUNT=$(($SUM_COUNT_LOW_COUNT+$count))
					SUM_COUNT_HIGH_COUNT=$(($SUM_COUNT_HIGH_COUNT+$count_HIGH))

					###  Calcul for all samples in APP
					SUM_NB_EXPECTED=$(($SUM_NB_EXPECTED+$NB_EXPECTED))
					
					# ALL
					SUM_NB_FOUND=$(($SUM_NB_FOUND+$NB_FOUND))
					SUM_NB_POSITIVE=$(($SUM_NB_POSITIVE+$NB_POSITIVE))
					SUM_NB_MISSING=$(($SUM_NB_MISSING+$NB_MISSING))
					SUM_NB_NOISE=$(($SUM_NB_NOISE+$NB_NOISE))
					
					rm -f $VCF_SAMPLE_filter $VCF_SAMPLE_filter.* $VCF_SAMPLE_HIGH_COV_FILTER $VCF_SAMPLE_HIGH_COV_FILTER.* $SAMPLE_VCF.gz* $VCF_FINAL.* $VCF_missing_variant $VCF_SAMPLE_filter_variant $VCF_SAMPLE_LOW_COV_FILTER $VCF_SAMPLE_LOW_COV_FILTER.* $RESULTS/results/$RUN/$SAMPLE_NAME/$SAMPLE_NAME.reports/*.sorted $DATA/SETS_STARK/$SET_NAME/*.sorted
				else
		
					echo "### Sample '$SAMPLE_NAME'" >> $VALIDATION_OUTPUT
					echo "# No VCF '$(basename $VCF_FINAL)' correctly generated" >> $VALIDATION_OUTPUT
									
				fi;

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
			echo -e "\n" >> $VALIDATION_OUTPUT
			echo -e "count low coverage variant = $SUM_COUNT_LOW_COUNT" >> $VALIDATION_OUTPUT
			echo -e "count high coverage variant = $SUM_COUNT_HIGH_COUNT" >> $VALIDATION_OUTPUT
			echo -e "######################################" >> $VALIDATION_OUTPUT
			echo -e "### VALIDATION '$VALIDATION'" >> $VALIDATION_OUTPUT
			echo -e "### SET '$SET_NAME' [$NB_SAMPLE SAMPLES]" >> $VALIDATION_OUTPUT
			echo -e "### APP '$APP'" >> $VALIDATION_OUTPUT
			echo -e "### COVERAGE '$COVERAGE X'" >> $VALIDATION_OUTPUT
			echo -e "######################################" >> $VALIDATION_OUTPUT
			echo -e "# EXPECTED:    $SUM_NB_EXPECTED" >> $VALIDATION_OUTPUT
			echo -e "#              ALL " >> $VALIDATION_OUTPUT
			echo -e "# FOUND:       $SUM_NB_FOUND" >> $VALIDATION_OUTPUT
			echo -e "# POSITIVE:    $SUM_NB_POSITIVE" >> $VALIDATION_OUTPUT
			echo -e "# MISSING:     $SUM_NB_MISSING" >> $VALIDATION_OUTPUT
			echo -e "# NOISE:       $SUM_NB_NOISE" >> $VALIDATION_OUTPUT
			echo -e "# SENSITIVITY: $SEN%" >> $VALIDATION_OUTPUT
			echo -e "# SPECIFICITY: $SPE% " >> $VALIDATION_OUTPUT
			echo -e "######################################" >> $VALIDATION_OUTPUT
			
			## copy to externel server
			if (($COPY_SERVER)); then
				mkdir -p $OUTPUT/$RUN
				rsync -r $RESULTS/$RUN $OUTPUT/
			fi;

		done; # APP
	fi; # SET
done;
