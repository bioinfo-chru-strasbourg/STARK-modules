#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARK_DEJAVU"
SCRIPT_DESCRIPTION="STARK DEJAVU ANNOVAR databases generation"
SCRIPT_RELEASE="0.11.0"
SCRIPT_DATE="29/09/2020"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS"
SCRIPT_LICENCE="GNU-AGPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-05/09/2017: Script creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1b-07/09/2017: Add generation of ANNOVAR generic database.\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.2b-02/11/2018: Use BCFTOOLS instead of VCFTOOLS.\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.10.0-12/08/2020: Many changes.\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.11.0-29/09/2020: Add STARK module json files.\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
ENV_CONFIG=$(find -L $SCRIPT_DIR/.. -name config.app)

source $ENV_CONFIG 1>/dev/null 2>/dev/null

ENV_TOOLS=$(find -L $SCRIPT_DIR/.. -name toold.app)

source $ENV_TOOLS 1>/dev/null 2>/dev/null


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
	echo "# USAGE: $(basename $0) [options...]";
	#echo "# -u/--update                 Update annovar option";
	echo "# -h/--help                    HELP option";
	echo "#";
}

# header
header;


####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "e:r:uvdnh" --long "update,app_folder:,application_folder:,repo_folder:,dejavu_folder:,dejavu_release:,dejavu_annotation:,tmp:,bcftools:,tabix:,bgzip:,annovar:,verbose,debug,release,help" -- "$@" 2> /dev/null)

eval set -- "$ARGS"
while true
do
	case "$1" in
		-u|--update)
			UPDATE=1;
			shift 1
			;;
		-e|--app_folder|--application_folder)
			APP_FOLDER="$2";
			shift 2
			;;
		-r|--repo_folder)
			REPO_FOLDER=$(echo "$2" | tr "," " ");
			shift 2
			;;
		--dejavu_folder)
			DEJAVU_FOLDER="$2";
			shift 2
			;;
		--dejavu_release)
			RELEASE="$2";
			shift 2
			;;
		--dejavu_annotation)
			DEJAVU_ANNOTATION="$2";
			shift 2
			;;
		--tmp)
			TMP_FOLDER_TMP="$2";
			shift 2
			;;
		--bcftools)
			BCFTOOLS="$2";
			shift 2
			;;
		--tabix)
			TABIX="$2";
			shift 2
			;;
		--bgzip)
			BGZIP="$2";
			shift 2
			;;
		--annovar)
			ANNOVAR="$2";
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
		*) 	echo "# Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done


## PARALETERS
##############

# APP FOLDER
if [ ! -z "$APP_FOLDER" ] && [ -d "$APP_FOLDER" ]; then
    STARK_FOLDER_APPS=$APP_FOLDER
fi;


# REPO_FOLDER
REPO_FOLDER=$(echo $REPO_FOLDER | tr "," " ")

# FOLDER DEJAVU
if [ -z "$DEJAVU_FOLDER" ]; then
    DEJAVU_FOLDER="."
fi;

# BCFTOOLS
if [ -z "$BCFTOOLS" ]; then
    BCFTOOLS="bcftools"
fi;

# TABIX
if [ -z "$TABIX" ]; then
    TABIX="tabix"
fi;

# BGZIP
if [ -z "$BGZIP" ]; then
    BGZIP="bgzip"
fi;

# DEJAVU ANNOTATION
#DEJAVU_ANNOTATION=ALL
# if [ -z "$DEJAVU_ANNOTATION" ]; then
#     DEJAVU_ANNOTATION="ALL"
# fi;

if [ -z "$RELEASE" ]; then
	RELEASE=$(date +%Y%m%d-%H%M%S)
fi;


# DEJAVU
DEJAVU=$DEJAVU_FOLDER
mkdir -p $DEJAVU/$RELEASE

DEJAVU_SUBFOLDER_ANNOVAR=annovar
DEJAVU_FOLDER_ANNOVAR=$DEJAVU/$RELEASE/$DEJAVU_SUBFOLDER_ANNOVAR
mkdir -p $DEJAVU_FOLDER_ANNOVAR

DEJAVU_SUBFOLDER_VCF=vcf
DEJAVU_FOLDER_VCF=$DEJAVU/$RELEASE/$DEJAVU_SUBFOLDER_VCF
mkdir -p $DEJAVU_FOLDER_VCF

DEJAVU_SUBFOLDER_LOG=log
DEJAVU_FOLDER_LOG=$DEJAVU/$RELEASE/$DEJAVU_SUBFOLDER_LOG
mkdir -p $DEJAVU_FOLDER_LOG

DEJAVU_SUBFOLDER_STATS=stats
#DEJAVU_FOLDER_STATS=$DEJAVU/$RELEASE/$DEJAVU_SUBFOLDER_STATS
#mkdir -p $DEJAVU_FOLDER_STATS



# DEJAVU FOLDER
#mkdir -p $DEJAVU_FOLDER_LOG
# mkdir -p $DEJAVU/$RELEASE/annovar
# mkdir -p $DEJAVU/$RELEASE/vcf
# mkdir -p $DEJAVU/$RELEASE/log


# ASSEMBLY PREFIX

ASSEMBLY_PREFIX_DEFAULT="hg19"


# VCF pattern

VCF_PATTERN=".final.vcf.gz"


# INFO SUFFIX (e.g. release): $PREFIX_dejavu.$GROUP.$PROJECT$SUFFIX.txt
#SUFFIX=".$RELEASE"
SUFFIX=""

# TMP
if [ "$TMP_FOLDER_TMP" == "" ]; then TMP_FOLDER_TMP=/tmp; fi;
TMP=$TMP_FOLDER_TMP/dejavu_$RANDOM$RANDOM$RANDOM$RANDOM
mkdir -p $TMP


# MK
MK=$DEJAVU_FOLDER_LOG/$RELEASE.mk
> $MK

# LOG
LOG=$DEJAVU_FOLDER_LOG/$RELEASE.log
> $LOG


#### INFO
(($VERBOSE)) && echo "#[INFO] RELEASE: $RELEASE"
(($VERBOSE)) && echo "#[INFO] STARK FOLDER APPLICATIONS: $STARK_FOLDER_APPS"
(($VERBOSE)) && echo "#[INFO] DEJAVU FOLDER: $DEJAVU_FOLDER"
(($VERBOSE)) && echo "#[INFO] REPOSITORY FOLDER: $REPO_FOLDER"

(($DEBUG)) && echo "#[INFO] TMP: $TMP"



### Find Group folders
########################

GP_FOLDER_LIST=""

if [ ! -z "$REPO_FOLDER" ]; then
	#GP_FOLDER_LIST=$(find -L $REPO_FOLDER -maxdepth 2 -mindepth 2 -type d 2>/dev/null)
	GP_FOLDER_LIST=$(ls $(ls $REPO_FOLDER/*/*/*/STARKCopyComplete.txt | xargs dirname | xargs dirname | sort -u | sed 's#$#/*/*/*'$VCF_PATTERN'#g') | xargs dirname | xargs dirname | xargs dirname | sort -u)
fi;



(($VERBOSE)) && echo "#"
(($VERBOSE)) && echo "#[INFO] DEJAVU database repository detection"


if [ -z "$GP_FOLDER_LIST" ]; then
    for ENV_DEF in $(find -L $STARK_FOLDER_APPS -name '*.app' -type f | sed s#$STARK_FOLDER_APPS/## | sort -f -t'/' -k2.3 -k2.2 -k2.1) $(find -L $STARK_FOLDER_APPS -name '*.plugapp' -type f | sed s#$STARK_FOLDER_APPS/## | sort -f -t'/' -k2.3 -k2.2 -k2.1); do
    
    	# APP INFO
        APP_FOLDER_ARCHIVES=$(source_app "$ENV_DEF"  2>/dev/null; echo $FOLDER_ARCHIVES)
        APP_FOLDER_REPOSITORY=$(source_app "$ENV_DEF"  2>/dev/null; echo $FOLDER_REPOSITORY)
        APP_GROUP=$(source_app "$ENV_DEF"  2>/dev/null; echo $APP_GROUP)
        APP_PROJECT=$(source_app "$ENV_DEF"  2>/dev/null; echo $APP_PROJECT)

        # REPOSITORY FOLDER
        if [ ! -z "$REPO_FOLDER" ] && [ -d "$REPO_FOLDER" ] && [ -d "$REPO_FOLDER/$APP_GROUP/$APP_PROJECT" ]; then
        	APP_FOLDER_ARCHIVES=$REPO_FOLDER
        fi;

        # GROUP FOLDER
        if [ -z "$APP_GROUP" ]; then
            APP_GROUP="UNKNOWN"
        fi;

        # PROJECT FOLDER
        if [ -z "$APP_PROJECT" ]; then
            APP_PROJECT="UNKNOWN"
        fi;
        
        # Add repository group project folder 
        if [ ! -z "$APP_FOLDER_ARCHIVES" ] && [ ! -z "$APP_GROUP" ] && [ ! -z "$APP_PROJECT" ] && [ -d "$APP_FOLDER_ARCHIVES/$APP_GROUP/$APP_PROJECT" ]; then
	        (($VERBOSE)) && echo "#[INFO] DEJAVU database '$APP_GROUP/$APP_PROJECT' repository '$APP_FOLDER_ARCHIVES' found"
	        GP_FOLDER_LIST="$GP_FOLDER_LIST\n$APP_FOLDER_ARCHIVES/$APP_GROUP/$APP_PROJECT"
	    fi;
	    if [ ! -z "$APP_FOLDER_REPOSITORY" ] && [ ! -z "$APP_GROUP" ] && [ ! -z "$APP_PROJECT" ] && [ -d "$APP_FOLDER_REPOSITORY/$APP_GROUP/$APP_PROJECT" ]; then
	        (($VERBOSE)) && echo "#[INFO] DEJAVU database '$APP_GROUP/$APP_PROJECT' repository '$APP_FOLDER_REPOSITORY' found"
	        GP_FOLDER_LIST="$GP_FOLDER_LIST\n$APP_FOLDER_REPOSITORY/$APP_GROUP/$APP_PROJECT"
	    fi;
        
        
    done;

else

	for GP_FOLDER in $GP_FOLDER_LIST; do
		if [ -d "$GP" ]; then
			REPO=$(dirname $(dirname "$GP_FOLDER"))
			GROUP=$(basename $(dirname "$GP_FOLDER"))
			PROJECT=$(basename "$GP_FOLDER")

	        (($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' repository '$REPO' "
	        #GP_FOLDER_LIST="$GP_FOLDER_LIST\n$APP_FOLDER_ARCHIVES/$APP_GROUP/$APP_PROJECT"
	    fi;
	done;

fi;

GP_FOLDER_LIST_UNIQ=$(echo -e $GP_FOLDER_LIST | sort -u)

GP_FOLDER_LIST_UNIQ_COUNT=$(echo $GP_FOLDER_LIST_UNIQ | wc -w)


#(($VERBOSE)) && echo "#[INFO] DEJAVU repository number: $GP_FOLDER_LIST_UNIQ_COUNT"

#echo $GP_FOLDER_LIST_UNIQ



### DEJAVU database copy file
##############################

(($VERBOSE)) && echo "#"
(($VERBOSE)) && echo "#[INFO] DEJAVU database file copy"

(($DEBUG)) && echo "#[INFO] DEJAVU database file pattern '$VCF_PATTERN'"



for GP_FOLDER in $GP_FOLDER_LIST_UNIQ; do

	REPO=$(dirname $(dirname $GP_FOLDER))
	GROUP=$(basename $(dirname $GP_FOLDER))
	PROJECT=$(basename $GP_FOLDER)

	# NB VARIANT
	#NB_VCF=$(find $GP_FOLDER/*/*/ -maxdepth 1 -name '*.final.vcf.gz' -a ! -name '*.*-*.final.vcf.gz' 2>/dev/null | wc -l)
	NB_VCF=$(find -L $GP_FOLDER/*/*/ -maxdepth 1 -name '*'$VCF_PATTERN -a ! -name '*.*-*'$VCF_PATTERN 2>/dev/null | wc -l)

	(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' repository '$REPO' $NB_VCF VCF found"

	if [ $NB_VCF -gt 0 ]; then
		
		> $MK.$GROUP.$PROJECT.log
		> $MK.$GROUP.$PROJECT.err

		# TMP folder creation
		mkdir -p $TMP/$GROUP/$PROJECT
		cp -f $(find -L $GP_FOLDER/*/*/ -maxdepth 1 -name '*'$VCF_PATTERN -a ! -name '*.*-*'$VCF_PATTERN) $TMP/$GROUP/$PROJECT/ 2>/dev/null
		NB_VCF_FOUND=$(ls $TMP/$GROUP/$PROJECT/*.vcf.gz | wc -w)
		(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' $NB_VCF_FOUND VCF considered files (some files may be found multiple times)"
		
	fi;

done;



### Database generation process

(($VERBOSE)) && echo "#"
(($VERBOSE)) && echo "#[INFO] DEJAVU database generation process"

for GP_FOLDER in $GP_FOLDER_LIST_UNIQ; do

	REPO=$(dirname $(dirname "$GP_FOLDER"))
	GROUP=$(basename $(dirname "$GP_FOLDER"))
	PROJECT=$(basename "$GP_FOLDER")

	#(($VERBOSE)) && echo "#"
	#(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' process"



	if [ ! -s "$DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.done" ]; then

		# MK files
		> $MK
		echo "%.vcf.gz.tbi: %.vcf.gz
			$TABIX $<

		" > $MK

		VCFGZ_LIST=""
		#for VCF in $(ls $TMP/$GROUP/$PROJECT/*); do
		VCFGZ_NB=0
		for VCF in $TMP/$GROUP/$PROJECT/*; do
			#echo "VCF: "$VCF
			if [ -s $VCF ] && (($(grep ^# -cv $VCF))); then
				echo "$VCF.simple.vcf.gz: $VCF
					mkdir $<.sort.
					$BCFTOOLS norm -m -any -c s --fasta-ref $GENOMES/current/$ASSEMBLY.fa $< | $BCFTOOLS norm --rm-dup=exact | $BCFTOOLS +fixploidy  -- -f 2 | $BCFTOOLS +setGT  -- -t . -n 0 | $BCFTOOLS sort -T $<.sort. 2>/dev/null | $BCFTOOLS annotate -x FILTER,QUAL,ID,INFO -o \$@ -O z 2>/dev/null;
					$TABIX \$@;
					rm -f $<.sort.*
				" >> $MK
				VCFGZ_LIST="$VCFGZ_LIST $VCF.simple.vcf.gz"
				((VCFGZ_NB++))
				#| sed s/ID=PL,Number=G/ID=PL,Number=./gi ,^FORMAT/GT
			fi;
		done

		
		echo "$TMP/$GROUP/$PROJECT/dejavu.simple.vcf: $VCFGZ_LIST
		if [ $VCFGZ_NB -gt 1 ]; then \
			$BCFTOOLS merge $TMP/$GROUP/$PROJECT/*.simple.vcf.gz | $BCFTOOLS norm -m -any | $BCFTOOLS norm --rm-dup=exact | $BCFTOOLS +setGT  -- -t . -n 0 | $BCFTOOLS +fill-tags -- -t AN,AC,AF,AC_Hemi,AC_Hom,AC_Het,ExcHet,HWE,MAF,NS > \$@; \
		else \
			$BCFTOOLS norm -m -any $VCFGZ_LIST | $BCFTOOLS +setGT  -- -t . -n 0 | $BCFTOOLS +fill-tags -- -t AN,AC,AF,AC_Hemi,AC_Hom,AC_Het,ExcHet,HWE,MAF,NS > \$@; \
		fi;
		" >> $MK

		echo "$TMP/$GROUP/$PROJECT/dejavu.vcf: $TMP/$GROUP/$PROJECT/dejavu.simple.vcf
			cp $TMP/$GROUP/$PROJECT/dejavu.simple.vcf $TMP/$GROUP/$PROJECT/dejavu.vcf
		" >> $MK

		if [ "$DEJAVU_ANNOTATION" != "" ]; then
			echo "$TMP/$GROUP/$PROJECT/dejavu.annotated.vcf: $TMP/$GROUP/$PROJECT/dejavu.simple.vcf
				+$HOWARD --input=$< --output=\$@.tmp.vcf --config=$HOWARD_CONFIG --config_annotation=$HOWARD_CONFIG_ANNOTATION  --annotation=$DEJAVU_ANNOTATION --calculation=$HOWARD_CALCULATION --nomen_fields=$HOWARD_NOMEN_FIELDS --annovar_folder=$ANNOVAR --annovar_databases=$ANNOVAR_DATABASES --snpeff_jar=$SNPEFF --snpeff_databases=$SNPEFF_DATABASES --multithreading --threads=$THREADS --snpeff_threads=$THREADS --tmp=$TMP_FOLDER_TMP --env=$CONFIG_TOOLS
				mkdir $<.sort.
				$BCFTOOLS sort -T $<.sort. \$@.tmp.vcf | $BCFTOOLS annotate -x INFO/AN,INFO/AC,INFO/AF  | $BCFTOOLS +fill-tags -- -t AN,AC,AF,AC_Hemi,AC_Hom,AC_Het,ExcHet,HWE,MAF,NS > \$@;
				rm \$@.tmp.vcf
				rm -f $<.sort.*
			" >> $MK
			# --multithreading --threads=$THREADS
			# --norm_options='--multiallelics=-any,--rm-dup=exact'

		else
			echo "$TMP/$GROUP/$PROJECT/dejavu.annotated.vcf: $TMP/$GROUP/$PROJECT/dejavu.simple.vcf
				cp $TMP/$GROUP/$PROJECT/dejavu.simple.vcf $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf
			" >> $MK
		fi;

		echo "$TMP/$GROUP/$PROJECT/VCF_LIST: $VCFGZ_LIST
		ls $^ > $TMP/$GROUP/$PROJECT/VCF_LIST

		" >> $MK

		echo "%.vcf.gz: %.vcf
			$BGZIP -c $< > \$@

		" >> $MK

		echo "%.gz.tbi: %.gz
			$TABIX $<

		" >> $MK

		echo "%.tsv: %.vcf
			+$HOWARD --input=$< --output=\$@ --env=$CONFIG_TOOLS --translation=TSV

		" >> $MK


		#echo "awk -F'\t' '{AF=\$\$6/($NB_VCF_FOUND*2)} {print \$\$1\"\\t\"\$\$2\"\\t\"\$\$3\"\\t\"\$\$4\"\\t\"AF}'"
		#exit 0

		#echo "$DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.txt: $TMP/$GROUP/$PROJECT/dejavu.vcf
		echo "$TMP/$GROUP/$PROJECT/dejavu.percent: $TMP/$GROUP/$PROJECT/dejavu.vcf
			$BCFTOOLS query -f'%CHROM\t%POS\t%REF\t%ALT\t%AF\t%NS\t%AN\t%AC\t%AC_Hom\t%AC_Het\n' $< > \$@

		" >> $MK
		#$BCFTOOLS query -f'%CHROM\t%POS\t%REF\t%ALT\t%AF\t%NS\t%AN\t%AC\t%AC_Hom\t%AF_Het\n' $< > \$@
		echo "$TMP/$GROUP/$PROJECT/dejavu.annovar: $TMP/$GROUP/$PROJECT/dejavu.vcf
			perl $ANNOVAR/convert2annovar.pl --format vcf4old  --allallele --outfile \$@.tmp $<
			cat \$@.tmp | cut -f1-5 > \$@
			rm -f \$@.tmp

		" >> $MK

		echo "$TMP/$GROUP/$PROJECT/dejavu.annovar.percent: $TMP/$GROUP/$PROJECT/dejavu.annovar $TMP/$GROUP/$PROJECT/dejavu.percent
			paste $^ | cut -f1-5,10-15 > \$@

		" >> $MK
		
		echo "$TMP/$GROUP/$PROJECT/dejavu.txt: $TMP/$GROUP/$PROJECT/dejavu.annovar.percent
			#perl $ANNOVAR/index_annovar.sh $< --outfile \$@
			perl $SCRIPT_DIR/index_annovar.pl $< --outfile \$@

		" >> $MK

		#cat $MK
		#(($VERBOSE)) && echo "#[INFO] DEJAVU database generation process..."
		make -j $THREADS -f $MK $TMP/$GROUP/$PROJECT/dejavu.txt $TMP/$GROUP/$PROJECT/dejavu.vcf $TMP/$GROUP/$PROJECT/dejavu.tsv $TMP/$GROUP/$PROJECT/dejavu.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.vcf.gz.tbi $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi $TMP/$GROUP/$PROJECT/dejavu.annotated.tsv 1>>$MK.$GROUP.$PROJECT.log 2>>$MK.$GROUP.$PROJECT.err
		# grep "\*\*\*" $MK.err
		(($DEBUG)) && grep "\*\*\*" $MK.$GROUP.$PROJECT.log -B30
		(($DEBUG)) && grep "\*\*\*" $MK.$GROUP.$PROJECT.err -B30
		
		# $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi 

		# echo
		if (($(cat $MK.$GROUP.$PROJECT.log $MK.$GROUP.$PROJECT.err | grep "\*\*\*" -c))); then
			echo "#[ERROR] File '$DEJAVU/$RELEASE/dejavu.$GROUP.$PROJECT.txt' generation..."
			(($DEBUG)) && cat $MK.$GROUP.$PROJECT.log $MK.$GROUP.$PROJECT.err | grep "\*\*\*" -B 80
			exit 1
		else
		
			# STATS
			> $TMP/$GROUP/$PROJECT/dejavu.stats.txt

			# SAMPLES
			NB_SAMPLES=$(echo $(grep "^#CHROM" $TMP/$GROUP/$PROJECT/dejavu.vcf | wc -w)" - 9" | bc)
			NB_VARIANTS=$(grep -cv "^#" $TMP/$GROUP/$PROJECT/dejavu.vcf)
			
			echo "### Number of samples: $NB_SAMPLES" >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
			echo "### Number of variants: $NB_VARIANTS" >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
			
			# VARIANTS
			echo "### variants frequency" >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
			# for p in $(seq 0 100); do
			#for p in 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do awk -F"\t" -v p=$p '$6>p{SUM++} {TOT++} END {PERC=((SUM+0)/TOT)*100; print "# "SUM+0" variants out of "TOT" ("PERC"%) found more than "(p*100)"% in the set"}' $TMP/$GROUP/$PROJECT/dejavu.txt; done >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
			for p in $(seq 1 100); do awk -F"\t" -v p=$( echo "scale=2; $p/100" | bc) '$6>p{SUM++} {TOT++} END {PERC=((SUM+0)/TOT)*100; print "# "SUM+0" variants out of "TOT" ("PERC"%) found more than "(p*100)"% in the set"}' $TMP/$GROUP/$PROJECT/dejavu.txt; done >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
			#cat $TMP/$GROUP/$PROJECT/dejavu.stats.txt

			# TSV
			echo "$NB_SAMPLES" >> $TMP/$GROUP/$PROJECT/dejavu.stats.nb_samples
			echo "$NB_VARIANTS" >> $TMP/$GROUP/$PROJECT/dejavu.stats.nb_variants

			echo -e "#NB_variant\tTOT_variant\tPercent\tmore_than" > $TMP/$GROUP/$PROJECT/dejavu.stats.tsv
			#for p in 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do awk -F"\t" -v p=$p '$6>p{SUM++} {TOT++} END {PERC=((SUM+0)/TOT)*100; print SUM+0"\t"TOT"\t"PERC"\t"(p*100)}' $TMP/$GROUP/$PROJECT/dejavu.txt; done > $TMP/$GROUP/$PROJECT/dejavu.stats.frequency.tsv
			for p in $(seq 1 100); do awk -F"\t" -v p=$( echo "scale=2; $p/100" | bc) '$6>p{SUM++} {TOT++} END {PERC=((SUM+0)/TOT)*100; print SUM+0"\t"TOT"\t"PERC"\t"(p*100)}' $TMP/$GROUP/$PROJECT/dejavu.txt; done >> $TMP/$GROUP/$PROJECT/dejavu.stats.tsv
			#cat $TMP/$GROUP/$PROJECT/dejavu.stats.tsv



			# find assembly prefix
			ASSEMBLY_PREFIX_VCF=$(grep "^##reference=" $TMP/$GROUP/$PROJECT/dejavu.vcf | cut -d= -f2 | xargs basename | sed 's/.fa$//' | sed 's/.fasta$//')
			if [ "$ASSEMBLY_PREFIX_VCF" == "" ]; then
				ASSEMBLY_PREFIX=$ASSEMBLY_PREFIX_DEFAULT
			else 
				ASSEMBLY_PREFIX=$ASSEMBLY_PREFIX_VCF
			fi;
			ASSEMBLY_PREFIX_ANNOVAR=$ASSEMBLY_PREFIX"_"
			
			# SUFFIX
			SUFFIX_ANNOVAR=$SUFFIX
			
			PATTERN=$ASSEMBLY_PREFIX_ANNOVAR"dejavu."$GROUP.$PROJECT$SUFFIX_ANNOVAR
			mkdir -p $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/

			cp $TMP/$GROUP/$PROJECT/dejavu.txt $DEJAVU_FOLDER_ANNOVAR/$PATTERN.txt
			cp $TMP/$GROUP/$PROJECT/dejavu.txt.idx $DEJAVU_FOLDER_ANNOVAR/$PATTERN.txt.idx

			cp $TMP/$GROUP/$PROJECT/dejavu.vcf $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/minimal.vcf
			cp $TMP/$GROUP/$PROJECT/dejavu.tsv $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/minimal.tsv
			cp $TMP/$GROUP/$PROJECT/dejavu.vcf.gz $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/minimal.vcf.gz
			cp $TMP/$GROUP/$PROJECT/dejavu.vcf.gz.tbi $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/minimal.vcf.gz.tbi
			cp $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/annotated.vcf
			cp $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/annotated.vcf.gz
			cp $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/annotated.vcf.gz.tbi
			cp $TMP/$GROUP/$PROJECT/dejavu.annotated.tsv $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/annotated.tsv
			cp $TMP/$GROUP/$PROJECT/dejavu.annotated.tsv $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/annotated.tsv

			mkdir -p $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/$DEJAVU_SUBFOLDER_STATS
			cp $TMP/$GROUP/$PROJECT/dejavu.stats* $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/$DEJAVU_SUBFOLDER_STATS/

			(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' generated with $NB_SAMPLES samples and $NB_VARIANTS variants"
			

			# end
			#(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' generated for release $RELEASE"
			echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' generated for release $RELEASE" > $DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.done
			echo "#[INFO] release=$RELEASE" >> $DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.done
			echo "#[INFO] samples=$NB_SAMPLES" >> $DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.done
			echo "#[INFO] variants=$NB_VARIANTS" >> $DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.done
			(($VERBOSE)) && echo "$RELEASE" > $DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.release
			

			# STARK module json

			# Database definition
			if [ ! -s $DEJAVU/STARK.database ]; then
				echo '

					{
						"code": "dejavu",
						"name": "DejaVu",
						"fullname": "STARK DejaVu databases",
						"website": "",
						"description": "STARK DejaVu databases is a compilation of all samples variants for each group/project, useful to calculate population frequencies"
					}
			
				' > $DEJAVU/STARK.database
			fi;

			# Release information
			echo '

					{
						"release": "'$RELEASE'",
						"date": "'$RELEASE'",
						"files": [ "release" ],
						"assembly": [ "'$ASSEMBLY'" ],
						"download": {
							"methode": "DejaVu Databases generation script ['$SCRIPT_RELEASE'-'$SCRIPT_DATE']",
							"date": "'.$(date).'"
						}
					}
		

		
			' > $DEJAVU/$RELEASE/STARK.database.release


			# Latest
			rm -f $DEJAVU/latest
			ln -s $RELEASE/ $DEJAVU/latest

			# CLeaning
			rm -rf $TMP/$GROUP/$PROJECT


		fi;

	# else

	# 	(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' already generated"

	fi;
	#fi;
done

# TODO: Add STARK.database and STARK.database.release



(($VERBOSE)) && echo "#"
(($VERBOSE)) && echo "#[INFO] DEJAVU database release '$RELEASE' done."
(($VERBOSE)) && echo "$RELEASE" > $DEJAVU/$RELEASE/release

rm -rf $TMP

exit 0;
