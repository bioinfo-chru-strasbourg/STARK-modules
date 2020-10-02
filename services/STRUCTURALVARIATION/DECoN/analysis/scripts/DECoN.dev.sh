

function INSTALL {

    # Based on STARK Image (docker)

    # PARAM
    THREADS=12


    # YUM INSTALL
    yum -y install gcc-gfortran zlib-devel #R


    # R INSTALL
    mkdir -p /STARK/tools/R
    cd /STARK/tools/R
    wget https://cran.r-project.org/src/base/R-3/R-3.1.2.tar.gz
    tar zxvf R-3.1.2.tar.gz; cd R-3.1.2/ 
    ./configure --with-readline=no --with-x=no --prefix=/STARK/tools/R/3.1.2
    make -j $THREADS;
    make -j $THREADS install
    export PATH=$PATH":/STARK/tools/R/3.1.2/bin"
    ln -s /STARK/tools/R/3.1.2/ /STARK/tools/R/current # ERROR

    # R INSTALL package
    wget https://cran.r-project.org/src/contrib/Rcpp_1.0.5.tar.gz
    R CMD INSTALL Rcpp_1.0.5.tar.gz



    # DECoN INSTALL
    mkdir -p /STARK/tools/decon
    cd /STARK/tools/decon
    wget https://github.com/RahmanTeam/DECoN/archive/v1.0.2.tar.gz
    tar -xvf v1.0.2.tar.gz
    mkdir 1.0.2
    ln -s 1.0.2/ current
    cp -prf DECoN-1.0.2/Linux/* 1.0.2/
    cd 1.0.2
    mkdir -p packrat/src/VGAM packrat/src/scales packrat/src/ggplot2
    wget https://cran.r-project.org/src/contrib/Archive/VGAM/VGAM_0.9-8.tar.gz -O packrat/src/VGAM/VGAM_0.9-8.tar.gz
    wget https://cran.r-project.org/src/contrib/Archive/scales/scales_0.2.4.tar.gz -O packrat/src/scales/scales_0.2.4.tar.gz
    wget https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_1.0.1.tar.gz -O packrat/src/ggplot2/ggplot2_1.0.1.tar.gz
    #./setup.sh
    [ -w ".Rprofile" ] && rm .Rprofile
    > setup.log
    #timeout 1000 tail -f setup.log &
    /STARK/tools/R/3.1.2/bin/Rscript sessionInfo.R --bootstrap-packrat > setup.log 2>&1
    cat setup.log
    cp packrat/packrat_source/.Rprofile ./
    # #TODO

    export PATH=$PATH":/STARK/tools/decon/1.0.2"

    # /STARK/tools/R/3.1.2/bin/Rscript /STARK/tools/decon/1.0.2/sessionInfo.R --bootstrap-packrat 

}

#INSTALL


### FUNCTIONS


function DECON_INIT {

    # ANALYSIS
    ANALYSIS=$(date +%Y%m%d-%H%M%S)

    # OUTPUT
    #OUTPUT_FOLDER=$INPUT_FOLDER/DECoN
    #OUTPUT_FOLDER=$( echo "$INPUT_FOLDER" | tr "*" "_")/DECoN/$ANALYSIS
    OUTPUT_FOLDER=$DECoN_main_results/$DECoN_sub_results/$( basename "$INPUT_FOLDER" | tr "*" "_")/$ANALYSIS
    LOG=$OUTPUT_FOLDER/DECoN.log
    ERR=$OUTPUT_FOLDER/DECoN.err
    BAMS_FILE=$OUTPUT_FOLDER/bams.file
    BED_FILE=$OUTPUT_FOLDER/design.bed
    PARAMS=$OUTPUT_FOLDER/DECoN.params
    OUTPUT_PREFIX=$OUTPUT_FOLDER/DECoN
    plotFolder=$OUTPUT_PREFIX"Plots"
    Rdata=$OUTPUT_PREFIX.RData
    

    # output folder
    if [ -d $OUTPUT_FOLDER ]; then
        mv $OUTPUT_FOLDER $OUTPUT_FOLDER"_moved_from_V"$ANALYSIS
    fi;
    mkdir -p $OUTPUT_FOLDER

    # LOG
    > $LOG
    > $ERR

    # PARAMS

    # Ressources
    echo "# DECoN Params" > $PARAMS
    echo "FASTA=$FASTA" >> $PARAMS
    echo "refGene=$refGene" >> $PARAMS

    # Data
    echo "design_ext=$design_ext" >> $PARAMS
    echo "BED_ADD_REFGENE=$BED_ADD_REFGENE" >> $PARAMS
    echo "bam_find_min_depth=$bam_find_min_depth" >> $PARAMS
    echo "bam_find_max_depth=$bam_find_max_depth" >> $PARAMS
    echo "bed_find_min_depth=$bed_find_min_depth" >> $PARAMS
    echo "bed_find_max_depth=$bed_find_max_depth" >> $PARAMS
    echo "BED_ADD_REFGENE=$BED_ADD_REFGENE" >> $PARAMS

    # PARAMS default
    echo "MINCORR=$MINCORR" >> $PARAMS
    echo "MINCOV=$MINCOV" >> $PARAMS
    echo "transProb=$transProb" >> $PARAMS
    echo "LOG_UP=$LOG_UP" >> $PARAMS
    echo "LOG_DOWN=$LOG_DOWN" >> $PARAMS
    echo "BF=$BF" >> $PARAMS
    echo "NB_EXON=$NB_EXON" >> $PARAMS
    echo "NB_COMP=$NB_COMP" >> $PARAMS

    

}

function DECON_ANALYSIS {

    (($VERBOSE)) && echo "#[INFO] OUTPUT: $OUTPUT_FOLDER"
    
    # Prepare bams.files and design

    # BAM
    (($VERBOSE)) && echo "#[INFO] BAMS list file generation" 

    # Exclude bam pattern
    bam_exclude_find_param="";
    for exclude_pattern in $bam_exclude_pattern; do
        #echo $exclude_pattern;
        bam_exclude_find_param=$bam_exclude_find_param" -and -not -name *$exclude_pattern*bam" ;
    done;
    #echo $bam_exclude_find_param

    find $INPUT_FOLDER -mindepth $bam_find_min_depth -maxdepth $bam_find_max_depth -name "*bam" $bam_exclude_find_param > $BAMS_FILE
    #find $INPUT_FOLDER -mindepth $bam_find_min_depth -maxdepth $bam_find_max_depth -name "*.validation.bam" $bam_exclude_find_param > $BAMS_FILE
    (($VERBOSE)) && echo "#[INFO] BAMS_FILE=$BAMS_FILE" && head $BAMS_FILE
    (($VERBOSE)) && echo "..."
    NB_BAM=$(cat $BAMS_FILE | wc -l)
    (($VERBOSE)) && echo "#[INFO] NB BAMS: " && echo $NB_BAM

    # BED
    (($VERBOSE)) && echo "#[INFO] BED generation" 
    #cat $(find $INPUT_FOLDER -name *bed) | bedtools sort | bedtools merge -c 4,5 -o distinct,distinct > $BED_FILE
    #cat $(find $INPUT_FOLDER -name *bed) | bedtools sort | bedtools merge -c 4,5 -o distinct,distinct | grep -v '^chrX\|^chrY\|^chrM' > $BED_FILE
    #cat $(find $INPUT_FOLDER -name *bed) | bedtools sort | bedtools merge -c 4,5 -o distinct,distinct | grep -v '^chrX\|^chrY\|^chrM' | awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $4}' > $BED_FILE
    #cat $(find $INPUT_FOLDER -name *bed) | bedtools sort | bedtools merge -c 4,5 -o distinct,distinct | grep -v '^chrX\|^chrY\|^chrM' | awk '{print $1 "\t" $2+1 "\t" $3 "\t" substr($5,1,20) "\t" $4}' | tr "," "_" > $BED_FILE
    #cat $(find $INPUT_FOLDER -name *genes  | grep -v ".list.genes$") | bedtools sort | bedtools merge -c 4,5 -o distinct,distinct | grep -v '^chrX\|^chrY\|^chrM' | awk '{print $1 "\t" $2+1 "\t" $3 "\t" substr($5,1,20) "\t" "0" "\t" $4}' | tr "," "_" > $BED_FILE
    # substr($2,1,length($2)-7)
    # bedtools intersect -wa -a /STARK/data/DECoN/datasets/refGene.hg19.ok.bed -b /STARK/data/DECoN/datasets/BBS_v1_Run4_062012/dataset.bed | bedtools sort | bedtools merge -c 4 -o distinct
    # bedtools intersect -wa -a /STARK/data/DECoN/datasets/BBS_v1_Run4_062012/dataset.bed -b /STARK/data/DECoN/datasets/refGene.hg19.ok.bed | bedtools sort | bedtools merge -c 4 -o distinct
    # cat /STARK/data/DECoN/datasets/BBS_v1_Run4_062012/dataset.bed | awk '{print $0 "\t" $1"_"$2"_"$3}' | bedtools intersect -wa -a - -b /STARK/data/DECoN/datasets/refGene.hg19.ok.bed | bedtools sort | bedtools merge -c 8 -o distinct
    # cat /STARK/data/DECoN/datasets/BBS_v1_Run4_062012/dataset.bed | awk '{print $0 "\t" $1"_"$2"_"$3}' | bedtools intersect -wa -wb -a - -b /STARK/data/DECoN/datasets/refGene.hg19.ok.bed | bedtools sort | bedtools merge -c 8 -o distinct



    if (($BED_ADD_REFGENE)); then
        (($VERBOSE)) && echo "#[INFO] Intersect with refGene" 
        cat $(find $INPUT_FOLDER -mindepth $bed_find_min_depth -maxdepth $bed_find_max_depth -name "*$design_ext" | grep -v ".list.genes$") | awk '{chr=$1} {if ($2>=$3) {start=$3-1} else {start=$2} } {stop=$3} {gene=$4} {print chr "\t" start "\t" stop}' | awk -f /tool/bin/bed_normalization.awk  | bedtools sort | bedtools merge -c 4 -o distinct | bedtools intersect -wa -wb -a - -b $refGene | bedtools sort | bedtools merge -c 4,5,6,7,8 -o distinct,distinct,distinct,distinct,distinct  | awk '{print $1 "\t" $2+1 "\t" $3 "\t" $8}'  | tr ",;-" "_" > $BED_FILE
    else
        cat $(find $INPUT_FOLDER -mindepth $bed_find_min_depth -maxdepth $bed_find_max_depth -name "*$design_ext" | grep -v ".list.genes$") | awk '{chr=$1} {if ($2>=$3) {start=$3-1} else {start=$2} } {stop=$3} {gene=$4} {print chr "\t" start "\t" stop "\t" gene}' | awk -f /tool/bin/bed_normalization.awk  | bedtools sort | bedtools merge -c 4 -o distinct | awk '{print $1 "\t" $2+1 "\t" $3 "\t" substr($4,1,20) "\t" "0" "\t" "+"}' | tr ",;-" "_" > $BED_FILE
    fi;

    (($VERBOSE)) && echo "#[INFO] BED_FILE=$BED_FILE" && head $BED_FILE && echo "..." && tail $BED_FILE
    (($VERBOSE)) && echo "#[INFO] Number of regions" && echo $(cat $BED_FILE | wc -l)


    # DECoN folder
    cd $DECoN_bin_folder;

    # ReadInBams
    (($VERBOSE)) && echo "#[INFO] DECoN ReadInBams" 
    $Rscript $DECoN_bin_folder/ReadInBams.R --bams $BAMS_FILE --bed $BED_FILE --fasta $FASTA --out $OUTPUT_PREFIX


    # IdentifyFailures
    (($VERBOSE)) && echo "#[INFO] DECoN IdentifyFailures" 
    > $OUTPUT_PREFIX"_Failures.txt"
    $Rscript $DECoN_bin_folder/IdentifyFailures.R --Rdata $Rdata --mincorr $MINCORR --mincov $MINCOV --out $OUTPUT_PREFIX #--exons customNumbers.file --custom FALSE --out output.prefix
    (($VERBOSE)) && echo "#[INFO] DECoN IdentifyFailures Failures" 
    (($VERBOSE)) && head $OUTPUT_PREFIX"_Failures.txt"
    (($VERBOSE)) && echo "..." 
    NB_Failures=$(grep ^CNV -v $OUTPUT_PREFIX"_Failures.txt" | wc -l)
    (($VERBOSE)) && echo "#[INFO] DECoN IdentifyFailures NB Failures: $NB_Failures"


    # makeCNVcalls
    (($VERBOSE)) && echo "#[INFO] DECoN makeCNVcalls" 
    $Rscript $DECoN_bin_folder/makeCNVcalls.R --Rdata $Rdata --transProb $transProb --out $OUTPUT_PREFIX --plot All --plotFolder $plotFolder  #--exons customNumbers.file --custom FALSE

    
    if [ ! -e $OUTPUT_PREFIX"_all.txt" ]; then

        (($VERBOSE)) && echo "#[ERROR] DECoN makeCNVcalls failed" 

    else 

        (($VERBOSE)) && echo "#[INFO] DECoN makeCNVcalls All CNV" 
        (($VERBOSE)) && head $OUTPUT_PREFIX"_all.txt"
        (($VERBOSE)) && echo "..." 
        NB_CNV=$(grep ^CNV -v $OUTPUT_PREFIX"_all.txt" | wc -l)
        (($VERBOSE)) && echo "#[INFO] DECoN makeCNVcalls NB CNV: $NB_CNV"

        # Filtration
        grep ^CNV $OUTPUT_PREFIX"_all.txt" > $OUTPUT_PREFIX"_filtered.txt"

        #grep ^CNV -v $OUTPUT_PREFIX"_all.txt" | awk '{if ( ($16>='$LOG_UP' || $16<='$LOG_DOWN' || $13>='$BF' ) && ($8>='$NB_EXON') && ($4>='$NB_COMP') ) {print $0}}' >> $OUTPUT_PREFIX"_filtered.txt"
        grep ^CNV -v $OUTPUT_PREFIX"_all.txt" | awk '{if ( ($16>='$LOG_UP' || $16<='$LOG_DOWN') && ($13>='$BF' ) && ($8>='$NB_EXON') && ($4>='$NB_COMP') ) {print $0}}' >> $OUTPUT_PREFIX"_filtered.txt"
        (($VERBOSE)) && head $OUTPUT_PREFIX"_filtered.txt"
        (($VERBOSE)) && echo "..." 
        NB_CNV_FILTERED=$(grep ^CNV -v $OUTPUT_PREFIX"_filtered.txt" | wc -l)
        (($VERBOSE)) && echo "#[INFO] DECoN makeCNVcalls NB CNV filtered: $NB_CNV_FILTERED"

        #grep ^CNV -v $OUTPUT_PREFIX"_all.txt" | awk '{if ($16>=1.2 || $16<=0.833333333) {print $1}}'



        # Convert to BED
        (($VERBOSE)) && echo "#[INFO] Concert to BED" 
        #echo -n "#" > $OUTPUT_PREFIX"_all.bed"
        #cat $OUTPUT_PREFIX"_all.txt" | awk '{print $11 "\t" $10 "\t" $9 "\t" $17 "\t" $13 "\t" "+" "\t" $2 }' >> $OUTPUT_PREFIX"_all.bed"
        cat $OUTPUT_PREFIX"_all.txt" | sed 1d | awk '{split($2,file_name,"."); sample_name=file_name[1]} {print $11 "\t" $10 "\t" $9 "\t" $17 "\t" $13 "\t" "+" "\t" sample_name }' > $OUTPUT_PREFIX"_all.bed"
        (($VERBOSE)) && echo "#[INFO] DECoN makeCNVcalls All CNV BED" 
        (($VERBOSE)) && head $OUTPUT_PREFIX"_all.bed"
        (($VERBOSE)) && echo "..." 
        NB_CNV_BED=$(grep ^CNV -v $OUTPUT_PREFIX"_all.bed" | wc -l)
        (($VERBOSE)) && echo "#[INFO] DECoN makeCNVcalls NB CNV (BED): $NB_CNV_BED"

        cat $OUTPUT_PREFIX"_filtered.txt" | sed 1d | awk '{split($2,file_name,"."); sample_name=file_name[1]} {print $11 "\t" $10 "\t" $9 "\t" $17 "\t" $13 "\t" "+" "\t" sample_name }' > $OUTPUT_PREFIX"_filtered.bed"
        (($VERBOSE)) && echo "#[INFO] DECoN makeCNVcalls Filtered CNV BED" 
        (($VERBOSE)) && head $OUTPUT_PREFIX"_filtered.bed"
        (($VERBOSE)) && echo "..." 
        NB_CNV_BED_filtered=$(grep ^CNV -v $OUTPUT_PREFIX"_filtered.bed" | wc -l)
        (($VERBOSE)) && echo "#[INFO] DECoN makeCNVcalls NB CNV filtered (BED): $NB_CNV_BED_filtered"


        # SAMPLE GENE List 
        cat $OUTPUT_PREFIX"_filtered.txt" | head -n1 | awk '{split($2,file_name,"."); sample_name=file_name[1]} {print sample_name "\t" $17 "\t" $7  "\t" $13 "\t" $16 }'  > $OUTPUT_PREFIX"_filtered.tsv"
        cat $OUTPUT_PREFIX"_filtered.txt" | sed 1d | awk '{split($2,file_name,"."); sample_name=file_name[1]} {print sample_name "\t" $17 "\t" $7  "\t" $13 "\t" $16 }' | sort -k1.1 -k2.1 >> $OUTPUT_PREFIX"_filtered.tsv"
        cat $OUTPUT_PREFIX"_filtered.tsv"



    fi;


}

function DECON_OUTPUT {
    echo "#[INFO] DECoN Output" 
    if [ -e $OUTPUT_PREFIX"_filtered.tsv" ]; then
        NB_CNV_TSV_filtered=$(grep ^Sample -v $OUTPUT_PREFIX"_filtered.tsv" | wc -l)
        (($VERBOSE)) && echo "#[INFO] DECoN NB CNV filtered: $NB_CNV_TSV_filtered"
        cat $OUTPUT_PREFIX"_filtered.tsv"
    else
        echo "#[INFO] NO Output" 
    fi;
}


function DECON {

    # RUN
    echo "#[INFO] DECoN Init" && DECON_INIT
    echo "#[INFO] DECoN Analysis..." && DECON_ANALYSIS 1>>$LOG 2>>$LOG
    DECON_OUTPUT

}




function DECON_PARAMS {

    # VERBOSE DEBUG
    VERBOSE=1


    # Ressources
    FASTA=/STARK/databases/genomes/current/hg19.fa 
    #refGene=/STARK/data/DECoN/datasets/refGene.hg19.ok.bed
    #refGene=/STARK/data/DECoN/refGene.hg19.ok.bed
    refGene=/STARK/databases/refGene/current/refGene.hg19.bed

    DECoN_main_results=/STARK/data/DECoN/results
    DECoN_sub_results=
    #DECoN_bin_folder=/STARK/tools/decon/1.0.2/DECoN-1.0.2/Linux
    #DECoN_bin_folder=/STARK/tools/decon/current
    DECoN_bin_folder=$(dirname $(whereis makeCNVcalls.R | cut -d' ' -f2))
    #Rscript=/STARK/tools/R/R-3.1.2/bin/Rscript
    Rscript=$(whereis Rscript | cut -d' ' -f2)
    #echo "Rscript: $Rscript"

    # R init
    #$Rscript $DECoN_bin_folder/sessionInfo.R --bootstrap-packrat 

    # Data
    design_ext=bed
    BED_ADD_REFGENE=1
    bam_find_min_depth=1
    bam_find_max_depth=2
    bed_find_min_depth=1
    bed_find_max_depth=1

    bam_exclude_pattern=" validation recal"
    #bam_exclude_pattern="  recal"

    # PARAMS default
    MINCORR=.98
    MINCOV=100
    transProb=0.01
    LOG_UP=1.2
    LOG_DOWN=$(echo "scale=6; 1/$LOG_UP" | bc)
    BF=10
    NB_EXON=1
    NB_COMP=1 #$(echo "$NB_BAM/2" | bc)

}



function DECON_DATASET {

    INPUT_FOLDER=$1
    #find $INPUT_FOLDER -mindepth 1 -maxdepth 1 -name *bed
    #find $INPUT_FOLDER -mindepth 3 -maxdepth 3 -name *bam
        
    # INIT
    NAME=$(basename $INPUT_FOLDER)
    (($VERBOSE)) && echo "#[INFO] ###### '$NAME' [$INPUT_FOLDER]"
    
    # PARAMS    
    echo "#[INFO] DECoN Parameters" && DECON_PARAMS
    DECoN_sub_results=datasets
    bam_find_min_depth=1
    bam_find_max_depth=2
    bed_find_min_depth=1
    bed_find_max_depth=1

    # LAUNCH
    echo "#[INFO] DECoN Launch" && DECON



}



function DECON_RUN {

    INPUT_FOLDER=$1
    #find $INPUT_FOLDER -mindepth 1 -maxdepth 1 -name *bed
    #find $INPUT_FOLDER -mindepth 3 -maxdepth 3 -name *bam
        
    # INIT
    NAME=$(basename $INPUT_FOLDER)
    (($VERBOSE)) && echo "#[INFO] ###### '$NAME' [$INPUT_FOLDER]"
    
    # PARAMS    
    echo "#[INFO] DECoN Parameters" && DECON_PARAMS
    DECoN_sub_results=custom
    bam_find_min_depth=3
    bam_find_max_depth=3
    bed_find_min_depth=3
    bed_find_max_depth=3

    # LAUNCH
    echo "#[INFO] DECoN Launch" && DECON



}


function DECON_PROJECT {

    INPUT_FOLDER=$1
    #find $INPUT_FOLDER -mindepth 1 -maxdepth 1 -name *bed
    #find $INPUT_FOLDER -mindepth 3 -maxdepth 3 -name *bam
        
    # INIT
    NAME=$(basename $INPUT_FOLDER)
    (($VERBOSE)) && echo "#[INFO] ###### '$NAME' [$INPUT_FOLDER]"
    
    # PARAMS    
    echo "#[INFO] DECoN Parameters" && DECON_PARAMS
    DECoN_sub_results=by_project
    bam_find_min_depth=4
    bam_find_max_depth=4
    bed_find_min_depth=4
    bed_find_max_depth=4

    # LAUNCH
    echo "#[INFO] DECoN Launch" && DECON


}


function DECON_GROUP {

    INPUT_FOLDER=$1
    #find $INPUT_FOLDER -mindepth 1 -maxdepth 1 -name *bed
    #find $INPUT_FOLDER -mindepth 3 -maxdepth 3 -name *bam
        
    # INIT
    NAME=$(basename $INPUT_FOLDER)
    (($VERBOSE)) && echo "#[INFO] ###### '$NAME' [$INPUT_FOLDER]"
    
    # PARAMS    
    echo "#[INFO] DECoN Parameters" && DECON_PARAMS
    DECoN_sub_results=by_group
    bam_find_min_depth=5
    bam_find_max_depth=5
    bed_find_min_depth=5
    bed_find_max_depth=5

    # LAUNCH
    echo "#[INFO] DECoN Launch" && DECON


}


function DECON_REPOSITORY {

    INPUT_FOLDER=$1
    #find $INPUT_FOLDER -mindepth 1 -maxdepth 1 -name *bed
    #find $INPUT_FOLDER -mindepth 3 -maxdepth 3 -name *bam
        
    # INIT
    NAME=$(basename $INPUT_FOLDER)
    (($VERBOSE)) && echo "#[INFO] ###### '$NAME' [$INPUT_FOLDER]"
    
    # PARAMS    
    echo "#[INFO] DECoN Parameters" && DECON_PARAMS
    DECoN_sub_results=repository
    bam_find_min_depth=6
    bam_find_max_depth=6
    bed_find_min_depth=6
    bed_find_max_depth=6

    # LAUNCH
    echo "#[INFO] DECoN Launch" && DECON


}





### RUN

#echo "#[INFO] Run '/STARK/data/DECoN/Repository/DIAG/BBS/BBS_1'"
#DECON_RUN /STARK/data/DECoN/Repository/DIAG/BBS/BBS_1

#RUN_FOLDER=/STARK/data/DECoN/Repository/DIAG/BBS/BBS_1
#echo "#[INFO] Run '$RUN_FOLDER'"
#DECON_RUN $RUN_FOLDER

# RUN_FOLDER=/STARK/data/DECoN/Repository/HUSHEMATO/TSOMYELOID/200814_M01658_0462_000000000-J7387
# echo "#[INFO] Run '$RUN_FOLDER'"
# DECON_RUN $RUN_FOLDER


## Datasets list by folders
if true; then

    VERBOSE=1

    # DATASETS
    if true; then
        (($VERBOSE)) && echo "#[INFO] Analyses by DATASETS"
        for INPUT_FOLDER in $(find /STARK/data/DECoN/DATASETS/ -mindepth 1 -maxdepth 1 -type d); do 
            time DECON_RUN $INPUT_FOLDER
        done;
    fi;


fi;


## datasets as run in repository
if false; then

    VERBOSE=1

    # RUNS
    if true; then
        (($VERBOSE)) && echo "#[INFO] Analyses by RUN"
        for INPUT_FOLDER in $(find /STARK/data/DECoN/Repository/ -mindepth 3 -maxdepth 3 -type d); do 
            time DECON_RUN $INPUT_FOLDER
        done;
    fi;

    # PROJECTS
    if true; then
        (($VERBOSE)) && echo "#[INFO] Analyses by PROJECTS"
        for INPUT_FOLDER in $(find /STARK/data/DECoN/Repository/ -mindepth 2 -maxdepth 2 -type d); do
            time DECON_PROJECT $INPUT_FOLDER
        done;
    fi;

    # GROUP
    if true; then
        (($VERBOSE)) && echo "#[INFO] Analyses by GROUP"
       for INPUT_FOLDER in $(find /STARK/data/DECoN/Repository/ -mindepth 1 -maxdepth 1 -type d); do
            time DECON_GROUP $INPUT_FOLDER
        done;
    fi;

    # REPOSITORY
    if true; then
        (($VERBOSE)) && echo "#[INFO] Analyses by REPOSITORY"
        for INPUT_FOLDER in $(find /STARK/data/DECoN/Repository/ -mindepth 0 -maxdepth 0 -type d); do
            time DECON_REPOSITORY $INPUT_FOLDER
        done;
    fi;


fi;



if false; then



    if false; then

        # ANALYSIS of a RUN

        #INPUT_FOLDER=/STARK/data/DECoN/Repository/DIAG/MYOPATHIE/MYOPATHIE_160928_NB551027_0056_AHGNJ5AFXX
        #INPUT_FOLDER=/STARK/data/DECoN/Repository/HUSHEMATO/TSOMYELOID/200814_M01658_0462_000000000-J7387
        INPUT_FOLDER=/STARK/data/DECoN/Repository/DIAG/BBS/BBS_1
        #find $INPUT_FOLDER -mindepth 1 -maxdepth 1 -name *bed
        #find $INPUT_FOLDER -mindepth 3 -maxdepth 3 -name *bam
            
        # INIT
        NAME=$(basename $INPUT_FOLDER)
        (($VERBOSE)) && echo "#[INFO] ###### '$NAME' [$INPUT_FOLDER]"
        
        # PARAMS    
        echo "#[INFO] DECoN Parameters" && DECON_PARAMS
        DECoN_sub_results=custom2
        bam_find_min_depth=3
        bam_find_max_depth=3
        bed_find_min_depth=3
        bed_find_max_depth=3

        # LAUNCH
            echo "#[INFO] DECoN Launch" && time DECON


    fi;

    # echo "#[INFO] Run '$RUN_FOLDER'"
    # DECON_RUN $RUN_FOLDER


    if true; then

        # REPOSITORY ANALYSIS by RUN

        # rm -rf $(find /STARK/data/DECoN/Repository  -mindepth 4 -maxdepth 4 -type d -name DECoN)

        for INPUT_FOLDER in $(find /STARK/data/DECoN/Repository/ -mindepth 3 -maxdepth 3 -type d); do
            
            # INIT
            NAME=$(basename $INPUT_FOLDER)
            (($VERBOSE)) && echo "#[INFO] ###### '$NAME' [$INPUT_FOLDER]"
            
            # PARAMS
            echo "#[INFO] DECoN Parameters" && DECON_PARAMS
            DECoN_sub_results=by_run
            bam_find_min_depth=3
            bam_find_max_depth=3
            bed_find_min_depth=3
            bed_find_max_depth=3

            # LAUNCH
            echo "#[INFO] DECoN Launch" && time DECON

        done;


    fi;


    if true; then

        # REPOSITORY ANALYSIS by PROJECT

        # rm -rf $(find /STARK/data/DECoN/Repository  -mindepth 4 -maxdepth 4 -type d -name DECoN)

        for INPUT_FOLDER in $(find /STARK/data/DECoN/Repository/ -mindepth 2 -maxdepth 2 -type d); do
            
            # INIT
            NAME=$(basename $INPUT_FOLDER)
            (($VERBOSE)) && echo "#[INFO] ###### '$NAME' [$INPUT_FOLDER]"
            
            # PARAMS
            DECoN_sub_results=by_project
            bam_find_min_depth=4
            bam_find_max_depth=4
            bed_find_min_depth=4
            bed_find_max_depth=4

            # LAUNCH
            echo "#[INFO] DECoN Launch" && time DECON

        done;


    fi;



    if true; then

        # REPOSITORY ANALYSIS by GROUP

        # rm -rf $(find /STARK/data/DECoN/Repository  -mindepth 4 -maxdepth 4 -type d -name DECoN)

        for INPUT_FOLDER in $(find /STARK/data/DECoN/Repository/ -mindepth 1 -maxdepth 1 -type d); do
            
            # INIT
            NAME=$(basename $INPUT_FOLDER)
            (($VERBOSE)) && echo "#[INFO] ###### '$NAME' [$INPUT_FOLDER]"
            
            # PARAMS
            DECoN_sub_results=by_group
            bam_find_min_depth=5
            bam_find_max_depth=5
            bed_find_min_depth=5
            bed_find_max_depth=5

            # LAUNCH
            echo "#[INFO] DECoN Launch" && time DECON

        done;


    fi;



    if true; then

        # REPOSITORY ANALYSIS FULL

        # rm -rf $(find /STARK/data/DECoN/Repository  -mindepth 4 -maxdepth 4 -type d -name DECoN)

        for INPUT_FOLDER in $(find /STARK/data/DECoN/Repository/ -mindepth 0 -maxdepth 0 -type d); do
            
            # INIT
            NAME=$(basename $INPUT_FOLDER)
            (($VERBOSE)) && echo "#[INFO] ###### '$NAME' [$INPUT_FOLDER]"
            
            # PARAMS
            DECoN_sub_results=all
            bam_find_min_depth=6
            bam_find_max_depth=6
            bed_find_min_depth=6
            bed_find_max_depth=6

            # LAUNCH
            echo "#[INFO] DECoN Launch" && time DECON

        done;


    fi;


fi;




if false; then


    if true; then

        # INPUT
        INPUT_FOLDER=/STARK/data/DECoN/datasets/RUN_TEST_TAG_LISTENER_MIN_TEST

        # PARAMS
        #MINCORR=.80
        #MINCOV=100
        #transProb=0.05

        # RUN
        (($VERBOSE)) && echo "#[INFO] ###### RUN $INPUT_FOLDER"
        echo "#[INFO] DECoN Init..." && DECON_INIT
        echo "#[INFO] DECoN Analysis..." && DECON 1>$LOG 2>$ERR
        DECON_OUTPUT

    fi

    if true; then

        # INPUT
        INPUT_FOLDER=/STARK/data/DECoN/datasets/191001_M01656_0336_000000000-CLH59

        # PARAMS
        bed_find_max_depth=2

        # RUN
        (($VERBOSE)) && echo "#[INFO] ###### RUN $INPUT_FOLDER"
        echo "#[INFO] DECoN Init..." && DECON_INIT
        echo "#[INFO] DECoN Analysis..." && DECON 1>$LOG 2>$ERR
        DECON_OUTPUT


    fi

    if true; then

        # INPUT
        INPUT_FOLDER=/STARK/data/DECoN/datasets/200720_NB551027_0749_AH2WL7AFX2

        # PARAMS
        bam_find_max_depth=3
        bed_find_max_depth=2

        # RUN
        (($VERBOSE)) && echo "#[INFO] ###### RUN $INPUT_FOLDER"
        echo "#[INFO] DECoN Init..." && DECON_INIT
        echo "#[INFO] DECoN Analysis..." && DECON 1>$LOG 2>$ERR
        DECON_OUTPUT


    fi;


    ### DATASETS

    if true; then

        # INPUT
        INPUT_FOLDER=/STARK/data/DECoN/datasets/BBS_1

        # RUN
        (($VERBOSE)) && echo "#[INFO] ###### RUN $INPUT_FOLDER"
        echo "#[INFO] DECoN Init..." && DECON_INIT
        echo "#[INFO] DECoN Analysis..." && DECON 1>$LOG 2>$ERR
        DECON_OUTPUT

    fi;


    if true; then

        # INPUT
        INPUT_FOLDER=/STARK/data/DECoN/datasets/BBS_v1_Run4_062012

        # RUN
        (($VERBOSE)) && echo "#[INFO] ###### RUN $INPUT_FOLDER"
        echo "#[INFO] DECoN Init..." && DECON_INIT
        echo "#[INFO] DECoN Analysis..." && DECON 1>$LOG 2>$ERR
        DECON_OUTPUT

    fi;


    if true; then

        # INPUT
        INPUT_FOLDER=/STARK/data/DECoN/datasets/CPS_AF206

        # RUN
        (($VERBOSE)) && echo "#[INFO] ###### RUN $INPUT_FOLDER"
        echo "#[INFO] DECoN Init..." && DECON_INIT
        echo "#[INFO] DECoN Analysis..." && DECON 1>$LOG 2>$ERR
        DECON_OUTPUT

    fi;



    if true; then

        # INPUT
        INPUT_FOLDER=/STARK/data/DECoN/datasets/BBS_*

        # RUN
        (($VERBOSE)) && echo "#[INFO] ###### RUN $INPUT_FOLDER"
        echo "#[INFO] DECoN Init..." && DECON_INIT
        echo "#[INFO] DECoN Analysis..." && DECON 1>$LOG 2>$ERR
        DECON_OUTPUT

    fi;


    if true; then

        # INPUT
        INPUT_FOLDER=/STARK/data/DECoN/datasets/*

        # PARAMS
        bam_find_max_depth=3
        bed_find_max_depth=2

        # RUN
        (($VERBOSE)) && echo "#[INFO] ###### RUN $INPUT_FOLDER"
        echo "#[INFO] DECoN Init..." && DECON_INIT
        echo "#[INFO] DECoN Analysis..." && DECON 1>$LOG 2>$ERR
        DECON_OUTPUT

    fi;


fi;

# Visu
#Rscript runShiny.R --Rdata $Rdata
