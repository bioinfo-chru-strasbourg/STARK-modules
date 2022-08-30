<?php


#############
### INFOS ###
#############

$APP_SECTION="Plugin JARVIS";



################
### INCLUDES ###
################


### CONFIG
############

include "config.inc.php";


### FUNCTIONS
###############

include "functions.inc.php";


### HEADER
############

#include "header.inc.php";


### VARIABLES
###############



### DATA
##########

#include "global_statistics.inc.php";
#include "activity.inc.php";



$DEBUG=$_REQUEST["DEBUG"];


# PARAMETERS
##############

$SAMPLE_PATH=$_GET["PATH"];
$SEARCH_LEVEL=$_GET["SEARCH_LEVEL"];
$SEARCH_EXT=$_GET["SEARCH_EXT"];
$PROCESS=$_GET["PROCESS"];
#$DEBUG=0;


# PROCESS
if ($PROCESS=="") {
	$PROCESS=0;
};

# SEARCH_EXT
if ($SEARCH_EXT=="") {
	$SEARCH_EXT="vcf.gz,vcf,gvcf,gvcf.gz,bcf,bcf.gz";
};



# FUNCTIONS
#############

function array_file_uniq ($FILES=array()) {
	# INPUT:
	# Array of files path
	# Array ([0]=>"/path/to/my/file1",[1]=>"/path/to/my/file2",...,[n]=>"/path/to/my/filen")
	# OUTPUT:
	# Array with uniq file base on filename
	# Array ([file1]=>"/path/to/my/file1",[file2]=>"/path/to/my/file2",...,[filen]=>"/path/to/my/filen")

	$return=array();

	foreach ($FILES as $FILE_key => $FILE_PATH) {
		$FILE_NAME=end(explode( "/", $FILE_PATH ));
		$return[$FILE_NAME]=$FILE_PATH;
	}
	return $return;

}


# MAIN URL
############

# JARVIS URL
$JARVIS_URL=$modules_obj_array["VARIANTBROWSER"]->{"services"}->{"VISION"}->{"href"};

# DATA URL
$DATA_URL=$modules_obj_array["STARK"]->{"services"}->{"DAS"}->{"href_inner"};



# VARIABLES
#############

# VCF_FILES
$VCF_FILES = "";



# SEARCH FOR TRACKS
#####################

#DEV
if ($DEBUG) {
	print_r($SAMPLE_PATH);
};

# TRACKS JSON ARRAY
$VCF_FILES_ARRAY=array();
$VCF_FILENAMES_ARRAY=array();

$VCF_FILES_KEY=0;

#echo $SEARCH_LEVEL;
$PATH_LEVEL="";
for ($i=0; $i<$SEARCH_LEVEL+0; $i++) {
	$PATH_LEVEL.="/*";
}


foreach ($SAMPLE_PATH as $key_sample => $ONE_SAMPLE_PATH) {

	# SAMPLE NAME
	$ONE_SAMPLE_NAME=end(explode( "/", $ONE_SAMPLE_PATH ));

	# Variants
	############

	# Search files
	$vcf=array();
	$vcf_root=glob("$ONE_SAMPLE_PATH$PATH_LEVEL/*{".$SEARCH_EXT."}",GLOB_BRACE);
	if ($CHECK_SUBFOLDER_DATA) {
		$vcf_data=glob("$ONE_SAMPLE_PATH$PATH_LEVEL/$RESULTS_SUBFOLDER_DATA/*{".$SEARCH_EXT."}",GLOB_BRACE);
	} else {
		$vcf_data=array();
	};
	$vcf=array_merge ($vcf_root,$vcf_data);

	#$vcf_merge=array_merge ($vcf_root,$vcf_data);
	#$vcf=array_file_uniq($vcf_merge);

	# Create tracks
	foreach ($vcf as $key_vcf => $VCF) {
		$VCF_FILES_KEY++;
		$VCF_NAME=end(explode( "/", $VCF ));
		$VCF_URL="$DATA_URL/".$VCF;
		$VCF_FILES_ARRAY[$VCF_FILES_KEY]="vcf_files[]=$VCF_URL";
		$VCF_FILENAMES_ARRAY[$VCF_FILES_KEY]="vcf_filenames[]=$VCF_NAME";
	}

};


# IMPLODE TRACKS
$VCF_FILES = implode ( "&" ,$VCF_FILES_ARRAY );
$VCF_FILENAMES = implode ( "&" ,$VCF_FILENAMES_ARRAY );


# DEV
if ($DEBUG) {
	echo "<br>";
	echo "<pre>";
	echo $VCF_FILES;
	echo "</pre>";
};



# EMBED FINAL URL
####################

# FINAL URL
$FINAL_URL="$JARVIS_URL?process=$PROCESS&$VCF_FILES";

# DEV
if ($DEBUG) {
	echo "<br>";
	echo "<pre>";
	echo $FINAL_URL;
	echo "</pre>";
};


#exit(0);
#echo $FINAL_URL;
# HEADER
header("Location: $FINAL_URL");
#exit;
#echo "<script>window.location.href = '".$FINAL_URL."';</script>";
#echo "<embed style='width:100%;height:100%' src='$FINAL_URL'>";

?>
