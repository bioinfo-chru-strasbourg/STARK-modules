<?php



#############
### INFOS ###
#############

$APP_SECTION="Plugin IGV";



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



# PARAMETERS
##############

$SAMPLE_PATH=$_GET["PATH"];
$CHECK_SUBFOLDER_DATA=$_GET["FULL_FILES"];
$SEARCH_LEVEL=$_GET["SEARCH_LEVEL"];
$DEBUG=$_REQUEST["DEBUG"];

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


### HEADER
############

echo '
<meta charset="UTF-8">
<meta http-equiv="X-UA-Compatible" content="IE=edge">
<meta name="generator" content="Mobirise v4.10.3, mobirise.com">
<meta name="viewport" content="width=device-width, initial-scale=1, minimum-scale=1">
<link rel="shortcut icon" href="assets/favicon.ico" type="image/x-icon">
<meta name="description" content="EDITH - '.$APP_SECTION.'">
<title>EDITH - '.$APP_SECTION.'</title>
<link rel="stylesheet" href="assets/web/assets/mobirise-icons/mobirise-icons.css">
<link rel="stylesheet" href="assets/tether/tether.min.css">
<link rel="stylesheet" href="assets/bootstrap/css/bootstrap.min.css">
<link rel="stylesheet" href="assets/bootstrap/css/bootstrap-grid.min.css">
<link rel="stylesheet" href="assets/bootstrap/css/bootstrap-reboot.min.css">
<link rel="stylesheet" href="assets/dropdown/css/style.css">
<link rel="stylesheet" href="assets/as-pie-progress/css/progress.min.css">
<link rel="stylesheet" href="assets/datatables/data-tables.bootstrap4.min.css">
<link rel="stylesheet" href="assets/theme/css/style.css">
<link rel="stylesheet" href="assets/gallery/style.css">
<link rel="stylesheet" href="assets/mobirise/css/mbr-additional.css" type="text/css">
<style>
	.div-wrapper {
		overflow: auto;
		max-height:400px;
		}
</style>
';



# MAIN URL
############

# IGV URL
if ($uri_igv == "") {
	if ($_ENV["URI_IGV"] != "") {
		$uri_igv=$_ENV["URI_IGV"];
	} elseif ($modules_obj_array["genomebrowser"]->{"services"}->{"igv"}->{"href"}!="") {
		$uri_igv=$modules_obj_array["genomebrowser"]->{"services"}->{"igv"}->{"href"};
	} else {
		$uri_igv="";
	};
};


# DATA URL

if ($uri_das == "") {
	if ($_ENV["URI_DAS"] != "") {
		$uri_das=$_ENV["URI_DAS"];
	} elseif ($modules_obj_array["stark"]->{"services"}->{"das"}->{"href"}!="") {
		$uri_das=$modules_obj_array["stark"]->{"services"}->{"das"}->{"href"};
	} else {
		$uri_das="";
	};
};



#DEV
if ($DEBUG) {
	echo "<pre>$uri_igv</pre>";
	echo "<pre>$uri_das</pre>";
};


# VARIABLES
#############

# DATE
$DATE=date('omd-his');

# JSON ID
$IGV_JSON_ID=$DATE.".".rand ( 1 , 1000000 );

# IGV_JSON_FILENAME
$IGV_JSON_FILENAME="$IGV_JSON_ID.json";

# JSON FILE
$IGV_JSON_FILE="".$folder_services."/GENOMEBROWSER/IGV/$IGV_JSON_FILENAME";

if ($DEBUG) {
	echo "<br>$IGV_JSON_FILE";
}


# JSON
$JSON = "";



# SEARCH FOR TRACKS
#####################

#DEV
if ($DEBUG) {
	print_r($SAMPLE_PATH);
};

# TRACKS JSON ARRAY
$TRACKS_JSON_ARRAY=array();

$TRACKS_JSON_KEY=0;

#echo $SEARCH_LEVEL;
$PATH_LEVEL="";
for ($i=0; $i<$SEARCH_LEVEL+0; $i++) {
	$PATH_LEVEL.="/*";
}

foreach ($SAMPLE_PATH as $key_sample => $ONE_SAMPLE_PATH) {

	# SAMPLE NAME
	$ONE_SAMPLE_NAME=end(explode( "/", $ONE_SAMPLE_PATH ));

	# Designs
	############

	# Search files
	$design=array();
	#$design_root=glob("$ONE_SAMPLE_PATH$PATH_LEVEL/*{bed,genes}",GLOB_BRACE);
	$design_root=glob("$ONE_SAMPLE_PATH$PATH_LEVEL/*{bed}",GLOB_BRACE);
	if ($CHECK_SUBFOLDER_DATA) {
		#$design_data=glob("$ONE_SAMPLE_PATH$PATH_LEVEL/$RESULTS_SUBFOLDER_DATA/*{bed,genes}",GLOB_BRACE);
		$design_data=glob("$ONE_SAMPLE_PATH$PATH_LEVEL/$RESULTS_SUBFOLDER_DATA/*{bed}",GLOB_BRACE);
	} else {
		$design_data=array();
	};
	$design_merge=array_merge ($design_root,$design_data);
	$design_merge_clean=array();

	# Test BED format
	foreach ($design_merge as $design_key=>$design_file) {
		$is_a_bed=1;
		foreach (file($design_file) as $line_nb=>$line_content) {
			$line_content_split=(explode( "\t", $line_content ));
			if ($line_content_split[0]!="browser" &&
				$line_content_split[0]!="track" &&
				substr(trim($line_content_split[0]),0,1)!="#" &&
				count($line_content_split)<3
			) {
				$is_a_bed=0;
				break;
			}
		};
		if ($is_a_bed) {
			$design_merge_clean[$design_key]=$design_file;
		};
	}

	$design=array_file_uniq($design_merge_clean);

	# Create tracks
	foreach ($design as $key_design => $DESIGN) {

		$TRACKS_JSON_KEY++;

		$DESIGN_NAME=end(explode( "/", $DESIGN ));
		$DESIGN_FORMAT=pathinfo($DESIGN,PATHINFO_EXTENSION);
		$TRACKS_URL=$uri_das."/".$DESIGN;

		$ALIGNMENT_JSON = ' {
			  "url": "'.$TRACKS_URL.'",
			  "indexURL": null,
			  "name": "'.$DESIGN_NAME.'",
			  "sourceType": "file",
			  "format": "'.$design_format.'",
			  "order": '.$TRACKS_JSON_KEY.',
			  "height": "40",
			  "type": "annotation"
			}
		';
		$TRACKS_JSON_ARRAY[$TRACKS_JSON_KEY]=$ALIGNMENT_JSON;

	}

	# Variants
	############

	# Search files
	$vcf=array();
	$vcf_root=glob("$ONE_SAMPLE_PATH$PATH_LEVEL/*{vcf.gz,gvcf.gz,bcf.gz}",GLOB_BRACE);
	if ($CHECK_SUBFOLDER_DATA) {
		$vcf_data=glob("$ONE_SAMPLE_PATH$PATH_LEVEL/$RESULTS_SUBFOLDER_DATA/*{vcf.gz,gvcf.gz,bcf.gz}",GLOB_BRACE);
	} else {
		$vcf_data=array();
	};
	$vcf_merge=array_merge ($vcf_root,$vcf_data);
	$vcf=array_file_uniq($vcf_merge);

	# Create tracks
	foreach ($vcf as $key_vcf => $VCF) {

		$TRACKS_JSON_KEY++;

		$VCF_NAME=end(explode( "/", $VCF ));
		$VCF_FORMAT="vcf";
		$TRACKS_URL=$uri_das."/".$VCF;

		$ALIGNMENT_JSON = ' {
			  "url": "'.$TRACKS_URL.'",
			  "indexURL": null,
			  "name": "'.$VCF_NAME.'",
			  "sourceType": "file",
			  "format": "'.$VCF_FORMAT.'",
			  "order": '.$TRACKS_JSON_KEY.',
			  "height": "40",
			  "type": "annotation"
			}
		';
		$TRACKS_JSON_ARRAY[$TRACKS_JSON_KEY]=$ALIGNMENT_JSON;

	}

	# Alignment
	##############

	# Search files
	$alignment=array();
	$alignment_root=glob("$ONE_SAMPLE_PATH$PATH_LEVEL/*{bam,cram}",GLOB_BRACE);
	if ($CHECK_SUBFOLDER_DATA) {
		$alignment_data=glob("$ONE_SAMPLE_PATH$PATH_LEVEL/$RESULTS_SUBFOLDER_DATA/*{bam,cram}",GLOB_BRACE);
	} else {
		$alignment_data=array();
	};
	$alignment_merge=array_merge ($alignment_root,$alignment_data );
	$alignment=array_file_uniq($alignment_merge);

	# Create tracks
	foreach ($alignment as $key_alignment => $ALIGNMENT) {

		$TRACKS_JSON_KEY++;

		$ALIGNMENT_NAME=end(explode( "/", $ALIGNMENT ));
		$ALIGNMENT_FORMAT=pathinfo($ALIGNMENT,PATHINFO_EXTENSION);
		$TRACKS_URL=$uri_das."/".$ALIGNMENT;
		$ALIGNMENT_INDEX="";

		# ALIGNMENT FILE
		if ($ALIGNMENT_FORMAT=="bam") {
			$ALIGNMENT_INDEX=glob( "$ALIGNMENT.bai")[0];
		} else if ($ALIGNMENT_FORMAT=="cram") {
			$ALIGNMENT_INDEX=glob( "$ALIGNMENT.crai")[0];
		};

		# INDEX
		if ($ALIGNMENT_INDEX == "") {
			$TRACKS_INDEXURL="null";
		} else {
			$TRACKS_INDEXURL=$uri_das."/".$ALIGNMENT_INDEX."";
		};

		# ALIGNMENT JSON
		$ALIGNMENT_JSON = ' {
			  "url": "'.$TRACKS_URL.'",
			  "indexURL": "'.$TRACKS_INDEXURL.'",
			  "name": "'.$ALIGNMENT_NAME.'",
			  "sourceType": "file",
			  "format": "'.$ALIGNMENT_FORMAT.'",
			  "order": '.$TRACKS_JSON_KEY.',
			  "type": "alignment"
			}
		';

		#"autoHeight": true,

		# ADD TRACKS
		if ($TRACKS_INDEXURL!="null") {
			$TRACKS_JSON_ARRAY[$TRACKS_JSON_KEY]=$ALIGNMENT_JSON;
		};

	};


	# IMPLODE TRACKS
	$TRACKS_JSON = implode ( " , " ,$TRACKS_JSON_ARRAY );

};


# JSON
$JSON='
{
	"genome": "hg19",
	"tracks": [
		'.$TRACKS_JSON.'
	]
}
';


# DEV
if ($DEBUG) {
	echo "<br>";
	echo $IGV_JSON_FILE;
	echo "<pre>";
	echo $JSON;
	echo "</pre>";
};


# JSON WRITE
###############

# JSON OPEN
$myfile = fopen($IGV_JSON_FILE, "w") or die("Unable to open file '$IGV_JSON_FILE'!");

# JSON WRITE
fwrite($myfile, $JSON);

# JSON CLOSE
fclose($myfile);


# EMBED FINAL URL
####################

# FINAL URL
$FINAL_URL="$uri_igv?file=$uri_das/services/GENOMEBROWSER/IGV/$IGV_JSON_FILENAME";

# DEV
if ($DEBUG) {
	echo "<br>";
	echo $FINAL_URL;
	echo "<pre>";
	#echo $JSON;
	echo "</pre>";
};


# EMBED
#echo "<embed style='width:100%;height:100%' src='$FINAL_URL'>";
echo "<embed style='width:100%;height:100%' src='$FINAL_URL'>";

?>
