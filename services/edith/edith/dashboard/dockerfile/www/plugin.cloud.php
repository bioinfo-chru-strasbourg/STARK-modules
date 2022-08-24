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

// echo '
// <meta charset="UTF-8">
// <meta http-equiv="X-UA-Compatible" content="IE=edge">
// <meta name="generator" content="Mobirise v4.10.3, mobirise.com">
// <meta name="viewport" content="width=device-width, initial-scale=1, minimum-scale=1">
// <link rel="shortcut icon" href="assets/favicon.ico" type="image/x-icon">
// <meta name="description" content="EDITH - '.$APP_SECTION.'">
// <title>EDITH - '.$APP_SECTION.'</title>
// <link rel="stylesheet" href="assets/web/assets/mobirise-icons/mobirise-icons.css">
// <link rel="stylesheet" href="assets/tether/tether.min.css">
// <link rel="stylesheet" href="assets/bootstrap/css/bootstrap.min.css">
// <link rel="stylesheet" href="assets/bootstrap/css/bootstrap-grid.min.css">
// <link rel="stylesheet" href="assets/bootstrap/css/bootstrap-reboot.min.css">
// <link rel="stylesheet" href="assets/dropdown/css/style.css">
// <link rel="stylesheet" href="assets/as-pie-progress/css/progress.min.css">
// <link rel="stylesheet" href="assets/datatables/data-tables.bootstrap4.min.css">
// <link rel="stylesheet" href="assets/theme/css/style.css">
// <link rel="stylesheet" href="assets/gallery/style.css">
// <link rel="stylesheet" href="assets/mobirise/css/mbr-additional.css" type="text/css">
// <style>
// 	.div-wrapper {
// 		overflow: auto;
// 		max-height:400px;
// 		}
// </style>
// ';



# MAIN URL
############


# CLOUD URL
if ($uri_cloud == "") {
	if ($_ENV["URI_CLOUD"] != "") {
		$uri_cloud=$_ENV["URI_CLOUD"];
	} elseif ($modules_obj_array["cloud"]->{"submodules"}->{"cloud"}->{"services"}->{"cloud"}->{"href"}!="") {
		$uri_cloud=$modules_obj_array["cloud"]->{"submodules"}->{"cloud"}->{"services"}->{"cloud"}->{"href"};
	} elseif ($modules_obj_array["cloud"]->{"services"}->{"cloud"}->{"href"}!="") {
		$uri_cloud=$modules_obj_array["cloud"]->{"services"}->{"cloud"}->{"href"};
	} else {
		$uri_cloud="";
	};
};


# PATH
$PATH=$SAMPLE_PATH[0];

# EMBED FINAL URL
####################

# FINAL URL
$FINAL_URL="$uri_cloud$PATH";

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
#echo "<embed style='width:100%;height:100%' src='$FINAL_URL'>";
header("Location: $FINAL_URL");

?>
