<?php


#############
### INFOS ###
#############

$APP_SECTION="Indexes";



################
### INCLUDES ###
################


### CONFIG
############

include "config.inc.php";


### FUNCTIONS
###############

#include "functions.inc.php";

### HEADER
############

#include "header.inc.php";


### VARIABLES
###############



### DATA
##########


### CONFIGURATION
###################



$indexes_old=60;
#$indexes_old=6;


### INDEX
$inputs_indexes_filemtime=filemtime("$folder_indexes/inputs.idx");
$now=time();
$inputs_indexes_old=($now - $inputs_indexes_filemtime);
#echo "$now - $inputs_indexes_filemtime = $inputs_indexes_old<br>";
$repositories_indexes_filemtime=filemtime("$folder_indexes/repositories.idx");
$now=time();
$repositories_indexes_old=($now - $repositories_indexes_filemtime);
#echo "$now - $repositories_indexes_filemtime = $repositories_indexes_old<br>";



if ( !is_file("$folder_indexes/inputs.idx") || ($inputs_indexes_old >= $indexes_old) || !is_file("$folder_indexes/repositories.idx") || ($repositories_indexes_old >= $indexes_old) ) {
	if (!is_file("$folder_indexes/inputs.idx.tmp") && !is_file("$folder_indexes/repositories.idx.tmp")) {

		#echo "create indexes<br>";

		# inputs filter
		$input_filter=array_keys($configuration["inputs"]["filter"]);
		if (empty($input_filter)) {
			$input_filter=array("*/*/*");
		}
		$input_filter_list=" $folder_inputs/".join(" $folder_inputs/",array_unique($input_filter));

		# repositories filter
		$repository_filter=array_keys($configuration["repositories"]["filter"]);
		if (empty($repository_filter)) {
			$repository_filter=array("*/*/*/*");
		}
		$repository_filter_list=" $folder_repositories/".join(" $folder_repositories/",array_unique($repository_filter));


		### INDEX

		# inputs
		$command="./indexes.sh '' '$input_filter_list' '0' '0' $folder_indexes/inputs.idx $folder_indexes/inputs.idx.tmp";
		$output = executeAsyncShellCommand($command);

		# inputs
		$command="./indexes.sh '' '$repository_filter_list' '0' '3' $folder_indexes/repositories.idx $folder_indexes/repositories.idx.tmp";
		$output = executeAsyncShellCommand($command);

	}
}


?>