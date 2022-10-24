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


# Index delay
$indexes_old=$configuration["parameters"]["data_indexes_delay"]+0;
$indexes_tmp_old=86400; # 1 day
#$indexes_old=6;
#echo "<br><br><br><br>indexes_old=$indexes_old";


### INDEX

# inputs - runs
$runs_indexes_filemtime=filemtime("$folder_indexes/runs.idx");
$now=time();
$runs_indexes_old=($now - $runs_indexes_filemtime);

# inputs - manifests
$inputs_indexes_filemtime=filemtime("$folder_indexes/manifests.idx");
$now=time();
$inputs_indexes_old=($now - $inputs_indexes_filemtime);

# outputs - repositories
$repositories_indexes_filemtime=filemtime("$folder_indexes/repositories.idx");
$now=time();
$repositories_indexes_old=($now - $repositories_indexes_filemtime);

# TMP

# inputs - runs
$runs_indexes_tmp_filemtime=filemtime("$folder_indexes/runs.idx.tmp");
$now=time();
$runs_indexes_tmp_old=($now - $runs_indexes_tmp_filemtime);

# inputs - manifests
$inputs_indexes_tmp_filemtime=filemtime("$folder_indexes/manifests.idx.tmp");
$now=time();
$inputs_indexes_tmp_old=($now - $inputs_indexes_tmp_filemtime);

# outputs - repositories
$repositories_indexes_tmp_filemtime=filemtime("$folder_indexes/repositories.idx.tmp");
$now=time();
$repositories_indexes_tmp_old=($now - $repositories_indexes_tmp_filemtime);

# fix if tmp too old (usually due to service hard down)
if (	($runs_indexes_tmp_old >= $indexes_tmp_old)
	||	($inputs_indexes_tmp_old >= $indexes_tmp_old)
	||	($repositories_indexes_tmp_old >= $indexes_tmp_old)
) {
	# tmp command
	$command="rm -f $folder_indexes/runs.idx.tmp $folder_indexes/manifests.idx.tmp $folder_indexes/repositories.idx.tmp";

	# inputs runs exec
	if ($indexes_tmp_old) {
		$output = executeAsyncShellCommand($command);
	} else {
		$output = shell_exec($command);
	}
}


if (   !is_file("$folder_indexes/runs.idx") || ($runs_indexes_old >= $indexes_old)
	|| !is_file("$folder_indexes/manifests.idx") || ($inputs_indexes_old >= $indexes_old)
	|| !is_file("$folder_indexes/manifests.idx") || ($repositories_indexes_old >= $indexes_old)
	) {

	if ( !is_file("$folder_indexes/runs.idx.tmp")
	&&   !is_file("$folder_indexes/manifests.idx.tmp")
	&&   !is_file("$folder_indexes/repositories.idx.tmp")
	) {

		### Runs

		# inputs runs filter
		$inputs_runs_filter=array_keys($configuration["inputs"]["runs"]["filter"]);
		if (empty($inputs_runs_filter)) {
			$inputs_runs_filter=array("*/runs/*");
		}

		# inputs runs parameters
		$inputs_runs_files_patterns=" -name SampleSheet.csv -o -name RTAComplete.txt -o -name *unParameters.xml ";
		$inputs_runs_mindepth="1";
		$inputs_runs_maxdepth="1";

		# inputs runs filter list
		$inputs_runs_filter_list=" $folder_inputs/".join(" $folder_inputs/",array_unique($inputs_runs_filter));

		# inputs runs command
		$command="./indexes.sh '.' '$inputs_runs_filter_list' '$inputs_runs_files_patterns' '$inputs_runs_mindepth' '$inputs_runs_maxdepth' $folder_indexes/runs.idx $folder_indexes/runs.idx.tmp";

		# inputs runs exec
		if ($indexes_old) {
			$output = executeAsyncShellCommand($command);
		} else {
			$output = shell_exec($command);
		}


		### Manifests

		# inputs manifests filter
		$inputs_manifests_filter=array_keys($configuration["inputs"]["manifests"]["filter"]);
		if (empty($inputs_manifests_filter)) {
			$inputs_manifests_filter=array("*/manifests/*");
		}

		# inputs manifests parameters
		$inputs_manifests_files_patterns="";
		$inputs_manifests_mindepth="0";
		$inputs_manifests_maxdepth="0";

		# inputs manifests filter list
		$inputs_manifests_filter_list=" $folder_inputs/".join(" $folder_inputs/",array_unique($inputs_manifests_filter));

		# inputs manifests command
		$command="./indexes.sh '.' '$inputs_manifests_filter_list' '$inputs_manifests_files_patterns' '$inputs_manifests_mindepth' '$inputs_manifests_maxdepth' $folder_indexes/manifests.idx $folder_indexes/manifests.idx.tmp";
		
		# inputs manifests exec
		if ($indexes_old) {
			$output = executeAsyncShellCommand($command);
		} else {
			$output = shell_exec($command);
		}


		### Repositories

		# repositories filter
		$repositories_filter=array_keys($configuration["repositories"]["filter"]);
		if (empty($repositories_filter)) {
			$repositories_filter=array("*/*/*/*");
		}

		# repositories parameters
		#$repositories_files_patterns=" -name *STARKCopyComplete.txt -o -name *stark.report.html -o -name *stark.report.pdf -o -name *tag -o -name *vcf.gz -o -name *vcf -o -name *tsv -o -name *cram -o -name *bed ";
		$repositories_files_patterns=" -name *STARKCopyComplete.txt -o -name *report.html -o -name *report.pdf -o -name *tag -o -name *vcf.gz -o -name *vcf -o -name *tsv -o -name *cram -o -name *bed ";
		$repositories_mindepth="2";
		$repositories_maxdepth="2";

		# repositories filter list
		$repositories_filter_list=" $folder_repositories/".join(" $folder_repositories/",array_unique($repositories_filter));

		# repositories command
		$command="./indexes.sh '.' '$repositories_filter_list' '$repositories_files_patterns' '$repositories_mindepth' '$repositories_maxdepth' $folder_indexes/repositories.idx $folder_indexes/repositories.idx.tmp";
		
		# repositories exec
		if ($indexes_old) {
			$output = executeAsyncShellCommand($command);
		} else {
			$output = shell_exec($command);
		}


	}
}


?>