<?php

# DEV
// Start the clock time in seconds 
$start_time = microtime(true); 

# DEV
#$_SESSION=[];



### INDEXES
###########

# Create indexes
include "indexes.php";




### INPUTS statistics
#####################


# scan (using index)
if (!isset($_SESSION["inputs"])) {
	$_SESSION["inputs"]["runs_inputs"]=preg_grep("#.*#", array_unique(file("$folder_indexes/inputs.idx")));
};


# cache
$runs_inputs=$_SESSION["inputs"]["runs_inputs"];


# process
foreach ($runs_inputs as $key => $input_path) {

	#echo "<br>$key => $input_path<br>";
	$run_split=explode ( "/" , trim($input_path) );
	$run_path_count=count($run_split);
	$input=$run_split[$run_path_count-3];
	$input_type=$run_split[$run_path_count-2];
	$input_file=$run_split[$run_path_count-1];
	$input_ext=pathinfo($input_path,PATHINFO_EXTENSION);
	#echo "<br>$input | $input_type | $input_path | $input_file | $input_ext";
	$input_list[$input]++;

	if ($input_type=="manifests") {
			#echo "<br>MANIFESTS";
		if ($input_ext=="genes") {
			$genes[$input][$input_file]++;
		} elseif ($input_ext=="transcripts") {
			$transcripts[$input][$input_file]++;
		} elseif ($input_ext=="manifest" || $input_ext=="txt") {
			$designs[$input][$input_file]++;
		} else {
			$designs[$input][$input_file]++;
		};
		
	} elseif ($input_type=="runs") {
		$runs[$input][$input_file]++;
		$total_runs_list[$input_file]++;
	};

};



### REPOSITORIES Statistics
###########################


# scan (using index)
if (!isset($_SESSION["repositories"])) {
	$_SESSION["repositories"]["runs_repositories"]=preg_grep("#^[^\/]*/[^\/]*/[^\/]*/[^\/]*/[^\/]*/[^\/]*/STARKCopyComplete.txt$#", array_unique(file("$folder_indexes/repositories.idx")));
}

# cache
$runs_repositories=$_SESSION["repositories"]["runs_repositories"];


# process
foreach ($runs_repositories as $key => $sample_path) {

	$sample_split=explode ( "/" , trim($sample_path) );
	$sample_path_count=count($sample_split);
	#print $sample_path_count;
	$repository=$sample_split[$sample_path_count-6];
	$group=$sample_split[$sample_path_count-5];
	$project=$sample_split[$sample_path_count-4];
	$run=$sample_split[$sample_path_count-3];
	$sample=$sample_split[$sample_path_count-2];

	$global[$repository]["groups"][$group]++;
	$global[$repository]["projects"][$project]++;
	$global[$repository]["runs"][$run]++;
	$global[$repository]["samples"][$sample]++;
	$global[$repository]["tree"][$group][$project][$run][$sample]++;

	$run_by_group[$repository][$group][$run]++;
	$sample_by_group[$repository][$group][$sample]++;
	$run_by_group_project[$repository][$group][$project][$run]++;
	$sample_by_group_project[$repository][$group][$project][$sample]++;

	$total_runs_list[$run]++;

};


### DEFAULT REPOSITORY
######################

# Default Repository 
if (isset($global["Repository"])) {
	$repository_default=$folder_repositories."/Repository";
} else {
	$repository_default=$folder_repositories."/$repository";
}


### DEV
// End the clock time in seconds 
$end_time = microtime(true); 

// Calculate the script execution time 
$execution_time = ($end_time - $start_time); 
  
if ($DEBUG) {
	echo "??? It takes ".$execution_time." seconds to execute the script";
};

?>