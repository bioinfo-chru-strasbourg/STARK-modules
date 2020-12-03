<?php

# DEV
// Start the clock time in seconds 
$start_time = microtime(true); 

# DEV
#$_SESSION=[];


### INPUTS statistics
#####################


$input_filter=array_keys($configuration["inputs"]["filter"]);
#print_r($input_filter);
if (empty($input_filter)) {
	$input_filter=array("*/*/*");
}
$input_filter_brace="{".join(",",array_unique($input_filter))."}";

#print_r($configuration);

if (!isset($_SESSION["inputs"])) {
	
	# scan
	$runs_inputs=glob("$folder_inputs/$input_filter_brace",GLOB_BRACE);
	$_SESSION["inputs"]["runs_inputs"]=$runs_inputs;
	
} else {

	# cache
	$runs_inputs=$_SESSION["inputs"]["runs_inputs"];

};

# process
foreach ($runs_inputs as $key => $input_path) {
	$run_split=explode ( "/" , $input_path );
	$run_path_count=count($run_split);
	$input=$run_split[$run_path_count-3];
	$input_type=$run_split[$run_path_count-2];
	$input_file=$run_split[$run_path_count-1];
	$input_ext=pathinfo($input_path,PATHINFO_EXTENSION);
	#echo "<br>$input | $input_type | $input_file | $input_ext";
	$input_list[$input]++;
	if ($input_type=="manifests" && is_file($input_path)) {
		if ($input_ext=="genes") {
			$genes[$input][$input_file]++;
		} elseif ($input_ext=="transcripts") {
			$transcripts[$input][$input_file]++;
		} elseif ($input_ext=="manifest" || $input_ext=="txt") {
			$designs[$input][$input_file]++;
		} else {
			$designs[$input][$input_file]++;
		};
	} elseif ($input_type=="runs" && is_dir($input_path)) {
		$runs[$input][$input_file]++;
		$total_runs_list[$input_file]++;
	};
};



### REPOSITORIES Statistics
###########################

#echo "<pre>";
#print_r($configuration);

#print_r($repository_filter);
$repository_filter=array_keys($configuration["repositories"]["filter"]);
#print_r($repository_filter);
if (empty($repository_filter)) {
	$repository_filter=array("*/*/*/*");
}
#echo "stop<br><br>";
#print_r($repository_filter);

#echo "</pre>";

#$repository_filter=array("*/*/*");
$repository_filter_brace="{".join(",",array_unique($repository_filter))."}";
#echo "repository_filter_brace=$repository_filter_brace";
#echo "scan";


if (!isset($_SESSION["repositories"])) {

	# scan
	
	# find run (with copy complete)
	$runs_repositories_folders=glob("$folder_repositories/$repository_filter_brace/STARKCopyComplete.txt",GLOB_BRACE|GLOB_NOSORT);
	#print_r(array_multisort(array_map('filemtime', $runs_repositories_folders), SORT_NUMERIC, SORT_ASC, $runs_repositories_folders));

	$run_folders=[];
	$project_folders=[];
	foreach ($runs_repositories_folders as $key=>$run_folder) {
		$run_folders[$key]=dirname($run_folder);
    	$project_folders[$key]=dirname(dirname($run_folder));
	};
	$run_folders_brace="{".join(",",array_unique($run_folders))."}";
	$project_folders_brace="{".join(",",array_unique($project_folders))."}";
	
	
	# fin samples (with copy complete)
	#$runs_repositories=glob("$run_folders_brace/*/STARKCopyComplete.txt",GLOB_BRACE|GLOB_NOSORT);
	$runs_repositories_tmp=glob("$project_folders_brace/*/*/STARKCopyComplete.txt",GLOB_BRACE|GLOB_NOSORT);
	$runs_repositories=[];
	foreach ($run_folders as $run_folder_filter) {
		$runs_repositories_filter=preg_grep('#^'.$run_folder_filter.'#', $runs_repositories_tmp);
		// echo "<br>COUNT="; echo count($runs_repositories_filter);
		// echo "<br>"; print_r($runs_repositories_filter);
		// echo "<br>";
		$runs_repositories=array_merge($runs_repositories,$runs_repositories_filter);
	};


	#$runs_repositories=glob("$run_folders_brace/*/STARKCopyComplete.txt",GLOB_BRACE);
	#arsort($runs_repositories);
	#print_r(array_multisort(array_map('filemtime', $runs_repositories), SORT_NUMERIC, SORT_DESC, $runs_repositories));
	#usort( $runs_repositories, function( $a, $b ) { return filemtime($a) - filemtime($b); } );

	# find samples (without copycomplete for run) - solution deprecated
	#$runs_repositories=glob("$folder_repositories/$repository_filter_brace/*/*/STARKCopyComplete.txt",GLOB_BRACE);

	# cache
	$_SESSION["repositories"]["runs_repositories"]=$runs_repositories;

} else {
	
	# cache
	$runs_repositories=$_SESSION["repositories"]["runs_repositories"];

};




# process
foreach ($runs_repositories as $key => $sample_path) {

	#print_r($sample_path); echo "<br>";

	$sample_split=explode ( "/" , $sample_path );
	$sample_path_count=count($sample_split);
	#print $sample_path_count;
	$repository=$sample_split[$sample_path_count-6];
	$group=$sample_split[$sample_path_count-5];
	$project=$sample_split[$sample_path_count-4];
	$run=$sample_split[$sample_path_count-3];
	$sample=$sample_split[$sample_path_count-2];

	#echo $sample;
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



# Default Repository 
if (isset($global["Repository"])) {
	$repository_default=$folder_repositories."/Repository";
} else {
	$repository_default=$folder_repositories."/$repository";
}


# DEV
// End the clock time in seconds 
$end_time = microtime(true); 

// Calculate the script execution time 
$execution_time = ($end_time - $start_time); 
  
if ($DEBUG) {
	echo "??? It takes ".$execution_time." seconds to execute the script";
};

?>