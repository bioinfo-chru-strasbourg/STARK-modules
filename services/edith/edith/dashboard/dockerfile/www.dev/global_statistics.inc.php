<?php


### INCLUDE
#############

include "scan.inc.php";


### INPUTS statistics


// #foreach (glob ($folder_inputs."/*/*/*" ) as $key => $input_path) {
// foreach ($runs_inputs as $key => $input_path) {
// 	$run_split=explode ( "/" , $input_path );
// 	$run_path_count=count($run_split);
// 	$input=$run_split[$run_path_count-3];
// 	$input_type=$run_split[$run_path_count-2];
// 	$input_file=$run_split[$run_path_count-1];
// 	$input_ext=pathinfo($input_path,PATHINFO_EXTENSION);
// 	#echo "<br>$input | $input_type | $input_file | $input_ext";
// 	$input_list[$input]++;
// 	if ($input_type=="manifests" && is_file($input_path)) {
// 		if ($input_ext=="genes") {
// 			$genes[$input][$input_file]++;
// 		} elseif ($input_ext=="transcripts") {
// 			$transcripts[$input][$input_file]++;
// 		} elseif ($input_ext=="manifest") {
// 			$designs[$input][$input_file]++;
// 		} elseif ($input_ext=="txt") {
// 			$designs[$input][$input_file]++;
// 		} else {
// 			#$designs[$input][$input_file]++;
// 		};
// 	} elseif ($input_type=="runs" && is_dir($input_path)) {
// 		$runs[$input][$input_file]++;
// 	};
// };


### REPOSITORIES Statistics

#STARKCopyComplete.txt

#foreach (glob ( $folder_repositories."/*/*/*/*/*/STARKCopyComplete.txt", GLOB_BRACE ) as $key => $sample_path) {
// foreach ($runs_repositories as $key => $sample_path) {

// 	print_r($sample_path); echo "<br>";

// 	$sample_split=explode ( "/" , $sample_path );
// 	$sample_path_count=count($sample_split);
// 	#print $sample_path_count;
// 	$repository=$sample_split[$sample_path_count-6];
// 	$group=$sample_split[$sample_path_count-5];
// 	$project=$sample_split[$sample_path_count-4];
// 	$run=$sample_split[$sample_path_count-3];
// 	$sample=$sample_split[$sample_path_count-2];

// 	#echo $sample;
// 	$global[$repository]["groups"][$group]++;
// 	$global[$repository]["projects"][$project]++;
// 	$global[$repository]["runs"][$run]++;
// 	$global[$repository]["samples"][$sample]++;
// 	$global[$repository]["tree"][$group][$project][$run][$sample]++;

// 	$run_by_group[$repository][$group][$run]++;
// 	$sample_by_group[$repository][$group][$sample]++;
// 	$run_by_group_project[$repository][$group][$project][$run]++;
// 	$sample_by_group_project[$repository][$group][$project][$sample]++;

// };


// foreach (glob ( $folder_repositories."/*/*/*/*/*", GLOB_ONLYDIR|GLOB_BRACE ) as $key => $sample_path) {
// 	$sample_split=explode ( "/" , $sample_path );
// 	$sample_path_count=count($sample_split);
// 	#print $sample_path_count;
// 	$repository=$sample_split[$sample_path_count-5];
// 	$group=$sample_split[$sample_path_count-4];
// 	$project=$sample_split[$sample_path_count-3];
// 	$run=$sample_split[$sample_path_count-2];
// 	$sample=$sample_split[$sample_path_count-1];
// 	#echo $sample;
// 	$global[$repository]["groups"][$group]++;
// 	$global[$repository]["projects"][$project]++;
// 	$global[$repository]["runs"][$run]++;
// 	$global[$repository]["samples"][$sample]++;
// 	$global[$repository]["tree"][$group][$project][$run][$sample]++;

// 	$run_by_group[$repository][$group][$run]++;
// 	$sample_by_group[$repository][$group][$sample]++;
// 	$run_by_group_project[$repository][$group][$project][$run]++;
// 	$sample_by_group_project[$repository][$group][$project][$sample]++;

// };

?>