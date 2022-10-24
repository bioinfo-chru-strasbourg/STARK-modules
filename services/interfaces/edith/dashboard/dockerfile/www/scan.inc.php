<?php



### INDEXES
###########

# Create indexes
include "indexes.php";


### INPUTS statistics
#####################


# scan (using index)
if (!isset($_SESSION["data"]["index"]["inputs"]["runs"])) {
	$_SESSION["data"]["index"]["inputs"]["runs"]=preg_grep("#.*#", array_unique(file("$folder_indexes/runs.idx")));
};

# cache
$inputs_runs=$_SESSION["data"]["index"]["inputs"]["runs"];

# process
foreach ($inputs_runs as $key => $input_run_file_path) {

	#echo "<br>$key => $input_run_file_path<br>";
	$input_run_file_path=trim($input_run_file_path);

	$input_run_file_name=basename($input_run_file_path);
	$input_run_path=dirname($input_run_file_path);
	$input_run_name=basename($input_run_path);
	$input_name=basename(dirname(dirname($input_run_path)));

	#echo "<br>$input_run_file_name | $input_run_path | $input_run_name | $input_name ";

	$runs[$input_name][$input_run_name][$input_run_file_path]=$input_run_file_name;
	$total_runs_list[$input_run_name][$input_run_file_path]=$input_run_file_name;


};

#print_r($runs);


# scan (using index)
if (!isset($_SESSION["data"]["index"]["inputs"]["manifests"])) {
	$_SESSION["data"]["index"]["inputs"]["manifests"]=preg_grep("#.*#", array_unique(file("$folder_indexes/manifests.idx")));
};

# cache
$inputs_manifests=$_SESSION["data"]["index"]["inputs"]["manifests"];
#print_r($inputs_manifests);


# process
foreach ($inputs_manifests as $key => $input_manifest_file_path) {

	#echo "<br>$key => $input_manifest_file_path<br>";
	$input_manifest_file_path=trim($input_manifest_file_path);

	$input_manifest_file_name=basename($input_manifest_file_path);
	$input_manifest_file_ext=pathinfo($input_manifest_file_name,PATHINFO_EXTENSION);
	$input_manifest_path=dirname($input_manifest_file_path);
	$input_manifest_name=basename($input_manifest_path);
	$input_name=basename(dirname(dirname($input_run_path)));

	#echo "<br>$input_manifest_file_path | $input_manifest_file_name | $input_manifest_file_ext | $input_name <br><br>";

	switch ($input_manifest_file_ext) {
		case "genes":
			$genes[$input_name][$input_manifest_file_path]=$input_manifest_file_name;
			break;
		case "transcripts":
			$transcripts[$input_name][$input_manifest_file_path]=$input_manifest_file_name;
			break;
		default:
			$designs[$input_name][$input_manifest_file_path]=$input_manifest_file_name;
			break;
	}


	// if ($input_manifest_file_ext=="genes") {
	// 	$genes[$input_name][$input_file]++;
	// } elseif ($input_manifest_file_ext=="transcripts") {
	// 	$transcripts[$input_name][$input_file]++;
	// } elseif ($input_manifest_file_ext=="manifest" || $input_ext=="txt") {
	// 	$designs[$input_name][$input_file]++;
	// } else {
	// 	$designs[$input_name][$input_file]++;
	// };


};

// echo "<pre>Manifests and designs";
// print_r($designs);
// print_r($genes);
// print_r($transcripts);
// echo "</pre>";



### REPOSITORIES Statistics
###########################


# scan (using index)
if (!isset($_SESSION["data"]["index"]["repositories"])) {
	#$_SESSION["data"]["index"]["repositories"]["runs_repositories"]=preg_grep("#^[^\/]*/[^\/]*/[^\/]*/[^\/]*/[^\/]*/[^\/]*/STARKCopyComplete.txt$#", array_unique(file("$folder_indexes/repositories.idx")));
	#$_SESSION["data"]["index"]["repositories"]["runs_repositories"]=preg_grep("#^[^\/]*/[^\/]*/[^\/]*/[^\/]*/[^\/]*/[^\/]*/[^\/]*/*stark.report.html$#", array_unique(file("$folder_indexes/repositories.idx")));	
	$_SESSION["data"]["index"]["repositories"]["runs_repositories"]=preg_grep("#^.*stark.report.(html|pdf)$#", array_unique(file("$folder_indexes/repositories.idx")));	
	$_SESSION["data"]["index"]["repositories"]["runs_repositories_samples_files"]=preg_grep("#^[^\/]*/[^\/]*/[^\/]*/[^\/]*/[^\/]*/[^\/]*/.*$#", array_unique(file("$folder_indexes/repositories.idx")));
	#$_SESSION["data"]["index"]["repositories"]["runs_repositories"]=preg_grep("#^[^\/]*/[^\/]*/[^\/]*/[^\/]*/[^\/]*/STARKCopyComplete.txt$#", array_unique(file("$folder_indexes/repositories.idx")));
	#$_SESSION["data"]["index"]["repositories"]["runs_repositories"]=preg_grep("#^[^\/]*/[^\/]*/[^\/]*/[^\/]*/[^\/]*/.*$#", array_unique(file("$folder_indexes/repositories.idx")));
	#$_SESSION["data"]["index"]["repositories"]["runs_repositories"]= array_unique(file("$folder_indexes/repositories.idx"));
	
}

#print_r($_SESSION["data"]["index"]["repositories"]["runs_repositories_samples_files"]);

# cache
$runs_repositories=$_SESSION["data"]["index"]["repositories"]["runs_repositories"];
$runs_repositories_samples_files=$_SESSION["data"]["index"]["repositories"]["runs_repositories_samples_files"];



// echo "<pre>";
// print_r($runs_repositories);
// #print_r(path_to_tree($runs_repositories));
// #print_r(paths_to_array_tree($runs_repositories));
// echo "</pre>";

# process
foreach ($runs_repositories as $key => $file_path) {

	$runs_repositories[$key]=trim($file_path);

	$file_path_split=explode("/",trim($file_path));
	$file_path_split_reverse=array_reverse($file_path_split);

	$file=$file_path_split_reverse[0];
	$sample=$file_path_split_reverse[1];
	$run=$file_path_split_reverse[2];
	$project=$file_path_split_reverse[3];
	$group=$file_path_split_reverse[4];

	$repository=$file_path_split[1];
	$repositories_root=$file_path_split[0];
	
	#echo "$repository | $group | $project | $run | $sample <br><br>";
	
	$global[$repository]["groups"][$group]++;
	$global[$repository]["projects"][$project]++;
	$global[$repository]["runs"][$run]++;
	$global[$repository]["samples"][$sample]++;
	$global[$repository]["tree"][$group][$project][$run][$sample]++;

	$run_by_group[$repository][$group][$run]++;
	$sample_by_group[$repository][$group][$sample]++;
	$run_by_group_project[$repository][$group][$project][$run]++;
	$sample_by_group_project[$repository][$group][$project][$sample]++;

	#$runs_repositories_corrected["$repositories_root/$repository/$group/$project/$run/$sample/$file"]=$file_path;
	$runs_repositories_corrected[$key]=trim("$repositories_root/$repository/$group/$project/$run/$sample/$file");

	$total_runs_list[$run][$runs_repositories_corrected[$key]]=$runs_repositories_corrected[$key];

};

# Clean/trim
foreach ($runs_repositories_samples_files as $key => $file_path) {
	$runs_repositories_samples_files[$key]=trim($file_path);
}

// echo "<pre>";
// print_r($runs_repositories_corrected);
// echo "</pre>";

### DEFAULT REPOSITORY
######################

# Default Repository 
if (isset($global["Repository"])) {
	$repository_default=$folder_repositories."/Repository";
} else {
	$repository_default=$folder_repositories."/$repository";
}


## DEV
// echo "<pre>";
// print_r($inputs_runs);
// print_r($inputs_manifests);
// print_r($runs_repositories);
// print_r($runs_repositories_samples_files);
// echo "</pre>";



?>