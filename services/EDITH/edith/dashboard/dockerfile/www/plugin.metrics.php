<?php



#############
### INFOS ###
#############

$APP_SECTION="Plugin Metrics";



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

include "header.inc.php";


### VARIABLES
###############



### DATA
##########

#include "global_statistics.inc.php";
#include "activity.inc.php";



# PARAMETERS
##############

#$PATH=$_GET["PATH"];
$PATH=$_GET["PATH"];
$SEARCH_LEVEL=$_GET["SEARCH_LEVEL"];
$SEARCH_EXT=$_GET["SEARCH_EXT"];
$PROCESS=$_GET["PROCESS"];
$DEBUG=0;


# PROCESS
if ($PROCESS=="") {
	$PROCESS=0;
};

# SEARCH_EXT
if ($SEARCH_EXT=="") {
	$SEARCH_EXT="report.*.metrics";
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

# DATA URL

if ($uri_das == "") {
	if ($_ENV["URI_DAS"] != "") {
		$uri_das=$_ENV["URI_DAS"];
	} elseif ($modules_obj_array["STARK"]->{"services"}->{"DAS"}->{"href"}!="") {
		$uri_das=$modules_obj_array["STARK"]->{"services"}->{"DAS"}->{"href"};
	} else {
		$uri_das="";
	};
};



#DEV
if ($DEBUG) {
	echo "<pre>$uri_das</pre>";
};


# VARIABLES
#############

# DATE
$DATE=date('omd-his');

# VCF_FILES
$VCF_FILES = "";




# SEARCH FOR METRICS
#####################

#DEV
if ($DEBUG) {
	print_r($PATH);
};

# TRACKS JSON ARRAY
$METRICS_FILES_ARRAY=array();
$METRICS_FILENAMES_ARRAY=array();

$METRICS_FILES_KEY=0;

#echo $SEARCH_LEVEL;
$PATH_LEVEL="";
for ($i=0; $i<$SEARCH_LEVEL+0; $i++) {
	$PATH_LEVEL.="/*";
}


foreach ($PATH as $key_sample => $ONE_SAMPLE_PATH) {

	# SAMPLE NAME
	$ONE_SAMPLE_NAME=end(explode( "/", $ONE_SAMPLE_PATH ));

	# METRICS
	############

	# Search files
	$metrics=array();
	$metrics=glob("$ONE_SAMPLE_PATH$PATH_LEVEL/*{".$SEARCH_EXT."}",GLOB_BRACE);

	# Create tracks
	foreach ($metrics as $key_metrics => $METRICS) {

		$METRICS_FILES_KEY++;
		$METRICS_split=explode( "/", $METRICS );
		# Find infos
		$root=$METRICS_split[0];
		$repository=$METRICS_split[1];
		$group=$METRICS_split[2];
		$project=$METRICS_split[3];
		$run=$METRICS_split[4];

		$METRICS_NAME=$METRICS_split[5];

		$METRICS_ID=explode( ".", $METRICS_NAME )[1];
		$METRICS_TYPE=explode( ".", $METRICS_NAME )[3];

		$METRICS_TYPE=pathinfo(preg_replace('/.metrics$/i','',$METRICS_NAME))['extension'];
		#$METRICS_BASENAME=pathinfo($METRICS_NAME)['basename'];

		#$METRICS_URL="http://$DATA_URL_EXTERNAL/".$METRICS;
		$METRICS_URL=$uri_das."/".$METRICS;
		#$METRICS_URL="http://$DATA_URL_EXTERNAL/".$METRICS;
		$METRICS_FILES_ARRAY[$METRICS_FILES_KEY]="metrics_files[]=$METRICS_URL";
		$METRICS_FILENAMES_ARRAY[$METRICS_FILES_KEY]="metrics_filenames[]=$METRICS_NAME";

		# Files
		$run_path=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run;
		$project_path=$root.'/'.$repository.'/'.$group.'/'.$project;
		$group_path=$root.'/'.$repository.'/'.$group;

		# Samples Reports LINKS
		$REPORTS_RUN_LINK="<a target='REPORT' href='index.reports.php?PATH=$run_path'>REPORT</a>";
		$REPORTS_PROJECT_LINK="<a target='REPORT' href='index.reports.php?PATH=$project_path'>REPORT</a>";
		$REPORTS_GROUP_LINK="<a target='REPORT' href='index.reports.php?PATH=$group_path'>REPORT</a>";

		# TBODY
		$tbody=$tbody.'
			<tr class="table-heads">
				<td class="head-item mbr-fonts-style display-7">
					<small title="'.$METRICS_NAME.'">
						<a href="'.$METRICS_URL.'" download><b>'.$METRICS_ID.'</b></a>
						<br>
						<b>'.$METRICS_TYPE.'</b>
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						<b>'.$run.'</b>
						<br>
						'.$IGV_RUN_LINK.'
						'.$JARVIS_RUN_LINK.'
						'.$REPORTS_RUN_LINK.'
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$project.'
						<br>
						'.$JARVIS_PROJECT_LINK.'
						'.$REPORTS_PROJECT_LINK.'
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$group.'
						<br>
						'.$JARVIS_GROUP_LINK.'
						'.$REPORTS_GROUP_LINK.'
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$repository.'
					</small>
				</td>
			</tr>
		';



	}

};


### SECTION
#############


### THEAD

$thead='
	<tr class="table-heads">
		<th class="head-item mbr-fonts-style display-7">
			Metrics
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Analysis/Run
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Project
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Group
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Repository
		</th>
	</tr>
';

echo '

<section class="section-table cid-ruiBoanIwc" id="METRICS">

  <div class="container container-table">

	  <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
		  Analysis/Run Metrics
	  </h2>

	  <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
		  Metrics files for runs
	  </p>

	  <div class="">

		<div class="container">
		  <div class="row search">
			<div class="col-md-6"></div>
			<div class="col-md-6">
				<div class="dataTables_filter">
				  <label class="searchInfo mbr-fonts-style display-7">Search:</label>
				  <input class="form-control input-sm" disabled="">
				</div>
			</div>
		  </div>
		</div>

		<div class="container scroll">
			<table class="table isSearch" cellspacing="0">
				<thead>
					'.$thead.'
				</thead>
				<tbody>
					'.$tbody.'
				</tbody>
			</table>
		</div>

		<div class="container table-info-container">
		  <div class="row info">
			<div class="col-md-6">
			  <div class="dataTables_info mbr-fonts-style display-7">
				<span class="infoBefore">Showing</span>
				<span class="inactive infoRows"></span>
				<span class="infoAfter">entries</span>
				<span class="infoFilteredBefore">(filtered from</span>
				<span class="inactive infoRows"></span>
				<span class="infoFilteredAfter"> total entries)</span>
			  </div>
			</div>
			<div class="col-md-6"></div>
		  </div>
		</div>

	  </div>

	</div>

</section>
';


##############
### FOOTER ###
##############

include "footer.inc.php";



?>
