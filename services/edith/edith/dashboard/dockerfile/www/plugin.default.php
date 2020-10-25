<?php


#############
### INFOS ###
#############

$APP_SECTION=$_REQUEST["APP_SECTION"];
if ($APP_SECTION=="") {
	$APP_SECTION="PLUGIN";	
};



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



$DEBUG=$_REQUEST["DEBUG"];


# PARAMETERS
##############

$PATH=$_GET["PATH"];
$SEARCH_LEVEL=$_GET["SEARCH_LEVEL"];



# SEARCH_EXT
$SEARCH_EXT=$_GET["SEARCH_EXT"];
if ($SEARCH_EXT=="") {
	$SEARCH_EXT="annotsv.tsv";
};


$RESULTS_SUBFOLDER_DATA=$_REQUEST["RESULTS_SUBFOLDER_DATA"];
if ($RESULTS_SUBFOLDER_DATA=="") {
	$RESULTS_SUBFOLDER_DATA="CANOE";	
};

$CHECK_SUBFOLDER_DATA=$_REQUEST["CHECK_SUBFOLDER_DATA"];
if ($CHECK_SUBFOLDER_DATA=="") {
	$CHECK_SUBFOLDER_DATA=0;	
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
#$JARVIS_URL=$modules_obj_array["VARIANTBROWSER"]->{"services"}->{"VISION"}->{"href"};

# DATA URL
$DATA_URL=$modules_obj_array["STARK"]->{"services"}->{"DAS"}->{"href_inner"};


# VARIABLES
#############

# VCF_FILES
$ANNOTSV_FILES = "";


$PATH_LEVEL="";
for ($i=0; $i<$SEARCH_LEVEL+0; $i++) {
	$PATH_LEVEL.="/*";
}

foreach ($PATH as $key_sample => $ONE_SAMPLE_PATH) {

	#echo "$key_sample => $ONE_SAMPLE_PATH";

	# SAMPLE NAME
	$ONE_SAMPLE_NAME=end(explode( "/", $ONE_SAMPLE_PATH ));

	# ANNOTSV
	############

	# Search files
	$result=array();
	#echo "$ONE_SAMPLE_PATH$PATH_LEVEL/*{".$SEARCH_EXT."}";
	$result_root=glob("$ONE_SAMPLE_PATH$PATH_LEVEL/*{".$SEARCH_EXT."}",GLOB_BRACE);
	if ($CHECK_SUBFOLDER_DATA) {
		#echo "$ONE_SAMPLE_PATH$PATH_LEVEL/$RESULTS_SUBFOLDER_DATA/*{".$SEARCH_EXT."}";
		$result_data=glob("$ONE_SAMPLE_PATH$PATH_LEVEL/$RESULTS_SUBFOLDER_DATA/*{".$SEARCH_EXT."}",GLOB_BRACE);
	} else {
		$result_data=array();
	};
	$result=array_merge ($result_root,$result_data);

}

// echo "<pre>";
// print_r($result);
// echo "</pre>";




### TABLE
###########


### THEAD

$thead='
	<tr class="table-heads">
		<th class="head-item mbr-fonts-style display-7">
			 File
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Sample
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Analysis / Run
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


### TBODY

$tbody="";



foreach ($result as $result_file) {

	#echo "<br><a href='$result_file' download>$result_file</a>";

	# Find infos
	$result_file_split=explode ( "/" , $result_file );
	$root=$result_file_split[0];
	$repository=$result_file_split[1];
	$group=$result_file_split[2];
	$project=$result_file_split[3];
	$run=$result_file_split[4];
	$sample=$result_file_split[5];
	$result_file_file=$result_file_split[6];
	$result_file_file_end=end(explode ( "/" , $result_file ));
	$result_file_file_ext_end=end(explode ( "." , $result_file_file_end ));

	#echo '<br>'.$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample.'/'.$sample;


	# TBODY
		$tbody=$tbody.'
			<tr class="table-heads">
				<td class="head-item mbr-fonts-style display-7" >
					<a href="'.$result_file.'" title="'.$result_file_file_end.'">
						<div class="media mb-2">
							<div class="card-img align-self-center">
								<span class="mbr-iconfont mbri-file" style="color: rgb(20, 157, 204); fill: rgb(20, 157, 204); "></span> 
								'.$result_file_file_ext_end.'
							</div>
						</div>
						
					</a>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						<b>'.$sample.'</b>
						<br>
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$run.'
						<br>
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$project.'
						<br>
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$group.'
						<br>
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$repository.'
						<br>
						
					</small>
				</td>
			</tr>
		';

}




### SECTION
#############

$CONTENT_SECTION_RESULT_FILES= '

<section class="section-table cid-ruiBoanIwc" id="REPORTS">

  <div class="container container-table">

	  <div class="table-wrapper">

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




###############
### CONTENT ###
###############

$CONTENT=$CONTENT_SECTION_RESULT_FILES;

echo $CONTENT;


##############
### FOOTER ###
##############

include "footer.inc.php";



?>
