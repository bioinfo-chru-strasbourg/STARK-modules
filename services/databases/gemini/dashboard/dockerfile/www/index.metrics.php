<?php



### INCLUDES
###############

include "config.php";




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




### MODULES/CONTAINERS
#############

$MODULES=array();

# SERVER JARVIS CONTAINER NAME
$DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME=$_ENV["DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME"];
if ($DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME!="") {
	$MODULES["JARVIS"]=$DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME;
};

# SERVER IGV CONTAINER NAME
$DOCKER_STARK_SERVICE_IGV_CONTAINER_NAME=$_ENV["DOCKER_STARK_SERVICE_IGV_CONTAINER_NAME"];
if ($DOCKER_STARK_SERVICE_IGV_CONTAINER_NAME!="") {
	$MODULES["IGV"]=$DOCKER_STARK_SERVICE_IGV_CONTAINER_NAME;
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


### HEADER
############

echo '<!-- Site made with Mobirise Website Builder v4.10.3, https://mobirise.com -->
<meta charset="UTF-8">
<meta http-equiv="X-UA-Compatible" content="IE=edge">
<meta name="generator" content="Mobirise v4.10.3, mobirise.com">
<meta name="viewport" content="width=device-width, initial-scale=1, minimum-scale=1">
<link rel="shortcut icon" href="assets/favicon.ico" type="image/x-icon">
<meta name="description" content="STARK Metrics">
<title>STARK Metrics</title>
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


echo '
<section class="menu cid-qTkzRZLJNu" once="menu" id="menu1-3">
  <nav class="navbar navbar-expand beta-menu navbar-dropdown align-items-center navbar-fixed-top collapsed bg-color transparent">
	  <button class="navbar-toggler navbar-toggler-right" type="button" data-toggle="collapse" data-target="#navbarSupportedContent" aria-controls="navbarSupportedContent" aria-expanded="false" aria-label="Toggle navigation">
		  <div class="hamburger">
			  <span></span>
			  <span></span>
			  <span></span>
			  <span></span>
		  </div>
	  </button>
	  <div class="menu-logo">
		  <div class="navbar-brand">
			  <span class="navbar-logo">
				  <a href="/">
					   <img src="assets/logo.png" alt="Mobirise" title="" style="height: 6rem;">
				  </a>
			  </span>
			  <span class="navbar-caption-wrap"><a class="navbar-caption text-secondary display-2" href="">Metrics</a></span>
		 </div>
	  </div>
	  <div class="collapse navbar-collapse align-center" id="navbarSupportedContent">
	  <!--
		  <p class="mbr-text pb-3 mbr-fonts-style display-5">
			  <img src="assets/logo.png" alt="STARK" title="" style="height: 24rem;">
		  </p>
		  <h1 class="mbr-section-title mbr-bold pb-3 mbr-fonts-style display-2">
			  $ENV_NAME
		  </h1>
		  <p class="mbr-text pb-3 mbr-fonts-style display-5">
			  $ENV_RELEASE - $ENV_DATE
			  <BR>
			  $ENV_DESCRIPTION
			  <BR>
			  © Copyright $ENV_COPYRIGHT - All Rights Reserved
			  <BR>
			  $ENV_LICENCE Licence
			  <BR>
			  <BR>
			  <BR>
		  </p>
		-->
	  </div>
  </nav>
</section>
';


# ENV VARIABLES
#################

# SERVER HOSTNAME/IP
$DOCKER_STARK_SERVER_NAME=$_SERVER["SERVER_NAME"];
if ($DOCKER_STARK_SERVER_NAME=="") { $DOCKER_STARK_SERVER_NAME="localhost"; };

# DOCKER_STARK_SERVICE_PORT_PATTERN
$DOCKER_STARK_SERVICE_PORT_PATTERN=$_ENV["DOCKER_STARK_SERVICE_PORT_PATTERN"];
if ($DOCKER_STARK_SERVICE_PORT_PATTERN=="") { $DOCKER_STARK_SERVICE_PORT_PATTERN="42"; };

# SERVER JARVIS PORT
$DOCKER_STARK_SERVICE_JARVIS_PORT=$_ENV["DOCKER_STARK_SERVICE_JARVIS_PORT"];
if ($DOCKER_STARK_SERVICE_JARVIS_PORT=="") { $DOCKER_STARK_SERVICE_JARVIS_PORT="92"; };

# SERVER DATA CONTAINER NAME
$DOCKER_STARK_SERVICE_DATA_CONTAINER_NAME=$_ENV["DOCKER_STARK_SERVICE_DATA_CONTAINER_NAME"];
if ($DOCKER_STARK_SERVICE_DATA_CONTAINER_NAME=="") { $DOCKER_STARK_SERVICE_DATA_CONTAINER_NAME="stark-service-data"; };

# DOCKER STARK SERVICE DATA PORT INNER
$DOCKER_STARK_SERVICE_DATA_PORT=$_ENV["DOCKER_STARK_SERVICE_DATA_PORT"];
if ($DOCKER_STARK_SERVICE_DATA_PORT=="") { $DOCKER_STARK_SERVICE_DATA_PORT="02"; };

# DOCKER STARK SERVICE DATA PORT INNER
$DOCKER_STARK_SERVICE_DATA_PORT_INNER=$_ENV["DOCKER_STARK_SERVICE_DATA_PORT_INNER"];
if ($DOCKER_STARK_SERVICE_DATA_PORT_INNER=="") { $DOCKER_STARK_SERVICE_DATA_PORT_INNER="5000"; };

# SERVER JARVIS CONTAINER NAME
$DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME=$_ENV["DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME"];
if ($DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME=="") { $DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME="stark-service-jarvis"; };

# SERVER DATA PUBLIC DIR
$DOCKER_STARK_SERVICE_DATA_PUBLIC_DIR=$_ENV["DOCKER_STARK_SERVICE_DATA_PUBLIC_DIR"];
if ($DOCKER_STARK_SERVICE_DATA_PUBLIC_DIR=="") { $DOCKER_STARK_SERVICE_DATA_PUBLIC_DIR="/static/data/public"; };

# DOCKER STARK INNER FOLDER OUTPUT LOG
$DOCKER_STARK_SERVICE_DATA_INNER_FOLDER_DATA=$_ENV["DOCKER_STARK_SERVICE_DATA_INNER_FOLDER_DATA"];
if ($DOCKER_STARK_SERVICE_DATA_INNER_FOLDER_DATA=="") { $DOCKER_STARK_SERVICE_DATA_INNER_FOLDER_DATA="/STARK/output/log"; };

# DOCKER STARK INNER FOLDER OUTPUT LOG
$DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_JARVIS=$_ENV["DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_JARVIS"];
if ($DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_JARVIS=="") { $DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_JARVIS="/STARK/output/log"; };



# MAIN URL
############

# JARVIS URL
#$JARVIS_URL=str_replace("//","/","$DOCKER_STARK_SERVER_NAME:$DOCKER_STARK_SERVICE_PORT_PATTERN$DOCKER_STARK_SERVICE_JARVIS_PORT")."/vision.php";
#$JARVIS_URL=str_replace("//","/","$DOCKER_STARK_SERVICE_JARVIS_CONTAINER_NAME:$DOCKER_STARK_SERVICE_JARVIS_PORT")."/vision.php";


# DATA URL
$DATA_URL_EXTERNAL=preg_replace("/\/$/", "", str_replace("//","/","$DOCKER_STARK_SERVER_NAME:$DOCKER_STARK_SERVICE_PORT_PATTERN$DOCKER_STARK_SERVICE_DATA_PORT/$DOCKER_STARK_SERVICE_DATA_PUBLIC_DIR"));
$DATA_URL=preg_replace("/\/$/", "", str_replace("//","/","$DOCKER_STARK_SERVICE_DATA_CONTAINER_NAME:$DOCKER_STARK_SERVICE_DATA_PORT_INNER/$DOCKER_STARK_SERVICE_DATA_PUBLIC_DIR"));

#DEV
if ($DEBUG) {
	echo "<pre>$JARVIS_URL</pre>";
	echo "<pre>$DATA_URL</pre>";
};


# VARIABLES
#############

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



echo '
	<section class="header1 cid-ru7OEConn1" id="header16-1k">
	    <div class="container">
			<h3 class="mbr-section-subtitle mbr-fonts-style align-center mbr-light display-2">
				&nbsp;
		  	</h3>
	    </div>
		<div class="container div-wrapper">

			<h3 class="mbr-section-subtitle mbr-fonts-style align-left mbr-light display-5">
				<small><small>
				<div class="media">

				</div>
				<br>

				<div class="div-wrapper">

				</div>

				</small></small>
		  	</h3>

	    </div>
	</section>

';




### TABLE
###########

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
		$METRICS_URL="http://$DATA_URL_EXTERNAL/".$METRICS;
		#$METRICS_URL="http://$DATA_URL_EXTERNAL/".$METRICS;
		$METRICS_FILES_ARRAY[$METRICS_FILES_KEY]="metrics_files[]=$METRICS_URL";
		$METRICS_FILENAMES_ARRAY[$METRICS_FILES_KEY]="metrics_filenames[]=$METRICS_NAME";

		# Files
		$run_path=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run;
		$project_path=$root.'/'.$repository.'/'.$group.'/'.$project;
		$group_path=$root.'/'.$repository.'/'.$group;

		# IGV LINKS
		if (isset($MODULES["IGV"])) {
			$IGV_RUN_LINK="<a target='IGV' href='index.igv.php?SEARCH_LEVEL=1&PROCESS=0&SAMPLE_PATH[]=$run_path'>IGV</a>";
		}

		# JARVIS LINKS
		if (isset($MODULES["JARVIS"])) {
			$JARVIS_RUN_LINK="<a target='JARVIS' href='index.jarvis.php?SEARCH_LEVEL=1&SEARCH_EXT=final.vcf.gz,final.bcf&PROCESS=0&SAMPLE_PATH[]=$run_path'>JARVIS</a>";
			$JARVIS_PROJECT_LINK="<a target='JARVIS' href='index.jarvis.php?SEARCH_LEVEL=2&SEARCH_EXT=final.vcf.gz,final.bcf&PROCESS=0&SAMPLE_PATH[]=$project_path'>JARVIS</a>";
			$JARVIS_GROUP_LINK="<a target='JARVIS' href='index.jarvis.php?SEARCH_LEVEL=3&SEARCH_EXT=final.vcf.gz,final.bcf&PROCESS=0&SAMPLE_PATH[]=$group_path'>JARVIS</a>";
		}

		# Samples Reports LINKS
		$REPORTS_RUN_LINK="<a target='REPORT' href='index.reports.php?PATH=$run_path'>REPORT</a>";
		$REPORTS_PROJECT_LINK="<a target='REPORT' href='index.reports.php?PATH=$project_path'>REPORT</a>";
		$REPORTS_GROUP_LINK="<a target='REPORT' href='index.reports.php?PATH=$group_path'>REPORT</a>";

		# TBODY
		$tbody=$tbody.'
			<tr class="table-heads">
				<td class="head-item mbr-fonts-style display-7">
					<small>
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

echo '

<section class="section-table cid-ruiBoanIwc" id="METRICS">

  <div class="container container-table">

	  <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
		  Run Metrics
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


### SCRIPTS
#############

echo '
	<script src="assets/web/assets/jquery/jquery.min.js"></script>
	<script src="assets/popper/popper.min.js"></script>
	<script src="assets/tether/tether.min.js"></script>
	<script src="assets/bootstrap/js/bootstrap.min.js"></script>
	<script src="assets/smoothscroll/smooth-scroll.js"></script>
	<script src="assets/dropdown/js/script.min.js"></script>
	<script src="assets/touchswipe/jquery.touch-swipe.min.js"></script>
	<script src="assets/viewportchecker/jquery.viewportchecker.js"></script>
	<script src="assets/as-pie-progress/jquery-as-pie-progress.min.js"></script>
	<script src="assets/vimeoplayer/jquery.mb.vimeo_player.js"></script>
	<script src="assets/datatables/jquery.data-tables.min.js"></script>
	<script src="assets/datatables/data-tables.bootstrap4.min.js"></script>
	<script src="assets/masonry/masonry.pkgd.min.js"></script>
	<script src="assets/imagesloaded/imagesloaded.pkgd.min.js"></script>
	<script src="assets/bootstrapcarouselswipe/bootstrap-carousel-swipe.js"></script>
	<script src="assets/mbr-switch-arrow/mbr-switch-arrow.js"></script>
	<script src="assets/theme/js/script.js"></script>
	<script src="assets/gallery/player.min.js"></script>
	<script src="assets/gallery/script.js"></script>
	<script src="assets/slidervideo/script.js"></script>
	<script src="assets/mbr-tabs/mbr-tabs.js"></script>
';

?>
