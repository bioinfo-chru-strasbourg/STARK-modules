<?php


	### INCLUDES
	###############

	include "config.php";


	### VARIABLES
	###############

	$DEBUG=1;



	### FUNCTIONS
	###############




	# ENV VARIABLES
	#################




	# DOCKER_STARK_SERVICE_DASHBOARD_INNER_FOLDER
	$DOCKER_STARK_SERVICE_DASHBOARD_INNER_FOLDER=$_SERVER["DOCKER_STARK_SERVICE_DASHBOARD_INNER_FOLDER"];
	if ($DOCKER_STARK_SERVICE_DASHBOARD_INNER_FOLDER=="") { $DOCKER_STARK_SERVICE_DASHBOARD_INNER_FOLDER="/app"; };


	// # DOCKER_STARK_SERVICE_DATA_INNER_FOLDER_DATA
	// $DOCKER_STARK_SERVICE_DATA_INNER_FOLDER_DATA=$_SERVER["DOCKER_STARK_SERVICE_DATA_INNER_FOLDER_DATA"];
	// if ($DOCKER_STARK_SERVICE_DATA_INNER_FOLDER_DATA=="") { $DOCKER_STARK_SERVICE_DATA_INNER_FOLDER_DATA="/app/igvjs/static/data/public"; };


	# DOCKER_STARK_SERVICE_DATA_SUBFOLDER_REPOSITORIES
	$DOCKER_STARK_SERVICE_DATA_SUBFOLDER_REPOSITORIES=$_ENV["DOCKER_STARK_SERVICE_DATA_SUBFOLDER_REPOSITORIES"];
	if ($DOCKER_STARK_SERVICE_DATA_SUBFOLDER_REPOSITORIES=="") { $DOCKER_STARK_SERVICE_DATA_SUBFOLDER_REPOSITORIES="repositories"; };


	# DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG
	$DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG=$_ENV["DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG"];
	#if ($DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG=="") { $DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG="analyses/stark-services/launcher"; };
	if ($DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG=="") { $DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG="services/stark/stark/api"; };


	# DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER
	$DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER=$_ENV["DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER"];
	if ($DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER=="") { $DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER="8000"; };


	# DOCKER_STARK_INNER_FOLDER_ANALYSES
	$DOCKER_STARK_INNER_FOLDER_ANALYSES=$_ENV["DOCKER_STARK_INNER_FOLDER_ANALYSES"];
	#if ($DOCKER_STARK_INNER_FOLDER_ANALYSES=="") { $DOCKER_STARK_INNER_FOLDER_ANALYSES="/STARK/analyses"; };
	if ($DOCKER_STARK_INNER_FOLDER_ANALYSES=="") { $DOCKER_STARK_INNER_FOLDER_ANALYSES="/STARK/services"; };


	# DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_LAUNCHER
	$DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_LAUNCHER=$_ENV["DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_LAUNCHER"];
	#if ($DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_LAUNCHER=="") { $DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_LAUNCHER="stark-services/launcher"; };
	if ($DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_LAUNCHER=="") { $DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_LAUNCHER="stark/stark/api"; };




	### VARIABLES
	###############


	$ANALYSIS_INPUT_FOLDER="$DOCKER_STARK_SERVICE_DASHBOARD_INNER_FOLDER/$DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG";
	$ANALYSIS_INPUT_FOLDER_STARK="$DOCKER_STARK_INNER_FOLDER_ANALYSES/$DOCKER_STARK_SERVICE_DATA_SUBFOLDER_SERVICES_LAUNCHER";
	$ANALYSIS_REPOSITORY_FOLDER="$DOCKER_STARK_SERVICE_DASHBOARD_INNER_FOLDER/$DOCKER_STARK_SERVICE_DATA_SUBFOLDER_REPOSITORIES/Analyses";

	# LAUNCHER
	$LAUNCHER="http://stark-service-launcher:$DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER/analysis";


	#$DATA_URL_EXTERNAL=preg_replace("/\/$/", "", str_replace("//","/","$DOCKER_STARK_SERVER_NAME:$DOCKER_STARK_SERVICE_PORT_PATTERN$DOCKER_STARK_SERVICE_DATA_PORT/$DOCKER_STARK_SERVICE_DATA_PUBLIC_DIR/$DOCKER_STARK_SERVICE_LAUNCHER_FOLDER_LOG"));


	### HEADER
	############

	echo '<!-- Site made with Mobirise Website Builder v4.10.3, https://mobirise.com -->
	<meta charset="UTF-8">
	<meta http-equiv="X-UA-Compatible" content="IE=edge">
	<meta name="generator" content="Mobirise v4.10.3, mobirise.com">
	<meta name="viewport" content="width=device-width, initial-scale=1, minimum-scale=1">
	<link rel="shortcut icon" href="assets/favicon.ico" type="image/x-icon">
	<meta name="description" content="STARK Analysis">
	<title>STARK Analysis</title>
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

	// <style>
	// 	.table-wrapper {
	// 		overflow: auto;
	// 		max-height:800px;
	// 		}
	// </style>


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
				  <span class="navbar-caption-wrap"><a class="navbar-caption text-secondary display-2" href="">Analysis</a></span>
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


	echo '
		<section class="header1 cid-ru7OEConn1" id="header16-1k">
		    <div class="container">
				<h3 class="mbr-section-subtitle mbr-fonts-style align-center mbr-light display-2">
					&nbsp;
			  	</h3>
		    </div>

		</section>

	';



	### Launch Analysis
	####################

	// echo "<pre>";
	// print_r($_POST);
	// print_r($_FILES);
	// echo "</pre>";


	$DATE=date("ymd");


	# FORM array
	$FORM_array=array();

	# JSON array
	$JSON_array=array();


	#echo "ANALYSIS=$ANALYSIS";

	// echo "<br>??? ANALYSIS_INPUT_FOLDER=$ANALYSIS_INPUT_FOLDER";
	// echo "<br>??? ANALYSIS_REPOSITORY_FOLDER=$ANALYSIS_REPOSITORY_FOLDER";

	if ($_POST["analysis"]=="Analysis") {
		$ANALYSIS=1;
		#$ANALYSIS_ID=$DATE."-".rand(0,100000);
		$ANALYSIS_RANDOM=rand(0,100000);
		$ANALYSIS_RANDOM=strtoupper(substr(md5(microtime()),rand(0,26),5));
		$ANALYSIS_NAME=$DATE."_STARK_ANALYSIS-".$ANALYSIS_RANDOM;
		$ANALYSIS_ID="ID-$ANALYSIS_NAME-NAME-$ANALYSIS_NAME";
		# Create ANALYSIS ID folder
		$ANALYSIS_ID_FOLDER="$ANALYSIS_INPUT_FOLDER/$ANALYSIS_ID";
		$ANALYSIS_ID_FOLDER_STARK="$ANALYSIS_INPUT_FOLDER_STARK/$ANALYSIS_ID";
		if (!mkdir($ANALYSIS_ID_FOLDER, 0777, true)) {
			die('Failed to create folders '.$ANALYSIS_ID_FOLDER.'...');
		}
		$JSON_array["analysis_name"]=$ANALYSIS_ID;
		#$JSON_array["repository"]=$ANALYSIS_REPOSITORY_FOLDER;
		$JSON_array["repository"]="$DOCKER_STARK_INNER_FOLDER_ANALYSES/stark";
		#$JSON_array["results"]="$DOCKER_STARK_INNER_FOLDER_ANALYSES/stark";

		# "$ANALYSIS_ID_FOLDER_STARK/$INPUT";
	} else {
		$ANALYSIS=0;
	};



	# APPS
	########


	$APPS=array();
	$APPS_FILE_ARRAY=file("config/apps.conf");

	foreach ($APPS_FILE_ARRAY as $key => $APPS_INFOS) {
		$APPS_INFOS_SPLIT=explode("\t",$APPS_INFOS);
		$APPS[$APPS_INFOS_SPLIT[0]]=trim($APPS_INFOS_SPLIT[1]);
	}

	#print_r($APPS);




	### CONTENT
	#############

	#$URL_LAUNCHER="http://$DOCKER_STARK_SERVICE_LAUNCHER_CONTAINER_NAME:$DOCKER_STARK_SERVICE_LAUNCHER_PORT_INNER/queue";

	#$TS_LIST=file($URL_LAUNCHER);




	# Copy all files into analysis folder
	foreach ($_FILES as $VARIABLE => $VARIABLE_FILE_INFOS) {
		$INPUT=$VARIABLE_FILE_INFOS["name"];
		$INPUT_tmp_name=$VARIABLE_FILE_INFOS["tmp_name"];
		if (is_array($INPUT)) {
			foreach ($INPUT as $INPUT_key => $INPUT_value) {
				$INPUT_value_tmp_name=$VARIABLE_FILE_INFOS["tmp_name"][$INPUT_key];
				if ($INPUT_value != "" && $INPUT_value_tmp_name != "") {
					if (!copy($INPUT_value_tmp_name, "$ANALYSIS_ID_FOLDER/$INPUT_value")) {
					    echo "failed to copy $INPUT_tmp_name...\n";
					}
				}
			}
		} else {
			if ($INPUT != "" && $INPUT_tmp_name != "") {
				if (!copy($INPUT_tmp_name, "$ANALYSIS_ID_FOLDER/$INPUT")) {
				    echo "failed to copy $INPUT_tmp_name...\n";
				}
			}
		}

	}





	// # Multiple files
	// ##################
	//
	// $VARIABLE="input_files[]";
	// $LABEL="Files";
	// $COMMENT="";
	//
	// if ($ANALYSIS) {
	// 	$INPUT=$_FILES[$VARIABLE]["name"];
	// 	#$INPUT_tmp_name=$_FILES[$VARIABLE]["tmp_name"];
	// 	if ($INPUT!="") {
	// 		$JSON_array[$VARIABLE]="$ANALYSIS_ID_FOLDER_STARK/$INPUT";
	// 	}
	// } else {
	// 	$INPUT='<input type="file" name="'.$VARIABLE.'" id="'.$VARIABLE.'" multiple></input>';
	// }
	//
	// $VARIABLE_array["LABEL"]=$LABEL;
	// $VARIABLE_array["INPUT"]=$INPUT;
	// $VARIABLE_array["COMMENT"]=$COMMENT;
	//
	// $FORM_array[$VARIABLE]=$VARIABLE_array;





	# Application
	##########

	$VARIABLE="application";
	$LABEL="Application";
	$COMMENT="STARK Application";

	if ($ANALYSIS) {
		$INPUT=$_POST[$VARIABLE];
		#$INPUT_tmp_name=$_FILES[$VARIABLE]["tmp_name"];
		if ($INPUT!="") {
			$JSON_array[$VARIABLE]="$INPUT";
		}
		$INPUT=$_POST[$VARIABLE]." - ".$APPS[$_POST[$VARIABLE]];
	} else {
		$OPTIONS="";
		foreach ($APPS as $APP => $APP_description) {
			$OPTIONS.="<option value='$APP' title='$APP_description'>$APP</option>";
		}
		$INPUT='<select name="'.$VARIABLE.'" id="'.$VARIABLE.'">'.$OPTIONS.'</select>';
	}

	$VARIABLE_array["LABEL"]=$LABEL;
	$VARIABLE_array["INPUT"]=$INPUT;
	$VARIABLE_array["COMMENT"]=$COMMENT;

	$FORM_array[$VARIABLE]=$VARIABLE_array;



	# READS1
	##########

	$VARIABLE="reads";
	$LABEL="Reads 1";
	$COMMENT="Either FASTQ, BAM, CRAM or RUN tarball";

	if ($ANALYSIS) {
		$INPUT=$_FILES[$VARIABLE]["name"];
		#$INPUT_tmp_name=$_FILES[$VARIABLE]["tmp_name"];
		if ($INPUT!="") {
			$JSON_array[$VARIABLE]="$ANALYSIS_ID_FOLDER_STARK/$INPUT";
		}
	} else {
		$INPUT='<input type="file" name="'.$VARIABLE.'" id="'.$VARIABLE.'"></input>';
	}

	$VARIABLE_array["LABEL"]=$LABEL;
	$VARIABLE_array["INPUT"]=$INPUT;
	$VARIABLE_array["COMMENT"]=$COMMENT;

	$FORM_array[$VARIABLE]=$VARIABLE_array;


	# READS2
	##########

	$VARIABLE="reads2";
	$LABEL="Reads 2";
	$COMMENT="Complementary Reads 2 FASTQ (only if Reads 1 as FASTQ format)";

	if ($ANALYSIS) {
		$INPUT=$_FILES[$VARIABLE]["name"];
		if ($INPUT!="") {
			$JSON_array[$VARIABLE]="$ANALYSIS_ID_FOLDER_STARK/$INPUT";
		}
	} else {
		$INPUT='<input type="file" name="'.$VARIABLE.'" id="'.$VARIABLE.'"></input>';
	}

	$VARIABLE_array["LABEL"]=$LABEL;
	$VARIABLE_array["INPUT"]=$INPUT;
	$VARIABLE_array["COMMENT"]=$COMMENT;

	$FORM_array[$VARIABLE]=$VARIABLE_array;


	# DESIGN
	##########

	$VARIABLE="design";
	$LABEL="Target Design";
	$COMMENT="Design in BED or Illumina Manifest format";

	if ($ANALYSIS) {
		$INPUT=$_FILES[$VARIABLE]["name"];
		if ($INPUT!="") {
			$JSON_array[$VARIABLE]="$ANALYSIS_ID_FOLDER_STARK/$INPUT";
		}
	} else {
		$INPUT='<input type="file" name="'.$VARIABLE.'" id="'.$VARIABLE.'"></input>';
	}

	$VARIABLE_array["LABEL"]=$LABEL;
	$VARIABLE_array["INPUT"]=$INPUT;
	$VARIABLE_array["COMMENT"]=$COMMENT;

	$FORM_array[$VARIABLE]=$VARIABLE_array;


	# GENES
	##########

	$VARIABLE="genes";
	$LABEL="Genes Panels";
	$COMMENT="Panels of genes in BED format";

	if ($ANALYSIS) {
		$INPUT_array=$_FILES[$VARIABLE]["name"];
		if (count($INPUT_array)) {
			$INPUT="$ANALYSIS_ID_FOLDER_STARK/".implode("+$ANALYSIS_ID_FOLDER_STARK/",$_FILES[$VARIABLE]["name"]);
			$JSON_array[$VARIABLE]="$INPUT";
		}
	} else {
		$INPUT='<input type="file" name="'.$VARIABLE.'[]" id="'.$VARIABLE.'[]" multiple></input>';
	}

	$VARIABLE_array["LABEL"]=$LABEL;
	$VARIABLE_array["INPUT"]=$INPUT;
	$VARIABLE_array["COMMENT"]=$COMMENT;

	$FORM_array[$VARIABLE]=$VARIABLE_array;



	# Submit
	##########

	$VARIABLE="analysis";
	$LABEL="";
	$COMMENT="";

	if ($ANALYSIS) {
		#$INPUT=$_POST[$VARIABLE];
		$INPUT="Analysis $ANALYSIS_ID";
	} else {
		$INPUT='<input type="submit" name="'.$VARIABLE.'" value="Analysis" id="'.$VARIABLE.'" class="nav-link"></input>';
	}

	$VARIABLE_array["LABEL"]=$LABEL;
	$VARIABLE_array["INPUT"]=$INPUT;
	$VARIABLE_array["COMMENT"]=$COMMENT;

	$FORM_array[$VARIABLE]=$VARIABLE_array;


	if ($ANALYSIS) {
		
		# JSON
		$JSON=json_encode($JSON_array)."\n";

		$JSON_file="$ANALYSIS_ID_FOLDER/analysis.json";
		file_put_contents($JSON_file, $JSON);

		$JSON_file_STARK="$ANALYSIS_ID_FOLDER_STARK/analysis.json";
		$JSON_analysis_array['analysis']=$JSON_file_STARK;

		$JSON_analysis=json_encode($JSON_analysis_array);
		
		# CMD="curl -s -X POST -H 'Content-Type: application/json' -d '{\"run\":\"$IFA\"}' $LAUNCHER"

		$CURL_OPTIONS=array(
			CURLOPT_POST => 1,
			CURLOPT_HEADER => 'Content-Type: application/json',
			CURLOPT_URL => $LAUNCHER,
			CURLOPT_FRESH_CONNECT => 1,
	        CURLOPT_RETURNTRANSFER => 1,
	        CURLOPT_FORBID_REUSE => 1,
	        CURLOPT_TIMEOUT => 4,
	        #CURLOPT_POSTFIELDS => http_build_query($JSON_analysis)
			#CURLOPT_POSTFIELDS => $JSON_analysis
			CURLOPT_POSTFIELDS => $JSON
		);

		// create a new cURL resource
		if (1) {
			$ch = curl_init();
			curl_setopt_array($ch, ($CURL_OPTIONS));
			if( ! $result = curl_exec($ch)) {
		    	trigger_error(curl_error($ch));
		   	}
		    curl_close($ch);
		    #echo $result;
		}

	}

	# /index.reports.php?PATH=repositories/Analyses/*/*/$ANALYSIS_ID/


	// echo "<pre>";
	// print_r($JSON_array);
	// print_r($CURL_OPTIONS);
	// echo $JSON;
	// echo "</pre>";





	echo '
	<section class="section-table cid-ruiBoanIwc" id="METRICS">

	  <div class="container container-table">

		  <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
			  STARK Analysis
		  </h2>

		  <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
		  	Select input files to launch STARK ANALYSIS<br>
			(Reads 1 mandatory)
		  </p>

		  <div class="div-wrapper" style="max-height:600px;">

			<div class="container ">

				<form method="POST" enctype="multipart/form-data">

				<table class="table" cellspacing="0">

					<tbody>';


	foreach ($FORM_array as $VARIABLE => $VARIABLE_array) {
		$LABEL=$VARIABLE_array["LABEL"];
		$COMMENT=$VARIABLE_array["COMMENT"];
		$INPUT=""; #$VARIABLE_array["INPUT"];

		if ($ANALYSIS) {
			foreach (explode("+",$VARIABLE_array["INPUT"]) as $key => $value) {
				$INPUT.=basename($value)." ";
			}
		} else {
			$INPUT=$VARIABLE_array["INPUT"];
		}


		echo '<tr class="table-heads">
				<td class="head-item mbr-fonts-style display-7" title="'.$COMMENT.'">
					<B>'.$LABEL.'</B>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					'.$INPUT.'
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small><I>'.$COMMENT.'</I></small>
				</td>
			</tr>
		';
	}



	echo '
					</tbody>
				</table>

				</form>

			</div>



		  </div>

		</div>

	</section>

	<br>
   	';




	# mbri-cust-feedback mbri-desktop

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
