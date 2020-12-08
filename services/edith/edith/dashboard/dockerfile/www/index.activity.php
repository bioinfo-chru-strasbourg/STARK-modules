<?php


#############
### INFOS ###
#############

$APP_SECTION="Activity";



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

#$DEBUG=0;
$DEBUG=$_REQUEST["DEBUG"];
$TS_SHOW=$_REQUEST["TS_SHOW"];


### DATA
##########

$data_read_file=$read_file_content;

#include "global_statistics.inc.php";
include "activity.inc.php";






### TASK SPOOLER
##################


### CONTENT


$thead_task_spooler='
	<tr class="table-heads">
		<th class="head-item mbr-fonts-style display-7">
			ID
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Status
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Time
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Analysis/Run & Command ID
		</th>
	</tr>
	';

$tbody__task_spooler="";

#foreach ($full_task["STARK"] as $ONE_TASK_RUN=>$ONE_TASK_RUN_INFO) {
foreach ($run_task as $ONE_TASK_RUN=>$ONE_TASK_RUN_INFO) {

	$ONE_TASK_ID=$ONE_TASK_RUN_INFO["task"]["id"];
	$ONE_TASK_STATE_COLOR=$ONE_TASK_RUN_INFO["task"]["color"];
	$ONE_TASK_STATE=$ONE_TASK_RUN_INFO["task"]["status"];

	$ONE_TASK_E_LEVEL=$ONE_TASK_RUN_INFO["task"]["error_level"];
	$ONE_TASK_E_LEVEL_HTML=$ONE_TASK_RUN_INFO["task"]["error_level_html"];

	$ONE_TASK_LOG_URL=$ONE_TASK_RUN_INFO["task"]["log_url"];
	$ONE_TASK_LOG_NAME=$ONE_TASK_RUN_INFO["task"]["log_name"];

	$ONE_TASK_TIME_COLOR=$ONE_TASK_RUN_INFO["task"]["time_color"];
	$ONE_TASK_TIME_FORMATED=$ONE_TASK_RUN_INFO["task"]["time_formated"];

	$ONE_TASK_KEY=$ONE_TASK_RUN_INFO["task"]["key"];
	$ONE_TASK_TYPE=$ONE_TASK_RUN_INFO["task"]["type"];

	$ONE_TASK_COMMAND_ID=$ONE_TASK_RUN_INFO["task"]["command_id"];

	$ONE_TASK_INFO_CONTENT_HTML=$ONE_TASK_RUN_INFO["task"]["info_content_html"];


	$ONE_TASK_RUN_HEAD="[$ONE_TASK_TYPE] $ONE_TASK_RUN";

	#if ($ONE_TASK_TYPE=="STARK") {
	if (1) {

		$tbody_task_spooler=$tbody_task_spooler.'
			<tr class="table-heads">
				<td class="head-item mbr-fonts-style display-7">
					<small>'.$ONE_TASK_ID.'</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small style="color:'.$ONE_TASK_STATE_COLOR.'">
						'.$ONE_TASK_STATE.'
						'.$ONE_TASK_E_LEVEL_HTML.'
						<br>
						<a href="'.$ONE_TASK_LOG_URL.'" download>'.$ONE_TASK_LOG_NAME.'</a>
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small style="color:'.$ONE_TASK_TIME_COLOR.'">'.$ONE_TASK_TIME_FORMATED.'</small>
				</td>
				<td class="head-item mbr-fonts-style display-7 ">
					'.$ONE_TASK_RUN_HEAD.'
					<br>
					<a>
					<small onclick=\'if(document.getElementById("info_'.$ONE_TASK_KEY.'").style.display=="none"){document.getElementById("info_'.$ONE_TASK_KEY.'").style.display="block"}else{document.getElementById("info_'.$ONE_TASK_KEY.'").style.display="none"};\'>'.$ONE_TASK_COMMAND_ID.'</small>
					</a>
					<div id="info_'.$ONE_TASK_KEY.'" style="display:none">
						<small><br>'.$ONE_TASK_INFO_CONTENT_HTML.'</small>
					</div>
				</td>

			</tr>

			';
	};
};


### SECTION

if ($TS_SHOW) {

	echo '

	<section class="section-table cid-ruiBoanIwc" id="METRICS">

	  <div class="container container-table">

		  <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
			  Task Spooler
		  </h2>

		  <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
			  Task Spooler queue up STARK analyses
			  <br>
			  <small>[<a href="?TS_SHOW=0">Hide Task Spooler activity</a>]</small>
		  </p>

		  <div class="div-wrapper" style="max-height:600px;">

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
						'.$thead_task_spooler.'
					</thead>
					<tbody>
						'.$tbody_task_spooler.'
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

	<br>

	';

}



### PROGRESS
##############


### CONTENT

$thead_progress='
<tr class="table-heads">
	<th class="head-item mbr-fonts-style display-7">
		Run
	</th>
	<th class="head-item mbr-fonts-style display-7">
		Sequencing
	</th>
	<th class="head-item mbr-fonts-style display-7">
		Analysis
	</th>
	<th class="head-item mbr-fonts-style display-7">
		Repository
	</th>
</tr>
';

$tbody_progress="";

foreach ($run_progress as $run=>$run_infos) {

	$tbody_progress=$tbody_progress.'
		<tr class="table-heads">
			<td class="head-item mbr-fonts-style display-7">
				<small><b>'.$run.'</b></small>
				<br>
				<small style="color:'.$run_color.'"></small>


			</td>
			<td class="head-item mbr-fonts-style display-7">
				<small style="color:'.$run_progress[$run]["sequencing"]["color"].'"><b>'.$run_progress[$run]["sequencing"]["status"].'</b></small>
				<br>
				<a>
				<small onclick=\'if(document.getElementById("info_sequencing_'.$run.'").style.display=="none"){document.getElementById("info_sequencing_'.$run.'").style.display="block"}else{document.getElementById("info_sequencing_'.$run.'").style.display="none"};\'>'.$run_progress[$run]["sequencing"]["message"].'</small>
				</a>
				<div id="info_sequencing_'.$run.'" style="display:none">
					<small>'.$run_progress[$run]["sequencing"]["message_plus"].'</small>
				</div>
			</td>
			<td class="head-item mbr-fonts-style display-7">
				<small style="color:'.$run_progress[$run]["analysis"]["color"].'"><b>'.$run_progress[$run]["analysis"]["status"].'</b></small>
				<br>
				<a>
				<small onclick=\'if(document.getElementById("info_analysis_'.$run.'").style.display=="none"){document.getElementById("info_analysis_'.$run.'").style.display="block"}else{document.getElementById("info_analysis_'.$run.'").style.display="none"};\'>'.$run_progress[$run]["analysis"]["message"].'</small>
				</a>
				<div id="info_analysis_'.$run.'" style="display:none">
					<small>'.$run_progress[$run]["analysis"]["message_plus"].'</small>
				</div>
			</td>
			<td class="head-item mbr-fonts-style display-7">
				<small style="color:'.$run_progress[$run]["repository"]["color"].'"><b>'.$run_progress[$run]["repository"]["status"].'</b></small>
				<br>
				<a>
				<small onclick=\'if(document.getElementById("info_repository_'.$run.'").style.display=="none"){document.getElementById("info_repository_'.$run.'").style.display="block"}else{document.getElementById("info_repository_'.$run.'").style.display="none"};\'>'.$run_progress[$run]["repository"]["message"].'</small>
				</a>
				<div id="info_repository_'.$run.'" style="display:none">
					<small>'.$run_progress[$run]["repository"]["message_plus"].'</small>
				</div>
			</td>

		</tr>

	';

};


### SECTION

echo '

<section class="section-table cid-ruiBoanIwc" id="METRICS">

  <div class="container container-table">

	  <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
		  Analyses/Runs Progress
	  </h2>

	  <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
		  Analyses and Runs process activity (sequencing, analysis, repository)
		  <br>
		  <small>[<a href="?TS_SHOW=1">Show Task Spooler activity</a>]</small>
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
					'.$thead_progress.'
				</thead>
				<tbody>
					'.$tbody_progress.'
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

<br><br><br>

';



// if ($DEBUG && 0) {
// 	echo "<br><br><br>RUNS_INFOS<br>";
// 	echo "<pre>";
// 	print_r($run_progress);
// 	print_r($runs_infos);
// 	echo "</pre>";
// };



##############
### FOOTER ###
##############

include "footer.inc.php";

?>
