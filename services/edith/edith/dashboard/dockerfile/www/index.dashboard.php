<?php


#############
### INFOS ###
#############

$APP_SECTION="Dashboard";


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

$data_read_file=false;

include "global_statistics.inc.php";
include "activity.inc.php";


### CONTENT
#############



#########################
### SECTION DASHBOARD ###
#########################


### DATA

$align="center";


# STATISTICS

$STATISTICS_DIV="";

#foreach ($input_list as $input=>$hash) {
foreach ($runs as $input=>$hash) {
	$STATISTICS_DIV.="
	<div class='p-2 col-4' style='float: left;'>
		<p class='mbr-text mbr-fonts-style display-7' style='text-align: $align;'>				
			<b>$input</b>
			<small>
			<br>
			<b>".count($runs[$input])."</b> runs
			<br>
			<b>".count($designs[$input])."</b> designs
			<br>
			<b>".count($genes[$input])." </b> panels
			<br>
			<b>".count($transcripts[$input])." </b> transcripts
			</small>
		</p>
	</div>
	";

};



foreach ($global as $repository=>$hash) {
	$nb_groups=count($hash["groups"]);
	$nb_projects=count($hash["projects"]);
	$nb_runs=count($hash["runs"]);
	$nb_samples=count($hash["samples"]);
	$STATISTICS_DIV.="
	<div class='p-2 col-4' style='float: left;'>
		<p class='mbr-text mbr-fonts-style display-7' style='text-align: $align;'>
			<b>$repository</b>
			<small>
			<br>
			<b>$nb_groups</b> groups
			<br>
			<b>$nb_projects</b> projects
			<br>
			<b>$nb_runs</b> analyses
			<br>
			<b>$nb_samples</b> samples
			</small>
		</p>
	</div>
	";
};




// echo "<pre>";
// #print_r($run_progress["RUN_TEST_TAG_LISTENER_MIN"]);
// print_r($run_progress);
// echo "</pre>";


$ACTIVITY_STATUS;
foreach ($run_progress as $run=>$run_infos) {

	foreach ($run_infos as $activity=>$activity_infos) {
		$activity_status=$activity_infos["status"];
		$ACTIVITY_STATUS[$activity][$activity_status]++;
		$ACTIVITY_STATUS_COLOR[$activity][$activity_status]=$activity_infos["color"];

	}

};

#print_r($ACTIVITY_STATUS);


$ACTIVITY_DIV="";
#$status_exclude=array("unavailable");
foreach ($ACTIVITY_STATUS as $activity=>$status_list) {
	$ACTIVITY_DIV.="
		<div class='p-2 col-4' style='float: left;'>
			<p class='mbr-text mbr-fonts-style display-7' style='text-align: $align;'>
				<b>".ucfirst($activity)."</b>
				<small>
	";
	foreach ($status_list as $status=>$nb) {
		if (!in_array($status,$status_exclude) ) {
			$ACTIVITY_DIV.="
					<br>
					<b>$nb</b> <span style='color:".$ACTIVITY_STATUS_COLOR[$activity][$status]."'>".$status."</span> 
			";
		};
	};
	$ACTIVITY_DIV.="
			</small>
			</p>
		</div>
	";
};



# List analysis/run
$runs_list="<select name='analysis_from_list' id='analysis_from_list' onchange='this.form.submit()' style='width: 330px;'>";
# onfocus='this.size=10;' onblur='this.size=1;'
#$runs_list="<select name='analysis' id='analysis' form='runs_list' onselect='this.form.submit()'>";
$runs_list.="<option value=''>Select an analysis...</option>";
$runs_list_array=[];
foreach ($global as $repository=>$repository_data) {
	$runs_list_array=array_merge($runs_list_array,$global[$repository]["runs"]);
}
#$runs_list_array=array_merge($runs_list_array,$global[$repository]["runs"]);
krsort($runs_list_array);
foreach ($runs_list_array as $run=>$nb_samples) {
	#echo "$run=>$nb_samples<br>";
	$runs_list.="<option value='$run'>".$run."</option>";
}
$runs_list.="</select>";




$CONTENT_SECTION_DASHBOARD='
<section class="features10 cid-ru7OEDbxhA" id="features10-1l">

	<div class="container ">
		<div class="row justify-content-center">

			<div class="card p-3 col-12 col-md-4 mb-4">
				<a href="index.statistics.php" class="navbar-caption text-secondary ">
					<div class="media mb-2">
						<div class="card-img align-self-center">
							<span class="mbr-iconfont mbri-growing-chart" style="color: rgb(20, 157, 204); fill: rgb(20, 157, 204);"></span>
						</div>
						<h4 class="card-title media-body py-3 mbr-fonts-style display-5">Statistics</h4>
					</div>
				</a>
				<div class="card-box">
					<p class="mbr-text mbr-fonts-style display-7">
						Statistics on runs/analyses, groups, projects and samples,
						stored in repository, archives and sequencing raw sequencing data.
						<br>
						'.$STATISTICS_DIV.'
					</p>
				</div>
			</div>

			<div class="card p-3 col-12 col-md-4 mb-4">
				<a href="index.activity.php" class="navbar-caption text-secondary ">
					<div class="media mb-2">
						<div class="card-img align-self-center">
							<span class="mbr-iconfont mbri-clock" style="color: rgb(20, 157, 204); fill: rgb(20, 157, 204);"></span>
						</div>
						<h4 class="card-title media-body py-3 mbr-fonts-style display-5">Activity</h4>
					</div>
				</a>
				<div class="card-box">
					<p class="mbr-text mbr-fonts-style display-7">
						Activity for runs/analyses sequenced, analysed, queued, running, archived and availabled in repository.
						<br>
						'.$ACTIVITY_DIV.'
					</p>
				</div>
				
			</div>

			<div class="card p-3 col-12 col-md-4 mb-4">
				<a href="index.reports.php?PATH='.$repository_default.'" class="navbar-caption text-secondary ">
					<div class="media mb-2">
						<div class="card-img align-self-center">
							<span class="mbr-iconfont mbri-file" style="color: rgb(20, 157, 204); fill: rgb(20, 157, 204);"></span>
						</div>
						<h4 class="card-title media-body py-3 mbr-fonts-style display-5">Analyses Reports</h4>
					</div>
				</a>
				<div class="card-box">
					<p class="mbr-text mbr-fonts-style display-7">
						Reports Browser provides a direct access to STARK reports and results availabled in reporitories.
						<br>
						<div class="p-2 col-12" style="float: left;">
						<p class="mbr-text mbr-fonts-style display-7" style="text-align: $align;">
							<b title="Search analyses or samples in repositories, using wildcard \'*\' like for analysis search by date/start (\'191009*\'), random id (\'*CYCYG\'), serial number (\'*_0335_*\'), or for sample search by name (\'sampleA\', \'*control*\') ">Search a analysis/sample in repositories</b>
							<form action="index.reports.php" method="POST" >
								<span class="head-item mbr-fonts-style display-7">
								<input name="analysis" value="" style="width:150px" class=""></input>
								<input name="sample" value="" style="width:100px" class=""></input>
								<input type="Submit" value="Search" class=""></input>
								<br>
								'.$runs_list.'
								</span>
							</form>
						

						</p>
						</div>
					</p>

				</div>
				
			</div>

		</div>
	</div>

</section>
';



#######################
### SECTION MODULES ###
#######################



$CONTENT_SECTION_MODULES_HEADER='
<section class="counters1 counters cid-ru7OEH6r7i" id="SECTION_ID" >

	<div class="container">

	   <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
		   STARK Modules & Services IHM
	   </h2>
	   <!--
	   <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
		   Statistics for inputs and repositories folders (duplication not considered)
	   </p>
	   -->
	</div>

</section>
';

		
$CONTENT_SECTION_MODULE_CONTENT='';

$MODULE_NUMBER=0;


### Check modules


foreach ($modules_obj_array as $module_name=>$module_obj) {
	
	
	# Variables
	$module_info_code=$module_obj->{'code'};
	$module_info_name=$module_obj->{'name'};
	$module_info_fullname=$module_obj->{'fullname'}!=""?$module_obj->{'fullname'}:$module_obj->{'name'};
	$module_info_release=$module_obj->{'release'}!=""?$module_obj->{'release'}:"unknown";
	$module_info_description=$module_obj->{'description'};
	$module_info_available=$module_obj->{'available'};

	# If module available
	if ($module_info_available) {

		$CONTENT_SECTION_MODULES_SERVICES_CONTENT="";
		$CONTENT_SECTION_MODULES_SERVICES_CONTENT_LI="";

		foreach ($module_obj->{'submodules'} as $submodule_name=>$submodule_obj) {

			# Variables
			$submodule_info_code=$submodule_obj->{'code'};
			$submodule_info_name=$submodule_obj->{'name'};
			$submodule_info_fullname=$submodule_obj->{'fullname'}!=""?$submodule_obj->{'fullname'}:$submodule_obj->{'name'};
			$submodule_info_release=$submodule_obj->{'release'}!=""?$submodule_obj->{'release'}:"unknown";
			$submodule_info_description=$submodule_obj->{'description'};
			$submodule_info_available=$submodule_obj->{'available'};


			foreach ($submodule_obj->{'services'} as $service=>$service_infos) {
				$service_code=$service_infos->{'code'};
				$service_name=$service_infos->{'name'};
				$service_fullname=$service_infos->{'fullname'};
				$service_description=$service_infos->{'description'};
				$service_release=($service_infos->{'release'}!="")?$service_infos->{'release'}:"unknown";
				$service_type=$service_infos->{'type'};
				$service_available=$service_infos->{'available'};
				$service_href='';
				$service_href_html='';
				
				if ( ($service_type=="IHM" || $service_type=="WEB")
					&& $service_infos->{'link'}!="" && $service_infos->{'link'}->{'available'}) {
					$service_href=$service_infos->{'href'};
					$service_href_html=" href='$service_href'";
				};

				if ( ($service_type=="IHM" || $service_type=="WEB")
					&& $service_available) {
					$CONTENT_SECTION_MODULES_SERVICES_CONTENT.="
					<p class='mbr-text mbr-fonts-style display-7'>
						+ <a $service_href_html $service_link_target title='[$service_type] $service_description'>
							$service_fullname
						</a> 
					</p>
					";
					$CONTENT_SECTION_MODULES_SERVICES_CONTENT_LI.="
					<li class='mbr-fonts-style'>
						
						<a $service_href_html $service_link_target title='[$submodule_name] [$service_type] [$service_release] $service_description'>
							<b>$service_fullname</b> 
						</a>
						
					</li>
					";
				};

			};

		};

		$CONTENT_SECTION_MODULE_CONTENT.='


					<div class="card p-3 col-12 col-md-4 mb-3">
						<a href="index.modules.php?module='.$module_info_code.'" class="navbar-caption text-secondary ">
							<div class="media mb-0">
								<div class="card-img align-self-center">
									<span class="mbr-iconfont mbri-extension" style="color: rgb(20, 157, 204); fill: rgb(20, 157, 204);"></span>
								</div>
								<h4 class="card-title media-body py-3 mbr-fonts-style display-5">
									'.$module_info_name.'
								</h4>
							</div>
						</a>

						<div class="card-box">
							<br>
							<p class="mbr-text mbr-fonts-style display-7">
								<b>'.$module_info_fullname.'</b>
								<br>
								'.$module_info_description.'
								
								<ul class="" role="tablist" style="list-style-type: circle;">
									'.$CONTENT_SECTION_MODULES_SERVICES_CONTENT_LI.'
								</ul>
							</p>
							
						</div>

					</div>

		';



	};


};



$CONTENT_SECTION_MODULES_CONTENT='
	<section class="features10 cid-ru7OEDbxhA" id="features10-1l">

		<div class="container ">

			<div class="row justify-content-left">

				'.$CONTENT_SECTION_MODULE_CONTENT.'

			</div>

		</div>
		
	</section>

	';



$CONTENT_SECTION_MODULES=$CONTENT_SECTION_MODULES_HEADER.$CONTENT_SECTION_MODULES_CONTENT;



###############
### CONTENT ###
###############

$CONTENT=$CONTENT_SECTION_DASHBOARD.$CONTENT_SECTION_MODULES.$CONTENT_OLD;

echo $CONTENT;


##############
### FOOTER ###
##############

include "footer.inc.php";



?>
