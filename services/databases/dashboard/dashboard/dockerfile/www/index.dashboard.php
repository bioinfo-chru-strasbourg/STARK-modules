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


$DEJAVU_RELEASE=$_REQUEST["release"];
if ($DEJAVU_RELEASE=="") {
	$DEJAVU_RELEASE="latest";
};


### DATA
##########

include "databases.inc.php";



### CONTENT
#############


#########################
### SECTION DASHBOARD ###
#########################






#########################
### SECTION DATABASES ###
#########################

if (1) {

	$CONTENT_SECTION_DATABASES_HEADER='
	<section class="counters1 counters cid-ru7OEH6r7i" id="SECTION_ID" >

		<div class="container">

		   <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
			   STARK Databases
		   </h2>
		   <!--
		   <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
			   Statistics for inputs and repositories folders (duplication not considered)
		   </p>
		   -->
		</div>

	</section>
	';

			
	$CONTENT_SECTION_DATABASE_CONTENT='';

	$DATABASE_NUMBER=0;


	### Check modules


	foreach ($databases as $database=>$database_obj) {
		
		// echo "<pre>";
		// print_r($database_obj);
		// echo "</pre>";

		# Variables
		$database_info_code=$database_obj->{"infos"}->{'code'};
		$database_info_name=$database_obj->{"infos"}->{'name'};
		$database_info_fullname=$database_obj->{"infos"}->{'fullname'}!=""?$database_obj->{"infos"}->{'fullname'}:$database_obj->{"infos"}->{'name'};
		$database_info_release=$database_obj->{"infos"}->{'release'}!=""?$database_obj->{"infos"}->{'release'}:"unknown";
		$database_info_description=$database_obj->{"infos"}->{'description'};
		$database_info_website=$database_obj->{"infos"}->{'website'};
		$database_info_available=$database_obj->{"infos"}->{'available'};


		# Website
		$database_info_website_url="";
		if ($database_info_website != "") {
			$database_info_website_url="<br><a href='$database_info_website' target='$database_info_code'>WebSite</a>";
		};


		
		# Current Release
		$database_info_current="";
		if ($database_obj->{'releases'}->{"current"}->{"release"}) {
			if (is_object($database_obj->{'releases'}->{"current"}->{"release"})) {
				$database_info_current="
					<br>Releases
					";
				foreach ($database_obj->{'releases'}->{"current"}->{"release"} as $key=>$value) {
					$database_info_current.="<br>&nbsp;&nbsp;&nbsp;$key - $value";
				};
			} elseif (is_array($database_obj->{'releases'}->{"current"}->{"release"})) {
				$database_info_current="
					<br>Releases"
					."["
					.implode($database_obj->{'releases'}->{"current"}->{"release"})
					."]"
					;
			} else {
				$database_info_current="
					<br>Release:
					<i>".$database_obj->{'releases'}->{'current'}->{'release'}."</i>"
				;
			}
		};

		# Current Date
		$database_info_current_date="";
		if ($database_obj->{'releases'}->{"current"}->{"date"}) {
			if (is_object($database_obj->{'releases'}->{"current"}->{"date"})) {
				$database_info_current_date="
					<br>Dates
					";
				foreach ($database_obj->{'releases'}->{"current"}->{"date"} as $key=>$value) {
					$database_info_current_date.="<br>&nbsp;&nbsp;&nbsp;$key - $value";
				};
			} elseif (is_array($database_obj->{'releases'}->{"current"}->{"date"})) {
				$database_info_current_date="
					<br>Dates "
					."["
					.implode($database_obj->{'releases'}->{"current"}->{"date"})
					."]"
					;
			} else {
				$database_info_current_date="
					<br>Date:
					<i>".$database_obj->{'releases'}->{'current'}->{'date'}."</i>"
				;
			}
		};



		# Current assembly
		$database_info_assembly="";
		if ($database_obj->{'releases'}->{"current"}->{"assembly"}) {
			if (is_object($database_obj->{'releases'}->{"current"}->{"assembly"})) {
				$database_info_assembly="
					<br>Assembly: 
					";
				foreach ($database_obj->{'releases'}->{"current"}->{"assembly"} as $key=>$value) {
					$database_info_assembly.="<br>&nbsp;&nbsp;&nbsp;$key - $value";
				};
			} elseif (is_array($database_obj->{'releases'}->{"current"}->{"assembly"})) {
				$database_info_assembly="
					<br>Assembly: "
					."["
					.implode($database_obj->{'releases'}->{"current"}->{"assembly"})
					."]"
				;

			} else {
				$database_info_assembly="
				<br>Assembly: 
				[".implode(", ",$database_obj->{'releases'}->{"current"}->{"assembly"})."]";
			};
		};

		# Current download
		$database_info_download="";
		if ($database_obj->{'releases'}->{"current"}->{"download"}) {
			if (is_object($database_obj->{'releases'}->{"current"}->{"download"})) {
				$database_info_download="
					<br>Download
					";
				foreach ($database_obj->{'releases'}->{"current"}->{"download"} as $key=>$value) {
					$database_info_download.="<br>&nbsp;&nbsp;&nbsp;$key - $value";
				};
			} elseif (is_array($database_obj->{'releases'}->{"current"}->{"download"})) {
				$database_info_download="
					<br>
					".implode(array_map(function($value, $key) {
				    return $key.' - '.$value.'';
				}, array_values($database_obj->{'releases'}->{"current"}->{"download"}), array_keys($database_obj->{'releases'}->{"current"}->{"download"})));
			} else {
				$database_info_download="
				<br>
				[".implode(", ",$database_obj->{'releases'}->{"current"}->{"download"})."]";
			};
		};

		$CONTENT_SECTION_DATABASE_CONTENT.='

			<div class="card p-3 col-12 col-md-4 mb-3">
				<a href="index.databases.php?database='.$database_info_code.'" class="navbar-caption text-secondary ">
					<div class="media mb-0">
						<div class="card-img align-self-center">
							<span class="mbr-iconfont mbri-database" style="color: rgb(20, 157, 204); fill: rgb(20, 157, 204);"></span>
						</div>
						<h4 class="card-title media-body py-3 mbr-fonts-style display-5">
							'.$database_info_name.'
						</h4>
					</div>
				</a>

				<div class="card-box">
					<br>
					<p class="mbr-text mbr-fonts-style display-7">
						<b>'.$database_info_fullname.'</b>
						'.$database_info_current.'
						'.$database_info_current_date.'
						'.$database_info_assembly.'
						<!--'.$database_info_download.'-->
						<br>'.$database_info_description.'
						'.$database_info_website_url.'
						
					</p>
					
				</div>

			</div>

		';

	};



	$CONTENT_SECTION_DATABASES_CONTENT='
		<section class="features10 cid-ru7OEDbxhA" id="features10-1l">

			<div class="container ">

				<div class="row justify-content-left">

					'.$CONTENT_SECTION_DATABASE_CONTENT.'

				</div>

			</div>
			
		</section>

		';



	$CONTENT_SECTION_DATABASES=$CONTENT_SECTION_DATABASES_HEADER.$CONTENT_SECTION_DATABASES_CONTENT;


};


if (0) {
	# 
	#echo $DEJAVU_RELEASE;

	$folder_dejavu=$folder_databases."/".$folder_databases_subfolder_dejavu;

	print_r(glob ($folder_dejavu."/$DEJAVU_RELEASE/*" ));

	foreach (glob ($folder_dejavu."/$DEJAVU_RELEASE/dejavu.*.done" ) as $key => $database_path) {

		echo "<pre>";
		#print_r($database_path);

		# Search DEJAVU databases
		
		$database_file_done=end(explode ( "/" , $database_path ));
		#echo "<br>".$database_file_done;

		# Search DEJAVU Gorup / Project

		$database_file_done_split=explode ( "." , str_replace("dejavu.","",str_replace(".done","",$database_file_done)) );
		$database_group=$database_file_done_split[0];
		$database_project=$database_file_done_split[1];
		echo "<br>Group: ".$database_group;
		echo "<br>Project: ".$database_project;
		echo "<br>".implode("",file($database_path));



		# Search DEJAVU assembly & stats
		$database_file_stats=glob ($folder_dejavu."/$DEJAVU_RELEASE/*dejavu.$database_group.$database_project.stats" )[0];
		$database_assembly=explode ( "_dejavu." , end(explode ( "/" , $database_file_stats )) )[0];
		echo "<br>Assembly ".$database_assembly;
		echo "<br>Statistics: ".$database_file_stats;
		echo "<br>".implode("",file($database_file_stats));

		# Search DEJAVU stats
		$database_file_stats=glob ($folder_dejavu."/$DEJAVU_RELEASE/*dejavu.$database_group.$database_project.txt" )[0];
		echo "<br>".$database_file_stats;
		#echo "<br>".implode("",file($database_file_stats));

		echo "</pre>";
	};


	### DATA

	$align="center";


	# STATISTICS

	$STATISTICS_DIV="";

	foreach ($input_list as $input=>$hash) {
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
					<b>$activity</b>
					<small>
		";
		foreach ($status_list as $status=>$nb) {
			if (!in_array($status,$status_exclude) ) {
				$ACTIVITY_DIV.="
						<br>
						<b>$nb</b> <span style='color:".$ACTIVITY_STATUS_COLOR[$activity][$status]."'>$status</span> 
				";
			};
		};
		$ACTIVITY_DIV.="
				</small>
				</p>
			</div>
		";
	};





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
					<a href="index.reports.php" class="navbar-caption text-secondary ">
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
								<b>Search a analysis/sample in repositories</b>
								<form action="index.reports.php" method="POST" >
									<span class="head-item mbr-fonts-style display-7">
									<input name="analysis" value="" style="width:150px" class=""></input>
									<input name="sample" value="" style="width:100px" class=""></input>
									<input type="Submit" value="Search" class=""></input>
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

			foreach ($module_obj->{'services'} as $service=>$service_infos) {
				$service_code=$service_infos->{'code'};
				$service_name=$service_infos->{'name'};
				$service_fullname=$service_infos->{'fullname'};
				$service_description=$service_infos->{'description'};
				$service_type=$service_infos->{'type'};
				$service_available=$service_infos->{'available'};
				$service_href='';
				
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
						
						<a $service_href_html $service_link_target title='[$service_type] $service_description'>
							<b>$service_fullname</b> 
						</a>
						
					</li>
					 ";
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


};#if

###############
### CONTENT ###
###############

$CONTENT=$CONTENT_SECTION_DASHBOARD.$CONTENT_SECTION_DATABASES;

echo $CONTENT;


##############
### FOOTER ###
##############

include "footer.inc.php";



?>
