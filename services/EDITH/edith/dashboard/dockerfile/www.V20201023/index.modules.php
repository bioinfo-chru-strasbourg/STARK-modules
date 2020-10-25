<?php


#############
### INFOS ###
#############

$APP_SECTION="Modules";



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

$module=$_REQUEST["module"];


### CONTENT
#############



#######################
### SECTION MODULES ###
#######################


# Init

$CONTENT_SECTION_MODULES_LIST='';
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

	if ($module_info_available) {
		// $CONTENT_SECTION_MODULES_LIST.='
		// 	<p class="mbr-text mbr-fonts-style display-7">
		// 		<big><a href="?module='.$module_info_code.'">'.$module_info_name.'</a></big> ['.$module_info_release.'] '.$module_info_fullname.'
		// 		<br>
		// 		'.$module_info_description.'
		// 	</p>
		// ';
		$CONTENT_SECTION_MODULES_LIST.='
		<li class="nav-item mbr-fonts-style">
			<a class="nav-link  display-7" role="tab" href="?module='.$module_info_code.'">
				<big>'.$module_info_name.'</big>
			</a>
		</li>
		';
		
	}

	# If module input
	if ($module_info_code==$module) {

		$CONTENT_SECTION_MODULES_SERVICES_CONTENT="";
		$CONTENT_SECTION_MODULES_SERVICES_CONTENT_LI="";

		foreach ($module_obj->{'services'} as $service=>$service_infos) {
			$service_code=$service_infos->{'code'};
			$service_name=$service_infos->{'name'};
			$service_fullname=$service_infos->{'fullname'};
			$service_description=$service_infos->{'description'};
			$service_release=$service_infos->{'release'};
			$service_type=$service_infos->{'type'};
			$service_available=$service_infos->{'available'};
			$service_href='';
			$service_href_html='';

			if ($service_type=="IHM" && $service_infos->{'link'}!="" && $service_infos->{'link'}->{'available'}) {
				$service_href=$service_infos->{'href'};
				$service_href_html=" href='$service_href'";
			};

			# Service Status
			
			$service_state=$docker_containers_obj_array[$service_infos->{'container'}]["obj_infos"]->{"State"};
			$service_status=$docker_containers_obj_array[$service_infos->{'container'}]["obj_infos"]->{"Status"};

			$service_state=$service_state==""?"unavailable":$service_state;
			$service_status=$service_status==""?"no information":$service_status;

			// echo "<pre>";
			// print_r($docker_containers_obj_array[$service_infos->{'container'}]["obj_infos"]);
			// echo "</pre>";

			$service_status_color="gray";
			switch ($service_state) {
			    case "running":
			        $service_status_color="green";
			        break;
			    case "restarting":
			        $service_status_color="orange";
			        break;
			    case "paused":
			        $service_status_color="orange";
			        break;
			    case "exited":
			        $service_status_color="red";
			        break;
			    case "dead":
			        $service_status_color="red";
			        break;
			    case "created":
			        $service_status_color="blue";
			        break;
			    case "unavailable":
			        $service_status_color="gray";
			        break;
			}

			$service_state_status="";
			if ($docker_containers_connected) {
				if ($service_state!="") { #  && $service_status!=""
					$service_state_status="<br><span style='color:$service_status_color'>$service_state</span> <small>[$service_status]</small>";
				};
			} else {
				$service_state_status="";
			};
			#$service_state_status="<br>".$service_state_status;

			# If service available
			if ($service_available) {
				$CONTENT_SECTION_MODULES_SERVICES_CONTENT.="
					<p class='mbr-text mbr-fonts-style display-7'>
						+ <a $service_href_html $service_link_target title='[$service_type] $service_description'>
							<b>$service_fullname</b>
						</a> [$service_release] - $service_type $service_state_status
						<br>
						$service_description
					</p> 
					";
				$CONTENT_SECTION_MODULES_SERVICES_CONTENT_LI.="
					<li class='mbr-text mbr-fonts-style display-7' style='color:gray'>
						<a $service_href_html $service_link_target title='[$service_type] $service_description'>
								<b>$service_fullname</b>
							</a> [$service_release] - $service_type 
							<br>
							$service_description
							$service_state_status
							<br><br>
					</li>
					";
			};

		};


		$CONTENT_SECTION_MODULE_CONTENT.='


					<div class="card p-3 col-12 col-md-12 mb-3">
						<a href="index.modules.php?module='.$module_info_code.'" class="navbar-caption text-secondary ">
							<div class="media mb-0">
								<div class="card-img align-self-center">
									<span class="mbr-iconfont mbri-extension" style="color: rgb(20, 157, 204); fill: rgb(20, 157, 204);"></span>
								</div>
								<h4 class="card-title media-body py-3 mbr-fonts-style display-5">
									'.$module_info_name.'
									<span class="mbr-text mbr-fonts-style display-7">['.$module_info_release.']</span>
								</h4>
							</div>
						</a>

						<div class="card-box">
							<p class="mbr-text mbr-fonts-style display-7">
								<br>
								<b>'.$module_info_fullname.'</b>
								<br>
								'.$module_info_description.'
								<br>
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
	<section class="features10 cid-rtQc0shhyd" id="features10-1l">

		<div class="container ">

			<h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
			   STARK Modules & Services
		   </h2>

			<div class="card p-3 col-12 col-md-12 mb-0">

				<div class="container-fluid col-md-12">
					<div class="row tabcont">
						<ul class="nav nav-tabs pt-0 mt-0" role="tablist">
							'.$CONTENT_SECTION_MODULES_LIST.' 
						</ul>
					</div>
				</div>

			</div>
		
		</div>
		
	</section>


	<section class=" cid-ru7OEDbxhA" id="truc"> 

		<div class="container">

			<div class="card p-3 col-12 col-md-12 mb-0 row justify-content-left">

				'.$CONTENT_SECTION_MODULE_CONTENT.'

			</div>

		</div>
		
	</section>

		

	';



$CONTENT_SECTION_MODULES=$CONTENT_SECTION_MODULES_CONTENT;




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
