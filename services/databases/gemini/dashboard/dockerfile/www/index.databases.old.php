<?php


#############
### INFOS ###
#############

$APP_SECTION="Database";



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

$db=$_REQUEST["database"];



### DATA
##########

include "databases.inc.php";




### CONTENT
#############



#########################
### SECTION DATABASES ###
#########################


# Init

$CONTENT_SECTION_DATABASES_LIST='';
$CONTENT_SECTION_DATABASE_CONTENT='';

$DATABASE_NUMBER=0;

// echo "<pre>";
// print_r($databases);
// echo "</pre>";


### Check databases

foreach ($databases as $database_name=>$database_obj) {

	$database_infos=$database_obj->{'infos'};
	$database_releases=$database_obj->{'releases'};


	# Variables
	$database_info_code=$database_infos->{'code'};
	$database_info_name=$database_infos->{'name'};
	$database_info_fullname=$database_infos->{'fullname'}!=""?$database_infos->{'fullname'}:$database_infos->{'name'};
	$database_info_release=$database_infos->{'release'}!=""?$database_infos->{'release'}:"unknown";
	$database_info_description=$database_infos->{'description'};
	$database_info_available=$database_infos->{'available'};
	$database_info_website=$database_infos->{'website'};
	
	# Current assembly
	$database_info_assembly="";
	if ($database_obj->{'releases'}->{"current"}->{"assembly"}) {
		if (is_object($database_obj->{'releases'}->{"current"}->{"assembly"})) {
			$database_info_assembly="
				<br>Releases
				";
			foreach ($database_obj->{'releases'}->{"current"}->{"assembly"} as $key=>$value) {
				$database_info_assembly.="<br>&nbsp;&nbsp;&nbsp;$key - $value";
			};
		} elseif (is_array($database_obj->{'releases'}->{"current"}->{"assembly"})) {
			$database_info_assembly="
				<br>"
				."["
				.implode($database_obj->{'releases'}->{"current"}->{"assembly"})
				."]"
			;

		} else {
			$database_info_assembly="
			<br>
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



	# Website
	$database_info_website_url="";
	if ($database_info_website != "") {
		$database_info_website_url="<br><a href='$database_info_website' target='$database_info_code'>WebSite</a>";
	};

	if (1) {
		// $CONTENT_SECTION_DATABASES_LIST.='
		// 	<p class="mbr-text mbr-fonts-style display-7">
		// 		<big><a href="?database='.$database_info_code.'">'.$database_info_name.'</a></big> ['.$database_info_release.'] '.$database_info_fullname.'
		// 		<br>
		// 		'.$database_info_description.'
		// 	</p>
		// ';
		$CONTENT_SECTION_DATABASES_LIST.='
		<li class="nav-item mbr-fonts-style">
			<a class="nav-link  display-7" role="tab" href="?database='.$database_info_code.'">
				<big>'.$database_info_name.'</big>
			</a>
		</li>
		';
		
	}

	# If database input
	if ($database_info_code==$db) {

		// echo "$database_name $database_info_code $database_info_name<br>" ;

		$CONTENT_SECTION_DATABASES_RELEASES_CONTENT="";
		$CONTENT_SECTION_DATABASES_RELEASES_CONTENT_LI_CURRENT="";
		$CONTENT_SECTION_DATABASES_RELEASES_CONTENT_LI_LATEST="";
		$CONTENT_SECTION_DATABASES_RELEASES_CONTENT_LI="";

		foreach ($database_releases as $release=>$release_release_infos) {

			$release_release=$release_release_infos->{'release'};
			$release_date=$release_release_infos->{'date'};


			# Assembly
			if (is_object($release_release_infos->{"assembly"})) {
				$database_info_assembly="
					<br>Assembly: 
					";
				foreach ($release_release_infos->{"assembly"} as $key=>$value) {
					$database_info_assembly.="<br>&nbsp;&nbsp;&nbsp;$key - $value";
				};
			} elseif (is_array($release_release_infos->{"assembly"})) {
				$database_info_assembly="
					<br>Assembly: "
					."["
					.implode($release_release_infos->{"assembly"})
					."]"
				;
	
			} else {
				$database_info_assembly="
				<br>Assembly: 
				[".implode(", ",$release_release_infos->{"assembly"})."]";
			};

			# Current download
			$database_info_download="";
			if ($release_release_infos->{"download"}) {
				if (is_object($release_release_infos->{"download"})) {
					$database_info_download="
						<br>Download:
						<ul>
						";
					foreach ($release_release_infos->{"download"} as $key=>$value) {
						$database_info_download.="<li>$key: $value</li>";
					};
					$database_info_download.="
						
						</ul>
						";
				} elseif (is_array($release_release_infos->{"download"})) {
					$database_info_download="
						<br>Download:
						".implode(array_map(function($value, $key) {
						return $key.' - '.$value.'';
					}, array_values($release_release_infos->{"download"}), array_keys($release_release_infos->{"download"})));
				} else {
					$database_info_download="
					<br>Download:
					[".implode(", ",$release_release_infos->{"download"})."]";
				};
			};

			# If release available
				
			$CONTENT_SECTION_DATABASES_RELEASES_CONTENT_LI_ONE="
				<li class='mbr-text mbr-fonts-style display-7' style='color:gray'>
					<a $release_href_html $release_link_target title='[$release_type] $release_description'>
							<b>$release_fullname</b>
						</a><b>$release</b> 
						<br>
						Release: <I>".$release_release."</I>
						<br>
						Date: <I>".$release_date."</I>
						".$database_info_assembly."
						".$database_info_download."
						<br>
				</li>
				";

			if ($release=="current") {
				$CONTENT_SECTION_DATABASES_RELEASES_CONTENT_LI_CURRENT=$CONTENT_SECTION_DATABASES_RELEASES_CONTENT_LI_ONE;
			} elseif ($release=="latest") {
				$CONTENT_SECTION_DATABASES_RELEASES_CONTENT_LI_LATEST=$CONTENT_SECTION_DATABASES_RELEASES_CONTENT_LI_ONE;
			} else {
				$CONTENT_SECTION_DATABASES_RELEASES_CONTENT_LI.=$CONTENT_SECTION_DATABASES_RELEASES_CONTENT_LI_ONE;
			};


		};


		$CONTENT_SECTION_DATABASE_CONTENT.='


					<div class="card p-3 col-12 col-md-12 mb-3">
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
							<p class="mbr-text mbr-fonts-style display-7">
								<br>
								<b>'.$database_info_fullname.'</b>
								<br>
								'.$database_info_description.'
								'.$database_info_website_url.'
								<br>
								<ul class="" role="tablist" style="list-style-type: square;">
									'.$CONTENT_SECTION_DATABASES_RELEASES_CONTENT_LI_CURRENT.'
									'.$CONTENT_SECTION_DATABASES_RELEASES_CONTENT_LI_LATEST.'
									'.$CONTENT_SECTION_DATABASES_RELEASES_CONTENT_LI.'
								</ul>
								
							</p>
							
						</div>

					</div>

		';

	};


};




$CONTENT_SECTION_DATABASES_CONTENT='
	<section class="features10 cid-rtQc0shhyd" id="features10-1l">

		<div class="container ">

			<h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
			   STARK databases
		   </h2>

			<div class="card p-3 col-12 col-md-12 mb-0">

				<div class="container-fluid col-md-12">
					<div class="row tabcont">
						<ul class="nav nav-tabs pt-0 mt-0" role="tablist">
							'.$CONTENT_SECTION_DATABASES_LIST.' 
						</ul>
					</div>
				</div>

			</div>
		
		</div>
		
	</section>


	<section class=" cid-ru7OEDbxhA" id="truc"> 

		<div class="container">

			<div class="card p-3 col-12 col-md-12 mb-0 row justify-content-left">

				'.$CONTENT_SECTION_DATABASE_CONTENT.'

			</div>

		</div>
		
	</section>

		

	';



$CONTENT_SECTION_DATABASES=$CONTENT_SECTION_DATABASES_CONTENT;




###############
### CONTENT ###
###############

$CONTENT=$CONTENT_SECTION_DASHBOARD.$CONTENT_SECTION_DATABASES.$CONTENT_OLD;

echo $CONTENT;


##############
### FOOTER ###
##############

include "footer.inc.php";



?>
