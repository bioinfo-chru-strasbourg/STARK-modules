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
	$DEJAVU_RELEASE="";
};
$default_release=$DEJAVU_RELEASE;

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

		#$default_release="";

		if ($default_release == "") {

			if (isset($database_obj->{'releases'}->{"current"})) {
				$default_release="current";
			} else if (isset($database_obj->{'releases'}->{"latest"})) {
				$default_release="latest";
			} else {
				foreach ($database_obj->{'releases'} as $database_obj_key=>$database_obj_value) {
					if ($default_release == "") {
						$default_release=$database_obj_key;
					};
				};
			};

		};


		$database_release_list="
			<!--
			<label for='release'>Select a release</label>
			<br>
			-->
			<select name='release' id='release' form='releaseform' onchange='this.form.submit()'>
		";
		foreach ($database_obj->{'releases'} as $database_obj_key=>$database_obj_value) {
				if ($database_obj_key==$default_release) {
					$database_release_selected=" selected ";
				} else {
					$database_release_selected="";
				};
			$database_release_list.="<option value='".$database_obj_key."' ".$database_release_selected.">".$database_obj_key."</option>";
		};
		$database_release_list.="
			</select>
		";	

		// echo "<pre>";
		// print_r($database_obj);
		// echo "</pre>";


		if ($database_obj->{'releases'}->{$default_release}->{"release"}) {
			if (is_object($database_obj->{'releases'}->{$default_release}->{"release"})) {
				$database_info_current="
					<br>Releases
					";
				foreach ($database_obj->{'releases'}->{$default_release}->{"release"} as $key=>$value) {
					$database_info_current.="<br>&nbsp;&nbsp;&nbsp;$key - $value";
				};
			} elseif (is_array($database_obj->{'releases'}->{$default_release}->{"release"})) {
				$database_info_current="
					<br>Releases"
					."["
					.implode($database_obj->{'releases'}->{$default_release}->{"release"})
					."]"
					;
			} else {
				$database_info_current="
					<br>Release:
					<i>".$database_obj->{'releases'}->{$default_release}->{'release'}."</i>"
				;
			}
		};

		# Current Date
		$database_info_current_date="";
		if ($database_obj->{'releases'}->{$default_release}->{"date"}) {
			if (is_object($database_obj->{'releases'}->{$default_release}->{"date"})) {
				$database_info_current_date="
					<br>Dates
					";
				foreach ($database_obj->{'releases'}->{$default_release}->{"date"} as $key=>$value) {
					$database_info_current_date.="<br>&nbsp;&nbsp;&nbsp;$key - $value";
				};
			} elseif (is_array($database_obj->{'releases'}->{$default_release}->{"date"})) {
				$database_info_current_date="
					<br>Dates "
					."["
					.implode($database_obj->{'releases'}->{$default_release}->{"date"})
					."]"
					;
			} else {
				$database_info_current_date="
					<br>Date:
					<i>".$database_obj->{'releases'}->{$default_release}->{'date'}."</i>"
				;
			}
		};



		# Current assembly
		$database_info_assembly="";
		if ($database_obj->{'releases'}->{$default_release}->{"assembly"}) {
			if (is_object($database_obj->{'releases'}->{$default_release}->{"assembly"})) {
				$database_info_assembly="
					<br>Assembly: 
					";
				foreach ($database_obj->{'releases'}->{$default_release}->{"assembly"} as $key=>$value) {
					$database_info_assembly.="<br>&nbsp;&nbsp;&nbsp;$key - $value";
				};
			} elseif (is_array($database_obj->{'releases'}->{$default_release}->{"assembly"})) {
				$database_info_assembly="
					<br>Assembly: "
					."["
					.implode($database_obj->{'releases'}->{$default_release}->{"assembly"})
					."]"
				;

			} else {
				$database_info_assembly="
				<br>Assembly: 
				[".implode(", ",$database_obj->{'releases'}->{$default_release}->{"assembly"})."]";
			};
		};

		# Current download
		$database_info_download="";
		if ($database_obj->{'releases'}->{$default_release}->{"download"}) {
			if (is_object($database_obj->{'releases'}->{$default_release}->{"download"})) {
				$database_info_download="
					<br>Download
					";
				foreach ($database_obj->{'releases'}->{$default_release}->{"download"} as $key=>$value) {
					$database_info_download.="<br>&nbsp;&nbsp;&nbsp;$key - $value";
				};
			} elseif (is_array($database_obj->{'releases'}->{$default_release}->{"download"})) {
				$database_info_download="
					<br>
					".implode(array_map(function($value, $key) {
				    return $key.' - '.$value.'';
				}, array_values($database_obj->{'releases'}->{$default_release}->{"download"}), array_keys($database_obj->{'releases'}->{$default_release}->{"download"})));
			} else {
				$database_info_download="
				<br>
				[".implode(", ",$database_obj->{'releases'}->{$default_release}->{"download"})."]";
			};
		};

		#echo $folder_databases."/".$folder_databases_subfolder_dejavu."/".$default_release."<br>";

		foreach (glob ($folder_databases."/".$folder_databases_subfolder_dejavu."/".$default_release."/vcf/*/*" ) as $key => $group_project) {

			# group/project
			$project=basename($group_project);
			$group=basename(dirname($group_project));

			# Links
			

			# Stats
			$stats_txt=glob($group_project."/stats/dejavu.stats.txt")[0];
			$stats_tsv=glob($group_project."/stats/dejavu.stats.tsv")[0];
			$link_download_stats='<a href="'.$stats_txt.'" download>TXT</a> <a href="'.$stats_txt.'" download>TSV</a>';

			# download links
			$link_download="";
			foreach (glob ($group_project."/*{tsv,vcf.gz}",GLOB_BRACE) as $key_file => $file) {
				$link_download.=' &nbsp;&nbsp;<a href="'.$file.'" download>'.basename($file).'</a>&nbsp;&nbsp; ';
			};

			# Counts
			$nb_variants=trim(join("",file(glob($group_project."/stats/dejavu.stats.nb_variants")[0])));
			$nb_samples=trim(join("",file(glob($group_project."/stats/dejavu.stats.nb_samples")[0])));


			# Plugins
			$PLUGINS_LINKS=array();
			foreach ($plugins_obj as $plugin=>$plugin_infos) {

				if (is_file($plugins_folder."/".$plugin_infos->{"script"})						# script exists
					&& is_dir($group_project."/".$plugin_infos->{"folder_data"})					# folder data exists
					&& glob($group_project."/{".$plugin_infos->{"files_data"}."}",GLOB_BRACE)		# files data exists (at least one)$plugin_infos->{"available"}
					&& $modules_obj_array[$plugin_infos->{"module"}]->{"available"}				# module available
					&& $modules_obj_array[$plugin_infos->{"module"}]->{"submodules"}->{$plugin_infos->{"submodule"}}->{"services"}->{$plugin_infos->{"service"}}->{"available"}		# service 
					&& is_file($plugins_folder."/".$plugin_infos->{"script"})
					) {
						foreach ($plugin_infos->{"level"} as $level=>$level_infos) {
							$url_full=$plugins_folder."/".$plugin_infos->{"script"}."?".$plugin_infos->{"parameters"};
							$url_full.=($level_infos->{"parameters"}!="")?"&".$level_infos->{"parameters"}:"";
							$url_full.="&PATH[]=".$group_project;
							$link="<a target='".$plugin_infos->{"target"}."' href='$url_full'>".$plugin_infos->{"label"}."</a>";
							$PLUGINS_LINKS[$level].=" ".$link." ";
						};
				}
				
			}



			$CONTENT_SECTION_DATABASE_CONTENT.='

				<div class="card p-3 col-12 col-md-4 mb-3">
					<a href="index.databases.php?database='.$database_info_code.'" class="navbar-caption text-secondary ">
						<div class="media mb-0">
							<div class="card-img align-self-center">
								<span class="mbr-iconfont mbri-database" style="color: rgb(20, 157, 204); fill: rgb(20, 157, 204);"></span>
							</div>
							<h3 class="card-title media-body py-3 mbr-fonts-style display-7">
								'.$group.'
								<br>
								'.$project.'
							</h3>
						</div>
					</a>

					<div class="card-box">
						<br>
						<p class="mbr-text mbr-fonts-style display-7">
						<b>'.$nb_variants.'</b> variants within <b>'.$nb_samples.'</b> samples
						<br>
						Statistics '.$link_download_stats.'
						<br>
						
						Plugins '.$PLUGINS_LINKS["project"].'
						<br>
						Files '.$link_download.'
							
						</p>
						
					</div>

				</div>

			';

			
		};


	};


	$CONTENT_SECTION_DATABASES_HEADER='
		<section class="counters1 counters cid-ru7OEH6r7i" id="SECTION_ID" >

			<div class="container">

			<h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
				'.$database_info_fullname.'
			</h2>
			
			<p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
				<span class="head-item mbr-fonts-style display-7">
					'.$database_release_list.'
					<br>
				</span>
				'.$database_info_current.'
				'.$database_info_current_date.'
				'.$database_info_assembly.'
				<br>
				
			</p>
			
			</div>
			<form action="" id="releaseform" method="POST">
			</form>
		</section>
		';




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
