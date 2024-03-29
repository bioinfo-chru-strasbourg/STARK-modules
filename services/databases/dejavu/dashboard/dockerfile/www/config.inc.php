<?php

### CONFIG
############


### APP
$APP_CODE="DejaVu";
$APP_NAME="DejaVu IHM";
$APP_RELEASE="1.0";
$APP_DATE="20201117";
$APP_DESCRIPTION="Dashboard for DejaVu databases tracking and exploration";
$APP_COPYRIGHT="HUS/CPS";
$APP_LICENCE="GPLA - GNU-GPL Affero";
$APP_AUTHORS="A. Le Béchec";

$APP_LOGO="assets/logo.png";
$APP_FAVICON="assets/favicon.ico";


### SERVER
$APP_SERVER_NAME=$_SERVER["SERVER_NAME"]!=""?$_SERVER["SERVER_NAME"]:"localhost";


### ENV

# DEBUG
$DEBUG=$_REQUEST["DEBUG"];





### PLUGINS

$plugins_folder=".";
$plugins_file="config/plugins.conf";
$plugins_json=implode("",file($plugins_file));
$plugins_obj=json_decode($plugins_json);




### FOLDERS


# SERVICES folder
if ($_ENV["FOLDER_SERVICES"] != "") {
	$folder_services=$_ENV["FOLDER_SERVICES"];
} else {
	$folder_services="services";
};

# DOCKER_STARK_PREFIX
if ($_ENV["DOCKER_STARK_PREFIX"] != "") {
	$docker_stark_prefix=$_ENV["DOCKER_STARK_PREFIX"];
} else {
	$docker_stark_prefix="stark";
};



# DATABASES folder
if ($_ENV["FOLDER_DATABASES"] != "") {
	$folder_databases=$_ENV["FOLDER_DATABASES"];
} else {
	$folder_databases="databases";
};


# FOLDER_DATABASES_SUBFOLDER_DEJAVU folder
if ($_ENV["FOLDER_DATABASES_SUBFOLDER_DEJAVU"] != "") {
	$folder_databases_subfolder_dejavu=$_ENV["FOLDER_DATABASES_SUBFOLDER_DEJAVU"];
} else {
	$folder_databases_subfolder_dejavu="dejavu";
};


# IGV
if ($_ENV["FOLDER_IGV"] != "") {
	$folder_igv=$_ENV["FOLDER_IGV"];
} else {
	$folder_igv=$folder_services."/genomebrowser/igv/igv";
};

# Create folder
if (!is_dir($folder_igv)) {
	if (!mkdir($folder_igv, 0777, true)) {
		die("Failed to create folders '$folder_igv'...");
	};
};



### MODULES

# Find modules configurations
$modules = array_filter(glob($folder_services.'/*/STARK*.module', GLOB_BRACE), 'is_file');
$submodules = array_filter(glob($folder_services.'/*/*/STARK*.module', GLOB_BRACE), 'is_file');


// echo "<pre>";

$module_folders=array();

foreach ($modules as $module_file) {
	
	$module_folder_name=basename(dirname($module_file));
	if (!array_key_exists($module_folder_name,$module_folders)) {
		$module_folders[$module_folder_name]=array();
	};
	
	$module_folders[$module_folder_name]=array_merge_recursive($module_folders[$module_folder_name],json_decode(implode("",file($module_file)), true));

}

foreach ($submodules as $module_file) {
	
	$module_folder_name=basename(dirname(dirname($module_file)));
	$submodule_folder_name=basename(dirname($module_file));
	if (!array_key_exists($module_folder_name,$module_folders)) {
		$module_folders[$module_folder_name]=array();
	};
	
	$module_folders[$module_folder_name]=array_merge_recursive($module_folders[$module_folder_name],json_decode(implode("",file($module_file)), true));

}

// print_r($module_folders);
// echo "</pre>";


# init
$modules_obj_array;

# Foeach modules
#foreach ($modules as $module_file) {
foreach ($module_folders as $module_folder=>$module_array) {
	
	# Read JSON file
	#$module_json=implode("",file($module_file));
	$module_json=json_encode($module_array);
	$module_obj=json_decode($module_json);
	
	# Module code
	$module_info_code=$module_obj->{'code'};

	$service_infos=null;

	foreach ($module_obj->{'submodules'} as $submodule_info_code=>$submodule_obj) {

		# Services array
		foreach ($submodule_obj->{'services'} as $service=>$service_infos) {
		
			# Service code
			$service_code=$service_infos->{'code'}!=""?$service_infos->{'code'}:$service;
			$service_infos->{'code'}=$service_code;

			# Name
			$service_name=$service_infos->{'name'}!=""?$service_infos->{'name'}:$service_infos->{'code'};
			$service_infos->{'name'}=$service_name;

			# Full Name
			$service_fullname=$service_infos->{'fullname'}!=""?$service_infos->{'fullname'}:$service_infos->{'name'};
			$service_infos->{'fullname'}=$service_fullname;

			# Description
			$service_description=$service_infos->{'description'}!=""?$service_infos->{'description'}:$service_infos->{'fullname'};
			$service_infos->{'description'}=$service_description;

			# Type
			$service_type=$service_infos->{'type'}!=""?$service_infos->{'type'}:"unknown";
			$service_infos->{'type'}=$service_type;

			# SubModule
			#$service_submodule=$service_infos->{'submodule'}!=""?$service_infos->{'submodule'}:"main";
			#$service_infos->{'submodule'}=$service_submodule;

			# Available
			$service_available=$service_infos->{'available'}!=""?$service_infos->{'available'}:false;
			$service_infos->{'available'}=$service_available;
		
			# Service container
			$service_container=$service_infos->{'container'}!=""?$service_infos->{'container'}:"$docker_stark_prefix-module-$module_info_code-submodule-$submodule_info_code-service-$service_code";
			$service_infos->{'container'}=$service_container;

			# Service HREF
			if ($service_infos->{'link'}!="") {
			
				# protocol
				$service_link_protocol=$service_infos->{'link'}->{'protocol'}!=""?$service_infos->{'link'}->{'protocol'}:"http";
				$service_infos->{'link'}->{'protocol'}=$service_link_protocol;
			
				# IP
				$service_link_ip=$service_infos->{'link'}->{'ip'}!=""?$service_infos->{'link'}->{'ip'}:$APP_SERVER_NAME;
				$service_infos->{'link'}->{'ip'}=$service_link_ip;
			
				# port Host
				$service_link_port=$service_infos->{'link'}->{'port'}!=""?":".$service_infos->{'link'}->{'port'}:"";
				$service_link_port=($service_link_port=="" && $docker_containers_obj_array[$service_container]["obj_infos"]->{"Ports"}[0]->{"PublicPort"}!="")?$docker_containers_obj_array[$service_container]["obj_infos"]->{"Ports"}[0]->{"PublicPort"}:$service_link_port;
				$service_infos->{'link'}->{'port'}=$service_link_port;

				# Port Inner
				$service_link_port_inner=$service_infos->{'link'}->{'port_inner'}!=""?":".$service_infos->{'link'}->{'port_inner'}:"";
				$service_link_port_inner=($service_link_port_inner=="" && $docker_containers_obj_array[$service_container]["obj_infos"]->{"Ports"}[0]->{"PrivatePort"}!="")?$docker_containers_obj_array[$service_container]["obj_infos"]->{"Ports"}[0]->{"PrivatePort"}:$service_link_port_inner;
				$service_link_port_inner=($service_link_port_inner=="")?$service_link_port:$service_link_port_inner;
				$service_infos->{'link'}->{'port_inner'}=$service_link_port_inner;
				
				# path
				$service_link_path=$service_infos->{'link'}->{'path'}!=""?$service_infos->{'link'}->{'path'}:"";
				$service_infos->{'link'}->{'path'}=$service_link_path;
				
				# target
				$service_link_target=$service_infos->{'link'}->{'target'}!=""?$service_infos->{'link'}->{'path'}:"";
				$service_infos->{'link'}->{'target'}=$service_link_target;
		
				# href
				$service_href="$service_link_protocol://$service_link_ip$service_link_port/$service_link_path";
				$service_href_inner="$service_link_protocol://$service_container$service_link_port_inner/$service_link_path";
				$service_infos->{'href'}=$service_href;
				$service_infos->{'href_inner'}=$service_href_inner;

			};

			# Module array
			#$module_obj->{'services'}->{$service}=$service_infos;
			$module_obj->{'submodules'}->{$submodule_info_code}->{'services'}->{$service}=$service_infos;

		};

	};

	# Module array
	$modules_obj_array[$module_info_code]=$module_obj;
};




### URI


# API
if ($_ENV["URINNER_API"] != "") {
	$urinner_api=$_ENV["URINNER_API"];
#} elseif ($modules_obj_array["STARK"]->{"services"}->{"API"}->{"href_inner"}!="") {
} elseif ($modules_obj_array["stark"]->{"submodules"}->{"stark"}->{"services"}->{"api"}->{"href_inner"}!="") {
	$urinner_api=$modules_obj_array["stark"]->{"submodules"}->{"stark"}->{"services"}->{"api"}->{"href_inner"};
} elseif ($modules_obj_array["stark"]->{"services"}->{"api"}->{"href_inner"}!="") {
	$urinner_api=$modules_obj_array["stark"]->{"services"}->{"api"}->{"href_inner"};
} else {
	$urinner_api="";
};

$urinner_api_list=$urinner_api."queue?list";
$urinner_api_info=$urinner_api."queue?info=";
$urinner_api_log=$urinner_api."queue?log=";


# DAS
if ($_ENV["URINNER_DAS"] != "") {
	$urinner_das=$_ENV["URINNER_DAS"];
} elseif ($modules_obj_array["stark"]->{"submodules"}->{"stark"}->{"services"}->{"das"}->{"href_inner"}!="") {
	$urinner_das=$modules_obj_array["stark"]->{"submodules"}->{"stark"}->{"services"}->{"das"}->{"href_inner"};
} elseif ($modules_obj_array["stark"]->{"services"}->{"das"}->{"href_inner"}!="") {
	$urinner_das=$modules_obj_array["stark"]->{"services"}->{"das"}->{"href_inner"};
} else {
	$urinner_das="";
};

# DAS
if ($_ENV["URI_DAS"] != "") {
	$uri_das=$_ENV["URI_DAS"];
} elseif ($modules_obj_array["stark"]->{"submodules"}->{"stark"}->{"services"}->{"das"}->{"href"}!="") {
	$uri_das=$modules_obj_array["stark"]->{"submodules"}->{"stark"}->{"services"}->{"das"}->{"href"};
} elseif ($modules_obj_array["stark"]->{"services"}->{"das"}->{"href"}!="") {
	$uri_das=$modules_obj_array["stark"]->{"services"}->{"das"}->{"href"};
} else {
	$uri_das="";
};


# IGV
if ($_ENV["URINNER_IGV"] != "") {
	$urinner_igv=$_ENV["URINNER_IGV"];
} elseif ($modules_obj_array["genomebrowser"]->{"submodules"}->{"igv"}->{"services"}->{"igv"}->{"href_inner"}!="") {
	$urinner_igv=$modules_obj_array["genomebrowser"]->{"submodules"}->{"igv"}->{"services"}->{"igv"}->{"href_inner"};
} elseif ($modules_obj_array["genomebrowser"]->{"services"}->{"igv"}->{"href_inner"}!="") {
	$urinner_igv=$modules_obj_array["genomebrowser"]->{"services"}->{"igv"}->{"href_inner"};
} else {
	$urinner_igv="";
};

# IGV
if ($_ENV["URI_IGV"] != "") {
	$uri_igv=$_ENV["URI_IGV"];
} elseif ($modules_obj_array["genomebrowser"]->{"submodules"}->{"igv"}->{"services"}->{"igv"}->{"href"}!="") {
	$uri_igv=$modules_obj_array["genomebrowser"]->{"submodules"}->{"igv"}->{"services"}->{"igv"}->{"href"};
} elseif ($modules_obj_array["genomebrowser"]->{"services"}->{"igv"}->{"href"}!="") {
	$uri_igv=$modules_obj_array["genomebrowser"]->{"services"}->{"igv"}->{"href"};
} else {
	$uri_igv="";
};

# VISION URL
if ($uri_vision == "") {
	if ($_ENV["URI_VISION"] != "") {
		$uri_vision=$_ENV["URI_VISION"];
	} elseif ($modules_obj_array["variantbrowser"]->{"submodules"}->{"jarvis"}->{"services"}->{"vision"}->{"href"}!="") {
		$uri_vision=$modules_obj_array["variantbrowser"]->{"submodules"}->{"jarvis"}->{"services"}->{"vision"}->{"href"};
	} elseif ($modules_obj_array["variantbrowser"]->{"services"}->{"vision"}->{"href"}!="") {
		$uri_vision=$modules_obj_array["variantbrowser"]->{"services"}->{"vision"}->{"href"};
	} else {
		$uri_vision="";
	};
};

# CLOUD URL
if ($uri_cloud == "") {
	if ($_ENV["URI_CLOUD"] != "") {
		$uri_cloud=$_ENV["URI_CLOUD"];
	} elseif ($modules_obj_array["cloud"]->{"submodules"}->{"cloud"}->{"services"}->{"cloud"}->{"href"}!="") {
		$uri_cloud=$modules_obj_array["cloud"]->{"submodules"}->{"cloud"}->{"services"}->{"cloud"}->{"href"};
	} elseif ($modules_obj_array["cloud"]->{"services"}->{"cloud"}->{"href"}!="") {
		$uri_cloud=$modules_obj_array["cloud"]->{"services"}->{"cloud"}->{"href"};
	} else {
		$uri_cloud="";
	};
};





?>
