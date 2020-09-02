<?php

### CONFIG
############


### APP
$APP_CODE="EDITH";
$APP_NAME="Excellent Dashboard and Interfaces To Help STARK experience";
$APP_RELEASE="1.0";
$APP_DATE="20200523";
$APP_DESCRIPTION="EDITH provides an IHM as a Dashboard for STARK analyses and results";
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


### FOLDERS

# SERVICES folder
if ($_ENV["FOLDER_SERVICES"] != "") {
	$folder_services=$_ENV["FOLDER_SERVICES"];
} else {
	$folder_services="services";
};


### DOCKER
# ensure /var/run/docker.sock is chmod o+rw

$docker_containers_obj="";
#$command="curl --unix-socket /var/run/docker.sock http://localhost/containers/json?all=1";
$command="curl --unix-socket /var/run/docker.sock http://localhost/containers/json?all=1";
#$command="curl -XGET --unix-socket /var/run/docker.sock -d '{\"all\":\"1\"}' -H 'Content-Type: application/json' http://localhost/containers/json";
#$command="curl --unix-socket /var/run/docker.sock http://sockguard_sockguard_1/containers/json?all=1";
#$command="curl --unix-socket /app/docker.sock http://localhost/containers/json?all=1";

$output_exec = shell_exec($command);
$docker_containers_obj=json_decode($output_exec);

if (is_array($docker_containers_obj)) {
	$docker_containers_connected=1;
} else {
	$docker_containers_connected=0;
};

// echo "<pre>";
// print_r($docker_containers_obj);
// echo "</pre>";


foreach ($docker_containers_obj as $docker_containers_obj_id=>$docker_containers_obj_infos) {

	# Container name
	$container_name=$docker_containers_obj_infos->{"Names"}[0];
	$container_name=ltrim($container_name, $container_name[0]);
	
	# docker_containers_obj_array
	$docker_containers_obj_array[$container_name]["obj_id"]=$docker_containers_obj_id;
	$docker_containers_obj_array[$container_name]["obj_infos"]=$docker_containers_obj_infos;

};



### MODULES

# Find modules configurations
$modules = array_filter(glob($folder_services.'/*/STARK.module'), 'is_file');

# init
$modules_obj_array;

# Foeach modules
foreach ($modules as $module_file) {

	# Read JSON file
	$module_json=implode("",file($module_file));
	$module_obj=json_decode($module_json);
	
	# Module code
	$module_info_code=$module_obj->{'code'};
	
	# Services array
	foreach ($module_obj->{'services'} as $service=>$service_infos) {
	
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

		# Available
		$service_available=$service_infos->{'available'}!=""?$service_infos->{'available'}:false;
		$service_infos->{'available'}=$service_available;
	
		# Service container
		$service_container=$service_infos->{'container'}!=""?$service_infos->{'container'}:"STARK-module-$module_info_code-service-$service_code";
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
		$module_obj->{'services'}->{$service}=$service_infos;

	};

	# Module array
	$modules_obj_array[$module_info_code]=$module_obj;
};


### PLUGINS

$plugins_folder=".";
$plugins_file="config/plugins.conf";
$plugins_json=implode("",file($plugins_file));
$plugins_obj=json_decode($plugins_json);



### FOLDERS

# INPUTS folder
if ($_ENV["FOLDER_INPUTS"] != "") {
	$folder_inputs=$_ENV["FOLDER_INPUTS"];
} else {
	$folder_inputs="inputs";
};

# REPOSITORIES folder
if ($_ENV["FOLDER_REPOSITORIES"] != "") {
	$folder_repositories=$_ENV["FOLDER_REPOSITORIES"];
} else {
	$folder_repositories="repositories";
};

# API
if ($_ENV["FOLDER_API"] != "") {
	$folder_api=$_ENV["FOLDER_API"];
} else {
	$folder_api=$folder_services."/STARK/API";
};

# listener
if ($_ENV["FOLDER_LISTENER"] != "") {
	$folder_listener=$_ENV["FOLDER_LISTENER"];
} else {
	$folder_listener=$folder_services."/STARK/listener";
};



### URI

# API
if ($_ENV["URINNER_API"] != "") {
	$urinner_api=$_ENV["URINNER_API"];
} elseif ($modules_obj_array["STARK"]->{"services"}->{"API"}->{"href_inner"}!="") {
	$urinner_api=$modules_obj_array["STARK"]->{"services"}->{"API"}->{"href_inner"};
} else {
	$urinner_api="";
};

$urinner_api_list=$urinner_api."queue?list";
$urinner_api_info=$urinner_api."queue?info=";
$urinner_api_log=$urinner_api."queue?log=";


# DAS
if ($_ENV["URINNER_DAS"] != "") {
	$urinner_das=$_ENV["URINNER_DAS"];
} elseif ($modules_obj_array["STARK"]->{"services"}->{"DAS"}->{"href_inner"}!="") {
	$urinner_das=$modules_obj_array["STARK"]->{"services"}->{"DAS"}->{"href_inner"};
} else {
	$urinner_das="";
};

# DAS
if ($_ENV["URI_DAS"] != "") {
	$uri_das=$_ENV["URI_DAS"];
} elseif ($modules_obj_array["STARK"]->{"services"}->{"DAS"}->{"href"}!="") {
	$uri_das=$modules_obj_array["STARK"]->{"services"}->{"DAS"}->{"href"};
} else {
	$uri_das="";
};



# IGV
if ($_ENV["URINNER_IGV"] != "") {
	$urinner_igv=$_ENV["URINNER_IGV"];
} elseif ($modules_obj_array["GENOMEBROWSER"]->{"services"}->{"IGV"}->{"href_inner"}!="") {
	$urinner_igv=$modules_obj_array["GENOMEBROWSER"]->{"services"}->{"IGV"}->{"href_inner"};
} else {
	$urinner_igv="";
};

# IGV
if ($_ENV["URI_IGV"] != "") {
	$uri_igv=$_ENV["URI_IGV"];
} elseif ($modules_obj_array["GENOMEBROWSER"]->{"services"}->{"IGV"}->{"href"}!="") {
	$uri_igv=$modules_obj_array["GENOMEBROWSER"]->{"services"}->{"IGV"}->{"href"};
} else {
	$uri_igv="";
};


$RESULTS_SUBFOLDER_DATA="STARK";


# DEBUG
// echo "<pre>";
// print_r($plugins_obj);
// echo "</pre>";


?>
