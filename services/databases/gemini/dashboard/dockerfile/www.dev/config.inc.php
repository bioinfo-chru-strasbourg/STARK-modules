<?php

### CONFIG
############


### APP
$APP_CODE="Gemini";
$APP_NAME="Gemini IHM";
$APP_RELEASE="1.0";
$APP_DATE="20200622";
$APP_DESCRIPTION="Dashboard for databases tracking and exploration";
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



?>
