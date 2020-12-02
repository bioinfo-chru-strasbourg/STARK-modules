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

#include "header.inc.php";


### VARIABLES
###############

# DEBUG
$SESSION=$_REQUEST["SESSION"];
$PHP_SELF=$_REQUEST["PHP_SELF"];
if ($SESSION=="unset") {
	#echo session_id();
	session_start();
	session_unset();
	session_destroy();
	#echo "unset";
	#echo $_SERVER["PHP_SELF"];
	header("Location: ".$PHP_SELF."");
} else {
	#echo "nothing";
};



?>