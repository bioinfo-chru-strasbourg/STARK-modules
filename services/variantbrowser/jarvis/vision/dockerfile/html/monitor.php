<?php
/*
 *
 * Released under the terms and conditions of the
 * GNU General Public License (http://www.gnu.org/licenses/gpl.txt)
 *
*/

#
# Versioning
##############

$release="0.9";
$date="23/03/2015";
$authors="Antony Le Béchec";

#
# HEADER
##########

if(!defined('e107_INIT'))
{
	require_once('../../class2.php');
}
$e107 = e107::getInstance();
$tp = e107::getParser();
$sql = e107::getDb();

require_once(HEADERF);
require("functions.inc.php");
require "connect.php"; 
require "config.php"; 

#
# PARAM
#########


#
# FUNCTIONS
#############



##########
## MAIN ##
##########

#
# INPUT
##########

#
# Main input
# 

$action_default="?";
$s=$_REQUEST["s"];
$q=$_REQUEST["q"];
$action_todo=$_REQUEST["action_todo"];
$input_run=$_REQUEST["run"];
$input_sample=$_REQUEST["sample"];

#
# VCF File
# 

$input_runs=$_REQUEST["runs"];
#if ($input_runs=="") {
#	$input_runs[]="";
#};#if


#
# Aligners
# 

$aligners=$_REQUEST["aligners"];
if ($aligners=="") {
	$aligners[]="";
};#if

#
# Callers
# 

$callers=$_REQUEST["callers"];
if ($callers=="") {
	$callers[]="";
};#if

#
# Annotators
# 

$annotators=$_REQUEST["annotators"];
if ($annotators=="") {
	$annotators[]="";
};#if

#
# Filters
# 

$filters=$_REQUEST["filters"];
if ($filters=="") {
	$filters[]="default";
};#if
$hardfiltering=$_REQUEST["hardfiltering"];

#
# Order by
# 

$orderby=$_REQUEST["orderby"];
if ($orderby=="") {
	$orderby="FilterScore";
};#if
$ascdesc=$_REQUEST["ascdesc"];
if ($ascdesc=="") {
	$ascdesc="DESC";
};#if

#
# Limit
# 

$limit=$_REQUEST["limit"];
if ($limit=="") {
	$limit="20";
};#if


#
# Global Filter
# 

$global_filter=$_REQUEST["global_filter"];
if ($global_filter=="") {
	$global_filter="";
};#if

#
# Gene Filter
# 

$gene_filter=$_REQUEST["gene_filter"];
if ($gene_filter=="") {
	$gene_filter="";
};#if


#
# GROUP / project_id
#
if ($project_id=="") {
	$project_id=1;
};#if



# Security
#if (!ADMIN) {
if (1) {
#	$text="<div style='text-align:center'>Please login</div>";
#} else {

	##
	## Style
	##

	
	$text="";
	$text.="<link rel='stylesheet' media='all' type='text/css' href='style.css' />";

	$text.="
	<SCRIPT>
	function show_hide(id) {
		if(document.getElementById(id).style.display=='none') {
			document.getElementById(id).style.display='inline';
		} else {
			document.getElementById(id).style.display='none';
		}
		return true;
	}
	</SCRIPT>
	";

	# Header
	##########

	include "header.php";

	# MENU
	########

	$command="ps -A | grep ' launch.sh$' | awk -F' ' '{print $1}'";
	$command_ID = shell_exec($command);
	$launched=0;
	if ($command_ID!="") {
		$launched=$command_ID;
	};#if;



	# Action
	##########

	$output_analysis="";
	if ($launched>0) {
		$output_analysis.="<p class='warning'>RUN Analysis process already launched: $launched</p>Please wait the end of the process to launch another analysis<BR><BR>";
	} else {
		$output_analysis.="RUN Analysis process unused<BR><BR>";
	};#if


	

	$command="$dir_bin/scripts/monitor.sh 1 details";
	$output_exec = shell_exec($command);
	$output_analysis.="<PRE>$output_exec</PRE>";

	#$command="top -b -n 1 | head -n 4";
	#$output_exec = shell_exec($command);
	#$output_analysis.="<PRE>$output_exec</PRE>";


	$text.=$output_analysis;

};#if

	$text.= "
			</div>
		</div>
	";


$title=" ";

$ns->tablerender($title, $text, "TRAKXS");



require_once(FOOTERF);

?>
