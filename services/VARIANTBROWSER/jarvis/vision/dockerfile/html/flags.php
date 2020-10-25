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
$date="20131023";
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

$text="";
$text.="<link rel='stylesheet' media='all' type='text/css' href='style.css' />";


#
# INPUT
##########

$dir_analysis="/media/IRC/V2/RES/ALL/";
$dir_bin="/NGS/bin";

#
# Filters
# 
$input_type=$_REQUEST["type"];
$input_run=$_REQUEST["run"];
$input_sample=$_REQUEST["sample"];
$input_variant=$_REQUEST["variant"];
$input_project_id=$_REQUEST["project_id"];
$input_sample_id=$_REQUEST["sample_id"];
$input_variant_id=$_REQUEST["variant_id"];
$input_flagvalue=$_REQUEST["flagvalue"];
$input_flag_id=$_REQUEST["flag_id"];


#$input_comment=strip_tags($input_comment);

#$input_comment=preg_replace( "/\r|\n/", "", $input_comment );
#$input_comment=preg_replace( "/<br>/", "\n", $input_comment );

print_r($_REQUEST);


if (1) { # Based on DB
	if ($input_type=="sample") {
		$query="SELECT sample.sample_id, sample.ref, sample.flag, sample.comment
			FROM sample
			INNER JOIN project ON (project.project_id=project.project_id)
			WHERE sample.run='$input_run'
			  AND sample.ref='$input_sample'
			  AND project.ref='$project'
			"; 
		print "$query<BR>";
		$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
		$summary= "<table class=' '>\n"; 
		while ($row = mysql_fetch_assoc($result)) {
			$sample_id=$row['sample_id'];
			$summary.= "<tr>"
				."<td>$project ".$row['sample_id']."</td>"
				."<td>".$row['ref']."</td>"
				."<td>".$row['flag']."</td>"
				."<td>".$row['comment']."</td>"
			."</tr>"; 
		};#while
		$summary.= "</table>\n"; 
		echo $summary;
	};#if
	if ($input_type=="variant_sample") {
		$query="SELECT VVS.flag_id, VVS.comment, VVS.history
			FROM validation_variant_sample AS VVS
			WHERE VVS.sample_id='$input_sample_id'
			  AND VVS.variant_id='$input_variant_id'
			"; 
		print "$query<BR>";
		$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
		$summary= "<table class=' '>\n"; 
		while ($row = mysql_fetch_assoc($result)) {
			$flag_id=$row['flag_id'];
			$flag_comment=$row['comment'];
			$flag_history=$row['history'];
			$summary.= "<tr>"
				."<td>$project ".$row['sample_id']."</td>"
				."<td>".$row['variant_id']."</td>"
				."<td>".$row['flag_id']."</td>"
			."</tr>"; 
		};#while
		$summary.= "</table>\n"; 
		echo $summary;
	};#if
};#if

## UPDATE SAMPLE
if ($input_type=="sample") {
$query_update="UPDATE sample SET flag=$input_flagvalue
		WHERE sample.sample_id='$sample_id'
		"; 
	print "<pre>$query_update</pre><BR>";
	$result = mysql_query($query_update,$mydbh) or die (mysql_error($mydbh));
};#if

## UPDATE VARIANT_SAMPLE
$date=date(DATE_ATOM);
if ($input_type=="variant_sample" and $input_flag_id!="") {
	$history="[$date]\nflag=$input_flag_id\ncomment=\n".$flag_history;
	if ($flag_id!="") {
		$query_update="UPDATE validation_variant_sample
		SET	flag_id=$input_flag_id,
			history='$history'
		WHERE validation_variant_sample.variant_id='$input_variant_id'
		  AND validation_variant_sample.sample_id='$input_sample_id'
		  AND validation_variant_sample.project_id='$input_project_id'
		"; 
		print "<pre>$query_update</pre><BR>";
		$result = mysql_query($query_update,$mydbh) or die (mysql_error($mydbh));
	} else {
		$query_insert="INSERT INTO validation_variant_sample
			VALUES ($input_variant_id,$input_sample_id,$input_project_id,
				$input_flag_id,'','$history','')
		"; 
		print "<pre>$query_insert</pre><BR>";
		$result = mysql_query($query_insert,$mydbh) or die (mysql_error($mydbh));
	};#if
};#if


?>
