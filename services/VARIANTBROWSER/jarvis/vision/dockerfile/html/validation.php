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

#$dir_analysis="/media/IRC/V2/RES/ALL/";
#$dir_bin="/NGS/bin";

#
# Filters
# 
$input_type=$_REQUEST["type"];

$input_run=$_REQUEST["run"];
$input_id=$_REQUEST["id"];
$input_project_id=$_REQUEST["project_id"];
$input_sample_id=$_REQUEST["sample_id"];
$input_variant_id=$_REQUEST["variant_id"];
$input_flag_id=$_REQUEST["flag_id"];
#if ($input_flag_id=="") {$input_flag_id="0";};#if
$input_comment=$_REQUEST["comment"];
$input_author=$_REQUEST["author"];

print_r($_REQUEST);

$date=date("Y-m-d H:i:s");

# Find Flag
$query="SELECT VS.flag_id, VS.comment, VS.history, flags.label
	FROM comment AS VS
	INNER JOIN flags ON (flags.id=VS.flag_id)
	WHERE VS.id='$input_id'
	  AND VS.type='$input_type'
	  
	"; 
#AND VS.project_id='$input_project_id'
print "<pre>$query</pre><BR>";
$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
while ($row = mysql_fetch_assoc($result)) {
	$flag=$row['flag'];
	$comment=$row['comment'];
	$history=$row['history']; 
	#$label=$row['label']; 
};#while

# Find Label
$query="SELECT label AS input_label, score AS input_score
	FROM flags
	WHERE type='$input_type'
	  AND id='$input_flag_id'
	"; 
print "<pre>$query</pre><BR>";
$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
while ($row = mysql_fetch_assoc($result)) {
	$input_label=$row['input_label'];
	$input_score=$row['input_score'];
};#while


# Prepare history
#$comment_ini_array="comment[]=".join("\ncomment[]=",split("\n",str_replace("'","\\\'",$input_comment))); # ERROR
$comment_ini_array="comment='".str_replace("'","\\\'",$input_comment)."'";
#$history="[$date]\nflag=$input_flag\nlabel=$input_label\ncomment=$input_comment\nauthor=$input_author\n".$history;
$history="[$date]\nflag=$input_flag_id\nlabel=$input_label\nscore=$input_score\n$comment_ini_array\nauthor=$input_author\n".$history;
# INSERT/UPDATE
$query="REPLACE comment
		(`id`,`type`,`flag_id`,`comment`,`history`,`author`,`date`)
	VALUES	($input_id,'$input_type',$input_flag_id,\"$input_comment\",\"$history\",\"$input_author\",'$date')
	";
print "<pre>$query</pre><BR>";
$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));




?>
