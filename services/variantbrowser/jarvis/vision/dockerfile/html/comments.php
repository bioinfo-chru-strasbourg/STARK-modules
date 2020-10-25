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
# FUNCTIONS
#############


#
# INPUT
##########

#$dir_analysis="/media/NGS/analysis";
$dir_analysis="/media/IRC/V2/RES/ALL/";
$dir_bin="/NGS/bin";

#
# Filters
# 
$input_run=$_REQUEST["run"];
$input_sample=$_REQUEST["sample"];
$input_comment=$_REQUEST["comment"];
$input_date=$_REQUEST["date"];

#$input_comment=strip_tags($input_comment);

#$input_comment=preg_replace( "/\r|\n/", "", $input_comment );
#$input_comment=preg_replace( "/<br>/", "\n", $input_comment );
$input_comment=preg_replace( "/\"/", "'", $input_comment );
#$input_comment=nl2br(htmlspecialchars($input_comment);

print_r($_REQUEST);


if (1) { # Based on DB
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
	#echo $summary;
};#if

## UPDATE SAMPLE
$query_update="UPDATE sample SET comment=\"$input_comment\"
		WHERE sample.sample_id='$sample_id'
		"; 
	print "<pre>$query_update</pre><BR>";
	$result = mysql_query($query_update,$mydbh) or die (mysql_error($mydbh));



?>
