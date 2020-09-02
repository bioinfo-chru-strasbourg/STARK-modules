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
$date="20141231";
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

#require_once(HEADERF);
require("functions.inc.php");
require "connect.php"; 
require "config.php"; 

$text="";
$text.="<link rel='stylesheet' media='all' type='text/css' href='style.css' /><BODY width='100%'>"; # style='background-color:white;'


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
$input_flag=$_REQUEST["flag"];
if ($input_flag=="") {$input_flag="0";};#if
$input_comment=$_REQUEST["comment"];
$input_author=$_REQUEST["author"];

#print_r($_REQUEST);
#print "TEST";
$date=date("Y-m-d H:i:s");

# Find Flag
$query="SELECT VS.flag_id AS flag, VS.comment, VS.history
	FROM comment AS VS
	WHERE VS.id='$input_id'
	  AND VS.type='$input_type'
	 # AND VS.project_id='$input_project_id'
	";
#print "<pre>$query</pre><BR>";
$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
while ($row = mysql_fetch_assoc($result)) {
	$flag_id=$row['flag'];
	$comment=$row['comment'];
	$history=$row['history']; 
};#while
# Prepare history
#$text.= "<BR>$history</BODY>"; # <TEXTAREA id='comment' name='comment' class='' style='border: none; width: 99%; height: 100%; -webkit-box-sizing: border-box; -moz-box-sizing: border-box; box-sizing: border-box; padding:0px;' rows=2>$history</TEXTAREA>
#$text.= "<TEXTAREA class='' disabled style='font-size:small; border: none; width: 100%; height: 100%; -webkit-box-sizing: border-box; -moz-box-sizing: border-box; box-sizing: border-box; padding:0px;' rows=2>$history</TEXTAREA></BODY>"; # 

$history_array=parse_ini_string($history,true);
#print_r($history_array);
$history_table="<TABLE  valign='top'>";
foreach ($history_array as $date=>$infos) {
	$history_table.="<TR><TD colspan='2' style='font-weight:;'>[$date]</TD></TR>";
	foreach ($infos as $variable=>$value) {
		if (is_array($value)) {
			$history_table.="<TR><TD valign='top'>$variable</TD><TD valign='top'>".join("<BR>",$value)."</TD></TR>";
		} else {
			$history_table.="<TR><TD valign='top'>$variable</TD><TD valign='top'>$value<TD></TR>";
		};#if
	};#foreach
};#foreach
$history_table.="</TABLE>";
$text.=$history_table;

#$history_html=str_replace("\n","<BR>",$history);

#$text.= "<small>$history_html</small>";
/*
$history_test="[section1]
value1=1
value2[]=2
value2[]=2b
value3=3
[section2]
value4=4";
$history_array=parse_ini_string($history,true);
print_r($history_array);
*/
#$history="[$date]\nflag=$input_flag\ncomment=$input_comment\nauthor=$input_author\n".$history;
# INSERT/UPDATE

echo $text;


?>
