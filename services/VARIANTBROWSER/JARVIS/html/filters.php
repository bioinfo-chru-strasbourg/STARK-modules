<?php
/*
 * e107 website system
 *
 * Copyright (C) 2008-2013 e107 Inc (e107.org)
 * Released under the terms and conditions of the
 * GNU General Public License (http://www.gnu.org/licenses/gpl.txt)
 *
 * Forum main page
 *
*/
#if(!defined('e107_INIT'))
#{
#	require_once('../../class2.php');
#}
#$e107 = e107::getInstance();
#$tp = e107::getParser();
#$sql = e107::getDb();

#require_once(HEADERF);
require("functions.inc.php");
#require "connect.php";
require "config.php";

#$text="";
#$text.="<link rel='stylesheet' media='all' type='text/css' href='style.css' />";


$text="";
$text.="<HEAD>";
$text.="<title>JARVIS/VISION</title>";
$text.="<link rel='stylesheet' media='all' type='text/css' href='style.css' />";
$text.="</HEAD>";

$text.="<BODY>";

# Header
##########

include "header.php";
#include "search.header.php";

$text.="

<TABLE class=' ' width='100%' border=0>
<TR width='100%'>
	<!--<TD class='suptitle ' width='50px'>
		$VISION_name
	</TD>
	<TD width='1px'></TD>-->
	<TD width='*' class=' subtitle'>
		<span class='' title=''><B><BIG><BIG>$VISION_name |</BIG></BIG> $VISION_description</B></span> [$VISION_release]
		<span style=''><BR>$VISION_comment<BR></span>
	</TD>
	<!--<TD width='1px'></TD>
	<TD width='80px' class=' subtitle' style='text-align:center;'>
		<span style=''><B>".$username."</B><BR>".date("d-m-Y")."<BR>".date("H:i:s")."</span>
	</TD>-->
</TR>
<TR width='100%'>
	<TD class=' ' heigth='1px'>

	</TD>
</TR>
</TABLE>
";


#
# INPUT
##########

$dir_analysis="/media/NGS/analysis";
$dir_bin="/NGS/bin";

#
# Filters
#
$filters=$_REQUEST["filters"];
if ($filters=="") {
	$filters[]="default";
};#if

#
# Pre computing
#################

#
# Filters
#



#print $filter_test;
$filter_table="";
#$filter_file="$dir_bin/scripts/filters.ini";
#$filter_file="$dir_howard/config.filter.ini";
#$filter_file="$dir_howard/config.filter.ini";
$filter_file=$config_filter_ini;
$filter_array=parse_ini_file($filter_file,TRUE);
#echo "Filter_array before:<BR>"; print_r($filter_array); echo "<BR>"; echo "<BR>";
$filter_array=basedon_filters($filter_array);
#echo "Filter_array after:<BR>"; print_r($filter_array); echo "<BR>"; echo "<BR>";
$filter_select="";
foreach ($filter_array as $filter_name=>$filter_criteria) {
	$filter_select.="<OPTION value='$filter_name' ".(in_array("$filter_name",$filters)?"selected":"").">$filter_name (".count($filter_criteria)." criteria)\n";
};#foreach
$filter_option_array=array();
$filter_option="";
$filter_option_array=load_filter($filters,$filter_array,$filter_option_array);
#print_r($filter_option_array);
if (!empty($filter_option_array)) {
	foreach ($filter_option_array as $key=>$filter_option2) {
		$filter_table.="<TR><TD>".join("</TD><TD>",split(":",$filter_option2))."</TD></TR>";
	};#foreach
};#if
$filter_table="
	<TABLE style='font-size:small;' border=1>
		<TR><TD>Annotation</TD><TD>Criterion</TD><TD>Ponderation</TD><TD>Comment</TD></TR>
		$filter_table
		</TABLE>
";

#
# Output
##########


#include "header.php";

$legend="
<BR>
	<UL>
		<LI>Annotation: annotation label (INFO/* or FORMAT/*) to filter</LI>

		<LI>Criterion: evaluation of the annotation value</LI>
		<LI>Ponderation : either a score (positive or negative number) or a flag (f for FILTER, p for force PASS), added in variant PZScore and PZFlag</LI>
		<LI>Comment: text to add in the variant comment (PZComment). 'VALUE' is the value of the annotation</LI>
	</UL>
<BR>
";




$text.="
	<div style='clear:left;'></DIV>
	<div style='text-align: left;' class=''>
		$legend
		<BR>
		Filter based on : ".join(", ",$filters)."<BR><BR>
		$filter_table
	</div>
	<BR>

	";





$title="";

print $text;

#$ns->tablerender($title, $text, "Filters");





#require_once(FOOTERF);

?>
