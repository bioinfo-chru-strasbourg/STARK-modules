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

$release="0.9.1";
$date="20131025";
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

#$dir_analysis="/media/NGS/analysis";
#$dir_miseq="/media/miseq/MSR";
$dir_analysis="/media/IRC/V2/RES/ALL";
$dir_miseq="/media/IRC/RAW/MSR";
$dir_bin="/NGS/bin";

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
$file_to_show=$_REQUEST["file_to_show"];

#
# VCF File
# 

$vcf_files=$_REQUEST["vcf_files"];
if ($vcf_files=="") {
	$vcf_files[]="";
};#if


#
# Aligners
# 

$aligners=$_REQUEST["aligners"];
if ($aligners=="") {
	$aligners[]="MiSeq-Aligner";
};#if

#
# Callers
# 

$callers=$_REQUEST["callers"];
if ($callers=="") {
	$callers[]="MiSeq-Caller";
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
# Stats
# 

$year=$_REQUEST["year"];
$month=$_REQUEST["month"];
$day=$_REQUEST["day"];


#
# Annotation List
#

$annotation_list=$_REQUEST["annotation_list"];
# Mandatory list
$annotation_list_mandatory=array("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FilterScore", "FilterFlag", "FilterComment", "geneSymbol", "location", "outcome", "GQ", "BQ", "DP", "AD", "VF", "AF", "FA", "dbSNP", "dbSNPNonFlagged");
#$annotation_list_mandatory=array("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FilterScore", "FilterFlag", "geneSymbol", "location", "outcome", "GQ", "BQ", "DP", "AD", "VF", "AF", "FA", "dbSNP", "dbSNPNonFlagged");
# Default list (combinason of a list of annotation and the madatory lists
$annotation_list_default="";
$annotation_list_default.=join(",",$annotation_list_mandatory);
if ($orderby != "" && !in_array($orderby,$annotation_list) && !in_array($orderby,$annotation_list_mandatory)) {
	$annotation_list[]=$orderby;
};#if	
# NO annotation if list is empty (or NO)
if (empty($annotation_list) || (in_array("NO",$annotation_list))) {# && count($annotation_list)<=1)) {
	$annotation_list[0]="NO";
	$annotation_list=array_merge($annotation_list,split(",",$annotation_list_default));
};#if

#$annotation_list=array_merge($annotation_list,$annotation_list_mandatory);
$annotation_list=array_merge($annotation_list_mandatory,$annotation_list);


##
## Style
##

#<div class='view-item'>
#width:32px;
#	filter:alpha(opacity=50);
#	opacity:0.5;
#	background-color:#DD55DD;
#	background-color:rgba(0, 0, 255, 0.5);
#	background-color:#55DD55;
#	background-color:rgba(128, 128, 255, 0.5);
#        -webkit-box-shadow: 0 1px 2px rgba(0,0,0,.05);
#           -moz-box-shadow: 0 1px 2px rgba(0,0,0,.05);
#                box-shadow: 0 1px 2px rgba(0,0,0,.05);
#	float:left;
##e5e5e5;
#	background-color:rgba(255, 165, 0, 0.5);
#	filter:progid:DXImageTransform.Microsoft.gradient(startColorstr=#7FFFA500,endColorstr=#7FFFA500);
#filter:progid:DXImageTransform.Microsoft.gradient(startColorstr=#7F9932cc,endColorstr=#7F9932cc); b22222
	
$text="";
$text.="<link rel='stylesheet' media='all' type='text/css' href='style.css' />";


# Header
##########

include "header.php";



# Security
if (!USER) {
	#$text.="<div style='text-align:center;background-color:'>Please<BR><A href='../../login.php'>Sign in</A></div>";
	include "welcome.php";
	#http://192.168.59.42/JARVISdev/login.php

} else {

include "search.header.php";

# Summary
###########

require "connect.php"; 


   $nblastruns=10;

if ($user_projects!="") {
	#$query="SELECT project.ref AS project, sample.run AS run, count(sample.sample_id) AS nb_sample FROM project, sample WHERE project.project_id=sample.project_id
	#	GROUP BY project.project_id, sample.run ORDER BY sample_id DESC LIMIT $nblastruns"; 
	$query="SELECT project.ref AS project, run.ref AS run, count(sample.id) AS nb_sample
		FROM sample
		INNER JOIN run ON (run.id=sample.run_id)
		INNER JOIN project ON (project.id=sample.project_id)
		WHERE project.id IN (".join(",",array_keys($user_projects)).")
		GROUP BY run.id
		ORDER BY run.ref DESC
		LIMIT $nblastruns
		 "; #$user_projects
	#echo "<pre>$query</pre>";
	$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
	#require "check.php"; 

	$lastruns= "<table class=' '>\n"; 
	#$summary.= "<tr><td>RUN</td><td>PROJECT</td><td>#SAMPLE</td></tr>\n"; 
	while ($row = mysql_fetch_assoc($result)) {
		$RUN=$row['run'];
		$RUNLink="run.php?action_todo=show_run&run=$RUN&q=$q";
		$RUNList.="<DIV style='clear: left;'>".ListSAMPLE($Directory.'/'.$RUN,$RUN,array("NOTHING"))."</DIV>";
		$lastruns.= "<tr>"
			."<td>".ListSAMPLE($Directory.'/'.$RUN,$RUN,array("NOTHING"))."</td>"
			#."<td>".$row['project']."</td>"
			#."<td>".$RUNList."</td>"
			#."<td>".$row['nb_sample']."</td>"
		  ."</tr>"; 
	}
   $lastruns.= "</table>\n"; 
   #echo $summary;

	if ($_REQUEST["all"] || 0) {

		$query="SELECT run.ref AS run, manifest.ref AS manifest, project.ref AS project, sample.ref AS sample
			FROM sample
			INNER JOIN run ON (run.id=sample.run_id)
			INNER JOIN project ON (project.id=sample.project_id)
			INNER JOIN manifest ON (manifest.id=sample.manifest_id)
			WHERE project.id IN (".join(",",array_keys($user_projects)).")
			ORDER BY manifest.ref ASC, run.ref ASC, sample.ref ASC";
		$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
		
		$allsamples_array=array();
		$fields=array("manifest","run");
		
		while ($row = mysql_fetch_assoc($result)) {
			$run=$row['run'];
			$manifest=$row['manifest'];
			$project=$row['project'];
			$sample=$row['sample'];
			$RUNLink="run.php?action_todo=show_run&run=$run";
			$SAMPLELink="sample.php?action_todo=show_sample&run=$run&sample=$sample";
			$RUNList.="<DIV style='clear: left;'>".ListSAMPLE($Directory.'/'.$RUN,$RUN)."</DIV>";
			
			$allsamples_array[$row[$fields[0]]][$row[$fields[1]]][]=$row['sample'];

			if ($manifest!=$manifest_prev) {
				$allsamples.= "<TR><TD colspan=10>$manifest</TD></TR>";
			};#if
			if ($manifest!=$manifest_prev||$run!=$run_prev) {
				$allsamples.= "<TR><TD width='20px'></TD><TD colspan=10><A href='$RUNLink'\">$run</A></TD></TR>";
			};#if
			
			$allsamples.= "<tr>"
				#."<td>".(($manifest!=$manifest_prev)?$manifest:"")."</td>"
				#."<td width='20px'></td>"
				#."<td>".(($manifest!=$manifest_prev||$run!=$run_prev)?"<A href='$RUNLink'\">$run</A>":"")."</td>"
				."<td></td>"
				."<TD width='20px'>"
				."<td><A href='$SAMPLELink'>$sample</A> ($project)</td>"
				#."<td>".$row['project']."</td>"
				#."<td>".$RUNList."</td>"
				#."<td>".$row['nb_sample']."</td>"
			  ."</tr>"; 

			$run_prev=$run;
			$manifest_prev=$manifest;
			$project_prev=$project;
			$sample_prev=$sample;


		};#while
		$allsamples="<table border=0 class='' width='400'>";
		foreach ($allsamples_array as $field1=>$infos1) {
			$allsamples.= "<TR><TD colspan=3>$field1</TD></TR>";
			foreach ($infos1 as $field2=>$infos) {
				#$RUNLink="run.php?action_todo=show_run&run=$field2";
				$allsamples.= "<TR><TD width='50px'></TD><TD colspan=2 class='brick run'>$field2</TD></TR>";
				$allsamples.= "<TR><TD width='50px'></TD><TD width='50px'></TD><TD colspan=1>";
				foreach ($infos as $key=>$val) {
					$SAMPLELink="sample.php?action_todo=show_sample&run=$field2&sample=$val";
					$allsamples.= "<A href='$SAMPLELink' class='brick sample'>$val</A> ";
				};#foreach
				$allsamples.= "</TD></TR>";
			};#foreach
		};#foreach
		
		$allsamples.="</table>";

		$allsamples="<table border=0 class='' width='90%'>";
		foreach ($allsamples_array as $field1=>$infos1) {
			$allsamples.= "<TR><TD colspan=10 height='32px'></TD></TR>";
			$allsamples.= "<TR><TD colspan=10 height=''><DIV class=' title' style='width:100%;height:24px;'>$field1</DIV></TD></TR>";
			foreach ($infos1 as $field2=>$infos) {
				#$RUNLink="run.php?action_todo=show_run&run=$field2";
				$allsamples.= "<TR><TD colspan=0 valign='top' width=''><DIV class='brick run' style='width:300px;text-align:center;'>$field2</DIV></TD><TD align='top'>";
				foreach ($infos as $key=>$val) {
					$SAMPLELink="sample.php?action_todo=show_sample&run=$field2&sample=$val";
					$allsamples.= "<span onclick=\"javascript:document.location='$SAMPLELink';\" class='brick sample'>$val</span> ";
				};#foreach
				$allsamples.= "</TD></TR>";
			};#foreach
		};#foreach
		
		$allsamples.="</table>";		


	};#if

};#if


 	#$query="SELECT project.ref AS project, count(sample.run_id) AS nb_run FROM project, sample WHERE project.project_id=sample.project_id
	#	GROUP BY project.project_id"; 
	#$year="15";
	#$month="";
	#$day="";
	$query="SELECT `group`.ref AS 'group', count(distinct project.id) AS nb_project, count(distinct manifest.id) AS nb_manifest, count(distinct run.id) AS nb_run, count(distinct sample.id) AS nb_sample
	FROM sample
	INNER JOIN run ON (run.id=sample.run_id)
	INNER JOIN manifest ON (manifest.id=sample.manifest_id)
	INNER JOIN project ON (project.id=sample.project_id)
	INNER JOIN `group` ON (`group`.id=project.group_id)
	WHERE `group`.ref NOT IN ('MISC')
	  AND `run`.ref NOT LIKE '%.%'
	  AND `run`.ref LIKE '$year$month$day%'
	GROUP BY `group`.id
	ORDER BY `group`.ref, project.ref";
	#print "$query";
   $result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));

   #require "check.php"; 

   $basicstats= "<table class=' ' cellpadding='10' >\n"; 
   #$basicstats.= "<tr  class=' ' style='background:gray;'   ><th>Group</th><th>Project</th><th>#Manifest</th><th width='20'><BR></th><th>#RUN</th><th>#SAMPLE</th></tr>\n"; 
   $basicstats.= "<tr  class=' ' style='background:gray;'><th>Group</th><th>#Project</th><th>#Manifest</th><th>#RUN</th><th>#SAMPLE</th></tr>\n";
 
   while ($row = mysql_fetch_assoc($result)) {
	$RUN=$row['run'];
	$RUNLink="run.php?action_todo=show_run&run=$RUN&q=$q";
	$RUNList.="<DIV style='clear: left;'>".ListSAMPLE($Directory.'/'.$RUN,$RUN)."</DIV>";
      $basicstats.= "<tr class=''  style='background:silver;'>"
		."<td>".$row['group']."</td>"
		."<td>".$row['nb_project']."</td>"
		."<td>".$row['nb_manifest']."</td>"
		#."<td>".$RUNList."</td>"
		."<td>".$row['nb_run']."</td>"
         	."<td>".$row['nb_sample']."</td>"
          ."</tr>"; 
	$data_run.=$row['nb_run'].",";
	$data_sample.=$row['nb_sample'].",";
	$data_manifest.=$row['nb_manifest'].",";
	$data_project.=$row['nb_project'].",";
	$group_list.=$row['group'].",";
   }

	#$url_pie_run="http://localhost/NGSV2/e107_plugins/NGSV2_plugin/jpgraph_stats.php?title=&subtitle=Run|by%20Group&data=$data_run&label=$group_list";
	#$url_pie_sample="http://localhost/NGSV2/e107_plugins/NGSV2_plugin/jpgraph_stats.php?title=&subtitle=Sample|by%20Group&data=$data_sample&label=$group_list";
	$url_pie_run="jpgraph_stats.php?title=&subtitle=Run|by%20Group&data=$data_run&label=$group_list";
	$url_pie_sample="jpgraph_stats.php?title=&subtitle=Sample|by%20Group&data=$data_sample&label=$group_list";

   $basicstats.= "</table>\n"; 


	# Stats for financial questions

	$query="SELECT run.ref AS run, `group`.ref AS 'group', count(sample.id) AS nb_sample
		FROM sample
		INNER JOIN run ON (run.id=sample.run_id)
		INNER JOIN project ON (project.id=sample.project_id)
		INNER JOIN `group` ON (`group`.id=project.group_id)
		WHERE `group`.ref NOT IN ('MISC')
		GROUP BY `run`.id , `group`.ref
		ORDER BY `run`.ref";
		#print "$query";
   $result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));

   #require "check.php"; 

   $basicstats2= "<table class=' ' cellpadding='10' >\n"; 
   #$basicstats.= "<tr  class=' ' style='background:gray;'   ><th>Group</th><th>Project</th><th>#Manifest</th><th width='20'><BR></th><th>#RUN</th><th>#SAMPLE</th></tr>\n"; 
   $basicstats2.= "<tr  class=' ' style='background:gray;'><th>RUN</th><th>GROUP</th><th>#SAMPLE</th></tr>\n";
 
   while ($row = mysql_fetch_assoc($result)) {
	#print_r($row); echo "<BR>";
	$RUN=$row['run'];
	$RUNLink="run.php?action_todo=show_run&run=$RUN&q=$q";
	$RUNList.="<DIV style='clear: left;'>".ListSAMPLE($Directory.'/'.$RUN,$RUN)."</DIV>";
      $basicstats2.= "<tr class=''  style='background:silver;'>"
		."<td>".$row['run']."</td>"
		."<td>".$row['group']."</td>"
		."<td>".$row['nb_sample']."</td>"
          ."</tr>"; 
	$run_distribution[$row['run']][$row['group']]=$row['nb_sample'];
	$run_nb_sample[$row['run']]+=$row['nb_sample'];
	
   }

	#$url_pie_run="http://localhost/NGSV2/e107_plugins/NGSV2_plugin/jpgraph_stats.php?title=&subtitle=Run%20by%20Group&data=$data_run&label=$group_list";
	#$url_pie_sample="http://localhost/NGSV2/e107_plugins/NGSV2_plugin/jpgraph_stats.php?title=&subtitle=Sample%20by%20Group&data=$data_sample&label=$group_list";

   $basicstats2.= "</table>\n"; 

	foreach ($run_distribution as $run=>$groups) {
		foreach ($groups as $group=>$nb_sample) {
			$nb_sample_in_run=$run_nb_sample[$run];
			$run_distribution_pond[$group]+=($nb_sample/$nb_sample_in_run);
		};#foreach
	};#foreach

	foreach ($run_distribution_pond as $group=>$pond) {
		$data_group.="$pond,";
		$group_list2.="$group,";
	};#foreach

	

	#$url_pie_activity="http://localhost/NGSV2/e107_plugins/NGSV2_plugin/jpgraph_stats.php?title=&subtitle=Activity|by group&data=$data_group&label=$group_list2";
	$url_pie_activity="jpgraph_stats.php?title=&subtitle=Activity|by group&data=$data_group&label=$group_list2";
	
	$text.= "<TABLE width='' border=0>
			<TR>
				<TD valign='top' width=''><p class=' title'>Last RUNs</p><BR>$lastruns</TD>
				<TD width='50'><BR></TD>
				<TD valign='top' width=''>
					<p class=' title'>Stats</p>
					| 
					<A HREF='?'>ALL</A> | 
					<A HREF='?year=13'>2013</A> | 
					<A HREF='?year=14'>2014</A> | 
					<A HREF='?year=15'>2015</A> | 
					<A HREF='?year=16'>2016</A> | 

					<BR>
					<CENTER>
						<A href='$url_pie_activity' target='pie'><IMG src='$url_pie_activity'></A></IMG>
						<BR>
						<A href='$url_pie_run' target='pie'><IMG src='$url_pie_run' width='200' height='200'></IMG></A>
						<A href='$url_pie_sample' target='pie'><IMG src='$url_pie_sample' width='200' height='200'></IMG></A>
					</CENTER>
					<BR><BR>
					<CENTER>$basicstats</CENTER>
					<BR>
				</TD>
			</TR>
		</TABLE>
		<BR>
		<!--<TABLE width='90%' border=0>
			<TR>
				<TD colspan='4' width='100%'><p class=' title'>ALL RUN/SAMPLE by MANIFEST</p><BR>
					$allsamples
				</TD>
			</TR>
		</TABLE>-->";


	
   mysql_free_result($result);


   mysql_close($mydbh); 
	

#DB Connexion TEST
#echo "Connexion to DB 'TRAKXS'<BR>";
#$mydb = new db();
#$mydb->db_Connect("localhost", "root", "", "trakxs");
#$topics = $mydb->db_Count("project_id", "project", "WHERE project_id>0");
#echo $topics;
#$mydb->db_Close();



};#if

	$text.= "
			</div>
		</div>
	";


$title="

";

$ns->tablerender($title, $text, "TRAKXS");



require_once(FOOTERF);

?>
