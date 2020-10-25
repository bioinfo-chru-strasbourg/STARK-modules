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

# Security
if (!USER && 0) {
	$text="<div style='text-align:center'>Please login</div>";
} else {

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
include "search.header.php";


	switch ($_REQUEST["action_todo"]) {

		
		# Show a RUN
		case "show_run":
		case "requeue":
		{
			#$text.="<H4>RUN '$input_run'</H4>";
			$text.=ListSAMPLE($dir_analysis."/".$input_run,$input_run,$input_sample);

			# SampleSheet
			$query="SELECT id, ref, samplesheet
				FROM run
				WHERE run.ref='$input_run'
				"; 
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
		   	

			$run_infos_patterns=array("Experiment Name","Date","Investigator Name","Description","Workflow","Application","Assay","Chemistry");
			while ($row = mysql_fetch_assoc($result)) {
				$samplesheet=$row["samplesheet"];
				#echo "<pre>"; echo str_replace("\\n","</TD></TR><TR></TD>",$samplesheet); echo "</pre>";
				$samplesheet_table="<TABLE><TR><TD>".str_replace(",","</TD><TD valign='top'>",str_replace("\n","</TD></TR><TR><TD valign='top'>",$samplesheet))."</TD></TR></TABLE>";
				foreach ($run_infos_patterns as $run_info_label) {
					
					$tmp=preg_grep("/^$run_info_label,.*$/",split("\n",$samplesheet));
					$run_infos[$run_info_label]=str_replace("$run_info_label,","",array_shift(array_values($tmp)));
					
				};#foreach
			};#while
			$run_infos_table="";
			foreach ($run_infos as $k=>$i) {
				$run_infos_table.="<TR><TD>$k</TD><TD width='20px'></TD><TD>$i</TD></TR>";
			};#foreach
			$comment="<TABLE>$run_infos_table</TABLE>";
			
			

			if ($comment!="") {
				$text.="<DIV class='commentsample' onclick=\"javascript:\"> 
						<DIV id='commentsample_header' class='commentsample_header'>Comment</DIV>
						$comment<BR><BR>
					</DIV>";
				$text.="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";
			};#if


			# QC
			# find all sample_id
			$query="SELECT sample.id AS sample_id, sample.ref, 1, ''
				FROM sample
				INNER JOIN project ON (project.id=sample.project_id)
				INNER JOIN run ON (run.id=sample.run_id)
				WHERE run.ref='$input_run'
				  #AND sample.ref='$input_sample'
				  #AND project.ref='$project'
				"; 
			#print "$query";
			
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
		   	$summary= "<table class=' '>\n"; 
			while ($row = mysql_fetch_assoc($result)) {
				#print_r($row);
				$sample_ids[$row['sample_id']]=$row['sample_id'];
				$sample_id=$row['sample_id'];
				$flag=$row['flag'];
				$comment=$row['comment'];
				$summary.= "<tr>"
					."<td>$sample_id</td>"
					."<td>$ref</td>"
					."<td>$flag</td>"
					."<td>$comment</td>"
				."</tr>"; 
			};#while
			$summary.= "</table>\n"; 
			#echo $summary.join(",",$sample_ids);

			$QC_infos=QC(join(",",$sample_ids),array("sample_interval_summary"));
			#$QC_infos=QC(join(",",$sample_ids));
			#echo "<pre>"; print_r($QC_infos); echo "</pre>"; 
			$QC=$QC_infos["QC"];
			$Files_sample_interval_summary_param_by_aligner=array();
			$Files_sample_interval_summary_param_legends_by_aligner=array();
			$depth_by_aligner=array();
			$depth_failed_by_aligner=array();
			$depth_warning_by_aligner=array();
			$depth_by_QC_id=array();
			#$depth_sum=0;
			#$depth_count=0;
			krsort($QC);
			foreach ($QC as $QC_id=>$infos) {
				foreach ($infos as $infos_k=>$info) {
					#print "$infos_k<BR>";
				};#foreach
				if ($infos["QC_metrics"]=="sample_interval_summary") { #QC_metrics
					#print "<pre>"; print_r($infos); echo "</pre><BR>"; 
					#print "<pre>"; print_r($infos["QC_QC"]); echo "</pre><BR>"; 
					$QC_file_tmp="tmp/QC.$QC_id.".$infos["QC_metrics"];
					file_put_contents($QC_file_tmp,$infos["QC_QC"]);
					$Files_sample_interval_summary_param_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']].="$QC_file_tmp,";
					$Files_sample_interval_summary_param_legends_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']].=$infos['sample_ref'].",";

					# Test QC
					$depth_failed_min=30;
					$depth_warning_min=100;
					
					foreach($infos["QCarray"] as $line => $split ) {
						if ($line>0 && trim($split[0])!="") {
							$chrompos=trim($split[0]);
							$split_chrompos = split(':',$chrompos);
							$chrom=$split_chrompos[0];
							$depth=trim($split[2]);
							$depth_sum+=$depth;
							$depth_count++;
							#print $infos['manifest_ref'].$infos['vcf_aligner'].$chrompos.$QC_id."=".$depth."<BR>";
							$depth_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']][$chrompos][$QC_id]=$depth;
							$depth_by_QC_id[$infos['manifest_ref']][$infos['vcf_aligner']][$QC_id]+=$depth;
							if ($depth<$depth_failed_min) {
								$depth_failed_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']][$chrompos][$QC_id]=$depth;	
							};#if
							if ($depth<$depth_warning_min) {
								$depth_warning_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']][$chrompos][$QC_id]=$depth;	
							};#if
						};#if
					};#foreach

				};#if
				$QCCoverage="jpgraph_coverage.php?coverage_file=".$Files_sample_interval_summary_param_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']]."&legends=".$Files_sample_interval_summary_param_legends_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']]."&title=Coverage of sequenced regions (Manifest=".$infos['manifest_ref'].", Aligner=".$infos['vcf_aligner'].")"; 
				$QCCoverage_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']]=$QCCoverage;
				if (join(",",array_keys($depth_failed_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']]))!="") {
					$QCCoverage_failed_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']]=$QCCoverage."&title=Failed regions&regions=".join(",",array_keys($depth_failed_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']]));
				};#if
				if (join(",",array_keys($depth_warning_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']]))!="") {
					$QCCoverage_warning_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']]=$QCCoverage."&title=Warning regions coverage&regions=".join(",",array_keys($depth_warning_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']]));
				};#if

				if (trim($QCCoverage_default)=="") {$QCCoverage_default=&$QCCoverage_by_aligner[$infos['manifest_ref']][$infos['vcf_aligner']];};
			};#foreach
			
			#DEPTH STATS
			$depth_stats="";
			$depth_stats_array=array();
			
			$depth_stats=$_REQUEST["depth_stats"];
			$depth_stats_zscore=($_REQUEST["depth_stats_zscore"]=="")?3:$_REQUEST["depth_stats_zscore"];
			$depth_stats_covar=($_REQUEST["depth_stats_covar"]=="")?0.6:$_REQUEST["depth_stats_covar"];
			$depth_stats_mindepth=($_REQUEST["depth_stats_mindepth"]=="")?100:$_REQUEST["depth_stats_mindepth"];
			$depth_stats_clusternb=($_REQUEST["depth_stats_clusternb"]=="")?3:$_REQUEST["depth_stats_clusternb"];
			$depth_stats_window=($_REQUEST["depth_stats_window"]=="")?10000:$_REQUEST["depth_stats_window"];


			if ($_REQUEST["depth_stats"]) {
			foreach($depth_by_aligner as $manifest => $split1 ) {	
				#echo "$manifest<BR> ";
				$depth_stats.= "<pre><B>$manifest</B></pre>";
				foreach($split1 as $aligner => $split2 ) {
					#echo "$aligner<BR> ";
					$depth_stats.= "<pre><B>$aligner</B></pre>";


					#print_r(array_search("chr1:43803397-43804497",array_keys($split2)));
					#print_r(array_search("chr1:43811994-43812655",array_keys($split2)));

					foreach($split2 as $chrompos => $split3 ) {
						#echo "$chrompos<BR> ";
						#print_r($split3); echo "<BR>";
						$depth_sum=array_sum($split3);
						$depth_count=count($split3);
						$depth_mean=$depth_sum/$depth_count;
						#$depth_variance=0.0;
						#foreach($split3 as $i ) { $depth_variance+=pow($i-$depth_mean,2); }; # foreach
						#$depth_variance /= $depth_count;
						#$depth_stddev= (float) sqrt($depth_variance/$depth_count);
						$depth_stddev=standard_deviation($split3);
						$depth_covar=$depth_stddev/$depth_mean;
						#$depth_stddev=stats_standard_deviation($split3);
						#echo "<pre>"; print "$chrompos\tdepth_count=$depth_count\tdepth_average=$depth_sum\tdepth_mean=$depth_mean\tdepth_stddev=$depth_stddev\tdepth_covar=$depth_covar<BR>"; echo "</pre>";
						foreach ($split3 as $QC_id=>$depth) {
							#echo "-$QC_id=>$depth<BR> ";
							$sample_ref=$QC_infos["QC"][$QC_id]["sample_ref"];
							#$depth_N=$depth*(1/$depth_count)/($depth/$depth_sum);
							$depth_QC=$depth_by_QC_id[$manifest][$aligner][$QC_id];
							$depth_sum_sum=array_sum($depth_by_QC_id[$manifest][$aligner]); # total depth average for all QC
							$depth_sum_sum_N=$depth_sum_sum-$depth; # total depth average for all QC - the depth of the QCID
							#echo "<pre>";print_r($depth_by_QC_id[$manifest][$aligner]); echo "</pre>";
							#$nb_s=count($sample_ids);
							$depth_N=$depth*(1/$depth_count)/($depth_QC/$depth_sum_sum); # $depth_by_QC_id[$infos['manifest_ref']][$infos['vcf_aligner']][$QC_id]
							#$depth_N=$depth*(1/($depth_count-1))/($depth_QC/$depth_sum_sum_N); # $depth_by_QC_id[$infos['manifest_ref']][$infos['vcf_aligner']][$QC_id]
							#print "$depth_N=$depth*(1/$depth_count)/($depth_QC/$depth_sum_sum)<BR>";
							$depth_mean_N=($depth_sum-$depth)/($depth_count-1);
							$depth_stddev_N=standard_deviation($split3,$QC_id);
							$depth_covar_N=$depth_stddev_N/$depth_mean_N;
							#if ($depth_stddev!=0 &&
							#	( abs($depth_mean-$depth_N)>(5*$depth_stddev) || abs($depth_mean-$depth_N)>(0.8*$depth_mean) )
							#) {
							#$zscore=(($depth_N-$depth_mean)/$depth_stddev);
							$zscore=(($depth_N-$depth_mean_N)/$depth_stddev_N);
							#$zscore=(($depth_N-$depth_mean)/($depth_stddev/sqrt($depth_sum)));
							#if ($depth_stddev!=0 && $depth_N!=0 &&
							#	( abs($depth_mean-$depth)>(0.8*$depth_mean) ) &&
							#	( $depth>100 )
							#) {
							if ($depth_mean_N>=$depth_stats_mindepth # $depth_stddev!=0 &&  $depth_N>=$depth_stats_mindepth && 
								&& abs($zscore)>=$depth_stats_zscore
								#&& abs($depth_mean-$depth)>(0.5*$depth_mean)
								&& $depth_covar_N<$depth_stats_covar
								
							) {
								#if ($depth_covar<0.5) {
								#if () {
									#echo "&& $depth_covar_N<$depth_stats_covar";
									$depth_stats_i= "<pre>";
									$depth_stats_i.= "$sample_ref\t$manifest\t$aligner\t$QC_id\t$chrompos<BR>";
									$depth_stats_i.= "$zscore\tDepth=$depth\tDNorm=$depth_N\tDMean=$depth_mean\tDMeanN=$depth_mean_N\tDSum=$depth_sum\tDStd=$depth_stddev\tDStdN=$depth_stddev_N\tDSCV=$depth_covar\tDSCVN=$depth_covar_N<BR>"; # \tDepthCount=$depth_count
									#$depth_stats_i.= "\n".join("\n",$split3);
									$depth_stats_i.= "</pre>";

									$depth_stats.=$depth_stats_i;
									$depth_stats_array[$sample_ref].=$depth_stats_i;
									$depth_stats_array2[$sample_ref][$chrompos][$manifest][$aligner]["text"].=$depth_stats_i;
									$depth_stats_array2[$sample_ref][$chrompos][$manifest][$aligner]["zscore"]=$zscore;

									#$cluster_i
									#if (!$in_cluster) { # new cluster
									#	$cluster_i++;
									#};#if
									#$depth_stats_array2[$sample_ref][$chrompos][$manifest][$aligner]["cluster"]=$cluster_i;

									#$in_cluster=true;
									

								#};#if
							} else {
								#$in_cluster=false;
							};#if
						};#foreach
					};#foreach
				};#foreach
			};#foreach
			};#if

			#print_r($QCCoverage_by_aligner);
			$coverage_table="<TABLE border=0>";
			$coverage_links="";
			ksort($QCCoverage_by_aligner);
			foreach ($QCCoverage_by_aligner as $manifest=>$QCCoverage_by_al) {
				$coverage_table.="<TR><TD valign='top' rowspan='' colspan='10'>$manifest&nbsp;&nbsp;&nbsp;</TD></TR>";
				#$coverage_table.="<TABLE border=1>";
				foreach ($QCCoverage_by_al as $aligner=>$QCCoverage) {
					#print $QCCoverage;
					$coverage_table.="<TR><TD valign='top' width='80'></TD>";
					$coverage_links.="&nbsp;&nbsp;&nbsp;<span class=' ' onclick=\"javascript: var newWin = window.open('$QCCoverage','cov','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>$manifest.$aligner</span>&nbsp;&nbsp;&nbsp;"; 
					$coverage_table.="<TD valign='top'>$aligner&nbsp;&nbsp;&nbsp;</TD>";
					if (trim($QCCoverage)!="") {
						$coverage_table.="<TD valign='top'><span class=' ' onclick=\"javascript: var newWin = window.open('$QCCoverage','cov','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>ALL</span>&nbsp;&nbsp;&nbsp;</TD>";
					} else {
						$coverage_table.="<TD valign='top'>noCOV&nbsp;&nbsp;&nbsp;</TD>";
					};#if
					if (trim($QCCoverage_failed_by_aligner[$manifest][$aligner])!="") {
						$coverage_table.="<TD valign='top'><span class=' ' onclick=\"javascript: var newWin = window.open('".$QCCoverage_failed_by_aligner[$manifest][$aligner]."','cov','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>FAIL</span>&nbsp;&nbsp;&nbsp;</TD>";
					} else {
						$coverage_table.="<TD valign='top'>noFAIL&nbsp;&nbsp;&nbsp;</TD>";
					};#if
					if (trim($QCCoverage_warning_by_aligner[$manifest][$aligner])!="") {
						$coverage_table.="<TD valign='top'><span class=' ' onclick=\"javascript: var newWin = window.open('".$QCCoverage_warning_by_aligner[$manifest][$aligner]."','cov','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>WARN</span>&nbsp;&nbsp;&nbsp;</TD>";
					
					} else {
						$coverage_table.="<TD valign='top'>noWARN&nbsp;&nbsp;&nbsp;</TD>";
					};#if
					$coverage_table.="</TD></TR>";
				};#foreach
				#$coverage_table.="</TABLE>";
				#$coverage_table.="</TD></TR>";
			};#foreach
			$coverage_table.="</TABLE>";

			$text.="<DIV class='commentsample' onclick=\"javascript:\"> 
					<DIV id='cov_header' class='commentsample_header'>Infos</DIV>
					$coverage_table<BR>
					<iframe name='cov' width='100%' height=650 frameborder=0 src='$QCCoverage_default'></iframe>
					<BR><BR>
				</DIV>";
			$text.="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";

			# Depth stats
			$depth_stats2="";
			$depth_stats_cluster="";
			# clusters
			$cluster_i=1;
			$in_cluster=false;
			$chrompos_i_prev=-10000000;
			$chr_prev=-100000000;
			$start_prev=-100000000;
			$stop_prev=-10000000000;
			$window=$depth_stats_window;
			$clusters=array();
			foreach ($depth_stats_array2 as $sample_ref=>$chromposs) {
				foreach ($chromposs as $chrompos=>$manifests) {
					
					$chrompos_i=array_search($chrompos,array_keys($split2));
					$chrompos_infos=preg_match('/^chr([0-9]*):([0-9]*)-([0-9]*)/',$chrompos,$matches);
					#echo "<pre>$chrompos "; print_r($matches); echo "</pre>";
					$chr=$matches[1];
					$start=$matches[2];
					$stop=$matches[3];
					
					#if (($chrompos_i-1)!=$chrompos_i_prev) { # same cluster
					#if (($chrompos_i-2)>$chrompos_i_prev) { # same cluster
					#print "<pre>$chrompos\t$chrompos_prev\t$start-$stop_prev(".($start-$stop_prev).")<=?$window || $start_prev-$stop(".($start_prev-$stop).")<=?$window</pre><BR>";
					if ($chr==$chr_prev
						&& ( 	   (abs($start-$stop_prev) <= $window)
							|| (abs($start_prev-$stop) <= $window)
						) 
					) { # same cluster
					} else {
						$cluster_i++;
					};#if
					$clusters[$sample_ref][$cluster_i][$chrompos]=$chrompos_i;
					$clusters_chrompos[$sample_ref][$chrompos]=$cluster_i;
					$chrompos_i_prev=$chrompos_i;
					$chrompos_prev=$chrompos;
					$chr_prev=$chr;
					$start_prev=$start;
					$stop_prev=$stop;
				};#foreach
			};#foreach
			foreach ($clusters as $sample_ref=>$clusters) {
				
				$depth_stats_cluster_by_sample="";
				$depth_stats_cluster_by_sample_nb_cluster=0;
				foreach ($clusters as $cluster=>$chromposs) {
					if (count($chromposs)>=$depth_stats_clusternb) {
						#print_r($chromposs);
						$depth_stats_cluster_by_sample_nb_cluster++;
						$depth_stats_cluster_by_sample.="<pre>\tCluster #$cluster\t".join(" ",array_keys($chromposs))."\n";
						$depth_stats_cluster_by_sample.="</pre>";
					};#if
				};#foreach
				if ($depth_stats_cluster_by_sample_nb_cluster>0) {
					$depth_stats_cluster.="<BR><pre>$sample_ref ($depth_stats_cluster_by_sample_nb_cluster clusters)</pre>$depth_stats_cluster_by_sample";
				};#if
			};#foreach

			foreach ($depth_stats_array2 as $sample_ref=>$chromposs) {
				$depth_stats2.="<BR><pre>$sample_ref</pre>";
				foreach ($chromposs as $chrompos=>$manifests) {
					$chrompos_i=array_search($chrompos,array_keys($split2));
					$depth_stats2.="<pre>\t$chrompos\n";
					foreach ($manifests as $manifest=>$aligners) {
						foreach ($aligners as $aligner=>$infos) {
							$zscore=$infos["zscore"];
							$cluster=$clusters_chrompos[$sample_ref][$chrompos];
							$infos_text=$infos["text"];
							$depth_stats2.="\t\t$zscore\t$cluster\t$manifest\t$aligner\t\n";
							#$depth_stats2.="\t\t$infos_text\t\n";
						};#foreach
					};#foreach
					$depth_stats2.="</pre>";
				};#foreach
			};#foreach

			if ($_REQUEST["depth_stats"]) {
			$text.="<DIV class='commentsample' onclick=\"javascript:\"> 
					<DIV id='depth_stats_header' class='commentsample_header'>Depth Stats</DIV>
					$depth_stats_cluster
					<BR><BR>
					$depth_stats2
					<!--".join("<BR><BR>",$depth_stats_array)."-->
					<!--$depth_stats-->
					<BR><BR>
				</DIV>";
			$text.="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";
			};#if


			if ($samplesheet!="") {
				$text.="<DIV class='commentsample' onclick=\"javascript:\"> 
						<DIV id='commentsample_header' class='commentsample_header'>SampleSheet</DIV>
						$samplesheet_table<BR><BR>
					</DIV>";
				$text.="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";
			};#if


			/*

			if ($file_to_show=="") {
				$file_to_show="SampleSheet.csv";
			};#if

			$file_to_show_selected=str_replace(".","_",$file_to_show)."_selected";
			$$file_to_show_selected="fileselected";
			#echo $file_to_show_selected;
			$files_to_provide=array("SampleSheet.csv","RunInfo.xml","runParameters.xml","CompletedJobInfo.xml","RunStatistics.xml","AnalysisLog.txt","AnalysisError.txt",);
			foreach ($files_to_provide as $file_to_provide) {
				if (is_file("$dir_miseq/$input_run/$file_to_provide")) {
					if ($file_to_provide==$file_to_show) { $class_file_to_provide="fileselected"; } else { $class_file_to_provide=""; };#if
					$path_parts = pathinfo($file_to_provide);
					$file_to_provide_show=$path_parts['filename']; #.".".$path_parts['extension'];
					$text.="<DIV class='brick file $class_file_to_provide' onclick=\"javascript:document.location='?action_todo=show_run&run=$input_run&file_to_show=$file_to_provide'\">$file_to_provide_show</DIV>";
				};#if
			};#foreach


			$text.="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";
			
			if ($file_to_show=="RunStatistics.xml") {
				foreach (glob("$dir_miseq/$input_run/*RunStatistics.xml") as $filename) {
					$path_parts = pathinfo($filename);
					$file_to_show=$path_parts['filename'].".".$path_parts['extension'];
				};#foreach
			};#if
			
			if (is_file("$dir_miseq/$input_run/$file_to_show")) {
				$file="$dir_miseq/$input_run/$file_to_show";
				if (is_file($file)) {
					$path_parts = pathinfo($file);
					$file_to_show_filename=$path_parts['filename'];
					
					$file_content=htmlspecialchars(file_get_contents($file));
					$text.="
						<!--<DIV class='brick file'>$file_to_show_filename</DIV><DIV style='margin-bottom:0px; clear:left;'><BR></DIV>-->
						<pre><small>$file_content</small></pre>
						";
				} else {
					$text.="No '$file' file...";
				};#if
			};#if

			*/

			break;
		}#case show_run	


	};#switch

};#if

	$text.= "</div></div>";

$title="

";

$ns->tablerender($title, $text, "TRAKXS");



require_once(FOOTERF);

?>
