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

require "connect.php"; 

	
	#$text.=ListRUN($dir_analysis,1,1,1,1);

	#$text.=" s=$s q=$q";

	if ( $_REQUEST["s"] == "Search" ) {
	

		# Init
		$output=""; 
		$RUNList="";
		$SAMPLEList="";
		$RUNNB=0;
		$SAMPLENB=0;
		$ListSample=1;
		#print "$dir_analysis";
		$Directory=$dir_analysis;
	  	$MyDirectory_RUN = opendir($Directory); # or die('Erreur');
		$RUN_SAMPLE_list=array();

		if (join(",",array_keys($user_projects)) == "") {
			$user_projects_ids="-1";
		} else {
			$user_projects_ids=join(",",array_keys($user_projects));
		}; # if

		$query="SELECT run.ref AS run, sample.ref AS sample
			FROM sample
			INNER JOIN project ON (project.id=sample.project_id)
			INNER JOIN run ON (run.id=sample.run_id)
			WHERE (run.ref LIKE '%$q%'
				OR sample.ref LIKE '%$q%'
				#OR run.samplesheet LIKE '%$q%'
				)
			   AND project.id IN ($user_projects_ids)
			"; 
		# AND project.id IN (".join(",",array_keys($user_projects)).")
		#print $query;
		$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));

		#$result_table= "<table class=' '>\n"; 
		#$result_table.= "<tr><td>RUN</td><td>SAMPLE</td><td>#SAMPLE</td></tr>\n"; 
		while ($row = mysql_fetch_assoc($result)) {
			#echo "<pre>"; print_r($row); echo "</pre>";
			$RUN=$row['run'];
			$SAMPLE=$row['sample'];
			$RUN_SAMPLE_list[$RUN][]=$SAMPLE;
			$SAMPLENB++;
		};#while

		# List of all RUNs/SAMPLEs
		foreach ($RUN_SAMPLE_list as $RUN=>$SAMPLES) {
			#echo "RUN=$RUN ";
			$output.="<DIV style='clear: left; margin: 20px;' class=''>".ListSAMPLE($Directory.'/'.$RUN,$RUN,$SAMPLES)."</DIV>";
			$RUNNB++;
			$RUNLink="run.php?action_todo=show_run&run=$RUN&q=$q";
			$SAMPLELink="sample.php?action_todo=show_sample&run=$RUN&sample=$SAMPLE&q=$q";
		};#foreach
		
		#echo "TEST $SAMPLENB $RUNNB";; die();

		# If 1 or 0 SAMPLE
		if (1) {
			if ($SAMPLENB==1) {
				header("Location: $SAMPLELink");
			} elseif ($RUNNB==1) {
				header("Location: $RUNLink");
			} elseif ($SAMPLENB==0 && $RUNNB==0) {
				#$output.="<DIV class='html'>No result. Please contact your favorite bioinformatician to perform the analysis<BR><BR></DIV>";
			};#if
		};#if

		#echo "$SAMPLELink<BR>";
		#die();
		#print "TEST";

		# If no RUN or SAMPLE
		if (1) {


		if (($SAMPLENB+$RUNNB)==0) {

			$q=str_replace("*","%",$q);
			$pattern = '/&&/'; #"|<[^>]+>(.*)</[^>]+>|U"
			$keywords=preg_split($pattern, $q);
			#$keywords = preg_split("/&/", "langage&hypertexte, programmation"); $keywords = preg_split("/[\s,]+/", "langage hypertexte, programmation");
			#echo "<pre>"; print_r($keywords); echo "</pre>";
	
			$filter="";
			if (count($keywords)>=1) {
				$filter=" ";
				foreach ($keywords as $keyword_i=>$keyword) {
					$keyword=trim($keyword);
					#print "$keyword<BR>";
					$pattern = '/^(.*):(.*)-(.*)$/'; #"|<[^>]+>(.*)</[^>]+>|U"
					#$pattern = "/^(\d+):(\d+)-(\d+)$/"; #"|<[^>]+>(.*)</[^>]+>|U"
					preg_match_all($pattern, $keyword, $matches, PREG_PATTERN_ORDER);
					$find_region=false;
					if ($keyword != "" && $matches[0][0]==$keyword && count($matches)==4) {
						#echo "<pre>"; print_r($matches); echo "</pre>";
						$chr=str_replace("chr","",$matches[1][0]);
						$start=$matches[2][0];
						$stop=$matches[3][0];
						$find_region=true;
						$filter_on_region="AND (variant.pos >= '$start' AND variant.pos <= '$stop' AND variant.chr = '$chr')";
						#$filter_on_region="OR (variant.pos >= '$start' AND variant.pos <= '$stop' AND variant.chr = '$chr')";
						#echo $filter_on_region;

					} else {
						/*$filter.=" AND (
						concat(variant.chr,':',variant.pos,variant.ref,'>',variant.alt) LIKE '%$q%'
							OR concat('chr',variant.chr,':',variant.pos,variant.ref,'>',variant.alt) LIKE '%$q%'
						
							$filter_on_region
							)
						";*/
						$filter.=" AND (
							concat(variant.chr,':',variant.pos,variant.ref,'>',variant.alt) LIKE '$keyword'
							OR concat('chr',variant.chr,':',variant.pos,variant.ref,'>',variant.alt) LIKE '$keyword'
							OR concat(variant.chr,':',variant.pos,'|',variant.ref,'>',variant.alt) LIKE '$keyword'
							OR concat('chr',variant.chr,':',variant.pos,'|',variant.ref,'>',variant.alt) LIKE '$keyword'
							OR concat('chr',variant.chr,':',variant.pos,'|',variant.ref,'|',variant.alt) LIKE '$keyword'
							OR concat(variant.chr,':',variant.pos,'|',variant.ref,'|',variant.alt) LIKE '$keyword'
							OR chr='$keyword'
							OR pos='$keyword'
							OR ref LIKE '$keyword'
							OR alt LIKE '$keyword'
							OR variant_annotation.value LIKE '%$keyword%'
							)
						"; # $filter_on_region
					};#if
					
				};#foreach
				#$filter_plus=" )";
				/*OR variant.hgvs LIKE '$keyword'
				OR variant.gene LIKE '$keyword'
				OR variant.genesymbol LIKE '$keyword'
				OR variant.snpid LIKE '$keyword'
				OR variant.location LIKE '$keyword'
				OR variant.outcome LIKE '$keyword'
				*/
			#$filter="AND variant.chr=12";
			
			/*
			$query="SELECT variant.id AS variant_id, variant.chr, variant.pos, variant.ref, variant.alt,
					variant_annotation.value
					
				FROM variant
				INNER JOIN variant_annotation ON (variant_annotation.variant_id=variant.id)
				WHERE variant.ref!='REF' AND variant.alt!='ALT' $filter $filter_on_region
				GROUP BY variant.id
				"; 
			*/
			$query="SELECT variant.id AS variant_id, variant.chr, variant.pos, variant.ref, variant.alt,
					variant_annotation.value
					
				FROM variant
				INNER JOIN variant_annotation ON (variant_annotation.variant_id=variant.id)
				WHERE ref!='REF' AND alt!='ALT' $filter $filter_on_region
				GROUP BY variant.id
				"; 

			#$query="SELECT variant.id AS variant_id, variant.chr, variant.pos, variant.ref, variant.alt, 'value'
			#	FROM variant
			#	WHERE variant.ref!='REF' AND variant.alt!='ALT' AND variant.chr=12
			#	GROUP BY variant.id
			#	LIMIT 100
			#	"; 


			#$q = "7:100-60000000";
			#echo "<pre>"; print($query); echo "</pre>"; 
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			## Mutation
			/*
			$query="SELECT variant.id AS variant_id, variant.chr, variant.pos, variant.ref, variant.alt,
					count(distinct variant_sample.sample_id) AS nb_sample, annotation.value
					
				FROM variant
				INNER JOIN annotation ON (annotation.variant_id=variant.id)
				
				INNER JOIN variant_sample ON (variant.id=variant_sample.variant_id)
				WHERE 1=1 $filter $filter_on_region
				GROUP BY variant.id
				"; 
			*/
			/*
			ann_symbol.value AS genesymbol, ann_location.value AS location, ann_outcome.value AS outcome, ann_hgvs.value AS hgvs, ann_outcome.value AS outcome
			INNER JOIN annotation AS ann_symbol ON (ann_symbol.source='symbol' AND ann_symbol.variant_id=variant.id)
				INNER JOIN annotation AS ann_location ON (ann_location.source='location' AND ann_location.variant_id=variant.id)
				INNER JOIN annotation AS ann_outcome ON (ann_outcome.source='outcome' AND ann_outcome.variant_id=variant.id)
				INNER JOIN annotation AS ann_hgvs ON (ann_hgvs.source='hgvs' AND ann_hgvs.variant_id=variant.id)
				INNER JOIN annotation AS ann_snpid ON (ann_snpid.source='snpid' AND ann_snpid.variant_id=variant.id)
			*/
		
			#echo "<pre>"; print($query); echo "</pre>"; die();
			#$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));

			$result_table= "<table border=1 class=' '>\n"; 
			$result_table.= "<tr>
						<td>VARIANT_ID</td><td>CHR</td><td>POS</td><td>REF</td><td>ALT</td><td>HGVS</td>
						<td>GENESYMBOL</td><td>GENE</td><td>SNPID</td><td>LOCATION</td><td>OUTCOME</td><td>#SAMPLE</td>
					</tr>\n"; 
			$variant_list=array();
			while ($row = mysql_fetch_assoc($result)) {
				$result_array[$row['variant_id']][]=$row;
				$result_table.= "<tr>"
				."<td>".$row['variant_id']."</td>"
				."<td>".$row['chr']."</td>"
				."<td>".$row['pos']."</td>"
				."<td>".$row['ref']."</td>"
			 	."<td>".$row['alt']."</td>"
			 	#."<td>".$row['hgvs']."</td>"
			 	."<td>".str_replace(",","<BR>",$row['hgvs'])."</td>"
			 	."<td>".$row['genesymbol']."</td>"
			 	."<td>".$row['gene']."</td>"
			 	."<td>".$row['snpid']."</td>"
			 	."<td>".$row['location']."</td>"
			 	."<td>".$row['outcome']."</td>"
			 	."<td>".$row['nb_sample']."</td>"
			 	#."<td>".$row['annotations']."</td>" 
			 	."</tr>"; 
			};#while
			$result_table.= "</table>";
			#$output=$result_table;
			#print_r($result_array);

			$variant_files=implode(",",array_keys($result_array));
			# Variant list
			$variant_files_option="";
			if ($variant_files!="") {
				$variant_files_option="--variant_files=$variant_files";
			};#if
			
			# Export from DB to VCF
			# EXPORT (V1 - Only DBexport.pl with muti vcf_id)
			$tmpfnamebase = tempnam("tmp/", "variant_file_");
			#$command="$dir_bin/scripts/export.pl --run=$input_run --samples=$input_sample $variant_files_option --vcf=$tmpfnamebase --annotation_list='ALL' ";
			$command="$dir_bin/scripts/DBexport.pl --config=$myconfig --variant_id=$variant_files --output=$tmpfnamebase.vcf --verbose"; #  --debug 
			#print "$command<BR>";
			$output_exec = shell_exec($command);
			#print "<pre>$output_exec</pre><BR>";
			#die();
			$t=file("$tmpfnamebase.vcf");
			#print "<pre>"; print_r($t); print "</pre>";
			#die();
			
			# Prioritization
			$filters_option=implode(",",$filters);
			$command="$dir_bin/scripts/VCFprioritization.pl --input=$tmpfnamebase.vcf --output=$tmpfnamebase.prioritized.vcf --filter=$filters_option --verbose"; # --debug
			#print "$command<BR>";
			#die();
			
			$output_exec = shell_exec($command);
			#print "<pre>$output_exec</pre><BR>";
		
			$Variants_VCFContent=file_get_contents("$tmpfnamebase.prioritized.vcf");
			#print "<pre>$Variants_VCFContent</pre><BR>";
			$Variants_HTMLContent=VCFtoHTML($Variants_VCFContent,$sample_id,$annotation_list,array("hardfiltering"=>$hardfiltering,"vcf_ids"=>$variant_files,"limit"=>$limit,"orderby"=>$orderby,"ascdesc"=>$ascdesc,"gene_filter"=>$gene_filter,"global_filter"=>$global_filter,"report"=>0));
			
			$text.="<div style='text-align: left; table-border:1; ' class=''>
					$Variants_HTMLContent
					</div>
			";


			};#if	

		};#if
		};#if

		$text .= $output;

	};# Search

	switch ($_REQUEST["action_todo"]) {

		# List of all RUN/SAMPLE
		case "list_run_sample":
		{
			$text.=ListRUN($dir_analysis,1);
			break;
		}#case list_run_sample
	
		# Show a RUN
		case "show_run":
		case "requeue":
		{
			#$text.="<H4>RUN '$input_run'</H4>";
			$text.=ListSAMPLE($dir_analysis."/".$input_run,$input_run,$input_sample);

			# Comment
			$comment="";
			$queued_file="$dir_analysis/$input_run/queued.txt";
			$analysislog_file="$dir_analysis/$input_run/analysis.log";
			$RTAComplete_file="$dir_miseq/$input_run/RTAComplete.txt";
			$MSRanalysislog_file="$dir_miseq/$input_run/AnalysisLog.txt";
			$CompletedJobInfo_file="$dir_miseq/$input_run/CompletedJobInfo.xml";
			$processing_file="$dir_analysis/$input_run/processing.txt";
			
			if (is_file($queued_file) || is_file($processing_file)) {
				$mod_date=date("Y-m-d H:i:s", filemtime($queued_file));
				$queued_file_content=trim(file_get_contents($queued_file));
				$comment.="RUN Analysis Processing...<BR>";
			};#if
			if (is_file($analysislog_file) && !is_file($processing_file)) {
				#$analysislog_file_content=trim(file_get_contents($analysislog_file));
				$mod_date=date("Y-m-d H:i:s", filemtime($analysislog_file));
				$comment.="[$mod_date] RUN Analysis Complete<BR>";
			};#if
			if (is_file($RTAComplete_file)) {
				$RTAComplete_file_content=trim(file_get_contents($RTAComplete_file));
				$mod_date=date("Y-m-d H:i:s", filemtime($RTAComplete_file));
				$comment.="<BR>[$mod_date] RTA Complete";
			} else {
				$comment.="<BR>[$mod_date] RTA Processing...";
			};#if
			if (is_file($MSRanalysislog_file)) {
				$MSRanalysislog_file_summary=implode("<BR>",preg_grep("/>>>/",file($MSRanalysislog_file)));
				$MSRanalysislog_file_steps="&nbsp;&nbsp;&nbsp;&nbsp".implode("<BR>&nbsp;&nbsp;&nbsp;&nbsp",preg_grep("/Step |Elapsed time/",file($MSRanalysislog_file)));
				#echo implode("<BR>",file($MSRanalysislog_file));
				$mod_date=date("Y-m-d H:i:s", filemtime($MSRanalysislog_file));
				if (is_file("$CompletedJobInfo_file")) {
					$comment.="<BR>[$mod_date] MSR Analysis Complete";
				} else {
					$comment.="<BR>[$mod_date] MSR Analysis processing...";
					$comment.="$MSRanalysislog_file_steps";
					
				};#if
			};#if
			if (is_file($analysislog_file)) {
				$analysislog_file_summary=str_replace(array("## $input_run -","ANALYSIS: "),array("#",""),implode("<BR>",preg_grep("/^## $input_run/",file($analysislog_file))));
				$mod_date=date("Y-m-d H:i:s", filemtime($analysislog_file));
				if (is_file($processing_file)) {
					$processing_file_content=trim(file_get_contents($processing_file));
					$comment.="<BR>[$mod_date] TRAKXS Analysis processing...<BR>&nbsp;&nbsp;&nbsp;&nbsp'$processing_file_content'";
				} else {
					$comment.="<BR>[$mod_date] TRAKXS Analysis Complete";
				};#if
				
			};#if
			
			
			

			if ($comment!="") {
				$text.="<DIV class='commentsample' onclick=\"javascript:\"> 
						<DIV id='commentsample_header' class='commentsample_header'>Comment</DIV>
						$comment<BR><BR>
					</DIV>";
				$text.="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";
			};#if

			if (ADMIN) {

				# requeue
				if ($_REQUEST["action_todo"]=="requeue") {
					$command_requeue="echo 'FORCE' > $queued_file";
					$output_exec = shell_exec($command_requeue);
					$comment_plus_message=" Requeued [$command_requeue]<BR>";
				};#requeue

				$text.="<DIV class='commentsample' onclick=\"javascript:\"> 
						<DIV id='commentsample_header' class='commentsample_header'>Admin Commands</DIV>
						<A href='?action_todo=requeue&run=$input_run'>Requeue</A><BR><BR>
						$comment_plus_message<BR>
					</DIV>";
				$text.="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";

			};#if



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

			/*
			$text.="<DIV class='brick file $SampleSheet_csv_selected' onclick=\"javascript:document.location='?action_todo=show_run&run=$input_run&file_to_show=SampleSheet.csv'\">SampleSheet</DIV>";
			$text.="<DIV class='brick file $RunInfo_xml_selected' onclick=\"javascript:document.location='?action_todo=show_run&run=$input_run&file_to_show=RunInfo.xml'\">RunInfo</DIV>";
			$text.="<DIV class='brick file $runParameters_xml_selected' onclick=\"javascript:document.location='?action_todo=show_run&run=$input_run&file_to_show=runParameters.xml'\">runParameters</DIV>";
			$text.="<DIV class='brick file $CompletedJobInfo_xml_selected' onclick=\"javascript:document.location='?action_todo=show_run&run=$input_run&file_to_show=CompletedJobInfo.xml'\">CompletedJobInfo</DIV>";
			$text.="<DIV class='brick file $RunStatistics_xml_selected' onclick=\"javascript:document.location='?action_todo=show_run&run=$input_run&file_to_show=RunStatistics.xml'\">RunStatistics</DIV>";
			$text.="<DIV class='brick file $AnalysisLog_txt_selected' onclick=\"javascript:document.location='?action_todo=show_run&run=$input_run&file_to_show=AnalysisLog.txt'\">AnalysisLog</DIV>";
			$text.="<DIV class='brick file $AnalysisError_txt_selected' onclick=\"javascript:document.location='?action_todo=show_run&run=$input_run&file_to_show=AnalysisError.txt'\">AnalysisError</DIV>";
			*/

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

			break;
		}#case show_run	

		# Show a SAMPLE
		case "show_sample":
		{


			
			# List of filters
			# Construct filter from file
			$filter_file="$dir_bin/scripts/filters.ini";
			$filter_array=parse_ini_file($filter_file,TRUE);
			#echo "Filter_array before:<BR>"; print_r($filter_array); echo "<BR>"; echo "<BR>";
			$filter_array=basedon_filters($filter_array);
			#echo "Filter_array after:<BR>"; print_r($filter_array); echo "<BR>"; echo "<BR>";
			$filter_select="";
			foreach ($filter_array as $filter_name=>$filter_criteria) {
				$filter_description="";
				if ($filter_criteria['description']) {
					$filter_description=$filter_criteria['description'];
				};#if
				$filter_basedon="";
				if ($filter_criteria['basedon']) {
					$filter_basedon="[based on ".$filter_criteria['basedon']." filter]";
				};#if
				$title=trim("$filter_basedon $filter_description");
				$filter_select.="<OPTION value='$filter_name' ".(in_array("$filter_name",$filters)?"selected":"")." title='$title'>$filter_name (".count($filter_criteria)." criteria)\n";
			};#foreach
			#print $vcf_files[0];
			$variants_output=ShowSAMPLE_VariantList($dir_analysis,$input_run,$input_sample,$aligners[0],$callers[0],$annotators[0],$vcf_files[0]);


			#print "$tabfileforannoationlist<BR>";
			$f = fopen($tabfileforannoationlist, 'r'); $line = fgets($f); fclose($f);
			$line_split=split("\t",$line);
			#print "$line<BR>"; print "<BR>";
			#print_r($line_split); print "<BR>";
			$annotation_list_select="";
			$annotation_list_for_order_by_select="";
			$annotation_list_nb=0;
			

			# Annotation list (show and order by) for select
			foreach (split("\t",$line) as $annotation) {
				if (!in_array($annotation,$annotation_list_mandatory)
					&& strpos($annotation, ":") === false) {
					$annotation_list_select.="<OPTION value='$annotation' ".((in_array($annotation,$annotation_list)&&!in_array("ALL",$annotation_list))?"selected":"").">$annotation";
					if (in_array($annotation,$annotation_list)) {
						$annotation_list_nb++;
					}#if
				};#if
				if (!$already_incuded_in_order_by[$annotation]) {
					$annotation_list_for_order_by_select.="<OPTION value='$annotation' ".($orderby==$annotation?"selected":"").">$annotation";
					
				};#if
				
				$already_incuded_in_order_by[$annotation]=1;
			};#foreach
			# case of 'NO' annotation to show
			if (($annotation_list_nb>0 && in_array("NO",$annotation_list)) || in_array("ALL",$annotation_list) ) {
				foreach (array_keys($annotation_list,"NO") as $key) {
					$annotation_list[$key]="";
				};#foreach
			};#if
			# case of FilterScore/FilterFlag/FilterComment not present (fucking bug!!!): just add FilterScore to order by
			#print_r($annotation);
			if (!in_array("FilterScore",$annotation) ) {
				$annotation_list_for_order_by_select="<OPTION value='FilterScore' ".($orderby=="FilerScore"?"selected":"").">FilterScore\n".$annotation_list_for_order_by_select;
			};#if
			

			# List of VCF Files
			$Directory="$dir_analysis/$input_run/$input_sample/";
			$VCFList=array();
			$MyDirectory = opendir($Directory); # or die('Erreur');
			while($Entry = @readdir($MyDirectory)) {
				if(is_file($Directory.'/'.$Entry)) {
					$VCFFile=$Directory.'/'.$Entry;
					$file_parts = pathinfo($VCFFile);
					if (strtolower($file_parts['extension'])=="vcf") {
						$VCFList[$VCFFile]=1;
						#print $VCFFile;
					};#if


				};#if
			}
			closedir($MyDirectory);

			$vcf_files_select="";
			foreach ($VCFList as $vcf_file=>$val) {
				$file_parts = pathinfo($vcf_file);
				$vcf_files_select.="<OPTION value='$vcf_file' ".(in_array($vcf_file,$vcf_files)?"selected":"").">".$file_parts['filename']."</OPTION>\n";
			};#foreach		




			# List of Aligners/Callers
			$Directory="$dir_analysis/$input_run/$input_sample/DATA/";
			$AlignerCallerList=array();
			$MyDirectory = opendir($Directory); # or die('Erreur');
			while($Entry = @readdir($MyDirectory)) {
				if(is_dir($Directory.'/'.$Entry) && $Entry != '.' && $Entry != '..') {
			    		$Aligner=$Entry;
					$Directory2="$dir_analysis/$input_run/$input_sample/DATA/$Aligner/";
					$MyDirectory2 = opendir($Directory2); # or die('Erreur');
					while($Entry2 = @readdir($MyDirectory2)) {
						if(is_dir($Directory2.'/'.$Entry2) && $Entry2 != '.' && $Entry2 != '..') {
					    		$Caller=$Entry2;
							$AlignerCallerList[$Aligner][$Caller]=1;
						};#if
					}
					closedir($MyDirectory2);
				};#if
			}
			closedir($MyDirectory);

			$Aligner_list_authorized=array("MiSeq-Aligner"); # Restriction to MiSeq-Aligner
			#$Aligner_list_authorized=array();
			$Caller_list_authorized=array();
			$Aligner_list_select="";
			$Caller_list_select="";
			#print_r($AlignerCallerList); print "<BR>";
			foreach ($AlignerCallerList as $Aligner=>$Callers_List) {
				if (in_array($Aligner,$Aligner_list_authorized) || empty($Aligner_list_authorized)) { 
					$Aligner_list_select.="<OPTION value='$Aligner' ".(in_array($Aligner,$aligners)?"selected":"").">$Aligner</OPTION>\n";
					foreach ($Callers_List as $Caller=>$value) {
						if (in_array($Caller,$Caller_list_authorized) || empty($Caller_list_authorized)) {
							$Caller_list_select.="<OPTION value='$Caller' ".(in_array($Caller,$callers)?"selected":"").">$Caller</OPTION>\n";
						};#if
					};#foreach
				};#if	
			};#foreach		



			# list of annotation
			#$command="$dir_bin/scripts/VCF.pl --input_file=$vcf_file --output_file=tmp/all.txt --annotation_list='ALL' --output_format=TAB --verbose --debug";
			#echo $command;
			#$output_exec = shell_exec($command);

			#$text.="<H4>RUN '$input_run'</H4>";
			#$text.="<H4>SAMPLE '$input_sample' <small>[<A HREF='?action_todo=show_run&run=$input_run'>$input_run</A>]</small></li></H4>";
			$SAMPLELink="?action_todo=show_sample&run=$input_run&sample=$input_sample";
			$text.="<DIV class='brick run' onclick=\"javascript:document.location='?action_todo=show_run&run=$input_run'\">$input_run</DIV>";
			$text.="<DIV class='brick sample' onclick=\"javascript:document.location='$SAMPLELink'\" title='RUN: $input_run'>$input_sample</DIV>";
			$text.="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";
			$text.="";

			# Comment
			$comment_file="$dir_analysis/$input_run/$input_sample/DATA/comments.txt";
			#$comment_content=nl2br(htmlspecialchars(file_get_contents($comment_file)));
			$comment_array=parse_ini_file($comment_file,TRUE);
			krsort($comment_array);
			#print_r($comment_array); echo "<BR>";
			$comment_array_current=current ($comment_array);
			$comment_array_current_key=key ($comment_array);
			#print_r($comment_array_current); echo "<BR>";
			$comment_current_nl=nl2br(htmlspecialchars($comment_array_current["comment"]));
			#$comment_current=htmlspecialchars($comment_array_current["comment"]);
			#$comment_current=htmlentities($comment_array_current["comment"]);
			
			$comment_current=$comment_array_current["comment"];
			#$input_comment=str_replace("\r\n","<br>",$input_comment);
			
			$comment_header="Comment [$comment_array_current_key]";
			if (trim($comment_current)=="") {
				$comment_current="Click to add a comment (and save it!)	";
			};#if
			if ($comment_array_current_key=="") {
				$comment_header="Comment [no comment yet...]";
			};##if


			if (trim($comment_current)=="rr") {
				$text.="<DIV class='commentsample' onclick=\"javascript:alert('nothing...')\">
						<DIV class='commentsample_header' onclick=\"javascript:\">Click to add a comment [TODO]</DIV>
					</DIV>";
			 	$text.="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";
			} else {
				$text.="

						<DIV class='commentsample' onclick=\"javascript:\"> 
						<DIV id='commentsample_header' class='commentsample_header'>$comment_header
						</DIV>
						<script>
						$(function(){
						   $('#save').click(function () {
							var mysave = $('#textBox').html();
							$('#comment').val(mysave);
							var d = new Date();
							var year = d.getUTCFullYear();
							var month = (d.getMonth()+1); if(month<10) month = '0'+month;
							var day = (d.getDate()); if(day<10) day = '0'+day;
							var hours = d.getUTCHours(); if(hours<10) hours = '0'+hours;
							var minutes = d.getUTCMinutes(); if(minutes<10) minutes = '0'+minutes;
							var seconds = d.getUTCSeconds(); if(seconds<10) seconds = '0'+seconds;
							var dateformatted = year+'-'+month+'-'+day+' '+hours+':'+minutes+':'+seconds;
							document.getElementById('commentsample_header').innerHTML = 'Comment ['+dateformatted+']';
							document.getElementById('date').value = dateformatted;
						    });
						});
						</script>
						<form action='comments.php' target='comments_script' METHOD='POST'>
							<input type='submit' id='save' name='save' value='Save comment' class='btn button search' onclick=\"javascript:var newWin = window.open('comments.php','comments_script','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes');\"/>
							<span style='visibility:hidden;float:right;display:none;'>
								<textarea id='comment' name='comment' style='visibility:hidden;float:right;display:none;'></textarea>
								<input type='hidden' id='date' name='date' style='visibility:hidden;float:right;display:none;'></input>
								<input type='hidden' id='run' name='run' style='visibility:hidden;float:right;display:none;' value='$input_run'></input>
								<input type='hidden' id='sample' name='sample' style='visibility:hidden;float:right;display:none;' value='$input_sample'></input>
								<iframe name='comments_script' style='visibility:hidden;float:right;display:none;'></iframe>
							</span>
							<pre class='prestyle' id='textBox' contenteditable='true' name='textBox'>$comment_current</pre>
						</form>
					</DIV>";
			 	$text.="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";






				#write_php_ini($comment_array, $comment_file.".3.txt");
			};#if
		

$AlignerCallerForm="
<!--<TR>
					<TD >Pipeline</TD>
					<TD>
						<SELECT name='aligners[]'> <!--  multiple size=1 -->
							<!--<OPTION value='MiSeq-Aligner' ".(in_array("MiSeq-Aligner",$aligners)?"selected":"").">MiSeq-Aligner-->
							$Aligner_list_select
						</SELECT>
						<SELECT name='callers[]'> <!--  multiple size=1 -->
							<!--<OPTION value='MiSeq-Caller' ".(in_array("MiSeq-Caller",$callers)?"selected":"").">MiSeq-Caller-->
							<!--<OPTION value='GATK' ".(in_array("GATK",$callers)?"selected":"").">GATK-->
							$Caller_list_select
							<OPTION value='union' ".(in_array("union",$callers)?"selected":"").">Union
							
						</SELECT>
					</TD>
				</TR>-->
";
			# Form filter display
			$text.="
				
				
				<div style='text-align: left; table-border:1; clear:left;' class='form'>
				<form action='?' method='get' class='' name='main' id='main'>
				<INPUT type='hidden' name='action_todo' value='$action_todo'>
				<INPUT type='hidden' name='sample' value='$input_sample'>
				<INPUT type='hidden' name='run' value='$input_run'>
				<INPUT type='hidden' name='q' value='$q'>
				
				<TABLE>
				<TR>
					<TD >VCF Files</TD>
					<TD>
						<SELECT name='vcf_files[]'> <!--  multiple size=1 -->
							<OPTION value=''>Select a VCF File</OPTION>
							$vcf_files_select
						</SELECT>
						
					</TD>
				</TR>
				
				<TR>
					<TD>Filter</TD>
					<TD>
						<SELECT name='filters[]' multiple size=".(count($filter_array)>5?5:count($filter_array)).">
							$filter_select
						</SELECT>
						<INPUT type='checkbox' name='hardfiltering' ".( $hardfiltering ? "checked" : "")."> Hard Filtering (hide 'FILTERED' variants)
						
						<BR><INPUT TYPE='submit' value='Display Filter' name='display_filter' class='btn button search' data-original-title=''  onclick=\"javascript:var newWin = window.open('filters.php','popup','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes');document.getElementById('main').action='filters.php';document.getElementById('main').target='popup';document.getElementById('main').submit();newWin.focus();document.getElementById('main').action='$action_default';document.getElementById('main').target='_self';return false;\"/>
						
				<!-- onclick=\"javascript:\" -->		
					</TD>
				</TR>
				<TR>
					<TD>Order by</TD>
					<TD>
						<SELECT name='orderby'>
							<!--<OPTION value='FilterScore' ".($orderby=="FilterFlag"?"selected":"").">FilterScore	-->
							<!--<OPTION value='FilterFlag' ".($orderby=="FilterFlag"?"selected":"").">FilterFlag	-->
							<!--<OPTION value='QUAL' ".($orderby=="QUAL"?"selected":"").">QUAL	-->
							$annotation_list_for_order_by_select
						</SELECT> 
						<SELECT name='ascdesc'>
							<OPTION value='DESC' ".($ascdesc=="DESC"?"selected":"").">DESC
							<OPTION value='ASC' ".($ascdesc=="ASC"?"selected":"").">ASC		
						</SELECT>
						<BR>
					</TD>
				</TR>
				<TR>
					<TD>Limit</TD>
					<TD>
						<INPUT type='text' name='limit' id='limit' value='$limit' size='10'>
						<BR>
					</TD>
				</TR>
				<TR>
					<TD>Annotation</TD>
					<TD>
						<SELECT name='annotation_list[]' multiple size=".(count($annotation_list)>3?5:count($filter_array)-2).">	
							<OPTION value='NO' ".(in_array("NO",$annotation_list)?"selected":"").">No more annotations 
							<OPTION value='ALL' ".(in_array("ALL",$annotation_list)?"selected":"").">ALL annotations	
							$annotation_list_select	
						</SELECT> 
						<BR>
					</TD>
				</TR>
				<TR>
					<TD></TD>
					<TD align='right'><INPUT TYPE='submit' value='Process' name='process' class='btn button search' data-original-title=''/></TD>
				</TR>
				</TABLE>
				
				
				</div>
				<BR>
				<div style='text-align: left; table-border:1; ' class=''>
				$variants_output
				</form>
				</div>
";

			break;
		}#case show_sample

	};#switch

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
