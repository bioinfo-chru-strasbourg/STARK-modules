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

if ($_REQUEST["action_todo"]!="show_sample_report" || 0) {
	require_once(HEADERF);
};#if
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


#
# Annotation List
#

$annotation_list=$_REQUEST["annotation_list"];
# Mandatory list
$annotation_list_mandatory=array("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FilterScore", "FilterFlag", "FilterComment", "geneSymbol", "location", "outcome", "GQ", "BQ", "DP", "AD", "VF", "AF", "FA", "dbSNP", "dbSNPNonFlagged","hgvs","outcomeinfos","variant_id","AlleleFrequency");
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
$text.="<link rel='stylesheet' media='print' type='text/css' href='print.css' />";
/*<link rel="stylesheet"
	type="text/css"
	media="print" href="print.css" />
*/
#$text.="<link rel='stylesheet' media='all' type='text/css' href='../../e107_web/css/e107.css?0' />";
#$text.="<link rel='stylesheet' media='all' type='text/css' href='../../e107_web/js/bootstrap/css/bootstrap.min.css?0' />";
#$text.="<link rel='stylesheet' media='all' type='text/css' href='../../e107_web/js/bootstrap/css/tooltips.css?0' />";
#$text.="<BODY style='font-size:6;'>";

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
if ($_REQUEST["action_todo"]!="show_sample_report"){
include "search.header.php";
};#


	switch ($_REQUEST["action_todo"]) {

			
		# Show a SAMPLE
		case "show_sample":
		case "show_sample_report":
		{

			# Sample info
			if (1) { # Based on DB
			#$query="SELECT sample.id, sample.ref, sample.flag, sample.comment
			#	FROM sample
			#	INNER JOIN project ON (project.project_id=project.project_id)
			#	WHERE sample.run='$input_run'
			#	  AND sample.ref='$input_sample'
			#	  AND project.ref='$project'
			#	"; 
			$query="SELECT sample.id AS sample_id, sample.ref, 1, ''
				FROM sample
				INNER JOIN project ON (project.id=sample.project_id)
				INNER JOIN run ON (run.id=sample.run_id)
				WHERE project.id IN (".join(",",array_keys($user_projects)).")
				  AND run.ref='$input_run'
				  AND sample.ref='$input_sample'
				  #AND project.ref='$project'
				"; 
			
			#echo "<pre>"; print "$query"; echo "</pre>";
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
		   	$summary= "<table class=' '>\n"; 
			while ($row = mysql_fetch_assoc($result)) {
				#print_r($row);
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
			#echo $summary;
			};#if
			#print "SAMPLE ID:".$sample_id;
			if ($sample_id=="") {
				$text.="<div class='brick run' style='text-align:'>$input_run</div>";
				$text.="<div class='brick sample' style='text-align:'>$input_sample</div>";
				$text.="<div class='brick warning' style='text-align:'>Access denied!</div>";
				break;
			};#if
			# die();
			# MANIFEST
			$query="SELECT sample.id AS sample_id, manifest.ref AS manifest_ref, manifest.manifest AS manifest
				FROM sample
				INNER JOIN manifest ON (manifest.id=sample.manifest_id)
				WHERE sample.id=$sample_id
				"; 
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
		   	$summary= "<table class=' '>\n"; 
			while ($row = mysql_fetch_assoc($result)) {
				#print_r($row);
				$manifest_ref=$row['manifest_ref'];
				$manifest=$row['manifest'];
				
			};#while

			#
			#print_r($test);

			#print "$manifest_ref<BR>";
			#print "<pre>$manifest</pre><BR>";

			# List of filters
			# Construct filter from file
			#$filter_file="$dir_bin/scripts/filters.ini";
			$filter_file="$dir_howard/config.filter.ini";
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


			# Variant List
			$tmpfnamebase = tempnam("tmp/", "variant_file_");
			
			# Exportation
			$variant_files=implode(",",$vcf_files);

			# flags
				$query="SELECT flags.id, flags.label
					FROM flags
					WHERE flags.type='variant_sample'
					  and score=0
					"; 
				#print "<pre>$query</pre><BR>";
				$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			   	$summary= "<table class=' '>\n"; 
				$nb_total_variant_id=0;
				$flag_0_label="unknown";
				$flag_0_id="0";
				while ($row = mysql_fetch_assoc($result)) {
					$flag_0_label=$row['label'];
					$flag_0_id=$row['id'];
				};#while
			#echo "$flag_0_label;$flag_0_id;<BR>";

			if ($_REQUEST["action_todo"]=="show_sample_report" && 1) {

				# Validation
				$query="SELECT variant.chr, variant.pos, variant.ref, variant.alt, variant_sample.variant_id AS variant_id, variant_annotation.value, flags.label, flags.score, flags.id, comment.comment, annotation.source, annotation.release
					FROM variant_sample
					LEFT OUTER JOIN comment ON (comment.id=variant_sample.id AND comment.type='variant_sample')
					LEFT OUTER JOIN variant ON (variant.id=variant_sample.variant_id)
					LEFT OUTER JOIN variant_annotation ON (variant_annotation.variant_id=variant.id)
					LEFT OUTER JOIN annotation ON (annotation.id=variant_annotation.annotation_id)
					LEFT OUTER JOIN flags ON (flags.id=comment.flag_id AND flags.type=comment.type)
					WHERE variant_sample.sample_id=$sample_id
					ORDER BY flags.score DESC
					"; 
				#print "<pre>$query</pre><BR>";
				# LEFT OUTER JOIN annotation ON (annotation.id=variant_annotation.annotation_id AND annotation.source='hgvs')
				$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			   	$summary= "<table class=' '>\n"; 
				$value="";
				while ($row = mysql_fetch_assoc($result)) {
					#print "<pre>"; print_r($row); print "</pre>"; 
					$variant_id=$row['variant_id'];
					$variant_chr=$row['chr'];
					$variant_pos=$row['pos'];
					$variant_ref=$row['ref'];
					$variant_alt=$row['alt'];
					$label=($row['label']=="")?"$flag_0_label":$row['label'];
					$score=($row['score']=="")?"0":$row['score'];
					$flag_id=($row['id']=="")?"$flag_0_id":$row['id'];
					#$flag_id=$row['id']+0;
					#echo "$flag_id; $label; $score;<BR>";
					$flag_label[$flag_id]=$label;
					$flag_score[$flag_id]=$score;
					$source=str_replace(",","<BR>",$row['source']);
					if (trim(strtolower($source))=="hgvs") {
						#$value=str_replace(",","<BR>",$row['value']);
						$value=$row['value'];
						$validation_table[$flag_id][$variant_id]["hgvs"]=$value;
						
					};#if
					$comment=$row['comment'];
					$validation_table[$flag_id][$variant_id]["comment"]=$comment;
					$validation_table[$flag_id][$variant_id]["genomic"]="chr$variant_chr:$variant_pos$variant_ref>$variant_alt";
					#$validation_table_variant_id[$variant_id]=$value;
					
				};#while
				krsort($validation_table);
				#print "<pre>"; print_r($validation_table); print "</pre><BR>";
				$variant_checked_summary_tmp="";
				foreach ($validation_table as $flag_id=>$variant_ids) {
					$label=$flag_label[$flag_id];
					$summary="<TR><TD>$label</TD><TD width='10px'>&nbsp;</TD><TD>".count($variant_ids)."</TD></TR>";
					
					if ($flag_score[$flag_id]>0) {
						$variant_checked_tmp.="<TR  valign='top' halign='left' style='font-weight:bold;'>
									<TD colspan='10'>$label</TD>
									</TR>";
						$variant_checked_summary_tmp=$summary.$variant_checked_summary_tmp;
						foreach ($variant_ids as $variant_id=>$infos) {
							$genomic=str_replace(",","<BR>",$infos['genomic']);
							#$hgvs=str_replace(",","<BR>",$infos['hgvs']);
							$customHGVS=find_hgvs($infos['hgvs'],"",$user_groups,$customNM_infos);
							$customHGVS=str_replace(",","<BR>",$customHGVS);
							#print $infos['hgvs']."";
							#print_r($user_groups); print_r($customNM_infos);
							$hgvs=$customHGVS;
							$comment=str_replace("\n","<BR>",$infos["comment"]);
							$variant_checked_tmp.="<TR  valign='top' halign='left' style='font-weight:;'>
									
									<TD width='10px'>#&nbsp;</TD>
									<TD width=''>$genomic<BR>$hgvs</TD>
									<!--<TD width=''>$hgvs<BR>$genomic</TD>-->
									<TD width='10px'>&nbsp;</TD>
									<TD colspan=1 width=''></TD>									
									<TD width='10px'>&nbsp;</TD>
									<TD colspan=1 width='99%'>$comment</TD>									
									
									";
						};#foreach
					} else {
						$variant_checked_summary_tmp.=$summary;
					}#if
				};#foreach

				#echo "<pre>"; print_r($validation_table); echo "</pre>";
				/*$variant_checked_header="<TR  valign='top' halign='left' style='font-weight:bold'>
									<TD width='' colspan='3'>#VID / HGVS</TD>
									<TD width='10px'>&nbsp;</TD>
									<TD width=99%>Validation Comment</TD>
									
									</TR>";
				*/
				$variant_checked_summary_header="<TR  valign='top' halign='left' style='font-weight:bold'>
									<TD width=''>Validation</TD>
									<TD width='10px'>&nbsp;</TD>
									<TD width=''>Number of variants</TD>
									</TR>";
				$variant_checked="<TR>
									<TD width='20px'>&nbsp;</TD>
									<TD width='99%' valign='top'>
										<TABLE valign='left' halign='top' border=0 width='100%'>$variant_checked_header$variant_checked_tmp</TABLE>
										<BR>
										<TABLE valign='left' halign='top'>$variant_checked_summary_header$variant_checked_summary_tmp</TABLE>
										<BR>
									</TD>
									<TD valign='top'>
										
									</TD>
									<TD width='20px'>&nbsp;</TD>
								</TR>";
				$Variants_HTMLContent="<DIV class=''>
				<TABLE width='95%' >
					<TR>
						<TD valign='top' class='validation' width='20%'>
						
							<TABLE width='100%' border=0 class='' style='padding:10px'>
								<TR>
									<TD width='20px'>&nbsp;</TD>
									<TD colspan=4 width='300px' height='40px' valign='middle'>
										<B>Variant Validation</B>
									</TD>
									<TD width='20px'>&nbsp;</TD>
								</TR>
								$variant_checked
							</TABLE>
							</DIV>
						</TD>

					</TR>
				</TABLE>";

				

			} elseif ($variant_files!="") {

				#print $variant_files;
				$variant_files_option="";
				if ($variant_files!="") {
					$variant_files_option="--variant_files=$variant_files";
				};#if
			
				# Export from DB to VCF
				#$command="$dir_bin/scripts/export.pl --run=$input_run --samples=$input_sample $variant_files_option --vcf=$tmpfnamebase --annotation_list='ALL' ";

			
				# EXPORT (V1 - Only DBexport.pl with muti vcf_id)
				$timestart=microtime(true);
				if (1) {
					#$command="$dir_bin/scripts/DBexport.pl --config=$myconfig --vcf_id=$variant_files --output=$tmpfnamebase.vcf --verbose"; #  --debug --config_annotation=$myAnnotationConfig 
					$command="$dir_pepper/DBexport.pl --config=$myconfig --vcf_id=$variant_files --output=$tmpfnamebase.vcf --verbose"; #  --debug --config_annotation=$myAnnotationConfig 
					#print "$command<BR>";
					$output_exec = shell_exec($command);
					#print "<pre>$output_exec</pre><BR>";
					#$t=file_get_contents("$tmpfnamebase.vcf");
					$t=file("$tmpfnamebase.vcf");
					#print "<pre>"; print_r($t); print "</pre>"; 
				};#if

				# EXPORT (V2 - DBexport.pl with ONE vcf_id each and VCFTOOLS VCF-MERGE)
				if (0) {
					foreach (split(",",$variant_files) as $k=>$vcf_id) {
						#print "<pre>$k=>$vcf_id</pre>";
						$tmpfnamebase_tmp=tempnam("tmp/", "variant_file_");
						$tmpfnamebase_array[$vcf_id] = $tmpfnamebase_tmp.".vcf";
						$command="$dir_pepper/DBexport.pl --config=$myconfig --vcf_id=$vcf_id --output=$tmpfnamebase_tmp.vcf --verbose"; #  --debug --config_annotation=$myAnnotationConfig 
						#print "$command<BR>";
						$output_exec = shell_exec($command);
						#print "<pre>$output_exec</pre><BR>";
						#$t=file_get_contents("$tmpfnamebase_tmp.vcf");
						$t=file("$tmpfnamebase_tmp.vcf");
						#print "<pre>"; print_r($t); print "</pre>";
					};#foreach
					$tmpfnamebase_merge=tempnam("tmp/", "variant_file_");
					$command="$dir_pepper/VCFmerge.sh \"".join(" ",$tmpfnamebase_array)."\" $tmpfnamebase_merge.vcf 1>$tmpfnamebase_merge.vcf.log 2>$tmpfnamebase_merge.vcf.err "; #  --debug --config_annotation=$myAnnotationConfig 
					#print "$command<BR>";
					$output_exec = shell_exec($command);
					#print "<pre>$output_exec</pre><BR>";
					#$t=file_get_contents("$tmpfnamebase_tmp.vcf");
					$tlog=file("$tmpfnamebase_merge.vcf.log");
					#print "FINAL<pre>"; print_r($tlog); print "</pre>"; 
					$terr=file("$tmpfnamebase_merge.vcf.err");
					#print "FINAL<pre>"; print_r($terr); print "</pre>"; 
					$t=file("$tmpfnamebase_merge.vcf");
					#print "FINAL<pre>"; print_r($t); print "</pre>"; 
					$tmpfnamebase=$tmpfnamebase_merge;
				};#if

				$timestop=microtime(true);
				$timeDBexport=$timestop-$timestart;
				#echo "Time DBexport: $timeDBexport<BR>";
				
				#echo "<pre>"; print_r(file("$tmpfnamebase.vcf")); echo "</pre>";
				# Prioritization
				$filters_option=implode(",",$filters); 
				$command="$dir_howard/VCFprioritization.pl --input=$tmpfnamebase.vcf --output=$tmpfnamebase.prioritized.vcf --filter=$filters_option --verbose --debug "; # --debug
				#print "$command<BR>";
				$output_exec = shell_exec($command);
				#print "<pre>$output_exec</pre><BR>";
				
			
				$Variants_VCFContent=file_get_contents("$tmpfnamebase.prioritized.vcf");
				#print "hardfiltering=>$hardfiltering,vcf_id=>$variant_files,limit=>$limit,orderby=>$orderby,ascdesc=>$ascdesc,gene_filter=>$gene_filter,global_filter=>$global_filter,report=>0<BR>";
				
				#print "<pre>$Variants_VCFContent</pre><BR>";
				$Variants_HTMLContent=VCFtoHTML($Variants_VCFContent,$sample_id,$annotation_list,array("hardfiltering"=>$hardfiltering,"vcf_ids"=>$variant_files,"limit"=>$limit,"orderby"=>$orderby,"ascdesc"=>$ascdesc,"gene_filter"=>$gene_filter,"global_filter"=>$global_filter,"report"=>0));

			};#if
				
			#$Variants_HTMLContent="";
			

			# Generate Files
			#$variantList=ShowSAMPLE_VariantList($dir_analysis,$input_run,$input_sample,$aligners[0],$callers[0],$annotators[0],$tmpfnamebase);
			#$variants_output=$variantList[0];
			#$variants_output=ShowSAMPLE_VariantList($dir_analysis,$input_run,$input_sample,$aligners[0],$callers[0],$annotators[0],$vcf_files[0]);
			

			# VCF to HTML
			
			#$tmpfnamebase=$variantList[1];
			#$VCFFile="$tmpfnamebase.vcf";
			#print $VCFFile;
			#$Variants_VCFContent=file_get_contents($VCFFile);
			#$f = file($VCFFile, 'r');
			#print_r($f); print "<BR>";
			#print "<pre>$Variants_VCFContent</pre>";
			#$Variants_HTMLContent=VCFtoHTML($Variants_VCFContent,$sample_id,$annotation_list);
			

			# Export files
			#$ExportFiles="".$variantList[1]." | ".$tmpfnamebase."<BR>";
			$ExportFiles.="<DIV class='brick file' onclick=\"javascript:var newWin = window.open('tmp/".basename($tmpfnamebase).".prioritized.vcf','EXPORT_VCF','menubar=yes, status=no, scrollbars=no, menubar=yes, width=800, height=600, visible=yes');newWin.focus()\">VCF</DIV> ";
			$ExportFiles.="<DIV class='brick file' onclick=\"javascript:var newWin = window.open('tmp/".basename($tmpfnamebase).".txt','EXPORT_TRAKXS','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes');newWin.focus()\">TRAKXS</DIV> ";
			$ExportFiles.="<DIV class='brick file' onclick=\"javascript:var newWin = window.open('tmp/".basename($tmpfnamebase).".html','EXPORT_HTML','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes');newWin.focus()\">HTML</DIV> ";
			#$ExportFiles.=$variantList[2];


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
			#print_r($annotation_list);

			# List of VCF Files
			if (0) { # based on files
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
			};#if
			

			#sample_interval_statistics
			# List of QC Files
			if (1) { # based on DB

				#$FASTQC_metrics_name="fastqc_data.txt";
				$FASTQC_metrics_name=array("FASTQC","fastqc_data.txt");

				$QC_metrics_array=array("sample_interval_summary","HsMetrics","FASTQC","fastqc_data.txt");
				#$QC_metrics_array=array("sample_interval_summary","FASTQC");
				$QC_infos=QC($sample_id,$QC_metrics_array,$FASTQC_metrics_name);
				$QC=$QC_infos["QC"];
				$QCFASTQC_id=$QC_infos["QCFASTQC_id"];
				$QCMetrics_array=$QC_infos["QCMetrics_array"];
				#echo "<PRE>"; print_r($QC_infos); echo "</PRE>";
				#echo "<PRE>"; print_r($QCMetrics_array); echo "</PRE>";
				$warnings=$QC_infos["warnings"];
				$module_for_warnings=$QC_infos["module_for_warnings"];
				#print_r($QC);
				#print_r($QCFASTQC_id);
				
				#print "$summary<BR>";
				#die();	
				$tmpfnamebase_QC = tempnam("tmp/", "QC_".$input_run."_".$input_sample."_");
				#print "$QCFASTQC_id";
				#print "<pre>"; print_r($QC); print "</pre>";
				$Directory="$dir_analysis/$input_run/$input_sample";
				


				#Q30
				#$QCMetrics_array[$BAMFile_name][$variable_name]=$variable_value;
				#$sequence_total=$QCFASTQFile_array[6][1];
				#print_r($QC[$QCFASTQC_id]["QCarray"]); # $QC[$QC_id]["QCarray"]
				#print_r($QC[$QCFASTQC_id]["QCarray"][6]); # $QC[$QC_id]["QCarray"]

				#FASTQC Modules
				$FASTQC_modules=array();
				foreach ($QC[$QCFASTQC_id]["QCarray"] as $line=>$cols) {
					# Line match
					$start_module_match = !preg_grep("/^>>.*/i",array($cols[0]),PREG_GREP_INVERT);
					$end_module_match = !preg_grep("/^>>END_MODULE/i",array($cols[0]),PREG_GREP_INVERT);
					$comment_match = !preg_grep("/^#/i",array($cols[0]),PREG_GREP_INVERT);
					# Module Name
					if ($start_module_match && !$end_module_match) {
						$module_name=$cols[0];
						$module_name=preg_replace(array("/^>>/"), array(""), $cols[0]);
					};#if
					foreach ($cols as $col=>$value) {
						# Module Hearder
						if (!$start_module_match && !$end_module_match && $comment_match) {
							#echo "$line $col $value<BR>";
							$value=preg_replace(array("/^#/"), array(""), $value);
							$FASTQC_modules_header[$module_name][$col]=$value;
						};#if
						# Module Values
						if (!$start_module_match && !$end_module_match && !$comment_match) {
							#echo "$line $col $value<BR>";
							$FASTQC_modules[$module_name][$line][$FASTQC_modules_header[$module_name][$col]]=$value;
						};#if
					};#foreach
				};#foreach

				#Per sequence quality scores
				$sequence_quality_scores_pattern="Per sequence quality scores";
				$Quality_total=0;
				$Quality_30=0;
				$Quality_cumul=array();
				$Quality_below=array();
				foreach ($FASTQC_modules[$sequence_quality_scores_pattern] as $line=>$values) {
					$Quality_below[$values["Quality"]]=$Quality_total;
					$Quality_total+=$values["Count"];
					$Quality_cumul[$values["Quality"]]=$Quality_total;
				};#foreach

				#print_r($FASTQC_modules);
				#print_r($Quality_total);
				#print_r($Quality30_cumul);
				$Quality_Scores=array();
				foreach ($Quality_cumul as $quality=>$cumul) {
					$Quality_Scores[$quality]=($Quality_total-$Quality_below[$quality])/$Quality_total;
				};#foreach
				#print_r($Quality_Scores);

				# Q30
				$sequence_q30=$Quality_Scores["30"];
				
				# QALL
				$sequence_q="# Quality Scores\n";
				foreach ($Quality_Scores as $quality=>$score) {
					$sequence_q.="Q$quality=".round($score*100,2)."\n";
				};#foreach
				#echo "Q30=$sequence_q30<BR>";
				
				#$sequence_total=$QC[$QCFASTQC_id]["QCarray"][6][1];
				#$sequence_q30_count=0;
				#for ($i=82; $i<=90; $i++) {
				#	$sequence_q30_count+=$QC[$QCFASTQC_id]["QCarray"][$i][1];
				#};#if
				#$sequence_q30=($sequence_q30_count/$sequence_total);
				# Warning QC
				if ($sequence_q30<0.95) {
					$warnings["Q30"][]="WARN Q30";
				};#if
				if ($sequence_q30<0.90) {
					$warnings["Q30"][]="FAIL Q30";
				};#if

				# QCData
				$QCData="<TABLE class=''>
						<TR>
							
							<TD title='$sequence_q'><!--<span class=' ' onclick=\"javascript:var newWin = window.open('$QCFASTQLink','QC_FASTQ','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=1200, height=800, visible=yes');newWin.focus()\" style='cursor:pointer;'>[+]</span>--> Q30=".round($sequence_q30*100,2)."%</TD>
							<TD></TD>
							<TD width='20px'>&nbsp;</TD>
							<TD width='20px'>&nbsp;</TD>
							</TR>
						
					</TABLE>
				";
				
				# Summary 
				$QCHsMetrics_array=array();
				$QCSummary_array=array();
				#$QCHsMetrics_variable_array=array("GENOME_SIZE","TOTAL_READS","PF_READS","PF_UQ_READS_ALIGNED","PCT_PF_UQ_READS_ALIGNED");
				$QCSummary="";

				# Per Sequence Quality link
				$per_sequence_quality_scores_file="";
				$per_sequence_quality_scores_file_legends="";
				foreach ($QC as $QC_id=>$infos) {
					#print_r($infos["QC_metrics"]); echo "<BR>";
					if ($infos["QC_metrics"]==$FASTQC_metrics_name) { #QC_metrics
						#print_r($infos["QCmodule"]); echo "<BR>"; 
						$QC_per_sequence_quality_scores_file_tmp="tmp/QC.$QC_id.per_sequence_quality_scores.".$infos["QC_metrics"];
						file_put_contents($QC_per_sequence_quality_scores_file_tmp,$infos["QCmodule"]["Per sequence quality scores"]);
						$per_sequence_quality_scores_file.="$QC_per_sequence_quality_scores_file_tmp,";
						$per_sequence_quality_scores_file_legends.="Per sequence quality scores".",";
					};#if
				};#foreach
				$QC_per_sequence_quality_scores="jpgraph_per_sequence_quality_scores.php?per_sequence_quality_scores_file=$per_sequence_quality_scores_file&legends=$per_sequence_quality_scores_file_legends";
				#echo $QC_per_sequence_quality_scores;

				# Per Base Quality link
				$per_base_quality_file="";
				$per_base_quality_file_legends="";
				foreach ($QC as $QC_id=>$infos) {
					if ($infos["QC_metrics"]==$FASTQC_metrics_name) { #QC_metrics
						#print_r($infos); echo "<BR>"; 
						$QC_per_base_quality_file_tmp="tmp/QC.$QC_id.per_base_quality.".$infos["QC_metrics"];
						file_put_contents($QC_per_base_quality_file_tmp,$infos["QCmodule"]["Per base sequence quality"]);
						$per_base_quality_file.="$QC_per_base_quality_file_tmp,";
						$per_base_quality_file_legends.="Per base quality".",";
					};#if
				};#foreach
				$QC_per_base_quality="jpgraph_per_base_quality.php?per_base_quality_file=$per_base_quality_file&legends=$per_base_quality_file_legends";

				# Coverage link
				$Files_sample_interval_summary_param="";
				$Files_sample_interval_summary_param_legends="";
				$depth_failed=array();
				$depth_warning=array();
				$depth_failed_manifest_region=array();
				foreach ($QC as $QC_id=>$infos) {
					if ($infos["QC_metrics"]=="sample_interval_summary") { #QC_metrics
						#print "<pre>"; print_r($infos); echo "</pre><BR>"; 
						#print "<pre>"; print_r($infos["QC_QC"]); echo "</pre><BR>"; 
						$QC_file_tmp="tmp/QC.$QC_id.".$infos["QC_metrics"];
						file_put_contents($QC_file_tmp,$infos["QC_QC"]);
						$Files_sample_interval_summary_param.="$QC_file_tmp,";
						$Files_sample_interval_summary_param_legends.=$infos['vcf_aligner'].",";

						# Test QC
						$depth_failed_min=30;
						$depth_warning_min=100;
						foreach( $infos["QCarray"] as $line => $split ) {
							#echo "<pre>".$infos['vcf_aligner']."<BR>"; print_r($split); echo "</pre>";
							if ($line>0 && trim($split[0])!="") {
								$chrompos=trim($split[0]);
								$split_chrompos = split('[:-]',$chrompos);
								$chrom=$split_chrompos[0];
								$start=$split_chrompos[1];
								$stop=$split_chrompos[2];
								
								$depth=trim($split[2]);
								
								if ($depth<max($depth_failed_min,$depth_warning_min)) {
									/*
									#Manifest region find
									$manifest_region=preg_grep("/$chrom.*$start.*$stop/i",split("\n",$manifest));
									foreach ($manifest_region as $manifest_region_k=>$manifest_region_line) {
										$manifest_region_split=split("\t",$manifest_region_line);
										$manifest_region_name=$manifest_region_split[0];
									};#foreach
									#print_r($manifest_region_name); print "<BR>";
									$depth_failed_manifest_region[$chrompos]=$manifest_region_name;
									*/
									if ($depth<$depth_failed_min) {
										$depth_failed[$chrompos][$QC_id]["depth"]=$depth;
									};#if
									if ($depth<$depth_warning_min) {
										$depth_warning[$chrompos][$QC_id]["depth"]=$depth;
									};#if
								};#if
							};#if
						};#foreach

					};#if
					
				};#foreach

				$option_region_translate="";
				foreach ($depth_failed_manifest_region as $chrompos=>$manifest_region_name) {
					#$option_region_translate.="$chrompos|$manifest_region_name,"; # DEBUG
				};#foreach
				#print "$option_region_translate<br>";
				
				$QCCoverage="jpgraph_coverage.php?coverage_file=$Files_sample_interval_summary_param&legends=$Files_sample_interval_summary_param_legends"; 
				if (join("",array_keys($depth_failed))!="") {
					$QCCoverage_failed=$QCCoverage."&title=Failed regions&region_translate=$option_region_translate&regions=".join(",",array_keys($depth_failed));
				};#if
				if (join("",array_keys($depth_warning))!="") {
					$QCCoverage_warning=$QCCoverage."&title=Warning regions coverage&region_translate=$option_region_translate&regions=".join(",",array_keys($depth_warning));
				};#if

				

				#$depth_failed=array();
				# Fail and Warning
				if ($depth_failed!=array()) {
					$warnings["Depth"][]="FAIL Depth (".count($depth_failed)." regions <span class=' ' onclick=\"javascript: var newWin = window.open('$QCCoverage_failed','QC_FASTQ','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=650, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>[+]</span>)";
					#foreach ($depth_failed as $chrompos=>$QC_ids) {
					#};#foreach
				};#if
				if ($depth_warning!=array()) {
					$warnings["Depth"][]="WARN Depth (".count($depth_warning)." regions <span class=' ' onclick=\"javascript: var newWin = window.open('$QCCoverage_warning','QC_FASTQ','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=650, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>[+]</span>)";
					#foreach ($depth_failed as $chrompos=>$QC_ids) {
					#};#foreach
				};#if

				
				
				#print $QCCoverage;
				
				# Show Metrics
				foreach ($QCMetrics_array as $BAMFile=>$values) {
					#echo "<PRE>";print_r($values); echo "</PRE><BR>";
					#$QCSummary_array["Reads Total"][$BAMFile]=$QCMetrics_array[$BAMFile]["TOTAL_READS"];
					#$QCSummary_array["Reads (PF|Aligned)"][$BAMFile]=$QCMetrics_array[$BAMFile]["TOTAL_READS"]." (".round($QCMetrics_array[$BAMFile]["PCT_PF_UQ_READS"]*100,2)."%"."|".round($QCMetrics_array[$BAMFile]["PCT_PF_UQ_READS_ALIGNED"]*100,2)."%)";
					$QCSummary_array["<span class=' ' onclick=\"javascript:var newWin = window.open('$QC_per_sequence_quality_scores','QC_FASTQ','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>[+]</span> Reads"][$BAMFile]=$QCMetrics_array[$BAMFile]["TOTAL_READS"];
					$QCSummary_array["&nbsp;&nbsp;PF/Aligned"][$BAMFile]="".round($QCMetrics_array[$BAMFile]["PCT_PF_UQ_READS"]*100,2)."%"."/".round($QCMetrics_array[$BAMFile]["PCT_PF_UQ_READS_ALIGNED"]*100,2)."%";
					# Coverage Warning
					if (($QCMetrics_array[$BAMFile]["PCT_PF_UQ_READS"]*100)<95) {
						$warnings["Reads"][]="WARN Passing Filter";
					};#if
					if (($QCMetrics_array[$BAMFile]["PCT_PF_UQ_READS"]*100)<90) {
						$warnings["Reads"][]="FAIL Passing Filter";
					};#if
					if (($QCMetrics_array[$BAMFile]["PCT_PF_UQ_READS_ALIGNED"]*100)<95) {
						$warnings["Reads"][]="WARN Alignment";
					};#if
					if (($QCMetrics_array[$BAMFile]["PCT_PF_UQ_READS_ALIGNED"]*100)<90) {
						$warnings["Reads"][]="FAIL Alignment";
					};#if


					#$QCSummary_array["Bases (On|Off)"][$BAMFile]=$QCMetrics_array[$BAMFile]["PF_UQ_BASES_ALIGNED"]." (".round($QCMetrics_array[$BAMFile]["PCT_SELECTED_BASES"]*100,2)."%"."|".round($QCMetrics_array[$BAMFile]["PCT_OFF_BAIT"]*100,2)."%)";
					$QCSummary_array["<span class=' ' onclick=\"javascript:var newWin = window.open('$QC_per_base_quality','QC_FASTQ','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>[+]</span> Bases"][$BAMFile]=$QCMetrics_array[$BAMFile]["PF_UQ_BASES_ALIGNED"];
					$QCSummary_array["&nbsp;&nbsp;On/Off"][$BAMFile]="".round($QCMetrics_array[$BAMFile]["PCT_SELECTED_BASES"]*100,2)."%"."/".round($QCMetrics_array[$BAMFile]["PCT_OFF_BAIT"]*100,2)."%";
					# Off Target Warning
					if (($QCMetrics_array[$BAMFile]["PCT_OFF_BAIT"]*100)>5) {
						$warnings["Coverage"][]="WARN Off Target";
					};#if
					if (($QCMetrics_array[$BAMFile]["PCT_OFF_BAIT"]*100)>10) {
						$warnings["Coverage"][]="FAIL Off Target very";
					};#if
					

					# Coverage
					$QCSummary_array["<span class=' ' onclick=\"javascript: var newWin = window.open('$QCCoverage','QC_FASTQ','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=650, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>[+]</span> Coverage"][$BAMFile]=round($QCMetrics_array[$BAMFile]["MEAN_TARGET_COVERAGE"],0);
					#$QCSummary_array["&nbsp;&nbsp;> 10X"][$BAMFile]=round($QCMetrics_array[$BAMFile]["PCT_TARGET_BASES_10X"]*100,2)."%";
					$QCSummary_array["&nbsp;&nbsp;> 30X"][$BAMFile]=round($QCMetrics_array[$BAMFile]["PCT_TARGET_BASES_30X"]*100,2)."%";
					#$QCSummary_array["&nbsp;&nbsp;> 50X"][$BAMFile]=round($QCMetrics_array[$BAMFile]["PCT_TARGET_BASES_50X"]*100,2)."%";
					$QCSummary_array["&nbsp;&nbsp;> 100X"][$BAMFile]=round($QCMetrics_array[$BAMFile]["PCT_TARGET_BASES_100X"]*100,2)."%";
					# Coverage Warning
					if (($QCMetrics_array[$BAMFile]["PCT_TARGET_BASES_100X"]*100)<95) {
						$warnings["Coverage"][]="WARN 100X";
					};#if
					if (($QCMetrics_array[$BAMFile]["PCT_TARGET_BASES_30X"]*100)<95) {
						$warnings["Coverage"][]="WARN 30X";
					};#if
					if (($QCMetrics_array[$BAMFile]["PCT_TARGET_BASES_100X"]*100)<90) {
						$warnings["Coverage"][]="FAIL 100X";
					};#if
					if (($QCMetrics_array[$BAMFile]["PCT_TARGET_BASES_30X"]*100)<90) {
						$warnings["Coverage"][]="FAIL 30X";
					};#if
					
				};#foreach


				$QCSummary="<TABLE border=0 style='border-color:white;'>";
				$QCSummary.="<TR><TD></TD><TD width='20px'>&nbsp;</TD>";
				foreach ($QCMetrics_array as $BAMFile=>$QCLink) {
					$QCFile_show=$BAMFile;
					$QCSummary.="<TD colspan=2 valign='top'><B>$QCFile_show</B></TD>";
				};#foreach
				$QCSummary.="</TR>";
				foreach ($QCSummary_array as $variable_name=>$variable_values) {
					$QCSummary.="<TR>";
					$QCSummary.="<TD valign='top'>$variable_name</TD><TD width='20px'>&nbsp;</TD>";
					foreach ($QCMetrics_array as $BAMFile=>$QCLink) {
						$QCSummary.="<TD>".$variable_values[$BAMFile]."</TD><TD width='20px'>&nbsp;</TD>";
					};#foreach
					$QCSummary.="</TR>";				
				};#foreach
			
				$QCSummary.="</TABLE>";
				
				

			};#if
			#print_r($QCList);
			
			# QC Files List
			$QCFiles_list_links="";
			$QCFiles_list_links.="";
			#foreach ($QCList_byExt as $extension=>$QCFiles_list) { #
			foreach ($QCList as $QCFile=>$QCLink) { #
				#foreach ($QCList as $QCFile=>$QCLink) { #
					$QCFiles_list_links.="<span class=' ' onclick=\"javascript:var newWin = window.open('".$QCLink["link"]."','QC_BAM','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes');newWin.focus()\" style='cursor:pointer;'>$QCFile</span><BR>";
				};#foreach
			#};#foreach
			$QCFiles_list_links.="";

			/*
			print "Depth failed ".count($depth_failed)."<pre>"; print_r($depth_failed); print "</pre><br>";
			print "Depth warning ".count($depth_warning)."<pre>"; print_r($depth_warning); print "</pre><br>";
			print "$QCCoverage_failed<BR>";
			print "$QCCoverage_warning<BR>";
			print "Warnings ".count($warnings)."<pre>"; print_r($warnings); print "</pre><br>";
			*/				

			# Warning
			$QCWarnings_line="";
			$QCWarnings_to_skip=array("Per base N content");
			foreach ($warnings as $class=>$warning_list) {
				$QCWarnings_warning_line="";
				foreach ($warning_list as $warning) {
					$warning_array=split(" ",$warning,2);
					if (!in_array($warning_array[1],$QCWarnings_to_skip)) {
						if (preg_grep("/FAIL/",array($warning))) {
							$warning_color="red";
						} elseif (preg_grep("/WARN/",array($warning))) {
							$warning_color="orange";
						} else {
							$warning_color="";
						};#if
						$QCWarnings_warning_line.="<span style='color:$warning_color'>$warning</span><BR>";
					};#if
				};#foreach
				if ($QCWarnings_warning_line != "") {
					$QCWarnings_line.="<TR class='' style='color:'><TD valign='top'><B>$class</B></TD><TD width='20px'>&nbsp;</TD><TD width='99%'><span class=''>$QCWarnings_warning_line</span></TD></TR>";
				};#if
				#$QCWarnings_line.=join("<BR> ",$warning_list);
				#$QCWarnings_line.="</span></TD></TR>";
			};#foreach
			if ($QCWarnings_line!="") {
				$QCWarnings="<TABLE width='100%'>$QCWarnings_line</TABLE>";
			} else {
				$QCWarnings="<TABLE width='100%'><TR><TD width='99%'>No QC Warning</TD></TR></TABLE>";
			};#if


			# QC files
			#$ExportFiles="".$variantList[1]." | ".$tmpfnamebase."<BR>";
			$QCFiles="

			<TABLE width='100%' border=0 class='' style='padding:10px'>
				<TR>
					<TD width='20px'>&nbsp;</TD>
					<TD colspan=1 height='40px;' valign='middle' style='min-width=400px'>
						<B>Quality Control</B>
					</TD>
					<TD width='40px'>&nbsp;&nbsp;&nbsp;&nbsp;</TD>
					<TD colspan=1 height='40px;' valign='middle' width='*'>
						<B>Warnings</B>
					</TD>
					<TD width='20px'>&nbsp;</TD>
				</TR>
				<TR>	
					<TD width='20px'>&nbsp;</TD>
					<TD colspan=1 valign='top'>
						$QCData
						$QCSummary
						<BR>
					</TD>
					<TD width='20px' valign='top'>&nbsp;</TD>
					<TD colspan=1 valign='top'>
						$QCWarnings
					</TD>
					<TD width='20px' valign='top'>&nbsp;</TD>
				</TR>
			</TABLE>
			";
			if ($_REQUEST["action_todo"]=="show_sample_report" && 1) {
			$QCFiles="

			<TABLE width='100%' border=0 class='' style='padding:10px'>
				<TR>
					<TD width='20px'>&nbsp;</TD>
					<TD colspan=1 height='40px;' valign='middle' width=''>
						<B>Quality Control</B>
					</TD>
				</TR>
				
				<TR>	
					<TD width='20px'>&nbsp;</TD>
					<TD colspan=1 valign='top'>
						$QCData
						$QCSummary
						<BR>
					</TD>
					
					<TD width='20px' valign='top'>&nbsp;</TD>
				</TR>
				<TR>
					<TD width='20px'>&nbsp;&nbsp;&nbsp;&nbsp;</TD>
					<TD colspan=1 height='40px;' valign='middle' width='99%'>
						<B>Warnings</B>
					</TD>
					<TD width='20px'>&nbsp;</TD>
				</TR><TR>	
					<TD width='20px' valign='top'>&nbsp;</TD>
					<TD colspan=1 valign='top' width=''>
						$QCWarnings
					</TD>
					<TD width='20px' valign='top'>&nbsp;</TD>
				</TR>
			</TABLE>
			";
			};#if
		
			#$QCFiles.="<DIV class=' ' onclick=\"javascript:var newWin = window.open('$QCFASTQLink','QC_FASTQ','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=1200, height=800, visible=yes');newWin.focus()\">UnAligned BAM/FASTQ</DIV>";
			


			/*
			$QCFiles.="<DIV class=' ' onclick=\"javascript:var newWin = window.open('$QCFASTQLink','QC_FASTQ','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=1200, height=800, visible=yes');newWin.focus()\">UnAligned BAM/FASTQ</DIV>";
			$QCfile_nb=0;
			foreach ($QCList as $QCFile=>$QCLink) {
				$QCfile_nb++;
				$QCFiles.="<DIV class=' ' onclick=\"javascript:var newWin = window.open('$QCLink','QC_BAM','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes');newWin.focus()\" style='cursor:pointer;'>$QCFile</DIV>";
			};#foreach
			foreach ($QCHsMetricsList as $QCFile=>$QCLink) {
				$QCfile_nb++;
				$QCFiles.="<DIV class=' ' onclick=\"javascript:var newWin = window.open('$QCLink','QC_BAM','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes');newWin.focus()\" style='cursor:pointer;'>$QCFile</DIV>";
			};#foreach
			*/
			#$QCFiles.="</TD></TR></TABLE>";
			
			

			if (1) { # Based on DB
				$query="SELECT vcf.id, vcf.ref
					FROM vcf
					INNER JOIN sample ON (sample.id=vcf.sample_id)
					INNER JOIN run ON (run.id=sample.run_id)
					WHERE run.ref='$input_run'
					  AND sample.ref='$input_sample'
					"; 
				#print "$query<BR>";
				$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			   	#$summary= "<table class=' '>\n"; 
				while ($row = mysql_fetch_assoc($result)) {
					#$VCFList[$row['ref']]=1;
					$VCFList[$row['id']]=$row['ref'];
					#$summary.= "<tr>"
					#	."<td>".$row['variant_file_id']."</td>"
					#	."<td>".$row['name']."</td>"
					#."</tr>"; 
				};#while
				#$summary.= "</table>\n"; 
				#echo $summary;
			};#if


			$vcf_files_select="";
			#print_r($vcf_files); print "S=".((count($vcf_files)==1 && $vcf_files[0]=="")?"selected":"")."<BR>";
			#print_r($vcf_files);
			#$vcf_files_select="<OPTION value='' ".((count($vcf_files)==1 && $vcf_files[0]=="")?"selected":"").">All VCF files</OPTION>";
			$vcf_files_select="<OPTION value='' ".((count($vcf_files)==1 && $vcf_files[0]=="")?"selected":"").">No pipeline</OPTION>";
			#foreach ($VCFList as $vcf_file=>$val) {
			foreach ($VCFList as $vcf_id=>$vcf_ref) {
				#$file_parts = pathinfo($vcf_file);
				#$vcf_files_select.="<OPTION value='$vcf_file' ".(in_array($vcf_file,$vcf_files)?"selected":"").">".$file_parts['filename']."</OPTION>\n";
				$vcf_files_select.="<OPTION value='$vcf_id' ".(in_array($vcf_id,$vcf_files)?"selected":"").">$vcf_ref</OPTION>\n";
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
			$RUNLink="run.php?action_todo=show_run&run=$input_run";
			$SAMPLELink="sample.php?action_todo=show_sample&run=$input_run&sample=$input_sample";
			$text.="<DIV class='brick run' onclick=\"javascript:document.location='$RUNLink'\">$input_run</DIV>";
			$text.="<DIV class='brick sample' onclick=\"javascript:document.location='$SAMPLELink'\" title='RUN: $input_run'>$input_sample</DIV>";
			$text.="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";
			$text.="";


			if (0) { # comment old with file ini
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

			};#if

			$comment_header="Comment";
			$comment_current=$comment;

			## VALIDATION

			if (1) { # Based on DB

				# variant flag
				$flag="0";
				$comment="";
				$variant_id=$variant_array["VARIANT_ID"];
				$query="SELECT VS.flag_id, VS.comment, VS.history
					FROM comment AS VS
					WHERE VS.id='$sample_id'
					  AND VS.type='sample'
					 # AND VS.project_id='$project_id'
					"; 
				#$query="SELECT sample_id FROM sample LIMIT 5";
				#print "$query<BR>";
				$result = mysql_query($query,$mydbh) or die ("OUPS! ".mysql_error($mydbh)." <BR><PRE>$query</PRE>");
				$summary= "<table class=' '>\n"; 
				while ($row = mysql_fetch_assoc($result)) {
					$comment_flag_id=$row['flag_id'];
					$comment=$row['comment'];
					$summary.= "<tr>"
						."<td>".$row['flag']."</td>"
						."<td>".$row['comment']."</td>"
					."</tr>"; 
				};#while
				$summary.= "</table>\n"; 
				#echo "<PRE>$query<BR>$summary</PRE>";
				
				$query="SELECT id, type, rank, score, label, code, definition
					FROM flags
					WHERE type='sample'
					ORDER BY abs(rank) ASC
					"; 
				#print "$query<BR>";
				$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			   	$summary= ""; 
				$flags_options="";
				while ($row = mysql_fetch_assoc($result)) {
					$flag_id=$row['id'];
					$flag_value=$row['flag'];
					#$flag_description=$row['label']." [".$row['level']."]";
					$flag_label=$row['label']; #." [".$row['level']."]";
					$summary.= "<tr>";
					$summary_head="<tr>";
					foreach ($row as $col=>$val) {
						$summary.= "<td>$val</td>";
						$summary_head.= "<td>$col</td>";
					};#foreach
					$summary.= "</tr>"; 
					$summary_head.= "</tr>"; 
					#$flags_array[$row['flag_value']]=$row['flag_label'];
					$flags_options.="<OPTION value='$flag_id' ".(($flag_id==$comment_flag_id)?"selected":"").">$flag_label</OPTION>\n";
					if ($flag_value==$flag) {
						$flag_description_show=$flag_description;
					};#if
				};#while
				$summary= "<table class=' '>$summary_head $summary</table>\n"; 
				#echo $summary."<textarea>$flags_options</textarea>";
			};#if

			$URI=$_SERVER["REQUEST_URI"];
			#$ReportLink="sample.php?action_todo=show_sample_report&run=$input_run&sample=$input_sample";
			$ReportLink="'$URI&action_todo=show_sample_report&rand='+Math.floor((Math.random()*10000000)+1)";
		
			$validation="
				<DIV class=''>
				<FORM  action='validation.php' name='validation' target='validation' METHOD='POST'>
				<iframe name='validation' style='visibility:hidden;display:none;'></iframe>
				<input type='hidden' id='type' name='type' style='visibility:hidden;float:right;display:none;' value='sample'></input>
				<input type='hidden' id='project_id' name='project_id' style='visibility:hidden;float:right;display:none;' value='$project_id'></input>
				<!--<input type='hidden' id='run' name='run' style='visibility:hidden;float:right;display:none;' value='$input_run'></input>-->
				<!--<input type='hidden' id='sample' name='sample' style='visibility:hidden;float:right;display:none;' value='$input_sample'></input>-->
				<!--<input type='hidden' id='sample_id' name='sample_id' style='visibility:hidden;float:right;display:none;' value='$sample_id'></input>-->
				<input type='hidden' id='id' name='id' style='visibility:hidden;float:right;display:none;' value='$sample_id'></input>
				<TABLE width='100%' border=0 class='' style='padding:10px'>
					<TR>
						<TD width='20px'>&nbsp;</TD>
						<TD colspan=3 width='300px' height='40px' valign='middle'  class=''> <!-- commentsample_header -->
							<B>Validation</B>
						</TD>
						<TD width='20px'>&nbsp;</TD>
					</TR>
					<TR>
						<TD width='20px'>&nbsp;</TD>
						<TD width='99%' valign='top'>
							<SELECT name='flag_id' onchange=\"javascript:;\" style='width:100%;'>
								$flags_options
							</SELECT>
						</TD>
						<TD valign='top'>
							<input type='submit' id='validation' name='validation' value='Save' class='btn button search' onclick=\"javascript\">
						</TD>
						<TD width='20px'>&nbsp;</TD>
					</TR>
					<TR>
						<TD width='20px'>&nbsp;</TD>
						<TD colspan=3 align='center' valign='top'>
							<TEXTAREA id='comment' name='comment' class='' style='border: none; width: 99%; height: 100%; -webkit-box-sizing: border-box; -moz-box-sizing: border-box; box-sizing: border-box; padding:10px;' rows=6>$comment</TEXTAREA>
						</TD>
						<TD width='20px'>&nbsp;</TD>
					</TR>
					<TR>
						<TD width='20px'>&nbsp;</TD>
						<TD colspan=3 align='center' valign='top'>
							<span onclick=\"javascript:var newWin = window.open($ReportLink,'Report','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=800, visible=yes, resizable=yes');newWin.focus()\" style='cursor:pointer;'>Report</span><!-- alert('Still under development!!!');  -->
						</TD>
						<TD width='20px'>&nbsp;</TD>
					</TR>
				</TABLE>
				</FORM>
				</DIV>
			";

			if ($_REQUEST["action_todo"]=="show_sample_report") {
				$validation="
				<DIV class=''>
				<TABLE width='100%' border=0 style='padding:10px' >
					<TR>
						<TD width='20px'>&nbsp;</TD>
						<TD colspan=3 width='300px' height='40px' valign='middle'>
							<B>Validation</B>
						</TD>
						<TD width='20px'>&nbsp;</TD>
					</TR>
					<TR>
						<TD width='20px'>&nbsp;</TD>
						<TD width='99%' valign='top' colspan=3>
							<B>$flag_description_show</B>
						</TD>
						<TD width='20px'>&nbsp;</TD>
					</TR>
					<TR>
						<TD width='20px'>&nbsp;</TD>
						<TD colspan=3 align='' valign='top'>
							".str_replace("\n","<BR>",$comment)."
						</TD>
						<TD width='20px'>&nbsp;</TD>
					</TR>
				</TABLE>
				</FORM>
				</DIV>
			";
			};#if
			#<p><pre class='' style='border:0;background-color:transparent;'>$comment</pre></p>


			$text.="
				<TABLE width='95%' >
					<TR>
						<TD valign='top' class='validation' style='width:250px;min-width:250px;max-width:400px;'>
							$validation
						</TD>
						<TD width='20px'>&nbsp;</TD>
						<TD valign='top' class='validation' width='*'>
							$QCFiles
						</TD>
					</TR>
				</TABLE>

				";
			
			$flags_target_frame="flags";
			$text.="<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>";


			# Form filter display
			$select_width="150"; # auto, 150...
			if ($_REQUEST["action_todo"]=="show_sample") {
				$filter_form="
				<TABLE width='95%' border=0 class='filter'>
					<TR><TD>
				
				<div style='text-align: left; table-border:1; clear:left;' class='' >
				<form action='?' method='get' class='' name='main' id='main'>
				<INPUT type='hidden' name='action_todo' value='$action_todo'>
				<INPUT type='hidden' name='sample' value='$input_sample'>
				<INPUT type='hidden' name='run' value='$input_run'>
				<INPUT type='hidden' name='q' value='$q'>
				

				<TABLE  border=0  width='100%' >
					<TR>
						<TD width='20px'>&nbsp;</TD>
						<TD colspan=6 height='40px;' valign='middle'>
							<B>Filters</B>
						</TD>
						<TD width='20px'>&nbsp;</TD>
					</TR>
					<TR>
						<TD width='20px;'>&nbsp;</TD>
						<TD valign='top' align='left' width='200px'>
							Workflows
						</TD>
						<TD width='20px;'>&nbsp;</TD>
						<TD valign='top' halign='center' width='200px'>
							Filters
						</TD>
						<!--<TD valign='top' halign='center' width='200px'>
							Annotations
						</TD>-->
						<TD width='20px;'>&nbsp;</TD>
						<TD valign='top' halign='center' width=''>
							Display
						</TD>
						<TD width='20px'>&nbsp;</TD>
					</TR>
					<TR>
						<TD width='20px;'>&nbsp;</TD>
						<TD valign='top' halign='center'>
							<SELECT name='vcf_files[]' multiple size=".(count($VCFList)+1>5?5:count($VCFList)+1)." style='width:$select_width;'> <!--  multiple size=1 -->
								$vcf_files_select
							</SELECT>
						</TD>
						<TD width='20px;'>&nbsp;</TD>
						<TD valign='top' halign='center'>
							<TABLE width=''>
								<TR><TD colspan=2>
									<SELECT name='filters[]' multiple size=".(count($filter_array)>5?5:count($filter_array))." style='width:$select_width;'>
										$filter_select
									</SELECT>
								</TD></TR>
								<TR>
									<TD halign=''>
										<INPUT type='checkbox' name='hardfiltering' ".( $hardfiltering ? "checked" : "")." title=\" (hide 'FILTERED' variants)\"> Hard Filtering</INPUT>
									</TD>
									<TD>
										[<span onclick=\"javascript:var newWin = window.open('filters.php','popup','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes');document.getElementById('main').action='filters.php';document.getElementById('main').target='popup';document.getElementById('main').submit();newWin.focus();document.getElementById('main').action='$action_default';document.getElementById('main').target='_self';return false;\" style='cursor:pointer;'>details</span>]
									</TD>
								</TR>
							</TABLE>

							
						</TD>
						<!--<TD valign='top' halign='center'>
							<SELECT name='annotation_list[]' multiple size=".(count($annotation_list)>3?5:count($annotation_array)-2)." style='width:$select_width;'>	
								<OPTION value='NO' ".(in_array("NO",$annotation_list)?"selected":"").">No more annotations 
								<OPTION value='ALL' ".(in_array("ALL",$annotation_list)?"selected":"").">ALL annotations	
								$annotation_list_select	
							</SELECT> 
						</TD>-->
						<TD width='20px;'>&nbsp;</TD>
						<TD valign='top' halign='center'>


							<TABLE width=''>
								<TR>
									<TD halign=''>
										Filter
									</TD>
									<TD width='20'>
										
									</TD>
									<TD>
										<INPUT type='text' name='global_filter' id='global_filter' value='$global_filter' style='width:120px;'>
									</TD>
								</TR>
								<TR>
									<TD halign=''>
										Gene/HGVS
									</TD>
									<TD width='20'>
										
									</TD>
									<TD>
										<INPUT type='text' name='gene_filter' id='gene_filter' value='$gene_filter' style='width:80px;'>
									</TD>
								</TR>
								<TR>
									<TD halign=''>
										Order By
									</TD>
									<TD width='20'>
										
									</TD>
									<TD>
										<SELECT name='orderby' style='width:120px;'>
											<!--$annotation_list_for_order_by_select-->
											<OPTION value='PZScore' ".($orderby=="PZScore"?"selected":"").">PZScore
											<OPTION value='VDScore' ".($orderby=="VDScore"?"selected":"").">VDScore
											<OPTION value='Symbol' ".($orderby=="Symbol"?"selected":"").">Symbol
										</SELECT> 
										<SELECT name='ascdesc' style='width:auto;'>
											<OPTION value='DESC' ".($ascdesc=="DESC"?"selected":"").">DESC
											<OPTION value='ASC' ".($ascdesc=="ASC"?"selected":"").">ASC		
										</SELECT>
									</TD>
								</TR>
								<TR>
									<TD halign=''>
										Limit
									</TD>
									<TD width='20'>
										
									</TD>
									<TD>
										<INPUT type='text' name='limit' id='limit' value='$limit' style='width:30px;'>
									</TD>
								</TR>
								
							</TABLE>

						</TD>
						<TD width='20px;'>&nbsp;</TD>
					</TR>
					<TR>
						<TD width='20px;'>&nbsp;</TD>
						<TD>
							<INPUT TYPE='submit' value='Process' name='process' class='btn button search' data-original-title=''/>
						</TD>
						<TD width='20px;'>&nbsp;</TD>
					</TR>
				</TABLE>
				</div>
				</form>
				
				</TD>
					</TR>
				</TABLE>
				
				
				";
			};#if

			if ($_REQUEST["action_todo"]=="show_sample") {
				$variant_list="
					<div style='text-align: left; table-border:1; ' class=''>
					<!--$ExportFiles <DIV style='margin-bottom:0px; clear:left;'><BR></DIV>-->
					$Variants_HTMLContent
					</div>
				";
			} elseif ($_REQUEST["action_todo"]=="show_sample_report") {
				$variant_list="
					<div style='text-align: left; table-border:1; ' class=''>
					$Variants_HTMLContent
					</div>
				";
			};#if
				$text.="
					$filter_form
					<BR>
					$variant_list

				";
				#$variants_output
			break;
		}#case show_sample

	};#switch

};#if

	#$text.= "</div></div>";


$title=" ";

$ns->tablerender($title, $text, "NGS");



require_once(FOOTERF);

?>
