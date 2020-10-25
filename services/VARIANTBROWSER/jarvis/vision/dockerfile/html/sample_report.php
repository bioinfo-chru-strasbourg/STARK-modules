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
$date="20140624";
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
$text.="<link rel='stylesheet' media='all' type='text/css' href='style.css' />
<STYLE>
.navbar {
   display:none;
}
</STYLE>

";

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
#include "search.header.php";


	switch ($_REQUEST["action_todo"]) {


		# Show a SAMPLE
		case "show_sample_report":
		{

			# Sample info
			if (1) { # Based on DB
			$query="SELECT sample.sample_id, sample.ref, sample.flag, sample.comment
				FROM sample
				INNER JOIN project ON (project.project_id=project.project_id)
				WHERE sample.run='$input_run'
				  AND sample.ref='$input_sample'
				  AND project.ref='$project'
				"; 
			#print "$query<BR>";
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
		   	$summary= "<table class=' '>\n"; 
			while ($row = mysql_fetch_assoc($result)) {
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

			# Sample info
			if (1) { # Based on DB
			$variant_ids="";
			$query="SELECT VS.variant_id, VS.sample_id, flags.label, flags.definition
				FROM variant_sample AS VS
				INNER JOIN project ON (project.project_id=project.project_id)
				INNER JOIN sample ON (sample.sample_id=VS.sample_id)
				INNER JOIN variant ON (variant.variant_id=VS.variant_id)
				INNER JOIN validation_variant_sample ON (validation_variant_sample.variant_id=VS.variant_id
				   AND validation_variant_sample.sample_id=VS.sample_id)
				INNER JOIN flags ON (flags.level=validation_variant_sample.flag)
				WHERE sample.run='$input_run'
				  AND sample.ref='$input_sample'
				  AND project.ref='$project'
				  AND flags.type='variant_sample'
				  #AND flags.label='validated'
				  AND flags.level>0
				GROUP BY VS.variant_id 
				ORDER BY flags.level ASC
				"; 

			#print "$query<BR>";
			$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
		   	$summary= "<table class=' '>\n"; 
			while ($row = mysql_fetch_assoc($result)) {
				#print_r($row);
				$variant_id=$row['variant_id'];
				$sample_id=$row['sample_id'];
				$flag_label=$row['flag_label'];
				$flag_comment=$row['flag_comment'];
				#$comment=$row['comment'];
				$summary.= "<tr>"
					."<td>$variant_id</td>"
					."<td>$sample_id</td>"
					."<td>$flag_label</td>"
					."<td>$flag_comment</td>"
				."</tr>"; 
				$variant_ids.="$variant_id,";
			};#while
			$summary.= "</table>\n"; 
			#echo $summary;
			};#if


			
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


			# Variant List
			$tmpfnamebase = tempnam("tmp/", "variant_file_");
			
			# Exportation
			$variant_files=implode(",",$vcf_files);
			$variant_files_option="";
			if ($variant_files!="") {
				$variant_files_option="--variant_files=$variant_files";
			};#if
			##if ($variant_ids!="") {
			#	$variant_ids_option="--variant_ids=$variant_ids";
			#};#if
			$command="$dir_bin/scripts/export.pl --run=$input_run --samples=$input_sample $variant_files_option --variant_ids=$variant_ids --vcf=$tmpfnamebase --annotation_list='ALL' ";
			#print "$command<BR>";
			$output_exec = shell_exec($command);

			# Generate Files
			

			$variantList=ShowSAMPLE_VariantList($dir_analysis,$input_run,$input_sample,$aligners[0],$callers[0],$annotators[0],$tmpfnamebase);
			$variants_output=$variantList[0];
			#$variants_output=ShowSAMPLE_VariantList($dir_analysis,$input_run,$input_sample,$aligners[0],$callers[0],$annotators[0],$vcf_files[0]);
			

			# VCF to HTML
			
			$tmpfnamebase=$variantList[1];
			$VCFFile="$tmpfnamebase.vcf";
			#print $VCFFile;
			$Variants_VCFContent=file_get_contents($VCFFile);
			#$f = file($VCFFile, 'r');
			#print_r($f); print "<BR>";
			#print "<pre>$Variants_VCFContent</pre>";
			
			$Variants_HTMLContent=VCFtoHTML($Variants_VCFContent,$sample_id,$annotation_list);
			

			# Export files
			#$ExportFiles="".$variantList[1]." | ".$tmpfnamebase."<BR>";
			$ExportFiles.="<DIV class='brick file' onclick=\"javascript:var newWin = window.open('tmp/".basename($tmpfnamebase).".vcf','EXPORT_VCF','menubar=yes, status=no, scrollbars=no, menubar=yes, width=800, height=600, visible=yes');newWin.focus()\">VCF</DIV> ";
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
			if (1) { # based on files
				$tmpfnamebase_QC = tempnam("tmp/", "QC_".$input_run."_".$input_sample."_");
				$Directory="$dir_analysis/$input_run/$input_sample";
				# FASTQ
				$QC_FAST_src="$Directory/$input_sample.fastqc";
				$QC_FAST_dst=$tmpfnamebase_QC."_DIR";
				@mkdir($QC_FAST_dest); 
				#print "$QC_FAST_src => $QC_FAST_dest<BR>\n";
				recurse_copy($QC_FAST_src,$QC_FAST_dst);
				$QCFASTQLink="tmp/".basename($QC_FAST_dst)."/".$input_sample.".unaligned_fastqc/fastqc_report.html";
				#P1261.unaligned_fastqc
				#if (copy($Directory.'/'.$Entry.'/'.$Entry2,$tmpfnamebase_QC.'_'.$Entry.'_'.$Entry2)) {
				#	#print "copy $Directory/$Entry/$Entry2 => $tmpfnamebase_QC.'_'.$Entry.'_'.$Entry2 ok<BR>\n";
				#	$QCFASTQFile=$tmpfnamebase_QC.".unaligned_fastqc";
				#	$QCFASTQLink=$tmpfnamebase_QC.".unaligned_fastqc";
				#};#if
				# BAM
				$QCList=array();
				$MyDirectory = opendir($Directory); # or die('Erreur');
				while($Entry = @readdir($MyDirectory)) {
					if(is_dir($Directory.'/'.$Entry) && $Entry != "." && $Entry != "..") {
						#print $Directory.'/'.$Entry."<BR>\n";
						$MyDirectory2 = opendir($Directory.'/'.$Entry); # or die('Erreur');
						while($Entry2 = @readdir($MyDirectory2)) {
							#print $Directory.'/'.$Entry.'/'.$Entry2."<BR>\n";
							if(is_file($Directory.'/'.$Entry.'/'.$Entry2)) {
								$QCFile=$Entry.'/'.$Entry2;
								$file_parts = pathinfo($QCFile);
								if (strtolower($file_parts['extension'])=="sample_interval_summary") {
									$QCFile_tmp=$tmpfnamebase_QC.'_'.$Entry.'_'.$Entry2;
									#if (copy($Directory.'/'.$Entry.'/'.$Entry2,$tmpfnamebase_QC.'_'.$Entry.'_'.$Entry2)) {
										#print "copy $Directory/$Entry/$Entry2 => $tmpfnamebase_QC.'_'.$Entry.'_'.$Entry2 ok<BR>\n";
										
									#};#if
									$QCList[$Entry2]="tmp/".basename($QCFile_tmp);
									#print $QCFile;
								};#if
							};#if
						};#while
						closedir($MyDirectory2);	


						


					};#if
				}
				closedir($MyDirectory);	
			};#if
			#print_r($QCList);
			

			# QC files
			#$ExportFiles="".$variantList[1]." | ".$tmpfnamebase."<BR>";
			$QCFiles="<DIV class='brick file' onclick=\"javascript:alert('Click on a file to open it!');\">Quality Control</DIV>";
			$QCFiles.="<DIV class='brick file' onclick=\"javascript:var newWin = window.open('$QCFASTQLink','QC_FASTQ','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=1200, height=800, visible=yes');newWin.focus()\">UnAligned BAM/FASTQ</DIV>";
			$QCfile_nb=0;
			foreach ($QCList as $QCFile=>$QCLink) {
				$QCfile_nb++;
				$QCFiles.="<DIV class='brick file' onclick=\"javascript:var newWin = window.open('$QCLink','QC_BAM','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes');newWin.focus()\">$QCFile</DIV> ";
			};#foreach

			
			

			if (1) { # Based on DB
				$query="SELECT variant_file_id, name
					FROM variant_file
					INNER JOIN sample ON (sample.sample_id=variant_file.sample_id)
					WHERE run='$input_run'
					  AND ref='$input_sample'
					"; 
				$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			   	#$summary= "<table class=' '>\n"; 
				while ($row = mysql_fetch_assoc($result)) {
					$VCFList[$row['name']]=1;
					#$summary.= "<tr>"
					#	."<td>".$row['variant_file_id']."</td>"
					#	."<td>".$row['name']."</td>"
					#."</tr>"; 
				};#while
				#$summary.= "</table>\n"; 
				#echo $summary;
			};#if


			$vcf_files_select="";
			#print_r($vcf_files);
			$vcf_files_select="<OPTION value='' ".((count($vcf_files)==1 && $vcf_files[0]=="")?"selected":"").">All VCF files</OPTION>";
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

			##FLAG

			if (1) { # Based on DB
				$query="SELECT level, label, type, definition
					FROM flags
					WHERE type='sample'
					"; 
				#print "$query<BR>";
				$result = mysql_query($query,$mydbh) or die (mysql_error($mydbh));
			   	$summary= ""; 
				$flags_options="";
				while ($row = mysql_fetch_assoc($result)) {
					$summary.= "<tr>";
					$summary_head="<tr>";
					foreach ($row as $col=>$val) {
						$summary.= "<td>$val</td>";
						$summary_head.= "<td>$col</td>";
					};#foreach
					$summary.= "</tr>"; 
					$summary_head.= "</tr>"; 
					$flags_arra[$row['flag_value']]=$row['flag_label'];
					$flags_options.="<OPTION value='".$row['flag_value']."' ".(($row['flag_value']==$flag)?"selected":"").">".$row['flag_label']."</OPTION>\n";
					if ($row['flag_value']==$flag) {
						$flag_label=$row['flag_label'];
					};#if
					#$flag_label=(($row['flag_value']==$flag)?$row['flag_label']:"");
				};#while
				$summary= "<table class=' '>$summary_head $summary</table>\n"; 
				#echo $summary;
			};#if

			$flags_target_frame="flags";

			# Flag and Comment
			$text.="<DIV class='commentsample' onclick=\"javascript:\"> 
					<pre class='prestyle'>Sample '$flag_label'<BR>$comment_current</pre>
				</DIV>";
		
			# QC
			$QC_infos="TODO";
			$text.="<DIV class='commentsample' onclick=\"javascript:\"> 
					<pre class='prestyle'>QC information<BR>$QC_infos</pre>
				</DIV>
				<DIV style='margin-bottom:0px; clear:left;'><BR></DIV>
				";

			# Form filter display
			$text.="
				
				
				
				<div style='text-align: left; table-border:1; ' class=''>
				
				$Variants_HTMLContent
				</div>
				
";
				#$variants_output
			break;
		}#case show_sample

	};#switch

};#if



$title="

";

$ns->tablerender($title, $text, "TRAKXS");



#require_once(FOOTERF);

?>
