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
$date="02/01/2015";
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

#
# VCF File
# 

$input_runs=$_REQUEST["runs"];
#if ($input_runs=="") {
#	$input_runs[]="";
#};#if


#
# Aligners
# 

$aligners=$_REQUEST["aligners"];
if ($aligners=="") {
	$aligners[]="";
};#if

#
# Callers
# 

$callers=$_REQUEST["callers"];
if ($callers=="") {
	$callers[]="";
};#if

#
# Annotators
# 

$annotators=$_REQUEST["annotators"];
if ($annotators=="") {
	$annotators[]="";
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



# Security
if (!ADMIN) {
	$text="<div style='text-align:center'>Please login</div>";
} else {

	##
	## Style
	##

	
	$text="";
	$text.="<link rel='stylesheet' media='all' type='text/css' href='style.css' />";

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

	# MENU
	########

	# List RUN available for Analysis
	# Directory exists
	
	$select_run="<select name='runs[]' id='runs' style='width:350px;' multiple size=5> ";
	$select_run.="<option value=''>No RUN selected</option> ";
	$Directory=$dir_miseq;

	$dir_list=preg_grep('/^([^.])/',scandir($Directory));
	$run_already_analyzed_text=""; #"[Already Analyzed]";
	$run_NOTalready_analyzed_text="[NOT Analyzed]";
	$nb_option=0;
	$select_run_show_hide="";
	foreach($dir_list as $Entry) {
		if(is_dir($Directory.'/'.$Entry) && $Entry != '.' && $Entry != '..') {
			#echo "$dir_analysis/$Entry<BR>";
			$nb_option++;
			if (glob("$dir_analysis/$Entry/*.release")) { #is_dir($dir_analysis/$Entry)
				$run_analyzed=$run_already_analyzed_text;
				$select_run_show_hide.="window.document.getElementById('runs').options[$nb_option].hidden=!window.document.getElementById('runs').options[$nb_option].hidden;";
			} else {
				$run_analyzed=$run_NOTalready_analyzed_text;				
			};#if
            		$select_run.= "<option value='$Entry' ".((in_array("$Entry",$input_runs))?"selected":"")." id='run_$Entry'>$Entry $run_analyzed</option>"; # 
			
		};#if
	};#foreach

	$select_run.="</select> ";
	$select_run.="
		<SCRIPT>
		function show_hide_options() {
			$select_run_show_hide
		}
		</SCRIPT>
		";
	
	# List Aligners
	$aligners_available=array(
		""=>"Default (bwamem bwamemUnclipped)",
		"bwamem"=>"BWA MEM",
		"bwamemUnclipped"=>"BWA MEM Unclipped",
		"bwasw"=>"BWA SW (Smith-Waterman)",
		"bwaswUnclipped"=>"BWA SW (Smith-Waterman) Unclipped",
		"bwaaln"=>"BWA ALN",
		"bwaalnUnclipped"=>"BWA ALN Unclipped",
		);
	$size_aligners=(count($aligners_available)<5)?count($aligners_available):"5";
	$select_aligners="<select name='aligners[]' id='aligners[]' multiple style='width:350px;' size=$size_aligners> ";
	foreach ($aligners_available as $aligner_id=>$aligner_label) {
		$select_aligners.="<option value='$aligner_id' ".((in_array("$aligner_id",$aligners))?"selected":"").">$aligner_label</option> ";
	};#foreach
	$select_aligners.="</select> ";
	
	# List Callers
	$callers_available=array(
		""=>"Default (gatkHC gatkUG)",
		"gatkHC"=>"GATK Haplotype Caller",
		"gatkUG"=>"GATK Unifeid Genotyper",
		);
	$size_callers=(count($callers_available)<5)?count($callers_available):"5";
	$select_callers="<select name='callers[]' id='callers[]' multiple style='width:350px;' size=$size_callers> ";
	foreach ($callers_available as $caller_id=>$caller_label) {
		$select_callers.="<option value='$caller_id' ".((in_array("$caller_id",$callers))?"selected":"").">$caller_label</option> ";
	};#foreach
	$select_callers.="</select> ";
	
	# List Annotators
	$annotators_available=array(
		""=>"Default (VAP)",
		"vap"=>"VAP",
		);
	$size_annotators=(count($annotators_available)<5)?count($annotators_available):"5";
	$select_annotators="<select name='annotators[]' id='annotators[]' multiple style='width:350px;' size=$size_annotators> ";
	foreach ($annotators_available as $annotator_id=>$annotator_label) {
		$select_annotators.="<option value='$annotator_id' ".((in_array("$annotator_id",$annotators))?"selected":"").">$annotator_label</option> ";
	};#foreach
	$select_annotators.="</select> ";

	
	$command="ps -A | grep ' launch.sh$' | awk -F' ' '{print $1}'";
	$command_ID = shell_exec($command);
	$launched=0;
	if ($command_ID!="") {
		$launched=$command_ID;
	};#if;


	if ($launched<1) {
	$text.= "<TABLE>
		<TR>
			<TD valign='top'>
				<p class=' title'>RUN Analysis process</p>
				<BR>
				<FORM action='?' method='get' class='' name='main' id='main'>
					<INPUT type='hidden' name='action_todo' name='action_todo' value='analysis'>
					<p class=''>Please check if analysis isn't already launched</p><p class='warning'>".(($launched>0)?"RUN Analysis process already launched: $command_ID<BR>":"")."</p><BR>
					RUN: [<span onclick='show_hide_options();' style='cursor:pointer;'>Show/hide already analyzed RUNs</span>]<BR>$select_run<BR>
					ALIGNERS:<BR>$select_aligners<BR>
					CALLERS:<BR>$select_callers<BR>
					ANNOTATORS:<BR>$select_annotators<BR>
					<BR><input type='submit' value='Start Analysis' name='submit' id='submit' class='btn button search' data-original-title=''/>
				</FORM>

				
				<BR>
			</TD>
			<TD width='50'><BR></TD>
		</TR>
	</TABLE>";
	};#if


	# Action
	##########

	$output_analysis="";
	if ($launched>0) {
		$output_analysis.="RUN Analysis process already launched: $launched<BR>Please wait the end of the process";
	};#if


	switch ($_REQUEST["action_todo"]) {

		
		# Analysis
		case "analysis":
		{

			if ($launched>0) {
				$output_analysis.="[ERROR] RUN Analysis process already launched: $launched<BR>Please wait the end of the process";
			} elseif (trim(join("",$input_runs))=="") {
				$output_analysis.="[ERROR] No RUN selected!";
			} else {
				# Command
				$command="$dir_bin/scripts/launch.sh \"".join(" ",$input_runs)."\" \"".join(" ",$aligners)."\" \"".join(" ",$callers)."\" \"".join(" ",$annotators)."\" > /dev/null &";
				$output_exec = shell_exec($command);
				$output_analysis.= "<TABLE>
					<TR>
						<TD valign='top'>
							
							<BR>
						
								RUN: ".((join(",",$input_runs)!="")?join(",",$input_runs):"[ERROR]")."<BR>
								ALIGNERS: ".((join(",",$aligners)!="")?join(",",$aligners):$aligners_available[""])."<BR>
								CALLERS: ".((join(",",$callers)!="")?join(",",$callers):$callers_available[""])."<BR>
								ANNOTATORS: ".((join(",",$annotators)!="")?join(",",$annotators):$annotators_available[""])."<BR>
							
								COMMAND: $command<BR>
			
								OUTPUT: <PRE>$output_exec</PRE>

							<BR>
						</TD>
						<TD width='50'><BR></TD>
					</TR>
				</TABLE>";
			};#if
			
			break;
		};#case analysis

	};#switch


	$command="$dir_bin/scripts/monitor.sh 1 details";
	$output_exec = shell_exec($command);
	$output_analysis.="<PRE>$output_exec</PRE>";

	$text.=$output_analysis;

};#if

	$text.= "
			</div>
		</div>
	";


$title=" ";

$ns->tablerender($title, $text, "TRAKXS");



require_once(FOOTERF);

?>
