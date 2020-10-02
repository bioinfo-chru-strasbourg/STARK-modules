<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">

<!--
<link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/handsontable@latest/dist/handsontable.full.min.css">
<link rel="stylesheet" type="text/css" href="https://handsontable.com/static/css/main.css">
<script src="https://cdn.jsdelivr.net/npm/handsontable@latest/dist/handsontable.full.min.js"></script>
<link rel='stylesheet' media='all' type='text/css' href='style.css' />
-->
<link rel="stylesheet" type="text/css" href="handsontable/dist/handsontable.full.min.css">
<link rel="stylesheet" type="text/css" href="handsontable/static/css/main.css">


<script src="handsontable/dist/handsontable.full.min.js"></script>


   
<?php

if ( ! session_id() ) @ session_start();


#
# Versioning
##############

$release="1.0.1";
$date="20191203";
$authors="Antony Le Béchec";

#
# HEADER
##########





#if(!defined('e107_INIT'))
#{
#   require_once('../../class2.php');
#}
#$e107 = e107::getInstance();
#$tp = e107::getParser();
#$sql = e107::getDb();

#require_once(HEADERF);
require("functions.inc.php");


//load the function
require_once 'tsv-to-array.inc.php';


#error_reporting(0);
#
# PARAM
#########

#$dir_analysis="/media/NGS/analysis";
#$dir_miseq="/media/miseq/MSR";
$dir_analysis="/media/IRC/RES/ALL";
$dir_analysisV1="/media/IRC/RES/ALL";
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
$process=$_REQUEST["process"];

if ($process=="") {
    $process=1;
}

#
# VCF File
#

$vcf_files=$_REQUEST["vcf_files"];
$vcf_filenames=$_REQUEST["vcf_files"];
#if ($vcf_files=="") {
#   $vcf_files[]="";
#};#if


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
    #$filters[]="default";
    $filters[]="";
};#if
$hardfiltering=$_REQUEST["hardfiltering"];

#
# Order by
#

$orderby=$_REQUEST["orderby"];
if ($orderby=="") {
    $orderby="";
} else {
    $orderby_split=explode(":",$orderby);
    $orderby=$orderby_split[0];
    $orderby_type=$orderby_split[1];
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
# URL file
#

$gene_filter=$_REQUEST["vcf"];
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
    $annotation_list=array_merge($annotation_list,explode(",",$annotation_list_default));
};#if

#$annotation_list=array_merge($annotation_list,$annotation_list_mandatory);
$annotation_list=array_merge($annotation_list_mandatory,$annotation_list);

#
# Security
#
if (!USER && 0) {
    $text="<div style='text-align:center'>Please login</div>";
} else {


#
# JARVIS temporary folder
#
$jarvis_tmp_folder="tmp/jarvis/";

#
# Environment
#
$env="env.sh";

##
## Style
##

#<div class='view-item'>
#width:32px;
#   filter:alpha(opacity=50);
#   opacity:0.5;
#   background-color:#DD55DD;
#   background-color:rgba(0, 0, 255, 0.5);
#   background-color:#55DD55;
#   background-color:rgba(128, 128, 255, 0.5);
#        -webkit-box-shadow: 0 1px 2px rgba(0,0,0,.05);
#           -moz-box-shadow: 0 1px 2px rgba(0,0,0,.05);
#                box-shadow: 0 1px 2px rgba(0,0,0,.05);
#   float:left;
##e5e5e5;
#   background-color:rgba(255, 165, 0, 0.5);
#   filter:progid:DXImageTransform.Microsoft.gradient(startColorstr=#7FFFA500,endColorstr=#7FFFA500);
#filter:progid:DXImageTransform.Microsoft.gradient(startColorstr=#7F9932cc,endColorstr=#7F9932cc); b22222

#$header="";
#$header.="<HEAD>";

#$header.="<link rel='stylesheet' media='all' type='text/css' href='style.css' />";
#$header.="</HEAD>";

#$header.="<BODY>";

# Header
##########

include "header.php";
#include "search.header.php";





$FRIDAY_name="FRIDAY";
$FRIDAY_description="Full Row Interface to Deal with Annotation dictionnarY";

$FRIDAY_release="V1.0b";
$FRIDAY_release_date="06/05/2020";
$FRIDAY_comment="$FRIDAY_name is an spreadsheet like interface providing genetic variants proritization and visualization in VCF format";




$header="


<TABLE class=' ' width='100%' border=0 height='20px'>
    <TR width='100%'>
        <TD class=' ' width='50px'>
            <img src=favicon.ico>
            <small>$NGS_name</small>
        </TD>
        <TD width='20px'></TD>
        <TD width='*' class=' '>
           <span class='' title=''><B><BIG>$FRIDAY_name |</BIG> $FRIDAY_description</B></span> [$FRIDAY_release]
            <span style=''><BR>$FRIDAY_comment<BR></span>
        </TD>
        <TD width='1px'></TD>
        <TD width='80px' class=' ' style='text-align:center;'>
            <span style=''><B>".$username."</B><BR>".date("d-m-Y")."<BR>".date("H:i:s")."</span>
        </TD>
    </TR>
</TABLE>

";




# cleanup temporary folder
if (1) {
    $command="find $jarvis_tmp_folder -atime 1 -delete ";
    $output_exec = shell_exec($command);
};

## PATCH !
#$dir_howard="$dir_tools/howard/0.9.13b/bin";


# List of filters
# Construct filter from file
#$filter_file="$dir_bin/scripts/filters.ini";
$filter_file=$config_filter_ini; #"$dir_howard/config.filter.ini";
#print $filter_file;
$filter_array=parse_ini_file($filter_file,TRUE);
#echo "Filter_array before:<BR>"; print_r($filter_array); echo "<BR>"; echo "<BR>";
$filter_array=basedon_filters($filter_array);
#echo "Filter_array after:<BR>"; print_r($filter_array); echo "<BR>"; echo "<BR>";
$filter_select="";
#$filter_select.="<OPTION value='' ".(in_array("$filter_name",$filters)?"selected":"")." title='NO filter'>NO filter\n";
$filter_select.="<OPTION value='' ".(in_array("$filter_name",$filters)?"selected":"")." title='NO (re)filter'>NO (re)filter \n";
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


# VCF Files through URL
if (count($vcf_files)>0) {
    #echo "VCF Files through URL";

    #DEV
    #http://localhost:4200/repositories/Archives/SOMATIC/HEMATOLOGY/RUN_TEST_TAG_LISTENER_MIN/P1439/P1439.final.vcf.gz
    #http://192.168.1.14:4200/repositories/Archives/SOMATIC/HEMATOLOGY/RUN_TEST_TAG_LISTENER_MIN/P1439/P1439.final.vcf.gz
    #http://localhost:4200/repositories/Archives/SOMATIC/SOLIDTUMOR/RUN_TEST_TAG_LISTENER_MIN/P1335/P1335.final.vcf.gz
    #http://192.168.1.14:4200/repositories/Archives/SOMATIC/SOLIDTUMOR/RUN_TEST_TAG_LISTENER_MIN/P1335/P1335.final.vcf.gz
    #print_r($vcf_files);
    #print_r(file($vcf_files[0]));

    #$file = $vcf_files[0];

    $analysis=rand(1,1000000);
    $analysis_tmp="$jarvis_tmp_folder/$analysis";
    $command=" mkdir -p $analysis_tmp ";
    $output_exec = shell_exec($command);
    #print "<pre>$output_exec</pre><BR>";

    foreach ($vcf_files as $key => $vcf_file) {
        $new_vcf_file[$key] = "$analysis_tmp/vcf_file_".rand(1,1000000);

        if ( copy($vcf_file, $new_vcf_file[$key]) ) {
            #echo "<br>Copy success! $file to ".$new_vcf_file[$key];
        } else{
            echo "<br>Copy failed.  $file to ".$new_vcf_file[$key];
        }

    };

    $input_files=implode(",",$new_vcf_file);
    $newfile_name="vcf_file_".rand(1,1000000);
    $newfile="$analysis_tmp/$newfile_name";

    #$newfile=$output_file;
    #print_r(file($input_files));

    #$command=" $dir_howard/HOWARD --env=$dir_howard/env.sh --input=$input_files --output=$newfile.uploaded --verbose  > $newfile.uploaded.log 2>$newfile.uploaded.err ;";
    $command=" $dir_howard/HOWARD --env=$env --input=$input_files --output=$newfile.uploaded --verbose ";
    #print "<pre>$command</pre><BR>";

    $output_exec = shell_exec($command);
    // print "<pre>$output_exec</pre><BR>";
    //
    // print "<pre>";
    // print_r(file($output_file));
    // print "</pre><BR>";
    //
    // print "<pre>";
    // print_r(file("$newfile.uploaded.err"));
    // print "</pre><BR>";

    #die();


    #$command="$dir_bcftools/bcftools view $newfile > $newfile.uploaded 2>$newfile.uploaded.err"; # --debug
    #print "<pre>$command</pre><BR>";

    #print "<pre>$output_exec</pre><BR>";
    #$output_exec = shell_exec($command);
    #print "<pre>$output_exec</pre><BR>";
    #$command="rm -f $newfile"; # --debug
    #print "<pre>$command</pre><BR>";
    #$output_exec = shell_exec($command);
    #print "<pre>$output_exec</pre><BR>";
    #die();

    #$_SESSION["loaded_vcf"]=$_FILES["vcf"];
    $t="$newfile.uploaded";
    $t_message="uploaded file ";
    $t_filename="$newfile_name"; # <a href='$t' download>$newfile_name</a>
    $input_run="";
    $input_sample="";


} elseif (is_file($_FILES["vcf"]["tmp_name"])) {
    copy($_FILES["vcf"]["tmp_name"],$_FILES["vcf"]["tmp_name"].".uploaded.tmp");
    #$dir_bcftools
    $command="$dir_bcftools/bcftools view ".$_FILES["vcf"]["tmp_name"].".uploaded.tmp > ".$_FILES["vcf"]["tmp_name"].".uploaded 2>".$_FILES["vcf"]["tmp_name"].".uploaded.err"; # --debug
    #print "<pre>$command</pre><BR>";

    #print "<pre>$output_exec</pre><BR>";
    $output_exec = shell_exec($command);
    #print "<pre>$output_exec</pre><BR>";
    $command="rm -f ".$_FILES["vcf"]["tmp_name"].".uploaded.tmp"; # --debug
    #print "<pre>$command</pre><BR>";
    $output_exec = shell_exec($command);
    #print "<pre>$output_exec</pre><BR>";
    #die();

    #$_SESSION["loaded_vcf"]=$_FILES["vcf"];
    $t=$_FILES["vcf"]["tmp_name"].".uploaded";
    $t_message="uploaded file ";
    $t_filename="".$_FILES["vcf"]["name"]."";
    $input_run="";
    $input_sample="";

    # copy file
} elseif ($input_run != "" && $input_sample !="") {
    #echo "TEST";
    $V1file="";
    $V1file_loc="";
    if (is_file("$dir_analysisV1/$input_run/$input_sample/DATA/MiSeq-Aligner-Alignment/MiSeq-Caller/trakxs/$input_sample.annotated.vcf")) {
        #echo "TEST2";
        $V1file="$dir_analysisV1/$input_run/$input_sample/DATA/MiSeq-Aligner-Alignment/MiSeq-Caller/trakxs/$input_sample.annotated.vcf";
        $V1file_loc="V1/RES";
    } else { #if (is_file("$dir_miseq/$input_run/$input_sample/Data/Intensities/BaseCalls/Alignment/$input_sample.vcf")) {
        #echo "TEST3 $dir_miseq/$input_run/Data/Intensities/BaseCalls/Alignment/".$input_sample."_S*.vcf";
        foreach (glob("$dir_miseq/$input_run/Data/Intensities/BaseCalls/Alignment/".$input_sample."_S*.vcf") as $V1file) {
            #echo "TEST3b $V1file";
            $V1file_loc="V1/RAW";
        };#foreach
        #$V1file="$dir_analysisV1/$input_run/$input_sample/DATA/MiSeq-Aligner-Alignment/MiSeq-Caller/trakxs/$input_sample.annotated.vcf";
    };#if
    if ($V1file!="") {
        echo "V1 filename = $V1file";
        copy($V1file,"/tmp/$input_run.$input_sample.vcf.uploaded");
        $t="/tmp/$input_run.$input_sample.vcf.uploaded";
        $t_message="Found in ";
        $t_filename="$V1file_loc/$input_run/$input_sample";
    } else {
        $V1_VCF_message="No V1 File";
    };#if

} else {
    if (is_file($_REQUEST["vcf_loaded"])) {
        $t=$_REQUEST["vcf_loaded"];
        $t_message="File ";
        $t_filename=$_REQUEST["vcf_loaded_filename"];
    } else {
        # no file
    };#if
};#if
#echo $t;
#print_r(file($_REQUEST["vcf_loaded"]));



if ( !file_exists($t) || filesize($t) == 0 ) {
    if ($_FILES["vcf"]["tmp_name"] == "") {
        $t_filename="";
        $t_message="No file loaded";
    } else {
        $t_filename="";
        $t_message="<span class=error>Error in loading file ".$_FILES["vcf"]["name"]."</span>";
    }; #if

};


# Translate
$command=" $dir_howard/HOWARD --env=$env --input=$t --output=$t.tsv --translation=tsv --verbose ";
#print "<pre>$command</pre><BR>";
$output_exec = shell_exec($command);



$vcf_header.="

    <TABLE  border=0  width='100%' style='margin:0px;' >
                    <TR>
                        <TD width='20px'>&nbsp;</TD>
                        <TD colspan=6 height='40px;' valign='middle'>
                            <B>VCF file</B> ".(($t_filename!="")?"[loaded: $t_message $t_filename] [<a href='$t' download>VCF</a>]  [<a href='$t.tsv' download>TSV</a>]":"[$t_message]")."
                        </TD>
                        <TD width='20px'>&nbsp;</TD>
                    </TR>
                    <TR>
                        <TD width='20px'>&nbsp;</TD>
                        <TD colspan=6 height='40px;' valign='middle'>

                            <TABLE>
                                <TR><TD colspan=3>Upload a VCF file</TD><TD width='20px'>&nbsp;</TD>
                                    <TD><input type='file' name='vcf' class='tbox search' data-original-title=''/>
                                            <input type='hidden' value='$t' name='vcf_loaded' data-original-title=''/>
                                        <input type='hidden' value='$t_filename' name='vcf_loaded_filename' data-original-title=''/>
                                    </TD>
                                </TR>
                            </TABLE>




                        </TD>
                        <TD width='20px'>&nbsp;</TD>
                    </TR>
                    <!--<TR>
                        <TD width='20px'>&nbsp;</TD>
                        <TD colspan=6 height='40px;' valign='middle'>

                            <TABLE>
                                <TR><TD colspan=3>OR try to found a VCF file on V1 file server (under development)</TD></TR>
                                <TR><TD width='20px'>&nbsp;</TD><TD>RUN:</TD><TD><input type='text' name='run' value='$input_run' class='tbox search' data-original-title=''/></TD></TR>
                                <TR><TD width='20px'>&nbsp;</TD><TD>SAMPLE:</TD><TD><input type='text' name='sample' value='$input_sample'class='tbox search' data-original-title=''/></TD></TR>
                            </TABLE>
                        $V1_VCF_message
                        </TD>
                        <TD width='20px'>&nbsp;</TD>
                    </TR>-->
    </TABLE>


    ";

$filter_form="
                <TABLE width='100%' border=0 ><TR><TD><TABLE width='100%' border=0 class=''>
                    <TR><TD>

                <div style='text-align: left; table-border:1; clear:left;' class='' >
                <form action='?' method='post' class='form-inline' enctype='multipart/form-data' name='main' id='main'>
                <INPUT type='hidden' name='action_todo' value='$action_todo'>
                <INPUT type='hidden' name='sample' value='$input_sample'>
                <INPUT type='hidden' name='run' value='$input_run'>
                <INPUT type='hidden' name='q' value='$q'>

                $vcf_header

                <TABLE  border=0  width='100%' >
                    <TR>
                        <TD width='20px'>&nbsp;</TD>
                        <TD colspan=6 height='40px;' valign='middle' width='200px'>
                            <B>Prioritization/Filters</B>
                        </TD>
                        <TD width='20px'>&nbsp;</TD>
                        <TD valign='top' halign='center' width=''>

                        </TD>
                        <TD width='20px'>&nbsp;</TD>
                    </TR>
                    <!--<TR>
                        <TD width='20px;'>&nbsp;</TD>
                        <TD valign='top' halign='center' width='200px'>

                        </TD>
                        <TD valign='top' halign='center' width='200px'>
                            Annotations
                        </TD>
                        <TD width='20px;'>&nbsp;</TD>
                        <TD valign='top' halign='center' width=''>

                        </TD>
                        <TD width='20px'>&nbsp;</TD>
                    </TR>-->
                    <TR>

                        <TD width='20px;'>&nbsp;</TD>
                        <TD valign='top' halign='center' width='200px'>
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
                                        [<SPAN onclick=\"javascript:var newWin = window.open('filters.php','popup','menubar=yes, status=no, scrollbars=yes, menubar=yes, width=800, height=600, visible=yes');document.getElementById('main').action='filters.php';document.getElementById('main').target='popup';document.getElementById('main').submit();newWin.focus();document.getElementById('main').action='$action_default';document.getElementById('main').target='_self'; return false;\" style='cursor:pointer;'>details</SPAN>]
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

                            <!--
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
                                            <OPTION value='' ".($orderby==""?"selected":"").">NO order
                                            <OPTION value='PZScore:n' ".($orderby=="PZScore"?"selected":"").">PZScore
                                            <OPTION value='symbol:' ".($orderby=="symbol"?"selected":"").">symbol
                                            <OPTION value='CHROM:' ".($orderby=="CHROM"?"selected":"").">CHROM
                                            <OPTION value='VARTYPE:' ".($orderby=="VARTYPE"?"selected":"").">VARTYPE
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
                            -->

                        </TD>
                        <TD width='20px;'>&nbsp;</TD>
                    </TR>
                    <TR>
                        <TD width='20px;'>&nbsp;</TD>
                        <TD>
                            <INPUT TYPE='submit' value='Process' name='process' class='btn button search btn size-medium bg-blue text-white shadow hover-moveup' style='margin-right: 5px;' data-original-title=''/>

                        </TD>
                        <TD width='20px;'>&nbsp;</TD>
                         
                        
                    </TR>
                </TABLE>
                </div>
                </form>

                </TD>
                    </TR>
                </TABLE>

                </TD>
                    </TR>
                </TABLE>


                ";


$text.=$filter_form."<BR>";


#<OPTION value='VDScore' ".($orderby=="VDScore"?"selected":"").">VDScore



#if (1) {
#if ($t != "") {
if ( file_exists($t) && filesize($t) != 0 && $process) {

    #print "<pre>"; print_r($t); print "</pre>";
    #$t_content=file_get_contents ( $t );
    #print "<pre>$t_content</pre>";
    #die();

    # Prioritization
    $filters_option=implode(",",$filters);
    # $dir_howard

    #$gene_filter_grep_symbol=($gene_filter=="")?"":"symbol=$gene_filter";
    #$gene_filter_grep_hgvs=($gene_filter=="")?"":"hgvs=$gene_filter";
    # hgvsEnsembl
    $gene_filter_grep="";
    $gene_filter_array=explode(" ",$gene_filter);
    foreach ($gene_filter_array as $gene_filter_one) {
        if ($gene_filter_one != "") {
            #print "<pre>'$gene_filter_one'</pre><BR>";;
            #$gene_filter_grep.=" -e 'symbol=$gene_filter_one' -e 'hgvs=$gene_filter_one'  -e 'NOMEN=$gene_filter_one' ";
            #$gene_filter_grep.=" -e '[;'$'\t'']symbol=[^;'$'\t'']*$gene_filter_one[^;'$'\t'']*[;'$'\t'']' -e '[;'$'\t'']hgvs=[^;'$'\t'']*$gene_filter_one[^;'$'\t'']*[;'$'\t'']'  -e '[;'$'\t'']NOMEN=[^;'$'\t'']*$gene_filter_one[^;'$'\t'']*[;'$'\t'']' -e '[;'$'\t'']hgvsEnsembl=[^;'$'\t'']*$gene_filter_one[^;'$'\t'']*[;'$'\t'']' ";
            $gene_filter_grep.=" -e '[; ]symbol=[^; ]*$gene_filter_one"."[^;    ]*[;    ]' ";
            $gene_filter_grep.=" -e '[; ]hgvs=[^;   ]*$gene_filter_one"."[^;    ]*[;    ]' ";
            $gene_filter_grep.=" -e '[; ][^;    ]*NOMEN=[^; ]*$gene_filter_one"."[^;    ]*[;    ]' ";
            $gene_filter_grep.=" -e '[; ]hgvsEnsembl=[^;    ]*$gene_filter_one"."[^;    ]*[;    ]' ";

        };#if
        # -e 'GNOMEN=$gene_filter_one'
    };#foreach
    if ($gene_filter_grep == "") {
        $gene_filter_grep=" -e '' ";
    };


    $global_filter_grep="";
    $global_filter_array=explode(" ",$global_filter);
    foreach ($global_filter_array as $global_filter_one) {
        if ($global_filter_one != "") {
            #print "<pre>$global_filter_one</pre><BR>";;
            #$global_filter_grep.=" -e '$global_filter_one'";
            $global_filter_grep.=" -e '[;   ][^;    ]*=[^;  ]*$global_filter_one"."[^;  ]*[;    ]' ";
        };#if
    };#foreach
    if ($global_filter_grep == "") {
        $global_filter_grep=" -e '' ";
    };

    $command_plus=" --split=10000 --tmp=$dir_tmp ";
    if ($hardfiltering) {
        $command_plus=" --hard ";
        $hard_filter_grep=" | grep -i -E -e '^#' -e 'PZFlag=PASS' ";
    };

    if ($orderby != "") {
        #$orderby_options=" --translation=vcf --sort_by=$orderby --order_by=$ascdesc ";
        $orderby_options=" --translation=vcf --sort='$orderby:$orderby_type:$ascdesc' ";
    };#if

    #config prioritization $filter_file


    $calculation=" --calculation='VAF,NOMEN,VARTYPE' " ;

    #$fields_default="location,outcome,CLINVAR,COSMIC,ALL";
    $fields_default="location,outcome,CLINVAR,COSMIC,snpeff_impact,snpeff_annotation,dbSNP,popfreq,INTERVAR,interpro_domain,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,gnomAD";

    #$fields=" --fields='NOMEN,VAF,VAF_average,PZFlag,location,outcome'  ";
    #$fields=" --fields='NOMEN,VAF,VAF_average,PZFlag,location,outcome,ANN'  ";
    #$fields=" --fields='NOMEN,VAF_average,location,outcome,ALL'  ";
    $fields=" --fields='NOMEN,PZFlag,VAF_average,$fields_default'  ";

    #$calculation="  " ;

    #$annotation=" --annotation='snpeff_split' " ;

    #print "<pre>$global_filter</pre><BR>";
    #print "<pre>";
    #print_r($global_filter_array);
    #print "</pre><BR>";
    #print "<pre>$global_filter_grep</pre><BR>";

    if ($filters_option != "") {

        $command=" cat $t |  grep -i -E -e '^#' $gene_filter_grep | grep -i -E -e '^#' $global_filter_grep > $t.prioritized.vcf.tmp ; $dir_howard/HOWARD --env=$env --config_prioritization=$filter_file --input=$t.prioritized.vcf.tmp --output=$t.prioritized.vcf --prioritization=$filters_option --pzfields='PZScore,PZFlag,PZComment,PZInfos' $orderby_options $calculation $annotation --verbose --force $command_plus 1>$t.prioritized.vcf.log 2>$t.prioritized.vcf.err ;";

    } else {
        #$command="cp $t $t.prioritized.vcf"; # --debug
        #$command=" cat $t |  grep -i -E -e '^#' $gene_filter_grep | grep -i -E -e '^#' $global_filter_grep > $t.prioritized.vcf ";
        $command=" cat $t |  grep -i -E -e '^#' $gene_filter_grep | grep -i -E -e '^#' $global_filter_grep $hard_filter_grep > $t.prioritized.vcf.tmp ; ";
        if ($orderby_options != "") {

            $command.=" $dir_howard/HOWARD --env=$env --config_prioritization=$filter_file --input=$t.prioritized.vcf.tmp --output=$t.prioritized.vcf $orderby_options $calculation $annotation --verbose $command_plus > $t.prioritized.vcf.log 2>$t.prioritized.vcf.err ;";
            # $command.=" $dir_howard/VCFtranslation.pl --input=$t.prioritized.vcf.tmp --output=$t.prioritized.vcf --translation=vcf --sort_by=$orderby --order_by=$ascdesc --verbose $command_plus > $t.prioritized.vcf.log ;";
        } else {
            $command.=" mv $t.prioritized.vcf.tmp $t.prioritized.vcf > $t.prioritized.vcf.log 2>$t.prioritized.vcf.err ;";
        };#if
    }; #if


    #$command="bcftools";

    #print "$command<BR>";
    #die();
    $output_exec = shell_exec($command);
    #print "<pre>$output_exec</pre><BR>";
    #$vcf=file_get_contents("$t.prioritized.vcf");
    #$log=file_get_contents("$t.prioritized.vcf.log");
    #$err=file_get_contents("$t.prioritized.vcf.err");
    #print "<pre>$vcf</pre><BR>";
    #print "<pre>$log</pre><BR>";
    #print "<pre>$err</pre><BR>";
    #die();



    # Version avec Content
    #$Variants_VCFContent=file_get_contents("$t.prioritized.vcf.log");
    #print "<pre>$Variants_VCFContent</pre><BR>";
    #exit();
    #$Variants_HTMLContent=VCFtoHTML($Variants_VCFContent,$sample_id,$annotation_list,array("hardfiltering"=>$hardfiltering,"vcf_ids"=>$variant_files,"limit"=>100000,"orderby"=>$orderby,"ascdesc"=>$ascdesc,"gene_filter"=>$gene_filter,"global_filter"=>$global_filter,"report"=>0));
    #$Variants_HTMLContent=VCFtoHTML($Variants_VCFContent,$sample_id,$annotation_list,array("hardfiltering"=>$hardfiltering,"vcf_ids"=>$variant_files,"limit"=>$limit,"orderby"=>$orderby,"ascdesc"=>$ascdesc,"gene_filter"=>$gene_filter,"global_filter"=>$global_filter,"report"=>0));


    # Version avec fichier
    $Variants_HTMLContent=VCFFiletoHTML("$t.prioritized.vcf",$sample_id,$annotation_list,array("hardfiltering"=>$hardfiltering,"vcf_ids"=>$variant_files,"limit"=>$limit,"orderby"=>$orderby,"ascdesc"=>$ascdesc,"gene_filter"=>$gene_filter,"global_filter"=>$global_filter,"report"=>0));



    # Translate
    $command=" $dir_howard/HOWARD --env=$env --input=$t.prioritized.vcf --output=$t.prioritized.tsv --translation=tsv $fields --verbose ";
    #print "<pre>$command</pre><BR>";
    $output_exec = shell_exec($command);

    
    $vcf=file_get_contents("$t.prioritized.vcf");
    $tsv=file_get_contents("$t.prioritized.tsv");
    //print "<pre>$vcf</pre><BR>";
    //print "<pre>$tsv</pre><BR>";
    

    //open up your file and convert it to an array
    
    $file = "$t.prioritized.tsv";
    #echo $file;
    $data = tsv_to_array($file,array('header_row'=>true,'remove_header_row'=>true));



    $vcf_array=file("$t.prioritized.vcf");
    #print_r($vcf_array);

    foreach ($vcf_array as $key=>$line) {
        #print($line[0].$line[0]."<br>");
        if ($line[0]=="#" & $line[1]=="#") {

            #print($line."<br>");
            preg_match('/##(.*)=<ID=(.*),Number=(.*),Type=(.*),Description=(.*)>/', $line, $matches_HEADER, PREG_OFFSET_CAPTURE);
            
            
            if (count($matches_HEADER)>0) {
                #print_r($matches_HEADER);
                $vcf_header_line_section=$matches_HEADER[1][0];
                $vcf_header_line_ID=$matches_HEADER[2][0];
                $vcf_header_line_Number=$matches_HEADER[3][0];
                $vcf_header_line_Type=$matches_HEADER[4][0];
                $vcf_header_line_Description=$matches_HEADER[5][0];

                $field=preg_replace("/[^a-zA-Z0-9]/", "", $vcf_header_line_ID);
                if (is_numeric($field[0])) {
                    $field="_".$field;
                };

                $vcf_header_array[$vcf_header_line_section][$$vcf_header_line_ID]=array(
                    "ID"=>$vcf_header_line_ID,
                    "Number"=>$vcf_header_line_Number,
                    "Type"=>$vcf_header_line_Type,
                    "Description"=>$vcf_header_line_Description
                );
                $vcf_header_array["ALL"][$field]=$vcf_header_array[$vcf_header_line_section][$$vcf_header_line_ID];
               
            }

        };
    };
    // echo "<pre>";
    // #print_r($vcf_header_array["ALL"]);
    // print_r($vcf_header_array);
    // echo "</pre>";

    $VCF_header_Type_translate=array(
        "String"=>"text",
        "Integer"=>"numeric",
        "Float"=>"numeric",

    );


    # Translate
    #$command=" grep ^## $t.prioritized.vcf  ";
    #print "<pre>$command</pre><BR>";
    #$output_exec = shell_exec($command);
    #print "<pre>$output_exec</pre><BR>";

    #preg_match('/##(.*)=(.*)/', $output_exec, $matches_HEADER, PREG_OFFSET_CAPTURE);
    #print_r($matches_HEADER);
    //this will print out the array to the screen
    #echo '<pre>'.print_r($data[2],true).'</pre>';
    
    /*
    echo "<pre>
      {
        id: 1,
        flag: 'EUR',
        currencyCode: 'EUR \\ntest2',
        currency: 'Euro',
        level: 0.9033,
        units: 'EUR / USD',
        asOf: '08/19/2019',
        onemore: 0.0026,
        onedChng: 0.0026
      }

      </pre>
      ";
    */

    # default
    $maxcharcol_default=50;
    $maxwidthcol_default=200;

      foreach ($data as $key => $value) {
        
        if ($value['#CHROM']!="") {


            // $location=preg_replace("/[,&]/", " ", $value['location']);
            // $outcome=preg_replace("/[,&]/", " ", $value['outcome']);
            // $CLINVAR=preg_replace("/[,&]/", " ", $value['CLINVAR']);
            // $COSMIC=preg_replace("/[,&]/", " ", $value['COSMIC']);

            # NOMEN
            $columns_fields_array["NOMEN"]["type"]='text';
            $columns_fields_array["NOMEN"]["width"]='300';
            $fixedColumnsLeft_array["NOMEN"]=1;

            # Comment
            $columns_fields_array["Comment"]["type"]='text';
            $columns_fields_array["Comment"]["width"]='100';
            $fixedColumnsLeft_array["Comment"]=1;


            # PZFlag and rest
            $PZFlag_ALL="";
            $fields_rest="";
            $default_field_type="text";
            $i=0;

            $columns_fields_array_PZFlag=array();
            $columns_fields_array_PZComment=array();
            $columns_fields_array_PZrest=array();

            #$fixedColumnsLeft_array=array();

            foreach ($value as $key2 => $value2) {

                $i++;
                $maxwidthcol="";

                #if ($i>=186 & $i<=186) {
                if ($i>7) {
                #if (1) {

                if ($key2 != "NOMEN") {


                #$field=preg_replace("/[-#]/", array("","_"), $key2);
                #$field="_".preg_replace("/[-#.+]/", "", $key2);
                #$field=preg_replace("/[^a-zA-Z0-9]/", "_", $key2);
                $field=preg_replace("/[^a-zA-Z0-9]/", "", $key2);
                if (is_numeric($field[0])) {
                    $field="_".$field;
                };
                #$field_value=$value2;
                $field_value=preg_replace("/[']/", "", $value2);
                #$field_value=preg_replace("/[,]/", " ", $field_value);
            
                $field_type=$VCF_header_Type_translate[$vcf_header_array["ALL"][$field]["Type"]];
                if ($field_type=="") {
                    $field_type=$default_field_type;
                };

                if (strlen($field_value)>=$maxcharcol_default) {
                    $maxwidthcol=$maxwidthcol_default;
                } 

                #echo "<pre>$key2 => $value2</pre>";

                    preg_match('/PZFlag(.*)/', $key2, $matches_PZFlag, PREG_OFFSET_CAPTURE);
                    preg_match('/PZComment(.*)/', $key2, $matches_PZComment, PREG_OFFSET_CAPTURE);
                    if (count($matches_PZFlag)>0 & $key2!="PZFlag-default" & $key2!="PZFlag") {
                        $PZFlag_ALL.=$field.": '".$field_value."', ";
                        $columns_fields_array_PZFlag[$field]["type"]='text';
                        #$fixedColumnsLeft_array[$field]=1;
                    // } elseif (count($matches_PZComment)>0 & $key2!="PZComment-default" & $key2!="PZComment") {
                    //     $PZComment_ALL.=$field.": '".$field_value."', ";
                    //     $columns_fields_array_PZComment[$field]["type"]='text';
                    //     #$fixedColumnsLeft_array[$field]=1;
                    } else {
                        $fields_rest.=$field.": '".$field_value."', ";
                        $columns_fields_array_rest[$field]["type"]=$field_type;
                    };

                    #$columns_fields_array_rest[$field]["autoWrapCol"]='true';
                    if ($maxwidthcol!="") {
                        $columns_fields_array_rest[$field]["width"]=$maxwidthcol; # editor: false
                    };
                    
                    #$columns_fields_array[$field]["editor"]='false'; # editor: false

                };

                };

            };
            
            /*
            echo "<br><br>line<br>";
            echo "<pre>";
            print_r($columns_fields_array_PZFlag);
            echo "</pre>";
            echo "<pre>";
            print_r($columns_fields_array_PZComment);
            echo "</pre>";
            
            echo "<pre>";
            print_r($columns_fields_array_rest);
            echo "</pre>";
            */

            $columns_fields_array=array_merge($columns_fields_array,$columns_fields_array_PZFlag,$columns_fields_array_PZComment,$columns_fields_array_rest);


            $dataObject_array[]="

                {
                    NOMEN: '".$value["NOMEN"]."',
                    $PZFlag_ALL
                    $PZComment_ALL
                    $fields_rest
                }

            ";
            
                    #PZFlag: '".$value["PZFlag"]."',
                    #VAF: ".$value["VAF_average"].",
                    #location: '".$location."',
                    #outcome: '".$outcome."',
                    #CLINVAR: '".$CLINVAR."',
                    #COSMIC: '".$COSMIC."',

            
            /*
            $columns_fields_array["VAF"]["type"]='numeric';
            $columns_fields_array["VAF"]["pattern"]='0.0000%';

            $columns_fields_array["location"]["type"]='text';

            $columns_fields_array["outcome"]["type"]='text';

            $columns_fields_array["CLINVAR"]["type"]='text';

            $columns_fields_array["COSMIC"]["type"]='text';
            */

        }

      }


/*
echo "<pre>";
print_r($dataObject_array);
echo "</pre>";
*/

// echo "<pre>";
// print_r($columns_fields_array);
// echo "</pre>";



$columns_fields_list="";
$colHeaders_fields_list="";
foreach ($columns_fields_array as $data => $value) {
    # code...
    $columns_fields_list_array="data: '$data'";
    
    foreach ($value as $key => $value2) {
        $columns_fields_list_array.=", $key: '$value2'";
    };

      $columns_fields_list.="{ ". $columns_fields_list_array." }, ";
      $colHeaders_fields_list.="'$data',";

};

#$columns_fields_list=str_replace('/,$/',"",$columns_fields_list);
#$colHeaders_fields_list=str_replace('/,$/',"",$colHeaders_fields_list);

$columns_fields_list=rtrim($columns_fields_list, ",");
$colHeaders_fields_list=rtrim($colHeaders_fields_list, ",");
 

$fixedColumnsLeft=count($fixedColumnsLeft_array);

#echo $columns_fields_list; #die(0);
#echo $colHeaders_fields_list; #die(0);



    $dataObject_variants=implode(",", $dataObject_array);
    #print "<pre>$dataObject_variants</pre>";


    $test_download='
    <div style="text-align: left; table-border:1; " class="">
    <table border="0" width="100%" class="  filter">
                <tbody><tr>
                <td width="20px"></td>
                <td>
                <br>
                Download Prioritized Files
                [<a href="'."$t.prioritized.vcf".'" download>VCF</a>]
                [<a href="'."$t.prioritized.tsv".'" download>TSV</a>]
                <br><br></td></tr>
                </tbody>
    </table>
    </div>
    <br>
    ';


    $text.=$test_download;

    $text.="<div style='text-align: left; table-border:1; ' class=''>
            $Variants_HTMLContent
            </div>
    ";


    # Remove files
    #$command="rm -f $t.prioritized.vcf"; #
    #$output_exec = shell_exec($command);
    #print "<pre>$output_exec</pre><BR>";


} else {

    #$text.="<div style='text-align: left; table-border:1; ' class='error'>
    #       Error in loading file...
    #       </div>
    #";


};#if


};#if

    $text.= "
            </div>
        </div>
    ";





$title="

";



$dataObject_list=$dataObject_variants;

$dataObject="
var dataObject = [
 $dataObject_list
];
";


$columns="[
    $columns_fields_list
]";


$colHeaders="[
    $colHeaders_fields_list
]";


#print "<pre>$dataObject</pre>";

#exit(0);


?>




 


<style>

#green {
    width: 100%;
    min-width: 500px;
    height: 20%;
    min-height: 400px;
    border: 0px solid green;
    /* float: left; */
    margin: 5px;
    padding: 5px;
}
#red {
    width: 100%;
    height: 90%;
    border: 0px solid red;    
    /*float: left;*/
}
#blue {
    width: 310px;
    height: 100px;
    border: 0px solid blue;
    margin: 5px;
    clear: both;
    padding: 5px;
}
#hidden {
    display: none;
}


</style>

<title><?php print "$NGS_name/$FRIDAY_name"; ?></title>

</head>




<body>


    

    <div id="green">

        <?php print $header; ?>
        <hr>
        <?php print $filter_form; ?> 

        <div id="export-buttons" class="visible">
            
            <!--<button id="export-string" class="btn size-medium bg-green text-white shadow hover-moveup">Export to string (see the browser console)</button>-->
        </div>

    </div>
    
    <div id="red">
            <button id="export-csv" class="btn size-medium bg-blue text-white " style="margin-right: 5px;">CSV</button>
          <div id="hot"></div>
    </div>

    <!--
    <div id="blue">
        
         
    </div>
    -->

<div id="#hidden">

<script>

<?php print $dataObject; ?>

var currencyCodes = ['EUR', 'EURR', 'JPY', 'GBP', 'CHF', 'CAD', 'AUD', 'NZD', 'SEK', 'NOK', 'BRL', 'CNY', 'RUB', 'INR', 'TRY', 'THB', 'IDR', 'MYR', 'MXN', 'ARS', 'DKK', 'ILS', 'PHP'];
var flagRenderer = function (instance, td, row, col, prop, value, cellProperties) {
  var currencyCode = value;
  while (td.firstChild) {
    td.removeChild(td.firstChild);
  }
  if (currencyCodes.indexOf(currencyCode) > -1) {
    var flagElement = document.createElement('DIV');
    flagElement.className = 'flag ' + currencyCode.toLowerCase();
    td.appendChild(flagElement);
  } else {
    var textNode = document.createTextNode(value === null ? '' : value);

    td.appendChild(textNode);
  }
};
var hotElement = document.querySelector('#hot');
var hotElementContainer = hotElement.parentNode;
var hotSettings = {
  data: dataObject,
  columns: <?php print $columns ?>,
  stretchH: 'all',
  //width: '800',
  autoWrapRow: true,
  height: '100%',
  //height: 800,
  maxRows: 1000,
  rowHeaders: true,
  colHeaders: <?php echo $colHeaders; ?>,
  mergeCells: true,
  contextMenu: true,
  columnSorting: {
    indicator: true
  },
  autoColumnSize: {
    samplingRatio: 23
  },
  manualRowResize: true,
  manualColumnResize: true,
  manualRowMove: true,
  manualColumnMove: true,
  filters: true,
  dropdownMenu: true,
  fixedRowsTop: 0,
  fixedRowsBottom: 0,
  fixedColumnsLeft: <?php echo $fixedColumnsLeft; ?>,
  multiColumnSorting: {
    indicator: true
  },
  exportFile: true,
  trimRows: [
    0,
    0,
    0,
    0
  ],
  licenseKey: 'non-commercial-and-evaluation'
};
var hot = new Handsontable(hotElement, hotSettings);
document.getElementById("export-csv").addEventListener("click", function(event) { hot.getPlugin("exportFile").downloadFile("csv", {filename: "CSV"});})
document.getElementById("export-string").addEventListener("click", function(event) {console.log(hot.getPlugin("exportFile").exportAsString("csv"));})


</script>

</div>

</body>
</html>                            
