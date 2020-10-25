<?PHP

$NGS_name="JARVIS";
$NGS_description="Just A Rather Marvellous Variant Interface System";
#$NGS_description_title="Just an Annotation and Reporting Variants Interpretation System";
#$NGS_description_title="Just A Rather Variant Interpretation System";
#$NGS_description_title2="Just A Rather cool Variant Interpretation System";

$NGS_release="V2.5b";
$NGS_release_date="07/11/2015";
$NGS_comment="$NGS_name is an efficient software tool functionally annotating genetic variants in order to detected causative cancer mutations";


$VISION_name="VISION";
$VISION_description="Vcf vIewer and Simplified prioritizatION";

$VISION_release="V1.0.2b";
$VISION_release_date="03/12/2019";
$VISION_comment="$VISION_name is an simplyfied interface providing genetic variants proritization and visualization in VCF format";




# Directories
$dir_analysis="/media/IRCV2/RES/ALL";
$dir_miseq="/media/IRCV2/RAW/MSR";
$dir_bin="/media/IRCV2/NGSEnv/bin";
$dir_scripts="/media/IRCV2/NGSEnv/scripts";
#$dir_tools="/media/IRCV2/NGSEnv/tools";
#$dir_tools=getcwd()."/tools";
$dir_tools="/home/TOOLS/tools";
#$dir_howard="$dir_tools/howard/current";
#$dir_howard="$dir_tools/howard/0.9.13b/bin";
$dir_howard="$dir_tools/howard/current/bin";
#$dir_bcftools="$dir_tools/bcftools/current/bin";
$dir_bcftools="$dir_tools/bcftools/current/bin";
$dir_pepper="$dir_tools/pepper/current";
#$customNM="$dir_bin/scripts/CustomNM.txt";
$customNM="config/config.customNM.txt";
#$dir_tmp="/media/IRCV2/NGSEnv/tmp";
$dir_tmp=getcwd()."/tmp";
$config_filter_ini="config/config.prioritization.ini";

# INFOS


# CustomNM
#print "customNM".$customNM;
$customNM_infos=array();
#$customNM_infos=array("???");
if (file_exists($customNM)) {
	$customNM_array=file($customNM);
	#print "<PRE>"; print_r($customNM_array);  print "</PRE>";
	foreach ($customNM_array as $k=>$line) {
		$line_split=explode("\t",$line);
		if ($k==0) {
			foreach ($line_split as $col=>$info) {

				$customNM_colname[$col]=str_replace("#","",$info);
				$first_col=($col==0)?$info:$first_col;
			};#foreach
		} else {
			if (trim($line_split[0])!="") {
				$customNM_infos[trim($line_split[0])][trim($line_split[1])][trim($line_split[2])]=trim($line_split[2]);
				#print_r($line_split);
			};#if
		};#if
	};#foreach
};#if
#print "<PRE>"; print_r($customNM_infos);  print "</PRE>";



?>
