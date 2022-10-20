<?php


#############
### INFOS ###
#############

$APP_SECTION="Reports";



################
### INCLUDES ###
################


### CONFIG
############

include "config.inc.php";


### FUNCTIONS
###############

include "functions.inc.php";


### HEADER
############

include "header.inc.php";


### VARIABLES
###############



### DATA
##########

include "scan.inc.php";


// Script start time
$rustart = getrusage();



### VARIABLES
###############


### PATH

# HOME

$HOME=$_REQUEST["HOME"];
if ($HOME=="" && $_ENV["FOLDER_REPOSITORIES"]!="") {
	$HOME=$_ENV["FOLDER_REPOSITORIES"];
} elseif ($HOME=="" && $_ENV["DOCKER_STARK_SERVICE_EDITH_DASHBOARD_SUBFOLDER_REPOSITORIES"]!="") {
	$HOME=$_ENV["DOCKER_STARK_SERVICE_EDITH_DASHBOARD_SUBFOLDER_REPOSITORIES"];
} elseif ($HOME=="" || !is_dir($HOME)) {
	$HOME="repositories";
};


# PATH & SEARCH input
$PATH=$_REQUEST["PATH"];
$SEARCH=$_REQUEST["search"];

# PATH from SEARCH
$INPUT_SAMPLE=$_REQUEST["sample"]!=""?$_REQUEST["sample"]:"*";
$INPUT_ANALYSIS=$_REQUEST["analysis"]!=""?$_REQUEST["analysis"]:"*";
$INPUT_ANALYSIS_FROM_LIST=$_REQUEST["analysis_from_list"]!=""?$_REQUEST["analysis_from_list"]:"";
if ($INPUT_ANALYSIS_FROM_LIST != "") {
	$INPUT_ANALYSIS=$INPUT_ANALYSIS_FROM_LIST;
}

if ($INPUT_SAMPLE!="*" || $INPUT_ANALYSIS!="*") {
	$PATH=$HOME."/*/*/*/".$INPUT_ANALYSIS."/".$INPUT_SAMPLE;
}


# DEFAULT PATH
if ($PATH=="") {
	$PATH=$HOME;
}


# Scan

#$PATH="repositories/Repository/HUS*/*/*/*/";
$PATH_FULL=path_full($PATH);
$PATH_FULL_array=explode("/",$PATH_FULL);
#$PATH_FULL_PREG="repositories\/Repository\/.*\/.*\/.*\/.*\/";
#$PATH_FULL_PREG="repositories/Repository/HUS.*/.*/.*/.*/";
$PATH_FULL_PREG=str_replace("*",".*",$PATH_FULL);
#$samples_paths_original=preg_grep("#$PATH_FULL_PREG#", $runs_repositories);
$samples_paths_corrected=preg_grep("#$PATH_FULL_PREG#", $runs_repositories_corrected);
#$samples_paths=$samples_paths_corrected; #preg_grep("#$PATH_FULL_PREG#", $runs_repositories_corrected);

// echo "PATH=$PATH";
// echo "PATH_FULL_PREG=$PATH_FULL_PREG";
// echo "<pre>"; print_r($runs_repositories_corrected); echo "</pre>";

$samples_paths=null;
foreach ($samples_paths_corrected as $key=>$path) {
	#$samples_paths[$key]=$runs_repositories[$key];
	$samples_paths[$key]=$runs_repositories_corrected[$key];
}

// echo "<pre>"; print_r($samples_paths); echo "</pre>";



# PATH HTML & LINKS
// $PATH_HTML=path_html($PATH,$samples_paths);
// $PATH_LINKS=path_links($PATH,$samples_paths);
$PATH_HTML=path_html($PATH,$samples_paths_corrected);
$PATH_LINKS=path_links($PATH,$samples_paths_corrected);

// echo "#### test";
// echo '<pre>';
// echo "PATH=$PATH<br>";
// echo "PATH_FULL=$PATH_FULL<br>";
// echo "<br>";
// echo "PATH_FULL_array<br>";
// print_r($PATH_FULL_array);
// echo "<br>";
// echo "samples_paths_original<br>";
// print_r($samples_paths_original);
// echo "<br>";
// echo "samples_paths_corrected<br>";
// print_r($samples_paths_corrected);
// echo "<br>";
// echo "samples_paths<br>";
// print_r($samples_paths);
// echo "<br>";
// echo "path_to_tree<br>";
// print_r(path_to_tree($samples_paths));
// echo '</pre>';


### MAX REPORTS
$MAX_REPORTS=$_REQUEST["MAX_REPORTS"];
if ($MAX_REPORTS=="") {
	$MAX_REPORTS=$configuration["parameters"]["max_reports"]+0;
};
if ($MAX_REPORTS=="" || $MAX_REPORTS<=0) {
	$MAX_REPORTS=1000;
};

#echo "<br><br><br>MAX_REPORTS=$MAX_REPORTS";


# DEV
$VERBOSE=0;


# DEBUG
if (0) {
	echo "<pre>";
	echo $HOME;
	echo $PATH;
	echo $PATH_HTML;
	echo $PATH_LINKS;
	echo "</pre>";
};




###############
### CONTENT ###
###############


### HEADER PATH
###############


$HEADER_PATH= '
	<section class="header1 cid-ru7OEConn1" id="header16-1k">

		<div class="container div-wrapper">

			<h3 class="mbr-section-subtitle mbr-fonts-style align-left mbr-light display-5">
				<small><small>
				<div class="media">
					'.$PATH_HTML.'
				</div>
				<br>

				<div class="div-wrapper">

				</div>

				</small></small>
		  	</h3>
	    </div>
	</section>

';

echo $HEADER_PATH;


### TABLE
###########


### THEAD

$thead='
	<tr class="table-heads">
		<th class="head-item mbr-fonts-style display-7">
			 
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Sample
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Analysis / Run
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Project
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Group
		</th>
		<th class="head-item mbr-fonts-style display-7">
			Repository
		</th>
	</tr>
';




### TBODY


$tbody="";


#if (count($reports)<=$MAX_REPORTS) {
if (count($samples_paths)<=0) {

	$CONTENT_SECTION_REPORTS= '

	<section class="section-table cid-ruiBoanIwc" id="REPORTS">

	  <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
	  	No reports found. Please browse folders...
	  </p>

	</section>
	';


} elseif (count($samples_paths)>$MAX_REPORTS) {

	$CONTENT_SECTION_REPORTS= '

	<section class="section-table cid-ruiBoanIwc" id="REPORTS">

	  <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
	  	Too many reports selected. Please browse folders...
	  </p>

	</section>
	';



} else {

	$samples_path_pattern=null;
	#echo "<pre>"; print_r($samples_paths); echo "</pre>";

	foreach ($samples_paths as $key=>$samples_path) {
		#dirname($samples_path);
		$samples_path_pattern[]=dirname($samples_path);
		#$samples_path_pattern[]=dirname($samples_paths_original[$key]);
	};
	#echo "<pre>"; print_r($samples_path_pattern); echo "</pre>";
	#$samples_path_pattern_brace="{".join(",",array_unique($samples_path_pattern))."}";
	$samples_path_pattern_brace="(".join("|",array_unique($samples_path_pattern)).")";
	#echo "<br>samples_path_pattern_brace=$samples_path_pattern_brace";
	#echo "<br>samples_path_pattern_brace2=$samples_path_pattern_brace2";

	#$samples_files=array_unique(glob($samples_path_pattern_brace."/{*stark.report.html,*tag,*analysis.tag}",GLOB_BRACE));

	#print_r($samples_files);
	#echo "<br><br>";

	#echo "<pre>"; print_r($runs_repositories_samples_files); echo "</pre>";
	#$samples_files=preg_grep("#(stark.report.html|tag)$#", $runs_repositories_samples_files);
	#$samples_files=preg_grep("#$samples_path_pattern_brace/.*(stark.report.html|tag)$#", $runs_repositories_samples_files);
	#$samples_files=preg_grep("#$samples_path_pattern_brace/.*(stark.report.html|stark.report.pdf|tag)$#", $runs_repositories_samples_files);
	$samples_files=preg_grep("#$samples_path_pattern_brace/.*(stark.report.html|stark.report.pdf|tag)$#", $runs_repositories_samples_files);

	#echo "samples_path_pattern_brace=$samples_path_pattern_brace";
	#echo "<pre>"; print_r($samples_files);echo "</pre>"; 
	#echo "<br><br>";

	$reports=null;
	$tags_sample_files=null;
	$tags_analysis_files=null;

	#echo "<br><br>";
	foreach ($samples_files as $key=>$sample_file) {

		$sample_file=trim($sample_file);

		#echo "basename ".basename($sample_file)."<br>";
		#switch (basename($sample_file)) {
		switch (true) {
			case (preg_match('/stark.report.html$/', basename($sample_file), $matches, PREG_OFFSET_CAPTURE)):
				$reports[]=$sample_file;
				break;
			case (preg_match('/stark.report.pdf$/', basename($sample_file), $matches, PREG_OFFSET_CAPTURE)):
				$reports[]=$sample_file;
				break;
			case (preg_match('/analysis.tag$/', basename($sample_file), $matches, PREG_OFFSET_CAPTURE)):
				$tags_analysis_files[]=$sample_file;
				break;
			case (preg_match('/tag$/', basename($sample_file), $matches, PREG_OFFSET_CAPTURE)):
				$tags_sample_files[]=$sample_file;
				break;
			default:
				break;
		}
	};

	#echo "<pre>"; print_r($tags_sample_files);echo "</pre>"; 


	# RESULTS_SUBFOLDER_DATA


	foreach ($reports as $key => $report_html) {

		$report_html_split=explode("/",trim($report_html));
		$report_html_split_reverse=array_reverse($report_html_split);
	
		// echo "<pre>";
		// echo "<br>";
		// echo "count: ".count($report_html_split)."<br>";
		// print_r($report_html_split);
		// echo "<br>";
		// echo "count: ".count($report_html_split_reverse)."<br>";
		// print_r($report_html_split_reverse);
		// echo "</pre>";

		$report_html_file=$report_html_split_reverse[0];
		$report_html_file_split=explode ( "." , $report_html_file );
		$report_html_id=$report_html_file_split[1];

		$sample=$report_html_split_reverse[1];
		$run=$report_html_split_reverse[2];
		$project=$report_html_split_reverse[3];
		$group=$report_html_split_reverse[4];
	
		$repository=$report_html_split[1];
		$root=$report_html_split[0];

		$sample_path=dirname($report_html);
		$run_path=dirname($sample_path);
		$project_path=dirname($run_path);
		$group_path=dirname($project_path);
		$repository_path=$root.'/'.$repository;

		# TSV
		// $report_tsv="$sample_path/$sample.final.tsv";
		// $report_tsv_link="";
		// #if (isset($runs_repositories_samples_files[$report_tsv])) {
		// if (in_array($report_tsv,$runs_repositories_samples_files)) {
		// 		$report_tsv_link='<a href="'.$report_tsv.'" download title="'.basename($report_tsv).'">TSV</a>&nbsp;';
		// }
		$report_tsv_link="";
		$report_tsv_pattern="$sample_path/$sample.*.tsv";
		$report_tsv_array=preg_grep("#$report_tsv_pattern$#", $runs_repositories_samples_files);
		foreach ($report_tsv_array as $report_tsv) {
			$report_tsv_link.='<a href="'.$report_tsv.'" download title="'.basename($report_tsv).'">TSV</a>&nbsp;';
		}
		
		# VCF
		// $report_vcf_gz="$sample_path/$sample.final.vcf.gz";
		// $report_vcf_link="";
		// #if (isset($runs_repositories_samples_files[$report_vcf])) {
		// if (in_array($report_vcf_gz,$runs_repositories_samples_files)) {
		// 	$report_vcf_link='<a href="'.$report_vcf_gz.'" download title="'.basename($report_vcf_gz).'">VCF</a>&nbsp;';
		// }
		$report_vcf_link="";
		$report_vcf_pattern="$sample_path/$sample.*.vcf.gz";
		$report_vcf_array=preg_grep("#$report_vcf_pattern$#", $runs_repositories_samples_files);
		foreach ($report_vcf_array as $report_vcf) {
			$report_vcf_link.='<a href="'.$report_vcf.'" download title="'.basename($report_vcf).'">VCF</a>&nbsp;';
		}

		# DESIGN BED
		$report_bed_link="";
		$report_bed="$sample_path/$sample.bed";
		if (in_array($report_bed,$runs_repositories_samples_files)) {
				$report_bed_link.='<a href="'.$report_bed.'" download title="'.basename($report_bed).'">Design</a>&nbsp;';
		}
		$report_bed_pattern="$sample_path/$sample.*Design.bed";
		$report_bed_array=preg_grep("#$report_bed_pattern$#", $runs_repositories_samples_files);
		foreach ($report_bed_array as $report_bed) {
			$report_bed_link.='<a href="'.$report_bed.'" download title="'.basename($report_bed).'">Design</a>&nbsp;';
		}

		# PANEL BED
		$report_genes_link="";
		$report_genes="$sample_path/$sample.genes.bed";
		if (in_array($report_genes,$runs_repositories_samples_files)) {
				$report_genes_link.='<a href="'.$report_genes.'" download title="'.basename($report_genes).'">Panel</a>&nbsp;';
		}
		$report_genes="$sample_path/$sample.genes";
		if (in_array($report_genes,$runs_repositories_samples_files)) {
				$report_genes_link.='<a href="'.$report_genes.'" download title="'.basename($report_genes).'">Panel</a>&nbsp;';
		}
		$report_genes_pattern="$sample_path/$sample.*Panel.*genes.bed";
		$report_genes_array=preg_grep("#$report_genes_pattern$#", $runs_repositories_samples_files);
		foreach ($report_genes_array as $report_genes) {
			$report_genes_link.='<a href="'.$report_genes.'" download title="'.basename($report_genes).'">Panel</a>&nbsp;';
		};


		# Find infos
		// $report_html_split=explode ( "/" , $report_html );
		// $root=$report_html_split[0];
		// $repository=$report_html_split[1];
		// $group=$report_html_split[2];
		// $project=$report_html_split[3];
		// $run=$report_html_split[4];
		// $sample=$report_html_split[5];
		// $report_html_file=$report_html_split[6];
		// $report_html_file_split=explode ( "." , $report_html_file );
		// $report_html_id=$report_html_file_split[1];

		# Files
		// $sample_path=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample;
		// $run_path=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run;
		// $project_path=$root.'/'.$repository.'/'.$group.'/'.$project;
		// $group_path=$root.'/'.$repository.'/'.$group;
		// $repository_path=$root.'/'.$repository;
		// $report_tsv=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample.'/'.$sample.'.final.tsv';
		// $report_vcf_gz=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample.'/'.$sample.'.final.vcf.gz';
		// $report_bed=$root.'/'.$repository.'/'.$group.'/'.$project.'/'.$run.'/'.$sample.'/'.$sample.'.bed';

		# paths level
		$level_path["sample"]=$sample_path;
		$level_path["analysis"]=$run_path;
		$level_path["project"]=$project_path;
		$level_path["group"]=$group_path;
		$level_path["repository"]=$repository_path;


		# TAGS

		# Deprecated
		#$tags_file=glob ( "".$sample_path."/".$sample.".tag" )[0];
		#$analysis_tags_file=glob ( "".$sample_path."/".$sample.".analysis.tag" )[0];
		#print_r($_SESSION["data"]["file_content"]);

		$tags_file=implode("",preg_grep('#'.$sample_path."/".$sample.".tag".'#', $tags_sample_files));
		$analysis_tags_file=implode("",preg_grep('#'.$sample_path."/".$sample.".analysis.tag".'#', $tags_analysis_files));
		
		#echo "$tags_file <br>";
		#echo "read_file_content=$read_file_content <br>";
		#$read_file_content=1;
		$data_read_file=false;

		$tags="";
		$tags_analysis="";

		if ($read_file_content) {

			if ($tags_file != "") {
				if (!isset($_SESSION["data"]["file_content"][$tags_file])) {
					$myfile = fopen($tags_file, "r") or die("Unable to open file '$tags_file'!");
					$_SESSION["data"]["file_content"][$tags_file]=trim(fread($myfile,filesize($tags_file)));
					fclose($myfile);
				};
				$tags=str_replace(" ","&nbsp;&nbsp; ",tags_extract($_SESSION["data"]["file_content"][$tags_file]));
			};
			if ($analysis_tags_file != "") {
				if (!isset($_SESSION["data"]["file_content"][$analysis_tags_file])) {
					$myfile = fopen($analysis_tags_file, "r") or die("Unable to open file '$analysis_tags_file'!");
					$_SESSION["data"]["file_content"][$analysis_tags_file]=trim(fread($myfile,filesize($analysis_tags_file)));
					fclose($myfile);
				};
				$tags_analysis=str_replace(" ","&nbsp;&nbsp; ",tags_extract($_SESSION["data"]["file_content"][$analysis_tags_file]));
				
			};

		};

		if ($tags != "") {
			$tags="<br><small><span style='color:gray'>$tags</span></small>";
		};
		if ($tags_analysis != "") {
			$tags_analysis="<br><small><span style='color:gray'>$tags_analysis</span></small>";
		};

		#echo "<pre>"; print_r($plugins_obj); echo "</pre>";


		# PLUGINS
		$PLUGINS_LINKS=array();
		foreach ($plugins_obj as $plugin=>$plugin_infos) {

			if (is_file($plugins_folder."/".$plugin_infos->{"script"})						# script exists
				&& is_dir($sample_path."/".$plugin_infos->{"folder_data"})					# folder data exists
				&& glob($sample_path."/{".$plugin_infos->{"files_data"}."}",GLOB_BRACE)		# files data exists (at least one)
				&& $plugin_infos->{"available"}												# plugin available
				&& $modules_obj_array[$plugin_infos->{"module"}]->{"available"}				# module available
				#&& $modules_obj_array[$plugin_infos->{"module"}]->{"services"}->{$plugin_infos->{"service"}}->{"available"}		# service available
				&& $modules_obj_array[$plugin_infos->{"module"}]->{"submodules"}->{$plugin_infos->{"submodule"}}->{"services"}->{$plugin_infos->{"service"}}->{"available"}		# service available
			) {

				#$uri=$modules_obj_array[$plugin_infos->{"module"}]->{"services"}->{$plugin_infos->{"service"}}->{"href"};
				#if ($uri!="") {
					foreach ($plugin_infos->{"level"} as $level=>$level_infos) {
						$url_full=$plugins_folder."/".$plugin_infos->{"script"}."?".$plugin_infos->{"parameters"};
						$url_full.=($level_infos->{"parameters"}!="")?"&".$level_infos->{"parameters"}:"";
						$url_full.="&PATH[]=".$level_path[$level];
						$link="<a target='".$plugin_infos->{"target"}."' href='$url_full'>".$plugin_infos->{"label"}."</a>";
						$PLUGINS_LINKS[$level].=" ".$link."&nbsp; ";
					};
				#};

			}
			
		}

		# TBODY
		$tbody=$tbody.'
			<tr class="table-heads">
				<td class="head-item mbr-fonts-style display-7" >
					<a href="'.$report_html.'" title="Report '.$report_html_id.'">
						 <img src="'.$APP_LOGO.'" style="height: 2rem;">
					</a>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						<b title="'.$report_html_id.'">'.$sample.'</b>
						'.$tags.'
						<br>
						'.$report_tsv_link.'&nbsp;
						'.$report_vcf_link.'&nbsp;
						'.$report_bed_link.'&nbsp;
						'.$report_genes_link.'&nbsp;
						
						'.$PLUGINS_LINKS["sample"].'&nbsp;
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$run.'
						'.$tags_analysis.'
						<br>
						'.$PLUGINS_LINKS["analysis"].'
						'.$METRICS_RUN_LINK.'
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$project.'
						<br>
						'.$PLUGINS_LINKS["project"].'
						'.$METRICS_PROJECT_LINK.'
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$group.'
						<br>
						'.$PLUGINS_LINKS["group"].'
						'.$METRICS_GROUP_LINK.'
					</small>
				</td>
				<td class="head-item mbr-fonts-style display-7">
					<small>
						'.$repository.'
						<br>
						'.$PLUGINS_LINKS["repository"].'
					</small>
				</td>
			</tr>
		';
	};


	### SECTION
	#############

	$CONTENT_SECTION_REPORTS= '

	<section class="section-table cid-ruiBoanIwc" id="REPORTS">

	  <div class="container container-table">

		  <div class="table-wrapper">

			<div class="container">
			  <div class="row search">
				<div class="col-md-6"></div>
				<div class="col-md-6">
					<div class="dataTables_filter">
					  <label class="searchInfo mbr-fonts-style display-7">Search:</label>
					  <input class="form-control input-sm" disabled="">
					</div>
				</div>
			  </div>
			</div>

			<div class="container scroll">
				<table class="table isSearch" cellspacing="0">
					<thead>
						'.$thead.'
					</thead>
					<tbody>
						'.$tbody.'
					</tbody>
				</table>
			</div>

			<div class="container table-info-container">
			  <div class="row info">
				<div class="col-md-6">
				  <div class="dataTables_info mbr-fonts-style display-7">
					<span class="infoBefore">Showing</span>
					<span class="inactive infoRows"></span>
					<span class="infoAfter">entries</span>
					<span class="infoFilteredBefore">(filtered from</span>
					<span class="inactive infoRows"></span>
					<span class="infoFilteredAfter"> total entries)</span>
				  </div>
				</div>
				<div class="col-md-6"></div>
			  </div>
			</div>

		  </div>

		</div>

	</section>
	';

};




###############
### CONTENT ###
###############

$CONTENT=$CONTENT_SECTION_REPORTS;

echo $CONTENT;


##############
### FOOTER ###
##############

include "footer.inc.php";



# Script end time

if ($VERBOSE) {

	$ru = getrusage();
	echo "This process used " . rutime($ru, $rustart, "utime") .
	    " ms for its computations\n";
	echo "It spent " . rutime($ru, $rustart, "stime") .
	    " ms in system calls\n";

};




?>
