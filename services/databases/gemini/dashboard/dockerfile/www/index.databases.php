<?php


#############
### INFOS ###
#############

$APP_SECTION="Database";



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

$db=$_REQUEST["database"];



### DATA
##########

include "databases.inc.php";




### CONTENT
#############



#########################
### SECTION DATABASES ###
#########################


# Init

$CONTENT_SECTION_DATABASES_LIST='';
$CONTENT_SECTION_DATABASE_CONTENT='';

$DATABASE_NUMBER=0;

// echo "<pre>";
// print_r($databases);
// echo "</pre>";


### Check databases

#$databases=glob("/gemini/db/*db");
// echo "db=$db ";
#$databases=glob("databases/gemini/current/db/*db");

$databases=[];
$databases[0]=$db;

// print_r($databases);
// echo "test";

foreach ($databases as $database_key=>$database_path) {


	$database_basename=basename($database_path);

	$database_basename_split=explode(".",$database_basename);
	$database_group=$database_basename_split[0];
	$database_project=$database_basename_split[1];


	if (0) {
		// $CONTENT_SECTION_DATABASES_LIST.='
		// 	<p class="mbr-text mbr-fonts-style display-7">
		// 		<big><a href="?database='.$database_info_code.'">'.$database_info_name.'</a></big> ['.$database_info_release.'] '.$database_info_fullname.'
		// 		<br>
		// 		'.$database_info_description.'
		// 	</p>
		// ';
		$CONTENT_SECTION_DATABASES_LIST.='
		<li class="nav-item mbr-fonts-style">
			<a class="nav-link  display-7" role="tab" href="?database='.$database_path.'">
				'.$database_group.'<br>'.$database_project.'
			</a>
		</li>
		';
		
	}

	# If database input
	if ($database_path==$db) {


		# NB variants
		$command="gemini query -q 'select count(variant_id) from variants' $database_path";
		$output = shell_exec($command);
		$database_nb_variants=trim($output);
		
		# NB samples
		$command="gemini query -q 'select count(sample_id) from samples' $database_path";
		$output = shell_exec($command);
		$database_nb_samples=trim($output);

		# NB genes
		$command="gemini query -q 'select count(gene) from variants where gene is not null' $database_path";
		$output = shell_exec($command);
		$database_nb_genes=trim($output);

		# clinvar_disease_name
		// $command="gemini query -q 'select clinvar_sig, count(variant_id) from variants group by clinvar_disease_name' $database_path";
		// $output = shell_exec($command);
		// $output_split=explode("\n",trim($output));
		// $clinvar_disease_name_html="";
		// foreach ($output_split as $output_line) {
		// 	$output_line_split=explode("\t",trim($output_line));
		// 	$clinvar_disease_name_html.="&nbsp;&nbsp;&nbsp; ".$output_line_split[1]." ".$output_line_split[0]."<br>";
		// };


		$query=$_REQUEST["query"];
		# query 
		if ($query!="") {
			$command="gemini query -q '".str_replace("\n","",$query)."' --header $database_path";
			$output = shell_exec($command);



		# clinvar_disease_name
		// $command="gemini query -q 'select clinvar_sig, count(variant_id) from variants group by clinvar_disease_name' $database_path";
		#$command="gemini query -q 'select chrom, start, end, ref, alt, gene from variants' --header $database_path";
		#$command="gemini query -q '$query' --header $database_path";
		// $command="gemini query -q 'select * from variants' --header $database_path";
		$output = shell_exec($command);
		$output_split=explode("\n",trim($output));
		#print_r($output_split);
		#$database_nb_samples=trim($output);
		

		$thead_content="";
		$tbody_content="";
		foreach ($output_split as $output_key=>$output_line) {
			#echo "output_line=$output_line<br>";
			$output_line_split=explode("\t",trim($output_line));
			#echo "output_line=$output_line";
			if ($output_key == 0) {
				#$thead.="&nbsp;&nbsp;&nbsp; ".$output_line_split[1]." ".$output_line_split[0]."<br>";
				#$thead_content.=str_replace("\t","</th><th class='head-item mbr-fonts-style display-7'>",$output_line);
				$thead_content.=implode("</th><th class='head-item mbr-fonts-style display-7'>",$output_line_split);
			} else {
				$tbody_content=$tbody_content.'
					<tr class="table-heads">
						<td class="head-item mbr-fonts-style display-7">
						'.implode("</td><td class='head-item mbr-fonts-style display-7'>",$output_line_split).'
						</td>
					</tr>
				';
				# '.str_replace("\t","</td><td class='head-item mbr-fonts-style display-7'>",$output_line).'
			};
		};
		$tbody=$tbody_content;


		$thead='
		<tr class="table-heads">
			<th class="head-item mbr-fonts-style display-7">
			'.$thead_content.'
			</th>
		</tr>
		';

	}
// foreach ($_SERVER as $k=>$i) {
// 	echo "$k=>$i<br>";
// }

		$CONTENT_SECTION_DATABASE_CONTENT.='


				<div class="row justify-content-left">

					<div class="card p-3 col-12 col-md-4 mb-3">
						<a href="index.databases.php?database='.$database_path.'" class="navbar-caption text-secondary ">
							<div class="media mb-0">
								<div class="card-img align-self-center">
									<span class="mbr-iconfont mbri-extension" style="color: rgb(20, 157, 204); fill: rgb(20, 157, 204);"></span>
								</div>
								<h4 class="card-title media-body py-3 mbr-fonts-style display-5">
								'.$database_group.'
								<br>
								'.$database_project.'
								</h4>
							</div>
						</a>

						<div class="card-box">
							<br>
							<p class="mbr-text mbr-fonts-style display-7">
							<b>'.$database_basename.'</b>
							<br>
							'.$database_nb_variants.' variants
							from
							'.$database_nb_genes.' genes
							within 
							'.$database_nb_samples.' samples

							<br>
							ClinVar Disease
							<br>
							'.$clinvar_disease_name_html.'
							</p>
							
						</div>

					</div>
					<div class="card p-3 col-12 col-md-8 mb-3">
						<a href="index.databases.php?database='.$module_info_code.'" class="navbar-caption text-secondary ">
							<div class="media mb-0">
								<div class="card-img align-self-center">
									<span class="mbr-iconfont mbri-extension" style="color: rgb(20, 157, 204); fill: rgb(20, 157, 204);"></span>
								</div>
								<h4 class="card-title media-body py-3 mbr-fonts-style display-5">
									Query
								</h4>
							</div>
						</a>

						<div class="card-box">
							<br>
							<p class="mbr-text mbr-fonts-style display-7">
							<form action="">
								<input id="database" name="database" type="hidden" value="'.$database_path.'">
								Example:
								<br>
								"select chrom, start, end, ref, alt from variants"
								<br>
								"select chrom, start, end, ref, alt, (gts).(*) from variants"
								<br>
								"select type, count(type) from variants group by type"
								<textarea id="query" name="query" rows="3" cols="80">'.$query.'</textarea>
								<br>
								<input id="Sumbit" name="Submit" type="submit">
							</form>
							</p>
							
						</div>

					</div>

				</div>
		';

	};


};




$CONTENT_SECTION_DATABASES_CONTENT='
	

		<!--
		<section class="features10 cid-rtQc0shhyd" id="features10-1l">

		<div class="container ">

			<h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
			   STARK Gemini database
		   </h2>
			
			<div class="card p-3 col-12 col-md-12 mb-0">

				<div class="container-fluid col-md-12">
					<div class="row tabcont">
						<ul class="nav nav-tabs pt-0 mt-0" role="tablist">
							'.$CONTENT_SECTION_DATABASES_LIST.' 
						</ul>
					</div>
				</div>

			</div>
			
		</div>
		
	</section>
	-->


	<section class=" cid-ru7OEDbxhA" id="truc"> 

		<div class="container">

			<h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
				STARK Gemini database
			</h2>			

				'.$CONTENT_SECTION_DATABASE_CONTENT.'


		</div>
		
	</section>

		

	';




### SECTION
#############

if ($query!="") {

	$CONTENT_SECTION_QUERY= '

	<section class="section-table cid-ruiBoanIwc" id="QUERY">

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

$CONTENT_SECTION_DATABASES=$CONTENT_SECTION_DATABASES_CONTENT.$CONTENT_SECTION_QUERY;




###############
### CONTENT ###
###############

$CONTENT=$CONTENT_SECTION_DASHBOARD.$CONTENT_SECTION_DATABASES.$CONTENT_OLD;

echo $CONTENT;


##############
### FOOTER ###
##############

include "footer.inc.php";



?>
