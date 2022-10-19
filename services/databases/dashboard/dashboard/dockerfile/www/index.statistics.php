<?php


#############
### INFOS ###
#############

$APP_SECTION="Statistics";



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

$count="  ";


### DATA
###########


include "global_statistics.inc.php";



### CONTENT
#############



### GLOBAL Statistics

$GLOBAL_HEADER='


	   <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
		   Global Statistics
	   </h2>

	   <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
		   Statistics for inputs and repositories folders (duplication not considered)
	   </p>
';


# GLOBAL Inputs

$GLOBAL_INPUT="";

foreach ($input_list as $input=>$hash) {

	$GLOBAL_INPUT.='

	   <!--
	   <h3 class="mbr-section-subtitle mbr-fonts-style align-center mbr-light display-2">
	   	'.$input.'
	   </h3>
	   -->

	   <div class="container pt-4 mt-2">
	   	<div class="media-container-row">

			<div class="card p-3 align-center col-12 col-md-6 col-lg-3">
	   			<div class="panel-item">
	   				<div class="card-img pb-3">
	   					<span class="mbr-iconfont mbri-info"></span>
	   				</div>

	   				<div class="card-text">
	   					<h4 class="mbr-section-subtitle align-center mbr-light pt-3 pb-3 mbr-fonts-style display-5">
	   						'.$input.'
	   					</h4>
	   				</div>
	   			</div>
	   		</div>


	   		<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
	   			<div class="panel-item">
	   				<div class="card-img pb-3">
	   					<span class="mbr-iconfont mbri-align-center"></span>
	   				</div>

	   				<div class="card-text">
	   					<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
	   						'.count($runs[$input]).'
	   					</h3>
	   					<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
	   						Runs
	   					</h4>
	   				</div>
	   			</div>
	   		</div>

			<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
	   			<div class="panel-item">
	   				<div class="card-img pb-3">
	   					<span class="mbr-iconfont mbri-target"></span>
	   				</div>

	   				<div class="card-text">
	   					<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
	   						'.count($designs[$input]).'
	   					</h3>
	   					<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
	   						Target Designs
	   					</h4>
	   				</div>
	   			</div>
	   		</div>

			<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
	   			<div class="panel-item">
	   				<div class="card-img pb-3">
	   					<span class="mbr-iconfont mbri-menu"></span>
	   				</div>

	   				<div class="card-text">
	   					<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
	   						'.count($genes[$input]).'
	   					</h3>
	   					<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
	   						Gene Panels
	   					</h4>
	   				</div>
	   			</div>
	   		</div>

			<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
	   			<div class="panel-item">
	   				<div class="card-img pb-3">
	   					<span class="mbr-iconfont mbri-bookmark"></span>
	   				</div>

	   				<div class="card-text">
	   					<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
	   						'.count($transcripts[$input]).'
	   					</h3>
	   					<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
	   						Transcript custom lists
	   					</h4>
	   				</div>
	   			</div>
	   		</div>

	   	</div>
	   </div>

	   <p class="mbr-content-text mbr-fonts-style display-7">
	       &nbsp;
	  </p>

	';

	# add-submenu align-justify bookmark
	# transcripts : bookmark
	# SAPLE : mbri-code

};

# GLOBAL Repositories

$GLOBAL_REPOSITORY="";

foreach ($global as $repository=>$hash) {

	$nb_groups=count($hash["groups"]);
	$nb_projects=count($hash["projects"]);
	$nb_runs=count($hash["runs"]);
	$nb_samples=count($hash["samples"]);

	$GLOBAL_REPOSITORY.='

			<!--
		   <h3 class="mbr-section-subtitle mbr-fonts-style align-center mbr-light display-2">
		   	'.$repository.'
		   </h3>
		   -->

		   <div class="container pt-4 mt-2">
		   	<div class="media-container-row">

		   		<div class="card p-3 align-center col-12 col-md-6 col-lg-3">
		   			<div class="panel-item">
		   				<div class="card-img pb-3">
		   					<span class="mbr-iconfont mbri-info"></span>
		   				</div>

		   				<div class="card-text">
		   					<h4 class="mbr-section-subtitle align-center mbr-light pt-3 pb-3 mbr-fonts-style display-5">
		   						'.$repository.'
		   					</h4>
		   				</div>
		   			</div>
		   		</div>

		   		<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
		   			<div class="panel-item">
		   				<div class="card-img pb-3">
		   					<span class="mbr-iconfont mbri-folder"></span>
		   				</div>

		   				<div class="card-text">
		   					<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
		   						'.$nb_groups.'
		   					</h3>
		   					<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
		   						Groups
		   					</h4>
		   				</div>
		   			</div>
		   		</div>

				<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
		   			<div class="panel-item">
		   				<div class="card-img pb-3">
		   					<span class="mbr-iconfont mbri-folder"></span>
		   				</div>

		   				<div class="card-text">
		   					<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
		   						'.$nb_projects.'
		   					</h3>
		   					<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
		   						Projects
		   					</h4>
		   				</div>
		   			</div>
		   		</div>

				<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
		   			<div class="panel-item">
		   				<div class="card-img pb-3">
		   					<span class="mbr-iconfont mbri-align-center"></span>
		   				</div>

		   				<div class="card-text">
		   					<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
		   						'.$nb_runs.'
		   					</h3>
		   					<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
		   						Analyses / Runs
		   					</h4>
		   				</div>
		   			</div>
		   		</div>

				<div class="card p-3 align-center col-12 col-md-6 col-lg-2">
		   			<div class="panel-item">
		   				<div class="card-img pb-3">
		   					<span class="mbr-iconfont mbri-code"></span>
		   				</div>

		   				<div class="card-text">
		   					<h3 class="'.$count.' pt-3 pb-3 mbr-fonts-style display-2">
		   						'.$nb_samples.'
		   					</h3>
		   					<h4 class="mbr-content-title mbr-bold mbr-fonts-style display-7">
		   						Samples
		   					</h4>
		   				</div>
		   			</div>
		   		</div>

		   	</div>
		   </div>

		   <p class="mbr-content-text mbr-fonts-style display-7">
		       &nbsp;
		  </p>

	';

};

# GLOBAL Statistics

$GLOBAL_STATISTICS='
	<section class="counters1 counters cid-ru7OEH6r7i" id="SECTION_ID" >

		<div class="container">

			'.$GLOBAL_HEADER.'
			'.$GLOBAL_INPUT.'
			'.$GLOBAL_REPOSITORY.'

		</div>
		
	</section>
';

echo $GLOBAL_STATISTICS;





# GROUP Statistics

$links_section="";

foreach ($global as $repository=>$hash) {

	# Analyses/Runs and Samples by group
	$run_by_group_nb=array();
	$sample_by_group_nb=array();
	$run_by_group_project_nb=array();
	$sample_by_group_project_nb=array();
	foreach ($run_by_group[$repository] as $group=>$hash2) {
		$run_by_group_nb[$group]=count($hash2);
		foreach ($run_by_group_project[$repository][$group] as $project=>$hash3) {
			$run_by_group_project_nb[$group][$project]=count($hash3);
		};
	};
	foreach ($sample_by_group[$repository] as $group=>$hash2) {
		$sample_by_group_nb[$group]=count($hash2);
		foreach ($sample_by_group_project[$repository][$group] as $project=>$hash3) {
			$sample_by_group_project_nb[$group][$project]=count($hash3);
		};
	};
	$group_list="'".implode("','",array_keys($run_by_group_nb))."'";
	$run_by_group_list=implode(",",$run_by_group_nb);
	$sample_by_group_list=implode(",",$sample_by_group_nb);


	$run_by_group_html = "
	<div class='section_div'>
		<div id='".$repository."run_by_group_list'>
			<div id='".$repository."run_by_group_list_figure'>
				<div class='figure' id='".$repository."plot_run_by_group_list' style='height:400px;width:500px;'></div>
			</div>

			<script type='text/javascript'>
				var data = [{
					  values: [".$run_by_group_list."],
					  labels: [".$group_list."],
					  type: 'pie'
					}];
				var layout={title:'Number of Analyses/Runs by group ', xaxis:{title:'duplication level'}, yaxis:{title:'Read percent (%) & GC ratio'}};
				Plotly.newPlot('".$repository."plot_run_by_group_list', data, layout);
			</script>

		</div>
	</div>
	";

	$sample_by_group_html = "
		<div class='section_div'>
			<div id='".$repository."sample_by_group_list'>
				<div id='".$repository."sample_by_group_list_figure'>
					<div class='figure' id='".$repository."plot_sample_by_group_list' style='height:400px;width:500px;'></div>
				</div>

				<script type='text/javascript'>
					var data = [{
						  values: [".$sample_by_group_list."],
						  labels: [".$group_list."],
						  type: 'pie'
						}];
					var layout={title:'Number of Samples by group', xaxis:{title:'duplication level'}, yaxis:{title:'Read percent (%) & GC ratio'}};
					Plotly.newPlot('".$repository."plot_sample_by_group_list', data, layout);
				</script>
			</div>
		</div>
	";

	$project_id=$repository."project_id";

	$project_link='
		<li class="nav-item mbr-fonts-style">
			<a class="nav-link active display-7" role="tab" data-toggle="tab" href="#tab'.$project_id.'">
				ALL Groups
			</a>
		</li>
	';

	$project_content='
		<div id="tab'.$project_id.'" class="tab-pane in active mbr-table" role="tabpanel">
			<div class="row tab-content-row">
				<div class="container pt-0 mt-0">
					<div class="media-container-row">
						'.$run_by_group_html.'
						'.$sample_by_group_html .'
					</div>
				</div>
			</div>
		</div>
	';


	foreach ($run_by_group_project_nb as $group=>$projects) {

		$project_id=str_replace(".","_",$repository.$group."project_id");

		$project_list="'".implode("','",array_keys($projects))."'";
		$run_by_group_project_list=implode(",",$projects);
		$sample_by_group_project_list=implode(",",$sample_by_group_project_nb[$group]);

		$run_by_group_html = "
			<div class='section_div'>
				<div id='".$project_id."run_by_group_list'>
					<div id='".$project_id."run_by_group_list_figure'>
						<div class='figure' id='".$project_id."plot_run_by_group_list' style='height:400px;width:500px;'></div>
					</div>

					<script type='text/javascript'>
						var data = [{
							  values: [".$run_by_group_project_list."],
							  labels: [".$project_list."],
							  type: 'pie'
							}];
						var layout={title:'Number of Analyses/Runs by projects', xaxis:{title:'duplication level'}, yaxis:{title:'Read percent (%) & GC ratio'}};
						Plotly.newPlot('".$project_id."plot_run_by_group_list', data, layout);
					</script>

				</div>
			</div>
		";

		$sample_by_group_html = "
			<div class='section_div'>
				<div id='".$project_id."sample_by_group_list'>
					<div id='".$project_id."sample_by_group_list_figure'>
						<div class='figure' id='".$project_id."plot_sample_by_group_list' style='height:400px;width:500px;'></div>
					</div>

					<script type='text/javascript'>
						var data = [{
							  values: [".$sample_by_group_project_list."],
							  labels: [".$project_list."],
							  type: 'pie'
							}];
						var layout={title:'Number of Samples by projects', xaxis:{title:'duplication level'}, yaxis:{title:'Read percent (%) & GC ratio'}};
						Plotly.newPlot('".$project_id."plot_sample_by_group_list', data, layout);
					</script>
				</div>
			</div>
		";


		# foreach group
		$project_link.='
			<li class="nav-item mbr-fonts-style">
				<a class="nav-link  display-7" role="tab" data-toggle="tab" href="#tab'.$project_id.'">
					'.$group.'
				</a>
			</li>
		';

		$project_content.='
			<div id="tab'.$project_id.'" class="tab-pane in mbr-table" role="tabpanel">
				<div class="row tab-content-row">

					<div class="container pt-0 mt-0">
			 		   	<div class="media-container-row">

							'.$run_by_group_html.'

							'.$sample_by_group_html .'

						</div>
					</div>


				</div>
			</div>
		';

	};

	$links_section.='

		<h3 class="mbr-section-subtitle mbr-fonts-style align-center mbr-light display-2">
		 '.$repository.'
		</h3>

		<div class="container-fluid col-md-12">
			<div class="row tabcont">
				<ul class="nav nav-tabs pt-0 mt-0" role="tablist">
					'.$project_link.'
				</ul>
			</div>
		</div>

		<div class="container">
			<div class="row px-1">
				<div class="tab-content">
					'.$project_content.'
				</div>
			</div>
		</div>

	';

};


$GLOBAL_STATISTICS='
<section class=" cid-rtQc0shhyd" id="truc">

	<div class="container">

	   <h2 class="mbr-section-title pb-3 align-center mbr-fonts-style display-2">
		   Groups Statistics
	   </h2>

	   <p class="mbr-section-subtitle mbr-fonts-style display-6 align-center">
		   Statistics on analyses/runs and samples for each projects in groups (duplication not considered)
	   </p>

	   '.$links_section.'

	</div>
</section>
';

echo $GLOBAL_STATISTICS;



##############
### FOOTER ###
##############

include "footer.inc.php";

?>
