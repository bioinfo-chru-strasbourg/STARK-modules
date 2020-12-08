<?php


### INCLUDE
#############






### Task Spooler List (as Array)
#################################

$TS_LIST=file($urinner_api_list);


### TASKS
foreach ($TS_LIST as $TASK_KEY => $ONE_TASK) {

	# SAMPLE NAME
	$ONE_TASK_split=explode( " ", preg_replace('/\s+/', ' ',($ONE_TASK)) ); # preg_replace('/\s+/', ' ',$row['message']);


	if ($TASK_KEY!=0) {

		$ONE_TASK_ID=$ONE_TASK_split[0];
		$ONE_TASK_STATE=$ONE_TASK_split[1];
		$ONE_TASK_OUTPUT=$ONE_TASK_split[2];
		if ($ONE_TASK_STATE=="running" || $ONE_TASK_STATE=="queued" ) {
			$ONE_TASK_E_LEVEL="";
			$ONE_TASK_TIMES="";
			$ONE_TASK_COMMAND=$ONE_TASK_split[3];
		} else {
			$ONE_TASK_E_LEVEL=$ONE_TASK_split[3];
			$ONE_TASK_TIMES=$ONE_TASK_split[4];
			$ONE_TASK_COMMAND=$ONE_TASK_split[5];
		}

		# INFO
		#echo $urinner_api_info."$ONE_TASK_ID"."<br>";

		$ONE_TASK_INFO_ARRAY=file($urinner_api_info."$ONE_TASK_ID");
		$ONE_TASK_INFO_CONTENT=implode("\n",$ONE_TASK_INFO_ARRAY);
		$ONE_TASK_INFO_CONTENT_HTML=implode("<br><br>",$ONE_TASK_INFO_ARRAY);

		if ($ONE_TASK_STATE=="running") {
			$ONE_TASK_TIME_RUNNING=str_replace('s','',explode(":",trim($ONE_TASK_INFO_ARRAY[4]))[1]);
		} else {
			$ONE_TASK_TIME_RUNNING="";
		}

		# LOG
		$ONE_TASK_LOG_ARRAY=file($urinner_api_log."$ONE_TASK_ID");
		$ONE_TASK_LOG_CONTENT=implode("",$ONE_TASK_LOG_ARRAY);

		$ONE_TASK_BASENAME=basename($ONE_TASK_OUTPUT);

		preg_match('#\[(.*?)\]#', $ONE_TASK_COMMAND, $match);
		$ONE_TASK_COMMAND_ID=$match[1];

		#$ONE_TASK_RUN=explode('.',$ONE_TASK_COMMAND_ID)[2];
		$ONE_TASK_RUN_TYPE=explode('.',$ONE_TASK_COMMAND_ID)[0];
		$ONE_TASK_RUN=explode('-',explode('.',$ONE_TASK_COMMAND_ID)[2])[3];
		preg_match("/(.*)\.(.*)\.ID-(.*)-NAME-(.*)/i", $ONE_TASK_COMMAND_ID,$ONE_TASK_COMMAND_ID_matches);
		#print_r($ONE_TASK_COMMAND_ID_matches);
		$ONE_TASK_RUN_TYPE=$ONE_TASK_COMMAND_ID_matches[1];
		$ONE_TASK_RUN=$ONE_TASK_COMMAND_ID_matches[4];

		if ($ONE_TASK_RUN=="") {
			preg_match("/(.*)\.(.*)\.(.*)/i", $ONE_TASK_COMMAND_ID,$ONE_TASK_COMMAND_ID_matches);
			$ONE_TASK_RUN_TYPE=$ONE_TASK_COMMAND_ID_matches[1];
			#$ONE_TASK_ANALYSIS_NAME=$ONE_TASK_COMMAND_ID_matches[3];
			$ONE_TASK_RUN=$ONE_TASK_COMMAND_ID_matches[3];
		}

		#echo "<br>$ONE_TASK_RUN_TYPE $ONE_TASK_RUN";


		if ($ONE_TASK_OUTPUT!="(file)") {
			$ONE_TASK_LOG_NAME="output";
			#$ONE_TASK_LOG_URL="http://$DATA_URL_EXTERNAL/$ONE_TASK_BASENAME";
			$ONE_TASK_LOG_URL="$folder_api/$ONE_TASK_BASENAME";
		} else {
			$ONE_TASK_LOG_NAME="";
			$ONE_TASK_LOG_URL="";
		};

		if ($ONE_TASK_E_LEVEL!="" && $ONE_TASK_E_LEVEL!="0" ) {
			$ONE_TASK_E_LEVEL_HTML="<br>Error '$ONE_TASK_E_LEVEL'";
		} else {
			$ONE_TASK_E_LEVEL_HTML="";
		};

		# Times
		$ONE_TASK_TIME="";
		if ($ONE_TASK_TIMES!="") {
			$ONE_TASK_TIME=$ONE_TASK_TIMES;
		} else if ($ONE_TASK_TIME_RUNNING!="") {
			$ONE_TASK_TIME=$ONE_TASK_TIME_RUNNING;
			#$ONE_TASK_TIME=str_replace(array("=","/",":","+"),array("d ","h ","m ","s "),gmdate("z=H/i:s+", explode("/",$ONE_TASK_TIME_RUNNING)[0]));
		}
		$ONE_TASK_TIME_FORMATED=str_replace(array("=","/",":","+"),array("d ","h ","m ","s "),gmdate("z=H/i:s+", explode("/",$ONE_TASK_TIME)[0]));

		# state finished/green queued/black running/orange other/red
		#echo "<BR>$ONE_TASK_STATE && $ONE_TASK_E_LEVEL";

		$ONE_TASK_STATE=($ONE_TASK_STATE=="")?"unavailable":$ONE_TASK_STATE;

		if ($ONE_TASK_STATE=="finished" && $ONE_TASK_E_LEVEL==0) {
			$ONE_TASK_STATE_COLOR="green";
			$ONE_TASK_TIME_COLOR="";
		} else if ($ONE_TASK_STATE=="finished" && $ONE_TASK_E_LEVEL!=0) {
			$ONE_TASK_STATE="error";
			$ONE_TASK_STATE_COLOR="red";
			$ONE_TASK_TIME_COLOR="red";
		} else if ($ONE_TASK_STATE=="queued") {
			$ONE_TASK_STATE_COLOR="gray";
			$ONE_TASK_TIME_COLOR="";
		} else if ($ONE_TASK_STATE=="running") {
			$ONE_TASK_STATE_COLOR="orange";
			$ONE_TASK_TIME_COLOR="orange";
		} else {
			$ONE_TASK_STATE_COLOR="red";
			$ONE_TASK_TIME_COLOR="";
		};

		#if (!isset($run_task[$ONE_TASK_RUN])) {
		if (!isset($run_task[$ONE_TASK_RUN]) || $ONE_TASK_STATE!="error" ) {

			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["status"]=$ONE_TASK_STATE;
			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["color"]=$ONE_TASK_STATE_COLOR;
			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["time"]=$ONE_TASK_TIME;
			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["time_formated"]=$ONE_TASK_TIME_FORMATED;
			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["time_color"]=$ONE_TASK_TIME_COLOR;
			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["error_level"]=$ONE_TASK_E_LEVEL;
			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["error_level_html"]=$ONE_TASK_E_LEVEL_HTML;
			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["id"]=$ONE_TASK_ID;
			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["log_url"]=$ONE_TASK_LOG_URL;
			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["log_name"]=$ONE_TASK_LOG_NAME;
			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["key"]=$TASK_KEY;
			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["type"]=$ONE_TASK_RUN_TYPE;
			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["command_id"]=$ONE_TASK_COMMAND_ID;
			$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN]["task"]["info_content_html"]=$ONE_TASK_INFO_CONTENT_HTML;


			# Last task (by ID)
			if ($ONE_TASK_RUN_TYPE=="STARK" || $TS_SHOW=="FULL") {
			
				if (!isset($run_task[$ONE_TASK_RUN]["task"]["id"]) || $ONE_TASK_ID>$run_task[$ONE_TASK_RUN]["task"]["id"] ) {

					$run_task[$ONE_TASK_RUN]=$full_task[$ONE_TASK_RUN_TYPE][$ONE_TASK_RUN];

				};

			};
		};

	}



};


############
### RUNS ###
############



# INPUTS
##########

$runs_infos=array();


# RUN folders
#array_multisort(array_map('filemtime', $inputs_runs), SORT_NUMERIC, SORT_DESC, $inputs_runs);

#print_r($runs);


foreach ($runs as $run_folder=>$runs_list) {
	foreach ($runs_list as $run_name=>$runs_number) {
		#echo "$run_folder=>$run_name<br>";
		#$run_path_split=explode("/",$run_path);
		#$run=$run_path_split[count($run_path_split)-1];
		$runs_infos[$run_name]["inputs"]["run_path"][$run_name]=$run_name;
	};
};

# Runs specific files
#print_r($runs);
$inputs_runs_pattern=null;
foreach ($runs as $input=>$runs_list) {
	foreach ($runs_list as $run=>$runs_nb) {
		#echo "$input/$run<br>";
		$inputs_runs_pattern[]="$input/runs/$run";
	};
};
$inputs_runs_pattern_brace="{".join(",",array_unique($inputs_runs_pattern))."}";
#echo "inputs_runs_pattern_brace=$inputs_runs_pattern_brace";

#$runs_inputs_files=glob($folder_inputs."/*/runs/*/{RTAComplete.txt,SampleSheet.csv}",GLOB_BRACE);
$runs_inputs_files=glob($folder_inputs."/$inputs_runs_pattern_brace/{RTAComplete.txt,SampleSheet.csv}",GLOB_BRACE);

#print_r($runs_inputs_files);

#print_r(preg_grep("#.*\[RTAComplete.txt\|SampleSheet.csv\]$#", $inputs_runs));
#print_r(array_merge(preg_grep("#RTAComplete.txt$#", $inputs_runs),preg_grep("#SampleSheet.csv$#", $inputs_runs)));
#print_r(preg_grep("#(RTAComplete.txt|SampleSheet.csv)$#", $inputs_runs));

$runs_inputs_files=preg_grep("#(RTAComplete.txt|SampleSheet.csv)$#", $inputs_runs);

#print_r($runs_inputs_files);

#array_multisort(array_map('filemtime', $runs_inputs_files), SORT_NUMERIC, SORT_DESC, $runs_inputs_files);


foreach ($runs_inputs_files as $runs_inputs_key=>$run_file_path) {

	$run_path=dirname(trim($run_file_path));
	$run_file_path_split=explode("/",trim($run_file_path));
	$run=$run_file_path_split[count($run_file_path_split)-2];
	$run_file=$run_file_path_split[count($run_file_path_split)-1];
	$runs_infos[$run]["inputs"]["run_files"][$run_file]["file"]=$run_file_path;

	#echo "$run_path | $run | $run_file <br>";

	#print_r($_SESSION["data"]["content_file"]);

	#echo $run_file;
	### READ FILE
	if ($data_read_file) {

		if ($run_file=="SampleSheet.csv") {

			#echo "SampleSheet read: $run_file_path<br>";

			# find INFOS within log file
			if (!isset($_SESSION["data"]["content_file"][trim($run_file_path)])) {
				$_SESSION["data"]["content_file"][trim($run_file_path)]=file(trim($run_file_path));
			}
			

			#$run_file_path_array=file(trim($run_file_path));
			$run_file_path_array=$_SESSION["data"]["content_file"][trim($run_file_path)];

			#print_r($run_file_path_array);

			# Explode listener file infos
			#$listener_infos=array();
			$ini_file_content="";
			foreach ($run_file_path_array as $run_file_path_array_key=>$run_file_path_array_info) {

				if (substr(trim($run_file_path_array_info),0,1)!="[") {
					$ini_file_content.=rand(1,10000)."=\"".trim($run_file_path_array_info)."\"\n";
				} else {
					$ini_file_content.=trim($run_file_path_array_info)."\n";
				}

				# ADD
				#$listener_infos[$runs_listener_log_file_array_info_matches[1]]=$runs_listener_log_file_array_info_matches[2];
			};
			#echo "<pre>$ini_file_content</pre>";
			$ini_filename="/tmp/".rand(1,100000);
			$ini_file = fopen($ini_filename, "w") or die("Unable to open file!");
			fwrite($ini_file, $ini_file_content);
			fclose($ini_file);
			$ini_array = parse_ini_file($ini_filename, true);
			unlink($ini_file);

			$runs_infos[$run]["inputs"]["run_files"][$run_file]["file_array"]=$ini_array;
			#echo "<pre>";
			#print_r($ini_array["Data"]);

			foreach ($ini_array["Data"] as $ini_array_key=>$ini_array_sample) {
				$ini_array_sample_split=explode(",",$ini_array_sample);
				#print_r($ini_array_sample_split);
				if ($ini_array_sample_split[0] != "Sample_ID") {
					if ($ini_array_sample_split[1]!="") {
						$runs_infos[$run]["inputs"]["run_files"][$run_file]["samples"][$ini_array_sample_split[1]]=$ini_array_sample_split[1];
						$runs_infos[$run]["inputs"]["run_files"][$run_file]["samples_infos"][$ini_array_sample_split[1]]=$ini_array_sample_split;
					} else if ($ini_array_sample_split[0]!="") {
						$runs_infos[$run]["inputs"]["run_files"][$run_file]["samples"][$ini_array_sample_split[0]]=$ini_array_sample_split[0];
						$runs_infos[$run]["inputs"]["run_files"][$run_file]["samples_infos"][$ini_array_sample_split[0]]=$ini_array_sample_split;
					}
				}
			}
			#echo "</pre>";
		}

	}

};

#echo "<pre>"; print_r($runs_infos); echo "</pre>";

# REPOSITORIES

#print_r($runs_repositories);


#$runs_repositories=glob($folder_repositories."/*/*/*/*/*",GLOB_ONLYDIR);
#array_multisort(array_map('filemtime', $runs_repositories), SORT_NUMERIC, SORT_DESC, $runs_repositories);
foreach ($runs_repositories as $runs_inputs_key=>$sample_path) {
	$sample_path_split=explode("/",$sample_path);
	$sample=basename(dirname($sample_path));
	#$run=$sample_path_split[count($sample_path_split)-2];
	$run_path=dirname(dirname($sample_path));
	$run=basename($run_path);
	#echo "<br>run:$run - sample=$sample";
	$runs_infos[$run]["repositories"]["run_path"][$run_path]=$run_path;
	$runs_infos[$run]["repositories"]["samples"][$sample]=$sample;
};


# LISTENER LOG
$LISTENER_LOG_PATTERN="{ID-*-NAME-*.log}";
#$runs_listener_log=glob("analyses/stark-services/listener/$LISTENER_LOG_PATTERN",GLOB_BRACE);
#$runs_listener_log=glob($folder_services."/".$STARK_LISTENER_SERVICE_FOLDER."/$LISTENER_LOG_PATTERN",GLOB_BRACE);

#echo "folder_listener=$folder_listener<br>";

$runs_listener_log=glob($folder_listener."/$LISTENER_LOG_PATTERN",GLOB_BRACE);
array_multisort(array_map('filemtime', $runs_listener_log), SORT_NUMERIC, SORT_DESC, $runs_listener_log);
foreach ($runs_listener_log as $runs_listener_log_key=>$runs_listener_log_file) {
	$runs_listener_log_file_split=explode("/",$runs_listener_log_file);
	$run_listener_log=$runs_listener_log_file_split[count($runs_listener_log_file_split)-1];

	#print_r("$runs_listener_log_key=>$runs_listener_log_file<br>");

	if (1) {

		# find INFOS within log file
		$runs_listener_log_file_array=file($runs_listener_log_file);

		# Explode listener file infos
		$listener_infos=array();
		foreach ($runs_listener_log_file_array as $runs_listener_log_file_array_key=>$runs_listener_log_file_array_info) {
			#echo "<br>$runs_listener_log_file_array_info";
			preg_match("/(.*): (.*)/i", $runs_listener_log_file_array_info,$runs_listener_log_file_array_info_matches);
			#print_r($runs_listener_log_file_array_info_matches);
			# ADD
			$listener_infos[$runs_listener_log_file_array_info_matches[1]]=$runs_listener_log_file_array_info_matches[2];
		};

	}

	# RUN
	$run=$listener_infos["RUN"];

	# ADD
	$runs_infos[$run]["listener"][$run_listener_log]["log"]=$run_listener_log;
	$runs_infos[$run]["listener"][$run_listener_log]["infos"]=$listener_infos;

};

# LAUNCHER LOG
#analyses/stark-services/igv/20191204-080425.432691.json
#$launcher_log=glob ( "igvjs/static/data/public/stark-services/launcher/*", GLOB_ONLYDIR );
#$LAUNCHER_LOG_PATTERN="{*.json,*.log,*.err,*.info,ts-out.*}";
#$LAUNCHER_LOG_PATTERN="{*.json,*.log,*.err,*.info,*.output,ts-out.*}";
$LAUNCHER_LOG_PATTERN="{*.json,*.log,*.err,*.info,*.output}";
#$runs_launcher_log=glob("analyses/stark-services/launcher/$LAUNCHER_LOG_PATTERN",GLOB_BRACE);

#$runs_launcher_log=glob($folder_services."/".$STARK_API_SERVICE_FOLDER."/$LAUNCHER_LOG_PATTERN",GLOB_BRACE);
$runs_launcher_log=glob($folder_api."/$LAUNCHER_LOG_PATTERN",GLOB_BRACE);
#echo "LAUNCHER API "; print_r($runs_launcher_log);
array_multisort(array_map('filemtime', $runs_launcher_log), SORT_NUMERIC, SORT_DESC, $runs_launcher_log);
foreach ($runs_launcher_log as $runs_launcher_log_key=>$runs_launcher_log_file) {
	#echo "<br><br>";
	#print_r($runs_launcher_log_file);
	$runs_listener_log_file_basename=basename($runs_launcher_log_file);
	preg_match("/(.*)\.(.*)\.ID-(.*)-NAME-(.*)\.(.*)/i", $runs_listener_log_file_basename,$runs_launcher_log_file_matches);
	$run_listener_log=$runs_listener_log_file_basename;
	$ext=$runs_launcher_log_file_matches[5];
	$run=$runs_launcher_log_file_matches[4];

	#echo "run=$run<br> ";
	#print_r($runs_list);
	#if ($run!="") {
	if ($run!="") {

		# Full list
		$runs_infos[$run]["launcher_list"][$ext][]=$runs_launcher_log_file;

		# Final output (test exit code)
		$run_launcher_info_split=file($runs_launcher_log_file);
		foreach ($run_launcher_info_split as $key => $value) {
			preg_match("/Exit status: died with exit code (.*)/i", $value,$run_launcher_info_split_matches);
			if (isset($run_launcher_info_split_matches[1])) {
				# If exit without error
				if ($run_launcher_info_split_matches[1]==0) {
					$runs_infos[$run]["launcher"][$ext]=$runs_launcher_log_file;
				# If exit with error and no other output
				} elseif (!isset($runs_infos[$run]["launcher"][$ext])) {
					$runs_infos[$run]["launcher"][$ext]=$runs_launcher_log_file;
				}
			}
		};
		#$runs_infos[$run]["launcher"][$ext]=$runs_launcher_log_file;

	};

};




# TABLE
# RUN / Sequencing / Analysis / Results


$tbody="";

foreach ($runs_infos as $run=>$run_infos) {

	if (isset($total_runs_list[$run])) {

		if ($DEBUG) {
			echo "<br><br>";
			echo "run: $run";
			echo "<pre>";
			print_r($run_infos);
			echo "</pre>";
		};

		$run_progress[$run]=array();

		$run_status="unknown";

		# Sequencing
		############

		$sequencing_status_code=isset($run_infos["inputs"])+isset($run_infos["inputs"]["run_files"]["RTAComplete.txt"]);

		switch ($sequencing_status_code) {
		case 0:
			$sequencing_status="unavailable";
			$sequencing_color="gray";
			$sequencing_message="No run sequencing folder";
			break;
		case 1:
			$sequencing_status="in progress";
			$sequencing_color="orange";
			$sequencing_message="Run sequencing seems to be in progress";
			break;
		case 2:
			#$sequencing_status="complete";
			if (!isset($run_infos["inputs"]["run_files"]["SampleSheet.csv"])) {
				$sequencing_status="ready";
				$sequencing_color="dodgerblue";
				$sequencing_message="Ready for analysis (SampleSheet missing)";
			} else {
				$sequencing_status="complete";
				$sequencing_color="green";
				$sequencing_message="Complete for analysis ";
			}
			break;
		default:
			$sequencing_status="unknown";
			$sequencing_color="gray";
			$sequencing_message="No information";
			break;
		}
		
		$sequencing_message_plus="";
		if (isset($runs_infos[$run]["inputs"]["run_files"]["RTAComplete.txt"])) {
			if ($data_read_file) {
				if (!isset($_SESSION["data"]["content_file"][trim($runs_infos[$run]["inputs"]["run_files"]["RTAComplete.txt"]["file"])])) {
					$_SESSION["data"]["content_file"][trim($runs_infos[$run]["inputs"]["run_files"]["RTAComplete.txt"]["file"])]=implode(file(trim($runs_infos[$run]["inputs"]["run_files"]["RTAComplete.txt"]["file"])));
				}
				$sequencing_message_plus.="<br>".$_SESSION["data"]["content_file"][trim($runs_infos[$run]["inputs"]["run_files"]["RTAComplete.txt"]["file"])]."<br>";
			}
			
		}
		if (count($runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]["samples"])) {
			$sequencing_message_plus.="<br><a target='SampleSheet' href='".$runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]["file"]."'>SampleSheet</a>";
			$sequencing_message_plus.="<br><span style='color:gray'>".implode("&nbsp;&nbsp;&nbsp; ",array_keys($runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]["samples"]))."</span>";
		}

		// echo "<pre>";
		// print_r($runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]);
		// echo "</pre>";


		$run_progress[$run]["sequencing"]["status"]=$sequencing_status;
		$run_progress[$run]["sequencing"]["color"]=$sequencing_color;
		$run_progress[$run]["sequencing"]["message"]=$sequencing_message;
		$run_progress[$run]["sequencing"]["message_plus"]=$sequencing_message_plus;

		if ($sequencing_status!="unavailable") {
			$run_status="Sequencing ".$sequencing_status;
			$run_color=$sequencing_color;
		}



		# Analysis
		#############

		// echo "<pre>";
		// echo "RUN: $run<br>";
		// print_r($run_infos);
		// echo "</pre>";

		// echo "RUN: $run<br>";
		// echo "<pre>";
		// print_r($run_infos["launcher"]);

		$analysis_infos=array();
		if (isset($run_task[$run]["task"])) {
			$analysis_status=$run_task[$run]["task"]["status"];
			$analysis_color=$run_task[$run]["task"]["color"];
			$analysis_message="";
			if ($run_task[$run]["task"]["error_level"]!=0) {
				$analysis_color="red";
				$run_status="error";
				$analysis_message.="<span style='color:red'>Error '".$run_task[$run]["task"]["error_level"]."' </span><br>";
			}
			$analysis_message.="<span style='color:".$run_task[$run]["task"]["time_color"]."'>".$run_task[$run]["task"]["time"]."</span>";

			# ANALYSIS INFOS
			$analysis_infos["error_level"]=$run_task[$run]["task"]["error_level"];
			$analysis_infos["status"]=$run_task[$run]["task"]["status"];
			$analysis_infos["color"]=$run_task[$run]["task"]["color"];
			$analysis_infos["time"]=$run_task[$run]["task"]["time"];
			#$analysis_infos["time_color"]=$run_task[$run]["task"]["time_color"];
			$analysis_infos["time_formated"]=$run_task[$run]["task"]["time_formated"];
			#$analysis_infos["time_color"]=$run_task[$run]["task"]["time_color"];

			#echo "<br><br>$run_task [$run][task]{status];".$analysis_infos["status"]."<br><br>";
			

		} elseif (isset($run_infos["launcher"]["info"])) {

				$run_launcher_info_split=file($run_infos["launcher"]["info"]);

				foreach ($run_launcher_info_split as $key => $value) {
					# ERROR_LEVEL
					preg_match("/Exit status: died with exit code (.*)/i", $value,$run_launcher_info_split_matches);
					#print_r($run_launcher_info_split_matches);
					if (isset($run_launcher_info_split_matches[1])) {
						$analysis_infos["error_level"]=$run_launcher_info_split_matches[1];
						if ($run_launcher_info_split_matches[1]!=0) {
							$analysis_infos["status"]="error";
							$analysis_infos["color"]="red";
						} else {
							$analysis_infos["status"]="finished";
							$analysis_infos["color"]="green";
						}
					}
					# TIME
					preg_match("/Time run: (.*)s/i", $value,$run_launcher_info_split_matches);
					#print_r($run_launcher_info_split_matches);
					if (isset($run_launcher_info_split_matches[1])) {
						$analysis_infos["time"]=$run_launcher_info_split_matches[1];
						$analysis_infos["time_formated"]=str_replace(array("=","/",":","+"),array("d ","h ","m ","s "),gmdate("z=H/i:s+", explode("/",$run_launcher_info_split_matches[1])[0]));
					}
				}

		} else {
			# TODO - Check on log file !!!
			$analysis_status="no information";
			$analysis_color="red";
			$analysis_message="Analysis launched but no Task Spooler information";
		};


		# Check analysis status code

		$analysis_status_code=isset($run_infos["listener"])+isset($run_infos["launcher"]);
		
		switch ($analysis_status_code) {
		case 0:
			$analysis_status="unavailable";
			$analysis_color="gray";
			$analysis_message="No analysis information available";
			break;
		case 1:
			#$analysis_status="detected";
			#$analysis_color="orange";
			if (isset($analysis_infos["status"])) {
				$analysis_status=$analysis_infos["status"];
			} else {
				$analysis_status="detected";
			}
			if (isset($analysis_infos["color"])) {
				$analysis_color=$analysis_infos["color"];
			} else {
				$analysis_color="orange";
			}
			$analysis_message="";
			if ( isset($analysis_infos["time_formated"]) ) {
				$analysis_message.="<span style='color:".$analysis_infos["time_color"]."'>".$analysis_infos["time_formated"];
			};

			#$analysis_message="Run sequencing ready detected";
			break;
		case 2:
			$analysis_status=$analysis_infos["status"];
			$analysis_color=$analysis_infos["color"];
			#$analysis_message="Run analysis launched";
			$analysis_message="";
			if ( isset($analysis_infos["time_formated"]) ) {
				$analysis_message.="<span style='color:".$analysis_infos["time_color"]."'>".$analysis_infos["time_formated"];
			};
			break;
		default:
			$analysis_status="unavailable";
			$analysis_color="gray";
			$analysis_message="No analysis information available";
			break;
		};

		# default
		$analysis_status=($analysis_status=="")?"unavailable":$analysis_status;
		$analysis_color=($analysis_color=="")?"gray":$analysis_color;
		$analysis_message=($analysis_message=="")?"No analysis information available":$analysis_message;
		$analysis_message_plus=($analysis_message_plus=="")?"":$analysis_message_plus;

		$analysis_message_plus="";
		// if (isset($runs_infos[$run]["launcher"]["output"])) {
		// 	$analysis_message_plus.="<br><br><a target='Output' href='".$runs_infos[$run]["launcher"]["output"]."' download>Output</a>";
		// };
		if (isset($runs_infos[$run]["launcher_list"]["output"])) {
			$analysis_message_plus.="";
			foreach ($runs_infos[$run]["launcher_list"]["output"] as $output_nb=>$output_file) {
				#print_r(file($runs_infos[$run]["launcher_list"]["info"][$output_nb]));

				# Final output (test exit code)
				$output_exit_msg="";
				$run_launcher_info_split=file($runs_infos[$run]["launcher_list"]["info"][$output_nb]);
				foreach ($run_launcher_info_split as $key => $value) {
					preg_match("/Exit status: died with exit code (.*)/i", $value,$run_launcher_info_split_matches);
					if (isset($run_launcher_info_split_matches[1])) {
						# If exit without error
						if ($run_launcher_info_split_matches[1]==0) {
							$output_exit_msg=" [<span style='color:green'>finished</span>]";
						# If exit with error and no other output
						} else {
							$output_exit_msg=" [<span style='color:red'>error</span>]";
						}
					}
				};


				$analysis_message_plus="<br><a target='Output' href='".$output_file."' download>Output".$output_exit_msg."</a>".$analysis_message_plus;
			};
			#$analysis_message_plus.="<br><br><a target='Output' href='".$runs_infos[$run]["launcher"]["output"]."' download>Output</a>";
		};

		$run_progress[$run]["analysis"]["status"]=$analysis_status;
		$run_progress[$run]["analysis"]["color"]=$analysis_color;
		$run_progress[$run]["analysis"]["message"]=$analysis_message;
		$run_progress[$run]["analysis"]["message_plus"]=$analysis_message_plus;

		if ($analysis_status!="unavailable") {
			$run_status="Analysis ".$analysis_status;
			$run_color=$analysis_color;
		};


		# Repository
		#############

		$repository_status_code=count($run_infos["repositories"]["run_path"]);
		if ($repository_status_code==0) {

			$repository_status="unavailable";
			$repository_color="gray";
			$repository_message="No results available";
			$repository_message_plus="";


		} else {

			$repository_status="available";
			$repository_color="green";
			$repository_message="Run results available in folders";
			$repository_message_plus="";

			# Check number of samples !
			if (isset($runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]) && isset($runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]["samples"])) {
				#$runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]["samples"]["truc"]="truc";
				$result1=array_diff($runs_infos[$run]["repositories"]["samples"],$runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]["samples"]);
				$result_missing=array_diff($runs_infos[$run]["inputs"]["run_files"]["SampleSheet.csv"]["samples"],$runs_infos[$run]["repositories"]["samples"]);
				$result=array_merge($result1,$result_missing);
				if (count($result_missing)) {
					if ($run_progress[$run]["analysis"]["color"]=="green") {
						$missing_samples_color="red";
					} else {
						$missing_samples_color=$run_task[$run]["task"]["color"];
					}
					$repository_message_plus.="<br><span style='color:$missing_samples_color'>Missing Samples (".implode(" ",$result_missing).")</span>";
				}
			}

			foreach ($run_infos["repositories"]["run_path"] as $repository_path=>$repository_path_infos) {
				# http://localhost:42002/index.reports.php?PATH=repositories/Archives/SOMATIC/HEMATOLOGY/RUN_TEST_TAG/
				#$REPORTS_RUN_LINK="<a target='METRICS' href='index.reports.php?PATH=$run_path'>REPORT</a>";
				$repository_path_split=explode("/",$repository_path);
				$repository_message_plus.="<br><a target='REPORT' href='index.reports.php?PATH=$repository_path'>".$repository_path_split[1]."/".$repository_path_split[2]."/".$repository_path_split[3]."</a>";
			}

		};

		$run_progress[$run]["repository"]["status"]=$repository_status;
		$run_progress[$run]["repository"]["color"]=$repository_color;
		$run_progress[$run]["repository"]["message"]=$repository_message;
		$run_progress[$run]["repository"]["message_plus"]=$repository_message_plus;

		if ($repository_status!="unavailable") {
			$run_status="Results ".$repository_status;
			$run_color=$repository_color;
		}

	}

}

// echo "<pre>";
// print_r($run_task);
// echo "</pre>";

// echo "<pre>";
// print_r($full_task);
// echo "</pre>";

// echo "<pre>";
// print_r($run_progress);
// echo "</pre>";



?>