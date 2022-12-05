<?php
error_reporting(E_ALL);
session_start();


if (session_status() == PHP_SESSION_NONE || !isset($_SESSION['user_name'])) {
	header("HTTP/1.1 401 Unauthorized");
	#print("<p>Your not login.</p>");
} else {
	// default to loading the slurm user token
	if (!isset($_SESSION['user_token']) || $_SESSION['user_token'] == "") {
		$_SESSION['user_token'] = file_get_contents("/auth/slurm");
	}

	header("X-SLURM-USER-NAME: ".$_SESSION['user_name']);
	header("X-SLURM-USER-TOKEN: ".$_SESSION['user_token']);
}


# REQUESTS

$OUTPUT_FORMAT = $_GET['format'];


# Connexion params
$user=$_SESSION['user_name']; #$user="root";

if ($user == "") {
	$user = "root";
};

$token=$_SESSION['user_token']; #$token="eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJleHAiOjE2ODAxNjY2ODksImlhdCI6MTY3MDE2NjY5MCwic3VuIjoic2x1cm0ifQ.r5LV-z4EWjvHDJUgExQS1EHdMx-nPXvNowGmXFuSAJg";
$rest_hostname=trim(file_get_contents("/auth/rest_hostname")); #$rest_url="rest";
$openapi_release=trim(file_get_contents("/auth/openapi_release")); #$openapi_release="v0.0.38";
// # curl -H "X-SLURM-USER-NAME:$user" -H "X-SLURM-USER-TOKEN:$token" http://$rest_url/slurmdb/$openapi_release/qos
// # curl -H "X-SLURM-USER-NAME:$user" -H "X-SLURM-USER-TOKEN:$token" http://$rest_url/slurm/$openapi_release/jobs

#print("curl -X GET -H 'X-SLURM-USER-NAME:$user' -H 'X-SLURM-USER-TOKEN:$token' -H 'Content-Type: application/json' -H 'Accept: application/json' -d '{\"start_time\":\"1970-01-01\"}' http://$rest_hostname/slurm/$openapi_release/jobs");


# DEV
$ch = curl_init();
try {
	curl_setopt($ch, CURLOPT_URL, "http://$rest_hostname/slurmdb/$openapi_release/jobs?start_time=1970-01-01");
	#curl_setopt($ch, CURLOPT_URL, "http://$rest_hostname/slurmdb/$openapi_release/jobs");
	curl_setopt($ch, CURLOPT_HTTPHEADER, array(
		"X-SLURM-USER-NAME: $user",
		"X-SLURM-USER-TOKEN: $token",
		// "Content-Type: application/json",
		// "Accept: application/json"
	));

// 	-H "Accept: application/json" \
//   -d '{"param0":"pradeep"}' \
	// curl_setopt($ch, CURLOPT_HEADER, false);
	curl_setopt($ch, CURLOPT_RETURNTRANSFER, true);
	// curl_setopt($ch, CURLOPT_CONNECTTIMEOUT, 5);   
	// curl_setopt($ch, CURLOPT_TIMEOUT, 5);         
	// curl_setopt($ch, CURLOPT_FOLLOWLOCATION, true);
	// curl_setopt($ch, CURLOPT_MAXREDIRS, 2);
	// curl_setopt($ch, CURLOPT_NOBODY, true);
	// $post = [
	// 	'start_time' => '1970-01-01'
	// ];
	// curl_setopt($ch, CURLOPT_POSTFIELDS, $post);
	
	$response = curl_exec($ch);

	if (curl_errno($ch)) {
		echo curl_error($ch);
		die();
	}
	
	$http_code = curl_getinfo($ch, CURLINFO_HTTP_CODE);
	if ($http_code == intval(200)){
		// echo "Ressource valide";
	}
	else{
		echo "No connexion: " . $http_code;
	}
} catch (\Throwable $th) {
	throw $th;
} finally {
	curl_close($ch);
}


# SACCT
#$sacct_json = file_get_contents("/service/mgmtnode.sacct.json"); # old way
$sacct_json = $response;
if ($sacct_json === false) {
    // deal with error...
}

$sacct = json_decode($sacct_json, true);
if ($sacct === null) {
    // deal with error...
}

# JSON output
# curl -H "X-SLURM-USER-NAME:$user" -H "X-SLURM-USER-TOKEN:$token" http://proxy:8080/dashboard/?format=json
if ($OUTPUT_FORMAT == "json") {

	echo $sacct_json;
	exit;

};


# JOBS
$jobs=$sacct["jobs"];

krsort($jobs, $flags = SORT_REGULAR);


if ( count($jobs) > 0 ) {

	

	$nb = 0;
	$cell_style = "";

	foreach ($jobs as $job_id => $job_infos) {

		// if ( $nb > 2000 ) {
		// 	continue;
		// }

		# JOB ID
		$job_id = $job_infos['job_id'];

		# Time
		$seconds = $job_infos['time']['elapsed'];
		$time=$output = sprintf('%02d:%02d:%02d', ($seconds/ 3600),($seconds/ 60 % 60), $seconds% 60);

		# Requested
		$requested = "";
		foreach ($job_infos["tres"]["requested"] as $key => $value) {
			if ($value["type"] == "cpu") {
				$requested .= $value["count"]."CPUs ";
			};
			if ($value["type"] == "mem") {
				$requested .= round($value["count"]/1024)."Go ";
			};
		};

		# Allocated
		$allocated = "";
		foreach ($job_infos["tres"]["allocated"] as $key => $value) {
			if ($value["type"] == "cpu") {
				$allocated .= $value["count"]."CPUs ";
			};
			if ($value["type"] == "mem") {
				$allocated .= round($value["count"]/1024)."Go ";
			};
		};
		if ($allocated == "") {
			$allocated = "-";
		}

		// # Time Submission
		// if ($job_infos['time']['submission']) {
		// 	$submission = date("Y/m/d H:m:s",$job_infos['time']['submission']);
		// } else {
		// 	$submission = "";
		// };

		# Time Start
		if ($job_infos['time']['start']) {
			$start = date("Y/m/d H:m:s",$job_infos['time']['start']);
		} else {
			$start = "";
		};
		
		# Time End
		if ($job_infos['time']['end']) {
			$end = date("Y/m/d H:m:s",$job_infos['time']['end']);
		} else {
			$end = "";
		};

		# State
		$state = $job_infos['state']['current'];

		# Name
		$bad_characters = array("'");
		$good_characters = array("\'");
		$name = str_replace($bad_characters, $good_characters, $job_infos['name']);

		# Account
		$account = $job_infos['account'];

		# User
		$user_name = $job_infos['association']['user'];


		// print("<pre>");
		// print_r($job_infos['association']);
		// print("</pre>");

		# Row
		$dataObject_array[]="
				{
					id: '".$job_id."',
					name: '".$name."',
					state: '".$state."',
					
					user: '".$user_name."',
					account: '".$account."',
					qos: '".$job_infos['qos']."',
					priority: '".$job_infos['priority']."',
					time: '".$time."',
					start: '".$start."',
					end: '".$end."',
					allocated: '".$requested." / ".$allocated."'
				}

			";

		# Cell Style Renderer

		# ID renderer
		$cell_style .= "{row: ".$nb.", col: 0, renderer: 'GRAY_StylesRenderer'},";

		# State renderer
		$cell_style .= "{row: ".$nb.", col: 1, renderer: '".$state."_StylesRenderer'},";

		# Name renderer
		$cell_style .= "{row: ".$nb.", col: 2, renderer: 'BOLD_StylesRenderer'},";

		$nb++;

	};

	# Implode Row
	$dataObject_jobs=implode(",", $dataObject_array);
	$dataObject="
	var dataObject = [
		$dataObject_jobs
	];
	";

	# Columns
	$columns="[
		{ data: 'id', title: 'ID', type: 'numeric', readOnly: true, width: 50 },
		{ data: 'state', title: 'State', type: 'text', readOnly: true, width: 100, renderer: 'html' },
		{ data: 'name', title: 'Name', type: 'text', readOnly: true, width: 620, wordWrap: true },
		{ data: 'user', title: 'User', type: 'text', readOnly: true, width: 80 },
		{ data: 'account', title: 'Account', type: 'text', readOnly: true, width: 80 },
		{ data: 'qos', title: 'QOS', type: 'text', readOnly: true, width: 80 },
		{ data: 'priority', title: 'Priority', type: 'text', readOnly: true, width: 80 },
		{ data: 'time', title: 'Time', type: 'text', readOnly: true, width: 70 },
		{ data: 'start', title: 'Start', type: 'text', readOnly: true, width: 130 },
		{ data: 'end', title: 'End', type: 'text', readOnly: true, width: 130 },
		{ data: 'allocated', title: 'Requested / Allocated', type: 'text', readOnly: true }
	]";
	$colHeaders="[ 'State', 'ID', 'Name', 'User', 'Account', 'QOS', 'Priority', 'Time', 'Start', 'End', 'Resources' ]";
	$fixedColumnsLeft=2;


	?>
	<!DOCTYPE html>
	<html lang="en">

		<head>
			<title>List of Jobs</title>
			<link rel="stylesheet" type="text/css" href="handsontable/dist/handsontable.full.min.css">
			<link rel="stylesheet" type="text/css" href="handsontable/static/css/main.css">
			<script src="handsontable/dist/handsontable.full.min.js"></script>
		</head>

		<body>

			<div id="hot"></div>

			<script>

				<?php print $dataObject; ?>

				Handsontable.renderers.registerRenderer('COMPLETED_StylesRenderer', (hotInstance, TD, ...rest) => {
					Handsontable.renderers.TextRenderer(hotInstance, TD, ...rest);

					TD.style.fontWeight = 'bold';
					TD.style.color = 'white';
					TD.style.background = 'green';
				});

				Handsontable.renderers.registerRenderer('RUNNING_StylesRenderer', (hotInstance, TD, ...rest) => {
					Handsontable.renderers.TextRenderer(hotInstance, TD, ...rest);

					TD.style.fontWeight = 'bold';
					TD.style.color = 'white';
					TD.style.background = 'blue';
				});

				Handsontable.renderers.registerRenderer('PENDING_StylesRenderer', (hotInstance, TD, ...rest) => {
					Handsontable.renderers.TextRenderer(hotInstance, TD, ...rest);

					TD.style.fontWeight = 'bold';
					TD.style.color = 'white';
					TD.style.background = 'gray';
				});

				Handsontable.renderers.registerRenderer('FAILED_StylesRenderer', (hotInstance, TD, ...rest) => {
					Handsontable.renderers.TextRenderer(hotInstance, TD, ...rest);

					TD.style.fontWeight = 'bold';
					TD.style.color = 'white';
					TD.style.background = 'red';
				});

				Handsontable.renderers.registerRenderer('CANCELLED_StylesRenderer', (hotInstance, TD, ...rest) => {
					Handsontable.renderers.TextRenderer(hotInstance, TD, ...rest);

					TD.style.fontWeight = 'bold';
					TD.style.color = 'white';
					TD.style.background = 'red';
				});

				Handsontable.renderers.registerRenderer('BOLD_StylesRenderer', (hotInstance, TD, ...rest) => {
					Handsontable.renderers.TextRenderer(hotInstance, TD, ...rest);

					TD.style.fontWeight = 'bold';
					// TD.style.color = 'white';
					// TD.style.background = 'red';
				});

				Handsontable.renderers.registerRenderer('GRAY_StylesRenderer', (hotInstance, TD, ...rest) => {
					Handsontable.renderers.TextRenderer(hotInstance, TD, ...rest);

					TD.style.fontWeight = 'bold';
					TD.style.color = 'white';
					TD.style.background = 'gray';
				});


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
				fillHandle: false,
				//overflow: auto,
				//autoRowSize: {syncLimit: 10},
				//RowSize: 10,
				trimRows: false,
				columns: <?php print $columns ?>,
				stretchH: 'last',
				//width: '800',
				colWidths: 100,
				autoWrapCol: true,
				autoWrapRow: false,
				autoRowSize: false,
				height: '100%',
				//height: 800,
				// maxRows: 2000,
				rowHeaders: true,
				rowHeights: '1px',
				search: true,
				colHeaders: <?php echo $colHeaders; ?>,
				mergeCells: true,
				contextMenu: true,
				// columnSorting: {
				// 	indicator: true
				// },
				columnHeaders: true,
				// exportHiddenColumns: true,
				// exportHiddenRows: true,
				autoColumnSize: {
					samplingRatio: 23
				},
				manualRowResize: true,
				manualRowMove: true,
				manualColumnResize: true,
				manualColumnMove: true,
				manualColumnFreeze: true,
				filters: true,
				dropdownMenu: true,
				dragToScroll: true,
				fixedRowsTop: 0,
				fixedRowsBottom: 0,
				fixedColumnsLeft: <?php echo $fixedColumnsLeft; ?>,
				noWordWrapClassName: 'is-noWrapCell',
				undo: true,
				multiColumnSorting: {
					indicator: true
				},
				// hiddenColumns: {
				// 	columns: [5, 6, 7, 8, 9],
				// 	indicators: true
				// },
				// hiddenRows: {
				// 	// set rows that are hidden by default
				// 	rows: [0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13],
				// 	// when copying or pasting data, take hidden rows into account
				// 	// copyPasteEnabled: true,
				// 	// show where hidden rows are
				// 	indicators: true
				// },
				exportFile: true,
				// trimRows: [
				// 	100,
				// 	0
				// ],
				cell: [ <?php echo $cell_style; ?>],
				licenseKey: 'non-commercial-and-evaluation'
				};
				var hot = new Handsontable(hotElement, hotSettings);

			</script>

		</body>

	</html>

<?php
} else {
?>
<!DOCTYPE html>
<html lang="en">

	<head>
		<title>List of Jobs</title>
		<link rel="stylesheet" type="text/css" href="handsontable/dist/handsontable.full.min.css">
		<link rel="stylesheet" type="text/css" href="handsontable/static/css/main.css">
		<script src="handsontable/dist/handsontable.full.min.js"></script>
	</head>

	<body>
		<div style="text-align:center">
			<p>No jobs in SLURM</p>
		</div>
	</body>

</html>
<?php
};
?>