<?php
error_reporting(E_ALL);
session_start();

if (session_status() == PHP_SESSION_NONE || !isset($_SESSION['user_name'])) {
	header("HTTP/1.1 401 Unauthorized");
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

# SACCT
$sacct_json = file_get_contents("/service/mgmtnode.sacct.json");
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

// echo "<pre>";
foreach ($jobs as $job_id => $job_infos) {
    $job_list[$job_id]['id'] = $job_infos['job_id'];
    $job_list[$job_id]['name'] = $job_infos['name'];
    $job_list[$job_id]['time'] = $job_infos['time']['elapsed'];
	#echo $job_list[$job_id]['id']." ".$job_list[$job_id]['name']." ".$job_list[$job_id]['time']."<br>";
	#$job_infos["tres"];
	$seconds = $job_infos['time']['elapsed'];
	$time=$output = sprintf('%02d:%02d:%02d', ($seconds/ 3600),($seconds/ 60 % 60), $seconds% 60);
	#print_r($job_infos['steps'][0]["tres"]["allocated"]);
	$allocated = "";
	foreach ($job_infos['steps'][0]["tres"]["allocated"] as $key => $value) {
		if ($value["type"] == "cpu") {
			$allocated .= $value["count"]."CPUs ";
		};
		if ($value["type"] == "mem") {
			$allocated .= round($value["count"]/1024)."Go ";
		};
		#print($value["type"]."=".$value["count"]);
		// if ($value["count"] != "") {
		// 	$allocated .= $value["type"]."=".$value["count"]." ";
		// };
	};

	$status = $job_infos['exit_code']['status'];
	// if ( $job_infos['exit_code']['return_code'] > 0 ) {
	// 	$status .= "[".$job_infos['exit_code']['return_code']."]";
	// }

	$dataObject_array[]="
			{
				id: '".$job_infos['job_id']."',
				name: '".$job_infos['name']."',
				state: '".$job_infos['state']['current']."',
				status: '".$status."',
				
				account: '".$job_infos['account']."',
				qos: '".$job_infos['qos']."',
				priority: '".$job_infos['priority']."',
				time: '".$time."',
				allocated: '".$allocated."'
			}

		";
		# account: '".$job_infos['account']."',
		#alloc: '".$job_infos['required']['CPUs']." CPUs / ".round($job_infos['required']['memory']/1024)."Go'

};
$dataObject_jobs=implode(",", $dataObject_array);

$columns_fields_list="
	{ data: 'id', title: 'ID', type: 'numeric', readOnly: true, width: 10 },
	{ data: 'name', title: 'Name', type: 'text', readOnly: true, width: 100, wordWrap: true, noWordWrapClassName: 'is-noWrapCell' },
	{ data: 'state', title: 'State', type: 'text', readOnly: true, width: 20 },
	
	{ data: 'account', title: 'Account', type: 'text', readOnly: true, width: 20 },
	{ data: 'qos', title: 'Quality Of Service', type: 'text', readOnly: true, width: 30 },
	{ data: 'priority', title: 'Priority', type: 'text', readOnly: true, width: 15 },
	{ data: 'time', title: 'Time Elapsed', type: 'text', readOnly: true, width: 25 },
	{ data: 'allocated', title: 'Resources Allocated', type: 'text', readOnly: true }

	";
	# 	{ data: 'status', title: 'Status', type: 'text', readOnly: true, width: 30 },


	#{ data: 'account', title: 'Job Account', type: 'text', readOnly: true, width: 30 },
	#{ data: 'alloc', title: 'Job Resources Allocation', type: 'text', readOnly: true, width: 50 },
$colHeaders_fields_list=" 'ID', 'Name', 'State', 'Account', 'QOS', 'Priority', 'Time', 'Allocated' ";

$fixedColumnsLeft=2;

$dataObject="
var dataObject = [
 $dataObject_jobs
];
";

$columns="[
    $columns_fields_list
]";


$colHeaders="[
    $colHeaders_fields_list
]";

// print_r($dataObject_array);
// print_r($dataObject_jobs);
# print_r($dataObject);

// echo "</pre>";




// echo "<pre>";
// print_r($jobs);
// echo "</pre>";


?>
<!DOCTYPE html>
<html lang="en">

	<head>
		<title></title>

		<link rel="stylesheet" type="text/css" href="handsontable/dist/handsontable.full.min.css">
		<link rel="stylesheet" type="text/css" href="handsontable/static/css/main.css">

		<script src="handsontable/dist/handsontable.full.min.js"></script>
	</head>

	<body>

		<div id="hot"></div>

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
			fillHandle: true,
			//overflow: auto,
			//autoRowSize: {syncLimit: 10},
			//RowSize: 10,
			trimRows: true,
			columns: <?php print $columns ?>,
			stretchH: 'all',
			//width: '800',
			colWidths: 100,
			autoWrapCol: true,
			autoWrapRow: true,
			autoRowSize: true,
			height: '100%',
			//height: 800,
			maxRows: 2000,
			rowHeaders: true,
			rowHeights: '10px',
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
			exportFile: true,
			// trimRows: [
			// 	0,
			// 	0,
			// 	0,
			// 	0
			// ],
			licenseKey: 'non-commercial-and-evaluation'
			};
			var hot = new Handsontable(hotElement, hotSettings);
			// document.getElementById("export-tsv").addEventListener("click", function(event) { hot.getPlugin("exportFile").downloadFile("csv", {
			// 	filename: "<?php echo $t_filename; ?>_export_[YYYY][MM][DD]",
			// 	fileExtension: "tsv",
			// 	columnDelimiter: "\t",
			// 	exportHiddenColumns: true,
			// 	columnHeaders: true
			// });})
			// document.getElementById("export-string").addEventListener("click", function(event) {console.log(hot.getPlugin("exportFile").exportAsString("csv"));})


		</script>




		
		


	</body>

</html>
