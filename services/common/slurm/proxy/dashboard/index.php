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


# Connexion params
$user=$_SESSION['user_name'];

if ($user == "") {
	$user = "root";
};

#$token=$_SESSION['user_token']; #$token="eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJleHAiOjE2ODAxNjY2ODksImlhdCI6MTY3MDE2NjY5MCwic3VuIjoic2x1cm0ifQ.r5LV-z4EWjvHDJUgExQS1EHdMx-nPXvNowGmXFuSAJg";
$token=trim(file_get_contents("/auth/slurm")); #$rest_url="rest";
$rest_hostname=trim(file_get_contents("/auth/rest_hostname")); #$rest_url="rest";
$openapi_release=trim(file_get_contents("/auth/openapi_release")); #$openapi_release="v0.0.38";
// # curl -H "X-SLURM-USER-NAME:$user" -H "X-SLURM-USER-TOKEN:$token" http://$rest_url/slurmdb/$openapi_release/qos
// # curl -H "X-SLURM-USER-NAME:$user" -H "X-SLURM-USER-TOKEN:$token" http://$rest_url/slurm/$openapi_release/jobs

#print("curl -X GET -H 'X-SLURM-USER-NAME:$user' -H 'X-SLURM-USER-TOKEN:$token' -H 'Content-Type: application/json' -H 'Accept: application/json' -d '{\"start_time\":\"1970-01-01\"}' http://$rest_hostname/slurm/$openapi_release/jobs");


# DEV
$ch = curl_init();
try {
	curl_setopt($ch, CURLOPT_URL, "http://$rest_hostname/slurmdb/$openapi_release/jobs?start_time=1970-01-01");
	curl_setopt($ch, CURLOPT_HTTPHEADER, array(
		"X-SLURM-USER-NAME: $user",
		"X-SLURM-USER-TOKEN: $token",
	));

	curl_setopt($ch, CURLOPT_RETURNTRANSFER, true);
	
	$response = curl_exec($ch);

	if (curl_errno($ch)) {
		echo curl_error($ch);
		die();
	}
	
	$http_code = curl_getinfo($ch, CURLINFO_HTTP_CODE);
	if ($http_code == intval(200)){
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
$sacct_json = $response;
if ($sacct_json === false) {
    // deal with error...
}

$sacct = json_decode($sacct_json, true);
if ($sacct === null) {
    // deal with error...
}

# JSON output
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
	$hiddenRows = [];

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
		$requested_cpu = "";
		$requested_mem = "";
		foreach ($job_infos["tres"]["requested"] as $key => $value) {
			if ($value["type"] == "cpu") {
				$requested .= $value["count"]."CPUs ";
				$requested_cpu = $value["count"];
			};
			if ($value["type"] == "mem") {
				$requested .= round($value["count"]/1024)."Go ";
				$requested_mem = round($value["count"]/1024)."";
			};
		};

		# Allocated
		$allocated = "";
		$allocated_cpu = "";
		$allocated_mem = "";
		foreach ($job_infos["tres"]["allocated"] as $key => $value) {
			if ($value["type"] == "cpu") {
				$allocated .= $value["count"]."CPUs ";
				$allocated_cpu = $value["count"];
			};
			if ($value["type"] == "mem") {
				$allocated .= round($value["count"]/1024)."Go ";
				$allocated_mem = round($value["count"]/1024)."";
			};
		};
		if ($allocated == "") {
			$allocated = "-";
		}

		# Time Submission
		if ($job_infos['time']['submission']) {
			$submission = date("Y/m/d H:m:s",$job_infos['time']['submission']);
		} else {
			$submission = "";
		};

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

		if ($start == "" && $end == "") {
			$time = "";
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

		
		$dataObject_table_tbody[]="
			<tr>
				<td></td>
				<td>".$state."</td>
				<td>".$name."</td>
				<td>".$job_id."</td>
				<td>".$user_name."</td>
				<td>".$account."</td>
				<td>".$job_infos['qos']."</td>
				<td>".$job_infos['priority']."</td>
				<td>".$time."</td>
				<td>".$submission."</td>
				<td>".$start."</td>
				<td>".$end."</td>
				<td>".$requested_cpu."</td>
				<td>".$requested_mem."</td>
				<td>".$allocated_cpu."</td>
				<td>".$allocated_mem."</td>
				

			</tr>
		";

		$dataObject_table_theader="
			<tr>
                <th colspan='3' rowspan='1'></th>
                <th colspan='3'>Launch Infos</th>
                <th colspan='2'>Service Quality</th>
                <th colspan='4'>Execution Time</th>
                <th colspan='2'>Requested</th>
                <th colspan='2'>Allocated</th>
            </tr>
			<tr>
				<th></th>
				<th>State</th>
				<th>Job name</th>
				<th>ID</th>
				<th>User</th>
				<th>Account</th>
				<th>QOS</th>
				<th>Priority</th>
				<th>Time</th>
				<th>Submission</th>
				<th>Start</th>
				<th>End</th>
				<th>CPU</th>
				<th>MEM</th>
				<th>CPU</th>
				<th>MEM</th>

			</tr>
		";

		$dataObject_table_tfooter="
			<tr>
				<th></th>
				<th>State</th>
				<th>Job name</th>
				<th>ID</th>
				<th>User</th>
				<th>Account</th>
				<th>QOS</th>
				<th>Priority</th>
				<th>Time</th>
				<th>Submission</th>
				<th>Start</th>
				<th>End</th>
				<th>CPU</th>
				<th>MEM</th>
				<th>CPU</th>
				<th>MEM</th>

			</tr>
		";


		$nb++;

	};

	$dataObject_table_jobs=implode("", $dataObject_table_tbody);

	?>
	<!DOCTYPE html>
	<html lang="en">

		<head>
			<title>SLURM Jobs' list</title>

			<!-- HandsOnTable -->
			<link rel="stylesheet" type="text/css" href="handsontable/dist/handsontable.full.min.css">
			<link rel="stylesheet" type="text/css" href="handsontable/static/css/main.css">
			<script src="handsontable/dist/handsontable.full.min.js"></script>

			<!-- DataTables -->
			<link rel="stylesheet" type="text/css" href="DataTables/datatables.min.css"/>
			<script type="text/javascript" src="DataTables/datatables.min.js"></script>


		</head>

		<body>

		<!-- <div>
        	Toggle column: <a class="toggle-vis" data-column="3">ID</a>
			- <a class="toggle-vis" data-column="4">User</a>
			- <a class="toggle-vis" data-column="5">Account</a>
			- <a class="toggle-vis" data-column="6">QOS</a>
			- <a class="toggle-vis" data-column="7">Priority</a>
			- <a class="toggle-vis" data-column="8">Time</a>
			- <a class="toggle-vis" data-column="9">Submission</a>
			- <a class="toggle-vis" data-column="10">Start</a>
			- <a class="toggle-vis" data-column="11">End</a>
			- <a class="toggle-vis" data-column="12">CPU</a>
			- <a class="toggle-vis" data-column="13">MEM</a>
    	</div> -->

			<table id="jobs" class="display nowrap table table-striped" style="width:100%">
				<thead>
					<?php print($dataObject_table_theader); ?>
				</thead>
				<tbody>
					<?php print($dataObject_table_jobs); ?>
				</tbody>
				<tfoot>
					<?php //print($dataObject_table_tfooter); ?>
				</tfoot>
			</table>

			<script>

			/* Formatting function for row details - modify as you need */
			function format(d) {
				// `d` is the original data object for the row
				var keys = Object.keys(d);

				// State,ID,Name,User,Account,QOS,Priority,Time,Start,End,Submission,CPU,MEM

				return (
					'<table cellpadding="5" cellspacing="0" border="0" style="padding-left:50px;">' +

					'<tr>' + '<td colspan="2">Job Name: </td><td width="20px"></td><td><b>' + d.Name + '</td></tr>' +
					'<tr>' + '<td colspan="2">State: </td><td width="20px"></td><td>' + d.State + '</td></tr>' +

					'<tr>' + '<td colspan="4">Launch Infos</td></tr>' +
					'<tr>' + '<td width="20px"></td><td>Job ID: </td><td width="20px"></td><td>' + d.ID + '</td></tr>' +
					'<tr>' + '<td width="20px"></td><td>User: </td><td width="20px"></td><td>' + d.User + '</td></tr>' +
					'<tr>' + '<td width="20px"></td><td>Account: </td><td width="20px"></td><td>' + d.Account + '</td></tr>' +

					'<tr>' + '<td colspan="4">Service Quality</td></tr>' +
					'<tr>' + '<td width="20px"><td>Quality Of Service: </td><td width="20px"></td><td>' + d.QOS + '</td></tr>' +
					'<tr>' + '<td width="20px"><td>Priority: </td><td width="20px"></td><td>' + d.Priority + '</td></tr>' +

					'<tr>' + '<td colspan="4">Execution Time</td></tr>' +
					'<tr>' + '<td width="20px"><td>Time: </td><td width="20px"></td><td>' + d.Time + '</td></tr>' +
					'<tr>' + '<td width="20px"><td>Submission: </td><td width="20px"></td><td>' + d.Submission + '</td></tr>' +
					'<tr>' + '<td width="20px"><td>Start: </td><td width="20px"></td><td>' + d.Start + '</td></tr>' +
					'<tr>' + '<td width="20px"><td>End: </td><td width="20px"></td><td>' + d.End + '</td></tr>' +
					
					'<tr>' + '<td colspan="4">Resources</td></tr>' +
					'<tr>' + '<td width="20px"><td>Requested CPU: </td><td width="20px"></td><td>' + d.requestedCPU + '</td></tr>' +
					'<tr>' + '<td width="20px"><td>Requested Memory (Go): </td><td width="20px"></td><td>' + d.requestedMEM + '</td></tr>' +
					'<tr>' + '<td width="20px"><td>Allocated CPU: </td><td width="20px"></td><td>' + d.allocatedCPU + '</td></tr>' +
					'<tr>' + '<td width="20px"><td>Allocated Memory (Go): </td><td width="20px"></td><td>' + d.allocatedMEM + '</td></tr>' +


					'</table>'
				);
			}


			$(document).ready(function () {

				// Setup - add a text input to each footer cell
				// $('#jobs tfoot th').each(function () {
				// 	var title = $(this).text();
				// 	$(this).html('<input type="text" placeholder="Search ' + title + '" />');
				// });

				// $('#jobs tbody').on('mouseenter', 'td', function () {
				// 	var colIdx = table.cell(this).index().column;
			
				// 	$(table.cells().nodes()).removeClass('highlight');
				// 	$(table.column(colIdx).nodes()).addClass('highlight');
				// });

				// Add event listener for opening and closing details
				$('#jobs tbody').on('click', 'td.dt-control', function () {
					var tr = $(this).closest('tr');
					var row = table.row(tr);
			
					if (row.child.isShown()) {
						// This row is already open - close it
						row.child.hide();
						// tr.removeClass('shown');
					} else {
						// Open this row
						row.child(format(row.data())).show();
						// tr.addClass('shown');
					}
				});

				$('#jobs').on('requestChild.dt', function (e, row) {
					row.child(format(row.data())).show();
				});

				$('a.toggle-vis').on('click', function (e) {
					e.preventDefault();
			
					// Get the column API object
					var column = table.column($(this).attr('data-column'));
			
					// Toggle the visibility
					column.visible(!column.visible());
				});

				// DataTable
				var table = $('#jobs').DataTable({
					paging: true,
					pagingType: 'full_numbers',
					order: [[3, 'desc']],
					lengthMenu: [
						[10, 25, 50, -1],
						[10, 25, 50, 'All'],
					],
					rowId: 'id',
					stateSave: true,
					// columnDefs: [
					// 	{
					// 		target: ,
					// 		visible: true,
					// 		searchable: true,
					// 	},
					// 	{
					// 		target: 3,
					// 		visible: true,
					// 	},
					// ],
					columns: [
						{
							className: 'dt-control',
							orderable: false,
							data: null,
							defaultContent: '',
						},
						// { data: '' },
						{
							data: 'State',
							render: function (data, type) {
								
								if (type === 'display') {
									let color = 'gray';
									if (data == "COMPLETED") {
										color = 'green';
									} else if (data == "RUNNING") {
										color = 'blue';
									} else if (data == "PENDING") {
										color = 'gray';
									} else if (data == "FAILED") {
										color = 'red';
									} else if (data == "CANCELLED") {
										color = 'red';
									}
			
									return '<span style="font-weight:bold;color:' + color + '">' + data + '</span>';
								}
			
								return data;
							},
						},
						{ data: 'Name' },
						{ data: 'ID' },
						{ data: 'User' },
						{ data: 'Account' },
						{ data: 'QOS' },
						{ data: 'Priority' },
						{ data: 'Time' },
						{ data: 'Submission' },
						{ data: 'Start' },
						{ data: 'End' },
						{ 
							data: 'requestedCPU',
							render: function (data, type) {
								
								if (type === 'display') {
									$unit = "";
									if (data != "") {
										$unit = " CPU";
									}
									return data + $unit;
								}
			
								return data;
							},
						},
						{ 
							data: 'requestedMEM',
							render: function (data, type) {
								
								if (type === 'display') {
									$unit = "";
									if (data != "") {
										$unit = " Go";
									}
									return data + $unit;
								}
			
								return data;
							},
						},
						{ 
							data: 'allocatedCPU',
							render: function (data, type) {
								
								if (type === 'display') {
									$unit = "";
									if (data != "") {
										$unit = " CPU";
									}
									return data + $unit;
								}
			
								return data;
							},
						},
						{ 
							data: 'allocatedMEM',
							render: function (data, type) {
								
								if (type === 'display') {
									$unit = "";
									if (data != "") {
										$unit = " Go";
									}
									return data + $unit;
								}
			
								return data;
							},
						},
					],
					// initComplete: function () {
					// 	// Apply the search
					// 	this.api()
					// 		.columns()
					// 		.every(function () {
					// 			var that = this;
			
					// 			$('input', this.footer()).on('keyup change clear', function () {
					// 				if (that.search() !== this.value) {
					// 					that.search(this.value).draw();
					// 				}
					// 			});
					// 		});
					// },
				});
			});
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