<?php

// TSV to Array Function
// Copyright (c) 2015, Ink Plant
// https://inkplant.com/code/
// this version was last updated June 23, 2015

function tsv_to_array($file,$args=array()) {

	//key => default
	$fields = array(
		'header_row'=>true,
		'remove_header_row'=>true,
		'trim_headers'=>true, //trim whitespace around header row values
		'trim_values'=>true, //trim whitespace around all non-header row values
		'debug'=>false, //set to true while testing if you run into troubles
		'lb'=>"\n", //line break character
		'tab'=>"\t", //tab character
	);
	foreach ($fields as $key => $default) {
		if (array_key_exists($key,$args)) { $$key = $args[$key]; }
		else { $$key = $default; }
	}

	if (!file_exists($file)) {
		if ($debug) { $error = 'File does not exist: '.htmlspecialchars($file).'.'; }
		else { $error = 'File does not exist.'; }
		#custom_die($error);
		die($error);
	}

    if ($debug) { echo '<p>Opening '.htmlspecialchars($file).'&hellip;</p>'; }
    $data = array();

    if (($handle = fopen($file,'r')) !== false) {
		$contents = fread($handle, filesize($file));
		fclose($handle);
    } else {
        custom_die('There was an error opening the file.');
    }

	$lines = explode($lb,$contents);
	if ($debug) { echo '<p>Reading '.count($lines).' lines&hellip;</p>'; }

	$row = 0;
	foreach ($lines as $line) {
		$row++;
		if (($header_row) && ($row == 1)) { $data['headers'] = array(); }
		else { $data[$row] = array(); }
		$values = explode($tab,$line);
		foreach ($values as $c => $value) {
			if (($header_row) && ($row == 1)) { //if this is part of the header row
				if (in_array($value,$data['headers'])) { custom_die('There are duplicate values in the header row: '.htmlspecialchars($value).'.'); }
				else {
					if ($trim_headers) { $value = trim($value); }
					$data['headers'][$c] = $value.''; //the .'' makes sure it's a string
				}
			} elseif ($header_row) { //if this isn't part of the header row, but there is a header row
				$key = $data['headers'][$c];
				if ($trim_values) { $value = trim($value); }
				$data[$row][$key] = $value;
			} else { //if there's not a header row at all
				$data[$row][$c] = $value;
			}
		}
	}

	if ($remove_header_row) {
		unset($data['headers']);
	}

	if ($debug) { echo '<pre>'.print_r($data,true).'</pre>'; }
	return $data;
}

?>