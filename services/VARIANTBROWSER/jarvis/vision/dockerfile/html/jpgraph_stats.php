<?php // content="text/plain; charset=utf-8"

$jpgraphlib="../jpgraph/src";

require_once ($jpgraphlib.'/jpgraph.php');
require_once ($jpgraphlib.'/jpgraph_line.php');
require_once ($jpgraphlib.'/jpgraph_bar.php');
require_once ($jpgraphlib.'/jpgraph_pie.php'); 
require_once ($jpgraphlib.'/jpgraph_pie3d.php'); 
 

// Some data
#$data = array(50,28,25,27,31,20);
$input_title=$_REQUEST["title"];
$title=$input_title;
$input_subtitle=$_REQUEST["subtitle"];
$subtitle=str_replace("|","\n",$input_subtitle);
$input_data=$_REQUEST["data"];
$explode_data=explode(",",$input_data);
$input_label=$_REQUEST["label"];
$label=explode(",",$input_label);
#$lbl = array("adam\n%.1f%%","bertil\n%.1f%%","johan\n%.1f%%",
#         "peter\n%.1f%%","daniel\n%.1f%%","erik\n%.1f%%");
#print $lbl;
foreach ($label as $key=>$value) {
	if (trim($explode_data[$key])!="") {
		#$lbl[]= "$value\n%.2f%%";
		$lbl[]= "$value\n%.2f%%";
		$data[]=$explode_data[$key];
	};#if
};#foreach
#print $lbl;

// A new pie graph
$graph = new PieGraph(400,400,'auto');
 
// Don't display the border
$graph->SetFrame(false);
 
// Uncomment this line to add a drop shadow to the border
// $graph->SetShadow();
 
// Setup title
$graph->title->Set($title);
#$graph->title->SetFont(FF_ARIAL,FS_BOLD,18);
$graph->title->SetMargin(8); // Add a little bit more margin from the top
 
// Create the pie plot
$p1 = new PiePlotC($data);
 
// Set size of pie
$p1->SetSize(0.40);
 
// Label font and color setup
#$p1->value->SetFont(FF_ARIAL,FS_BOLD,12);
$p1->value->SetColor('white');
 
$p1->value->Show();
 
// Setup the title on the center circle
$p1->midtitle->Set($subtitle);
#$p1->midtitle->SetFont(FF_ARIAL,FS_NORMAL,14);
 
// Set color for mid circle
$p1->SetMidColor('white');
 
// Use percentage values in the legends values (This is also the default)
$p1->SetLabelType(PIE_VALUE_PER);
 
// The label array values may have printf() formatting in them. The argument to the
// form,at string will be the value of the slice (either the percetage or absolute
// depending on what was specified in the SetLabelType() above.
#$lbl = array("adam\n%.1f%%","bertil\n%.1f%%","johan\n%.1f%%",
#         "peter\n%.1f%%","daniel\n%.1f%%","erik\n%.1f%%");
$p1->SetLabels($lbl);
 
// Uncomment this line to remove the borders around the slices
// $p1->ShowBorder(false);
 
// Add drop shadow to slices
$p1->SetShadow();
 
// Explode all slices 15 pixels
$p1->ExplodeAll(10);

// Theme
$p1->SetTheme("earth");
 
// Add plot to pie graph
$graph->Add($p1);

 
// .. and send the image on it's marry way to the browser
$graph->Stroke();











die;

function readsunspotdata($aFile, &$aYears, &$aSunspots) {
    $lines = @file($aFile,FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES);
    if( $lines === false ) {
        throw new JpGraphException('Can not read sunspot data file.');
    }
    foreach( $lines as $line => $datarow ) {
        $split = preg_split('/[\s]+/',$datarow);
        $aYears[] = substr(trim($split[0]),0,4);
        $aSunspots[] = trim($split[1]);
    }
}




#Test
if (0) {
$datay=array(2,3,5,8,12,6,3);
$datax=array('Jan','Feb','Mar','Apr','May','Jun','Jul');
 
// Size of graph
$width=400;
$height=500;
 
// Set the basic parameters of the graph
$graph = new Graph($width,$height,'auto');
$graph->SetScale('textlin');
 
// Rotate graph 90 degrees and set margin
$graph->Set90AndMargin(50,20,50,30);
 
// Nice shadow
$graph->SetShadow();
 
// Setup title
$graph->title->Set('Horizontal bar graph ex 1');
#$graph->title->SetFont(FF_VERDANA,FS_BOLD,14);
 
// Setup X-axis
$graph->xaxis->SetTickLabels($datax);
#$graph->xaxis->SetFont(FF_VERDANA,FS_NORMAL,12);
 
// Some extra margin looks nicer
$graph->xaxis->SetLabelMargin(10);
 
// Label align for X-axis
$graph->xaxis->SetLabelAlign('right','center');
 
// Add some grace to y-axis so the bars doesn't go
// all the way to the end of the plot area
$graph->yaxis->scale->SetGrace(20);
 
// We don't want to display Y-axis
$graph->yaxis->Hide();
 
// Now create a bar pot
$bplot = new BarPlot($datay);
$bplot->SetFillColor('orange');
$bplot->SetShadow();
 
//You can change the width of the bars if you like
//$bplot->SetWidth(0.5);
 
// We want to display the value of each bar at the top
$bplot->value->Show();
$bplot->value->SetFont(FF_ARIAL,FS_BOLD,12);
$bplot->value->SetAlign('left','center');
$bplot->value->SetColor('black','darkred');
$bplot->value->SetFormat('%.1f mkr');
 
// Add the bar to the graph
$graph->Add($bplot);
 
// .. and stroke the graph
$graph->Stroke();
die();

};#if



$dir="/media/IRC/V2/RES/ALL/140522_M01656_0009_000000000-D02GU/P974/P974.bwamem.bam.metrics";
$dir="/media/IRC/V2/RES/ALL/140627_M01656_0013_000000000-A947V/14C0001/14C0001.bwamem.bam.metrics";
#$file="yearssn.txt"; 
#$file="$dir/P974.bwamem.sample_interval_summary";
$file="$dir/14C0001.bwamem.sample_interval_summary";

$files=$_REQUEST["coverage_file"];
$legends=$_REQUEST["legends"];

$files_split = preg_split('/,/',$files);
$legends_split = preg_split('/,/',$legends);

foreach ($files_split as $k=>$file) {
	if ($file!="") {
		$file=trim($file);
		#print("$file<BR>");

		$year = array();
		$ydata = array();
		readsunspotdata2($file,$year,$ydata);

		#print_r($year);print_r($ydata);
		$nb_values=count($ydata)-20;


		// Just keep the last 20 values in the arrays
		$year = array_slice($year, 0);
		$ydata = array_slice($ydata, 0);
		foreach ($ydata as $v) {
			$ydata2[]=$v/2;
		};#foreach

		 // Width and height of the graph
		$width = 800;
		$height = 600;
		#$height = ($nb_values*30);
		#if ($height<800) { $height=800; };#if
		$width = ($nb_values*40);
		if ($width<600) { $width=600; };#if
		 
		// Create a graph instance

		
		
		
		// Create the bar plot
		$barplot=new BarPlot($ydata); #BarPlot
		$barplot->value->Show();
		$barplot->value->SetColor('black','darkred');

		#$barplot2->SetFillColor('orange');
		$barplot->SetLegend($legends_split[$k]);
		#$barplot->SetColor("blue");
		#$barplot->SetFillColor('orange');
		$barplot->SetShadow();
		#$barplot->SetWeight(0);
		#$barplot->SetFillGradient('red','yellow');
		#$barplot->legend->SetPos(0.5,0.98,'center','bottom');
		$GroupBar[]=$barplot;
		#$GroupBar[]=$barplot2;
	};#if
};#foreach

$graph = new Graph($width,$height);
$graph->SetMargin(150,50,50,50);
#$graph->Set90AndMargin(0,0,0,0);
// Specify what scale we want to use,
// text = txt scale for the X-axis
// int = integer scale for the Y-axis
$graph->SetScale('textint'); #textint
$graph->xaxis->SetLabelAngle(45);
$graph->yaxis->scale->SetGrace(20);

// Setup a title for the graph
$graph->title->Set('Coverage of sequenced regions');
 
// Setup titles and X-axis labels
#$graph->xaxis->SetFont(FS_VERDANA,FS_NORMAL,10);
#$graph->yaxis->SetFont(FF_VERDANA,FS_NORMAL,10);
$graph->xaxis->title->Set(''); #(sequenced regions)
$graph->xaxis->SetTickLabels($year);


// Setup Y-axis title
$graph->yaxis->title->Set(''); #(coverage)

$gbarplot = new GroupBarPlot($GroupBar);

// Add the plot to the graph
#$graph->Add($barplot);
#$graph->Add($barplot2);
$graph->Add($gbarplot);
$graph->legend->SetPos(0,0,'right','top');


// Display the graph
$graph->Stroke();


#echo "END JPGRAPH test<BR>";


?>
