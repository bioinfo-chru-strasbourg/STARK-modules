<?php // content="text/plain; charset=utf-8"

$jpgraphlib="../jpgraph/src";

require_once ($jpgraphlib.'/jpgraph.php');
require_once ($jpgraphlib.'/jpgraph_line.php');
require_once ($jpgraphlib.'/jpgraph_bar.php');

function readdata($aFile, &$aYears, &$aSunspots) {
	# Read File
	$lines = @file($aFile,FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES);
	# If no data
	if( $lines === false ) {
		throw new JpGraphException('Can not read sunspot data file .');
	};#if
	# Read lines
	foreach( $lines as $line => $datarow ) {
		if ($line!=0) { # Skip first line
			$split = preg_split('/[\s]+/',$datarow);
			$aYears[] = trim($split[0]);
			$aSunspots[] = trim($split[1]);
		};#if
	};#foreach
}#function readdata

# Input parameters
$files=$_REQUEST["per_sequence_quality_scores_file"];
$legends=$_REQUEST["legends"];

# Split input parameters
$files_split = preg_split('/,/',$files);
$legends_split = preg_split('/,/',$legends);

#print_r($legends_split);

$width = 750;
$height = 550;


$graph = new Graph($width,$height);

# Foreach coverage files
foreach ($files_split as $k=>$file) {
	#print "$file";
	if ($file!="") {
		$file=trim($file);

		# Read data
		$year = array();
		$ydata = array();
		readdata($file,$year,$ydata);

		$nb_values=count($ydata)-20;

		$year = array_slice($year, 0);
		$ydata = array_slice($ydata, 0);
		foreach ($ydata as $v) {
			$ydata2[]=$v/2;
		};#foreach

		$width = 750;
		$height = 550;
		$width = ($nb_values*40);
		if ($width<600) { $width=600; };#if
		#print_r($ydata);
		# Create the bar plot
		#$barplot=new BarPlot($ydata); #BarPlot
		$lineplot=new LinePlot($ydata); #BarPlot
		# Legend
		$lineplot->SetLegend($legends_split[$k]);
		#$barplot->SetShadow();
		$lineplot->SetWeight( 2 );   // Two pixel wide
		// Add the plot to the graph
		$graph->Add($lineplot);
		#$GroupBar[]=$lineplot;


		# TEST
		/*$lineplot=new LinePlot($ydata2); #BarPlot
		$lineplot->SetLegend($legends_split[$k]);
		#$barplot->SetShadow();
		$lineplot->SetWeight( 2 );   // Two pixel wide
		// Add the plot to the graph
		$graph->Add($lineplot);
		*/
	};#if
};#foreach


$graph->SetMargin(150,50,50,50);
// Specify what scale we want to use,
// text = txt scale for the X-axis
// int = integer scale for the Y-axis
#$graph->SetScale('textint'); #textint intlin
$graph->SetScale('intlin'); #textint intlin
$graph->xaxis->SetLabelAngle(45);
$graph->yaxis->scale->SetGrace(20);

// Setup a title for the graph
$graph->title->Set('Per sequence quality scores');
 
// Setup titles and X-axis labels
$graph->xaxis->title->Set(''); #(sequenced regions)
$graph->xaxis->SetTickLabels($year);

// Setup Y-axis title
$graph->yaxis->title->Set(''); #(coverage)

#$gbarplot = new GroupBarPlot($GroupBar);
#$gbarplot = new GroupLinePlot($GroupBar);

// Add the plot to the graph
#$graph->Add($gbarplot);
#$graph->legend->SetPos(0,0,'right','top');

//File
$coverage_file="tmp/coverage_file_".rand().".png";

// Display the graph
$graph->Stroke($coverage_file);
print "<IMG src='$coverage_file'>";


?>
