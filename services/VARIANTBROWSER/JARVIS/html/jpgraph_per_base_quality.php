<?php // content="text/plain; charset=utf-8"

$jpgraphlib="../jpgraph/src";

require_once ($jpgraphlib.'/jpgraph.php');
require_once ($jpgraphlib.'/jpgraph_line.php');
require_once ($jpgraphlib.'/jpgraph_bar.php');
require_once ($jpgraphlib.'/jpgraph_stock.php');

function readdata($aFile, &$aYears, &$aSunspots, &$mean, &$median, &$top, &$middle, &$bottom) {
	# Read File
	$lines = @file($aFile,FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES);
	#print_r($lines);
	# If no data
	if( $lines === false ) {
		throw new JpGraphException('Can not read sunspot data file .');
	};#if
	# Read lines
	foreach( $lines as $line => $datarow ) {
		if ($line!=0) { # Skip first line
			$split = preg_split('/[\s]+/',$datarow);
			$aYears[] = trim($split[0]);
			$top[] = 40;
			$middle[] = 30;
			$bottom[] = 20;
			$mean[] = trim($split[1]);
			$median[] = trim($split[2]);
			$aSunspots[] = trim($split[3]);
			$aSunspots[] = trim($split[4]);
			$aSunspots[] = trim($split[5]);
			$aSunspots[] = trim($split[6]);
		};#if
	};#foreach
}#function readdata

# Input parameters
$files=$_REQUEST["per_base_quality_file"];
$legends=$_REQUEST["legends"];

# Split input parameters
$files_split = preg_split('/,/',$files);
$legends_split = preg_split('/,/',$legends);

#print_r($files_split);
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
		$year = array("");
		$top = array(40);
		$middle = array(30);
		$bottom = array(20);
		$mean = array("");
		$median = array("");
		$ydata = array(0,0,0,0);
		
		readdata($file,$year,$ydata, $mean, $median, $top, $middle, $bottom);
		$top[]="40"; $middle[]="30"; $bottom[]="20";
		$mean[]=""; $median[]="";
		$ydata[] = "";$ydata[] = "";$ydata[] = "";$ydata[] = "";
		$year[]="";
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
		
		#print_r($mean);
		
		
		# Surface
		# top
		$lineplot=new LinePlot($top); #BarPlot
		$lineplot->SetFillColor("lightgreen");
		$graph->Add($lineplot);
		# middle
		$lineplot=new LinePlot($middle); #BarPlot
		$lineplot->SetFillColor("lightyellow");
		$graph->Add($lineplot);
		# bottom
		$lineplot=new LinePlot($bottom); #BarPlot
		$lineplot->SetFillColor("lightred");
		$graph->Add($lineplot);
		

		# Mean
		$lineplot=new LinePlot($mean); #BarPlot
		#$lineplot->SetFillColor("orange");
		$graph->Add($lineplot);
		# Median
		$lineplot=new LinePlot($median); #BarPlot
		#$lineplot->SetFillColor("orange");
		$graph->Add($lineplot);

		# Create the bar plot
		#$barplot=new BarPlot($ydata); #BarPlot
		$StockPlot=new StockPlot($ydata); #BarPlot
		# Legend
		$StockPlot->SetLegend($legends_split[$k]);
		#$lineplot->SetShadow();
		$StockPlot->SetWeight( 2 );   // Two pixel wide
		// Add the plot to the graph
		$graph->Add($StockPlot);
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
#die();

$graph->SetMargin(150,50,50,50);
// Specify what scale we want to use,
// text = txt scale for the X-axis
// int = integer scale for the Y-axis
#$graph->SetScale('textint'); #textint intlin
$graph->SetScale('textlin'); #textint intlin
$graph->xaxis->SetLabelAngle(45);
$graph->yaxis->scale->SetGrace(20);
$graph->SetMarginColor('lightblue');

// Setup a title for the graph
$graph->title->Set('Per base quality');
#$graph->xaxis->SetPos("max");
// Setup titles and X-axis labels
$graph->yscale->SetGrace(0);
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
$coverage_file="tmp/per_base_quality_file_".rand().".png";

// Display the graph
$graph->Stroke($coverage_file);
print "<IMG src='$coverage_file'>";


?>
