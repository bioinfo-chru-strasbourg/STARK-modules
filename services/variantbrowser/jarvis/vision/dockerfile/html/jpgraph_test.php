<?php // content="text/plain; charset=utf-8"

$jpgraphlib="../jpgraph/src";

require_once ($jpgraphlib.'/jpgraph.php');
require_once ($jpgraphlib.'/jpgraph_line.php');
require_once ($jpgraphlib.'/jpgraph_bar.php');


#echo "BEGIN JPGRAPH test<BR>";


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

function readsunspotdata2($aFile, &$aYears, &$aSunspots) {
    $lines = @file($aFile,FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES);
    if( $lines === false ) {
        throw new JpGraphException('Can not read sunspot data file.');
    }
    foreach( $lines as $line => $datarow ) {
	if ($line!=0) {
        $split = preg_split('/[\s]+/',$datarow);
        #$split = preg_split('/[\t]+/',$datarow);
        #$aYears[] = substr(trim($split[0]),0,4);
        #$aYears[] = substr(trim($split[0]),0,30);
        #$aYears[] = str_replace(array(":","-","1"),"\n",substr(trim($split[0]),0,50));
        $aYears[] = $line;
        $aSunspots[] = trim($split[1]);
	};#if
    }
}


$dir="/media/IRC/V2/RES/ALL/140522_M01656_0009_000000000-D02GU/P974/P974.bwamem.bam.metrics";
$dir="/media/IRC/V2/RES/ALL/140627_M01656_0013_000000000-A947V/14C0001/14C0001.bwamem.bam.metrics";
#$file="yearssn.txt"; 
#$file="$dir/P974.bwamem.sample_interval_summary";
$file="$dir/14C0001.bwamem.sample_interval_summary";

$year = array();
$ydata = array();
readsunspotdata2($file,$year,$ydata);

#print_r($year);print_r($ydata);
$nb_values=count($ydata)-20;


// Just keep the last 20 values in the arrays
$year = array_slice($year, 0);
$ydata = array_slice($ydata, 0);

 // Width and height of the graph
$width = 1200; $height = 600;
 
// Create a graph instance
$graph = new Graph($width,$height);
 
// Specify what scale we want to use,
// text = txt scale for the X-axis
// int = integer scale for the Y-axis
$graph->SetScale('textint');
 
// Setup a title for the graph
$graph->title->Set('Coverage of sequenced regions');
 
// Setup titles and X-axis labels
$graph->xaxis->title->Set('(sequenced regions)');
$graph->xaxis->SetTickLabels($year);

// Setup Y-axis title
$graph->yaxis->title->Set('(coverage)');
 
// Create the bar plot
$barplot=new BarPlot($ydata);
   
// Add the plot to the graph
$graph->Add($barplot);
 
// Display the graph
$graph->Stroke();
 

#echo "END JPGRAPH test<BR>";


?>
