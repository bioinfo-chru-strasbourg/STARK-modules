<?php // content="text/plain; charset=utf-8"

$jpgraphlib="../jpgraph/src";

require_once ($jpgraphlib.'/jpgraph.php');
require_once ($jpgraphlib.'/jpgraph_line.php');
require_once ($jpgraphlib.'/jpgraph_bar.php');

function readdata($aFile, &$aYears, &$aSunspots, $chrs=array(), $regions=array(), $region_translate=array()) {
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
			$chrompos=trim($split[0]);
			$split_chrompos = split(':',$chrompos);
			$chrom=$split_chrompos[0];
			
			#$chrs=array("chr10");
			#print(join("",$chrs)." ".join("",$regions)."<BR>");
			if ( ( join("",$chrs)=="" || in_array($chrom,$chrs) )
				&& ( join("",$regions)=="" || in_array($chrompos,$regions) )
				) {

				$region_name=$chrompos;
				if ($region_translate[$chrompos]!="") {
					$region_translate_chrompos_split=split("[ _]",$region_translate[$chrompos]);
					$region_name=$region_translate_chrompos_split[0];
					#$aYears[] = str_replace(array("_"," "),"\n",$region_translate[$chrompos]); #."\ntest test t tset set setsrtsetest ";
					#$aYears[] = $region_translate_chrompos_split; #."\ntest test t tset set setsrtsetest ";
				} else {
					$region_name = $chrompos;
				};#if
				$aYears[] = $region_name;
				
				#$aYears[] = "A";
				$aSunspots[] = trim($split[2]);
				#$aSunspots[] = rand(0,10);
			};#if
		};#if
	};#foreach
}#function readdata

# Input parameters
$files=$_REQUEST["coverage_file"];
$legends=$_REQUEST["legends"];
$chrs=split(',',$_REQUEST["chrs"]);
#$_REQUEST["regions"]="chr1:43803495-43804416";
$regions=split(',',$_REQUEST["regions"]);
$title="Coverage of sequenced regions";
if (trim($_REQUEST["title"])!="") {
	$title=$_REQUEST["title"];
};#if
$region_translate_option_split=split(',',$_REQUEST["region_translate"]);
$region_translate=array();
foreach ($region_translate_option_split as $option) {
	if (trim($option)!="") {
		$option_split=split("\|",$option);
		$region_translate[$option_split[0]]=$option_split[1];
		#foreach ( as $chrompos=>$manifest_region_name) {
		#	print_r("$chrompos=>$manifest_region_name");
		#	$region_translate[$chrompos]=$manifest_region_name;
		#};#foreach
	};#if
};#foreach
#print_r($region_translate);

# Split input parameters
$files_split = preg_split('/,/',$files);
$legends_split = preg_split('/,/',$legends);

# Foreach coverage files
foreach ($files_split as $k=>$file) {
	if ($file!="") {
		$file=trim($file);
		
		# Read data
		$year = array();
		$ydata = array();
		readdata($file,$year,$ydata,$chrs,$regions,$region_translate);

		if (count($ydata)>0) {
			$nb_values=count($ydata)-0;

			$year = array_slice($year, 0);
			$ydata = array_slice($ydata, 0);
			#print $legends_split[$k]."<BR>";
			#print "<pre>"; print_r($year); print "</pre>";
			#print "<pre>"; print_r($ydata); print "</pre>";
			#foreach ($ydata as $v) {
			#	$ydata2[]=$v/2;
			#};#foreach

			$width = 750;
			$height = 550;
			$width = ($nb_values*40);
			if ($width<600) { $width=600; };#if
			if ($width>16000) { $width=16000; };#if
			#print " W:$width H:$height <BR>";
			# Create the bar plot
			$barplot=new BarPlot($ydata); #BarPlot
			$barplot->value->Show();
			$barplot->value->SetColor('black','darkred');

			# Legend
			$barplot->SetLegend($legends_split[$k]);
			$barplot->SetShadow();
			$GroupBar[]=$barplot;
		};#if

	};#if
};#foreach

# Header

echo "<B>$title</B><BR>";

$URI=$_SERVER["REQUEST_URI"];
$chr_links="Chromosomes: <A href='$URI&chrs='>ALL</A> ";
for ($chr=1; $chr<=22; $chr++) {
	$chr_links.=" <A href='$URI&chrs=chr$chr'>$chr</A>";
};#for
foreach (array("X","Y","M") as $chr) {
	$chr_links.=" <A href='$URI&chrs=chr$chr'>$chr</A>";
};#foreach
print $chr_links;
echo "<BR>";

if ($GroupBar!="") {

	$graph = new Graph($width,$height);
	$graph->SetMargin(150,50,50,50);
	// Specify what scale we want to use,
	// text = txt scale for the X-axis
	// int = integer scale for the Y-axis
	$graph->SetScale('textint'); #textint
	$graph->xaxis->SetLabelAngle(45);
	$graph->yaxis->scale->SetGrace(20);

	// Setup a title for the graph
	$graph->title->Set($title);
	 
	// Setup titles and X-axis labels
	$graph->xaxis->title->Set(''); #(sequenced regions)
	$graph->xaxis->SetTickLabels($year);

	// Setup Y-axis title
	$graph->yaxis->title->Set(''); #(coverage)

	$gbarplot = new GroupBarPlot($GroupBar);

	// Add the plot to the graph
	$graph->Add($gbarplot);
	$graph->legend->SetPos(0,0,'right','top');

	//File
	$coverage_file="tmp/coverage_file_".rand().".png";
	#print "$coverage_file";
	// Display the graph
	$graph->Stroke($coverage_file);
	#print "CHR<BR>";
	print "<IMG src='$coverage_file'>";

} else {

	print "<BR>No data...";

};#if

?>
